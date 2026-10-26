# %% v38（事後・amendment B-2）: 浸水図はどのセンサーで作られたか — レーダーのみ・光学のみの AOI で名前の増分を見る
"""ラベルの誤差は撮像の種類で違う（レーダーは市街地・植生下の水を取りこぼしやすく、光学は雲で欠ける）。
主分割の各 AOI について、EMS 製品の source 表（第 1–3 世代の全製品にある: source_nam・eventphase・
sensor_gsd / sensor_res）から洪水後の画像のセンサー名を読み、レーダーのみ・光学のみ・混在・不明に分ける。
そのうえで、確証解析（world_v2、主分割）の held-out 予測を AOI の群ごとに集計し直し、名前の場所内の増分
ΔBS(M1c → M3n) と R（群内の ΔBS(M0 → M1) に対する比）を記述する。再学習はしない。主判定は変えない。

2 段階:
- ``--extract``（クラスタ、NAS 上の製品 zip を読む）: data/verification/v38_ems_sensors_raw.json
- 引数なし（ローカル、上を引き上げた後）: 分類・集計 → stdout の表と docs/verification/v38_ems_sensors.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import zipfile
from pathlib import Path

import numpy as np

# コンテナ（/workspace）には CLAUDE.md が無いので、見つからなければ backend の 1 つ上をルートとみなす
ROOT = next((p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file()),
            Path(__file__).resolve().parents[3])
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR, RUNS_DIR

RAW = "verification/v38_ems_sensors_raw.json"  # DATA_DIR からの相対
OUT = ROOT / "docs" / "verification" / "v38_ems_sensors.json"
B = 4000
SEED = 20261026
# 洪水後の画像のセンサー名 → 種類（綴りの揺れは小文字・記号除去で吸収）
SAR = (
    "sentinel1",
    "terrasar",
    "tandem",
    "cosmo",
    "skymed",
    "radarsat",
    "rcm",
    "alos",
    "palsar",
    "iceye",
    "capella",
    "paz",
    "saocom",
    "kompsat5",
    "gaofen3",
    "gf3",
    "risat",
    "umbra",
)
OPTICAL = (
    "sentinel2",
    "landsat",
    "spot",
    "pleiades",
    "worldview",
    "geoeye",
    "kompsat3",
    "deimos",
    "superview",
    "planet",
    "rapideye",
    "triplesat",
    "gaofen1",
    "gaofen2",
    "gaofen6",
    "gf1",
    "gf2",
    "gf6",
    "cbers",
    "resourcesat",
    "cartosat",
    "ikonos",
    "quickbird",
    "formosat",
    "kanopus",
    "resurs",
    "jilin",
    "beijing",
    "aster",
    "skysat",
    "vision1",
    "ziyuan",
    "zy3",
    "dmc",
    "aerial",
    "ortho",
    "uav",
    "drone",
    "sentinel3",
    "modis",
    "viirs",
    "vrss",
    "alsat",
    "nigeriasat",
    "legion",
)


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (s or "").lower().replace("é", "e"))


def kind(name: str) -> str:
    n = norm(name)
    if any(k in n for k in SAR):
        return "sar"
    if any(k in n for k in OPTICAL):
        return "optical"
    return "other"


# %% 抽出（クラスタ）
def extract() -> None:
    import pyogrio

    aois = sorted(
        json.loads((RUNS_DIR / "world_v2" / "aoi_table_primary.json").read_text())
    )
    rows = {}
    for a in aois:
        meta = json.loads(
            (DATA_DIR / "global" / "labels" / a / "meta.json").read_text()
        )
        code = meta["activation"]["code"]
        zp = DATA_DIR / "global" / "ems" / "zips" / code / meta["product"]["zip_name"]
        row = {
            "schema_generation": meta.get("schema_generation"),
            "zip": zp.name,
            "sources": [],
        }
        try:
            names = [
                n
                for n in zipfile.ZipFile(zp).namelist()
                if "source" in n.lower() and n.lower().endswith(".dbf")
            ]
        except (OSError, zipfile.BadZipFile) as e:
            row["error"] = repr(e)
            rows[a] = row
            continue
        for n in names:
            meta_l, _, _, fields = pyogrio.raw.read(
                f"/vsizip/{zp}/{n}", read_geometry=False
            )
            cols = [str(c) for c in meta_l["fields"]]
            for vals in zip(*fields, strict=True):
                rec = {
                    c: (v.item() if hasattr(v, "item") else v)
                    for c, v in zip(cols, vals, strict=True)
                }
                rows_src = {
                    "name": rec.get("source_nam", ""),
                    "phase": rec.get("eventphase", ""),
                    "gsd": rec.get("sensor_gsd", rec.get("sensor_res", "")),
                    "date": rec.get("src_date", ""),
                }
                row["sources"].append({k: str(v) for k, v in rows_src.items()})
        rows[a] = row
    dst = DATA_DIR / RAW
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"{len(rows)} AOIs -> {dst}")


# %% 分類と集計（ローカル）
def classify(raw: dict) -> dict[str, dict]:
    out = {}
    for a, r in raw.items():
        post = [s for s in r["sources"] if "post" in s["phase"].lower()]
        kinds = sorted({kind(s["name"]) for s in post} - {"other"})
        gsd = [
            float(m.group(1))
            for s in post
            if (m := re.search(r"(\d+(?:\.\d+)?)\s*m", s["gsd"]))
            and kind(s["name"]) != "other"
        ]
        cls = (
            "unknown"
            if not kinds
            else "sar"
            if kinds == ["sar"]
            else "optical"
            if kinds == ["optical"]
            else "mixed"
        )
        out[a] = {
            "class": cls,
            "post_sensors": sorted({s["name"] for s in post}),
            "min_gsd_m": min(gsd) if gsd else None,
            "schema_generation": r.get("schema_generation"),
        }
    return out


def load_primary() -> dict[str, np.ndarray]:
    parts = [
        np.load(f)
        for f in sorted(
            (RUNS_DIR / "world_v2" / "predictions" / "primary").glob("fold*.npz")
        )
    ]
    keys = ("y", "aoi", "unit", "p_M0", "p_M1", "p_M1c", "p_M3n", "p_M3np")
    return {k: np.concatenate([p[k] for p in parts]) for k in keys}


def group_increment(d: dict[str, np.ndarray], mask: np.ndarray) -> dict:
    """群内の ΔBS(M1c → M3n)・ΔBS(M0 → M1)・R と、場所グループの cluster bootstrap 区間。"""
    y = d["y"][mask].astype(np.float64)
    err = {
        a: (d[f"p_{a}"][mask].astype(np.float64) - y) ** 2
        for a in ("M0", "M1", "M1c", "M3n")
    }
    units, inv = np.unique(d["unit"][mask], return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    s = {a: np.bincount(inv, weights=v) for a, v in err.items()}
    bs = {a: v.sum() / n.sum() for a, v in s.items()}
    dn, d01 = bs["M1c"] - bs["M3n"], bs["M0"] - bs["M1"]
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(n), np.full(len(n), 1 / len(n)), size=B).astype(np.float64)
    bb = {a: (w @ v) / (w @ n) for a, v in s.items()}
    bd = bb["M1c"] - bb["M3n"]
    return {
        "n_cells": int(mask.sum()),
        "n_units": len(units),
        "wet_rate": float(y.mean()),
        "delta_bs_names": float(dn),
        "ci90_delta": [float(np.quantile(bd, q)) for q in (0.05, 0.95)],
        "delta_bs_M0_M1": float(d01),
        "R": float(dn / d01) if d01 else None,
        "ci90_R": [float(np.quantile(bd, q) / d01) for q in (0.05, 0.95)]
        if d01
        else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--extract", action="store_true")
    args = ap.parse_args()
    if args.extract:
        extract()
        return
    raw = json.loads((DATA_DIR / RAW).read_text(encoding="utf-8"))
    cls = classify(raw)
    d = load_primary()
    counts = {
        c: sum(v["class"] == c for v in cls.values())
        for c in ("sar", "optical", "mixed", "unknown")
    }
    by_gen = {}
    for v in cls.values():
        by_gen.setdefault(str(v["schema_generation"]), {}).setdefault(v["class"], 0)
        by_gen[str(v["schema_generation"])][v["class"]] += 1
    aoi_cls = np.array([cls.get(a, {"class": "unknown"})["class"] for a in d["aoi"]])
    res = {"all": group_increment(d, np.ones(len(d["y"]), dtype=bool))}
    for c in ("sar", "optical", "mixed", "unknown"):
        m = aoi_cls == c
        if m.sum() and d["y"][m].sum() >= 20:
            res[c] = group_increment(d, m)
    print("AOI classes:", counts, "by schema generation:", by_gen)
    print(
        f"{'group':8s} {'cells':>9s} {'units':>5s} {'wet':>6s} {'dBS names':>10s} {'R':>8s}  90% CI of R"
    )
    for k, v in res.items():
        r = f"{v['R']:+.4f}" if v["R"] is not None else "   —"
        ci = f"[{v['ci90_R'][0]:+.4f}, {v['ci90_R'][1]:+.4f}]" if v["ci90_R"] else ""
        print(
            f"{k:8s} {v['n_cells']:9,d} {v['n_units']:5d} {v['wet_rate']:6.3f} "
            f"{v['delta_bs_names']:+10.2e} {r}  {ci}"
        )
    sensors = {}
    for v in cls.values():
        for s in v["post_sensors"]:
            sensors[s] = sensors.get(s, 0) + 1
    OUT.write_text(
        json.dumps(
            {
                "aoi_class_counts": counts,
                "by_schema_generation": by_gen,
                "post_event_sensor_counts": dict(
                    sorted(sensors.items(), key=lambda x: -x[1])
                ),
                "increments": res,
                "aoi_classes": cls,
            },
            indent=1,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
