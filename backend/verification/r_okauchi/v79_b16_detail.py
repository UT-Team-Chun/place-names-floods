# %% v79（事後・amendment B-16 の C5）: EMS の詳細図と概観図のラベル
"""確証判定は変えない。B-16 の追記 1 の 15 で計算の前に固定した量。

解析した AOI のうち、旧い配布で同じ AOI 番号の DETAIL の製品（MONIT を除く、DETAIL の番号ごとに最新の版）があり、
ラベルが DETAIL 以外の製品から作られた AOI。詳細図は全ての DETAIL の製品の層を合わせて一つの製品とし、AOI の格子に
同じ規則で落とす。比べるセルは両方で評価できるセル。

    uv run python verification/r_okauchi/v79_b16_detail.py --fetch   # 詳細図の製品を取得（data/global/ems/detail_zips/）
    uv run python verification/r_okauchi/v79_b16_detail.py --build   # AOI の格子に落とす（data/verification/v79_detail/）
    uv run python verification/r_okauchi/v79_b16_detail.py --score   # 比べる（docs/verification/v79_b16_detail.json）
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR
from topsem.world.b16 import aoi_grid, detail_fractions
from topsem.world.dataset import evaluable, riverine_labels
from topsem.world.ems import _detail_no, download

LABELS = DATA_DIR / "global" / "labels"
ZIPS = DATA_DIR / "global" / "ems" / "detail_zips"
CELLS = DATA_DIR / "verification" / "v79_detail"
OUT = ROOT / "docs" / "verification" / "v79_b16_detail.json"
INVENTORY = DATA_DIR / "global" / "ems" / "inventory.json"
SEED = 20261105
NEAR_M = 500.0
WHOLE_AOI = 0.9  # 追記 2
F_COLS = ("f_riverine", "f_flash", "f_coastal", "f_other", "f_trace", "f_flood_any", "f_aoi", "f_observed",
          "f_perm_water")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def plan() -> dict[str, list[dict]]:
    """AOI → 使う詳細図の製品（DETAIL の番号ごとに最新の版）。ラベルが DETAIL の製品から作られた AOI は除く。"""
    v77 = _load("v77_b16_local")
    aois = sorted(v77._superclusters()["aoi_sc"])
    acts = {a["code"]: a for a in _json(INVENTORY)["activations"]}
    out: dict[str, list[dict]] = {}
    for tag in aois:
        prod = _json(LABELS / tag / "meta.json").get("product") or {}
        if str(prod.get("scale", "")).upper() == "DETAIL":
            continue
        act, num = tag.split("_AOI")
        dets = [p for p in (acts.get(act) or {}).get("del_products") or []
                if p.get("kind") == "PRODUCT" and str(p.get("scale", "")).upper() == "DETAIL" and p.get("aoi") == num]
        latest: dict[int, dict] = {}
        for p in dets:
            k = _detail_no(p["url"])
            if k not in latest or int(p.get("version", 0)) > int(latest[k].get("version", 0)):
                latest[k] = p
        if latest:
            out[tag] = [latest[k] for k in sorted(latest)]
    return out


def _zip_path(url: str) -> Path:
    return ZIPS / url.rsplit("/", 1)[-1]


def fetch() -> dict:
    pl = plan()
    urls = sorted({p["url"] for ps in pl.values() for p in ps})
    status: dict[str, str] = {}

    def one(u: str) -> tuple[str, str]:
        try:
            return u, download(u, _zip_path(u))
        except Exception as e:  # noqa: BLE001
            return u, f"error: {type(e).__name__}: {e}"

    with ThreadPoolExecutor(8) as ex:
        for u, s in ex.map(one, urls):
            status[u] = s
    n_err = sum(s.startswith("error") for s in status.values())
    print(f"{len(pl)} AOIs, {len(urls)} products, {n_err} errors", flush=True)
    return {"n_aoi": len(pl), "n_products": len(urls), "errors": {u: s for u, s in status.items()
                                                                     if s.startswith("error")}}


def build() -> dict:
    pl = plan()
    v77 = _load("v77_b16_local")
    untyped = v77._superclusters()["aoi_untyped"]
    CELLS.mkdir(parents=True, exist_ok=True)
    info: dict = {}
    for tag, prods in pl.items():
        zips = [_zip_path(p["url"]) for p in prods if _zip_path(p["url"]).is_file()]
        if not zips:
            info[tag] = "no product downloaded"
            continue
        meta = _json(LABELS / tag / "meta.json")
        grid = aoi_grid(meta)
        fr = detail_fractions(zips, grid)
        c = pq.read_table(LABELS / tag / "cells.parquet").to_pandas()
        flat = (c["iy"].to_numpy() * grid.nx + c["ix"].to_numpy()).astype(np.int64)
        d = pd.DataFrame({k: np.asarray(fr[k])[flat] for k in F_COLS})
        flag = bool(untyped.get(tag, False))
        y_o, keep_o = riverine_labels(c, unspecified_as_riverine=flag)
        y_d, keep_d = riverine_labels(d, unspecified_as_riverine=flag)
        out = pd.DataFrame({"ix": c["ix"].to_numpy(), "iy": c["iy"].to_numpy(),
                            "y_over": y_o, "ok_over": evaluable(c) & keep_o,
                            "y_det": y_d, "ok_det": evaluable(d) & keep_d, "det_aoi": d["f_aoi"].to_numpy()})
        out.to_parquet(CELLS / f"{tag}.parquet")
        info[tag] = {"n_products": len(zips), "generations": fr["generations"].tolist(),
                     "n_cells_in_detail": int((out["det_aoi"] >= 0.5).sum()),
                     "n_both_evaluable": int((out["ok_over"] & out["ok_det"]).sum())}
    return info


def _agree(e: np.ndarray, g: np.ndarray) -> dict:
    """e = 詳細図（基準）、g = 概観図。"""
    v77 = _load("v77_b16_local")
    return v77._agree(e, g)


def score() -> dict:
    v77 = _load("v77_b16_local")
    v39 = v77.v39
    t41 = pd.read_parquet(DATA_DIR / "verification" / "v41_table.parquet")
    t37 = pd.read_parquet(v39.TABLE)
    if not (np.array_equal(t37["aoi"].to_numpy(), t41["aoi"].to_numpy())
            and np.array_equal(t37["y"].to_numpy(), t41["y"].to_numpy())):
        return {"note": "v37 and v41 tables are not aligned"}
    t = t37.assign(ix=t41["ix"].to_numpy(), iy=t41["iy"].to_numpy(), p_M0=t41["p_M0"].to_numpy())
    files = sorted(CELLS.glob("*.parquet"))
    if not files:
        return {"note": "run --build first"}
    # 追記 2: 範囲の層を合わせたものが AOI のセルの 0.9 以上を覆う AOI（範囲の層が AOI 全体を指す）は除く
    whole, parts = [], []
    for f in files:
        x = pd.read_parquet(f)
        if float((x["det_aoi"].to_numpy() >= 0.5).mean()) >= WHOLE_AOI:
            whole.append(f.stem)
            continue
        parts.append(x.assign(aoi=f.stem))
    det = pd.concat(parts, ignore_index=True)
    m = t.merge(det, on=["aoi", "ix", "iy"], how="inner", validate="one_to_one")
    both = m[m["ok_over"] & m["ok_det"]].reset_index(drop=True)
    check_same = bool((both["y_over"].to_numpy() == both["y"].to_numpy()).all())
    wet = both["lex_wetside_n1"].to_numpy() > 0
    band = v77._band_index(both["P_d1"].to_numpy())
    e, g = both["y_det"].to_numpy(), both["y_over"].to_numpy()
    out: dict = {"removed_whole_aoi_layer": whole, "n_aoi": int(both["aoi"].nunique()), "n_cells": len(both),
                 "overview_labels_match_analysis": check_same,
                 "agreement": {"all": _agree(e, g)}}
    for i, bn in enumerate(("b500", "b1k", "bfar")):
        for w, wn in ((True, "wet"), (False, "other")):
            sel = (band == i) & (wet == w)
            out["agreement"][f"{bn}_{wn}"] = _agree(e[sel], g[sel])
    rng = np.random.default_rng(SEED)
    tt = both.assign(y_d=both["y_det"].astype(np.float64), y_o=both["y_over"].astype(np.float64))
    ones = np.ones(len(tt), dtype=bool)
    out["scoring"] = {"detail_labels": v77._r_fixed(tt, ones, rng, ycol="y_d"),
                      "overview_labels": v77._r_fixed(tt, ones, rng, ycol="y_o")}
    out["eq2"] = {}
    for key, ycol in (("detail_labels", "y_det"), ("overview_labels", "y_over")):
        s = both.assign(y=both[ycol].astype(np.float64))
        rate = s.groupby("aoi")["y"].transform("mean")
        s = s[(rate > 0) & (rate < 1)].reset_index(drop=True)
        if s.empty:
            out["eq2"][key] = {"n_cells": 0}
            continue
        z, names = v39.design(s, ("wetside",), True, banded=True)
        j = names.index("wetside:le500m")
        near = (s["P_d1"].to_numpy() <= NEAR_M) & (s["lex_wetside_n1"].to_numpy() > 0) & (s["y"].to_numpy() == 1)
        try:
            r = v39.fit(s["y"].to_numpy(np.float64), z, s["aoi"].to_numpy(), s["unit"].to_numpy(str))
            b, se = float(r["beta"][j]), float(r["se"][j])
            out["eq2"][key] = {"odds_ratio": float(np.exp(b)), "ci95": [float(np.exp(b - 1.96 * se)),
                                                                         float(np.exp(b + 1.96 * se))],
                               "n_cells": len(s), "n_aoi": int(s["aoi"].nunique()), "n_flooded_near": int(near.sum())}
        except np.linalg.LinAlgError as ex:
            out["eq2"][key] = {"error": str(ex), "n_flooded_near": int(near.sum())}
    out["rule"] = ("B-16 C5 (addendum 1, item 15): DETAIL products of the old distribution with the AOI's number (latest "
                   "version of each DETAIL number), layers taken together and rasterized onto the AOI's grid with the "
                   "same rules; cells of the analysed (capped) table evaluable on both; detail labels as the reference; "
                   "confirmatory predictions (models fixed), 4,000 draws of location groups; Eq. (2) on the v37 table")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()
    res = _json(OUT) if OUT.is_file() else {}
    for key, fn, on in (("fetch", fetch, a.fetch), ("build", build, a.build), ("score", score, a.score)):
        if not on:
            continue
        t0 = time.time()
        res[key] = fn()
        res[f"{key}_seconds"] = round(time.time() - t0)
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        print(f"{key} ({res[f'{key}_seconds']} s) {json.dumps(res[key], default=float)[:1500]}", flush=True)


if __name__ == "__main__":
    main()
