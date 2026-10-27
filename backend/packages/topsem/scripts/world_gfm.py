"""事後 B-3: 解析した AOI の独立の浸水ラベル（Copernicus GFM、`topsem.world.gfm`）。

AOI ごとに ``labels/<tag>/gfm.parquet``（cells.parquet と同じ行順）と ``gfm_meta.json``（窓・撮像・項目の
ID と URL・処理版）を書く。EMS の浸水ラベルは読まない（meta.json の格子と、v38 が製品から読んだ出典表の
事後画像の日付だけを使う）。

``--summary`` は全シャードの完了後に、解析した AOI ごとの採用（B-3: EMS と GFM の両方で評価可能なセルが
400 以上、GFM で浸水 20 以上、浸水率 0.8 以下）と事後画像の種類を ``gfm_summary.json`` にまとめる。

    uv run python packages/topsem/scripts/world_gfm.py --workers 4 --shard 0/4
    uv run python packages/topsem/scripts/world_gfm.py --summary
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from topsem.config import DATA_DIR
from topsem.world.dataset import evaluable, riverine_labels
from topsem.world.gfm import (
    gfm_include,
    gfm_labels,
    gfm_labels_v1,
    post_event_window,
    sensor_group,
)

ROOT = DATA_DIR / "global"
SOURCES = DATA_DIR / "verification" / "v38_ems_sensors_raw.json"
SUMMARY = ROOT / "gfm_summary.json"


def one(d: Path, sources: list[dict], *, force: bool = False) -> dict:
    out = d / "gfm.parquet"
    if not force and out.is_file():
        return {"tag": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        ix, iy = cells["ix"].to_numpy(), cells["iy"].to_numpy()
        cols, prov = gfm_labels(meta, sources, ix, iy)
        prov["built_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        (d / "gfm_meta.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
        tmp = out.with_suffix(".parquet.tmp")  # parquet があれば完了（途中で止まっても半端な表を残さない）
        pq.write_table(pa.table({"ix": ix, "iy": iy, **cols}), tmp, compression="zstd")
        tmp.replace(out)
        return {"tag": d.name, "status": "ok", "reason": prov["reason"],
                "n_scenes": prov.get("n_scenes", 0), "n_eval": int(cols["gfm_eval"].sum()),
                "n": len(cells), "sec": round(time.time() - t0, 1)}
    except Exception as e:  # noqa: BLE001
        return {"tag": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1200:]}


def one_v1(d: Path, sources: list[dict], *, force: bool = False) -> dict:
    """事後 B-6 の 4: 層の成分と変種 V1 を ``gfm_v1.parquet`` に（B-3 のラベルも作り直して照合する）。"""
    out = d / "gfm_v1.parquet"
    if not force and out.is_file():
        return {"tag": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        ix, iy = cells["ix"].to_numpy(), cells["iy"].to_numpy()
        cols, prov = gfm_labels_v1(meta, sources, ix, iy)
        old = pq.read_table(d / "gfm.parquet", columns=["gfm_eval", "gfm_y"]).to_pandas()
        same = bool((old["gfm_eval"].to_numpy() == cols["gfm_eval"]).all()
                    and (old["gfm_y"].to_numpy() == cols["gfm_y"]).all())
        prov |= {"built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "b3_labels_reproduced": same,
                 "rule": "amendment B-6 item 4"}
        (d / "gfm_v1_meta.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
        tmp = out.with_suffix(".parquet.tmp")
        pq.write_table(pa.table({"ix": ix, "iy": iy, **cols}), tmp, compression="zstd")
        tmp.replace(out)
        return {"tag": d.name, "status": "ok", "b3_reproduced": same, "n_v1_eval": int(cols["v1_eval"].sum()),
                "n": len(cells), "sec": round(time.time() - t0, 1)}
    except Exception as e:  # noqa: BLE001
        return {"tag": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1200:]}


def build_v1(args: argparse.Namespace) -> None:
    """B-6 の 4: B-3 で採用した AOI だけ。"""
    raw = json.loads(SOURCES.read_text(encoding="utf-8"))
    summ = json.loads(SUMMARY.read_text(encoding="utf-8"))
    dirs = sorted(ROOT / "labels" / t for t, v in summ["aois"].items() if v.get("include"))
    i, n = (int(v) for v in args.shard.split("/"))
    dirs = dirs[i::n]
    print(f"B-3 included AOIs: {len(dirs)} (shard {args.shard}) workers={args.workers}", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(one_v1, d, raw[d.name]["sources"], force=args.force): d for d in dirs}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    (ROOT / "labels" / f"_gfm_v1_{stamp}_{i}of{n}.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"→ {sum(r['status'] in ('ok', 'cached') for r in rows)}/{len(rows)} ok; "
          f"B-3 reproduced {sum(bool(r.get('b3_reproduced')) for r in rows)}")


def build(args: argparse.Namespace) -> None:
    raw = json.loads(SOURCES.read_text(encoding="utf-8"))
    dirs = sorted(ROOT / "labels" / t for t in raw
                  if post_event_window(raw[t].get("sources") or []) is not None)
    i, n = (int(v) for v in args.shard.split("/"))
    dirs = dirs[i::n]
    print(f"AOIs with a post-event date: {len(dirs)} (shard {args.shard}) "
          f"workers={args.workers}", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(one, d, raw[d.name]["sources"], force=args.force): d for d in dirs}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    (ROOT / "labels" / f"_gfm_{stamp}_{i}of{n}.json").write_text(
        json.dumps(rows, indent=1), encoding="utf-8")
    print(f"→ {sum(r['status'] in ('ok', 'cached') for r in rows)}/{len(rows)} ok")


def summary() -> None:
    import world_eval as ev  # 確証解析と同じ AOI の採用（select_aois）と、その内陸判定を使う

    raw = json.loads(SOURCES.read_text(encoding="utf-8"))
    gn = pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")
    dirs, _, infos = ev.select_aois(gn, ev.Sensitivity())
    out: dict[str, dict] = {}
    for d in dirs:
        sens = sensor_group(raw.get(d.name, {}).get("sources") or [])
        if not (d / "gfm.parquet").is_file():
            out[d.name] = {"include": False, "why": "no-gfm-table", "sensor": sens}
            continue
        cells = pq.read_table(d / "cells.parquet").to_pandas()
        _, keep = riverine_labels(cells,
                                  unspecified_as_riverine=infos[d.name].unspecified_as_riverine)
        g = pq.read_table(d / "gfm.parquet").to_pandas()
        if not (g["ix"].to_numpy() == cells["ix"].to_numpy()).all():
            msg = f"{d.name}: gfm.parquet rows do not match cells.parquet"
            raise SystemExit(msg)
        ok, info = gfm_include(evaluable(cells) & keep, g["gfm_eval"].to_numpy(),
                               g["gfm_y"].to_numpy())
        meta = json.loads((d / "gfm_meta.json").read_text(encoding="utf-8"))
        out[d.name] = {"include": ok, **info, "sensor": sens, "reason": meta.get("reason"),
                       "n_scenes": meta.get("n_scenes", 0)}
    why = Counter(v["why"] for v in out.values())
    doc = {"rule": "amendment B-3", "n_aoi": len(out),
           "n_include": sum(v["include"] for v in out.values()), "why": dict(why),
           "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "aois": out}
    SUMMARY.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in doc.items() if k != "aois"}, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--shard", default="0/1", help="i/n: AOI を n 分割した i 番目だけ処理")
    ap.add_argument("--summary", action="store_true", help="全シャードの後: AOI の採用をまとめる")
    ap.add_argument("--v1", action="store_true", help="事後 B-6 の 4: 層の成分と季節水域を戻した変種 V1")
    args = ap.parse_args()
    if args.summary:
        summary()
    elif args.v1:
        build_v1(args)
    else:
        build(args)


if __name__ == "__main__":
    main()
