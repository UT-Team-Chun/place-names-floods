"""EMS 台帳 → AOI ごとの 250 m 観測浸水ラベル（v2 全球検定の M2）。

``data/global/ems/inventory.json``（scripts/ems_inventory.py）を読み、activation ごとに
AOI 1 製品を選んで取得・ラスタ化し、``data/global/labels/<EMSR>_AOI<nn>/`` に
``cells.parquet`` と ``meta.json``（製品 URL・zip sha256・版・層・属性の内訳）を書く。

ラベルとモデル特徴の関係はここでは一切計算しない（protocol v2 凍結前の工程）。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/world_labels.py --codes EMSR517,EMSR271
    uv run python packages/topsem/scripts/world_labels.py --workers 8
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyarrow as pa
import pyarrow.parquet as pq

from topsem.config import DATA_DIR
from topsem.world.ems import (
    build_aoi_labels,
    cells_digest,
    download,
    select_products,
    write_meta,
)

ROOT = DATA_DIR / "global"


def one(act: dict, prod: dict, *, force: bool = False) -> dict:
    """1 AOI を処理して要約行を返す（失敗も行として返す）。force は既存を上書き再構築。"""
    tag = f"{act['code']}_AOI{prod['aoi']}"
    out = ROOT / "labels" / tag
    if not force and (out / "cells.parquet").is_file():
        return {"tag": tag, "status": "cached"}
    t0 = time.time()
    try:
        name = prod["url"].split("/")[-1].split("?")[0]
        zpath = ROOT / "ems" / "zips" / act["code"] / name
        sha = download(prod["url"], zpath)
        cols, meta = build_aoi_labels(zpath)
        out.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.table(cols), out / "cells.parquet", compression="zstd")
        write_meta(
            out / "meta.json",
            {
                "tag": tag,
                "activation": {k: act.get(k) for k in (
                    "code", "name", "countries", "activation_time", "event_time",
                    "sub_category", "continent", "system")},
                "product": {**prod, "zip_sha256": sha, "zip_name": name},
                "cells_sha256": cells_digest(out / "cells.parquet"),
                **meta,
                "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        return {"tag": tag, "status": "ok", "n_cells": meta["n_cells"],
                "sec": round(time.time() - t0, 1)}
    except Exception as e:  # noqa: BLE001
        return {"tag": tag, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-800:]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--codes", default="", help="カンマ区切りの EMSR コードに限定")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true", help="既存ラベルを上書き再構築")
    ap.add_argument("--shard", default="0/1", help="i/n: AOI ジョブを n 分割した i 番目だけ処理")
    ap.add_argument("--inventory", default="",
                    help="EMS 台帳 JSON（既定 data/global/ems/inventory.json。クラスタでは "
                         "docs/data/ems_inventory_<日付>.json の凍結版を渡す）")
    args = ap.parse_args()

    inv_path = Path(args.inventory) if args.inventory else ROOT / "ems" / "inventory.json"
    inv = json.loads(inv_path.read_text(encoding="utf-8"))
    print(f"inventory: {inv_path} (retrieved {inv.get('retrieved_utc')})", flush=True)
    acts = inv["activations"]
    if args.codes:
        want = {c.strip() for c in args.codes.split(",") if c.strip()}
        acts = [a for a in acts if a["code"] in want]
    jobs = [(a, p) for a in acts for p in select_products(a["del_products"])]
    if args.limit:
        jobs = jobs[: args.limit]
    si, sn = (int(v) for v in args.shard.split("/"))
    jobs = jobs[si::sn]
    print(f"AOI jobs: {len(jobs)}  workers={args.workers}", flush=True)
    rows = []
    if args.workers <= 1:
        for a, p in jobs:
            rows.append(one(a, p, force=args.force))
            print(json.dumps({k: v for k, v in rows[-1].items() if k != "trace"}), flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(one, a, p, force=args.force) for a, p in jobs]
            for f in as_completed(futs):
                rows.append(f.result())
                print(json.dumps({k: v for k, v in rows[-1].items() if k != "trace"}),
                      flush=True)
    log = ROOT / "labels" / (
        f"_build_{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}_{si}of{sn}.json")
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    n_ok = sum(r["status"] in ("ok", "cached") for r in rows)
    print(f"→ {n_ok}/{len(rows)} ok; log {log}")


if __name__ == "__main__":
    main()
