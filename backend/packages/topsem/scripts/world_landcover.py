"""v2 全球検定の土地被覆特徴（洪水の前年版の ESA CCI / C3S、`topsem.world.landcover`）。

AOI ごとに ``labels/<tag>/landcover.parquet``（cells.parquet と同じ行順）と
``landcover_meta.json``（版の年・タイル URL・群の定義）を書く。ラベルは読まない
（meta.json の格子と activation の事象時刻だけを使う）。

    uv run python packages/topsem/scripts/world_landcover.py --workers 6 --shard 0/4
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

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from topsem.config import DATA_DIR
from topsem.world.landcover import cci_features

ROOT = DATA_DIR / "global"


def one(d: Path, *, force: bool = False) -> dict:
    out = d / "landcover.parquet"
    if not force and out.is_file():
        return {"tag": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        act = meta.get("activation", {})
        event_time = str(act.get("event_time") or act.get("activation_time") or "")
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        g = meta["grid"]
        flat = (cells["iy"].to_numpy() * g["nx"] + cells["ix"].to_numpy()).astype(np.int64)
        feats, prov = cci_features(meta, event_time)
        pq.write_table(pa.table({k: v[flat] for k, v in feats.items()}), out,
                       compression="zstd")
        prov["built_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        (d / "landcover_meta.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
        return {"tag": d.name, "status": "ok", "year": prov["year"], "n": len(cells),
                "sec": round(time.time() - t0, 1)}
    except Exception as e:  # noqa: BLE001
        return {"tag": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1200:]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--shard", default="0/1", help="i/n: AOI を n 分割した i 番目だけ処理")
    args = ap.parse_args()
    dirs = sorted(p for p in (ROOT / "labels").glob("EMSR*") if (p / "cells.parquet").is_file())
    i, n = (int(v) for v in args.shard.split("/"))
    dirs = dirs[i::n]
    print(f"AOIs: {len(dirs)} (shard {args.shard}) workers={args.workers}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(one, d, force=args.force): d for d in dirs}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    (ROOT / "labels" / f"_landcover_{stamp}_{i}of{n}.json").write_text(
        json.dumps(rows, indent=1), encoding="utf-8")
    print(f"→ {sum(r['status'] in ('ok', 'cached') for r in rows)}/{len(rows)} ok")


if __name__ == "__main__":
    main()
