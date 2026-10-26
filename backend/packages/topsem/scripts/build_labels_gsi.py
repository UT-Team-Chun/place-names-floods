"""国土地理院「浸水推定段彩図/浸水推定図 輪郭線」→ 250 m ラベルグリッド。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/build_labels_gsi.py \
        --extracted ../data/labels/<event>/raw/extracted \
        --event 2024_mogami \
        --source-url https://www.gsi.go.jp/common/000260109.zip \
        --version-date 20240726 \
        --role robustness_check_only

role は protocol §3 の使途（primary / robustness_check_only / snapshot_only）。
地理院の免責（推定混じり・時点情報）を metadata に必ず残す。
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.config import LABELS_DIR
from topsem.labels.grid import rasterize_wet_frac, rasterize_window
from topsem.labels.sources import load_gsi_shinsui_dir
from topsem.runs import write_json

# 凍結値（protocol §2）。結果を見て変更しない。
WET_TH = 0.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extracted", type=Path, required=True)
    parser.add_argument("--event", required=True, help="data/labels/<event>/ の event 名")
    parser.add_argument("--source-url", required=True)
    parser.add_argument("--version-date", required=True, help="地理院作成日 (YYYYMMDD)")
    parser.add_argument("--role", default="robustness_check_only")
    parser.add_argument("--size-m", type=float, default=250.0)
    args = parser.parse_args()

    retrieved = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    extents = load_gsi_shinsui_dir(
        args.extracted,
        version_date=args.version_date,
        source_url=args.source_url,
        retrieved_utc=retrieved,
    )
    out_root = LABELS_DIR / args.event / "derived"
    n_total = 0
    for ext in extents:
        cells = rasterize_wet_frac(ext.union_geom, size_m=args.size_m)
        window = rasterize_window(ext.union_geom, size_m=args.size_m)
        river = ext.river
        river_dir = out_root / river.replace("/", "_")
        river_dir.mkdir(parents=True, exist_ok=True)
        with (river_dir / "grid.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ix", "iy", "wet_frac"])
            w.writerows(cells)
        with (river_dir / "grid_window.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ix", "iy", "wet_frac"])
            w.writerows(window)
        wet_cells = sum(1 for _, _, fr in cells if fr >= WET_TH)
        cell_km2 = (args.size_m / 1000.0) ** 2
        minx, miny, maxx, maxy = ext.union_geom.bounds
        write_json(
            river_dir / "metadata.json",
            {
                "event": args.event,
                "river": river,
                "grid_size_m": args.size_m,
                "wet_threshold": WET_TH,
                "role": args.role,
                "lat_ref": (miny + maxy) / 2.0,
                "bbox": [minx, miny, maxx, maxy],
                "n_cells_touched": len(cells),
                "n_window_cells": len(window),
                "n_cells_wet": wet_cells,
                "wet_area_km2_grid": sum(fr for _, _, fr in cells) * cell_km2,
                "source": ext.meta,
            },
        )
        n_total += len(cells)
        print(f"{river:36s} cells={len(cells):6d} wet>={WET_TH}: {wet_cells:6d}")
    print(f"-- total {n_total} cells -> {out_root}")


if __name__ == "__main__":
    main()
