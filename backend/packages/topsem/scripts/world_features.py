"""AOI ごとのラベルに地形・土地被覆・河道距離・地名割当を接合（v2 全球検定の M3 前段）。

前提: ``world_labels.py`` で ``data/global/labels/<tag>/{cells.parquet,meta.json}`` があること。
初回は ``--prepare`` で GeoNames（全件 → AOI 範囲の部分集合）と HydroRIVERS を用意する。

出力: ``data/global/labels/<tag>/features.parquet``（cells と同じ行順）。ラベルと特徴の
関係（モデル fit・相関）はここでは計算しない（protocol v2 凍結前の工程）。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/world_features.py --prepare
    uv run python packages/topsem/scripts/world_features.py --workers 8 [--tags EMSR772_AOI01]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyproj import CRS

from topsem.config import DATA_DIR
from topsem.world.hydro import ensure_downloaded, river_distances
from topsem.world.names import assign_names, geonames_to_parquet, subset_for_bounds
from topsem.world.terrain import grid_lonlat_bounds, provenance, terrain_features

ROOT = DATA_DIR / "global"
GN_ZIP_URL = "https://download.geonames.org/export/dump/allCountries.zip"
_GN: pd.DataFrame | None = None


def _label_dirs(tags: str = "") -> list[Path]:
    dirs = sorted(p for p in (ROOT / "labels").glob("EMSR*") if (p / "cells.parquet").is_file())
    if tags:
        want = {t.strip() for t in tags.split(",") if t.strip()}
        dirs = [d for d in dirs if d.name in want]
    return dirs


def _download(url: str, dst: Path) -> None:
    """既定 UA を拒否するサーバがあるので研究用 UA を明示して取得。"""
    req = urllib.request.Request(  # noqa: S310
        url, headers={"User-Agent": "toponym-semantics research (ryoutaokauchi@g.ecc.u-tokyo.ac.jp)"})
    tmp = dst.with_suffix(dst.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as f:  # noqa: S310
        while chunk := resp.read(1 << 20):
            f.write(chunk)
    tmp.rename(dst)


def prepare() -> None:
    gdir = ROOT / "geonames"
    gdir.mkdir(parents=True, exist_ok=True)
    zp = gdir / "allCountries.zip"
    if not zp.is_file():
        _download(GN_ZIP_URL, zp)
    full = gdir / "allCountries.parquet"
    if not full.is_file():
        n = geonames_to_parquet(zp, full)
        print(f"geonames rows: {n:,}")
    gn = pd.read_parquet(full)
    bounds = [grid_lonlat_bounds(json.loads((d / "meta.json").read_text())) for d in _label_dirs()]
    sub = subset_for_bounds(gn, bounds)
    sub.to_parquet(gdir / "aoi_subset.parquet", compression="zstd", index=False)
    print(f"geonames subset for {len(bounds)} AOIs: {len(sub):,} rows")
    ensure_downloaded(ROOT / "hydrorivers")
    print("hydrorivers ready")
    coast = ROOT / "coast" / "ne_10m_coastline.zip"  # §2.4 内陸判定（種別未記載の製品）
    if not coast.is_file():
        coast.parent.mkdir(parents=True, exist_ok=True)
        _download("https://naciscdn.org/naturalearth/10m/physical/ne_10m_coastline.zip", coast)
    print("coastline ready")
    ci = gdir / "countryInfo.txt"  # 国 → 大陸（場所グループの fold 割当に使う）
    if not ci.is_file():
        _download("https://download.geonames.org/export/dump/countryInfo.txt", ci)
    print("countryInfo ready")


def _init_worker() -> None:
    global _GN  # noqa: PLW0603
    _GN = pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")


def _gn() -> pd.DataFrame:
    if _GN is None:
        msg = "worker の GeoNames 部分集合が未ロード"
        raise RuntimeError(msg)
    return _GN


def one(d: Path, *, force: bool = False) -> dict:
    out = d / "features.parquet"
    if not force and out.is_file():
        return {"tag": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet").to_pandas()
        g = meta["grid"]
        flat = (cells["iy"].to_numpy() * g["nx"] + cells["ix"].to_numpy()).astype(np.int64)
        terr = terrain_features(meta)
        bounds = grid_lonlat_bounds(meta)
        zips = sorted((ROOT / "hydrorivers").glob("HydroRIVERS_v10_*_shp.zip"))
        hyd = river_distances(meta, zips, bounds)
        names = assign_names((cells["x"].to_numpy(dtype=float), cells["y"].to_numpy(dtype=float)),
                             CRS.from_proj4(meta["crs_proj4"]), _gn(), bounds)
        cols = {k: v[flat] for k, v in {**terr, **hyd}.items()}
        cols |= names
        pq.write_table(pa.table(cols), out, compression="zstd")
        prov = {"terrain": provenance(), "hydro": "HydroRIVERS v1.0 (all regions, bbox read)",
                "names": "GeoNames allCountries (see data/global/geonames/*.sha256)",
                "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        (d / "features_meta.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
        return {"tag": d.name, "status": "ok", "n": len(cells), "sec": round(time.time() - t0, 1)}
    except Exception as e:  # noqa: BLE001
        return {"tag": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1200:]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--tags", default="")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--shard", default="0/1", help="i/n: AOI を n 分割した i 番目だけ処理")
    args = ap.parse_args()
    if args.prepare:
        prepare()
        return
    dirs = _label_dirs(args.tags)
    i, n = (int(v) for v in args.shard.split("/"))
    dirs = dirs[i::n]
    print(f"AOIs: {len(dirs)} (shard {args.shard}) workers={args.workers}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers), initializer=_init_worker) as ex:
        futs = {ex.submit(one, d, force=args.force): d for d in dirs}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    (ROOT / "labels" / f"_features_{stamp}_{i}of{n}.json").write_text(
        json.dumps(rows, indent=1), encoding="utf-8")
    print(f"→ {sum(r['status'] in ('ok', 'cached') for r in rows)}/{len(rows)} ok")


if __name__ == "__main__":
    main()
