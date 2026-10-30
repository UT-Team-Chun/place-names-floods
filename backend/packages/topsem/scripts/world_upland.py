"""amendment B-16 の C1: 最寄りの河道の上流の面積（features_upland.parquet）を、解析する AOI ごとに作る。

    python world_upland.py           # 解析する全ての AOI（HydroRIVERS は大陸の zip ごとに一度だけ読む）
    python world_upland.py --check   # 書いたファイルの行数と、値の得られたセルの割合

出力: 各 AOI のラベルの場所に ``features_upland.parquet``（cells.parquet と同じ行の順。lupl = 最寄りの河道の頂点の
UPLAND_SKM の log10、dist_upl_m = その頂点までの距離）と ``features_upland_meta.json``。頂点は AOI の範囲 ± 0.2°
（world_fabdem.py の river_vertices と同じ範囲）。範囲に頂点が無いセルは NaN。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

import world_eval as ev
from world_fabdem import _padded, included_dirs

from topsem.world.b9 import cell_lonlat
from topsem.world.b16 import in_box, merge_nearest, nearest_upland
from topsem.world.terrain import grid_lonlat_bounds

OUT_NAME = "features_upland.parquet"


def continent_vertices(z: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """大陸の zip の全ての河道の頂点の経度・緯度・上流の面積（km²）。"""
    import pyogrio.raw
    import shapely

    with zipfile.ZipFile(z) as zf:
        shp = next(x for x in zf.namelist() if x.endswith(".shp"))
    _, _, wkb, fields = pyogrio.raw.read(f"/vsizip/{z}/{shp}", columns=["UPLAND_SKM"])
    xy, ix = shapely.get_coordinates(shapely.from_wkb(wkb), return_index=True)
    return xy[:, 0], xy[:, 1], np.asarray(fields[0], dtype=np.float64)[ix]


def _one_degree_index(vlon: np.ndarray, vlat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    key = ((np.floor(vlat) + 90) * 360 + (np.floor(vlon) + 180)).astype(np.int64)
    order = np.argsort(key, kind="stable")
    return order, key[order]


def _candidates(order: np.ndarray, key_s: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    w, s, e, n = box
    lo_lon, hi_lon = int(np.floor(max(w, -180.0))) + 180, int(np.floor(min(e, 179.999999))) + 180
    parts = []
    for la in range(int(np.floor(max(s, -90.0))), int(np.floor(min(n, 89.999999))) + 1):
        k0, k1 = (la + 90) * 360 + lo_lon, (la + 90) * 360 + hi_lon
        a, b = np.searchsorted(key_s, [k0, k1 + 1])
        parts.append(order[a:b])
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)


def build() -> None:
    t0 = time.time()
    aois = []
    for d in included_dirs():
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        flat = (cells["iy"].to_numpy() * meta["grid"]["nx"] + cells["ix"].to_numpy()).astype(np.int64)
        lon, lat = cell_lonlat(meta)
        aois.append({"d": d, "box": _padded(grid_lonlat_bounds(meta)), "lon": lon[flat], "lat": lat[flat],
                     "best_d": np.full(len(cells), np.inf), "best_v": np.full(len(cells), np.nan), "n_vert": 0})
    print(f"{len(aois)} AOIs ({time.time() - t0:.0f}s)", flush=True)
    zips = sorted((ev.ROOT / "hydrorivers").glob("HydroRIVERS_v10_*_shp.zip"))
    for z in zips:
        vlon, vlat, vup = continent_vertices(z)
        order, key_s = _one_degree_index(vlon, vlat)
        hit = 0
        for a in aois:
            c = _candidates(order, key_s, a["box"])
            c = c[in_box(vlon[c], vlat[c], a["box"])]
            if len(c) == 0:
                continue
            lupl, dist = nearest_upland(a["lon"], a["lat"], vlon[c], vlat[c], vup[c])
            merge_nearest(a["best_d"], a["best_v"], dist.astype(np.float64), lupl.astype(np.float64))
            a["n_vert"] += len(c)
            hit += 1
        print(f"{z.name}: {len(vlon):,} vertices, {hit} AOIs ({time.time() - t0:.0f}s)", flush=True)
    for a in aois:
        dist = np.where(np.isfinite(a["best_d"]), a["best_d"], np.nan)
        pq.write_table(pa.table({"lupl": a["best_v"].astype(np.float32), "dist_upl_m": dist.astype(np.float32)}),
                       a["d"] / OUT_NAME, compression="zstd")
        info = {"aoi": a["d"].name, "n_cells": len(dist), "n_river_vertices": a["n_vert"],
                "valid_share": float(np.isfinite(a["best_v"]).mean()), "zips": [z.name for z in zips],
                "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        (a["d"] / "features_upland_meta.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    print(f"written {len(aois)} files ({time.time() - t0:.0f}s)", flush=True)


def check() -> dict:
    out = {"n_aoi": 0, "row_mismatch": [], "missing": [], "valid_share_min": 1.0, "cells": 0, "valid_cells": 0}
    for d in included_dirs():
        f = d / OUT_NAME
        if not f.is_file():
            out["missing"].append(d.name)
            continue
        n = pq.read_metadata(d / "cells.parquet").num_rows
        t = pq.read_table(f, columns=["lupl"]).to_pandas()
        if len(t) != n:
            out["row_mismatch"].append(d.name)
        ok = np.isfinite(t["lupl"].to_numpy(np.float64))
        out["n_aoi"] += 1
        out["cells"] += len(t)
        out["valid_cells"] += int(ok.sum())
        out["valid_share_min"] = min(out["valid_share_min"], float(ok.mean()) if len(t) else 1.0)
    out["valid_share"] = out["valid_cells"] / max(out["cells"], 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if args.check:
        print(json.dumps(check(), indent=1), flush=True)
        return
    build()


if __name__ == "__main__":
    main()
