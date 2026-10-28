# %% v58（事後・amendment B-8 の C4）: セルごとの繰り返す水（JRC Global Surface Water v1.4）
"""解析する表（確証と同じ規則で作る上限つきの表）の各 250 m のセルについて、セルの中の 5 × 5 の点で 30 m の画素を
読み、JRC GSW の occurrence（水があった時間の割合、%）・seasonality（水のあった月の数）の平均と、max extent（一度でも
水だった画素）の割合を出す。浸水のラベルは使わない。タイルの読み方は v49（B-6 の 3(b)）と同じ。

出力: ``data/global/gsw/cell_gsw.parquet``（aoi, ix, iy, occ, seas, maxext, n_valid）と ``cell_gsw_meta.json``。
クラスタの script 段で実行する（ネットワークに出る）。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import CRS, Transformer

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import v49_gsw_settlements as v49  # noqa: E402
import world_eval as ev  # noqa: E402

from topsem.config import DATA_DIR  # noqa: E402

OUT = DATA_DIR / "global" / "gsw"
OFFS = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])  # セルの幅に対する点の位置（5 × 5）
PAD_DEG = 0.01


def cell_points(sub: pd.DataFrame, meta: dict) -> tuple[np.ndarray, np.ndarray]:
    """各セルの 25 点の経緯度（セル数 × 25）。セルの中心は AOI の格子（x0・y1・cell_m）から。"""
    g = meta["grid"]
    cm = float(g["cell_m"])
    x = g["x0"] + (sub["ix"].to_numpy(np.float64) + 0.5) * cm
    y = g["y1"] - (sub["iy"].to_numpy(np.float64) + 0.5) * cm
    dx, dy = np.meshgrid(OFFS * cm, OFFS * cm)
    px = (x[:, None] + dx.ravel()[None, :])
    py = (y[:, None] + dy.ravel()[None, :])
    tr = Transformer.from_crs(CRS.from_proj4(meta["crs_proj4"]), "EPSG:4326", always_xy=True)
    lon, lat = tr.transform(px.ravel(), py.ravel())
    return np.asarray(lon).reshape(px.shape), np.asarray(lat).reshape(px.shape)


def sample(arr: np.ndarray, origin: tuple[float, float], lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    res = 0.00025
    w, n = origin
    c = np.floor((lon - w) / res).astype(np.int64)
    r = np.floor((n - lat) / res).astype(np.int64)
    ok = (c >= 0) & (c < arr.shape[1]) & (r >= 0) & (r < arr.shape[0])
    out = np.full(lon.shape, v49.NODATA, dtype=np.uint8)
    out[ok] = arr[r[ok], c[ok]]
    return out


def main() -> None:
    t0 = time.time()
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    sens = ev.Sensitivity()
    dirs, _, infos = ev.select_aois(gn, sens)
    df = ev.build_table(dirs, seed=20261025, infos=infos, sens=sens, extra_cols=("ix", "iy"))
    parts = []
    for k, d in enumerate(dirs):
        sub = df.loc[df["aoi"] == d.name, ["aoi", "ix", "iy"]].reset_index(drop=True)
        if not len(sub):
            continue
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        lon, lat = cell_points(sub, meta)
        box = (float(lon.min()) - PAD_DEG, float(lat.min()) - PAD_DEG, float(lon.max()) + PAD_DEG,
               float(lat.max()) + PAD_DEG)
        vals = {}
        for layer in ("occurrence", "seasonality", "extent"):
            arr, origin = v49.read_layer(layer, *box)
            vals[layer] = sample(arr, origin, lon, lat)
        ok = (vals["occurrence"] != v49.NODATA) & (vals["seasonality"] != v49.NODATA) & (vals["extent"] != v49.NODATA)
        n_ok = ok.sum(axis=1)
        with np.errstate(invalid="ignore"):
            occ = np.where(ok, vals["occurrence"], 0).sum(axis=1) / n_ok
            seas = np.where(ok, vals["seasonality"], 0).sum(axis=1) / n_ok
            ext = np.where(ok, vals["extent"] == 1, False).sum(axis=1) / n_ok
        parts.append(sub.assign(occ=occ.astype(np.float32), seas=seas.astype(np.float32),
                                maxext=ext.astype(np.float32), n_valid=n_ok.astype(np.int16)))
        if k % 25 == 0:
            print(f"{k + 1}/{len(dirs)} {d.name}: {len(sub)} cells ({time.time() - t0:.0f}s)", flush=True)
    out = pd.concat(parts, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT / "cell_gsw.parquet", index=False)
    meta = {"rule": "amendment B-8 item C4", "source": "JRC Global Surface Water v1.4 (1984-2021)",
            "url_pattern": v49.URL, "points_per_cell": int(len(OFFS) ** 2), "n_cells": len(out),
            "share_with_values": float((out["n_valid"] > 0).mean()),
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sec": round(time.time() - t0)}
    (OUT / "cell_gsw_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
