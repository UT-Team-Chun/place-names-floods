"""HydroRIVERS v1.0 → AOI 格子セルの河道距離（m）。

HydroSHEDS（15 秒格子由来、上流面積 ≥ 10 km² の河道）を全球一様の水文共変量として使う。
EMS 製品同梱の hydrography 層は製品ごとに取得元・密度が違うので使わない。

セル中心から (a) 全河道、(b) 上流面積 ≥ 1000 km² の主要河道 までの距離。
AOI 外の河道も効くように、格子を 20 km 余白で広げた 100 m 格子に線を焼いて
ユークリッド距離変換をかけ、セル中心で読む（余白を超える距離は余白長で打ち切る）。
"""

from __future__ import annotations

import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import shapely
from pyogrio import raw
from pyproj import CRS, Transformer
from rasterio import features
from rasterio.transform import from_origin
from scipy.ndimage import distance_transform_edt

REGIONS = ("af", "ar", "as", "au", "eu", "gr", "na", "sa", "si")
URL = "https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_{r}_shp.zip"
PAD_M = 20_000.0
RES_M = 100.0
MAJOR_UPLAND_KM2 = 1000.0
UA = {"User-Agent": "toponym-semantics research (ryoutaokauchi@g.ecc.u-tokyo.ac.jp)"}


def ensure_downloaded(root: Path) -> list[Path]:
    """全地域の zip を取得（キャッシュ）。"""
    out = []
    for r in REGIONS:
        dst = root / f"HydroRIVERS_v10_{r}_shp.zip"
        if not dst.is_file():
            dst.parent.mkdir(parents=True, exist_ok=True)
            tmp = dst.with_suffix(".part")
            # data.hydrosheds.org は既定の "Python-urllib" UA に 403 を返す（2026-10-25 実測）
            req = urllib.request.Request(URL.format(r=r), headers=UA)  # noqa: S310
            with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as f:  # noqa: S310
                while chunk := resp.read(1 << 20):
                    f.write(chunk)
            tmp.rename(dst)
        out.append(dst)
    return out


def _shp_member(zip_path: Path) -> str:
    with zipfile.ZipFile(zip_path) as zf:
        return next(n for n in zf.namelist() if n.lower().endswith(".shp"))


def rivers_in_bounds(zips: list[Path], bounds: tuple[float, float, float, float],
                     ) -> tuple[np.ndarray, np.ndarray]:
    """lon/lat 範囲の河道線（WGS84 shapely 配列）と上流面積（km²）。"""
    geoms, upland = [], []
    for z in zips:
        vsi = f"/vsizip/{z}/{_shp_member(z)}"
        _meta, _, g, f = raw.read(vsi, bbox=bounds, columns=["UPLAND_SKM"])
        if g is None or len(g) == 0:
            continue
        geoms.append(shapely.from_wkb(g))
        upland.append(np.asarray(f[0], dtype=float))
    if not geoms:
        return np.array([], dtype=object), np.array([], dtype=float)
    return np.concatenate(geoms), np.concatenate(upland)


def river_distances(meta: dict, zips: list[Path], bounds_ll: tuple[float, float, float, float],
                    ) -> dict[str, np.ndarray]:
    """AOI 格子全セル（ny*nx、行 0 = 北）の河道距離。"""
    g = meta["grid"]
    crs = CRS.from_proj4(meta["crs_proj4"])
    w, s, e, n = bounds_ll
    pad_deg = PAD_M / 111_000.0 * 1.5
    lines, upland = rivers_in_bounds(zips, (w - pad_deg, s - pad_deg, e + pad_deg, n + pad_deg))
    x0, y1 = g["x0"] - PAD_M, g["y1"] + PAD_M
    width = g["nx"] * g["cell_m"] + 2 * PAD_M
    height = g["ny"] * g["cell_m"] + 2 * PAD_M
    shape = (int(np.ceil(height / RES_M)), int(np.ceil(width / RES_M)))
    tr = from_origin(x0, y1, RES_M, RES_M)
    to_local = Transformer.from_crs(4326, crs, always_xy=True)
    # セル中心 → 100 m 格子の行列
    ix, iy = np.meshgrid(np.arange(g["nx"]), np.arange(g["ny"]))
    cx = g["x0"] + (ix.ravel() + 0.5) * g["cell_m"]
    cy = g["y1"] - (iy.ravel() + 0.5) * g["cell_m"]
    col = np.clip(((cx - x0) / RES_M).astype(int), 0, shape[1] - 1)
    row = np.clip(((y1 - cy) / RES_M).astype(int), 0, shape[0] - 1)
    out: dict[str, np.ndarray] = {}
    for tag, thr in (("all", 0.0), ("major", MAJOR_UPLAND_KM2)):
        sel = lines[upland >= thr] if len(lines) else lines
        if len(sel) == 0:
            out[f"dist_river_{tag}"] = np.full(len(cx), PAD_M, dtype=np.float32)
            continue
        loc = shapely.transform(
            sel, lambda xy: np.column_stack(to_local.transform(xy[:, 0], xy[:, 1])))
        burn = features.rasterize(((geom, 1) for geom in loc), out_shape=shape, transform=tr,
                                  fill=0, dtype="uint8", all_touched=True)
        dist = distance_transform_edt(burn == 0) * RES_M
        out[f"dist_river_{tag}"] = np.minimum(dist[row, col], PAD_M).astype(np.float32)
    return out
