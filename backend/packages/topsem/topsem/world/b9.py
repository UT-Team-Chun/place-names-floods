"""amendment B-9（初見の査読 第 5 回、事後）の実装。

- C3: 形 AH の集落ごとの重み w_s ~ 指数分布（平均 1）。集落の gid と回の塩で
  決まる（δ によらず同じ重み）
- C5: 各セルの JRC GSW の occurrence の平均（v58 の 25 点）が閾値以上のセルを
  永久水域として除く
- C7: FABDEM v1-2（建物・樹冠を除いた GLO-30、1° タイル、配布は 10° の zip）から、
  GLO-30 と同じ標高の特徴と、集水 100・1000 km² 以上の HydroRIVERS の河道の
  最寄りの頂点からの高さ（B-8 の B7 と同じ作り方）
"""

from __future__ import annotations

import math
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS, Transformer
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.ndimage import minimum_filter
from scipy.spatial import cKDTree

from topsem.world.b6 import hashed_uniform
from topsem.world.terrain import GDAL_ENV, _block, _box_mean, _neigh_min, grid_lonlat_bounds

FABDEM_BASE = "https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn"
FABDEM_NODATA = -9999.0
UPLAND_CUTS = (100.0, 1000.0)  # B-8 の B7 と同じ（km²）
RIVER_PAD_DEG = 0.2  # 河道の頂点を探す範囲（AOI ± 0.2°、B-8 の B7 と同じ）
RIVER_WIN_PX = 9  # 河道の標高: 頂点の周り 9 × 9 画素（約 270 m）の最小
EARTH_M = 6_371_008.8
LON_WRAP = 180  # 経度 180 の名前（W180）


# %% C3: 集落ごとの重み
def settlement_weights(gids: np.ndarray, rep: int) -> np.ndarray:
    """形 AH の重み w_s = −log(1 − u_s)（指数分布、平均 1）。

    u_s は集落の gid と回で決まる一様乱数。
    """
    salt = zlib.crc32(f"B-9/AH/weights/{rep}".encode())
    u = hashed_uniform(np.asarray(gids).astype(np.int64).astype(np.uint64), salt)
    return -np.log1p(-np.minimum(u, 1 - 1e-12))


# %% C5: 共通の永久水域
def gsw_perm_mask(df: pd.DataFrame, gsw: pd.DataFrame, threshold: float) -> np.ndarray:
    """セル（aoi, ix, iy）の GSW の occurrence の平均（%）が閾値以上なら True。

    GSW の無いセルは False。
    """
    key = ["aoi", "ix", "iy"]
    m = df[key].merge(gsw[[*key, "occ"]], on=key, how="left")
    occ = m["occ"].to_numpy(np.float64)
    return np.isfinite(occ) & (np.nan_to_num(occ, nan=-1.0) >= threshold)


# %% C7: FABDEM
def _ns(v: int) -> str:
    return f"N{v:02d}" if v >= 0 else f"S{-v:02d}"


def _ew(v: int) -> str:
    if v >= LON_WRAP:  # 公開の名前の規則（例: N00E170-N10W180）
        return "W180"
    return f"E{v:03d}" if v >= 0 else f"W{-v:03d}"


def fabdem_zip_name(lat: int, lon: int) -> str:
    """1° タイル（南西角 lat, lon）を含む 10° の zip の名前。"""
    la, lo = math.floor(lat / 10) * 10, math.floor(lon / 10) * 10
    return f"{_ns(la)}{_ew(lo)}-{_ns(la + 10)}{_ew(lo + 10)}_FABDEM_V1-2.zip"


def fabdem_tile_name(lat: int, lon: int) -> str:
    return f"{_ns(lat)}{_ew(lon)}_FABDEM_V1-2.tif"


def fabdem_path(lat: int, lon: int, cache: Path | None) -> str:
    """GDAL のパス。

    cache（NAS の fabdem/）に取り出したタイル（tiles/）があればそれ、
    無いと記録したタイル（tiles/*.absent、海域）は存在しないパス（読むと RasterioIOError）、
    10° の zip（zips/）があれば /vsizip/、どれも無ければ /vsizip//vsicurl/。
    """
    z, t = fabdem_zip_name(lat, lon), fabdem_tile_name(lat, lon)
    if cache is not None:
        tile = cache / "tiles" / t
        if tile.is_file() or (cache / "tiles" / f"{t}.absent").is_file():
            return str(tile)
        if (cache / "zips" / z).is_file():
            return f"/vsizip/{cache / 'zips' / z}/{t}"
    return f"/vsizip//vsicurl/{FABDEM_BASE}/{z}/{t}"


def tiles_for(bounds: tuple[float, float, float, float]) -> list[tuple[int, int]]:
    w, s, e, n = bounds
    return [(la, lo) for la in range(math.floor(s), math.floor(n) + 1)
            for lo in range(math.floor(w), math.floor(e) + 1)]


def _env() -> rasterio.Env:
    return rasterio.Env(**(GDAL_ENV | {"CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".zip,.tif"}))


def warp_fabdem(meta: dict, cache: Path | None) -> tuple[np.ndarray, int, int]:
    """FABDEM を AOI の副格子へ（GLO-30 と同じ双一次）。

    戻り値: 配列、読めたタイルの数、無かったタイルの数。
    """
    g = meta["grid"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shape = (g["ny"] * sub, g["nx"] * sub)
    dst = np.full(shape, np.nan, dtype=np.float32)
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    dst_crs = CRS.from_proj4(meta["crs_proj4"])
    ok_n, miss = 0, 0
    with _env():
        for lat, lon in tiles_for(grid_lonlat_bounds(meta)):
            try:
                src = rasterio.open(fabdem_path(lat, lon, cache))
            except rasterio.errors.RasterioIOError:
                miss += 1  # 海域などタイル無し
                continue
            with src:
                tmp = np.full(shape, np.nan, dtype=np.float32)
                reproject(source=rasterio.band(src, 1), destination=tmp, dst_transform=dst_tr,
                          dst_crs=dst_crs, src_nodata=FABDEM_NODATA, dst_nodata=np.nan,
                          resampling=Resampling.bilinear)
                ok = np.isfinite(tmp)
                dst[ok] = tmp[ok]
                ok_n += 1
    return dst, ok_n, miss


def fabdem_elevation_features(dem: np.ndarray, meta: dict) -> dict[str, np.ndarray]:
    """GLO-30 の `terrain_features` と同じ定義の標高の特徴（格子全体、行 0 = 北、ravel 済み）。"""
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    res = g["fine_res_m"]
    gy, gx = np.gradient(np.where(np.isfinite(dem), dem, np.nan), res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy))).astype(np.float32)
    with np.errstate(all="ignore"):
        db = _block(dem, sub, ny, nx)
        e_mean = np.nanmean(db, axis=2)
        out = {"elev_1km": _box_mean(e_mean, 2), "elev_mean": e_mean.astype(np.float32),
               "depress": (e_mean - _neigh_min(np.where(np.isfinite(e_mean), e_mean, np.inf)))
               .astype(np.float32),
               "elev_min": np.nanmin(db, axis=2).astype(np.float32),
               "elev_std": np.nanstd(db, axis=2).astype(np.float32),
               "slope_mean": np.nanmean(_block(slope, sub, ny, nx), axis=2).astype(np.float32),
               "dem_valid_frac": np.mean(np.isfinite(db), axis=2).astype(np.float32)}
    return {k: v.ravel() for k, v in out.items()}


def _unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def river_elevations(vert: pd.DataFrame, cache: Path | None) -> np.ndarray:
    """河道の頂点の標高: FABDEM の頂点の周り RIVER_WIN_PX 四方の画素の最小（無ければ NaN）。"""
    lon, lat = vert["lon"].to_numpy(np.float64), vert["lat"].to_numpy(np.float64)
    out = np.full(len(vert), np.nan)
    la0, lo0 = np.floor(lat).astype(int), np.floor(lon).astype(int)
    with _env():
        for la, lo in sorted(set(zip(la0.tolist(), lo0.tolist(), strict=True))):
            idx = np.flatnonzero((la0 == la) & (lo0 == lo))
            try:
                src = rasterio.open(fabdem_path(la, lo, cache))
            except rasterio.errors.RasterioIOError:
                continue
            with src:
                a = src.read(1).astype(np.float32)
                a[a == FABDEM_NODATA] = np.nan
                m = minimum_filter(np.where(np.isfinite(a), a, np.inf), size=RIVER_WIN_PX,
                                   mode="nearest")
                r, c = rasterio.transform.rowcol(src.transform, lon[idx], lat[idx])
                r = np.clip(np.asarray(r), 0, a.shape[0] - 1)
                c = np.clip(np.asarray(c), 0, a.shape[1] - 1)
                v = m[r, c]
                out[idx] = np.where(np.isfinite(v), v, np.nan)
    return out


def cell_lonlat(meta: dict) -> tuple[np.ndarray, np.ndarray]:
    """格子全体（ny*nx、行 0 = 北）のセルの中心の経緯度。"""
    g = meta["grid"]
    cm = float(g["cell_m"])
    xs = g["x0"] + (np.arange(g["nx"]) + 0.5) * cm
    ys = g["y1"] - (np.arange(g["ny"]) + 0.5) * cm
    xx, yy = np.meshgrid(xs, ys)
    tr = Transformer.from_crs(CRS.from_proj4(meta["crs_proj4"]), 4326, always_xy=True)
    lon, lat = tr.transform(xx.ravel(), yy.ravel())
    return np.asarray(lon), np.asarray(lat)


def river_heights(elev_mean: np.ndarray, meta: dict, vert: pd.DataFrame,
                  cache: Path | None) -> dict[str, np.ndarray]:
    """各セルの、集水 cut km² 以上の最寄りの河道の頂点からの高さと距離。

    高さ = セルの平均標高 − 頂点の標高。

    vert は AOI の範囲 ± RIVER_PAD_DEG の頂点（lon, lat, upland）。
    """
    lon, lat = cell_lonlat(meta)
    xyz = _unit_xyz(lon, lat)
    out: dict[str, np.ndarray] = {}
    for cut in UPLAND_CUTS:
        sub = vert.loc[vert["upland"] >= cut].reset_index(drop=True)
        name = f"h_above_{int(cut)}"
        if sub.empty:
            out[name] = np.full(len(lon), np.nan, dtype=np.float32)
            out[f"dist_{int(cut)}_m"] = np.full(len(lon), np.nan, dtype=np.float32)
            continue
        el_r = river_elevations(sub, cache)
        ok = np.isfinite(el_r)
        if not ok.any():
            out[name] = np.full(len(lon), np.nan, dtype=np.float32)
            out[f"dist_{int(cut)}_m"] = np.full(len(lon), np.nan, dtype=np.float32)
            continue
        sub, el_r = sub.loc[ok].reset_index(drop=True), el_r[ok]
        tree = cKDTree(_unit_xyz(sub["lon"].to_numpy(), sub["lat"].to_numpy()))
        d, j = tree.query(xyz)
        out[name] = (elev_mean - el_r[j]).astype(np.float32)
        out[f"dist_{int(cut)}_m"] = (d * EARTH_M).astype(np.float32)
    return out
