"""GLO-30 DEM・ASF GLO-30 HAND・ESA WorldCover を AOI の副格子へ投影し 250 m に集約。

ラベル（`topsem.world.ems`）と同じ LAEA 副格子（25 m、大 AOI は 50 m）に載せる。
出所（AWS Open Data、匿名 COG 読み出し）:

- GLO-30: ``copernicus-dem-30m`` — DSM（建物・樹冠込み）。1°タイル、名前は南西角
- HAND: ``glo-30-hand`` v1/2021 — GLO-30 から ASF が算出した排水網からの比高
- WorldCover: ``esa-worldcover`` v200/2021 — 10 m 土地被覆、3°タイル

地形情報の梯子（protocol v2 の RQ3）に要る列をここで一度に出す:
T1 = 1 km 平均標高、T2 = 250 m 平均標高 + 相対低地度（paper_1 と同定義）、
T3 = セル内 30 m 統計（最小・標準偏差・傾斜）、T4 = HAND。

WorldCover は建物・農地・樹林・草地・裸地の被覆率のみ特徴に使う。水域（80）・
湿地（90）は洪水年の撮像が混入し得るので**特徴にしない**（QC 用に保存のみ）。
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import rasterio
from pyproj import CRS, Transformer
from rasterio.enums import Resampling
from rasterio.warp import reproject

GLO30 = (
    "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)
HAND = (
    "/vsicurl/https://glo-30-hand.s3.amazonaws.com/v1/2021/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_HAND.tif"
)
WORLDCOVER = (
    "/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
    "ESA_WorldCover_10m_2021_v200_{ns}{lat:02d}{ew}{lon:03d}_Map.tif"
)
GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "GDAL_HTTP_MULTIRANGE": "YES",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "VSI_CACHE": "TRUE",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "2",
}
WC_FEATURE_CLASSES = {"built": 50, "crop": 40, "tree": 10, "shrub": 20, "grass": 30, "bare": 60}
WC_QC_CLASSES = {"water": 80, "wetland": 90}


def _tiles(bounds: tuple[float, float, float, float], step: int) -> list[tuple[int, int]]:
    """lon/lat 範囲を覆うタイルの南西角 (lat, lon)。"""
    w, s, e, n = bounds
    lats = range(math.floor(s / step) * step, math.floor(n / step) * step + step, step)
    lons = range(math.floor(w / step) * step, math.floor(e / step) * step + step, step)
    return [(la, lo) for la in lats for lo in lons]


def _fmt(tpl: str, lat: int, lon: int) -> str:
    return tpl.format(ns="N" if lat >= 0 else "S", lat=abs(lat),
                      ew="E" if lon >= 0 else "W", lon=abs(lon))


def grid_lonlat_bounds(meta: dict, pad_m: float = 0.0) -> tuple[float, float, float, float]:
    """AOI 格子（LAEA）の外接 lon/lat 範囲（辺上を密に変換して取る）。"""
    g = meta["grid"]
    x0, y1 = g["x0"] - pad_m, g["y1"] + pad_m
    x1 = g["x0"] + g["nx"] * g["cell_m"] + pad_m
    y0 = g["y1"] - g["ny"] * g["cell_m"] - pad_m
    t = np.linspace(0, 1, 101)
    xs = np.concatenate([x0 + (x1 - x0) * t, np.full_like(t, x1), x1 - (x1 - x0) * t,
                         np.full_like(t, x0)])
    ys = np.concatenate([np.full_like(t, y1), y1 - (y1 - y0) * t, np.full_like(t, y0),
                         y0 + (y1 - y0) * t])
    tr = Transformer.from_crs(CRS.from_proj4(meta["crs_proj4"]), 4326, always_xy=True)
    lon, lat = tr.transform(xs, ys)
    return float(np.min(lon)), float(np.min(lat)), float(np.max(lon)), float(np.max(lat))


def warp_to_fine(  # noqa: PLR0913 — タイル様式・格子・再標本化・型・nodata は独立に指定する
    tpl: str, step: int, meta: dict, resampling: Resampling, *,
    dtype: str = "float32", nodata: float = np.nan,
) -> np.ndarray:
    """タイル群を AOI の副格子（ny*sub, nx*sub）へ再投影（無いタイルは nodata のまま）。"""
    g = meta["grid"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shape = (g["ny"] * sub, g["nx"] * sub)
    dst = np.full(shape, nodata, dtype=dtype)
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    dst_crs = CRS.from_proj4(meta["crs_proj4"])
    with rasterio.Env(**GDAL_ENV):
        for lat, lon in _tiles(grid_lonlat_bounds(meta), step):
            url = _fmt(tpl, lat, lon)
            try:
                src = rasterio.open(url)
            except rasterio.errors.RasterioIOError:
                continue  # 海域などタイル無し
            with src:
                tmp = np.full(shape, nodata, dtype=dtype)
                reproject(
                    source=rasterio.band(src, 1),
                    destination=tmp,
                    dst_transform=dst_tr,
                    dst_crs=dst_crs,
                    dst_nodata=nodata,
                    resampling=resampling,
                )
                is_float = np.issubdtype(np.dtype(dtype), np.floating)
                ok = ~np.isnan(tmp) if is_float else tmp != nodata
                dst[ok] = tmp[ok]
    return dst


def _block(a: np.ndarray, sub: int, ny: int, nx: int) -> np.ndarray:
    return a.reshape(ny, sub, nx, sub).transpose(0, 2, 1, 3).reshape(ny, nx, sub * sub)


def _neigh_min(a: np.ndarray) -> np.ndarray:
    g = np.pad(a, 1, constant_values=np.inf)
    m = np.full_like(a, np.inf)
    for dy in range(3):
        for dx in range(3):
            m = np.minimum(m, g[dy:dy + a.shape[0], dx:dx + a.shape[1]])
    return np.where(np.isfinite(m), m, a)


def _box_mean(a: np.ndarray, r: int) -> np.ndarray:
    """(2r+1)² 近傍の nan 無視平均（1 km 平均標高 = r=2 で 5×5 セル ≈ 1.25 km）。"""
    v = np.where(np.isfinite(a), a, 0.0)
    w = np.isfinite(a).astype(np.float64)
    pv = np.pad(v, r)
    pw = np.pad(w, r)
    s = np.zeros_like(v, dtype=np.float64)
    c = np.zeros_like(v, dtype=np.float64)
    for dy in range(2 * r + 1):
        for dx in range(2 * r + 1):
            s += pv[dy:dy + a.shape[0], dx:dx + a.shape[1]]
            c += pw[dy:dy + a.shape[0], dx:dx + a.shape[1]]
    with np.errstate(invalid="ignore", divide="ignore"):
        return (s / c).astype(np.float32)


def terrain_features(meta: dict) -> dict[str, np.ndarray]:
    """1 AOI の全セル（格子全体 ny*nx、行 0 = 北）の地形・土地被覆特徴。"""
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    res = g["fine_res_m"]

    dem = warp_to_fine(GLO30, 1, meta, Resampling.bilinear)
    hand = warp_to_fine(HAND, 1, meta, Resampling.bilinear)
    wc = warp_to_fine(WORLDCOVER, 3, meta, Resampling.mode, dtype="uint8", nodata=0)

    gy, gx = np.gradient(np.where(np.isfinite(dem), dem, np.nan), res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy))).astype(np.float32)

    with np.errstate(all="ignore"):
        db = _block(dem, sub, ny, nx)
        e_mean = np.nanmean(db, axis=2)
        e_min = np.nanmin(db, axis=2)
        e_std = np.nanstd(db, axis=2)
        s_mean = np.nanmean(_block(slope, sub, ny, nx), axis=2)
        hb = _block(hand, sub, ny, nx)
        h_mean = np.nanmean(hb, axis=2)
        h_min = np.nanmin(hb, axis=2)
        dem_valid = np.mean(np.isfinite(db), axis=2)
    out: dict[str, np.ndarray] = {
        "elev_1km": _box_mean(e_mean, 2),
        "elev_mean": e_mean.astype(np.float32),
        "depress": (e_mean - _neigh_min(np.where(np.isfinite(e_mean), e_mean, np.inf)))
        .astype(np.float32),
        "elev_min": e_min.astype(np.float32),
        "elev_std": e_std.astype(np.float32),
        "slope_mean": s_mean.astype(np.float32),
        "hand_mean": h_mean.astype(np.float32),
        "hand_min": h_min.astype(np.float32),
        "dem_valid_frac": dem_valid.astype(np.float32),
    }
    wb = _block(wc, sub, ny, nx)
    for k, code in {**WC_FEATURE_CLASSES, **WC_QC_CLASSES}.items():
        out[f"lc_{k}"] = np.mean(wb == code, axis=2).astype(np.float32)
    out["lc_valid_frac"] = np.mean(wb != 0, axis=2).astype(np.float32)
    return {k: v.ravel() for k, v in out.items()}


def provenance() -> dict[str, Any]:
    return {
        "glo30": GLO30.replace("/vsicurl/", ""),
        "hand": HAND.replace("/vsicurl/", ""),
        "worldcover": WORLDCOVER.replace("/vsicurl/", ""),
        "worldcover_feature_classes": WC_FEATURE_CLASSES,
        "worldcover_qc_only": WC_QC_CLASSES,
    }
