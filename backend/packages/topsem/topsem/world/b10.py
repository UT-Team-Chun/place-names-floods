"""amendment B-10（初見の査読 第 6 回、事後）の実装。

- C1: 距離帯（500 m 以内・500 m〜1 km・1 km 超）と大陸の指標、名前の特徴との交互作用
- C2: 二段目の L2 の強さ C を、訓練側の場所グループで分けた入れ子の交差検証（Brier）で選ぶ
- C3: 国（セルの最寄りの集落の国）の指標
- C4: 形 AC（形 A の対象に、反転を集落の地面の中で相関させて注入）。
  u_i = Φ(√ρ z_s + √(1−ρ) e_i)（z_s は集落ごと、e_i はセルごとの標準正規）。
  周辺は一様なので反転の確率は形 A と同じ
- C5: JRC の全球河川洪水ハザード図 v2.1.2（CEMS-GloFAS、3 秒角、10° タイル）の
  再現期間 10・100・500 年の浸水深から、各セルの浸水の割合と 100 年の平均の深さ
  （恒久水域の画素は分子・分母から除く）
"""

from __future__ import annotations

import zlib
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.special import ndtr, ndtri
from shapely.geometry import box, shape

from topsem.world.b6 import hashed_uniform
from topsem.world.terrain import GDAL_ENV, _block

BANDS_M = (500.0, 1000.0)  # 距離帯の境（m）: ≤ 500、500〜1000、> 1000
BAND_NAMES = ("b500", "b1k", "bfar")
JRC_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard"
JRC_RPS = (10, 100, 500)
JRC_NODATA = -9999.0
PW_NODATA = 255
C_GRID = tuple(float(x) for x in np.logspace(-4, 2, 13))  # C2: 10^-4〜10^2 の 13 点


# %% C1・C3: 指標と交互作用
def band_indicators(d1: np.ndarray) -> np.ndarray:
    """最寄りの集落までの距離の帯（列: ≤ 500 m、500 m〜1 km、> 1 km）。"""
    d = np.asarray(d1, dtype=np.float64)
    near = d <= BANDS_M[0]
    mid = (d > BANDS_M[0]) & (d <= BANDS_M[1])
    return np.column_stack([near, mid, ~(near | mid)]).astype(np.float32)


def indicators(values: Iterable[str], levels: list[str]) -> np.ndarray:
    """カテゴリの指標（levels の順の列）。levels に無い値は全ての列が 0。"""
    v = np.asarray(list(values), dtype=object)
    return np.column_stack([v == lv for lv in levels]).astype(np.float32) if levels \
        else np.empty((len(v), 0), dtype=np.float32)


def interactions(x: np.ndarray, ind: np.ndarray) -> np.ndarray:
    """特徴 × 指標の全ての組（列の順: 特徴ごとに指標を並べる）。"""
    x = np.asarray(x, dtype=np.float32)
    return (x[:, :, None] * ind[:, None, :]).reshape(len(x), -1)


# %% C2: 二段目の罰則の入れ子の交差検証
def choose_c(a: np.ndarray, y: np.ndarray, inner: np.ndarray, grid: tuple[float, ...] = C_GRID,
             ) -> tuple[float, list[float]]:
    """標準化済みの二段目の入力 a で、内側の fold（inner）の Brier が最小の C。

    inner は場所グループで分けた内側の fold の番号。戻り値: 選んだ C と、格子の各点の
    内側の Brier（全ての内側の評価セルで平均）。
    """
    from sklearn.linear_model import LogisticRegression  # noqa: PLC0415

    folds = [k for k in np.unique(inner) if k >= 0]
    sse = np.zeros(len(grid))
    n = 0
    for k in folds:
        te = inner == k
        tr = ~te
        if tr.sum() == 0 or te.sum() == 0 or len(np.unique(y[tr])) < 2:  # noqa: PLR2004
            continue
        n += int(te.sum())
        for j, c in enumerate(grid):
            m = LogisticRegression(C=c, max_iter=1000)
            m.fit(a[tr], y[tr])
            p = m.predict_proba(a[te])[:, 1]
            sse[j] += float(((p - y[te]) ** 2).sum())
    bs = sse / max(n, 1)
    return float(grid[int(np.argmin(bs))]), [float(v) for v in bs]


# %% C4: 空間にまとまる注入
def settlement_normal(gids: np.ndarray, rep: int) -> np.ndarray:
    """集落ごとの標準正規 z_s（gid と回で決まる。δ によらず同じ）。"""
    salt = zlib.crc32(f"B-10/AC/settlement/{rep}".encode())
    u = hashed_uniform(np.asarray(gids).astype(np.int64).astype(np.uint64), salt)
    return ndtri(np.clip(u, 1e-12, 1 - 1e-12))


def correlated_uniform(u_cell: np.ndarray, gids: np.ndarray, rep: int, rho: float) -> np.ndarray:
    """セルの一様乱数 u_cell を、同じ集落の中で相関させた一様乱数に変える（周辺は一様のまま）。

    集落の無いセル（gid < 0）は u_cell をそのまま返す。
    """
    e = ndtri(np.clip(np.asarray(u_cell, dtype=np.float64), 1e-12, 1 - 1e-12))
    z = settlement_normal(gids, rep)
    out = ndtr(np.sqrt(rho) * z + np.sqrt(1 - rho) * e)
    return np.where(np.asarray(gids) >= 0, out, u_cell)


# %% C5: JRC の洪水ハザード図
def jrc_tiles(index: dict, bounds: tuple[float, float, float, float]) -> list[tuple[int, str]]:
    """tile_extents.geojson（読んだ dict）のうち、経緯度の範囲 bounds にかかるタイル。

    戻り値は (id, name) の並び。
    """
    q = box(*bounds)
    return sorted((int(f["properties"]["id"]), str(f["properties"]["name"]))
                  for f in index["features"] if shape(f["geometry"]).intersects(q))


def jrc_file(tile: tuple[int, str], rp: int | None) -> str:
    """タイルのファイルの相対パス（rp=None は恒久水域）。"""
    i, name = tile
    if rp is None:
        return f"Permanent_WaterBodies/ID{i}_{name}_permanent_water.tif"
    return f"RP{rp}/ID{i}_{name}_RP{rp}_depth.tif"


def _warp_nearest(path: Path, meta: dict, shape_: tuple[int, int],
                  src_nodata: float | None) -> np.ndarray:
    """タイルを AOI の副格子へ最近傍で（タイルの外は NaN）。

    src_nodata が None ならタイルの nodata も NaN に、どの画素にも無い値（_PASS）なら
    nodata を値のまま通す。
    """
    g = meta["grid"]
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    tmp = np.full(shape_, np.nan, dtype=np.float32)
    with rasterio.open(path) as src:
        reproject(source=rasterio.band(src, 1), destination=tmp, dst_transform=dst_tr,
                  dst_crs=CRS.from_proj4(meta["crs_proj4"]), src_nodata=src_nodata,
                  dst_nodata=np.nan, init_dest_nodata=True, resampling=Resampling.nearest)
    return tmp


# 元の nodata を値として通すための、どの画素にも無い値
_PASS = -3.0e38


def warp_jrc(meta: dict, files: dict[int | None, list[Path]]) -> dict[int | None, np.ndarray]:
    """再現期間ごと（と恒久水域 None）のタイルを AOI の副格子へ。タイルの無い画素は NaN。

    深さは元の nodata（−9999、乾いた陸と海）を −9999 のまま残す。恒久水域は 1 を残し、
    nodata（255）は NaN。
    """
    g = meta["grid"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shp = (g["ny"] * sub, g["nx"] * sub)
    out = {}
    with rasterio.Env(**GDAL_ENV):
        for rp, paths in files.items():
            dst = np.full(shp, np.nan, dtype=np.float32)
            for p in paths:
                tmp = _warp_nearest(p, meta, shp, _PASS if rp is not None else None)
                ok = np.isfinite(tmp)
                dst[ok] = tmp[ok]
            out[rp] = dst
    return out


def jrc_features(warped: dict[int | None, np.ndarray], meta: dict) -> dict[str, np.ndarray]:
    """セルの特徴（格子全体、行 0 = 北、ravel 済み）。

    hz_f{rp}: 恒久水域でない画素のうち深さ > 0 の割合。hz_d100: 恒久水域でない画素の
    100 年の平均の深さ（乾いた画素は 0）。hz_cover: 100 年のタイルのある画素の割合。
    恒久水域のタイルが無い画素は恒久水域でないとする。
    """
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    pw = warped.get(None)
    perm = np.zeros_like(warped[JRC_RPS[0]], dtype=bool) if pw is None else (pw == 1)
    out: dict[str, np.ndarray] = {}
    with np.errstate(all="ignore"):
        for rp in JRC_RPS:
            d = warped[rp]
            have = np.isfinite(d) & ~perm
            wet = have & (d > 0)
            nb = _block(have.astype(np.float32), sub, ny, nx).sum(axis=2)
            wb = _block(wet.astype(np.float32), sub, ny, nx).sum(axis=2)
            out[f"hz_f{rp}"] = np.where(nb > 0, wb / nb, np.nan).astype(np.float32)
            if rp == 100:  # noqa: PLR2004
                depth = np.where(wet, d, 0.0).astype(np.float32)
                db = _block(np.where(have, depth, 0.0), sub, ny, nx).sum(axis=2)
                out["hz_d100"] = np.where(nb > 0, db / nb, np.nan).astype(np.float32)
                cov = _block(np.isfinite(d).astype(np.float32), sub, ny, nx).mean(axis=2)
                out["hz_cover"] = cov.astype(np.float32)
    return {k: v.ravel() for k, v in out.items()}


def fill_hazard(df: pd.DataFrame, cols: list[str]) -> dict:
    """ハザードの特徴の欠け（タイルの無いセル・全て恒久水域のセル）を 0 にし、割合を返す。"""
    info = {}
    for c in cols:
        v = df[c].to_numpy(np.float64)
        miss = ~np.isfinite(v)
        info[f"share_missing_{c}"] = float(miss.mean())
        df[c] = np.where(miss, 0.0, v).astype(np.float32)
    return info
