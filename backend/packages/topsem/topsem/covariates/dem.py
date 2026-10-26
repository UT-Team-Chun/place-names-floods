"""SRTM 30 m 標高（terrarium png・AWS Open Data）の取得と格子サンプリング。

版・出所（docs/data/licenses.md B2 に対応）:
- タイル: https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png
- 変換: height = R*256 + G + B/256 - 32768（terrarium 符号付きエンコーディング）
- 内容: NASA SRTM 1 秒（30 m）由来の公開 terrain tiles。
- 時点: SRTM は 2000 年取得が基本。以降の盛土・河道改修は未反映
  （時点台帳 B2 にその旨記録済み）。
- 主解析では国土地理院 DEM 5A に置換する（キャッシュ方式は同一、取得に要登録）。
"""

from __future__ import annotations

import math
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from topsem.config import COVARIATES_DIR

TILE_ROOT = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium"
Z = 13
TILE_PX = 256
_N_TILE = 1 << Z


def tile_xy(lon: float, lat: float, z: int = Z) -> tuple[int, int]:
    """Web Mercator のタイル番号 (x, y)。"""
    n = 1 << z
    x = math.floor((lon + 180.0) / 360.0 * n)
    lat_r = math.radians(lat)
    y = math.floor((1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r)) / math.pi) / 2.0 * n)
    return x, y


def _fractile(lon: float, lat: float) -> tuple[float, int, float, int]:
    """(世界 x 小数, タイル x, 世界 y 小数, タイル y)。"""
    fx = (lon + 180.0) / 360.0 * _N_TILE
    lat_r = math.radians(lat)
    fy = (1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r)) / math.pi) / 2.0 * _N_TILE
    return fx, math.floor(fx), fy, math.floor(fy)


def _tile(x: int, y: int, cache_dir: Path) -> np.ndarray:
    """タイル標高 [256,256]（行 0 = 北）。無効タイルは 0 埋め。"""
    p = cache_dir / str(x) / f"{y}.npy"
    if p.is_file():
        return np.load(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    png = cache_dir / str(x) / f"{y}.png"
    ok = True
    try:
        urlretrieve(f"{TILE_ROOT}/{Z}/{x}/{y}.png", png)  # noqa: S310
        arr = _decode_png_height(png).astype(np.float32)
    except Exception:  # noqa: BLE001
        # 失敗タイルを 0 埋めで**キャッシュしない**（恒久的な低地偽装は M1 を汚し、
        # 地名特徴がその偽低地を「当てる」見かけの利得を作る — レビュー C-4）
        arr = np.zeros((TILE_PX, TILE_PX), dtype=np.float32)
        ok = False
    finally:
        png.unlink(missing_ok=True)
    if ok:
        np.save(p, arr)
    return arr


def _decode_png_height(png: Path) -> np.ndarray:
    # plots 追加依存をコアから遅延させる（CLAUDE.md の import 規則の明示的例外・要約に記録）
    from matplotlib.image import imread  # noqa: PLC0415

    img = (imread(str(png)) * 255.0).astype(np.uint16)
    return (
        img[:, :, 0].astype(np.float64) * 256.0
        + img[:, :, 1].astype(np.float64)
        + img[:, :, 2].astype(np.float64) / 256.0
        - 32768.0
    ).astype(np.float32)


def sample_grid(
    ix: np.ndarray,
    iy: np.ndarray,
    cell_lon_deg: float,
    cell_lat_deg: float,
    cache_dir: Path | None = None,
) -> np.ndarray:
    """格子中心の標高を最近傍サンプリング。海・欠測は 0。

    `labels.grid` と同一の「地域基準緯度で一定経度幅」の格子（cell_*_deg を
    その基準緯度で計算して渡す）を仮定する。
    """
    cache = cache_dir or COVARIATES_DIR / "terrarium" / str(Z)
    lon = (np.asarray(ix, dtype=float) + 0.5) * cell_lon_deg
    lat = (np.asarray(iy, dtype=float) + 0.5) * cell_lat_deg
    out = np.zeros(len(lon), dtype=np.float32)
    loaded: dict[tuple[int, int], np.ndarray] = {}
    for j, (lo, la) in enumerate(zip(lon, lat, strict=True)):
        ffx, tx, ffy, ty = _fractile(float(lo), float(la))
        if (tx, ty) not in loaded:
            loaded[(tx, ty)] = _tile(tx, ty, cache)
        arr = loaded[(tx, ty)]
        px = min(TILE_PX - 1, max(0, int((ffx - tx) * TILE_PX)))
        py = min(TILE_PX - 1, max(0, int((ffy - ty) * TILE_PX)))
        v = arr[py, px]  # png 行 0 = 北端。mercator fy も北端基準で反転不要
        out[j] = 0.0 if not np.isfinite(v) else float(v)
    return out


def relief_depression(elev_grid: np.ndarray) -> np.ndarray:
    """3x3 近傍最小との差（>= 0 の相対低地度）。grid は [iy, ix] 2 次元。"""
    g = np.pad(elev_grid, 1, constant_values=np.inf)
    m = np.full_like(elev_grid, np.inf)
    for dy in range(3):
        for dx in range(3):
            m = np.minimum(
                m,
                g[dy : dy + elev_grid.shape[0], dx : dx + elev_grid.shape[1]],
            )
    m = np.where(np.isfinite(m), m, elev_grid)
    return (elev_grid - m).astype(np.float32)
