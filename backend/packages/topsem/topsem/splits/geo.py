"""窓の地理的排他（split v5・amendment A-7 / r4' F3）。

watershed-loro-4 まで watershed_holdout は「水系キー」でしか両側を分けて
おらず、隣接水系の窓が物理地面を共有していた（久慈×那珂 14.7×4.4 km、
久慈 wet 113 中 92 の中心が那珂訓練窓内）。CLAUDE.md 漏洩規則 4 の
「同一地点・隣接セルを両側に置かない」を幾何面で達成する。

格納 ix,iy は unit の lat_ref に依存するローカル指標なので、**共通メート
ル座標**へ変換してから近接判定する（固定基準緯度 36.5° の等距離円筒）。
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

_DEG_PER_M = 1.0 / 111_320.0
_REF_LAT = 36.5  # プロジェクト中央緯度（全窓で同一値を使うこと。意味は近似のみ）
_SIZE_M = 250.0


def cell_meters(
    ix: np.ndarray, iy: np.ndarray, lat_ref: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """ix,iy（unit 固有）→ 共通メートル座標。cell_lon は metadata の lat_ref で復元する。"""
    cell_lon = _SIZE_M * _DEG_PER_M / np.cos(np.radians(lat_ref))
    lon = (np.asarray(ix, dtype=float) + 0.5) * cell_lon
    lat = (np.asarray(iy, dtype=float) + 0.5) * _SIZE_M * _DEG_PER_M
    xm = lon * 111_320.0 * math.cos(math.radians(_REF_LAT))
    ym = lat * 111_320.0
    return xm, ym


def exclude_near_holdout(
    tr: pd.DataFrame,
    holdout: pd.DataFrame,
    buffer_m: float,
) -> tuple[pd.DataFrame, int]:
    """訓練側から、holdout セル中心から buffer_m 以内のセルを除去する。

    返り値: (縮小した tr, 除去件数)。tr/holdout はいずれも ix,iy,lat_ref 列を
    持つ DataFrame（build_dataset 由来）。
    """
    if len(holdout) == 0 or buffer_m <= 0:
        return tr, 0
    hx, hy = cell_meters(
        holdout["ix"].to_numpy(), holdout["iy"].to_numpy(), holdout["lat_ref"].to_numpy()
    )
    tx, ty = cell_meters(tr["ix"].to_numpy(), tr["iy"].to_numpy(), tr["lat_ref"].to_numpy())
    tree = cKDTree(np.column_stack([hx, hy]))
    d, _ = tree.query(np.column_stack([tx, ty]), k=1)
    keep = d > buffer_m
    return tr[np.asarray(keep)], int((~keep).sum())
