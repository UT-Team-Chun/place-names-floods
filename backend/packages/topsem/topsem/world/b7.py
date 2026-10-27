"""事後 B-7（初見の査読 第 3 回への対応）の部品。

計算と規則は amendment B-7 で計算の前に固定したもの。

- 1: 置換ごとに地域構成（AOI 平均）を置き換えた名前から計算し直す・AOI 内で中心化する
- 4: 別の EMS の事象で同じ場所が浸水していた割合（繰り返し）と、集落の点から 500 m 以内の平均
- 5: 氾濫原の地面を指す読み（凍結した辞書の河川・湿地の区分だけ）
- 6: 部分的に浸水したセル（y = 0 で浸水被覆が 0.05 を超える）
"""

from __future__ import annotations

from datetime import UTC, datetime

import numpy as np
from pyproj import Transformer

FLOODPLAIN_CATEGORIES = ("riverine", "wet")  # 5
RECUR_MIN_DAYS = 30.0  # 4: これ以下の日付の差の事象は同じ洪水とみなして使わない
PARTIAL_COVER = 0.05  # 6
PLASMODE_REFIT_DELTAS = (0.5, 1.0, 1.5)  # 2: OR 1.65・2.72・4.48（B-6 の 2 の格子の上の点）
PLASMODE_REFIT_REPS = 3


# %% 1: 地域構成を計算し直す・中心化
def aoi_means(x: np.ndarray, aoi: np.ndarray) -> np.ndarray:
    """各列の AOI 内の平均（セルで重みづけ）を、各セルへ配った配列（x と同じ形）。"""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        return aoi_means(x[:, None], aoi)[:, 0]
    _, inv = np.unique(np.asarray(aoi), return_inverse=True)
    sums = np.zeros((inv.max() + 1, x.shape[1]))
    np.add.at(sums, inv, x)
    cnt = np.bincount(inv).astype(np.float64)
    return (sums / cnt[:, None])[inv]


def centre_within(x: np.ndarray, aoi: np.ndarray) -> np.ndarray:
    """AOI 内で中心化した特徴（x − AOI 平均）。"""
    return np.asarray(x, dtype=np.float64) - aoi_means(x, aoi)


def unit_sse(y: np.ndarray, p: np.ndarray, unit_codes: np.ndarray, n_units: int) -> np.ndarray:
    """場所ごとの二乗誤差の和（置換の分布を場所の単位で残す）。"""
    return np.bincount(unit_codes, weights=(np.asarray(p, np.float64) - y) ** 2, minlength=n_units)


# %% 4: 繰り返す浸水
def _utc(t: str) -> datetime:
    d = datetime.fromisoformat(t)
    return d.replace(tzinfo=UTC) if d.tzinfo is None else d.astimezone(UTC)


def days_apart(t1: str, t2: str) -> float:
    """二つの事象の日付（ISO の文字列、時差の無いものは UTC）の差の日数。読めなければ無限大。"""
    try:
        a, b = _utc(t1), _utc(t2)
    except ValueError:
        return float("inf")
    return abs((a - b).total_seconds()) / 86400.0


def recurrence(lonlat: np.ndarray,
               others: list[tuple[str, dict, np.ndarray, np.ndarray, np.ndarray]],
               ) -> tuple[np.ndarray, np.ndarray]:
    """別の事象で浸水していた割合と、そのセルを評価した事象の数。

    lonlat: 対象のセルの経緯度（n × 2）。others: [(crs_proj4, grid, ix, iy, y)]、y ∈ {0, 1} は
    別の事象の評価できてラベルのあるセル（確証と同じ規則）。セルの中心を別の AOI の格子に写し、
    その格子のセルのラベルを読む。評価した事象が無いセルの割合は NaN。
    """
    lonlat = np.asarray(lonlat, dtype=np.float64)
    flooded = np.zeros(len(lonlat))
    n = np.zeros(len(lonlat), dtype=np.int64)
    for crs, grid, ix, iy, y in others:
        tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        x, yy = tr.transform(lonlat[:, 0], lonlat[:, 1])
        cx = np.floor((np.asarray(x) - grid["x0"]) / grid["cell_m"]).astype(np.int64)
        cy = np.floor((grid["y1"] - np.asarray(yy)) / grid["cell_m"]).astype(np.int64)
        lab = np.full((int(grid["ny"]), int(grid["nx"])), -1, dtype=np.int8)
        lab[np.asarray(iy, np.int64), np.asarray(ix, np.int64)] = np.asarray(y, np.int8)
        inside = (cx >= 0) & (cx < grid["nx"]) & (cy >= 0) & (cy < grid["ny"])
        v = np.full(len(lonlat), -1, dtype=np.int8)
        v[inside] = lab[cy[inside], cx[inside]]
        ok = v >= 0
        flooded[ok] += v[ok]
        n[ok] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(n > 0, flooded / np.maximum(n, 1), np.nan)
    return rate, n


# %% 6: 部分的に浸水したセル
def partial_mask(y: np.ndarray, f_flood_any: np.ndarray, lo: float = PARTIAL_COVER) -> np.ndarray:
    """y = 0 なのに浸水被覆が lo を超える（部分的に浸水した）セル。"""
    return (np.asarray(y) == 0) & (np.asarray(f_flood_any, np.float64) > lo)
