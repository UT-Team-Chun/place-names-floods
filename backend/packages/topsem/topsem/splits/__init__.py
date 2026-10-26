"""凍結分割（watershed / event / temporal / strict）。

ルール（CLAUDE.md「漏洩防止」）: 分割は data/splits/ の凍結 manifest から読む。
その場で train_test_split しない。格子単位は 250 m メッシュを主、100/500 m を
感度分析に使う。**格子の定義関数の正本はここ**（labels.grid はこれを使う）。
"""

from __future__ import annotations

import math

_DEG_PER_M = 1.0 / 111_320.0


def cell_of_point(
    lon: float,
    lat: float,
    size_m: float = 250.0,
    lat_ref: float | None = None,
) -> tuple[int, int]:
    """点 (lon, lat) を地域基準緯度 lat_ref の size_m 格子に対応付ける。

    lat_ref=None の後方互換モードでは点自身の纬度を使う（非一貫に注意 —
    リージョン解析では必ず lat_ref を渡す）。正式版は統計グリッド code へ
    置き換え、置き換えたら split version を上げる。
    """
    cell_lat = size_m * _DEG_PER_M
    iy = math.floor(lat / cell_lat)
    ref = lat if lat_ref is None else lat_ref
    cell_lon = size_m * _DEG_PER_M / math.cos(math.radians(ref))
    ix = math.floor(lon / cell_lon)
    return ix, iy


def cell_center(
    ix: int,
    iy: int,
    size_m: float = 250.0,
    lat_ref: float = 0.0,
) -> tuple[float, float]:
    """(lon, lat) 格子中心。cell_of_point と相互に整合する。"""
    cell_lat = size_m * _DEG_PER_M
    lat = (iy + 0.5) * cell_lat
    cell_lon = size_m * _DEG_PER_M / math.cos(math.radians(lat_ref))
    return (ix + 0.5) * cell_lon, lat
