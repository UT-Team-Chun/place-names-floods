"""ポリゴン → メッシュ（label side）。

格子は**地域基準緯度 lat_ref での一定経度幅**の一貫した格子（暫定実装、
統計グリッド code へ置き換え予定。置き換えたら split version を上げる）。
境界とセル矩形の intersection 面積比 `wet_frac` を出す。ラベル化閾値
（wet_frac >= t、既定 0.5）は結果を見て変えない（protocol §2）。
窓外は情報なし（UNOBSERVED）なので行を作らない。
"""

from __future__ import annotations

import math

from shapely.geometry import box
from shapely.geometry.base import BaseGeometry

from topsem.splits import cell_of_point  # 格子の点対応は splits が正本（再エクスポート）

__all__ = ["cell_of_point", "cell_spans", "rasterize_wet_frac", "rasterize_window"]

_DEG_PER_M = 1.0 / 111_320.0


def cell_spans(
    geom: BaseGeometry,
    size_m: float = 250.0,
    buffer_m: float = 0.0,
    lat_ref: float | None = None,
) -> tuple[float, float, int, int, int, int]:
    """(cell_lat, cell_lon, iy0, iy1, ix0, ix1): geom+buffer を覆う格子範囲。"""
    minx, miny, maxx, maxy = geom.bounds
    ref = lat_ref if lat_ref is not None else (miny + maxy) / 2.0
    pad = buffer_m * _DEG_PER_M
    cell_lat = size_m * _DEG_PER_M
    cell_lon = size_m * _DEG_PER_M / math.cos(math.radians(ref))
    iy0 = math.floor((miny - pad) / cell_lat)
    iy1 = math.ceil((maxy + pad) / cell_lat)
    ix0 = math.floor((minx - pad) / cell_lon)
    ix1 = math.ceil((maxx + pad) / cell_lon)
    return cell_lat, cell_lon, iy0, iy1, ix0, ix1


def rasterize_wet_frac(
    geom: BaseGeometry,
    size_m: float = 250.0,
    buffer_m: float = 2000.0,
    lat_ref: float | None = None,
) -> list[tuple[int, int, float]]:
    """geom を含み buffer_m 拡張した窓を size_m 格子に切って (ix, iy, wet_frac) を返す。

    wet_frac は 0 より大きいセルのみ出力（辺共有だけの接触セルは除外）。
    """
    cell_lat, cell_lon, iy0, iy1, ix0, ix1 = cell_spans(geom, size_m, buffer_m, lat_ref)
    out: list[tuple[int, int, float]] = []
    for iy in range(iy0, iy1 + 1):
        y0 = iy * cell_lat
        for ix in range(ix0, ix1 + 1):
            x0 = ix * cell_lon
            cell = box(x0, y0, x0 + cell_lon, y0 + cell_lat)
            inter = cell.intersection(geom)
            if not inter.is_empty and inter.area > 0.0:
                out.append((ix, iy, inter.area / cell.area))
    return out


def rasterize_window(
    geom: BaseGeometry,
    size_m: float = 250.0,
    buffer_m: float = 2000.0,
    lat_ref: float | None = None,
) -> list[tuple[int, int, float]]:
    """窓内全セルを wet_frac 付きで返す（湿 0.0 = この製品では未浸水扱い）。

    注意: GSI 浸水推定図では着色外にも実際の浸水があり得る（免責）。
    パイロットでは dry として扱い、妥当性チェックで感度分析する。
    """
    cell_lat, cell_lon, iy0, iy1, ix0, ix1 = cell_spans(geom, size_m, buffer_m, lat_ref)
    out: list[tuple[int, int, float]] = []
    for iy in range(iy0, iy1 + 1):
        y0 = iy * cell_lat
        for ix in range(ix0, ix1 + 1):
            x0 = ix * cell_lon
            cell = box(x0, y0, x0 + cell_lon, y0 + cell_lat)
            frac = 0.0
            if cell.intersects(geom):
                inter = cell.intersection(geom)
                frac = inter.area / cell.area
            out.append((ix, iy, frac))
    return out
