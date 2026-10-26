import math

from shapely.geometry import box
from shapely.ops import unary_union

from topsem.labels.grid import _DEG_PER_M, cell_spans, rasterize_wet_frac

CELL = 250.0
LAT_REF = 36.5
CELL_LAT = CELL * _DEG_PER_M
CELL_LON = CELL * _DEG_PER_M / math.cos(math.radians(LAT_REF))


def _cells_box(iy: int, ix: int, ny: int, nx: int) -> box:
    return box(ix * CELL_LON, iy * CELL_LAT, (ix + nx) * CELL_LON, (iy + ny) * CELL_LAT)


def test_aligned_square_maps_to_full_cells() -> None:
    sq = _cells_box(16000, 55000, 2, 2)
    cells = rasterize_wet_frac(sq, size_m=CELL, buffer_m=0.0, lat_ref=LAT_REF)
    assert len(cells) == 4
    assert all(fr > 0.999 for *_, fr in cells)


def test_half_shift_splits_frac() -> None:
    shifted = box(
        55000 * CELL_LON + CELL_LON / 2,
        16000 * CELL_LAT,
        55000 * CELL_LON + 1.5 * CELL_LON,
        16000 * CELL_LAT + CELL_LAT,
    )
    cells = rasterize_wet_frac(shifted, size_m=CELL, buffer_m=0.0, lat_ref=LAT_REF)
    assert len(cells) == 2
    assert all(abs(fr - 0.5) < 0.02 for *_, fr in cells)


def test_multipolygon_union_span() -> None:
    g1 = _cells_box(16000, 55000, 1, 1)
    g2 = _cells_box(16000, 55004, 1, 1)
    cells = rasterize_wet_frac(unary_union([g1, g2]), size_m=CELL, buffer_m=0.0, lat_ref=LAT_REF)
    full = [c for c in cells if c[2] > 0.99]
    assert len(full) == 2


def test_buffer_widens_scan_window_not_output() -> None:
    sq = _cells_box(16000, 55000, 1, 1)
    _, _, iy0a, iy1a, ix0a, ix1a = cell_spans(sq, CELL, 0.0, LAT_REF)
    _, _, iy0b, iy1b, ix0b, ix1b = cell_spans(sq, CELL, CELL * 1.2, LAT_REF)
    assert (iy1b - iy0b) > (iy1a - iy0a)
    assert (ix1b - ix0b) > (ix1a - ix0a)
    # 出力は wet_frac>0 のセルのみ（乾セルは後段の解析域定義で作る）
    assert len(rasterize_wet_frac(sq, size_m=CELL, buffer_m=CELL * 1.2, lat_ref=LAT_REF)) == 1
