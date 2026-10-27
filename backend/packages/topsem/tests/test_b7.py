"""事後 B-7 の部品（topsem.world.b7）のテスト。"""

from __future__ import annotations

import numpy as np
import pytest
from pyproj import Transformer

from topsem.world.b7 import (
    aoi_means,
    centre_within,
    days_apart,
    partial_mask,
    recurrence,
    unit_sse,
)


def test_aoi_means_weights_by_cell_and_broadcasts() -> None:
    x = np.array([[1.0, 0.0], [3.0, 2.0], [10.0, 4.0]])
    aoi = np.array(["a", "a", "b"])
    m = aoi_means(x, aoi)
    np.testing.assert_allclose(m, [[2.0, 1.0], [2.0, 1.0], [10.0, 4.0]])
    np.testing.assert_allclose(aoi_means(np.array([1.0, 3.0, 5.0]), aoi), [2.0, 2.0, 5.0])


def test_centre_within_sums_to_zero_in_each_aoi() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(50, 3))
    aoi = rng.choice(["a", "b", "c"], size=50)
    c = centre_within(x, aoi)
    for a in "abc":
        np.testing.assert_allclose(c[aoi == a].sum(axis=0), 0.0, atol=1e-12)


def test_unit_sse_sums_squared_errors_by_place() -> None:
    y = np.array([0.0, 1.0, 1.0, 0.0])
    p = np.array([0.5, 0.5, 1.0, 0.0])
    s = unit_sse(y, p, np.array([0, 0, 1, 2]), 4)
    np.testing.assert_allclose(s, [0.5, 0.0, 0.0, 0.0])


def test_days_apart_and_unreadable_dates() -> None:
    assert days_apart("2021-07-14T00:00:00Z", "2021-07-24T12:00:00") == pytest.approx(10.5)
    assert days_apart("", "2021-07-24") == float("inf")


def _grid_lonlat(crs: str, grid: dict, ix: np.ndarray, iy: np.ndarray) -> np.ndarray:
    x = grid["x0"] + (ix + 0.5) * grid["cell_m"]
    y = grid["y1"] - (iy + 0.5) * grid["cell_m"]
    lon, lat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(x, y)
    return np.column_stack([lon, lat])


def test_recurrence_reads_other_events_on_their_own_grids() -> None:
    crs_b = "+proj=laea +lat_0=50 +lon_0=10 +datum=WGS84 +units=m +no_defs"
    crs_c = "+proj=laea +lat_0=50.0005 +lon_0=10.0008 +datum=WGS84 +units=m +no_defs"
    grid = {"x0": -500.0, "y1": 500.0, "nx": 4, "ny": 4, "cell_m": 250.0}
    # 対象の点: b の格子の (0,0)・(1,0)・(3,3) の中心と、どちらの格子の外の点
    pts = np.vstack([_grid_lonlat(crs_b, grid, np.array([0, 1, 3]), np.array([0, 0, 3])),
                     [[11.0, 51.0]]])
    # b: (0,0) 浸水、(1,0) 浸水なし、(3,3) は評価できない（ラベル無し）
    other_b = (crs_b, grid, np.array([0, 1]), np.array([0, 0]), np.array([1, 0]))
    # c: 対象の点を c の格子に写した位置に、(0,0) の点を浸水なし、(1,0) の点を浸水として置く
    tr = Transformer.from_crs("EPSG:4326", crs_c, always_xy=True)
    x, y = tr.transform(pts[:2, 0], pts[:2, 1])
    cx = np.floor((x - grid["x0"]) / grid["cell_m"]).astype(int)
    cy = np.floor((grid["y1"] - y) / grid["cell_m"]).astype(int)
    assert ((cx >= 0) & (cx < 4) & (cy >= 0) & (cy < 4)).all()
    other_c = (crs_c, grid, cx, cy, np.array([0, 1]))
    rate, n = recurrence(pts, [other_b, other_c])
    np.testing.assert_array_equal(n, [2, 2, 0, 0])
    np.testing.assert_allclose(rate[:2], [0.5, 0.5])
    assert np.isnan(rate[2:]).all()


def test_partial_mask_only_flags_unflooded_cells_with_cover() -> None:
    y = np.array([0, 0, 0, 1])
    f = np.array([0.0, 0.05, 0.3, 0.9])
    np.testing.assert_array_equal(partial_mask(y, f), [False, False, True, False])
