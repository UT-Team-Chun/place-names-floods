"""amendment B-16 の実装（topsem.world.b16、層別の置換）のテスト。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from topsem.world.b16 import aoi_grid, detail_fractions, in_box, merge_nearest, nearest_upland
from topsem.world.channels import (
    placebo_lookups_strata,
    placebo_perms,
    placebo_perms_strata,
    site_strata,
)
from topsem.world.ems import build_aoi_labels

from .test_world import _make_product


def _frame(n: int) -> pd.DataFrame:
    return pd.DataFrame({"aoi": ["X"] * n, "P_gid1": list(range(1, n + 1)),
                         "P_gid2": [-1] * n, "P_gid3": [-1] * n})


def test_site_strata_terciles_ties_low_and_missing() -> None:
    gids = np.array([1, 2, 3, 4, 5, 6, 7])
    value = {1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0, 6: 6.0}  # 7 は値が無い
    lab = site_strata(gids, value)
    assert lab.tolist() == [0, 0, 1, 1, 2, 2, -1]
    # 同じ値は低い層へ（分位点と等しい値）
    lab2 = site_strata(np.array([1, 2, 3]), {1: 1.0, 2: 1.0, 3: 1.0})
    assert lab2.tolist() == [0, 0, 0]


def test_placebo_perms_strata_stay_within_strata() -> None:
    df = _frame(6)
    value = {1: 0.0, 2: 0.1, 3: 5.0, 4: 5.1, 5: 9.0, 6: 9.1}
    gids = np.arange(1, 7)
    lab = dict(zip(gids.tolist(), site_strata(gids, value).tolist(), strict=True))
    for seed in range(20):
        m = placebo_perms_strata(df, seed, value)["X"]
        assert sorted(m) == sorted(m.values()) == [1, 2, 3, 4, 5, 6]
        assert all(lab[g] == lab[p] for g, p in m.items())
    # 層が一つだけなら、層の無い置換と同じ集合の上の置換（乱数が違うので並べ方の一致は求めない）
    one = placebo_perms_strata(df, 3, dict.fromkeys(range(1, 7), 1.0))["X"]
    assert sorted(one.values()) == sorted(placebo_perms(df, 3)["X"].values())


def test_placebo_lookups_strata_maps_names() -> None:
    df = _frame(4)
    look = {i: (f"n{i}", "IT") for i in range(1, 5)}
    value = {1: 0.0, 2: 0.0, 3: 9.0, 4: 9.0}
    lk = placebo_lookups_strata(df, look, 5, value)["X"]
    perm = placebo_perms_strata(df, 5, value)["X"]
    assert lk == {g: look[p] for g, p in perm.items()}


def test_nearest_upland_and_merge() -> None:
    lon, lat = np.array([0.0, 1.0]), np.array([0.0, 0.0])
    vlon, vlat, vup = np.array([0.0, 1.01]), np.array([0.001, 0.0]), np.array([100.0, 1e4])
    lupl, dist = nearest_upland(lon, lat, vlon, vlat, vup)
    assert np.allclose(lupl, [2.0, 4.0])
    assert 100 < dist[0] < 120  # 0.001° ≈ 111 m
    empty_l, empty_d = nearest_upland(lon, lat, np.array([]), np.array([]), np.array([]))
    assert np.isnan(empty_l).all()
    assert np.isnan(empty_d).all()
    best_d, best_v = np.array([np.inf, 50.0]), np.array([np.nan, 1.0])
    merge_nearest(best_d, best_v, np.array([10.0, 60.0]), np.array([3.0, 9.0]))
    assert best_d.tolist() == [10.0, 50.0]
    assert best_v.tolist() == [3.0, 1.0]
    box = (-1.0, -1.0, 1.0, 1.0)
    assert in_box(np.array([0.0, 2.0]), np.array([0.0, 0.0]), box).tolist() == [True, False]


def test_detail_fractions_reproduce_the_overview_rules(tmp_path: Path) -> None:
    """C5: 同じ製品を詳細図として落とすと、概観図のラベル（build_aoi_labels）と同じ被覆率。"""
    z = _make_product(tmp_path, footprint=True)
    cols, meta = build_aoi_labels(z)
    out = detail_fractions([z], aoi_grid(meta))
    for k in ("f_aoi", "f_observed", "f_riverine", "f_perm_water", "f_trace", "f_flood_any"):
        np.testing.assert_allclose(out[k], cols[k], atol=1e-9)
    assert out["generations"].tolist() == [meta["schema_generation"]]


def test_detail_fractions_merge_products_counting_overlaps_once(tmp_path: Path) -> None:
    """C5: 複数の詳細図の製品は一つとして合わせ、重なった浸水は一度だけ数える。"""
    z1 = _make_product(tmp_path, footprint=True)
    _, meta = build_aoi_labels(z1)
    d2 = tmp_path / "second"
    d2.mkdir()
    z2 = _make_product(d2, footprint=True)  # 同じ浸水の多角形（西半分）をもう一度
    one = detail_fractions([z1], aoi_grid(meta))
    two = detail_fractions([z1, z2], aoi_grid(meta))
    np.testing.assert_allclose(two["f_riverine"], one["f_riverine"], atol=1e-9)
    assert two["f_riverine"].max() <= 1.0
    assert two["generations"].tolist() == [3, 3]
