"""amendment B-10 の実装の試験。

距離帯・指標・交互作用、二段目の C の選び方、集落の中で相関させた一様乱数、
JRC のハザード図のタイルと集計。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from topsem.world.b6 import cell_keys, hashed_uniform
from topsem.world.b10 import (
    BAND_NAMES,
    C_GRID,
    JRC_RPS,
    band_indicators,
    choose_c,
    correlated_uniform,
    fill_hazard,
    indicators,
    interactions,
    jrc_features,
    jrc_file,
    jrc_tiles,
)


def test_band_indicators_partition() -> None:
    d = np.array([0.0, 500.0, 500.1, 1000.0, 1000.1, 5e4])
    b = band_indicators(d)
    assert b.shape == (6, len(BAND_NAMES))
    assert (b.sum(axis=1) == 1).all()
    assert b[:, 0].tolist() == [1, 1, 0, 0, 0, 0]
    assert b[:, 1].tolist() == [0, 0, 1, 1, 0, 0]
    assert b[:, 2].tolist() == [0, 0, 0, 0, 1, 1]


def test_indicators_unknown_level_is_zero() -> None:
    v = indicators(["DE", "FR", "", "XX"], ["DE", "FR"])
    assert v.tolist() == [[1, 0], [0, 1], [0, 0], [0, 0]]
    assert indicators(["DE"], []).shape == (1, 0)


def test_interactions_order() -> None:
    x = np.array([[1.0, 2.0], [3.0, 4.0]])
    ind = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    out = interactions(x, ind)
    # 特徴ごとに指標を並べる: f0×b0, f0×b1, f0×b2, f1×b0, ...
    assert out.tolist() == [[1, 0, 0, 2, 0, 0], [0, 0, 3, 0, 0, 4]]


def test_choose_c_prefers_weak_penalty_with_signal() -> None:
    rng = np.random.default_rng(3)
    n = 4000
    x = rng.normal(size=(n, 2))
    y = (rng.random(n) < 1 / (1 + np.exp(-(2.5 * x[:, 0])))).astype(int)
    inner = np.arange(n) % 5
    c, curve = choose_c(x, y, inner)
    assert len(curve) == len(C_GRID)
    assert c >= 1e-2  # 強い信号では 10^-4 のような強い罰則は選ばない
    assert min(curve) == curve[C_GRID.index(c)]


def test_correlated_uniform_margins_and_clustering() -> None:
    n_s, per = 2000, 20
    gids = np.repeat(np.arange(n_s), per)
    ix = np.arange(len(gids))
    keys = cell_keys(["A"] * len(gids), ix % 60000, ix // 60000)
    u = hashed_uniform(keys, 7)
    v = correlated_uniform(u, gids, rep=0, rho=0.5)
    # 周辺は一様（平均 0.5、分散 1/12、10 分位の割合がほぼ等しい）
    assert abs(v.mean() - 0.5) < 0.01
    assert abs(v.var() - 1 / 12) < 0.004
    h = np.histogram(v, bins=10, range=(0, 1))[0] / len(v)
    assert np.abs(h - 0.1).max() < 0.012
    # 集落の中で相関する（集落の平均の分散が独立のときより大きい）
    m = pd.Series(v).groupby(gids).mean().to_numpy()
    assert m.var() > 5 * (1 / 12) / per
    # 集落の無いセルはそのまま、δ によらず同じ集落には同じ z_s（回で決まる）
    w = correlated_uniform(u, np.where(gids < 10, -1, gids), rep=0, rho=0.5)
    assert np.array_equal(w[gids < 10], u[gids < 10])
    assert np.array_equal(correlated_uniform(u, gids, 0, 0.5), v)
    assert not np.array_equal(correlated_uniform(u, gids, 1, 0.5), v)


def test_correlated_uniform_rho_zero_is_identity() -> None:
    u = np.linspace(0.001, 0.999, 101)
    v = correlated_uniform(u, np.zeros(101, dtype=int), rep=0, rho=0.0)
    assert np.allclose(v, u, atol=1e-9)


def _index() -> dict:
    def poly(w: float, s: float, e: float, n: float) -> dict:
        return {"type": "Polygon", "coordinates": [[[w, n], [e, n], [e, s], [w, s], [w, n]]]}
    return {"features": [
        {"properties": {"id": 121, "name": "N60_E10"}, "geometry": poly(9.98, 50.0, 19.98, 60.0)},
        {"properties": {"id": 114, "name": "N50_W0"}, "geometry": poly(-0.02, 40.0, 9.98, 50.0)},
    ]}


def test_jrc_tiles_and_files() -> None:
    idx = _index()
    assert jrc_tiles(idx, (12.8, 51.5, 13.1, 51.7)) == [(121, "N60_E10")]
    assert jrc_tiles(idx, (9.0, 49.5, 11.0, 50.5)) == [(114, "N50_W0"), (121, "N60_E10")]
    assert jrc_tiles(idx, (30.0, 0.0, 31.0, 1.0)) == []
    assert jrc_file((121, "N60_E10"), 100) == "RP100/ID121_N60_E10_RP100_depth.tif"
    pw = "Permanent_WaterBodies/ID121_N60_E10_permanent_water.tif"
    assert jrc_file((121, "N60_E10"), None) == pw


def test_jrc_features_fraction_and_depth() -> None:
    # 1 セル = 2 × 2 画素
    meta = {"grid": {"nx": 2, "ny": 1, "cell_m": 250.0, "fine_res_m": 125.0}}
    nan = np.nan
    d100 = np.array([[2.0, -9999.0, nan, nan],
                     [1.0, -9999.0, nan, nan]], dtype=np.float32)
    pw = np.array([[nan, nan, nan, nan],
                   [1.0, nan, nan, nan]], dtype=np.float32)
    warped = {rp: d100.copy() for rp in JRC_RPS} | {None: pw}
    f = jrc_features(warped, meta)
    # セル 0: 恒久水域の 1 画素を除く 3 画素のうち 1 画素が浸水（深さ 2）
    # セル 1: タイルの外
    assert np.isclose(f["hz_f100"][0], 1 / 3)
    assert np.isclose(f["hz_d100"][0], 2 / 3)
    assert np.isclose(f["hz_cover"][0], 1.0)
    assert np.isnan(f["hz_f100"][1])
    assert f["hz_cover"][1] == 0
    df = pd.DataFrame({"hz_f100": f["hz_f100"], "hz_d100": f["hz_d100"]})
    info = fill_hazard(df, ["hz_f100", "hz_d100"])
    assert info["share_missing_hz_f100"] == 0.5
    assert df["hz_f100"].tolist()[1] == 0.0
