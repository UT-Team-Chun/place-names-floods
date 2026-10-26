"""cell_features（格子集約・設計行列）の単体テスト。データ不要な合成 frame で。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from topsem.analysis.cell_features import _aggregate_to_grid, build_xy
from topsem.analysis.models import arm_names
from topsem.features.route_b import DIM, embed_name


def _grid(cells: list[tuple[int, int]]) -> pd.DataFrame:
    ix, iy = zip(*cells, strict=True)
    return pd.DataFrame(
        {"ix": ix, "iy": iy, "wet_frac": [0.0] * len(cells)},
        columns=["ix", "iy", "wet_frac"],
    )


def test_aggregate_fills_zero_without_rows() -> None:
    g = _aggregate_to_grid(_grid([(0, 0), (1, 1)]), [], pd.DataFrame({"name": []}), {})
    assert g["a_flood"].eq(0.0).all()
    assert g["n_addr"].eq(0.0).all()
    assert np.vstack(g["emb"].to_numpy()).shape == (2, DIM)
    assert np.vstack(g["emb"].to_numpy()).sum() == 0.0


def test_aggregate_max_and_count_per_cell() -> None:
    # セル (1,2) に住所 2 件、セル (5,6) に 1 件、無住所セル (9,9)
    grid = _grid([(1, 2), (5, 6), (9, 9)])
    rows = [(1, 2, 0.0, 1.0, 0.0), (1, 2, 2.0, 1.0, 1.0), (5, 6, 1.0, 3.0, 0.0)]
    sub = pd.DataFrame({"name": ["a", "b", "c"]})
    g = _aggregate_to_grid(grid, rows, sub, {})
    assert g.loc[0, "a_flood"] == 2.0  # max（同じセルの洪水語は打ち消さない）
    assert g.loc[0, "a_any"] == 2.0  # sum
    assert g.loc[0, "a_sub"] == 1.0
    assert g.loc[0, "n_addr"] == 2.0
    assert g.loc[1, "n_addr"] == 1.0
    assert g.loc[2, ["a_flood", "a_any", "a_sub", "n_addr"]].eq(0.0).all()


def test_aggregate_emb_is_mean_of_member_names() -> None:
    grid = _grid([(1, 2)])
    rows = [(1, 2, 0.0, 1.0, 0.0), (1, 2, 0.0, 1.0, 0.0)]
    sub = pd.DataFrame({"name": ["水田", "河原"]})
    g = _aggregate_to_grid(grid, rows, sub, {})
    want = (np.asarray(embed_name("水田")) + np.asarray(embed_name("河原"))) / 2
    np.testing.assert_allclose(np.vstack(g["emb"].to_numpy()), [want], atol=1e-6)


def test_build_xy_arm_columns() -> None:
    df = pd.DataFrame(
        {
            "ix": [0, 1],
            "iy": [2, 3],
            "elev": [10.0, 20.0],
            "depress": [1.0, 2.0],
            "a_flood": [0.0, 1.0],
            "a_any": [1.0, 2.0],
            "a_sub": [0.0, 0.0],
            "emb": [np.zeros(DIM, dtype=np.float32), np.ones(DIM, dtype=np.float32)],
        }
    )
    assert build_xy(df, "M0").shape == (2, 2)
    assert build_xy(df, "M1").shape == (2, 4)
    assert build_xy(df, "M2").shape == (2, 7)
    # protocol §1: M3 = M1 + 埋め込み（辞書特征は M2 の別枝で累進ではない）
    assert build_xy(df, "M3").shape == (2, 4 + DIM)


def test_build_xy_drop_names_ablation() -> None:
    df = pd.DataFrame(
        {
            "ix": [0],
            "iy": [0],
            "elev": [5.0],
            "depress": [0.5],
            "a_flood": [1.0],
            "a_any": [1.0],
            "a_sub": [0.0],
            "emb": [np.ones(DIM, dtype=np.float32)],
        }
    )
    np.testing.assert_array_equal(
        build_xy(df, "M3", drop_names=True), build_xy(df, "M1")
    )


def test_fit_arms_ensemble_matches_single_seed() -> None:
    """--seeds 1 個は fit_arms と完全一致、複数個は平均した別物になる。"""
    import numpy as np  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415

    from topsem.analysis.models import (  # noqa: PLC0415
        fit_arms,
        fit_arms_ensemble,
    )

    rng = np.random.default_rng(7)
    n = 400
    df = pd.DataFrame(
        {
            "ix": rng.integers(0, 20, n),
            "iy": rng.integers(0, 20, n),
            "elev": rng.normal(30, 10, n),
            "depress": rng.uniform(0, 5, n),
            "a_flood": rng.integers(0, 2, n).astype(float),
            "a_any": rng.integers(0, 3, n).astype(float),
            "a_sub": rng.integers(0, 2, n).astype(float),
            "y": rng.integers(0, 2, n).astype(float),
        }
    )
    df["emb"] = [rng.normal(size=DIM).astype(np.float32) for _ in range(n)]
    tr, te = df.iloc[:300], df.iloc[300:]
    a1 = fit_arms(tr, te, 42)
    e1 = fit_arms_ensemble(tr, te, [42])
    for arm in ("M0", "M1", "M2", "M3"):
        np.testing.assert_allclose(a1[arm], e1[arm])
    e3 = fit_arms_ensemble(tr, te, [42, 43, 44])
    assert e3["M3"].shape == a1["M3"].shape


def test_build_xy_m1_river_flag_columns():
    """A-4: dist_river 列があれば M1 は +1 列、M1r は base+距離のみ。無ければ旧来と同一。"""
    df = pd.DataFrame(
        {
            "ix": [0, 1, 2], "iy": [0, 1, 2],
            "elev": [1.0, 2.0, 3.0], "depress": [0.5, 0.5, 0.5],
            "dist_river": [10.0, 200.0, 5000.0],
            "a_flood": [0.0, 1.0, 0.0], "a_any": [0.0, 2.0, 0.0], "a_sub": [0.0, 1.0, 0.0],
            "emb": [np.zeros(64, dtype=np.float32)] * 3,
        }
    )
    x1 = build_xy(df, "M1")
    assert x1.shape[1] == 5  # ix,iy,elev,depress,log1p(dist)
    xr = build_xy(df, "M1r")
    assert xr.shape[1] == 3
    assert np.isclose(xr[0, 2], np.log1p(10.0))
    xrn = build_xy(df, "M1rN")
    assert xrn.shape[1] == 3 + 64  # base+距離+埋め込み（地形なし・r3-F4）
    assert arm_names(df) == ("M0", "M1", "M2", "M3", "M1r", "M1rN")
    assert arm_names(df.drop(columns=["dist_river"])) == ("M0", "M1", "M2", "M3")
    legacy = build_xy(df.drop(columns=["dist_river"]), "M1")
    assert np.array_equal(legacy[:, :4], x1[:, :4])  # 列追加のみ・既存は不変
    # drop_names 消融去は腕名に依らず M1 同一設計（r3 で修正した分岐の回帰テスト）
    assert np.array_equal(build_xy(df, "M3", drop_names=True), x1)
