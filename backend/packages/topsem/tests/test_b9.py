"""amendment B-9 の実装の試験。

注入の向き・配列の δ・集落ごとの重み・永久水域の印・FABDEM の名前と読む場所。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from topsem.world.b6 import hashed_uniform, inject_labels
from topsem.world.b9 import (
    FABDEM_BASE,
    fabdem_path,
    fabdem_tile_name,
    fabdem_zip_name,
    gsw_perm_mask,
    settlement_weights,
    tiles_for,
)


def _old_inject(y: np.ndarray, p: np.ndarray, target: np.ndarray, delta: float,
                u: np.ndarray) -> np.ndarray:  # B-6 の実装（δ > 0 の値が変わらないことを確かめる）
    p = np.clip(p.astype(np.float64), 1e-6, 1 - 1e-6)
    e = np.exp(delta)
    p2 = p * e / (1 - p + p * e)
    q = (p2 - p) / (1 - p)
    flip = target & (y == 0) & (u < q)
    return np.where(flip, 1, y).astype(np.int8)


def test_positive_delta_unchanged_from_b6() -> None:
    rng = np.random.default_rng(1)
    n = 50_000
    p = rng.uniform(0.001, 0.6, n)
    y = (rng.random(n) < p).astype(np.int8)
    target = rng.random(n) < 0.3
    u = hashed_uniform(np.arange(n, dtype=np.uint64), 5)
    for d in (0.05, 0.25, 1.0, 3.0):
        assert np.array_equal(inject_labels(y, p, target, d, u), _old_inject(y, p, target, d, u))


def test_negative_delta_multiplies_odds() -> None:
    rng = np.random.default_rng(0)
    n = 400_000
    p = np.full(n, 0.2)
    y = (rng.random(n) < p).astype(np.int8)
    target = np.zeros(n, dtype=bool)
    target[: n // 2] = True
    y2 = inject_labels(y, p, target, np.log(0.5), hashed_uniform(np.arange(n, dtype=np.uint64), 2))
    p2 = 0.2 * 0.5 / (1 - 0.2 + 0.2 * 0.5)
    assert abs(y2[target].mean() - p2) < 0.004
    assert np.array_equal(y2[~target], y[~target])
    assert (y2 <= y).all()  # 0 を 1 にはしない


def test_array_delta_and_zero() -> None:
    n = 10_000
    p = np.full(n, 0.3)
    y = np.zeros(n, dtype=np.int8)
    y[::2] = 1
    target = np.ones(n, dtype=bool)
    d = np.where(np.arange(n) < n // 2, 0.0, 1.0)
    y2 = inject_labels(y, p, target, d, hashed_uniform(np.arange(n, dtype=np.uint64), 3))
    assert np.array_equal(y2[: n // 2], y[: n // 2])  # δ = 0 は変えない
    assert (y2[n // 2:] >= y[n // 2:]).all()
    assert (y2[n // 2:] != y[n // 2:]).any()


def test_settlement_weights_exponential_and_fixed() -> None:
    g = np.arange(200_000, dtype=np.int64)
    w = settlement_weights(g, 7)
    assert abs(w.mean() - 1.0) < 0.01
    assert abs(w.std() - 1.0) < 0.02
    assert (w >= 0).all()
    assert np.array_equal(w, settlement_weights(g, 7))
    assert not np.array_equal(w, settlement_weights(g, 8))
    # 同じ集落は同じ重み（セルの並びによらない）
    gg = np.array([5, 9, 5, 9], dtype=np.int64)
    ww = settlement_weights(gg, 7)
    assert ww[0] == ww[2]
    assert ww[1] == ww[3]


def test_gsw_perm_mask() -> None:
    df = pd.DataFrame({"aoi": ["a", "a", "b", "b"], "ix": [0, 1, 0, 1], "iy": [0, 0, 0, 0]})
    gsw = pd.DataFrame({"aoi": ["a", "a", "b"], "ix": [0, 1, 0], "iy": [0, 0, 0],
                        "occ": [80.0, 10.0, 50.0]})
    assert gsw_perm_mask(df, gsw, 50.0).tolist() == [True, False, True, False]


def test_fabdem_names() -> None:
    assert fabdem_zip_name(51, 6) == "N50E000-N60E010_FABDEM_V1-2.zip"
    assert fabdem_tile_name(51, 6) == "N51E006_FABDEM_V1-2.tif"
    assert fabdem_zip_name(5, -5) == "N00W010-N10E000_FABDEM_V1-2.zip"
    assert fabdem_zip_name(5, 175) == "N00E170-N10W180_FABDEM_V1-2.zip"
    assert fabdem_zip_name(-3, -45) == "S10W050-N00W040_FABDEM_V1-2.zip"
    assert fabdem_tile_name(-3, -45) == "S03W045_FABDEM_V1-2.tif"
    want = [(50, 5), (50, 6), (50, 7), (51, 5), (51, 6), (51, 7)]
    assert tiles_for((5.5, 50.2, 7.1, 51.9)) == want


def test_fabdem_path_prefers_tiles(tmp_path: Path) -> None:
    t, z = fabdem_tile_name(51, 6), fabdem_zip_name(51, 6)
    assert fabdem_path(51, 6, None) == f"/vsizip//vsicurl/{FABDEM_BASE}/{z}/{t}"
    assert fabdem_path(51, 6, tmp_path).startswith("/vsizip//vsicurl/")
    (tmp_path / "zips").mkdir()
    (tmp_path / "zips" / z).write_bytes(b"")
    assert fabdem_path(51, 6, tmp_path) == f"/vsizip/{tmp_path / 'zips' / z}/{t}"
    (tmp_path / "tiles").mkdir()
    (tmp_path / "tiles" / t).write_bytes(b"")
    assert fabdem_path(51, 6, tmp_path) == str(tmp_path / "tiles" / t)
    sea = fabdem_tile_name(55, 4)
    (tmp_path / "tiles" / f"{sea}.absent").write_text(z)
    assert fabdem_path(55, 4, tmp_path) == str(tmp_path / "tiles" / sea)  # 読むと失敗するパス
    assert not (tmp_path / "tiles" / sea).exists()
