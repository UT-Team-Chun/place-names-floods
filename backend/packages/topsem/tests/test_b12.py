"""amendment B-12 の実装（topsem.world.b12）のテスト。"""

from __future__ import annotations

import numpy as np

from topsem.world.b12 import PLASMODE_1KM_M, PLASMODE_LOCAL_M, plasmode_target


def test_shapes_select_the_documented_distances() -> None:
    t = np.array([True, True, True, False, True])
    d1 = np.array([100.0, 500.0, 900.0, 100.0, 1500.0])
    assert plasmode_target("A", t, d1).tolist() == t.tolist()
    assert plasmode_target("B", t, d1).tolist() == [True, True, False, False, False]
    assert plasmode_target("CB", t, d1).tolist() == [True, True, False, False, False]
    assert plasmode_target("BK", t, d1).tolist() == [True, True, True, False, False]
    assert plasmode_target("AH", t, d1).tolist() == t.tolist()
    assert plasmode_target("AC", t, d1).tolist() == t.tolist()


def test_bk_contains_b() -> None:
    rng = np.random.default_rng(0)
    t = rng.random(1000) < 0.3
    d1 = rng.uniform(0, 3000, 1000)
    b, bk = plasmode_target("B", t, d1), plasmode_target("BK", t, d1)
    assert np.all(bk[b])
    assert bk.sum() >= b.sum()
    assert PLASMODE_1KM_M == 2 * PLASMODE_LOCAL_M
