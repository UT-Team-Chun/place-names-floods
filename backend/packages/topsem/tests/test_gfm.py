"""事後 B-3 の独立ラベル（GFM）の規則の単体テスト。合成データのみ（ネットワークに出ない）。"""

from __future__ import annotations

import datetime as dt
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from topsem.world.gfm import (
    cell_fractions,
    combine_scenes,
    gfm_include,
    group_scenes,
    item_provenance,
    post_event_window,
    sensor_group,
    valid_and_wet,
)

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):  # noqa: ANN202
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_post_event_window_uses_post_event_dates_only() -> None:
    src = [{"phase": "Pre-event", "date": "01/05/2023", "name": "Sentinel-2"},
           {"phase": "Post-event", "date": "24/05/2023", "name": "Sentinel-1"},
           {"phase": "Post-event", "date": "2023/05/26", "name": "COSMO-SkyMed"},
           {"phase": "Post-event", "date": "Not Applicable or dd/mm/aaaa", "name": "OSM"}]
    assert post_event_window(src) == (dt.date(2023, 5, 23), dt.date(2023, 5, 27))
    assert post_event_window(src[:1]) is None
    assert post_event_window([]) is None


def test_sensor_group() -> None:
    def s(*names: str) -> list[dict]:
        return [{"phase": "Post-event", "name": n} for n in names]

    assert sensor_group(s("Sentinel-1", "COSMO-SkyMed")) == "s1"
    assert sensor_group(s("COSMO-SkyMed Second Generation")) == "sar"
    assert sensor_group(s("RCM")) == "sar"
    assert sensor_group(s("Sentinel-2", "WorldView-2")) == "optical"
    assert sensor_group([{"phase": "Pre-event", "name": "Sentinel-1"}]) == "none"


def test_group_scenes_by_acquisition_time() -> None:
    items = [{"id": "a", "properties": {"datetime": "2023-05-24T05:03:27Z"}},
             {"id": "b", "properties": {"datetime": "2023-05-24T05:03:27Z"}},
             {"id": "c", "properties": {"datetime": "2023-05-25T17:10:00Z"}}]
    g = group_scenes(items)
    assert [len(v) for v in g.values()] == [2, 1]


def test_item_provenance_tolerates_missing_asset() -> None:
    it = {"id": "X", "properties": {"processing:version": "1"},
          "assets": {"ensemble_flood_extent": {"href": "a"},
                     "reference_water_mask": {"href": "b"}}}
    pr = item_provenance(it)
    assert pr["id"] == "X"
    assert pr["version"] == "1"
    assert pr["hrefs"]["exclusion_mask"] is None
    assert pr["hrefs"]["ensemble_flood_extent"] == "a"


def test_valid_and_wet_masks() -> None:
    flood = np.array([0, 1, 1, 255, 1, 0], dtype=np.uint8)
    excl = np.array([0, 0, 1, 0, 0, 2], dtype=np.uint8)
    refw = np.array([0, 0, 0, 0, 1, 0], dtype=np.uint8)
    valid, wet = valid_and_wet(flood, excl, refw)
    assert valid.tolist() == [True, True, False, False, False, False]
    assert wet.tolist() == [False, True, False, False, False, False]


def test_cell_fractions_per_250m_cell() -> None:
    # 1 行 2 セル、副格子 2×2: 左セルは 4 画素とも判定可能で 1 画素浸水、
    # 右セルは 1 画素だけ判定可能で浸水
    valid = np.array([[1, 1, 0, 0], [1, 1, 1, 0]], dtype=bool)
    wet = np.array([[0, 1, 0, 0], [0, 0, 1, 0]], dtype=bool)
    v, f = cell_fractions(valid, wet, sub=2, ny=1, nx=2)
    np.testing.assert_allclose(v, [1.0, 0.25])
    np.testing.assert_allclose(f, [0.25, 1.0])


def test_combine_scenes_rule() -> None:
    # セル 0: 撮像 1 で評価可能・非浸水、撮像 2 で評価可能・浸水 → 浸水
    # セル 1: どちらも判定可能な割合 < 0.5 → 評価不能（−1）
    # セル 2: 撮像 1 だけ評価可能、浸水の割合 0.4 → 非浸水
    #         （撮像 2 の浸水は v < 0.5 なので数えない）
    vs = [np.array([0.9, 0.4, 0.8]), np.array([0.6, 0.3, 0.2])]
    fs = [np.array([0.1, 1.0, 0.4]), np.array([0.7, 1.0, 1.0])]
    vmax, ev, y, n = combine_scenes(vs, fs)
    np.testing.assert_allclose(vmax, [0.9, 0.4, 0.8])
    assert ev.tolist() == [True, False, True]
    assert y.tolist() == [1, -1, 0]
    assert n.tolist() == [2, 0, 1]


def test_gfm_include_thresholds() -> None:
    n = 1000
    ems = np.ones(n, dtype=bool)
    ev = np.zeros(n, dtype=bool)
    ev[:500] = True
    y = np.full(n, -1, dtype=np.int8)
    y[:500] = 0
    y[:30] = 1
    ok, info = gfm_include(ems, ev, y)
    assert ok
    assert info == {"n_eval": 500, "n_wet": 30, "wet_frac": 0.06, "why": "ok"}
    assert gfm_include(ems, ev, np.where(y == 1, 0, y).astype(np.int8))[1]["why"] == "few-flooded"
    ev_small = ev.copy()
    ev_small[399:] = False
    assert gfm_include(ems, ev_small, y)[1]["why"] == "few-evaluable"
    assert gfm_include(ems, ev, np.where(ev, 1, -1).astype(np.int8))[1]["why"] == "too-flooded"
    # EMS で評価不能なセルは数えない
    assert gfm_include(np.zeros(n, dtype=bool), ev, y)[1]["n_eval"] == 0


def test_apply_gfm_labels_filters_and_replaces_without_reordering() -> None:
    ev = _load("world_eval")
    df = pd.DataFrame({"aoi": ["A", "A", "B", "A", "C"], "ix": [0, 1, 0, 2, 0],
                       "iy": [0, 0, 0, 0, 0], "y": [1, 1, 1, 0, 0],
                       "fold": [3, 3, 1, 3, 0], "unit": ["g1", "g1", "g2", "g1", "g3"]})
    tables = {"A": pd.DataFrame({"ix": [0, 1, 2], "iy": [0, 0, 0],
                                 "gfm_eval": [True, False, True], "gfm_y": [0, -1, 1]}),
              "B": pd.DataFrame({"ix": [0], "iy": [0], "gfm_eval": [True], "gfm_y": [1]})}
    summary = {"aois": {"A": {"include": True}, "B": {"include": False}, "C": {"include": True}}}
    out, info = ev.apply_gfm_labels(df, summary, tables)
    assert out["aoi"].tolist() == ["A", "A"]  # B は不採用、C は表が無い、A の ix=1 は評価不能
    assert out["ix"].tolist() == [0, 2]
    assert out["y"].tolist() == [0, 1]  # EMS の y（1, 0）を GFM のラベル（0, 1）に置き換える
    assert out["fold"].tolist() == [3, 3]
    assert info["n_cells"] == 2
    assert info["n_cells_before"] == 5
