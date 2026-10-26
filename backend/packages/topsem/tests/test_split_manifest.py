from pathlib import Path

import pytest

from topsem.splits.manifest import (
    FoldDef,
    SplitSpec,
    SplitUnit,
    load_split,
    make_watershed_holdout,
    save_split,
)


def test_make_and_save_load_roundtrip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("topsem.splits.manifest.SPLITS_DIR", tmp_path)
    units = [("2019", "a"), ("2019", "b"), ("2024", "b")]
    spec = make_watershed_holdout(units, "t1")
    path = save_split(spec)
    assert path.is_file()
    back = load_split("t1")
    assert [f.name for f in back.folds] == ["hold_a", "hold_b"]


def test_validate_flags_missing_and_extra_units() -> None:
    spec = SplitSpec(
        version="t2",
        kind="watershed_holdout",
        folds=[FoldDef(name="hold_a", holdout=[SplitUnit(event="e", watershed="a")])],
    )
    ok = spec.validate(["e", "e"], ["a", "b"])
    assert not ok["ok"]
    assert any("never held out" in p for p in ok["problems"])
    spec2 = SplitSpec(
        version="t2",
        kind="watershed_holdout",
        folds=[
            FoldDef(name="h1", holdout=[SplitUnit(event="e", watershed="a")]),
            FoldDef(name="h2", holdout=[SplitUnit(event="e", watershed="b")]),
        ],
    )
    assert spec2.validate(["e", "e"], ["a", "b"])["ok"]


def test_same_key_multi_units_in_one_fold_is_ok() -> None:
    """同流域に束ねた sub-basin unit（同一 (event,watershed) が 1 fold に複数）は正常。

    v3 で筑後本川+矢部+大牟田が同一 watershed キーになり、entry 数で数えた
    validator が誤検知した回帰（2026-09-22）。
    """
    spec = SplitSpec(
        version="t5b",
        kind="watershed_holdout",
        folds=[
            FoldDef(
                name="hold_chikugo",
                holdout=[SplitUnit(event="2020", watershed="chikugo")] * 3,
            ),
        ],
    )
    res = spec.validate(["2020", "2020", "2020"], ["chikugo"] * 3)
    assert res["ok"], res["problems"]
    assert res["n_units"] == 1


def test_same_key_in_two_folds_still_flagged() -> None:
    spec = SplitSpec(
        version="t5c",
        kind="watershed_holdout",
        folds=[
            FoldDef(name="f1", holdout=[SplitUnit(event="e", watershed="w")]),
            FoldDef(name="f2", holdout=[SplitUnit(event="e", watershed="w")]),
        ],
    )
    res = spec.validate(["e"], ["w"])
    assert not res["ok"]
    assert any(">1 fold" in p for p in res["problems"])


def test_temporal_train_mask_excludes_future_events() -> None:
    """time_2020 fold の訓練側は 2024 を含んではならない（レビュー A-1 の回帰）。

    holdout の素朴補集合（旧 run_split 実装）は未来イベントを訓練に入れて
    temporal 系を event 系と bit 一致にしていた。train_mask が両側除外する。
    """
    spec = SplitSpec(
        version="t3",
        kind="temporal_holdout",
        event_dates={"prior": "1990-01-01", "mid": "2000-01-01", "future": "2020-01-01"},
        folds=[FoldDef(name="time_mid", holdout=[SplitUnit(event="mid", watershed="w1")])],
    )
    events = ["prior", "prior", "mid", "future"]
    wss = ["wA", "wB", "w1", "wA"]
    tr = spec.train_mask(spec.folds[0], events, wss)
    assert tr == [True, True, False, False], "future は訓練からも除外、mid は test 側"
    te = spec.fold_mask(spec.folds[0], events, wss)
    assert te == [False, False, True, False]
    assert sum(tr) == 2


def test_temporal_multi_event_fold_flagged() -> None:
    spec = SplitSpec(
        version="t3b",
        kind="temporal_holdout",
        event_dates={"a": "2000-01-01", "b": "2010-01-01"},
        folds=[
            FoldDef(
                name="bad",
                holdout=[SplitUnit(event="a", watershed="w"), SplitUnit(event="b", watershed="w")],
            ),
        ],
    )
    res = spec.validate(["a", "b"], ["w", "w"])
    assert not res["ok"]
    assert any(">1 event" in p for p in res["problems"])


def test_event_axis_holdout() -> None:
    spec = make_watershed_holdout(
        [("a", "w1"), ("b", "w1")], "t4", kind="event_holdout"
    )
    assert {f.name for f in spec.folds} == {"hold_a", "hold_b"}


def test_missing_manifest_raises(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("topsem.splits.manifest.SPLITS_DIR", tmp_path)
    with pytest.raises(FileNotFoundError):
        load_split("nope")


def test_temporal_unevaluable_units_are_not_missing() -> None:
    """temporal で前期なしイベントの unit は holdout 不要（評価不能は正当）。"""
    spec = SplitSpec(
        version="t5",
        kind="temporal_holdout",
        event_dates={"old": "2000-01-01", "new": "2010-01-01"},
        folds=[FoldDef(name="h_new", holdout=[SplitUnit(event="new", watershed="w1")])],
    )
    res = spec.validate(["old", "new"], ["w1", "w1"])
    assert res["ok"], res["problems"]
    assert res["n_units"] == 2
    assert res["n_units_evaluable"] == 1


def test_watershed_split_still_requires_full_coverage() -> None:
    spec = SplitSpec(
        version="t6",
        kind="watershed_holdout",
        folds=[FoldDef(name="h1", holdout=[SplitUnit(event="e", watershed="a")])],
    )
    res = spec.validate(["e", "e"], ["a", "b"])
    assert not res["ok"]
    assert "1 units never held out" in res["problems"]
