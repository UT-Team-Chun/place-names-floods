"""v2 評価ハーネス（dataset / channels / 集計判定）の単体テスト。合成データのみ。"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topsem.world.channels import (
    NREAD_FEATURES,
    nlearn_scores,
    nread_block,
    nread_lookup,
    placebo_lookups,
    placebo_perms,
)
from topsem.world.dataset import (
    AoiInfo,
    assign_folds,
    group_quotas,
    include,
    location_groups,
    merge_by_activation,
    riverine_labels,
)
from topsem.world.lexicon import Lexicon
from topsem.world.nread import NREAD_CATEGORIES, parse, prompt_sha256, user_prompt

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):  # noqa: ANN202
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclass は定義時に sys.modules から自モジュールを引く
    spec.loader.exec_module(mod)
    return mod


def _info(**kw) -> AoiInfo:  # noqa: ANN003
    base = {"tag": "EMSR1_AOI01", "activation": "EMSR1", "countries": ("Italy",),
            "event_time": "2020-01-01", "bounds": (10.0, 45.0, 10.1, 45.1), "n_eval": 1000,
            "n_wet": 100, "riverine_share": 1.0, "footprint_implicit": False}
    return AoiInfo(**{**base, **kw})


@pytest.mark.parametrize(
    ("kw", "names", "ok", "why"),
    [
        ({}, 3, True, "ok"),
        ({"n_eval": 399}, 3, False, "eval<400"),
        ({"n_wet": 19}, 3, False, "wet<20"),
        ({"n_wet": 900}, 3, False, "wetfrac>0.8"),
        ({"riverine_share": 0.4}, 3, False, "not-riverine"),
        # 種別未記載が過半: 内陸なら河川とみなし、沿岸なら主解析から外す
        ({"riverine_share": 0.0, "unspecified_share": 1.0, "coast_km": 50.0}, 3, True, "ok"),
        ({"riverine_share": 0.0, "unspecified_share": 1.0, "coast_km": 3.0}, 3, False,
         "not-riverine"),
        ({}, 0, False, "no_P_name"),
    ],
)
def test_include_rules(kw: dict, names: int, *, ok: bool, why: str) -> None:
    assert include(_info(**kw), n_names_inside=names) == (ok, why)


def test_location_groups_merge_within_buffer_only() -> None:
    b = {
        "A": (10.00, 45.00, 10.10, 45.10),
        "B": (10.15, 45.00, 10.25, 45.10),  # A から約 4 km → 5 km 緩衝で同じ組
        "C": (11.00, 45.00, 11.10, 45.10),  # 70 km 離れる → 別の組
    }
    g = location_groups(b, buffer_km=5.0)
    assert g["A"] == g["B"]
    assert g["A"] != g["C"]


def test_assign_folds_is_deterministic_round_robin() -> None:
    unit = {"a1": "u1", "a2": "u1", "a3": "u2", "a4": "u3"}
    sort = {"u1": ("EU", -10, "u1"), "u2": ("EU", -5, "u2"), "u3": ("AF", -1, "u3")}
    f1 = assign_folds(unit, sort, n_folds=2)
    f2 = assign_folds(unit, sort, n_folds=2)
    assert f1 == f2
    assert f1["a1"] == f1["a2"]  # 同じ単位は同じ fold
    assert f1 == {"a4": 0, "a1": 1, "a2": 1, "a3": 0}  # AF→EU(大)→EU(小) の順に巡回


def test_group_quotas_cap_each_place_and_allocate_proportionally() -> None:
    n = {"A": 9000, "B": 3000, "C": 2000, "D": 500}
    g = {"A": 1, "B": 1, "C": 2, "D": 2}
    q = group_quotas(n, g, cap=6000)
    assert q["A"] + q["B"] == 6000  # 大きな場所は合計 cap に
    assert (q["A"], q["B"]) == (4500, 1500)  # 評価セル数に比例
    assert (q["C"], q["D"]) == (2000, 500)  # cap 以下の場所はそのまま
    odd = group_quotas({"x": 1, "y": 1, "z": 1}, {"x": 0, "y": 0, "z": 0}, cap=2)
    assert sum(odd.values()) == 2  # 最大剰余法で合計がちょうど cap
    assert odd == group_quotas({"x": 1, "y": 1, "z": 1}, {"x": 0, "y": 0, "z": 0}, cap=2)


def test_placebo_permutes_within_aoi_and_is_reproducible() -> None:
    df = pd.DataFrame({
        "aoi": ["X"] * 3 + ["Y"] * 2,
        "P_gid1": [1, 2, 3, 4, 5], "P_gid2": [2, 3, 1, 5, 4], "P_gid3": [3, 1, 2, -1, -1],
    })
    look = {i: (f"n{i}", "IT") for i in range(1, 6)}
    a = placebo_lookups(df, look, seed=7)
    b = placebo_lookups(df, look, seed=7)
    assert a == b
    assert sorted(v[0] for v in a["X"].values()) == ["n1", "n2", "n3"]  # 名前集合は不変
    assert set(a["Y"]) == {4, 5}  # 他 AOI の名前は混ざらない


def test_nlearn_scores_shapes_and_signal() -> None:
    rng = np.random.default_rng(0)
    wet = ["riverside", "marshton", "fordham", "brookby"]
    dry = ["hilltop", "mountview", "highcrest", "rockridge"]
    names = np.array([rng.choice(wet if i % 2 else dry) for i in range(400)], dtype=object)
    y = np.array([1 if n in wet else 0 for n in names], dtype=np.int8)
    groups = np.array([f"g{i % 10}" for i in range(400)])
    te = np.array(["riverside", "hilltop", ""], dtype=object)
    oof, sc = nlearn_scores(names, y, groups, te)
    assert oof.shape == (400,)
    assert sc.shape == (3,)
    assert sc[0] > sc[1]  # 学習した名前の得点順序


@pytest.mark.parametrize(
    ("r", "c95", "c90", "kw", "label", "equiv"),
    [
        (0.08, [0.01, 0.15], [0.02, 0.14], {}, "positive-meaningful", False),
        (0.03, [0.001, 0.08], [0.005, 0.07], {}, "positive-small", False),
        # 0 より大きいが SESOI 未満と同等性の両方が成立 = 統計的には正・実務的には無視できる
        (0.02, [0.001, 0.04], [0.003, 0.035], {}, "positive-negligible", True),
        (0.08, [0.01, 0.15], [0.02, 0.14], {"placebo_ci95": [-0.01, 0.05]}, "not-meaning",
         False),
        (0.08, [0.01, 0.15], [0.02, 0.14], {"placebo_ci95": [0.01, 0.05]},
         "positive-meaningful", False),
        (-0.02, [-0.04, -0.001], [-0.035, -0.004], {}, "negative", True),
        (0.0, [-0.04, 0.045], [-0.03, 0.04], {}, "equivalent-null", True),
        (0.0, [-0.2, 0.2], [-0.15, 0.15], {}, "inconclusive", False),
        (0.08, [0.01, 0.15], [0.02, 0.14], {"den_ci95": [-1e-4, 3e-3]},
         "undefined-denominator", False),
    ],
)
def test_decide_table(r: float, c95: list, c90: list, kw: dict, label: str, *,
                      equiv: bool) -> None:
    agg = _load("world_aggregate")
    d = agg.decide(r, c95, c90, **kw)
    assert d["label"] == label
    assert d["equivalent"] is equiv
    assert d["positive"] is label.startswith("positive")


def test_time_folds_exclude_earlier_floods_at_evaluation_places() -> None:
    ev = _load("world_eval")
    infos = {a: _info(tag=a, event_time=t) for a, t in
             (("A", "2019-05-01"), ("B", "2024-03-01"), ("C", "2020-01-01"), ("D", "2025-06-01"))}
    unit = {"A": "g1", "B": "g1", "C": "g2", "D": "g3"}  # A は B と同じ場所の以前の洪水
    f = ev.time_folds(sorted(infos), infos, "2024-01-01", unit)
    assert f == {"A": -1, "B": 0, "C": 1, "D": 0}
    g = ev.time_folds(sorted(infos), infos, ev.LLM_CUT, unit)
    assert g == {"A": 1, "B": 1, "C": 1, "D": 0}


def test_folds_for_uses_the_given_group_buffer() -> None:
    """§7-6: 緩衝幅は引数で効く（DataFrame.attrs に頼らない。amendment B-1）。"""
    ev = _load("world_eval")
    # A と B は約 16 km 離れる: 5 km 緩衝（10 km 以内で結合）では別の場所、
    # 10 km 緩衝（20 km 以内で結合）では同じ場所
    infos = {"A": _info(tag="A", bounds=(10.00, 45.00, 10.10, 45.10)),
             "B": _info(tag="B", bounds=(10.30, 45.00, 10.40, 45.10)),
             "C": _info(tag="C", bounds=(12.00, 45.00, 12.10, 45.10))}
    df = pd.DataFrame({"aoi": ["A", "A", "B", "C"], "cc": ["IT"] * 4})
    df.attrs["buffer_km"] = 5.0  # 古い経路の値が残っていても使われないこと
    ev.folds_for("primary", df, infos, {"IT": "EU"}, buffer_km=10.0)
    u = df.groupby("aoi")["unit"].first()
    assert u["A"] == u["B"] != u["C"]
    ev.folds_for("primary", df, infos, {"IT": "EU"}, buffer_km=5.0)
    u = df.groupby("aoi")["unit"].first()
    assert len({u["A"], u["B"], u["C"]}) == 3


def test_assign_folds_seed_keeps_balance_and_changes_partition() -> None:
    units = {f"g{k:02d}": f"g{k:02d}" for k in range(23)}
    sort = {u: (u,) for u in units}
    base = assign_folds(units, sort)
    a = assign_folds(units, sort, seed=1)
    b = assign_folds(units, sort, seed=1)
    assert a == b  # 固定シードで決定的
    assert a != base  # 分け方は変わる
    for start in range(0, 23, 10):  # 各区画（大きさの帯）から各 fold へ高々 1 つ
        block = [f"g{k:02d}" for k in range(start, min(start + 10, 23))]
        assert len({a[u] for u in block}) == len(block)


def test_merge_by_activation_joins_groups_of_the_same_event() -> None:
    grp = {"E1_A": 0, "E1_B": 3, "E2_A": 1, "E3_A": 2, "E3_B": 2}
    act = {"E1_A": "E1", "E1_B": "E1", "E2_A": "E2", "E3_A": "E3", "E3_B": "E3"}
    m = merge_by_activation(grp, act)
    assert m["E1_A"] == m["E1_B"]  # 同じ洪水の別の場所は同じグループに
    assert len({m["E1_A"], m["E2_A"], m["E3_A"]}) == 3
    assert m["E3_A"] == m["E3_B"]


def test_recalibration_only_stacks_exist() -> None:
    ev = _load("world_eval")
    assert ev.STACK["M0s"] == ("M0", [])
    assert ev.STACK["T0s"] == ("T0", [])


def test_country_units_keep_cross_border_places_together() -> None:
    ev = _load("world_eval")
    aois = ["a", "b", "c"]
    grp = {"a": 0, "b": 0, "c": 1}  # a と b は国境をまたぐ同じ場所
    cc = {"a": "DE", "b": "AT", "c": "AT"}
    n = {"a": 300, "b": 100, "c": 50}
    assert ev.country_units(aois, grp, cc, n) == {"a": "DE", "b": "DE", "c": "AT"}


def test_placebo_perms_match_placebo_lookups() -> None:
    df = pd.DataFrame({
        "aoi": ["X"] * 3 + ["Y"] * 2,
        "P_gid1": [1, 2, 3, 4, 5], "P_gid2": [2, 3, 1, 5, 4], "P_gid3": [3, 1, 2, -1, -1],
    })
    look = {i: (f"n{i}", "IT") for i in range(1, 6)}
    perms = placebo_perms(df, seed=7)
    lk = placebo_lookups(df, look, seed=7)
    assert {a: {g: look[p] for g, p in m.items()} for a, m in perms.items()} == lk


def test_nread_block_features_and_placebo() -> None:
    df = pd.DataFrame({
        "aoi": ["X", "X"], "P_gid1": [1, 2], "P_gid2": [2, 1], "P_gid3": [-1, -1],
        "P_d1": [0.0, 1000.0], "P_d2": [1000.0, 3000.0], "P_d3": [np.inf, np.inf],
    })
    tab = pd.DataFrame({"gid": [1, 2, 9], "category": ["crossing", "highland", "water"],
                        "confidence": [0.9, 0.6, 0.5], "parsed_ok": [True, True, False]})
    nr = nread_lookup(tab)
    assert set(nr) == {1, 2}  # 解釈できなかった応答は登録しない
    f = nread_block(df, nr)
    assert list(f.columns) == NREAD_FEATURES
    assert f.loc[0, "nr_crossing_n1"] == 1.0
    assert f.loc[0, "nr_conf_n1"] == pytest.approx(0.9)
    assert f.loc[0, "nr_highland_k"] == pytest.approx(0.5)  # 2 番目の名前は 1 km 先で重み 1/2
    assert f.loc[1, "nr_crossing_k"] == pytest.approx(0.25)  # 3 km 先で重み 1/4
    swap = {"X": {1: 2, 2: 1}}  # 置換すると最近傍名の判定が入れ替わる
    fp = nread_block(df, nr, perms=swap, prefix="p_")
    assert fp.loc[0, "p_nr_highland_n1"] == 1.0
    assert fp.loc[1, "p_nr_crossing_n1"] == 1.0


def test_bootstrap_bs_recovers_pooled_mean() -> None:
    agg = _load("world_aggregate")
    s = np.array([[1.0, 2.0], [3.0, 2.0], [2.0, 2.0]])
    n = np.array([10.0, 10.0, 10.0])
    boot = agg.bootstrap_bs(s, n, b=2000, seed=1)
    assert boot.shape == (2000, 2)
    assert boot[:, 0].mean() == pytest.approx(0.2, abs=0.01)
    assert np.allclose(boot[:, 1], 0.2)


# ---- N-read（プロンプトと応答の解釈） ----------------------------------------------


def test_nread_prompt_and_parse() -> None:
    p = user_prompt("Frankfurt", "Frankfurt", "Germany")
    assert '"Frankfurt"' in p
    assert "local spelling" not in p  # 現地表記が同じなら併記しない
    assert "local spelling" in user_prompt("Muenchen", "München", "Germany")
    assert len(prompt_sha256()) == 64
    assert parse('{"category": "crossing", "confidence": 0.9}') == ("crossing", 0.9, True)
    assert parse('Sure: {"category": "Wetland", "confidence": 2}') == ("wetland", 1.0, True)
    assert parse("no json here") == ("unknown", 0.0, False)
    assert parse('{"category": "volcano", "confidence": 0.5}') == ("unknown", 0.0, False)
    assert "unknown" in NREAD_CATEGORIES


def test_vectorised_lexicon_features_match_reference_loop() -> None:
    """一意 gid で配る実装が、セル × K の素朴なループと同一の値を返す。"""
    from topsem.config import REPO_ROOT  # noqa: PLC0415
    from topsem.world.lexicon import (  # noqa: PLC0415
        CATEGORIES,
        Lexicon,
        _match_cached,
        cell_lexicon_features,
    )

    lex = Lexicon(REPO_ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    look = {1: ("Frankfurt", "DE"), 2: ("Oberhausen", "DE"), 3: ("Oxford", "GB"),
            4: ("Riverside", "US"), 5: ("Zzyx", "XX")}
    rng = np.random.default_rng(3)
    gids = rng.choice([-1, 1, 2, 3, 4, 5, 9], size=(200, 3))
    dists = rng.uniform(0, 5000, size=(200, 3))
    got = cell_lexicon_features(gids, dists, look, lex)
    cat = np.zeros((200, 3, len(CATEGORIES)), dtype=np.float32)
    for i in range(200):
        for j in range(3):
            g = int(gids[i, j])
            if g >= 0:
                name, cc = look.get(g, ("", ""))
                cat[i, j] = _match_cached(lex, name, cc)
    w = 1.0 / (1.0 + dists / 1000.0)
    for c_i, c in enumerate(CATEGORIES):
        assert np.array_equal(got[f"lex_{c}_n1"], cat[:, 0, c_i])
        assert np.allclose(got[f"lex_{c}_k"], (cat[:, :, c_i] * w).sum(axis=1))


def test_nlearn_score_on_unique_names_is_identical() -> None:
    from topsem.world.channels import _VEC, _fit_name_model, _score  # noqa: PLC0415

    names = np.array(["riverside", "hilltop", "fordham", "riverside", ""] * 40, dtype=object)
    y = np.array([1, 0, 1, 1, 0] * 40, dtype=np.int8)
    model = _fit_name_model(names, y)
    ref = model.decision_function(_VEC.transform(names)).astype(np.float32)
    assert np.array_equal(_score(model, names), ref)


def test_stage2_is_nested_recalibration_plus_names() -> None:
    """二段目: 名前が空なら基礎予測の単調な再較正、情報のある名前特徴を足せば評価側で改善する。"""
    ev = _load("world_eval")
    rng = np.random.default_rng(5)
    n = 20_000
    base = rng.normal(size=n)  # 基礎モデルのロジット
    name = (rng.random(n) < 0.1).astype(float)  # 情報のある名前特徴
    p = 1 / (1 + np.exp(-(-2.5 + base + 1.5 * name)))
    y = (rng.random(n) < p).astype(np.int8)
    tr, te = slice(0, n // 2), slice(n // 2, n)
    recal = ev.stage2(base[tr], base[te], np.empty((n // 2, 0)), np.empty((n // 2, 0)), y[tr])
    order = np.argsort(base[te])
    assert np.all(np.diff(recal[order]) >= -1e-7)  # 単調（順位は基礎予測のまま）
    named = ev.stage2(base[tr], base[te], name[tr, None], name[te, None], y[tr])
    bs = lambda q: float(np.mean((q - y[te]) ** 2))  # noqa: E731
    assert bs(named) < bs(recal)


def test_crossfit_base_partitions_training_side() -> None:
    ev = _load("world_eval")
    rng = np.random.default_rng(1)
    n = 3000
    df = pd.DataFrame({"x": rng.normal(size=n), "y": (rng.random(n) < 0.3).astype(np.int8),
                       "g": rng.integers(0, 40, size=n)})
    tr, te = df.iloc[:2500], df.iloc[2500:]
    from topsem.world.channels import inner_folds  # noqa: PLC0415

    folds = inner_folds(tr["g"].to_numpy())
    assert set(np.unique(folds)) == {0, 1, 2, 3, 4}
    # 同じ場所グループは同じ内側 fold（場所の暗記を内側でも防ぐ）
    assert tr.assign(f=folds).groupby("g")["f"].nunique().max() == 1
    oof, p_te = ev.crossfit_base(tr, te, ["x"], [42], folds)
    assert oof.shape == (2500,)
    assert p_te.shape == (500,)
    assert np.all((oof > 0) & (oof < 1))


def test_evaluate_and_diagnostics_on_random_predictions(tmp_path: Path) -> None:
    """乱数の予測で集計の評価経路（判定・梯子・層別・二方向・診断図）が最後まで通る。"""
    agg = _load("world_aggregate")
    rng = np.random.default_rng(11)
    n = 6000
    y = (rng.random(n) < 0.06).astype(np.int8)
    base = np.clip(0.06 + 0.05 * (y - 0.06) + rng.normal(0, 0.02, n), 1e-3, 0.999)
    arms = ["M0", "M1", "M1c", "M3n", "M3np", "M3nx", "M3nl", "M0c", "M0n", "M1s", "M3s", "M3ps",
            "T0", "M3", "M3p", *(f"T{k}c" for k in range(4)), *(f"T{k}n" for k in range(4))]
    d = {"y": y, "unit": np.array([f"g{i % 40:03d}" for i in range(n)]),
         "aoi": np.array([f"EMSR{i % 25:03d}_AOI01" for i in range(n)]),
         "cc": np.array(["DE", "FR", "KE"])[rng.integers(0, 3, n)],
         "s_lex_any": rng.integers(0, 2, n).astype(np.int8),
         "s_wetside": rng.integers(0, 2, n).astype(np.int8),
         "s_slope": rng.random(n).astype(np.float32), "s_built": rng.random(n).astype(np.float32),
         "s_hand": rng.random(n).astype(np.float32), "lon": rng.random(n).astype(np.float32),
         "lat": rng.random(n).astype(np.float32)}
    for a in arms:
        d[f"p_{a}"] = np.clip(base + rng.normal(0, 0.005, n), 1e-3, 0.999).astype(np.float32)
    d["p_M0"] = np.full(n, 0.06, dtype=np.float32)  # 位置だけ = ほぼ定数
    res = agg.evaluate(d, {"DE": "EU", "FR": "EU", "KE": "AF"})
    assert res["primary"]["comparison"] == "M1c->M3n"
    assert res["total"]["comparison"] == "M1s->M3s"
    assert res["flexible"]["comparison"] == "M1->M3"
    assert res["ladder"]["rungs"][-1] == "M1c->M3n"
    # 固定マージン法: R の区間は ΔBS の区間を地形の増分の点推定で割ったもの
    p = res["primary"]
    assert p["ci95_R"][0] == pytest.approx(p["ci95_delta"][0] / p["delta_bs_M0_M1"])
    assert "nearest_name_in_lexicon" in res["H5_strata"]
    preds = {k[2:]: v for k, v in d.items() if k.startswith("p_")}
    out = tmp_path / "diag.png"
    agg.diagnostics(d, preds, ("M1c", "M3n"), out)
    assert out.is_file()
    assert out.stat().st_size > 10_000


def test_composition_is_aoi_mean_of_name_features() -> None:
    ev = _load("world_eval")
    cols = {f: [0.0, 1.0, 1.0, 0.0] for f in ev.NAMES}
    tr = pd.DataFrame({"aoi": ["A", "A", "B", "B"], **cols, "y": [1, 0, 0, 0]})
    te = pd.DataFrame({"aoi": ["C", "C"], **{f: [1.0, 1.0] for f in ev.NAMES}})
    ev.add_composition(tr, te)
    assert tr["c_nlearn"].tolist() == [0.5, 0.5, 0.5, 0.5]
    assert te["c_lex_water_n1"].tolist() == [1.0, 1.0]
    assert set(ev.COMP) <= set(tr.columns)


def test_demeaned_name_features_have_zero_aoi_mean_and_no_composition() -> None:
    ev = _load("world_eval")
    feats = {f: [0.0, 1.0, 3.0, 1.0] for f in (*ev.NAMES, *ev.NAMES_P)}
    tr = pd.DataFrame({"aoi": ["A", "A", "B", "B"], **feats, "y": [1, 0, 0, 0]})
    te = pd.DataFrame({"aoi": ["C", "C"], **{f: [2.0, 2.0] for f in (*ev.NAMES, *ev.NAMES_P)}})
    ev.add_composition(tr, te)
    assert tr["d_nlearn"].tolist() == [-0.5, 0.5, 1.0, -1.0]
    assert te["d_lex_water_n1"].tolist() == [0.0, 0.0]  # AOI 内で一定なら情報は残らない
    assert set(ev.DEMEANED) | set(ev.DEMEANED_P) <= set(tr.columns)
    assert ev.STACK["M3d"] == ("M1", ev.DEMEANED)
    assert not set(ev.STACK["M3d"][1]) & set(ev.COMP)


def test_local_synthetic_effect_only_near_settlements() -> None:
    ev = _load("world_eval")
    rng = np.random.default_rng(0)
    n = 4000
    df = pd.DataFrame({"aoi": np.repeat(["A", "B"], n // 2), "y": rng.integers(0, 2, n),
                       "hand_min": rng.random(n), "slope_mean": rng.random(n),
                       "ldr_all": rng.random(n), "cci_crop": rng.random(n),
                       "lon": rng.random(n), "lat": rng.random(n),
                       "lex_wetside_n1": rng.integers(0, 2, n).astype(float),
                       "P_d1": np.where(np.arange(n) % 2 == 0, 100.0, 5000.0)})
    base = ev.synth_labels(df, 0.0, 0)
    near = ev.synth_labels(df, 3.0, 0, local_m=500.0)
    far = df["P_d1"].to_numpy() > 500
    wet = df["lex_wetside_n1"].to_numpy() > 0
    # 近い水系側のセルだけが浸水しやすくなり、遠いセルでは名前による差が出ない
    assert near[~far & wet].mean() > base[~far & wet].mean() + 0.05
    assert abs(near[far & wet].mean() - near[far & ~wet].mean()) < 0.1


def test_all_types_labels_and_inclusion() -> None:
    from topsem.world.dataset import riverine_labels  # noqa: PLC0415

    df = pd.DataFrame({"f_riverine": [0.9, 0.0, 0.0], "f_other": [0.0, 0.0, 0.0],
                       "f_flood_any": [0.9, 0.7, 0.1]})
    y, keep = riverine_labels(df)
    assert y.tolist() == [1, 0, 0]
    assert keep.tolist() == [True, False, True]  # 他種別だけの浸水は主解析から除外
    y, keep = riverine_labels(df, all_types=True)
    assert y.tolist() == [1, 1, 0]
    assert keep.all()
    coastal = _info(riverine_share=0.1, all_types=True)
    assert include(coastal, n_names_inside=3) == (True, "ok")  # §7-1 は種別の基準を外す


def test_quotas_for_sensitivity_modes() -> None:
    ev = _load("world_eval")
    dirs = [Path("A"), Path("B"), Path("C")]
    infos = {"A": _info(tag="A", n_eval=9000, bounds=(0.0, 0.0, 0.1, 0.1)),
             "B": _info(tag="B", n_eval=3000, bounds=(0.05, 0.05, 0.15, 0.15)),
             "C": _info(tag="C", n_eval=500, bounds=(5.0, 5.0, 5.1, 5.1))}
    q = ev.quotas_for(dirs, infos, ev.Sensitivity(group_cap=6000))
    assert q["A"] + q["B"] == 6000  # A と B は同じ場所グループ
    assert q["C"] == 500
    assert ev.quotas_for(dirs, infos, ev.Sensitivity(group_cap=0)) == {"A": 9000, "B": 3000,
                                                                       "C": 500}
    assert ev.quotas_for(dirs, infos, ev.Sensitivity(aoi_cap=2000)) == {"A": 2000, "B": 2000,
                                                                        "C": 500}


def test_real_label_evaluation_requires_frozen_protocol(tmp_path: Path) -> None:
    ev = _load("world_eval")
    draft = tmp_path / "p.md"
    draft.write_text("# x\n\n状態: **DRAFT**（凍結時にこの行を「FROZEN <UTC> / commit <sha>」"
                     "に置き換え）\n", encoding="utf-8")
    assert not ev.protocol_frozen(draft)  # 説明文中の FROZEN では凍結とみなさない
    frozen = tmp_path / "f.md"
    frozen.write_text("# x\n\n状態: **FROZEN** 2026-10-26T01:02:03Z / commit abc1234\n",
                      encoding="utf-8")
    assert ev.protocol_frozen(frozen)


def test_run_fold_all_arms_on_random_table(monkeypatch: pytest.MonkeyPatch) -> None:
    """本番の経路（柔軟腕 + 全ての入れ子腕 + N-read）を乱数の小さな表で 1 fold 回す。

    経路の確認が目的なので、HistGB の反復は 5 回に差し替える（本番の設定は変えない）。
    """
    ev = _load("world_eval")
    hgb = ev.HistGradientBoostingClassifier
    monkeypatch.setattr(ev, "HistGradientBoostingClassifier",
                        lambda **kw: hgb(**{**kw, "max_iter": 5}))
    rng = np.random.default_rng(7)
    n = 1200
    feats = {c for a in ev.ARMS.values() for c in a} | {
        c for _, extra in ev.STACK.values() for c in extra if not c.startswith("c_")}
    df = pd.DataFrame({c: rng.normal(size=n).astype(np.float32) for c in sorted(feats)})
    df["aoi"] = [f"EMSR{i % 12:03d}_AOI01" for i in range(n)]
    df["unit"] = [f"g{i % 12:03d}" for i in range(n)]
    df["fold"] = [i % 12 % 3 for i in range(n)]
    df["y"] = (rng.random(n) < 0.2).astype(np.int8)
    names = np.array(["riverside", "hilltop", "fordham", "marshby"], dtype=object)
    df["name_real"] = names[rng.integers(0, 4, n)]
    df["name_plac"] = names[rng.integers(0, 4, n)]
    arms = [a for a in ev.ARMS if a not in ("M3r", "M3rp")] + ["M3r", "M3rp"]
    r = ev.run_fold(df, 0, arms, [42], list(ev.STACK))
    assert set(r["preds"]) >= set(arms) | set(ev.STACK)
    for v in r["preds"].values():
        assert v.shape == (int((df["fold"] == 0).sum()),)
        assert np.all((v > 0) & (v < 1))


def test_lexicon_country_keys_are_iso_codes_and_norway_is_read() -> None:
    """B-4: 引用符なしの NO（ノルウェー）が YAML 1.1 の false にならず、意図した言語に対応。"""
    repo = Path(__file__).resolve().parents[4]
    path = repo / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
    lex = Lexicon(path)
    assert all(len(cc) == 2 and cc.isalpha() and cc.isupper() for cc in lex.country_langs)
    assert "FALSE" not in lex.country_langs
    assert lex.country_langs["NO"] == ("no",)
    assert lex.coverage("NO")
    assert sum(lex.match("Kloftefoss", "NO").values()) > 0


def test_label_variants_default_unchanged_and_variants_shift_labels() -> None:
    """B-5 の 7: 既定は凍結時と同じ。notrace は痕跡を除き、nonwater は水域を除いた面積で測る。"""
    df = pd.DataFrame({"f_riverine": np.array([0.6, 0.6, 0.3, 0.0], dtype=np.float32),
                       "f_flood_any": np.array([0.6, 0.6, 0.3, 0.6], dtype=np.float32),
                       "f_other": np.zeros(4, dtype=np.float32),
                       "f_trace": np.array([0.2, 0.0, 0.0, 0.0], dtype=np.float32),
                       "f_perm_water": np.array([0.0, 0.0, 0.45, 0.0], dtype=np.float32)})
    y, keep = riverine_labels(df)
    assert y.tolist() == [1, 1, 0, 0]
    assert keep.tolist() == [True, True, True, False]  # 河川以外の浸水だけのセルは除外
    y, _ = riverine_labels(df, variant="notrace")
    assert y.tolist() == [0, 1, 0, 0]  # 痕跡 0.2 を除くと 0.4 < 0.5
    y, _ = riverine_labels(df, variant="nonwater")
    assert y.tolist() == [1, 1, 1, 0]  # 0.3 / (1 − 0.45) ≈ 0.55 ≥ 0.5
    with pytest.raises(ValueError, match="unknown label variant"):
        riverine_labels(df, variant="x")


def test_attach_oracle_maps_site_hand_and_fills_missing() -> None:
    ev = _load("world_eval")
    df = pd.DataFrame({"aoi": ["A", "A", "A", "B"], "P_gid1": [1, 2, 3, 4]})
    ev.attach_oracle(df, {1: 0.5, 2: 1.5})
    assert df["o_site_hand"].tolist() == [0.5, 1.5, 1.0, 1.0]  # 3 は A の中央値、4 は全体の中央値


def test_stage2_returns_standardised_coefficients() -> None:
    ev = _load("world_eval")
    rng = np.random.default_rng(0)
    z = rng.normal(size=400)
    x = rng.normal(size=(400, 2))
    y = (rng.random(400) < 1 / (1 + np.exp(-(z + 0.8 * x[:, 0])))).astype(int)
    cf: list = []
    p = ev.stage2(z, z, x, x, y, coefs=cf)
    assert p.shape == (400,)
    assert len(cf) == 4  # 切片、基礎予測のロジット、名前特徴 2 つ
    assert cf[2] > abs(cf[3])  # 効果を入れた特徴の係数が大きい
