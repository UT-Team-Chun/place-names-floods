"""amendment B-13 の実装（topsem.world.b13）のテスト。"""

from __future__ import annotations

import numpy as np
import pytest
import yaml

from topsem.config import REPO_ROOT
from topsem.world.b13 import (
    QUALIFIERS,
    correct_lexicon_doc,
    event_before,
    fieller_interval,
    grouped_n_iter,
    l_pm,
    match_pairs,
    nadeau_bengio_var,
    nested_weights,
    qualifier_side,
    split_qualifier,
    tost_from_replicates,
)
from topsem.world.lexicon import Lexicon, _Yaml12Loader

LEX = REPO_ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"


def _doc() -> dict:
    return yaml.load(LEX.read_bytes(), Loader=_Yaml12Loader)  # noqa: S506


def test_corrected_lexicon_applies_each_rule(tmp_path) -> None:  # noqa: ANN001
    out, rec = correct_lexicon_doc(_doc())
    p = tmp_path / "lex.yaml"
    p.write_text(yaml.safe_dump(out, allow_unicode=True, sort_keys=False), encoding="utf-8")
    lex0, lex1 = Lexicon(LEX), Lexicon(p)
    # (b) -ach は -bach だけ、mont- の頭は外し語の mont は残す、kacha・Ried
    assert lex0.match("Birkach", "DE")["water"] == 1
    assert lex1.match("Birkach", "DE")["water"] == 0
    assert lex1.match("Steinbach", "DE")["water"] == 1
    assert lex0.match("Montreuil", "FR")["high"] == 1
    assert lex1.match("Montreuil", "FR")["high"] == 0
    assert lex1.match("Mont Saint Michel", "FR")["high"] == 1
    assert lex1.match("Kacha Pakka Mor", "PK")["riverine"] == 0
    assert lex1.match("Oberried", "DE")["wet"] == 0
    # (d) 相対の修飾語は低地・高地から外す、(e) 泉・井戸、(f) 早瀬・滝
    assert lex0.match("Niederhausen", "DE")["low"] == 1
    assert lex1.match("Niederhausen", "DE")["low"] == 0
    assert lex1.match("Ala-Kittila", "FI")["low"] == 0
    assert lex0.match("Chah Bahar", "IR")["water"] == 1
    assert lex1.match("Chah Bahar", "IR")["water"] == 0
    assert lex1.match("Storfors", "SE")["water"] == 0
    # (g) 足した形、(h) ベルギーとイタリアのドイツ語のブロック
    assert lex1.match("Bar-sur-Aube", "FR")["riverine"] == 1
    assert lex1.match("Usti nad Labem", "CZ")["riverine"] == 1
    assert lex1.match("Siofok", "HU")["water"] == 1
    assert lex0.match("Chateau", "BE")["riverine"] == 1
    assert lex1.match("Chateau", "BE")["riverine"] == 0
    assert lex1.country_langs["BE"] == ("nl", "fr")
    assert lex1.country_langs["IT"] == ("it",)
    assert lex1.country_langs["DE"] == lex0.country_langs["DE"]
    rules = {d["rule"] for d in rec["dropped"]}
    assert rules == {"a_v57_clean", "b_homograph", "c_riverine_named", "d_qualifier",
                     "e_spring_well", "f_rapids"}


def test_corrected_lexicon_strict_catches_missing_forms() -> None:
    entry = {"form": ["bach"], "pos": "suffix", "cat": "water"}
    doc = {"version": "x", "languages": {"de": {"entries": [entry]}}, "countries": {"DE": ["de"]}}
    with pytest.raises(ValueError, match="not found"):
        correct_lexicon_doc(doc)
    out, _ = correct_lexicon_doc(doc, strict=False)
    forms = [f for e in out["languages"]["de"]["entries"] for f in e["form"]]
    assert "bach" in forms
    assert "hamm" in forms


def test_split_qualifier_word_and_prefix_positions() -> None:
    side = qualifier_side(QUALIFIERS, ("de",))
    pre = ("unter", "nieder", "ober")
    assert split_qualifier("Oberhausen", side, pre) == ("upper", "hausen")
    assert split_qualifier("Niederhausen", side, pre) == ("lower", "hausen")
    assert split_qualifier("Unter Ach", side, pre) == ("lower", "ach")
    assert split_qualifier("Ober-Unter", side, pre) is None  # 修飾語が二つ
    assert split_qualifier("Obere", side, pre) is None  # 頭の後が 2 文字未満
    assert split_qualifier("Ober", side, pre) is None  # 基が無い
    en = qualifier_side(QUALIFIERS, ("en",))
    assert split_qualifier("Upper Heyford", en, ()) == ("upper", "heyford")
    assert split_qualifier("Low Bentham", en, ()) == ("lower", "bentham")


def test_match_pairs_is_greedy_by_distance_and_one_to_one() -> None:
    lo = np.array([[0.0, 0.0], [5000.0, 0.0]])
    hi = np.array([[100.0, 0.0], [20000.0, 0.0]])
    assert match_pairs(lo, hi, 10_000) == [(0, 0, 100.0)]
    assert match_pairs(lo, np.empty((0, 2)), 10_000) == []


def test_event_before_needs_more_than_the_gap() -> None:
    assert event_before("2020-01-01", "2020-03-01", 30)
    assert not event_before("2020-02-15", "2020-03-01", 30)
    assert not event_before("2020-05-01", "2020-03-01", 30)
    assert not event_before("", "2020-03-01", 30)


def test_grouped_n_iter_picks_a_count_within_range() -> None:
    rng = np.random.default_rng(0)
    n = 4000
    x = rng.normal(size=(n, 3)).astype(np.float32)
    y = (x[:, 0] + rng.normal(size=n) > 0).astype(int)
    groups = np.repeat(np.arange(40), n // 40).astype(str)
    k, info = grouped_n_iter(x, y, groups, seed=1, max_iter=50)
    assert 1 <= k <= 50
    assert info["n_val_groups"] == 4


def test_tost_and_fieller_agree_on_simple_cases() -> None:
    rng = np.random.default_rng(1)
    d = 1.0 + 0.01 * rng.normal(size=4000)
    e = 0.0 + 0.005 * rng.normal(size=4000)
    lp, lm = l_pm(e, d)
    t = tost_from_replicates(lp, lm)
    assert t["equivalent"]
    assert t["l_plus_upper"] < 0 < t["l_minus_lower"]
    ci = fieller_interval(0.0, 1.0, (0.005**2, 0.01**2, 0.0))
    assert ci is not None
    assert ci[0] < 0 < ci[1]
    assert abs(ci[0] + ci[1]) < 1e-9
    assert fieller_interval(0.0, 0.01, (1e-6, 1.0, 0.0)) is None  # 分母が 0 と区別できない
    e2 = 0.06 + 0.005 * rng.normal(size=4000)
    lp2, lm2 = l_pm(e2, d)
    assert not tost_from_replicates(lp2, lm2)["equivalent"]


def test_nested_weights_resample_activations_within_drawn_countries() -> None:
    rng = np.random.default_rng(2)
    country_of_act = np.array([0, 0, 0, 1, 2, 2])
    for _ in range(50):
        w = nested_weights(rng, country_of_act)
        assert w.shape == (6,)
        assert (w >= 0).all()
        for c in range(3):  # 国の中の重みの和は、その国の activation の数 × 国が引かれた回数
            s = w[country_of_act == c].sum()
            assert s % (country_of_act == c).sum() == 0


def test_nested_weights_draws_only_countries_with_activations() -> None:
    """B-13 の追記 6: activation の無い国の符号（ここでは 1 と 3）があっても止まらない。"""
    rng = np.random.default_rng(3)
    country_of_act = np.array([0, 0, 2, 4, 4])
    for _ in range(200):
        w = nested_weights(rng, country_of_act)
        # 引いた国の数は activation を持つ国の数（3）。
        # 国ごとの引かれた回数 = 重みの和 ÷ 国の activation の数
        times = sum(w[country_of_act == c].sum() / (country_of_act == c).sum() for c in (0, 2, 4))
        assert times == pytest.approx(3)


def test_nadeau_bengio_variance() -> None:
    v = np.array([0.0, 1.0, 2.0, 3.0])
    assert nadeau_bengio_var(v) == pytest.approx((1 / 4 + 1 / 9) * np.var(v, ddof=1))
