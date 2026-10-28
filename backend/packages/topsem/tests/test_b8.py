"""事後 B-8 の部品（topsem.world.b8）のテスト。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from topsem.world.b8 import (
    clean_lexicon_doc,
    contains_hydronym,
    correct_misclassification,
    frequent_tokens,
    group_weights,
    groups_needed,
    hydronym_cores,
    jackknife_se,
    name_tokens,
    stage2_weighted,
)


def test_groups_needed_for_equivalence_and_positive() -> None:
    # 半幅 0.05、点推定 0.0、余白 0.05 → 半幅を 0.05 以下に: 今の数で足りる
    assert groups_needed(0.0, -0.05, 0.05, 100, margin=0.05) == pytest.approx(100.0)
    # 点推定 0.03 → 余地 0.02、半幅 0.08 → (0.08/0.02)^2 = 16 倍
    assert groups_needed(0.03, -0.05, 0.11, 195, margin=0.05) == pytest.approx(195 * 16)
    # 点推定が余白の外なら届かない
    assert math.isinf(groups_needed(0.06, 0.0, 0.12, 195, margin=0.05))
    # 正: 下限 > 0 には (半幅 / 点推定)^2
    assert groups_needed(0.02, -0.02, 0.06, 195) == pytest.approx(195 * 4)
    assert math.isinf(groups_needed(-0.01, -0.03, 0.01, 195))


def test_jackknife_se_matches_formula() -> None:
    assert jackknife_se([1.0, 2.0, 3.0]) == pytest.approx(math.sqrt(2 / 3 * 2.0))
    assert math.isnan(jackknife_se([1.0]))


def test_correct_misclassification_identity_and_example() -> None:
    np.testing.assert_allclose(correct_misclassification(1.3, 1.0, 0.0), 1.3)
    # 精度 0.76・他の名の 0.18 で 1.01 は 1 の近く、0.77 と 1.32 は外へ広がる
    v = correct_misclassification(np.array([0.77, 1.01, 1.32]), 0.76, 0.18)
    assert v[0] < 0.77
    assert v[2] > 1.32
    assert v[1] == pytest.approx(1.0173, abs=1e-3)
    # 誤分類で説明しきれない値は NaN
    assert np.isnan(correct_misclassification(0.2, 0.76, 0.18))


def test_clean_lexicon_doc_drops_short_affixes_and_english_suffixes_only() -> None:
    doc = {"version": "1.0", "countries": {"DE": ["de"], "GB": ["en"]}, "languages": {
        "de": {"name": "German", "entries": [
            {"form": ["au", "aue", "werder"], "pos": "suffix", "cat": "riverine"},
            {"form": ["aue", "au"], "pos": ["word", "suffix"], "cat": "riverine"}]},
        "en": {"name": "English", "entries": [
            {"form": ["burn", "brook", "well"], "pos": "suffix", "cat": "water"},
            {"form": ["burn", "well"], "pos": "word", "cat": "water"}]}}}
    out, dropped = clean_lexicon_doc(doc)
    de = out["languages"]["de"]["entries"]
    assert {"form": ["aue", "werder"], "pos": "suffix", "cat": "riverine"} in de
    assert {"form": ["aue", "au"], "pos": "word", "cat": "riverine"} in de  # 語の位置は残す
    assert {"form": ["aue"], "pos": "suffix", "cat": "riverine"} in de
    en = out["languages"]["en"]["entries"]
    assert {"form": ["brook"], "pos": "suffix", "cat": "water"} in en
    assert {"form": ["burn", "well"], "pos": "word", "cat": "water"} in en
    assert ("de", "au", "suffix") in dropped
    assert ("en", "burn", "suffix") in dropped
    assert doc["languages"]["de"]["entries"][0]["form"] == ["au", "aue", "werder"]  # 元のまま


def test_hydronym_cores_and_matching() -> None:
    assert name_tokens("Frankfurt (Oder)") == ["frankfurt", "oder"]
    generic = frequent_tokens(["River Main", "Neckar River", "Lake Inari", "Oder River"], top=1)
    assert generic == {"river"}
    assert hydronym_cores("River Main", generic, set()) == {"main"}
    assert hydronym_cores("Lake Inari", generic | {"lake"}, set()) == {"inari"}
    # 国の集落の名で多いトークンは除く
    assert hydronym_cores("New River", {"river"}, {"new"}) == set()
    assert contains_hydronym("Frankfurt am Main", {"main"})
    assert contains_hydronym("Neckarsulm", {"neckar"})
    assert not contains_hydronym("Mainz", {"main"})  # 残り 1 文字
    assert contains_hydronym("Ounasjoki", {"ounasjoki"})
    assert not contains_hydronym("Ennerdale Bridge", {"ehen"})


def test_group_weights_sum() -> None:
    w = group_weights(195, 20261028)
    assert w.sum() == 195
    assert len(w) == 195
    assert (w >= 0).all()


def test_stage2_weighted_reproduces_unweighted_and_replication() -> None:
    rng = np.random.default_rng(1)
    n = 600
    z = rng.normal(size=n)
    x = rng.normal(size=(n, 3))
    y = (rng.random(n) < 1 / (1 + np.exp(-(z + 0.5 * x[:, 0])))).astype(int)
    p0, _ = stage2_weighted(z[:400], z[400:], x[:400], x[400:], y[:400])
    # 整数の重みは行を複製したのと同じ
    w = rng.integers(0, 3, size=400).astype(float)
    rep = np.repeat(np.arange(400), w.astype(int))
    p1, coef = stage2_weighted(z[:400], z[400:], x[:400], x[400:], y[:400], w_tr=w)
    p2, _ = stage2_weighted(z[rep], z[400:], x[rep], x[400:], y[rep])
    np.testing.assert_allclose(p1, p2, atol=1e-6)
    # 初期値から始めても同じ収束先
    p3, _ = stage2_weighted(z[:400], z[400:], x[:400], x[400:], y[:400], w_tr=w, init=coef)
    np.testing.assert_allclose(p1, p3, atol=1e-5)
    assert p0.shape == (200,)
