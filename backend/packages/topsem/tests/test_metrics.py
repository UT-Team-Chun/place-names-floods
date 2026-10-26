import numpy as np
import pytest

from topsem.evaluation.metrics import brier_score, log_loss, pr_auc, top_area_capture


def test_brier_known_values() -> None:
    y = np.array([1.0, 0.0, 1.0, 0.0])
    p = np.array([0.5, 0.5, 0.5, 0.5])
    assert brier_score(y, p) == pytest.approx(0.25)
    assert brier_score(y, y) == pytest.approx(0.0)


def test_log_loss_penalizes_confident_wrong() -> None:
    y = np.array([1.0, 0.0])
    assert log_loss(y, np.array([0.9, 0.9])) > log_loss(y, np.array([0.6, 0.4]))


def test_pr_auc_bounds() -> None:
    rng = np.random.default_rng(42)
    y = (rng.random(200) < 0.1).astype(float)
    assert pr_auc(y, y) == pytest.approx(1.0)
    shuffled = pr_auc(y, rng.permutation(y))
    assert 0.0 <= shuffled <= 1.0


def test_top_area_capture_perfect_and_floor() -> None:
    y = np.array([1.0] * 10 + [0.0] * 990)
    assert top_area_capture(y, y, frac=0.05) == pytest.approx(1.0)
    inverse = np.where(y > 0, 0.0, 1.0)
    assert top_area_capture(y, inverse, frac=0.01) == pytest.approx(0.0)


def test_top_area_capture_undefined_without_positives() -> None:
    with pytest.raises(ValueError, match="未定義"):
        top_area_capture(np.zeros(10), np.random.default_rng(0).random(10))


def test_top_area_capture_random_tiebreak_expectation() -> None:
    """全同値 → capture = k/n（ランダムタイブレークの期待値）。ソート順に依存しない。"""
    y = np.array([1.0, 1.0, 0.0, 0.0])
    p = np.array([0.1, 0.1, 0.1, 0.1])
    assert top_area_capture(y, p, frac=0.5) == pytest.approx(0.5)
    # 入替えても不変（permutation invariant）
    y2 = np.array([0.0, 1.0, 0.0, 1.0])
    assert top_area_capture(y2, p, frac=0.5) == pytest.approx(0.5)


def test_top_area_capture_no_ties_includes_kth_cell() -> None:
    """無 ties のとき k 番目のセルを含む（レビュー A-5 の回帰・2026-09-23 修正）。"""
    y = np.array([0.0, 1.0, 0.0, 0.0])
    p = np.array([0.9, 0.5, 0.3, 0.1])
    # k=2: 上位 {0.9(dry), 0.5(wet)} → capture = 1/1 = 1.0（旧実装は 0.0）
    assert top_area_capture(y, p, frac=0.5) == pytest.approx(1.0)
    # wet が上位 k 外なら 0.0
    y2 = np.array([0.0, 0.0, 1.0, 0.0])
    assert top_area_capture(y2, p, frac=0.5) == pytest.approx(0.0)


def test_top_area_capture_partial_tie() -> None:
    y = np.array([1.0, 0.0, 1.0, 0.0, 0.0])
    p = np.array([0.9, 0.5, 0.5, 0.5, 0.1])
    # k=2: 上位確定 1 件(0.9) + tie 3 件中 1 スロット → 期待 = (1 + 1*1/3)/2
    assert top_area_capture(y, p, frac=0.4) == pytest.approx((1 + 1 / 3) / 2)
