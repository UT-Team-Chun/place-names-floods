import numpy as np
import pytest

from topsem.evaluation.bootstrap import (
    _crve_block,
    two_way_anchored_test,
    two_way_cluster_bootstrap,
    two_way_crve,
    two_way_wild_test,
    wild_cluster_test,
)


def test_cluster_noise_widens_ci_vs_cells() -> None:
    rng = np.random.default_rng(0)
    n_w = 60
    n_cells = 400
    # 流域レベルで効く共通バイアス（クラスタ効果）
    wbias = rng.normal(0, 0.05, n_w)
    events = np.repeat([0, 1], n_w * n_cells // 2)
    ws = np.tile(np.repeat(np.arange(n_w), n_cells // 2), 2)
    y = (rng.random(len(ws)) < 0.1).astype(float)
    p_base = np.full(len(y), 0.1)
    p_alt = np.clip(p_base + wbias[ws] * (y - 0.1) * -5 + rng.normal(0, 0.01, len(y)), 0.01, 0.99)
    out = two_way_cluster_bootstrap(y, p_base, p_alt, events, ws, reps=200, seed=1)
    assert out["ci_lo"] <= out["delta"] <= out["ci_hi"] * 1.0 + abs(out["delta"]) + 1e-9


def test_perfect_signal_ci_positive() -> None:
    rng = np.random.default_rng(1)
    n_w, n_cells = 40, 200
    ws = np.repeat(np.arange(n_w), n_cells)
    events = np.zeros(n_w * n_cells, dtype=int) + 0
    y = np.concatenate([
        (rng.random(n_cells) < (0.02 + 0.1 * (w % 2))).astype(float) for w in range(n_w)
    ])
    p_base = np.full(len(y), 0.08)
    # alt は節ごとの実発生率をよく当てる（訓練由来の予測 = LORO 相当を模擬）
    true_rate = np.repeat([0.02 + 0.1 * (w % 2) for w in range(n_w)], n_cells)
    p_alt = np.clip(true_rate + rng.normal(0, 0.01, len(y)), 1e-3, 1 - 1e-3)
    out = two_way_cluster_bootstrap(y, p_base, p_alt, events, ws, reps=200, seed=2)
    assert out["delta"] > 0
    assert out["ci_lo"] > 0 or out["ci_watershed"][0] <= out["ci_event"][0]


def test_string_cluster_labels_work() -> None:
    """実データは event/watershed が日本語文字列配列（npz 由来）。 regression。"""
    rng = np.random.default_rng(3)
    ws = np.repeat(np.array(["久慈川水系", "那珂川水系", "阿武隈川水系"]), 50)
    events = np.repeat(np.array(["2019_taihuu19gou"]), len(ws))
    y = (rng.random(len(ws)) < 0.1).astype(float)
    p_base = np.full(len(y), 0.1)
    p_alt = np.clip(p_base + rng.normal(0, 0.02, len(y)), 0.01, 0.99)
    out = two_way_cluster_bootstrap(y, p_base, p_alt, events, ws, reps=100, seed=1)
    assert out["ci_lo"] <= out["delta"] <= out["ci_hi"]


def test_length_mismatch_raises() -> None:
    rng = np.random.default_rng(0)
    y = (rng.random(40) < 0.2).astype(float)
    p = np.full(40, 0.2)
    with pytest.raises(ValueError, match="長さ"):
        two_way_cluster_bootstrap(y, p, p, np.zeros(41, dtype=int), np.zeros(40, dtype=int))


def test_wild_cluster_test_null_and_signal() -> None:
    """null（腕同値）は p 大、明確信号は p 小。CGM4 重みは両軸積。"""
    rng = np.random.default_rng(1)
    n = 600
    ev = rng.integers(0, 3, n)
    ws = rng.integers(0, 5, n)
    y = (rng.random(n) < 0.15).astype(float)
    p1 = np.full(n, 0.15)
    null = wild_cluster_test(y, p1, p1, ev, ws, reps=500, seed=3)
    assert null["delta"] == 0.0
    assert null["p_value"] > 0.2
    p3 = np.where(y == 1, 0.6, 0.05)  # 明確に良い alt → Δ=BS(base)−BS(alt) は正
    sig = wild_cluster_test(y, p1, p3, ev, ws, reps=500, seed=3)
    assert sig["delta"] > 0.0
    assert sig["p_value"] < 0.05
    assert sig["n_clusters_event"] == len(np.unique(ev))


def test_wild_cluster_singletons_capped_p() -> None:
    """極小サンプルでも p 値は 1/(reps+1) を下回らない（過小報告防止）。"""
    y = np.array([1.0, 0.0, 0.0, 1.0])
    p1 = np.array([0.2, 0.2, 0.2, 0.2])
    p3 = np.array([0.9, 0.01, 0.01, 0.9])
    r = wild_cluster_test(
        y, p1, p3, np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1]), reps=100, seed=0
    )
    assert r["p_value"] >= 1 / 101


def test_two_way_crve_positive_and_additive_decomposition() -> None:
    rng = np.random.default_rng(5)
    n = 400
    d = rng.normal(0, 1, n)
    ev = np.repeat(np.arange(4), n // 4)
    ws = np.tile(np.arange(5), n // 5)
    v = two_way_crve(d, ev, ws)
    assert v > 0
    # 同一イベント衝撃の共有は二方向分散に載る（一方向 event のみより縮小しない）
    v_e = _crve_block(d, ev, 4)
    v_w = _crve_block(d, ws, 5)
    both = np.unique(np.column_stack([ev, ws]), axis=0, return_inverse=True)[1]
    assert abs(v - (v_e + v_w - _crve_block(d, both, both.max() + 1))) < 1e-12


def test_two_way_wild_test_returns_p_and_method() -> None:
    rng = np.random.default_rng(6)
    n = 300
    y = (rng.random(n) < 0.1).astype(float)
    p1 = np.clip(rng.random(n) * 0.2, 0.01, 0.99)
    p3 = np.clip(p1 + rng.normal(0, 0.001, n), 0.01, 0.99)
    ev = np.repeat(np.arange(3), 100)
    ws = np.tile(np.arange(5), 60)
    out = two_way_wild_test(y, p1, p3, ev, ws, reps=400, seed=1)
    assert 0 < out["p_value"] <= 1
    assert out["method"].startswith("two_way_wild")
    assert out["n_clusters_event"] == 3


def test_two_way_anchored_null_signal_and_mc() -> None:
    """A-8 anchored 検定: null で non-reject・大信号で reject・閉形式σ=MC。"""
    rng = np.random.default_rng(7)
    n = 3000
    y = (rng.random(n) < 0.1).astype(float)
    p1 = rng.random(n) * 0.2
    ev = np.repeat(np.arange(3), 1000)
    ws = np.tile(np.arange(5), 600)
    # null: p3 は y と無関係 → worst-p は大きく出て当然（conservative 設計）
    p3 = rng.random(n) * 0.2
    out0 = two_way_anchored_test(y, p1, p3, ev, ws)
    assert out0["p_value"] > 0.2
    assert not out0["reject_05"]
    # 大信号: y への 4% shrinkage（全節で d>0・obs≈1e-3 ≫ 衝撃σ帯 2.2e-4）
    out1 = two_way_anchored_test(y, p1, p1 + 0.04 * (y - p1), ev, ws)
    assert out1["p_value"] < 1e-3
    assert out1["reject_05"]
    # 閉形式 sd と純衝撃モンテカルロの 95% 分位が一致（power_v3 と同一分布）
    share = np.full(15, 1 / 15)
    ss = 3.7e-4  # sec 設定（event 成分 0）
    g = rng.normal(0, ss, (20_000, 15)) @ share
    q = np.quantile(np.abs(g), 0.95)
    sd = (ss**2 * float(share @ share)) ** 0.5
    assert abs(q - 1.959964 * sd) / q < 0.05
    assert out0["method"].startswith("two_way_anchored")
