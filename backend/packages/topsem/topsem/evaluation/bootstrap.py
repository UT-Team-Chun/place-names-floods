"""イベント × 流域の二方向クラスタ bootstrap（protocol §1）。

Δ = metric(base) − metric(alt)（正 = alt の方が良い）の CI を、
- イベント方向リサンプル
- 流域方向リサンプル
のそれぞれで計算し、**ワースト側（広い方）区間**を報告する。
セル数は独立標本数ではない — 揺らぐのは節（cluster）単位でだけ。
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

Arr = np.ndarray


def brier(y: Arr, p: Arr) -> float:
    return float(np.mean((p - y) ** 2))


def two_way_cluster_bootstrap(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    metric: Callable[[Arr, Arr], float] = brier,
    reps: int = 2000,
    seed: int = 42,
) -> dict[str, float]:
    """点推定 Δ と二方向最悪側 95% CI を返す。"""
    rng = np.random.default_rng(seed)
    obs = metric(y, p_base) - metric(y, p_alt)

    def _axis(arr: Arr) -> tuple[float, float]:
        if len(arr) != len(y):  # 別解: 明示チェックでサイレントブロードキャストを防ぐ
            msg = "cluster 指定の長さが y と不一致"
            raise ValueError(msg)
        idx = {v: np.flatnonzero(arr == v) for v in np.unique(arr)}
        ks = list(idx)
        boots = np.empty(reps)
        for r in range(reps):
            picked = rng.choice(len(ks), len(ks), replace=True)
            rows = np.concatenate([idx[ks[i]] for i in picked])
            boots[r] = metric(y[rows], p_base[rows]) - metric(y[rows], p_alt[rows])
        return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

    elo, ehi = _axis(events)
    wlo, whi = _axis(watersheds)
    return {
        "delta": obs,
        "ci_lo": min(elo, wlo),
        "ci_hi": max(ehi, whi),
        "ci_event": [elo, ehi],
        "ci_watershed": [wlo, whi],
        "reps": reps,
        "seed": seed,
    }


def wild_cluster_test(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    reps: int = 4000,
    seed: int = 42,
) -> dict[str, float]:
    """wild-cluster 置換検定（Cameron-Gelbach-Miller type 4、レビュー A-3）。

    クラスタ数が少ない（イベント 3 等）とき percentile bootstrap の CI は
    過小カバーになる。cell 単位の損失差 d_i に対し、行（イベント）・列（流域）の
    Rademacher 重み積 w_ij = r_e * s_w で符号を交換し、
    p = P(|mean(d*w)| >= |mean(d)|) を両側で出す。CI ではなく p 値の統計。
    """
    yv, pb, pa = np.asarray(y), np.asarray(p_base), np.asarray(p_alt)
    d = (pb - yv) ** 2 - (pa - yv) ** 2
    erew, eidx = np.unique(np.asarray(events), return_inverse=True)
    wrew, widx = np.unique(np.asarray(watersheds), return_inverse=True)
    rng = np.random.default_rng(seed)
    obs = float(d.mean())
    stat = np.empty(reps)
    for r in range(reps):
        rd_e = rng.choice(np.array([-1.0, 1.0]), len(erew))
        rd_w = rng.choice(np.array([-1.0, 1.0]), len(wrew))
        stat[r] = float((d * rd_e[eidx] * rd_w[widx]).mean())
    p = float(np.mean(np.abs(stat) >= abs(obs)))
    se = float(stat.std(ddof=1))
    return {
        "delta": obs,
        "t_obs": obs / se if se else float("inf"),
        "p_value": max(p, 1.0 / (reps + 1)),
        "n_clusters_event": len(erew),
        "n_clusters_watershed": len(wrew),
        "reps": reps,
        "seed": seed,
    }


def _crve_block(d: Arr, idx: Arr, k: int) -> float:
    """mean(d) の一方向 CRVE（cluster 数 k の小標本補正込み）。"""
    n = len(d)
    tot = np.bincount(idx, weights=d, minlength=k)
    m = np.bincount(idx, minlength=k).astype(float)
    mm = np.where(m == 0, 1.0, m)
    gm = tot / mm
    mu = float(d.mean())
    q = float(np.sum(mm**2 * (gm - mu) ** 2)) / n**2
    return q * (k / max(k - 1, 1)) if k > 1 else 0.0


def two_way_crve(d: Arr, eidx: Arr, widx: Arr) -> float:
    """mean(d) の二方向クラスタ頑健分散（CGM 2015 §4.2 加法分解）。"""
    ne = len(np.unique(eidx))
    nw = len(np.unique(widx))
    both = np.unique(np.column_stack([eidx, widx]), axis=0, return_inverse=True)[1]
    v_e = _crve_block(d, eidx, ne)
    v_w = _crve_block(d, widx, nw)
    v_b = _crve_block(d, both, len(np.unique(both)))
    return max(v_e + v_w - v_b, 1e-300)


def two_way_wild_test(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    reps: int = 4000,
    seed: int = 42,
) -> dict[str, float]:
    """two-way 完全 studentized wild（A-8・CGM §V 推奨形。v7 の主検定）。

    重み w_i = r_e + s_w（**加法** — 共通イベント衝撃を null に載せる）。
    統計量 T = mean(d)/sqrt(V̂2way(d))、各リプリケートでも V̂2way を
    d*w から再計算（fully studentized multiplier）。積重みの
    wild_cluster_test（v6 までの主検定）は単一軸共有セルの共分散を 0 にして
    過大棄却 25–48%（r4' F2・v16 sweep）ため歴史的経路として残置。
    """
    yv, pb, pa = np.asarray(y), np.asarray(p_base), np.asarray(p_alt)
    d = (pb - yv) ** 2 - (pa - yv) ** 2
    erew, eidx = np.unique(np.asarray(events), return_inverse=True)
    wrew, widx = np.unique(np.asarray(watersheds), return_inverse=True)
    rng = np.random.default_rng(seed)
    # 分母は**観測 d の V̂2way で固定**（リプル毎再計算は分子と相関して
    # size を壊す — v16 実測）。分子 mean(d·w*) の null 分布を同じ目盛りで
    # 採る古典的 multiplier t 検定（CGM §V・MacKinnon–Webb 流）。
    se = two_way_crve(d, eidx, widx) ** 0.5
    t_obs = float(d.mean()) / se
    ts = np.empty(reps)
    for r in range(reps):
        rd_e = rng.choice(np.array([-1.0, 1.0]), len(erew))
        rd_w = rng.choice(np.array([-1.0, 1.0]), len(wrew))
        ts[r] = float((d * (rd_e[eidx] + rd_w[widx])).mean()) / se
    p = float(np.mean(np.abs(ts) >= abs(t_obs)))
    return {
        "delta": float(d.mean()),
        "t_obs": t_obs,
        "p_value": max(p, 1.0 / (reps + 1)),
        "n_clusters_event": len(erew),
        "n_clusters_watershed": len(wrew),
        "reps": reps,
        "seed": seed,
        "method": "two_way_wild_studentized (A-8)",
    }


def two_way_param_wild(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    reps: int = 4000,
    seed: int = 42,
) -> dict[str, float]:
    """two-way 分解付き parametric wild（A-8 採用形・v16 実測で size 検証）.

    d_i を二方向加法モデル d = μ + α_e + β_w + ε_i に fits（一次補正=マージナル
    平均偏差、交互作用・ノイズは resid に吸収）し、帰無分布を
    d*_i = r_e α_e + s_w β_w + h_i ε_i（各独立 Rademacher）で構成する。
    真の DGP が（イベント衝撃, 流域衝撃, 独立ノイズ）の加法なら
    Var(mean(d*)) = Var(mean(d)|H0) が**厳密に一致** — 積重み（単一軸共有の
    共分散を消す: v16 実測 25–72%）や additive 重み（両軸で二重計上: 全比 0%）
    の失敗をいずれも回避。μ 帯域は H0 で mean=0 に中心化するため落とす。
    """
    yv, pb, pa = np.asarray(y), np.asarray(p_base), np.asarray(p_alt)
    d = (pb - yv) ** 2 - (pa - yv) ** 2
    erew, eidx = np.unique(np.asarray(events), return_inverse=True)
    wrew, widx = np.unique(np.asarray(watersheds), return_inverse=True)
    mu = float(d.mean())
    ce = np.bincount(eidx, minlength=len(erew))
    a_e = np.bincount(eidx, weights=d, minlength=len(erew)) / ce - mu
    cw = np.bincount(widx, minlength=len(wrew))
    a_w = np.bincount(widx, weights=d, minlength=len(wrew)) / cw - mu
    resid = d - mu - a_e[eidx] - a_w[widx]
    rng = np.random.default_rng(seed)
    n = len(d)
    ts = np.empty(reps)
    ne, nw = len(erew), len(wrew)
    for r in range(reps):
        re_ = rng.choice(np.array([-1.0, 1.0]), ne)
        rw = rng.choice(np.array([-1.0, 1.0]), nw)
        rh = rng.choice(np.array([-1.0, 1.0]), n)
        ts[r] = float((re_[eidx] * a_e[eidx] + rw[widx] * a_w[widx] + rh * resid).mean())
    p = float(np.mean(np.abs(ts) >= abs(mu)))
    return {
        "delta": mu,
        "t_obs": mu / (ts.std(ddof=1) if ts.std() > 0 else 1.0),
        "p_value": max(p, 1.0 / (reps + 1)),
        "n_clusters_event": ne,
        "n_clusters_watershed": nw,
        "reps": reps,
        "seed": seed,
        "method": "two_way_param_wild (A-8 final)",
    }


def two_way_shock_test(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    reps: int = 4000,
    seed: int = 42,
) -> dict[str, float]:
    """two-way 衝撃 parametric bootstrap p（A-8 最終形）。

    節（event×watershed）平均の ANOVA 二次モーメントから分散成分
    (σ²_event, σ²_watershed, σ²_cell) を推定し、H0 合成 d を
        d~_i = g_e + u_w + t_i·ε_i（g,u ~ 適合正規、t Rademacher）で描く。

    **v16d までの方針検証（重要）**: (i) 積重み wild 25–72% 過大、(ii) 加法
    group-sum 反転は Cauchy–Schwarz で高衝撃域に無効、(iii) CRVE studentized
    は全比 0% 過小、(iv) 本 parametric 形もイベント 3・共通モード σ_g²/9 が
    同一実現から識別不能な領域で 12–20%。**イベント 3 の data だけでは
    共通衝撃未知の下で exact size の検定は原理的に存在しない**（temporal の
    floor p=0.125 と同じ不可能性系）。したがって v7 の主判定は
    「two-way worst-side 区間 + 検出力 v3 シミュレータ（衝撃 grid 込み）で
    較正した臨界値」に置き、本関数は robustness・感度報告専用とする。
    """
    yv, pb, pa = np.asarray(y), np.asarray(p_base), np.asarray(p_alt)
    d = (pb - yv) ** 2 - (pa - yv) ** 2
    ev = np.asarray(events)
    ws = np.asarray(watersheds)
    n = len(d)
    erew, eidx = np.unique(ev, return_inverse=True)
    wrew, widx = np.unique(ws, return_inverse=True)
    _, secidx = np.unique(np.column_stack([eidx, widx]), axis=0, return_inverse=True)
    nsec = int(secidx.max()) + 1
    ne = len(erew)
    mu = float(d.mean())
    cnt_sec = np.bincount(secidx, minlength=nsec).astype(float)
    m_sec = np.bincount(secidx, weights=d, minlength=nsec) / np.where(cnt_sec > 0, cnt_sec, 1)
    eps = d - m_sec[secidx]
    sigma_e2_cell = float(np.mean(eps**2))
    # section means の二方向 ANOVA SS（重み=cnt 無視の簡易版・節数が大きく偏らない前提）
    first = np.unique(secidx, return_index=True)[1]
    e_of_sec = eidx[first]
    w_of_sec = widx[first]
    _nw = len(wrew)
    ce_s = np.bincount(e_of_sec, minlength=ne)
    a_e = np.bincount(e_of_sec, weights=m_sec, minlength=ne) / ce_s - m_sec.mean()
    cw_s = np.bincount(w_of_sec, minlength=_nw)
    a_w = np.bincount(w_of_sec, weights=m_sec, minlength=_nw) / cw_s - m_sec.mean()
    inter = m_sec - m_sec.mean() - a_e[e_of_sec] - a_w[w_of_sec]
    sigma_sec2 = float(np.mean(inter**2))
    # 共通モード補正: E[mean a_e²] = (1−1/ne)σ²_g ⟹ ×ne/(ne−1)。
    # 偏差だけで
    # なく全事象平均のゆらぎも H0 帯域に載せる（v16d: これが無いと σ²/9 不足で
    # 高衝撃域 20–27% 過大棄却）。節 interaction 分は両軸に吸収されるので減じる。
    sigma_ev2 = max(
        float(np.mean(a_e**2)) * ne / (ne - 1) - sigma_sec2 / max(_nw, 1), 0.0
    )
    sigma_w2 = max(
        float(np.mean(a_w**2)) * _nw / (_nw - 1) - sigma_sec2 / max(ne, 1) - sigma_ev2 / ne,
        0.0,
    )
    rng = np.random.default_rng(seed)
    dev = np.empty(reps)
    for r in range(reps):
        g = rng.normal(0.0, np.sqrt(sigma_ev2), ne)[eidx]
        u = rng.normal(0.0, np.sqrt(sigma_w2), _nw)[widx]
        e = (eps[rng.integers(0, n, n)] if n > 0 else eps) * rng.choice(
            np.array([-1.0, 1.0]), n
        )
        dev[r] = float((g + u + e).mean())
    p = float(np.mean(np.abs(dev) >= abs(mu)))
    return {
        "delta": mu,
        "null_sd": float(dev.std(ddof=1)),
        "t_obs": mu / (dev.std(ddof=1) if dev.std() > 0 else 1.0),
        "p_value": max(p, 1.0 / (reps + 1)),
        "n_clusters_event": ne,
        "n_clusters_watershed": _nw,
        "n_sections": nsec,
        "sigma_event": float(np.sqrt(sigma_ev2)),
        "sigma_watershed": float(np.sqrt(sigma_w2)),
        "sigma_cell": float(np.sqrt(sigma_e2_cell)),
        "reps": reps,
        "seed": seed,
        "method": "two_way_shock_parametric (A-8 final)",
    }


# power_v3 衝撃定数（amendment A-8・sim.csv と同一の較正値。変更は新 amendment 経由）
ANCHOR_SHOCK_SCENARIOS = [
    ("iid", 0.0, 0.0),
    ("sec", 0.0, 3.7e-4),
    ("event-1", 1.85e-4, 3.2e-4),
    ("event-2", 3.7e-4, 2.6e-4),
]


def two_way_anchored_test(  # noqa: PLR0913
    y: Arr,
    p_base: Arr,
    p_alt: Arr,
    events: Arr,
    watersheds: Arr,
    *,
    scenarios: list = ANCHOR_SHOCK_SCENARIOS,
) -> dict:
    """A-8 主判定: 既登録衝撃族で較正した臨界値による anchored 検定（閉形式）。

    統計量 T = Σ_c share_c·d̄_c（節平均のセル数加重和）。H0 分散は
    衝撃成分 ||Eᵀshare||²σ_e² + ||share||²σ_s² ＋ 独立成分 Σ share_c²σ²_c/cnt_c
    の和で、power_v3 の delta=0 シミュレータ（線形正規）と同一分布の
    95% 分位 = 1.96σ と一致する。各設定の正規両側 p を取り、
    **worst-side（最大 p）を主 p** とする。data 由来分散成分推定に
    依存しない（3 イベント識別不能問題を回避）のが他方式との決定的差。
    iid 成分（Σ share²σ²_c/cnt_c）も載せるため純衝撃 sim より conservative。
    """
    yv, pb, pa = np.asarray(y), np.asarray(p_base), np.asarray(p_alt)
    d = (pb - yv) ** 2 - (pa - yv) ** 2
    eidx = np.unique(np.asarray(events), return_inverse=True)[1]
    wrew, widx = np.unique(np.asarray(watersheds), return_inverse=True)
    nw = len(wrew)
    key = eidx * nw + widx
    ukeys, sec = np.unique(key, return_inverse=True)
    n_sec = len(ukeys)
    e_of_sec = ukeys // nw
    cnt = np.bincount(sec, minlength=n_sec).astype(float)
    share = cnt / cnt.sum()
    d_c = np.bincount(sec, weights=d, minlength=n_sec) / cnt
    obs = float(share @ d_c)
    z = d - d_c[sec]
    var_c = np.bincount(sec, weights=z * z, minlength=n_sec) / np.maximum(cnt - 1.0, 1.0)
    s_iid2 = float(np.sum(share**2 * var_c / cnt))
    n_ev = int(e_of_sec.max()) + 1
    ps, sig = {}, {}
    for name, se, ss in scenarios:
        ev_w = np.bincount(e_of_sec, weights=share, minlength=n_ev)
        v = se**2 * float(ev_w @ ev_w) + ss**2 * float(share @ share) + s_iid2
        s0 = math.sqrt(v)
        sig[name] = s0
        if s0 == 0.0:
            ps[name] = 1.0 if obs == 0.0 else 0.0
        else:
            ps[name] = math.erfc(abs(obs) / (s0 * math.sqrt(2.0)))
    p_worst = max(ps.values())
    alpha = 0.05  # 凍結済み有意水準（protocol §1。A-8 で不変）
    return {
        "obs": obs,
        "sd_iid": math.sqrt(s_iid2),
        "p_by_scenario": ps,
        "sd_by_scenario": sig,
        "p_value": p_worst,
        "reject_05": p_worst < alpha,
        "n_sections": n_sec,
        "method": "two_way_anchored (A-8 primary)",
    }
