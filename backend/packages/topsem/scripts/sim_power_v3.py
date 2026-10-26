# %% 検出力 v3（A-8 主判定基盤・r4' F4 処置）: 実出現表 + two-way 衝撃族 + 方式別 size。
# v2 からの変更:
#   (1) 出現表は data/labels/*/derived/*/metadata.json の実（event×watershed×n_cells）
#   (2) 節ゆらぎを event 共通衝撃 σ_e と節固有 σ_s に分解（grid を回す。v8 FOLD_SD は
#       σ_s 相当の上限として内部に保持）
#   (3) 判定方式を並列評価: product wild / basin wild / parametric shock /
#       simulator-anchored（H0 分布の 95% 分位を臨界値に使う＝方式ではなく較正）
# 使い方（backend/ から）::
#     uv run python packages/topsem/scripts/sim_power_v3.py [--trials 400]
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.config import LABELS_DIR
from topsem.runs import start_run, write_json

ALPHA = 0.05
SHOCK_SCENARIOS = [  # (σ_e, σ_s) — 実データ較正値と、より強い adversarial 設定
    ("iid", 0.0, 0.0),
    ("sec", 0.0, 3.7e-4),
    ("event-1", 1.85e-4, 3.2e-4),
    ("event-2", 3.7e-4, 2.6e-4),
]
DELTAS = [0.0, 5e-5, 1e-4, 2e-4, 5e-4]


def real_incidence() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ev_l, n_l = [], []
    for f in sorted(LABELS_DIR.glob("*/derived/*/metadata.json")):
        m = json.loads(f.read_text(encoding="utf-8"))
        if "n_window_cells" not in m:
            continue
        ev_l.append((m["event"], m["river"]))
        n_l.append(int(m["n_window_cells"]))
    events = sorted({e for e, _ in ev_l})
    basins = sorted({w for _, w in ev_l})
    eidx = np.array([events.index(e) for e, _ in ev_l])
    widx = np.array([basins.index(w) for _, w in ev_l])
    return eidx, widx, np.array(n_l, dtype=float), (len(events), len(basins))


def trial_rejections(  # noqa: PLR0913, PLR0917

    rng: np.random.Generator,
    eidx: np.ndarray,
    widx: np.ndarray,
    cnt: np.ndarray,
    delta: float,
    sig_e: float,
    sig_s: float,
    *,
    reps: int = 800,
) -> dict[str, bool]:
    """1 データ描画に対する 4 判定の棄却真偽（section-level 等価形）。"""
    n_s = len(cnt)
    share = cnt / cnt.sum()
    d_c = delta + rng.normal(0, sig_e, eidx.max() + 1)[eidx] + rng.normal(0, sig_s, n_s)
    obs = float(d_c @ share)
    # (a) product wild（v6 主）
    r = rng.choice([-1.0, 1.0], (reps, int(eidx.max()) + 1))
    s = rng.choice([-1.0, 1.0], (reps, int(widx.max()) + 1))
    stat_p = (r[:, eidx] * s[:, widx]) @ (d_c * share)
    # (b) basin wild
    stat_b = (s[:, widx] * (d_c * share)).sum(1)
    # (c) parametric shock（production two_way_shock_test 相当）:
    #     event raw 平均の r-flip + 節残差（interaction+basin+cell）の t-flip
    m_e = np.array([d_c[eidx == k].mean() if (eidx == k).any() else 0.0
                    for k in range(int(eidx.max()) + 1)])
    sh_e = np.array([share[eidx == k].sum() for k in range(len(m_e))])
    z = d_c - m_e[eidx]
    t = rng.choice([-1.0, 1.0], (reps, n_s))
    stat_a = ((rng.choice([-1.0, 1.0], (reps, len(m_e))) * (sh_e * m_e)).sum(1)
              + (t * (share * z)).sum(1))
    # (d) simulator-anchored: H0（delta=0）の同衝撃帯域は main 側で事前計算済み临界值を使う
    return {
        "prod_wild": float((np.abs(stat_p) >= abs(obs)).mean()) < ALPHA,
        "basin_wild": float((np.abs(stat_b) >= abs(obs)).mean()) < ALPHA,
        "param_shock": float((np.abs(stat_a) >= abs(obs)).mean()) < ALPHA,
        "_obs": obs,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--tag", default="power_v3")
    args = ap.parse_args()
    run = start_run(args.tag)
    eidx, widx, cnt, dims = real_incidence()
    print(f"incidence: events={dims[0]} basins={dims[1]} sections={len(cnt)} "
          f"cells={int(cnt.sum())}")
    # simulator-anchored 臨界値: 各衝撃設定で delta=0 の |obs| 95% 分位
    rows = []
    for name, se, ss in SHOCK_SCENARIOS:
        rng = np.random.default_rng(hash(name) % 2**31)
        d0 = (rng.normal(0, se, (args.trials, dims[0]))[:, eidx]
              + rng.normal(0, ss, (args.trials, len(cnt))))
        crit = float(np.quantile(np.abs(d0 @ (cnt / cnt.sum())), 1 - ALPHA))
        for delta in DELTAS:
            rej = {"prod_wild": 0, "basin_wild": 0, "param_shock": 0, "anchored": 0}
            r2 = np.random.default_rng(int(delta * 1e9) + 7)
            for _ in range(args.trials):
                out = trial_rejections(r2, eidx, widx, cnt, delta, se, ss)
                for k in ("prod_wild", "basin_wild", "param_shock"):
                    rej[k] += out[k]
                rej["anchored"] += abs(out["_obs"]) > crit
            row = {
                "shock": name, "sigma_event": se, "sigma_section": ss, "delta": delta,
                "crit_anchored": crit,
                **{k: round(v / args.trials, 3) for k, v in rej.items()},
                "trials": args.trials,
            }
            rows.append(row)
            print(f"{name:9s} δ={delta:.1e} crit={crit:.2e} " + " ".join(f"{k}={row[k]:.2f}" for k in ("prod_wild", "basin_wild", "param_shock", "anchored")))
    with (run.path("sim.csv")).open("w", encoding="utf-8") as f:
        import csv as _csv
        wtr = _csv.DictWriter(f, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)
    write_json(run.path("summary.json"), {
        **run.provenance(),
        "incidence": {"events": dims[0], "basins": dims[1], "sections": len(cnt),
                      "cells": [int(x) for x in cnt]},
        "note": "A-8 simulator-anchored 主判定。delta=0 行が size、>0 が power。",
    })
    print(f"run -> {run.path('')}")


if __name__ == "__main__":
    main()
