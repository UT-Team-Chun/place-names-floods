# %% 検出力 v2（ROADMAP 設計確定分）: 疎な出現表・wild p 判定・モデル揺らぎ込み。
# 節レベル CGM4 の正確な同値形（セルでなく節の損失差 Δ_c に重み w_c を掛けた和）を
# ベクトル化シミュレーションで回す。cell 単位で回すと本地マシンで重すぎるための
# 等価変形 — 断り: 節内サイズ均一仮定を捨てるため w_c = n_c/N を明示的に使う。
# 使い方（backend/ から）: uv run python packages/topsem/scripts/sim_power_v2.py
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.runs import start_run, write_json

FOLD_SD = 3.7e-4  # v8: fold あたり ΔBS のシード間 SD（モデル当てはめ揺らぎの代理）
ALPHA = 0.05  # 検定水準（protocol §1 と同じ 5% 両側）
MIN_SECTION_SHARE = 0.5  # sparse 出現表で節が生まれる確率
MIN_PAIRS = 3


def simulate(  # noqa: PLR0913
    rng: np.random.Generator,
    n_events: int,
    n_basins: int,
    delta: float,
    mode: str,
    *,
    trials: int = 400,
    reps: int = 1000,
    cells: int = 20_000,
) -> tuple[float, int]:
    if mode == "share":  # 全イベントが全水系に現れる（2019 型）
        pairs = [(e, w) for e in range(n_events) for w in range(n_basins)]
    else:  # sparse — 実データ（3×10 で 13 節）のような薄い出現表
        pairs = [(e, w) for e in range(n_events) for w in range(n_basins)
                 if rng.random() < MIN_SECTION_SHARE]
    if len(pairs) < MIN_PAIRS:
        pairs = pairs or [(0, 0)]
    ev = np.array([e for e, _ in pairs])
    ws = np.array([w for _, w in pairs])
    nc = rng.integers(cells // 4, 4 * cells, len(pairs))  # 節サイズの実データ風散布
    share = nc / nc.sum()
    rej = 0
    for _ in range(trials):  # 外: データ描画（Δ_c）、内: 同データ上の重み null
        d_c = delta + rng.normal(0, FOLD_SD, len(pairs))
        obs = float(d_c @ share)
        r = rng.choice([-1.0, 1.0], size=(reps, n_events))
        s = rng.choice([-1.0, 1.0], size=(reps, n_basins))
        stat = (d_c * (r[:, ev] * s[:, ws]) * share).sum(1)
        p = (np.abs(stat) >= abs(obs)).mean()
        rej += p < ALPHA
    return rej / trials, len(pairs)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    run = start_run(args.tag or "power_v2")
    deltas = [0.0, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3]
    grids = [(3, 9, "share"), (3, 15, "share"), (5, 10, "sparse"), (5, 20, "share"),
             (8, 20, "sparse"), (10, 30, "share")]
    rows = []
    for e, w, mode in grids:
        for d in deltas:
            rng = np.random.default_rng(11 * e + 101 * w)
            power, npair = simulate(rng, e, w, d, mode, trials=args.trials)
            rows.append({"events": e, "basins": w, "incidence": mode, "n_sections": npair,
                         "delta": d, "power": round(power, 3), "trials": args.trials})
            print(f"E={e:2d} W={w:2d} {mode:6s} 節={npair:3d} δ={d:.1e} power={power:.2f}")
    with (run.path("sim.csv")).open("w", newline="", encoding="utf-8") as f:
        wtr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)
    write_json(run.path("summary.json"), {**run.provenance(), "fold_sd": FOLD_SD, "results": rows})
    print(f"→ {run.run_dir}")


if __name__ == "__main__":
    main()
