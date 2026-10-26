"""Pilot: 2019 東日本台風 GSI ラベル × 全国地名、leave-one-river-out。

protocol §1 の主比較（M3 vs M1）ではなく、**Gate 2 準備の探査的試験**として
M0（公開位置情報だけ）も対照に ΔBS を見る。ラベルは robustness 級
（浸水推定段彩図 — 地理院免責あり）である点を summary.json / ログに明記する。
データ組立は本体解析と共通（`topsem.analysis.cell_features`）。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/run_pilot_2019.py --tag pilot2019_v3
"""

from __future__ import annotations

import argparse
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from topsem.analysis.cell_features import (
    discover_label_units,
    load_population,
    river_cells,
)
from topsem.analysis.models import ARMS, fit_arms
from topsem.corpus.snapshot import load_snapshot
from topsem.evaluation.metrics import brier_score, log_loss, pr_auc, top_area_capture
from topsem.features.route_a import ElementIndex
from topsem.runs import start_run, write_json
from topsem.splits.manifest import load_split

EVENT = "2019_taihuu19gou"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--split-version", default="pilot-loro-2019")
    args = ap.parse_args()

    run = start_run(args.tag or datetime.now(UTC).strftime("pilot2019_%Y%m%dT%H%M%SZ"))
    snap = load_snapshot()
    index = ElementIndex.from_elements(snap.elements)
    units = [u for u in discover_label_units() if u["event"] == EVENT]
    if not units:
        msg = f"data/labels/{EVENT}/derived がない。make labels-2019 を先に実行"
        raise FileNotFoundError(msg)
    pop = load_population()
    feat_cache: dict = {}
    emb_cache: dict = {}
    data = pd.concat(
        [river_cells(u, pop, index, feat_cache, emb_cache) for u in units],
        ignore_index=True,
    )

    split = load_split(args.split_version)
    leak = split.validate(data["event"].tolist(), data["watershed"].tolist())
    if not leak["ok"]:
        msg = f"分割 manifest とデータが不一致: {leak['problems']}"
        raise RuntimeError(msg)

    preds_dir = run.path("predictions")
    preds_dir.mkdir(exist_ok=True)
    fold_rows = []
    for fold in split.folds:
        held = fold.holdout[0].watershed
        mask = np.array(split.fold_mask(fold, data["event"].tolist(), data["watershed"].tolist()))
        tr, te = data[~mask], data[mask]
        arms = fit_arms(tr, te, args.seed)
        y = te["y"].to_numpy()
        row = {"held_river": held, "n_cells": len(te), "wet_rate": float(y.mean())}
        for arm in ARMS:
            row[f"bs_{arm}"] = brier_score(y, arms[arm])
            row[f"ll_{arm}"] = log_loss(y, arms[arm])
            if y.sum() > 0:
                row[f"pr_auc_{arm}"] = pr_auc(y, arms[arm])
                row[f"cap10_{arm}"] = top_area_capture(y, arms[arm])
        row["delta_bs_M1_M0"] = row["bs_M0"] - row["bs_M1"]
        row["delta_bs_M2_M0"] = row["bs_M0"] - row["bs_M2"]
        row["delta_bs_M3_M0"] = row["bs_M0"] - row["bs_M3"]
        row["delta_bs_M3_M1"] = row["bs_M1"] - row["bs_M3"]
        fold_rows.append(row)
        np.savez_compressed(
            preds_dir / f"{held}.npz", **{f"p_{k}": v for k, v in arms.items()}, y=y
        )
        print(
            f"{held:34s} ΔBS(M3-M1)={row['delta_bs_M3_M1']:+.6f} "
            f"ΔBS(M1-M0)={row['delta_bs_M1_M0']:+.6f} n={row['n_cells']}"
        )

    deltas = np.array([r["delta_bs_M3_M1"] for r in fold_rows])
    deltas0 = np.array([r["delta_bs_M3_M0"] for r in fold_rows])
    summary = {
        **run.provenance(),
        "corpus_snapshot": snap.manifest(),
        "protocol_version": "v1-draft",
        "splits_version": args.split_version,
        "leakage_checks": leak,
        "label_role": "robustness_check_only (GSI 浸水推定段彩図; 免責: 推定混じり)",
        "arms": {
            "M0": "ix,iy",
            "M1": "+DEM(SRTM30) elev, 相対低地度（地形のみ）",
            "M2": "+route_a 辞書",
            "M3": "+route_b n-gram hashed embedding",
        },
        "seeds": [args.seed],
        "event": EVENT,
        "n_rivers": len(units),
        "n_cells_total": len(data),
        "delta_bs_M3_M1_by_river": {r["held_river"]: r["delta_bs_M3_M1"] for r in fold_rows},
        "delta_bs_M3_M0_by_river": {r["held_river"]: r["delta_bs_M3_M0"] for r in fold_rows},
        "delta_bs_mean": float(deltas.mean()),
        "delta_bs_mean_vs_m0": float(deltas0.mean()),
        "sign_positive": int((deltas > 0).sum()),
    }
    write_json(run.path("metrics.json"), {"folds": fold_rows})
    write_json(run.path("summary.json"), summary)
    print(
        f"\nΔBS(M3-M1) mean={deltas.mean():+.6f} positive {int((deltas > 0).sum())}/{len(deltas)}"
        f" | ΔBS(M3-M0) mean={deltas0.mean():+.6f}"
    )
    print(f"→ {run.run_dir}")


if __name__ == "__main__":
    main()
