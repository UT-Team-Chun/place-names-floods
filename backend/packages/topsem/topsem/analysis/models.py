"""腕（arm）梯子の HistGB 訓練 — pilot/本体解析で共通の単一実装。

HistGradientBoostingClassifier(max_iter=200, lr=0.1) は protocol §6.1 の
凍結済み既定値。seed だけ実行時に渡す。
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier

from topsem.analysis.cell_features import build_xy

ARMS = ("M0", "M1", "M2", "M3")


def arm_names(tr: pd.DataFrame) -> tuple[str, ...]:
    """データが dist_river を持つ（=A-4 発効）とき消融去腕 M1r と M1rN を追加。"""
    return ARMS + (("M1r", "M1rN") if "dist_river" in tr.columns else ())


def fit_arms(tr: pd.DataFrame, te: pd.DataFrame, seed: int) -> dict[str, np.ndarray]:
    """fold 内での 4 腕訓練＋ホールドアウト予測（対応比較可能にするため全部保存）。"""
    arms: dict[str, np.ndarray] = {}
    y_tr = tr["y"].to_numpy()
    for arm in arm_names(tr):
        model = HistGradientBoostingClassifier(random_state=seed, max_iter=200, learning_rate=0.1)
        model.fit(build_xy(tr, arm), y_tr)
        arms[arm] = model.predict_proba(build_xy(te, arm))[:, 1]
    return arms


def fit_arms_ensemble(
    tr: pd.DataFrame, te: pd.DataFrame, seeds: list[int]
) -> dict[str, np.ndarray]:
    """複数シードで予測確率を平均した腕（アルゴリズム揺らぎの抑制・protocol §7 起案中）。"""
    acc: dict[str, np.ndarray] = {}
    for seed in seeds:
        for arm, p in fit_arms(tr, te, seed).items():
            acc[arm] = acc.get(arm, np.zeros_like(p)) + p / len(seeds)
    return acc
