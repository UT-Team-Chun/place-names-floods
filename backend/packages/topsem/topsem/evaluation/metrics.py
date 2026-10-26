"""主評価量（事前登録 docs/protocol/ に固定して使う指標群）。

y は {0,1}（観測マスクで used=True のセルのみで呼ぶこと）、p は確率予測。
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

_EPS = 1e-12


def brier_score(y: npt.NDArray[np.float64], p: npt.NDArray[np.float64]) -> float:
    return float(np.mean((np.asarray(p) - np.asarray(y)) ** 2))


def log_loss(y: npt.NDArray[np.float64], p: npt.NDArray[np.float64]) -> float:
    pc = np.clip(np.asarray(p), _EPS, 1.0 - _EPS)
    yv = np.asarray(y)
    return float(-np.mean(yv * np.log(pc) + (1.0 - yv) * np.log(1.0 - pc)))


def pr_auc(y: npt.NDArray[np.float64], p: npt.NDArray[np.float64]) -> float:
    """Average-precision 型の PR-AUC（ties は順位で解消）。"""
    yv = np.asarray(y)
    order = np.argsort(-np.asarray(p), kind="stable")
    ys = yv[order]
    tp = np.cumsum(ys)
    precision = tp / np.arange(1, ys.size + 1)
    recall = tp / max(float(tp[-1]), 1.0)
    delta_r = np.diff(np.concatenate(([0.0], recall)))
    return float(np.sum(precision * delta_r))


def top_area_capture(
    y: npt.NDArray[np.float64],
    p: npt.NDArray[np.float64],
    frac: float = 0.10,
) -> float:
    """予測上位 frac 面積（セル数）が実浸水面積をどれだけ捕捉するか（H3）。

    カットオフ上で予測値が同値（ties）なセルは、**ランダムタイブレーク下の
    期待値**として寄与を按分する。ソート順（=格子行順≈空間順）依存を残す
    旧実装の definition fix で、凍結時に protocol へ注記する。
    """
    yv = np.asarray(y)
    pv = np.asarray(p)
    k = max(1, round(frac * yv.size))
    flooded = float(yv.sum())
    if flooded == 0.0:
        msg = "浸水セルが 0 の集合では capture は未定義"
        raise ValueError(msg)
    t = np.partition(pv, pv.size - k)[pv.size - k]
    strict = pv > t
    wet = float(yv[strict].sum())
    ties = pv == t
    n_tie = int(ties.sum())
    # 同値（境界セル自身を含む n_tie 個）へ残り枠 (k - strict) を按分。
    # 無 ties の通常ケースは n_tie==1・strict==k-1 なので境界セル y をそのまま加算する
    # （v5 の定義修正で n_tie>1 ガードが k 番目を落としていたバグを 2026-09-23 解消）。
    if n_tie:
        wet += float(yv[ties].sum()) * (k - int(strict.sum())) / n_tie
    return wet / flooded
