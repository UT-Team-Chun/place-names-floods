"""amendment B-12（初見の査読 第 8 回、事後）の実装。

- C1: plasmode の形 BK。形 B（水系側の名が最寄りで、その集落から 500 m 以内の
  セル）の距離を 1,000 m にした形。注入の対象の選び方を形ごとに一つの関数に
  まとめる（形 A・B・CA・CB・AH・AC は従来と同じ対象）
"""

from __future__ import annotations

import numpy as np

PLASMODE_LOCAL_M = 500.0  # 形 B・CB（B-6）
PLASMODE_1KM_M = 1000.0  # 形 BK（B-12 の C1）


def plasmode_target(form: str, t: np.ndarray, d1: np.ndarray) -> np.ndarray:
    """注入の対象のセル。

    t は水系側（形 C* は読み違いを入れた潜在の印）の最寄りの名、d1 は最寄りの
    集落までの距離（m）。形の名の最後が B なら 500 m 以内、BK なら 1,000 m 以内に
    限る。それ以外（A・CA・AH・AC）は t のまま。
    """
    t = np.asarray(t, dtype=bool)
    if form == "BK":
        return t & (np.asarray(d1) <= PLASMODE_1KM_M)
    if form.endswith("B"):
        return t & (np.asarray(d1) <= PLASMODE_LOCAL_M)
    return t
