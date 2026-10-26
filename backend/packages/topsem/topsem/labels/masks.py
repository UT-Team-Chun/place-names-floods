"""観測マスク契約 — 未観測 ≠ 非浸水（計画 §6.3 / CLAUDE.md）。

浸水ラベルは各セル 3 状態。解析側は必ず used で絞ってから y を読むこと。

- FLOODED: 観測あり・浸水（航空写真・被災地図・GCP 等での判読）
- DRY:     観測あり・非浸水（判読可能領域で水没を確認しなかった）
- UNOBSERVED: 雲・樹冠・未調査（loss / metric から除外）
"""

from __future__ import annotations

from enum import IntEnum


class ObsState(IntEnum):
    FLOODED = 1
    DRY = 0
    UNOBSERVED = -1


def used_mask(states: list[ObsState]) -> list[bool]:
    return [s is not ObsState.UNOBSERVED for s in states]
