"""route_a: 字面・読み・辞書特徴（elements.yaml 由来）。

名称と（あれば）カタカナ読みから、要素の surfaces / readings /
substitution_of / specificity / hazard_types を引く。
辞書データそのものは**評価イベントより前の版**を使う（時点台帳 — 現時点では
app 側の更新履歴が浅く、pin した SHA がそのまま版になる）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

_FLOOD_TYPES = frozenset({"flood", "inundation"})
_MIN_READING_LEN = 2  # 1 音の読みは誤爆が多すぎる


@dataclass(frozen=True)
class ElementIndex:
    """要素辞書の検索インデックス。"""

    surfaces: list[tuple[str, dict[str, Any]]]
    readings: list[tuple[str, dict[str, Any]]]

    @classmethod
    def from_elements(cls, elements: tuple[dict[str, Any], ...]) -> ElementIndex:
        surf = [(str(s), e) for e in elements for s in (e.get("surfaces") or [])]
        read = [
            (str(r), e)
            for e in elements
            for r in (e.get("readings") or [])
            if len(str(r)) >= _MIN_READING_LEN
        ]
        return cls(surfaces=surf, readings=read)


def name_features(name: str, name_kana: str | None, index: ElementIndex) -> dict[str, float]:
    """1 名称 → route_a 特徴ベクトル（dict 稀疏）。

    重みなしカウントに加え、specificity 加重和（低頻度・高特異の要素を大きく）を
    持つ。洪水特化（flood_*）と全災害（any_*）を分離して、H3 の解釈に使えるようにする。
    """
    hits: list[dict[str, Any]] = []
    seen: set[str] = set()
    for s, e in index.surfaces:
        if s and s in name and e.get("id") not in seen:
            hits.append(e)
            seen.add(str(e.get("id")))
    if name_kana:
        for r, e in index.readings:
            if r in name_kana and e.get("id") not in seen:
                hits.append(e)
                seen.add(str(e.get("id")))
    is_flood = [bool(_FLOOD_TYPES & set(h.get("hazard_types", []) or [])) for h in hits]
    specs = [float(h.get("specificity", 0.0)) for h in hits]
    sub = sum(1 for h in hits if h.get("substitution_of"))
    return {
        "a_n_elements": float(len(hits)),
        "a_n_flood_elements": float(sum(is_flood)),
        "a_spec_sum": float(sum(specs)),
        "a_spec_max": float(max(specs, default=0.0)),
        "a_flood_spec_sum": float(sum(s for s, fl in zip(specs, is_flood, strict=True) if fl)),
        "a_n_substituted": float(sub),
        "a_any_hit": float(bool(hits)),
        "a_flood_hit": float(any(is_flood)),
    }
