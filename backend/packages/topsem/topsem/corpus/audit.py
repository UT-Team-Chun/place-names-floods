"""Gate 1 監査: コーパスの成立性を数える（モデルより先にこれを回す）。

見るべきは「件数」ではなく、河川洪水ラベル付けに使えるか:

- 洪水（flood / inundation）系の hazard_types を持つ件数
- 記録内の証拠の最高 level の分布（level 3 がどれだけ残るか）
- 反証（stance: disputes）を持つ件数
- location（lat/lon）保持率と precision の内訳
"""

from __future__ import annotations

import json
from collections import Counter
from typing import Any

from topsem.corpus.snapshot import CorpusSnapshot

_FLOOD_TYPES = frozenset({"flood", "inundation"})


def _max_level(toponym: dict[str, Any]) -> int:
    levels = [int(e.get("level", 0)) for e in toponym.get("evidence", [])]
    return max(levels, default=0)


def audit_corpus(snapshot: CorpusSnapshot) -> dict[str, Any]:
    """監査サマリ（JSON 化可能な dict）を返す。"""
    topos = list(snapshot.toponyms)
    hazards: Counter[str] = Counter()
    for t in topos:
        for h in t.get("hazard_types", []):
            hazards[str(h)] += 1
    flood_n = sum(1 for t in topos if _FLOOD_TYPES & set(t.get("hazard_types", [])))
    levels = Counter(_max_level(t) for t in topos)
    disputed = sum(
        1 for t in topos if any(e.get("stance") == "disputes" for e in t.get("evidence", []))
    )
    with_loc = sum(1 for t in topos if "location" in t)
    precision: Counter[str] = Counter(
        str(t["location"].get("precision", "unknown")) for t in topos if "location" in t
    )
    status: Counter[str] = Counter(str(t.get("status", "unknown")) for t in topos)

    elems = list(snapshot.elements)
    elem_hazards: Counter[str] = Counter()
    for e in elems:
        for h in e.get("hazard_types", []):
            elem_hazards[str(h)] += 1
    flood_elems = [
        str(e.get("id", "?")) for e in elems if _FLOOD_TYPES & set(e.get("hazard_types", []))
    ]
    substituted = sum(1 for e in elems if e.get("substitution_of"))

    return {
        "n_toponyms": len(topos),
        "toponyms_by_status": dict(status),
        "toponyms_by_hazard_type": dict(hazards),
        "n_toponyms_flood_any": flood_n,
        "toponym_evidence_max_level_hist": {str(k): v for k, v in sorted(levels.items())},
        "n_toponyms_with_dispute": disputed,
        "toponyms_with_location": with_loc,
        "location_precision": dict(precision),
        "n_elements": len(elems),
        "elements_by_hazard_type": dict(elem_hazards),
        "element_ids_flood_any": sorted(flood_elems),
        "n_elements_substitution": substituted,
    }


def audit_built_geojson(snapshot: CorpusSnapshot) -> dict[str, Any] | None:
    """アプリ組み立て済み geojson があれば件数を読む（無ければ None）。

    組み立て成果物なので git 管理外・要 `make build-data`。監査の主対象は
    curated YAML の方であり、これは「ビルド版との一致感」の確認用。
    """
    path = snapshot.app_root / "frontend" / "public" / "data" / "toponyms.geojson"
    if not path.is_file():
        return None
    gj = json.loads(path.read_text(encoding="utf-8"))
    feats = (gj or {}).get("features", [])
    # 注: ビルド版の properties は camelCase（hazardTypes / evidenceLevel）
    hazards: Counter[str] = Counter()
    levels: Counter[int] = Counter()
    flood_n = 0
    disputed = 0
    no_hazard = 0
    for f in feats:
        p = f.get("properties") or {}
        hts = {str(h) for h in p.get("hazardTypes") or []}
        for h in hts:
            hazards[h] += 1
        if not hts:
            no_hazard += 1
        if hts & _FLOOD_TYPES:
            flood_n += 1
        levels[int(p.get("evidenceLevel", 0))] += 1
        disputed += int(bool(p.get("disputed")))
    return {
        "n_features": len(feats),
        "by_hazard_type": dict(hazards),
        "n_flood_any": flood_n,
        "n_without_hazard_type": no_hazard,
        "evidence_level_hist": {str(k): v for k, v in sorted(levels.items())},
        "n_disputed": disputed,
    }
