"""解析母集団（一般地名を含む）の構築（計画 §4.3）。

候補スクリーニングは警鐘地名「だけ」だと、地名を全く持たない一般セルが
母集団から消えて選択バイアスになる。そのため Geolonia 住所テーブル相当の
**全国地名を母集団**とし、キュレーション済みコーパスと辞書要素をそこに
「印を付けた」形で作る（app 側の candidates screen が落とした負例の復元）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd

from topsem.corpus.snapshot import CorpusSnapshot

_FLOOD_TYPES = frozenset({"flood", "inundation"})


@dataclass(frozen=True)
class NameIndex:
    """母集団行の文字列に付ける印の辞書（部分一致）。"""

    curated: dict[str, str]  # surface -> 地名レコード name
    hazard: dict[str, str]  # surface -> element id
    flood: frozenset[str]  # hazard のうち flood / inundation を含む id

    @classmethod
    def from_snapshot(cls, snapshot: CorpusSnapshot) -> NameIndex:
        curated: dict[str, str] = {}
        for t in snapshot.toponyms:
            name = str(t.get("name", ""))
            surfaces = [name, *map(str, t.get("variants", []) or [])]
            for s in surfaces:
                if s:
                    curated[s] = name
        hazard: dict[str, str] = {}
        flood_ids: set[str] = set()
        for e in snapshot.elements:
            eid = str(e.get("id", "?"))
            if _FLOOD_TYPES & set(e.get("hazard_types", []) or []):
                flood_ids.add(eid)
            for s in map(str, e.get("surfaces", []) or []):
                if s:
                    hazard[s] = eid
        return NameIndex(curated=curated, hazard=hazard, flood=frozenset(flood_ids))


def tag_rows(
    df: pd.DataFrame,
    name_col: str,
    index: NameIndex,
) -> pd.DataFrame:
    """name_col を持つ行に is_curated / hazard_elements / is_flood_word を付す。

    母集団の正式定義（同一地名の重複行、大字/小字の扱い）は
    scripts/build_population.py 側の docs に従う。ここは素朴な部分一致版。
    """
    haz_items = list(index.hazard.items())
    curated_surfaces = list(index.curated)
    elem_hits: list[str] = []
    curated_hits: list[bool] = []
    flood_hits: list[bool] = []
    for value in df[name_col]:
        s = "" if pd.isna(value) else str(value)
        matched_ids = sorted({eid for surf, eid in haz_items if surf in s})
        elem_hits.append(",".join(matched_ids))
        curated_hits.append(any(k in s for k in curated_surfaces))
        flood_hits.append(any(eid in index.flood for eid in matched_ids))
    out = df.copy()
    out["is_curated"] = curated_hits
    out["hazard_elements"] = elem_hits
    out["is_flood_word"] = flood_hits
    return out


def population_stats(tagged: pd.DataFrame) -> dict[str, Any]:
    n = len(tagged)
    return {
        "n_rows": n,
        "frac_any_hazard_word": float((tagged["hazard_elements"] != "").mean()) if n else 0.0,
        "frac_flood_word": float(tagged["is_flood_word"].mean()) if n else 0.0,
        "frac_curated": float(tagged["is_curated"].mean()) if n else 0.0,
    }


def infer_name_column(df: pd.DataFrame) -> str:
    """Geolonia 住所 CSV 等から地名らしき列を推定する（無ければ分かりやすい失敗）。"""
    for cand in (
        "大字町丁目名",  # Geolonia 住所（2026 時点の列名）
        "大字・丁目名",
        "name",
        "oaza",
        "chome",
        "大字",
        "町名",
    ):
        if cand in df.columns:
            return cand
    for col in df.columns:
        if "name" in str(col).lower():
            return str(col)
    msg = "地名列を推定できない。--name-col で明示してください。列: " + ", ".join(
        map(str, df.columns)
    )
    raise KeyError(msg)
