"""N-read: 言語モデルによる地名の語義判定（preregistration_v2.md §3・H4）。

プロンプトと応答の解釈だけを持つ純関数モジュール（vLLM には依存しない）。推論は
``scripts/nread_infer.py`` が Miyabi-G で行う。プロンプトは凍結対象（文字列の sha256 を記録）。

暗記の統制（H4）: プロンプトは「名前の要素の字義（語源）だけで判定し、その場所の出来事の
知識は使わない」と指示するが、モデルが従う保証はないので、評価は**モデルの学習データ
打切り日より後に起きた洪水**に限る（Gemma-4-31B-it: 2025-01、モデルカード記載）。
"""

from __future__ import annotations

import hashlib
import json
import re

MODEL = "google/gemma-4-31B-it"
# HF の snapshot（コミット）を固定する。Miyabi のオフラインキャッシュには refs/main が無いので、
# 名前だけでは解決できない（2026-10-25 の試走 3700861 で確認）。凍結対象
REVISION = "842da3794eaa0b77d5f08bae87a17459d91ff475"
NREAD_CATEGORIES = (
    "water", "wetland", "riverine", "crossing", "lowland", "highland",
    "other_environmental", "non_environmental", "unknown",
)
SYSTEM = "You are an expert in toponymy (the etymology of place names) and historical linguistics."
_USER = """Place name: "{name}"{native}
Country: {country}

Classify the LITERAL meaning of this place name's elements (its etymology). Use only linguistic
knowledge of the words in the name; do not use any knowledge of events, floods or other facts
about this particular place. Which category best describes what the name refers to?
- water: river, stream, lake, pond, spring, well
- wetland: marsh, swamp, bog, fen, reeds, mud
- riverine: river island, floodplain meadow, riverbank, sandbank, levee
- crossing: ford, bridge, ferry, landing place
- lowland: valley, hollow, low or lower ground
- highland: hill, mountain, height, ridge, rock, upper ground
- other_environmental: forest, field, plant, animal, soil or other nature not related to water
- non_environmental: person, saint, religion, building, administration, event, or other
- unknown: the meaning is not known
Answer with only a JSON object:
{{"category": "<one category>", "confidence": <number from 0 to 1>}}"""
_JSON = re.compile(r"\{[^{}]*\}")


def user_prompt(asciiname: str, native: str, country: str) -> str:
    nat = f' (local spelling: "{native}")' if native and native != asciiname else ""
    return _USER.format(name=asciiname, native=nat, country=country)


def prompt_sha256() -> str:
    return hashlib.sha256((SYSTEM + "\n" + _USER).encode("utf-8")).hexdigest()


def parse(text: str) -> tuple[str, float, bool]:
    """応答 → (category, confidence, parsed_ok)。壊れた応答は unknown / 0。"""
    m = _JSON.search(text or "")
    if not m:
        return "unknown", 0.0, False
    try:
        obj = json.loads(m.group(0))
    except json.JSONDecodeError:
        return "unknown", 0.0, False
    cat = str(obj.get("category", "")).strip().lower()
    if cat not in NREAD_CATEGORIES:
        return "unknown", 0.0, False
    try:
        conf = float(obj.get("confidence", 0.0))
    except (TypeError, ValueError):
        conf = 0.0
    return cat, min(max(conf, 0.0), 1.0), True
