"""多言語の水文・地形形態素辞書（N-lex）の照合。

辞書本体は ``docs/protocol/v2/lexicon_hydro_v1.yaml``（凍結対象・sha256 を run に記録）。
辞書は**地名・ラベル・イベントを一切見ずに**作成した一般名詞の形態素だけを持つ
（固有名・河川名は入れない）ので、場所ごとの洪水履歴を暗記できない。

照合は GeoNames の ``asciiname`` に対して行う:
小文字化 → 英字以外で分割 → 地名の国（``cc``）に対応する言語の項目だけを適用。
位置規則:

- ``word``   : トークン全体が一致
- ``prefix`` : トークンが form で始まり、残りが 2 文字以上
- ``suffix`` : トークンが form で終わり、前が 2 文字以上
- ``infix``  : トークン内のどこか（form が 4 文字以上のときだけ許す）

1 項目に form・pos のリストを書ける
（例: ``{form: [berg, berge], pos: [word, suffix]}``）。

カテゴリ（6 種）: water / wet / riverine / crossing / low / high。
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml

CATEGORIES = ("water", "wet", "riverine", "crossing", "low", "high")
WET_SIDE = ("water", "wet", "riverine", "crossing", "low")
_SPLIT = re.compile(r"[^a-z]+")
MIN_REST = 2
MIN_INFIX = 4


class _Yaml12Loader(yaml.SafeLoader):
    """YAML 1.2 と同じく true / false だけを真偽値にする（yes / no / on / off は文字列のまま）。

    PyYAML の既定（YAML 1.1）では引用符なしの国コード ``NO``（ノルウェー）が false になり、
    ノルウェーの地名に辞書が当たらなかった（amendment B-4）。凍結した辞書のファイルは変えず、
    読み方を直す。
    """


_Yaml12Loader.yaml_implicit_resolvers = {
    k: [(tag, rx) for tag, rx in v if tag != "tag:yaml.org,2002:bool"]
    for k, v in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
_Yaml12Loader.add_implicit_resolver(
    "tag:yaml.org,2002:bool", re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"), list("tTfF"))


@dataclass(frozen=True)
class Entry:
    form: str
    pos: str
    cat: str
    lang: str


class Lexicon:
    """YAML 辞書 → 国別の照合器。"""

    def __init__(self, path: Path) -> None:
        raw = path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        doc = yaml.load(raw, Loader=_Yaml12Loader)  # noqa: S506 — SafeLoader の派生（真偽値の解釈だけ変える）
        self.version = str(doc.get("version", ""))
        self.by_lang: dict[str, tuple[Entry, ...]] = {}
        for lang, spec in (doc.get("languages") or {}).items():
            ents = []
            for e in spec.get("entries") or []:
                cat = str(e["cat"])
                poss = e.get("pos", "word")
                poss = [poss] if isinstance(poss, str) else list(poss)
                forms = e["form"] if isinstance(e["form"], list) else [e["form"]]
                if cat not in CATEGORIES:
                    msg = f"unknown category {cat!r} in {lang}"
                    raise ValueError(msg)
                for pos in poss:
                    if pos not in ("word", "prefix", "suffix", "infix"):
                        msg = f"unknown pos {pos!r} in {lang}"
                        raise ValueError(msg)
                    for f in forms:
                        form = str(f).lower()
                        if pos == "infix" and len(form) < MIN_INFIX:
                            msg = f"infix too short {form!r} in {lang}"
                            raise ValueError(msg)
                        ents.append(Entry(form=form, pos=pos, cat=cat, lang=lang))
            self.by_lang[lang] = tuple(ents)
        self.country_langs: dict[str, tuple[str, ...]] = {
            str(cc).upper(): tuple(v) for cc, v in (doc.get("countries") or {}).items()
        }

    def entries_for(self, cc: str) -> tuple[Entry, ...]:
        langs = self.country_langs.get(str(cc).upper(), ())
        return tuple(e for lang in langs for e in self.by_lang.get(lang, ()))

    def match(self, asciiname: str, cc: str) -> dict[str, int]:
        return dict(zip(CATEGORIES, _match_cached(self, str(asciiname), str(cc)), strict=True))

    def coverage(self, cc: str) -> bool:
        return bool(self.entries_for(cc))


def _hit(tok: str, e: Entry) -> bool:
    if e.pos == "word":
        return tok == e.form
    if e.pos == "prefix":
        return tok.startswith(e.form) and len(tok) - len(e.form) >= MIN_REST
    if e.pos == "suffix":
        return tok.endswith(e.form) and len(tok) - len(e.form) >= MIN_REST
    return e.form in tok


@lru_cache(maxsize=1_000_000)
def _match_cached(lex: Lexicon, asciiname: str, cc: str) -> tuple[int, ...]:
    toks = [t for t in _SPLIT.split(asciiname.lower()) if t]
    counts = dict.fromkeys(CATEGORIES, 0)
    for e in lex.entries_for(cc):
        if any(_hit(t, e) for t in toks):
            counts[e.cat] += 1
    return tuple(min(counts[c], 1) for c in CATEGORIES)


def cell_lexicon_features(
    gids: np.ndarray,
    dists: np.ndarray,
    gid_lookup: dict[int, tuple[str, str]],
    lex: Lexicon,
) -> dict[str, np.ndarray]:
    """最近傍 K 名（gids/dists: [n, K]）→ セルの辞書特徴（最近傍 1 名と距離重み付き K 名）。"""
    n, k = gids.shape
    # 一意な gid ごとに 1 回だけ照合して配る（数百万セル × K の Python ループを避ける）
    uniq, inv = np.unique(np.asarray(gids, dtype=np.int64).ravel(), return_inverse=True)
    cat_u = np.zeros((len(uniq), len(CATEGORIES)), dtype=np.float32)
    cov_u = np.zeros(len(uniq), dtype=np.float32)
    for i, g in enumerate(uniq.tolist()):
        if g < 0:
            continue
        name, cc = gid_lookup.get(g, ("", ""))
        cat_u[i] = _match_cached(lex, name, cc)
        cov_u[i] = float(lex.coverage(cc))
    cat = cat_u[inv].reshape(n, k, len(CATEGORIES))
    covered = cov_u[inv].reshape(n, k)
    w = 1.0 / (1.0 + np.where(np.isfinite(dists), dists, 1e9) / 1000.0)
    out: dict[str, np.ndarray] = {}
    for c_i, c in enumerate(CATEGORIES):
        out[f"lex_{c}_n1"] = cat[:, 0, c_i]
        out[f"lex_{c}_k"] = (cat[:, :, c_i] * w).sum(axis=1).astype(np.float32)
    out["lex_wetside_n1"] = np.clip(cat[:, 0, [CATEGORIES.index(c) for c in WET_SIDE]].sum(axis=1),
                                    0, 1).astype(np.float32)
    out["lex_covered_n1"] = covered[:, 0]
    return out
