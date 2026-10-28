# %% v57a（事後・amendment B-8 の B6(d)・C5）: 掃除した辞書（派生のファイル）
"""凍結した辞書（docs/protocol/v2/lexicon_hydro_v1.yaml）から、2 文字以下の語頭・語末の項目と、英語の接尾辞
-burn・-bourne・-borne・-pool・-well・-wells を機械的に除いた派生の辞書を書く。凍結したファイルは変えない。

出力: ``docs/verification/v57_lexicon_clean.yaml`` と、除いた (言語, 形, 位置) の一覧
``docs/verification/v57_lexicon_clean_dropped.json``。
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import yaml

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.world.b8 import clean_lexicon_doc  # noqa: E402
from topsem.world.lexicon import Lexicon, _Yaml12Loader  # noqa: E402

SRC = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
DST = ROOT / "docs" / "verification" / "v57_lexicon_clean.yaml"
DROPPED = ROOT / "docs" / "verification" / "v57_lexicon_clean_dropped.json"


def main() -> None:
    raw = SRC.read_bytes()
    doc = yaml.load(raw, Loader=_Yaml12Loader)  # noqa: S506 — SafeLoader の派生（国コード NO を文字列のまま読む）
    out, dropped = clean_lexicon_doc(doc)
    head = ("# 派生の辞書（amendment B-8 の B6(d)・C5）。凍結した lexicon_hydro_v1.yaml"
            f"（sha256 {hashlib.sha256(raw).hexdigest()}）から\n"
            "# 2 文字以下の語頭・語末の項目と、英語の接尾辞 -burn・-bourne・-borne・-pool・-well・-wells を\n"
            "# backend/verification/r_okauchi/v57_lexicon_clean.py で機械的に除いたもの。手で編集しない。\n")
    DST.write_text(head + yaml.safe_dump(out, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")
    lex0, lex1 = Lexicon(SRC), Lexicon(DST)
    n0 = sum(len(v) for v in lex0.by_lang.values())
    n1 = sum(len(v) for v in lex1.by_lang.values())
    assert lex1.country_langs == lex0.country_langs, "country map changed"
    rec = {"rule": "amendment B-8 items B6(d) and C5", "src_sha256": lex0.sha256, "dst_sha256": lex1.sha256,
           "entries_src": n0, "entries_dst": n1, "dropped": [list(x) for x in dropped]}
    DROPPED.write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in rec.items() if k != "dropped"}), len(dropped), "dropped")
    print(sorted({(a, b, c) for a, b, c in dropped}))


if __name__ == "__main__":
    main()
