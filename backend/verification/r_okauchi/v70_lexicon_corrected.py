# %% v70a（事後・amendment B-13 の C2）: 誤りを直した辞書（派生のファイル）
"""凍結した辞書（docs/protocol/v2/lexicon_hydro_v1.yaml）から、B-13 の C2 の (a)〜(h) を機械的に当てた派生の辞書を
書く。凍結したファイルは変えない。規則は topsem.world.b13 の定数（HOMOGRAPHS・RIVERINE_NAMED・QUALIFIERS・SPRINGS・
RAPIDS・ADDITIONS・BLOCK_REMOVALS）で、B-13 の文のとおり。

出力: ``docs/verification/v70_lexicon_corrected.yaml`` と、外した・足した (言語, 形, 位置, 区分, 規則) と国の
ブロックの変更の一覧 ``docs/verification/v70_lexicon_corrected_record.json``。
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.world.b13 import correct_lexicon_doc  # noqa: E402
from topsem.world.lexicon import Lexicon, _Yaml12Loader  # noqa: E402

SRC = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
DST = ROOT / "docs" / "verification" / "v70_lexicon_corrected.yaml"
RECORD = ROOT / "docs" / "verification" / "v70_lexicon_corrected_record.json"


def main() -> None:
    raw = SRC.read_bytes()
    doc = yaml.load(raw, Loader=_Yaml12Loader)  # noqa: S506 — SafeLoader の派生（国コード NO を文字列のまま読む）
    out, rec = correct_lexicon_doc(doc)
    head = ("# 派生の辞書（amendment B-13 の C2）。凍結した lexicon_hydro_v1.yaml"
            f"（sha256 {hashlib.sha256(raw).hexdigest()}）から\n"
            "# B-13 の C2 の (a)〜(h)（v57 の掃除、同形異義、名指された川沿いの項目、相対の修飾語、泉・井戸、早瀬・滝、\n"
            "# 名指された欠け、国のブロック）を backend/verification/r_okauchi/v70_lexicon_corrected.py で機械的に\n"
            "# 当てたもの。手で編集しない。\n")
    DST.write_text(head + yaml.safe_dump(out, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")
    lex0, lex1 = Lexicon(SRC), Lexicon(DST)
    n0 = sum(len(v) for v in lex0.by_lang.values())
    n1 = sum(len(v) for v in lex1.by_lang.values())
    summary = {"rule": "amendment B-13 item C2 (a)-(h)", "src_sha256": lex0.sha256, "dst_sha256": lex1.sha256,
               "entries_src": n0, "entries_dst": n1,
               "dropped_by_rule": dict(Counter(d["rule"] for d in rec["dropped"])),
               "added": len(rec["added"]), "countries": rec["countries"]}
    RECORD.write_text(json.dumps({**summary, **rec}, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
