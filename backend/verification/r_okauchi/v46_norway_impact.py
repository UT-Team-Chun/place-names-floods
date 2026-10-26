# %% v46（amendment B-4 の記述の再計算）: ノルウェーの辞書の欠陥が触れたセルの数
"""B-4 に書いた影響の大きさ（主分割の標本で最寄りの集落がノルウェーのセル、直した読み方なら辞書に一致する
セル）を、論文の数値の出どころとして成果物から計算し直す。浸水ラベルは使わない。

出力: stdout と docs/verification/v46_norway_impact.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR
from topsem.world.lexicon import CATEGORIES, Lexicon

LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
OUT = ROOT / "docs" / "verification" / "v46_norway_impact.json"


def main() -> None:
    t = pd.read_parquet(DATA_DIR / "verification" / "v37_table.parquet", columns=["aoi", "P_gid1"])
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    cc = t["P_gid1"].map(gn["cc"])
    no = (cc == "NO").to_numpy()
    lex = Lexicon(LEX)  # B-4 の読み方
    names = t.loc[no, "P_gid1"].map(gn["asciiname"]).astype(str)
    hit = [any(lex.match(n, "NO")[c] for c in CATEGORIES) for n in names]
    res = {"n_cells": len(t), "n_cells_norway": int(no.sum()),
           "n_aoi_norway": int(t.loc[no, "aoi"].nunique()),
           "n_cells_norway_matched": int(sum(hit)),
           "share_matched_of_norway": sum(hit) / max(1, int(no.sum())),
           "share_matched_of_all": sum(hit) / len(t),
           "lexicon_covers_norway_after_fix": lex.coverage("NO"),
           # 登録した層「辞書の被覆」（B-2）: 直した読みでは、最寄りの集落の国のどれも辞書の言語を持つ
           "n_countries_nearest": int(cc.nunique()),
           "uncovered_countries_after_fix": sorted(str(c) for c in cc.dropna().unique() if not lex.coverage(str(c)))}
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
