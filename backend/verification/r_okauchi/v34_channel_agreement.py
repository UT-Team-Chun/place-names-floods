# %% v34: 名前チャネルの一致度 — 多言語辞書（N-lex）と言語モデル（N-read）は同じ名前を同じ意味に読むか
"""ラベルを使わない照合（凍結前に許される）。採用 AOI に割り当てられた一意な P 名 186,638 件について、
N-lex の 6 区分（water / wet / riverine / crossing / low / high）と N-read の 9 類を突き合わせる。

- 辞書が区分を付けた名前のうち、N-read が対応する類（water↔water、wet↔wetland、riverine↔riverine、
  crossing↔crossing、low↔lowland、high↔highland）を返した割合（辞書から見た一致率）
- N-read が水辺系（water / wetland / riverine / crossing / lowland）と読んだ名前のうち、辞書も水辺系の
  区分を付けた割合（N-read から見た辞書の被覆）
- 大陸別の割合（国 → 大陸は GeoNames countryInfo）

入力: data/global/nread/names.parquet・nread.parquet、docs/protocol/v2/lexicon_hydro_v1.yaml
出力: stdout の表と docs/verification/v34_channel_agreement.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

import numpy as np
import pandas as pd

from topsem.world.lexicon import CATEGORIES, Lexicon, _match_cached

PAIRS = {"water": "water", "wet": "wetland", "riverine": "riverine", "crossing": "crossing",
         "low": "lowland", "high": "highland"}
WET_LEX = ("water", "wet", "riverine", "crossing", "low")
WET_NR = ("water", "wetland", "riverine", "crossing", "lowland")


# %%
def main() -> None:
    g = ROOT / "data" / "global"
    names = pd.read_parquet(g / "nread" / "names.parquet")
    nr = pd.read_parquet(g / "nread" / "nread.parquet")[["gid", "category", "parsed_ok"]]
    d = names.merge(nr, on="gid", how="left")
    lex = Lexicon(ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    m = np.array([_match_cached(lex, a, c) for a, c in
                  zip(d["asciiname"].astype(str), d["cc"].astype(str), strict=True)])
    cont = {}
    for line in (g / "geonames" / "countryInfo.txt").read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            cont[f[0]] = f[8]
    d["continent"] = d["cc"].map(cont).fillna("??")
    out: dict = {"n_names": len(d), "lexicon_sha256": lex.sha256,
                 "nread_parsed_ok": float(d["parsed_ok"].mean())}
    agree = {}
    for c_lex, c_nr in PAIRS.items():
        hit = m[:, CATEGORIES.index(c_lex)] > 0
        agree[c_lex] = {"n_lexicon": int(hit.sum()),
                        "nread_same_class": float((d.loc[hit, "category"] == c_nr).mean())
                        if hit.any() else None,
                        "nread_class_counts": dict(Counter(d.loc[hit, "category"]).most_common(4))}
    out["lexicon_to_nread"] = agree
    lex_wet = m[:, [CATEGORIES.index(c) for c in WET_LEX]].max(axis=1) > 0
    nr_wet = d["category"].isin(WET_NR).to_numpy()
    out["wetside"] = {"lexicon_share": float(lex_wet.mean()), "nread_share": float(nr_wet.mean()),
                      "nread_given_lexicon": float(nr_wet[lex_wet].mean()),
                      "lexicon_given_nread": float(lex_wet[nr_wet].mean())}
    by = {}
    for k, sub in d.groupby("continent"):
        idx = sub.index.to_numpy()
        by[k] = {"n": len(sub), "lexicon_wetside": float(lex_wet[idx].mean()),
                 "nread_wetside": float(nr_wet[idx].mean())}
    out["by_continent"] = by
    (ROOT / "docs" / "verification" / "v34_channel_agreement.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"names {len(d):,}  lexicon wetside {out['wetside']['lexicon_share']:.3f}  "
          f"N-read wetside {out['wetside']['nread_share']:.3f}")
    print(f"N-read wetside | lexicon wetside = {out['wetside']['nread_given_lexicon']:.3f}; "
          f"lexicon wetside | N-read wetside = {out['wetside']['lexicon_given_nread']:.3f}")
    for c, v in agree.items():
        same = v["nread_same_class"]
        print(f"  lexicon {c:9s} n={v['n_lexicon']:6d}  N-read same class "
              f"{'' if same is None else f'{same:.3f}'}  {v['nread_class_counts']}")
    for k, v in sorted(by.items(), key=lambda kv: -kv[1]["n"]):
        print(f"  {k}: n={v['n']:6d} lexicon {v['lexicon_wetside']:.3f} N-read {v['nread_wetside']:.3f}")


if __name__ == "__main__":
    main()
