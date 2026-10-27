# %% v51（事後・amendment B-6 の 5）: 名前の読みの文書化と、厳密な読み
"""確証判定は変えない。B-6 の 5 で計算の前に固定した量を、辞書・解析した表・v44 の判定から出す。

- 辞書の表: 言語の符号ごとに、含む言語（YAML の name）、区分ごとの形の数、その符号を当てる国、照合の位置。
  衝突のため除いた形（YAML の冒頭の注記）
- 一致の統計: 解析したセルの最寄りの名前が、いずれかの区分・水系側に一致した割合（言語の符号別・大陸別）。
  水系側で一致の多い形の上位 20（形・位置・セル数・集落数）。読みは B-4 の読み方（直したコード）
- 一致の位置別の精度: v44 の判定（由来あり）のうち辞書が水系側と読んだ名について、水系側の形が一致した位置
  （語・接頭・接尾・語中）別の精度（重み無し、Wilson の区間）。位置が複数の名はそれぞれに数える
- 厳密な定義: 辞書の厳密な読み（水は泉・井戸・水場の形を除く・湿地・河川）と、記述の厳密な由来（水・湿地・河川、
  ただし根拠の引用が泉・井戸の語だけを含む水は除く）。厳密な読みの精度・再現率（v44 と同じ層の重み、層ごとの
  bootstrap 2,000 回）、記述の区分別の精度・再現率（重み無しの件数つき）
- 名前の例: 辞書の水系側の区分ごとに、番号の若い順に記述の由来が区分と一致した名を 2 つ、一致しなかった名を 1 つ。
  図 3 の例の AOI の集落名と読み
- 人による照合の準備: v44 の判定から層・言語で層別した 100 件の無作為標本（seed 20261028）と手順書

語の照合の操作的な定義（この版で決めた）: 引用を小文字にして英字（アクセント付きを含む）の連なりに分け、表の語と
一致するか、語に s / es を足した形と一致すれば「現れる」とする。

出力: stdout と docs/verification/v51_readings.json、表は docs/verification/v51_*.csv、手順書は
docs/verification/v51_audit_guide.md
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b6 import SPRING_FORMS, STRICT_CATEGORIES, strict_lexicon
from topsem.world.dataset import continent_map
from topsem.world.lexicon import (
    _SPLIT,
    CATEGORIES,
    WET_SIDE,
    Lexicon,
    _hit,
    _Yaml12Loader,
)

DOCS = ROOT / "docs" / "verification"
OUT = DOCS / "v51_readings.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
EXAMPLE_AOI = "EMSR511_AOI01"  # 図 3（集落を P 類に限って選び直した例、make_figures の pick_example）
SEED = 20261028
N_AUDIT = 100
N_BOOT = 2000
TOP = 20
SPRING_WORDS = (*sorted({f for v in SPRING_FORMS.values() for f in v}),
                "spring", "well", "fountain", "source", "quelle", "brunnen", "fontaine", "puits", "pozzo", "fonte",
                "fuente", "manantial", "studnia", "studna", "studenec", "izvir", "izvor", "kút", "forrás", "bunar",
                "ain", "bir")
RIVER_WORDS = ("river", "stream", "brook", "creek", "lake", "pond", "marsh", "fluss", "bach", "see", "rivière",
               "ruisseau", "lac", "étang", "rzeka", "potok", "jezioro", "reka", "jezero", "folyó", "patak", "tó")
_WORD = re.compile(r"[^\W\d_]+", re.UNICODE)


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v44 = _load("v44_wikipedia_etymology")


def matched_entries(lex: Lexicon, name: str, cc: str) -> list:
    toks = [t for t in _SPLIT.split(str(name).lower()) if t]
    return [e for e in lex.entries_for(cc) if any(_hit(t, e) for t in toks)]


# %% 1. 辞書の表
def lexicon_table(lex: Lexicon) -> tuple[pd.DataFrame, list[str]]:
    doc = yaml.load(LEX.read_bytes(), Loader=_Yaml12Loader)
    rows = []
    for code, spec in doc["languages"].items():
        ents = lex.by_lang[code]
        rows.append({"code": code, "languages": spec.get("name", ""),
                     "countries": " ".join(sorted(cc for cc, ls in lex.country_langs.items() if code in ls)),
                     **{f"n_{c}": len({e.form for e in ents if e.cat == c}) for c in CATEGORIES},
                     "positions": " ".join(p for p in ("word", "prefix", "suffix", "infix")
                                           if any(e.pos == p for e in ents)),
                     "n_spring_forms_strict": len(SPRING_FORMS.get(code, ()))})
    head = []
    for line in LEX.read_text(encoding="utf-8").splitlines():
        if not line.startswith("#"):
            break
        head.append(line[1:].strip())
    i = next(k for k, s in enumerate(head) if s.startswith("方針"))
    return pd.DataFrame(rows), head[i:i + 2]


# %% 2. 一致の統計（解析したセル）
def match_stats(lex: Lexicon) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    t = pd.read_parquet(V37T, columns=["aoi", "P_gid1"])
    gn = pd.read_parquet(GN, columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    gid_cells = t["P_gid1"].value_counts()
    rows, ent_rows = [], []
    for g, n_cells in gid_cells.items():
        if g < 0 or g not in gn.index:
            rows.append({"gid": int(g), "n_cells": int(n_cells), "cc": "", "covered": False})
            continue
        name, cc = str(gn.at[g, "asciiname"]), str(gn.at[g, "cc"])
        ents = matched_entries(lex, name, cc)
        rows.append({"gid": int(g), "n_cells": int(n_cells), "cc": cc, "covered": lex.coverage(cc),
                     "any": bool(ents), "wetside": any(e.cat in WET_SIDE for e in ents),
                     "codes": tuple(lex.country_langs.get(cc.upper(), ())),
                     "any_by_code": {code: any(e.lang == code for e in ents)
                                     for code in lex.country_langs.get(cc.upper(), ())},
                     "wet_by_code": {code: any(e.lang == code and e.cat in WET_SIDE for e in ents)
                                     for code in lex.country_langs.get(cc.upper(), ())}})
        for e in {(e.lang, e.form, e.pos, e.cat) for e in ents if e.cat in WET_SIDE}:
            ent_rows.append({"gid": int(g), "n_cells": int(n_cells), "lang": e[0], "form": e[1], "pos": e[2],
                             "cat": e[3]})
    s = pd.DataFrame(rows)
    s["continent"] = s["cc"].map(lambda c: cont.get(c, "") if c else "")
    w = s["n_cells"].to_numpy(dtype=float)
    tot = w.sum()
    out = {"n_cells": int(tot), "n_settlements": len(s),
           "share_cells_covered": float(w[s["covered"].to_numpy(dtype=bool)].sum() / tot),
           "share_cells_any": float(w[s["any"].fillna(value=False).to_numpy(dtype=bool)].sum() / tot),
           "share_cells_wetside": float(w[s["wetside"].fillna(value=False).to_numpy(dtype=bool)].sum() / tot)}
    by_code: dict = defaultdict(lambda: {"n_cells": 0, "any": 0, "wet": 0, "n_settlements": 0})
    for r in s.itertuples():
        if not r.covered:
            continue
        for code in r.codes:
            b = by_code[code]
            b["n_cells"] += r.n_cells
            b["n_settlements"] += 1
            b["any"] += r.n_cells * r.any_by_code[code]
            b["wet"] += r.n_cells * r.wet_by_code[code]
    out["by_code"] = {k: {"n_cells": v["n_cells"], "n_settlements": v["n_settlements"],
                          "share_any": v["any"] / v["n_cells"], "share_wetside": v["wet"] / v["n_cells"]}
                      for k, v in sorted(by_code.items(), key=lambda kv: -kv[1]["n_cells"])}
    out["by_continent"] = {}
    for k, d in s.groupby("continent"):
        wd = d["n_cells"].to_numpy(dtype=float)
        out["by_continent"][k or "unknown"] = {
            "n_cells": int(wd.sum()), "share_covered": float(wd[d["covered"].to_numpy(dtype=bool)].sum() / wd.sum()),
            "share_any": float(wd[d["any"].fillna(value=False).to_numpy(dtype=bool)].sum() / wd.sum()),
            "share_wetside": float(wd[d["wetside"].fillna(value=False).to_numpy(dtype=bool)].sum() / wd.sum())}
    e = pd.DataFrame(ent_rows)
    top = (e.groupby(["lang", "form", "pos", "cat"]).agg(n_cells=("n_cells", "sum"), n_settlements=("gid", "nunique"))
           .reset_index().sort_values(["n_cells", "n_settlements"], ascending=False).head(TOP))
    # 確かめ: 一致の読みが、直したコードの実行（world_v2_lexfix）の水系側の印と同じ
    files = sorted((RUNS_DIR / "world_v2_lexfix" / "predictions" / "primary").glob("fold[0-9].npz"),
                   key=lambda f: int(f.stem[4:]))
    flag = np.concatenate([np.load(f)["s_wetside"] for f in files]) > 0
    mine = t["P_gid1"].map(s.set_index("gid")["wetside"].fillna(value=False)).to_numpy(dtype=bool)
    out["reproduces_lexfix_wetside"] = bool(len(flag) == len(mine) and (flag == mine).all())
    return out, top, s


# %% 3〜5. v44 の判定
def coded() -> tuple[pd.DataFrame, dict]:
    samp = json.loads(v44.SAMPLE.read_text(encoding="utf-8"))
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    codes = pd.read_csv(v44.CODES)
    codes["gid"] = codes["id"].map(key)
    d = codes.join(meta, on="gid")
    d["cats_set"] = d["cats"].fillna("").map(lambda s: frozenset(c for c in s.split(";") if c))
    d["truth"] = d["cats_set"].map(lambda s: bool(s & set(WET_SIDE)))
    return d, samp["stratum_sizes"]


def _words(text: str) -> set[str]:
    return set(_WORD.findall(str(text).lower()))


def _has(words: set[str], vocab: tuple[str, ...]) -> bool:
    return any(w in words or f"{w}s" in words or f"{w}es" in words for w in vocab)


def spring_only(quote: str) -> bool:
    w = _words(quote)
    return _has(w, SPRING_WORDS) and not _has(w, RIVER_WORDS)


def by_position(d: pd.DataFrame, lex: Lexicon) -> dict:
    c = d[(d["status"] == "coded") & d["lex_wetside"]].copy()
    c["pos_set"] = [frozenset(e.pos for e in matched_entries(lex, r.name, r.cc) if e.cat in WET_SIDE)
                    for r in c.itertuples()]
    out = {}
    for p in ("word", "prefix", "suffix", "infix"):
        m = c["pos_set"].map(lambda s, p=p: p in s)
        only = c["pos_set"].map(lambda s, p=p: s == frozenset({p}))
        out[p] = {"any": v44._rate(int(c.loc[m, "truth"].sum()), int(m.sum())),
                  "only_this_position": v44._rate(int(c.loc[only, "truth"].sum()), int(only.sum()))}
    out["n_names"] = len(c)
    out["n_multiple_positions"] = int((c["pos_set"].map(len) > 1).sum())
    return out


def strict(d: pd.DataFrame, sizes: dict, lex: Lexicon) -> dict:
    ls = strict_lexicon(lex)
    c = d[d["status"] == "coded"].copy()
    c["strict_flag"] = [any(ls.match(r.name, r.cc)[k] for k in STRICT_CATEGORIES) for r in c.itertuples()]
    c["spring_only"] = [("water" in r.cats_set) and spring_only(r.quote) for r in c.itertuples()]
    c["strict_truth"] = [bool((r.cats_set - ({"water"} if r.spring_only else set())) & set(STRICT_CATEGORIES))
                         for r in c.itertuples()]
    e = c.assign(truth=c["strict_truth"])
    point = v44.weighted(e, sizes, "strict_flag")
    rng = np.random.default_rng(SEED)
    parts = {s: e[e["stratum"] == s] for s in sizes}
    boots = [v44.weighted(pd.concat([p.iloc[rng.integers(0, len(p), len(p))] for p in parts.values()]),
                          sizes, "strict_flag") for _ in range(N_BOOT)]
    res = {k: {"est": point[k], "ci95": [float(np.nanpercentile([b[k] for b in boots], 2.5)),
                                         float(np.nanpercentile([b[k] for b in boots], 97.5))]}
           for k in ("precision", "recall")}
    res |= {"n_coded": len(c), "n_strict_flag": int(c["strict_flag"].sum()),
            "n_strict_truth": int(c["strict_truth"].sum()), "n_water_spring_only": int(c["spring_only"].sum()),
            "ids_water_spring_only": sorted(int(i) for i in c.loc[c["spring_only"], "id"])}
    per = {}
    for k in CATEGORIES:
        lf = c["lex_cats"].map(lambda v, k=k: k in v)
        tr = c["cats_set"].map(lambda s, k=k: k in s)
        per[k] = {"precision": v44._rate(int((lf & tr).sum()), int(lf.sum())),
                  "recall_unweighted": v44._rate(int((lf & tr).sum()), int(tr.sum()))}
    res["by_description_category"] = per
    return res


def examples(d: pd.DataFrame) -> pd.DataFrame:
    c = d[d["status"] == "coded"].sort_values("id")
    rows = []
    for k in WET_SIDE:
        sub = c[c["lex_cats"].map(lambda v, k=k: k in v)]
        hit = sub[sub["cats_set"].map(lambda s, k=k: k in s)].head(2)
        miss = sub[~sub["cats_set"].map(lambda s, k=k: k in s)].head(1)
        for kind, part in (("agrees", hit), ("disagrees", miss)):
            rows += [{"lexicon_category": k, "kind": kind, "id": int(r.id), "name": r.name, "country": r.cc,
                      "reading": ";".join(r.lex_cats), "derivation_codes": r.cats, "quote": r.quote, "url": r.url}
                     for r in part.itertuples()]
    return pd.DataFrame(rows)


def example_aoi(lex: Lexicon) -> list[dict]:
    from pyproj import Transformer
    meta = json.loads((DATA_DIR / "global" / "labels" / EXAMPLE_AOI / "meta.json").read_text(encoding="utf-8"))
    gn = pd.read_parquet(GN, columns=["gid", "asciiname", "cc", "lat", "lon", "fclass"]).drop_duplicates("gid")
    gn = gn[gn["fclass"] == "P"]
    tr = Transformer.from_proj("+proj=longlat +datum=WGS84", meta["crs_proj4"], always_xy=True)
    gr = meta["grid"]
    x, y = tr.transform(gn["lon"].to_numpy(), gn["lat"].to_numpy())
    col, row = (x - gr["x0"]) / gr["cell_m"], (gr["y1"] - y) / gr["cell_m"]
    ok = (col >= 0) & (col < gr["nx"]) & (row >= 0) & (row < gr["ny"])
    out = []
    for r in gn.loc[ok].sort_values("asciiname").itertuples():
        m = lex.match(r.asciiname, r.cc)
        out.append({"name": r.asciiname, "cc": r.cc, "reading": [k for k in CATEGORIES if m[k]],
                    "wetside": any(m[k] for k in WET_SIDE),
                    "forms": sorted({f"{e.form} ({e.pos}, {e.cat})" for e in matched_entries(lex, r.asciiname, r.cc)})})
    return out


# %% 6. 人による照合の準備
def audit_sample(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    sec = {json.loads(x)["id"]: json.loads(x)["sections"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    cells = d.groupby(["stratum", "lang"]).size()
    raw = cells / cells.sum() * N_AUDIT  # 層 × 言語の件数に比例（最大剰余で丸める）
    alloc = np.floor(raw).astype(int)
    rest = N_AUDIT - int(alloc.sum())
    alloc[(raw - alloc).sort_values(ascending=False).index[:rest]] += 1
    rng = np.random.default_rng(SEED)
    picks = []
    for (s, lang), k in alloc.items():
        if k == 0:
            continue
        sub = d[(d["stratum"] == s) & (d["lang"] == lang)].sort_values("id")
        picks.append(sub.iloc[np.sort(rng.choice(len(sub), size=k, replace=False))])
    p = pd.concat(picks).reset_index(drop=True)
    p["audit_id"] = rng.permutation(len(p)) + 1
    p = p.sort_values("audit_id")
    sheet = pd.DataFrame({"audit_id": p["audit_id"], "name": p["name"], "country": p["cc"],
                          "articles": [" ".join(dict.fromkeys(s["url"] for s in sec[i])) for i in p["id"]],
                          "status": "", "categories": "", "quote": "", "comment": ""})
    key = pd.DataFrame({"audit_id": p["audit_id"], "v44_id": p["id"], "gid": p["gid"], "stratum": p["stratum"],
                        "lang": p["lang"], "v44_status": p["status"], "v44_categories": p["cats"]})
    return sheet, key, {"n": len(p), "allocation": {f"{s}/{lang}": int(k) for (s, lang), k in alloc.items() if k}}


GUIDE = """# Etymology audit: coding guide (amendment B-6, item 5)

Purpose: an independent check of the place-name readings used in the study. You will read the Wikipedia
article(s) listed for each settlement and record what the article says about the origin of the name.
Please do not consult the study's code repository, the paper or the study's own codes
(`v44_etym_codes.csv`) until your coding is complete; the codes to be compared are kept in a separate file.

For each row of `v51_audit_sample.csv`:

1. Open the article(s) in the `articles` column and find the section on the name (etymology, toponymy, name).
2. `status`:
   - `coded` — the article gives a derivation of the name (even if tentative);
   - `no-derivation` — the article has no statement on where the name comes from;
   - `unclear` — the article discusses the name but the derivation cannot be classified.
3. `categories` (only for `coded`; several allowed, separated by `;`). Code what the recorded derivation
   refers to, not what the modern place is like:
   - `water` — a river, stream, brook, lake, pond, spring or well;
   - `wet` — a marsh, fen, bog, swamp, mud, reeds, wet meadow;
   - `riverine` — an island in a river, a river meadow, floodplain, riverbank, sandbank, dyke or polder;
   - `crossing` — a ford, bridge, ferry, landing place or river port;
   - `low` — a valley, lowland, hollow, basin, or "lower";
   - `high` — a hill, mountain, height, rock, or "upper";
   - `other` — anything else (a person, saint, tribe, building, plant other than reeds, animal, etc.).
   A name taken from the name of a river counts as `water`.
4. `quote`: copy the shortest sentence(s) from the article that support the codes.
5. `comment`: optional (for example, competing derivations or a derivation you consider doubtful).

The sample is 100 names drawn at random within strata of the study's design (seed 20261028); the strata and
the study's codes are in `v51_audit_key.csv` and are used only after your coding is complete. Agreement will
be reported as Cohen's kappa for `status` and for each category.
"""


def main() -> None:
    lex = Lexicon(LEX)
    table, policy = lexicon_table(lex)
    table.to_csv(DOCS / "v51_lexicon_table.csv", index=False)
    stats, top, _ = match_stats(lex)
    top.to_csv(DOCS / "v51_top_forms.csv", index=False)
    d, sizes = coded()
    ex = examples(d)
    ex.to_csv(DOCS / "v51_examples.csv", index=False)
    sheet, key, audit = audit_sample(d)
    sheet.to_csv(DOCS / "v51_audit_sample.csv", index=False)
    key.to_csv(DOCS / "v51_audit_key.csv", index=False)
    (DOCS / "v51_audit_guide.md").write_text(GUIDE, encoding="utf-8")
    res = {"rule": "amendment B-6 item 5", "seed": SEED, "lexicon_sha256": lex.sha256,
           "lexicon": {"n_codes": len(table), "policy_note": policy,
                       "n_forms_by_category": {c: int(table[f"n_{c}"].sum()) for c in CATEGORIES}},
           "match": stats, "top_forms": top.to_dict(orient="records"),
           "precision_by_position": by_position(d, lex), "strict": strict(d, sizes, lex),
           "examples": ex.to_dict(orient="records"), "example_aoi": {EXAMPLE_AOI: example_aoi(lex)},
           "audit": audit}
    OUT.write_text(json.dumps(res, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: res[k] for k in ("lexicon",)}, ensure_ascii=False)[:600])
    m = res["match"]
    print("cells", m["n_cells"], "covered", round(m["share_cells_covered"], 3), "any", round(m["share_cells_any"], 3),
          "wetside", round(m["share_cells_wetside"], 3), "reproduces lexfix", m["reproduces_lexfix_wetside"])
    print(top.to_string(index=False))
    print("by position", {k: (v["any"]["rate"], v["any"]["n"]) for k, v in res["precision_by_position"].items()
                          if isinstance(v, dict)})
    s = res["strict"]
    print("strict", s["precision"], s["recall"], "spring-only water codes", s["n_water_spring_only"])
    print("audit", audit)


if __name__ == "__main__":
    main()
