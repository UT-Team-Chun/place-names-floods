# %% v44（事後・amendment B-5 の 10、追記 2）: 名前の読みを Wikipedia の語源の記述と照合する
"""確証判定は変えない。専門家による語源の照合の代わりに、人が書いた語源の記述で辞書と言語モデルの読みを確かめる。

手順（B-5 追記 2 で計算の前に固定）:

1. ``--fetch``: 解析した表の最寄りの集落のうち、Wikidata の項目（v43 のキャッシュ）と、英語版か国の辞書の最初の
   言語の版の記事があるものを母集団にし、辞書が水系側と読んだ層と読まなかった層から 300 ずつ無作為に選ぶ
   （seed 20261026）。語源・名前の由来の見出しの節を取り、通し番号以外を伏せた判定用のファイルを書く
2. ``--show START N``: 判定用のファイルを番号順に表示する（読み・層・国は出さない）
3. 判定は docs/verification/v44_etym_codes.csv に書く（id, status, cats, quote, url）。status は
   coded / no-derivation / unclear、cats は water・wet・riverine・crossing・low・high・other の集合（; 区切り）
4. ``--score``: 判定と読みを結び、層の重みを付けた精度・再現率を出す

出力: docs/verification/v44_wikipedia_etymology.json と stdout の表。取得した記述は data/verification/ に置く。
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR
from topsem.world.lexicon import WET_SIDE, Lexicon

UA = "toponym-semantics research (https://github.com/UT-Team-Chun/place-names-floods)"
V43 = DATA_DIR / "verification" / "v43_wikidata_raw.json"
VDIR = DATA_DIR / "verification"
SITELINKS = VDIR / "v44_sitelinks.json"
SAMPLE = VDIR / "v44_sample.json"  # 層・読みを含む（判定中は開かない）
TO_CODE = VDIR / "v44_to_code.jsonl"  # 伏せた判定用
CODES = ROOT / "docs" / "verification" / "v44_etym_codes.csv"
OUT = ROOT / "docs" / "verification" / "v44_wikipedia_etymology.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
SEED = 20261026
PER_STRATUM = 300  # 追記 2（最初の取得）
SEED2 = 20261027  # 追記 4: 標本を広げる（水系側の層は全数、他の層は 1,500）
SIZES = {"wet": None, "not": 1500}
ROUND1 = VDIR / "v44_sample_round1.json"
ROUND1_CODE = VDIR / "v44_to_code_round1.jsonl"
MAX_CHARS = 2000
NREAD_WETSIDE = ("water", "wetland", "riverine", "crossing", "lowland")
CODE_CATS = ("water", "wet", "riverine", "crossing", "low", "high", "other")
MIN_GROUP = 30
N_BOOT = 2000

# 見出しの語（小文字）。見出しと一致するか、見出しがこの語と区切りで始まれば語源の節とみなす
HEAD = {
    "en": ["etymology", "toponymy", "name", "names", "naming", "origin of the name", "place name",
           "place-name", "name origin", "the name", "name and etymology"],
    "de": ["name", "namensherkunft", "etymologie", "ortsname", "toponymie", "namensgebung",
           "herkunft des namens", "namensdeutung", "ortsnamen"],
    "fr": ["toponymie", "étymologie", "nom", "origine du nom", "dénomination"],
    "es": ["toponimia", "etimología", "nombre", "origen del nombre", "denominación"],
    "it": ["toponimo", "toponomastica", "etimologia", "origine del nome", "nome"],
    "pt": ["toponímia", "etimologia", "origem do nome", "nome", "topónimo", "topônimo"],
    "nl": ["etymologie", "naam", "toponymie", "naamgeving", "herkomst van de naam",
           "naamsverklaring"],
    "pl": ["nazwa", "etymologia", "toponimia", "historia nazwy", "pochodzenie nazwy"],
    "cs": ["název", "etymologie", "původ názvu", "jméno"],
    "sh": ["ime", "naziv", "etimologija", "toponimija", "porijeklo imena", "poreklo imena",
           "porijeklo naziva", "poreklo naziva", "toponim"],
    "sl": ["ime", "imenoslovje", "etimologija", "izvor imena", "toponim", "poimenovanje"],
    "bg": ["име", "название", "етимология", "произход на името", "топоним", "наименование"],
    "uk": ["назва", "етимологія", "походження назви", "топонімія", "топонім"],
    "ro": ["toponimie", "toponimia", "etimologie", "nume", "denumire", "originea numelui"],
    "hu": ["neve", "nevének eredete", "etimológia", "név", "a név eredete", "helynév"],
    "sv": ["namnet", "etymologi", "namn", "ortnamnet", "namnets ursprung"],
    "no": ["navn", "navnet", "etymologi", "stedsnavnet", "navnets opprinnelse"],
    "da": ["navn", "navnet", "etymologi", "bynavnet", "stednavnet"],
    "fi": ["nimi", "nimen alkuperä", "etymologia", "nimistö"],
    "et": ["nimi", "nime päritolu", "etümoloogia"],
    "lv": ["nosaukums", "etimoloģija", "nosaukuma izcelsme"],
    "lt": ["vardas", "pavadinimas", "etimologija", "pavadinimo kilmė", "vardo kilmė"],
    "el": ["ονομασία", "ετυμολογία", "τοπωνύμιο", "όνομα", "προέλευση του ονόματος",
           "τοπωνυμία"],
    "tr": ["etimoloji", "adı", "isim", "adının kökeni", "ad"],
    "id": ["etimologi", "asal-usul nama", "toponimi", "nama", "asal usul nama", "asal nama"],
    "vi": ["từ nguyên", "tên gọi", "địa danh", "nguồn gốc tên gọi"],
    "tl": ["etimolohiya", "pangalan", "pinagmulan ng pangalan"],
    "sw": ["jina", "asili ya jina", "etimolojia"],
    "mg": ["anarana"],
    "so": ["magaca", "magac"],
    "sq": ["emri", "etimologjia", "toponimia", "origjina e emrit"],
    "am": ["ስያሜ", "ስም"],
    "ar": ["التسمية", "أصل التسمية", "الاسم", "سبب التسمية", "أصل الاسم", "تسمية"],
    "fa": ["وجه تسمیه", "وجه‌تسمیه", "نام", "ریشه نام", "ریشه‌شناسی"],
    "hi": ["नाम", "व्युत्पत्ति", "नामकरण", "नाम की उत्पत्ति"],
    "ka": ["სახელწოდება", "ეტიმოლოგია", "ტოპონიმი"],
    "kk": ["атауы", "этимологиясы", "атау"],
    "km": ["ឈ្មោះ"],
    "my": ["အမည်"],
    "th": ["ที่มาของชื่อ", "ชื่อ", "นิรุกติศาสตร์"],
    "zh": ["名称", "名稱", "地名由来", "地名由來", "名称由来", "名稱由來", "词源", "詞源", "地名"],
}
# 見出しのどこにあっても語源の節とみなす語の一部
STRONG = ("etym", "etim", "toponym", "toponim", "этимол", "етимол", "ετυμολ", "τοπων",
          "топоним", "топонім", "namensherkunft", "herkunft des namens", "origin of the name",
          "origine du nom", "origen del nombre", "origine del nome", "origem do nome",
          "pochodzenie nazwy", "původ názvu", "nimen alkuperä")
_HDR = re.compile(r"^(={2,6})\s*(.*?)\s*\1\s*$")
# 見出しが語の後に続けてよいもの（「Name and etymology」など。「Names of notable people」は取らない）
_JOIN = r"and |und |et |y |e |en |i |og |och |ja |és |ve |dan |và |na |,|/|–|-|\(|:"


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except OSError:
            if attempt == 5:
                raise
            time.sleep(5 * (attempt + 1))
    return {}


def is_etym_heading(head: str, lang: str) -> bool:
    h = head.casefold().strip(" :.")
    if any(s in h for s in STRONG):
        return True
    return any(h == kw or re.match(rf"{re.escape(kw)}\s*({_JOIN})", h)
               for kw in HEAD.get(lang, []))


def etym_sections(text: str, lang: str) -> list[tuple[str, str]]:
    """プレーンテキストの記事（== 見出し == 付き）→ 語源の見出しの節（下位の節を含む）の (見出し, 本文)。"""
    out: list[tuple[str, str]] = []
    cur: tuple[str, int] | None = None
    buf: list[str] = []
    for line in text.splitlines():
        m = _HDR.match(line)
        if m:
            level = len(m.group(1))
            if cur is not None and level <= cur[1]:
                out.append((cur[0], "\n".join(buf).strip()))
                cur, buf = None, []
            if cur is None and is_etym_heading(m.group(2), lang):
                cur, buf = (m.group(2), level), []
            elif cur is not None:
                buf.append(m.group(2))
            continue
        if cur is not None:
            buf.append(line)
    if cur is not None:
        out.append((cur[0], "\n".join(buf).strip()))
    return [(h, b) for h, b in out if b]


def sitelinks(qids: list[str]) -> dict[str, dict[str, str]]:
    if SITELINKS.is_file():
        return json.loads(SITELINKS.read_text(encoding="utf-8"))
    res: dict[str, dict[str, str]] = {}
    for i in range(0, len(qids), 50):
        q = urllib.parse.urlencode({"action": "wbgetentities", "ids": "|".join(qids[i:i + 50]),
                                    "props": "sitelinks", "format": "json"})
        d = _get(f"https://www.wikidata.org/w/api.php?{q}")
        for qid, ent in (d.get("entities") or {}).items():
            res[qid] = {k: v["title"] for k, v in (ent.get("sitelinks") or {}).items()
                        if k.endswith("wiki") and not k.startswith(("commons", "species"))}
        time.sleep(0.3)
        if i % 2000 == 0:
            print(f"sitelinks {i}/{len(qids)}", flush=True)
    SITELINKS.write_text(json.dumps(res), encoding="utf-8")
    return res


def extract(lang: str, title: str) -> str:
    q = urllib.parse.urlencode({"action": "query", "prop": "extracts", "explaintext": 1,
                                "exsectionformat": "wiki", "redirects": 1, "titles": title,
                                "format": "json", "formatversion": 2})
    d = _get(f"https://{lang}.wikipedia.org/w/api.php?{q}")
    pages = (d.get("query") or {}).get("pages") or []
    return str(pages[0].get("extract") or "") if pages else ""


def population(lex: Lexicon) -> pd.DataFrame:
    t = pd.read_parquet(DATA_DIR / "verification" / "v37_table.parquet", columns=["P_gid1"])
    gids = sorted(int(g) for g in t["P_gid1"].unique() if g >= 0)
    raw = json.loads(V43.read_text(encoding="utf-8"))
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    nr = pd.read_parquet(DATA_DIR / "global" / "nread" / "nread.parquet", columns=["gid", "category"])
    nr = nr.drop_duplicates("gid").set_index("gid")["category"]
    rows = []
    for g in gids:
        qid = raw["items"].get(str(g))
        if qid is None or g not in gn.index:
            continue
        name, cc = str(gn.at[g, "asciiname"]), str(gn.at[g, "cc"])
        m = lex.match(name, cc)
        langs = lex.country_langs.get(cc.upper(), ())
        rows.append({"gid": g, "qid": qid, "name": name, "cc": cc,
                     "lang": langs[0] if langs else "",
                     "lex_wetside": any(m[c] for c in WET_SIDE),
                     "lex_cats": [c for c in m if m[c]],
                     "nread_wetside": str(nr.get(g, "")) in NREAD_WETSIDE,
                     "nread_cat": str(nr.get(g, ""))})
    return pd.DataFrame(rows)


def fetch() -> None:
    import datetime as dt
    import shutil
    lex = Lexicon(LEX)  # B-4 の読み方
    pop = population(lex)
    cache: dict[int, list] = {}
    if SAMPLE.is_file():  # 最初の取得（追記 2）の抜粋を使い回す。最初の標本の記録は残す
        old = json.loads(SAMPLE.read_text(encoding="utf-8"))
        if old.get("round", 1) == 1 and not ROUND1.is_file():
            shutil.copyfile(SAMPLE, ROUND1)
            shutil.copyfile(TO_CODE, ROUND1_CODE)
        cache = {int(r["gid"]): r["sections"] for r in old["records"]}
    links = sitelinks(sorted(pop["qid"].unique()))
    pop["en_title"] = [links.get(q, {}).get("enwiki") for q in pop["qid"]]
    pop["local_title"] = [links.get(q, {}).get(f"{lg}wiki") if lg and lg != "en" else None
                          for q, lg in zip(pop["qid"], pop["lang"], strict=True)]
    pop = pop[pop["en_title"].notna() | pop["local_title"].notna()].sort_values("gid")
    rng = np.random.default_rng(SEED2)
    strata = {"wet": pop[pop["lex_wetside"]], "not": pop[~pop["lex_wetside"]]}
    picks = []
    for s, d in strata.items():
        n = SIZES[s]
        idx = np.arange(len(d)) if n is None or n >= len(d) else rng.choice(len(d), size=n, replace=False)
        picks.append(d.iloc[np.sort(idx)].assign(stratum=s))
    samp = pd.concat(picks, ignore_index=True)
    recs = []
    for i, r in samp.iterrows():
        if int(r["gid"]) in cache:
            recs.append({**{k: r[k] for k in ("gid", "qid", "name", "cc", "lang", "stratum",
                                               "lex_wetside", "nread_wetside", "nread_cat")},
                         "lex_cats": list(r["lex_cats"]), "sections": cache[int(r["gid"])]})
            continue
        secs = []
        for lang, title in (("en", r["en_title"]), (r["lang"], r["local_title"])):
            if not title:
                continue
            for head, body in etym_sections(extract(lang, title), lang):
                secs.append({"lang": lang, "title": title, "heading": head,
                             "url": f"https://{lang}.wikipedia.org/wiki/"
                                    f"{urllib.parse.quote(title.replace(' ', '_'))}",
                             "text": body[:MAX_CHARS]})
            time.sleep(0.2)
        recs.append({**{k: r[k] for k in ("gid", "qid", "name", "cc", "lang", "stratum",
                                           "lex_wetside", "nread_wetside", "nread_cat")},
                     "lex_cats": list(r["lex_cats"]), "sections": secs})
        if i % 50 == 0:
            print(f"extracts {i}/{len(samp)}", flush=True)
    sizes = {s: len(d) for s, d in strata.items()}
    now = dt.datetime.now(dt.UTC).strftime("%Y-%m-%d")
    SAMPLE.write_text(json.dumps({"round": 2, "stratum_sizes": sizes, "n_population": len(pop),
                                  "sample_sizes": {s: int((samp["stratum"] == s).sum()) for s in strata},
                                  "retrieved_utc": now, "records": recs}, default=bool), encoding="utf-8")
    have = [r for r in recs if r["sections"]]
    order = np.random.default_rng(SEED2 + 1).permutation(len(have))
    with TO_CODE.open("w", encoding="utf-8") as f:
        for k, j in enumerate(order, start=1):
            f.write(json.dumps({"id": k, "key": int(have[j]["gid"]),
                                "sections": [{x: s[x] for x in ("url", "heading", "text")}
                                             for s in have[j]["sections"]]},
                               ensure_ascii=False) + "\n")
    print(json.dumps({"population": len(pop), "strata": sizes, "sampled": len(samp),
                      "from_cache": int(sum(int(g) in cache for g in samp["gid"])), "with_section": len(have)}))


def show(start: int, n: int) -> None:
    for line in TO_CODE.read_text(encoding="utf-8").splitlines()[start - 1:start - 1 + n]:
        r = json.loads(line)
        print(f"### {r['id']}")
        for s in r["sections"]:
            print(f"[{s['url']} | {s['heading']}]\n{s['text']}\n")


def _rate(k: int, n: int) -> dict:
    if n == 0:
        return {"k": 0, "n": 0, "rate": None, "ci95": None}
    z, p = 1.96, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return {"k": k, "n": n, "rate": p, "ci95": [float(max(0.0, c - h)), float(min(1.0, c + h))]}


def weighted(d: pd.DataFrame, sizes: dict[str, int], flag: str) -> dict[str, float]:
    """層の大きさで重みを付けた、flag の精度と再現率（正解 = 記述が水系側）。"""
    w = {s: sizes[s] / max(1, int((d["stratum"] == s).sum())) for s in sizes}
    wt = d["stratum"].map(w).to_numpy(dtype=float)
    y, f = d["truth"].to_numpy(dtype=bool), d[flag].to_numpy(dtype=bool)
    tp = float((wt * (y & f)).sum())
    return {"precision": tp / float((wt * f).sum()) if f.any() else float("nan"),
            "recall": tp / float((wt * y).sum()) if y.any() else float("nan")}


def score() -> None:
    samp = json.loads(SAMPLE.read_text(encoding="utf-8"))
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"]
           for x in TO_CODE.read_text(encoding="utf-8").splitlines()}
    codes = pd.read_csv(CODES)
    codes["gid"] = codes["id"].map(key)
    d = codes.join(meta, on="gid")
    counts = d["status"].value_counts().to_dict()
    d = d[d["status"] == "coded"].copy()
    d["cats_set"] = d["cats"].fillna("").map(lambda s: {c for c in s.split(";") if c})
    bad = {c for s in d["cats_set"] for c in s} - set(CODE_CATS)
    if bad:
        msg = f"unknown code categories {bad}"
        raise SystemExit(msg)
    d["truth"] = d["cats_set"].map(lambda s: bool(s & set(WET_SIDE)))
    res: dict = {"n_population": samp["n_population"], "stratum_sizes": sizes,
                 "n_sampled": len(samp["records"]), "sample_sizes": samp.get("sample_sizes"),
                 "retrieved_utc": samp.get("retrieved_utc"), "round": samp.get("round", 1),
                 "n_with_section": len(key), "status_counts": counts,
                 "truth_rate_in_wet": _rate(int(d.loc[d["stratum"] == "wet", "truth"].sum()),
                                            int((d["stratum"] == "wet").sum())),
                 "truth_rate_in_not": _rate(int(d.loc[d["stratum"] == "not", "truth"].sum()),
                                            int((d["stratum"] == "not").sum()))}
    if ROUND1.is_file():  # 追記 4: 最初の標本の件数（広げた理由として報告する。判定はしていない）
        r1 = json.loads(ROUND1.read_text(encoding="utf-8"))["records"]
        res["round1"] = {"n_sampled": len(r1), "n_with_section": sum(1 for r in r1 if r.get("sections"))}
    for flag in ("lex_wetside", "nread_wetside"):
        point = weighted(d, sizes, flag)
        rng = np.random.default_rng(SEED)
        boots = []
        parts = {s: d[d["stratum"] == s] for s in sizes}
        for _ in range(N_BOOT):
            bb = pd.concat([p.iloc[rng.integers(0, len(p), len(p))] for p in parts.values()])
            boots.append(weighted(bb, sizes, flag))
        res[flag] = {k: {"est": point[k],
                         "ci95": [float(np.nanpercentile([b[k] for b in boots], 2.5)),
                                  float(np.nanpercentile([b[k] for b in boots], 97.5))]}
                     for k in ("precision", "recall")}
    res["by_lexicon_category"] = {}
    for c in (*WET_SIDE, "high"):  # 追記 2 の 4: 辞書の区分ごと（水系側でない high も）
        sub = d[d["lex_cats"].map(lambda v, c=c: c in v)]
        res["by_lexicon_category"][c] = {
            "same_category": _rate(int(sub["cats_set"].map(lambda s, c=c: c in s).sum()),
                                   len(sub)),
            "any_wetside": _rate(int(sub["truth"].sum()), len(sub))}
    res["by_language"] = {}
    for lang, s in d.groupby("lang"):
        if len(s) >= MIN_GROUP:
            res["by_language"][lang] = {
                "n": len(s),
                "precision_lex": _rate(int((s["truth"] & s["lex_wetside"]).sum()),
                                       int(s["lex_wetside"].sum())),
                "truth_rate_not_flagged": _rate(int((s["truth"] & ~s["lex_wetside"]).sum()),
                                                int((~s["lex_wetside"]).sum()))}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps(res, indent=1, default=float))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--show", nargs=2, type=int, metavar=("START", "N"))
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()
    if a.fetch:
        fetch()
    if a.show:
        show(*a.show)
    if a.score:
        score()


if __name__ == "__main__":
    main()
