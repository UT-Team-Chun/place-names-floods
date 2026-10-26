# %% v43（事後・amendment B-5 の 10）: 名前の読みを、人が記録した名前の由来（Wikidata P138）と照合する
"""確証判定は変えない。専門家による語源の照合の代わりに、人が記録した構造化データで辞書と言語モデルの読みを確かめる。

解析した標本の最寄りの集落（GeoNames の gid）を Wikidata の GeoNames ID（P1566）で引き、「名前の由来」（P138）が
何かを、その項目の分類（P31/P279*）で判定する:

- 水に関わる由来（hydro）: 水域・水路（body of water・watercourse・river・stream・spring・waterfall・water well・
  estuary）、湿地（wetland・marsh・swamp・bog・fen・floodplain）、島（island・river island）、渡り（ford・bridge）
- 記念・転用（commemorative）: 人（human）または他の集落（human settlement）で、水に関わる由来を持たないもの

報告: 由来の記録がある集落の数、水に関わる由来の集落のうち辞書・言語モデルが水系側と読んだ割合（再現率の目安）、
水系側と読んだ名前のうち由来が人・他の集落だった割合（読み違いの目安）。国別にも。人の記録がある名前（主に
著名な集落）に限られる偏りがあり、専門家による語源の照合ではない。LLM は使わない（構造化データだけ）。

出力: stdout の表と docs/verification/v43_wikidata_etymology.json（問い合わせの結果は
data/verification/v43_wikidata_raw.json にキャッシュ）
"""

from __future__ import annotations

import json
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

ENDPOINT = "https://query.wikidata.org/sparql"
UA = "toponym-semantics research (https://github.com/UT-Team-Chun/place-names-floods)"
CACHE = DATA_DIR / "verification" / "v43_wikidata_raw.json"
OUT = ROOT / "docs" / "verification" / "v43_wikidata_etymology.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
ROOTS = {
    "water": ["Q15324", "Q355304", "Q4022", "Q47521", "Q124714", "Q34038", "Q43483", "Q47053"],
    "wetland": ["Q170321", "Q30198", "Q166735", "Q1681353", "Q514050", "Q193110"],
    "island": ["Q23442", "Q162602"],
    "crossing": ["Q12743", "Q12280"],
    "human": ["Q5"],
    "settlement": ["Q486972"],
}
HYDRO = ("water", "wetland", "island", "crossing")
NREAD_WETSIDE = ("water", "wetland", "riverine", "crossing", "lowland")
BATCH = 300


def sparql(q: str) -> list[dict]:
    url = f"{ENDPOINT}?{urllib.parse.urlencode({'query': q, 'format': 'json'})}"
    req = urllib.request.Request(url, headers={"User-Agent": UA})  # noqa: S310 — 固定の https
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:  # noqa: S310
                return json.load(r)["results"]["bindings"]
        except OSError:
            if attempt == 5:
                raise
            time.sleep(10 * (attempt + 1))
    return []


def fetch(gids: list[int]) -> dict:
    if CACHE.is_file():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    named: dict[str, list[str]] = {}
    items: dict[str, str] = {}
    for i in range(0, len(gids), BATCH):
        vals = " ".join(f'"{g}"' for g in gids[i:i + BATCH])
        rows = sparql(f"SELECT ?gn ?item ?na WHERE {{ VALUES ?gn {{ {vals} }} ?item wdt:P1566 ?gn . "
                      f"OPTIONAL {{ ?item wdt:P138 ?na }} }}")
        for b in rows:
            gn = b["gn"]["value"]
            items[gn] = b["item"]["value"].rsplit("/", 1)[1]
            if "na" in b:
                named.setdefault(gn, []).append(b["na"]["value"].rsplit("/", 1)[1])
        time.sleep(1.0)
        print(f"items {i + BATCH}/{len(gids)}: {len(items)} items, {len(named)} with P138", flush=True)
    targets = sorted({t for v in named.values() for t in v})
    cats: dict[str, list[str]] = {t: [] for t in targets}
    for cat, roots in ROOTS.items():
        rv = " ".join(f"wd:{r}" for r in roots)
        for i in range(0, len(targets), BATCH):
            tv = " ".join(f"wd:{t}" for t in targets[i:i + BATCH])
            rows = sparql(f"SELECT DISTINCT ?t WHERE {{ VALUES ?t {{ {tv} }} VALUES ?root {{ {rv} }} "
                          f"?t wdt:P31/wdt:P279* ?root }}")
            for b in rows:
                cats[b["t"]["value"].rsplit("/", 1)[1]].append(cat)
            time.sleep(1.0)
        print(f"classified {cat}", flush=True)
    raw = {"items": items, "named_after": named, "target_categories": cats}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(raw), encoding="utf-8")
    return raw


def wilson(k: int, n: int) -> list[float] | None:
    if n == 0:
        return None
    z, p = 1.96, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [float(max(0.0, c - h)), float(min(1.0, c + h))]


def rate(mask: np.ndarray, sub: np.ndarray) -> dict:
    n, k = int(sub.sum()), int((mask & sub).sum())
    return {"k": k, "n": n, "rate": k / n if n else None, "ci95": wilson(k, n)}


def main() -> None:
    t = pd.read_parquet(DATA_DIR / "verification" / "v37_table.parquet", columns=["P_gid1"])
    gids = sorted(int(g) for g in t["P_gid1"].unique() if g >= 0)
    raw = fetch(gids)
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    lex = Lexicon(LEX)  # B-4 の読み方（ノルウェーを含む）
    nr = pd.read_parquet(DATA_DIR / "global" / "nread" / "nread.parquet", columns=["gid", "category"])
    nr = nr.drop_duplicates("gid").set_index("gid")["category"]
    rows = []
    for g in gids:
        if g not in gn.index:
            continue
        name, cc = str(gn.at[g, "asciiname"]), str(gn.at[g, "cc"])
        m = lex.match(name, cc)
        cats = {c for tgt in raw["named_after"].get(str(g), []) for c in raw["target_categories"].get(tgt, [])}
        rows.append({"gid": g, "cc": cc, "has_item": str(g) in raw["items"],
                     "has_p138": str(g) in raw["named_after"],
                     "hydro": bool(cats & set(HYDRO)),
                     "commemorative": bool(cats & {"human", "settlement"}) and not (cats & set(HYDRO)),
                     **{f"na_{c}": c in cats for c in ROOTS},
                     "lex_wetside": any(m[c] for c in WET_SIDE),
                     "nread_wetside": str(nr.get(g, "")) in NREAD_WETSIDE})
    d = pd.DataFrame(rows)
    hy, cm, p138 = d["hydro"].to_numpy(), d["commemorative"].to_numpy(), d["has_p138"].to_numpy()
    lw, nw = d["lex_wetside"].to_numpy(), d["nread_wetside"].to_numpy()
    res = {"n_settlements": len(d), "n_with_item": int(d["has_item"].sum()), "n_with_p138": int(p138.sum()),
           "n_hydro": int(hy.sum()), "n_commemorative": int(cm.sum()),
           "recall_lexicon": rate(lw, hy), "recall_nread": rate(nw, hy),
           "commemorative_among_lexicon_wetside": rate(cm, lw & p138),
           "commemorative_among_nread_wetside": rate(cm, nw & p138),
           "hydro_among_lexicon_wetside": rate(hy, lw & p138),
           "hydro_among_nread_wetside": rate(hy, nw & p138),
           "hydro_among_not_wetside_lexicon": rate(hy, ~lw & p138),
           "by_country": {}}
    for cc, s in d.groupby("cc"):
        h = s["hydro"].to_numpy()
        if h.sum() >= 30:
            res["by_country"][cc] = {"n_hydro": int(h.sum()),
                                     "recall_lexicon": rate(s["lex_wetside"].to_numpy(), h),
                                     "recall_nread": rate(s["nread_wetside"].to_numpy(), h)}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "by_country"}, indent=1, default=float))
    for cc, v in sorted(res["by_country"].items(), key=lambda kv: -kv[1]["n_hydro"]):
        print(cc, v["n_hydro"], f"lex {v['recall_lexicon']['rate']:.2f}", f"nread {v['recall_nread']['rate']:.2f}")


if __name__ == "__main__":
    main()
