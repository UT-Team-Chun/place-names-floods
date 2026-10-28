# %% v57b（事後・amendment B-8 の B13）: 二人目の符号化者（別の開発元の言語モデル）
"""B-5 の Wikipedia の抜粋（名の節がある集落、v44_to_code.jsonl。辞書の読み・層・最初の符号は含まない）を、辞書の
作成にも最初の符号化にも関わっていない別の開発元のモデル（OpenAI の Codex の既定のモデル、codex-companion の task）に、
人の符号化者のための手引き（v51_audit_guide.md）と同じ規則で符号化させる。手引きに足す符号: 名が他の場所や
この集落に無い地物から取られた・記念の名か（transferred）、名の指すものがこの場所にあるか（at_site）、水の要素が
泉・井戸だけか（spring_or_well_only、厳密な読みの採点に使う）。

抜粋は 20 ずつ固定の prompt で渡し（道具・ファイル・ウェブを使わないよう指示し、空の作業フォルダで走らせる）、
出力は JSON の配列。回ごとの結果は data/verification/v57b/batch_<k>.json に残して使い回す。

    uv run python verification/r_okauchi/v57b_second_coder.py --code   # 符号化（4 本並行）
    uv run python verification/r_okauchi/v57b_second_coder.py --score  # 最初の符号と比べる

出力: docs/verification/v57b_second_coder.json と docs/verification/v57b_codes_coder2.csv
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR
from topsem.world.b6 import STRICT_CATEGORIES, strict_lexicon
from topsem.world.lexicon import WET_SIDE, Lexicon

TO_CODE = DATA_DIR / "verification" / "v44_to_code.jsonl"
CACHE = DATA_DIR / "verification" / "v57b"
OUT = ROOT / "docs" / "verification" / "v57b_second_coder.json"
CODES2 = ROOT / "docs" / "verification" / "v57b_codes_coder2.csv"
COMPANION = Path.home() / ".claude" / "plugins" / "cache" / "openai-codex" / "codex" / "1.0.6" / "scripts" / \
    "codex-companion.mjs"
BATCH = 20
MAX_CHARS = 2000  # 最初の符号化と同じ（v44: 記事あたり 2000 文字まで）
CATS = ("water", "wet", "riverine", "crossing", "low", "high", "other")
STRICT = ("water", "wet", "riverine")

PROMPT = """You are coding what Wikipedia says about the ORIGIN of place names. Do not use any tool, file, command or web \
access; use only the excerpts below.
For each item, return one JSON object with keys:
- "id": the item id;
- "status": "coded" if the excerpt gives a derivation of the name (even tentative), "no-derivation" if it says nothing \
about where the name comes from, "unclear" if it discusses the name but the derivation cannot be classified;
- "categories": for "coded" only, a list of what the recorded derivation refers to (several allowed): "water" (river, \
stream, brook, lake, pond, spring or well; a name taken from a river's name counts as water), "wet" (marsh, fen, bog, \
swamp, mud, reeds, wet meadow), "riverine" (island in a river, river meadow, floodplain, riverbank, sandbank, dyke or \
polder), "crossing" (ford, bridge, ferry, landing place, river port), "low" (valley, lowland, hollow, basin, or \
"lower"), "high" (hill, mountain, height, rock, or "upper"), "other" (anything else: person, saint, tribe, building, \
plant other than reeds, animal, etc.). Code an element only when the excerpt states its meaning. Code what the \
derivation refers to, not what the modern place is like;
- "spring_or_well_only": "yes" if "water" is coded and the only water element is a spring or a well, otherwise "no";
- "transferred": "yes" if the excerpt says the settlement took its name from another place, or from a feature that is \
not at the settlement (for example a town named after a lake some distance away), or that the name is commemorative; \
"no" if the name describes the settlement's own site or a feature at it; "unknown" otherwise;
- "at_site": "yes" if the thing the name refers to (the river, marsh, ford, valley ...) is at or beside the settlement, \
"no" if the excerpt places it elsewhere, "unknown" otherwise (also "unknown" when nothing in the name refers to a place \
feature);
- "quote": the shortest sentence(s) from the excerpt supporting the codes.
Return ONLY a JSON array of these objects, nothing else.

Items:
"""


def batches() -> list[list[dict]]:
    rows = [json.loads(x) for x in TO_CODE.read_text(encoding="utf-8").splitlines() if x.strip()]
    items = [{"id": r["id"], "sections": [{"heading": s.get("heading"), "text": s["text"][:MAX_CHARS]}
                                          for s in r["sections"]]} for r in rows]
    return [items[i:i + BATCH] for i in range(0, len(items), BATCH)]


def _parse(text: str, ids: list[int]) -> list[dict]:
    m = re.search(r"\[\s*\{.*\}\s*\]", text, flags=re.S)
    if not m:
        msg = "no JSON array in the output"
        raise ValueError(msg)
    arr = json.loads(m.group(0))
    got = sorted(int(x["id"]) for x in arr)
    if got != sorted(ids):
        msg = f"ids differ: {got} vs {sorted(ids)}"
        raise ValueError(msg)
    return arr


def code_batch(k: int, items: list[dict]) -> str:
    dst = CACHE / f"batch_{k:03d}.json"
    if dst.is_file():
        return f"{k}: cached"
    prompt = PROMPT + json.dumps(items, ensure_ascii=False, indent=0)
    ids = [int(x["id"]) for x in items]
    for attempt in range(3):
        with tempfile.TemporaryDirectory() as tmp:  # 空の作業フォルダ（研究のファイルを見えなくする）
            p = subprocess.run(["node", str(COMPANION), "task", "--effort", "medium", prompt],  # noqa: S603, S607
                               cwd=tmp, capture_output=True, text=True, check=False, timeout=1800)
        try:
            arr = _parse(p.stdout, ids)
        except (ValueError, json.JSONDecodeError) as e:
            print(f"{k}: attempt {attempt} failed ({e})", flush=True)
            continue
        dst.write_text(json.dumps({"batch": k, "ids": ids, "codes": arr, "attempt": attempt,
                                   "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1,
                                  ensure_ascii=False), encoding="utf-8")
        return f"{k}: ok (attempt {attempt})"
    return f"{k}: FAILED"


def code_all(workers: int) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    bs = batches()
    with ThreadPoolExecutor(workers) as ex:
        for r in ex.map(lambda kb: code_batch(*kb), enumerate(bs)):
            print(r, flush=True)


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _kappa(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a), np.asarray(b)
    labs = np.union1d(a, b)
    po = float((a == b).mean())
    pe = float(sum((a == v).mean() * (b == v).mean() for v in labs))
    return (po - pe) / (1 - pe) if pe < 1 else float("nan")


def score() -> None:
    v44 = _load("v44_wikipedia_etymology")
    c2 = []
    for f in sorted(CACHE.glob("batch_*.json")):
        c2 += json.loads(f.read_text(encoding="utf-8"))["codes"]
    c2 = pd.DataFrame(c2)
    c2["id"] = c2["id"].astype(int)
    c2["cats2"] = c2["categories"].map(lambda v: set(v) if isinstance(v, list) else set())
    c2.assign(categories=c2["cats2"].map(lambda s: ";".join(sorted(s))))[
        ["id", "status", "categories", "spring_or_well_only", "transferred", "at_site", "quote"]].to_csv(
        CODES2, index=False)
    samp = json.loads(v44.SAMPLE.read_text(encoding="utf-8"))
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in TO_CODE.read_text(encoding="utf-8").splitlines()}
    c1 = pd.read_csv(v44.CODES)
    c1["cats1"] = c1["cats"].fillna("").map(lambda s: {c for c in s.split(";") if c})
    d = c1.merge(c2, on="id", how="inner", suffixes=("_1", "_2"))
    d["gid"] = d["id"].map(key)
    d = d.join(meta, on="gid")
    both = d[(d["status_1"] == "coded") & (d["status_2"] == "coded")].copy()
    out: dict = {"n_items": len(d), "n_coder2": len(c2), "status_kappa": _kappa(d["status_1"], d["status_2"]),
                 "status_agreement": float((d["status_1"] == d["status_2"]).mean()),
                 "n_both_coded": len(both)}
    out["category_kappa"] = {c: _kappa(both["cats1"].map(lambda s, c=c: c in s), both["cats2"].map(lambda s, c=c: c in s))
                             for c in CATS}
    w1 = both["cats1"].map(lambda s: bool(s & set(WET_SIDE)))
    w2 = both["cats2"].map(lambda s: bool(s & set(WET_SIDE)))
    out["water_related_kappa"] = _kappa(w1, w2)
    out["water_related_agreement"] = float((w1 == w2).mean())
    # 二人目の符号での精度・再現率（最初の符号と同じ重み、coded の行だけ）
    d2 = d[d["status_2"] == "coded"].copy()
    d2["truth"] = d2["cats2"].map(lambda s: bool(s & set(WET_SIDE)))
    # 厳密: v51 と同じく、厳密な読みの辞書の印と、泉・井戸だけの水を除いた水域・湿地・河川の派生
    ls = strict_lexicon(Lexicon(ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"))
    d2["strict_flag"] = [any(ls.match(r.name, r.cc)[k] for k in STRICT_CATEGORIES) for r in d2.itertuples()]
    d2["truth_strict"] = [bool((s - ({"water"} if sw == "yes" else set())) & set(STRICT_CATEGORIES))
                          for s, sw in zip(d2["cats2"], d2["spring_or_well_only"], strict=True)]
    res = {"broad": {flag: v44.weighted(d2, sizes, flag) for flag in ("lex_wetside", "nread_wetside")},
           "strict": {"strict_flag": v44.weighted(d2.assign(truth=d2["truth_strict"]), sizes, "strict_flag")}}
    out["coder2_precision_recall"] = res
    d1 = d[d["status_1"] == "coded"].copy()
    d1["truth"] = d1["cats1"].map(lambda s: bool(s & set(WET_SIDE)))
    out["coder1_precision_recall_same_rows"] = {flag: v44.weighted(d1, sizes, flag)
                                                for flag in ("lex_wetside", "nread_wetside")}
    agree = d2[d2["lex_wetside"] & d2["truth"]]
    out["transferred_among_lexicon_wetside_true"] = {"yes": int((agree["transferred"] == "yes").sum()),
                                                     "n": len(agree)}
    out["at_site_no_among_lexicon_wetside_true"] = {"no": int((agree["at_site"] == "no").sum()), "n": len(agree)}
    p1 = out["coder1_precision_recall_same_rows"]["lex_wetside"]
    p2 = res["broad"]["lex_wetside"]
    out["rule_depends_on_coder"] = bool(abs(p1["precision"] - p2["precision"]) >= 0.10
                                        or abs(p1["recall"] - p2["recall"]) >= 0.10)
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out, indent=1, default=float))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.code:
        code_all(a.workers)
    if a.score:
        score()


if __name__ == "__main__":
    main()
