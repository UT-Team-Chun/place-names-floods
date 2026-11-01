# %% v82（事後・amendment B-17 の追記 2）: 第 13 回の報告の残りの求めに応じた記述の確かめ
"""確証判定は変えない。追記 2 で計算の前に固定した量。

- L: 地域構成の漏れを直に見る（文字の模型の点の AOI 平均と、自分・訓練側の兄弟の浸水の割合の Spearman）
- S: GFM の画像の枚数と、GFM の浸水の割合・EMS の浸水を GFM が拾った割合の Spearman
- K: 上限を掛けた後に浸水（乾き）のセルが 0 の AOI の数
- H: fa の ab・aab の照合（語・接尾辞、国ごと、パキスタンの語）
- E: ぼかした語源の数と、それを除いた precision・recall
- C: 言語の被覆の二つの量と「別の言語の地域」の定義（既にある計算を並べる）
- G: 場所グループで早期停止した基礎モデル（B-13 の C3、既にある計算を並べる）

出力: stdout と docs/verification/v82_b17_addendum2.json（部分ごとに足す。`--only L,S` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR

DOCV = ROOT / "docs" / "verification"
OUT = DOCV / "v82_b17_addendum2.json"
SEED = 20261106
N_BOOT = 1000
_MODS: dict = {}


def _load(name: str):
    if name in _MODS:
        return _MODS[name]
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    _MODS[name] = mod
    return mod


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr

    return float(spearmanr(a, b).statistic) if len(a) > 2 else float("nan")


def _boot_spearman(df: pd.DataFrame, x: str, y: str, rng: np.random.Generator) -> dict:
    """activation の再標本（復元抽出）で Spearman の 95% 区間。"""
    acts = df["activation"].unique()
    groups = {a: df.index[df["activation"] == a].to_numpy() for a in acts}
    vals = []
    for _ in range(N_BOOT):
        pick = rng.choice(acts, size=len(acts), replace=True)
        idx = np.concatenate([groups[a] for a in pick])
        s = df.loc[idx]
        vals.append(_spearman(s[x].to_numpy(), s[y].to_numpy()))
    v = np.array(vals)
    v = v[np.isfinite(v)]
    return {"rho": _spearman(df[x].to_numpy(), df[y].to_numpy()), "n_aoi": len(df), "n_activations": len(acts),
            "ci95": np.quantile(v, [0.025, 0.975]).tolist() if len(v) else None}


# %% L: 地域構成の漏れ
def l_leakage() -> dict:
    v77 = _load("v77_b16_local")
    folds = v77._dump_folds()
    if not folds:
        return {"note": "B1 dump of B-16 not found"}
    act = v77._superclusters()["aoi_activation"]
    rows = []
    for k, f in enumerate(folds):
        col = list(map(str, f["comp_cols"])).index("c_nlearn")
        te = pd.DataFrame({"aoi": f["aoi_te"].astype(str), "y": f["y_te"].astype(np.float64),
                           "score": f["xc_te"][:, col].astype(np.float64)})
        tr = pd.DataFrame({"aoi": f["aoi_tr"].astype(str), "y": f["y_tr"].astype(np.float64)})
        tr["activation"] = tr["aoi"].map(act)
        a_te = te.groupby("aoi").agg(score=("score", "mean"), share=("y", "mean"), n=("y", "size"))
        tr_act = tr.groupby("activation").agg(sib_share=("y", "mean"), sib_cells=("y", "size"),
                                              sib_aois=("aoi", "nunique"))
        for aoi, r in a_te.iterrows():
            a = act.get(aoi)
            sib = tr_act.loc[a] if a in tr_act.index else None
            rows.append({"aoi": aoi, "fold": k, "activation": a, "score": r["score"], "share": r["share"],
                         "n_cells": int(r["n"]), "has_sib": sib is not None,
                         "sib_share": float(sib["sib_share"]) if sib is not None else np.nan,
                         "sib_aois": int(sib["sib_aois"]) if sib is not None else 0})
    d = pd.DataFrame(rows).reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    w, wo = d[d["has_sib"]].reset_index(drop=True), d[~d["has_sib"]].reset_index(drop=True)
    out = {"n_aoi": len(d), "n_with_training_side_siblings": len(w), "n_without": len(wo),
           "score_vs_own_share_with_siblings": _boot_spearman(w, "score", "share", rng),
           "score_vs_own_share_without_siblings": _boot_spearman(wo, "score", "share", rng),
           "score_vs_sibling_share": _boot_spearman(w, "score", "sib_share", rng),
           "own_vs_sibling_share": _boot_spearman(w, "share", "sib_share", rng),
           "constant_score_within_aoi": bool(all(
               np.allclose(f["xc_te"][f["aoi_te"].astype(str) == a][:, list(map(str, f["comp_cols"])).index("c_nlearn")],
                           f["xc_te"][f["aoi_te"].astype(str) == a][0, list(map(str, f["comp_cols"])).index("c_nlearn")])
               for f in folds[:1] for a in np.unique(f["aoi_te"].astype(str))[:20]))}
    out["rule"] = ("B-17 addendum 2, item 1: held-out AOIs of the B1 dump of B-16 (primary split, frozen reading); score = AOI "
                   "mean of the N-learn column of M1c's regional composition; training-side siblings = other analysed AOIs "
                   "of the same activation on the training side of the fold (cell-weighted flooded share); Spearman across "
                   "AOIs, 95% interval from 1,000 resamples of activations (seed 20261106)")
    return out


# %% S: GFM の画像の枚数
def s_scenes() -> dict:
    v80 = _load("v80_b17_local")
    v36 = _load("v77_b16_local").v36
    gs = _json(DATA_DIR / "global" / "gfm_summary.json")["aois"]
    inc = {a: r for a, r in gs.items() if r.get("include")}
    sc = pd.DataFrame({"aoi": list(inc), "n_scenes": [inc[a].get("n_scenes") for a in inc],
                       "wet_frac": [inc[a].get("wet_frac") for a in inc]}).dropna()
    e = v36.load("world_v2", "primary")
    gf = v36.load("world_v2_gfm", "primary")
    ke, kg = v80._key(e["aoi"], e["lon"], e["lat"]), v80._key(gf["aoi"], gf["lon"], gf["lat"])
    pos = pd.Series(np.arange(len(ke)), index=ke).reindex(kg).to_numpy()
    ok = ~np.isnan(pos)
    ye = np.asarray(e["y"])[pos[ok].astype(int)]
    yg = np.asarray(gf["y"])[ok]
    j = pd.DataFrame({"aoi": np.asarray(gf["aoi"]).astype(str)[ok], "ems": ye, "gfm": yg})
    hit = j[j["ems"] == 1].groupby("aoi")["gfm"].mean().rename("hit")
    sc = sc.merge(hit, left_on="aoi", right_index=True, how="left")
    act = _load("v77_b16_local")._superclusters()["aoi_activation"]
    sc["activation"] = sc["aoi"].map(act)
    sc = sc.reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    q = sc["n_scenes"].quantile([0, 0.25, 0.5, 0.75, 1]).tolist()
    sh = sc.dropna(subset=["hit"]).reset_index(drop=True)
    return {"n_aoi": len(sc), "n_scenes_quantiles_0_25_50_75_100": q,
            "share_one_scene": float((sc["n_scenes"] == 1).mean()),
            "n_cells_joined": int(ok.sum()), "n_cells_gfm": len(kg),
            "scenes_vs_gfm_flooded_share": _boot_spearman(sc, "n_scenes", "wet_frac", rng),
            "n_aoi_with_ems_flooded_on_gfm_cells": len(sh),
            "scenes_vs_hit_rate": _boot_spearman(sh, "n_scenes", "hit", rng),
            "hit_rate_pooled": float(j.loc[j["ems"] == 1, "gfm"].mean()),
            "rule": "B-17 addendum 2, item 2: n_scenes and wet_frac of gfm_summary.json (192 AOIs); hit rate = share of "
                    "cells mapped flooded by EMS that GFM labels flooded (capped cells assessable on both); Spearman "
                    "across AOIs, 95% interval from 1,000 resamples of activations"}


# %% K: 上限で浸水のセルを失う AOI
def k_cap() -> dict:
    t = pd.read_parquet(DATA_DIR / "verification" / "v41_table.parquet", columns=["aoi", "y"])
    g = t.groupby("aoi")["y"].agg(["sum", "size"])
    return {"n_aoi": len(g), "n_no_flooded_after_cap": int((g["sum"] == 0).sum()),
            "n_no_dry_after_cap": int((g["sum"] == g["size"]).sum()),
            "n_both_classes": int(((g["sum"] > 0) & (g["sum"] < g["size"])).sum()),
            "aois_no_flooded": sorted(g.index[g["sum"] == 0].astype(str).tolist()),
            "rule": "B-17 addendum 2, item 3: capped analysed table (v41_table.parquet, EMS labels)"}


# %% H: ab の照合
def h_ab() -> dict:
    v73 = _load("v73_b14_local")
    lex = v73.Lexicon(v73.LEX)
    u = v73._analysed()
    gn = v73._gn()
    s = gn.reindex(np.unique(u["P_gid1"].to_numpy()))
    s = s[s["asciiname"].notna()]
    forms = {"ab", "aab"}
    by_pos: Counter = Counter()
    by_cc: Counter = Counter()
    by_cc_pos: Counter = Counter()
    pk_words: Counter = Counter()
    n = 0
    for r in s.itertuples():
        es = [e for e in v73._matches(lex, r.asciiname, r.cc) if e.form in forms and e.lang == "fa"]
        if not es:
            continue
        n += 1
        poss = sorted({e.pos for e in es})
        key = "+".join(poss)
        by_pos[key] += 1
        by_cc[str(r.cc)] += 1
        by_cc_pos[f"{r.cc}|{key}"] += 1
        if str(r.cc) == "PK":
            toks = [x for x in v73._TOK.split(str(r.asciiname).lower()) if x]
            for tk in toks:
                if any(v73._hit(tk, e) for e in es):
                    pk_words[tk] += 1
    return {"n_settlements_matched": n, "by_position": dict(by_pos), "by_country": dict(by_cc.most_common()),
            "by_country_position": dict(by_cc_pos.most_common()), "pk_words_top30": pk_words.most_common(30),
            "n_settlements_examined": len(s),
            "rule": "B-17 addendum 2, item 4: nearest settlements of the analysed cells; entries ab and aab of block fa "
                    "of the frozen lexicon (YAML 1.2 reading, as v73); a settlement counts once per position class"}


# %% E: ぼかした語源
HEDGES = ("perhaps", "probably", "possibly", "may", "might", "likely", "uncertain", "unclear", "disputed", "believed",
          "said to", "supposedly", "legend", "legendary", "folk etymology", "tradition", "traditionally", "according to")
_HEDGE_RX = re.compile(r"\b(?:" + "|".join(re.escape(h) for h in HEDGES) + r")\b", re.IGNORECASE)


def e_hedges() -> dict:
    v44 = _load("v44_wikipedia_etymology")
    v57b = _load("v57b_second_coder")
    samp = _json(v44.SAMPLE)
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    c1 = pd.read_csv(v44.CODES)
    c1["hedged"] = c1["quote"].fillna("").map(lambda q: bool(_HEDGE_RX.search(str(q))))
    c1["cats1"] = c1["cats"].fillna("").map(lambda s: {c for c in s.split(";") if c})
    c2 = pd.read_csv(v57b.CODES2)
    c2["cats2"] = c2["categories"].fillna("").map(lambda s: {c for c in str(s).split(";") if c})
    d = c1.merge(c2[["id", "status", "cats2"]], on="id", how="left", suffixes=("_1", "_2"))
    d["gid"] = d["id"].map(key)
    d = d.join(meta, on="gid")
    wet = set(v44.WET_SIDE)
    out: dict = {"terms": list(HEDGES), "n_quotes": len(d),
                 "n_hedged_all": int(d["hedged"].sum()),
                 "n_hedged_coded_coder1": int((d["hedged"] & (d["status_1"] == "coded")).sum()),
                 "n_coded_coder1": int((d["status_1"] == "coded").sum())}
    rng = np.random.default_rng(SEED)
    for who, status, cats in (("coder1", "status_1", "cats1"), ("coder2", "status_2", "cats2")):
        base = d[d[status] == "coded"].copy()
        base["truth"] = base[cats].map(lambda s: bool(set(s) & wet) if isinstance(s, set) else False)
        res = {}
        for lab, x in (("all", base), ("without_hedged", base[~base["hedged"]])):
            point = v44.weighted(x, sizes, "lex_wetside")
            parts = {s: x[x["stratum"] == s] for s in sizes}
            boots = []
            for _ in range(N_BOOT * 2):
                bb = pd.concat([p.iloc[rng.integers(0, len(p), len(p))] for p in parts.values() if len(p)])
                boots.append(v44.weighted(bb, sizes, "lex_wetside"))
            res[lab] = {"n": len(x), "n_hedged": int(x["hedged"].sum()),
                        **{k: {"est": point[k], "ci95": [float(np.nanpercentile([b[k] for b in boots], 2.5)),
                                                         float(np.nanpercentile([b[k] for b in boots], 97.5))]}
                           for k in ("precision", "recall")}}
        out[who] = res
    out["rule"] = ("B-17 addendum 2, item 5: hedged = coder 1's quote (v44_etym_codes.csv) matches any term at word "
                   "boundaries, case-insensitive; precision and recall of the lexicon's water-related flag with the "
                   "stratum weights of v44, 2,000 stratified bootstrap resamples (seed 20261106)")
    return out


# %% C: 言語の被覆の二つの量
def c_coverage() -> dict:
    v69 = _json(DOCV / "v69_b12_local.json").get("B3") or {}
    sh = (_json(DOCV / "v77_b16_local.json").get("B9c") or {}).get("share_cells") or {}
    v66 = (_json(DOCV / "v66_b10_local.json").get("B14") or {}).get("regional") or {}
    v66mod = _load("v66_b10_local")
    return {"first_language_outside_block_groups": v69.get("country_level_share"),
            "first_language_outside_block_groups_countries": {c: r["n_cells"] for c, r in
                                                              (v69.get("country_level") or {}).items()},
            "block_reads_first_language": sh.get("block_reads_first_language"),
            "first_language_without_own_block": (1 - sh["block_reads_first_language"])
            if sh.get("block_reads_first_language") is not None else None,
            "regions_of_another_language_share": v69.get("regional_share_v66"),
            "regions": {k: {"n_cells": r["n_cells"], "n_aoi": r["n_aoi"]} for k, r in v66.items()},
            "region_definitions": {k: list(v[:2]) + [sorted(v[2]) if v[2] else None]
                                   for k, v in v66mod.REGIONAL.items()},
            "source": "Natural Earth 1:10m admin-1 states and provinces (iso_a2 and name, region or geonunit)",
            "rule": "B-17 addendum 2, item 6: B3 of v69 (block language groups of LANG_COVER), B9 (c) of v77 (first "
                    "GeoNames language among the country's block names, Norway excluded in the frozen reading), B14 of "
                    "v66 (13 regions)"}


# %% G: 場所グループの早期停止（B-13 の C3）
def g_grouped_es() -> dict:
    c = _json(DOCV / "v71_b13_cluster.json").get("C3_grouped_es") or {}
    nm = c.get("names") or {}
    return {"run": c.get("run"), "yardstick": c.get("yardstick"), "ci95_yardstick": c.get("ci95_yardstick"),
            "R": nm.get("R"), "ci90_R": nm.get("ci90_R"), "ci95_R": nm.get("ci95_R"), "decision": nm.get("decision"),
            "R_minus_Rp": (nm.get("placebo") or {}).get("R_minus_Rp"),
            "ci95_R_minus_Rp": (nm.get("placebo") or {}).get("ci95_R_minus_Rp"),
            "brier": c.get("brier"), "brier_training_constant": c.get("brier_training_constant"),
            "iterations": c.get("grouped_es_iterations"),
            "rule": "B-17 addendum 2, item 7: C3 of B-13 (v71), copied"}


PARTS = {"L": l_leakage, "S": s_scenes, "K": k_cap, "H": h_ab, "E": e_hedges, "C": c_coverage, "G": g_grouped_es}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    want = [p for p in (a.only.split(",") if a.only else PARTS) if p]
    for k in want:
        t0 = time.time()
        r = PARTS[k]()
        r["seconds"] = round(time.time() - t0)
        res = _json(OUT) if OUT.is_file() else {}
        res[k] = r
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        print(f"{k} ({r['seconds']} s) {json.dumps(r, ensure_ascii=False, default=float)[:2000]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
