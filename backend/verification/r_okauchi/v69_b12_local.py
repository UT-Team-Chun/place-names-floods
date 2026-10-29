# %% v69（事後・amendment B-12 の B1〜B21）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-12 で計算の前に固定した量を、既存の表と held-out の予測から計算する。

- B1: 判定の規則の動作特性（固定した余白の規則と比の規則、入れ子の場所グループの bootstrap）
- B2: 上限なしの集落の近くの標本の構成（国・activation）と、パキスタン・バングラデシュを除いた OR
- B3: 言語の割り当ての数え直し（ブロックの言語の一覧で、国の第一言語を読むか）
- B4: 包含の 5 つの基準を順に当てたときの除外の数
- B5: 帰無の模擬の外れた 13 回の、ΔBS に最も寄与した場所グループ
- B6: BH で残る推定（高地、GFM）と、高地の名の近くの GFM の評価域
- B7r: 再起の記録（227 AOI）を全体の尺度に
- B8: 表 2 の柔軟な「辞書だけ」「学習した得点だけ」の固定した余白の区間
- B9: 二人目の符号化者の precision・recall の区間
- B10: 主分割の基準のモデル（文字モデルの地域構成を除いた M1c）と −0.74 の区間
- B11: 名前の係数（と地域構成）だけに罰則をかけ、格子を内側の最適まで広げた二段目
- B12: GFM の時間窓（画像の時刻と EMS の画像の日の差）
- B13: 浸水の種別（内陸の規則、種別のある非河川の浸水、flood trace の割合）
- B14: 浸水とも乾きとも決まらず除いたセル
- B15: HAND だけの地図の CSI（同点は無作為の順）
- B16: JRC の図そのものの CSI（手元に特徴があれば）
- B17: 各分割の M0 と定数の予測の Brier
- B18: 市街地の 10 m の割合（ESA WorldCover）を足した距離帯の対比と水系側の OR
- B19: 集落の近くの水系側の OR の全ての版（目録から）
- B20: p0（読まなかった名のうち語源が水に関わる割合）
- B21: 照合の標本の設計の改め（読まなかった名の層を足す。答えは別のファイル、公開しない）

出力: stdout と docs/verification/v69_b12_local.json（部分ごとに足す。`--only B1,B4` で一部だけ）。
"""

from __future__ import annotations

import csv
import datetime as dt
import importlib.util
import json
import math
import os
import re
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import stage2_weighted
from topsem.world.b10 import C_GRID
from topsem.world.dataset import evaluable, riverine_labels
from topsem.world.ems import flood_class
from topsem.world.lexicon import WET_SIDE

OUT = ROOT / "docs" / "verification" / "v69_b12_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
LABELS = DATA_DIR / "global" / "labels"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
COUNTRY_INFO = DATA_DIR / "global" / "geonames" / "countryInfo.txt"
INCIDENCE = DATA_DIR / "global" / "incidence.json"
GFM_SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
SEED = 20261101
SESOI_R = 0.05
B_OUT, B_IN = 2000, 1000
B_BOOT = 4000
B_PAIRS = 2000
NEAR_M = 500.0
N_WORKERS = int(os.environ.get("V69_WORKERS", "5"))
K_LOGIT = 1000.0  # B11: 基礎の logit の列をこの倍にして、罰則が名前と地域構成の列だけにかかるようにする
AUDIT_NOT_N = 100  # B21: 読まなかった名の層の名の数


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v36 = _load("v36_yardstick_robustness")
v39 = _load("v39_within_aoi_odds")
v55 = _load("v55_b7_local")
v68 = _load("v68_b11_local")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _group_sums(d: dict, arms: tuple[str, ...]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    y = d["y"].astype(np.float64)
    units, uc = np.unique(np.asarray(d["unit"]).astype(str), return_inverse=True)
    g = len(units)
    n = np.bincount(uc, minlength=g).astype(np.float64)
    s = {a: np.bincount(uc, weights=(d[f"p_{a}"].astype(np.float64) - y) ** 2, minlength=g) for a in arms}
    return n, s, units


# %% B1: 判定の規則の動作特性
R_GRID = (-0.06, -0.05, -0.04, -0.02, 0.0, 0.02, 0.04, 0.05, 0.06)


def b1_operating() -> dict:
    """真の R を r0 に置いた母集団（確証の実行の予測を場所グループごとに直したもの）から、入れ子の bootstrap で
    二つの規則が同等（90% 区間が ±0.05 の内）と言う割合。母集団の直し方は二つ:

    - proportional（B-12 の B1 の文）: グループごとの名前の増分を eg − Ê n + r0 dg にする（重みを付けた標本では
      E_w − Ê + r0 D_w。名前の効果が各グループの基礎の利得に比例する）
    - constant（実装で並べた別の形）: eg − Ê n + r0 D̂ n（名前の効果がセルに一様、E_w − Ê + r0 D̂）

    外の回は「新しい標本」で、その標本の点の推定と、その標本の中の内の bootstrap の 90% 区間で規則を当てる。分母の
    95% 区間が 0 を含む回は同等と言わない（判定表の undefined-denominator）。
    """
    d = v36.load("world_v2", "primary")
    n, s, units = _group_sums(d, ("M0", "M1", "M1c", "M3n"))
    g = len(units)
    dg, eg = s["M0"] - s["M1"], s["M1c"] - s["M3n"]
    d_hat, e_hat = float(dg.sum() / n.sum()), float(eg.sum() / n.sum())
    rng = np.random.default_rng(SEED)
    cnt = {sh: {r0: {"fixed": 0, "ratio": 0, "undefined": 0} for r0 in R_GRID} for sh in ("proportional", "constant")}
    outer_e, outer_d = [], []
    t0 = time.time()
    for _ in range(B_OUT):
        w = rng.multinomial(g, np.full(g, 1.0 / g)).astype(np.float64)
        nw = w @ n
        d_b, e_b = (w @ dg) / nw, (w @ eg) / nw
        outer_e.append(e_b)
        outer_d.append(d_b)
        v = rng.multinomial(g, w / g, size=B_IN).astype(np.float64)
        nv = v @ n
        dv, ev = (v @ dg) / nv, (v @ eg) / nv
        d_lo = float(np.quantile(dv, 0.025))
        for r0 in R_GRID:
            for sh in ("proportional", "constant"):
                if sh == "proportional":
                    e_out, e_in = e_b - e_hat + r0 * d_b, ev - e_hat + r0 * dv
                else:
                    e_out, e_in = e_b - e_hat + r0 * d_hat, ev - e_hat + r0 * d_hat
                c = cnt[sh][r0]
                if d_lo <= 0:
                    c["undefined"] += 1
                    continue
                lo, hi = np.quantile(e_in, [0.05, 0.95]) / d_b
                if -SESOI_R < lo and hi < SESOI_R:
                    c["fixed"] += 1
                rlo, rhi = np.quantile(e_in / dv, [0.05, 0.95])
                if -SESOI_R < rlo and rhi < SESOI_R:
                    c["ratio"] += 1
                del e_out
    oe, od = np.array(outer_e), np.array(outer_d)
    res = {sh: {f"{r0:+.2f}": {k: v / B_OUT for k, v in c.items()} for r0, c in cs.items()} for sh, cs in cnt.items()}
    return {"n_units": g, "replicates_outer": B_OUT, "replicates_inner": B_IN, "seed": SEED,
            "delta_bs_hat": e_hat, "yardstick_hat": d_hat, "R_hat": e_hat / d_hat,
            "se_delta_bs": float(oe.std(ddof=1)), "se_yardstick": float(od.std(ddof=1)),
            "corr_delta_yardstick": float(np.corrcoef(oe, od)[0, 1]),
            "share_equivalent": res, "seconds": round(time.time() - t0),
            "rule": "B-12 B1: share of outer replicates (new samples) in which each rule declares equivalence "
                    "(90% interval inside ±0.05; denominator's 95% interval above zero) with the true R set to r0"}


# %% B2: 上限なしの標本の構成と、パキスタン・バングラデシュを除いた OR
def _activation_countries() -> dict[str, list[str]]:
    out = {}
    for r in _json(INCIDENCE)["aois"]:
        out[r["tag"]] = list(r.get("countries") or [])
    return out


def b2_composition() -> dict:
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    rate = u.groupby("aoi")["y"].transform("mean").to_numpy()
    t = u.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)  # odds() と同じ AOI の採否
    gn = pd.read_parquet(GN, columns=["gid", "cc"]).drop_duplicates("gid").set_index("gid")["cc"]
    t["cc"] = gn.reindex(t["P_gid1"].astype(np.int64).to_numpy()).fillna("").to_numpy(str)
    t["activation"] = t["aoi"].str.split("_").str[0]
    near = (t["lex_wetside_n1"].to_numpy() > 0) & (t["P_d1"].to_numpy() <= NEAR_M)
    fl = near & (t["y"].to_numpy() == 1)
    t["near"], t["flooded_near"] = near, fl
    tot = int(fl.sum())

    def summarise(key: str, top: int = 15) -> list[dict]:
        gg = t.groupby(key).agg(cells_near=("near", "sum"), flooded_near=("flooded_near", "sum"),
                                aois=("aoi", "nunique")).reset_index()
        gg = gg.loc[gg["cells_near"] > 0].sort_values(["flooded_near", "cells_near"], ascending=False)
        return [{key: str(r[key]), "flooded_near": int(r["flooded_near"]), "cells_near": int(r["cells_near"]),
                 "share_of_flooded_near": float(r["flooded_near"] / tot) if tot else 0.0, "aois": int(r["aois"])}
                for _, r in gg.head(top).iterrows()]

    by_cc, by_act = summarise("cc"), summarise("activation")
    out: dict = {"flooded_near_total": tot, "cells_near_total": int(near.sum()),
                 "n_countries_with_flooded_near": int(t.loc[fl, "cc"].nunique()),
                 "n_activations_with_flooded_near": int(t.loc[fl, "activation"].nunique()),
                 "share_top5_countries": float(sum(r["flooded_near"] for r in by_cc[:5]) / tot) if tot else None,
                 "share_pk_bd": float(sum(r["flooded_near"] for r in by_cc if r["cc"] in {"PK", "BD"}) / tot)
                 if tot else None,
                 "by_country": by_cc, "by_activation": by_act}
    ac = _activation_countries()
    pkbd = u["aoi"].map(lambda a: any(c in ac.get(a, []) for c in ("Pakistan", "Bangladesh"))).to_numpy()
    v = u.loc[~pkbd].reset_index(drop=True)
    out["n_aoi_pk_bd"] = int(u.loc[pkbd, "aoi"].nunique())
    out["without_pk_bd"] = v68.pairs_or(v, "b12_nopkbd_uncapped")
    out["rule"] = ("B-12 B2: uncapped sample of B-11 B2 (corrected readings); countries of the nearest settlement; "
                   "Pakistan and Bangladesh removed by the activation's countries")
    return out


# %% B3: 言語の割り当ての数え直し
LANG_COVER = {
    "sh": {"sh", "sr", "hr", "bs", "cnr"}, "bg": {"bg", "mk"}, "cs": {"cs", "sk"}, "uk": {"uk", "ru", "be"},
    "no": {"no", "nb", "nn", "is"}, "fa": {"fa", "ur", "ps", "sd"}, "hi": {"hi", "bn", "ne"}, "id": {"id", "ms"},
    "tl": {"tl", "fil"}, "qu": {"qu", "ay"}, "th": {"th", "lo"}, "kk": {"kk", "ky", "uz", "tk", "mn"},
    "am": {"am", "ti"}, "ng": {"ha", "yo", "ig"}, "ro": {"ro", "mo"},
    "bantu": {"ny", "bem", "sn", "nd", "zu", "xh", "rw", "rn", "lg", "ln", "kg", "ts", "tn", "st", "ss", "ve", "nr",
              "ki", "yao", "tum", "loz", "ng", "hz", "sw", "mg"},
}


def b3_language() -> dict:
    lex = v68.lexicon()
    first = v36.first_languages(COUNTRY_INFO)
    d = v36.load("world_v2", "primary")
    cc = np.asarray(d["cc"]).astype(str)
    n_all = len(cc)
    rows, rows_old = {}, {}
    old = _load("v66_b10_local").LANG_GROUPS
    for c in sorted(set(cc) - {""}):
        langs = set(lex.country_langs.get(c, ()))
        fl = first.get(c, "")
        cover = set().union(*(LANG_COVER.get(x, {x}) for x in langs)) if langs else set()
        cover_old = set().union(*(old.get(x, {x}) for x in langs)) if langs else set()
        m = int((cc == c).sum())
        if fl and fl not in cover:
            rows[c] = {"first_language": fl, "lexicon_blocks": sorted(langs), "n_cells": m}
        if fl and fl not in cover_old:
            rows_old[c] = {"first_language": fl, "lexicon_blocks": sorted(langs), "n_cells": m}
    v66 = _json(DOCV / "v66_b10_local.json").get("B14") or {}
    reg = v66.get("regional") or {}
    reg_cells = int(sum(r["n_cells"] for r in reg.values()))
    n_country = int(sum(r["n_cells"] for r in rows.values()))
    n_old = int(sum(r["n_cells"] for r in rows_old.values()))
    return {"n_cells": n_all, "country_level": rows, "country_level_share": n_country / n_all,
            "country_level_share_v66": n_old / n_all, "removed_by_fix": sorted(set(rows_old) - set(rows)),
            "regional_share_v66": reg_cells / n_all, "regional": {k: v["n_cells"] for k, v in reg.items()},
            "rule": "B-12 B3: a country is flagged when the first GeoNames language of the nearest settlement's country "
                    "is not among the languages of the blocks applied to it (block languages from the block names); "
                    "regional flags of v66 B14 unchanged and counted separately"}


# %% B4: 包含の順の数
def b4_inclusion() -> dict:
    rows = _json(INCIDENCE)["aois"]
    steps = [("evaluable cells >= 400", lambda r: r["n_eval"] >= 400),
             ("flooded cells >= 20", lambda r: r["n_wet"] >= 20),
             ("flooded share <= 80 %", lambda r: r["n_wet"] <= 0.8 * r["n_eval"]),
             ("riverine (or untyped inland) share >= 0.5",
              lambda r: (r["riverine_share"] + (r["unspecified_share"] if r.get("unspecified_as_riverine") else 0))
              >= 0.5),
             ("a settlement inside the grid", lambda r: r.get("n_P_inside", 0) >= 1)]
    left = list(rows)
    out = []
    for name, ok in steps:
        keep = [r for r in left if ok(r)]
        out.append({"criterion": name, "excluded": len(left) - len(keep), "remaining": len(keep)})
        left = keep
    inc = sum(1 for r in rows if r.get("included"))
    return {"n_labelled": len(rows), "steps": out, "n_included_recorded": inc,
            "matches_recorded": len(left) == inc,
            "rule": "B-12 B4: the five inclusion criteria applied in the order of the protocol"}


# %% B5: 帰無の模擬の外れた回
def b5_null_outliers() -> dict:
    v55j = _json(DOCV / "v55_b7_local.json")["null_outliers"]
    sim = RUNS_DIR / "power_v3_null" / "sim"
    prim = v36.load("world_v2", "primary")
    pu = np.asarray(prim["unit"]).astype(str)
    pcc = np.asarray(prim["cc"]).astype(str)
    unit_cc = pd.Series(pcc).groupby(pu).agg(lambda s: s.value_counts().index[0]).to_dict()
    out = []
    for o in v55j["outlying"]:
        z = np.load(sim / f"{o['rep']}.npz", allow_pickle=False)
        y = z["y"].astype(np.float64)
        units, uc = np.unique(z["unit"].astype(str), return_inverse=True)
        g = len(units)
        n = np.bincount(uc, minlength=g).astype(np.float64)
        e1 = np.bincount(uc, weights=(z["M1c"].astype(np.float64) - y) ** 2, minlength=g)
        e3 = np.bincount(uc, weights=(z["M3n"].astype(np.float64) - y) ** 2, minlength=g)
        d0 = np.bincount(uc, weights=(z["M0"].astype(np.float64) - y) ** 2, minlength=g)
        d1 = np.bincount(uc, weights=(z["M1"].astype(np.float64) - y) ** 2, minlength=g)
        yard = (d0.sum() - d1.sum()) / n.sum()
        contrib = (e1 - e3) / n.sum() / yard
        i = int(np.argmax(np.abs(contrib)))
        ysum = np.bincount(uc, weights=y, minlength=g)
        out.append({"rep": o["rep"], "R": float(contrib.sum()), "top_unit": str(units[i]),
                    "top_unit_R": float(contrib[i]), "top_unit_share_of_R": float(contrib[i] / contrib.sum())
                    if contrib.sum() != 0 else None,
                    "top_unit_country": unit_cc.get(str(units[i]), ""),
                    "top_unit_flooded_share": float(ysum[i] / n[i]), "top_unit_cells": int(n[i])})
    return {"outlying": out, "n": len(out),
            "rule": "B-12 B5: the location group with the largest absolute contribution to R in each outlying null "
                    "replicate (v55), with its country (majority of cells) and flooded share"}


# %% B6: 高地の区分の推定と GFM の評価域
def b6_highland() -> dict:
    v41 = _json(DOCV / "v41_gfm_independent_labels.json")
    sens = v41.get("sensitivity_independent_sensor") or {}

    def term(lab: str) -> dict | None:
        x = (((sens.get(f"odds_{lab}") or {}).get("categories_adjusted_M1") or {}).get("terms") or {}).get(
            "high:le500m")
        return x

    t = v55.table()
    summ = _json(GFM_SUMMARY)
    keep = sorted(a for a, x in summ["aois"].items() if x.get("include"))
    g = pd.concat([pd.read_parquet(LABELS / a / "gfm.parquet", columns=["ix", "iy", "gfm_eval"]).assign(aoi=a)
                   for a in keep], ignore_index=True)
    m = t.loc[t["aoi"].isin(keep)].merge(g, on=["aoi", "ix", "iy"], how="left", validate="many_to_one")
    ev = m["gfm_eval"].astype("boolean").fillna(value=False).to_numpy(dtype=bool)
    near = m["P_d1"].to_numpy() <= NEAR_M
    hi = m["lex_high_n1"].to_numpy() > 0
    ws = m["lex_wetside_n1"].to_numpy() > 0
    out: dict = {"gfm_not_s1": term("gfm"), "ems_not_s1": term("ems"), "n_aoi_not_s1": sens.get("n_aoi"),
                 "n_gfm_aoi": len(keep)}
    out["share_assessable_le500m"] = {
        "highland": float(ev[near & hi].mean()) if (near & hi).any() else None,
        "water_side": float(ev[near & ws].mean()) if (near & ws).any() else None,
        "other": float(ev[near & ~hi & ~ws].mean()) if (near & ~hi & ~ws).any() else None,
        "n_cells": {"highland": int((near & hi).sum()), "water_side": int((near & ws).sum()),
                    "other": int((near & ~hi & ~ws).sum())}}
    if "hand_min" in m:
        hand = m["hand_min"].to_numpy(np.float64)
        out["median_hand_min_le500m"] = {
            "highland": float(np.nanmedian(hand[near & hi])), "water_side": float(np.nanmedian(hand[near & ws])),
            "other": float(np.nanmedian(hand[near & ~hi & ~ws])),
            "highland_assessable": float(np.nanmedian(hand[near & hi & ev])),
            "highland_not_assessable": float(np.nanmedian(hand[near & hi & ~ev]))}
    out["rule"] = "B-12 B6: the estimate surviving BH (v41) and GFM assessability within 500 m by class"
    return out


# %% B7r: 再起の記録を全体の尺度に
def b7r_recurrence() -> dict:
    r = _json(DOCV / "v56_b7_cluster.json")["recur"]
    s = r["subset_aois_with_other_events"]
    rs = s["recur_settlement"]
    share = s["n_cells"] / r["n_cells"]
    yard = r["yardstick"]
    f = share / yard
    return {"subset_R": rs["R"], "subset_yardstick": s["yardstick"], "subset_cells": s["n_cells"],
            "all_cells": r["n_cells"], "pooled_yardstick": yard,
            "pooled_scale_R": rs["delta_bs_M1cE_M3oE"] * f,
            "pooled_scale_ci90": [rs["ci90_delta"][0] * f, rs["ci90_delta"][1] * f],
            "rule": "B-12 B7': the subset's ΔBS times the subset's share of cells, divided by the baseline gain of "
                    "all cells of the same run"}


# %% B8: 表 2 の柔軟な辞書だけ・学習した得点だけ
def b8_flexible_rows() -> dict:
    d = v36.load("world_v2", "primary")
    n, s, units = _group_sums(d, ("M0", "M1", "M3x", "M3l"))
    g = len(units)
    yard = float((s["M0"].sum() - s["M1"].sum()) / n.sum())
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_BOOT).astype(np.float64)
    nw = w @ n
    out: dict = {"yardstick": yard}
    for arm in ("M3x", "M3l"):
        e = s["M1"] - s[arm]
        pt = float(e.sum() / n.sum())
        b = (w @ e) / nw
        lo, hi = np.quantile(b, [0.025, 0.975])
        out[arm] = {"delta_bs": pt, "R": pt / yard, "ci95_R_fixed": [float(lo / yard), float(hi / yard)],
                    "ci95_delta_bs": [float(lo), float(hi)]}
    out["rule"] = "B-12 B8: location-group bootstrap (4,000) of ΔBS(M1→arm), divided by the point baseline gain"
    return out


# %% B9: 二人目の符号化者の区間
def b9_coder2() -> dict:
    v44 = _load("v44_wikipedia_etymology")
    v57b = _load("v57b_second_coder")
    c2 = []
    for f in sorted(v57b.CACHE.glob("batch_*.json")):
        c2 += _json(f)["codes"]
    if not c2:
        return {"note": "coder-2 batches not found locally"}
    c2 = pd.DataFrame(c2)
    c2["id"] = c2["id"].astype(int)
    c2["cats2"] = c2["categories"].map(lambda v: set(v) if isinstance(v, list) else set())
    samp = _json(v44.SAMPLE)
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v57b.TO_CODE.read_text(encoding="utf-8").splitlines()}
    c1 = pd.read_csv(v44.CODES)
    d = c1.merge(c2, on="id", how="inner", suffixes=("_1", "_2"))
    d["gid"] = d["id"].map(key)
    d = d.join(meta, on="gid")
    d2 = d[d["status_2"] == "coded"].copy()
    d2["truth"] = d2["cats2"].map(lambda s: bool(s & set(WET_SIDE)))
    out: dict = {}
    for flag in ("lex_wetside", "nread_wetside"):
        point = v44.weighted(d2, sizes, flag)
        rng = np.random.default_rng(v44.SEED)
        parts = {s: d2[d2["stratum"] == s] for s in sizes}
        boots = []
        for _ in range(v44.N_BOOT):
            bb = pd.concat([p.iloc[rng.integers(0, len(p), len(p))] for p in parts.values()])
            boots.append(v44.weighted(bb, sizes, flag))
        out[flag] = {k: {"est": point[k], "ci95": [float(np.nanpercentile([b[k] for b in boots], 2.5)),
                                                    float(np.nanpercentile([b[k] for b in boots], 97.5))]}
                     for k in ("precision", "recall")}
    out["rule"] = "B-12 B9: the stratified bootstrap of v44 (same seed and replicates) on coder 2's codes"
    return out


# %% B10: 主分割の基準のモデルと −0.74 の区間
def _dump_folds() -> list[dict]:
    src = RUNS_DIR / "world_v2_b8stage2_local" / "predictions" / "primary"
    return [dict(np.load(src / f"fold{k}_stage2.npz", allow_pickle=False)) for k in range(10)]


def b10_reference() -> dict:
    folds = _dump_folds()
    sse = {k: 0.0 for k in ("M1s", "M1c", "M1c_no_nlearn", "M0", "M1")}
    n = 0
    for z in folds:
        comp = [str(c) for c in z["comp_cols"]]
        j = comp.index("c_nlearn")
        keep = [i for i in range(len(comp)) if i != j]
        y = z["y_te"].astype(np.float64)
        empty_tr, empty_te = np.zeros((len(z["z_tr"]), 0)), np.zeros((len(z["z_te"]), 0))
        p_s, _ = stage2_weighted(z["z_tr"], z["z_te"], empty_tr, empty_te, z["y_tr"])
        p_c, _ = stage2_weighted(z["z_tr"], z["z_te"], z["xc_tr"], z["xc_te"], z["y_tr"])
        p_cn, _ = stage2_weighted(z["z_tr"], z["z_te"], z["xc_tr"][:, keep], z["xc_te"][:, keep], z["y_tr"])
        for k_, p in (("M1s", p_s), ("M1c", p_c), ("M1c_no_nlearn", p_cn), ("M0", z["p_M0_te"]),
                      ("M1", z["p_M1_te"])):
            sse[k_] += float(((p.astype(np.float64) - y) ** 2).sum())
        n += len(y)
    bs = {k: v / n for k, v in sse.items()}
    yard = bs["M0"] - bs["M1"]
    v55j = _json(DOCV / "v55_b7_local.json")["residual_between_aois"]["same_activation"]["0-50km"]
    r, m = float(v55j["corr"]), int(v55j["n_pairs"])
    zf, se = math.atanh(r), 1 / math.sqrt(m - 3)
    return {"brier": bs, "yardstick": yard,
            "gain_M1c_over_M1s_share_of_yardstick": (bs["M1s"] - bs["M1c"]) / yard,
            "gain_M1c_no_nlearn_over_M1s_share_of_yardstick": (bs["M1s"] - bs["M1c_no_nlearn"]) / yard,
            "residual_corr": {"r": r, "n_pairs": m,
                              "ci95_fisher": [math.tanh(zf - 1.96 * se), math.tanh(zf + 1.96 * se)]},
            "rule": "B-12 B10: second stage refitted on the dumps of B-8 C2 (corrected code, primary split): M1s (no "
                    "composition), M1c, and M1c without the AOI mean of the character model's score"}


# %% B11: 名前の係数（と地域構成）だけに罰則
def _inner_folds(units: np.ndarray, k: int = 5) -> np.ndarray:
    u, inv = np.unique(units.astype(str), return_inverse=True)
    lab = np.array([zlib.crc32(f"B-12/B11/{x}".encode()) % k for x in u])
    return lab[inv]


def _choose_c_scaled(a: np.ndarray, y: np.ndarray, inner: np.ndarray) -> tuple[float, dict]:
    from topsem.world.b10 import choose_c

    grid = list(C_GRID)
    hist = {}
    for _ in range(4):  # 最適が格子の端なら一桁ずつ広げる（最大 4 回）
        c, bs = choose_c(a, y, inner, tuple(grid))
        hist.update(dict(zip(grid, bs, strict=True)))
        if c == grid[0]:
            grid = [grid[0] / 10, *grid]
        elif c == grid[-1]:
            grid = [*grid, grid[-1] * 10]
        else:
            return c, {"grid": grid, "inner_brier": [hist[x] for x in grid], "interior": True}
    return c, {"grid": grid, "inner_brier": [hist.get(x) for x in grid], "interior": False}


def _b11_fold(k: int) -> dict:
    _threads1()
    from sklearn.linear_model import LogisticRegression

    z = _dump_folds_one(k)
    inner = _inner_folds(z["unit_tr"])
    res = {"fold": k}
    for arm, x_tr, x_te in (("M1c", z["xc_tr"], z["xc_te"]),
                            ("M3n", np.column_stack([z["xc_tr"], z["xn_tr"]]),
                             np.column_stack([z["xc_te"], z["xn_te"]]))):
        a = np.column_stack([z["z_tr"], x_tr]).astype(np.float64)
        b = np.column_stack([z["z_te"], x_te]).astype(np.float64)
        mu, sd = a.mean(axis=0), a.std(axis=0)
        sd[sd == 0] = 1.0
        a, b = (a - mu) / sd, (b - mu) / sd
        a[:, 0] *= K_LOGIT
        b[:, 0] *= K_LOGIT
        c, info = _choose_c_scaled(a, z["y_tr"], inner)
        m = LogisticRegression(C=c, max_iter=2000)
        m.fit(a, z["y_tr"])
        p = m.predict_proba(b)[:, 1]
        res[arm] = {"C": c, **info, "sse": float(((p - z["y_te"]) ** 2).sum())}
    y = z["y_te"].astype(np.float64)
    res["sse_M0"] = float(((z["p_M0_te"] - y) ** 2).sum())
    res["sse_M1"] = float(((z["p_M1_te"] - y) ** 2).sum())
    res["n"] = len(y)
    return res


def _dump_folds_one(k: int) -> dict:
    src = RUNS_DIR / "world_v2_b8stage2_local" / "predictions" / "primary"
    return dict(np.load(src / f"fold{k}_stage2.npz", allow_pickle=False))


def b11_names_only_penalty() -> dict:
    _threads1()
    with ProcessPoolExecutor(N_WORKERS) as ex:
        folds = list(ex.map(_b11_fold, range(10)))
    n = sum(f["n"] for f in folds)
    bs = {k: sum(f[k]["sse"] for f in folds) / n for k in ("M1c", "M3n")}
    yard = (sum(f["sse_M0"] for f in folds) - sum(f["sse_M1"] for f in folds)) / n
    return {"R": (bs["M1c"] - bs["M3n"]) / yard, "brier": bs, "yardstick": yard,
            "C_chosen": {arm: [f[arm]["C"] for f in folds] for arm in ("M1c", "M3n")},
            "interior": {arm: [f[arm]["interior"] for f in folds] for arm in ("M1c", "M3n")},
            "folds": folds, "k_logit": K_LOGIT,
            "rule": "B-12 B11: base-logit column scaled by K so that the L2 penalty acts on composition and name "
                    "features only; C chosen per fold and arm by five inner folds of location groups (inner Brier), "
                    "grid widened until the optimum is interior"}


# %% B12: GFM の時間窓
def b12_gfm_window() -> dict:
    summ = _json(GFM_SUMMARY).get("aois") or {}
    used = sorted(a for a, x in summ.items() if (LABELS / a / "gfm.parquet").is_file() and x.get("n_scenes", 1))
    n_aoi, multi_day = 0, 0
    day_offsets: list[int] = []
    aoi_same_day, aoi_only_adjacent = 0, 0
    frames = []
    for a in used:
        p = LABELS / a / "gfm_meta.json"
        if not p.is_file():
            continue
        m = _json(p)
        if not m.get("window") or not m.get("scenes"):
            continue
        n_aoi += 1
        d0 = dt.date.fromisoformat(m["window"][0]) + dt.timedelta(days=1)
        d1 = dt.date.fromisoformat(m["window"][1]) - dt.timedelta(days=1)
        if d1 > d0:
            multi_day += 1
        ems_days = [d0 + dt.timedelta(days=i) for i in range((d1 - d0).days + 1)]
        days = sorted({dt.datetime.fromisoformat(x["datetime"].replace("Z", "+00:00")).date() for x in m["scenes"]})
        frames.append(len(m["scenes"]))
        offs = [min((dd - e).days for e in ems_days if (dd - e).days >= 0) if dd >= ems_days[0]
                else (dd - ems_days[0]).days for dd in days]
        day_offsets += offs
        if 0 in offs:
            aoi_same_day += 1
        else:
            aoi_only_adjacent += 1
    o = np.array(day_offsets)
    return {"n_aoi": n_aoi, "aoi_with_several_ems_days": multi_day,
            "acquisition_days_by_offset": {str(k): int((o == k).sum()) for k in sorted(set(o.tolist()))},
            "aoi_with_an_image_on_an_ems_day": aoi_same_day, "aoi_with_images_only_on_adjacent_days": aoi_only_adjacent,
            "frames_per_aoi_median": float(np.median(frames)) if frames else None,
            "rule": "B-12 B12: AOIs with GFM labels; window = from 00:00 of the day before the earliest EMS image date "
                    "to 24:00 of the day after the latest; distinct acquisition days of the GFM frames found in the "
                    "window, offset from the nearest EMS image day"}


# %% B13: 浸水の種別
def b13_typing() -> dict:
    rows = {r["tag"]: r for r in _json(INCIDENCE)["aois"] if r.get("included")}
    sens = {a: s.get("sensor", "none") for a, s in (_json(GFM_SUMMARY).get("aois") or {}).items()}
    typed_other, typed_flash, typed_coastal, untyped_any = 0, 0, 0, 0
    typed_rows: list[dict] = []
    other_desc: dict[str, int] = {}
    trace = {}
    coast = []
    for a, r in rows.items():
        m = _json(LABELS / a / "meta.json")
        classes = {}
        for desc, k in (m.get("obj_desc") or {}).items():
            cl = "untyped" if not desc else flood_class(desc)
            classes[cl] = classes.get(cl, 0) + int(k)
            if cl == "other" and desc:
                other_desc[desc] = other_desc.get(desc, 0) + int(k)
        typed = [dd for dd in (m.get("obj_desc") or {}) if dd and flood_class(dd) == "other"]
        if typed:
            typed_other += 1
            typed_rows.append({"aoi": a, "descriptions": typed, "inland_rule": bool(r.get("unspecified_as_riverine")),
                               "riverine_share": r.get("riverine_share"), "unspecified_share": r.get("unspecified_share")})
        typed_flash += int(classes.get("flash", 0) > 0)
        typed_coastal += int(classes.get("coastal", 0) > 0)
        untyped_any += int(classes.get("untyped", 0) > 0)
        c = pq.read_table(LABELS / a / "cells.parquet",
                          columns=["f_trace", "f_flood_any", "f_aoi", "f_observed", "f_perm_water"]).to_pandas()
        ev = evaluable(c)
        s = sens.get(a, "none")
        tr = trace.setdefault(s, [0.0, 0.0])
        tr[0] += float(c.loc[ev, "f_trace"].sum())
        tr[1] += float(c.loc[ev, "f_flood_any"].sum())
        coast.append(float(r.get("coast_km", np.nan)))
    co = np.array(coast)
    return {"n_aoi": len(rows),
            "riverine_by_inland_rule": int(sum(1 for r in rows.values() if r.get("unspecified_as_riverine"))),
            "aoi_with_typed_other": typed_other, "typed_other_aois": typed_rows, "aoi_with_flash": typed_flash, "aoi_with_coastal": typed_coastal,
            "aoi_with_untyped": untyped_any,
            "other_descriptions_top": sorted(other_desc.items(), key=lambda kv: -kv[1])[:15],
            "trace_share_by_sensor": {s: (v[0] / v[1] if v[1] else None) for s, v in trace.items()},
            "coast_km": {"min": float(np.nanmin(co)), "q10": float(np.nanquantile(co, 0.1)),
                         "median": float(np.nanmedian(co)), "n_le_20km": int((co <= 20).sum())},
            "rule": "B-12 B13: polygon descriptions classed by topsem.world.ems.flood_class; trace share = area of "
                    "flood-trace polygons over all flood polygons in evaluable cells, by the sensor of the EMS map"}


# %% B14: どちらでもないセル
def b14_neither() -> dict:
    rows = {r["tag"]: r for r in _json(INCIDENCE)["aois"] if r.get("included")}
    n_eval = n_drop = 0
    for a, r in rows.items():
        c = pq.read_table(LABELS / a / "cells.parquet").to_pandas()
        ev = evaluable(c)
        _, keep = riverine_labels(c, unspecified_as_riverine=bool(r.get("unspecified_as_riverine")))
        n_eval += int(ev.sum())
        n_drop += int((ev & ~keep).sum())
    return {"n_evaluable": n_eval, "n_neither": n_drop, "share": n_drop / n_eval if n_eval else None,
            "rule": "B-12 B14: evaluable cells of the 500 analysed AOIs (before the cap) that are neither flooded "
                    "nor dry under the riverine rule (flood cover >= 0.5 but riverine cover < 0.5)"}


# %% B15: HAND だけの地図の CSI
def _csi_by_rank(y: np.ndarray, score: np.ndarray, aoi: np.ndarray, rng: np.random.Generator) -> tuple:
    tp_all = fp_all = fn_all = 0
    per = []
    for a in np.unique(aoi):
        m = aoi == a
        yy, ss = y[m], score[m]
        k = int(yy.sum())
        if k == 0:
            continue
        tie = rng.random(len(ss))
        order = np.lexsort((tie, -ss))  # 得点の高い順、同点は無作為
        flag = np.zeros(len(ss), dtype=bool)
        flag[order[:k]] = True
        tp = int((flag & (yy == 1)).sum())
        fp, fn = int((flag & (yy == 0)).sum()), int((~flag & (yy == 1)).sum())
        per.append(tp / (tp + fp + fn))
        tp_all, fp_all, fn_all = tp_all + tp, fp_all + fp, fn_all + fn
    return float(np.median(per)), tp_all / (tp_all + fp_all + fn_all), len(per)


def b15_hand_csi() -> dict:
    v66 = _json(DOCV / "v66_b10_local.json")["B12"]
    return {"HAND_mean": v66["HAND_mean"], "HAND_min_ties_by_mean": v66["HAND_min"], "M1": v66["M1"],
            "M1c": v66["M1c"], "M3n": v66["M3n"],
            "rule": "B-12 B15: one definition for the HAND-only map: mean HAND of the cell, lowest first, ties in a "
                    "random order (v66 B12, corrected run); the other tables take this value. The minimum-HAND map "
                    "(ties by mean HAND) is a different map and is labelled as such where shown"}


# %% B16: JRC の図そのものの技能
def b16_jrc() -> dict:
    hits = sorted(LABELS.glob("*/jrc*.parquet"))
    if not hits:
        return {"note": "JRC hazard features are not available locally (computed on the cluster); not computed"}
    return {"note": f"{len(hits)} files found; computation not implemented for this layout"}


# %% B17: 定数の予測との比較
def b17_constant() -> dict:
    out = {}
    for label, run, split in (("primary", "world_v2", "primary"), ("country", "world_v2", "country"),
                              ("time", "world_v2", "temporal"), ("llmcut", "world_v2", "llmcut"),
                              ("activation_frozen", "world_v2_actgrp", "primary"),
                              ("activation_corrected", "world_v2_b7actfix", "primary")):
        d = v36.load(run, split)
        if d is None:
            out[label] = {"note": f"predictions of {run}/{split} not found locally"}
            continue
        y = d["y"].astype(np.float64)
        pbar = float(y.mean())
        r = {"n_cells": len(y), "brier_M0": float(np.mean((d["p_M0"] - y) ** 2)),
             "brier_M1": float(np.mean((d["p_M1"] - y) ** 2)), "brier_best_constant": pbar * (1 - pbar)}
        if "fold" in d and len(np.unique(d["fold"])) > 1:
            f = np.asarray(d["fold"])
            sse = 0.0
            for k in np.unique(f):
                m = f == k
                sse += float(((y[m] - y[~m].mean()) ** 2).sum())
            r["brier_training_constant"] = sse / len(y)
        r["M0_worse_than_best_constant"] = r["brier_M0"] > r["brier_best_constant"]
        out[label] = r
    out["rule"] = ("B-12 B17: held-out Brier of M0 and M1 against the best constant (the test cells' flooded share) "
                   "and, for k-fold splits, the training-side constant of each fold")
    return out


# %% B18: 市街地の 10 m の割合
def _lc_built(t: pd.DataFrame) -> np.ndarray:
    out = np.full(len(t), np.nan)
    for a, ia in t.groupby("aoi").indices.items():
        f = pq.read_table(LABELS / a / "features.parquet", columns=["lc_built"]).to_pandas()
        c = pq.read_table(LABELS / a / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        key = pd.MultiIndex.from_arrays([c["ix"].to_numpy(), c["iy"].to_numpy()])
        s = pd.Series(f["lc_built"].to_numpy(np.float64), index=key)
        q = pd.MultiIndex.from_arrays([t["ix"].to_numpy()[ia].astype(np.int64), t["iy"].to_numpy()[ia].astype(np.int64)])
        out[ia] = s.reindex(q).to_numpy()
    return out


def _pairs_with_extra(s: pd.DataFrame, extra: tuple[str, ...]) -> dict:
    """水系側 × 500 m 以内の OR を、extra の列を足した形と足さない形で、同じ場所グループの引き方（種 SEED、B_PAIRS 回）の
    pairs bootstrap で（度数の重みの形、v68.fit_w）。s は AOI の採否（0 < 率 < 1）を済ませた表。"""
    y = s["y"].to_numpy(dtype=np.float64)
    z0, names = v39.design(s, ("wetside",), True, banded=True)
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    rng = np.random.default_rng(SEED)
    sels = [rng.integers(0, len(unit_lab), size=len(unit_lab)) for _ in range(B_PAIRS)]
    j = names.index("wetside:le500m")
    out = {}
    for label, z in (("without", z0), ("with", np.column_stack([z0, s[list(extra)].to_numpy(np.float64)]))):
        cache = VER / f"v69_pairs_b18_{label}.npy"
        if cache.is_file():
            b = np.load(cache)
        else:
            _threads1()
            with ProcessPoolExecutor(N_WORKERS, initializer=v68._pairs_init_w,
                                     initargs=(y, z, aoi, unit, len(unit_lab))) as ex:
                res = list(ex.map(v68._pairs_rep_w, sels, chunksize=8))
            b = np.array([r[j] if r is not None else np.nan for r in res])
            np.save(cache, b)
        ok = np.isfinite(b)
        out[label] = {"pairs_ci95": [float(np.exp(np.quantile(b[ok], 0.025))),
                                     float(np.exp(np.quantile(b[ok], 0.975)))], "reps_ok": int(ok.sum())}
    return out


def b18_builtup() -> dict:
    t = v55.table()
    t["lc_built10"] = _lc_built(t)
    miss = float(np.isnan(t["lc_built10"]).mean())
    t["lc_built10"] = t["lc_built10"].fillna(0.0)
    y = t["y"].to_numpy(np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    s = t.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
    yy = s["y"].to_numpy(np.float64)
    out: dict = {"share_missing_lc_built": miss}
    for label, extra in (("without", ()), ("with_built10", ("lc_built10",))):
        z, names = v39.design(s, ("wetside",), True, banded=True)
        if extra:
            z = np.column_stack([z, s[list(extra)].to_numpy(np.float64)])
            names = [*names, *extra]
        res = v39.fit(yy, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
        tab = {}
        for nm in ("band_500m_1km", "band_gt1km", "wetside:le500m", *extra):
            j = names.index(nm)
            b, se = float(res["beta"][j]), float(res["se"][j])
            tab[nm] = {"coef": b, "odds_ratio": math.exp(b), "ci95": [math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)]}
        out[label] = tab
    out["pairs"] = _pairs_with_extra(s, ("lc_built10",))
    out["rule"] = ("B-12 B18: the Sect. 5 model (AOI fixed effects, distance bands, logit of M1, water-side × band) "
                   "on the capped table, with and without the ESA WorldCover built-up share of the cell (10 m)")
    return out


# %% B19: 仕様の曲線
GFM_DOMAIN_FILES = ("v41", "v45", "v50")  # GFM が評価できるセルの上の推定を出すファイル（EMS の推定もその範囲）


def b19_labels(file: str, path: str) -> tuple[str, str]:
    """B-12 の追記 2: ラベル（GFM か EMS か）は道の名で、EMS の推定には評価の範囲（GFM が評価できるセルか全てか）。"""
    segs = path.split("/")
    if any(x.startswith("odds_gfm") or x in ("gfm", "odds_v1") for x in segs):
        return "GFM", "gfm_cells"
    return "EMS", "gfm_cells" if file.split("_")[0] in GFM_DOMAIN_FILES else "all"


# 追記 3: 道の名に le500m が無い第 5 節の推定量の版のうち、水系側の区分の全体の OR を推定するもの（ファイル, 道の名）
B19_KEY_VERSIONS = (("v55", r"^odds_(ems|gfm)/wetside_adjusted_cci$"), ("v55", r"^partial/odds_"), ("v57", r"^B6/"),
                    ("v57", r"^B7/odds_(adjusted_h_main|same_cells_unadjusted)$"), ("v61", r"^without_"),
                    ("v64", r"^C[57]_near$"))


def b19_takes(file: str, path: str) -> bool:
    if "wetside" in path and "le500m" in path:
        return True
    pre = file.split("_")[0]
    return any(pre == f and re.search(rx, path) for f, rx in B19_KEY_VERSIONS)


def b19_spec_curve() -> dict:
    rows, seen = [], set()

    def add(row: dict) -> None:
        key = (round(row["estimate"], 3), round(row["ci95"][0], 3), round(row["ci95"][1], 3))
        if key not in seen:  # 追記 3: 点と区間が小数 3 桁で同じ行は 1 回
            seen.add(key)
            rows.append(row)

    with open(DOCV / "v69_estimates.csv", encoding="utf-8") as f:  # 追記 1 で直した目録（追記 2）
        for r in csv.DictReader(f):
            if (r["kind"] != "odds_ratio" or r["near_settlement"] != "True" or r["duplicate"] == "True"
                    or r["valid"] != "True"):
                continue
            p = r["path"]
            if not b19_takes(r["file"], p):
                continue
            lab, dom = b19_labels(r["file"], p)
            add({"file": r["file"], "path": p, "estimate": float(r["estimate"]),
                 "ci95": [float(r["ci95_lo"]), float(r["ci95_hi"])], "amendment": r["amendment"],
                 "labels": lab, "domain": dom, "gfm": lab == "GFM", "uncapped": "uncapped" in p})
    b2 = (_json(DOCV / "v68_b11_local.json").get("B2") or {})
    for k, lab, dom in (("ems_all", "EMS", "all"), ("gfm", "GFM", "gfm_cells"), ("ems_same_cells", "EMS", "gfm_cells"),
                        ("gfm_not_s1", "GFM", "gfm_cells"), ("ems_not_s1", "EMS", "gfm_cells")):
        x = b2.get(k) or {}
        if x.get("pairs_ci95"):
            add({"file": "v68_b11_local.json", "path": f"B2/{k}", "estimate": float(x["odds_ratio"]),
                 "ci95": [float(v) for v in x["pairs_ci95"]], "amendment": "11", "labels": lab,
                 "domain": dom, "gfm": lab == "GFM", "uncapped": True})
    rng = {}
    for lab, dom in (("GFM", "gfm_cells"), ("EMS", "gfm_cells"), ("EMS", "all")):
        e = [r["estimate"] for r in rows if r["labels"] == lab and r["domain"] == dom]
        if e:
            rng[f"{lab}_{dom}"] = {"n": len(e), "min": min(e), "max": max(e)}
    return {"rows": rows, "n": len(rows), "ranges": rng,
            "rule": "B-12 B19 (addenda 2 and 3): every version of the water-side odds ratio within 500 m in the "
                    "corrected census (v69_estimates.csv), including the Sect. 5 estimator versions flagged by keys, "
                    "without strata, single-activation leave-outs or other classes, with the uncapped estimates of B-11 "
                    "B2 (pairs intervals); equal values counted once; labels from the path, and for EMS the domain "
                    "(cells GFM can assess, or all cells)"}


# %% B20: p0
def b20_p0() -> dict:
    v44 = _json(DOCV / "v44_wikipedia_etymology.json")
    return {"p1": v44["truth_rate_in_wet"], "p0": v44["truth_rate_in_not"],
            "rule": "B-12 B20: water-related derivations among names read as water-side (p1) and among the others (p0)"}


# %% B21: 照合の標本の設計の改め
AUDIT_KEY_OLD = VER / "v68_audit_key.csv"
AUDIT_SAMPLE_OLD = DOCV / "v68_audit_sample.csv"
AUDIT_SAMPLE = DOCV / "v69_audit_sample.csv"
AUDIT_KEY = VER / "v69_audit_key.csv"


def b21_audit() -> dict:
    old = pd.read_csv(AUDIT_SAMPLE_OLD)
    key_old = pd.read_csv(AUDIT_KEY_OLD)
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    near = u["P_d1"].to_numpy() <= NEAR_M
    ws = u["lex_wetside_n1"].to_numpy() > 0
    fl = near & ws & (u["y"].to_numpy() == 1)
    carry = set(u.loc[fl, "aoi"])  # 集落の近くの推定を出す AOI（水系側の名の 500 m 以内に浸水セル）
    gn = pd.read_parquet(GN, columns=["gid", "name", "asciiname", "cc", "lat", "lon"]).drop_duplicates("gid")
    gn = gn.set_index("gid")
    cand = u.loc[near & ~ws & u["aoi"].isin(carry), "P_gid1"].astype(np.int64)
    cand = pd.Series(sorted(set(cand.tolist()) - set(old["gid"].astype(np.int64))))
    cand = cand[cand.isin(gn.index)]
    cc = gn.loc[cand.to_numpy(), "cc"].astype(str).to_numpy()
    # 国の層: 水系側の名の近くの浸水セルの国の割合に比例（最大の剰余で丸める）
    t_cc = gn["cc"].reindex(u.loc[fl, "P_gid1"].astype(np.int64).to_numpy()).fillna("").to_numpy(str)
    share = pd.Series(t_cc).value_counts(normalize=True)
    share = share[share.index.isin(set(cc))]
    share = share / share.sum()
    raw = share * AUDIT_NOT_N
    alloc = np.floor(raw).astype(int)
    for c in (raw - alloc).sort_values(ascending=False).index[: AUDIT_NOT_N - int(alloc.sum())]:
        alloc[c] += 1
    rng = np.random.default_rng(SEED)
    pick = []
    for c, k in alloc.items():
        pool = cand.to_numpy()[cc == c]
        k = min(int(k), len(pool))
        if k:
            pick += rng.choice(pool, size=k, replace=False).tolist()
    new = gn.loc[pick].reset_index()[["gid", "name", "asciiname", "cc", "lat", "lon"]]
    allr = pd.concat([old[["gid", "name", "asciiname", "cc", "lat", "lon"]], new], ignore_index=True)
    allr = allr.sample(frac=1.0, random_state=int(SEED % (2**31))).reset_index(drop=True)
    allr.insert(0, "item", np.arange(1, len(allr) + 1))
    allr.to_csv(AUDIT_SAMPLE, index=False)
    kn = pd.DataFrame({"gid": new["gid"], "stratum_form": "", "stratum_cc": new["cc"], "stratum_weight": np.nan,
                       "read_as_water_side": 0})
    ko = key_old.assign(read_as_water_side=1)
    key = pd.concat([ko, kn], ignore_index=True)
    key.to_csv(AUDIT_KEY, index=False)
    return {"n_names": len(allr), "n_read_water_side": len(ko), "n_not_read": len(kn),
            "not_read_by_country": {str(k): int(v) for k, v in alloc.items() if v},
            "n_carrying_aoi": len(carry),
            "rule": "B-12 B21: names not read as water-side, nearest to cells within 500 m in the AOIs that hold "
                    "flooded cells near water-side names, allocated to countries in proportion to those flooded "
                    "cells; mixed with the B-11 sample in a random order; the key (read or not) is kept private"}


# %% B8c（追記 1 の 1・2）: 目録の「集落の近く」の印の直しと BY の q
NEAR_KEYS = ("n_group_flooded", "n_group_cells", "n_flooded_near")


def _walk_fixed(node, path: list[str], rows: list[dict], file: str, amend: str) -> None:
    """v68._walk と同じ規則で推定を集め、第 5 節の推定量の出力（NEAR_KEYS を持つ節）にも「集落の近く」の印を付ける。"""
    n0 = len(rows)
    if isinstance(node, dict):
        flat = {k: v for k, v in node.items()
                if not isinstance(v, (dict, list)) or (isinstance(v, list) and all(isinstance(x, (int, float)) for x in v))}
        v68._walk(flat, path, rows, file, amend)
        if len(rows) > n0 and rows[-1]["kind"] == "odds_ratio" and any(k in node for k in NEAR_KEYS):
            rows[-1]["near_settlement"] = True
        for k, v in node.items():
            if isinstance(v, (dict, list)):
                _walk_fixed(v, [*path, str(k)], rows, file, amend)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            if isinstance(v, (dict, list)):
                _walk_fixed(v, [*path, str(i)], rows, file, amend)


def _bh_by(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = len(p)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    c = float(np.sum(1.0 / np.arange(1, m + 1)))
    return np.minimum(q, 1.0), np.minimum(q * c, 1.0)


def b8c_census() -> dict:
    import re

    rows: list[dict] = []
    for f in sorted(DOCV.glob("v[3-6][0-9]*_*.json")):
        m = re.match(r"v(\d+)", f.name)
        v = int(m.group(1)) if m else 0
        if not 36 <= v <= 67:
            continue
        _walk_fixed(_json(f), [], rows, f.name, v68._amend(v))
    df = pd.DataFrame(rows)
    df["key"] = (df["kind"] + "|" + df["estimate"].round(6).astype(str) + "|" + df["ci95_lo"].round(6).astype(str)
                 + "|" + df["ci95_hi"].round(6).astype(str))
    df["duplicate"] = df.duplicated("key")
    sep = (df["kind"] == "odds_ratio") & ((np.abs(np.log(df["estimate"].clip(lower=1e-300))) > 10)
                                          | (df["ci95_hi"] / df["ci95_lo"].clip(lower=1e-300) > 1e6))
    df["valid"] = np.isfinite(df["ci95_lo"]) & np.isfinite(df["ci95_hi"]) & (df["ci95_lo"] <= df["ci95_hi"]) & ~sep
    df["family"] = [v68._family(p, k) for p, k in zip(df["path"], df["kind"], strict=True)]
    base = df.loc[~df["duplicate"] & df["valid"] & df["name_effect"]].copy()
    fl = base.loc[(base["family"] == "flood") & base["near_settlement"] & base["p_approx"].notna()].copy()
    q_bh, q_by = _bh_by(fl["p_approx"].to_numpy())
    fl["q_bh"], fl["q_by"] = q_bh, q_by
    df = df.merge(fl[["file", "path", "q_bh", "q_by"]], on=["file", "path"], how="left")
    df.drop(columns="key").to_csv(DOCV / "v69_estimates.csv", index=False)
    old = pd.read_csv(DOCV / "v68_estimates.csv")
    out: dict = {"n_rows": len(df), "n_unique": int((~df["duplicate"]).sum())}
    for fam in ("flood", "site", "nuisance"):
        u = base.loc[base["family"] == fam]
        out[fam] = {"n": len(u), "n_excluding_null": int(u["excludes_null"].sum()),
                    "share_excluding_null": float(u["excludes_null"].mean()) if len(u) else None,
                    "by_kind": {k: {"n": len(g), "excl": int(g["excludes_null"].sum())} for k, g in u.groupby("kind")}}
    top = fl.sort_values("q_bh").head(10)
    out["near_settlement_or"] = {
        "n": len(fl), "n_excluding_one": int(fl["excludes_null"].sum()),
        "n_bh_q_below_0_05": int((q_bh < 0.05).sum()), "n_bh_q_below_0_10": int((q_bh < 0.10).sum()),
        "n_by_q_below_0_05": int((q_by < 0.05).sum()), "n_by_q_below_0_10": int((q_by < 0.10).sum()),
        "smallest_q": [{"file": r.file, "path": r.path, "estimate": r.estimate, "ci95": [r.ci95_lo, r.ci95_hi],
                        "p": r.p_approx, "q_bh": r.q_bh, "q_by": r.q_by} for r in top.itertuples()],
        "n_flagged_by_keys_only": int(len(fl) - int(((old["near_settlement"]) & (old["family"] == "flood")
                                                     & (~old["duplicate"]) & (old["valid"]) & (old["name_effect"])
                                                     & old["p_approx"].notna()).sum()))}
    out["rule"] = ("B-12 addendum 1: the census of B-11 B8 with odds ratios flagged near settlements also when the node "
                   "carries " + ", ".join(NEAR_KEYS) + "; BH and BY q within the near-settlement family")
    return out


# %% B21b（追記 1 の 3）: 照合の標本に上限なしの標本の形を足す
AUDIT_ADD_N = 50
AUDIT_SAMPLE2 = DOCV / "v69_audit_sample.csv"


def b21b_audit_add() -> dict:
    from topsem.world.lexicon import _SPLIT, _hit

    key_old = pd.read_csv(AUDIT_KEY_OLD)
    have = set(zip(key_old["stratum_form"].astype(str), key_old["stratum_cc"].astype(str), strict=True))
    lex = v68.lexicon()
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    gn = pd.read_parquet(GN, columns=["gid", "name", "asciiname", "cc", "lat", "lon"]).drop_duplicates("gid")
    gn = gn.set_index("gid")
    gids = u["P_gid1"].to_numpy(np.int64)
    ws = u["lex_wetside_n1"].to_numpy() > 0
    fl = ws & (u["P_d1"].to_numpy() <= NEAR_M) & (u["y"].to_numpy() == 1)
    per = pd.DataFrame({"gid": gids[ws], "fl": fl[ws]}).groupby("gid").agg(cells=("fl", "size"), flooded=("fl", "sum"))
    rows = []
    for gid, r in per.iterrows():
        if gid not in gn.index:
            continue
        name, cc = str(gn.at[gid, "asciiname"]), str(gn.at[gid, "cc"])
        toks = [x for x in _SPLIT.split(name.lower()) if x]
        forms = sorted({e.form for e in lex.entries_for(cc) if e.cat in WET_SIDE for tok in toks if _hit(tok, e)})
        for f_ in forms or ["?"]:
            rows.append({"gid": int(gid), "cc": cc, "form": f_, "cells": int(r["cells"]), "flooded": int(r["flooded"])})
    d = pd.DataFrame(rows)
    st = d.groupby(["form", "cc"]).agg(cells=("cells", "sum"), flooded=("flooded", "sum"), n=("gid", "nunique"))
    st["w"] = st["cells"] / st["cells"].sum() + (st["flooded"] / st["flooded"].sum() if st["flooded"].sum() else 0)
    st = st.sort_values("w", ascending=False)
    cur = pd.read_csv(AUDIT_SAMPLE)
    taken = set(cur["gid"].astype(np.int64))
    rng = np.random.default_rng(SEED + 1)
    picked, key = [], []
    for (form, cc), r in st.iterrows():
        if len(picked) >= AUDIT_ADD_N:
            break
        if (str(form), str(cc)) in have:
            continue
        cand = [g for g in d.loc[(d["form"] == form) & (d["cc"] == cc), "gid"].unique() if g not in taken]
        take = rng.choice(cand, size=min(3, len(cand), AUDIT_ADD_N - len(picked)), replace=False) if cand else []
        for g in take:
            g = int(g)
            taken.add(g)
            picked.append({"gid": g, "name": str(gn.at[g, "name"]), "asciiname": str(gn.at[g, "asciiname"]), "cc": cc,
                           "lat": float(gn.at[g, "lat"]), "lon": float(gn.at[g, "lon"])})
            key.append({"gid": g, "stratum_form": form, "stratum_cc": cc, "stratum_weight": float(r["w"]),
                        "read_as_water_side": 1})
    allr = pd.concat([cur.drop(columns="item"), pd.DataFrame(picked)], ignore_index=True)
    allr = allr.sample(frac=1.0, random_state=int((SEED + 1) % (2**31))).reset_index(drop=True)
    allr.insert(0, "item", np.arange(1, len(allr) + 1))
    allr.to_csv(AUDIT_SAMPLE2, index=False)
    k = pd.concat([pd.read_csv(AUDIT_KEY), pd.DataFrame(key)], ignore_index=True)
    k.to_csv(AUDIT_KEY, index=False)
    strata = sorted({(r["stratum_form"], r["stratum_cc"]) for r in key})
    w_cov = float(st.loc[strata, "w"].sum()) / float(st["w"].sum())
    return {"n_added": len(picked), "n_names": len(allr),
            "added_by_country": pd.Series([p["cc"] for p in picked]).value_counts().to_dict(),
            "added_strata": sorted({f"{r['stratum_form']}|{r['stratum_cc']}" for r in key}),
            "share_uncapped_weight_of_added_strata": w_cov,
            "rule": "B-12 addendum 1 item 3: strata (form x country) of the uncapped sample absent from the B-11 sample, "
                    "heaviest first, up to three names each, 50 names; mixed into the B21 sample in a random order"}


# %% 実行
ITEMS = (("B1", b1_operating), ("B3", b3_language), ("B4", b4_inclusion), ("B5", b5_null_outliers),
         ("B7r", b7r_recurrence), ("B8", b8_flexible_rows), ("B9", b9_coder2), ("B10", b10_reference),
         ("B12", b12_gfm_window), ("B13", b13_typing), ("B14", b14_neither), ("B15", b15_hand_csi),
         ("B16", b16_jrc), ("B17", b17_constant), ("B19", b19_spec_curve), ("B20", b20_p0), ("B6", b6_highland),
         ("B18", b18_builtup), ("B21", b21_audit), ("B8c", b8c_census), ("B21b", b21b_audit_add),
         ("B2", b2_composition), ("B11", b11_names_only_penalty))


def main() -> None:
    global N_WORKERS
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else set()
    if "--workers" in sys.argv:
        N_WORKERS = int(sys.argv[sys.argv.index("--workers") + 1])
    for k, fn in ITEMS:
        if only and k not in only:
            continue
        t0 = time.time()
        r = fn()
        print(k, f"({time.time() - t0:.0f} s)", json.dumps(r, default=float)[:800], flush=True)
        cur = _json(OUT) if OUT.is_file() else {}
        cur |= {"rule": "amendment B-12 items B1-B21 (local; post hoc; confirmatory verdict unchanged)",
                "seed": SEED, k: r}
        OUT.write_text(json.dumps(cur, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
