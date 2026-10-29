# %% v73（事後・amendment B-14 の B1〜B15）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-14（追記 2）で計算の前に固定した量を、既存の held-out の予測・二段目の保存・表から計算する。

- B1: 直した検定の表（(a)〜(d) の L± の限界を百分位と正規近似で、TOST の p、片側 95% の上限を R・Brier・注入の OR で、
  百分位の限界のモンテカルロ誤差、改訂 11 の二方向の区間）
- B2: 重い裾への頑健性（jackknife、studentized bootstrap、G*、上位 2 を除いた R と重みを揃えた R）
- B3: 国の中の activation の二段の再抽出での動作特性（L 規則と上限だけの規則）
- B4: pigeonhole の再抽出（場所グループ × activation、二段目を当て直す）
- B5: 基準の分解と組成を含まない比較
- B6: 保たれる割合（集落を通じて ÷ セルに直接）
- B7: 名の特徴の held-out の係数（EMS・GFM）
- B8: 地名学の表
- B9: 修飾語の対の標高と川の大きさ
- B10: 目録の向き
- B11: 0.84 の再標本の分布
- B12: 余白を Brier skill で
- B13: 集落の近くの推定の点検（上限ありの表）
- B14: 位置で結ぶ表の整合
- B15: AOI ごとの分布と activation を一つずつ除いた R

出力: stdout と docs/verification/v73_b14_local.json（部分ごとに足す。`--only B1,B5` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import stage2_weighted
from topsem.world.b13 import SESOI_R, fieller_interval, l_pm
from topsem.world.lexicon import WET_SIDE, Lexicon, _hit

OUT = ROOT / "docs" / "verification" / "v73_b14_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
LEX_CORR = DOCV / "v70_lexicon_corrected.yaml"
SEED = 20261103
B_FIXED = 4000
B_REFIT = 500
B_OUT, B_IN = 1000, 500
B_MC = 1000
B_PAIR_BOOT = 2000
Z95 = 1.6448536269514722  # 片側 95%
CELL_KM2 = 0.0625
N_WORKERS = int(os.environ.get("V73_WORKERS", "3"))
CORR = ("world_v2_lexfix", "primary")
ARMS4 = ("M0", "M1", "M1c", "M3n")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v70 = _load("v70_b13_local")
v69, v36, v68 = v70.v69, v70.v36, v70.v68
v39 = v69.v39


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _act(aoi: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(aoi).astype(str)).str.split("_").str[0].to_numpy()


def _cell_terms(d: dict) -> tuple[np.ndarray, np.ndarray]:
    """セルごとの名の増分 e = (M1c の二乗誤差) − (M3n の) と物差し d = (M0 の) − (M1 の)。"""
    y = d["y"].astype(np.float64)
    sq = {a: (d[f"p_{a}"].astype(np.float64) - y) ** 2 for a in ARMS4}
    return sq["M1c"] - sq["M3n"], sq["M0"] - sq["M1"]


def _gsum(key: np.ndarray, *vals: np.ndarray) -> tuple[np.ndarray, ...]:
    u, inv = np.unique(np.asarray(key).astype(str), return_inverse=True)
    g = len(u)
    return (u, np.bincount(inv, minlength=g).astype(np.float64),
            *(np.bincount(inv, weights=v, minlength=g) for v in vals))


def _boot_w(rng: np.random.Generator, g: int, b: int = B_FIXED) -> np.ndarray:
    return rng.multinomial(g, np.full(g, 1.0 / g), size=b).astype(np.float64)


# %% B1: 直した検定の表
def _levels_detail(e0: float, d0: float, eb: np.ndarray, db: np.ndarray, rng: np.random.Generator) -> dict:
    """百分位と正規近似の L± の限界、TOST の p、Fieller の区間、百分位の限界のモンテカルロ誤差。"""
    ok = np.isfinite(eb) & np.isfinite(db)
    eb, db = eb[ok], db[ok]
    lp, lm = l_pm(eb, db)
    lp0, lm0 = (float(x) for x in l_pm(e0, d0))
    se_p, se_m = float(np.std(lp, ddof=1)), float(np.std(lm, ddof=1))
    c = np.cov(eb, db)
    fi = fieller_interval(e0, d0, (float(c[0, 0]), float(c[1, 1]), float(c[0, 1])))
    idx = rng.integers(0, len(lp), size=(B_MC, len(lp)))
    mc_up = np.quantile(lp[idx], 0.95, axis=1)
    mc_lo = np.quantile(lm[idx], 0.05, axis=1)
    pct = {"l_plus_upper": float(np.quantile(lp, 0.95)), "l_minus_lower": float(np.quantile(lm, 0.05))}
    nrm = {"l_plus_upper": lp0 + Z95 * se_p, "l_minus_lower": lm0 - Z95 * se_m}
    return {"R": e0 / d0, "delta_bs": e0, "yardstick": d0, "L_plus": lp0, "L_minus": lm0, "reps": len(lp),
            "percentile": pct | {"equivalent": bool(pct["l_plus_upper"] < 0 < pct["l_minus_lower"]),
                                 "mc_sd_l_plus_upper": float(np.std(mc_up, ddof=1)),
                                 "mc_sd_l_minus_lower": float(np.std(mc_lo, ddof=1))},
            "normal": nrm | {"equivalent": bool(nrm["l_plus_upper"] < 0 < nrm["l_minus_lower"]),
                             "se_l_plus": se_p, "se_l_minus": se_m},
            "tost_p": float(max(stats.norm.cdf(lp0 / se_p), stats.norm.sf(lm0 / se_m))),
            "fieller90": list(fi) if fi else None, "upper_R": float(fi[1]) if fi else None,
            "upper_brier": float(fi[1] * d0) if fi else None,
            "var": [float(c[0, 0]), float(c[1, 1]), float(c[0, 1])]}


def _or_conversions(uppers: dict[str, float]) -> dict:
    """各上限（R）に届く注入の OR（B-13 の B2 と同じ曲線・基準・換算）。"""
    v64 = v70._load("v64_b9_aggregate")
    v59 = v64.v59
    grids, gains, _ = v64.combined_curves()
    refs = v59.references()
    yard = float(_json(RUNS_DIR / "world_v2_lexfix" / "metrics_primary.json")["primary"]["delta_bs_M0_M1"])
    _, conv = v70._conv_factory()
    rng = np.random.default_rng(SEED)
    out: dict = {"yardstick": yard, "references": refs, "bounds": {}}
    for kind, ur in uppers.items():
        out["bounds"][kind] = {}
        for ref, dnull in refs.items():
            grid = np.concatenate([[0.0], grids["A"]])
            e = v59._entry(gains["A"], grid, "bound", ur * yard - dnull, float(np.exp(grids["A"][-1])), rng)
            e["estimator_scale"] = conv("A", e.get("odds_ratio"))
            out["bounds"][kind][ref] = e
        ors = [b["odds_ratio"] for b in out["bounds"][kind].values() if b.get("odds_ratio")]
        est = [b["estimator_scale"]["estimator"] for b in out["bounds"][kind].values()
               if b.get("estimator_scale") and b["estimator_scale"].get("estimator")]
        out["bounds"][kind]["range"] = {"nominal": [min(ors), max(ors)] if ors else None,
                                        "estimator": [min(est), max(est)] if est else None}
    return out


def b1_table() -> dict:
    v70j = _json(DOCV / "v70_b13_local.json")
    rng = np.random.default_rng(SEED)
    out: dict = {"levels": {}}
    # (a): B-13 の B1(a) の重みを同じ種・同じ順で作り直す（確証の実行の分を先に引く）
    ra = np.random.default_rng(v70.SEED)
    conf = v36.load("world_v2", "primary")
    n_c, _, _ = v69._group_sums(conf, ("M0",))
    ra.multinomial(len(n_c), np.full(len(n_c), 1.0 / len(n_c)), size=v70.B_FIXED)
    lex = v36.load(*CORR)
    n, s, units = v69._group_sums(lex, ARMS4)
    g = len(units)
    w = ra.multinomial(g, np.full(g, 1.0 / g), size=v70.B_FIXED).astype(np.float64)
    eg, dg = s["M1c"] - s["M3n"], s["M0"] - s["M1"]
    big_n = float(n.sum())
    e0, d0 = float(eg.sum() / big_n), float(dg.sum() / big_n)
    nw = w @ n
    la = _levels_detail(e0, d0, (w @ eg) / nw, (w @ dg) / nw, rng)
    saved = v70j["B1a"]["corrected"]["fieller90"]
    la["replay_matches_b13"] = bool(np.allclose(la["fieller90"], saved, rtol=0, atol=1e-12))
    out["levels"]["a_fixed_groups"] = la
    # (b)(c): 保存した再標本
    lx = v70j["B1bc"]["b_groups"]["refit"]
    for k, scheme in (("b_refit_groups", "groups"), ("c_refit_nested", "nested")):
        r = np.load(VER / f"v70_b1_refit_reps_{scheme}.npy")
        wt = r[:, 6]
        out["levels"][k] = _levels_detail(lx["delta_bs"], lx["yardstick"], (r[:, 0] - r[:, 1]) / wt,
                                          (r[:, 4] - r[:, 5]) / wt, rng)
    # (d): 正規近似だけ
    d = v70j["B1d"]
    c = v70j["B1bc"]["c_nested"]["refit"]
    ve, vd, vc = c["var"]
    nb = d["nb_var"]
    se_p = math.sqrt(ve + SESOI_R**2 * vd - 2 * SESOI_R * vc + nb["L_plus"])
    se_m = math.sqrt(ve + SESOI_R**2 * vd + 2 * SESOI_R * vc + nb["L_minus"])
    out["levels"]["d_refit_nested_folds"] = {
        "R": c["R"], "delta_bs": c["delta_bs"], "yardstick": c["yardstick"],
        "L_plus": d["L_plus"], "L_minus": d["L_minus"], "percentile": None,
        "normal": {"l_plus_upper": d["l_plus_upper"], "l_minus_lower": d["l_minus_lower"],
                   "equivalent": d["equivalent"], "se_l_plus": se_p, "se_l_minus": se_m},
        "tost_p": float(max(stats.norm.cdf(d["L_plus"] / se_p), stats.norm.sf(d["L_minus"] / se_m))),
        "fieller90": d["fieller90"], "upper_R": d["fieller90"][1], "upper_brier": d["fieller90"][1] * c["yardstick"]}
    tw = _json(DOCV / "v68_b11_local.json")["B7"]["two_way"]
    out["b11_two_way_refit"] = {"ratio_pct90": tw["ci90_R_refit_ratio"], "upper_R": tw["ci90_R_refit_ratio"][1],
                                "upper_brier": tw["ci90_R_refit_ratio"][1] * c["yardstick"],
                                "note": "B-11 B7 saved only the percentile interval of the ratio (no replicates)"}
    uppers = {k: v["upper_R"] for k, v in out["levels"].items() if v.get("upper_R") is not None}
    uppers["b11_two_way_refit"] = out["b11_two_way_refit"]["upper_R"]
    out["injected_or"] = _or_conversions(uppers)
    out["rule"] = ("B-14 B1 (addendum 2): (a) weights of B-13 B1(a) redrawn and checked; (b)(c) saved replicates; (d) "
                   "normal approximation only; TOST p = max(Phi(L+/se+), 1 - Phi(L-/se-)); one-sided 95% upper limit "
                   "of R = upper end of the 90% Fieller interval; Monte Carlo SD of percentile limits from 1,000 "
                   "resamples of the replicates; injected OR by the curves of B-13 B2")
    return out


# %% B2: 重い裾への頑健性
def _jackknife(key: np.ndarray, ec: np.ndarray, dc: np.ndarray) -> dict:
    units, n, e, d = _gsum(key, ec, dc)
    big_n, big_e, big_d = n.sum(), e.sum(), d.sum()
    g = len(units)
    te, td = (big_e - e) / (big_n - n), (big_d - d) / (big_n - n)
    lp, lm = l_pm(te, td)
    e0, d0 = big_e / big_n, big_d / big_n
    lp0, lm0 = (float(x) for x in l_pm(e0, d0))

    def jv(x: np.ndarray) -> float:
        return float((g - 1) / g * np.sum((x - x.mean()) ** 2))

    cov = float((g - 1) / g * np.sum((te - te.mean()) * (td - td.mean())))
    fi = fieller_interval(e0, d0, (jv(te), jv(td), cov))
    r_minus = (big_e - e) / (big_d - d)
    psi = e0 / d0 - r_minus
    top = np.argsort(-np.abs(psi))[:5]
    return {"n_clusters": g, "R": e0 / d0, "se_delta_bs": math.sqrt(jv(te)), "se_yardstick": math.sqrt(jv(td)),
            "l_plus_upper": lp0 + Z95 * math.sqrt(jv(lp)), "l_minus_lower": lm0 - Z95 * math.sqrt(jv(lm)),
            "equivalent": bool(lp0 + Z95 * math.sqrt(jv(lp)) < 0 < lm0 - Z95 * math.sqrt(jv(lm))),
            "fieller90": list(fi) if fi else None,
            "top5": [{"cluster": str(units[i]), "R_without": float(r_minus[i]), "psi": float(psi[i]),
                      "share_cells": float(n[i] / big_n)} for i in top]}


def _studentized(key: np.ndarray, ec: np.ndarray, dc: np.ndarray, rng: np.random.Generator) -> dict:
    units, n, e, d = _gsum(key, ec, dc)
    g = len(units)
    w = _boot_w(rng, g)
    out: dict = {"n_clusters": g, "replicates": B_FIXED}
    for side, sgn in (("plus", -1.0), ("minus", 1.0)):
        lg = e + sgn * SESOI_R * d
        big_n = n.sum()
        lhat = float(lg.sum() / big_n)
        se_hat = math.sqrt(g / (g - 1) * np.sum((lg - lhat * n) ** 2)) / big_n
        nw = w @ n
        ls = (w @ lg) / nw
        ss = (w @ lg**2) - 2 * ls * (w @ (lg * n)) + ls**2 * (w @ n**2)
        se_s = np.sqrt(g / (g - 1) * np.clip(ss, 0, None)) / nw
        t = (ls - lhat) / se_s
        if side == "plus":
            lim = lhat - float(np.quantile(t, 0.05)) * se_hat
            out["l_plus_upper"] = lim
        else:
            lim = lhat - float(np.quantile(t, 0.95)) * se_hat
            out["l_minus_lower"] = lim
        out[f"L_{side}"] = lhat
        out[f"se_L_{side}"] = se_hat
        out[f"t_quantiles_{side}"] = np.quantile(t, [0.05, 0.5, 0.95]).tolist()
    out["equivalent"] = bool(out["l_plus_upper"] < 0 < out["l_minus_lower"])
    return out


def _g_star(n: np.ndarray, e: np.ndarray) -> dict:
    big_n = n.sum()
    dhat = e.sum() / big_n
    gam = ((e - dhat * n) / big_n) ** 2
    g = len(n)
    big_g = float(np.mean((gam - gam.mean()) ** 2) / gam.mean() ** 2)
    return {"G": g, "Gamma": big_g, "G_star": g / (1 + big_g),
            "share_variance_top5": float(np.sort(gam)[::-1][:5].sum() / gam.sum())}


def _weighted_r(lex: dict, omega: np.ndarray, rng: np.random.Generator, keep: np.ndarray | None = None) -> dict:
    ec, dc = _cell_terms(lex)
    unit = np.asarray(lex["unit"]).astype(str)
    if keep is not None:
        ec, dc, unit, omega = ec[keep], dc[keep], unit[keep], omega[keep]
    _, nw_, ew, dw = _gsum(unit, ec * omega, dc * omega)
    nw_ = np.bincount(np.unique(unit, return_inverse=True)[1], weights=omega)
    w = _boot_w(rng, len(nw_))
    den = w @ nw_
    return v70._summ(float(ew.sum() / nw_.sum()), float(dw.sum() / nw_.sum()), (w @ ew) / den, (w @ dw) / den)


def b2_heavy_tails() -> dict:
    rng = np.random.default_rng(SEED)
    lex = v36.load(*CORR)
    ec, dc = _cell_terms(lex)
    unit = np.asarray(lex["unit"]).astype(str)
    aoi = np.asarray(lex["aoi"]).astype(str)
    act = _act(aoi)
    out: dict = {"jackknife": {"groups": _jackknife(unit, ec, dc), "activations": _jackknife(act, ec, dc)},
                 "studentized_groups": _studentized(unit, ec, dc, rng)}
    gs: dict = {}
    for label, run, split in (("corrected_primary", *CORR), ("confirmatory_primary", "world_v2", "primary"),
                              ("confirmatory_actgrp", "world_v2_actgrp", "primary"),
                              ("confirmatory_country", "world_v2", "country")):
        d = v36.load(run, split)
        if d is None:
            gs[label] = {"note": f"{run}/{split} not found"}
            continue
        e_c, d_c = _cell_terms(d)
        _, n, e, dd = _gsum(np.asarray(d["unit"]).astype(str), e_c, d_c)
        gs[label] = {"run": run, "split": split, "delta_bs": _g_star(n, e), "yardstick": _g_star(n, dd)}
    out["g_star"] = gs
    b6 = _json(DOCV / "v70_b13_local2.json")["B6"]["corrected"]["top10"][:2]
    top2 = [x["unit"] for x in b6]
    keep = ~np.isin(unit, top2)
    out["without_top2"] = {"units": top2, "countries": [x["country"] for x in b6],
                           **_weighted_r(lex, np.ones(len(unit)), rng, keep)}
    for lab, key in (("groups", unit), ("aois", aoi), ("activations", act)):
        _, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        out[f"equal_weight_{lab}"] = _weighted_r(lex, 1.0 / cnt[inv], rng)
    out["rule"] = ("B-14 B2 (addendum 2): corrected run, fixed models; (i) jackknife over location groups and over "
                   "activations; (ii) studentized bootstrap over location groups with a robust variance per replicate; "
                   "(iii) G* of Carter, Schnepel and Steigerwald (2017) for the location-group contributions; (iv) R "
                   "without the two most influential groups of B-13 B6 and with equal weights per group, AOI and "
                   "activation (Fieller intervals from 4,000 location-group replicates)")
    return out


# %% B3: 動作特性（国の中の activation の二段）
R0S = (-0.05, 0.0, 0.05)


def b3_operating() -> dict:
    lex = v36.load(*CORR)
    ec, dc = _cell_terms(lex)
    act = _act(np.asarray(lex["aoi"]))
    cc = np.asarray(lex["cc"]).astype(str)
    acts, n_a, e_a, d_a = _gsum(act, ec, dc)
    coa = (pd.DataFrame({"a": act, "c": cc}).groupby("a")["c"].agg(lambda s: s.value_counts().index[0])
           .reindex(acts).to_numpy())
    countries, ic = np.unique(coa, return_inverse=True)
    members = [np.flatnonzero(ic == k) for k in range(len(countries))]
    n_c = len(countries)
    e_hat, d_hat = float(e_a.sum() / n_a.sum()), float(d_a.sum() / n_a.sum())
    rng = np.random.default_rng(SEED)
    cnt = {sh: {r0: {"L_equivalent": 0, "upper_nonsuperior": 0} for r0 in R0S} for sh in ("proportional", "uniform")}
    gate = 0
    t0 = time.time()
    for b in range(B_OUT):
        kc = rng.integers(0, n_c, size=n_c)
        copies = [members[k][rng.integers(0, len(members[k]), size=len(members[k]))] for k in kc]
        m = rng.multinomial(n_c, np.full(n_c, 1.0 / n_c), size=B_IN)
        big_e, big_d, big_n = np.zeros(B_IN), np.zeros(B_IN), np.zeros(B_IN)
        for j, arr in enumerate(copies):
            k = len(arr)
            c = rng.multinomial(m[:, j] * k, np.full(k, 1.0 / k)).astype(np.float64)
            big_e += c @ e_a[arr]
            big_d += c @ d_a[arr]
            big_n += c @ n_a[arr]
        ev, dv = big_e / big_n, big_d / big_n
        gate += int(np.quantile(dv, 0.025) <= 0)
        for r0 in R0S:
            for sh in ("proportional", "uniform"):
                e_in = ev - e_hat + r0 * (dv if sh == "proportional" else d_hat)
                lp, lm = l_pm(e_in, dv)
                up, lo = np.quantile(lp, 0.95), np.quantile(lm, 0.05)
                cnt[sh][r0]["L_equivalent"] += int(up < 0 < lo)
                cnt[sh][r0]["upper_nonsuperior"] += int(up < 0)
        if b % 100 == 0:
            print(f"  B3 outer {b} ({time.time() - t0:.0f} s)", flush=True)

    def rate(k: int) -> dict:
        p = k / B_OUT
        return {"share": p, "mc_se": math.sqrt(p * (1 - p) / B_OUT)}

    res = {sh: {f"{r0:+.2f}": {k: rate(v) for k, v in c.items()} for r0, c in cs.items()} for sh, cs in cnt.items()}
    return {"n_activations": len(acts), "n_countries_with_activations": n_c, "replicates_outer": B_OUT,
            "replicates_inner": B_IN, "seed": SEED, "R_hat": e_hat / d_hat, "share_outer_denominator_lo_le0":
            gate / B_OUT, "rates": res, "seconds": round(time.time() - t0),
            "rule": "B-14 B3 (addendum 2): corrected run, fixed models; countries with activations, then activations "
                    "within each drawn copy of a country (outer and inner alike); L rule = one-sided 95% upper limit "
                    "of L+ below 0 and lower limit of L- above 0; upper-only rule = upper limit of L+ below 0"}


# %% B4: pigeonhole（二段目を当て直す）
_B: dict = {}


def _b4_init() -> None:
    v68._b7_init()
    _B.update(v68._B7)


def _b4_rep(b: int) -> list[float]:
    rng = np.random.default_rng(SEED + 200_000 + b)
    wu = rng.multinomial(_B["nu"], np.full(_B["nu"], 1 / _B["nu"])).astype(np.float64)
    wa = rng.multinomial(_B["na"], np.full(_B["na"], 1 / _B["na"])).astype(np.float64)
    acc = np.zeros(7)
    for f in _B["folds"]:
        w_tr, w_te = wu[f["u_tr"]] * wa[f["a_tr"]], wu[f["u_te"]] * wa[f["a_te"]]
        if w_te.sum() == 0 or w_tr.sum() == 0:
            continue
        pc, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xc_tr"], f["xc_te"], f["y_tr"], w_tr=w_tr, init=f["init_c"])
        pn, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xn_tr"], f["xn_te"], f["y_tr"], w_tr=w_tr, init=f["init_n"])
        y = f["y_te"].astype(np.float64)
        for i, p in enumerate((pc, pn, f["p_M1c_te"], f["p_M3n_te"], f["p_M0_te"], f["p_M1_te"])):
            acc[i] += float((w_te * (p - y) ** 2).sum())
        acc[6] += float(w_te.sum())
    return acc.tolist()


def b4_pigeonhole() -> dict:
    f = VER / "v73_b4_pigeonhole_reps.npy"
    if f.is_file():
        r = np.load(f)
    else:
        v69._threads1()
        t0 = time.time()
        with ProcessPoolExecutor(N_WORKERS, initializer=_b4_init) as ex:
            r = np.array(list(ex.map(_b4_rep, range(B_REFIT), chunksize=2)), dtype=np.float64)
        np.save(f, r)
        print(f"  B4: {len(r)} replicates ({time.time() - t0:.0f} s)", flush=True)
    c = _json(DOCV / "v70_b13_local.json")["B1bc"]["b_groups"]["refit"]
    w = r[:, 6]
    rng = np.random.default_rng(SEED)
    out = {"replicates": len(r), "seed_rule": "SEED + 200000 + b",
           "refit": _levels_detail(c["delta_bs"], c["yardstick"], (r[:, 0] - r[:, 1]) / w, (r[:, 4] - r[:, 5]) / w, rng),
           "fixed": _levels_detail(c["delta_bs"], c["yardstick"], (r[:, 2] - r[:, 3]) / w, (r[:, 4] - r[:, 5]) / w, rng)}
    out["injected_or"] = _or_conversions({"pigeonhole_refit": out["refit"]["upper_R"]})
    out["rule"] = ("B-14 B4 (addendum 2): pigeonhole bootstrap (Owen 2007): multinomial weights of location groups and "
                   "of activations drawn independently, cell weight = product; stage 2 refitted on the dumps of B-8 C2 "
                   "(design and signs of B-13 B1(b)(c))")
    return out


# %% B4d（追記 4、B4 の結果を見た後の記述）: 場所グループと activation の重なりと、外れた再標本
def b4d_describe() -> dict:
    lex = v36.load(*CORR)
    df = pd.DataFrame({"u": np.asarray(lex["unit"]).astype(str), "a": _act(np.asarray(lex["aoi"]))})
    ga, gu = df.groupby("a")["u"].nunique(), df.groupby("u")["a"].nunique()
    a1, u1 = df["a"].map(ga).eq(1), df["u"].map(gu).eq(1)
    r = np.load(VER / "v73_b4_pigeonhole_reps.npy")
    e = (r[:, 0] - r[:, 1]) / r[:, 6]
    mad = 1.4826 * np.median(np.abs(e - np.median(e)))
    z = (e - np.median(e)) / mad
    keep = np.abs(z) <= 4
    return {"n_activations": len(ga), "share_activations_in_one_group": float((ga == 1).mean()),
            "n_groups": len(gu), "share_groups_with_one_activation": float((gu == 1).mean()),
            "n_pairs": int(df.groupby(["u", "a"]).ngroups), "share_cells_activation_in_one_group": float(a1.mean()),
            "share_cells_group_with_one_activation": float(u1.mean()), "share_cells_both": float((a1 & u1).mean()),
            "n_reps": len(e), "n_reps_robust_z_gt4": int((~keep).sum()), "max_abs_robust_z": float(np.abs(z).max()),
            "sd_delta_bs_refit": float(np.std(e, ddof=1)), "sd_delta_bs_refit_without": float(np.std(e[keep], ddof=1)),
            "rule": "B-14 addendum 4 (after the B4 result): overlap of location groups and activations (corrected run, "
                    "primary split) and the refitted B4 resamples with a robust z of delta BS above 4; descriptive"}


# %% B5: 基準の分解と組成を含まない比較
def _brier(d: dict, arms: tuple[str, ...]) -> dict:
    y = d["y"].astype(np.float64)
    return {a: float(np.mean((d[f"p_{a}"].astype(np.float64) - y) ** 2)) for a in arms if f"p_{a}" in d}


def _test_a(d: dict, base: str, model: str, rng: np.random.Generator) -> dict:
    n, s, units = v69._group_sums(d, ("M0", "M1", base, model))
    eg, dg = s[base] - s[model], s["M0"] - s["M1"]
    big_n = float(n.sum())
    w = _boot_w(rng, len(units))
    nw = w @ n
    return v70._summ(float(eg.sum() / big_n), float(dg.sum() / big_n), (w @ eg) / nw, (w @ dg) / nw)


def b5_decomposition() -> dict:
    rng = np.random.default_rng(SEED)
    arms = ("M0", "M1", "M1s", "M1c", "M3n", "M3d", "M3s")
    specs = (("corrected", "world_v2_lexfix", "primary", None), ("corrected", "world_v2_lexfix", "llmcut", None),
             ("confirmatory", "world_v2", "primary", "world_v2_dm"),
             ("confirmatory", "world_v2_actgrp", "primary", "world_v2_dm_actgrp"),
             ("confirmatory", "world_v2", "country", "world_v2_dm"),
             ("confirmatory", "world_v2", "temporal", "world_v2_dm"),
             ("confirmatory", "world_v2", "llmcut", "world_v2_dm"))
    out: dict = {}
    for code, run, split, dm in specs:
        d = v36.load(run, split)
        if d is None:
            continue
        note = None
        if "p_M3d" not in d and dm:
            x = v36.load(dm, split)
            if x is not None and np.array_equal(x["y"], d["y"]) and np.array_equal(x["aoi"], d["aoi"]):
                d["p_M3d"] = x["p_M3d"]
                note = f"M3d from {dm}"
        bs = _brier(d, arms)
        yard = bs["M0"] - bs["M1"]
        row: dict = {"code": code, "run": run, "split": split, "brier": bs, "yardstick": yard, "note": note}
        inc = {"composition_M1s_M1c": ("M1s", "M1c"), "names_given_composition_M1c_M3n": ("M1c", "M3n"),
               "within_aoi_M1s_M3d": ("M1s", "M3d"), "total_M1s_M3s": ("M1s", "M3s")}
        row["increments"] = {k: {"delta_bs": bs[a] - bs[b], "R": (bs[a] - bs[b]) / yard}
                             for k, (a, b) in inc.items() if a in bs and b in bs}
        for k, (a, b) in (("test_a_within_aoi_M1s_M3d", ("M1s", "M3d")), ("test_a_total_M1s_M3s", ("M1s", "M3s"))):
            if a in bs and b in bs:
                row[k] = _test_a(d, a, b, rng)
        out[f"{code}_{run}_{split}"] = row
    out["rule"] = ("B-14 B5 (addendum 2): Brier scores of the reference arms and increments over the yardstick; "
                   "corrected test (a) (4,000 location-group replicates) for M1s -> M3d and M1s -> M3s; corrected run's "
                   "splits as specified, confirmatory code's splits for reference")
    return out


# %% B6: 保たれる割合
RECORDS = (("gsw_from_2022", "world_v2_b13gsw22", ("M1cG", "M3oG"), "world_v2_b13gsw22", ("M1cW", "M3oW")),
           ("gsw_from_2022_nonpermanent", "world_v2_b13gsw22", ("M1cGn", "M3oGn"), "world_v2_b13gsw22",
            ("M1cW", "M3oW")),
           ("ems_earlier_events", "world_v2_b13recprior", ("M1cE", "M3oE"), "world_v2_b13recprior", ("M1cEc", "M3oEc")),
           ("jrc_100yr", "world_v2_b13jrcset", ("M1cJ", "M3oJ"), "world_v2_b13jrcset", ("M1cJc", "M3oJc")),
           ("ems_events_before_after", "world_v2_b7recur", ("M1cE", "M3oE"), "world_v2_b7recur", ("M1cEc", "M3oEc")),
           ("gsw_all_events", "world_v2_b6gsw", ("M1cG", "M3oG"), "world_v2_b8gswcell_local", ("M1cW", "M3oW")),
           ("gsw_all_events_nonpermanent", "world_v2_b6gsw", ("M1cGn", "M3oGn"), "world_v2_b8gswcell_local",
            ("M1cW", "M3oW")))


def b6_retained() -> dict:
    rng = np.random.default_rng(SEED)
    out: dict = {}
    cache: dict = {}
    for key, run_s, arms_s, run_c, arms_c in RECORDS:
        for r in (run_s, run_c):
            if r not in cache:
                cache[r] = v36.load(r, "primary")
        ds, dcl = cache[run_s], cache[run_c]
        if run_s != run_c and not (np.array_equal(ds["y"], dcl["y"]) and np.array_equal(ds["aoi"], dcl["aoi"])):
            out[key] = {"note": "rows differ between the two runs"}
            continue
        y = ds["y"].astype(np.float64)

        def gain(d: dict, a: tuple[str, str], y: np.ndarray = y) -> np.ndarray:
            return (d[f"p_{a[0]}"].astype(np.float64) - y) ** 2 - (d[f"p_{a[1]}"].astype(np.float64) - y) ** 2

        yard = (ds["p_M0"].astype(np.float64) - y) ** 2 - (ds["p_M1"].astype(np.float64) - y) ** 2
        units, _, es, ecl, dd = _gsum(np.asarray(ds["unit"]), gain(ds, arms_s), gain(dcl, arms_c), yard)
        w = _boot_w(rng, len(units))
        bs, bc = w @ es, w @ ecl
        pos = bc > 0
        share = bs[pos] / bc[pos]
        out[key] = {"runs": [run_s, run_c], "arms_settlement": list(arms_s), "arms_cell": list(arms_c),
                    "R_settlement": float(es.sum() / dd.sum()), "R_cell": float(ecl.sum() / dd.sum()),
                    "retained_share": float(es.sum() / ecl.sum()),
                    "share_pct90": np.quantile(share, [0.05, 0.95]).tolist() if pos.any() else None,
                    "share_reps_denominator_le0": float(1 - pos.mean()), "n_cells": len(y), "n_units": len(units)}
    vals = [v["retained_share"] for k, v in out.items() if "retained_share" in v and "nonpermanent" not in k]
    out["range_retained_share"] = [min(vals), max(vals)] if vals else None
    out["rule"] = ("B-14 B6 (addendum 2): gain of the record through settlements over its gain directly in the cell, "
                   "same location-group replicates (4,000); 90% percentile interval over replicates with a positive "
                   "denominator")
    return out


# %% B7: 名の特徴の held-out の係数
def _std_fit(y: np.ndarray, cols: list[tuple[str, np.ndarray]], unit: np.ndarray, n_names: int) -> dict:
    names = [c for c, _ in cols]
    x = np.column_stack([v for _, v in cols]).astype(np.float64)
    sd = x.std(axis=0)
    keep = sd > 0
    x = (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]
    kept = [nm for nm, k in zip(names, keep, strict=True) if k]
    res = v39.fit(y, x, np.zeros(len(y), dtype=np.int64), unit)
    k_n = int(keep[-n_names:].sum())
    b, se = res["beta"][-k_n:], res["se"][-k_n:]
    v = res["cov"][-k_n:, -k_n:]
    w = float(b @ np.linalg.solve(v, b))
    return {"n_cells": len(y), "n_units": int(res["n_groups"]),
            "coefficients": {nm: {"beta": float(bb), "ci95": [float(bb - 1.96 * s), float(bb + 1.96 * s)],
                                  "favours_flooding": bool(bb - 1.96 * s > 0), "against": bool(bb + 1.96 * s < 0)}
                             for nm, bb, s in zip(kept[-k_n:], b, se, strict=True)},
            "n_positive": int((b > 0).sum()), "n_negative": int((b < 0).sum()),
            "wald": w, "df": k_n, "p": float(stats.chi2.sf(w, k_n))}


def _logit(p: np.ndarray) -> np.ndarray:
    q = np.clip(p.astype(np.float64), 1e-6, 1 - 1e-6)
    return np.log(q / (1 - q))


def b7_direction() -> dict:
    folds = v69._dump_folds()
    names = [str(c) for c in folds[0]["name_cols"]]
    comp = [str(c) for c in folds[0]["comp_cols"]]
    cat = {k: np.concatenate([f[k] for f in folds]) for k in ("z_te", "xc_te", "xn_te", "y_te", "unit_te",
                                                                 "p_M1c_te")}
    y = cat["y_te"].astype(np.float64)
    unit = cat["unit_te"].astype(str)
    out: dict = {}
    cols = [("base_logit", cat["z_te"])] + [(c, cat["xc_te"][:, i]) for i, c in enumerate(comp)] + \
        [(c, cat["xn_te"][:, i]) for i, c in enumerate(names)]
    out["ems_b8_estimator"] = _std_fit(y, cols, unit, len(names))
    cols2 = [("logit_M1c", _logit(cat["p_M1c_te"]))] + [(c, cat["xn_te"][:, i]) for i, c in enumerate(names)]
    out["ems_m1c_base"] = _std_fit(y, cols2, unit, len(names))
    # GFM: 計画した GFM の表、基礎は GFM で学習し採点した実行の M1c
    v41 = v70._load("v41_gfm_independent_labels")
    t41 = v41.build().copy()
    if not (np.array_equal(t41["unit"].to_numpy(str), unit) and np.array_equal(t41["y"].to_numpy(int), y.astype(int))):
        out["gfm"] = {"note": "v41 table is not aligned with the dumps"}
        return out
    t41["_row"] = np.arange(len(t41))
    ll = pd.read_parquet(VER / "v37_table.parquet", columns=["aoi", "lon", "lat"])
    t41["lon"], t41["lat"] = ll["lon"].to_numpy(np.float32), ll["lat"].to_numpy(np.float32)  # 予測の座標は float32
    g = v41.attach_gfm(t41, _json(v69.GFM_SUMMARY))
    gp = v36.load("world_v2_gfm", "primary")
    gpd = pd.DataFrame({"aoi": np.asarray(gp["aoi"]).astype(str), "lon": gp["lon"].astype(np.float32),
                        "lat": gp["lat"].astype(np.float32),
                        "y_run": gp["y"].astype(int), "p_M1c_gfm": gp["p_M1c"].astype(np.float64)})
    m = g.merge(gpd, on=["aoi", "lon", "lat"], how="left", validate="one_to_one")
    ok = m["p_M1c_gfm"].notna().to_numpy()
    agree = float((m.loc[ok, "y_run"].to_numpy() == m.loc[ok, "y_gfm"].to_numpy()).mean())
    m = m.loc[ok].reset_index(drop=True)
    rows = m["_row"].to_numpy()
    cols3 = [("logit_M1c_gfm", _logit(m["p_M1c_gfm"].to_numpy()))] + \
        [(c, cat["xn_te"][rows, i]) for i, c in enumerate(names)]
    out["gfm_m1c_base"] = _std_fit(m["y_gfm"].to_numpy(np.float64), cols3, m["unit"].to_numpy(str), len(names)) | {
        "n_cells_table": len(g), "n_joined": int(ok.sum()), "label_agreement_with_run": agree}
    out["rule"] = ("B-14 B7 (addendum 2): logistic regression with a single intercept, standardised columns, "
                   "cluster-robust by location group (v39.fit); EMS on the held-out cells of the stage-2 dumps of B-8 "
                   "C2 (B-13 B8 estimator, and with the logit of the corrected run's M1c as the base); GFM on the "
                   "planned GFM table with the logit of the held-out M1c of world_v2_gfm as the base")
    return out


# %% B8: 地名学の表
_TOK = re.compile(r"[^a-z]+")
WIKIDATA_CACHE = VER / "v73_wikidata_named_after_river.json"
SPARQL = """SELECT DISTINCT ?gn WHERE {
  ?item wdt:P1566 ?gn ; wdt:P138 ?x .
  ?x wdt:P31/wdt:P279* wd:Q4022 .
}"""


def _analysed() -> pd.DataFrame:
    u = pd.read_parquet(VER / "v68_uncapped.parquet", columns=["aoi", "P_gid1"])
    u = u.loc[u["P_gid1"].notna()]
    u["P_gid1"] = u["P_gid1"].astype(np.int64)
    return u.loc[u["P_gid1"] >= 0]


def _gn() -> pd.DataFrame:
    return (pd.read_parquet(v69.GN, columns=["gid", "name", "asciiname", "lat", "lon", "fcode", "cc", "population"])
            .drop_duplicates("gid").set_index("gid"))


def _wikidata() -> dict:
    if WIKIDATA_CACHE.is_file():
        return _json(WIKIDATA_CACHE)
    cmd = ["curl", "-sS", "--max-time", "300", "-G", "https://query.wikidata.org/sparql",
           "--data-urlencode", f"query={SPARQL}", "-H", "Accept: application/sparql-results+json",
           "-H", "User-Agent: toponym-semantics-research/1.0 (verification v73)"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=330, check=False)
        js = json.loads(r.stdout)
        gids = sorted({int(b["gn"]["value"]) for b in js["results"]["bindings"] if b["gn"]["value"].isdigit()})
        res = {"ok": True, "n": len(gids), "gids": gids, "query": SPARQL,
               "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    except (json.JSONDecodeError, KeyError, subprocess.TimeoutExpired) as exc:
        res = {"ok": False, "error": f"{type(exc).__name__}: {str(exc)[:200]}", "query": SPARQL}
    WIKIDATA_CACHE.write_text(json.dumps(res), encoding="utf-8")
    return res


def _wilson(k: int, n: int) -> list[float] | None:
    return v70._wilson(k, n) if n else None


def _matches(lex: Lexicon, name: str, cc: str) -> list:
    toks = [x for x in _TOK.split(str(name).lower()) if x]
    return [e for e in lex.entries_for(cc) if any(_hit(tk, e) for tk in toks)]


def b8_toponymy() -> dict:
    lex, lexc = Lexicon(LEX), Lexicon(LEX_CORR)
    rec = _json(DOCV / "v70_lexicon_corrected_record.json")
    cat_of = {(e.lang, e.form, e.pos): e.cat for es in lex.by_lang.values() for e in es}
    rows = [{"change": "removed", "lang": d["lang"], "form": d["form"], "pos": d["pos"],
             "class": cat_of.get((d["lang"], d["form"], d["pos"]), "?"), "rule": d["rule"]} for d in rec["dropped"]]
    rows += [{"change": "added", "lang": d["lang"], "form": d["form"], "pos": d["pos"], "class": d["cat"],
              "rule": d["rule"]} for d in rec["added"]]
    ed = pd.DataFrame(rows).sort_values(["change", "lang", "pos", "form"])
    ed.to_csv(DOCV / "v73_lexicon_edits.csv", index=False)
    out: dict = {"i_edits": {k: {str(kk): int(vv) for kk, vv in ed.groupby(["change", k]).size().items()}
                             for k in ("rule", "lang", "pos", "class")} | {"n_removed": int((ed["change"] == "removed")
                                                                                            .sum()),
                                                                           "n_added": int((ed["change"] == "added")
                                                                                          .sum()),
                                                                           "file": "docs/verification/v73_lexicon_edits.csv",
                                                                           "countries": rec["countries"]}}
    out["ii_blocks"] = {"frozen": {cc: list(v) for cc, v in sorted(lex.country_langs.items())},
                        "corrected": {cc: list(v) for cc, v in sorted(lexc.country_langs.items())}}
    u = _analysed()
    gn = _gn()
    gids = np.unique(u["P_gid1"].to_numpy())
    s = gn.reindex(gids)
    s = s[s["asciiname"].notna()]
    words = s["asciiname"].str.lower().map(lambda x: [w for w in _TOK.split(x) if w])
    flags = {"was_suffix": words.map(lambda ws: any(w.endswith("was") and len(w) > 3 for w in ws)),
             "ham_suffix": words.map(lambda ws: any(w.endswith(("ham", "hamm")) and len(w) > 4
                                                    for w in ws)),
             "bruhl": words.map(lambda ws: any("bruhl" in w or "bruehl" in w for w in ws)),
             "terp_wierde": words.map(lambda ws: any(w.endswith(("terp", "wierde")) for w in ws)),
             "sur_word": words.map(lambda ws: "sur" in ws), "nad_word": words.map(lambda ws: "nad" in ws)}
    wd = _wikidata()
    if wd.get("ok"):
        flags["wikidata_named_after_river"] = pd.Series(s.index.isin(set(wd["gids"])), index=s.index)
    out["iii_omitted_elements"] = {
        "n_settlements": len(s), "wikidata": {k: wd.get(k) for k in ("ok", "n", "error", "retrieved_utc")},
        **{k: {"n": int(v.sum()), "by_country": {str(c): int(x) for c, x in s.loc[v.to_numpy(bool), "cc"]
                                                 .value_counts().head(15).items()}} for k, v in flags.items()}}
    forms = ("char", "au", "ab")
    iv: dict = {}
    for lab, lx in (("frozen", lex), ("corrected", lexc)):
        c: dict = {f: {} for f in forms}
        for r in s.itertuples():
            hit = {e.form for e in _matches(lx, r.asciiname, r.cc) if e.form in forms}
            for f in hit:
                c[f][str(r.cc)] = c[f].get(str(r.cc), 0) + 1
        iv[lab] = {f: {"n": sum(v.values()), "by_country": dict(sorted(v.items(), key=lambda kv: -kv[1]))}
                   for f, v in c.items()}
    out["iv_char_au_ab"] = iv
    tc = pd.read_csv(DOCV / "v70_two_coders.csv")
    tc["block"] = tc["cc"].map(lambda c: "+".join(lex.country_langs.get(str(c), ("none",))))
    tc["forms"] = [sorted({e.form for e in _matches(lex, n, c) if e.cat in WET_SIDE})
                   for n, c in zip(tc["name"].astype(str), tc["cc"].astype(str), strict=True)]
    v: dict = {"n_items": len(tc), "by_country_stratum": {f"{c}|{st}": int(n) for (c, st), n in
                                                          tc.groupby(["cc", "stratum"]).size().items()}}
    for who in ("coder1", "coder2"):
        d = tc[tc[f"status_{who}"] == "coded"].copy()
        d["truth"] = d[f"codes_{who}"].fillna("").map(lambda x: bool(set(str(x).split(";")) & set(WET_SIDE)))

        def rate(x: pd.DataFrame) -> dict:
            k, n = int(x["truth"].sum()), len(x)
            return {"n": n, "k": k, "rate": k / n if n else None, "ci95": _wilson(k, n)}

        wet, nott = d[d["stratum"] == "wet"], d[d["stratum"] == "not"]
        fx = wet.explode("forms")
        top = fx["forms"].value_counts().head(20).index
        v[who] = {"precision_all": rate(wet), "miss_all": rate(nott),
                  "precision_by_block": {b: rate(x) for b, x in wet.groupby("block")},
                  "miss_by_block": {b: rate(x) for b, x in nott.groupby("block")},
                  "precision_by_form_top20": {f: rate(fx[fx["forms"] == f]) for f in top}}
    out["v_audit_sample"] = v
    u2 = u.assign(cc=gn["cc"].reindex(u["P_gid1"].to_numpy()).to_numpy())
    area = u2.groupby("P_gid1").size() * CELL_KM2
    per_aoi = u2.groupby("aoi").agg(n=("P_gid1", "nunique"), cc=("cc", lambda x: x.value_counts().index[0]))
    lat60, lon60 = s["lat"] * 60, s["lon"] * 60
    v37 = v70._load("v37_positive_control")
    rounded = ((lat60 - lat60.round()).abs() < v37.MINUTE_TOL) & ((lon60 - lon60.round()).abs() < v37.MINUTE_TOL)
    s2 = s.assign(area=area.reindex(s.index).to_numpy(), rounded=rounded.to_numpy())
    vi: dict = {}
    for c, x in s2.groupby("cc"):
        fc = x["fcode"].value_counts(normalize=True)
        vi[str(c)] = {"n_settlements": len(x), "fcode_shares": {str(k): float(q) for k, q in fc.head(6).items()},
                      "median_settlements_per_aoi": float(per_aoi.loc[per_aoi["cc"] == c, "n"].median())
                      if (per_aoi["cc"] == c).any() else None,
                      "n_aoi": int((per_aoi["cc"] == c).sum()),
                      "median_nearest_ground_km2": float(x["area"].median()),
                      "share_whole_minutes": float(x["rounded"].mean())}
    out["vi_gazetteer_by_country"] = vi
    out["rule"] = ("B-14 B8 (addendum 2): lexicon edits with the class of the original entry; blocks; settlements nearest "
                   "to an evaluable cell of the uncapped table; omitted elements by word; char, au and ab matches by "
                   "country; precision and misses of the published two coders by block and form (Wilson 95%); "
                   "gazetteer by country")
    return out


# %% B9: 修飾語の対の標高と川の大きさ
def _pair_stats(x: np.ndarray, units: np.ndarray, rng: np.random.Generator) -> dict:
    neg, pos = int((x < 0).sum()), int((x > 0).sum())
    uu, inv = np.unique(units, return_inverse=True)
    med, share = [], []
    for _ in range(B_PAIR_BOOT):
        w = rng.multinomial(len(uu), np.full(len(uu), 1 / len(uu)))[inv]
        xs = np.repeat(x, w)
        if len(xs):
            med.append(float(np.median(xs)))
            share.append(float((xs < 0).mean()))
    return {"n": len(x), "median": float(np.median(x)), "n_negative": neg, "n_positive": pos,
            "share_lower_member_lower": neg / max(1, neg + pos), "n_units": len(uu),
            "sign_test_p": float(stats.binomtest(neg, neg + pos).pvalue) if neg + pos else None,
            "median_ci90": np.quantile(med, [0.05, 0.95]).tolist(), "share_ci90": np.quantile(share, [0.05, 0.95])
            .tolist()}


def b9_qualifiers() -> dict:
    rng = np.random.default_rng(SEED)
    p = pd.read_csv(VER / "v70_b15_pairs.csv")
    sr = pd.read_parquet(VER / "v57_settlement_rivers.parquet", columns=["gid", "elev_site", "upland_nearest"])
    by = sr.groupby("gid")[["elev_site", "upland_nearest"]].median()
    lo, hi = by.reindex(p["lower"].to_numpy()), by.reindex(p["upper"].to_numpy())
    p["elev_diff"] = lo["elev_site"].to_numpy() - hi["elev_site"].to_numpy()
    p["log_upland_ratio"] = np.log(lo["upland_nearest"].to_numpy()) - np.log(hi["upland_nearest"].to_numpy())
    out: dict = {"n_pairs": len(p)}
    for col in ("elev_diff", "log_upland_ratio"):
        q = p[np.isfinite(p[col].to_numpy(np.float64))]
        out[col] = _pair_stats(q[col].to_numpy(np.float64), q["unit"].to_numpy(str), rng)
    out["log_upland_ratio"]["median_ratio"] = float(np.exp(out["log_upland_ratio"]["median"]))
    out["rule"] = ("B-14 B9 (addendum 2): pairs of B-13 B15; elevation of the cell holding the settlement point (GLO-30 "
                   "mean) and drainage area of the nearest reach (v57); lower minus upper (log ratio for the area); "
                   "sign test and location-group bootstrap (2,000)")
    return out


# %% B10: 目録の向き
WET_CLASSES = ("wetside", "water", "wet", "riverine", "crossing", "low", "floodplain", "any", "lex")
WET_VARS = ("low_hand", "low_half", "occ", "seas", "flood", "ems", "gfm", "recur", "y")
DRY_VARS = ("hand_mean", "hand_pct", "hand_min", "log_river", "elev", "h_above")


def _class_of(path: str) -> str | None:
    segs = re.split(r"[/:]", path)
    if "high" in segs:
        return "high"
    for c in WET_CLASSES:
        if c in segs:
            return "wet"
    return None


def _direction(row: pd.Series) -> str:
    if row["family"] == "nuisance":
        return "not_applicable"
    k, est, path = row["kind"], float(row["estimate"]), str(row["path"])
    null = 1.0 if k == "odds_ratio" else 0.0
    up = est > null
    if row["family"] == "flood" and k == "R":
        return "favours" if up else "against"
    cls = _class_of(path)
    if cls is None:
        return "unclassified"
    sign = 1 if cls == "wet" else -1
    if row["family"] == "flood":
        return "favours" if (up if sign > 0 else not up) else "against"
    last = path.split("/")[-1]
    if any(last.startswith(v) for v in DRY_VARS):
        want_up = sign < 0
    elif any(last.startswith(v) for v in WET_VARS):
        want_up = sign > 0
    else:
        return "unclassified"
    return "favours" if up == want_up else "against"


def b10_census() -> dict:
    df = pd.read_csv(DOCV / "v69_estimates.csv")
    b = df.loc[~df["duplicate"] & df["valid"] & df["name_effect"] & df["excludes_null"]].copy()
    b["direction"] = b.apply(_direction, axis=1)
    out: dict = {}
    for fam, x in b.groupby("family"):
        vc = x["direction"].value_counts()
        out[str(fam)] = {"n_excluding_null": len(x), **{str(k): int(v) for k, v in vc.items()},
                         "unclassified_paths": [f"{r.file}:{r.path}" for r in x[x["direction"] == "unclassified"]
                                                .itertuples()][:40]}
    hl = df[(df["file"] == "v41_gfm_independent_labels.json")
            & (df["path"] == "sensitivity_independent_sensor/odds_gfm/categories_adjusted_M1/terms/high:le500m")]
    out["highland_gfm_not_s1"] = hl[["file", "path", "estimate", "ci95_lo", "ci95_hi", "family", "near_settlement",
                                     "duplicate", "valid", "p_approx", "q_bh", "q_by"]].to_dict("records")
    out["rule"] = ("B-14 B10 (addendum 2): final census of B-12 (v69_estimates.csv), unique valid name-effect rows that "
                   "exclude the null; direction by family, class and quantity as in addendum 2")
    return out


# %% B11・B12
def b11_pairs_084() -> dict:
    x = np.load(VER / "v68_pairs_v71_c2_uncapped.npy")
    x = x[np.isfinite(x)]
    c2 = _json(DOCV / "v71_b13_cluster.json")["C2_near"]["ems_uncapped"]
    return {"n": len(x), "point_or": c2["odds_ratio"], "median_or": float(np.exp(np.median(x))),
            "mean_log_or": float(np.mean(x)), "sd_log_or": float(np.std(x, ddof=1)),
            "skewness_log_or": float(stats.skew(x)), "share_or_above_1": float((x > 0).mean()),
            "q025_q975_or": np.exp(np.quantile(x, [0.025, 0.975])).tolist(),
            "point_minus_median_log": float(math.log(c2["odds_ratio"]) - np.median(x)),
            "rule": "B-14 B11: the 2,000 pairs replicates (log odds ratio) of B-13 C2, uncapped EMS"}


def b12_margin(b1: dict) -> dict:
    k = _json(DOCV / "v69_b12_local.json")["B17"]["primary"]
    clim = k.get("brier_training_constant", k["brier_best_constant"])
    lv = b1["levels"]["a_fixed_groups"]
    d0 = lv["yardstick"]
    out = {"climatology_brier": clim, "yardstick": d0, "margin_brier": SESOI_R * d0,
           "margin_bss": SESOI_R * d0 / clim, "yardstick_bss": d0 / clim, "upper": {}}
    for key, v in b1["levels"].items():
        if v.get("upper_R") is not None:
            out["upper"][key] = {"upper_R": v["upper_R"], "upper_brier": v["upper_R"] * d0,
                                 "upper_bss": v["upper_R"] * d0 / clim}
    out["rule"] = "B-14 B12: Brier skill increments relative to the training-side constant of B-12 B17 (primary split)"
    return out


# %% B13: 集落の近くの推定の点検
def _gfm_jackknife() -> dict:
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    y = g["y_gfm"].to_numpy(np.float64)
    rate = pd.Series(y).groupby(g["aoi"].to_numpy()).transform("mean").to_numpy()
    s = g.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
    ys = s["y_gfm"].to_numpy(np.float64)
    z, names = v39.design(s, ("wetside",), True, banded=True)
    j = names.index("wetside:le500m")
    unit, aoi = s["unit"].to_numpy(str), s["aoi"].to_numpy()
    full = v39.fit(ys, z, aoi, unit)
    b0 = float(full["beta"][j])
    uu = np.unique(unit)
    bs = []
    for u in uu:
        m = unit != u
        bs.append(float(v39.fit(ys[m], z[m], aoi[m], unit[m])["beta"][j]))
    bs = np.array(bs)
    se = math.sqrt((len(uu) - 1) / len(uu) * np.sum((bs - bs.mean()) ** 2))
    return {"odds_ratio": math.exp(b0), "se_cluster_robust": float(full["se"][j]), "se_jackknife": se,
            "ci95_jackknife": [math.exp(b0 - 1.96 * se), math.exp(b0 + 1.96 * se)],
            "ci95_cluster_robust": [math.exp(b0 - 1.96 * full["se"][j]), math.exp(b0 + 1.96 * full["se"][j])],
            "n_units": len(uu), "range_or_leave_one_out": [float(np.exp(bs.min())), float(np.exp(bs.max()))],
            "most_influential": [{"unit": str(uu[i]), "or_without": float(np.exp(bs[i]))}
                                 for i in np.argsort(-np.abs(bs - b0))[:5]]}


def b13_near_checks() -> dict:
    v55 = v70._load("v55_b7_local")
    t = v55.table()
    _, infos, _ = v55._infos()
    untyped = t["aoi"].map(lambda a: bool(infos[a].unspecified_as_riverine)).to_numpy()
    out: dict = {"main_capped": v55.odds(t)}
    riv = t["f_riverine"].to_numpy() + np.where(untyped, t["f_other"].to_numpy(), 0.0)
    out["i_fractional_cover"] = v55.odds(t.assign(y_frac=np.clip(riv, 0.0, 1.0)), "y_frac")
    out["ii_typed_only"] = v55.odds(t, mask=~untyped) | {"n_aoi_untyped_dropped": int(t.loc[untyped, "aoi"].nunique())}
    gn = _gn()
    pop = gn["population"].reindex(t["P_gid1"].to_numpy(np.int64)).to_numpy(np.float64)
    built = v69._lc_built(t)
    tt = t.assign(log_pop=np.log1p(pop), lc_built=built)
    ok = np.isfinite(tt["log_pop"].to_numpy()) & np.isfinite(built)
    out["iii_population_builtup"] = {"adjusted": v55.odds(tt, extra=("log_pop", "lc_built"), mask=ok),
                                     "same_rows_unadjusted": v55.odds(tt, mask=ok), "n_rows_dropped": int((~ok).sum())}
    out["iv_gfm_jackknife"] = _gfm_jackknife()
    out["rule"] = ("B-14 B13 (addendum 2): capped table and estimator of v55 (frozen lexicon reading); (i) riverine cover "
                   "as a fractional response; (ii) typed AOIs only; (iii) log(1 + population) of the nearest settlement "
                   "and the 10 m built-up share; (iv) planned GFM estimate with a leave-one-location-group-out "
                   "jackknife")
    return out


# %% B14: 位置で結ぶ表の整合
def b14_join_check() -> dict:
    from pyproj import CRS, Transformer

    v55 = v70._load("v55_b7_local")
    _, _, dirs = v55._infos()
    gn = _gn()
    rows = []
    for d in dirs:
        c = pq.read_table(d / "cells.parquet", columns=["x", "y"]).to_pandas()
        f = pq.read_table(d / "features.parquet", columns=["P_gid1", "P_d1"]).to_pandas()
        n_lc = pq.read_metadata(d / "landcover.parquet").num_rows
        rec = {"aoi": d.name, "n_cells": len(c), "n_features": len(f), "n_landcover": n_lc}
        if len(c) == len(f):
            meta = _json(d / "meta.json")
            gid = f["P_gid1"].to_numpy(np.int64)
            has = gid >= 0
            lat = gn["lat"].reindex(gid[has]).to_numpy(np.float64)
            lon = gn["lon"].reindex(gid[has]).to_numpy(np.float64)
            tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(meta["crs_proj4"]), always_xy=True)
            px, py = tr.transform(lon, lat)
            dist = np.hypot(np.asarray(px) - c["x"].to_numpy()[has], np.asarray(py) - c["y"].to_numpy()[has])
            diff = np.abs(dist - f["P_d1"].to_numpy(np.float64)[has])
            rec |= {"n_with_settlement": int(has.sum()), "n_gid_missing": int(np.isnan(lat).sum()),
                    "n_mismatch_gt1m": int(np.nansum(diff > 1.0)), "max_abs_diff_m": float(np.nanmax(diff))
                    if has.any() else 0.0}
        rows.append(rec)
    r = pd.DataFrame(rows)
    bad_n = r[(r["n_cells"] != r["n_features"]) | (r["n_cells"] != r["n_landcover"])]
    bad_d = r[r.get("n_mismatch_gt1m", pd.Series(0, index=r.index)).fillna(0) > 0]
    return {"n_aoi": len(r), "n_aoi_row_counts_differ": len(bad_n), "aoi_row_counts_differ": bad_n["aoi"].tolist(),
            "n_aoi_distance_mismatch": len(bad_d), "aoi_distance_mismatch": bad_d["aoi"].tolist()[:50],
            "total_rows_checked": int(r["n_with_settlement"].sum()), "total_gid_missing": int(r["n_gid_missing"].sum()),
            "max_abs_diff_m": float(r["max_abs_diff_m"].max()),
            "rule": "B-14 B14 (addendum 2): row counts of cells, features and landcover; P_d1 recomputed from the cell "
                    "centre and the GeoNames coordinates of P_gid1 in the AOI's LAEA (agreement within 1 m)"}


# %% B15: AOI ごとの分布
def b15_per_aoi() -> dict:
    v55 = v70._load("v55_b7_local")
    lex = v36.load(*CORR)
    t = v55.table()
    aoi = np.asarray(lex["aoi"]).astype(str)
    if not np.array_equal(aoi, t["aoi"].to_numpy(str)):
        return {"note": "corrected predictions are not aligned with the v55 table"}
    _, infos, _ = v55._infos()
    y = lex["y"].astype(np.int64)
    rows = []
    cont = t.groupby("aoi")["continent"].agg(lambda s: s.value_counts().index[0])
    for a in np.unique(aoi):
        m = aoi == a
        ya = y[m]
        pc, pn = lex["p_M1c"][m].astype(np.float64), lex["p_M3n"][m].astype(np.float64)
        dbs = float(np.mean((pc - ya) ** 2) - np.mean((pn - ya) ** 2))
        dauc = v70._auc(ya, pn) - v70._auc(ya, pc) if 0 < ya.mean() < 1 else np.nan
        rows.append({"aoi": a, "n": int(m.sum()), "delta_bs": dbs, "delta_auc": dauc, "continent": cont[a],
                     "flood_type": "untyped" if infos[a].unspecified_as_riverine else "typed"})
    r = pd.DataFrame(rows)
    r.to_csv(VER / "v73_b15_per_aoi.csv", index=False)

    def summ(x: pd.DataFrame) -> dict:
        o: dict = {"n_aoi": len(x)}
        for c in ("delta_bs", "delta_auc"):
            v = x[c].dropna().to_numpy()
            o[c] = {"n": len(v), "median": float(np.median(v)) if len(v) else None,
                    "q25_q75": np.quantile(v, [0.25, 0.75]).tolist() if len(v) else None,
                    "share_positive": float((v > 0).mean()) if len(v) else None}
        return o

    b2 = (_json(OUT).get("B2") or {}) if OUT.is_file() else {}
    return {"all": summ(r), "by_continent": {str(k): summ(x) for k, x in r.groupby("continent")},
            "by_flood_type": {str(k): summ(x) for k, x in r.groupby("flood_type")},
            "leave_one_activation_top5": ((b2.get("jackknife") or {}).get("activations") or {}).get("top5"),
            "rule": "B-14 B15 (addendum 2): per AOI of the corrected run, delta BS = BS(M1c) - BS(M3n) and delta AUC = "
                    "AUC(M3n) - AUC(M1c); R without each of the five most influential activations from B2(i)"}


# %% 実行
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    want = {x for x in args.only.split(",") if x}

    def run(key: str, fn: callable, *a: object) -> None:
        if want and key not in want:
            return
        t0 = time.time()
        val = fn(*a)
        res = _json(OUT) if OUT.is_file() else {}
        res |= {"rule": "amendment B-14 (post hoc), local computations; confirmatory verdict unchanged",
                "seed": SEED, key: val}
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        print(f"{key} ({time.time() - t0:.0f} s) {json.dumps(val, default=float)[:700]}", flush=True)

    run("B1", b1_table)
    run("B2", b2_heavy_tails)
    run("B5", b5_decomposition)
    run("B6", b6_retained)
    run("B9", b9_qualifiers)
    run("B10", b10_census)
    run("B11", b11_pairs_084)
    if OUT.is_file() and "B1" in _json(OUT):
        run("B12", b12_margin, _json(OUT)["B1"])
    run("B15", b15_per_aoi)
    run("B7", b7_direction)
    run("B13", b13_near_checks)
    run("B14", b14_join_check)
    run("B8", b8_toponymy)
    run("B3", b3_operating)
    run("B4", b4_pigeonhole)
    if (VER / "v73_b4_pigeonhole_reps.npy").is_file():
        run("B4d", b4d_describe)


if __name__ == "__main__":
    main()
