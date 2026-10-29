# %% v70b（事後・amendment B-13 の B6・B8〜B12・B14・B15）: 既存の成果物から手元で計算するもの（第二部）
"""確証判定は変えない。B-13 の追記 4 で計算の前に固定した量。第一部（v70_b13_local.py）と同時に走るので、出力は
別のファイル docs/verification/v70_b13_local2.json（部分ごとに足す。`--only B6,B8` で一部だけ）。

- B6: 主の R の場所グループの影響（jackknife）、実効のクラスタ数、影響の大きい 10 を除いた値（固定と当て直し）
- B8: 名の 14 の係数の同時 Wald 検定（二段目の M3n を全てのセルで、場所グループでクラスタ頑健）
- B9: 集落の近くの点検（座標の丸め・activation を一つずつ・世代・置換・2×2・センサー・地物の符号）
- B10: 辞書の点検（同形異義・ブロック・照合の位置・Unicode の名）
- B11: N-read の区分と unknown の割合（ブロックと国ごと）
- B12: 二つの符号化者の一致と、公開の符号のファイル（照合の標本の名は除く）
- B14: GFM での当て直しの崩れの診断（AOI ごとの浸水の割合と基礎の利得）
- B15: 相対の修飾語の対
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import sys
import time
import unicodedata
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import stage2_weighted
from topsem.world.b13 import (
    PREFIX_QUALIFIERS,
    QUALIFIERS,
    match_pairs,
    qualifier_side,
    split_qualifier,
)
from topsem.world.lexicon import WET_SIDE, Lexicon, _hit

OUT = ROOT / "docs" / "verification" / "v70_b13_local2.json"
TWO_CODERS = ROOT / "docs" / "verification" / "v70_two_coders.csv"
DOCV = ROOT / "docs" / "verification"
VER = DATA_DIR / "verification"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
SEED = 20261102
B_PERM = 200
NEAR_M = 500.0
SITE_M = 177.0
GEN_EXCLUDE = (0, 1)
FCODE_EXCLUDE = ("PPLX", "PPLH", "PPLQ", "PPLW")
PAIR_MAX_M = 10_000.0
B_PAIR_BOOT = 2000
N_WORKERS = 3
AUDIT_SAMPLES = ("v51_audit_sample.csv", "v68_audit_sample.csv", "v69_audit_sample.csv")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v69 = _load("v69_b12_local")
v36, v68, v39 = v69.v36, v69.v68, v69.v39
v37 = _load("v37_positive_control")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    v69._threads1()


_C: dict = {}


def _codes(a: np.ndarray) -> np.ndarray:
    return np.unique(a, return_inverse=True)[1]


def _mode(a: np.ndarray) -> str:
    v, c = np.unique(a, return_counts=True)
    return str(v[np.argmax(c)])


def _gn() -> pd.DataFrame:
    if "gn" not in _C:
        _C["gn"] = (pd.read_parquet(v69.GN, columns=["gid", "name", "asciiname", "lat", "lon", "fcode", "cc"])
                    .drop_duplicates("gid").set_index("gid"))
    return _C["gn"]


def _uncapped() -> pd.DataFrame:
    if "u" not in _C:
        u = v68.uncapped_full()
        u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
        _C["u"] = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    return _C["u"]


def _planned_gfm() -> pd.DataFrame:
    """計画した GFM の推定の表（v41 の表に GFM のラベル。最寄りの集落は上限なしの表から (aoi, ix, iy) で結ぶ）。"""
    if "g" not in _C:
        v41 = _load("v41_gfm_independent_labels")
        g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
        u = v68.uncapped_full()[["aoi", "ix", "iy", "P_gid1"]].dropna()
        u = u.astype({"ix": np.int64, "iy": np.int64})
        g = g.astype({"ix": np.int64, "iy": np.int64}).merge(u, on=["aoi", "ix", "iy"], how="left",
                                                              validate="many_to_one")
        g["P_gid1"] = g["P_gid1"].fillna(-1).astype(np.int64)
        _C["g"] = g
    return _C["g"]


def _odds_planned(g: pd.DataFrame) -> dict:
    v41 = _load("v41_gfm_independent_labels")
    return v41.odds(g, "y_gfm")["wetside_adjusted_M1"]["terms"]["wetside:le500m"]


def _rounded(gids: np.ndarray) -> np.ndarray:
    """最寄りの集落の座標が整数の分か（v37 の規則: 緯度と経度の 60 倍が整数から 2e-3 以内）。"""
    gn = _gn()
    lat = gn["lat"].reindex(gids).to_numpy(np.float64) * 60
    lon = gn["lon"].reindex(gids).to_numpy(np.float64) * 60
    return (np.abs(lat - np.round(lat)) < v37.MINUTE_TOL) & (np.abs(lon - np.round(lon)) < v37.MINUTE_TOL)


def _variation(t: pd.DataFrame, ycol: str = "y") -> np.ndarray:
    y = t[ycol].to_numpy(dtype=np.float64)
    ia = _codes(t["aoi"].to_numpy())
    r = (np.bincount(ia, weights=y) / np.bincount(ia))[ia]
    return (r > 0) & (r < 1)


# %% B6: 影響
_B: dict = {}


def _b6_init() -> None:
    v68._b7_init()
    _B.update(v68._B7)


def _b6_rep(g: int) -> list[float]:
    acc = np.zeros(7)
    for f in _B["folds"]:
        w_tr = (f["u_tr"] != g).astype(np.float64)
        w_te = (f["u_te"] != g).astype(np.float64)
        pc, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xc_tr"], f["xc_te"], f["y_tr"], w_tr=w_tr, init=f["init_c"])
        pn, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xn_tr"], f["xn_te"], f["y_tr"], w_tr=w_tr, init=f["init_n"])
        y = f["y_te"].astype(np.float64)
        for i, p in enumerate((pc, pn, f["p_M1c_te"], f["p_M3n_te"], f["p_M0_te"], f["p_M1_te"])):
            acc[i] += float((w_te * (p - y) ** 2).sum())
        acc[6] += float(w_te.sum())
    return acc.tolist()


def b6_influence() -> dict:
    out: dict = {}
    for label, run in (("confirmatory", "world_v2"), ("corrected", "world_v2_lexfix")):
        d = v36.load(run, "primary")
        n, s, units = v69._group_sums(d, ("M0", "M1", "M1c", "M3n"))
        eg, dg = s["M1c"] - s["M3n"], s["M0"] - s["M1"]
        e, dd = float(eg.sum()), float(dg.sum())
        r_minus = (e - eg) / (dd - dg)
        psi = e / dd - r_minus
        top = np.argsort(-np.abs(psi))[:10]
        cc, uu = np.asarray(d["cc"]).astype(str), np.asarray(d["unit"]).astype(str)
        out[label] = {"R": e / dd, "n_units": len(units), "n_eff": float((psi**2).sum() ** 2 / (psi**4).sum()),
                      "top10": [{"unit": str(units[i]), "unit_index": int(i), "country": _mode(cc[uu == units[i]]),
                                 "R_without": float(r_minus[i]), "psi": float(psi[i]),
                                 "share_cells": float(n[i] / n.sum())} for i in top],
                      "range_R_without_top10": [float(r_minus[top].min()), float(r_minus[top].max())]}
    idx = [x["unit_index"] for x in out["corrected"]["top10"]]
    _threads1()
    with ProcessPoolExecutor(2, initializer=_b6_init) as ex:
        res = np.array(list(ex.map(_b6_rep, idx)), dtype=np.float64)
    rr = (res[:, 0] - res[:, 1]) / (res[:, 4] - res[:, 5])
    out["refit_corrected"] = {"R_without": [{"unit": out["corrected"]["top10"][k]["unit"], "R_without": float(rr[k])}
                                            for k in range(len(rr))], "range": [float(rr.min()), float(rr.max())]}
    out["rule"] = ("B-13 B6: jackknife over location groups (fixed models); n_eff = (sum psi^2)^2 / sum psi^4; the "
                   "10 most influential groups of the corrected run removed one at a time with the second stage "
                   "refitted on the design of B-8 C2")
    return out


# %% B8: 名の係数の同時検定
def b8_joint_wald() -> dict:
    folds = v69._dump_folds()
    z = np.concatenate([f["z_te"] for f in folds]).astype(np.float64)
    xc = np.concatenate([f["xc_te"] for f in folds]).astype(np.float64)
    xn = np.concatenate([f["xn_te"] for f in folds]).astype(np.float64)
    y = np.concatenate([f["y_te"] for f in folds]).astype(np.float64)
    unit = np.concatenate([f["unit_te"] for f in folds]).astype(str)
    names = [str(c) for c in folds[0]["name_cols"]]
    x = np.column_stack([z, xc, xn])
    sd = x.std(axis=0)
    keep = sd > 0
    x = (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]
    kn = keep[1 + xc.shape[1]:]
    k_n = int(kn.sum())
    res = v39.fit(y, x, np.zeros(len(y), dtype=np.int64), unit)
    b, v = res["beta"][-k_n:], res["cov"][-k_n:, -k_n:]
    w = float(b @ np.linalg.solve(v, b))
    v55j = _json(DOCV / "v55_b7_local.json")
    return {"n_cells": len(y), "n_units": int(res["n_groups"]), "n_name_coefs": k_n,
            "name_cols_kept": [nm for nm, k in zip(names, kn, strict=True) if k],
            "wald": w, "df": k_n, "p": float(stats.chi2.sf(w, k_n)),
            "v55_heldout_coef_test": v55j.get("heldout_coef_test"), "v55_wald_tests": v55j.get("wald_tests"),
            "rule": "B-13 B8: logistic regression of y on standardised [base logit, composition, name features] over "
                    "the held-out cells of all folds (stage-2 dumps of B-8 C2), cluster-robust by location group "
                    "(v39.fit with a single intercept); joint Wald test of the name coefficients (constant columns "
                    "dropped)"}


# %% B9: 集落の近くの点検
def _loo_init(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, act: np.ndarray) -> None:
    v68._pairs_init_w(y, z, aoi, np.zeros(len(y), dtype=np.int64), 1)
    _C["act"] = act


def _loo_rep(a: str) -> np.ndarray | None:
    w = v68._WW
    idx = np.flatnonzero(_C["act"] != a)
    try:
        return v68.fit_w(w["y"][idx], w["z"][idx], w["aoi"][idx], np.ones(len(idx)))
    except np.linalg.LinAlgError:
        return None


def b9_near() -> dict:
    u = _uncapped()
    g = _planned_gfm()
    out: dict = {"n_cells_uncapped": len(u), "n_cells_gfm": len(g)}
    rd, grd = _rounded(u["P_gid1"].to_numpy(np.int64)), _rounded(g["P_gid1"].to_numpy(np.int64))
    out["a_unrounded"] = {"share_cells_rounded": float(rd.mean()),
                          "ems_uncapped": v68.pairs_or(u.loc[~rd].reset_index(drop=True), "v70_b9a_unrounded"),
                          "gfm_planned": _odds_planned(g.loc[~grd].reset_index(drop=True))}
    print(f"  B9(a) {json.dumps(out['a_unrounded'], default=float)[:300]}", flush=True)
    act = np.array([a.split("_")[0] for a in u["aoi"].to_numpy(dtype=str)])
    near = (u["P_d1"].to_numpy() <= NEAR_M) & (u["lex_wetside_n1"].to_numpy() > 0)
    acts = sorted(set(act[near]))
    y, z, names, aoi, _, _, _ = v68._pairs_arrays(u)
    j = names.index("wetside:le500m")
    act_k = act[_variation(u)]
    _threads1()
    with ProcessPoolExecutor(N_WORKERS, initializer=_loo_init, initargs=(y, z, aoi, act_k)) as ex:
        fs = list(ex.map(_loo_rep, acts, chunksize=1))
    loo = {a: (float(np.exp(b[j])) if b is not None else None) for a, b in zip(acts, fs, strict=True)}
    vals = [x for x in loo.values() if x is not None]
    out["b_leave_one_activation"] = {"n_activations": len(acts), "min": min(vals), "max": max(vals),
                                     "by_activation": loo}
    keep2 = ~np.isin(act, ("EMSR052", "EMSR130"))
    out["b_without_emsr052_130"] = v68.pairs_or(u.loc[keep2].reset_index(drop=True), "v70_b9b_no052_130")
    gen = {a: int(_json(v69.LABELS / a / "meta.json").get("schema_generation", -1)) for a in u["aoi"].unique()}
    ug, gg = u["aoi"].map(gen).to_numpy(), g["aoi"].map(gen).to_numpy()
    vals_g, cnt_g = np.unique(list(gen.values()), return_counts=True)
    out["c_generation"] = {"n_aoi_by_generation": {str(k): int(c) for k, c in zip(vals_g, cnt_g, strict=True)},
                           "ems_uncapped": v68.pairs_or(u.loc[~np.isin(ug, GEN_EXCLUDE)].reset_index(drop=True),
                                                        "v70_b9c_gen2plus"),
                           "gfm_planned": _odds_planned(g.loc[~np.isin(gg, GEN_EXCLUDE)].reset_index(drop=True))}
    fc = _gn()["fcode"].reindex(u["P_gid1"].to_numpy(np.int64)).astype(str).to_numpy()
    gfc = _gn()["fcode"].reindex(g["P_gid1"].to_numpy(np.int64)).astype(str).to_numpy()
    out["g_fcodes"] = {"share_cells_excluded": float(np.isin(fc, FCODE_EXCLUDE).mean()),
                       "ems_uncapped": v68.pairs_or(u.loc[~np.isin(fc, FCODE_EXCLUDE)].reset_index(drop=True),
                                                    "v70_b9g_fcodes"),
                       "gfm_planned": _odds_planned(g.loc[~np.isin(gfc, FCODE_EXCLUDE)].reset_index(drop=True))}
    out["rule"] = ("B-13 B9 (a)(b)(c)(g): Sect. 5 estimator on the uncapped table (corrected readings; pairs cluster "
                   "bootstrap of B-11 B2) and the planned GFM estimate (v41, cluster-robust); (b) points by "
                   "frequency-weighted IRLS without each activation that has water-related cells within 500 m")
    return out


def _perm_index(t: pd.DataFrame) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    """(AOI, 最寄りの集落) の組の番号、AOI ごとの組の番号の一覧、組ごとの水系側の印。"""
    ia = _codes(t["aoi"].to_numpy())
    ig = _codes(t["P_gid1"].to_numpy(np.int64))
    key = _codes(ia.astype(np.int64) * (int(ig.max()) + 1) + ig)
    first = np.full(int(key.max()) + 1, -1)
    first[key[::-1]] = np.arange(len(key))[::-1]
    s_aoi = ia[first]
    groups = [np.flatnonzero(s_aoi == a) for a in np.unique(s_aoi)]
    return key, groups, t["lex_wetside_n1"].to_numpy(np.float64)[first] > 0


def _permute(idx: tuple, rng: np.random.Generator) -> np.ndarray:
    """集落ごとの水系側の印を、AOI の中で集落の間で入れ替える（セルは最寄りの集落の印をとる）。"""
    key, groups, s_flag = idx
    perm = np.empty(len(s_flag), dtype=bool)
    for m in groups:
        perm[m] = s_flag[rng.permutation(m)]
    return perm[key].astype(np.float64)


def b9d_permutation() -> dict:
    u = _uncapped()
    s = u.loc[_variation(u)].reset_index(drop=True)
    g = _planned_gfm()
    rng = np.random.default_rng(SEED)
    y = s["y"].to_numpy(np.float64)
    aoi = _codes(s["aoi"].to_numpy())
    z0, names = v39.design(s, ("wetside",), True, banded=True)
    j = names.index("wetside:le500m")
    obs_u = float(v68.fit_w(y, z0, aoi, np.ones(len(y)))[j])
    obs_g = float(np.log(_odds_planned(g)["odds_ratio"]))
    pu, pg = [], []
    iu, ig = _perm_index(s), _perm_index(g)
    t0 = time.time()
    for k in range(B_PERM):
        z, _ = v39.design(s.assign(lex_wetside_n1=_permute(iu, rng)), ("wetside",), True, banded=True)
        pu.append(float(v68.fit_w(y, z, aoi, np.ones(len(y)))[j]))
        pg.append(float(np.log(_odds_planned(g.assign(lex_wetside_n1=_permute(ig, rng)))["odds_ratio"])))
        if k % 20 == 0:
            print(f"  B9(d) permutation {k} ({time.time() - t0:.0f} s)", flush=True)
    pu_a, pg_a = np.array(pu), np.array(pg)
    return {"observed_or": {"ems_uncapped": math.exp(obs_u), "gfm_planned": math.exp(obs_g)},
            "p_two_sided": {"ems_uncapped": float((1 + (np.abs(pu_a) >= abs(obs_u)).sum()) / (1 + B_PERM)),
                            "gfm_planned": float((1 + (np.abs(pg_a) >= abs(obs_g)).sum()) / (1 + B_PERM))},
            "perm_or_q025_q50_q975": {"ems_uncapped": np.exp(np.quantile(pu_a, [0.025, 0.5, 0.975])).tolist(),
                                      "gfm_planned": np.exp(np.quantile(pg_a, [0.025, 0.5, 0.975])).tolist()},
            "n_perm": B_PERM,
            "rule": "B-13 B9 (d): water-related flags permuted among the settlements of each AOI (cells take the flag "
                    "of their nearest settlement); point estimates; two-sided p on |log OR|"}


def b9ef_labels() -> dict:
    g = _planned_gfm()
    built = v69._lc_built(g)
    d1 = g["P_d1"].to_numpy()
    band = np.where(d1 <= NEAR_M, "le500m", np.where(d1 <= 1000, "500m_1km", "gt1km"))
    cuts = np.nanquantile(built, [1 / 3, 2 / 3])
    bt = np.where(built <= cuts[0], "low", np.where(built <= cuts[1], "middle", "high"))
    e, f = g["y_ems"].to_numpy(), g["y_gfm"].to_numpy()
    tab: dict = {}
    for lab, key in (("band", band), ("built_tercile", bt)):
        tab[lab] = {}
        for k in np.unique(key):
            m = key == k
            tab[lab][str(k)] = {"n": int(m.sum()), "ems1_gfm1": int(((e == 1) & (f == 1) & m).sum()),
                                "ems1_gfm0": int(((e == 1) & (f == 0) & m).sum()),
                                "ems0_gfm1": int(((e == 0) & (f == 1) & m).sum()),
                                "ems0_gfm0": int(((e == 0) & (f == 0) & m).sum())}
    out: dict = {"e_two_by_two": tab, "built_cuts": cuts.tolist()}
    v41 = _load("v41_gfm_independent_labels")
    t = v41.build()
    sens = {a: s.get("sensor", "none") for a, s in _json(v69.GFM_SUMMARY)["aois"].items()}
    sen = t["aoi"].map(lambda a: sens.get(a, "none")).to_numpy()
    s = t.loc[np.isin(sen, ("s1", "sar", "optical"))].reset_index(drop=True)
    s = s.loc[_variation(s)].reset_index(drop=True)
    opt = (s["aoi"].map(lambda a: sens.get(a, "none")).to_numpy() == "optical").astype(np.float64)
    z, names = v39.design(s, ("wetside",), True, banded=True)
    near = ((s["P_d1"].to_numpy() <= NEAR_M) & (s["lex_wetside_n1"].to_numpy() > 0)).astype(np.float64)
    z = np.column_stack([z, near * opt])
    names = [*names, "wetside:le500m:optical"]
    res = v39.fit(s["y"].to_numpy(np.float64), z, s["aoi"].to_numpy(), s["unit"].to_numpy())
    jj, j0 = names.index("wetside:le500m:optical"), names.index("wetside:le500m")
    b, se = float(res["beta"][jj]), float(res["se"][jj])
    out["f_sensor"] = {"n_aoi": int(s["aoi"].nunique()), "n_aoi_optical": int(s.loc[opt > 0, "aoi"].nunique()),
                       "or_radar": float(np.exp(res["beta"][j0])), "or_optical": float(np.exp(res["beta"][j0] + b)),
                       "ratio_optical_to_radar": float(np.exp(b)),
                       "ratio_ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                       "p": float(2 * stats.norm.sf(abs(b / se)))}
    out["rule"] = ("B-13 B9 (e): counts on the planned GFM table by distance band and tercile of the 10 m built-up "
                   "share; (f): capped EMS table (v41), AOIs with a known sensor (Sentinel-1 and other SAR = radar, "
                   "optical), term wetside within 500 m x optical, cluster-robust by location group")
    return out


# %% B10: 辞書の点検
_TOK = re.compile(r"[^a-z]+")
_FOLD = str.maketrans({"ä": "a", "ö": "o", "ü": "u", "å": "a", "ø": "o", "Ä": "A", "Ö": "O", "Ü": "U", "Å": "A",
                       "Ø": "O"})


def _norm(s: str) -> str:
    s = str(s).translate(_FOLD).replace("æ", "ae").replace("Æ", "Ae").replace("ß", "ss")
    return "".join(ch for ch in unicodedata.normalize("NFKD", s) if not unicodedata.combining(ch))


def _analysed_settlements() -> pd.DataFrame:
    u = _uncapped()
    gids = np.unique(u["P_gid1"].to_numpy(np.int64))
    gn = _gn().reindex(gids)
    return gn[gn["asciiname"].notna()]


def b10_lexicon() -> dict:
    lex = Lexicon(LEX)
    gn = _analysed_settlements()
    c = dict.fromkeys(("ach_not_bach", "ach_with_bach", "ore_ora_all", "ore_ora_store", "mont_prefix",
                       "mont_prefix_monastery", "ried", "kacha", "taung", "be_eau_read_as_au", "it_read_by_de",
                       "ir_ab_word", "ir_ab_suffix", "fi_ala", "fi_yli_or_yla"), 0)
    pos: dict = {}
    uni: dict = {}
    for r in gn.itertuples():
        name, cc = str(r.asciiname), str(r.cc)
        toks = [x for x in _TOK.split(name.lower()) if x]
        hits = [e for e in lex.entries_for(cc) if any(_hit(tk, e) for tk in toks)]
        for e in hits:
            pos.setdefault(e.lang, {}).setdefault(e.pos, 0)
            pos[e.lang][e.pos] += 1
        f = {(e.lang, e.form, e.pos) for e in hits}
        if ("de", "ach", "suffix") in f:
            c["ach_with_bach" if any(tk.endswith("bach") for tk in toks) else "ach_not_bach"] += 1
        if any(x[:2] in (("no", "ore"), ("no", "ora")) for x in f):
            c["ore_ora_all"] += 1
            c["ore_ora_store"] += int(any(tk.startswith(("store", "stora")) for tk in toks))
        if ("fr", "mont", "prefix") in f:
            c["mont_prefix"] += 1
            c["mont_prefix_monastery"] += int(any(tk.startswith(("montreuil", "montier", "moutier")) for tk in toks))
        c["ried"] += int(any(x[:2] == ("de", "ried") for x in f))
        c["kacha"] += int(("fa", "kacha", "word") in f)
        c["taung"] += int(("my", "taung", "word") in f)
        c["be_eau_read_as_au"] += int(cc == "BE" and ("de", "au", "suffix") in f and any(tk.endswith("eau")
                                                                                         for tk in toks))
        c["it_read_by_de"] += int(cc == "IT" and any(x[0] == "de" for x in f))
        if cc == "IR":
            c["ir_ab_word"] += int(("fa", "ab", "word") in f)
            c["ir_ab_suffix"] += int(("fa", "ab", "suffix") in f)
        if cc == "FI":
            c["fi_ala"] += int("ala" in toks)
            c["fi_yli_or_yla"] += int(any(tk in ("yli", "yla") for tk in toks))
        a = any(lex.match(name, cc)[k] for k in WET_SIDE)
        b = any(lex.match(_norm(r.name), cc)[k] for k in WET_SIDE)
        if a != b:
            lang = (lex.country_langs.get(cc) or ("none",))[0]
            d = uni.setdefault(lang, {"gained": 0, "lost": 0})
            d["gained" if b else "lost"] += 1
    return {"n_settlements": len(gn), "counts": c, "positions_by_language": pos, "unicode_vs_ascii_wetside": uni,
            "blocks_by_country": {cc: list(v) for cc, v in sorted(lex.country_langs.items())},
            "rule": "B-13 B10: settlements nearest to at least one evaluable cell (uncapped table), frozen lexicon "
                    "read with the corrected loader; positions counted per matched entry; Unicode names folded "
                    "(a/o/u umlauts, a-ring, o-slash, ae, ss, combining marks removed)"}


# %% B11: N-read
def b11_nread() -> dict:
    lex = Lexicon(LEX)
    gn = _analysed_settlements()
    nr = pd.read_parquet(DATA_DIR / "global" / "nread" / "nread.parquet", columns=["gid", "category"])
    nr = nr.drop_duplicates("gid").set_index("gid")["category"]
    t = gn.assign(category=nr.reindex(gn.index).fillna("missing").to_numpy())
    t["block"] = t["cc"].map(lambda cc: "+".join(lex.country_langs.get(str(cc), ("none",))))
    out: dict = {}
    for key in ("block", "cc"):
        rows = {}
        for k, s in t.groupby(key):
            vc = s["category"].value_counts(normalize=True)
            rows[str(k)] = {"n": len(s), "unknown": float(vc.get("unknown", 0.0)),
                            "non_environmental": float(vc.get("non_environmental", 0.0)),
                            "water_related": float(sum(vc.get(x, 0.0) for x in WET_SIDE))}
        out[f"by_{key}"] = rows
    out["overall"] = t["category"].value_counts(normalize=True).to_dict()
    out["rule"] = ("B-13 B11: N-read categories of the settlements nearest to at least one evaluable cell, by the "
                   "lexicon blocks of the country and by country (shares of settlements)")
    return out


# %% B12: 二つの符号化者
def b12_coders() -> dict:
    v44 = _load("v44_wikipedia_etymology")
    samp = json.loads(v44.SAMPLE.read_text(encoding="utf-8"))
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"]
           for x in (DATA_DIR / "verification" / "v44_to_code.jsonl").read_text(encoding="utf-8").splitlines()}
    c1 = pd.read_csv(DOCV / "v44_etym_codes.csv")[["id", "status", "cats"]].rename(
        columns={"status": "status_coder1", "cats": "codes_coder1"})
    c2 = pd.read_csv(DOCV / "v57b_codes_coder2.csv")[["id", "status", "categories"]].rename(
        columns={"status": "status_coder2", "categories": "codes_coder2"})
    d = c1.merge(c2, on="id", how="outer")
    d["gid"] = d["id"].map(key)
    d = d.join(meta[["name", "cc", "stratum"]], on="gid")
    audit: set[int] = set()
    audit_names: set[str] = set()
    for f in AUDIT_SAMPLES:
        p = DOCV / f
        if p.is_file():
            a = pd.read_csv(p)
            if "gid" in a.columns:
                audit |= set(a["gid"].astype(np.int64).tolist())
            else:  # v51 の標本は ID を持たないので名で除く（同名の集落もまとめて除く）
                audit_names |= set(a["name"].astype(str).str.casefold().tolist())
    excl = d["gid"].isin(audit) | d["name"].astype(str).str.casefold().isin(audit_names)
    pub = d.loc[~excl, ["gid", "name", "cc", "stratum", "status_coder1", "codes_coder1", "status_coder2",
                        "codes_coder2"]].sort_values("gid")
    pub.to_csv(TWO_CODERS, index=False, quoting=csv.QUOTE_MINIMAL)
    s2 = _json(DOCV / "v57b_second_coder.json")
    return {"n_items": len(d), "n_published": len(pub), "n_excluded_audit_overlap": int(excl.sum()),
            "water_related_kappa": s2.get("water_related_kappa"), "status_kappa": s2.get("status_kappa"),
            "category_kappa": s2.get("category_kappa"), "file": str(TWO_CODERS.relative_to(ROOT)),
            "rule": "B-13 B12: both coders' codes with GeoNames ID, country and stratum (names in any audit sample "
                    "left out); agreement from v57b"}


# %% B14: GFM での当て直しの崩れ
def b14_gfm() -> dict:
    g = _planned_gfm()
    by = g.groupby("aoi").agg(n=("y_ems", "size"), ems=("y_ems", "mean"), gfm=("y_gfm", "mean"))
    d = v36.load("world_v2_gfm", "primary")
    out: dict = {"n_aoi": len(by), "median_ems_share": float(by["ems"].median()),
                 "median_gfm_share": float(by["gfm"].median()),
                 "share_aoi_gfm_below_ems": float((by["gfm"] < by["ems"]).mean()),
                 "corr_ems_gfm_share": float(np.corrcoef(by["ems"], by["gfm"])[0, 1])}
    if d is not None:
        y = d["y"].astype(np.float64)
        aoi = np.asarray(d["aoi"]).astype(str)
        rows = []
        for a in np.unique(aoi):
            m = aoi == a
            ya = y[m]
            gain = float(np.mean((d["p_M0"][m] - ya) ** 2) - np.mean((d["p_M1"][m] - ya) ** 2))
            rows.append((a, gain, float(ya.mean()), int(m.sum())))
        r = pd.DataFrame(rows, columns=["aoi", "gain", "gfm_share", "n"])
        tot = float((r["gain"] * r["n"]).sum() / r["n"].sum())
        out |= {"gfm_run_n_aoi": len(r), "yardstick_gfm": tot, "median_gain_per_aoi": float(r["gain"].median()),
                "share_aoi_negative_gain": float((r["gain"] < 0).mean()),
                "share_of_gain_top10_aoi": float(r.nlargest(10, "gain").eval("gain * n").sum()
                                                 / (r["gain"] * r["n"]).sum()) if (r["gain"] * r["n"]).sum() else None,
                "yardstick_ems_primary": float(_json(RUNS_DIR / "world_v2" / "metrics_primary.json")["primary"]
                                               ["delta_bs_M0_M1"])}
    out["rule"] = ("B-13 B14: flooded shares of EMS and GFM by AOI on the planned GFM table; terrain gain (Brier of "
                   "M0 minus M1) by AOI in the run trained and scored on GFM labels (world_v2_gfm)")
    return out


# %% B15: 相対の修飾語の対
def _site_table() -> pd.DataFrame:
    """集落ごと: 国、座標、地点の HAND（177 m 以内のセルの hand_min の平均）、500 m 以内の EMS の浸水割合、場所グループ。"""
    t = pd.read_parquet(VER / "v37_table.parquet", columns=["aoi", "unit", "y", "P_gid1", "P_d1", "hand_min",
                                                             "set_lat", "set_lon", "cc"])
    t = t[t["P_gid1"] >= 0]
    site = t[t["P_d1"] <= SITE_M].groupby("P_gid1")["hand_min"].mean()
    near = t[t["P_d1"] <= NEAR_M].groupby("P_gid1")["y"].agg(["mean", "size"])
    first = t.groupby("P_gid1").agg(cc=("cc", "first"), lat=("set_lat", "first"), lon=("set_lon", "first"),
                                    unit=("unit", lambda s: _mode(s.to_numpy(dtype=str))))
    s = first.join(site.rename("hand_site")).join(near.rename(columns={"mean": "ems_le500", "size": "n_le500"}))
    g = _planned_gfm()
    gn_ = g[(g["P_d1"] <= NEAR_M) & (g["P_gid1"] >= 0)].groupby("P_gid1")["y_gfm"].mean().rename("gfm_le500")
    s = s.join(gn_)
    lat60, lon60 = s["lat"] * 60, s["lon"] * 60
    s["rounded"] = ((lat60 - lat60.round()).abs() < v37.MINUTE_TOL) & ((lon60 - lon60.round()).abs() < v37.MINUTE_TOL)
    s["asciiname"] = _gn()["asciiname"].reindex(s.index).to_numpy()
    return s


def b15_pairs() -> dict:
    lex = Lexicon(LEX)
    s = _site_table()
    s = s[~s["rounded"] & s["asciiname"].notna()]
    rows = []
    for cc, part in s.groupby("cc"):
        langs = lex.country_langs.get(str(cc), ())
        side_of = qualifier_side(QUALIFIERS, tuple(langs))
        pre = tuple(f for lang in langs for f in PREFIX_QUALIFIERS.get(lang, ()))
        if not side_of:
            continue
        sq = [(gid, split_qualifier(r.asciiname, side_of, pre)) for gid, r in part.iterrows()]
        df = pd.DataFrame([(gid, x[0], x[1]) for gid, x in sq if x], columns=["gid", "side", "base"])
        for _, grp in df.groupby("base"):
            lo, hi = grp[grp["side"] == "lower"]["gid"].to_numpy(), grp[grp["side"] == "upper"]["gid"].to_numpy()
            if not len(lo) or not len(hi):
                continue
            lat0 = math.radians(float(part.loc[np.concatenate([lo, hi]), "lat"].mean()))

            def xy(g: np.ndarray, lat0: float = lat0, part: pd.DataFrame = part) -> np.ndarray:
                return np.column_stack([part.loc[g, "lon"].to_numpy() * 111_320 * math.cos(lat0),
                                        part.loc[g, "lat"].to_numpy() * 110_540])

            for i, k, dist in match_pairs(xy(lo), xy(hi), PAIR_MAX_M):
                a, b = part.loc[lo[i]], part.loc[hi[k]]
                rows.append({"cc": str(cc), "lower": int(lo[i]), "upper": int(hi[k]), "distance_m": dist,
                             "unit": a["unit"], "hand_diff": a["hand_site"] - b["hand_site"],
                             "ems_diff": a["ems_le500"] - b["ems_le500"], "gfm_diff": a["gfm_le500"] - b["gfm_le500"],
                             "lower_name": a["asciiname"], "upper_name": b["asciiname"]})
    p = pd.DataFrame(rows)
    out: dict = {"n_pairs": len(p), "by_country": p["cc"].value_counts().to_dict() if len(p) else {}}
    rng = np.random.default_rng(SEED)
    for col in ("hand_diff", "ems_diff", "gfm_diff"):
        q = p[np.isfinite(p[col].to_numpy(np.float64))] if len(p) else p
        if len(q) < 5:
            out[col] = {"n": len(q)}
            continue
        x = q[col].to_numpy(np.float64)
        neg, pos = int((x < 0).sum()), int((x > 0).sum())
        units, inv = np.unique(q["unit"].to_numpy(dtype=str), return_inverse=True)
        med, share = [], []
        for _ in range(B_PAIR_BOOT):
            w = rng.multinomial(len(units), np.full(len(units), 1 / len(units)))[inv]
            xs = np.repeat(x, w)
            if len(xs):
                med.append(float(np.median(xs)))
                share.append(float((xs < 0).mean()))
        out[col] = {"n": len(x), "median": float(np.median(x)), "share_lower_member_lower": neg / max(1, neg + pos),
                    "n_negative": neg, "n_positive": pos, "n_units": len(units),
                    "sign_test_p": float(stats.binomtest(neg, neg + pos).pvalue) if neg + pos else None,
                    "median_ci90": np.quantile(med, [0.05, 0.95]).tolist(),
                    "share_ci90": np.quantile(share, [0.05, 0.95]).tolist()}
    p.to_csv(VER / "v70_b15_pairs.csv", index=False)
    out["rule"] = ("B-13 B15: pairs of settlements in the same country with the same base name once one relative "
                   "qualifier (C2 (d)) is removed, lower and upper member within 10 km (greedy by distance), both with "
                   "coordinates not on whole arc-minutes; differences lower minus upper in HAND at the site, and in "
                   "the flooded share of cells within 500 m (EMS, planned GFM table); sign test and location-group "
                   "bootstrap (2,000)")
    return out


# %% 実行
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    want = {x for x in args.only.split(",") if x}

    def run(key: str, fn: callable) -> None:
        if want and key not in want:
            return
        t0 = time.time()
        val = fn()
        res = _json(OUT) if OUT.is_file() else {}
        res |= {"rule": "amendment B-13 (post hoc), local computations, part 2", "seed": SEED, key: val}
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        print(f"{key} ({time.time() - t0:.0f} s) {json.dumps(val, default=float)[:700]}", flush=True)

    for key, fn in (("B12", b12_coders), ("B10", b10_lexicon), ("B11", b11_nread), ("B14", b14_gfm),
                    ("B15", b15_pairs), ("B8", b8_joint_wald), ("B9ef", b9ef_labels), ("B6", b6_influence),
                    ("B9", b9_near), ("B9d", b9d_permutation)):
        run(key, fn)


if __name__ == "__main__":
    main()
