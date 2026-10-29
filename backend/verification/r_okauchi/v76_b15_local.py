# %% v76（事後・amendment B-15 の B1〜B11）: 初見の査読 第 11 回に応じた点検
"""確証判定は変えない。B-15 で計算の前に固定した量を、既存の held-out の予測・二段目の保存・表から計算する。

- B1: 事前に定めた全ての比較の表（確証の実行の metrics から。計算はしない）
- B2: 少数クラスタの参照（臨界値 z・t(G−1)・t(ν) で Fieller の片側 95% 上限）
- B3: super-cluster の再解析（活性化を分けない分割の二段目の保存、super-cluster を引いて二段目を当て直す）
- B4: 名と記録の対の差（集落を通じて割り当てた記録と、同じセル・同じ再標本で）
- B5: 集落の輪の観測のしやすさ（市街地・樹木・GFM で評価できない画素・部分的な恒久水）
- B6: GFM の診断（EMS で学習したモデルの GFM のラベルでの物差し、AOI 内の AUC・CSI、同じ日の場面で分けた地域内の結果）
- B7: 読みの表（ブロック、地名辞典、N-read、類ごとの精度、ローマ字化と規則、Wikidata の名の由来）
- B8: 目録の拡張（v36〜v76 の全ての出力、族ごとの BH・BY と向き）
- B9: ノルウェーの直しがパキスタンの組に届いた仕組み（B3 の保存で凍結した読みを再現）
- B10: 集落の近くの推定の推測（wild cluster bootstrap、ラベルを AOI の中で入れ替えた OR）
- B11: 点検（表 S9 の fa の層、大きさの三分位、48% と 49% など）

出力: stdout と docs/verification/v76_b15_local.json（部分ごとに足す。`--only B1,B2` で一部だけ）。
"""

from __future__ import annotations

import argparse
import datetime as dt
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
from topsem.world.b13 import SESOI_R, fieller_interval
from topsem.world.lexicon import WET_SIDE, Lexicon

OUT = ROOT / "docs" / "verification" / "v76_b15_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
SEED = 20261104
B_FIXED = 4000
B_SC = 500
B_WILD = 9999
B_PERM = 200
Z90 = float(stats.norm.ppf(0.95))  # 両側 90% = 片側 95%
N_WORKERS = int(os.environ.get("V76_WORKERS", "5"))
SC_RUN = "world_v2_b15sc_stage2_local"
SPLITS = ("primary", "country", "temporal", "llmcut")
WIKIDATA_CACHE = VER / "v76_wikidata_named_after.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v73 = _load("v73_b14_local")
v70, v69, v68, v36, v39 = v73.v70, v73.v69, v73.v68, v73.v36, v73.v39


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


# %% B1: 事前に定めた全ての比較の表
B1_KEYS = ("primary", "placebo", "total", "total_placebo", "flexible", "flexible_placebo", "H1", "H1_flexible", "H4",
           "H4_placebo", "H4_total", "H4_total_placebo", "H4_flexible", "H4_flexible_placebo")


def _scalars(b: dict) -> dict:
    return {k: v for k, v in b.items() if isinstance(v, (int, float, str, bool)) or v is None
            or (isinstance(v, list) and len(v) == 2 and all(isinstance(x, (int, float)) for x in v))}


def b1_prespecified() -> dict:
    out: dict = {}
    for sp in SPLITS:
        f = RUNS_DIR / "world_v2" / f"metrics_{sp}.json"
        if not f.is_file():
            continue
        m = _json(f)
        row = {k: _scalars(m[k]) for k in B1_KEYS if isinstance(m.get(k), dict)}
        ch = m.get("channels") or {}
        for k in ("M3nx", "M3nl"):
            if isinstance(ch.get(k), dict):
                row[f"channel_{k}"] = _scalars(ch[k])
        row["n_cells"], row["n_units"] = m.get("n_cells"), m.get("n_units")
        out[sp] = row
    out["rule"] = ("B-15 B1: every pre-specified comparison of the confirmatory run (world_v2) on the four splits, copied "
                   "from its metrics (nested, channels, total, flexible, H1, H4 and their placebos); no new computation")
    return out


# %% B2: 少数クラスタの参照
def _nu(x: np.ndarray, n: np.ndarray) -> float:
    """単位ごとの寄与 x（和）と行の数 n からの Satterthwaite の有効自由度 (Σγ)² / Σγ²。"""
    big_n = float(n.sum())
    xhat = float(x.sum()) / big_n
    g = ((x - xhat * n) / big_n) ** 2
    return float(g.sum() ** 2 / (g ** 2).sum())


def _unit_contribs() -> pd.DataFrame:
    """主分割の直したコードの二段目の保存（B-8 の C2）の評価側のセルごとの寄与と、場所グループ・activation・国の符号。"""
    c = np.load(v68.B7_CODES)
    parts = []
    for k in range(10):
        z = np.load(v68.B7_SRC / f"fold{k}_stage2.npz", allow_pickle=False)
        y = z["y_te"].astype(np.float64)
        sq = {a: (z[f"p_{a}_te"].astype(np.float64) - y) ** 2 for a in ("M0", "M1", "M1c", "M3n")}
        parts.append(pd.DataFrame({"e": sq["M1c"] - sq["M3n"], "d": sq["M0"] - sq["M1"], "u": c[f"u_te{k}"],
                                   "a": c[f"a_te{k}"], "c": c[f"c_te{k}"]}))
    return pd.concat(parts, ignore_index=True)


def _nus(df: pd.DataFrame, key: str) -> dict:
    g = df.groupby(key).agg(n=("e", "size"), e=("e", "sum"), d=("d", "sum"))
    n, e, d = g["n"].to_numpy(np.float64), g["e"].to_numpy(), g["d"].to_numpy()
    return {"G": len(g), "nu_delta": _nu(e, n), "nu_yardstick": _nu(d, n), "nu_l_plus": _nu(e - SESOI_R * d, n),
            "nu_l_minus": _nu(e + SESOI_R * d, n)}


def _crit_rows(e0: float, d0: float, var: tuple[float, float, float], g: int, nu: float) -> dict:
    out = {}
    for lab, c in (("z", Z90), ("t_G_minus_1", float(stats.t.ppf(0.95, g - 1))), ("t_nu", float(stats.t.ppf(0.95, nu)))):
        fi = fieller_interval(e0, d0, var, z=c)
        out[lab] = {"critical_value": c, "fieller90": list(fi) if fi else None, "upper_R": fi[1] if fi else None}
    return out


def _cov3(eb: np.ndarray, db: np.ndarray) -> tuple[float, float, float]:
    ok = np.isfinite(eb) & np.isfinite(db)
    c = np.cov(eb[ok], db[ok])
    return float(c[0, 0]), float(c[1, 1]), float(c[0, 1])


def b2_few_clusters() -> dict:
    df = _unit_contribs()
    coa, _ = v70._country_of_activation()
    df["country_of_act"] = coa[df["a"].to_numpy()]
    nus = {"location_groups": _nus(df, "u"), "activations": _nus(df, "a"), "countries": _nus(df, "country_of_act")}
    v70j = _json(DOCV / "v70_b13_local.json")
    out: dict = {"effective_df": nus, "constructions": {}}
    # (a): B-13 の B1(a) の重みを同じ種・同じ順で作り直す（確証の実行の分が先）
    ra = np.random.default_rng(v70.SEED)
    conf = v36.load("world_v2", "primary")
    n_c, s_c, _ = v69._group_sums(conf, v73.ARMS4)
    wc = ra.multinomial(len(n_c), np.full(len(n_c), 1.0 / len(n_c)), size=v70.B_FIXED).astype(np.float64)
    lex = v36.load(*v73.CORR)
    n, s, _ = v69._group_sums(lex, v73.ARMS4)
    w = ra.multinomial(len(n), np.full(len(n), 1.0 / len(n)), size=v70.B_FIXED).astype(np.float64)
    gl = nus["location_groups"]
    for key, nn, ss, ww in (("protocol_run_fixed_groups", n_c, s_c, wc), ("a_fixed_groups", n, s, w)):
        eg, dg = ss["M1c"] - ss["M3n"], ss["M0"] - ss["M1"]
        e0, d0 = float(eg.sum() / nn.sum()), float(dg.sum() / nn.sum())
        nw = ww @ nn
        var = _cov3((ww @ eg) / nw, (ww @ dg) / nw)
        out["constructions"][key] = {"R": e0 / d0, "var": list(var), "units": "location groups", "G": len(nn),
                                     "nu": gl["nu_l_plus"], **_crit_rows(e0, d0, var, len(nn), gl["nu_l_plus"])}
    lx = v70j["B1bc"]["b_groups"]["refit"]
    for key, scheme, units, nu_key in (("b_refit_groups", "groups", "location groups", "location_groups"),
                                       ("c_refit_nested", "nested", "countries, then activations", "countries")):
        r = np.load(VER / f"v70_b1_refit_reps_{scheme}.npy")
        wt = r[:, 6]
        var = _cov3((r[:, 0] - r[:, 1]) / wt, (r[:, 4] - r[:, 5]) / wt)
        g = nus[nu_key]["G"]
        out["constructions"][key] = {"R": lx["delta_bs"] / lx["yardstick"], "var": list(var), "units": units, "G": g,
                                     "nu": nus[nu_key]["nu_l_plus"],
                                     **_crit_rows(lx["delta_bs"], lx["yardstick"], var, g, nus[nu_key]["nu_l_plus"])}
    c = v70j["B1bc"]["c_nested"]["refit"]
    nb, nbc = v70j["B1d"]["nb_var"], float(v70j["B1d"]["nb_cov"])
    var_d = (c["var"][0] + nb["delta_bs"], c["var"][1] + nb["yardstick"], c["var"][2] + nbc)
    g = nus["countries"]["G"]
    out["constructions"]["d_refit_nested_folds"] = {
        "R": c["R"], "var": list(var_d), "units": "countries, then activations; fold assignments", "G": g,
        "nu": nus["countries"]["nu_l_plus"],
        "check_z_matches_b13": fieller_interval(c["delta_bs"], c["yardstick"], var_d, z=Z90),
        **_crit_rows(c["delta_bs"], c["yardstick"], var_d, g, nus["countries"]["nu_l_plus"])}
    r = np.load(VER / "v73_b4_pigeonhole_reps.npy")
    pv = _json(DOCV / "v73_b14_local.json")["B4"]["refit"]
    var_p = _cov3((r[:, 0] - r[:, 1]) / r[:, 6], (r[:, 4] - r[:, 5]) / r[:, 6])
    nu_p = min(nus["location_groups"]["nu_l_plus"], nus["activations"]["nu_l_plus"])
    g_p = min(nus["location_groups"]["G"], nus["activations"]["G"])
    out["constructions"]["pigeonhole_refit"] = {
        "R": pv["R"], "var": list(var_p), "units": "location groups x activations", "G": g_p, "nu": nu_p,
        **_crit_rows(pv["delta_bs"], pv["yardstick"], var_p, g_p, nu_p)}
    tw = _json(DOCV / "v68_b11_local.json")["B7"]["two_way"]
    out["constructions"]["b11_two_way_refit"] = {"ratio_pct90": tw["ci90_R_refit_ratio"],
                                                 "note": "no replicates saved; only the percentile interval"}
    ups = {k: v["z"]["upper_R"] for k, v in out["constructions"].items() if isinstance(v.get("z"), dict)}
    out["range_upper_z_refit"] = [min(v for k, v in ups.items() if "refit" in k),
                                  max(v for k, v in ups.items() if "refit" in k)]
    tn = {k: v["t_nu"]["upper_R"] for k, v in out["constructions"].items() if isinstance(v.get("t_nu"), dict)}
    out["max_upper_t_nu_refit"] = max((v for k, v in tn.items() if "refit" in k and v is not None), default=None)
    out["rule"] = ("B-15 B2: one-sided 95% upper limit of R (upper end of the 90% Fieller interval) with the resample "
                   "covariance of each construction and critical values z, t(G-1) and t(nu); nu = Satterthwaite "
                   "effective df (sum gamma)^2 / sum gamma^2 of the unit contributions to L+ (location groups for (a) and "
                   "(b), countries for (c) and (d), the smaller of location groups and activations for the pigeonhole)")
    return out


# %% B3・B9: super-cluster の再解析（活性化を分けない分割の二段目の保存）
_S: dict = {}


def _sc_folds() -> list[dict]:
    src = RUNS_DIR / SC_RUN / "predictions" / "primary"
    out = []
    for k in range(10):
        f = src / f"fold{k}_stage2.npz"
        if not f.is_file():
            continue
        z = np.load(f, allow_pickle=False)
        out.append({kk: z[kk] for kk in z.files})
    return out


def _demean(x: np.ndarray, aoi: np.ndarray) -> np.ndarray:
    m = pd.DataFrame(x.astype(np.float64)).groupby(pd.Series(aoi).to_numpy()).transform("mean").to_numpy()
    return x.astype(np.float64) - m


def _composition(xn: np.ndarray, aoi: np.ndarray) -> np.ndarray:
    return pd.DataFrame(xn.astype(np.float64)).groupby(pd.Series(aoi).to_numpy()).transform("mean").to_numpy()


def _designs(f: dict, *, frozen: bool = False) -> dict:
    """腕ごとの二段目の説明変数（学習側, 評価側）。frozen=True はノルウェーの名の辞書の特徴を 0 にした凍結の読み。"""
    names = [str(c) for c in f["name_cols"]]
    lexcols = np.array([c.startswith("lex_") for c in names])
    xn_tr, xn_te = f["xn_tr"].astype(np.float64).copy(), f["xn_te"].astype(np.float64).copy()
    xc_tr, xc_te = f["xc_tr"].astype(np.float64), f["xc_te"].astype(np.float64)
    if frozen:
        for x, cc in ((xn_tr, f["cc_tr"]), (xn_te, f["cc_te"])):
            m = np.asarray(cc).astype(str) == "NO"
            x[np.ix_(m, lexcols)] = 0.0
        xc_tr, xc_te = _composition(xn_tr, f["aoi_tr"]), _composition(xn_te, f["aoi_te"])
    empty_tr, empty_te = np.zeros((len(xn_tr), 0)), np.zeros((len(xn_te), 0))
    return {"M1c": (xc_tr, xc_te), "M3n": (np.column_stack([xc_tr, xn_tr]), np.column_stack([xc_te, xn_te])),
            "M1s": (empty_tr, empty_te), "M3d": (_demean(xn_tr, f["aoi_tr"]), _demean(xn_te, f["aoi_te"]))}


ARMS_SC = ("M1c", "M3n", "M1s", "M3d")


def _sc_init() -> None:
    _threads1()
    folds = _sc_folds()
    units = np.unique(np.concatenate([np.concatenate([f["unit_tr"], f["unit_te"]]) for f in folds]).astype(str))
    idx = {u: i for i, u in enumerate(units)}
    _S["nu"] = len(units)
    _S["folds"] = []
    for f in folds:
        d = {"z_tr": f["z_tr"], "z_te": f["z_te"], "y_tr": f["y_tr"], "y_te": f["y_te"],
             "u_tr": np.array([idx[u] for u in f["unit_tr"].astype(str)]),
             "u_te": np.array([idx[u] for u in f["unit_te"].astype(str)]),
             "p_M0_te": f["p_M0_te"], "p_M1_te": f["p_M1_te"], "x": _designs(f), "init": {}}
        for a in ARMS_SC:
            xtr, xte = d["x"][a]
            _, d["init"][a] = stage2_weighted(d["z_tr"], d["z_te"], xtr, xte, d["y_tr"])
        _S["folds"].append(d)


def _sc_rep(b: int) -> list[float]:
    rng = np.random.default_rng(SEED + b)
    wu = rng.multinomial(_S["nu"], np.full(_S["nu"], 1.0 / _S["nu"])).astype(np.float64)
    acc = np.zeros(len(ARMS_SC) + 3)
    for f in _S["folds"]:
        w_tr, w_te = wu[f["u_tr"]], wu[f["u_te"]]
        if w_tr.sum() == 0 or w_te.sum() == 0:
            continue
        y = f["y_te"].astype(np.float64)
        for i, a in enumerate(ARMS_SC):
            xtr, xte = f["x"][a]
            p, _ = stage2_weighted(f["z_tr"], f["z_te"], xtr, xte, f["y_tr"], w_tr=w_tr, init=f["init"][a])
            acc[i] += float((w_te * (p - y) ** 2).sum())
        acc[len(ARMS_SC)] += float((w_te * (f["p_M0_te"].astype(np.float64) - y) ** 2).sum())
        acc[len(ARMS_SC) + 1] += float((w_te * (f["p_M1_te"].astype(np.float64) - y) ** 2).sum())
        acc[len(ARMS_SC) + 2] += float(w_te.sum())
    return acc.tolist()


def _point_sc(folds: list[dict], *, frozen: bool = False) -> dict:
    """重みの無い当て直しの点の推定（腕ごとの SSE、M0・M1、場所グループごとの寄与）。"""
    sse: dict = {a: 0.0 for a in (*ARMS_SC, "M0", "M1")}
    unit_rows = []
    coefs = {}
    n = 0
    for k, f in enumerate(folds):
        y = f["y_te"].astype(np.float64)
        x = _designs(f, frozen=frozen)
        pr = {}
        for a in ARMS_SC:
            xtr, xte = x[a]
            pr[a], coefs[f"{a}_fold{k}"] = stage2_weighted(f["z_tr"], f["z_te"], xtr, xte, f["y_tr"])
            sse[a] += float(((pr[a] - y) ** 2).sum())
        for a in ("M0", "M1"):
            sse[a] += float(((f[f"p_{a}_te"].astype(np.float64) - y) ** 2).sum())
        n += len(y)
        unit_rows.append(pd.DataFrame({"unit": f["unit_te"].astype(str), "aoi": f["aoi_te"].astype(str),
                                       "cc": f["cc_te"].astype(str),
                                       "e": (pr["M1c"] - y) ** 2 - (pr["M3n"] - y) ** 2,
                                       "d": (f["p_M0_te"] - y) ** 2 - (f["p_M1_te"] - y) ** 2}))
    bs = {a: v / n for a, v in sse.items()}
    return {"brier": bs, "n": n, "units": pd.concat(unit_rows, ignore_index=True), "coefs": coefs}


def b3_supercluster() -> dict:
    folds = _sc_folds()
    if len(folds) < 10:
        return {"note": f"stage-2 dump of {SC_RUN} incomplete ({len(folds)} folds)"}
    units = np.unique(np.concatenate([f["unit_te"] for f in folds]).astype(str))
    pt = _point_sc(folds)
    bs, yard = pt["brier"], pt["brier"]["M0"] - pt["brier"]["M1"]
    cache = VER / "v76_b3_supercluster_reps.npy"
    if cache.is_file():
        r = np.load(cache)
    else:
        _threads1()
        t0 = time.time()
        with ProcessPoolExecutor(N_WORKERS, initializer=_sc_init) as ex:
            r = np.array(list(ex.map(_sc_rep, range(B_SC), chunksize=2)), dtype=np.float64)
        np.save(cache, r)
        print(f"  B3: {len(r)} replicates ({time.time() - t0:.0f} s)", flush=True)
    wt = r[:, -1]
    db = (r[:, len(ARMS_SC)] - r[:, len(ARMS_SC) + 1]) / wt
    u = pt["units"]
    g = u.groupby("unit").agg(n=("e", "size"), e=("e", "sum"), d=("d", "sum"))
    nu = _nu(g["e"].to_numpy() - SESOI_R * g["d"].to_numpy(), g["n"].to_numpy(np.float64))
    out: dict = {"n_superclusters": len(units), "n_cells": pt["n"], "brier": bs, "yardstick": yard,
                 "nu_l_plus": nu, "comparisons": {}}
    rng = np.random.default_rng(SEED)
    for key, (a, b) in (("M1c_M3n", ("M1c", "M3n")), ("M1s_M3d", ("M1s", "M3d"))):
        ia, ib = ARMS_SC.index(a), ARMS_SC.index(b)
        eb = (r[:, ia] - r[:, ib]) / wt
        e0 = bs[a] - bs[b]
        det = v73._levels_detail(e0, yard, eb, db, rng)
        var = _cov3(eb, db)
        out["comparisons"][key] = {"delta_bs": e0, "R": e0 / yard, "detail": det,
                                   **_crit_rows(e0, yard, var, len(units), nu)}
    b1 = _json(DOCV / "v68_b11_local.json").get("B1") or {}
    out["cluster_runs_R"] = {"corrected_world_v2_b7actfix": b1.get("R_corrected"),
                             "frozen_world_v2_actgrp": b1.get("R_frozen")}
    out["rule"] = ("B-15 B3: split with activations kept whole (super-clusters of location groups linked by a shared "
                   "activation), stage-2 design saved locally (world_v2_b15sc_stage2_local); 500 draws of super-clusters "
                   "(seed 20261104 + b), stage 2 refitted with weights in each draw; M1c -> M3n and M1s -> M3d")
    return out


def b9_norway() -> dict:
    folds = _sc_folds()
    if len(folds) < 10:
        return {"note": f"stage-2 dump of {SC_RUN} incomplete ({len(folds)} folds)"}
    cor, fro = _point_sc(folds), _point_sc(folds, frozen=True)
    out: dict = {}
    for lab, pt in (("corrected", cor), ("frozen_emulated", fro)):
        bs = pt["brier"]
        out[lab] = {"brier": bs, "R": (bs["M1c"] - bs["M3n"]) / (bs["M0"] - bs["M1"])}
    gc = cor["units"].groupby("unit").agg(n=("e", "size"), e=("e", "sum"), share_no=("cc", lambda s: float((s == "NO")
                                                                                                    .mean())))
    gf = fro["units"].groupby("unit").agg(e=("e", "sum"))
    big_n = float(gc["n"].sum())
    yard = out["corrected"]["brier"]["M0"] - out["corrected"]["brier"]["M1"]
    dchange = ((gc["e"] - gf["e"]) / big_n / yard).sort_values()
    top = dchange.abs().sort_values(ascending=False).head(5).index
    cc_of = cor["units"].groupby("unit")["cc"].agg(lambda s: s.value_counts().index[0])
    out["change_R_total"] = out["corrected"]["R"] - out["frozen_emulated"]["R"]
    b1 = _json(DOCV / "v68_b11_local.json").get("B1") or {}
    out["cluster_runs"] = {k: b1.get(k) for k in ("R_frozen", "R_corrected", "change_R", "change_R_norwegian_units")}
    out["top_units"] = [{"unit": str(uu), "country": str(cc_of[uu]), "change_R": float(dchange[uu]),
                         "share_cells_norwegian": float(gc.loc[uu, "share_no"]), "n_cells": int(gc.loc[uu, "n"])}
                        for uu in top]
    out["change_R_norwegian_units"] = float(dchange[gc["share_no"] > 0].sum())
    names = [str(c) for c in folds[0]["name_cols"]]
    comp = [f"c_{c}" for c in names]
    cols = [*comp, *names]
    dif = []
    for k in range(10):
        cc_, cf = np.asarray(cor["coefs"][f"M3n_fold{k}"][0]), np.asarray(fro["coefs"][f"M3n_fold{k}"][0])
        dif.append(cc_ - cf)
    md = np.mean(dif, axis=0)  # 標準化した尺度、[基礎の logit, 組成..., 名...]
    out["mean_coef_change_M3n"] = dict(zip(["base_logit", *cols], (float(x) for x in md), strict=True))
    out["rule"] = ("B-15 B9: stage 2 refitted on the dump of B3 with the corrected reading and with the frozen reading "
                   "emulated (lexicon features of cells whose nearest settlement is in Norway set to 0, composition "
                   "recomputed); change in R by super-cluster; mean change of the standardised M3n coefficients")
    return out


# %% B4: 名と記録の対の差
B4_RECORDS = ("gsw_from_2022", "ems_earlier_events", "jrc_100yr")


def b4_paired() -> dict:
    rng = np.random.default_rng(v73.SEED)  # B-14 の B6 と同じ種・同じ順で重みを引く
    out: dict = {}
    cache: dict = {}
    for key, run_s, arms_s, run_c, _ in v73.RECORDS:
        for r_ in (run_s, run_c):
            if r_ not in cache:
                cache[r_] = v36.load(r_, "primary")
        ds, dcl = cache[run_s], cache[run_c]
        if run_s != run_c and not (np.array_equal(ds["y"], dcl["y"]) and np.array_equal(ds["aoi"], dcl["aoi"])):
            continue  # B-14 の B6 もこの記録では重みを引かない（引く順をそろえる）
        y = ds["y"].astype(np.float64)
        units = np.unique(np.asarray(ds["unit"]))
        w = v73._boot_w(rng, len(units))
        if key not in B4_RECORDS or "p_M3n" not in ds:
            continue

        def gain(a: str, b: str, d: dict = ds, y: np.ndarray = y) -> np.ndarray:
            return (d[f"p_{a}"].astype(np.float64) - y) ** 2 - (d[f"p_{b}"].astype(np.float64) - y) ** 2

        yard = (ds["p_M0"].astype(np.float64) - y) ** 2 - (ds["p_M1"].astype(np.float64) - y) ** 2
        _, _, en, er, dd = v73._gsum(np.asarray(ds["unit"]), gain("M1c", "M3n"), gain(*arms_s), yard)
        rn, rr = (w @ en) / (w @ dd), (w @ er) / (w @ dd)
        dif = rn - rr
        out[key] = {"run": run_s, "arms_record": list(arms_s), "R_names": float(en.sum() / dd.sum()),
                    "R_record": float(er.sum() / dd.sum()), "difference": float((en.sum() - er.sum()) / dd.sum()),
                    "difference_pct90": np.quantile(dif, [0.05, 0.95]).tolist(),
                    "share_reps_names_below_record": float((dif < 0).mean()), "n_cells": len(y), "n_units": len(units)}
    # 集落の地点の HAND（陽性対照、M1 に）: 主の実行と同じ行。B-15 の種で新しく引く
    o, lx = v36.load("world_v2_b14oracle", "primary"), v36.load(*v73.CORR)
    if o is not None and np.array_equal(o["y"], lx["y"]) and np.array_equal(o["aoi"], lx["aoi"]):
        y = lx["y"].astype(np.float64)

        def gain2(d: dict, a: str, b: str) -> np.ndarray:
            return (d[f"p_{a}"].astype(np.float64) - y) ** 2 - (d[f"p_{b}"].astype(np.float64) - y) ** 2

        yard = (lx["p_M0"].astype(np.float64) - y) ** 2 - (lx["p_M1"].astype(np.float64) - y) ** 2
        units, _, en, er, dd = v73._gsum(np.asarray(lx["unit"]), gain2(lx, "M1c", "M3n"), gain2(o, "M1co", "M3o"), yard)
        w = v73._boot_w(np.random.default_rng(SEED), len(units))
        dif = (w @ en) / (w @ dd) - (w @ er) / (w @ dd)
        out["site_hand_to_M1"] = {"R_names": float(en.sum() / dd.sum()), "R_record": float(er.sum() / dd.sum()),
                                  "difference": float((en.sum() - er.sum()) / dd.sum()),
                                  "difference_pct90": np.quantile(dif, [0.05, 0.95]).tolist(),
                                  "share_reps_names_below_record": float((dif < 0).mean()), "n_cells": len(y),
                                  "n_units": len(units), "draws": "new draws, seed 20261104"}
    else:
        out["site_hand_to_M1"] = {"note": "rows of world_v2_b14oracle and world_v2_lexfix differ"}
    out["rule"] = ("B-15 B4: R(names, M1c -> M3n) minus R(record assigned through settlements) on the record's cells and "
                   "the location-group draws of B-14 B6 (models fixed); site HAND against the corrected run with new "
                   "draws; 90% percentile intervals")
    return out


# %% B5: 集落の輪の観測のしやすさ
def _cell_join(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for a, ia in t.groupby("aoi").indices.items():
        d = v69.LABELS / a
        f = pq.read_table(d / "features.parquet", columns=["lc_built", "lc_tree"]).to_pandas()
        c = pq.read_table(d / "cells.parquet", columns=["ix", "iy", "f_perm_water"]).to_pandas()
        x = pd.concat([c, f], axis=1)
        gf = d / "gfm.parquet"
        if gf.is_file():
            g = pq.read_table(gf, columns=["ix", "iy", "gfm_v"]).to_pandas()
            x = x.merge(g, on=["ix", "iy"], how="left", validate="one_to_one")
        else:
            x["gfm_v"] = np.nan
        q = pd.DataFrame({"ix": t["ix"].to_numpy()[ia].astype(np.int64), "iy": t["iy"].to_numpy()[ia].astype(np.int64)})
        m = q.merge(x, on=["ix", "iy"], how="left", validate="many_to_one")
        m.index = ia
        rows.append(m)
    return pd.concat(rows).sort_index()


def _fe_diff(v: np.ndarray, wet: np.ndarray, aoi: np.ndarray, unit: np.ndarray) -> dict:
    ok = np.isfinite(v)
    v, wet, aoi, unit = v[ok], wet[ok].astype(np.float64), aoi[ok], unit[ok]
    s = pd.DataFrame({"v": v, "w": wet, "a": aoi})
    vm = v - s.groupby("a")["v"].transform("mean").to_numpy()
    wm = wet - s.groupby("a")["w"].transform("mean").to_numpy()
    den = float((wm ** 2).sum())
    if den == 0:
        return {"n": int(ok.sum()), "note": "no within-AOI variation"}
    b = float((wm * vm).sum() / den)
    res = vm - b * wm
    sc = pd.Series(wm * res).groupby(pd.Series(unit).to_numpy()).sum().to_numpy()
    g = len(sc)
    se = math.sqrt(float((sc ** 2).sum()) * g / (g - 1)) / den
    return {"difference": b, "ci95": [b - 1.96 * se, b + 1.96 * se], "n": int(ok.sum()), "n_units": g,
            "mean_water_related": float(v[wet > 0].mean()), "mean_other": float(v[wet == 0].mean())}


def b5_rings() -> dict:
    t = v69.v55.table()
    near = t["P_d1"].to_numpy() <= 500.0
    s = t.loc[near].reset_index(drop=True)
    j = _cell_join(s)
    wet = (s["lex_wetside_n1"].to_numpy() > 0)
    aoi, unit = s["aoi"].to_numpy(), s["unit"].to_numpy()
    vals = {"lc_built": j["lc_built"].to_numpy(np.float64), "lc_tree": j["lc_tree"].to_numpy(np.float64),
            "gfm_not_assessable": 1.0 - j["gfm_v"].to_numpy(np.float64),
            "partly_permanent_water": j["f_perm_water"].to_numpy(np.float64)}
    out: dict = {"n_cells_within_500m": len(s), "n_water_related": int(wet.sum())}
    for k, v in vals.items():
        out[k] = _fe_diff(v, wet, aoi, unit)
    diffs = [out[k] for k in vals if "difference" in out[k]]
    out["all_within_0_05"] = bool(diffs and all(abs(d["difference"]) < 0.05 and -0.05 < d["ci95"][0]
                                                and d["ci95"][1] < 0.05 for d in diffs))
    out["rule"] = ("B-15 B5: capped EMS table, cells within 500 m of the nearest settlement; difference between cells near "
                   "water-related names and other cells in each quantity, AOI fixed effects, 95% intervals clustered by "
                   "location group; reporting rule: bounds wording only if every |difference| < 0.05 with its interval "
                   "within +-0.05")
    return out


# %% B6: GFM の診断
def _auc(y: np.ndarray, p: np.ndarray) -> float | None:
    return v70._auc(y, p) if 0 < y.sum() < len(y) else None


def _csi_matched(y: np.ndarray, p: np.ndarray) -> tuple[int, int]:
    k = int(y.sum())
    if k == 0:
        return 0, 0
    top = np.argsort(-p, kind="stable")[:k]
    tp = int(y[top].sum())
    return tp, 2 * k - tp


def _within_scores(y: np.ndarray, p: np.ndarray, aoi: np.ndarray) -> dict:
    aucs, csis, tp_all, den_all = [], [], 0, 0
    for ia in pd.Series(np.arange(len(y))).groupby(pd.Series(aoi).to_numpy()).groups.values():
        ia = np.asarray(ia)
        yy, pp = y[ia], p[ia]
        a = _auc(yy, pp)
        if a is not None:
            aucs.append(a)
            tp, den = _csi_matched(yy, pp)
            csis.append(tp / den)
            tp_all, den_all = tp_all + tp, den_all + den
    return {"n_aoi": len(aucs), "median_auc": float(np.median(aucs)) if aucs else None,
            "median_csi": float(np.median(csis)) if csis else None,
            "pooled_csi": tp_all / den_all if den_all else None}


def _gfm_day_class() -> dict[str, bool]:
    summ = _json(v69.GFM_SUMMARY).get("aois") or {}
    out = {}
    for a in summ:
        p = v69.LABELS / a / "gfm_meta.json"
        if not p.is_file():
            continue
        m = _json(p)
        if not m.get("window") or not m.get("scenes"):
            continue
        d0 = dt.date.fromisoformat(m["window"][0]) + dt.timedelta(days=1)
        d1 = dt.date.fromisoformat(m["window"][1]) - dt.timedelta(days=1)
        ems_days = {d0 + dt.timedelta(days=i) for i in range((d1 - d0).days + 1)}
        days = {dt.datetime.fromisoformat(x["datetime"].replace("Z", "+00:00")).date() for x in m["scenes"]}
        out[a] = bool(days & ems_days)
    return out


def _within_subset(run: str, seed: int, keep: set[str]) -> dict:
    """v66 の within_placebo（AOI の中だけの評価、判定表）を AOI の部分集合で。"""
    v66 = v70._load("v66_b10_local")
    v63 = v66.v63
    d = v36.load(run, "primary")
    m = np.isin(np.asarray(d["aoi"]).astype(str), list(keep))
    y = np.asarray(d["y"], dtype=np.float64)[m]
    aoi = pd.factorize(pd.Series(np.asarray(d["aoi"])[m]))[0]
    n_a = np.bincount(aoi).astype(np.float64)
    rate = np.bincount(aoi, weights=y) / n_a
    inner = (rate > 0) & (rate < 1)
    inv, cnt = v63._boot_counts(np.asarray(d["unit"])[m], np.random.default_rng(seed))
    g = cnt.shape[1]
    sse = {}
    for a in ("M0", "M1", "M1c", "M3n", "M3np"):
        p = np.clip(np.asarray(d[f"p_{a}"], dtype=np.float64)[m], 1e-7, 1 - 1e-7)
        q = v63._shift_to_rate(np.log(p / (1 - p)), aoi, rate, n_a)
        q = np.where(inner[aoi], q, rate[aoi])
        sse[a] = np.bincount(inv, weights=(q - y) ** 2, minlength=g)
    w = np.vstack([np.ones(g), cnt])
    den = w @ sse["M0"] - w @ sse["M1"]
    r = (w @ sse["M1c"] - w @ sse["M3n"]) / den
    rp = (w @ sse["M1c"] - w @ sse["M3np"]) / den
    dif = r - rp
    return {"n_aoi": len(np.unique(np.asarray(d["aoi"])[m])), "n_units": g, "R": float(r[0]),
            "ci95_R": np.quantile(r[1:], [0.025, 0.975]).tolist(), "R_placebo": float(rp[0]),
            "R_minus_placebo": float(dif[0]), "ci95_R_minus_placebo": np.quantile(dif[1:], [0.025, 0.975]).tolist(),
            "yardstick_ci95": np.quantile(den[1:] / (w[1:] @ np.bincount(inv, minlength=g).astype(np.float64)),
                                          [0.025, 0.975]).tolist()}


def b6_gfm() -> dict:
    v41 = v70._load("v41_gfm_independent_labels")
    t = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    has_c2 = v41.attach_c2(t)
    y = t["y_gfm"].to_numpy(np.float64)
    unit = t["unit"].to_numpy()
    out: dict = {"n_cells": len(t), "n_aoi": int(t["aoi"].nunique())}
    sq = {a: (t[f"p_{a}"].to_numpy(np.float64) - y) ** 2 for a in ("M0", "M1")}
    units, n, e0, e1 = v73._gsum(unit, sq["M0"], sq["M1"])
    w = v73._boot_w(np.random.default_rng(SEED), len(units))
    db = (w @ e0 - w @ e1) / (w @ n)
    out["a_ems_trained_on_gfm_labels"] = {"yardstick": float((e0.sum() - e1.sum()) / n.sum()),
                                          "ci95": np.quantile(db, [0.025, 0.975]).tolist(),
                                          "brier_M0": float(e0.sum() / n.sum()), "brier_M1": float(e1.sum() / n.sum())}
    aoi = t["aoi"].to_numpy(dtype=str)
    out["b_within_aoi"] = {"ems_trained_M1_on_gfm": _within_scores(y, t["p_M1"].to_numpy(np.float64), aoi),
                           "ems_trained_M1_on_ems": _within_scores(t["y"].to_numpy(np.float64),
                                                                    t["p_M1"].to_numpy(np.float64), aoi)}
    if has_c2:
        out["b_within_aoi"]["gfm_trained_M1_on_gfm"] = _within_scores(y, t["p_c2_M1"].to_numpy(np.float64), aoi)
    day = _gfm_day_class()
    same = {a for a, v in day.items() if v}
    other = {a for a, v in day.items() if not v}
    out["c_within_by_scene_day"] = {"same_day": _within_subset("world_v2_gfm", SEED, same),
                                    "other": _within_subset("world_v2_gfm", SEED, other),
                                    "n_aoi_same_day": len(same), "n_aoi_other": len(other)}
    out["rule"] = ("B-15 B6: (a) yardstick of the EMS-trained confirmatory M0 and M1 on the GFM labels of the cells GFM "
                   "can assess (4,000 location-group draws); (b) median within-AOI AUC and CSI with the mapped area "
                   "matched to the flooded area; (c) the within-AOI GFM evaluation of B-10 B10 split by whether a GFM "
                   "scene fell on a day of the EMS images")
    return out


# %% B7: 読みの表
_TOK = re.compile(r"[^a-z0-9]+")


def _settlements() -> tuple[pd.DataFrame, pd.Series]:
    u = v73._analysed()
    gn = v73._gn()
    gids = np.unique(u["P_gid1"].to_numpy())
    s = gn.reindex(gids)
    s = s[s["asciiname"].notna()].copy()
    cells = u.groupby("P_gid1").size()
    s["cells"] = cells.reindex(s.index).fillna(0).astype(int).to_numpy()
    from topsem.world.dataset import continent_map

    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    s["continent"] = s["cc"].map(lambda c: cont.get(str(c), "??"))
    return s, cells


def _wikidata_named_after(gids: list[int]) -> dict:
    """解析した集落（GeoNames の ID）のうち P138（名の由来）を持つ項目の、由来の項目の類。3,000 件ずつ照会する。"""
    if WIKIDATA_CACHE.is_file() and _json(WIKIDATA_CACHE).get("ok"):
        return _json(WIKIDATA_CACHE)
    cls = {"Q5": "person", "Q4022": "water", "Q47521": "water", "Q355304": "water", "Q23397": "water",
           "Q486972": "place", "Q515": "place", "Q3957": "place", "Q532": "place", "Q5084": "place",
           "Q1549591": "place", "Q15284": "place", "Q8502": "mountain", "Q54050": "mountain"}
    vals = " ".join(f"wd:{q}" for q in cls)
    rows = []
    for i in range(0, len(gids), 3000):
        part = " ".join(f'"{g}"' for g in gids[i:i + 3000])
        q = ("SELECT DISTINCT ?gn ?cls WHERE { VALUES ?gn { " + part + " } ?item wdt:P1566 ?gn ; wdt:P138 ?t . "
             "?t wdt:P31 ?cls . VALUES ?cls { " + vals + " } }")
        cmd = ["curl", "-sS", "--max-time", "120", "-X", "POST", "https://query.wikidata.org/sparql",
               "--data-urlencode", f"query={q}", "-H", "Accept: application/sparql-results+json",
               "-H", "User-Agent: toponym-semantics-research/1.0 (verification v76)"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=150, check=False)
            js = json.loads(r.stdout, strict=False)
            rows += [{"gn": b["gn"]["value"], "kind": cls[b["cls"]["value"].rsplit("/", 1)[-1]]}
                     for b in js["results"]["bindings"]]
        except (json.JSONDecodeError, KeyError, subprocess.TimeoutExpired) as exc:
            res = {"ok": False, "error": f"batch {i}: {type(exc).__name__}: {str(exc)[:200]}"}
            WIKIDATA_CACHE.write_text(json.dumps(res), encoding="utf-8")
            return res
        time.sleep(1.0)
    res = {"ok": True, "n_rows": len(rows), "rows": rows, "n_gids_queried": len(gids), "batch": 3000,
           "classes": cls, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    WIKIDATA_CACHE.write_text(json.dumps(res), encoding="utf-8")
    return res


PERSON = {"Q5"}
WATER = {"Q4022", "Q47521", "Q355304", "Q23397"}
PLACE = {"Q486972", "Q515", "Q3957", "Q532", "Q5084", "Q1549591", "Q15284"}
MOUNTAIN = {"Q8502", "Q54050"}


def b7_readings() -> dict:
    lex = Lexicon(v73.LEX)
    s, _ = _settlements()
    blocks_of = {cc: list(v) for cc, v in lex.country_langs.items()}
    hits = []
    for r in s.itertuples():
        ms = v73._matches(lex, r.asciiname, r.cc)
        hits.append(ms)
    s["forms"] = [sorted({(e.lang, e.form, e.cat, e.pos) for e in ms}) for ms in hits]
    s["wet"] = [any(e.cat in WET_SIDE for e in ms) for ms in hits]
    s["any"] = [bool(ms) for ms in hits]
    out: dict = {"n_settlements": len(s), "n_cells": int(s["cells"].sum())}
    # (i) ブロックごと
    bl: dict = {}
    for lang, es in sorted(lex.by_lang.items()):
        countries = sorted(cc for cc, ls in blocks_of.items() if lang in ls)
        ins = s["cc"].isin(countries).to_numpy()
        x = s.loc[ins]
        fx = x.assign(f=x["forms"]).explode("f")
        fx = fx[fx["f"].notna()]
        fx = fx.loc[np.array([f[0] == lang for f in fx["f"]], dtype=bool)]
        wetf = fx.loc[np.array([f[2] in WET_SIDE for f in fx["f"]], dtype=bool)]
        top = wetf.groupby(wetf["f"].map(lambda f: f[1]))["cells"].sum().sort_values(ascending=False).head(5)
        bl[lang] = {"countries": countries, "n_entries": len(es),
                    "entries_by_class": pd.Series([e.cat for e in es]).value_counts().to_dict(),
                    "entries_by_position": pd.Series([e.pos for e in es]).value_counts().to_dict(),
                    "n_forms_matched": int(fx["f"].map(lambda f: f[1]).nunique()) if len(fx) else 0,
                    "n_settlements": int(ins.sum()), "n_cells": int(x["cells"].sum()),
                    "n_settlements_matched": int(fx.index.nunique()) if len(fx) else 0,
                    "n_settlements_water_related": int(wetf.index.nunique()) if len(wetf) else 0,
                    "cells_water_related": int(x.loc[x.index.isin(wetf.index), "cells"].sum()) if len(wetf) else 0,
                    "top5_water_forms_by_cells": {str(k): int(v) for k, v in top.items()}}
    out["i_blocks"] = bl
    # (ii) 地名辞典（国ごと、name と asciiname が違う割合）
    s["non_ascii"] = (s["name"].astype(str) != s["asciiname"].astype(str))
    out["ii_name_differs_by_country"] = {str(c): {"n": len(x), "share_name_differs": float(x["non_ascii"].mean()),
                                                  "share_water_related": float(x["wet"].mean())}
                                         for c, x in s.groupby("cc") if len(x) >= 200}
    # (iii) N-read の類
    nr = pd.read_parquet(DATA_DIR / "global" / "nread" / "nread.parquet", columns=["gid", "category"]).set_index("gid")
    s["nread"] = nr["category"].reindex(s.index).fillna("missing").to_numpy()
    s["block"] = s["cc"].map(lambda c: "+".join(blocks_of.get(str(c), ("none",))))

    def shares(x: pd.DataFrame) -> dict:
        w = x.groupby("nread")["cells"].sum()
        return {"n_settlements": len(x), "n_cells": int(x["cells"].sum()),
                **{str(k): float(v / w.sum()) for k, v in w.items()}} if w.sum() else {"n_settlements": len(x)}
    out["iii_nread_by_continent"] = {str(c): shares(x) for c, x in s.groupby("continent")}
    out["iii_nread_by_block"] = {str(b): shares(x) for b, x in s.groupby("block") if x["cells"].sum() >= 5000}
    # (iv) 二人の符号化者の類ごとの精度と見落とし
    tc = pd.read_csv(DOCV / "v70_two_coders.csv")
    tc["cats"] = [sorted({e.cat for e in v73._matches(lex, n, c)}) for n, c in
                  zip(tc["name"].astype(str), tc["cc"].astype(str), strict=True)]
    iv: dict = {}
    for who in ("coder1", "coder2"):
        d = tc[tc[f"status_{who}"] == "coded"].copy()
        d["codes"] = d[f"codes_{who}"].fillna("").map(lambda x: set(str(x).split(";")) - {""})
        rows = {}
        for c in sorted({c for cs in d["cats"] for c in cs} | set(WET_SIDE)):
            rd = d[d["cats"].map(lambda cs, c=c: c in cs)]
            k, n = int(rd["codes"].map(lambda cs, c=c: c in cs).sum()), len(rd)
            nott = d[d["stratum"] == "not"]
            km, nm = int(nott["codes"].map(lambda cs, c=c: c in cs).sum()), len(nott)
            rows[c] = {"read_as_class_n": n, "recorded_same_class": k, "precision": k / n if n else None,
                       "precision_ci95": v73._wilson(k, n), "not_read_n": nm, "recorded_class_among_not_read": km,
                       "miss_rate": km / nm if nm else None}
        iv[who] = rows
    out["iv_class_precision"] = iv
    # (v) ローマ字化と規則の点検
    def words(x: pd.Series) -> pd.Series:
        return x.str.lower().map(lambda v: [w for w in _TOK.split(v) if w])
    s["words"] = words(s["asciiname"].astype(str))
    chk: dict = {}
    for cc, toks in (("KH", ("phnum", "stoeng", "phumi", "phnom", "stung", "steung", "phum", "kampong", "kompong",
                             "kaoh", "koh", "prek", "preaek")),
                     ("IR", ("pol", "darreh", "rudkhaneh", "darra", "dara", "pul", "ab", "chah", "cheshmeh", "rud")),
                     ("VN", ("cau", "con", "song", "ben", "bai", "dam")),
                     ("PK", ("chah", "ab", "pul", "nala", "nullah")),
                     ("MM", ("myit", "kyun", "chaung", "ywa"))):
        x = s[s["cc"] == cc]
        chk[cc] = {"n_settlements": len(x), **{t: {"n": int(x["words"].map(lambda ws, t=t: t in ws).sum()),
                                                   "n_read_water_related": int(x.loc[x["words"].map(
                                                       lambda ws, t=t: t in ws), "wet"].sum())} for t in toks}}
    fi = s[s["cc"] == "FI"]
    chk["FI"] = {"n_settlements": len(fi),
                 "ala_prefix": int(fi["words"].map(lambda ws: any(w.startswith("ala") for w in ws)).sum()),
                 "yli_prefix": int(fi["words"].map(lambda ws: any(w.startswith("yli") for w in ws)).sum())}
    ham = s[s["words"].map(lambda ws: any(w.endswith(("ham", "hamm")) and len(w) > 4 for w in ws))]
    chk["ham_hamm_by_country"] = {str(k): int(v) for k, v in ham["cc"].value_counts().head(12).items()}
    sur = s[s["words"].map(lambda ws: "sur" in ws)]
    chk["sur_word_by_country"] = {str(k): int(v) for k, v in sur["cc"].value_counts().head(12).items()}
    out["v_checks"] = chk
    # (vi) Wikidata の名の由来
    wd = _wikidata_named_after(sorted(int(g) for g in s.index))
    if wd.get("ok"):
        rw = pd.DataFrame(wd["rows"])
        rw["gid"] = pd.to_numeric(rw["gn"], errors="coerce")
        rw = rw[rw["gid"].isin(s.index)]
        order = {"person": 0, "water": 1, "place": 2, "mountain": 3}
        kind = rw.groupby("gid")["kind"].agg(lambda c: min(c, key=lambda k: order[k]))
        x = s.loc[kind.index].assign(kind=kind.to_numpy())
        out["vi_named_after"] = {"ok": True, "retrieved_utc": wd.get("retrieved_utc"),
                                 "by_kind": {str(k): {"n": len(g), "n_read_water_related": int(g["wet"].sum()),
                                                      "cells_water_related": int(g.loc[g["wet"], "cells"].sum())}
                                             for k, g in x.groupby("kind")}}
    else:
        out["vi_named_after"] = {"ok": False, "error": wd.get("error")}
    out["rule"] = ("B-15 B7: settlements nearest to an evaluable cell of the uncapped table (cells counted per "
                   "settlement); frozen lexicon; blocks by country; N-read classes weighted by cells; class precision "
                   "and misses of the two published coders; romanisation and rule checks by word; Wikidata P138 targets "
                   "by class (person, water feature, place, mountain)")
    return out


# %% B8: 目録の拡張
def _amend_all(v: int, name: str) -> str:
    if v <= 67:
        return v68._amend(v)
    if v == 68:
        return "11"
    if v == 69:
        return "12"
    if v in (70, 71, 72):
        return "13"
    if v in (73, 74, 75):
        return "14"
    return "15"


def b8_census_all() -> dict:
    rows: list[dict] = []
    for f in sorted(DOCV.glob("v[3-7][0-9]*_*.json")):
        m = re.match(r"v(\d+)", f.name)
        v = int(m.group(1)) if m else 0
        if not 36 <= v <= 76 or "estimates" in f.name:
            continue
        v69._walk_fixed(_json(f), [], rows, f.name, _amend_all(v, f.name))
    df = pd.DataFrame(rows)
    df["key"] = (df["kind"] + "|" + df["estimate"].round(6).astype(str) + "|" + df["ci95_lo"].round(6).astype(str)
                 + "|" + df["ci95_hi"].round(6).astype(str))
    df["duplicate"] = df.duplicated("key")
    sep = (df["kind"] == "odds_ratio") & ((np.abs(np.log(df["estimate"].clip(lower=1e-300))) > 10)
                                          | (df["ci95_hi"] / df["ci95_lo"].clip(lower=1e-300) > 1e6))
    df["valid"] = np.isfinite(df["ci95_lo"]) & np.isfinite(df["ci95_hi"]) & (df["ci95_lo"] <= df["ci95_hi"]) & ~sep
    df["family"] = [v68._family(p, k) for p, k in zip(df["path"], df["kind"], strict=True)]
    base = df.loc[~df["duplicate"] & df["valid"] & df["name_effect"]].copy()
    out: dict = {"n_rows": len(df), "n_unique_valid_name_effect": len(base), "families": {}}
    qcols = []
    for fam, x in base.groupby("family"):
        xp = x[x["p_approx"].notna()].copy()
        qbh, qby = v69._bh_by(xp["p_approx"].to_numpy()) if len(xp) else (np.array([]), np.array([]))
        xp["q_bh_family"], xp["q_by_family"] = qbh, qby
        qcols.append(xp[["file", "path", "q_bh_family", "q_by_family"]])
        ex = x[x["excludes_null"]].copy()
        ex["direction"] = ex.apply(v73._direction, axis=1)
        out["families"][str(fam)] = {
            "n": len(x), "n_excluding_null": len(ex), "by_amendment": {str(k): int(v) for k, v in
                                                                      x["amendment"].astype(str).value_counts().items()},
            "n_bh_q_below_0_05": int((qbh < 0.05).sum()), "n_by_q_below_0_05": int((qby < 0.05).sum()),
            "n_by_q_below_0_10": int((qby < 0.10).sum()),
            "direction_of_null_excluding": {str(k): int(v) for k, v in ex["direction"].value_counts().items()}}
    df = df.merge(pd.concat(qcols), on=["file", "path"], how="left") if qcols else df
    df.drop(columns="key").to_csv(DOCV / "v76_estimates.csv", index=False)
    nb = base[base["near_settlement"] & (base["family"] == "flood")]
    out["near_settlement_or"] = {"n": len(nb), "n_excluding_one": int(nb["excludes_null"].sum())}
    out["rule"] = ("B-15 B8: the rules of the v69 census applied to every output v36 to v76 (amendments 2 to 15); "
                   "duplicates by estimate and interval; BH and BY within each family; direction of the estimates "
                   "that exclude the null by the rules of B-14 B10; file docs/verification/v76_estimates.csv")
    return out


# %% B10: 集落の近くの推定の推測
def _fit_scores(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, groups: np.ndarray) -> dict:
    """v39.fit と同じ推定（AOI 固定効果 + z）で、係数・Schur 補行列の逆・場所グループごとのスコアの和を返す。"""
    from scipy import sparse

    ua, ia = np.unique(aoi, return_inverse=True)
    n, k, a = len(y), z.shape[1], len(ua)
    d = sparse.csr_matrix((np.ones(n), (np.arange(n), ia)), shape=(n, a))
    rate = np.bincount(ia, weights=y) / np.bincount(ia)
    alpha = np.log(rate / (1 - rate))
    beta = np.zeros(k)
    for _ in range(50):
        p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta)))
        w = p * (1 - p)
        r = y - p
        haa = np.bincount(ia, weights=w)
        haz = d.T @ (z * w[:, None])
        hzz = z.T @ (z * w[:, None])
        ga, gz = np.bincount(ia, weights=r), z.T @ r
        s = hzz - haz.T @ (haz / haa[:, None])
        db = np.linalg.solve(s, gz - haz.T @ (ga / haa))
        da = (ga - haz @ db) / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max()) < 1e-10:
            break
    p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta)))
    w = p * (1 - p)
    haa = np.bincount(ia, weights=w)
    haz = d.T @ (z * w[:, None])
    hzz = z.T @ (z * w[:, None])
    sinv = np.linalg.inv(hzz - haz.T @ (haz / haa[:, None]))
    zc = z - (haz / haa[:, None])[ia]
    ug, ig = np.unique(groups, return_inverse=True)
    u = np.zeros((len(ug), k))
    np.add.at(u, ig, zc * (y - p)[:, None])
    return {"beta": beta, "sinv": sinv, "u": u, "G": len(ug)}


WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5), math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def _wild(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, groups: np.ndarray, j: int, rng: np.random.Generator) -> dict:
    f = _fit_scores(y, z, aoi, groups)
    c = f["u"] @ f["sinv"][j]
    v = WEBB[rng.integers(0, 6, size=(B_WILD, f["G"]))]
    bstar = f["beta"][j] + v @ c
    g = f["G"]
    se_cr = math.sqrt(float(c @ c) * g / (g - 1))
    tq = float(stats.t.ppf(0.975, g - 1))
    return {"odds_ratio": math.exp(f["beta"][j]), "G": g,
            "wild_ci95": np.exp(np.quantile(bstar, [0.025, 0.975])).tolist(),
            "cr_normal_ci95": [math.exp(f["beta"][j] - 1.96 * se_cr), math.exp(f["beta"][j] + 1.96 * se_cr)],
            "cr_t_ci95": [math.exp(f["beta"][j] - tq * se_cr), math.exp(f["beta"][j] + tq * se_cr)]}


_P: dict = {}


def _perm_init(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, groups: np.ndarray, j: int) -> None:
    _threads1()
    _P.update(y=y, z=z, aoi=aoi, groups=groups, j=j)


def _perm_rep(b: int) -> float:
    rng = np.random.default_rng(SEED + 500_000 + b)
    y, aoi = _P["y"], _P["aoi"]
    s = pd.Series(np.arange(len(y)))
    yp = y.copy()
    for ix in s.groupby(aoi).groups.values():
        ix = np.asarray(ix)
        yp[ix] = y[rng.permutation(ix)]
    return float(v39.fit(yp, _P["z"], aoi, _P["groups"])["beta"][_P["j"]])


def b10_near_inference() -> dict:
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    out: dict = {}
    for key, ycol in (("planned_gfm", "y_gfm"), ("ems_same_cells", "y")):
        y = g[ycol].to_numpy(np.float64)
        rate = pd.Series(y).groupby(g["aoi"].to_numpy()).transform("mean").to_numpy()
        s = g.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
        ys = s[ycol].to_numpy(np.float64)
        z, names = v39.design(s, ("wetside",), True, banded=True)
        j = names.index("wetside:le500m")
        aoi, unit = s["aoi"].to_numpy(), s["unit"].to_numpy(str)
        res = _wild(ys, z, aoi, unit, j, np.random.default_rng(SEED))
        cache = VER / f"v76_b10_perm_{key}.npy"
        if cache.is_file():
            bp = np.load(cache)
        else:
            t0 = time.time()
            with ProcessPoolExecutor(N_WORKERS, initializer=_perm_init, initargs=(ys, z, aoi, unit, j)) as ex:
                bp = np.array(list(ex.map(_perm_rep, range(B_PERM), chunksize=2)), dtype=np.float64)
            np.save(cache, bp)
            print(f"  B10 {key}: {len(bp)} permutations ({time.time() - t0:.0f} s)", flush=True)
        res["label_permutation"] = {"n": len(bp), "median_or": float(np.exp(np.median(bp))),
                                    "q025_q975_or": np.exp(np.quantile(bp, [0.025, 0.975])).tolist(),
                                    "share_below_observed": float((bp < math.log(res["odds_ratio"])).mean())}
        res["n_cells"], res["n_flooded_near"] = len(s), int(((s["P_d1"].to_numpy() <= 500)
                                                            & (s["lex_wetside_n1"].to_numpy() > 0) & (ys == 1)).sum())
        out[key] = res
    out["rule"] = ("B-15 B10: planned GFM estimate and the EMS estimate on the same cells (v41 table, frozen reading, AOI "
                   "fixed effects, distance bands, logit of M1); (a) unrestricted wild cluster score bootstrap over "
                   "location groups with Webb six-point weights (9,999 draws) and CR intervals with normal and t(G-1); "
                   "(b) 200 permutations of the labels within AOIs with names fixed")
    return out


# %% B11: 点検
def b11_checks() -> dict:
    out: dict = {}
    st = ((_json(DOCV / "v36_yardstick_robustness.json").get("splits") or {}).get("primary") or {}).get("strata") or {}
    out["language_strata"] = {k: _scalars(v) for k, v in st.items() if k.startswith("language_")
                              and (v.get("n_cells") or 0) >= 20000}
    hits = []
    for f in sorted(DOCV.glob("v[3-7][0-9]*_*.json")):
        txt = f.read_text(encoding="utf-8")
        hits += [f"{f.name}:{m.group(0)}" for m in re.finditer(r'"[^"]*(?:tercile|magnitude)[^"]*"', txt, re.IGNORECASE)][:5]
    out["magnitude_tercile_keys"] = hits
    pa = pd.read_csv(VER / "v73_b15_per_aoi.csv")
    x = pa["delta_auc"].dropna()
    out["auc_share_up"] = {"n_with_both_classes": len(x), "share_positive": float((x > 0).mean()),
                           "share_nonnegative": float((x >= 0).mean()), "n_aoi_total": len(pa),
                           "share_positive_of_all_aoi": float((pa["delta_auc"] > 0).sum() / len(pa))}
    out["rule"] = "B-15 B11: checks of numbers named by the eleventh critiques"
    return out


# %%
PARTS = {"B1": b1_prespecified, "B2": b2_few_clusters, "B3": b3_supercluster, "B4": b4_paired, "B5": b5_rings,
         "B6": b6_gfm, "B7": b7_readings, "B8": b8_census_all, "B9": b9_norway, "B10": b10_near_inference,
         "B11": b11_checks}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="例: B1,B2")
    a = ap.parse_args()
    want = [p for p in (a.only.split(",") if a.only else PARTS) if p]
    res = _json(OUT) if OUT.is_file() else {}
    for k in want:
        t0 = time.time()
        print(f"{k} ...", flush=True)
        res[k] = PARTS[k]()
        res[k]["seconds"] = round(time.time() - t0)
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        print(f"{k} ({res[k]['seconds']} s) {json.dumps(res[k], ensure_ascii=False, default=float)[:1500]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
