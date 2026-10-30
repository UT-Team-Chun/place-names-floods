# %% v77（事後・amendment B-16 の B1〜B10 と C2）: 初見の査読 第 12 回に応じた点検
"""確証判定は変えない。B-16 と追記 1 で計算の前に固定した量を、二段目の保存・held-out の予測・表から計算する。

- B1: 指定した事後の構成。主分割の場所グループを同じ activation で束ねた super-cluster を 500 回引き、二段目を
  重みつきで当て直す（凍結した読み）。C2 の帯の模型（M1cd → M3nd、M1sd → M3db）と B9 (e) の水名の腕も同じ引きで
- B2: Δ の分散への寄与の大きい五つの場所グループを除いた当て直し
- B3: 少数クラスタの参照（θ_U と ν を交互に収束させ t(ν − 1)）を、保存した再標本のある全ての構成に
- B4: 上限の注入の OR
- B5: 集落の近くの推定の、帰無を課した wild cluster score bootstrap と検定の反転の区間
- B6: ラベルの層・目的・恒久水・距離帯（直した実行、モデル固定）
- B7: GFM と EMS の一致
- B8: 集落の近くの OR の誤分類の偏りの解析
- B9: 読みの点検（B-7 の置換の placebo、ブロック・国・上位の場所グループ、部分集合、符号化者の κ、水名は B1）
- C2: B-10 の帯の模型（world_v2_b10stk）の帯ごとの R

出力: stdout と docs/verification/v77_b16_local.json（部分ごとに足す。`--only B1,B2` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import os
import pickle
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse, stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import stage2_weighted
from topsem.world.b10 import band_indicators, interactions
from topsem.world.b13 import SESOI_R, fieller_interval
from topsem.world.lexicon import CATEGORIES, WET_SIDE, Lexicon

OUT = ROOT / "docs" / "verification" / "v77_b16_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
LABELS = DATA_DIR / "global" / "labels"
SEED = 20261105
B_SC = 500
B_FIXED = 4000
B_WILD = 9999
Z90 = float(stats.norm.ppf(0.95))  # 両側 90% = 片側 95%
N_WORKERS = int(os.environ.get("V77_WORKERS", "6"))
DUMP_RUN = "world_v2_b16_stage2_local"
CORR_RUN = "world_v2_lexfix"
B10_RUN = "world_v2_b10stk"
PERM_RUN = "world_v2_b7perm"
BANDS = ("all", "b500", "b1k", "bfar")
NEAR_M = 500.0
FROZEN_LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
GRID = np.round(np.arange(-1.6, 1.6 + 1e-9, 0.01), 2)  # B5: log OR の格子
SE_GRID = (0.6, 0.8, 1.0)  # B8
SP_GRID = (0.99, 0.995, 1.0)


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v76 = _load("v76_b15_local")
v73, v70, v69, v68, v36, v39 = v76.v73, v76.v70, v76.v69, v76.v68, v76.v36, v76.v39


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[k] = "1"


def _band_index(d1: np.ndarray) -> np.ndarray:
    """0: 500 m 以内、1: 0.5〜1 km、2: 1 km 超（B-10 の帯と同じ境）。"""
    return band_indicators(np.asarray(d1, dtype=np.float64)).argmax(axis=1).astype(np.int8)


# %% 共通: B3 の少数クラスタの参照の規則（追記 1 の 5）
def _nu(x: np.ndarray, n: np.ndarray) -> float:
    """単位の寄与 x（和）と行の数 n からの Satterthwaite の有効自由度 (Σγ)² / Σγ²（v76 の _nu と同じ）。"""
    return v76._nu(np.asarray(x, dtype=np.float64), np.asarray(n, dtype=np.float64))


def _nu_same(xe: np.ndarray, xd: np.ndarray, n: np.ndarray):
    """L+(θ) = Δ − θD（同じセルの平均）の単位の寄与 xe − θ xd からの ν(θ)。"""
    return lambda th: _nu(xe - th * xd, n)


def _nu_mixed(key: np.ndarray, e_cell: np.ndarray, d_cell: np.ndarray, m: np.ndarray):
    """帯の Δ（帯のセルの平均）と全体の D（全てのセルの平均）の L+(θ) の、単位の影響の二乗からの ν(θ)。"""
    _, inv = np.unique(np.asarray(key), return_inverse=True)
    g = int(inv.max()) + 1
    eg = np.bincount(inv, weights=np.where(m, e_cell, 0.0), minlength=g)
    nb = np.bincount(inv, weights=m.astype(np.float64), minlength=g)
    dg = np.bincount(inv, weights=d_cell, minlength=g)
    ng = np.bincount(inv, minlength=g).astype(np.float64)
    big_nb, big_n = nb.sum(), ng.sum()
    e_hat, d_hat = eg.sum() / big_nb, dg.sum() / big_n

    def f(th: float) -> float:
        gam = ((eg - e_hat * nb) / big_nb - th * (dg - d_hat * ng) / big_n) ** 2
        return float(gam.sum() ** 2 / (gam ** 2).sum())

    return f


def _nu_rule(e0: float, d0: float, var: tuple[float, float, float], nu_of, g: int) -> dict:
    """θ_U と ν(θ_U) を交互に収束させた t(ν − 1) の片側 95% 上限（Fieller の上端）と、同じ ν の z・t(ν)。"""
    out: dict = {"G": g, "G_star_nu_delta": nu_of(0.0)}
    fz = fieller_interval(e0, d0, var, z=Z90)
    out["z"] = {"fieller90": list(fz) if fz else None, "upper_R": fz[1] if fz else None}
    if fz is None:
        out["t_nu_minus_1"] = {"upper_R": None, "reason": "Fieller interval unbounded with z"}
        return out
    theta, nu, it = fz[1], float("nan"), 0
    reason = "converged"
    for it in range(1, 101):
        nu = nu_of(theta)
        if nu - 1 < 1:
            reason = "nu - 1 < 1"
            break
        fi = fieller_interval(e0, d0, var, z=float(stats.t.ppf(0.95, nu - 1)))
        if fi is None:
            reason = "Fieller interval unbounded"
            break
        done = abs(fi[1] - theta) < 1e-6
        theta = fi[1]
        if done:
            break
    else:
        reason = "not converged in 100 iterations"
    out["nu_l_plus_at_upper"] = nu
    out["iterations"] = it
    if reason in ("converged", "not converged in 100 iterations"):
        out["t_nu_minus_1"] = {"critical_value": float(stats.t.ppf(0.95, nu - 1)), "upper_R": theta,
                               "reason": reason}
    else:
        out["t_nu_minus_1"] = {"upper_R": None, "reason": reason}
    if nu > 0:
        ft = fieller_interval(e0, d0, var, z=float(stats.t.ppf(0.95, nu)))
        out["t_nu"] = {"critical_value": float(stats.t.ppf(0.95, nu)), "upper_R": ft[1] if ft else None}
    return out


def _cov3(eb: np.ndarray, db: np.ndarray) -> tuple[float, float, float]:
    ok = np.isfinite(eb) & np.isfinite(db)
    c = np.cov(eb[ok], db[ok])
    return float(c[0, 0]), float(c[1, 1]), float(c[0, 1])


# %% B1: 指定した事後の構成（super-cluster、二段目を当て直す）
ARMS_B1 = ("M1c", "M3n", "M1s", "M3d", "M1cd", "M3nd", "M1sd", "M3db", "M1cH", "M3nH", "M3dH")
COMPS_B1 = {"M1c_M3n": ("M1c", "M3n"), "M1s_M3d": ("M1s", "M3d"), "M1cd_M3nd": ("M1cd", "M3nd"),
            "M1sd_M3db": ("M1sd", "M3db"), "M1cH_M3nH": ("M1cH", "M3nH"), "M1s_M3dH": ("M1s", "M3dH")}
SC_CACHE = VER / "v77_superclusters.json"
POINT_CACHE = VER / "v77_b1_point.npz"
INIT_CACHE = VER / "v77_b1_init.pkl"
REPS_CACHE = VER / "v77_b1_reps.npy"


def _dump_folds() -> list[dict]:
    src = RUNS_DIR / DUMP_RUN / "predictions" / "primary"
    out = []
    for k in range(10):
        f = src / f"fold{k}_stage2.npz"
        if not f.is_file():
            return []
        with np.load(f, allow_pickle=False) as z:
            out.append({kk: z[kk] for kk in z.files})
    return out


def _superclusters() -> dict:
    """解析した AOI の場所グループ（主分割、緩衝 5 km）と、同じ activation で束ねた super-cluster。"""
    if SC_CACHE.is_file():
        return _json(SC_CACHE)
    import world_eval as ev

    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    dirs, _, infos = ev.select_aois(gn)
    names = sorted(d.name for d in dirs)
    grp = ev.location_groups({a: infos[a].bounds for a in names}, ev.GROUP_BUFFER_KM)
    sc = ev.merge_by_activation(grp, {a: infos[a].activation for a in names})
    out = {"aoi_group": {a: int(grp[a]) for a in names}, "aoi_sc": {a: int(sc[a]) for a in names},
           "aoi_activation": {a: infos[a].activation for a in names},
           "aoi_event": {a: infos[a].event_time for a in names},
           "aoi_untyped": {a: bool(infos[a].unspecified_as_riverine) for a in names}}
    SC_CACHE.write_text(json.dumps(out), encoding="utf-8")
    return out


_HYD: dict = {}


def _hydronym(gids: np.ndarray) -> np.ndarray:
    if "s" not in _HYD:
        t = pd.read_parquet(VER / "v57_hydronym.parquet")
        _HYD["s"] = pd.Series(t["hydronym"].to_numpy(np.float64), index=t["gid"].to_numpy(np.int64))
    return _HYD["s"].reindex(np.asarray(gids, dtype=np.int64)).fillna(0.0).to_numpy()


def _base16(f: dict, *, frozen: bool) -> dict:
    """腕の設計の元（地域構成・名の特徴・中心化した特徴・帯・水名）。凍結した読みはノルウェーの辞書の特徴を 0 に。"""
    x = v76._designs(f, frozen=frozen)
    k = x["M1c"][0].shape[1]
    h = {s: _hydronym(f[f"gid_{s}"])[:, None] for s in ("tr", "te")}
    ch = {s: v76._composition(h[s], f[f"aoi_{s}"]) for s in ("tr", "te")}
    out = {}
    for i, s in enumerate(("tr", "te")):
        out[s] = {"xc": x["M1c"][i].astype(np.float32), "xn": x["M3n"][i][:, k:].astype(np.float32),
                  "dm": x["M3d"][i].astype(np.float32), "bd": band_indicators(f[f"d1_{s}"]),
                  "h": h[s].astype(np.float32), "ch": ch[s].astype(np.float32)}
    return out


def _x(b: dict, arm: str) -> np.ndarray:
    """腕の二段目の説明変数（_base16 の片側から）。"""
    if arm == "M1s":
        return np.zeros((len(b["xc"]), 0), dtype=np.float32)
    parts = {"M1c": ("xc",), "M3n": ("xc", "xn"), "M3d": ("dm",), "M1cd": ("xc", "bd"), "M1sd": ("bd",),
             "M1cH": ("xc", "ch"), "M3nH": ("xc", "ch", "h")}
    if arm in parts:
        return np.column_stack([b[p] for p in parts[arm]])
    if arm == "M3nd":
        return np.column_stack([b["xc"], b["bd"], interactions(b["xn"], b["bd"])])
    if arm == "M3db":
        return np.column_stack([b["bd"], interactions(b["dm"], b["bd"])])
    if arm == "M3dH":
        return b["h"] - b["ch"]
    msg = f"unknown arm {arm}"
    raise KeyError(msg)


def _sc_codes(aoi: np.ndarray, sc: dict) -> np.ndarray:
    m = sc["aoi_sc"]
    return np.array([m[a] for a in np.asarray(aoi).astype(str)], dtype=np.int32)


def _b1_point(folds: list[dict], sc: dict, *, frozen: bool) -> dict:
    """重みの無い当て直し（点）。セルごとの予測と、腕ごとの重みの無い当ての係数（引きの初期値）。"""
    preds = {a: [] for a in ARMS_B1}
    init: list[dict] = []
    keep = {k: [] for k in ("y", "p_M0", "p_M1", "unit", "aoi", "cc", "d1", "gid", "ix", "iy")}
    for f in folds:
        b = _base16(f, frozen=frozen)
        ini = {}
        for a in ARMS_B1:
            p, ini[a] = stage2_weighted(f["z_tr"], f["z_te"], _x(b["tr"], a), _x(b["te"], a), f["y_tr"])
            preds[a].append(p.astype(np.float32))
        init.append(ini)
        for k, src in (("y", "y_te"), ("p_M0", "p_M0_te"), ("p_M1", "p_M1_te"), ("unit", "unit_te"),
                       ("aoi", "aoi_te"), ("cc", "cc_te"), ("d1", "d1_te"), ("gid", "gid_te"), ("ix", "ix_te"),
                       ("iy", "iy_te")):
            keep[k].append(f[src])
    out = {k: np.concatenate(v) for k, v in keep.items()} | {f"p_{a}": np.concatenate(v) for a, v in preds.items()}
    out["sc"] = _sc_codes(out["aoi"], sc)
    out["band"] = _band_index(out["d1"])
    return {"arrays": out, "init": init}


_S: dict = {}


def _b1_init() -> None:
    _threads1()
    folds = _dump_folds()
    sc = _superclusters()
    with INIT_CACHE.open("rb") as fh:
        init = pickle.load(fh)
    labels = sorted(set(sc["aoi_sc"].values()))
    _S["S"] = len(labels)
    pos = {v: i for i, v in enumerate(labels)}
    _S["folds"] = []
    for f, ini in zip(folds, init, strict=True):
        b = _base16(f, frozen=True)
        _S["folds"].append({
            "z_tr": f["z_tr"], "z_te": f["z_te"], "y_tr": f["y_tr"], "y_te": f["y_te"].astype(np.float64),
            "s_tr": np.array([pos[v] for v in _sc_codes(f["aoi_tr"], sc)], dtype=np.int32),
            "s_te": np.array([pos[v] for v in _sc_codes(f["aoi_te"], sc)], dtype=np.int32),
            "band": _band_index(f["d1_te"]), "p0": f["p_M0_te"].astype(np.float64),
            "p1": f["p_M1_te"].astype(np.float64), "b": b, "init": ini})


def _b1_rep(r: int) -> list[float]:
    """一回の引き: 腕・M0・M1 の帯ごと（全て・3 帯）の重みつき二乗誤差の和と、帯ごとの重みの和。"""
    rng = np.random.default_rng(SEED + r)
    w = rng.multinomial(_S["S"], np.full(_S["S"], 1.0 / _S["S"])).astype(np.float64)
    nm = len(ARMS_B1) + 2
    acc = np.zeros((nm + 1, len(BANDS)))
    for f in _S["folds"]:
        w_tr, w_te = w[f["s_tr"]], w[f["s_te"]]
        if w_tr.sum() == 0 or w_te.sum() == 0:
            continue
        y = f["y_te"]
        sq = []
        for a in ARMS_B1:
            p, _ = stage2_weighted(f["z_tr"], f["z_te"], _x(f["b"]["tr"], a), _x(f["b"]["te"], a), f["y_tr"],
                                   w_tr=w_tr, init=f["init"][a])
            sq.append((p - y) ** 2)
        sq += [(f["p0"] - y) ** 2, (f["p1"] - y) ** 2, np.ones(len(y))]
        for i, v in enumerate(sq):
            wv = w_te * v
            acc[i, 0] += wv.sum()
            acc[i, 1:] += np.bincount(f["band"], weights=wv, minlength=3)
    return acc.ravel().tolist()


def _b1_reps() -> np.ndarray:
    if REPS_CACHE.is_file():
        return np.load(REPS_CACHE)
    t0 = time.time()
    with ProcessPoolExecutor(N_WORKERS, initializer=_b1_init) as ex:
        r = np.array(list(ex.map(_b1_rep, range(B_SC), chunksize=1)), dtype=np.float64)
    np.save(REPS_CACHE, r)
    print(f"  B1: {len(r)} draws ({time.time() - t0:.0f} s)", flush=True)
    return r


def _unit_sums(key: np.ndarray, e: np.ndarray, d: np.ndarray, mask: np.ndarray | None = None) -> tuple:
    m = np.ones(len(e), dtype=bool) if mask is None else mask
    _, inv = np.unique(np.asarray(key)[m], return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    return np.bincount(inv, weights=e[m]), np.bincount(inv, weights=d[m]), n


def b1_designated() -> dict:
    folds = _dump_folds()
    if len(folds) < 10:
        return {"note": f"stage-2 dump of {DUMP_RUN} incomplete ({len(folds)} folds)"}
    sc = _superclusters()
    if POINT_CACHE.is_file() and INIT_CACHE.is_file():
        with np.load(POINT_CACHE, allow_pickle=False) as z:
            arr = {k: z[k] for k in z.files}
    else:
        pt = _b1_point(folds, sc, frozen=True)
        arr = pt["arrays"]
        np.savez_compressed(POINT_CACHE, **arr)
        with INIT_CACHE.open("wb") as fh:
            pickle.dump(pt["init"], fh)
    y = arr["y"].astype(np.float64)
    # 実装の確かめ: 単位の符号の一致、確証の実行（world_v2）の予測との差
    grp = sc["aoi_group"]
    unit_ok = bool(all(f"g{grp[a]:05d}" == u for a, u in zip(arr["aoi"].astype(str), arr["unit"].astype(str),
                                                             strict=True)))
    conf = v36.load("world_v2", "primary")
    check: dict = {"unit_codes_match_location_groups": unit_ok}
    if conf is not None and len(conf["y"]) == len(y) and np.array_equal(conf["y"].astype(np.float64), y):
        for a in ("M0", "M1", "M1c", "M3n"):
            pa = arr[f"p_{a}"].astype(np.float64)
            check[f"max_abs_diff_p_{a}"] = float(np.max(np.abs(pa - conf[f"p_{a}"].astype(np.float64))))
            check[f"brier_{a}_refit"] = float(np.mean((pa - y) ** 2))
            check[f"brier_{a}_confirmatory"] = float(np.mean((conf[f"p_{a}"].astype(np.float64) - y) ** 2))
    else:
        check["confirmatory_rows"] = "differ or missing"
    # 直した辞書の読み（点だけ）
    corr = _b1_point(folds, sc, frozen=False)["arrays"]
    bs_corr = {a: float(np.mean((corr[f"p_{a}"].astype(np.float64) - y) ** 2)) for a in ARMS_B1}
    d_all = float(np.mean((arr["p_M0"] - y) ** 2 - (arr["p_M1"] - y) ** 2))
    corrected = {k: {"delta_bs": bs_corr[a] - bs_corr[b], "R": (bs_corr[a] - bs_corr[b]) / d_all}
                 for k, (a, b) in COMPS_B1.items()}
    r = _b1_reps()
    nm = len(ARMS_B1) + 2
    reps = r.reshape(len(r), nm + 1, len(BANDS))
    idx = {a: i for i, a in enumerate((*ARMS_B1, "M0", "M1"))}
    wt = reps[:, nm, :]
    labels, sc_n = np.unique(arr["sc"], return_counts=True)
    out: dict = {"n_superclusters": len(labels), "largest_supercluster_share": float(sc_n.max() / sc_n.sum()),
                 "n_location_groups": len(np.unique(arr["unit"])), "n_cells": len(y), "draws": len(r),
                 "check": check, "reading": "frozen (Norwegian lexicon features set to 0)",
                 "corrected_reading_points": corrected, "comparisons": {}}
    sq = {a: (arr[f"p_{a}"].astype(np.float64) - y) ** 2 for a in (*ARMS_B1, "M0", "M1")}
    d_cell = sq["M0"] - sq["M1"]
    band = arr["band"]
    rng = np.random.default_rng(SEED)
    for key, (a, b) in COMPS_B1.items():
        e_cell = sq[a] - sq[b]
        res = {}
        for j, bn in enumerate(BANDS):
            m = np.ones(len(y), dtype=bool) if bn == "all" else band == (j - 1)
            e0, d0 = float(e_cell[m].mean()), float(d_cell[m].mean())
            eb = (reps[:, idx[a], j] - reps[:, idx[b], j]) / wt[:, j]
            db = (reps[:, idx["M0"], j] - reps[:, idx["M1"], j]) / wt[:, j]
            det = v73._levels_detail(e0, d0, eb, db, rng)
            xe, xd, n = _unit_sums(arr["sc"], e_cell, d_cell, m)
            row = {"n_cells": int(m.sum()), "band_yardstick": det,
                   "few_cluster": _nu_rule(e0, d0, _cov3(eb, db), _nu_same(xe, xd, n), len(n))}
            if bn != "all":  # 全体の物差し（全てのセルの D）で割った帯の R
                d_all_b = (reps[:, idx["M0"], 0] - reps[:, idx["M1"], 0]) / wt[:, 0]
                var2 = _cov3(eb, d_all_b)
                fi = fieller_interval(e0, d_all, var2, z=Z90)
                row["overall_yardstick"] = {
                    "R": e0 / d_all, "fieller90": list(fi) if fi else None,
                    "few_cluster": _nu_rule(e0, d_all, var2, _nu_mixed(arr["sc"], e_cell, d_cell, m),
                                            len(np.unique(arr["sc"])))}
            res[bn] = row
        out["comparisons"][key] = res
    out["rule"] = ("B-16 B1 (addendum 1, items 1-3): stage-2 design of the primary split re-dumped locally "
                   f"({DUMP_RUN}); frozen reading; super-clusters = location groups (5 km) merged by shared "
                   f"activation; {B_SC} multinomial draws of super-clusters (seed {SEED} + b) weighting training "
                   "and evaluation cells, stage 2 refitted with weights from the unweighted fit; base predictions "
                   "fixed; Fieller 90% interval from the draws' covariance; few-cluster rule of B3 with "
                   "super-clusters as units")
    return out


# %% B2: 分散への寄与の大きい五つの場所グループを除いた当て直し
def _point_arrays() -> dict | None:
    if not POINT_CACHE.is_file():
        return None
    with np.load(POINT_CACHE, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def b2_leaveout() -> dict:
    folds = _dump_folds()
    arr = _point_arrays()
    if len(folds) < 10 or arr is None or not INIT_CACHE.is_file():
        return {"note": "B1 point estimate not available (run B1 first)"}
    with INIT_CACHE.open("rb") as fh:
        init = pickle.load(fh)
    y = arr["y"].astype(np.float64)
    e = (arr["p_M1c"] - y) ** 2 - (arr["p_M3n"] - y) ** 2
    u, inv = np.unique(arr["unit"].astype(str), return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    eg = np.bincount(inv, weights=e)
    big_n = n.sum()
    gam = ((eg - eg.sum() / big_n * n) / big_n) ** 2
    top = [str(x) for x in u[np.argsort(-gam)[:5]]]
    cc = pd.Series(arr["cc"].astype(str)).groupby(inv).agg(lambda s: s.mode().iloc[0]).to_numpy()
    info = [{"unit": str(u[i]), "majority_country": str(cc[i]), "n_cells": int(n[i]),
             "share_of_variance": float(gam[i] / gam.sum()), "contribution_delta": float(eg[i] / big_n)}
            for i in np.argsort(-gam)[:5]]
    out: dict = {"top_five": info, "removed": {}}
    for drop in [[t] for t in top] + [top]:
        sse = dict.fromkeys(("M1c", "M3n", "M1s", "M3d", "M0", "M1"), 0.0)
        nn = 0
        for f, ini in zip(folds, init, strict=True):
            b = _base16(f, frozen=True)
            ktr = ~np.isin(f["unit_tr"].astype(str), drop)
            kte = ~np.isin(f["unit_te"].astype(str), drop)
            yt = f["y_te"].astype(np.float64)[kte]
            for a in ("M1c", "M3n", "M1s", "M3d"):
                p, _ = stage2_weighted(f["z_tr"][ktr], f["z_te"][kte], _x(b["tr"], a)[ktr], _x(b["te"], a)[kte],
                                       f["y_tr"][ktr], init=ini[a])
                sse[a] += float(((p - yt) ** 2).sum())
            sse["M0"] += float(((f["p_M0_te"][kte] - yt) ** 2).sum())
            sse["M1"] += float(((f["p_M1_te"][kte] - yt) ** 2).sum())
            nn += int(kte.sum())
        dd = (sse["M0"] - sse["M1"]) / nn
        out["removed"]["+".join(drop)] = {
            "n_cells": nn, "yardstick": dd,
            **{k: {"delta_bs": (sse[a] - sse[b]) / nn, "R": (sse[a] - sse[b]) / nn / dd}
               for k, (a, b) in (("M1c_M3n", ("M1c", "M3n")), ("M1s_M3d", ("M1s", "M3d")))}}
    out["rule"] = ("B-16 B2: the five location groups with the largest gamma_g = ((e_g - Delta n_g)/N)^2 for M1c -> M3n "
                   "(frozen reading, B1 point refit) removed one at a time and all five together from training and "
                   "evaluation; stage 2 refitted (unweighted, started from the full fit); base predictions fixed")
    return out


# %% B3: 少数クラスタの参照の規則を全ての構成に
def b3_few_cluster() -> dict:
    df = v76._unit_contribs()
    coa, _ = v70._country_of_activation()
    df["country_of_act"] = coa[df["a"].to_numpy()]
    v70j = _json(DOCV / "v70_b13_local.json")
    out: dict = {"constructions": {}}

    def sums(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        g = df.groupby(key).agg(n=("e", "size"), e=("e", "sum"), d=("d", "sum"))
        return g["e"].to_numpy(), g["d"].to_numpy(), g["n"].to_numpy(np.float64)

    # 確証の実行と (a): 場所グループ、モデル固定（B-13 の B1(a) と同じ重み）。ν はそれぞれの予測の単位の寄与から
    ra = np.random.default_rng(v70.SEED)
    conf = v36.load("world_v2", "primary")
    n_c, s_c, _ = v69._group_sums(conf, v73.ARMS4)
    wc = ra.multinomial(len(n_c), np.full(len(n_c), 1.0 / len(n_c)), size=v70.B_FIXED).astype(np.float64)
    lex = v36.load(*v73.CORR)
    n, s, _ = v69._group_sums(lex, v73.ARMS4)
    w = ra.multinomial(len(n), np.full(len(n), 1.0 / len(n)), size=v70.B_FIXED).astype(np.float64)
    for key, nn, ss, ww in (("protocol_run_fixed_groups", n_c, s_c, wc), ("a_fixed_groups", n, s, w)):
        eg, dg = ss["M1c"] - ss["M3n"], ss["M0"] - ss["M1"]
        e0, d0 = float(eg.sum() / nn.sum()), float(dg.sum() / nn.sum())
        nw = ww @ nn
        var = _cov3((ww @ eg) / nw, (ww @ dg) / nw)
        out["constructions"][key] = {"R": e0 / d0, "delta_bs": e0, "yardstick": d0, "units": "location groups",
                                     **_nu_rule(e0, d0, var, _nu_same(eg, dg, nn), len(nn))}
    lx = v70j["B1bc"]["b_groups"]["refit"]
    for key, scheme, units, by in (("b_refit_groups", "groups", "location groups", "u"),
                                   ("c_refit_nested", "nested", "countries", "country_of_act")):
        r = np.load(VER / f"v70_b1_refit_reps_{scheme}.npy")
        wt = r[:, 6]
        var = _cov3((r[:, 0] - r[:, 1]) / wt, (r[:, 4] - r[:, 5]) / wt)
        xe, xd, nu_n = sums(by)
        out["constructions"][key] = {"R": lx["delta_bs"] / lx["yardstick"], "delta_bs": lx["delta_bs"],
                                     "yardstick": lx["yardstick"], "units": units,
                                     **_nu_rule(lx["delta_bs"], lx["yardstick"], var, _nu_same(xe, xd, nu_n),
                                                len(nu_n))}
    c = v70j["B1bc"]["c_nested"]["refit"]
    nb, nbc = v70j["B1d"]["nb_var"], float(v70j["B1d"]["nb_cov"])
    var_d = (c["var"][0] + nb["delta_bs"], c["var"][1] + nb["yardstick"], c["var"][2] + nbc)
    xe, xd, nu_n = sums("country_of_act")
    out["constructions"]["d_refit_nested_folds"] = {
        "R": c["R"], "delta_bs": c["delta_bs"], "yardstick": c["yardstick"], "units": "countries",
        **_nu_rule(c["delta_bs"], c["yardstick"], var_d, _nu_same(xe, xd, nu_n), len(nu_n))}
    r = np.load(VER / "v73_b4_pigeonhole_reps.npy")
    pv = _json(DOCV / "v73_b14_local.json")["B4"]["refit"]
    var_p = _cov3((r[:, 0] - r[:, 1]) / r[:, 6], (r[:, 4] - r[:, 5]) / r[:, 6])
    fu, fa = _nu_same(*sums("u")), _nu_same(*sums("a"))
    g_p = min(len(sums("u")[2]), len(sums("a")[2]))
    out["constructions"]["pigeonhole_refit"] = {
        "R": pv["R"], "delta_bs": pv["delta_bs"], "yardstick": pv["yardstick"],
        "units": "smaller of location groups and activations",
        **_nu_rule(pv["delta_bs"], pv["yardstick"], var_p, lambda th: min(fu(th), fa(th)), g_p)}
    # B-15 の B3（活性化を分けない分割、super-cluster）: 比較ごとのセルの差から ν を求め直す（B-15 の誤りの直し）
    b3 = (_json(DOCV / "v76_b15_local.json").get("B3") or {}).get("comparisons") or {}
    sc_folds = v76._sc_folds()
    if b3 and len(sc_folds) == 10:
        rows = []
        for f in sc_folds:
            x = v76._designs(f)
            yt = f["y_te"].astype(np.float64)
            p = {a: stage2_weighted(f["z_tr"], f["z_te"], x[a][0], x[a][1], f["y_tr"])[0]
                 for a in ("M1c", "M3n", "M1s", "M3d")}
            rows.append(pd.DataFrame({
                "unit": f["unit_te"].astype(str), "d": (f["p_M0_te"] - yt) ** 2 - (f["p_M1_te"] - yt) ** 2,
                "e_M1c_M3n": (p["M1c"] - yt) ** 2 - (p["M3n"] - yt) ** 2,
                "e_M1s_M3d": (p["M1s"] - yt) ** 2 - (p["M3d"] - yt) ** 2}))
        u = pd.concat(rows, ignore_index=True)
        for k in ("M1c_M3n", "M1s_M3d"):
            det = b3[k]["detail"]
            g = u.groupby("unit").agg(n=("d", "size"), e=(f"e_{k}", "sum"), d=("d", "sum"))
            out["constructions"][f"b15_supercluster_split_{k}"] = {
                "R": det["R"], "delta_bs": det["delta_bs"], "yardstick": det["yardstick"],
                "units": "super-clusters of the activation-merged split",
                **_nu_rule(det["delta_bs"], det["yardstick"], tuple(det["var"]),
                           _nu_same(g["e"].to_numpy(), g["d"].to_numpy(), g["n"].to_numpy(np.float64)), len(g))}
    b1 = (_json(OUT).get("B1") or {}) if OUT.is_file() else {}
    for k in ("M1c_M3n", "M1s_M3d"):
        row = ((b1.get("comparisons") or {}).get(k) or {}).get("all")
        if row:
            out["constructions"][f"b16_b1_{k}"] = {"R": row["band_yardstick"]["R"],
                                                   "units": "super-clusters of the primary split",
                                                   **row["few_cluster"]}
    flagged = []
    for k, v in out["constructions"].items():
        up = (v.get("t_nu_minus_1") or {}).get("upper_R")
        if up is None or up > SESOI_R:
            flagged.append(k)
    out["constructions_without_robust_bound_below_margin"] = flagged
    out["reporting_rule_triggered"] = bool(flagged)
    out["rule"] = ("B-16 B3 (addendum 1, item 5): nu(theta) = (sum gamma)^2 / sum gamma^2 from the unit contributions to "
                   "L+(theta); theta_U and nu iterated with t(nu - 1) until the change < 1e-6; no upper limit if "
                   "nu - 1 < 1 or the Fieller interval is unbounded; z and t(nu) with the same nu; G* = nu at theta 0")
    return out


# %% B4: 上限の注入の OR
def b4_odds() -> dict:
    res = _json(OUT) if OUT.is_file() else {}
    uppers: dict[str, float] = {}
    unbounded: list[str] = []
    for part, path in (("B3", "constructions"),):
        for k, v in ((res.get(part) or {}).get(path) or {}).items():
            for crit in ("z", "t_nu", "t_nu_minus_1"):
                up = (v.get(crit) or {}).get("upper_R")
                lab = f"{k}:{crit}"
                if up is None:
                    unbounded.append(lab)
                else:
                    uppers[lab] = float(up)
    for k, v in ((res.get("B1") or {}).get("comparisons") or {}).items():
        for bn in ("all", "b500"):
            fc = ((v.get(bn) or {}).get("few_cluster")) or {}
            for crit in ("z", "t_nu", "t_nu_minus_1"):
                up = (fc.get(crit) or {}).get("upper_R")
                lab = f"B1_{k}_{bn}:{crit}"
                if up is None:
                    unbounded.append(lab)
                elif bn == "all":
                    uppers[lab] = float(up)
    if not uppers:
        return {"note": "run B1 and B3 first"}
    conv = v73._or_conversions(uppers)
    return {"conversions": conv, "unbounded": unbounded, "margin_remove_flooding": _margin_remove(),
            "rule": "B-16 B4: _or_conversions of v73 (curves of B2 of B-13, four no-effect references); bands other "
                    "than all ground are not converted (the curves are for effects on all ground). For A3, the margin "
                    "(0.05 D) is also translated to an effect that removes flooding, with the curves of C2 of B-9 "
                    "(v64.c2_negative, the same interpolation and references)"}


def _margin_remove() -> dict:
    """A3: 余白（0.05 D）に届く、洪水を減らす注入の OR（B-9 の C2 の曲線と同じ換算、四つの無効果の基準）。"""
    v64 = v70._load("v64_b9_aggregate")
    r = v64.run_gains(v64.NEG, ("A",))
    if r is None:
        return {"note": "curves of B-9 C2 not found"}
    deltas, gains, _ = r
    o = np.argsort(np.abs(deltas))
    ad = np.abs(deltas)[o]
    grid = np.concatenate([[0.0], ad])
    yard = float(_json(RUNS_DIR / CORR_RUN / "metrics_primary.json")["primary"]["delta_bs_M0_M1"])
    rng = np.random.default_rng(SEED)
    out: dict = {"yardstick": yard, "references": {}}
    for ref, dnull in v64.v59.references().items():
        e = v64.v59._entry(gains["A"][o], grid, "bound", SESOI_R * yard - dnull, float(np.exp(ad[-1])), rng)
        x = e.get("delta")
        out["references"][ref] = {"odds_ratio": float(np.exp(-x)) if x is not None else None,
                                  "beyond_grid_below": float(np.exp(-ad[-1])) if x is None else None}
    ors = [v["odds_ratio"] for v in out["references"].values() if v["odds_ratio"] is not None]
    out["range"] = [min(ors), max(ors)] if ors else None
    return out


# %% B5: 帰無を課した wild cluster score bootstrap（Kline and Santos 2012 の形）
WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5), math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def _fe_fit(y: np.ndarray, zr: np.ndarray, ia: np.ndarray, na: int, off: np.ndarray,
            start: tuple[np.ndarray, np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """AOI 固定効果 + zr のロジスティック回帰（オフセットつき、v39.fit と同じ IRLS）。(alpha, beta) を返す。"""
    n, k = zr.shape
    d = sparse.csr_matrix((np.ones(n), (np.arange(n), ia)), shape=(n, na))
    if start is None:
        rate = np.bincount(ia, weights=y, minlength=na) / np.bincount(ia, minlength=na)
        alpha, beta = np.log(rate / (1 - rate)), np.zeros(k)
    else:
        alpha, beta = start[0].copy(), start[1].copy()
    for _ in range(100):
        p = 1 / (1 + np.exp(-(alpha[ia] + zr @ beta + off)))
        w = p * (1 - p)
        r = y - p
        haa = np.bincount(ia, weights=w, minlength=na)
        ga = np.bincount(ia, weights=r, minlength=na)
        if k:
            haz = d.T @ (zr * w[:, None])
            hzz = zr.T @ (zr * w[:, None])
            s = hzz - haz.T @ (haz / haa[:, None])
            db = np.linalg.solve(s, zr.T @ r - haz.T @ (ga / haa))
            da = (ga - haz @ db) / haa
        else:
            db, da = np.zeros(0), ga / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max() if k else 0.0) < 1e-10:
            break
    return alpha, beta


def _eff_scores(y: np.ndarray, z: np.ndarray, ia: np.ndarray, na: int, ig: np.ndarray, ng: int, j: int,
                p: np.ndarray) -> np.ndarray:
    """帰無の当ての下での列 j の効率的スコアの場所グループの和（AOI を除き、残りの係数で調整）。"""
    n = len(y)
    d = sparse.csr_matrix((np.ones(n), (np.arange(n), ia)), shape=(n, na))
    w = p * (1 - p)
    haa = np.bincount(ia, weights=w, minlength=na)
    haz = d.T @ (z * w[:, None])
    s = z.T @ (z * w[:, None]) - haz.T @ (haz / haa[:, None])
    zc = z - (haz / haa[:, None])[ia]
    u = np.zeros((ng, z.shape[1]))
    np.add.at(u, ig, zc * (y - p)[:, None])
    o = [i for i in range(z.shape[1]) if i != j]
    coef = np.linalg.solve(s[np.ix_(o, o)], s[o, j]) if o else np.zeros(0)
    return u[:, j] - (u[:, o] @ coef if o else 0.0)


def _restricted_wild(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, groups: np.ndarray, j: int) -> dict:
    ua, ia = np.unique(aoi, return_inverse=True)
    ug, ig = np.unique(groups, return_inverse=True)
    na, ng = len(ua), len(ug)
    zr = np.delete(z, j, axis=1)
    rng = np.random.default_rng(SEED)
    v = WEBB[rng.integers(0, 6, size=(B_WILD, ng))]
    v2 = v ** 2
    full = v39.fit(y, z, aoi, groups)
    b_hat = float(full["beta"][j])
    pvals = {}
    order = np.argsort(np.abs(GRID - b_hat))  # 推定に近い δ₀ から順に（前の当てを初期値に）
    start = None
    for i in order:
        d0 = float(GRID[i])
        alpha, beta = _fe_fit(y, zr, ia, na, d0 * z[:, j], start)
        start = (alpha, beta)
        p = 1 / (1 + np.exp(-(alpha[ia] + zr @ beta + d0 * z[:, j])))
        se = _eff_scores(y, z, ia, na, ig, ng, j, p)
        w_obs = float(se.sum() ** 2 / (se ** 2).sum())
        w_star = (v @ se) ** 2 / (v2 @ (se ** 2))
        pvals[d0] = float((1 + np.sum(w_star >= w_obs)) / (1 + B_WILD))
    acc = sorted(d for d, pv in pvals.items() if pv > 0.05)
    contiguous = bool(acc) and all(abs(b - a - 0.01) < 1e-9 for a, b in itertools.pairwise(acc))
    lo, hi = (acc[0], acc[-1]) if acc else (None, None)
    return {"odds_ratio": math.exp(b_hat), "G": ng, "n_cells": len(y), "n_aoi": na,
            "cr_normal_ci95": [math.exp(b_hat - 1.96 * full["se"][j]), math.exp(b_hat + 1.96 * full["se"][j])],
            "p_restricted_or_1": pvals.get(0.0),
            "ci95_restricted": [math.exp(lo), math.exp(hi)] if acc else None,
            "ci_at_grid_edge": bool(acc and (lo <= GRID[0] + 1e-9 or hi >= GRID[-1] - 1e-9)),
            "ci_contiguous": contiguous}


def b5_near_inference() -> dict:
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    out: dict = {}
    for key, ycol in (("planned_gfm", "y_gfm"), ("ems_same_cells", "y")):
        y = g[ycol].to_numpy(np.float64)
        rate = pd.Series(y).groupby(g["aoi"].to_numpy()).transform("mean").to_numpy()
        s = g.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
        z, names = v39.design(s, ("wetside",), True, banded=True)
        out[key] = _restricted_wild(s[ycol].to_numpy(np.float64), z, s["aoi"].to_numpy(),
                                    s["unit"].to_numpy(str), names.index("wetside:le500m"))
        print(f"  B5 {key}: {out[key]}", flush=True)
    t = pd.read_parquet(v39.TABLE)
    rate = t.groupby("aoi")["y"].transform("mean")
    t = t[(rate > 0) & (rate < 1)].reset_index(drop=True)
    z, names = v39.design(t, ("wetside",), True, banded=False)
    out["all_ground_ems_capped"] = _restricted_wild(t["y"].to_numpy(np.float64), z, t["aoi"].to_numpy(),
                                                    t["unit"].to_numpy(str), names.index("wetside:all"))
    prev = (_json(DOCV / "v76_b15_local.json").get("B10") or {})
    out["b15_unrestricted"] = {k: {kk: prev[k].get(kk) for kk in ("odds_ratio", "wild_ci95", "cr_t_ci95")}
                               for k in ("planned_gfm", "ems_same_cells") if isinstance(prev.get(k), dict)}
    out["without_EMSR052_EMSR130"] = "existing value (Table S29: 1.10, 0.83-1.41); no new computation"
    out["rule"] = ("B-16 B5 (addendum 1, item 7): restricted wild cluster score bootstrap over location groups; "
                   f"efficient score of the term adjusted for the other coefficients with AOI effects profiled out; "
                   f"W = (sum s)^2 / sum s^2; Webb six-point weights, {B_WILD} draws, seed {SEED}, same weights at "
                   "every delta0; p = (1 + #{W* >= W}) / (1 + B); 95% interval = delta0 on a grid of log odds ratios "
                   "(-1.6 to 1.6, step 0.01) with p > 0.05")
    return out


# %% B6: ラベルの層・目的・恒久水・帯（直した実行、モデル固定）
CELLS_CACHE = VER / "v77_corr_cells.parquet"


def _join_cells(aoi: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> pd.DataFrame:
    """予測のセルを AOI と経緯度（float32 にして 10⁻⁶ 度に丸める）で cells.parquet・features.parquet に結ぶ。"""
    key = pd.DataFrame({"aoi": np.asarray(aoi).astype(str),
                        "klon": np.round(np.asarray(lon, dtype=np.float32).astype(np.float64), 6),
                        "klat": np.round(np.asarray(lat, dtype=np.float32).astype(np.float64), 6)})
    key["row"] = np.arange(len(key))
    parts = []
    for a, sub in key.groupby("aoi", sort=False):
        c = pq.read_table(LABELS / a / "cells.parquet", columns=["ix", "iy", "lon", "lat", "f_riverine"]).to_pandas()
        fe = pq.read_table(LABELS / a / "features.parquet", columns=["P_d1", "P_gid1"]).to_pandas()
        c["P_d1"], c["P_gid1"] = fe["P_d1"].to_numpy(), fe["P_gid1"].to_numpy()
        c["klon"] = np.round(c["lon"].to_numpy(np.float32).astype(np.float64), 6)
        c["klat"] = np.round(c["lat"].to_numpy(np.float32).astype(np.float64), 6)
        m = sub.merge(c.drop(columns=["lon", "lat"]), on=["klon", "klat"], how="left", validate="one_to_one")
        parts.append(m)
    out = pd.concat(parts, ignore_index=True).sort_values("row").reset_index(drop=True)
    if out["ix"].isna().any():
        msg = f"{int(out['ix'].isna().sum())} cells not joined by AOI and coordinates"
        raise SystemExit(msg)
    return out


def _corr_table() -> pd.DataFrame:
    if CELLS_CACHE.is_file():
        return pd.read_parquet(CELLS_CACHE)
    d = v36.load(CORR_RUN, "primary")
    t = pd.DataFrame({k: d[k] for k in ("y", "aoi", "unit", "lon", "lat", "p_M0", "p_M1", "p_M1c", "p_M3n")})
    j = _join_cells(t["aoi"].to_numpy(), t["lon"].to_numpy(), t["lat"].to_numpy())
    for c in ("ix", "iy", "P_d1", "P_gid1", "f_riverine"):
        t[c] = j[c].to_numpy()
    t.to_parquet(CELLS_CACHE)
    return t


def _r_fixed(t: pd.DataFrame, mask: np.ndarray, rng: np.random.Generator, *, ycol: str = "y",
             arms: tuple[str, str, str, str] = ("M0", "M1", "M1c", "M3n"), overall: bool = False) -> dict:
    """モデル固定の Δ・R（部分集合の物差し）と、overall なら全体の物差しの R（場所グループの 4,000 回）。"""
    yy = t[ycol].to_numpy(np.float64)
    a0, a1, ac, an = arms
    sq = {a: (t[f"p_{a}"].to_numpy(np.float64) - yy) ** 2 for a in arms}
    u, inv = np.unique(t["unit"].to_numpy(str), return_inverse=True)
    g = len(u)
    m = np.asarray(mask, dtype=bool)
    if m.sum() == 0:
        return {"n_cells": 0}
    nb = np.bincount(inv, weights=m.astype(np.float64), minlength=g)
    eg = np.bincount(inv, weights=np.where(m, sq[ac] - sq[an], 0.0), minlength=g)
    dg = np.bincount(inv, weights=np.where(m, sq[a0] - sq[a1], 0.0), minlength=g)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)
    nw = w @ nb
    ok = nw > 0
    e0, d0 = float(eg.sum() / nb.sum()), float(dg.sum() / nb.sum())
    det = v73._levels_detail(e0, d0, (w @ eg)[ok] / nw[ok], (w @ dg)[ok] / nw[ok], rng)
    out = {"n_cells": int(m.sum()), "n_units": int((nb > 0).sum()), "flooded_share": float(yy[m].mean()),
           "subset_yardstick": det}
    if overall:
        na = np.bincount(inv, minlength=g).astype(np.float64)
        da = np.bincount(inv, weights=sq[a0] - sq[a1], minlength=g)
        d_all = float(da.sum() / na.sum())
        var = _cov3((w @ eg)[ok] / nw[ok], (w @ da)[ok] / (w @ na)[ok])
        fi = fieller_interval(e0, d_all, var, z=Z90)
        out["overall_yardstick"] = {"R": e0 / d_all, "fieller90": list(fi) if fi else None}
    return out


def _aoi_meta(aois: list[str]) -> pd.DataFrame:
    summ = _json(DATA_DIR / "global" / "gfm_summary.json")["aois"]
    sens = {"s1": "radar", "sar": "radar", "optical": "optical"}
    v38 = _json(VER / "v38_ems_sensors_raw.json")
    inv = _json(DATA_DIR / "global" / "ems" / "inventory.json")
    img: dict[str, list[str]] = {}
    for a in inv["activations"]:
        for p in a.get("del_products") or []:
            if "MONIT" in str(p.get("kind", "")) or not p.get("image_times"):
                continue
            img.setdefault(f"{a['code']}_AOI{p['aoi']}", []).extend(str(x) for x in p["image_times"])
    sc = _superclusters()
    rows = []
    for a in aois:
        meta = _json(LABELS / a / "meta.json")
        post = []
        for s in (v38.get(a) or {}).get("sources") or []:
            if "post" in str(s.get("phase", "")).lower():
                try:
                    post.append(pd.to_datetime(str(s.get("date", "")), format="%d/%m/%Y"))
                except (ValueError, TypeError):
                    continue
        first = min(post) if post else (pd.to_datetime(min(img[a])[:10]) if a in img else None)
        ev = pd.to_datetime(str(sc["aoi_event"].get(a, ""))[:10], errors="coerce")
        delay = float((first - ev).days) if first is not None and not pd.isna(ev) else float("nan")
        rows.append({"aoi": a, "sensor": sens.get(str(summ.get(a, {}).get("sensor", "none")), "unknown"),
                     "generation": meta.get("schema_generation"),
                     "hydro_layer": bool((meta.get("layers") or {}).get("hydro")),
                     "untyped": bool(sc["aoi_untyped"].get(a, False)), "delay_days": delay})
    return pd.DataFrame(rows)


def _capped_or(aois: set[str]) -> dict:
    """Eq. (2) の上限ありの EMS の推定（v37 の表、凍結した読み）を AOI の部分集合で。"""
    t = pd.read_parquet(v39.TABLE)
    t = t[t["aoi"].isin(aois)]
    rate = t.groupby("aoi")["y"].transform("mean")
    t = t[(rate > 0) & (rate < 1)].reset_index(drop=True)
    if t.empty:
        return {"n_cells": 0}
    z, names = v39.design(t, ("wetside",), True, banded=True)
    j = names.index("wetside:le500m")
    near = (t["P_d1"].to_numpy() <= NEAR_M) & (t["lex_wetside_n1"].to_numpy() > 0) & (t["y"].to_numpy() == 1)
    try:
        r = v39.fit(t["y"].to_numpy(np.float64), z, t["aoi"].to_numpy(), t["unit"].to_numpy(str))
    except np.linalg.LinAlgError as e:
        return {"error": str(e), "n_flooded_near": int(near.sum())}
    b, se = float(r["beta"][j]), float(r["se"][j])
    return {"odds_ratio": math.exp(b), "ci95": [math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)],
            "n_cells": len(t), "n_aoi": int(t["aoi"].nunique()), "n_flooded_near": int(near.sum())}


def b6_labels() -> dict:
    t = _corr_table()
    rng = np.random.default_rng(SEED)
    aois = sorted(t["aoi"].unique())
    meta = _aoi_meta(aois).set_index("aoi")
    known = meta["delay_days"].dropna()
    q1, q2 = np.quantile(known.to_numpy(), [1 / 3, 2 / 3]) if len(known) else (np.nan, np.nan)
    meta["delay"] = np.where(meta["delay_days"].isna(), "unknown",
                             np.where(meta["delay_days"] <= q1, "t1", np.where(meta["delay_days"] <= q2, "t2", "t3")))
    out: dict = {"overall": _r_fixed(t, np.ones(len(t), dtype=bool), rng),
                 "delay_tercile_cuts_days": [float(q1), float(q2)], "strata": {}}
    for col in ("sensor", "generation", "hydro_layer", "untyped", "delay"):
        res = {}
        for lev in sorted(meta[col].dropna().unique(), key=str):
            sel = set(meta.index[meta[col] == lev])
            m = t["aoi"].isin(sel).to_numpy()
            res[str(lev)] = {"n_aoi": len(sel), **_r_fixed(t, m, rng), "capped_ems_or": _capped_or(sel)}
        out["strata"][col] = res
    fr = t["f_riverine"].to_numpy(np.float64)
    out["targets"] = {"fraction": _r_fixed(t.assign(yf=fr), np.ones(len(t), dtype=bool), rng, ycol="yf")}
    for thr in (0.10, 0.25):
        out["targets"][f"threshold_{thr:.2f}"] = _r_fixed(t.assign(yt=(fr >= thr).astype(np.float64)),
                                                          np.ones(len(t), dtype=bool), rng, ycol="yt")
    gsw = pd.read_parquet(DATA_DIR / "global" / "gsw" / "cell_gsw.parquet", columns=["aoi", "ix", "iy", "occ"])
    occ = t[["aoi", "ix", "iy"]].merge(gsw, on=["aoi", "ix", "iy"], how="left", validate="many_to_one")["occ"]
    keep = ~(occ.to_numpy(np.float64) >= 50)
    out["gsw_mask"] = {"n_removed": int((~keep).sum()), **_r_fixed(t, keep, rng)}
    band = _band_index(t["P_d1"].to_numpy())
    out["bands"] = {bn: _r_fixed(t, band == i, rng, overall=True) for i, bn in enumerate(("b500", "b1k", "bfar"))}
    out["rule"] = ("B-16 B6 (addendum 1, item 8): held-out predictions of world_v2_lexfix (primary, models fixed); "
                   "cells joined to cells.parquet by AOI and float32 coordinates; R with the stratum's or target's "
                   "own yardstick BS(M0) - BS(M1); 4,000 draws of location groups (seed 20261105); capped EMS "
                   "estimate of Eq. (2) on the v37 table (frozen reading) restricted to the stratum's AOIs")
    return out


# %% B7: GFM と EMS の一致
def _agree(e: np.ndarray, g: np.ndarray) -> dict:
    a, b = int(((g == 1) & (e == 1)).sum()), int(((g == 1) & (e == 0)).sum())
    c, d = int(((g == 0) & (e == 1)).sum()), int(((g == 0) & (e == 0)).sum())
    n = a + b + c + d
    po = (a + d) / n if n else float("nan")
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2 if n else float("nan")
    return {"n": n, "both_flooded": a, "gfm_only": b, "ems_only": c, "neither": d,
            "agreement": po, "kappa": (po - pe) / (1 - pe) if n and pe < 1 else None,
            "hit_rate": a / (a + c) if a + c else None, "miss_rate": c / (a + c) if a + c else None,
            "false_alarm_ratio": b / (a + b) if a + b else None, "csi": a / (a + b + c) if a + b + c else None}


def b7_agreement() -> dict:
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    e, f = g["y_ems"].to_numpy(), g["y_gfm"].to_numpy()
    band = _band_index(g["P_d1"].to_numpy())
    sens = g["sensor"].map({"s1": "radar", "sar": "radar", "optical": "optical"}).fillna("unknown").to_numpy()
    out = {"all": _agree(e, f), "n_aoi": int(g["aoi"].nunique()),
           "bands": {bn: _agree(e[band == i], f[band == i]) for i, bn in enumerate(("b500", "b1k", "bfar"))},
           "ems_sensor": {s: _agree(e[sens == s], f[sens == s]) for s in ("radar", "optical", "unknown")}}
    out["rule"] = ("B-16 B7: cells of the v41 table assessable by GFM; EMS as the reference; hit rate P(GFM=1|EMS=1), "
                   "false-alarm ratio P(EMS=0|GFM=1), CSI a/(a+b+c), Cohen's kappa")
    return out


# %% B8: 集落の近くの OR の誤分類の偏りの解析
def _bias_range(y: np.ndarray, wet: np.ndarray, near: np.ndarray, or_adj: float) -> dict:
    a_w, n_w = float(y[near & wet].sum()), float((near & wet).sum())
    a_o, n_o = float(y[near & ~wet].sum()), float((near & ~wet).sum())
    crude = (a_w / (n_w - a_w)) / (a_o / (n_o - a_o))

    def fix(a: float, n: float, se: float, sp: float) -> float | None:
        v = (a - (1 - sp) * n) / (se + sp - 1)
        return v if 0 < v < n else None

    rows = []
    for sew in SE_GRID:
        for spw in SP_GRID:
            for seo in SE_GRID:
                for spo in SP_GRID:
                    aw, ao = fix(a_w, n_w, sew, spw), fix(a_o, n_o, seo, spo)
                    if aw is None or ao is None:
                        rows.append({"se_wet": sew, "sp_wet": spw, "se_other": seo, "sp_other": spo, "or": None})
                        continue
                    c = (aw / (n_w - aw)) / (ao / (n_o - ao))
                    rows.append({"se_wet": sew, "sp_wet": spw, "se_other": seo, "sp_other": spo,
                                 "or": or_adj * c / crude})
    ok = [r["or"] for r in rows if r["or"] is not None]
    same = [r["or"] for r in rows if r["or"] is not None and r["se_wet"] == r["se_other"]
            and r["sp_wet"] == r["sp_other"]]
    return {"counts": {"wet_flooded": a_w, "wet_cells": n_w, "other_flooded": a_o, "other_cells": n_o},
            "crude_or": crude, "adjusted_or": or_adj, "n_combinations": len(rows), "n_impossible": len(rows) - len(ok),
            "range_all": [min(ok), max(ok)] if ok else None,
            "range_nondifferential": [min(same), max(same)] if same else None, "combinations": rows}


def b8_bias() -> dict:
    out: dict = {}
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    y = g["y_gfm"].to_numpy(np.float64)
    rate = pd.Series(y).groupby(g["aoi"].to_numpy()).transform("mean").to_numpy()
    s = g.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
    b10 = (_json(DOCV / "v76_b15_local.json").get("B10") or {}).get("planned_gfm") or {}
    out["planned_gfm"] = _bias_range(s["y_gfm"].to_numpy(np.float64), s["lex_wetside_n1"].to_numpy() > 0,
                                     s["P_d1"].to_numpy() <= NEAR_M, float(b10["odds_ratio"]))
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    yu = u["y"].to_numpy(np.float64)
    ru = pd.Series(yu).groupby(u["aoi"].to_numpy()).transform("mean").to_numpy()
    su = u.loc[(ru > 0) & (ru < 1)].reset_index(drop=True)
    b2 = (_json(DOCV / "v68_b11_local.json").get("B2") or {}).get("ems_all") or {}
    out["uncapped_ems"] = _bias_range(su["y"].to_numpy(np.float64), su["lex_wetside_n1"].to_numpy() > 0,
                                      su["P_d1"].to_numpy() <= NEAR_M, float(b2["odds_ratio"]))
    out["rule"] = ("B-16 B8 (addendum 1, item 10): within 500 m, Se in {0.6, 0.8, 1.0} and Sp in {0.99, 0.995, 1.0} for "
                   "the ground of water-related names and other ground separately (81 combinations); corrected count "
                   "a = (a_obs - (1 - Sp) N) / (Se + Sp - 1), impossible if <= 0 or >= N; corrected OR = adjusted OR x "
                   "corrected crude OR / observed crude OR; rows of each estimator (AOIs with variation)")
    return out


# %% B9: 読みの点検
def _perm_table() -> dict | None:
    """B-7 の置換（地域構成を作り直す、200 回）の場所グループごとの二乗誤差の和を fold をまたいで並べる。"""
    src = RUNS_DIR / PERM_RUN / "predictions" / "primary"
    files: dict[int, list[tuple[int, Path]]] = {}
    for f in src.glob("fold*_permsrc_*.npz"):
        files.setdefault(int(f.name.split("_")[0][4:]), []).append((int(f.stem.split("_")[-1]), f))
    if len(files) < 10:
        return None
    parts: dict[str, list] = {k: [] for k in ("units", "n", "real_M1c", "real_M3n", "real_M1s", "real_M3d", "ref",
                                              "names", "cent")}
    for k in sorted(files):
        zs = [np.load(f, allow_pickle=False) for _, f in sorted(files[k])]
        z0 = zs[0]
        parts["units"].append(z0["units"].astype(str))
        parts["n"].append(z0["n_unit"].astype(np.float64))
        for a in ("M1c", "M3n", "M1s", "M3d"):
            parts[f"real_{a}"].append(z0[f"real_{a}"])
        for key, src_key in (("ref", "sse_ref"), ("names", "sse_names"), ("cent", "sse_cent")):
            parts[key].append(np.concatenate([z[src_key] for z in zs], axis=0))
    out = {k: np.concatenate(v, axis=-1) for k, v in parts.items()}
    out["n_perms"] = out["ref"].shape[0]
    return out


def b9a_placebo() -> dict:
    import world_aggregate as wa

    p = _perm_table()
    if p is None:
        return {"note": f"{PERM_RUN} not found locally"}
    units = p["units"]
    # 物差しの腕（M0・M1）の場所グループの和: 実の M1c・M3n が一致する実行から
    pick, sums = None, {}
    for run in (CORR_RUN, "world_v2"):
        d = v36.load(run, "primary")
        n, s, u = v69._group_sums(d, ("M0", "M1", "M1c", "M3n"))
        pos = pd.Series(np.arange(len(u)), index=u.astype(str)).reindex(units).to_numpy()
        if np.isnan(pos).any():
            continue
        pos = pos.astype(int)
        diff = max(float(np.max(np.abs(s[a][pos] - p[f"real_{a}"]) / np.maximum(p[f"real_{a}"], 1e-12)))
                   for a in ("M1c", "M3n"))
        sums[run] = {"max_rel_diff_real": diff, "m0": s["M0"][pos], "m1": s["M1"][pos]}
        if pick is None or diff < sums[pick]["max_rel_diff_real"]:
            pick = run
    if pick is None:
        return {"note": "no run with the same location groups"}
    m0, m1, n = sums[pick]["m0"], sums[pick]["m1"], p["n"]
    dg = m0 - m1
    rng = np.random.default_rng(SEED)
    g = len(units)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)
    db = w @ dg
    den = db / (w @ n)
    out: dict = {"yardstick_run": pick, "match": {k: v["max_rel_diff_real"] for k, v in sums.items()},
                 "n_perms": int(p["n_perms"]), "n_units": g}
    for key, real_e, perm_e in (("nested", p["real_M1c"] - p["real_M3n"], p["ref"] - p["names"]),
                                ("centred", p["real_M1s"] - p["real_M3d"], p["real_M1s"][None, :] - p["cent"])):
        r_hat = float(real_e.sum() / dg.sum())
        rb = (w @ real_e) / db
        rp_hat = perm_e.sum(axis=1) / dg.sum()
        diffb = rb[:, None] - (w @ perm_e.T) / db[:, None]
        ci_r95, ci_r90 = wa.ci(rb, 0.95), wa.ci(rb, 0.90)
        labels, lows, highs = [], [], []
        for j in range(perm_e.shape[0]):
            ci_d = wa.ci(diffb[:, j], 0.95)
            lows.append(ci_d[0])
            highs.append(ci_d[1])
            labels.append(wa.decide(r_hat, ci_r95, ci_r90, placebo_ci95=ci_d, den_ci95=wa.ci(den, 0.95))["label"])
        base = wa.decide(r_hat, ci_r95, ci_r90, den_ci95=wa.ci(den, 0.95))["label"]
        out[key] = {"R": r_hat, "ci95_R": ci_r95, "ci90_R": ci_r90, "verdict_without_placebo": base,
                    "Rp_range": [float(rp_hat.min()), float(rp_hat.max())],
                    "R_minus_Rp_ci95_lower_range": [float(min(lows)), float(max(lows))],
                    "R_minus_Rp_ci95_upper_range": [float(min(highs)), float(max(highs))],
                    "share_R_minus_Rp_lower_above_0": float(np.mean(np.array(lows) > 0)),
                    "verdicts": {lab: labels.count(lab) for lab in sorted(set(labels))},
                    "share_verdict_changes": float(np.mean(np.array(labels) != base))}
    out["rule"] = ("B-16 B9 (a): each of the 200 reassignments of B-7 (composition rebuilt) as the placebo; "
                   "Rp_j = [SSE_ref,j - SSE_names,j] / (N D) (centred: SSE(M1s) - SSE_cent,j); 95% interval of "
                   "R - Rp_j from 4,000 draws of location groups (models fixed, seed 20261105, same draws for all j); "
                   "verdict from decide of world_aggregate with that interval as the placebo interval")
    return out


def _wet_forms(name: str, cc: str, lex: Lexicon) -> list[str]:
    from topsem.world.lexicon import _SPLIT, _hit

    toks = [t for t in _SPLIT.split(str(name).lower()) if t]
    return sorted({f"{e.form}({e.pos})" for e in lex.entries_for(cc) if e.cat in WET_SIDE
                   and any(_hit(t, e) for t in toks)})


def b9b_blocks() -> dict:
    arr = _point_arrays()
    if arr is None:
        return {"note": "run B1 first"}
    y = arr["y"].astype(np.float64)
    big_n = len(y)
    sq = {a: (arr[f"p_{a}"].astype(np.float64) - y) ** 2 for a in ("M0", "M1", "M1c", "M3n", "M1s", "M3d")}
    e = {"M1c_M3n": sq["M1c"] - sq["M3n"], "M1s_M3d": sq["M1s"] - sq["M3d"]}
    d = sq["M0"] - sq["M1"]
    lex = Lexicon(FROZEN_LEX)
    cc = arr["cc"].astype(str)
    block = np.array(["+".join(lex.country_langs.get(c.upper(), ())) or "none" for c in cc])
    out: dict = {"delta_total": {k: float(v.sum() / big_n) for k, v in e.items()}, "yardstick": float(d.mean())}
    for key, lab in (("by_block", block), ("by_country", cc)):
        df = pd.DataFrame({"lab": lab, "e1": e["M1c_M3n"], "e2": e["M1s_M3d"], "d": d, "y": y})
        g = df.groupby("lab").agg(n=("y", "size"), e1=("e1", "sum"), e2=("e2", "sum"), d=("d", "sum"))
        g = g.sort_values("e1")
        out[key] = {str(k): {"n_cells": int(r.n), "contribution_M1c_M3n": float(r.e1 / big_n),
                             "contribution_M1s_M3d": float(r.e2 / big_n), "contribution_D": float(r.d / big_n)}
                    for k, r in g.iterrows()}
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet", columns=["gid", "asciiname", "cc"])
    gn = gn.drop_duplicates("gid").set_index("gid")
    meta = _aoi_meta(sorted(set(arr["aoi"].astype(str)))).set_index("aoi")
    u, inv = np.unique(arr["unit"].astype(str), return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    eg = np.bincount(inv, weights=e["M1c_M3n"])
    gam = ((eg - eg.sum() / big_n * n) / big_n) ** 2
    rows = []
    for i in np.argsort(-gam)[:10]:
        m = inv == i
        aois = pd.Series(arr["aoi"][m].astype(str))
        top_aoi = aois.mode().iloc[0]
        gids = arr["gid"][m]
        names = gn.reindex(gids)
        forms: dict[str, int] = {}
        for nm, c in zip(names["asciiname"].fillna("").to_numpy(), cc[m], strict=True):
            if c.upper() == "NO":  # 凍結した読みではノルウェーの名は読まれない
                continue
            for f in _wet_forms(nm, c, lex):
                forms[f] = forms.get(f, 0) + 1
        mt = meta.reindex(aois.unique())
        rows.append({"unit": str(u[i]), "share_of_variance": float(gam[i] / gam.sum()),
                     "contribution_M1c_M3n": float(eg[i] / big_n),
                     "contribution_D": float(np.bincount(inv, weights=d)[i] / big_n),
                     "majority_country": str(pd.Series(cc[m]).mode().iloc[0]), "n_cells": int(m.sum()),
                     "n_aoi": int(aois.nunique()), "flooded_share": float(y[m].mean()),
                     "sensor": str(meta.loc[top_aoi, "sensor"]), "generation": meta.loc[top_aoi, "generation"],
                     "hydro_layer_share_aoi": float(mt["hydro_layer"].mean()),
                     "top_wet_forms": sorted(forms.items(), key=lambda kv: -kv[1])[:3]})
    out["top_ten_location_groups"] = rows
    out["rule"] = ("B-16 B9 (b): frozen-reading refit of B1; block = the languages the lexicon assigns to the country "
                   "of the nearest settlement; contributions sum(e)/N; top ten location groups by gamma_g as in B2")
    return out


MODDATE_CACHE = VER / "v77_gn_moddate.parquet"


def _moddates(gids: set[int]) -> pd.Series:
    if MODDATE_CACHE.is_file():
        t = pd.read_parquet(MODDATE_CACHE)
        return pd.Series(t["moddate"].to_numpy(), index=t["gid"].to_numpy())
    rows = []
    with zipfile.ZipFile(DATA_DIR / "global" / "geonames" / "allCountries.zip") as zf, \
            zf.open("allCountries.txt") as fh:
        for line in fh:
            i = line.find(b"\t")
            g = int(line[:i])
            if g in gids:
                rows.append((g, line.rstrip(b"\n").split(b"\t")[-1].decode()))
    t = pd.DataFrame(rows, columns=["gid", "moddate"])
    t.to_parquet(MODDATE_CACHE)
    return pd.Series(t["moddate"].to_numpy(), index=t["gid"].to_numpy())


def _first_language() -> dict[str, str]:
    out = {}
    for line in (DATA_DIR / "global" / "geonames" / "countryInfo.txt").read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        langs = f[15].split(",") if len(f) > 15 else []
        out[f[0]] = langs[0].split("-")[0].lower() if langs and langs[0] else ""
    return out


def b9c_subsets() -> dict:
    arr = _point_arrays()
    if arr is None:
        return {"note": "run B1 first"}
    t = pd.DataFrame({"y": arr["y"].astype(np.float64), "unit": arr["unit"].astype(str),
                      **{f"p_{a}": arr[f"p_{a}"] for a in ("M0", "M1", "M1c", "M3n", "M1s", "M3d")}})
    gid = arr["gid"].astype(np.int64)
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "lat", "lon", "fcode", "cc"]).drop_duplicates("gid").set_index("gid")
    info = gn.reindex(gid)
    lat60, lon60 = info["lat"].to_numpy(np.float64) * 60, info["lon"].to_numpy(np.float64) * 60
    whole = (np.abs(lat60 - np.round(lat60)) < 1e-4) & (np.abs(lon60 - np.round(lon60)) < 1e-4)
    md = _moddates({int(g) for g in np.unique(gid) if g >= 0}).reindex(gid)
    ev = pd.Series(arr["aoi"].astype(str)).map(_superclusters()["aoi_event"]).str[:10]
    before = pd.to_datetime(md.to_numpy(), errors="coerce") < pd.to_datetime(ev.to_numpy(), errors="coerce")
    lex = Lexicon(FROZEN_LEX)
    first = _first_language()
    cc = arr["cc"].astype(str)
    reads_first = np.array([first.get(c, "") in lex.country_langs.get(c.upper(), ()) and c.upper() != "NO"
                            for c in cc])
    rng = np.random.default_rng(SEED)
    out: dict = {}
    for key, m in (("not_whole_minutes", ~whole), ("ppl_only", (info["fcode"] == "PPL").to_numpy()),
                   ("modified_before_event", np.asarray(before, dtype=bool)),
                   ("block_reads_first_language", reads_first)):
        out[key] = {"nested": _r_fixed(t, m, rng, overall=True),
                    "centred": _r_fixed(t, m, rng, arms=("M0", "M1", "M1s", "M3d"), overall=True)}
    out["share_cells"] = {"whole_minutes": float(whole.mean()), "ppl": float((info["fcode"] == "PPL").mean()),
                          "modified_before_event": float(np.mean(before)),
                          "modification_date_missing": float(md.isna().mean()),
                          "block_reads_first_language": float(reads_first.mean())}
    out["rule"] = ("B-16 B9 (c): frozen-reading refit of B1 (models fixed); (i) nearest settlement not on whole "
                   "arc-minutes in both coordinates; (ii) feature code PPL; (iii) GeoNames modification date before "
                   "the AOI's event date; (iv) the block the lexicon assigns reads the first language of the country "
                   "in countryInfo.txt (Norway excluded, as in the frozen reading); 4,000 draws of location groups")
    return out


def b9d_coders() -> dict:
    v44 = v70._load("v44_wikipedia_etymology")
    v57b = v70._load("v57b_second_coder")
    samp = _json(v44.SAMPLE)
    recs = samp["records"]
    c1 = pd.read_csv(v44.CODES)
    to_code = [json.loads(x) for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines() if x.strip()]
    flow = {"drawn": len(recs), "with_name_section": int(sum(bool(r.get("sections")) for r in recs)),
            "passed_to_coder": len(to_code), "coder1_status": c1["status"].value_counts().to_dict()}
    c2 = []
    for f in sorted(v57b.CACHE.glob("batch_*.json")):
        c2 += _json(f)["codes"]
    out: dict = {"flow": flow}
    if not c2:
        out["kappa"] = "no codes of a second coder found locally"
        return out
    c2 = pd.DataFrame(c2)
    c2["id"] = c2["id"].astype(int)
    d = c1.merge(c2, on="id", how="inner", suffixes=("_1", "_2"))
    flow["coder2_status"] = d["status_2"].value_counts().to_dict()
    both = d[(d["status_1"] == "coded") & (d["status_2"] == "coded")].copy()
    s1 = both["cats"].fillna("").map(lambda s: {c for c in str(s).split(";") if c})
    s2 = both["categories"].map(lambda v: set(v) if isinstance(v, list) else set())
    out["n_coded_by_both"] = len(both)
    out["status_agreement"] = float((d["status_1"] == d["status_2"]).mean())
    kap = {}
    for c in (*CATEGORIES, "wetside"):
        a1 = np.array([bool(x & set(WET_SIDE)) if c == "wetside" else c in x for x in s1])
        a2 = np.array([bool(x & set(WET_SIDE)) if c == "wetside" else c in x for x in s2])
        po = float((a1 == a2).mean())
        pe = float(a1.mean() * a2.mean() + (1 - a1.mean()) * (1 - a2.mean()))
        kap[c] = {"n_coder1": int(a1.sum()), "n_coder2": int(a2.sum()), "agreement": po,
                  "kappa": (po - pe) / (1 - pe) if pe < 1 else None}
    out["kappa_by_class"] = kap
    out["rule"] = ("B-16 B9 (d): coder 1 = the codes of v44 (first coding), coder 2 = v57b (a model of another "
                   "developer); Cohen's kappa by class on the excerpts coded by both; flow from the v44 sample")
    return out


# %% C2: B-10 の帯の模型の帯ごとの R（world_v2_b10stk の保存した予測）
def c2_bands() -> dict:
    d = v36.load(B10_RUN, "primary")
    if d is None:
        return {"note": f"{B10_RUN} not found locally"}
    j = _join_cells(d["aoi"], d["lon"], d["lat"])
    t = pd.DataFrame({"y": d["y"].astype(np.float64), "unit": d["unit"].astype(str),
                      **{f"p_{a}": d[f"p_{a}"] for a in ("M0", "M1", "M1c", "M3n", "M1cd", "M3nd")}})
    band = _band_index(j["P_d1"].to_numpy())
    rng = np.random.default_rng(SEED)
    arms = ("M0", "M1", "M1cd", "M3nd")
    out = {"all": _r_fixed(t, np.ones(len(t), dtype=bool), rng, arms=arms),
           "bands": {bn: _r_fixed(t, band == i, rng, arms=arms, overall=True)
                     for i, bn in enumerate(("b500", "b1k", "bfar"))},
           "same_run_M1c_M3n_all": _r_fixed(t, np.ones(len(t), dtype=bool), rng)}
    out["rule"] = ("B-16 C2 (addendum 1, item 1(a)): saved held-out predictions of world_v2_b10stk; distance from P_d1 "
                   "of features.parquet after joining by AOI and coordinates; M1cd -> M3nd by band with the band's "
                   "yardstick and the overall one; 4,000 draws of location groups (models fixed, seed 20261105)")
    return out


# %%
PARTS = {"B1": b1_designated, "B2": b2_leaveout, "B3": b3_few_cluster, "B4": b4_odds, "B5": b5_near_inference,
         "B6": b6_labels, "B7": b7_agreement, "B8": b8_bias, "B9a": b9a_placebo, "B9b": b9b_blocks,
         "B9c": b9c_subsets, "B9d": b9d_coders, "C2": c2_bands}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="例: B1,B2")
    a = ap.parse_args()
    want = [p for p in (a.only.split(",") if a.only else PARTS) if p]
    for k in want:
        t0 = time.time()
        print(f"{k} ...", flush=True)
        r = PARTS[k]()
        r["seconds"] = round(time.time() - t0)
        res = _json(OUT) if OUT.is_file() else {}  # 並べて回す部分が互いの結果を消さないよう、書く直前に読み直す
        res[k] = r
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        print(f"{k} ({r['seconds']} s) {json.dumps(r, ensure_ascii=False, default=float)[:2000]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
