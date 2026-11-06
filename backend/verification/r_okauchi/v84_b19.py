# %% v84（事後・amendment B-19 の B1〜B5）: Natural Hazards への組み替えに伴う記述的な集計（utens のジョブで回す）
"""確証判定は変えない。B-19 で計算の前に固定した量を、既にある実行の held-out 予測とセルの表から計算する。
新しくモデルを当てることはしない（Eq. (2) の形のロジスティック回帰だけは当てる）。

- B1: 活性化（activation）ごとの Δ と D への寄与（主分割と活性化を保つ分割、凍結した読みと直した読み）。寄与の 90% 区間は
  その活性化の場所グループの再標本（他は固定、4,000 回）。一つずつ外した R、上位の寄与の割合、AOI の種別の推定の有無
- B2: 局所の効果の天井 R_max を半径（250 m・500 m・1 km・2 km）と、水に関わる名の集落か全ての集落かで。AOI の集落の
  密度の三分位ごとの寄与。浸水したセルと乾いたセルの分け（天井がラベルと参照の模型に依ることを示す）
- B3: 二値化で失われる浸水（25 m の河川の浸水の被覆のうち、乾いたラベルのセルにある割合）を距離帯ごとに、全ての評価できる
  セル（上限なし）と解析したセル（上限あり）で。集落の 500 m の輪に浸水のラベル・25 m の浸水があるかの割合
- B4: 観測できる地面（ESA WorldCover 10 m の市街地 < 0.10 かつ樹木 < 0.50）で、主の比較の R（モデル固定）、距離帯の OR、
  Eq. (2) の水に関わる名 × 500 m 以内の OR
- B5: 集落の近くのラベルの欠け（距離帯の OR）を地図の作り方の層（センサー・製品の世代・水系の層・撮影の遅れ・痕跡の有無）で

出力: stdout と docs/verification/v84_b19.json（部分ごとに足す。`--only B1`）。Mac では回さない（2026-10-31 の指示）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
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

from topsem.config import DATA_DIR
from topsem.world.b13 import fieller_interval

OUT = ROOT / "docs" / "verification" / "v84_b19.json"
LABELS = DATA_DIR / "global" / "labels"
VER = DATA_DIR / "verification"
SEED = 20261108
B_DRAWS = 4000
Z90 = float(stats.norm.ppf(0.95))
RADII = (250.0, 500.0, 1000.0, 2000.0)
BANDS = (("le500m", -np.inf, 500.0), ("500m_1km", 500.0, 1000.0), ("gt1km", 1000.0, np.inf))
BUILT_MAX, TREE_MAX = 0.10, 0.50  # 観測できる地面（B4）
_MODS: dict = {}


def _load(name: str):
    if name not in _MODS:
        spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        _MODS[name] = mod
    return _MODS[name]


def _q(a: np.ndarray, lo: float = 0.05, hi: float = 0.95) -> list[float] | None:
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    return [float(np.quantile(a, lo)), float(np.quantile(a, hi))] if len(a) else None


# %% B1: 活性化ごとの寄与
RUNS_B1 = (("primary_frozen", "world_v2"), ("primary_corrected", "world_v2_lexfix"),
           ("activations_whole_frozen", "world_v2_actgrp"), ("activations_whole_corrected", "world_v2_b7actfix"))
TOP_B1 = 12


def _location_group_by_aoi() -> dict[str, str]:
    """AOI の場所グループ（主分割の実行 world_v2 の予測の unit）。活性化を保つ分割の実行の unit は活性化の群なので、
    B1 の区間はこの対応で場所グループを引き直す（B-19 追記 1）。"""
    d = _load("v36_yardstick_robustness").load("world_v2", "primary")
    s = pd.DataFrame({"aoi": d["aoi"].astype(str), "unit": d["unit"].astype(str)}).drop_duplicates()
    if s["aoi"].duplicated().any():
        msg = "an AOI lies in more than one location group"
        raise SystemExit(msg)
    return dict(zip(s["aoi"], s["unit"], strict=True))


def _cc_by_aoi() -> dict[str, str]:
    """AOI の国（主の実行の予測の cc の最頻値）。活性化を保つ分割の実行に cc が無いときに使う。"""
    d = _load("v36_yardstick_robustness").load("world_v2", "primary")
    s = pd.DataFrame({"aoi": d["aoi"].astype(str), "cc": d["cc"].astype(str)})
    return s.groupby("aoi")["cc"].agg(lambda x: x.value_counts().index[0]).to_dict()


def b1_activations() -> dict:
    v36 = _load("v36_yardstick_robustness")
    v77 = _load("v77_b16_local")
    sc = v77._superclusters()
    untyped = sc["aoi_untyped"]
    event = sc["aoi_event"]
    ccmap = _cc_by_aoi()
    lgmap = _location_group_by_aoi()
    out: dict = {}
    for key, run in RUNS_B1:
        d = v36.load(run, "primary")
        if d is None:
            out[key] = {"note": f"{run} predictions missing"}
            continue
        rng = np.random.default_rng(SEED)
        y = d["y"].astype(np.float64)
        e = (d["p_M1c"].astype(np.float64) - y) ** 2 - (d["p_M3n"].astype(np.float64) - y) ** 2
        dd = (d["p_M0"].astype(np.float64) - y) ** 2 - (d["p_M1"].astype(np.float64) - y) ** 2
        aoi = d["aoi"].astype(str)
        act = np.array([a.split("_")[0] for a in aoi])
        unit = np.array([lgmap[a] for a in aoi])  # 場所グループ（B-19 追記 1）
        if key.startswith("primary") and not np.array_equal(unit, d["unit"].astype(str)):
            msg = f"{run}: location groups differ from the run's units"
            raise SystemExit(msg)
        se, sd = float(e.sum()), float(dd.sum())
        ua, ia = np.unique(act, return_inverse=True)
        ea = np.bincount(ia, weights=e, minlength=len(ua))
        da = np.bincount(ia, weights=dd, minlength=len(ua))
        order = np.argsort(-np.abs(ea))
        absum = float(np.abs(ea).sum())
        rows = []
        for r, j in enumerate(order):
            m = ia == j
            aois = sorted(set(aoi[m]))
            row = {"activation": str(ua[j]), "rank_abs_delta": r + 1,
                   "country": pd.Series([ccmap.get(a, "??") for a in aoi[m]]).value_counts().index[0],
                   "event": str(event.get(aois[0], ""))[:10], "n_aoi": len(aois),
                   "n_aoi_untyped": int(sum(bool(untyped.get(a, False)) for a in aois)),
                   "n_cells": int(m.sum()), "n_flooded": int(y[m].sum()), "flooded_share": float(y[m].mean()),
                   "contrib_R": float(ea[j] / sd), "share_D": float(da[j] / sd),
                   "share_abs_delta": float(abs(ea[j]) / absum) if absum > 0 else None,
                   "R_without": float((se - ea[j]) / (sd - da[j])) if sd != da[j] else None}
            if r < TOP_B1:  # その活性化の場所グループだけを引き直す（他は固定。分母 D は固定）
                uu, iu = np.unique(unit[m], return_inverse=True)
                row["n_units"] = len(uu)
                if len(uu) > 1:
                    eg = np.bincount(iu, weights=e[m], minlength=len(uu))
                    dg = np.bincount(iu, weights=dd[m], minlength=len(uu))
                    w = rng.multinomial(len(uu), np.full(len(uu), 1.0 / len(uu)), size=B_DRAWS).astype(np.float64)
                    row["contrib_R_ci90"] = _q((w @ eg) / sd)
                    row["share_D_ci90"] = _q((w @ dg) / sd)
                else:
                    row["contrib_R_ci90"] = row["share_D_ci90"] = None
            rows.append(row)
        top3 = order[:3]
        keep = ~np.isin(ia, top3)
        out[key] = {"run": run, "n_cells": len(y), "n_activations": len(ua), "R": se / sd,
                    "share_abs_delta_top1": float(np.abs(ea[order[:1]]).sum() / absum),
                    "share_abs_delta_top3": float(np.abs(ea[order[:3]]).sum() / absum),
                    "share_abs_delta_top10": float(np.abs(ea[order[:10]]).sum() / absum),
                    "R_without_top3": float(e[keep].sum() / dd[keep].sum()),
                    "R_without_range": [float(min(x["R_without"] for x in rows if x["R_without"] is not None)),
                                        float(max(x["R_without"] for x in rows if x["R_without"] is not None))],
                    "activations": rows}
    out["rule"] = ("B-19 B1: per-activation contributions to R (sum of the cell differences of the reference and names "
                   "models' squared errors over the activation, divided by the pooled sum of the position and baseline "
                   "models' differences) and shares of D; models fixed; 90% intervals from 4,000 draws of the "
                   "activation's own location groups with all else fixed (seed 20261108; none for one group); R without "
                   "each activation with predictions fixed; location groups are those of the primary split for every "
                   "run (addendum 1: the runs that keep activations whole carry activation groups as units)")
    return out


# %% セルの表（v41 の表 + 直した読みの予測 + 水に関わる名の印 + セルの被覆・土地被覆・最寄りの集落）
def cell_table() -> pd.DataFrame:
    v41 = _load("v41_gfm_independent_labels")
    t = v41.build()  # 凍結したコードの予測（fold 順）。p_M0, p_M1, p_M1c, p_M3n, p_M3np
    if not v41.attach_lexfix(t):
        msg = "world_v2_lexfix predictions missing"
        raise SystemExit(msg)
    from topsem.config import RUNS_DIR

    for tag, col in (("world_v2", "wetside_frozen"), ("world_v2_lexfix", "wetside_corrected")):
        files = sorted((RUNS_DIR / tag / "predictions" / "primary").glob("fold[0-9].npz"),
                       key=lambda f: int(f.stem[4:]))
        t[col] = np.concatenate([np.load(f)["s_wetside"] for f in files]) > 0
    if not np.array_equal(t["wetside_frozen"].to_numpy(), t["lex_wetside_n1"].to_numpy() > 0):
        msg = "frozen water-related flags differ from the v41 table"
        raise SystemExit(msg)
    add = {c: np.full(len(t), np.nan) for c in ("P_gid1", "f_riverine", "f_trace", "lc_built", "lc_tree")}
    for a, ia in t.groupby("aoi").indices.items():
        c = pq.read_table(LABELS / a / "cells.parquet", columns=["ix", "iy", "f_riverine", "f_trace"]).to_pandas()
        f = pq.read_table(LABELS / a / "features.parquet", columns=["P_gid1", "lc_built", "lc_tree"]).to_pandas()
        key = pd.MultiIndex.from_arrays([c["ix"].to_numpy(np.int64), c["iy"].to_numpy(np.int64)])
        s = pd.DataFrame({"P_gid1": f["P_gid1"].to_numpy(np.float64), "f_riverine": c["f_riverine"].to_numpy(),
                          "f_trace": c["f_trace"].to_numpy(), "lc_built": f["lc_built"].to_numpy(np.float64),
                          "lc_tree": f["lc_tree"].to_numpy(np.float64)}, index=key)
        q = pd.MultiIndex.from_arrays([t["ix"].to_numpy()[ia].astype(np.int64),
                                       t["iy"].to_numpy()[ia].astype(np.int64)])
        r = s.reindex(q)
        for k, arr in add.items():
            arr[ia] = r[k].to_numpy(np.float64)
    for k, v in add.items():
        t[k] = v
    if np.isnan(t["f_riverine"].to_numpy()).any():
        msg = "flood cover missing for analysed cells"
        raise SystemExit(msg)
    return t


# %% B2: 天井を半径と密度で
def _ceiling(t: pd.DataFrame, pre: str, wet: np.ndarray) -> dict:
    y = t["y"].to_numpy(np.float64)
    n = len(y)
    bs = {a: float(np.mean((t[f"{pre}{a}"].to_numpy(np.float64) - y) ** 2)) for a in ("M0", "M1")}
    yard = bs["M0"] - bs["M1"]
    p = t[f"{pre}M1c"].to_numpy(np.float64)
    err = (p - y) ** 2
    d1 = t["P_d1"].to_numpy(np.float64)
    out: dict = {"yardstick": yard, "n_cells": n}
    for grp, g in (("water_related", wet), ("any_settlement", np.ones(n, dtype=bool))):
        res = {}
        for r in RADII:
            m = g & (d1 <= r)
            res[f"{int(r)}m"] = {"n_cells": int(m.sum()), "share_cells": float(m.mean()),
                                 "n_flooded": int(y[m].sum()), "R_max": float(err[m].sum() / n / yard),
                                 "R_max_flooded_part": float(err[m & (y == 1)].sum() / n / yard),
                                 "R_max_dry_part": float(err[m & (y == 0)].sum() / n / yard)}
        out[grp] = res
    # AOI の集落の密度（AOI の解析したセルから最寄りの集落までの距離の中央値）の三分位（AOI を単位に）
    med = t.groupby("aoi")["P_d1"].median()
    q1, q2 = np.quantile(med.to_numpy(), [1 / 3, 2 / 3])
    terc = t["aoi"].map(lambda a: "dense" if med[a] <= q1 else ("middle" if med[a] <= q2 else "sparse")).to_numpy()
    near = wet & (d1 <= 500.0)
    dens = {}
    for lev in ("dense", "middle", "sparse"):
        m = terc == lev
        dens[lev] = {"n_aoi": int((med <= q1).sum() if lev == "dense" else
                                  ((med > q1) & (med <= q2)).sum() if lev == "middle" else (med > q2).sum()),
                     "share_cells": float(m.mean()), "share_cells_near500_water": float((m & near).sum() / max(m.sum(), 1)),
                     "R_max_500m_contribution": float(err[m & near].sum() / n / yard),
                     "median_d1_m": float(np.median(d1[m]))}
    out["density_terciles"] = {"cuts_median_d1_m": [float(q1), float(q2)], "levels": dens}
    return out


def b2_ceiling(t: pd.DataFrame) -> dict:
    out = {"corrected": _ceiling(t, "p_lf_", t["wetside_corrected"].to_numpy(bool)),
           "frozen": _ceiling(t, "p_", t["wetside_frozen"].to_numpy(bool))}
    out["rule"] = ("B-19 B2: R_max(S) = sum over S of the reference model's held-out squared error / N / D (models "
                   "fixed), S = cells whose nearest settlement lies within the radius and has a water-related name (or "
                   "any name); split into flooded and dry cells; AOI terciles of the median distance from the AOI's "
                   "analysed cells to their nearest settlement")
    return out


# %% B3: 二値化で失われる浸水と集落の輪
def _band(d1: np.ndarray) -> np.ndarray:
    out = np.full(len(d1), "", dtype=object)
    for name, lo, hi in BANDS:
        out[(d1 > lo) & (d1 <= hi)] = name
    return out


def _binarisation(y: np.ndarray, f: np.ndarray, d1: np.ndarray) -> dict:
    band = _band(d1)
    res = {"all": None}
    for name in ("all", *(b[0] for b in BANDS)):
        m = np.ones(len(y), dtype=bool) if name == "all" else band == name
        area = float(f[m].sum())
        res[name] = {"n_cells": int(m.sum()), "flooded_label_share": float(y[m].mean()) if m.any() else None,
                     "riverine_cover_cells": area,
                     "share_cover_in_dry_cells": float(f[m & (y == 0)].sum() / area) if area > 0 else None,
                     "share_cells_any_cover": float((f[m] > 0).mean()) if m.any() else None,
                     "share_dry_cells_with_cover_ge_0.10": float((f[m & (y == 0)] >= 0.10).mean())
                     if (m & (y == 0)).any() else None}
    return res


def _rings(gid: np.ndarray, y: np.ndarray, f: np.ndarray, d1: np.ndarray) -> dict:
    m = (d1 <= 500.0) & np.isfinite(gid) & (gid >= 0)
    s = pd.DataFrame({"g": gid[m].astype(np.int64), "y": y[m], "f": f[m]})
    a = s.groupby("g").agg(any_label=("y", "max"), any_cover=("f", "max"), cover=("f", "sum"), cells=("y", "size"))
    return {"n_settlements": len(a), "share_any_flood_label": float((a["any_label"] > 0).mean()),
            "share_any_riverine_cover": float((a["any_cover"] > 0).mean()),
            "share_cover_ge_0.10_in_a_cell": float((a["any_cover"] >= 0.10).mean()),
            "share_cover_ge_0.5_cell_in_ring": float((a["cover"] >= 0.5).mean()),
            "median_cells_in_ring": float(a["cells"].median())}


def b3_binarisation(t: pd.DataFrame) -> dict:
    out: dict = {"capped": {"cells": _binarisation(t["y"].to_numpy(np.float64), t["f_riverine"].to_numpy(np.float64),
                                                   t["P_d1"].to_numpy(np.float64)),
                            "rings": _rings(t["P_gid1"].to_numpy(np.float64), t["y"].to_numpy(np.float64),
                                            t["f_riverine"].to_numpy(np.float64), t["P_d1"].to_numpy(np.float64))}}
    up = VER / "v68_uncapped.parquet"
    if up.is_file():
        u = pd.read_parquet(up, columns=["aoi", "y", "P_d1", "P_gid1", "f_riverine"])
        ok = np.isfinite(u["f_riverine"].to_numpy(np.float64))
        u = u.loc[ok]
        out["uncapped"] = {"n_unmatched_dropped": int((~ok).sum()),
                           "cells": _binarisation(u["y"].to_numpy(np.float64), u["f_riverine"].to_numpy(np.float64),
                                                  u["P_d1"].to_numpy(np.float64)),
                           "rings": _rings(u["P_gid1"].to_numpy(np.float64), u["y"].to_numpy(np.float64),
                                           u["f_riverine"].to_numpy(np.float64), u["P_d1"].to_numpy(np.float64))}
    else:
        out["uncapped"] = {"note": "v68_uncapped.parquet missing"}
    out["rule"] = ("B-19 B3: share of the riverine flood cover on the 25 m sub-grid (fractions of the whole cell, "
                   "summed) that lies in cells labelled dry, by distance band to the nearest settlement; rings: "
                   "settlements nearest to at least one cell within 500 m, with any flood label or any riverine cover "
                   "among those cells (capped analysed cells and all evaluable cells of the uncapped run)")
    return out


# %% B4: 観測できる地面
def _band_or(s: pd.DataFrame, *, wet: str | None) -> dict:
    """Eq. (2) の形（AOI の固定効果、距離帯、M1 の logit、wet があれば 水に関わる名 × 距離帯）。場所でクラスタ頑健。"""
    v39 = _load("v39_within_aoi_odds")
    y = s["y"].to_numpy(np.float64)
    rate = pd.Series(y).groupby(s["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, y = s.loc[keep].reset_index(drop=True), y[keep]
    if wet is not None:
        s = s.assign(lex_wetside_n1=s[wet].astype(np.float64))
    z, names = v39.design(s, ("wetside",) if wet else (), True, banded=True)
    try:
        r = v39.fit(y, z, s["aoi"].to_numpy(), s["unit"].to_numpy(str))
    except np.linalg.LinAlgError as e:
        return {"error": str(e)}
    out = {"n_cells": len(s), "n_aoi": int(s["aoi"].nunique()), "n_units": int(r["n_groups"])}
    for nm in ("band_500m_1km", "band_gt1km", *(("wetside:le500m",) if wet else ())):
        j = names.index(nm)
        b, se = float(r["beta"][j]), float(r["se"][j])
        out[nm] = {"odds_ratio": math.exp(b), "ci95": [math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)]}
    if wet:
        near = (s["lex_wetside_n1"].to_numpy() > 0) & (s["P_d1"].to_numpy() <= 500.0)
        out["n_flooded_near_water_related"] = int(y[near].sum())
    return out


def _r_mask(t: pd.DataFrame, m: np.ndarray, pre: str, rng: np.random.Generator) -> dict:
    """モデル固定の R（部分集合の物差しと全体の物差し）、場所グループの 4,000 回、Fieller の 90%。"""
    y = t["y"].to_numpy(np.float64)
    sq = {a: (t[f"{pre}{a}"].to_numpy(np.float64) - y) ** 2 for a in ("M0", "M1", "M1c", "M3n")}
    u, inv = np.unique(t["unit"].to_numpy(str), return_inverse=True)
    g = len(u)
    nb = np.bincount(inv, weights=m.astype(np.float64), minlength=g)
    na = np.bincount(inv, minlength=g).astype(np.float64)
    eg = np.bincount(inv, weights=np.where(m, sq["M1c"] - sq["M3n"], 0.0), minlength=g)
    dg = np.bincount(inv, weights=np.where(m, sq["M0"] - sq["M1"], 0.0), minlength=g)
    dall = np.bincount(inv, weights=sq["M0"] - sq["M1"], minlength=g)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_DRAWS).astype(np.float64)
    ok = (w @ nb) > 0
    eb, db, dab = (w @ eg)[ok] / (w @ nb)[ok], (w @ dg)[ok] / (w @ nb)[ok], (w @ dall)[ok] / (w @ na)[ok]
    e0, d0, da0 = eg.sum() / nb.sum(), dg.sum() / nb.sum(), dall.sum() / na.sum()

    def fi(num0: float, den0: float, nb_: np.ndarray, db_: np.ndarray) -> list[float] | None:
        c = np.cov(nb_, db_)
        r = fieller_interval(float(num0), float(den0), (float(c[0, 0]), float(c[1, 1]), float(c[0, 1])), z=Z90)
        return list(r) if r else None

    return {"n_cells": int(m.sum()), "share_cells": float(m.mean()), "flooded_share": float(y[m].mean()),
            "R_subset_yardstick": float(e0 / d0), "fieller90_subset": fi(e0, d0, eb, db), "D_subset": float(d0),
            "R_overall_yardstick": float(e0 / da0), "fieller90_overall": fi(e0, da0, eb, dab)}


def b4_observable(t: pd.DataFrame) -> dict:
    built, tree = t["lc_built"].to_numpy(np.float64), t["lc_tree"].to_numpy(np.float64)
    known = np.isfinite(built) & np.isfinite(tree)
    obs = known & (built < BUILT_MAX) & (tree < TREE_MAX)
    d1 = t["P_d1"].to_numpy(np.float64)
    band = _band(d1)
    out: dict = {"definition": f"WorldCover built-up < {BUILT_MAX} and tree cover < {TREE_MAX} (10 m, 2021)",
                 "share_missing_landcover": float((~known).mean()), "share_observable": float(obs.mean()),
                 "share_observable_by_band": {b[0]: float(obs[band == b[0]].mean()) for b in BANDS}}
    rng = np.random.default_rng(SEED)
    out["R"] = {"frozen": _r_mask(t, obs, "p_", rng), "corrected": _r_mask(t, obs, "p_lf_", rng)}
    out["band_or"] = {"all_ground": _band_or(t, wet=None), "observable": _band_or(t.loc[obs], wet=None)}
    out["water_related_or"] = {"all_ground": _band_or(t, wet="wetside_frozen"),
                               "observable": _band_or(t.loc[obs], wet="wetside_frozen")}
    out["rule"] = ("B-19 B4: observable ground = cells with ESA WorldCover (2021, 10 m) built-up share < 0.10 and tree "
                   "share < 0.50; R of the protocol's comparison with models fixed on that ground (own and overall "
                   "yardstick; 4,000 draws of location groups, seed 20261108, 90% Fieller); distance-band odds ratios and "
                   "Eq. (2) (frozen reading) with AOI fixed effects and the logit of M1, clustered by location group")
    return out


# %% B5: 集落の近くのラベルの欠けを地図の作り方の層で
def b5_strata(t: pd.DataFrame) -> dict:
    v77 = _load("v77_b16_local")
    aois = sorted(t["aoi"].unique())
    meta = v77._aoi_meta(aois).set_index("aoi")
    known = meta["delay_days"].dropna()
    q1, q2 = np.quantile(known.to_numpy(), [1 / 3, 2 / 3])
    meta["delay"] = np.where(meta["delay_days"].isna(), "unknown",
                             np.where(meta["delay_days"] <= q1, "t1", np.where(meta["delay_days"] <= q2, "t2", "t3")))
    tr = t.groupby("aoi")["f_trace"].max()
    meta["traces"] = [bool(tr.get(a, 0.0) > 0) for a in meta.index]
    out: dict = {"delay_tercile_cuts_days": [float(q1), float(q2)], "strata": {}}
    for col in ("sensor", "generation", "hydro_layer", "delay", "traces"):
        res = {}
        for lev in sorted(meta[col].dropna().unique(), key=str):
            sel = set(meta.index[meta[col] == lev])
            res[str(lev)] = {"n_aoi_stratum": len(sel), **_band_or(t.loc[t["aoi"].isin(sel)], wet=None)}
        out["strata"][col] = res
    out["rule"] = ("B-19 B5: distance-band odds ratios (ground 0.5-1 km and beyond 1 km from the nearest settlement "
                   "against ground within 500 m), with AOI fixed effects and the logit of M1, clustered by location group, "
                   "within strata of the EMS maps (sensor, schema generation, hydrography layer, delay tercile, flood "
                   "traces present in the AOI's labels)")
    return out


# %% 実行
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="カンマ区切り（B1,B2,...）")
    args = ap.parse_args()
    only = {x.strip() for x in args.only.split(",") if x.strip()}
    res = json.loads(OUT.read_text(encoding="utf-8")) if OUT.is_file() else {}
    t0 = time.time()
    if not only or "B1" in only:
        res["B1"] = b1_activations()
        print("B1", {k: (round(v["R"], 4), round(v["share_abs_delta_top1"], 3), v["activations"][0]["activation"])
                     for k, v in res["B1"].items() if isinstance(v, dict) and "R" in v}, f"{time.time() - t0:.0f}s",
              flush=True)
    if not only or only & {"B2", "B3", "B4", "B5"}:
        t = cell_table()
        print("cells", len(t), f"{time.time() - t0:.0f}s", flush=True)
        for key, fn in (("B2", b2_ceiling), ("B3", b3_binarisation), ("B4", b4_observable), ("B5", b5_strata)):
            if not only or key in only:
                res[key] = fn(t)
                print(key, "done", f"{time.time() - t0:.0f}s", flush=True)
                OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    res["meta"] = {"script": "backend/verification/r_okauchi/v84_b19.py", "seed": SEED, "draws": B_DRAWS,
                   "amendment": "B-19", "seconds": round(time.time() - t0, 1)}
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    print("B2 ceiling 500m (corrected):", res.get("B2", {}).get("corrected", {}).get("water_related", {}).get("500m"))
    print("B4 R (frozen):", {k: res.get("B4", {}).get("R", {}).get("frozen", {}).get(k)
                             for k in ("R_subset_yardstick", "fieller90_subset", "share_cells")})


if __name__ == "__main__":
    main()
