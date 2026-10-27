# %% v55（事後・amendment B-7 の 5〜11 のうち、既存の成果物から手元で計算するもの）
"""確証判定は変えない。B-7 で計算の前に固定した量を計算する（クラスタの実行の集計は v56）。

- 5: 氾濫原の地面を指す読み（凍結した辞書の河川・湿地）の、集落の近く（500 m 以内）のオッズ比（EMS・GFM、地形を
  与えて）と集落の地点の HAND の差（B-2 の正の対照と同じ表・同じ推定量）
- 6: (a) y = 0 のセルの浸水被覆の分布、(c) 部分的に浸水したセルを除いたオッズ比と、河川の浸水被覆を割合の結果に
  したオッズ比（準二項。点推定は二項のスコア方程式、区間は場所でクラスタ頑健）
- 7: 事象の大きさ（浸水セルの割合 / HAND ≤ 2 m のセルの割合、解析したセルで）の三分位ごとのオッズ比
- 8: CCI の市街地・樹木の割合で調整したオッズ比（EMS・GFM）と集落の地点の HAND の差（AOI の固定効果の回帰で、
  共変量なし・ありを同じ標本で並べる）
- 9: GFM のラベルでの AOI 内の AUC と、観測の浸水数で切った CSI（M1・M1c・M3n。直した実行を主、凍結した実行も）
- 10: オッズ比（水系側・区分別、距離帯別）の場所の pairs cluster bootstrap（2000 回）
- 12〜14（追記 1）: 集落の近くのモデルの係数の同時 Wald 検定、二重字（ae・oe・ue）の読み替えで増える照合、語源の照合の
  言語別の再現率。15（追記 3）: held-out の予測での辞書の特徴 13 列の係数の同時 Wald 検定
- 11: (b) 置換の null の外れた反復とその物差し、(d) 除いた AOI の 25 m の浸水面積（基準別）、(e) 解析した集落の ASCII の
  名の二重字、(f) 語源の照合の言語別の精度（v44 の値の再掲）、(g) 同じ activation で別の場所の AOI の間の残差の相関、
  (h) 欧州の R と大陸別の浸水率、(i) H3 の傾きの bootstrap、(j) 固定の閾値での M1 の CSI
  （(a)(c)(k) は既存の JSON・表から補足の表を作るだけなので、ここでは計算しない）

オッズ比は v39 と同じ形（AOI の固定効果、凍結した実行の M1 の held-out 予測の logit、場所でクラスタ頑健な標準誤差）、
B-2〜B-6 のオッズ比と同じ表（v37・v41 の表、凍結した辞書の照合）。大陸はセルの国（最寄りの集落の国）の GeoNames の
大陸コードで、論文の大陸別の表と同じ。
出力: stdout と docs/verification/v55_b7_local.json（pairs bootstrap は data/verification/v55_pairs_<key>.json に残して
使い回す）
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import os
import re
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

import world_aggregate as wa

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b7 import FLOODPLAIN_CATEGORIES, PARTIAL_COVER
from topsem.world.dataset import continent_map, evaluable, include
from topsem.world.lexicon import CATEGORIES, WET_SIDE, Lexicon

OUT = ROOT / "docs" / "verification" / "v55_b7_local.json"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
LABELS = DATA_DIR / "global" / "labels"
SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
DIGRAPH = re.compile(r"([aou])e")
SEED = 20261028
B_PAIRS = 2000
B_SET = 2000
TERM = "le500m"
LOW_HAND_M = 2.0
CELL_KM2 = 0.0625
AOI_BANDS_KM = (0.0, 50.0, 100.0, 250.0, 500.0, np.inf)
N_WORKERS = 8


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v37 = _load("v37_positive_control")
v39 = _load("v39_within_aoi_odds")
v41 = _load("v41_gfm_independent_labels")
v50 = _load("v50_gfm_diagnostics")
v44 = _load("v44_wikipedia_etymology")
ev = v41.ev


# %% 表
def table() -> pd.DataFrame:
    """v50 の表（v41 の表 + HAND・CCI + 直した実行の予測 p_lf_*）に、v37 の列・浸水被覆・氾濫原の読みを足す。"""
    t = v50.base_table()
    s = pd.read_parquet(V37T, columns=["aoi", "y", "hand_min", "lon", "lat", "P_gid1", "cc"])
    if not (np.array_equal(t["aoi"].to_numpy(), s["aoi"].to_numpy())
            and np.array_equal(t["y"].to_numpy(), s["y"].to_numpy())):
        msg = "v37 and v41 tables are not aligned"
        raise SystemExit(msg)
    for c in ("hand_min", "lon", "lat", "P_gid1", "cc"):
        t[c] = s[c].to_numpy()
    cov = []
    for a, g in t.groupby("aoi", sort=False):
        c = pq.read_table(LABELS / a / "cells.parquet",
                          columns=["ix", "iy", "f_riverine", "f_other", "f_flood_any"]).to_pandas()
        cov.append(g[["ix", "iy"]].reset_index().merge(c, on=["ix", "iy"], how="left",
                                                        validate="many_to_one").set_index("index"))
    cov = pd.concat(cov).loc[t.index]
    for c in ("f_riverine", "f_other", "f_flood_any"):
        t[c] = cov[c].to_numpy(dtype=np.float64)
    if t[["f_riverine", "f_flood_any"]].isna().any().any():
        msg = "flood cover missing for analysed cells"
        raise SystemExit(msg)
    t["lex_floodplain_n1"] = np.logical_or.reduce(
        [t[f"lex_{c}_n1"].to_numpy() > 0 for c in FLOODPLAIN_CATEGORIES]).astype(np.float64)
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    t["continent"] = t["cc"].map(lambda c: cont.get(str(c), "??"))
    t["activation"] = t["aoi"].str.split("_").str[0]
    return t


def _infos() -> tuple[pd.DataFrame, dict, list]:
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    dirs, _, infos = ev.select_aois(gn, ev.Sensitivity())
    return gn, infos, dirs


# %% オッズ比（v39 の形）
def odds(t: pd.DataFrame, ycol: str = "y", cls: str = "wetside", *, extra: tuple[str, ...] = (),
         mask: np.ndarray | None = None) -> dict:
    """集落から 500 m 以内の区分の群のオッズ比（距離帯・M1 の logit・AOI の固定効果、場所でクラスタ頑健）。"""
    s = t if mask is None else t.loc[mask]
    y = s[ycol].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(s["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, y = s.loc[keep].reset_index(drop=True), y[keep]
    z, names = v39.design(s, (cls,), True, banded=True)
    if extra:
        z = np.column_stack([z, s[list(extra)].to_numpy(dtype=np.float64)])
        names = [*names, *extra]
    res = v39.fit(y, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
    j = names.index(f"{cls}:{TERM}")
    b, se = float(res["beta"][j]), float(res["se"][j])
    g = (s["lex_wetside_n1"] if cls == "wetside" else s[f"lex_{cls}_n1"]).to_numpy() > 0
    near = g & (s["P_d1"].to_numpy() <= 500)
    out = {"odds_ratio": float(np.exp(b)), "ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
           "n_cells": len(s), "n_aoi": int(s["aoi"].nunique()), "n_units": int(res["n_groups"]),
           "n_group_cells": int(near.sum()), "n_group_flooded": float(y[near].sum()),
           "n_group_units": int(s.loc[near, "unit"].nunique())}
    for c in extra:
        out[f"coef_{c}"] = float(res["beta"][names.index(c)])
    return out


def gfm_table(t: pd.DataFrame) -> pd.DataFrame:
    """GFM で評価できるセル（B-3 で採用した AOI）と GFM のラベル y_gfm（v41・v50 と同じ）。"""
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    return v41.attach_gfm(t, summary)


# %% 5 と 8: 集落の地点の尺度（v37 の正の対照の表）
def settlement_points() -> pd.DataFrame:
    """v37 の表と派生量（HAND ≤ 2 m など）で集落の点の表を作る（v37.settlement_table と同じ。CCI の列も平均する）。"""
    t = v37.derived(pd.read_parquet(V37T))
    fl = v37.flags(t)
    fl["floodplain"] = np.logical_or.reduce([fl[c] for c in FLOODPLAIN_CATEGORIES])
    m = (t["P_d1"].to_numpy() <= v37.POINT_M) & (t["P_gid1"].to_numpy() >= 0)
    s = t.loc[m, ["aoi", "unit", "P_gid1", *v37.PC_VARS, "hand_min", "cci_built", "cci_tree"]].copy()
    for g, v in fl.items():
        s[g] = v[m].astype(np.float64)
    return s.groupby(["aoi", "unit", "P_gid1"], as_index=False).mean()


def fe_contrast(s: pd.DataFrame, var: str, gname: str, covs: tuple[str, ...]) -> dict | None:
    """集落の点の尺度の群の差を AOI の固定効果の回帰で（v37 と同じ AOI の採否、場所の cluster bootstrap）。"""
    keep = [x for _, x in s.groupby("aoi", sort=True)
            if (x[gname] > 0.5).sum() >= v37.MIN_SETTLEMENTS and (x[gname] <= 0.5).sum() >= v37.MIN_SETTLEMENTS]
    if len(keep) < 5:
        return None
    d = pd.concat(keep, ignore_index=True).dropna(subset=[var, *covs])
    d["_g"] = (d[gname] > 0.5).astype(np.float64)
    cols = [var, "_g", *covs]
    dm = d[cols] - d.groupby("aoi")[cols].transform("mean")  # 重複して引いた場所でも AOI の平均は同じ
    x = dm[["_g", *covs]].to_numpy(dtype=np.float64)
    y = dm[var].to_numpy(dtype=np.float64)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    uniq, inv = np.unique(d["unit"].to_numpy(), return_inverse=True)
    rows = [np.flatnonzero(inv == u) for u in range(len(uniq))]
    rng = np.random.default_rng(SEED)
    reps = np.full(B_SET, np.nan)
    for r in range(B_SET):
        idx = np.concatenate([rows[u] for u in rng.integers(0, len(uniq), size=len(uniq))])
        reps[r] = np.linalg.lstsq(x[idx], y[idx], rcond=None)[0][0]
    ok = np.isfinite(reps)
    return {"estimate": float(beta[0]), "ci95": [float(np.quantile(reps[ok], 0.025)),
                                                 float(np.quantile(reps[ok], 0.975))],
            "covariates": list(covs), "coef_covariates": [float(b) for b in beta[1:]],
            "n_aoi": int(d["aoi"].nunique()), "n_settlements": len(d), "n_units": len(uniq)}


def settlement_results(s: pd.DataFrame) -> dict:
    out: dict = {"n_settlement_points": len(s)}
    out["floodplain_v37_estimator"] = {v: v37.settlement_contrast(s, v, "floodplain") for v in v37.PC_VARS}
    out["wetside_v37_estimator"] = {"hand_mean": v37.settlement_contrast(s, "hand_mean", "wetside")}
    for g in ("wetside", "floodplain"):
        out[f"{g}_fe"] = {"unadjusted": fe_contrast(s, "hand_mean", g, ()),
                          "adjusted_cci": fe_contrast(s, "hand_mean", g, ("cci_built", "cci_tree"))}
    return out


# %% 6: 部分的に浸水したセル
def partial(t: pd.DataFrame, infos: dict) -> dict:
    y0 = t["y"].to_numpy() == 0
    f = t["f_flood_any"].to_numpy()
    near = t["P_d1"].to_numpy() <= 500

    def dist(sel: np.ndarray) -> dict:
        v = f[sel]
        return {"n_cells": int(sel.sum()), "share_zero": float(np.mean(v == 0)),
                "share_0_0.05": float(np.mean((v > 0) & (v <= PARTIAL_COVER))),
                "share_0.05_0.25": float(np.mean((v > PARTIAL_COVER) & (v <= 0.25))),
                "share_0.25_0.5": float(np.mean((v > 0.25) & (v < 0.5))),
                "quantiles": {str(q): float(np.quantile(v, q)) for q in (0.5, 0.9, 0.95, 0.99)}}

    pm = y0 & (f > PARTIAL_COVER)
    untyped = t["aoi"].map({a: bool(infos[a].unspecified_as_riverine) for a in t["aoi"].unique()}).to_numpy()
    riv = t["f_riverine"].to_numpy() + np.where(untyped, t["f_other"].to_numpy(), 0.0)  # ラベルと同じ河川の被覆
    tf = t.assign(y_frac=np.clip(riv, 0.0, 1.0))
    return {"distribution_y0_all": dist(y0), "distribution_y0_le500m": dist(y0 & near),
            "n_partial": int(pm.sum()), "share_partial_all_cells": float(pm.mean()),
            "share_partial_le500m": float(pm[near].mean()),
            "odds_all": odds(t), "odds_without_partial": odds(t, mask=~pm),
            "odds_fractional_riverine_cover": odds(tf, "y_frac")}


# %% 7: 事象の大きさ
def severity(t: pd.DataFrame) -> dict:
    a = t.groupby("aoi").agg(flooded=("y", "mean"),
                             low=("hand_min", lambda v: float(np.mean(v.to_numpy() <= LOW_HAND_M))))
    a["magnitude"] = a["flooded"] / a["low"].clip(lower=1e-6)  # HAND ≤ 2 m のセルが無い AOI は最上位の層へ
    q = a["magnitude"].quantile([1 / 3, 2 / 3]).to_numpy()
    lab = pd.cut(a["magnitude"], [-np.inf, q[0], q[1], np.inf], labels=["low", "middle", "high"])
    out: dict = {"tercile_cuts": q.tolist(), "n_aoi_no_low_hand": int((a["low"] == 0).sum())}
    for k in ("low", "middle", "high"):
        aois = lab.index[lab == k]
        sel = t["aoi"].isin(aois).to_numpy()
        out[k] = odds(t, mask=sel) | {"n_aoi_tercile": len(aois),
                                       "median_magnitude": float(a.loc[aois, "magnitude"].median())}
    return out


# %% 9: GFM での AOI 内の判別
def _auc(y: np.ndarray, p: np.ndarray) -> float:
    r = pd.Series(p).rank().to_numpy()
    n1 = y.sum()
    n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def gfm_discrimination(g: pd.DataFrame, prefix: str) -> dict:
    arms = ("M1", "M1c", "M3n")
    auc: dict[str, list[float]] = {m: [] for m in arms}
    csi: dict[str, list[float]] = {m: [] for m in arms}
    for _, x in g.groupby("aoi"):
        y = x["y_gfm"].to_numpy()
        k = int(y.sum())
        if k in (0, len(y)):
            continue
        for m in arms:
            p = x[f"{prefix}{m}"].to_numpy(dtype=np.float64)
            auc[m].append(_auc(y, p))
            top = np.argsort(-p, kind="stable")[:k]
            tp = int(y[top].sum())
            csi[m].append(tp / (2 * k - tp))
    out: dict = {m: {"auc_median": float(np.median(auc[m])), "csi_median": float(np.median(csi[m]))} for m in arms}
    d_auc = np.array(auc["M3n"]) - np.array(auc["M1c"])
    d_csi = np.array(csi["M3n"]) - np.array(csi["M1c"])
    out["M3n_minus_M1c"] = {"auc_median": float(np.median(d_auc)), "auc_iqr": np.quantile(d_auc, [0.25, 0.75]).tolist(),
                            "csi_median": float(np.median(d_csi)), "csi_iqr": np.quantile(d_csi, [0.25, 0.75]).tolist(),
                            "share_auc_positive": float(np.mean(d_auc > 0)), "n_aoi": len(d_auc)}
    return out


# %% 10: pairs cluster bootstrap（場所を復元抽出し、重複した場所は別の場所・別の AOI として当て直す。v50 と同じ）
_W: dict = {}


def _pairs_init(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, rows: list[np.ndarray], na: int) -> None:
    _W.update(y=y, z=z, aoi=aoi, rows=rows, na=na)


def _pairs_rep(sel: np.ndarray) -> np.ndarray | None:
    rows = _W["rows"]
    idx = np.concatenate([rows[u] for u in sel])
    pos = np.concatenate([np.full(len(rows[u]), k, dtype=np.int64) for k, u in enumerate(sel)])
    try:
        return v39.fit(_W["y"][idx], _W["z"][idx], _W["aoi"][idx] + pos * _W["na"], pos)["beta"]
    except np.linalg.LinAlgError:
        return None


def pairs_bootstrap(t: pd.DataFrame, key: str, classes: tuple[str, ...], reps: int = B_PAIRS) -> dict:
    cache = DATA_DIR / "verification" / f"v55_pairs_{key}.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))
    y = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = t.loc[keep].reset_index(drop=True)
    y = y[keep]
    z, names = v39.design(s, classes, True, banded=True)
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    g, na = len(unit_lab), int(aoi.max()) + 1
    rows = [np.flatnonzero(unit == u) for u in range(g)]
    point = v39.fit(y, z, aoi, unit)
    rng = np.random.default_rng(SEED)
    sels = [rng.integers(0, g, size=g) for _ in range(reps)]
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"  # 子のプロセスの BLAS は 1 スレッド（子は numpy を読む前にこの環境を引き継ぐ）
    t0 = time.time()
    with ProcessPoolExecutor(N_WORKERS, initializer=_pairs_init, initargs=(y, z, aoi, rows, na)) as ex:
        res = list(ex.map(_pairs_rep, sels, chunksize=8))
    beta = np.array([r if r is not None else np.full(len(names), np.nan) for r in res])
    ok = np.isfinite(beta).all(axis=1)
    out: dict = {}
    for j, nm in enumerate(names):
        if ":" not in nm:
            continue
        b0, se = float(point["beta"][j]), float(point["se"][j])
        out[nm] = {"odds_ratio": float(np.exp(b0)),
                   "sandwich_ci95": [float(np.exp(b0 - 1.96 * se)), float(np.exp(b0 + 1.96 * se))],
                   "pairs_ci95": [float(np.exp(np.quantile(beta[ok, j], 0.025))),
                                  float(np.exp(np.quantile(beta[ok, j], 0.975)))],
                   "pairs_ci90": [float(np.exp(np.quantile(beta[ok, j], 0.05))),
                                  float(np.exp(np.quantile(beta[ok, j], 0.95)))],
                   "share_beta_below_minus5": float((beta[ok, j] < -5).mean())}
    out["_meta"] = {"reps": reps, "n_ok": int(ok.sum()), "n_units": g, "n_cells": len(y),
                    "sec": round(time.time() - t0), "workers": N_WORKERS}
    cache.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"  pairs bootstrap {key}: {out['_meta']['sec']} s", flush=True)
    return out


# %% 11
def null_outliers() -> dict:
    """(b) 合成ラベル δ = 0 の回の外れ（中央値から頑健な SD の 3 倍を超える回）と、その回の物差し。"""
    reps = json.loads((RUNS_DIR / "power_v3_null" / "power_summary.json").read_text(encoding="utf-8"))["reps"]
    r = np.array([x["R"] for x in reps])
    se = np.array([(x["ci90"][1] - x["ci90"][0]) / (2 * 1.645) for x in reps])
    yard = np.array([x["delta_bs_M0_M1"] for x in reps])
    med = float(np.median(r))
    rsd = float(1.4826 * np.median(np.abs(r - med)))
    out_idx = np.flatnonzero(np.abs(r - med) > 3 * rsd)
    order = out_idx[np.argsort(-np.abs(r[out_idx] - med))]
    rank_yard = pd.Series(yard).rank(pct=True).to_numpy()
    return {"n_reps": len(reps), "robust_sd_R": rsd, "median_R": med,
            "median_yardstick": float(np.median(yard)), "median_se": float(np.median(se)),
            "n_outlying_3rsd": len(out_idx),
            "outlying": [{"rep": reps[i]["rep"], "R": float(r[i]), "se": float(se[i]), "yardstick": float(yard[i]),
                          "yardstick_percentile": float(rank_yard[i]), "decision": reps[i]["decision"]}
                         for i in order],
            "sd_R_all": float(r.std(ddof=1)),
            "sd_R_without_outlying": float(np.delete(r, out_idx).std(ddof=1)),
            "spearman_absdev_vs_yardstick": float(pd.Series(np.abs(r - med)).corr(pd.Series(yard),
                                                                                    method="spearman")),
            "spearman_se_vs_yardstick": float(pd.Series(se).corr(pd.Series(yard), method="spearman"))}


def excluded_area(gn: pd.DataFrame, infos: dict, dirs: list) -> dict:
    """(d) 採否の基準ごとの、評価できるセルの 25 m の浸水面積（被覆 × 250 m のセル）の分布。採用した AOI も並べる。"""
    keep = {d.name for d in dirs}
    rows = []
    for tag, info in sorted(infos.items()):
        why = "included" if tag in keep else include(info, n_names_inside=ev._p_names_inside(gn, info.bounds))[1]
        c = pq.read_table(LABELS / tag / "cells.parquet",
                          columns=["f_riverine", "f_other", "f_flood_any", "f_aoi", "f_observed",
                                   "f_perm_water"]).to_pandas()
        e = evaluable(c)
        riv = c["f_riverine"].to_numpy() + (c["f_other"].to_numpy() if info.unspecified_as_riverine else 0.0)
        rows.append({"tag": tag, "reason": why, "riverine_km2": float(riv[e].sum()) * CELL_KM2,
                     "any_km2": float(c.loc[e, "f_flood_any"].sum()) * CELL_KM2, "n_wet": info.n_wet})
    d = pd.DataFrame(rows)
    out: dict = {}
    for why, x in d.groupby("reason"):
        v = x["riverine_km2"].to_numpy()
        out[why] = {"n_aoi": len(x), "median_riverine_km2": float(np.median(v)),
                    "q10_q90_riverine_km2": np.quantile(v, [0.1, 0.9]).tolist(),
                    "share_riverine_ge_1km2": float(np.mean(v >= 1)), "share_riverine_ge_5km2": float(np.mean(v >= 5)),
                    "median_any_km2": float(x["any_km2"].median()), "max_riverine_km2": float(v.max())}
    return out


def digraphs(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    """(e) 解析したセルの最寄りの集落（集落の単位）の ASCII の名に ae・oe・ue を含む割合（全体と国別）。"""
    gids = np.unique(t["P_gid1"].to_numpy())
    g = gn.drop_duplicates("gid").set_index("gid").reindex(gids)
    name = g["asciiname"].fillna("").str.lower()
    has = name.str.contains("ae|oe|ue", regex=True)
    by = pd.DataFrame({"cc": g["cc"], "h": has}).groupby("cc")["h"].agg(["mean", "size"])
    by = by[by["size"] >= 50].sort_values("mean", ascending=False).head(12)
    cell_share = float(t["P_gid1"].map(dict(zip(gids, has.to_numpy(), strict=True))).mean())
    return {"n_settlements": len(gids), "share_settlements": float(has.mean()), "share_cells": cell_share,
            "by_country_settlements": {k: {"share": float(v["mean"]), "n_settlements": int(v["size"])}
                                       for k, v in by.iterrows()}}


def block_precision() -> dict:
    """(f) 語源の照合（v44、符号は言語モデル）の、辞書の水系側の読みの精度の言語別（件数つき）の再掲。"""
    v44 = json.loads((ROOT / "docs" / "verification" / "v44_wikipedia_etymology.json").read_text(encoding="utf-8"))
    return {"by_language": v44["by_language"], "overall_lex_wetside": v44["lex_wetside"],
            "source": "docs/verification/v44_wikipedia_etymology.json"}


def residual_between(t: pd.DataFrame) -> dict:
    """(g) 同じ activation で別の場所の AOI の対の、M1（凍結した実行）の held-out 残差の AOI 平均の相関（中心間の距離帯別）。

    参照として、別の activation で別の場所の AOI の対の同じ距離帯の相関も並べる。
    """
    a = t.assign(res=t["y"] - t["p_M1"]).groupby("aoi").agg(
        res=("res", "mean"), lon=("lon", "mean"), lat=("lat", "mean"), unit=("unit", "first"),
        act=("activation", "first"))
    lon, lat = np.radians(a["lon"].to_numpy()), np.radians(a["lat"].to_numpy())
    i, j = np.triu_indices(len(a), k=1)
    h = (np.sin((lat[j] - lat[i]) / 2) ** 2
         + np.cos(lat[i]) * np.cos(lat[j]) * np.sin((lon[j] - lon[i]) / 2) ** 2)
    dist = 2 * 6371.0088 * np.arcsin(np.sqrt(h))
    unit, act, res = a["unit"].to_numpy(), a["act"].to_numpy(), a["res"].to_numpy()
    diff_place = unit[i] != unit[j]
    same_act = act[i] == act[j]
    mu = res.mean()

    def corr(sel: np.ndarray) -> dict:
        if sel.sum() < 5:
            return {"n_pairs": int(sel.sum()), "corr": None}
        x = np.concatenate([res[i[sel]], res[j[sel]]]) - mu  # 対の両方の順序（対称化した対の相関）
        z = np.concatenate([res[j[sel]], res[i[sel]]]) - mu
        return {"n_pairs": int(sel.sum()), "corr": float((x * z).sum() / np.sqrt((x * x).sum() * (z * z).sum()))}

    out: dict = {"n_aoi": len(a), "same_activation": {}, "different_activation": {}}
    for lo, hi in itertools.pairwise(AOI_BANDS_KM):
        b = (dist > lo) & (dist <= hi) & diff_place
        key = f"{lo:g}-{hi:g}km"
        out["same_activation"][key] = corr(b & same_act)
        out["different_activation"][key] = corr(b & ~same_act)
    out["same_activation"]["all"] = corr(diff_place & same_act)
    out["different_activation"]["all"] = corr(diff_place & ~same_act)
    return out


def _r_subset(s: pd.DataFrame, prefix: str) -> dict:
    preds = {a: s[f"{prefix}{a}"].to_numpy(dtype=np.float64) for a in ("M0", "M1", "M1c", "M3n", "M3np")}
    sums, n, arms = wa.unit_sums(s["y"].to_numpy(dtype=np.float64), preds, s["unit"].to_numpy())
    idx = {a: j for j, a in enumerate(arms)}
    hat = sums.sum(axis=0) / n.sum()
    delta = wa.make_delta(hat, wa.bootstrap_bs(sums, n), idx)
    d01, b01 = delta("M0", "M1")
    v, p = wa.verdict(delta, "M1c", "M3n", "M3np", d01=d01, b01=b01)
    return v | (p or {}) | {"n_cells": len(s), "n_units": len(n), "n_aoi": int(s["aoi"].nunique()),
                            "flooded_share": float(s["y"].mean()), "brier_M1": float(hat[idx["M1"]])}


def europe(t: pd.DataFrame) -> dict:
    """(h) 欧州（セルの国の大陸コード EU）に限った R（欧州の物差し・欧州の場所の bootstrap、直した実行）と大陸別の浸水率。"""
    out: dict = {"flooded_share_by_continent": {k: {"share": float(v.mean()), "n_cells": len(v)}
                                                for k, v in t.groupby("continent")["y"]}}
    eu = t[t["continent"] == "EU"]
    out["europe_lexfix"] = _r_subset(eu, "p_lf_")
    out["europe_frozen"] = _r_subset(eu, "p_")
    out["non_europe_lexfix"] = _r_subset(t[t["continent"] != "EU"], "p_lf_")
    return out


def h3_slope() -> dict:
    """(i) 梯子 T0〜T4 の R（段の増分 / 物差し）の段の番号への最小二乗の傾き（凍結した実行、場所の bootstrap 4000 回）。"""
    data = wa.load_split(RUNS_DIR / "world_v2", "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    s, n, arms = wa.unit_sums(y, preds, data["unit"])
    idx = {a: j for j, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    delta = wa.make_delta(hat, wa.bootstrap_bs(s, n), idx)
    d01, b01 = delta("M0", "M1")
    est = np.array([delta(a, b)[0] for a, b in wa.LADDER_S]) / d01
    reps = np.column_stack([delta(a, b)[1] for a, b in wa.LADDER_S]) / b01[:, None]
    k = np.arange(len(wa.LADDER_S), dtype=np.float64) - (len(wa.LADDER_S) - 1) / 2
    slope = float((est * k).sum() / (k * k).sum())
    sb = (reps * k).sum(axis=1) / (k * k).sum()
    return {"rungs": [f"{a}->{b}" for a, b in wa.LADDER_S], "R_by_rung": est.tolist(), "slope_R_per_rung": slope,
            "ci95": wa.ci(sb, 0.95), "ci90": wa.ci(sb, 0.90),
            "p_one_sided_slope_ge_0": float(np.mean(sb >= 0)), "B": len(sb),
            "note": "H3 predicted a larger name increment at lower rungs, i.e. a negative slope"}


def fixed_threshold_csi(t: pd.DataFrame) -> dict:
    """(j) M1（凍結した実行）の CSI: 確率の閾値 0.5 と、その fold の訓練側の浸水率（全セル・AOI の中央値）。"""
    y = t["y"].to_numpy()
    p = t["p_M1"].to_numpy(dtype=np.float64)
    fold = t["fold"].to_numpy()
    train = {int(k): float(y[fold != k].mean()) for k in np.unique(fold)}
    out: dict = {"train_rate_by_fold": train}
    for key, thr in (("p_ge_0.5", np.full(len(p), 0.5)), ("p_ge_train_rate", np.array([train[k] for k in fold]))):
        f = p >= thr
        tp, fp, fn = int((f & (y == 1)).sum()), int((f & (y == 0)).sum()), int((~f & (y == 1)).sum())
        per = []
        for _, g in t.assign(_f=f).groupby("aoi"):
            yy, ff = g["y"].to_numpy(), g["_f"].to_numpy()
            a, b, c = int((ff & (yy == 1)).sum()), int((ff & (yy == 0)).sum()), int((~ff & (yy == 1)).sum())
            if a + b + c:
                per.append(a / (a + b + c))
        out[key] = {"csi_pooled": tp / max(tp + fp + fn, 1), "csi_median_aoi": float(np.median(per)),
                    "share_flagged": float(f.mean()), "hit_rate": tp / max(tp + fn, 1),
                    "false_alarm_ratio": fp / max(tp + fp, 1)}
    return out


# %% 12〜14（追記 1）
def wald_tests(t: pd.DataFrame) -> dict:
    """12: EMS の集落の近くのモデルで、区分 × 距離帯の係数が全て 0 の同時 Wald 検定（場所でクラスタ頑健）。"""
    y = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, y = t.loc[keep].reset_index(drop=True), y[keep]
    out: dict = {}
    for key, classes in (("wetside", ("wetside",)), ("classes", tuple(CATEGORIES))):
        z, names = v39.design(s, classes, True, banded=True)
        res = v39.fit(y, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
        idx = [j for j, nm in enumerate(names) if ":" in nm]
        b, v = res["beta"][idx], res["cov"][np.ix_(idx, idx)]
        w, q, g = float(b @ np.linalg.solve(v, b)), len(idx), int(res["n_groups"])
        out[key] = {"terms": [names[j] for j in idx], "wald": w, "df": q, "n_groups": g,
                    "p_chi2": float(stats.chi2.sf(w, q)), "p_F": float(stats.f.sf(w / q, q, g - 1))}
    return out


def digraph_matches(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    """13: ae・oe・ue を a・o・u と読んだ名でも照合したとき（元の読みとの和）に新たに付く照合（直したコードの辞書の照合）。"""
    lex = Lexicon(LEX)
    gids, counts = np.unique(t["P_gid1"].to_numpy(), return_counts=True)
    g = gn.drop_duplicates("gid").set_index("gid").reindex(gids)
    ccs = g["cc"].fillna("").astype(str).to_numpy()
    rows = []
    for gid, name, cc, n in zip(gids, g["asciiname"].fillna("").astype(str), ccs, counts, strict=True):
        low = name.lower()
        alt = DIGRAPH.sub(r"\1", low)
        if alt == low:
            continue
        m0, m1 = lex.match(low, cc), lex.match(alt, cc)
        w0, w1 = any(m0[c] for c in WET_SIDE), any(m1[c] for c in WET_SIDE)
        rows.append({"gid": int(gid), "cc": cc, "n_cells": int(n),
                     "gain_any": any(m1[c] and not m0[c] for c in CATEGORIES), "gain_wet": (not w0) and w1,
                     "lose_any": any(m0[c] and not m1[c] for c in CATEGORIES), "lose_wet": w0 and not w1})
    d = pd.DataFrame(rows)
    n_set, n_cells = len(gids), int(counts.sum())
    out: dict = {"n_settlements": n_set, "n_with_digraph": len(d), "n_cells": n_cells}
    for k in ("gain_any", "gain_wet", "lose_any", "lose_wet"):
        sel = d[k].to_numpy(dtype=bool)
        out[k] = {"n_settlements": int(sel.sum()), "share_settlements": float(sel.sum() / n_set),
                  "share_cells": float(d.loc[sel, "n_cells"].sum() / n_cells)}
    tot = pd.Series(ccs).value_counts()
    by = d.groupby("cc")[["gain_any", "gain_wet", "lose_wet"]].sum()
    by = by.join(tot.rename("n_settlements"), how="left")
    by = by[(by["n_settlements"] >= 50) & (by["gain_any"] > 0)].sort_values("gain_wet", ascending=False).head(12)
    out["by_country"] = {k: {"n_settlements": int(v["n_settlements"]), "gain_any": int(v["gain_any"]),
                             "gain_wet": int(v["gain_wet"]), "lose_wet": int(v["lose_wet"]),
                             "share_gain_wet": float(v["gain_wet"] / v["n_settlements"])} for k, v in by.iterrows()}
    out["examples_gain_wet"] = d.loc[d["gain_wet"], "gid"].head(15).tolist()
    return out


def language_recall(reps: int = 2000) -> dict:
    """14: 語源の照合（v44 の符号）の、辞書の水系側の読みの言語別の再現率（標本全体の層の重みの部分集団の推定）。"""
    samp = json.loads(v44.SAMPLE.read_text(encoding="utf-8"))
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    codes = pd.read_csv(v44.CODES)
    codes["gid"] = codes["id"].map(key)
    d = codes.join(meta, on="gid")
    d = d[d["status"] == "coded"].reset_index(drop=True)
    d["truth"] = d["cats"].fillna("").map(lambda s: bool({c for c in s.split(";") if c} & set(WET_SIDE)))
    d["flag"] = d["lex_wetside"].astype(bool)
    langs = [lang for lang, x in d.groupby("lang") if len(x) >= v44.MIN_GROUP]

    def recall(x: pd.DataFrame) -> dict[str, float]:
        w = x["stratum"].map({s: sizes[s] / max(1, int((x["stratum"] == s).sum())) for s in sizes})
        out = {}
        for lang in langs:
            m = (x["lang"] == lang).to_numpy()
            y, f, ww = x["truth"].to_numpy()[m], x["flag"].to_numpy()[m], w.to_numpy()[m]
            out[lang] = float((ww * (y & f)).sum() / (ww * y).sum()) if y.any() else float("nan")
        return out

    point = recall(d)
    rng = np.random.default_rng(SEED)
    parts = [d[d["stratum"] == s] for s in sizes]
    boots = [recall(pd.concat([p.iloc[rng.integers(0, len(p), len(p))] for p in parts], ignore_index=True))
             for _ in range(reps)]
    out: dict = {}
    for lang in langs:
        x = d[d["lang"] == lang]
        b = np.array([bb[lang] for bb in boots])
        out[lang] = {"recall": point[lang], "ci95": [float(np.nanpercentile(b, 2.5)), float(np.nanpercentile(b, 97.5))],
                     "n_coded": len(x), "n_true": int(x["truth"].sum()), "n_true_flagged": int((x["truth"] & x["flag"]).sum())}
    return {"by_language": out, "codes": "language-model codes of v44 (not independent human codes)", "reps": reps}



def heldout_coef_test() -> dict:
    """15（追記 3）: 直したコードの表の辞書の特徴 13 列の係数の同時 Wald 検定（直した実行の M1c の held-out の logit を与え、
    AOI の固定効果、場所でクラスタ頑健）。表の行の並び・y・AOI が直した実行の予測と一致することを先に確かめる。"""
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    sens = ev.Sensitivity()
    dirs, _, infos = ev.select_aois(gn, sens)
    df = ev.build_table(dirs, seed=20261025, infos=infos, sens=sens)
    df, _, _ = ev.add_channels(df, gn, None)
    cont = continent_map(ev.ROOT / "geonames" / "countryInfo.txt")
    df["fold"] = ev.folds_for("primary", df, infos, cont, buffer_km=sens.buffer_km).to_numpy()
    cols = [*(f"lex_{c}_n1" for c in CATEGORIES), *(f"lex_{c}_k" for c in CATEGORIES), "lex_wetside_n1"]
    t = pd.concat([df[df["fold"] == k][["aoi", "unit", "y", *cols]] for k in sorted(df["fold"].unique()) if k >= 0],
                  ignore_index=True)
    pred = RUNS_DIR / "world_v2_lexfix" / "predictions" / "primary"
    parts = [np.load(f) for f in sorted(pred.glob("fold[0-9].npz"), key=lambda f: int(f.stem[4:]))]
    y = np.concatenate([q["y"] for q in parts])
    aoi = np.concatenate([q["aoi"] for q in parts])
    if not (np.array_equal(y, t["y"].to_numpy()) and np.array_equal(aoi, t["aoi"].to_numpy(dtype=str))):
        msg = "rebuilt table does not match the corrected run's predictions"
        raise SystemExit(msg)
    p = np.clip(np.concatenate([q["p_M1c"] for q in parts]).astype(np.float64), 1e-5, 1 - 1e-5)
    yv = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(yv).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    z = np.column_stack([np.log(p / (1 - p)), t[cols].to_numpy(dtype=np.float64)])[keep]
    res = v39.fit(yv[keep], z, t["aoi"].to_numpy()[keep], t["unit"].to_numpy()[keep])
    idx = list(range(1, z.shape[1]))
    b, v = res["beta"][idx], res["cov"][np.ix_(idx, idx)]
    w, q, g = float(b @ np.linalg.solve(v, b)), len(idx), int(res["n_groups"])
    return {"features": cols, "coef": [float(x) for x in b], "se": [float(x) for x in res["se"][idx]],
            "coef_logit_M1c": float(res["beta"][0]), "wald": w, "df": q, "n_groups": g,
            "n_cells": int(keep.sum()), "n_aoi": int(res["n_aoi"]),
            "p_chi2": float(stats.chi2.sf(w, q)), "p_F": float(stats.f.sf(w / q, q, g - 1)),
            "predictions": "world_v2_lexfix p_M1c (held out)", "table": "rebuilt with the corrected code"}

# %% main
def main() -> None:
    t0 = time.time()
    t = table()
    gn, infos, dirs = _infos()
    print(f"table {len(t):,} cells ({time.time() - t0:.0f}s)", flush=True)
    res: dict = {"rule": "amendment B-7 items 5-11 (local)", "seed": SEED, "B_pairs": B_PAIRS,
                 "B_settlement": B_SET}
    g = gfm_table(t)
    cci = ("cci_built", "cci_tree")
    res["odds_ems"] = {"wetside": odds(t), "floodplain": odds(t, cls="floodplain"),
                       "wetside_adjusted_cci": odds(t, extra=cci),
                       "floodplain_adjusted_cci": odds(t, cls="floodplain", extra=cci)}
    res["odds_gfm"] = {"wetside": odds(g, "y_gfm"), "floodplain": odds(g, "y_gfm", cls="floodplain"),
                       "wetside_adjusted_cci": odds(g, "y_gfm", extra=cci)}
    print(f"odds done ({time.time() - t0:.0f}s)", flush=True)
    res["settlement"] = settlement_results(settlement_points())
    res["partial"] = partial(t, infos)
    res["severity"] = severity(t)
    res["gfm_discrimination"] = {"lexfix": gfm_discrimination(g, "p_lf_"), "frozen": gfm_discrimination(g, "p_")}
    print(f"items 5-9 done ({time.time() - t0:.0f}s)", flush=True)
    res["null_outliers"] = null_outliers()
    res["excluded_area"] = excluded_area(gn, infos, dirs)
    res["digraphs"] = digraphs(t, gn)
    res["block_precision"] = block_precision()
    res["residual_between_aois"] = residual_between(t)
    res["europe"] = europe(t)
    res["h3_slope"] = h3_slope()
    res["fixed_threshold_csi"] = fixed_threshold_csi(t)
    print(f"item 11 done ({time.time() - t0:.0f}s)", flush=True)
    res["wald_tests"] = wald_tests(t)
    res["digraph_matches"] = digraph_matches(t, gn)
    res["language_recall"] = language_recall()
    res["heldout_coef_test"] = heldout_coef_test()
    print(f"items 12-15 done ({time.time() - t0:.0f}s)", flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    res["pairs_bootstrap"] = {"wetside": pairs_bootstrap(t, "wetside", ("wetside",)),
                              "classes": pairs_bootstrap(t, "classes", tuple(CATEGORIES))}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps(res, indent=1, default=float)[:12000])
    print(f"done ({time.time() - t0:.0f}s) -> {OUT}")


if __name__ == "__main__":
    main()
