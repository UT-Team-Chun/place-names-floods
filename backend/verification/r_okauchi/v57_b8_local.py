# %% v57（事後・amendment B-8 の B のうち、既存の成果物から手元で計算するもの）
"""確証判定は変えない。B-8 で計算の前に固定した量を計算する（B1 とクラスタの実行の集計は v59、B13 は v57b）。

- B2: 達成した同等性の境界、B3: 決着しなかった判定に要る場所グループの数、B4: 置換による検定
- B5: GFM（同じセルでのラベルの違いの pairs bootstrap、jackknife の区間、Sentinel-1 以外の AOI、氾濫原の読み、
  当て直しの失敗の診断、指定した比較の ΔBS）
- B6: 集落の近くのオッズ比の版（型のある AOI、水系の層の無い AOI を除く、入植社会、掃除した読み、上限の無い標本、
  誤分類の補正、区分ごとの調整、川の大きさ、割合の結果）
- B7: 本川からの高さと最寄りの河道の集水面積（正の対照）、B8: 較正、B9: 地形の基礎の技能、B11: 辞書の監査、
  B12: 言語に依らない水名の読み

オッズ比は v39・v55 と同じ形（AOI の固定効果、凍結した実行の M1 の held-out 予測の logit、距離帯、場所でクラスタ
頑健な標準誤差）。pairs bootstrap は場所グループの復元抽出（2,000 回、種 20261028、v55 と同じ）。
出力: stdout と docs/verification/v57_b8_local.json（重い中間の表は data/verification/v57_*.parquet・npy）
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyproj import CRS, Transformer
from scipy import stats
from scipy.spatial import cKDTree

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import (
    HYDRONYM_CODES,
    HYDRONYM_RADIUS_M,
    contains_hydronym,
    correct_misclassification,
    frequent_tokens,
    groups_needed,
    hydronym_cores,
    jackknife_se,
)
from topsem.world.lexicon import CATEGORIES, WET_SIDE, Lexicon, _hit

OUT = ROOT / "docs" / "verification" / "v57_b8_local.json"
VER = DATA_DIR / "verification"
LABELS = DATA_DIR / "global" / "labels"
SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
FROZEN_LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
CLEAN_LEX = ROOT / "docs" / "verification" / "v57_lexicon_clean.yaml"
SEED = 20261028
B_PAIRS = 2000
B_PBA = 20000
TERM = "le500m"
UPLAND_CUTS = (100.0, 1000.0)  # B6(h)・B7: 最寄りの河道の集水面積の層の境（km²）
N_WORKERS = int(os.environ.get("V57_WORKERS", "8"))  # 結果に影響しない（再標本は先に引く）。手元のメモリに合わせて減らせる
EARTH_M = 6_371_008.8


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v55 = _load("v55_b7_local")
v39, v37, v41, v50, ev = v55.v39, v55.v37, v55.v41, v55.v50, v55.ev


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# %% B2・B3・B4
def b2_achieved() -> dict:
    out: dict = {}
    for run in ("world_v2", "world_v2_lexfix"):
        p = _json(RUNS_DIR / run / "metrics_primary.json")["primary"]
        lo, hi = p["ci90_R"]
        out[run] = {"R": p["R"], "ci90_R": [lo, hi], "achieved_bound_R": max(abs(lo), abs(hi)),
                    "achieved_bound_brier": max(abs(x) for x in p["ci90_delta"])}
    v36 = _json(ROOT / "docs" / "verification" / "v36_yardstick_robustness.json")["splits"]["primary"]
    combos = {f"{b}/{y}": max(abs(v36[b]["primary"][y]["ci90_R"][0]), abs(v36[b]["primary"][y]["ci90_R"][1]))
              for b in ("cluster", "two_way") for y in ("Y0", "Y1", "Y2", "Y3", "Y4")}
    k = max(combos, key=combos.get)
    out["combos"] = combos
    out["max_over_combos"] = {"which": k, "achieved_bound_R": combos[k]}
    return out


def b3_needed() -> dict:
    m = _json(RUNS_DIR / "world_v2" / "metrics_primary.json")
    near = _json(ROOT / "docs" / "verification" / "v37_positive_control.json")["primary_by_distance_wetside"][TERM]
    rows = {"total_M1s_M3s": (m["total"], 195), "flexible_M1_M3": (m["flexible"], 195),
            "near_500m_wetside_not_planned": (near, int(near["n_units"]))}
    z = stats.norm.ppf(0.975) / stats.norm.ppf(0.95)
    out = {}
    for k, (d, n) in rows.items():
        est, (lo90, hi90) = float(d["R"]), d["ci90_R"]
        if "ci95_R" in d:
            lo95, hi95 = d["ci95_R"]
            how95 = "reported"
        else:  # 95% 区間が無ければ 90% 区間の半幅を正規の比で広げる
            h = (hi90 - lo90) / 2 * z
            lo95, hi95, how95 = est - h, est + h, "scaled from the 90% interval"
        ge, gp = groups_needed(est, lo90, hi90, n, margin=0.05), groups_needed(est, lo95, hi95, n)
        out[k] = {"R": est, "ci90_R": [lo90, hi90], "ci95_R": [lo95, hi95], "ci95_how": how95, "n_units": n,
                  "groups_equivalence": ge, "ratio_equivalence": ge / n, "groups_positive": gp,
                  "ratio_positive": gp / n}
    return out


def b4_perm() -> dict:
    d = RUNS_DIR / "world_v2_b7perm" / "predictions" / "primary"
    files = sorted(d.glob("fold*_permsrc_*.npz"))
    v56 = _json(ROOT / "docs" / "verification" / "v56_b7_cluster.json")["permutations"]
    yard = float(v56["yardstick_lexfix"])
    zs = [np.load(f) for f in files]
    n = float(sum(int(z["n_unit"].sum()) for z in zs)) / (len(files) / 10)  # 範囲ごとに同じセルを数える
    by_first: dict[str, list] = {}
    for f, z in zip(files, zs, strict=True):
        by_first.setdefault(f.stem.split("_")[-1], []).append(z)
    ref, nam, real_c, real_n = [], [], None, None
    for _, parts in sorted(by_first.items()):
        ref.append(sum(z["sse_ref"].sum(axis=1) for z in parts))
        nam.append(sum(z["sse_names"].sum(axis=1) for z in parts))
        if real_c is None:
            real_c = sum(float(z["real_M1c"].sum()) for z in parts)
            real_n = sum(float(z["real_M3n"].sum()) for z in parts)
    dq = (np.concatenate(ref) - np.concatenate(nam)) / n
    d_real = (real_c - real_n) / n
    r_q, r_real = dq / yard, d_real / yard
    k = int((dq >= d_real).sum())
    return {"n_perms": len(dq), "R_real": r_real, "R_perm_q95": float(np.quantile(r_q, 0.95)),
            "R_perm_mean": float(r_q.mean()), "p_one_sided": (1 + k) / (1 + len(dq)),
            "check_vs_v56_mean": float(v56["recomputed_composition"]["R_perm_mean"]),
            "label_with_permutation_placebo": "unchanged (equivalent-null): the placebo condition applies only "
                                              "when the 95% interval of R lies above zero"}


# %% 表（v55 の表 + 掃除した読み・水名の読み・集落の川）
def analysis_table() -> tuple[pd.DataFrame, pd.DataFrame, dict, list]:
    t = v55.table()
    gn, infos, dirs = v55._infos()  # noqa: SLF001
    return t, gn, infos, dirs


def settlement_lex(gids: np.ndarray, gn: pd.DataFrame, lex: Lexicon) -> pd.DataFrame:
    """集落ごとの辞書の区分（名前の最寄りの 1 つ、n1 と同じ照合）。"""
    g = gn.drop_duplicates("gid").set_index("gid").reindex(gids)
    rows = []
    for gid, name, cc in zip(gids, g["asciiname"].fillna("").to_numpy(), g["cc"].fillna("").to_numpy(), strict=True):
        m = lex.match(str(name), str(cc))
        rows.append([gid, *(int(m.get(c, 0) > 0) for c in CATEGORIES)])
    s = pd.DataFrame(rows, columns=["gid", *CATEGORIES])
    s["wetside"] = s[list(WET_SIDE)].max(axis=1)
    return s.set_index("gid")


# %% B5: GFM
_W: dict = {}


def _pairs2_init(y1: np.ndarray, y2: np.ndarray, z: np.ndarray, aoi: np.ndarray, rows: list, na: int, j: int) -> None:
    _W.update(y1=y1, y2=y2, z=z, aoi=aoi, rows=rows, na=na, j=j)


def _pairs2_rep(sel: np.ndarray) -> tuple[float, float] | None:
    rows = _W["rows"]
    idx = np.concatenate([rows[u] for u in sel])
    pos = np.concatenate([np.full(len(rows[u]), k, dtype=np.int64) for k, u in enumerate(sel)])
    a = _W["aoi"][idx] + pos * _W["na"]
    try:
        b1 = v39.fit(_W["y1"][idx], _W["z"][idx], a, pos)["beta"][_W["j"]]
        b2 = v39.fit(_W["y2"][idx], _W["z"][idx], a, pos)["beta"][_W["j"]]
    except np.linalg.LinAlgError:
        return None
    return float(b1), float(b2)


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def same_cell_difference(tg: pd.DataFrame, cls: str = "wetside") -> dict:
    """(a) 同じセルで EMS と GFM の水系側 × 500 m の係数の差（GFM − EMS）。AOI の採否は両方のラベルで変化のある AOI。"""
    cache = VER / f"v57_gfm_pairs_{cls}.json"
    if cache.is_file():
        return _json(cache)
    ye = tg["y_ems"].to_numpy(dtype=np.float64)
    yg = tg["y_gfm"].to_numpy(dtype=np.float64)
    re_ = pd.Series(ye).groupby(tg["aoi"].to_numpy()).transform("mean").to_numpy()
    rg = pd.Series(yg).groupby(tg["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (re_ > 0) & (re_ < 1) & (rg > 0) & (rg < 1)
    s = tg.loc[keep].reset_index(drop=True)
    ye, yg = ye[keep], yg[keep]
    z, names = v39.design(s, (cls,), True, banded=True)
    j = names.index(f"{cls}:{TERM}")
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    g, na = len(unit_lab), int(aoi.max()) + 1
    rows = [np.flatnonzero(unit == u) for u in range(g)]
    fe, fg = v39.fit(ye, z, aoi, unit), v39.fit(yg, z, aoi, unit)
    rng = np.random.default_rng(SEED)
    sels = [rng.integers(0, g, size=g) for _ in range(B_PAIRS)]
    _threads1()
    with ProcessPoolExecutor(N_WORKERS, initializer=_pairs2_init, initargs=(ye, yg, z, aoi, rows, na, j)) as ex:
        res = [r for r in ex.map(_pairs2_rep, sels, chunksize=8) if r is not None]
    b = np.array(res)
    d = b[:, 1] - b[:, 0]
    # (b) 場所グループを一つずつ除いた推定からの jackknife
    loo = np.full((g, 2), np.nan)
    for u in range(g):
        kp = unit != u
        try:
            loo[u, 0] = v39.fit(ye[kp], z[kp], aoi[kp], unit[kp])["beta"][j]
            loo[u, 1] = v39.fit(yg[kp], z[kp], aoi[kp], unit[kp])["beta"][j]
        except np.linalg.LinAlgError:
            continue
    ok = np.isfinite(loo).all(axis=1)
    tq = stats.t.ppf(0.975, int(ok.sum()) - 1)
    out = {"n_cells": len(s), "n_aoi": int(s["aoi"].nunique()), "n_units": g, "reps_ok": len(b)}
    for k, (lab, fit) in enumerate((("ems", fe), ("gfm", fg))):
        b0 = float(fit["beta"][j])
        se_j = jackknife_se(loo[ok, k])
        out[lab] = {"odds_ratio": float(np.exp(b0)),
                    "sandwich_ci95": [float(np.exp(b0 - 1.96 * fit["se"][j])), float(np.exp(b0 + 1.96 * fit["se"][j]))],
                    "pairs_ci95": [float(np.exp(np.quantile(b[:, k], 0.025))),
                                   float(np.exp(np.quantile(b[:, k], 0.975)))],
                    "jackknife_se": se_j, "jackknife_t_ci95": [float(np.exp(b0 - tq * se_j)),
                                                               float(np.exp(b0 + tq * se_j))]}
    d0 = float(fg["beta"][j] - fe["beta"][j])
    out["log_or_difference_gfm_minus_ems"] = {"estimate": d0, "ratio_of_odds_ratios": float(np.exp(d0)),
                                              "pairs_ci95": [float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))],
                                              "pairs_ci95_ratio": [float(np.exp(np.quantile(d, 0.025))),
                                                                   float(np.exp(np.quantile(d, 0.975)))]}
    cache.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def _auc(y: np.ndarray, p: np.ndarray) -> float:
    return v55._auc(y, p)  # noqa: SLF001


def _match_cells(a: pd.DataFrame, b: pd.DataFrame) -> np.ndarray:
    """a の各行に、同じ AOI で経緯度が最も近い b の行（1e-5° 以内、無ければ −1）。"""
    idx = np.full(len(a), -1, dtype=np.int64)
    for aoi, ia in a.groupby("aoi").indices.items():
        ib = np.flatnonzero(b["aoi"].to_numpy() == aoi)
        if not len(ib):
            continue
        tree = cKDTree(b.iloc[ib][["lon", "lat"]].to_numpy(np.float64))
        dist, j = tree.query(a.iloc[ia][["lon", "lat"]].to_numpy(np.float64))
        ok = dist < 1e-5
        idx[ia[ok]] = ib[j[ok]]
    return idx


def gfm_refit_diagnosis(tg: pd.DataFrame) -> dict:
    """(e) GFM で当て直したモデルと、EMS で学習した M1 を、同じセル・GFM のラベルで比べる。"""
    d = RUNS_DIR / "world_v2_gfm" / "predictions" / "primary"
    parts = [np.load(f, allow_pickle=False) for f in sorted(d.glob("fold[0-9].npz"))]
    cols = ("p_M0", "p_M1", "p_M1s", "p_M1c")
    r = pd.DataFrame({"aoi": np.concatenate([q["aoi"] for q in parts]).astype(str),
                      "lon": np.concatenate([q["lon"] for q in parts]).astype(np.float64),
                      "lat": np.concatenate([q["lat"] for q in parts]).astype(np.float64),
                      "y": np.concatenate([q["y"] for q in parts]).astype(np.float64),
                      **{c: np.concatenate([q[c] for q in parts]).astype(np.float64) for c in cols}})
    idx = _match_cells(r, tg)
    ok = idx >= 0
    r = r.loc[ok].reset_index(drop=True)
    r["p_EMS_M1"] = tg["p_M1"].to_numpy(np.float64)[idx[ok]]
    r["y_check"] = tg["y_gfm"].to_numpy(np.float64)[idx[ok]]
    out: dict = {"n_cells_refit": int(ok.size), "n_matched": int(ok.sum()),
                 "labels_equal": bool(np.array_equal(r["y"].to_numpy(), r["y_check"].to_numpy()))}
    y = r["y"].to_numpy()
    for m in (*cols, "p_EMS_M1"):
        p = r[m].to_numpy()
        aucs = []
        for _, x in r.groupby("aoi"):
            yy = x["y"].to_numpy()
            if 0 < yy.sum() < len(yy):
                aucs.append(_auc(yy, x[m].to_numpy()))
        g = r.assign(_p=p).groupby("aoi").agg(n=("y", "size"), ybar=("y", "mean"), pbar=("_p", "mean"))
        between = float((g["n"] * (g["pbar"] - g["ybar"]) ** 2).sum() / len(r))
        brier = float(np.mean((p - y) ** 2))
        out[m] = {"brier": brier, "between_aoi": between, "within_aoi": brier - between,
                  "auc_within_median": float(np.median(aucs)), "auc_within_q25_q75": np.quantile(aucs, [0.25, 0.75]).tolist(),
                  "n_aoi_auc": len(aucs)}
    a = r.groupby("aoi")["y"].mean()
    e = tg.groupby("aoi")["y_ems"].mean().reindex(a.index)
    out["aoi_flooded_share"] = {"gfm_quantiles": a.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).tolist(),
                                "ems_quantiles": e.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).tolist(),
                                "gfm_sd": float(a.std()), "ems_sd": float(e.std()), "n_aoi": len(a)}
    diff = out["p_M1"]["auc_within_median"] - out["p_EMS_M1"]["auc_within_median"]
    out["rule"] = ("learned the terrain ordering within AOIs; the AOI-level flooded share did not transfer"
                   if abs(diff) <= 0.05 else "could not learn the relation with terrain from the GFM labels")
    out["auc_difference_refit_minus_ems"] = diff
    return out


def b5_gfm(t: pd.DataFrame) -> dict:
    summary = _json(SUMMARY)
    tg = v41.attach_gfm(t, summary)
    out: dict = {"same_cells": same_cell_difference(tg)}
    out["floodplain_same_cells"] = {"ems": v55.odds(tg, "y_ems", "floodplain"), "gfm": v55.odds(tg, "y_gfm", "floodplain")}
    v41j = _json(ROOT / "docs" / "verification" / "v41_gfm_independent_labels.json")
    out["non_s1"] = v41j.get("sensitivity_independent_sensor")
    out["designated_c2"] = v41j.get("C2")
    out["refit_diagnosis"] = gfm_refit_diagnosis(tg)
    return out


# %% B6: 集落の近くの版
def pairs_draws(t: pd.DataFrame, cls: str, key: str) -> np.ndarray:
    """主の推定の pairs bootstrap の係数（水系側 × 500 m）を保存して使い回す（v55 と同じ再標本）。"""
    cache = VER / f"v57_pairs_draws_{key}.npy"
    if cache.is_file():
        return np.load(cache)
    y = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = t.loc[keep].reset_index(drop=True)
    y = y[keep]
    z, names = v39.design(s, (cls,), True, banded=True)
    j = names.index(f"{cls}:{TERM}")
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    g, na = len(unit_lab), int(aoi.max()) + 1
    rows = [np.flatnonzero(unit == u) for u in range(g)]
    rng = np.random.default_rng(SEED)
    sels = [rng.integers(0, g, size=g) for _ in range(B_PAIRS)]
    _threads1()
    with ProcessPoolExecutor(N_WORKERS, initializer=v55._pairs_init, initargs=(y, z, aoi, rows, na)) as ex:  # noqa: SLF001
        res = list(ex.map(v55._pairs_rep, sels, chunksize=8))  # noqa: SLF001
    b = np.array([r[j] if r is not None else np.nan for r in res])
    np.save(cache, b)
    return b


def _with_pairs(t: pd.DataFrame, cls: str, key: str, mask: np.ndarray | None = None) -> dict:
    s = t if mask is None else t.loc[mask].reset_index(drop=True)
    o = v55.odds(s, cls=cls)
    b = pairs_draws(s, cls, key)
    ok = np.isfinite(b)
    o["pairs_ci95"] = [float(np.exp(np.quantile(b[ok], 0.025))), float(np.exp(np.quantile(b[ok], 0.975)))]
    o["pairs_reps_ok"] = int(ok.sum())
    return o


def misclassification(t: pd.DataFrame) -> dict:
    v44 = _json(ROOT / "docs" / "verification" / "v44_wikipedia_etymology.json")
    w, o = v44["truth_rate_in_wet"], v44["truth_rate_in_not"]
    b = pairs_draws(t, "wetside", "main_wetside")
    b = b[np.isfinite(b)]
    rng = np.random.default_rng(SEED)
    orr = np.exp(rng.choice(b, size=B_PBA, replace=True))
    p1 = rng.beta(w["k"] + 1, w["n"] - w["k"] + 1, size=B_PBA)
    p0 = rng.beta(o["k"] + 1, o["n"] - o["k"] + 1, size=B_PBA)
    c = correct_misclassification(orr, p1, p0)
    ok = np.isfinite(c)
    point = float(correct_misclassification(np.exp(np.median(b)), w["rate"], o["rate"]))
    return {"p1": w, "p0": o, "point_at_rates": point, "median": float(np.median(c[ok])),
            "ci95": [float(np.quantile(c[ok], 0.025)), float(np.quantile(c[ok], 0.975))],
            "share_undefined": float(1 - ok.mean()), "draws": B_PBA,
            "note": "p1 and p0 come from the codes of the model that compiled the lexicon"}


def class_adjusted(t: pd.DataFrame) -> dict:
    """(g) 区分ごとのオッズ比（区分を一度に入れた v39 のモデル）、Holm の p と Bonferroni の 95% 区間。"""
    y = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, y = t.loc[keep].reset_index(drop=True), y[keep]
    z, names = v39.design(s, tuple(CATEGORIES), True, banded=True)
    res = v39.fit(y, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
    k = len(CATEGORIES)
    zb = stats.norm.ppf(1 - 0.05 / (2 * k))
    rows = []
    for c in CATEGORIES:
        j = names.index(f"{c}:{TERM}")
        b, se = float(res["beta"][j]), float(res["se"][j])
        rows.append((c, b, se, float(2 * stats.norm.sf(abs(b / se)))))
    order = np.argsort([r[3] for r in rows])
    holm = np.zeros(k)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (k - rank) * rows[i][3]))
        holm[i] = run
    return {c: {"odds_ratio": float(np.exp(b)), "ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                "bonferroni_ci95": [float(np.exp(b - zb * se)), float(np.exp(b + zb * se))], "p": p,
                "p_holm": float(holm[i])} for i, (c, b, se, p) in enumerate(rows)}


def fractional_all_bands(t: pd.DataFrame, infos: dict) -> dict:
    untyped = t["aoi"].map({a: bool(infos[a].unspecified_as_riverine) for a in t["aoi"].unique()}).to_numpy()
    riv = np.clip(t["f_riverine"].to_numpy() + np.where(untyped, t["f_other"].to_numpy(), 0.0), 0.0, 1.0)
    rate = pd.Series(riv).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, y = t.loc[keep].reset_index(drop=True), riv[keep]
    z, names = v39.design(s, ("wetside",), True, banded=True)
    res = v39.fit(y, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
    return {b: v39.table(res, names)[f"wetside:{b}"] for b, *_ in v39.BANDS}


def uncapped_table(t: pd.DataFrame) -> pd.DataFrame:
    """(e) 上限の無い標本: 上限なしの感度の実行の held-out の予測と、各 AOI の特徴の表の最寄りの集落の距離。"""
    cache = VER / "v57_uncapped.parquet"
    if cache.is_file():
        return pd.read_parquet(cache)
    d = RUNS_DIR / "world_v2_s2a_nocap" / "predictions" / "primary"
    parts = []
    for f in sorted(d.glob("fold[0-9].npz")):
        q = np.load(f, allow_pickle=False)
        parts.append(pd.DataFrame({"aoi": q["aoi"].astype(str), "unit": q["unit"].astype(str),
                                   "y": q["y"].astype(np.int8), "p_M1": q["p_M1"].astype(np.float64),
                                   "lex_wetside_n1": q["s_wetside"].astype(np.float64),
                                   "lon": q["lon"].astype(np.float64), "lat": q["lat"].astype(np.float64)}))
    u = pd.concat(parts, ignore_index=True)
    u["P_d1"] = np.nan
    for aoi, ia in u.groupby("aoi").indices.items():
        c = pq.read_table(LABELS / aoi / "cells.parquet", columns=["lon", "lat"]).to_pandas()
        f = pq.read_table(LABELS / aoi / "features.parquet", columns=["P_d1"]).to_pandas()
        tree = cKDTree(c[["lon", "lat"]].to_numpy(np.float64))
        dist, j = tree.query(u.iloc[ia][["lon", "lat"]].to_numpy(np.float64))
        v = f["P_d1"].to_numpy(np.float64)[j]
        v[dist > 1e-5] = np.nan
        u.loc[u.index[ia], "P_d1"] = v
    u = u.loc[u["P_d1"].notna()].reset_index(drop=True)
    u.to_parquet(cache, index=False)
    return u


def b6_variants(t: pd.DataFrame, gn: pd.DataFrame, infos: dict, dirs: list) -> dict:
    out: dict = {"main": _with_pairs(t, "wetside", "main_wetside")}
    untyped = t["aoi"].map({a: bool(infos[a].unspecified_as_riverine) for a in t["aoi"].unique()}).to_numpy()
    out["typed_aois"] = _with_pairs(t, "wetside", "typed", ~untyped)
    nohydro = {d.name for d in dirs if not (_json(d / "meta.json").get("layers") or {}).get("hydro")}
    m_nh = ~t["aoi"].isin(nohydro).to_numpy()
    out["without_no_hydro"] = _with_pairs(t, "wetside", "nohydro", m_nh) | {"n_aoi_dropped": len(nohydro)}
    settler = t["continent"].isin(["OC", "NA", "SA"]).to_numpy()
    out["without_oceania_americas"] = _with_pairs(t, "wetside", "nosettler", ~settler)
    out["without_australia"] = _with_pairs(t, "wetside", "noau", (t["cc"] != "AU").to_numpy())
    gids = np.unique(t["P_gid1"].to_numpy())
    clean = settlement_lex(gids, gn, Lexicon(CLEAN_LEX))
    t["lex_wetclean_n1"] = clean["wetside"].reindex(t["P_gid1"].to_numpy()).fillna(0).to_numpy(np.float64)
    out["cleaned_reading"] = _with_pairs(t, "wetclean", "clean") | {
        "share_cells_changed": float(((t["lex_wetclean_n1"] > 0) != (t["lex_wetside_n1"] > 0)).mean())}
    u = uncapped_table(t)
    out["uncapped"] = _with_pairs(u, "wetside", "uncapped") | {"note": "frozen-code readings of the uncapped run"}
    out["misclassification_corrected"] = misclassification(t)
    out["class_adjusted"] = class_adjusted(t)
    out["fractional_all_bands"] = fractional_all_bands(t, infos)
    e = out["uncapped"]["pairs_ci95"]
    same = all(np.sign(np.log(out[k]["odds_ratio"])) == np.sign(np.log(out["uncapped"]["odds_ratio"]))
               for k in ("typed_aois", "cleaned_reading"))
    out["rule"] = {"uncapped_excludes_one": bool(e[0] > 1 or e[1] < 1), "typed_and_cleaned_same_direction": bool(same),
                   "conclusion_changes": bool((e[0] > 1 or e[1] < 1) and same)}
    return out


# %% B7・B6(h): HydroRIVERS
def _unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def river_vertices(boxes: list[tuple[float, float, float, float]]) -> pd.DataFrame:
    """AOI の範囲（± 0.2°）にかかる HydroRIVERS の全ての河道の頂点と、その河道の集水面積（UPLAND_SKM）。"""
    cache = VER / "v57_river_vertices.parquet"
    if cache.is_file():
        return pd.read_parquet(cache)
    import pyogrio.raw
    import shapely

    b = np.array(boxes)
    pad = 0.2
    parts = []
    for z in sorted((DATA_DIR / "global" / "hydrorivers").glob("HydroRIVERS_v10_*_shp.zip")):
        with zipfile.ZipFile(z) as zf:
            shp = next(n for n in zf.namelist() if n.endswith(".shp"))
        _, _, wkb, fields = pyogrio.raw.read(f"/vsizip/{z}/{shp}", columns=["UPLAND_SKM"])
        geoms = shapely.from_wkb(wkb)
        bounds = shapely.bounds(geoms)
        hit = np.zeros(len(geoms), dtype=bool)
        for w, s, e, n in b:
            hit |= ((bounds[:, 2] >= w - pad) & (bounds[:, 0] <= e + pad) & (bounds[:, 3] >= s - pad)
                    & (bounds[:, 1] <= n + pad))
        xy, ix = shapely.get_coordinates(geoms[hit], return_index=True)
        up = np.asarray(fields[0], dtype=np.float64)[hit][ix]
        parts.append(np.column_stack([xy[:, 0], xy[:, 1], up]))
        print(f"  {z.name}: {int(hit.sum())} reaches, {len(xy)} vertices", flush=True)
    v = pd.DataFrame(np.concatenate(parts), columns=["lon", "lat", "upland"])
    v.to_parquet(cache, index=False)
    return v


def settlement_rivers(t: pd.DataFrame, gn: pd.DataFrame, dirs: list) -> pd.DataFrame:
    """集落（AOI ごと）の最寄りの河道の集水面積と、集水面積 ≥ 100・≥ 1000 km² の河道からの高さ。"""
    cache = VER / "v57_settlement_rivers.parquet"
    if cache.is_file():
        return pd.read_parquet(cache)
    dmap = {d.name: d for d in dirs}
    boxes = [tuple(v) for v in t.groupby("aoi").agg(w=("lon", "min"), s=("lat", "min"), e=("lon", "max"),
                                                       n=("lat", "max")).to_numpy()]
    rv = river_vertices(boxes)
    trees = {"any": (cKDTree(_unit_xyz(rv["lon"].to_numpy(), rv["lat"].to_numpy())), rv)}
    for cut in UPLAND_CUTS:
        sub = rv.loc[rv["upland"] >= cut].reset_index(drop=True)
        trees[str(int(cut))] = (cKDTree(_unit_xyz(sub["lon"].to_numpy(), sub["lat"].to_numpy())), sub)
    g = gn.drop_duplicates("gid").set_index("gid")
    rows = []
    for aoi, x in t.groupby("aoi"):
        meta = _json(dmap[aoi] / "meta.json")
        grid = meta["grid"]
        feat = pq.read_table(dmap[aoi] / "features.parquet", columns=["elev_mean", "elev_min"]).to_pandas()
        cells = pq.read_table(dmap[aoi] / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        pos = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(cells["ix"], cells["iy"], strict=True))}
        tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(meta["crs_proj4"]), always_xy=True)

        def cell_of(lon: np.ndarray, lat: np.ndarray, grid: dict = grid, tr: Transformer = tr,
                    pos: dict = pos) -> np.ndarray:
            xx, yy = tr.transform(lon, lat)
            cx = np.floor((np.asarray(xx) - grid["x0"]) / grid["cell_m"]).astype(np.int64)
            cy = np.floor((grid["y1"] - np.asarray(yy)) / grid["cell_m"]).astype(np.int64)
            return np.array([pos.get((int(a), int(b)), -1) for a, b in zip(cx, cy, strict=True)])

        gids = np.unique(x["P_gid1"].to_numpy())
        gids = gids[gids >= 0]
        lon = g["lon"].reindex(gids).to_numpy(np.float64)
        lat = g["lat"].reindex(gids).to_numpy(np.float64)
        xyz = _unit_xyz(lon, lat)
        own = cell_of(lon, lat)
        el_own = np.where(own >= 0, feat["elev_mean"].to_numpy(np.float64)[np.maximum(own, 0)], np.nan)
        rec = {"aoi": aoi, "gid": gids, "elev_site": el_own}
        d, j = trees["any"][0].query(xyz)
        rec["upland_nearest"] = trees["any"][1]["upland"].to_numpy()[j]
        rec["dist_nearest_m"] = d * EARTH_M
        for cut in UPLAND_CUTS:
            tree, sub = trees[str(int(cut))]
            d, j = tree.query(xyz)
            rc = cell_of(sub["lon"].to_numpy()[j], sub["lat"].to_numpy()[j])
            el_r = np.where(rc >= 0, feat["elev_min"].to_numpy(np.float64)[np.maximum(rc, 0)], np.nan)
            rec[f"h_above_{int(cut)}"] = el_own - el_r
            rec[f"dist_{int(cut)}_m"] = d * EARTH_M
        rows.append(pd.DataFrame(rec))
    s = pd.concat(rows, ignore_index=True)
    s.to_parquet(cache, index=False)
    return s


def b7_rivers(t: pd.DataFrame, gn: pd.DataFrame, dirs: list) -> dict:
    sr = settlement_rivers(t, gn, dirs)
    sp = v55.settlement_points()
    sp = sp.merge(sr, left_on=["aoi", "P_gid1"], right_on=["aoi", "gid"], how="left")
    sp["log_upland"] = np.log(sp["upland_nearest"])
    out: dict = {"n_settlement_points": len(sp)}
    for var in ("h_above_100", "h_above_1000", "log_upland"):
        out[var] = {g: v55.fe_contrast(sp, var, g, ()) for g in ("wetside", "floodplain")}
        out[var]["by_class"] = {c: v55.fe_contrast(sp, var, c, ()) for c in WET_SIDE}
    # B6(h): 最寄りの集落の最寄りの河道の大きさの層ごとのオッズ比、本川からの高さを共変量に
    key = sr.set_index(["aoi", "gid"])
    idx = pd.MultiIndex.from_arrays([t["aoi"].to_numpy(), t["P_gid1"].to_numpy()])
    up = key["upland_nearest"].reindex(idx).to_numpy()
    h = key["h_above_1000"].reindex(idx).to_numpy()
    strata = {"lt100": up < UPLAND_CUTS[0], "100_1000": (up >= UPLAND_CUTS[0]) & (up < UPLAND_CUTS[1]),
              "ge1000": up >= UPLAND_CUTS[1]}
    out["odds_by_river_size"] = {k: v55.odds(t, mask=m) for k, m in strata.items()}
    tt = t.assign(h_main=h)
    ok = np.isfinite(h)
    out["odds_adjusted_h_main"] = v55.odds(tt, extra=("h_main",), mask=ok)
    out["odds_same_cells_unadjusted"] = v55.odds(tt, mask=ok)
    return out


# %% B8・B9: 較正と技能
def _calibration(y: np.ndarray, p: np.ndarray) -> dict:
    q = np.clip(p, 1e-6, 1 - 1e-6)
    x = np.column_stack([np.ones_like(q), np.log(q / (1 - q))])
    beta = np.zeros(2)
    for _ in range(50):
        eta = x @ beta
        mu = 1 / (1 + np.exp(-eta))
        w = mu * (1 - mu)
        step = np.linalg.solve(x.T @ (x * w[:, None]), x.T @ (y - mu))
        beta += step
        if np.abs(step).max() < 1e-10:
            break
    edges = np.quantile(p, np.linspace(0, 1, 11))
    b = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, 9)
    nb = np.bincount(b, minlength=10).astype(np.float64)
    pb = np.bincount(b, weights=p, minlength=10) / np.maximum(nb, 1)
    yb = np.bincount(b, weights=y, minlength=10) / np.maximum(nb, 1)
    ybar = y.mean()
    return {"intercept": float(beta[0]), "slope": float(beta[1]), "brier": float(np.mean((p - y) ** 2)),
            "reliability": float((nb * (pb - yb) ** 2).sum() / len(y)),
            "resolution": float((nb * (yb - ybar) ** 2).sum() / len(y)), "uncertainty": float(ybar * (1 - ybar))}


def _preds(run: str, split: str) -> pd.DataFrame | None:
    d = RUNS_DIR / run / "predictions" / split
    files = sorted(d.glob("fold[0-9].npz"))
    if not files:
        return None
    parts = [np.load(f, allow_pickle=False) for f in files]
    cols = [c for c in ("p_M1", "p_M1s", "p_M1c", "p_M3n") if c in parts[0].files]
    return pd.DataFrame({"aoi": np.concatenate([q["aoi"] for q in parts]).astype(str),
                         "y": np.concatenate([q["y"] for q in parts]).astype(np.float64),
                         **{c: np.concatenate([q[c] for q in parts]).astype(np.float64) for c in cols}})


def b8_calibration() -> dict:
    runs = {"primary": ("world_v2", "primary"), "country": ("world_v2", "country"),
            "temporal": ("world_v2", "temporal"), "llmcut": ("world_v2", "llmcut"),
            "activation_frozen": ("world_v2_actgrp", "primary"), "activation_corrected": ("world_v2_b7actfix", "primary"),
            "gfm_refit": ("world_v2_gfm", "primary")}
    out: dict = {}
    for k, (run, split) in runs.items():
        df = _preds(run, split)
        if df is None:
            out[k] = "missing"
            continue
        res: dict = {"n_cells": len(df), "n_aoi": int(df["aoi"].nunique())}
        for m in (c for c in ("p_M1", "p_M1s", "p_M1c", "p_M3n") if c in df):
            res[m[2:]] = _calibration(df["y"].to_numpy(), df[m].to_numpy())
            aucs = [_auc(x["y"].to_numpy(), x[m].to_numpy()) for _, x in df.groupby("aoi")
                    if 0 < x["y"].sum() < len(x)]
            res[m[2:]]["auc_within_median"] = float(np.median(aucs))
        out[k] = res
    return out


def b9_skill(t: pd.DataFrame) -> dict:
    rows = {"M1": [], "HAND": []}
    for _, x in t.groupby("aoi"):
        y = x["y"].to_numpy()
        k = int(y.sum())
        if k in (0, len(y)):
            continue
        for name, score in (("M1", x["p_M1"].to_numpy(np.float64)), ("HAND", -x["hand_min"].to_numpy(np.float64))):
            top = np.argsort(-score, kind="stable")[:k]
            tp = int(y[top].sum())
            rows[name].append((tp / (2 * k - tp), tp / k, (k - tp) / k, _auc(y, score)))
    out = {}
    for name, r in rows.items():
        a = np.array(r)
        out[name] = {m: {"median": float(np.median(a[:, i])), "q25_q75": np.quantile(a[:, i], [0.25, 0.75]).tolist()}
                     for i, m in enumerate(("csi", "hit_rate", "false_alarm_ratio", "auc_within"))}
        out[name]["n_aoi"] = len(a)
    return out


# %% B11: 辞書の監査
def b11_audit(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    lex = Lexicon(FROZEN_LEX)
    dropped = {tuple(x) for x in _json(ROOT / "docs" / "verification" / "v57_lexicon_clean_dropped.json")["dropped"]}
    watch = set(dropped) | {("tupi", f, "word") for f in ("parana", "iguacu", "paranapanema")} | \
        {("id", "krueng", "word"), ("id", "way", "word"), ("hi", "char", "word")} | \
        {("fa", f, p) for f in ("tibba", "bela", "pattan", "kacha") for p in ("word", "prefix", "suffix")}
    gids, cnt = np.unique(t["P_gid1"].to_numpy(), return_counts=True)
    keep = gids >= 0
    gids, cnt = gids[keep], cnt[keep]
    g = gn.drop_duplicates("gid").set_index("gid").reindex(gids)
    names = g["asciiname"].fillna("").str.lower().to_numpy()
    ccs = g["cc"].fillna("").to_numpy()
    hits: dict = {}
    form_cc: dict = {}
    for gid, name, cc, n in zip(gids, names, ccs, cnt, strict=True):
        toks = [x for x in __import__("re").split(r"[^a-z]+", name) if x]
        for e in lex.entries_for(cc):
            if any(_hit(tok, e) for tok in toks):
                key = (e.lang, e.form, e.pos)
                fc = form_cc.setdefault((cc, e.form, e.cat), [0, 0])
                fc[0] += 1
                fc[1] += int(n)
                if key in watch:
                    h = hits.setdefault(key, {"settlements": 0, "cells": 0, "by_country": {}, "examples": {}})
                    h["settlements"] += 1
                    h["cells"] += int(n)
                    h["by_country"][cc] = h["by_country"].get(cc, 0) + 1
                    h["examples"][name] = h["examples"].get(name, 0) + int(n)
    out = {f"{a}:{b}:{c}": {**v, "examples": dict(sorted(v["examples"].items(), key=lambda kv: -kv[1])[:10])}
           for (a, b, c), v in sorted(hits.items())}
    fc = pd.DataFrame([(cc, f, cat, s, c) for (cc, f, cat), (s, c) in form_cc.items()],
                      columns=["cc", "form", "cat", "settlements", "cells"])
    top_cc = fc.groupby("cc")["cells"].sum().sort_values(ascending=False).head(15).index
    by_cc = {cc: fc.loc[fc["cc"] == cc].sort_values("cells", ascending=False).head(10)[
        ["form", "cat", "settlements", "cells"]].to_dict("records") for cc in top_cc}
    return {"watched_forms": out, "top_forms_by_country": by_cc, "n_settlements": len(gids)}


# %% B12: 水名の読み
def hydronym_table(t: pd.DataFrame) -> pd.DataFrame:
    cache = VER / "v57_hydronym.parquet"
    if cache.is_file():
        return pd.read_parquet(cache)
    allc = pd.read_parquet(DATA_DIR / "global" / "geonames" / "allCountries.parquet",
                           columns=["gid", "asciiname", "lat", "lon", "fclass", "fcode", "cc"])
    h_all = allc.loc[allc["fclass"] == "H"]
    generic = frequent_tokens(h_all["asciiname"].fillna("").to_numpy())
    h = h_all.loc[h_all["fcode"].fillna("").str.startswith(HYDRONYM_CODES)].reset_index(drop=True)
    gids = np.unique(t["P_gid1"].to_numpy())
    gids = gids[gids >= 0]
    p = allc.loc[allc["fclass"] == "P"]
    s = p.drop_duplicates("gid").set_index("gid").reindex(gids)
    local = {cc: frequent_tokens(x["asciiname"].fillna("").to_numpy())
             for cc, x in p.loc[p["cc"].isin(s["cc"].dropna().unique())].groupby("cc")}
    tree = cKDTree(_unit_xyz(h["lon"].to_numpy(np.float64), h["lat"].to_numpy(np.float64)))
    xyz = _unit_xyz(s["lon"].to_numpy(np.float64), s["lat"].to_numpy(np.float64))
    near = tree.query_ball_point(xyz, r=2 * np.sin(HYDRONYM_RADIUS_M / EARTH_M / 2))
    hn = h["asciiname"].fillna("").to_numpy()
    rows = []
    for gid, name, cc, idx in zip(gids, s["asciiname"].fillna("").to_numpy(), s["cc"].fillna("").to_numpy(), near,
                                  strict=True):
        loc = local.get(cc, set())
        cores: set[str] = set()
        for i in idx:
            cores |= hydronym_cores(hn[i], generic, loc)
        rows.append((gid, int(contains_hydronym(name, cores)), len(idx)))
    out = pd.DataFrame(rows, columns=["gid", "hydronym", "n_water_features_5km"])
    out.to_parquet(cache, index=False)
    return out


def b12_hydronym(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    ht = hydronym_table(t).set_index("gid")
    t["lex_hydronym_n1"] = ht["hydronym"].reindex(t["P_gid1"].to_numpy()).fillna(0).to_numpy(np.float64)
    gids = ht.index.to_numpy()
    lx = settlement_lex(gids, gn, Lexicon(FROZEN_LEX))
    both = pd.DataFrame({"hyd": ht["hydronym"], "wet": lx["wetside"].reindex(gids)})
    raw = _json(DATA_DIR / "verification" / "v43_wikidata_raw.json")
    cats = raw["target_categories"]
    hydro = {"water", "wetland", "island", "crossing"}
    water_ns = [int(g) for g, qs in raw["named_after"].items() if any(set(cats.get(q, [])) & hydro for q in qs)]
    water_ns = [g for g in water_ns if g in ht.index]
    v44 = v55.v44
    samp = _json(v44.SAMPLE)
    sizes = samp["stratum_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    codes = pd.read_csv(v44.CODES)
    codes["gid"] = codes["id"].map(key)
    d = codes.join(meta, on="gid")
    d = d[d["status"] == "coded"].copy()
    d["truth"] = d["cats"].fillna("").map(lambda s: bool({c for c in s.split(";") if c} & set(WET_SIDE)))
    d["hydronym"] = ht["hydronym"].reindex(d["gid"].to_numpy()).fillna(0).to_numpy() > 0
    hm = d["hydronym"].to_numpy()
    out: dict = {"n_settlements": len(ht), "share_hydronym": float(ht["hydronym"].mean()),
                 "share_cells": float(t["lex_hydronym_n1"].mean()),
                 "overlap": {"hydronym_and_wetside": int(((both["hyd"] > 0) & (both["wet"] > 0)).sum()),
                             "hydronym_only": int(((both["hyd"] > 0) & (both["wet"] == 0)).sum()),
                             "wetside_only": int(((both["hyd"] == 0) & (both["wet"] > 0)).sum())},
                 "recall_wikidata_water_namesakes": {"k": int(ht["hydronym"].reindex(water_ns).sum()),
                                                     "n": len(water_ns)}}
    out["wikipedia_sample"] = {"weighted": v44.weighted(d, sizes, "hydronym"),
                               "n_hydronym_coded": int(hm.sum()),
                               "k_water_related_among_hydronym": int(d.loc[hm, "truth"].sum())}
    out["near_500m"] = v55.odds(t, cls="hydronym")
    sp = v55.settlement_points()
    sp["hydronym"] = ht["hydronym"].reindex(sp["P_gid1"].to_numpy()).fillna(0).to_numpy(np.float64)
    out["positive_control_hand"] = v55.fe_contrast(sp, "hand_mean", "hydronym", ())
    return out


# %% main
def _save(res: dict) -> None:
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")


def main() -> None:
    t0 = time.time()
    res: dict = _json(OUT) if OUT.is_file() and "--resume" in sys.argv else {}
    res |= {"rule": "amendment B-8 part B (local)", "seed": SEED}
    steps = [("B2", b2_achieved), ("B3", b3_needed), ("B4", b4_perm)]
    for k, f in steps:
        if k not in res:
            res[k] = f()
    _save(res)
    t, gn, infos, dirs = analysis_table()
    later = [("B9", lambda: b9_skill(t)), ("B8", b8_calibration), ("B11", lambda: b11_audit(t, gn)),
             ("B12", lambda: b12_hydronym(t, gn)), ("B7", lambda: b7_rivers(t, gn, dirs)), ("B5", lambda: b5_gfm(t)),
             ("B6", lambda: b6_variants(t, gn, infos, dirs))]
    for k, f in later:
        if k in res:
            continue
        res[k] = f()
        res["sec"] = round(time.time() - t0)
        _save(res)
        print(k, "done", res["sec"], flush=True)


if __name__ == "__main__":
    main()
