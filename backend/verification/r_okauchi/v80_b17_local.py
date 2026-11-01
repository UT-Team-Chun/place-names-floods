# %% v80（事後・amendment B-17 の B1〜B15）: 初見の査読 第 13 回に応じた点検（手元）
"""確証判定は変えない。B-17 で計算の前に固定した量を、保存したシミュレーション・二段目の保存・held-out の予測から計算する。

- B1: 被覆の点検。プロトコルの帰無と余白のシミュレーション（各 100 回）で、プロトコルの区間・z の Fieller・B-16 の B3 の
  少数クラスタの規則が真の値を覆う割合
- B2: 余白を全標本の実現値に固定した規則の動作特性
- B3: 無作為化推論（置換の分布に対する実の名の位置。クラスタの実行を引き上げた後）
- B4: ラベルの別定義で二段目を当て直す（B-16 の B1 の保存）
- B5: GFM と EMS の同じ比較（予測を各 AOI の浸水の割合に移す）
- B6: GFM の利得と集落の属性（C7 の保存の後）
- B7: 文字モデルの n-gram
- B8: 言語ごとの係数
- B9: 希釈の分解
- B10: 注入の判定を判定表の順で
- B11: 名の特徴の相関と分散拡大係数
- B12: 由来の標本の流れと言語の被覆の定義
- B13: 読める名の部分集合
- B14: 目録を改訂 16 と 17 に
- B15: 多重性

出力: stdout と docs/verification/v80_b17_local.json（部分ごとに足す。`--only B1,B2` で一部だけ）。
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
from scipy import stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b13 import SESOI_R

OUT = ROOT / "docs" / "verification" / "v80_b17_local.json"
VER = DATA_DIR / "verification"
SEED = 20261106
Z90 = float(stats.norm.ppf(0.95))
B_COVER = 1000
B_FIXED = 4000
SIMS = {"null": "power_v3_null", "margin": "power_v3_sesoi"}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


_MODS: dict = {}


def v77():
    if "v77" not in _MODS:
        _MODS["v77"] = _load("v77_b16_local")
    return _MODS["v77"]


def _wilson(k: int, n: int) -> list[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    z = 1.959964
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [mid - half, mid + half]


def _unit_sums(y: np.ndarray, p: dict, unit: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """場所グループごとの行数 n、名の増分 e = Σ[(M1c − y)² − (M3n − y)²]、物差し d = Σ[(M0 − y)² − (M1 − y)²]。"""
    _, inv = np.unique(unit, return_inverse=True)
    g = int(inv.max()) + 1
    yy = y.astype(np.float64)

    def sq(a: str) -> np.ndarray:
        return (p[a].astype(np.float64) - yy) ** 2

    n = np.bincount(inv, minlength=g).astype(np.float64)
    e = np.bincount(inv, weights=sq("M1c") - sq("M3n"), minlength=g)
    d = np.bincount(inv, weights=sq("M0") - sq("M1"), minlength=g)
    return n, e, d


# %% B1: 被覆の点検
def _cover_rep(n: np.ndarray, e: np.ndarray, d: np.ndarray, rng: np.random.Generator) -> dict:
    """一回のシミュレーションの点の推定と、場所グループの bootstrap（モデル固定）からの区間と上限。"""
    big_n = n.sum()
    e0, d0 = float(e.sum() / big_n), float(d.sum() / big_n)
    g = len(n)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_COVER).astype(np.float64)
    nw = w @ n
    eb, db = (w @ e) / nw, (w @ d) / nw
    var = v77()._cov3(eb, db)
    out = {"R": e0 / d0, "delta": e0, "yard": d0}
    lo, hi = np.quantile(eb, [0.05, 0.95])
    out["protocol90"] = [float(lo / d0), float(hi / d0)]  # Δ の百分位を全標本の D̂ で割る（プロトコルの区間）
    rule = v77()._nu_rule(e0, d0, var, v77()._nu_same(e, d, n), g)
    out["z_fieller90"] = rule["z"]["fieller90"]
    out["z_upper"] = rule["z"]["upper_R"]
    out["t_upper"] = rule["t_nu_minus_1"]["upper_R"]
    out["t_reason"] = rule["t_nu_minus_1"]["reason"]
    out["nu"] = rule.get("nu_l_plus_at_upper")
    # Δ だけ（Brier の単位）: 百分位の 90% 区間と、t(G* − 1) の片側 95% 上限（G* = θ = 0 の ν）
    gstar = rule["G_star_nu_delta"]
    out["delta_pct90"] = [float(lo), float(hi)]
    sd = float(np.std(eb, ddof=1))
    out["delta_t_upper"] = (e0 + float(stats.t.ppf(0.95, gstar - 1)) * sd) if gstar - 1 >= 1 else None
    out["G_star"] = gstar
    return out


def b1_coverage() -> dict:
    out: dict = {}
    for key, tag in SIMS.items():
        files = sorted((RUNS_DIR / tag / "sim").glob("rep*.npz"))
        reps = []
        t0 = time.time()
        for i, f in enumerate(files):
            with np.load(f, allow_pickle=False) as z:
                p = {a: z[a] for a in ("M0", "M1", "M1c", "M3n")}
                n, e, d = _unit_sums(z["y"], p, z["unit"])
            reps.append(_cover_rep(n, e, d, np.random.default_rng(SEED + i)))
        r_true = float(np.mean([x["R"] for x in reps]))
        d_true = float(np.mean([x["delta"] for x in reps]))
        m = len(reps)

        def share(cond: list[bool], m: int = m) -> dict:
            k = int(sum(cond))
            return {"k": k, "n": m, "share": k / m, "wilson95": _wilson(k, m)}

        res = {"n_reps": m, "R_true_mean": r_true, "delta_true_mean": d_true,
               "R_sd_between": float(np.std([x["R"] for x in reps], ddof=1)),
               "protocol90_two_sided": share([x["protocol90"][0] <= r_true <= x["protocol90"][1] for x in reps]),
               "protocol90_upper": share([x["protocol90"][1] >= r_true for x in reps]),
               "z_fieller90_two_sided": share([bool(x["z_fieller90"]) and x["z_fieller90"][0] <= r_true
                                               <= x["z_fieller90"][1] for x in reps]),
               "z_upper": share([x["z_upper"] is None or x["z_upper"] >= r_true for x in reps]),
               "t_upper": share([x["t_upper"] is None or x["t_upper"] >= r_true for x in reps]),
               "t_no_bound": share([x["t_upper"] is None for x in reps]),
               "delta_pct90_two_sided": share([x["delta_pct90"][0] <= d_true <= x["delta_pct90"][1] for x in reps]),
               "delta_t_upper": share([x["delta_t_upper"] is None or x["delta_t_upper"] >= d_true for x in reps]),
               "delta_t_no_bound": share([x["delta_t_upper"] is None for x in reps]),
               "nu_quantiles": np.quantile([x["nu"] for x in reps if x["nu"] is not None and math.isfinite(x["nu"])],
                                           [0.1, 0.5, 0.9]).tolist(),
               "G_star_quantiles": np.quantile([x["G_star"] for x in reps], [0.1, 0.5, 0.9]).tolist(),
               "seconds": round(time.time() - t0)}
        out[key] = res | {"reps": reps}
        print(f"  B1 {key}: {m} reps, R_true {r_true:.4f}; protocol two-sided {res['protocol90_two_sided']['share']:.2f}"
              f", z upper {res['z_upper']['share']:.2f}, t upper {res['t_upper']['share']:.2f}", flush=True)
    # 報告の規則（B-17）: 少数クラスタの規則の片側 95% の上限が二つのシナリオの両方で 90% 以上の回で真の値を上回るか
    out["rule_checked"] = all(out[k]["t_upper"]["share"] >= 0.90 for k in SIMS)
    out["rule"] = ("B-17 B1: 1,000 location-group bootstrap draws per replicate (models fixed, seed 20261106 + replicate); "
                   "true value = mean point estimate over the 100 replicates; 'no bound' counts as covering")
    return out


# %% B2: 余白を実現値に固定した動作特性
R_GRID2 = (-0.05, 0.0, 0.05)


def b2_fixed_margin() -> dict:
    """v69 の b1_operating と同じ外 2,000 × 内 1,000 の bootstrap で、規則の余白を全標本の 0.05 D̂ に固定した版。"""
    v69 = v77().v69
    d = v77().v36.load("world_v2", "primary")
    n, s, units = v69._group_sums(d, ("M0", "M1", "M1c", "M3n"))
    g = len(units)
    dg, eg = s["M0"] - s["M1"], s["M1c"] - s["M3n"]
    d_hat, e_hat = float(dg.sum() / n.sum()), float(eg.sum() / n.sum())
    rng = np.random.default_rng(SEED)
    rules = ("realised_margin", "outer_margin")
    cnt = {sh: {r0: dict.fromkeys(rules, 0) | {"undefined": 0} for r0 in R_GRID2} for sh in ("proportional", "constant")}
    t0 = time.time()
    for _ in range(v69.B_OUT):
        w = rng.multinomial(g, np.full(g, 1.0 / g)).astype(np.float64)
        nw = w @ n
        d_b = (w @ dg) / nw
        v = rng.multinomial(g, w / g, size=v69.B_IN).astype(np.float64)
        nv = v @ n
        dv, ev = (v @ dg) / nv, (v @ eg) / nv
        d_lo = float(np.quantile(dv, 0.025))
        for r0 in R_GRID2:
            for sh in ("proportional", "constant"):
                e_in = ev - e_hat + r0 * (dv if sh == "proportional" else d_hat)
                c = cnt[sh][r0]
                if d_lo <= 0:
                    c["undefined"] += 1
                    continue
                q = np.quantile(e_in, [0.05, 0.95])
                lo, hi = q / d_hat  # 余白を全標本の実現値に固定
                if -SESOI_R < lo and hi < SESOI_R:
                    c["realised_margin"] += 1
                lo2, hi2 = q / d_b  # v69 と同じ（外の標本の D で割る）
                if -SESOI_R < lo2 and hi2 < SESOI_R:
                    c["outer_margin"] += 1
    b = v69.B_OUT
    res = {sh: {f"{r0:+.2f}": {k: {"share": x / b, "mc_se": math.sqrt((x / b) * (1 - x / b) / b)} for k, x in c.items()}
                for r0, c in cs.items()} for sh, cs in cnt.items()}
    return {"n_units": g, "outer": b, "inner": v69.B_IN, "seed": SEED, "yardstick_hat": d_hat, "delta_hat": e_hat,
            "margin_brier": SESOI_R * d_hat, "share_equivalent": res, "seconds": round(time.time() - t0),
            "rule": "B-17 B2: v69 b1_operating with the margin fixed at 0.05 times the full-sample yardstick "
                    "(realised_margin) and, for comparison, divided by the outer resample's yardstick (outer_margin)"}


# %% 共通: B-16 の B1 の保存（主分割、凍結した読み）と、セルの割合
LABELS = DATA_DIR / "global" / "labels"
FRAC_COLS = ("f_riverine", "f_flood_any", "f_other", "f_trace", "f_perm_water", "f_observed")
_FR: dict = {}


def _fractions(aoi: np.ndarray, ix: np.ndarray, iy: np.ndarray) -> pd.DataFrame:
    """セルを AOI と格子の位置で cells.parquet に結び、浸水・痕跡・恒久水・観測の割合を返す（行の順は入力と同じ）。"""
    key = pd.DataFrame({"aoi": np.asarray(aoi).astype(str), "ix": np.asarray(ix, dtype=np.int64),
                        "iy": np.asarray(iy, dtype=np.int64)})
    key["row"] = np.arange(len(key))
    parts = []
    for a, sub in key.groupby("aoi", sort=False):
        if a not in _FR:
            c = pd.read_parquet(LABELS / a / "cells.parquet", columns=["ix", "iy", *FRAC_COLS])
            _FR[a] = c.astype({"ix": np.int64, "iy": np.int64})
        parts.append(sub.merge(_FR[a], on=["ix", "iy"], how="left", validate="many_to_one"))
    out = pd.concat(parts, ignore_index=True).sort_values("row").reset_index(drop=True)
    if out["f_riverine"].isna().any():
        msg = f"{int(out['f_riverine'].isna().sum())} cells not joined by AOI and grid position"
        raise SystemExit(msg)
    return out


_GSW: dict = {}


def _gsw_occ(aoi: np.ndarray, ix: np.ndarray, iy: np.ndarray) -> np.ndarray:
    if "t" not in _GSW:
        g = pd.read_parquet(DATA_DIR / "global" / "gsw" / "cell_gsw.parquet", columns=["aoi", "ix", "iy", "occ"])
        _GSW["t"] = g.astype({"aoi": str, "ix": np.int64, "iy": np.int64})
    key = pd.DataFrame({"aoi": np.asarray(aoi).astype(str), "ix": np.asarray(ix, dtype=np.int64),
                        "iy": np.asarray(iy, dtype=np.int64)})
    return key.merge(_GSW["t"], on=["aoi", "ix", "iy"], how="left", validate="many_to_one")["occ"].to_numpy(np.float64)


VARIANTS = ("t010", "t025", "notrace", "nonwater", "landobs", "gswmask", "fraction")


def _variant_labels(fr: pd.DataFrame, untyped: np.ndarray, occ: np.ndarray, v: str) -> tuple[np.ndarray, np.ndarray]:
    """B-17 の B4 のラベル（y、残すセル）。種別の無い内陸の製品は確証と同じく未記載の浸水を河川に数える。"""
    riv = fr["f_riverine"].to_numpy(np.float64) + np.where(untyped, fr["f_other"].to_numpy(np.float64), 0.0)
    anyf = fr["f_flood_any"].to_numpy(np.float64)
    th = 0.5
    keep = np.ones(len(fr), dtype=bool)
    if v in ("t010", "t025"):
        th = {"t010": 0.10, "t025": 0.25}[v]
    elif v == "notrace":
        tr = fr["f_trace"].to_numpy(np.float64)
        riv, anyf = np.clip(riv - tr, 0.0, None), np.clip(anyf - tr, 0.0, None)
    elif v == "nonwater":
        land = np.clip(1.0 - fr["f_perm_water"].to_numpy(np.float64), 1e-6, None)
        riv, anyf = riv / land, anyf / land
    elif v in ("landobs", "fraction"):
        lo = fr["f_observed"].to_numpy(np.float64) - fr["f_perm_water"].to_numpy(np.float64)
        keep = lo > 1e-6
        land = np.clip(lo, 1e-6, None)
        riv, anyf = riv / land, anyf / land
        if v == "fraction":
            return np.clip(riv, 0.0, 1.0), keep
    elif v == "gswmask":
        keep = ~(np.nan_to_num(occ, nan=0.0) >= 50.0)
    y1, y0 = riv >= th, anyf < th
    return y1.astype(np.float64), keep & (y1 | y0)


def _fit_stage2(f: dict, x: dict, arm: str, y_tr: np.ndarray, keep_tr: np.ndarray, *, frac: bool) -> np.ndarray:
    xtr, xte = x[arm]
    if not frac:
        p, _ = v77().stage2_weighted(f["z_tr"][keep_tr], f["z_te"], xtr[keep_tr], xte, y_tr[keep_tr].astype(np.int8))
        return p
    z2 = np.concatenate([f["z_tr"][keep_tr], f["z_tr"][keep_tr]])
    x2 = np.concatenate([xtr[keep_tr], xtr[keep_tr]])
    y2 = np.concatenate([np.ones(int(keep_tr.sum())), np.zeros(int(keep_tr.sum()))]).astype(np.int8)
    w2 = np.concatenate([y_tr[keep_tr], 1.0 - y_tr[keep_tr]])
    p, _ = v77().stage2_weighted(z2, f["z_te"], x2, xte, y2, w_tr=w2)
    return p


def _auc_change(y: np.ndarray, p_ref: np.ndarray, p_new: np.ndarray, aoi: np.ndarray) -> dict:
    from sklearn.metrics import roc_auc_score

    ch = []
    for a in np.unique(aoi):
        m = aoi == a
        if 0 < y[m].sum() < m.sum():
            ch.append(roc_auc_score(y[m], p_new[m]) - roc_auc_score(y[m], p_ref[m]))
    ch = np.array(ch)
    return {"n_aoi": len(ch), "median_change": float(np.median(ch)) if len(ch) else None,
            "share_higher": float((ch > 0).mean()) if len(ch) else None}


def _eq2_or(t: pd.DataFrame) -> dict:
    v39 = v77().v39
    rate = t.groupby("aoi")["y"].transform("mean")
    t = t[(rate > 0) & (rate < 1)].reset_index(drop=True)
    z, names = v39.design(t, ("wetside",), True, banded=True)
    j = names.index("wetside:le500m")
    try:
        r = v39.fit(t["y"].to_numpy(np.float64), z, t["aoi"].to_numpy(), t["unit"].to_numpy(str))
    except np.linalg.LinAlgError as e:
        return {"error": str(e)}
    b, se = float(r["beta"][j]), float(r["se"][j])
    near = (t["P_d1"].to_numpy() <= 500.0) & (t["lex_wetside_n1"].to_numpy() > 0) & (t["y"].to_numpy() == 1)
    return {"odds_ratio": math.exp(b), "ci95": [math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)],
            "n_cells": len(t), "n_flooded_near": int(near.sum())}


def _r_from_cells(y: np.ndarray, p: dict, unit: np.ndarray, rng: np.random.Generator,
                  comps: tuple[tuple[str, str], ...] = (("M1c", "M3n"), ("M1s", "M3d"))) -> dict:
    """モデル固定の Δ・D・R・90% Fieller（場所グループの 4,000 回）。"""
    u, inv = np.unique(unit, return_inverse=True)
    g = len(u)
    n = np.bincount(inv, minlength=g).astype(np.float64)
    sq = {a: np.bincount(inv, weights=(p[a] - y) ** 2, minlength=g) for a in p}
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)
    nw = w @ n
    dg = sq["M0"] - sq["M1"]
    out = {"n_cells": len(y), "n_units": g, "mean_label": float(y.mean()), "yardstick": float(dg.sum() / n.sum())}
    for ref, new in comps:
        eg = sq[ref] - sq[new]
        e0, d0 = float(eg.sum() / n.sum()), float(dg.sum() / n.sum())
        det = v77().v73._levels_detail(e0, d0, (w @ eg) / nw, (w @ dg) / nw, rng)
        out[f"{ref}_{new}"] = {"delta_bs": e0, "R": e0 / d0, "detail": det}
    return out


FR_CACHE = VER / "v80_fractions_fold{k}.npz"


def _b4_prep() -> None:
    """保存の各 fold の学習・評価のセルの割合・種別の無い内陸の印・GSW の頻度をキャッシュに。"""
    vv = v77()
    folds = vv._dump_folds()
    untyped_map = vv._superclusters()["aoi_untyped"]
    for k, f in enumerate(folds):
        path = Path(str(FR_CACHE).format(k=k))
        if path.is_file():
            continue
        arrs = {}
        for s in ("tr", "te"):
            fr = _fractions(f[f"aoi_{s}"], f[f"ix_{s}"], f[f"iy_{s}"])
            for c in FRAC_COLS:
                arrs[f"{c}_{s}"] = fr[c].to_numpy(np.float32)
            arrs[f"unt_{s}"] = np.array([untyped_map.get(a, False) for a in f[f"aoi_{s}"].astype(str)])
            arrs[f"occ_{s}"] = _gsw_occ(f[f"aoi_{s}"], f[f"ix_{s}"], f[f"iy_{s}"]).astype(np.float32)
        np.savez_compressed(path, **arrs)


def _sides(k: int) -> dict:
    with np.load(Path(str(FR_CACHE).format(k=k)), allow_pickle=False) as z:
        a = {kk: z[kk] for kk in z.files}
    return {s: {"fr": pd.DataFrame({c: a[f"{c}_{s}"].astype(np.float64) for c in FRAC_COLS}), "unt": a[f"unt_{s}"],
                "occ": a[f"occ_{s}"].astype(np.float64)} for s in ("tr", "te")}


def _b4_one(v: str) -> dict:
    vv = v77()
    v77()._threads1()
    folds = vv._dump_folds()
    names = [str(c) for c in folds[0]["name_cols"]]
    jw = names.index("lex_wetside_n1")
    prep = [(f, vv.v76._designs(f, frozen=True), _sides(k)) for k, f in enumerate(folds)]
    t0 = time.time()
    frac = v == "fraction"
    cols = {k: [] for k in ("y", "unit", "aoi", "d1", "wet", "M0", "M1", "M1c", "M3n", "M1s", "M3d")}
    for f, x, sides in prep:
        y_tr, k_tr = _variant_labels(sides["tr"]["fr"], sides["tr"]["unt"], sides["tr"]["occ"], v)
        y_te, k_te = _variant_labels(sides["te"]["fr"], sides["te"]["unt"], sides["te"]["occ"], v)
        for arm in ("M1c", "M3n", "M1s", "M3d"):
            cols[arm].append(_fit_stage2(f, x, arm, y_tr, k_tr, frac=frac)[k_te])
        cols["y"].append(y_te[k_te])
        cols["unit"].append(f["unit_te"][k_te])
        cols["aoi"].append(f["aoi_te"][k_te])
        cols["d1"].append(f["d1_te"][k_te])
        xn = x["M3n"][1][:, x["M1c"][1].shape[1]:]
        cols["wet"].append(xn[k_te, jw])
        cols["M0"].append(f["p_M0_te"].astype(np.float64)[k_te])
        cols["M1"].append(f["p_M1_te"].astype(np.float64)[k_te])
    c = {k: np.concatenate(vals) for k, vals in cols.items()}
    rng = np.random.default_rng(SEED)
    res = _r_from_cells(c["y"], {a: c[a] for a in ("M0", "M1", "M1c", "M3n", "M1s", "M3d")}, c["unit"], rng)
    if not frac:
        res["auc_within_aoi"] = _auc_change(c["y"], c["M1c"], c["M3n"], c["aoi"])
        t = pd.DataFrame({"y": c["y"], "aoi": c["aoi"], "unit": c["unit"], "P_d1": c["d1"],
                          "lex_wetside_n1": c["wet"], "p_M1": c["M1"]})
        res["eq2_capped_ems"] = _eq2_or(t)
    print(f"  B4 {v}: R {res['M1c_M3n']['R']:.4f} centred {res['M1s_M3d']['R']:.4f} "
          f"({time.time() - t0:.0f} s)", flush=True)
    return res


def b4_label_variants() -> dict:
    """B-16 の B1 の保存（凍結した読み）で、ラベルの別定義に二段目を当て直す（基礎の M0・M1 は固定）。"""
    from concurrent.futures import ProcessPoolExecutor

    t0 = time.time()
    _b4_prep()
    print(f"  B4: fractions cached ({time.time() - t0:.0f} s)", flush=True)
    with ProcessPoolExecutor(min(len(VARIANTS), 7)) as ex:
        res = dict(zip(VARIANTS, ex.map(_b4_one, VARIANTS), strict=True))
    out: dict = dict(res)
    out["rule"] = ("B-17 B4: dump of B1 of B-16 (primary split, frozen reading); fractions joined by AOI and grid "
                   "position; second stage (M1c, M3n, M1s, M3d) refitted on each label definition with base M0 and M1 "
                   "fixed; the fraction target is fitted by a binomial likelihood weighted by the fraction; R with the "
                   "definition's own yardstick BS(M0) - BS(M1); 4,000 draws of location groups (seed 20261106)")
    return out


# %% B11: 名の特徴の相関と分散拡大係数
def b11_correlation() -> dict:
    vv = v77()
    folds = vv._dump_folds()
    xs = []
    for f in folds:
        x = vv.v76._designs(f, frozen=True)
        xs.append(x["M3n"][1][:, x["M1c"][1].shape[1]:])
    xn = np.concatenate(xs)
    names = [str(c) for c in folds[0]["name_cols"]]
    sd = xn.std(axis=0)
    use = [i for i in range(xn.shape[1]) if sd[i] > 0]
    r = np.corrcoef(xn[:, use], rowvar=False)
    vif = np.diag(np.linalg.pinv(r))
    pairs = sorted(((names[use[i]], names[use[j]], float(r[i, j])) for i in range(len(use))
                    for j in range(i + 1, len(use))), key=lambda z: -abs(z[2]))
    return {"n_cells": len(xn), "features": [names[i] for i in use], "constant": [names[i] for i in range(len(names))
                                                                                    if i not in use],
            "vif": {names[use[i]]: float(vif[i]) for i in range(len(use))}, "top_pairs": pairs[:12],
            "rule": "B-17 B11: all held-out cells of the B1 dump (frozen reading); Pearson correlations and VIF = "
                    "diagonal of the inverse correlation matrix"}


# %% B13: 読める名の部分集合
def b13_read_subset() -> dict:
    vv = v77()
    pt = dict(np.load(vv.POINT_CACHE))
    folds = vv._dump_folds()
    names = [str(c) for c in folds[0]["name_cols"]]
    cls = [names.index(f"lex_{c}_n1") for c in ("water", "wet", "riverine", "crossing", "low", "high")]
    xs = []
    for f in folds:
        x = vv.v76._designs(f, frozen=True)
        xs.append(x["M3n"][1][:, x["M1c"][1].shape[1]:][:, cls])
    read = (np.concatenate(xs) > 0).any(axis=1)
    if len(read) != len(pt["y"]):
        return {"error": "row counts differ between the dump and the point refit"}
    y = pt["y"].astype(np.float64)
    p = {a: pt[f"p_{a}"].astype(np.float64) for a in ("M0", "M1", "M1c", "M3n", "M1s", "M3d")}
    unit, aoi = pt["unit"].astype(str), pt["aoi"].astype(str)
    rng = np.random.default_rng(SEED)
    out = {"share_read": float(read.mean()),
           "read_only": _r_from_cells(y[read], {a: v[read] for a, v in p.items()}, unit[read], rng)}
    share = pd.Series(read.astype(float)).groupby(aoi).transform("mean").to_numpy()
    by_aoi = pd.Series(read.astype(float)).groupby(aoi).mean()
    q1, q2 = np.quantile(by_aoi.to_numpy(), [1 / 3, 2 / 3])
    out["aoi_share_tercile_cuts"] = [float(q1), float(q2)]
    for lab, m in (("low", share <= q1), ("mid", (share > q1) & (share <= q2)), ("high", share > q2)):
        out[f"aoi_tercile_{lab}"] = _r_from_cells(y[m], {a: v[m] for a, v in p.items()}, unit[m], rng)
    out["rule"] = ("B-17 B13: predictions of the point refit of B1 of B-16 (frozen reading, models fixed); 'read' = the "
                   "nearest name matches an element of any of the six classes; AOI terciles of the share of read "
                   "cells (AOIs as units, unweighted); R with the subset's own yardstick; 4,000 draws of location groups")
    return out


# %% B5: GFM と EMS の同じ比較（予測を各 AOI の浸水の割合に移す）
def _within(y: np.ndarray, aoi_raw: np.ndarray, unit: np.ndarray, preds: dict, rng: np.random.Generator) -> dict:
    """B-10 の B10（v66 の within_placebo）と同じ採点: 予測を AOI ごとに浸水率へ移し、R と置換との差の区間。"""
    if "v63" not in _MODS:
        _MODS["v63"] = _load("v63_b9_local")
    v63 = _MODS["v63"]
    aoi = pd.factorize(pd.Series(aoi_raw))[0]
    n_a = np.bincount(aoi).astype(np.float64)
    rate = np.bincount(aoi, weights=y) / n_a
    inner = (rate > 0) & (rate < 1)
    u, inv = np.unique(unit, return_inverse=True)
    g = len(u)
    cnt = rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)
    sse = {}
    for m, pv in preds.items():
        pp = np.clip(pv.astype(np.float64), 1e-7, 1 - 1e-7)
        q = v63._shift_to_rate(np.log(pp / (1 - pp)), aoi, rate, n_a)
        q = np.where(inner[aoi], q, rate[aoi])
        sse[m] = np.bincount(inv, weights=(q - y) ** 2, minlength=g)
    w = np.vstack([np.ones(g), cnt])
    den = w @ sse["M0"] - w @ sse["M1"]
    r = (w @ sse["M1c"] - w @ sse["M3n"]) / den
    rp = (w @ sse["M1c"] - w @ sse["M3np"]) / den
    dif = r - rp
    q2 = (0.025, 0.975)
    return {"n_cells": len(y), "n_aoi_inner": int(inner.sum()), "n_units": g, "R": float(r[0]),
            "ci95_R": np.quantile(r[1:], q2).tolist(), "ci90_R": np.quantile(r[1:], (0.05, 0.95)).tolist(),
            "R_placebo": float(rp[0]), "R_minus_placebo": float(dif[0]), "ci95_R_minus_placebo": np.quantile(dif[1:], q2).tolist()}


def _key(aoi: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> pd.MultiIndex:
    return pd.MultiIndex.from_arrays([np.asarray(aoi).astype(str),
                                      np.round(np.asarray(lon, dtype=np.float32).astype(np.float64), 6),
                                      np.round(np.asarray(lat, dtype=np.float32).astype(np.float64), 6)])


def b5_gfm_ems() -> dict:
    v36 = v77().v36
    e = v36.load("world_v2", "primary")
    gf = v36.load("world_v2_gfm", "primary")
    ke, kg = _key(e["aoi"], e["lon"], e["lat"]), _key(gf["aoi"], gf["lon"], gf["lat"])
    on_gfm = ke.isin(kg)
    arms = ("M0", "M1", "M1c", "M3n", "M3np")
    y = np.asarray(e["y"], dtype=np.float64)
    out = {"n_gfm_cells": len(kg), "n_joined": int(on_gfm.sum())}
    for lab, m in (("ems_on_gfm_cells", on_gfm), ("ems_all_aois", np.ones(len(y), dtype=bool))):
        out[lab] = _within(y[m], np.asarray(e["aoi"])[m], np.asarray(e["unit"]).astype(str)[m],
                           {a: np.asarray(e[f"p_{a}"])[m] for a in arms}, np.random.default_rng(SEED))
    out["rule"] = ("B-17 B5: predictions of the confirmatory run world_v2 (trained on EMS, frozen reading, models fixed) "
                   "scored on EMS labels after shifting each AOI's predictions on the logit scale to its flooded share "
                   "(v63 _shift_to_rate; all-flooded or all-dry AOIs scored at their rate); R with the shifted baseline "
                   "gain; real minus reassigned names (M3np); 4,000 draws of location groups (seed 20261106); the GFM "
                   "cells are those of world_v2_gfm, joined by AOI and coordinates")
    return out


# %% B6: GFM の利得と集落の属性（追記 3）
GFM_DUMP = "world_v2_b17gfmdump2"  # 追記 4: utens で回し直した保存（手元の途中の保存は使わない）
FCODE_ADMIN = {"PPLA", "PPLA2", "PPLA3", "PPLA4", "PPLC"}


def _lc_cells(aoi: np.ndarray, ix: np.ndarray, iy: np.ndarray) -> pd.DataFrame:
    """セルの WorldCover の市街地と樹木の割合（features.parquet、cells.parquet と同じ行の順）。"""
    key = pd.DataFrame({"aoi": np.asarray(aoi).astype(str), "ix": np.asarray(ix, dtype=np.int64),
                        "iy": np.asarray(iy, dtype=np.int64)})
    key["row"] = np.arange(len(key))
    parts = []
    for a, sub in key.groupby("aoi", sort=False):
        c = pd.read_parquet(LABELS / a / "cells.parquet", columns=["ix", "iy"]).astype({"ix": np.int64, "iy": np.int64})
        f = pd.read_parquet(LABELS / a / "features.parquet", columns=["lc_built", "lc_tree"])
        parts.append(sub.merge(pd.concat([c, f], axis=1), on=["ix", "iy"], how="left", validate="many_to_one"))
    return pd.concat(parts, ignore_index=True).sort_values("row").reset_index(drop=True)


def _settlement_attrs(f: dict, gn: pd.DataFrame, first: dict[str, str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """一つの fold の学習と評価のセルの集落の属性（人口・地物の類・500 m の市街地と樹木・言語の地域）。"""
    gid = np.concatenate([f["gid_tr"], f["gid_te"]]).astype(np.int64)
    d1 = np.concatenate([f["d1_tr"], f["d1_te"]]).astype(np.float64)
    aoi = np.concatenate([f["aoi_tr"], f["aoi_te"]]).astype(str)
    ix = np.concatenate([f["ix_tr"], f["ix_te"]])
    iy = np.concatenate([f["iy_tr"], f["iy_te"]])
    lc = _lc_cells(aoi, ix, iy)
    near = d1 <= 500.0
    by = pd.DataFrame({"gid": gid[near], "b": lc["lc_built"].to_numpy(np.float64)[near],
                       "t": lc["lc_tree"].to_numpy(np.float64)[near]}).groupby("gid").mean()
    g = gn.reindex(gid)
    pop = np.log10(1.0 + g["population"].fillna(0).clip(lower=0).to_numpy(np.float64))
    fc = g["fcode"].fillna("").astype(str).to_numpy()
    cls = np.column_stack([np.isin(fc, list(FCODE_ADMIN)), fc == "PPLX",
                           ~np.isin(fc, [*FCODE_ADMIN, "PPL", "PPLX"])]).astype(np.float64)
    bt = by.reindex(gid)
    miss = bt["b"].isna().to_numpy()
    built, tree = bt["b"].fillna(0).to_numpy(np.float64), bt["t"].fillna(0).to_numpy(np.float64)
    cc = np.concatenate([f["cc_tr"], f["cc_te"]]).astype(str)
    blk = np.array([first.get(c, "none") for c in cc])
    n_tr = len(f["gid_tr"])
    vals, cnts = np.unique(blk[:n_tr], return_counts=True)
    areas = sorted(v for v, n in zip(vals, cnts, strict=True) if n >= 0.02 * n_tr and v != "none")
    lang = np.column_stack([blk == a for a in areas]).astype(np.float64) if areas else np.zeros((len(blk), 0))
    x = np.column_stack([pop, cls, built, tree, miss.astype(np.float64), lang])
    names = ["log_pop", "fc_admin", "fc_section", "fc_other", "built500", "tree500", "no_cells500",
             *[f"lang_{a}" for a in areas]]
    return x[:n_tr], x[n_tr:], names


def _b6_one(k: int) -> tuple[int, dict]:
    vv = v77()
    vv._threads1()
    src = RUNS_DIR / GFM_DUMP / "predictions" / "primary" / f"fold{k}_stage2.npz"
    with np.load(src, allow_pickle=False) as z:
        f = {kk: z[kk] for kk in z.files}
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "fcode", "population"]).drop_duplicates("gid").set_index("gid")
    lex = vv.Lexicon(vv.FROZEN_LEX)
    first = {cc: langs[0] for cc, langs in lex.country_langs.items() if langs}
    a_tr, a_te, names = _settlement_attrs(f, gn, first)
    preds = {}
    for arm, extra_tr, extra_te in (("M1c", None, None), ("M3n", f["xn_tr"], f["xn_te"]),
                                    ("M3np", f["xnp_tr"], f["xnp_te"])):
        for lab, add_tr, add_te in (("", None, None), ("A", a_tr, a_te)):
            xs_tr = [f["xc_tr"]] + ([add_tr] if add_tr is not None else []) + ([extra_tr] if extra_tr is not None else [])
            xs_te = [f["xc_te"]] + ([add_te] if add_te is not None else []) + ([extra_te] if extra_te is not None else [])
            preds[arm + lab] = _stage2_c(f["z_tr"], f["z_te"], np.column_stack(xs_tr), np.column_stack(xs_te),
                                         f["y_tr"], 1.0).astype(np.float32)
    return k, {"preds": preds, "names": names, "n_te": len(f["y_te"])}


def b6_gfm_attrs() -> dict:
    from concurrent.futures import ProcessPoolExecutor

    src = RUNS_DIR / GFM_DUMP / "predictions" / "primary"
    if len(list(src.glob("fold*_stage2.npz"))) < 10:
        return {"note": f"{GFM_DUMP} incomplete"}
    with ProcessPoolExecutor(5) as ex:
        res = dict(ex.map(_b6_one, range(10)))
    parts = {kk: [] for kk in ("y", "aoi", "unit", "p_M0", "p_M1")}
    arms: dict[str, list] = {}
    for k in range(10):
        with np.load(src / f"fold{k}_stage2.npz", allow_pickle=False) as z:
            parts["y"].append(z["y_te"].astype(np.float64))
            parts["aoi"].append(z["aoi_te"].astype(str))
            parts["unit"].append(z["unit_te"].astype(str))
            parts["p_M0"].append(z["p_M0_te"])
            parts["p_M1"].append(z["p_M1_te"])
        for a, pv in res[k]["preds"].items():
            arms.setdefault(a, []).append(pv)
    y = np.concatenate(parts["y"])
    aoi, unit = np.concatenate(parts["aoi"]), np.concatenate(parts["unit"])
    base = {"M0": np.concatenate(parts["p_M0"]), "M1": np.concatenate(parts["p_M1"])}
    cat = {a: np.concatenate(v) for a, v in arms.items()}
    out: dict = {"n_cells": len(y), "attributes_fold0": res[0]["names"]}
    for lab, suf in (("without_attributes", ""), ("with_attributes", "A")):
        preds = base | {"M1c": cat["M1c" + suf], "M3n": cat["M3n" + suf], "M3np": cat["M3np" + suf]}
        w = _within(y, aoi, unit, preds, np.random.default_rng(SEED))
        bs = {a: float(np.mean((preds[a].astype(np.float64) - y) ** 2)) for a in ("M1c", "M3n", "M3np")}
        out[lab] = w | {"delta_bs": bs["M1c"] - bs["M3n"], "brier": bs}
    out |= {k: out["with_attributes"][k] for k in ("R", "ci95_R", "R_minus_placebo", "ci95_R_minus_placebo")}
    # 確かめ: 属性なしの当て直しが保存した二段目の予測（world_eval の M1c・M3n）を再現するか
    out["refit_vs_saved_max_abs"] = {a: float(np.max(np.abs(cat[a].astype(np.float64) - np.concatenate(
        [np.load(src / f"fold{k}_stage2.npz")[f"p_{a}_te"] for k in range(10)]).astype(np.float64))))
        for a in ("M1c", "M3n")}
    out["rule"] = ("B-17 B6 with addenda 3 and 4: dump world_v2_b17gfmdump2 on utens (GFM labels, corrected code, "
                   "placebo names "
                   "saved); second stage refitted (C = 1, standardised) for M1c, M3n and M3np with and without the "
                   "settlement attributes (log10(1 + population), feature class PPLA-PPLA4/PPLC, PPLX, other against PPL, "
                   "mean WorldCover built-up and tree shares of the settlement's cells within 500 m with a no-cell "
                   "indicator, first lexicon block of the country with at least 2% of the fold's training cells); "
                   "predictions shifted to each AOI's flooded share; 4,000 draws of location groups (seed 20261106)")
    return out


def b6_check_cluster() -> dict:
    """追記 3・4: クラスタの C7 の保存（world_v2_b17gfmdump）の実の配列が、B6 に使う保存と一致するか（共通の鍵の最大の差）。"""
    loc = RUNS_DIR / GFM_DUMP / "predictions" / "primary"  # B6 に使う保存（utens）
    clu = RUNS_DIR / "world_v2_b17gfmdump" / "predictions" / "primary"
    if len(list(clu.glob("fold*_stage2.npz"))) < 10:
        return {"note": "world_v2_b17gfmdump incomplete locally"}
    out: dict = {}
    for k in range(10):
        with np.load(loc / f"fold{k}_stage2.npz") as a, np.load(clu / f"fold{k}_stage2.npz") as b:
            common = sorted(set(a.files) & set(b.files))
            diff = {}
            for key in common:
                x, y = a[key], b[key]
                if x.shape != y.shape:
                    diff[key] = "shape differs"
                elif np.issubdtype(x.dtype, np.number):
                    diff[key] = float(np.max(np.abs(x.astype(np.float64) - y.astype(np.float64)))) if x.size else 0.0
                else:
                    diff[key] = "equal" if np.array_equal(x, y) else "differs"
            out[f"fold{k}"] = {"only_local": sorted(set(a.files) - set(b.files)), "max_abs_diff": diff}
    num = [v for f in out.values() for v in f["max_abs_diff"].values() if isinstance(v, float)]
    bad = [f"{k}:{kk}" for k, f in out.items() for kk, v in f["max_abs_diff"].items() if isinstance(v, str) and v != "equal"]
    return {"max_abs_diff_all": max(num) if num else None, "non_numeric_mismatch": bad, "folds": out,
            "rule": "B-17 addenda 3 and 4: real arrays of the cluster dump C7 against the dump used for B6"}


# %% B9: 希釈の分解
DIL = {"gsw22": ("world_v2_b13gsw22", ("M1cW", "M3oW"), ("M1cG", "M3oG")),
       "recur_prior": ("world_v2_b13recprior", ("M1cEc", "M3oEc"), ("M1cE", "M3oE")),
       "jrc": ("world_v2_b13jrcset", ("M1cJc", "M3oJc"), ("M1cJ", "M3oJ"))}


def b9_dilution() -> dict:
    vv = v77()
    out: dict = {}
    for key, (run, cell, setv) in DIL.items():
        d = vv.v36.load(run, "primary")
        y = np.asarray(d["y"], dtype=np.float64)
        j = vv._join_cells(np.asarray(d["aoi"]), np.asarray(d["lon"]), np.asarray(d["lat"]))
        near = j["P_d1"].to_numpy(np.float64) <= 500.0
        u, inv = np.unique(np.asarray(d["unit"]).astype(str), return_inverse=True)
        g = len(u)
        rng = np.random.default_rng(SEED)
        w = np.vstack([np.ones(g), rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)])

        def gain(pair: tuple[str, str], m: np.ndarray, d=d, y=y, inv=inv, g=g, w=w) -> np.ndarray:
            a, b = (np.asarray(d[f"p_{x}"], dtype=np.float64) for x in pair)
            e = np.where(m, (a - y) ** 2 - (b - y) ** 2, 0.0)
            n = np.bincount(inv, weights=m.astype(np.float64), minlength=g)
            return (w @ np.bincount(inv, weights=e, minlength=g)) / (w @ n)

        res = {"run": run, "n_cells": len(y), "share_near": float(near.mean())}
        for ground, m in (("all", np.ones(len(y), dtype=bool)), ("near500", near)):
            gc, gs = gain(cell, m), gain(setv, m)
            ok = np.isfinite(gc) & np.isfinite(gs)
            var = vv._cov3(gs[1:][ok[1:]], gc[1:][ok[1:]])
            from topsem.world.b13 import fieller_interval

            fi = fieller_interval(float(gs[0]), float(gc[0]), var, z=Z90)
            res[ground] = {"delta_cell": float(gc[0]), "delta_set": float(gs[0]), "share_kept": float(gs[0] / gc[0]),
                           "fieller90_share": list(fi) if fi else None}
        out[key] = res
    out["rule"] = ("B-17 B9: predictions of the record runs of amendment 13 (models fixed); the gain of each version over "
                   "its reference on all ground and on ground within 500 m of the nearest settlement (P_d1 joined by AOI "
                   "and coordinates); share kept = through-settlement gain / per-cell gain, with a 90 % Fieller interval "
                   "from 4,000 draws of location groups (seed 20261106)")
    return out


# %% B10: 注入の判定を判定表の順で
def b10_injected_verdicts() -> dict:
    from collections import Counter

    c4 = _json(ROOT / "docs" / "verification" / "v64_b9_aggregate.json")["C4"]
    out = {}
    for b in c4["by_delta"]:
        out[f"{b['odds_ratio']:.2f}"] = {"reps": b["reps"], "verdicts": dict(Counter(r["decision"] for r in b["rows"])),
                                         "ci95_above_zero": b["n_positive"], "tost_equivalent": b["n_equivalent"]}
    out["rule"] = ("B-17 B10: verdicts of the decision table (positive verdicts before equivalent-null) of each replicate of "
                   "C4 of B-9 (every model refitted), from v64_b9_aggregate.json")
    return out


# %% B12: 由来の標本の流れと言語の被覆の定義
def b12_flow() -> dict:
    out: dict = {}
    v44 = ROOT / "docs" / "verification" / "v44_wikipedia_etymology.json"
    if v44.is_file():
        w = _json(v44)
        out["v44_keys"] = sorted(w.keys())[:40]
        for k in ("n_sampled", "n_with_section", "n_coded", "n_derivation", "n_no_derivation", "n_disputed",
                  "round1", "round2", "flow"):
            if k in w:
                out[k] = w[k]
    codes = ROOT / "docs" / "verification" / "v44_etym_codes.csv"
    if codes.is_file():
        c = pd.read_csv(codes)
        out["v44_codes_rows"] = len(c)
        out["v44_codes_columns"] = list(c.columns)[:30]
    two = ROOT / "docs" / "verification" / "v70_two_coders.csv"
    if two.is_file():
        c2 = pd.read_csv(two)
        out["v70_two_coders_rows"] = len(c2)
    aud = ROOT / "docs" / "verification" / "v69_audit_sample.csv"
    if aud.is_file():
        out["v69_audit_sample_rows"] = len(pd.read_csv(aud))
    out["rule"] = "B-17 B12: counts read from the published records (v44, v70, v69); the flow is assembled in the paper's table"
    return out


# %% B15: 多重性
def b15_multiplicity() -> dict:
    s = _json(RUNS_DIR / "power_v3_null" / "power_summary.json")
    dc = s["decision_counts"]
    n = s["n_reps"]
    non_eq = n - dc.get("equivalent-null", 0)
    return {"null_simulation_decisions": dc, "share_not_equivalent": non_eq / n,
            "rule": "B-17 B15: share of null-simulation replicates whose verdict is not equivalent-null; the number of "
                    "pre-specified verdicts is counted in the paper's tables (Tables 2, C1 and C2)"}


# %% B7: 文字モデルの n-gram
def _ngram_table(names: np.ndarray, y: np.ndarray, k: int = 20) -> dict:
    from sklearn.utils import murmurhash3_32

    from topsem.world import channels as ch

    model = ch._fit_name_model(names, y)
    if isinstance(model, float):
        return {"note": "constant model"}
    analyzer = ch._VEC.build_analyzer()
    cands: dict[int, str] = {}
    count: dict[str, int] = {}
    for s in pd.unique(names):
        if not s:
            continue
        for g in set(analyzer(str(s))):
            cands.setdefault(abs(murmurhash3_32(g, seed=0)) % ch._VEC.n_features, g)
            count[g] = count.get(g, 0) + 1
    coef = model.coef_.ravel()
    seen = np.array(sorted(cands))
    order = seen[np.argsort(-coef[seen])]

    def rows(idx: np.ndarray) -> list[dict]:
        return [{"ngram": cands[int(j)], "weight": float(coef[j]), "n_names": count.get(cands[int(j)], 0)}
                for j in idx]
    return {"n_names": len(set(names) - {""}), "largest": rows(order[:k]), "smallest": rows(order[::-1][:k])}


def b7_ngrams() -> dict:
    import world_eval as ev

    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    look = dict(zip(gn["gid"].astype(np.int64), gn["asciiname"].astype(str), strict=True))
    t = pd.read_parquet(v77().v39.TABLE, columns=["aoi", "lon", "lat", "y", "P_gid1"])
    names = np.array([look.get(int(g), "") if g >= 0 else "" for g in t["P_gid1"].to_numpy()], dtype=object)
    out = {"ems": _ngram_table(names, t["y"].to_numpy(np.float64))}
    gf = v77().v36.load("world_v2_gfm", "primary")
    pos = pd.Series(np.arange(len(t)), index=_key(t["aoi"], t["lon"], t["lat"]))
    idx = pos.reindex(_key(gf["aoi"], gf["lon"], gf["lat"])).to_numpy()
    ok = ~np.isnan(idx)
    out["gfm"] = _ngram_table(names[idx[ok].astype(np.int64)], np.asarray(gf["y"], dtype=np.float64)[ok])
    out["gfm"]["n_cells_joined"] = int(ok.sum())
    out["rule"] = ("B-17 B7: the character model of the protocol (char_wb 1-4-grams hashed to 2^20, names aggregated to "
                   "their flooded share, L2 logistic regression with C = 1) fitted once on all analysed cells with EMS "
                   "labels and on the cells GFM can assess with GFM labels; hash buckets read back to the n-grams of the "
                   "analysed names (collisions keep the first n-gram met)")
    return out


# %% B8: 言語ごとの係数
def _area_codes(cc: np.ndarray, areas: list[str], first_block: dict[str, str]) -> np.ndarray:
    blk = np.array([first_block.get(str(c), "none") for c in cc])
    return np.array([areas.index(b) if b in areas else len(areas) for b in blk])


def _stage2_c(z_tr, z_te, x_tr, x_te, y_tr, c: float) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression

    a = np.column_stack([z_tr, x_tr]).astype(np.float64)
    b = np.column_stack([z_te, x_te]).astype(np.float64)
    mu, sd = a.mean(axis=0), a.std(axis=0)
    sd[sd == 0] = 1.0
    m = LogisticRegression(C=c, max_iter=2000)
    m.fit((a - mu) / sd, y_tr)
    return m.predict_proba((b - mu) / sd)[:, 1]


def _b8_one(args: tuple[int, float]) -> tuple[int, float, np.ndarray]:
    k, c = args
    vv = v77()
    vv._threads1()
    f = vv._dump_folds()[k]
    x = vv.v76._designs(f, frozen=True)
    lex = vv.Lexicon(vv.FROZEN_LEX)
    first = {cc: langs[0] for cc, langs in lex.country_langs.items() if langs}
    tr_blk = np.array([first.get(str(cc), "none") for cc in f["cc_tr"]])
    vals, cnts = np.unique(tr_blk, return_counts=True)
    areas = sorted(v for v, n in zip(vals, cnts, strict=True) if n >= 0.02 * len(tr_blk) and v != "none")
    k1 = x["M1c"][0].shape[1]
    xn_tr, xn_te = x["M3n"][0][:, k1:], x["M3n"][1][:, k1:]
    a_tr, a_te = _area_codes(f["cc_tr"], areas, first), _area_codes(f["cc_te"], areas, first)
    ind_tr = np.eye(len(areas) + 1)[a_tr]
    ind_te = np.eye(len(areas) + 1)[a_te]
    inter_tr = np.concatenate([xn_tr * ind_tr[:, [j]] for j in range(len(areas) + 1)], axis=1)
    inter_te = np.concatenate([xn_te * ind_te[:, [j]] for j in range(len(areas) + 1)], axis=1)
    p = _stage2_c(f["z_tr"], f["z_te"], np.column_stack([x["M1c"][0], xn_tr, inter_tr]),
                  np.column_stack([x["M1c"][1], xn_te, inter_te]), f["y_tr"], c)
    return k, c, p.astype(np.float32)


C_GRID8 = (0.01, 0.1, 1.0)


def b8_language() -> dict:
    from concurrent.futures import ProcessPoolExecutor

    vv = v77()
    pt = dict(np.load(vv.POINT_CACHE))
    folds = vv._dump_folds()
    sizes = [len(f["y_te"]) for f in folds]
    jobs = [(k, c) for c in C_GRID8 for k in range(len(folds))]
    with ProcessPoolExecutor(6) as ex:
        res = list(ex.map(_b8_one, jobs))
    y = pt["y"].astype(np.float64)
    unit = pt["unit"].astype(str)
    out: dict = {}
    for c in C_GRID8:
        parts = sorted((k, p) for k, cc, p in res if cc == c)
        pm = np.concatenate([p for _, p in parts]).astype(np.float64)
        if len(pm) != len(y) or [len(p) for _, p in parts] != sizes:
            return {"error": "row order differs from the point refit"}
        preds = {a: pt[f"p_{a}"].astype(np.float64) for a in ("M0", "M1", "M1c")} | {"M3nL": pm}
        out[f"C={c:g}"] = _r_from_cells(y, preds, unit, np.random.default_rng(SEED), comps=(("M1c", "M3nL"),))
    out["rule"] = ("B-17 B8: B1 dump of B-16 (frozen reading); names model = reference M1c + the 15 name features + their "
                   "interactions with the language area (first lexicon block of the nearest settlement's country; blocks "
                   "with at least 2 % of the fold's training cells, others pooled), L2 penalty C on standardised features; "
                   "the reference is the point refit of B1 (C = 1); 4,000 draws of location groups, models fixed")
    return out


# %% B14: 目録を改訂 16 と 17 に
def b14_census() -> dict:
    """B-15 の B8（v76 の b8_census_all）と同じ規則を、v36〜v82 の全ての出力に（改訂 2〜17）。"""
    import re

    v76 = v77().v76
    v68, v69, v73 = v76.v68, v76.v69, v76.v73
    docv = ROOT / "docs" / "verification"

    def amend(v: int) -> str:
        if v <= 76:
            return v76._amend_all(v, "")
        return "16" if v <= 79 else "17"

    rows: list[dict] = []
    for f in sorted(docv.glob("v[3-8][0-9]*_*.json")):
        m = re.match(r"v(\d+)", f.name)
        v = int(m.group(1)) if m else 0
        if not 36 <= v <= 82 or "estimates" in f.name:  # v82 は追記 2（記述）
            continue
        v69._walk_fixed(_json(f), [], rows, f.name, amend(v))
    df = pd.DataFrame(rows)
    df["key"] = (df["kind"] + "|" + df["estimate"].round(6).astype(str) + "|" + df["ci95_lo"].round(6).astype(str)
                 + "|" + df["ci95_hi"].round(6).astype(str))
    df["duplicate"] = df.duplicated("key")
    sep = (df["kind"] == "odds_ratio") & ((np.abs(np.log(df["estimate"].clip(lower=1e-300))) > 10)
                                          | (df["ci95_hi"] / df["ci95_lo"].clip(lower=1e-300) > 1e6))
    df["valid"] = np.isfinite(df["ci95_lo"]) & np.isfinite(df["ci95_hi"]) & (df["ci95_lo"] <= df["ci95_hi"]) & ~sep
    df["family"] = [v68._family(pp, k) for pp, k in zip(df["path"], df["kind"], strict=True)]
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
            "n": len(x), "n_excluding_null": len(ex),
            "by_amendment": {str(k): int(v) for k, v in x["amendment"].astype(str).value_counts().items()},
            "n_bh_q_below_0_05": int((qbh < 0.05).sum()), "n_by_q_below_0_05": int((qby < 0.05).sum()),
            "n_by_q_below_0_10": int((qby < 0.10).sum()),
            "direction_of_null_excluding": {str(k): int(v) for k, v in ex["direction"].value_counts().items()}}
    df = df.merge(pd.concat(qcols), on=["file", "path"], how="left") if qcols else df
    df.drop(columns="key").to_csv(docv / "v80_estimates.csv", index=False)
    out["rule"] = ("B-17 B14: the rules of the census of B8 of B-15 applied to every output v36 to v82 (amendments 2 to "
                   "17); file docs/verification/v80_estimates.csv")
    return out


PARTS = {"B1": b1_coverage, "B2": b2_fixed_margin, "B4": b4_label_variants, "B5": b5_gfm_ems, "B6": b6_gfm_attrs,
         "B6c": b6_check_cluster,
         "B9": b9_dilution,
         "B10": b10_injected_verdicts, "B11": b11_correlation, "B12": b12_flow, "B13": b13_read_subset,
         "B15": b15_multiplicity, "B7": b7_ngrams, "B8": b8_language, "B14": b14_census}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    want = [p for p in (a.only.split(",") if a.only else PARTS) if p]
    for k in want:
        print(f"{k} ...", flush=True)
        t0 = time.time()
        r = PARTS[k]()
        r["seconds"] = round(time.time() - t0)
        res = _json(OUT) if OUT.is_file() else {}
        res[k] = r
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        short = {kk: vv for kk, vv in r.items() if kk not in ("reps",)}
        print(f"{k} ({r['seconds']} s) {json.dumps(short, ensure_ascii=False, default=float)[:2500]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
