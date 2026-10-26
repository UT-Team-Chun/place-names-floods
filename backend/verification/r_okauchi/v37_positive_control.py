# %% v37（事後・amendment B-2）: 名前は地形を記憶しているか（正の対照）と、名前の区分・距離帯・重み別の記述
"""確証判定は変えない。内部審査の指摘 3・8・9・13 に応じた記述。

A. 解析表（world_eval の主分割と同じ手順・同じ行）を作り直し、確証解析の held-out 予測と行をそろえる
   （data/verification/v37_table.parquet にキャッシュ。v39 も使う）
B. 正の対照（ラベル非依存）: AOI 内で、最寄りの集落名が水系側（WET_SIDE）・各区分のセルとそれ以外の
   セルの、HAND（セル内最小）≤ 2 m の割合・AOI 内の標高の下半分の割合・log(1 + 河道距離) の差。
   全セル・最寄り集落までの距離帯別・集落の点（中心が集落から 177 m 以内のセル = 集落の点の周りの 1–4 セル）の 3 尺度
C. 浸水の対比（探索的）: 区分別・距離帯別に、観測の浸水率の差と、観測 − M1c の残差の差。重みは
   v35（標本内の小さい群のセル数）と、上限をかける前の AOI の大きさに引き伸ばした版の 2 通り
D. 予測の上位 10% セルの M1c と M3n の重なり（全体・AOI 内）
E. 効果量の翻訳: 合成ラベル（§8 の形）の δ 格子での平均 ΔBS(M1c → M3n) から、余白と主分割の 90% 上限に
   当たる δ を単調補間し、オッズ比 e^δ（最寄りの集落名が水系側か否か、地形を固定）で示す
区間はすべて場所グループの cluster bootstrap（B = 4000）。

出力: stdout の表と docs/verification/v37_positive_control.json
"""

from __future__ import annotations

import json
import sys
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_eval as ev

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.dataset import continent_map
from topsem.world.lexicon import CATEGORIES

CACHE = DATA_DIR / "verification" / "v37_table.parquet"
OUT = ROOT / "docs" / "verification" / "v37_positive_control.json"
PRED = RUNS_DIR / "world_v2" / "predictions" / "primary"
B = 4000
SEED = 20261026
MIN_CELLS = 20
MIN_SETTLEMENTS = 2
SESOI = 0.05
MINUTE_TOL = 2e-3  # GeoNames の座標は小数 5 桁: 分の整数なら 60 倍が整数から 3e-4 以内
PC_VARS = ("low_hand", "low_half", "log_river", "hand_mean", "hand_pct")
LOW_HAND_M = 2.0
POINT_M = 177.0  # 250 m セルの中心から角まで: 集落の点を含むセルは必ず入る（点が辺・角の近くなら隣も、最大 4 セル）
BANDS = (
    ("le250m", 0, 250),
    ("250_500m", 250, 500),
    ("500m_1km", 500, 1000),
    ("gt1km", 1000, np.inf),
)
GROUPS = ("wetside", *CATEGORIES)
PRED_ARMS = ("M0c", "M1", "M1c", "M3n")


# %% A. 解析表
NEEDED = (
    "cci_tree",
    "set_lat",
    "set_lon",
)  # B-2 追記 2 で足した列（無い古いキャッシュは作り直す）


def build() -> pd.DataFrame:
    if CACHE.is_file():
        cached = pd.read_parquet(CACHE)
        if all(c in cached.columns for c in NEEDED):
            return cached
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    sens = ev.Sensitivity()
    dirs, _, infos = ev.select_aois(gn, sens)
    df = ev.build_table(dirs, seed=20261025, infos=infos, sens=sens)
    df, _, _ = ev.add_channels(df, gn, None)
    cont = continent_map(ev.ROOT / "geonames" / "countryInfo.txt")
    df["fold"] = ev.folds_for(
        "primary", df, infos, cont, buffer_km=sens.buffer_km
    ).to_numpy()
    df["n_eval_aoi"] = df["aoi"].map({a: infos[a].n_eval for a in df["aoi"].unique()})
    keep = [
        "aoi",
        "unit",
        "fold",
        "y",
        "lon",
        "lat",
        "P_gid1",
        "P_d1",
        "hand_min",
        "hand_mean",
        "elev_mean",
        "dist_river_all",
        "cci_built",
        "cci_tree",
        "cc",
        "n_eval_aoi",
        "lex_wetside_n1",
        *(f"lex_{c}_n1" for c in CATEGORIES),
    ]
    t = pd.concat(
        [df[df["fold"] == k][keep] for k in sorted(df["fold"].unique()) if k >= 0],
        ignore_index=True,
    )
    parts = [
        np.load(f)
        for f in sorted(PRED.glob("fold*.npz"), key=lambda p: int(p.stem[4:]))
    ]
    y = np.concatenate([p["y"] for p in parts])
    aoi = np.concatenate([p["aoi"] for p in parts])
    if not (
        np.array_equal(y, t["y"].to_numpy())
        and np.array_equal(aoi, t["aoi"].to_numpy(dtype=str))
    ):
        msg = "rebuilt table does not match the held-out predictions"
        raise SystemExit(msg)
    for a in PRED_ARMS:
        t[f"p_{a}"] = np.concatenate([p[f"p_{a}"] for p in parts])
    # 最寄りの集落の座標（GeoNames）: 角度の分の整数に丸められた座標の判定に使う
    coords = gn.drop_duplicates("gid").set_index("gid")
    t["set_lat"] = t["P_gid1"].map(coords["lat"]).astype(float)
    t["set_lon"] = t["P_gid1"].map(coords["lon"]).astype(float)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(CACHE)
    return t


def flags(t: pd.DataFrame) -> dict[str, np.ndarray]:
    out = {"wetside": t["lex_wetside_n1"].to_numpy() > 0}
    out |= {c: t[f"lex_{c}_n1"].to_numpy() > 0 for c in CATEGORIES}
    return out


def derived(t: pd.DataFrame) -> pd.DataFrame:
    t = t.copy()
    t["low_hand"] = (t["hand_min"] <= LOW_HAND_M).astype(np.float64)
    t["low_half"] = (
        t["elev_mean"] <= t.groupby("aoi")["elev_mean"].transform("median")
    ).astype(np.float64)
    t["log_river"] = np.log1p(t["dist_river_all"].clip(lower=0))
    t["resid_M1c"] = t["y"] - t["p_M1c"]
    # HAND の連続量: セル平均（m）と AOI 内の百分位（0–1、低いほど排水路に近い）
    t["hand_pct"] = t.groupby("aoi")["hand_mean"].rank(pct=True)
    lat60, lon60 = t["set_lat"] * 60, t["set_lon"] * 60
    t["minute_rounded"] = ((lat60 - lat60.round()).abs() < MINUTE_TOL) & (
        (lon60 - lon60.round()).abs() < MINUTE_TOL
    )
    return t


# %% 対比と bootstrap
def contrast(
    t: pd.DataFrame, var: str, g: np.ndarray, m: np.ndarray, weighting: str
) -> dict | None:
    """AOI ごとの (群 − それ以外) の var の平均差を重み付けで平均。区間は場所グループの bootstrap。"""
    sub = pd.DataFrame(
        {
            "aoi": t["aoi"].to_numpy()[m],
            "unit": t["unit"].to_numpy()[m],
            "v": t[var].to_numpy(dtype=np.float64)[m],
            "g": g[m],
            "n_eval": t["n_eval_aoi"].to_numpy()[m],
        }
    )
    sub = sub[np.isfinite(sub["v"].to_numpy())]  # 欠測（HAND が無いセルなど）は除く
    n_aoi = t.groupby("aoi").size()
    rows = []
    for a, s in sub.groupby("aoi", sort=True):
        gi = s["g"].to_numpy()
        ng, no = int(gi.sum()), int((~gi).sum())
        if ng < MIN_CELLS or no < MIN_CELLS:
            continue
        v = s["v"].to_numpy()
        w = float(min(ng, no))
        if weighting == "uncapped":  # 上限をかける前の AOI の大きさへ引き伸ばす
            w *= float(s["n_eval"].iloc[0]) / float(n_aoi[a])
        rows.append(
            (
                s["unit"].iloc[0],
                w,
                v[gi].mean() - v[~gi].mean(),
                v[gi].mean(),
                v[~gi].mean(),
            )
        )
    if len(rows) < 5:
        return None
    unit = np.array([r[0] for r in rows])
    w, diff, mg, mo = (np.array([r[k] for r in rows]) for k in (1, 2, 3, 4))
    est = float(np.average(diff, weights=w))
    uniq, inv = np.unique(unit, return_inverse=True)
    rng = np.random.default_rng(SEED)
    k = rng.multinomial(len(uniq), np.full(len(uniq), 1 / len(uniq)), size=B)[:, inv]
    ww = k * w
    reps = (ww @ diff) / ww.sum(axis=1)
    return {
        "estimate": est,
        "ci95": [float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))],
        "mean_group": float(np.average(mg, weights=w)),
        "mean_other": float(np.average(mo, weights=w)),
        "n_aoi": len(rows),
        "n_units": len(uniq),
    }


def settlement_table(t: pd.DataFrame, fl: dict[str, np.ndarray]) -> pd.DataFrame:
    """集落の点の尺度: 中心が集落から 177 m 以内のセル（集落の点の周りの 1–4 セル）を集落ごとに平均。"""
    m = (t["P_d1"].to_numpy() <= POINT_M) & (t["P_gid1"].to_numpy() >= 0)
    s = t.loc[m, ["aoi", "unit", "P_gid1", *PC_VARS, "hand_min"]].copy()
    for gname, gv in fl.items():
        s[gname] = gv[m]
    return s.groupby(["aoi", "unit", "P_gid1"], as_index=False).mean()


def settlement_contrast(s: pd.DataFrame, var: str, gname: str) -> dict | None:
    rows, n_set = [], 0
    for _, x in s.groupby("aoi", sort=True):
        gi = x[gname].to_numpy() > 0.5
        if gi.sum() < MIN_SETTLEMENTS or (~gi).sum() < MIN_SETTLEMENTS:
            continue
        v = x[var].to_numpy()
        n_set += len(x)
        rows.append(
            (
                x["unit"].iloc[0],
                float(min(gi.sum(), (~gi).sum())),
                v[gi].mean() - v[~gi].mean(),
                v[gi].mean(),
                v[~gi].mean(),
            )
        )
    if len(rows) < 5:
        return None
    unit = np.array([r[0] for r in rows])
    w, diff = np.array([r[1] for r in rows]), np.array([r[2] for r in rows])
    uniq, inv = np.unique(unit, return_inverse=True)
    rng = np.random.default_rng(SEED)
    k = rng.multinomial(len(uniq), np.full(len(uniq), 1 / len(uniq)), size=B)[:, inv]
    reps = ((k * w) @ diff) / (k * w).sum(axis=1)
    return {
        "estimate": float(np.average(diff, weights=w)),
        "ci95": [float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))],
        "mean_group": float(np.average([r[3] for r in rows], weights=w)),
        "mean_other": float(np.average([r[4] for r in rows], weights=w)),
        "n_aoi": len(rows),
        "n_settlements": n_set,
    }


# %% F. 集落からの距離帯ごとの主比較の増分
def attach_arms(
    t: pd.DataFrame, arms: tuple[str, ...] = ("M0", "M3np")
) -> pd.DataFrame:
    """キャッシュに無い腕の予測を held-out 予測から足す（キャッシュと同じ fold 順）。"""
    parts = [
        np.load(f)
        for f in sorted(PRED.glob("fold*.npz"), key=lambda p: int(p.stem[4:]))
    ]
    if not np.array_equal(np.concatenate([p["y"] for p in parts]), t["y"].to_numpy()):
        msg = "cached table and predictions are out of order"
        raise SystemExit(msg)
    t = t.copy()
    for a in arms:
        t[f"p_{a}"] = np.concatenate([p[f"p_{a}"] for p in parts])
    return t


def band_increment(t: pd.DataFrame, m: np.ndarray) -> dict:
    """層内の ΔBS(M1c → M3n)・placebo 差・層内の物差し ΔBS(M0 → M1) と、場所グループの bootstrap 区間。"""
    y = t["y"].to_numpy(dtype=np.float64)[m]
    uniq, inv = np.unique(t["unit"].to_numpy()[m], return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    s = {
        a: np.bincount(
            inv, weights=(t[f"p_{a}"].to_numpy(dtype=np.float64)[m] - y) ** 2
        )
        for a in ("M0", "M1", "M1c", "M3n", "M3np")
    }
    hat = {a: v.sum() / n.sum() for a, v in s.items()}
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(n), np.full(len(n), 1 / len(n)), size=B).astype(np.float64)
    bb = {a: (w @ v) / (w @ n) for a, v in s.items()}
    dn, y0 = hat["M1c"] - hat["M3n"], hat["M0"] - hat["M1"]
    bd = bb["M1c"] - bb["M3n"]
    c90 = [float(np.quantile(bd, 0.05)), float(np.quantile(bd, 0.95))]
    # B-2 追記 4: 95% 区間と、real − permuted の差（ΔBS の尺度、同じ bootstrap）
    dp = bb["M3np"] - bb["M3n"]
    return {
        "n_cells": int(m.sum()),
        "n_units": len(uniq),
        "wet_rate": float(y.mean()),
        "delta_bs_names": float(dn),
        "ci90_delta": c90,
        "ci95_delta": [float(np.quantile(bd, 0.025)), float(np.quantile(bd, 0.975))],
        "real_minus_permuted_delta_bs": float(hat["M3np"] - hat["M3n"]),
        "ci95_real_minus_permuted": [
            float(np.quantile(dp, 0.025)),
            float(np.quantile(dp, 0.975)),
        ],
        "yardstick_M0_M1": float(y0),
        "R": float(dn / y0) if y0 > 0 else None,
        "ci90_R": [c / y0 for c in c90] if y0 > 0 else None,
        "equivalent": bool(y0 > 0 and -SESOI * y0 < c90[0] and c90[1] < SESOI * y0),
        "R_minus_placebo": float((hat["M3np"] - hat["M3n"]) / y0) if y0 > 0 else None,
    }


# %% D. 上位 10% の重なり
def top_overlap(
    t: pd.DataFrame, a: str = "p_M1c", b: str = "p_M3n", share: float = 0.1
) -> dict:
    pa, pb = t[a].to_numpy(), t[b].to_numpy()
    ta, tb = pa >= np.quantile(pa, 1 - share), pb >= np.quantile(pb, 1 - share)
    within = []
    for idx in t.groupby("aoi").indices.values():
        if len(idx) < 50:
            continue
        xa, xb = pa[idx], pb[idx]
        sa, sb = xa >= np.quantile(xa, 1 - share), xb >= np.quantile(xb, 1 - share)
        within.append((sa & sb).sum() / max(sa.sum(), 1))
    return {
        "pooled": float((ta & tb).sum() / ta.sum()),
        "within_aoi_median": float(np.median(within)),
        "within_aoi_p10": float(np.quantile(within, 0.1)),
        "n_aoi": len(within),
    }


# %% E. 効果量の翻訳
GRID = {
    0.0: "power_v3_null",
    0.2: "power_v3_d020",
    0.4: "power_v3_d040",
    0.6: "power_v3_d060",
    0.8: "power_v3_d080",
    1.2: "power_v3_sesoi",
    1.6: "power_v3_d160",
}
GRID_LOCAL = {
    0.0: "power_v3_null",
    1.0: "power_v3_loc_d100",
    2.0: "power_v3_loc_d200",
    3.0: "power_v3_loc_d300",
}  # 効果を集落から 500 m 以内に限った合成ラベル（B-2 追記 2）


def power_curve(runs: dict[float, str] = GRID) -> dict:
    """δ ごとの平均 ΔBS(M1c → M3n)（合成ラベル）とモンテカルロの 95% 区間。

    power_summary.json に回ごとの ΔBS があればそれを使い、無ければ予測ファイル（sim/rep*.npz）から計算する。
    """
    out = {}
    for delta, run in runs.items():
        s = RUNS_DIR / run / "power_summary.json"
        summ = json.loads(s.read_text()) if s.is_file() else {}
        reps = summ.get("reps") or []
        if reps and "delta_bs" in reps[0]:
            dn = [r["delta_bs"] for r in reps]
            d01 = [r["delta_bs_M0_M1"] for r in reps]
        else:
            files = sorted((RUNS_DIR / run / "sim").glob("rep*.npz"))
            if not files:
                continue
            dn, d01 = [], []
            for f in files:
                with np.load(f) as z:
                    y = z["y"].astype(np.float64)
                    bs = {
                        a: float(np.mean((z[a].astype(np.float64) - y) ** 2))
                        for a in ("M0", "M1", "M1c", "M3n")
                    }
                dn.append(bs["M1c"] - bs["M3n"])
                d01.append(bs["M0"] - bs["M1"])
        mu, sd = float(np.mean(dn)), (float(np.std(dn, ddof=1)) if len(dn) > 1 else 0.0)
        se = sd / np.sqrt(len(dn))
        out[delta] = {
            "run": run,
            "n_reps": len(dn),
            "mean_delta_bs": mu,
            "sd_delta_bs": sd,
            "mc_ci95_mean": [mu - 1.96 * se, mu + 1.96 * se],
            "mean_yardstick": float(np.mean(d01)),
        }
    return out


def null_coverage() -> dict | None:
    """合成ラベル δ = 0 の各回の区間が、回の間の平均 R（設計の真値の目安）を含む割合（B-2 追記 2）。"""
    s = RUNS_DIR / "power_v3_null" / "power_summary.json"
    if not s.is_file():
        return None
    reps = json.loads(s.read_text())["reps"]
    r = np.array([x["R"] for x in reps])
    mu = float(r.mean())
    cov90 = float(np.mean([x["ci90"][0] <= mu <= x["ci90"][1] for x in reps]))
    cov95 = float(np.mean([x["ci95"][0] <= mu <= x["ci95"][1] for x in reps]))
    se = np.array([(x["ci90"][1] - x["ci90"][0]) / (2 * 1.645) for x in reps])
    se_boot = float(np.median(se))
    # B-2 追記 3: 0 を含む割合（名前の効果が無いという意味での真値）、頑健な散らばり、標準化した偏差の SD。
    # 回の間の SD / SE の中央値（2.1）は、SE の大きい少数の回に引っ張られる
    mad = float(1.4826 * np.median(np.abs(r - np.median(r))))
    z = (r - mu) / se
    return {
        "n_reps": len(reps),
        "mean_R": mu,
        "sd_R_between_reps": float(r.std(ddof=1)),
        "median_bootstrap_se_R": se_boot,
        "coverage90": cov90,
        "coverage95": cov95,
        "coverage90_zero": float(
            np.mean([x["ci90"][0] <= 0 <= x["ci90"][1] for x in reps])
        ),
        "coverage95_zero": float(
            np.mean([x["ci95"][0] <= 0 <= x["ci95"][1] for x in reps])
        ),
        "robust_sd_R": mad,
        "sd_standardised_deviation": float(z.std(ddof=1)),
    }


def invert(curve: dict, target: float) -> float | None:
    """平均 ΔBS が target になる δ（隣り合う格子点の線形補間。曲線が単調な区間で）。"""
    pts = sorted((d, v["mean_delta_bs"]) for d, v in curve.items())
    for (d0, v0), (d1, v1) in pairwise(pts):
        if v0 <= target <= v1 and v1 > v0:
            return d0 + (target - v0) * (d1 - d0) / (v1 - v0)
    return None


# %% 実行
def main() -> None:
    t = derived(build())
    fl = flags(t)
    d1 = t["P_d1"].to_numpy()
    allm = np.ones(len(t), dtype=bool)
    res: dict = {
        "n_cells": len(t),
        "share_by_group": {g: float(v.mean()) for g, v in fl.items()},
        "nearest_settlement_distance_m": {
            q: float(np.quantile(d1[np.isfinite(d1)], q)) for q in (0.25, 0.5, 0.75)
        },
        "share_within_1km": float(np.mean(d1 <= 1000)),
    }
    # B. 正の対照（ラベル非依存）
    pc: dict = {
        "cells": {},
        "bands": {},
        "settlements": {},
        "settlements_not_rounded": {},
    }
    for g in GROUPS:  # 全セル（B-2 で登録した尺度）
        pc["cells"][g] = {v: contrast(t, v, fl[g], allm, "capped") for v in PC_VARS}
    for name, lo, hi in BANDS:  # 距離帯（登録外、v37 を書く段階で追加）
        m = (d1 > lo) & (d1 <= hi)
        pc["bands"][name] = {
            v: contrast(t, v, fl["wetside"], m, "capped") for v in PC_VARS
        }
    st = settlement_table(t, fl)  # 集落の点（登録外、v37 を書く段階で追加）
    pc["n_settlement_points"] = len(st)
    for g in GROUPS:
        pc["settlements"][g] = {v: settlement_contrast(st, v, g) for v in PC_VARS}
    # 座標が角度の分の整数に丸められた集落（最大 ±0.9 km）を除いた集落の点（B-2 追記 2）
    rounded = t["minute_rounded"].to_numpy()
    st_nr = settlement_table(t[~rounded], {g: v[~rounded] for g, v in fl.items()})
    pc["n_settlement_points_not_rounded"] = len(st_nr)
    pc["share_cells_minute_rounded"] = float(rounded.mean())
    pc["share_settlement_points_minute_rounded"] = float(1 - len(st_nr) / len(st))
    for g in ("wetside", "high"):
        pc["settlements_not_rounded"][g] = {
            v: settlement_contrast(st_nr, v, g) for v in PC_VARS
        }
    res["positive_control"] = pc
    # C. 浸水の対比（探索的）
    fc: dict = {"by_group": {}, "by_band": {}}
    for wt in ("capped", "uncapped"):
        fc["by_group"][wt] = {
            g: {v: contrast(t, v, fl[g], allm, wt) for v in ("y", "resid_M1c")}
            for g in GROUPS
        }
        fc["by_band"][wt] = {}
        for name, lo, hi in BANDS:
            m = (d1 > lo) & (d1 <= hi)
            fc["by_band"][wt][name] = {
                v: contrast(t, v, fl["wetside"], m, wt) for v in ("y", "resid_M1c")
            }
    # 物差し: AOI 内で HAND が中央値以下のセル（v35 と同じ）
    low = (t["hand_min"] <= t.groupby("aoi")["hand_min"].transform("median")).to_numpy()
    fc["hand_yardstick"] = {
        wt: contrast(t, "y", low, allm, wt) for wt in ("capped", "uncapped")
    }
    res["flood_contrast"] = fc
    # D・E
    res["top_decile_overlap_M1c_M3n"] = top_overlap(t)
    ta = attach_arms(t)
    res["primary_by_distance"] = {
        name: band_increment(ta, (d1 > lo) & (d1 <= hi))
        for name, lo, hi in (*BANDS, ("le500m", 0, 500), ("le1km", 0, 1000))
    }
    res["primary_by_class"] = {
        g: band_increment(ta, fl[g]) for g in GROUPS
    }  # B-2 で登録
    # 集落から 500 m 以内で建物群（CCI built > 10%）や樹木（CCI tree > 50%）の多いセルを除いた対比（B-2 追記 2）
    clean = ((t["cci_built"] <= 0.1) & (t["cci_tree"] <= 0.5)).to_numpy()
    near = d1 <= 500
    res["near_settlement_clean"] = {
        "share_cells_removed_within_500m": float(1 - clean[near].mean()),
        "flood_contrast": {
            wt: {
                v: contrast(t, v, fl["wetside"], near & clean, wt)
                for v in ("y", "resid_M1c")
            }
            for wt in ("capped", "uncapped")
        },
        "primary": band_increment(ta, near & clean),
        "primary_wetside": band_increment(ta, near & clean & fl["wetside"]),
    }
    res["null_coverage"] = null_coverage()
    res["primary_by_distance_wetside"] = {
        name: band_increment(ta, (d1 > lo) & (d1 <= hi) & fl["wetside"])
        for name, lo, hi in (("le500m", 0, 500), ("le1km", 0, 1000))
    }
    curve = power_curve(GRID)
    curve_local = power_curve(GRID_LOCAL)
    m90 = json.loads((RUNS_DIR / "world_v2" / "metrics_primary.json").read_text())[
        "primary"
    ]
    targets = {
        "margin": m90["margin_delta_bs"],
        "achieved_upper90": m90["ci90_delta"][1],
    }
    v36p = ROOT / "docs" / "verification" / "v36_yardstick_robustness.json"
    if v36p.is_file():  # 二方向 bootstrap の 90% 上限も翻訳する（指摘 B5）
        tw = json.loads(v36p.read_text())["splits"]["primary"]["two_way"]["primary"]
        targets["two_way_upper90"] = tw["ci90_delta"][1]

    def translate(cv: dict) -> dict:
        """目標の ΔBS に当たる δ と、平均 ΔBS のモンテカルロ区間の端で補間した δ の幅。"""
        lo_c = {
            k: v | {"mean_delta_bs": v["mc_ci95_mean"][1]} for k, v in cv.items()
        }  # 平均が大きい側 → δ 小
        hi_c = {k: v | {"mean_delta_bs": v["mc_ci95_mean"][0]} for k, v in cv.items()}
        out = {}
        for k, v in targets.items():
            d = invert(cv, v)
            dr = [invert(lo_c, v), invert(hi_c, v)]
            top = max(cv.items(), key=lambda kv: kv[0]) if cv else None
            out[k] = {
                "delta": d,
                "odds_ratio": float(np.exp(d)) if d is not None else None,
                "delta_mc_range": dr,
                "odds_ratio_mc_range": [
                    float(np.exp(x)) if x is not None else None for x in dr
                ],
                # 格子の最大の δ でも目標に届かないとき（オッズ比は e^{δmax} より大きい、とだけ言える）
                "beyond_grid": bool(top and v > top[1]["mean_delta_bs"]),
                "delta_max": top[0] if top else None,
            }
        return out

    by_run = {}
    for run in (
        "world_v2",
        "world_v2_fs1",
        "world_v2_fs2",
        "world_v2_s6_buffer_b1",
        "world_v2_actgrp",
    ):
        mp = RUNS_DIR / run / "metrics_primary.json"
        if not mp.is_file():
            continue
        up = json.loads(mp.read_text())["primary"]["ci90_delta"][1]
        one = {"upper90": up}
        for kind, cv in (("homogeneous", curve), ("local_500m", curve_local)):
            d = invert(cv, up)
            top = max(cv.items(), key=lambda kv: kv[0]) if cv else None
            one[kind] = {
                "delta": d,
                "odds_ratio": float(np.exp(d)) if d is not None else None,
                "beyond_grid": bool(top and up > top[1]["mean_delta_bs"]),
                "below_grid": bool(
                    d is None and not (top and up > top[1]["mean_delta_bs"])
                ),
                "delta_max": top[0] if top else None,
            }
        by_run[run] = one

    res["effect_size_translation"] = {
        # B-2 追記 3: 観測のクラスタ 90% 上限の翻訳を、fold の引き方・束ね方の違う実行ごとに
        "by_run": by_run,
        "curve": {str(k): v for k, v in curve.items()},
        "curve_local": {str(k): v for k, v in curve_local.items()},
        "targets": targets,
        "homogeneous": translate(curve),
        "local_500m": translate(curve_local),
        # 以前の出力形（図 5d・マクロが参照）
        "delta": {k: invert(curve, v) for k, v in targets.items()},
        "odds_ratio": {
            k: (float(np.exp(x)) if (x := invert(curve, v)) is not None else None)
            for k, v in targets.items()
        },
        "share_cells_wetside": float(fl["wetside"].mean()),
        "share_wetside_cells_within_500m": float(
            (fl["wetside"] & (d1 <= 500)).sum() / fl["wetside"].sum()
        ),
    }
    show(res)
    show_round2(res)
    show_round3(res)
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"-> {OUT}")


def fmt(r: dict | None, scale: float = 100.0, unit: str = "pp") -> str:
    if r is None:
        return "   (too few AOIs)"
    return (
        f"{scale * r['estimate']:+6.2f} {unit} [{scale * r['ci95'][0]:+6.2f}, {scale * r['ci95'][1]:+6.2f}]"
        f" ({r['n_aoi']} AOIs)"
    )


def show(res: dict) -> None:
    print(
        f"cells {res['n_cells']:,}; nearest-settlement distance quartiles (m) "
        f"{res['nearest_settlement_distance_m']}; within 1 km {res['share_within_1km']:.2f}"
    )
    print(
        "share of cells by nearest-name group:",
        {k: round(v, 4) for k, v in res["share_by_group"].items()},
    )
    pc = res["positive_control"]
    print("\nB. positive control (label-free), group − other within AOI")
    print(
        f"{'group':9s} {'scale':11s} {'HAND≤2m share':>34s} {'lower-half elev share':>34s} {'log(1+river dist)':>30s}"
    )
    for g in GROUPS:
        for scale in ("cells", "settlements"):
            r = pc[scale][g]
            print(
                f"{g:9s} {scale:11s} {fmt(r['low_hand']):>34s} {fmt(r['low_half']):>34s} "
                f"{fmt(r['log_river'], 1.0, ''):>30s}"
            )
    print(f"settlement points: {pc['n_settlement_points']:,}")
    for name, r in pc["bands"].items():
        print(
            f"  wetside, cells {name:9s} HAND≤2m {fmt(r['low_hand'])}  log river {fmt(r['log_river'], 1.0, '')}"
        )
    fc = res["flood_contrast"]
    print(
        "\nC. flood contrast (exploratory), group − other within AOI: observed flooded share / residual vs M1c"
    )
    for wt in ("capped", "uncapped"):
        print(f"  weighting: {wt}")
        for g in GROUPS:
            r = fc["by_group"][wt][g]
            print(f"    {g:9s} observed {fmt(r['y'])}   residual {fmt(r['resid_M1c'])}")
        for name, r in fc["by_band"][wt].items():
            print(
                f"    wetside {name:9s} observed {fmt(r['y'])}   residual {fmt(r['resid_M1c'])}"
            )
        print(
            f"    HAND yardstick (low-HAND half − other) {fmt(fc['hand_yardstick'][wt])}"
        )
    print("\nD. top-decile overlap M1c vs M3n:", res["top_decile_overlap_M1c_M3n"])
    print(
        "\nF. primary increment by distance to the nearest settlement (R within band, 90% CI)"
    )
    for key in ("primary_by_distance", "primary_by_distance_wetside"):
        for name, r in res[key].items():
            ci = (
                f"[{r['ci90_R'][0]:+.4f}, {r['ci90_R'][1]:+.4f}]" if r["ci90_R"] else ""
            )
            print(
                f"  {key[-8:]:8s} {name:9s} cells {r['n_cells']:9,d} units {r['n_units']:3d} wet {r['wet_rate']:.3f}"
                f"  dBS {r['delta_bs_names']:+.2e}  R {r['R']:+.4f} {ci} eq={r['equivalent']}"
            )
    e = res["effect_size_translation"]
    print(
        "\nE. synthetic δ curves (mean ΔBS M1c→M3n, Monte-Carlo 95% interval of the mean)"
    )
    for kind, key in (("homogeneous", "curve"), ("local_500m", "curve_local")):
        print(
            f"   {kind}: "
            + ", ".join(
                f"δ={k} {v['mean_delta_bs']:.2e} [{v['mc_ci95_mean'][0]:.2e}, {v['mc_ci95_mean'][1]:.2e}] (n={v['n_reps']})"
                for k, v in e[key].items()
            )
        )
        for tgt, r in e[kind].items():
            rng = r["odds_ratio_mc_range"]
            print(
                f"     {tgt:17s} ΔBS {e['targets'][tgt]:.2e} → OR {r['odds_ratio']}  (MC range {rng})"
            )


def show_round3(res: dict) -> None:
    """B-2 追記 3 で足した量の表示。"""
    e = res["effect_size_translation"]
    print(
        "\nround 3: translation of the observed cluster 90% bound by run (homogeneous / within 500 m)"
    )
    for run, r in e["by_run"].items():
        h, loc = r["homogeneous"], r["local_500m"]
        f = lambda x: (
            f"OR {x['odds_ratio']:.2f}"
            if x["odds_ratio"]
            else ("beyond grid" if x["beyond_grid"] else "below grid")
        )
        print(f"  {run:24s} upper90 {r['upper90']:.2e}  {f(h)} / {f(loc)}")
    c = res["null_coverage"]
    print(
        f"  null coverage of zero: 90% {c['coverage90_zero']:.2f}, 95% {c['coverage95_zero']:.2f}; robust SD {c['robust_sd_R']:.4f}"
        f" vs median SE {c['median_bootstrap_se_R']:.4f}; SD of standardised deviations {c['sd_standardised_deviation']:.2f}"
    )


def show_round2(res: dict) -> None:
    """B-2 追記 2 で足した量の表示。"""
    pc = res["positive_control"]
    print("\nround 2: continuous measures and base rates (group mean / other mean)")
    for scale in ("cells", "settlements", "settlements_not_rounded"):
        for g in ("wetside", "high"):
            r = pc[scale].get(g)
            if not r:
                continue
            cells = []
            for v in PC_VARS:
                x = r.get(v)
                if x:
                    cells.append(
                        f"{v} {x['estimate']:+.3f} [{x['ci95'][0]:+.3f},{x['ci95'][1]:+.3f}]"
                        f" ({x['mean_group']:.3f}/{x['mean_other']:.3f})"
                    )
            print(f"  {scale:24s} {g:8s} " + " | ".join(cells))
    print(
        f"  minute-rounded: cells {pc['share_cells_minute_rounded']:.3f}, settlement points "
        f"{pc['share_settlement_points_minute_rounded']:.3f} (kept {pc['n_settlement_points_not_rounded']})"
    )
    print("  primary increment by name class (R within class, 90%):")
    for g, b in res["primary_by_class"].items():
        print(
            f"    {g:9s} cells {b['n_cells']:8,d} R {b['R']:+.4f} [{b['ci90_R'][0]:+.4f}, {b['ci90_R'][1]:+.4f}]"
            f" eq={b['equivalent']}"
        )
    nc = res["near_settlement_clean"]
    print(
        f"  within 500 m, built-up/tree cells removed ({nc['share_cells_removed_within_500m']:.2f} of cells):"
    )
    for wt, d in nc["flood_contrast"].items():
        print(f"    {wt:8s} observed {fmt(d['y'])}  residual {fmt(d['resid_M1c'])}")
    for k in ("primary", "primary_wetside"):
        b = nc[k]
        print(
            f"    {k:16s} R {b['R']:+.4f} [{b['ci90_R'][0]:+.4f}, {b['ci90_R'][1]:+.4f}] cells {b['n_cells']:,}"
        )
    print("  null coverage:", res["null_coverage"])


if __name__ == "__main__":
    main()
