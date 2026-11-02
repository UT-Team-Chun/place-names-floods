# %% v83（事後・amendment B-18 の B1・B2）: 初見の査読 第 14 回に応じた点検（utens のジョブで回す）
"""確証判定は変えない。B-18 で計算の前に固定した量を、既にある実行の予測と表から計算する。

- B1: 表 3 の全ての記録について、全ての地面で保った割合の 90% Fieller 区間（B-17 の B9 と同じ作り方、種 20261107）
- B2: 集落の近くの OR の δ（500 m 以内 × 水の名）に寄与する場所グループの実効の数と、寄与の大きい三つの群の割合
  （計画した GFM の推定、同じセルの EMS、上限ありの EMS の全てのセル、上限なしの EMS）

B3（読みの版の一覧）は計算を伴わないので、既にある値から論文の make_figures で表にする（追記 1）。
出力: stdout と docs/verification/v83_b18.json（部分ごとに足す。`--only B1`）。Mac では回さない（2026-10-31 の指示）。
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

from topsem.world.b13 import fieller_interval

OUT = ROOT / "docs" / "verification" / "v83_b18.json"
SEED = 20261107
B_DRAWS = 4000
Z90 = float(stats.norm.ppf(0.95))
_MODS: dict = {}


def _load(name: str):
    if name not in _MODS:
        spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        _MODS[name] = mod
    return _MODS[name]


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# %% B1: 保った割合の Fieller の区間（全ての地面）
# (鍵, 集落を通じた版の実行, 腕の組, セルごとの版の実行, 腕の組)。B-14 の B6（v73 の RECORDS）と同じ、恒久の水を除いた版は除く
RECORDS = (("gsw_from_2022", "world_v2_b13gsw22", ("M1cG", "M3oG"), "world_v2_b13gsw22", ("M1cW", "M3oW")),
           ("ems_earlier_events", "world_v2_b13recprior", ("M1cE", "M3oE"), "world_v2_b13recprior", ("M1cEc", "M3oEc")),
           ("jrc_100yr", "world_v2_b13jrcset", ("M1cJ", "M3oJ"), "world_v2_b13jrcset", ("M1cJc", "M3oJc")),
           ("ems_events_before_after", "world_v2_b7recur", ("M1cE", "M3oE"), "world_v2_b7recur", ("M1cEc", "M3oEc")),
           ("gsw_all_events", "world_v2_b6gsw", ("M1cG", "M3oG"), "world_v2_b8gswcell_local", ("M1cW", "M3oW")))


def b1_share_fieller() -> dict:
    v36 = _load("v36_yardstick_robustness")
    rng = np.random.default_rng(SEED)
    cache: dict = {}
    out: dict = {}
    w = None
    for key, run_s, arms_s, run_c, arms_c in RECORDS:
        for r in (run_s, run_c):
            if r not in cache:
                cache[r] = v36.load(r, "primary")
        ds, dc = cache[run_s], cache[run_c]
        if run_s != run_c and not (np.array_equal(ds["y"], dc["y"]) and np.array_equal(ds["aoi"], dc["aoi"])):
            out[key] = {"note": "rows differ between the two runs"}
            continue
        y = np.asarray(ds["y"], dtype=np.float64)
        u, inv = np.unique(np.asarray(ds["unit"]).astype(str), return_inverse=True)
        g = len(u)
        if w is None or w.shape[1] != g:  # 全ての記録の実行は同じ 195 の場所グループ（同じ重みを使う）
            w = np.vstack([np.ones(g), rng.multinomial(g, np.full(g, 1.0 / g), size=B_DRAWS).astype(np.float64)])

        def gsum(d: dict, a: tuple[str, str], y: np.ndarray = y, inv: np.ndarray = inv, g: int = g) -> np.ndarray:
            e = (np.asarray(d[f"p_{a[0]}"], dtype=np.float64) - y) ** 2 - (np.asarray(d[f"p_{a[1]}"], dtype=np.float64) - y) ** 2
            return np.bincount(inv, weights=e, minlength=g)

        es, ec = w @ gsum(ds, arms_s), w @ gsum(dc, arms_c)
        ok = np.isfinite(es[1:]) & np.isfinite(ec[1:])
        c = np.cov(es[1:][ok], ec[1:][ok])
        fi = fieller_interval(float(es[0]), float(ec[0]), (float(c[0, 0]), float(c[1, 1]), float(c[0, 1])), z=Z90)
        out[key] = {"runs": [run_s, run_c], "arms_settlement": list(arms_s), "arms_cell": list(arms_c),
                    "share_kept": float(es[0] / ec[0]), "fieller90": list(fi) if fi else None,
                    "n_cells": len(y), "n_units": g}
    out["rule"] = ("B-18 B1: share of the per-cell gain kept through settlements on all ground (gain over the record "
                   "run's own reference, models fixed); 90% Fieller interval from the covariance of 4,000 draws of "
                   "location groups (seed 20261107; None when unbounded), as in B9 of B-17")
    return out


# %% B2: δ に寄与する場所グループの実効の数
def _irls(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, iters: int = 50) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """AOI の切片つきのロジスティック回帰（v68 の fit_w と同じ IRLS、重み 1）。切片 α と係数 β と添字を返す。"""
    ua, ia = np.unique(aoi, return_inverse=True)
    a, k = len(ua), z.shape[1]
    rate = np.bincount(ia, weights=y, minlength=a) / np.bincount(ia, minlength=a)
    alpha, beta = np.log(rate / (1 - rate)), np.zeros(k)
    haz, hzz = np.empty((a, k)), np.empty((k, k))
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta)))
        wv, r = p * (1 - p), y - p
        haa = np.bincount(ia, weights=wv, minlength=a)
        for jj in range(k):
            zw = z[:, jj] * wv
            haz[:, jj] = np.bincount(ia, weights=zw, minlength=a)
            hzz[:, jj] = z.T @ zw
        ga, gz = np.bincount(ia, weights=r, minlength=a), z.T @ r
        s_ = hzz - haz.T @ (haz / haa[:, None])
        db = np.linalg.solve(s_, gz - haz.T @ (ga / haa))
        da = (ga - haz @ db) / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max()) < 1e-10:
            break
    return alpha, beta, ia


def _cluster_shares(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, unit: np.ndarray, j: int) -> dict:
    """δ = β_j の場所グループごとの寄与（影響関数 ψ_g = (S⁻¹ u_g)_j）と、分散の割合の実効の数。"""
    alpha, beta, ia = _irls(y, z, aoi)
    p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta)))
    wv = p * (1 - p)
    a, k = len(np.unique(ia)), z.shape[1]
    haa = np.bincount(ia, weights=wv, minlength=a)
    haz, hzz = np.empty((a, k)), np.empty((k, k))
    for jj in range(k):
        zw = z[:, jj] * wv
        haz[:, jj] = np.bincount(ia, weights=zw, minlength=a)
        hzz[:, jj] = z.T @ zw
    sinv = np.linalg.inv(hzz - haz.T @ (haz / haa[:, None]))
    zc = z - (haz / haa[:, None])[ia]
    ug, ig = np.unique(unit, return_inverse=True)
    sc = np.zeros((len(ug), k))
    np.add.at(sc, ig, zc * (y - p)[:, None])
    psi = sc @ sinv[j]
    v = psi ** 2
    share = v / v.sum()
    order = np.argsort(-share)
    g = len(ug)
    se = math.sqrt(float(v.sum()) * g / (g - 1))
    contrib = np.bincount(ig, weights=((np.abs(z[:, j]) > 0) & (y == 1)).astype(np.float64), minlength=g)
    return {"odds_ratio": math.exp(beta[j]), "cr_ci95": [math.exp(beta[j] - 1.96 * se), math.exp(beta[j] + 1.96 * se)],
            "G": g, "G_with_flooded_near": int((contrib > 0).sum()),
            "effective_groups": float(1.0 / (share ** 2).sum()),
            "top3": [{"unit": str(ug[i]), "variance_share": float(share[i]), "flooded_near": int(contrib[i])}
                     for i in order[:3]],
            "n_cells": len(y)}


def _near_arrays(s0: pd.DataFrame, ycol: str) -> tuple:
    v39 = _load("v39_within_aoi_odds")
    y0 = s0[ycol].to_numpy(dtype=np.float64)
    rate = pd.Series(y0).groupby(s0["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = s0.loc[keep].reset_index(drop=True)
    z, names = v39.design(s, ("wetside",), True, banded=True)
    return (y0[keep], z, pd.factorize(s["aoi"])[0].astype(np.int64), s["unit"].to_numpy(str),
            names.index("wetside:le500m"))


def b2_effective_clusters() -> dict:
    v41 = _load("v41_gfm_independent_labels")
    v68 = _load("v68_b11_local")
    v69 = _load("v69_b12_local")
    out: dict = {}
    t = v41.build()
    g = v41.attach_gfm(t, _json(v69.GFM_SUMMARY))
    for key, tab, ycol in (("planned_gfm", g, "y_gfm"), ("ems_same_cells", g, "y"), ("ems_capped_all", t, "y")):
        t0 = time.time()
        out[key] = _cluster_shares(*_near_arrays(tab, ycol))
        out[key]["seconds"] = round(time.time() - t0)
        print(f"  B2 {key}: OR {out[key]['odds_ratio']:.3f}, effective groups {out[key]['effective_groups']:.2f}",
              flush=True)
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    t0 = time.time()
    out["ems_uncapped"] = _cluster_shares(*_near_arrays(u, "y"))
    out["ems_uncapped"]["seconds"] = round(time.time() - t0)
    out["rule"] = ("B-18 B2: Eq. (2) estimator (AOI intercepts, distance bands, water-related x band, logit of M1; AOIs "
                   "with both classes) on the v41 table (frozen reading; GFM-assessable cells for the planned GFM and the "
                   "EMS same-cell estimates; all analysed cells for the capped EMS estimate) and on the uncapped table "
                   "(corrected reading, as v68 B2); influence of delta per location group psi_g = (S^-1 u_g)_j from the "
                   "AOI-profiled score; effective number = 1 / sum of squared variance shares; CR interval with G/(G-1)")
    return out


PARTS = {"B1": b1_share_fieller, "B2": b2_effective_clusters}


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
        print(f"{k} ({r['seconds']} s) {json.dumps(r, ensure_ascii=False, default=float)[:2000]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
