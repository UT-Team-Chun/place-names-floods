# %% v39（事後・amendment B-2）: 同じ AOI の中で、水系の名の集落の近くは何倍浸水しやすいか — 区分別・距離帯別のオッズ比
"""確証判定は変えない。内部審査の指摘 3・9・14（解釈できる効果量）に応じた記述。

AOI の固定効果つきロジスティック回帰を、主分割の全セル（確証解析と同じ 1,507,486 セル）で当てる:

    logit P(y = 1) = α_AOI + Σ_b η_b·[距離帯 b] + Σ_{c,b} γ_{c,b}·[最寄り名の区分 c]·[距離帯 b]   （素の関連）
    同上 + β·logit(p̂_M1)                                                                 （地形を固定）

p̂_M1 は確証解析の held-out 予測（地形・水文・土地被覆・集落配置の柔軟モデル）。e^γ は、同じ AOI・同じ
距離帯の中で、最寄りの集落名が区分 c のセルの浸水のオッズが、そうでないセルの何倍か。距離帯は最寄りの
集落の点まで ≤ 500 m・500 m–1 km・> 1 km。区間は場所グループのクラスタ頑健なサンドイッチ分散
（G/(G − 1) の補正）による 95% 区間。入力は v37 のキャッシュ（data/verification/v37_table.parquet）。

出力: stdout の表と docs/verification/v39_within_aoi_odds.json
`--bootstrap N` で、水系側の名のモデルの区間を場所グループの cluster bootstrap（N 回）でも出す（サンドイッチ
分散の確認）。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR
from topsem.world.lexicon import CATEGORIES

TABLE = DATA_DIR / "verification" / "v37_table.parquet"
OUT = ROOT / "docs" / "verification" / "v39_within_aoi_odds.json"
BANDS = (("le500m", 0, 500), ("500m_1km", 500, 1000), ("gt1km", 1000, np.inf))
EPS = 1e-5


# %% 固定効果つきロジスティック回帰（IRLS、AOI ダミーは対角ブロック）
def fit(
    y: np.ndarray, z: np.ndarray, aoi: np.ndarray, groups: np.ndarray, iters: int = 50
) -> dict:
    """AOI 固定効果 + 稠密な説明変数 z の最尤推定と、groups でのクラスタ頑健な分散。"""
    ua, ia = np.unique(aoi, return_inverse=True)
    n, k, a = len(y), z.shape[1], len(ua)
    d = sparse.csr_matrix((np.ones(n), (np.arange(n), ia)), shape=(n, a))
    rate = np.bincount(ia, weights=y) / np.bincount(ia)
    alpha = np.log(rate / (1 - rate))
    beta = np.zeros(k)
    for _ in range(iters):
        eta = alpha[ia] + z @ beta
        p = 1 / (1 + np.exp(-eta))
        w = p * (1 - p)
        r = y - p
        # ヘッセ行列のブロック: D'WD は対角、D'WZ、Z'WZ
        haa = np.bincount(ia, weights=w)
        haz = d.T @ (z * w[:, None])
        hzz = z.T @ (z * w[:, None])
        ga, gz = np.bincount(ia, weights=r), z.T @ r
        # Schur 補行列で z のブロックを解き、AOI の切片を戻す
        s = hzz - haz.T @ (haz / haa[:, None])
        db = np.linalg.solve(s, gz - haz.T @ (ga / haa))
        da = (ga - haz @ db) / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max()) < 1e-10:
            break
    eta = alpha[ia] + z @ beta
    p = 1 / (1 + np.exp(-eta))
    w = p * (1 - p)
    haa = np.bincount(ia, weights=w)
    haz = d.T @ (z * w[:, None])
    hzz = z.T @ (z * w[:, None])
    s = hzz - haz.T @ (haz / haa[:, None])
    # z のブロックの逆ヘッセ（Schur 補行列の逆）と、AOI を射影したスコアのクラスタ和
    sinv = np.linalg.inv(s)
    r = y - p
    # 切片を最適化した後の z のスコア寄与: x_i − (D'WZ)_a / (D'WD)_a で中心化（AOI は場所グループの中に入れ子）
    zc = z - (haz / haa[:, None])[ia]
    ug, ig = np.unique(groups, return_inverse=True)
    u = np.zeros((len(ug), k))
    np.add.at(u, ig, zc * r[:, None])
    g = len(ug)
    meat = (u.T @ u) * g / (g - 1)
    v = sinv @ meat @ sinv
    return {"beta": beta, "se": np.sqrt(np.diag(v)), "cov": v, "n": n, "n_aoi": a, "n_groups": g}


def design(
    t: pd.DataFrame, classes: tuple[str, ...], offset: bool, *, banded: bool = True
) -> tuple[np.ndarray, list[str]]:
    """banded=False は距離帯に分けない区分の指示変数（B-2 で登録した「全セルでの区分別オッズ比」）。"""
    d1 = t["P_d1"].to_numpy()
    cols, names = [], []
    if banded:
        for bname, lo, hi in BANDS[1:]:  # 基準の帯は ≤ 500 m
            cols.append(((d1 > lo) & (d1 <= hi)).astype(np.float64))
            names.append(f"band_{bname}")
    for c in classes:
        g = (t["lex_wetside_n1"] if c == "wetside" else t[f"lex_{c}_n1"]).to_numpy() > 0
        if not banded:
            cols.append(g.astype(np.float64))
            names.append(f"{c}:all")
            continue
        for bname, lo, hi in BANDS:
            cols.append((g & (d1 > lo) & (d1 <= hi)).astype(np.float64))
            names.append(f"{c}:{bname}")
    if offset:
        p = np.clip(t["p_M1"].to_numpy(dtype=np.float64), EPS, 1 - EPS)
        cols.append(np.log(p / (1 - p)))
        names.append("logit_M1")
    return np.column_stack(cols), names


def table(res: dict, names: list[str]) -> dict:
    out = {}
    for j, nm in enumerate(names):
        b, se = float(res["beta"][j]), float(res["se"][j])
        out[nm] = {
            "coef": b,
            "se": se,
            "odds_ratio": float(np.exp(b)),
            "ci95_odds_ratio": [
                float(np.exp(b - 1.96 * se)),
                float(np.exp(b + 1.96 * se)),
            ],
        }
    return out


# %% サンドイッチ分散の確認: 場所グループの cluster bootstrap（水系側の名のモデルだけ）
def bootstrap_check(t: pd.DataFrame, reps: int, seed: int = 20261026) -> dict:
    y = t["y"].to_numpy(dtype=np.float64)
    aoi, grp = t["aoi"].to_numpy(dtype=str), t["unit"].to_numpy(dtype=str)
    ug, ig = np.unique(grp, return_inverse=True)
    rows = {k: np.flatnonzero(ig == k) for k in range(len(ug))}
    rng = np.random.default_rng(seed)
    out: dict = {}
    for offset in (False, True):
        z, names = design(t, ("wetside",), offset)
        draws = []
        for _ in range(reps):
            pick = rng.integers(0, len(ug), len(ug))
            idx = np.concatenate([rows[k] for k in pick])
            # 同じ AOI が複数回選ばれたら別の AOI として扱う（固定効果を複製ごとに置く）
            rep_id = np.concatenate(
                [np.full(len(rows[k]), j) for j, k in enumerate(pick)]
            )
            a2 = np.char.add(np.char.add(aoi[idx], "#"), rep_id.astype(str))
            r = fit(y[idx], z[idx], a2, grp[idx], iters=30)
            draws.append(r["beta"])
        d = np.array(draws)
        key = "adjusted_M1" if offset else "raw"
        out[key] = {
            nm: {
                "odds_ratio_ci95_bootstrap": [
                    float(np.exp(np.quantile(d[:, j], 0.025))),
                    float(np.exp(np.quantile(d[:, j], 0.975))),
                ]
            }
            for j, nm in enumerate(names)
            if nm.startswith("wetside")
        }
    out["reps"] = reps
    return out


# %% 実行
def main() -> None:
    t = pd.read_parquet(TABLE)
    # 標本の中で浸水セルが 0（または全て）の AOI は固定効果が発散し、係数に寄与しないので除く
    rate = t.groupby("aoi")["y"].transform("mean")
    dropped = sorted(t.loc[(rate == 0) | (rate == 1), "aoi"].unique().tolist())
    t = t[(rate > 0) & (rate < 1)].reset_index(drop=True)
    y = t["y"].to_numpy(dtype=np.float64)
    aoi, grp = t["aoi"].to_numpy(dtype=str), t["unit"].to_numpy(dtype=str)
    res = {"n_cells": len(t), "aoi_dropped_no_variation": dropped}
    for label, classes in (
        ("wetside", ("wetside",)),
        ("categories", tuple(CATEGORIES)),
    ):
        for offset in (False, True):
            z, names = design(t, classes, offset)
            r = fit(y, z, aoi, grp)
            key = f"{label}_{'adjusted_M1' if offset else 'raw'}"
            res[key] = {
                "n_aoi": r["n_aoi"],
                "n_groups": r["n_groups"],
                "terms": table(r, names),
            }
    # 全セルでの区分別オッズ比（距離帯なし、B-2 で登録）
    for label, classes in (
        ("wetside_allcells", ("wetside",)),
        ("categories_allcells", tuple(CATEGORIES)),
    ):
        for offset in (False, True):
            z, names = design(t, classes, offset, banded=False)
            r = fit(y, z, aoi, grp)
            res[f"{label}_{'adjusted_M1' if offset else 'raw'}"] = {
                "n_aoi": r["n_aoi"],
                "n_groups": r["n_groups"],
                "terms": table(r, names),
            }
    # 建物群（CCI built > 10%）・樹木（CCI tree > 50%）の多いセルを除いた推定（B-2 追記 2）
    keep = ((t["cci_built"] <= 0.1) & (t["cci_tree"] <= 0.5)).to_numpy()
    tc = t[keep]
    rc = tc.groupby("aoi")["y"].transform("mean")
    tc = tc[(rc > 0) & (rc < 1)].reset_index(drop=True)
    res["clean_share_cells_kept"] = float(keep.mean())
    for offset in (False, True):
        z, names = design(tc, ("wetside",), offset)
        r = fit(
            tc["y"].to_numpy(dtype=np.float64),
            z,
            tc["aoi"].to_numpy(dtype=str),
            tc["unit"].to_numpy(dtype=str),
        )
        res[f"wetside_clean_{'adjusted_M1' if offset else 'raw'}"] = {
            "n_aoi": r["n_aoi"],
            "n_groups": r["n_groups"],
            "terms": table(r, names),
        }
    counts = {}
    d1 = t["P_d1"].to_numpy()
    for c in ("wetside", *CATEGORIES):
        g = (t["lex_wetside_n1"] if c == "wetside" else t[f"lex_{c}_n1"]).to_numpy() > 0
        counts[c] = {b: int((g & (d1 > lo) & (d1 <= hi)).sum()) for b, lo, hi in BANDS}
    res["cells_by_class_and_band"] = counts
    if "--bootstrap" in sys.argv:
        reps = int(sys.argv[sys.argv.index("--bootstrap") + 1])
        res["bootstrap_check_wetside"] = bootstrap_check(t, reps)
        print("bootstrap check:", json.dumps(res["bootstrap_check_wetside"], indent=0))
    show(res)
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"-> {OUT}")


def show(res: dict) -> None:
    for key in (
        "wetside_allcells_raw",
        "wetside_allcells_adjusted_M1",
        "categories_allcells_adjusted_M1",
        "wetside_clean_raw",
        "wetside_clean_adjusted_M1",
        "wetside_raw",
        "wetside_adjusted_M1",
        "categories_raw",
        "categories_adjusted_M1",
    ):
        r = res[key]
        print(f"\n== {key} ({r['n_aoi']} AOIs, {r['n_groups']} location groups)")
        for nm, v in r["terms"].items():
            lo, hi = v["ci95_odds_ratio"]
            print(f"  {nm:22s} OR {v['odds_ratio']:6.3f} [{lo:6.3f}, {hi:6.3f}]")
    print("\ncells by class and band:", res["cells_by_class_and_band"])


if __name__ == "__main__":
    main()
