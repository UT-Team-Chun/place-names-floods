# %% v50（事後・amendment B-6 の 4）: GFM のラベルの診断と、季節水域を戻した変種 V1
"""確証判定は変えない。B-6 の 4 で計算の前に固定した量を計算する。

- 層の成分（gfm_v1.parquet が要る）: EMS で評価可能なセルのうち GFM で評価できない割合（全体・水系側の名の集落から
  500 m 以内・他の集落から 500 m 以内）と、そのセルの成分（除外域・恒常水域・季節水域・データなし）の平均、EMS の
  浸水率。除外域が過半のセル（c_excl ≥ 0.5）の HAND（セル平均）≥ 10 m の割合と、CCI の市街地・樹木の割合（≥ 0.5 の
  割合と平均）。比較に GFM で評価できるセルの同じ量
- 変種 V1（gfm_v1.parquet が要る）: 季節水域を判定可能に戻したラベルで、B-3 の B（水系側の群のオッズ比）と C1
  （採点し直し）。同じセルの EMS も
- AOI の選択とセルの選択の分離: B-3 で採用した 192 AOI の EMS で評価可能な全セルでの EMS のオッズ比
- 寄与する単位の数: 500 m 以内の水系側の群で浸水セルを持つ AOI・場所の数と浸水セルの数
- 区間の頑健性: 500 m 以内の水系側の群のオッズ比（GFM・同じセルの EMS・全 AOI の EMS）の、場所グループのペアの
  クラスタ bootstrap（1,000 回、回帰を当て直す。重複して引いた場所は別の場所・別の AOI として扱う）の区間と、場所を
  1 つずつ除いた範囲
- C2（GFM で学習し直し、world_v2_gfm）の診断: 定数の予測（fold ごとの訓練側の浸水率、評価側の浸水率）と各腕の
  Brier、M1s → M3d の増分（地域構成の項を使わない比較）
- 一致の指標: CSI と条件つきの率（全体と事後画像の種類別）、Sentinel-1 を含まない AOI の採点し直し（v41 の値の再掲）

出力: stdout と docs/verification/v50_gfm_diagnostics.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR, RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v50_gfm_diagnostics.json"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
LABELS = DATA_DIR / "global" / "labels"
SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
C2 = RUNS_DIR / "world_v2_gfm" / "predictions" / "primary"
SEED = 20261028  # B-6
B_PAIRS = 1000
B = 4000
TERM = "wetside:le500m"
EXCL_MAJORITY = 0.5  # 除外域が過半のセル（計算の時点で決めた操作的な定義。結果を見る前）
HAND_HIGH_M = 10.0
CCI_HIGH = 0.5


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v39 = _load("v39_within_aoi_odds")
v41 = _load("v41_gfm_independent_labels")


def base_table() -> pd.DataFrame:
    t = v41.build()
    s = pd.read_parquet(V37T, columns=["aoi", "y", "hand_mean", "cci_built", "cci_tree"])
    if not (np.array_equal(t["aoi"].to_numpy(), s["aoi"].to_numpy())
            and np.array_equal(t["y"].to_numpy(), s["y"].to_numpy())):
        msg = "v37 and v41 tables are not aligned"
        raise SystemExit(msg)
    for c in ("hand_mean", "cci_built", "cci_tree"):
        t[c] = s[c].to_numpy()
    v41.attach_lexfix(t)
    return t


def attach_v1(t: pd.DataFrame, summary: dict) -> pd.DataFrame | None:
    """B-3 で採用した AOI のセルに gfm_v1.parquet（B-3 のラベルの作り直し・V1・層の成分）を結ぶ。"""
    keep = sorted(a for a, s in summary["aois"].items() if s.get("include"))
    paths = [LABELS / a / "gfm_v1.parquet" for a in keep]
    if not all(p.is_file() for p in paths):
        return None
    g = pd.concat([pd.read_parquet(p).assign(aoi=a) for a, p in zip(keep, paths, strict=True)], ignore_index=True)
    sub = t[t["aoi"].isin(keep)].reset_index(drop=True)
    m = sub.merge(g, on=["aoi", "ix", "iy"], how="left", validate="many_to_one")
    if len(m) != len(sub) or m["gfm_eval"].isna().any():
        msg = "gfm_v1 does not cover the analysed cells"
        raise SystemExit(msg)
    m["gfm_eval"] = m["gfm_eval"].astype(bool)
    m["v1_eval"] = m["v1_eval"].astype(bool)
    m["sensor"] = m["aoi"].map({a: s["sensor"] for a, s in summary["aois"].items()})
    return m


def reproduction(summary: dict) -> dict:
    """作り直した B-3 のラベル（gfm_v1.parquet の gfm_eval・gfm_y）と B-3 のラベル（gfm.parquet）の、AOI の全セルでの一致。"""
    keep = sorted(a for a, s in summary["aois"].items() if s.get("include"))
    n_same = n_cells = n_eval_diff = n_y_diff = n_skipped = 0
    for a in keep:
        d = LABELS / a
        meta = json.loads((d / "gfm_v1_meta.json").read_text(encoding="utf-8"))
        n_same += bool(meta.get("b3_labels_reproduced"))
        n_skipped += sum(len(sc.get("skipped_unreadable", [])) for sc in meta.get("scenes", []))
        old = pd.read_parquet(d / "gfm.parquet", columns=["gfm_eval", "gfm_y"])
        new = pd.read_parquet(d / "gfm_v1.parquet", columns=["gfm_eval", "gfm_y"])
        eo, en = old["gfm_eval"].to_numpy(bool), new["gfm_eval"].to_numpy(bool)
        n_cells += len(old)
        n_eval_diff += int((eo != en).sum())
        n_y_diff += int((old["gfm_y"].to_numpy()[eo & en] != new["gfm_y"].to_numpy()[eo & en]).sum())
    return {"n_aoi": len(keep), "n_aoi_identical": n_same, "n_cells_all": n_cells, "n_eval_differs": n_eval_diff,
            "n_label_differs": n_y_diff, "n_items_unreadable_on_rebuild": n_skipped}


# %% 層の成分
def components(m: pd.DataFrame) -> dict:
    d1 = m["P_d1"].to_numpy()
    wet = m["lex_wetside_n1"].to_numpy() > 0
    ev = m["gfm_eval"].to_numpy()
    comp = [f"c_{k}" for k in ("excl", "perm", "seas", "nodata")]
    out = {}
    for key, sel in (("all", np.ones(len(m), dtype=bool)), ("wetside_le500m", wet & (d1 <= 500)),
                     ("other_le500m", ~wet & (d1 <= 500))):
        s = m.loc[sel]
        ne = s.loc[~s["gfm_eval"]]
        out[key] = {"n_cells": int(sel.sum()), "share_not_gfm_evaluable": float((~s["gfm_eval"]).mean()),
                    "mean_components_not_evaluable": {c: float(ne[c].mean()) for c in comp},
                    "mean_components_evaluable": {c: float(s.loc[s["gfm_eval"], c].mean()) for c in comp},
                    "ems_rate_not_evaluable": float(ne["y"].mean()) if len(ne) else None,
                    "ems_rate_evaluable": float(s.loc[s["gfm_eval"], "y"].mean())}
    exc = (~ev) & (m["c_excl"].to_numpy() >= EXCL_MAJORITY)

    def prof(sel: np.ndarray) -> dict:
        s = m.loc[sel]
        return {"n_cells": int(sel.sum()), "share_hand_ge10m": float((s["hand_mean"] >= HAND_HIGH_M).mean()),
                "share_built_ge05": float((s["cci_built"] >= CCI_HIGH).mean()),
                "share_tree_ge05": float((s["cci_tree"] >= CCI_HIGH).mean()),
                "mean_built": float(s["cci_built"].mean()), "mean_tree": float(s["cci_tree"].mean()),
                "median_hand_m": float(s["hand_mean"].median()), "ems_rate": float(s["y"].mean())}

    out["excluded_majority"] = prof(exc)
    out["gfm_evaluable"] = prof(ev)
    out["share_not_evaluable_with_excl_majority"] = float(exc.sum() / max((~ev).sum(), 1))
    return out


# %% 寄与する単位の数
def contributing(t: pd.DataFrame, ycol: str) -> dict:
    g = (t["lex_wetside_n1"].to_numpy() > 0) & (t["P_d1"].to_numpy() <= 500)
    y = t[ycol].to_numpy()
    f = g & (y == 1)
    return {"n_cells_group": int(g.sum()), "n_flooded_group": int(f.sum()),
            "n_aoi_group": int(t.loc[g, "aoi"].nunique()), "n_aoi_with_flooded": int(t.loc[f, "aoi"].nunique()),
            "n_units_group": int(t.loc[g, "unit"].nunique()), "n_units_with_flooded": int(t.loc[f, "unit"].nunique())}


# %% 区間の頑健性
def _prep(t: pd.DataFrame, ycol: str) -> tuple:
    y = t[ycol].to_numpy(dtype=np.float64)
    rate = t.assign(_y=y).groupby("aoi")["_y"].transform("mean").to_numpy()
    m = (rate > 0) & (rate < 1)
    s = t.loc[m].reset_index(drop=True)
    z, names = v39.design(s, ("wetside",), True, banded=True)
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    return y[m], z, names.index(TERM), aoi, unit, unit_lab


def robustness(t: pd.DataFrame, ycol: str, key: str, reps: int = B_PAIRS) -> dict:
    """重い（全 AOI の EMS で約 30 分）ので、結果を data/verification/v50_robust_<key>.json に残して使い回す。"""
    cache = DATA_DIR / "verification" / f"v50_robust_{key}.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))
    out = _robustness(t, ycol, reps)
    cache.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def _robustness(t: pd.DataFrame, ycol: str, reps: int) -> dict:
    y, z, j, aoi, unit, unit_lab = _prep(t, ycol)
    g, na = len(unit_lab), int(aoi.max()) + 1
    rows = [np.flatnonzero(unit == u) for u in range(g)]
    point = v39.fit(y, z, aoi, unit)
    b0 = float(point["beta"][j])
    rng = np.random.default_rng(SEED)
    beta = np.full(reps, np.nan)
    t0 = time.time()
    for r in range(reps):
        sel = rng.integers(0, g, size=g)
        idx = np.concatenate([rows[u] for u in sel])
        pos = np.concatenate([np.full(len(rows[u]), k, dtype=np.int64) for k, u in enumerate(sel)])
        try:  # 重複して引いた場所は別の場所・別の AOI として当てる
            beta[r] = v39.fit(y[idx], z[idx], aoi[idx] + pos * na, pos)["beta"][j]
        except np.linalg.LinAlgError:
            continue
    loo = np.full(g, np.nan)
    for u in range(g):
        keep = unit != u
        try:
            loo[u] = v39.fit(y[keep], z[keep], aoi[keep], unit[keep])["beta"][j]
        except np.linalg.LinAlgError:
            continue
    ok = np.isfinite(beta)
    lo, hi = np.nanargmin(loo), np.nanargmax(loo)
    return {"odds_ratio": float(np.exp(b0)),
            "sandwich_ci95": [float(np.exp(b0 - 1.96 * point["se"][j])), float(np.exp(b0 + 1.96 * point["se"][j]))],
            "pairs_bootstrap": {"reps": reps, "n_ok": int(ok.sum()),
                                "ci95_odds_ratio": [float(np.exp(np.quantile(beta[ok], 0.025))),
                                                    float(np.exp(np.quantile(beta[ok], 0.975)))],
                                "ci90_odds_ratio": [float(np.exp(np.quantile(beta[ok], 0.05))),
                                                    float(np.exp(np.quantile(beta[ok], 0.95)))],
                                "share_beta_below_minus5": float((beta[ok] < -5).mean()),
                                "sec": round(time.time() - t0)},
            "leave_one_place_out": {"n_units": g, "min_odds_ratio": float(np.exp(loo[lo])),
                                    "min_drop_unit": str(unit_lab[lo]), "max_odds_ratio": float(np.exp(loo[hi])),
                                    "max_drop_unit": str(unit_lab[hi])}}


# %% C2 の診断
def c2_diag() -> dict:
    files = sorted(C2.glob("fold[0-9].npz"), key=lambda f: int(f.stem[4:]))
    if len(files) < 10:
        return {}
    parts = [np.load(f) for f in files]
    meta = [json.loads(f.with_suffix(".json").read_text(encoding="utf-8")) for f in files]
    y = np.concatenate([q["y"] for q in parts]).astype(np.float64)
    fold = np.concatenate([np.full(len(q["y"]), k) for k, q in enumerate(parts)])
    unit = np.concatenate([q["unit"] for q in parts])
    n_all = len(y)
    train_rate = np.array([y[fold != k].mean() for k in range(len(parts))])
    train_is_rest = all(mt["n_train"] == n_all - mt["n_test"] for mt in meta)
    arms = ("M0", "M0s", "M1", "M1s", "M1c", "M3n", "M3s", "M3d")
    pred = {a: np.concatenate([q[f"p_{a}"] for q in parts]) for a in arms}
    pred["const_train_rate"] = train_rate[fold]
    pred["const_eval_rate_pooled"] = np.full(n_all, y.mean())
    pred["const_eval_rate_by_fold"] = np.array([y[fold == k].mean() for k in range(len(parts))])[fold]
    uniq, inv = np.unique(unit, return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    sse = {a: np.bincount(inv, weights=(p - y) ** 2) for a, p in pred.items()}
    brier = {a: float(v.sum() / n.sum()) for a, v in sse.items()}
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(n), np.full(len(n), 1 / len(n)), size=B).astype(np.float64)
    d = (w @ sse["M1s"] - w @ sse["M3d"]) / (w @ n)
    return {"n_cells": n_all, "n_units": len(uniq), "wet_rate": float(y.mean()),
            "train_rate_by_fold": [float(x) for x in train_rate], "train_is_all_other_folds": train_is_rest,
            "brier": brier, "delta_bs_M1s_M3d": brier["M1s"] - brier["M3d"],
            "ci95_delta_M1s_M3d": [float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))],
            "ci90_delta_M1s_M3d": [float(np.quantile(d, 0.05)), float(np.quantile(d, 0.95))],
            "yardstick_M0s_M1s": brier["M0s"] - brier["M1s"], "yardstick_M0_M1": brier["M0"] - brier["M1"]}


# %% 一致の指標
def agreement(t: pd.DataFrame, gcol: str = "y_gfm") -> dict:
    e, g = t["y_ems"].to_numpy() == 1, t[gcol].to_numpy() == 1
    both, either = int((e & g).sum()), int((e | g).sum())
    return {"n_cells": len(t), "csi": both / either if either else None,
            "p_gfm1_given_ems0": float(g[~e].mean()), "p_ems1_given_gfm0": float(e[~g].mean()),
            "p_gfm1_given_ems1": float(g[e].mean()), "p_ems1_given_gfm1": float(e[g].mean())}


def main() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    t0 = base_table()
    t = v41.attach_gfm(t0, summary)  # GFM で評価できるセル（B-3 の解析のセル）
    keep = {a for a, s in summary["aois"].items() if s.get("include")}
    t192 = t0[t0["aoi"].isin(keep)].reset_index(drop=True)
    res: dict = {"rule": "amendment B-6 item 4", "seed": SEED, "B_pairs": B_PAIRS, "B": B,
                 "definitions": {"excluded_majority": f"not GFM-evaluable and c_excl >= {EXCL_MAJORITY}",
                                 "hand_high_m": HAND_HIGH_M, "cci_high": CCI_HIGH}}
    res["agreement"] = {"all": agreement(t)} | {f"sensor_{k}": agreement(s) for k, s in t.groupby("sensor")}
    v41_json = json.loads((ROOT / "docs" / "verification" / "v41_gfm_independent_labels.json").read_text())
    res["non_s1_rescoring_v41"] = {k: v41_json["sensitivity_independent_sensor"][k] for k in
                                   ("n_aoi", "C1_gfm", "C1_ems")}
    res["aoi_vs_cell_selection"] = {"ems_192aoi_all_cells": v41.odds(t192, "y")["wetside_adjusted_M1"]["terms"][TERM]
                                    | {"n_cells": len(t192), "n_aoi": int(t192["aoi"].nunique())},
                                    "ems_gfm_cells": v41.odds(t, "y_ems")["wetside_adjusted_M1"]["terms"][TERM],
                                    "gfm_gfm_cells": v41.odds(t, "y_gfm")["wetside_adjusted_M1"]["terms"][TERM]}
    res["contributing"] = {"gfm_gfm_cells": contributing(t, "y_gfm"), "ems_gfm_cells": contributing(t, "y_ems"),
                           "ems_192aoi_all_cells": contributing(t192, "y"), "ems_all_aoi": contributing(t0, "y")}
    res["c2"] = c2_diag()
    m = attach_v1(t0, summary)
    if m is not None:
        rep = m.loc[m["gfm_eval"]].reset_index(drop=True)
        same = t.merge(rep[["aoi", "ix", "iy", "gfm_y"]], on=["aoi", "ix", "iy"], how="left")
        res["v1_reproduces_b3"] = {"n_cells_b3": len(t), "n_cells_rebuilt": len(rep),
                                   "labels_equal": bool(len(rep) == len(t)
                                                        and (same["gfm_y"].to_numpy() == t["y_gfm"].to_numpy()).all())}
        res["v1_reproduces_b3"] |= reproduction(summary)
        res["components"] = components(m)
        tv = m.loc[m["v1_eval"]].reset_index(drop=True)
        tv["y_ems"] = tv["y"].to_numpy(dtype=np.int8)
        tv["y_v1"] = tv["v1_y"].to_numpy(dtype=np.int8)
        res["v1"] = {"n_cells": len(tv), "n_aoi": int(tv["aoi"].nunique()), "rate_v1": float(tv["y_v1"].mean()),
                     "agreement": agreement(tv, "y_v1"),
                     "odds_v1": v41.odds(tv, "y_v1"), "odds_ems": v41.odds(tv, "y_ems"),
                     "C1_v1": v41.scored(tv, "y_v1"), "C1_ems": v41.scored(tv, "y_ems"),
                     "contributing_v1": contributing(tv, "y_v1"), "contributing_ems": contributing(tv, "y_ems")}
        if "p_lf_M1c" in tv:
            res["v1"]["C1_v1_lexfix"] = v41.scored(tv, "y_v1", prefix="p_lf_")
        res["robustness_v1"] = {"v1_v1_cells": robustness(tv, "y_v1", "v1_v1_cells"),
                                "ems_v1_cells": robustness(tv, "y_ems", "ems_v1_cells")}
    res["robustness"] = {"gfm_gfm_cells": robustness(t, "y_gfm", "gfm_gfm_cells"),
                         "ems_gfm_cells": robustness(t, "y_ems", "ems_gfm_cells"),
                         "ems_all_aoi": robustness(t0, "y", "ems_all_aoi")}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    show(res)


def show(res: dict) -> None:
    print("agreement", {k: round(v["csi"], 3) for k, v in res["agreement"].items()})
    for k, v in res["aoi_vs_cell_selection"].items():
        print(k, round(v["odds_ratio"], 3), [round(x, 3) for x in v["ci95_odds_ratio"]])
    for k, v in res["contributing"].items():
        print(k, v)
    c = res["c2"]
    if c:
        print("C2 brier", {k: round(v, 5) for k, v in c["brier"].items()})
        print("C2 M1s->M3d", c["delta_bs_M1s_M3d"], c["ci95_delta_M1s_M3d"])
    for k, v in res["robustness"].items():
        print(k, round(v["odds_ratio"], 3), "pairs95", [round(x, 3) for x in v["pairs_bootstrap"]["ci95_odds_ratio"]],
              "LOO", round(v["leave_one_place_out"]["min_odds_ratio"], 3),
              round(v["leave_one_place_out"]["max_odds_ratio"], 3))
    if "v1" in res:
        o = res["v1"]["odds_v1"]["wetside_adjusted_M1"]["terms"][TERM]
        e = res["v1"]["odds_ems"]["wetside_adjusted_M1"]["terms"][TERM]
        print("V1 OR", round(o["odds_ratio"], 3), o["ci95_odds_ratio"], "EMS same cells", round(e["odds_ratio"], 3))
        print("V1 C1", res["v1"]["C1_v1"].get("R"), res["v1"]["C1_v1"].get("decision"))


if __name__ == "__main__":
    main()
