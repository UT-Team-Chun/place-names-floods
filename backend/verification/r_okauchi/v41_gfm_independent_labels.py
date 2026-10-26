# %% v41（事後・amendment B-3）: 独立の浸水ラベル（Copernicus GFM）による確認 — ラベルの点検・集落の規模・主比較の採点し直し
"""確証判定は変えない。B-3 で計算の前に登録した量を計算する（C2 の学習し直しは world_eval --labels gfm）。

入力: 確証解析の主分割の標本と held-out 予測（world_v2）、GFM のラベル（labels/<AOI>/gfm.parquet、
gfm_summary.json の採用 AOI）。解析のセル = 採用 AOI の標本のうち GFM で評価可能なセル（fold と場所グループは
確証解析のまま）。

- A0: EMS と GFM の浸水率・一致率・κ・条件つきの率（全体と事後画像の種類別）
- A1: AOI 内で ≤ 500 m のセルと > 1 km のセルの [EMS − GFM] の差（v37 の対比、小さい群の大きさの重み）
- A2: AOI 内で最寄りの集落名が水系側のセルとそれ以外の [EMS − GFM] の差（≤ 500 m と全セル）
- A3: v39 の回帰（AOI の切片・距離帯・M1 の held-out 予測のロジット）の > 1 km 対 ≤ 500 m のオッズ比（GFM と EMS）
- B: v39 の回帰の水系側の群と 6 区分のオッズ比（≤ 500 m は距離帯つきの模型の項、全セルは帯なし）。GFM と EMS
- C1: world_v2 の held-out 予測を GFM のラベルで採点（ΔBS(M1c → M3n)、物差し ΔBS(M0 → M1)、R、placebo、
  §6 の判定表）。EMS のラベルでも同じセルで
- 感度（B・C1）: EMS の事後画像に Sentinel-1 を含まない AOI だけ
- 感度（B）: 調整の M1 を C2 で GFM のラベルで学習し直した M1（world_v2_gfm）の held-out 予測に替えて

区間は場所グループのクラスタ bootstrap（B = 4,000、シード 20261027）と、回帰はクラスタ頑健なサンドイッチ分散。
出力: stdout の表と docs/verification/v41_gfm_independent_labels.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa
import world_eval as ev

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.dataset import continent_map
from topsem.world.lexicon import CATEGORIES


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v37 = _load("v37_positive_control")
v39 = _load("v39_within_aoi_odds")

CACHE = DATA_DIR / "verification" / "v41_table.parquet"
OUT = ROOT / "docs" / "verification" / "v41_gfm_independent_labels.json"
PRED = RUNS_DIR / "world_v2" / "predictions" / "primary"
GFM_SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
LABELS = DATA_DIR / "global" / "labels"
B = 4000
SEED = 20261027  # B-3
SESOI = 0.05
ARMS = ("M0", "M1", "M1c", "M3n", "M3np")
v37.SEED = SEED  # 対比の bootstrap も B-3 のシードで
v37.B = B


# %% A. 解析表（確証解析の主分割の標本を fold 順に作り直し、held-out 予測と突き合わせる）
def build() -> pd.DataFrame:
    if CACHE.is_file():
        return pd.read_parquet(CACHE)
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    sens = ev.Sensitivity()
    dirs, _, infos = ev.select_aois(gn, sens)
    df = ev.build_table(dirs, seed=20261025, infos=infos, sens=sens, extra_cols=("ix", "iy"))
    df, _, _ = ev.add_channels(df, gn, None)
    cont = continent_map(ev.ROOT / "geonames" / "countryInfo.txt")
    df["fold"] = ev.folds_for("primary", df, infos, cont, buffer_km=sens.buffer_km).to_numpy()
    keep = ["aoi", "unit", "fold", "y", "ix", "iy", "P_d1", "lex_wetside_n1",
            *(f"lex_{c}_n1" for c in CATEGORIES)]
    t = pd.concat([df[df["fold"] == k][keep] for k in sorted(df["fold"].unique()) if k >= 0],
                  ignore_index=True)
    parts = [np.load(f) for f in sorted(PRED.glob("fold*.npz"), key=lambda p: int(p.stem[4:]))]
    y = np.concatenate([p["y"] for p in parts])
    aoi = np.concatenate([p["aoi"] for p in parts])
    if not (np.array_equal(y, t["y"].to_numpy()) and np.array_equal(aoi, t["aoi"].to_numpy(dtype=str))):
        msg = "rebuilt table does not match the held-out predictions"
        raise SystemExit(msg)
    for a in ARMS:
        t[f"p_{a}"] = np.concatenate([p[f"p_{a}"] for p in parts])
    t["n_eval_aoi"] = t["aoi"].map({a: infos[a].n_eval for a in t["aoi"].unique()})
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(CACHE)
    return t


def attach_lexfix(t: pd.DataFrame) -> bool:
    """B-4: 辞書を直した再実行（world_v2_lexfix）の held-out 予測を、同じ行順で p_lf_* として足す。"""
    pred = RUNS_DIR / "world_v2_lexfix" / "predictions" / "primary"
    files = sorted(pred.glob("fold[0-9].npz"), key=lambda f: int(f.stem[4:])) if pred.is_dir() else []
    if len(files) < 10:
        return False
    parts = [np.load(f) for f in files]
    y = np.concatenate([q["y"] for q in parts])
    aoi = np.concatenate([q["aoi"] for q in parts])
    if not (np.array_equal(y, t["y"].to_numpy()) and np.array_equal(aoi, t["aoi"].to_numpy(dtype=str))):
        msg = "world_v2_lexfix predictions are not aligned with the table"
        raise SystemExit(msg)
    for a in ARMS:
        t[f"p_lf_{a}"] = np.concatenate([q[f"p_{a}"] for q in parts])
    return True


def attach_c2(t: pd.DataFrame) -> bool:
    """B の感度（B-3 で登録）: C2 の GFM で学習し直した M1（world_v2_gfm）の held-out 予測を p_c2_M1 として足す。"""
    pred = RUNS_DIR / "world_v2_gfm" / "predictions" / "primary"
    files = sorted(pred.glob("fold[0-9].npz"), key=lambda f: int(f.stem[4:])) if pred.is_dir() else []
    if len(files) < 10:
        return False
    parts = [np.load(f) for f in files]
    y = np.concatenate([q["y"] for q in parts])
    aoi = np.concatenate([q["aoi"] for q in parts])
    if not (np.array_equal(y, t["y_gfm"].to_numpy()) and np.array_equal(aoi, t["aoi"].to_numpy(dtype=str))):
        msg = "world_v2_gfm predictions are not aligned with the GFM table"
        raise SystemExit(msg)
    t["p_c2_M1"] = np.concatenate([q["p_M1"] for q in parts])
    return True


def attach_gfm(t: pd.DataFrame, summary: dict) -> pd.DataFrame:
    """採用 AOI の GFM のラベルを (aoi, ix, iy) で結び、GFM で評価可能なセルだけを返す（並べ替えない）。"""
    keep = sorted(a for a, s in summary["aois"].items() if s.get("include"))
    g = pd.concat([pd.read_parquet(LABELS / a / "gfm.parquet", columns=["ix", "iy", "gfm_eval", "gfm_y"])
                   .assign(aoi=a) for a in keep], ignore_index=True)
    m = t.merge(g, on=["aoi", "ix", "iy"], how="left", validate="many_to_one")
    ok = m["gfm_eval"].astype("boolean").fillna(value=False).to_numpy(dtype=bool)
    out = t.loc[ok].reset_index(drop=True)
    out["y_ems"] = out["y"].to_numpy(dtype=np.int8)
    out["y_gfm"] = m.loc[ok, "gfm_y"].to_numpy(dtype=np.int8)
    out["sensor"] = out["aoi"].map({a: s["sensor"] for a, s in summary["aois"].items()})
    return out


# %% A0. 一致
def agreement(t: pd.DataFrame) -> dict:
    e, g = t["y_ems"].to_numpy(), t["y_gfm"].to_numpy()
    n = len(e)
    po = float((e == g).mean())
    pe = float(e.mean() * g.mean() + (1 - e.mean()) * (1 - g.mean()))
    return {"n_cells": n, "n_aoi": int(t["aoi"].nunique()), "n_units": int(t["unit"].nunique()),
            "rate_ems": float(e.mean()), "rate_gfm": float(g.mean()), "agreement": po,
            "kappa": (po - pe) / (1 - pe) if pe < 1 else None,
            "p_gfm1_given_ems0": float(g[e == 0].mean()) if (e == 0).any() else None,
            "p_ems1_given_gfm0": float(e[g == 0].mean()) if (g == 0).any() else None,
            "p_gfm1_given_ems1": float(g[e == 1].mean()) if (e == 1).any() else None}


# %% A1・A2. AOI 内の [EMS − GFM] の対比（v37.contrast）
def diff_contrasts(t: pd.DataFrame) -> dict:
    t = t.copy()
    t["ems_minus_gfm"] = t["y_ems"].astype(np.float64) - t["y_gfm"].astype(np.float64)
    d1 = t["P_d1"].to_numpy()
    near, far = d1 <= 500, d1 > 1000
    wet = t["lex_wetside_n1"].to_numpy() > 0
    out = {"A1_near_minus_far": v37.contrast(t, "ems_minus_gfm", near, near | far, "capped"),
           "A2_wetside_le500m": v37.contrast(t, "ems_minus_gfm", wet, near, "capped"),
           "A2_wetside_all": v37.contrast(t, "ems_minus_gfm", wet, np.ones(len(t), dtype=bool), "capped")}
    for lab, col in (("ems", "y_ems"), ("gfm", "y_gfm")):  # 記述: それぞれのラベルでの同じ対比
        t[f"{lab}_y"] = t[col].astype(np.float64)
        out[f"near_minus_far_{lab}"] = v37.contrast(t, f"{lab}_y", near, near | far, "capped")
        out[f"wetside_le500m_{lab}"] = v37.contrast(t, f"{lab}_y", wet, near, "capped")
    return out


# %% A3・B. v39 の回帰（AOI の切片、変動の無い AOI は除く）
def odds(t: pd.DataFrame, ycol: str) -> dict:
    y = t[ycol].to_numpy(dtype=np.float64)
    rate = t.assign(_y=y).groupby("aoi")["_y"].transform("mean").to_numpy()
    m = (rate > 0) & (rate < 1)
    s = t.loc[m].reset_index(drop=True)
    ys = y[m]
    out = {"n_cells": int(m.sum()), "n_aoi": int(s["aoi"].nunique()),
           "aoi_dropped_no_variation": int(t["aoi"].nunique() - s["aoi"].nunique())}
    specs = {"bands_adjusted_M1": ((), True, True), "wetside_adjusted_M1": (("wetside",), True, True),
             "wetside_raw": (("wetside",), False, True),
             "categories_adjusted_M1": (tuple(CATEGORIES), True, True),
             "wetside_allcells_adjusted_M1": (("wetside",), True, False),
             "categories_allcells_adjusted_M1": (tuple(CATEGORIES), True, False)}
    for key, (classes, offset, banded) in specs.items():
        z, names = v39.design(s.assign(p_M1=s["p_M1"]), classes, offset, banded=banded)
        res = v39.fit(ys, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
        out[key] = {"terms": v39.table(res, names), "n_groups": res["n_groups"]}
    return out


# %% C1. held-out 予測の採点し直し
def scored(t: pd.DataFrame, ycol: str, prefix: str = "p_") -> dict:
    y = t[ycol].to_numpy(dtype=np.float64)
    uniq, inv = np.unique(t["unit"].to_numpy(), return_inverse=True)
    n = np.bincount(inv).astype(np.float64)
    s = {a: np.bincount(inv, weights=(t[f"{prefix}{a}"].to_numpy(dtype=np.float64) - y) ** 2) for a in ARMS}
    hat = {a: v.sum() / n.sum() for a, v in s.items()}
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(n), np.full(len(n), 1 / len(n)), size=B).astype(np.float64)
    bb = {a: (w @ v) / (w @ n) for a, v in s.items()}
    dn, y0 = hat["M1c"] - hat["M3n"], hat["M0"] - hat["M1"]
    bd, b0 = bb["M1c"] - bb["M3n"], bb["M0"] - bb["M1"]
    bp = bb["M3np"] - bb["M3n"]  # real − permuted（ΔBS の尺度）
    q = lambda x, a, b: [float(np.quantile(x, a)), float(np.quantile(x, b))]
    out = {"n_cells": len(y), "n_units": len(uniq), "wet_rate": float(y.mean()),
           "brier": {a: float(v) for a, v in hat.items()},
           "delta_bs_names": float(dn), "ci95_delta": q(bd, 0.025, 0.975), "ci90_delta": q(bd, 0.05, 0.95),
           "yardstick_M0_M1": float(y0), "ci95_yardstick": q(b0, 0.025, 0.975),
           "real_minus_permuted_delta_bs": float(hat["M3np"] - hat["M3n"]),
           "ci95_real_minus_permuted": q(bp, 0.025, 0.975)}
    if y0 > 0:
        r95 = [c / y0 for c in out["ci95_delta"]]
        r90 = [c / y0 for c in out["ci90_delta"]]
        rp95 = [c / y0 for c in out["ci95_real_minus_permuted"]]
        out |= {"R": float(dn / y0), "ci95_R": r95, "ci90_R": r90, "margin_delta_bs": SESOI * float(y0),
                "equivalent_TOST": bool(-SESOI < r90[0] and r90[1] < SESOI),
                "decision": wa.decide(float(dn / y0), r95, r90, placebo_ci95=rp95,
                                      den_ci95=out["ci95_yardstick"])}
    else:
        out |= {"R": None, "decision": "undefined-denominator"}
    return out


def main() -> None:
    summary = json.loads(GFM_SUMMARY.read_text(encoding="utf-8"))
    t0 = build()
    lexfix = attach_lexfix(t0)
    t = attach_gfm(t0, summary)
    res: dict = {"rule": "amendment B-3", "seed": SEED, "B": B,
                 "aoi_inclusion": {"n_aoi": summary["n_aoi"], "n_include": summary["n_include"],
                                   "why": summary["why"]},
                 "sensor_counts_included": t.groupby("sensor")["aoi"].nunique().to_dict()}
    res["A0"] = {"all": agreement(t)} | {f"sensor_{k}": agreement(s) for k, s in t.groupby("sensor")}
    res["A1_A2"] = diff_contrasts(t)
    res["odds_gfm"] = odds(t, "y_gfm")
    res["odds_ems"] = odds(t, "y_ems")
    if attach_c2(t):  # B の感度: 調整に C2 の GFM で学習した M1 を使う
        tc = t.assign(p_M1=t["p_c2_M1"])
        res["odds_gfm_c2_M1"] = odds(tc, "y_gfm")
        res["odds_ems_c2_M1"] = odds(tc, "y_ems")
    res["C1_gfm"] = scored(t, "y_gfm")
    res["C1_ems"] = scored(t, "y_ems")
    if lexfix:  # B-4: 辞書を直した再実行の予測でも
        res["C1_gfm_lexfix"] = scored(t, "y_gfm", prefix="p_lf_")
        res["C1_ems_lexfix"] = scored(t, "y_ems", prefix="p_lf_")
    ind = t[t["sensor"].isin(["sar", "optical"])].reset_index(drop=True)  # 感度: Sentinel-1 を含まない
    res["sensitivity_independent_sensor"] = {
        "n_aoi": int(ind["aoi"].nunique()), "odds_gfm": odds(ind, "y_gfm"), "odds_ems": odds(ind, "y_ems"),
        "C1_gfm": scored(ind, "y_gfm"), "C1_ems": scored(ind, "y_ems"), "A0": agreement(ind)}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    show(res)


def show(res: dict) -> None:
    a = res["A0"]["all"]
    print(f"cells {a['n_cells']:,}  AOIs {a['n_aoi']}  places {a['n_units']}  "
          f"EMS {a['rate_ems']:.3f}  GFM {a['rate_gfm']:.3f}  agreement {a['agreement']:.3f}  kappa {a['kappa']:.3f}")
    for k, v in res["A1_A2"].items():
        if v:
            print(f"{k:24s} {100 * v['estimate']:+.2f} pp [{100 * v['ci95'][0]:+.2f}, {100 * v['ci95'][1]:+.2f}] "
                  f"AOIs {v['n_aoi']}")
    for lab in ("gfm", "ems"):
        o = res[f"odds_{lab}"]
        t = o["bands_adjusted_M1"]["terms"]["band_gt1km"]
        w = o["wetside_adjusted_M1"]["terms"]["wetside:le500m"]
        print(f"{lab}: >1 km vs <=500 m OR {t['odds_ratio']:.2f} {t['ci95_odds_ratio']}; "
              f"water-side <=500 m OR {w['odds_ratio']:.2f} {w['ci95_odds_ratio']}")
    for lab in ("gfm", "ems"):
        c = res[f"C1_{lab}"]
        print(f"C1 {lab}: dBS {c['delta_bs_names']:.2e} R {c.get('R')} 90% {c.get('ci90_R')} {c['decision']}")


if __name__ == "__main__":
    main()
