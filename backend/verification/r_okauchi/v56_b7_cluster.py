# %% v56（事後・amendment B-7 の 1〜6(b)）: クラスタで回した B-7 の実行の集計
"""確証判定は変えない。B-7 で計算の前に固定した比較を、各実行の held-out 予測から計算する。

比較はどれも確証の集計（world_aggregate）と同じ固定余白の形（場所グループのクラスタ bootstrap 4,000 回、
シード 20261025、R = ΔBS / ΔBS(M0 → M1)、§6 の判定表）。

- 1（world_v2_b7perm）: 地域構成を置き換えた名前から計算し直す置換（基準 M1 + c(q) → 名前 M1 + c(q) + q）と
  AOI 内で中心化した置換（M1s → M1s + q − c(q)）の分布と実の増分（M1c → M3n、M1s → M3d）の百分位。B-5 の置換の
  分布（world_v2_b5perm）での確証の placebo（M3np）の百分位も
- 2（world_v2_b7plasref）: plasmode の形 A を注入したラベルで基礎モデルから当て直した利得と、基礎を固定した利得
  （同じ δ・同じ回）、その比、当て直した曲線から読む余白のオッズ比 δ_M と上限 δ_A（v52 と同じ補間）
- 3（world_v2_b7actfix）: activation で束ねた分割を直したコードで（凍結したコードの world_v2_actgrp と並べる）
- 4（world_v2_b7recur）: 繰り返す浸水の対照 (a) セルに直接 M1cEc → M3oEc、(b) 集落を通じて M1cE → M3oE、
  別の事象のある AOI に限った R（その部分集合の物差し・場所の bootstrap）
- 5（world_v2_b7fp）: 氾濫原の地面を指す読み M1cF → M3nF（placebo M3nFp）
- 6(b)（world_v2_b7partial）: 部分的に浸水した y = 0 のセルを除いた確証の比較

出力: stdout と docs/verification/v56_b7_cluster.json（無い入力は "missing"）
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v56_b7_cluster.json"
V47 = ROOT / "docs" / "verification" / "v47_b5_cluster.json"
V52 = ROOT / "docs" / "verification" / "v52_b6_cluster.json"
LEXFIX = RUNS_DIR / "world_v2_lexfix"
N_PERMS = 200
OR_GRID_MAX = float(np.exp(1.5))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v52 = _load("v52_b6_cluster")

# 実行 → [(key, 基準, 名前あり, placebo)]
RUNS = {
    "actfix": ("world_v2_b7actfix", [("primary", "M1c", "M3n", "M3np"), ("total", "M1s", "M3s", "M3ps")]),
    "recur": ("world_v2_b7recur", [("recur_cell", "M1cEc", "M3oEc", None), ("recur_settlement", "M1cE", "M3oE", None),
                                   ("primary", "M1c", "M3n", "M3np")]),
    "floodplain": ("world_v2_b7fp", [("floodplain", "M1cF", "M3nF", "M3nFp"), ("primary", "M1c", "M3n", "M3np")]),
    "partial": ("world_v2_b7partial", [("primary", "M1c", "M3n", "M3np")]),
}
SAME_AS_LEXFIX = {"recur", "floodplain"}  # 同じ表・種 → M1c・M3n が直した実行と一致するはず


def _complete(run: Path, pattern: str = "fold[0-9].npz") -> bool:
    d = run / "predictions" / "primary"
    return d.is_dir() and len(list(d.glob(pattern))) >= 10


def _verdicts(y: np.ndarray, preds: dict[str, np.ndarray], unit: np.ndarray, pairs: list) -> dict:
    s, n, arms = wa.unit_sums(y, preds, unit)
    idx = {a: j for j, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    delta = wa.make_delta(hat, wa.bootstrap_bs(s, n), idx)
    d01, b01 = delta("M0", "M1")
    out: dict = {"n_cells": len(y), "n_units": len(n), "wet_rate": float(y.mean()),
                 "brier": {a: float(hat[idx[a]]) for a in arms}, "yardstick": d01, "ci95_yardstick": wa.ci(b01, 0.95)}
    for name, base, alt, plac in pairs:
        if base in idx and alt in idx:
            v, p = wa.verdict(delta, base, alt, plac if plac in idx else None, d01=d01, b01=b01)
            out[name] = v | ({"placebo": p} if p else {})
    return out


def compare(key: str) -> dict | str:
    tag, pairs = RUNS[key]
    run = RUNS_DIR / tag
    if not _complete(run):
        return "missing"
    data = wa.load_split(run, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    out = {"run": tag, "n_aoi": len(np.unique(data["aoi"]))} | _verdicts(y, preds, data["unit"], pairs)
    rec = json.loads((run / "predictions" / "primary" / "fold0.json").read_text(encoding="utf-8"))
    if rec.get("b6"):
        out["record"] = rec["b6"]
    if key in SAME_AS_LEXFIX and _complete(LEXFIX):
        ref = wa.load_split(LEXFIX, "primary")
        same_y = len(ref["y"]) == len(y) and bool(np.array_equal(ref["y"], data["y"]))
        out["check_vs_lexfix"] = {"same_cells": same_y} | (
            {f"max_abs_diff_{a}": float(np.max(np.abs(ref[f"p_{a}"].astype(float) - preds[a].astype(float))))
             for a in ("M1c", "M3n") if a in preds} if same_y else {})
    if key == "recur" and rec.get("b6", {}).get("recurrence", {}).get("aois_with_other_events"):
        keep = np.isin(data["aoi"], rec["b6"]["recurrence"]["aois_with_other_events"])
        sub = {a: v[keep] for a, v in preds.items()}
        out["subset_aois_with_other_events"] = {"n_aoi": len(np.unique(data["aoi"][keep]))} | _verdicts(
            y[keep], sub, data["unit"][keep], pairs[:2])
    if key == "actfix":
        out["frozen_code_run"] = actgrp_frozen()
    return out


def actgrp_frozen() -> dict | str:
    """3: 凍結したコードの activation で束ねた分割（B-2）の同じ比較（Brier の変化を記述するため）。"""
    run = RUNS_DIR / "world_v2_actgrp"
    if not _complete(run):
        return "missing"
    data = wa.load_split(run, "primary")
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    return _verdicts(data["y"].astype(np.float64), preds, data["unit"], RUNS["actfix"][1])


# %% 1. 地域構成を計算し直した置換
def permutations() -> dict | str:
    d = RUNS_DIR / "world_v2_b7perm" / "predictions" / "primary"
    files = sorted(d.glob("fold[0-9]_permsrc_*.npz")) if d.is_dir() else []
    met = json.loads((LEXFIX / "metrics_primary.json").read_text(encoding="utf-8"))
    yard = float(met["primary"]["delta_bs_M0_M1"])
    folds: dict[int, list] = {}
    for f in files:
        folds.setdefault(int(f.name[4]), []).append(f)
    have = {k: len(v) for k, v in folds.items()}
    if len(folds) < 10:
        return f"missing (folds {sorted(folds)})"
    real = {a: [] for a in ("M1c", "M3n", "M3np", "M1s", "M3d")}
    sse = {k: [] for k in ("ref", "names", "cent")}
    n_unit = []
    for k in sorted(folds):
        parts = [np.load(f) for f in sorted(folds[k])]  # 種の範囲の順（ファイル名の先頭の種）
        units = parts[0]["units"]
        if not all(np.array_equal(p["units"], units) for p in parts):
            msg = f"fold {k}: unit order differs between permutation ranges"
            raise SystemExit(msg)
        n_unit.append(parts[0]["n_unit"].astype(np.float64))
        for a, acc in real.items():
            acc.append(parts[0][f"real_{a}"])
        for s, acc in sse.items():
            acc.append(np.concatenate([p[f"sse_{s}"] for p in parts], axis=0))
    n_perm = {s: {v.shape[0] for v in acc} for s, acc in sse.items()}
    if any(len(v) != 1 or next(iter(v)) != N_PERMS for v in n_perm.values()):
        return f"incomplete ({have})"
    n = float(np.concatenate(n_unit).sum())
    tot = {a: float(np.concatenate(v).sum()) for a, v in real.items()}
    ssum = {s: np.concatenate(v, axis=1).sum(axis=1) for s, v in sse.items()}  # [置換]
    d_recomp = (ssum["ref"] - ssum["names"]) / n
    d_real = (tot["M1c"] - tot["M3n"]) / n
    d_cent = (tot["M1s"] - ssum["cent"]) / n
    d_real_c = (tot["M1s"] - tot["M3d"]) / n

    def dist(d_perm: np.ndarray, d_r: float) -> dict:
        r = d_perm / yard
        return {"R_real": d_r / yard, "R_perm_mean": float(r.mean()), "R_perm_q025": float(np.quantile(r, 0.025)),
                "R_perm_q975": float(np.quantile(r, 0.975)), "share_perm_ge_real": float(np.mean(d_perm >= d_r)),
                "percentile_of_real": float(100 * np.mean(d_perm < d_r))}

    out = {"n_perms": N_PERMS, "n_cells": int(n), "yardstick_lexfix": yard,
           "check_real_brier_vs_lexfix": {a: [tot[a] / n, float(met["brier"][a])] for a in ("M1c", "M3n")},
           "recomputed_composition": dist(d_recomp, d_real), "centred": dist(d_cent, d_real_c),
           "reference_brier_perm_mean": float((ssum["ref"] / n).mean()), "brier_M1c_real": tot["M1c"] / n}
    out["b5_single_placebo"] = b5_placebo_percentile(yard)
    return out


def b5_placebo_percentile(yard: float) -> dict | str:
    """1: 確証の placebo（M3np）の、B-5 の置換の分布（world_v2_b5perm）での百分位。"""
    d = RUNS_DIR / "world_v2_b5perm" / "predictions" / "primary"
    files = sorted(d.glob("fold*_perms.npz")) if d.is_dir() else []
    if len(files) < 10:
        return "missing"
    y, m1c, m3np, perms = [], [], [], []
    for f in files:
        with np.load(f) as z:
            y.append(z["y"].astype(np.float64))
            m1c.append(z["p_M1c"].astype(np.float64))
            m3np.append(z["p_M3np"].astype(np.float64))
            perms.append(z["p_perms"].astype(np.float64))
    y, m1c, m3np = np.concatenate(y), np.concatenate(m1c), np.concatenate(m3np)
    pp = np.concatenate(perms, axis=1)
    bs_c = float(np.mean((m1c - y) ** 2))
    d_plac = bs_c - float(np.mean((m3np - y) ** 2))
    d_perm = bs_c - np.mean((pp - y[None, :]) ** 2, axis=1)
    return {"R_placebo": d_plac / yard, "percentile_of_single_placebo": float(100 * np.mean(d_perm < d_plac)),
            "R_perm_mean": float(d_perm.mean() / yard)}


# %% 2. 基礎モデルも当て直した plasmode
def plasmode_refit() -> dict | str:
    d = RUNS_DIR / "world_v2_b7plasref" / "predictions" / "primary"
    files = sorted(d.glob("fold[0-9]_plasref.npz")) if d.is_dir() else []
    if len(files) < 10:
        return f"missing ({len(files)} of 10 folds)"
    lex = json.loads((LEXFIX / "metrics_primary.json").read_text(encoding="utf-8"))["primary"]
    p47 = json.loads(V47.read_text(encoding="utf-8"))["permutations"]
    yard = float(lex["delta_bs_M0_M1"])
    u90 = float(lex["ci90_delta"][1])
    d_null = float(p47["R_perm_mean"]) * float(p47["yardstick"])
    margin = wa.SESOI * yard
    zs = [np.load(f) for f in files]
    models = [str(x) for x in zs[0]["models"]]
    mi = {m: j for j, m in enumerate(models)}
    deltas = np.asarray(zs[0]["deltas"], dtype=np.float64)
    n = float(sum(int(z["n_te"]) for z in zs))
    base = sum(z["base"] for z in zs) / n  # [model]: 注入しない実行（その実行の基礎）
    refit = sum(z["refit"] for z in zs) / n  # [δ, 回, model]
    fixed = sum(z["fixed"] for z in zs) / n  # [δ, 回, (M1c, M3n)]
    flips = sum(z["flips"] for z in zs) / n
    d0 = float(base[mi["M1c"]] - base[mi["M3n"]])
    inc_refit = refit[:, :, mi["M1c"]] - refit[:, :, mi["M3n"]]
    inc_fixed = fixed[:, :, 0] - fixed[:, :, 1]
    gain_refit, gain_fixed = inc_refit.mean(axis=1) - d0, inc_fixed.mean(axis=1) - d0
    grid = np.concatenate([[0.0], deltas])
    g = np.concatenate([[0.0], gain_refit])
    out: dict = {"n_cells": int(n), "deltas": deltas.tolist(), "odds_ratios": np.exp(deltas).tolist(),
                 "reps": int(refit.shape[1]), "delta_bs_0": d0, "margin": margin, "U90": u90,
                 "delta_bs_null_perm_mean": d_null, "yardstick_lexfix": yard,
                 "gain_refit": gain_refit.tolist(), "gain_fixed": gain_fixed.tolist(),
                 "ratio_refit_to_fixed": (gain_refit / gain_fixed).tolist(),
                 "gain_refit_by_rep": (inc_refit - d0).tolist(), "gain_fixed_by_rep": (inc_fixed - d0).tolist(),
                 "yardstick_injected_refit": (refit[:, :, mi["M0"]] - refit[:, :, mi["M1"]]).mean(axis=1).tolist(),
                 "label_flip_share": flips.mean(axis=1).tolist(), "monotone_refit": bool(np.all(np.diff(g) >= 0))}
    for name, tgt in (("delta_M", margin), ("delta_M_prime", margin - d_null)):
        x, how = v52._first_reach(g, tgt, grid)
        out[name] = {"delta": x, "odds_ratio": float(np.exp(x)) if x is not None else None, "how": how,
                     "odds_ratio_text": f"{np.exp(x):.2f}" if x is not None else f"> {OR_GRID_MAX:.2f}"}
    for name, f in (("delta_A", d_null + g - u90), ("delta_B", g - u90)):
        x, how = v52._cross(f, grid)
        out[name] = {"delta": x, "odds_ratio": float(np.exp(x)) if x is not None else None, "how": how,
                     "odds_ratio_text": (f"{np.exp(x):.2f}" if x is not None else
                                         f"> {OR_GRID_MAX:.2f}" if how == "above grid" else "undefined")}
    if V52.is_file():  # B-6 の基礎を固定した plasmode（形 A、10 回）の同じ δ の利得と、ここでの回 0〜2 の一致
        p52 = json.loads(V52.read_text(encoding="utf-8")).get("plasmode")
        if isinstance(p52, dict):
            d52 = list(p52["deltas"])
            out["b6_fixed_gain_10reps"] = [p52["forms"]["A"]["gain"][d52.index(float(x))] for x in deltas]
            out["b6_delta_M"] = p52["forms"]["A"]["delta_M"]
    return out


def main() -> None:
    res: dict = {"rule": "amendment B-7 items 1, 2, 3, 4, 5, 6(b)"}
    for key in RUNS:
        res[key] = compare(key)
    res["permutations"] = permutations()
    res["plasmode_refit"] = plasmode_refit()
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for key in RUNS:
        v = res[key]
        if not isinstance(v, dict):
            print(f"{key:11s} {v}")
            continue
        for name, *_ in RUNS[key][1]:
            c = v.get(name)
            if c:
                print(f"{key:11s} {name:17s} R {c['R']:+.4f} 90% [{c['ci90_R'][0]:+.4f}, {c['ci90_R'][1]:+.4f}] "
                      f"{c['decision']:22s} cells {v['n_cells']:,} AOIs {v['n_aoi']}")
        for extra in ("check_vs_lexfix", "subset_aois_with_other_events"):
            if extra in v:
                print(f"{'':11s} {extra}: {json.dumps(v[extra], default=float)[:300]}")
    print("permutations", json.dumps(res["permutations"], default=float)[:800])
    print("plasmode_refit", json.dumps(res["plasmode_refit"], default=float)[:1200])


if __name__ == "__main__":
    main()
