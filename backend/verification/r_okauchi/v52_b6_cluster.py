# %% v52（事後・amendment B-6 の 2・3・5・6）: クラスタで回した B-6 の実行の集計
"""確証判定は変えない。B-6 で計算の前に固定した比較を、各実行の held-out 予測から計算する。

比較はどれも確証の集計（world_aggregate）と同じ固定余白の形（場所グループのクラスタ bootstrap 4,000 回、
シード 20261025、R = ΔBS / ΔBS(M0 → M1)、§6 の判定表）。

- 2（world_v2_b6plasmode）: 注入した効果 δ ごとの増分の利得 gain(δ) = 平均_回 ΔBS(δ) − ΔBS(0) と、上限 δ_A
  （ΔBS_null + gain(δ) ≤ U90）・δ_B（gain(δ) ≤ U90）・U95 の変種・余白のオッズ比 δ_M（gain = 余白）・δ_M′
  （ΔBS_null + gain = 余白）。ΔBS_null は B-5 の 1 の置換 200 回の増分の平均（v47）、U90・U95 は直した実行
  （world_v2_lexfix）の確証の増分の 90%・95% 区間の上限。格子の間は δ（log オッズ比）で線形補間
- 3（world_v2_b6lab, world_v2_b6gsw）: 事象のラベルの割合（M1cL → M3oL）、繰り返す水（M1cG → M3oG）と恒常水域を
  除いた変種（M1cGn → M3oGn）
- 5（world_v2_b6strict）: 厳密な読み（M1cS → M3nS、placebo M3nSp）
- 6: (a) world_v2_b6gn・(b) world_v2_b6t0x（T0xc → T0xn と T0c → T0n）・(c) world_v2_b6nohydro・
  (d) world_v2_b6typed・(e) world_v2_b6nread（M1rc → M3rn）・(f) world_v2_b6actperm（置換 200 回の分布での実の
  増分の順位）

出力: stdout と docs/verification/v52_b6_cluster.json（無い入力は "missing"）
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v52_b6_cluster.json"
V47 = ROOT / "docs" / "verification" / "v47_b5_cluster.json"
LEXFIX = RUNS_DIR / "world_v2_lexfix"
OR_MAX = float(np.exp(3.0))

# 実行 → [(key, 基準, 名前あり, placebo)]。check は同じ表・同じ種の直した実行と予測が一致するかの照合
RUNS = {
    "strict": ("world_v2_b6strict", [("strict", "M1cS", "M3nS", "M3nSp"), ("primary", "M1c", "M3n", "M3np")]),
    "t0x": ("world_v2_b6t0xb", [  # 最初の実行（world_v2_b6t0x）は T0x の腕が回らなかった（f516c6a で修正）
           ("T0x", "T0xc", "T0xn", None), ("T0", "T0c", "T0n", None),
                               ("primary", "M1c", "M3n", None)]),
    "gn": ("world_v2_b6gn", [("primary", "M1c", "M3n", "M3np"), ("total", "M1s", "M3s", "M3ps")]),
    "nohydro": ("world_v2_b6nohydro", [("primary", "M1c", "M3n", "M3np"), ("total", "M1s", "M3s", "M3ps")]),
    "typed": ("world_v2_b6typed", [("primary", "M1c", "M3n", "M3np"), ("total", "M1s", "M3s", "M3ps")]),
    "nread": ("world_v2_b6nread", [("H4", "M1rc", "M3rn", "M3rnp"), ("primary", "M1c", "M3n", "M3np")]),
    "label_oracle": ("world_v2_b6lab", [("oracle_label", "M1cL", "M3oL", None), ("primary", "M1c", "M3n", None)]),
    "gsw_oracle": ("world_v2_b6gsw", [("oracle_gsw", "M1cG", "M3oG", None),
                                      ("oracle_gsw_nonperm", "M1cGn", "M3oGn", None),
                                      ("primary", "M1c", "M3n", None)]),
}
SAME_AS_LEXFIX = {"strict", "t0x", "nread", "label_oracle", "gsw_oracle"}  # 同じ表・種 → M1c・M3n が一致するはず


def _complete(run: Path, pattern: str = "fold[0-9].npz") -> bool:
    d = run / "predictions" / "primary"
    return d.is_dir() and len(list(d.glob(pattern))) >= 10


def compare(key: str) -> dict | str:
    tag, pairs = RUNS[key]
    run = RUNS_DIR / tag
    if not _complete(run):
        return "missing"
    data = wa.load_split(run, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    s, n, arms = wa.unit_sums(y, preds, data["unit"])
    idx = {a: j for j, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    delta = wa.make_delta(hat, wa.bootstrap_bs(s, n), idx)
    d01, b01 = delta("M0", "M1")
    out: dict = {"run": tag, "n_cells": len(y), "n_units": len(n), "n_aoi": len(np.unique(data["aoi"])),
                 "wet_rate": float(y.mean()), "brier": {a: float(hat[idx[a]]) for a in arms},
                 "yardstick": d01, "ci95_yardstick": wa.ci(b01, 0.95)}
    rec = json.loads((run / "predictions" / "primary" / "fold0.json").read_text(encoding="utf-8"))
    if rec.get("b6"):  # 実行の記録（地名辞典の絞り込みの件数など）
        out["b6_record"] = rec["b6"]
    for name, base, alt, plac in pairs:
        if base in idx and alt in idx:
            v, p = wa.verdict(delta, base, alt, plac if plac in idx else None, d01=d01, b01=b01)
            out[name] = v | ({"placebo": p} if p else {})
    if key in SAME_AS_LEXFIX and _complete(LEXFIX):
        ref = wa.load_split(LEXFIX, "primary")
        same_y = len(ref["y"]) == len(y) and bool(np.array_equal(ref["y"], data["y"]))
        out["check_vs_lexfix"] = {"same_cells": same_y} | (
            {f"max_abs_diff_{a}": float(np.max(np.abs(ref[f"p_{a}"].astype(float) - preds[a].astype(float))))
             for a in ("M1c", "M3n") if a in preds} if same_y else {})
    return out


# %% 6(f) activation で束ねた分割の置換分布
def actperm() -> dict | str:
    run = RUNS_DIR / "world_v2_b6actperm"
    d = run / "predictions" / "primary"
    files = sorted(d.glob("fold*_perms.npz")) if d.is_dir() else []
    if len(files) < 10:
        return f"missing ({len(files)} of 10 folds)"
    y, m1c, m3n, m3np, unit, perms = [], [], [], [], [], []
    for f in files:
        with np.load(f) as z:
            y.append(z["y"].astype(np.float64))
            m1c.append(z["p_M1c"].astype(np.float64))
            m3n.append(z["p_M3n"].astype(np.float64))
            m3np.append(z["p_M3np"].astype(np.float64))
            unit.append(z["unit"])
            perms.append(z["p_perms"].astype(np.float64))
    y, m1c, m3n, m3np = (np.concatenate(v) for v in (y, m1c, m3n, m3np))
    pp = np.concatenate(perms, axis=1)
    out: dict = {"n_perms": int(pp.shape[0]), "n_cells": len(y)}
    ref = RUNS_DIR / "world_v2_actgrp" / "metrics_primary.json"  # B-2 の activation で束ねた分割の確証の形
    if ref.is_file():
        m = json.loads(ref.read_text(encoding="utf-8"))
        out["actgrp_reference"] = {k: m["primary"][k] for k in ("R", "ci95_R", "ci90_R", "decision",
                                                                 "delta_bs_M0_M1")}
        yard = float(m["primary"]["delta_bs_M0_M1"])
    else:
        yard = None
    bs_c = float(np.mean((m1c - y) ** 2))
    d_real = bs_c - float(np.mean((m3n - y) ** 2))
    d_plac = bs_c - float(np.mean((m3np - y) ** 2))
    d_perm = bs_c - np.mean((pp - y[None, :]) ** 2, axis=1)
    out |= {"delta_bs_real": d_real, "delta_bs_placebo_single": d_plac,
            "delta_bs_perm_mean": float(d_perm.mean()), "delta_bs_perm_q025": float(np.quantile(d_perm, 0.025)),
            "delta_bs_perm_q975": float(np.quantile(d_perm, 0.975)),
            "share_perm_ge_real": float(np.mean(d_perm >= d_real)),
            "percentile_of_real": float(100 * np.mean(d_perm < d_real)),
            "percentile_of_single_placebo": float(100 * np.mean(d_perm < d_plac))}
    if yard:
        out |= {"yardstick": yard, "R_real": d_real / yard, "R_perm_mean": float(d_perm.mean() / yard)}
    return out


# %% 2. plasmode の上限
def _cross(f: np.ndarray, grid: np.ndarray) -> tuple[float | None, str]:
    """f(δ) ≤ 0 を満たす最大の δ（格子の間は線形補間）。f は δ = 0 を含む格子の上の値。"""
    ok = np.flatnonzero(f <= 0)
    if len(ok) == 0:
        return None, "below grid (f(0) > 0)"
    j = int(ok.max())
    if j == len(grid) - 1:
        return None, "above grid"
    x0, x1, y0, y1 = grid[j], grid[j + 1], f[j], f[j + 1]
    return float(x0 + (0 - y0) * (x1 - x0) / (y1 - y0)), "interpolated"


def _first_reach(g: np.ndarray, target: float, grid: np.ndarray) -> tuple[float | None, str]:
    """g(δ) が target に初めて届く δ（線形補間）。"""
    hit = np.flatnonzero(g >= target)
    if len(hit) == 0:
        return None, "above grid"
    j = int(hit.min())
    if j == 0:
        return 0.0, "at zero"
    x0, x1, y0, y1 = grid[j - 1], grid[j], g[j - 1], g[j]
    return float(x0 + (target - y0) * (x1 - x0) / (y1 - y0)), "interpolated"


def plasmode() -> dict | str:
    d = RUNS_DIR / "world_v2_b6plasmode" / "predictions" / "primary"
    files = sorted(d.glob("fold[0-9]_plasmode.npz")) if d.is_dir() else []
    if len(files) < 10:
        return f"missing ({len(files)} of 10 folds)"
    lex = json.loads((LEXFIX / "metrics_primary.json").read_text(encoding="utf-8"))["primary"]
    v47 = json.loads(V47.read_text(encoding="utf-8"))["permutations"]
    yard = float(lex["delta_bs_M0_M1"])
    u90, u95 = float(lex["ci90_delta"][1]), float(lex["ci95_delta"][1])
    d_null = float(v47["R_perm_mean"]) * float(v47["yardstick"])
    margin = wa.SESOI * yard
    se = (float(lex["ci95_delta"][1]) - float(lex["ci95_delta"][0])) / 3.92  # 追記 1 の 8: 正規近似の標準誤差
    z80 = 1.959964 + 0.841621  # 検出力 0.8: (ΔBS_null + gain) / SE − 1.96 = Φ⁻¹(0.8)
    zs = [np.load(f) for f in files]
    forms = [str(x) for x in zs[0]["forms"]]
    deltas = np.asarray(zs[0]["deltas"], dtype=np.float64)
    n = float(sum(int(z["n_te"]) for z in zs))
    sse_c = sum(z["sse_c"] for z in zs)
    sse_n = sum(z["sse_n"] for z in zs)
    flips = sum(z["flips_te"] for z in zs)
    d0 = float(sum(float(z["base_c"]) - float(z["base_n"]) for z in zs) / n)
    dbs = (sse_c - sse_n) / n  # [form, δ, 回]
    grid = np.concatenate([[0.0], deltas])
    out: dict = {"n_cells": int(n), "reps": int(dbs.shape[2]), "deltas": deltas.tolist(),
                 "yardstick_lexfix": yard, "U90": u90, "U95": u95, "delta_bs_null_perm_mean": d_null,
                 "delta_bs_0": d0, "delta_bs_0_lexfix": float(lex["delta_bs_M1c_M3n"]), "margin": margin,
                 "se_normal": se, "forms": {}}
    for i, form in enumerate(forms):
        gain = dbs[i].mean(axis=1) - d0
        gain_se = dbs[i].std(axis=1, ddof=1) / np.sqrt(dbs.shape[2])  # 回の間（正規近似の se とは別）
        g = np.concatenate([[0.0], gain])
        res = {"gain": gain.tolist(), "gain_se": gain_se.tolist(), "R_gain": (gain / yard).tolist(),
               "monotone": bool(np.all(np.diff(g) >= 0)),
               "label_flip_share": (flips[i].mean(axis=1) / n).tolist()}
        for name, f in (("delta_A", d_null + g - u90), ("delta_B", g - u90),
                        ("delta_A_U95", d_null + g - u95), ("delta_B_U95", g - u95)):
            x, how = _cross(f, grid)
            res[name] = {"delta": x, "odds_ratio": float(np.exp(x)) if x is not None else None, "how": how,
                         "odds_ratio_text": (f"{np.exp(x):.2f}" if x is not None else
                                             f"> {OR_MAX:.0f}" if how == "above grid" else "undefined")}
        for name, tgt in (("delta_M", margin), ("delta_M_prime", margin - d_null), ("mde80", z80 * se - d_null)):
            x, how = _first_reach(g, tgt, grid)
            res[name] = {"delta": x, "odds_ratio": float(np.exp(x)) if x is not None else None, "how": how,
                         "odds_ratio_text": (f"{np.exp(x):.2f}" if x is not None else f"> {OR_MAX:.0f}")}
        out["forms"][form] = res
    return out


def main() -> None:
    res: dict = {"rule": "amendment B-6 items 2, 3, 5, 6"}
    for key in RUNS:
        res[key] = compare(key)
    res["actperm"] = actperm()
    res["plasmode"] = plasmode()
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for key in RUNS:
        v = res[key]
        if not isinstance(v, dict):
            print(f"{key:13s} {v}")
            continue
        for name, *_ in RUNS[key][1]:
            c = v.get(name)
            if c:
                print(f"{key:13s} {name:18s} R {c['R']:+.4f} 90% [{c['ci90_R'][0]:+.4f}, {c['ci90_R'][1]:+.4f}] "
                      f"{c['decision']:22s} cells {v['n_cells']:,} AOIs {v['n_aoi']}")
        if "check_vs_lexfix" in v:
            print(f"{'':13s} check vs lexfix {v['check_vs_lexfix']}")
    a = res["actperm"]
    print("actperm", a if not isinstance(a, dict) else {k: a.get(k) for k in ("percentile_of_real", "R_real",
                                                                               "share_perm_ge_real")})
    p = res["plasmode"]
    if isinstance(p, dict):
        print("plasmode dBS(0)", p["delta_bs_0"], "lexfix", p["delta_bs_0_lexfix"], "null", p["delta_bs_null_perm_mean"],
              "U90", p["U90"])
        for form, r in p["forms"].items():
            print(form, "A", r["delta_A"]["odds_ratio_text"], "B", r["delta_B"]["odds_ratio_text"],
                  "M", r["delta_M"]["odds_ratio_text"], "M'", r["delta_M_prime"]["odds_ratio_text"],
                  "gain/yard", [round(x, 4) for x in r["R_gain"]], "monotone", r["monotone"])
    else:
        print("plasmode", p)


if __name__ == "__main__":
    main()
