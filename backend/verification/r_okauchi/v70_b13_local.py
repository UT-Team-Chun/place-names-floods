# %% v70（事後・amendment B-13 の B1〜B15）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-13（追記 1〜）で計算の前に固定した量を、既存の表と held-out の予測から計算する。

- B1: 比の余白に対する Fieller 型の同等性の検定（L± = ΔBS ∓ 0.05·Δ_base）。(a) 場所グループ、モデル固定、
  (b) 場所グループ、二段目を当て直す、(c) 国の中の activation を二段で再抽出、二段目を当て直す、(d) (c) に fold の
  割り当ての分散（Nadeau–Bengio）。主比較・H4 の入れ子と全体・国の分割。動作特性（L 規則）
- B2: 各区間の上限に対応する注入の OR（名目と推定量の尺度、無効果の 4 つの基準、モンテカルロの範囲）
- B3: 符号化者 2 の p0・p1 と、その率での誤読の補正
- B4: 集落の近くの規模の層（v55 の三分位）
- B5: 文字モデルの得点の AOI 平均を除いた M1c と M3n の R
- B7: 分母の代替（Y3）・ΔBS・気候値に対する Brier skill の増分
- B13: AOI で釣り合わせた指標（ΔAUC・ΔPR-AUC の中央値、AOI で重みを揃えた Brier skill の増分）

出力: stdout と docs/verification/v70_b13_local.json（部分ごとに足す。`--only B1,B4` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import correct_misclassification, stage2_weighted
from topsem.world.b13 import (
    SESOI_R,
    fieller_interval,
    l_pm,
    nadeau_bengio_var,
    nested_weights,
    tost_from_replicates,
)

OUT = ROOT / "docs" / "verification" / "v70_b13_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
SEED = 20261102
B_FIXED = 4000
B_REFIT = 500
B_OUT, B_IN = 2000, 1000
B_PBA = 20000
N_WORKERS = int(os.environ.get("V70_WORKERS", "3"))
FOLD_RUNS = ("world_v2_lexfix", *(f"world_v2_fs{k}" for k in ("1b", "2b", 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)))
Z95 = 1.6448536269514722  # 片側 95%


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v69 = _load("v69_b12_local")
v36, v68 = v69.v36, v69.v68


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


# %% B1(a): モデル固定の場所グループの再抽出
COMPARISONS = (
    ("confirmatory", "world_v2", "primary", ("M0", "M1", "M1c", "M3n")),
    ("corrected", "world_v2_lexfix", "primary", ("M0", "M1", "M1c", "M3n")),
    ("H4_nested", "world_v2", "llmcut", ("M0", "M1", "M1rc", "M3rn")),
    ("H4_total", "world_v2", "llmcut", ("M0", "M1", "M1s", "M3rs")),
    ("country", "world_v2", "country", ("M0", "M1", "M1c", "M3n")),
)


def _point(e: float, d: float) -> dict:
    lp, lm = l_pm(e, d)
    return {"R": e / d, "delta_bs": e, "yardstick": d, "L_plus": float(lp), "L_minus": float(lm)}


def _summ(e: float, d: float, eb: np.ndarray, db: np.ndarray) -> dict:
    """点の推定と再標本 (eb, db) から TOST・Fieller の区間・比の百分位の区間。"""
    lp, lm = l_pm(eb, db)
    ok = np.isfinite(eb) & np.isfinite(db)
    c = np.cov(eb[ok], db[ok])
    fi = fieller_interval(e, d, (float(c[0, 0]), float(c[1, 1]), float(c[0, 1])))
    r = eb[ok] / db[ok]
    return (_point(e, d) | tost_from_replicates(lp, lm)
            | {"fieller90": list(fi) if fi else None, "fieller_within_margin": bool(fi and -SESOI_R < fi[0]
                                                                                    and fi[1] < SESOI_R),
               "ratio_pct90": np.quantile(r, [0.05, 0.95]).tolist(), "sd_delta_bs": float(np.std(eb[ok], ddof=1)),
               "sd_yardstick": float(np.std(db[ok], ddof=1)), "corr": float(np.corrcoef(eb[ok], db[ok])[0, 1]),
               "var": [float(c[0, 0]), float(c[1, 1]), float(c[0, 1])]})


def b1a_fixed() -> dict:
    rng = np.random.default_rng(SEED)
    out: dict = {"replicates": B_FIXED, "seed": SEED}
    for label, run, split, (a0, a1, ac, an) in COMPARISONS:
        d = v36.load(run, split)
        if d is None or any(f"p_{a}" not in d for a in (a0, a1, ac, an)):
            out[label] = {"note": f"{run}/{split} not found locally"}
            continue
        n, s, units = v69._group_sums(d, (a0, a1, ac, an))
        eg, dg = s[ac] - s[an], s[a0] - s[a1]
        big_n = float(n.sum())
        g = len(units)
        w = rng.multinomial(g, np.full(g, 1.0 / g), size=B_FIXED).astype(np.float64)
        nw = w @ n
        out[label] = {"run": run, "split": split, "arms": [a0, a1, ac, an], "n_units": g,
                      **_summ(float(eg.sum() / big_n), float(dg.sum() / big_n), (w @ eg) / nw, (w @ dg) / nw)}
    return out


# %% B1(b)(c): 二段目を当て直す再抽出（B-11 の B7 と同じ設計の保存と符号）
_B: dict = {}


def _b1_init(country_of_act: np.ndarray) -> None:
    v68._b7_init()
    _B.update(v68._B7)
    _B["country_of_act"] = country_of_act


def _b1_rep(job: tuple[str, int]) -> list[float]:
    """1 回の再標本: [当て直した M1c, M3n, 固定の M1c, M3n, M0, M1 の重みつき SSE, 評価側の重みの和]。"""
    scheme, b = job
    if scheme == "groups":
        rng = np.random.default_rng(SEED + b)
        wu = rng.multinomial(_B["nu"], np.full(_B["nu"], 1 / _B["nu"])).astype(np.float64)
    else:
        rng = np.random.default_rng(SEED + 100_000 + b)
        wa = nested_weights(rng, _B["country_of_act"])
    acc = np.zeros(7)
    for f in _B["folds"]:
        if scheme == "groups":
            w_tr, w_te = wu[f["u_tr"]], wu[f["u_te"]]
        else:
            w_tr, w_te = wa[f["a_tr"]], wa[f["a_te"]]
        if w_te.sum() == 0 or w_tr.sum() == 0:
            continue
        pc, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xc_tr"], f["xc_te"], f["y_tr"], w_tr=w_tr, init=f["init_c"])
        pn, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xn_tr"], f["xn_te"], f["y_tr"], w_tr=w_tr, init=f["init_n"])
        y = f["y_te"].astype(np.float64)
        for i, p in enumerate((pc, pn, f["p_M1c_te"], f["p_M3n_te"], f["p_M0_te"], f["p_M1_te"])):
            acc[i] += float((w_te * (p - y) ** 2).sum())
        acc[6] += float(w_te.sum())
    return acc.tolist()


def _country_of_activation() -> tuple[np.ndarray, dict]:
    """activation → その activation のセルが最も多い国（二段の再抽出の入れ子に使う）。"""
    c = np.load(v68.B7_CODES)
    na, nc = int(c["na"]), int(c["nc"])
    cnt = np.zeros((na, nc), dtype=np.int64)
    for k in range(10):
        np.add.at(cnt, (c[f"a_te{k}"], c[f"c_te{k}"]), 1)
    coa = cnt.argmax(axis=1)
    multi = int(((cnt > 0).sum(axis=1) > 1).sum())
    return coa, {"n_activations": na, "n_countries": nc, "n_activations_in_several_countries": multi}


def b1bc_refit() -> dict:
    coa, info = _country_of_activation()
    info["n_countries_with_activations"] = len(np.unique(coa))  # 追記 6: 引くのはこの国だけ
    jobs = [(s, b) for s in ("groups", "nested") for b in range(B_REFIT)]
    _threads1()
    parts = []
    for scheme in ("groups", "nested"):  # 追記 6: 方式ごとに保存し、保存があれば読む
        f = VER / f"v70_b1_refit_reps_{scheme}.npy"
        if f.is_file():
            parts.append(np.load(f))
            continue
        t0 = time.time()
        with ProcessPoolExecutor(N_WORKERS, initializer=_b1_init, initargs=(coa,)) as ex:
            r = np.array(list(ex.map(_b1_rep, [j for j in jobs if j[0] == scheme], chunksize=2)), dtype=np.float64)
        np.save(f, r)
        print(f"  B1({scheme}): {len(r)} replicates ({time.time() - t0:.0f} s)", flush=True)
        parts.append(r)
    res = np.concatenate(parts)
    # 点の推定: 重みの無い当て直し（= 直した実行の二段目）
    lex = v36.load("world_v2_lexfix", "primary")
    y = lex["y"].astype(np.float64)
    e0 = float(np.mean((lex["p_M1c"] - y) ** 2) - np.mean((lex["p_M3n"] - y) ** 2))
    d0 = float(np.mean((lex["p_M0"] - y) ** 2) - np.mean((lex["p_M1"] - y) ** 2))
    out: dict = {"replicates": B_REFIT, "seed": SEED, "nesting": info,
                 "rule": "B-13 B1 (b)(c): stage 2 refitted in each replicate on the stage-2 design of the corrected "
                         "code (world_v2_b8stage2_local, as B-11 B7); (b) location groups, (c) countries then "
                         "activations within each drawn country (activation assigned to its majority country)"}
    for k, scheme in (("b_groups", "groups"), ("c_nested", "nested")):
        r = res[[s == scheme for s, _ in jobs]]
        w = r[:, 6]
        eb_ref, eb_fix = (r[:, 0] - r[:, 1]) / w, (r[:, 2] - r[:, 3]) / w
        db = (r[:, 4] - r[:, 5]) / w
        out[k] = {"refit": _summ(e0, d0, eb_ref, db), "fixed": _summ(e0, d0, eb_fix, db)}
    return out


# %% B1(d): fold の割り当ての分散を足す
def b1d_folds(b1bc: dict) -> dict:
    vals = []
    for run in FOLD_RUNS:
        d = v36.load(run, "primary")
        if d is None:
            return {"note": f"{run} predictions not found locally"}
        y = d["y"].astype(np.float64)
        e = float(np.mean((d["p_M1c"] - y) ** 2) - np.mean((d["p_M3n"] - y) ** 2))
        dd = float(np.mean((d["p_M0"] - y) ** 2) - np.mean((d["p_M1"] - y) ** 2))
        vals.append((e, dd))
    v = np.array(vals)
    lp, lm = l_pm(v[:, 0], v[:, 1])
    nb = {"L_plus": nadeau_bengio_var(lp), "L_minus": nadeau_bengio_var(lm), "delta_bs": nadeau_bengio_var(v[:, 0]),
          "yardstick": nadeau_bengio_var(v[:, 1])}
    cov_nb = float((1 / len(v) + 1 / 9) * np.cov(v[:, 0], v[:, 1])[0, 1])
    c = b1bc["c_nested"]["refit"]
    e0, d0 = c["delta_bs"], c["yardstick"]
    ve, vd, vc = c["var"]
    lp0, lm0 = l_pm(e0, d0)
    var_lp = ve + SESOI_R**2 * vd - 2 * SESOI_R * vc + nb["L_plus"]
    var_lm = ve + SESOI_R**2 * vd + 2 * SESOI_R * vc + nb["L_minus"]
    up, lo = float(lp0 + Z95 * math.sqrt(var_lp)), float(lm0 - Z95 * math.sqrt(var_lm))
    fi = fieller_interval(e0, d0, (ve + nb["delta_bs"], vd + nb["yardstick"], vc + cov_nb))
    return {"runs": list(FOLD_RUNS), "R_by_run": (v[:, 0] / v[:, 1]).tolist(), "nb_var": nb, "nb_cov": cov_nb,
            "L_plus": float(lp0), "L_minus": float(lm0), "l_plus_upper": up, "l_minus_lower": lo,
            "equivalent": bool(up < 0 < lo), "fieller90": list(fi) if fi else None,
            "fieller_within_margin": bool(fi and -SESOI_R < fi[0] and fi[1] < SESOI_R),
            "rule": "B-13 B1 (d): variance of (c) (nested, refit) plus the Nadeau-Bengio variance over the 13 fold "
                    "assignments of v47; normal approximation, one-sided 95%"}


# %% B1: 動作特性（L 規則）
R_GRID = (-0.06, -0.05, -0.04, -0.02, 0.0, 0.02, 0.04, 0.05, 0.06)


def b1_operating() -> dict:
    """B-12 の B1 と同じ二重の再抽出（確証の実行、モデル固定）で、固定の規則・比の規則・L 規則が同等と言う割合。"""
    d = v36.load("world_v2", "primary")
    n, s, units = v69._group_sums(d, ("M0", "M1", "M1c", "M3n"))
    g = len(units)
    dg, eg = s["M0"] - s["M1"], s["M1c"] - s["M3n"]
    d_hat, e_hat = float(dg.sum() / n.sum()), float(eg.sum() / n.sum())
    rng = np.random.default_rng(SEED)
    rules = ("fixed", "ratio", "L")
    cnt = {sh: {r0: dict.fromkeys((*rules, "undefined"), 0) for r0 in R_GRID} for sh in ("proportional", "constant")}
    t0 = time.time()
    for _ in range(B_OUT):
        w = rng.multinomial(g, np.full(g, 1.0 / g)).astype(np.float64)
        nw = w @ n
        d_b = (w @ dg) / nw
        e_b = (w @ eg) / nw
        v = rng.multinomial(g, w / g, size=B_IN).astype(np.float64)
        nv = v @ n
        dv, ev = (v @ dg) / nv, (v @ eg) / nv
        d_lo = float(np.quantile(dv, 0.025))
        for r0 in R_GRID:
            for sh in ("proportional", "constant"):
                e_in = ev - e_hat + r0 * (dv if sh == "proportional" else d_hat)
                c = cnt[sh][r0]
                if d_lo <= 0:
                    c["undefined"] += 1
                    continue
                lo, hi = np.quantile(e_in, [0.05, 0.95]) / d_b
                c["fixed"] += int(-SESOI_R < lo and hi < SESOI_R)
                rlo, rhi = np.quantile(e_in / dv, [0.05, 0.95])
                c["ratio"] += int(-SESOI_R < rlo and rhi < SESOI_R)
                lp, lm = l_pm(e_in, dv)
                c["L"] += int(np.quantile(lp, 0.95) < 0 < np.quantile(lm, 0.05))
        del e_b
    res = {sh: {f"{r0:+.2f}": {k: v / B_OUT for k, v in c.items()} for r0, c in cs.items()} for sh, cs in cnt.items()}
    return {"n_units": g, "replicates_outer": B_OUT, "replicates_inner": B_IN, "seed": SEED,
            "share_equivalent": res, "seconds": round(time.time() - t0),
            "rule": "B-13 B1: as B-12 B1 (true R set to r0 by proportional or constant shifts), adding the L rule "
                    "(one-sided 95% upper limit of L+ below 0 and lower limit of L- above 0)"}


# %% B2: 上限の換算
def _conv_factory() -> tuple[dict, callable]:
    v66 = _json(DOCV / "v66_b10_local.json")
    rows = v66["B11"]["rows"]
    curves = {}
    for form, col in (("A", "mean_or_ratio_to_base_all"), ("B", "mean_or_ratio_to_base_le500m")):
        rr = sorted((r["nominal_or"], r[col]) for r in rows if r["form"] == form)
        curves[form] = (np.log(np.array([1.0] + [a for a, _ in rr])), np.log(np.array([1.0] + [b for _, b in rr])))

    def conv(form: str, nominal: float | None) -> dict | None:
        if nominal is None or not np.isfinite(nominal):
            return None
        x, yv = curves[form]
        if np.log(nominal) > x[-1]:
            return {"nominal": nominal, "estimator": None, "text": f"> {np.exp(yv[-1]):.2f}"}
        e = float(np.exp(np.interp(np.log(nominal), x, yv)))
        return {"nominal": nominal, "estimator": e, "text": f"{e:.2f}"}

    return curves, conv


def b2_bounds(b1: dict) -> dict:
    v64 = _load("v64_b9_aggregate")
    v59 = v64.v59
    grids, gains, _ = v64.combined_curves()
    refs = v59.references()
    lex = _json(RUNS_DIR / "world_v2_lexfix" / "metrics_primary.json")["primary"]
    yard = float(lex["delta_bs_M0_M1"])
    c2 = _json(DOCV / "v59_b8_aggregate.json")["C2"]["corrected"]
    b7 = _json(DOCV / "v68_b11_local.json")["B7"]
    fs = _json(DOCV / "v47_b5_cluster.json")["fold_seeds"]["fixed_code"]
    uppers = {"fixed": float(lex["ci90_delta"][1]) / yard, "refit_b8_c2": float(c2["refit"]["ci90_R"][1]),
              "groups_refit_b11_b7": float(b7["groups"]["ci90_R_refit_ratio"][1]),
              "two_way_b11_b7": float(b7["two_way"]["ci90_R_refit_ratio"][1]),
              "fold_assignments_v47": float(fs["combined_ci90"][1])}
    for k, v in (("b1_a_corrected", ((b1.get("B1a") or {}).get("corrected") or {}).get("fieller90")),
                 ("b1_b_groups_refit", (((b1.get("B1bc") or {}).get("b_groups") or {}).get("refit") or {})
                  .get("fieller90")),
                 ("b1_c_nested_refit", (((b1.get("B1bc") or {}).get("c_nested") or {}).get("refit") or {})
                  .get("fieller90")),
                 ("b1_d_with_folds", (b1.get("B1d") or {}).get("fieller90"))):
        if v:
            uppers[k] = float(v[1])
    _, conv = _conv_factory()
    rng = np.random.default_rng(SEED)
    out: dict = {"yardstick": yard, "upper_R": uppers, "references": refs, "bounds": {}}
    for kind, ur in uppers.items():
        out["bounds"][kind] = {}
        for ref, dnull in refs.items():
            out["bounds"][kind][ref] = {}
            for f in ("A", "B"):
                grid = np.concatenate([[0.0], grids[f]])
                e = v59._entry(gains[f], grid, "bound", ur * yard - dnull, float(np.exp(grids[f][-1])), rng)
                e["estimator_scale"] = conv(f, e.get("odds_ratio"))
                out["bounds"][kind][ref][f] = e
    out["rule"] = ("B-13 B2: injected odds ratio whose mean gain reaches each interval's upper limit (in Brier units, "
                   "upper R times the corrected run's yardstick) minus each no-effect reference, from the curves of "
                   "v64; Monte Carlo ranges as in v59; estimator scale by the B-10 B11 curves (v68 B3 rule)")
    return out


# %% B3: 符号化者 2 の率
def _wilson(k: int, n: int, z: float = 1.959964) -> list[float]:
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [c - h, c + h]


def _coder2_frame() -> object:
    import pandas as pd

    v44 = _load("v44_wikipedia_etymology")
    samp = json.loads(v44.SAMPLE.read_text(encoding="utf-8"))
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    to_code = DATA_DIR / "verification" / "v44_to_code.jsonl"
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in to_code.read_text(encoding="utf-8").splitlines()}
    c2 = pd.read_csv(DOCV / "v57b_codes_coder2.csv")
    c2["cats2"] = c2["categories"].fillna("").map(lambda s: {c for c in str(s).split(";") if c})
    c2["gid"] = c2["id"].map(key)
    return c2.join(meta, on="gid"), samp["stratum_sizes"]


def b3_coder2() -> dict:
    from topsem.world.lexicon import WET_SIDE

    d, sizes = _coder2_frame()
    d = d[d["status"] == "coded"].copy()
    d["truth"] = d["cats2"].map(lambda s: bool(s & set(WET_SIDE)))
    wet, nott = d[d["lex_wetside"].astype(bool)], d[~d["lex_wetside"].astype(bool)]
    p1 = {"k": int(wet["truth"].sum()), "n": len(wet)}
    p0 = {"k": int(nott["truth"].sum()), "n": len(nott)}
    for p in (p1, p0):
        p["rate"] = p["k"] / p["n"]
        p["ci95"] = _wilson(p["k"], p["n"])
    v44 = _json(DOCV / "v44_wikipedia_etymology.json")
    b = np.load(VER / "v57_pairs_draws_main_wetside.npy")
    b = b[np.isfinite(b)]
    rng = np.random.default_rng(SEED)
    orr = np.exp(rng.choice(b, size=B_PBA, replace=True))
    out: dict = {"coder2": {"p1": p1, "p0": p0}, "coder1": {"p1": v44["truth_rate_in_wet"],
                                                              "p0": v44["truth_rate_in_not"]},
                 "stratum_sizes": sizes}
    for who in ("coder1", "coder2"):
        a, z = out[who]["p1"], out[who]["p0"]
        q1 = rng.beta(a["k"] + 1, a["n"] - a["k"] + 1, size=B_PBA)
        q0 = rng.beta(z["k"] + 1, z["n"] - z["k"] + 1, size=B_PBA)
        c = correct_misclassification(orr, q1, q0)
        ok = np.isfinite(c)
        out[who]["misclassification"] = {
            "point_at_rates": float(correct_misclassification(np.exp(np.median(b)), a["rate"], z["rate"])),
            "median": float(np.median(c[ok])) if ok.any() else None,
            "ci95": [float(np.quantile(c[ok], 0.025)), float(np.quantile(c[ok], 0.975))] if ok.any() else None,
            "share_undefined": float(1 - ok.mean()), "draws": B_PBA}
        # 層の重みで直した再現率（p0 を層の大きさで重みづけ）
        nw, nn = sizes["wet"] if "wet" in sizes else sizes.get("wetside"), sizes.get("not")
        if nw and nn:
            out[who]["recall_implied"] = a["rate"] * nw / (a["rate"] * nw + z["rate"] * nn)
    out["rule"] = ("B-13 B3: p1, p0 from coder 2's codes (coded rows; water-related = any of the five water-side "
                   "classes), Wilson 95% intervals; probabilistic bias analysis of the main within-500 m odds ratio "
                   "(cached pairs draws of v57) with Beta(k+1, n-k+1) rates for each coder")
    return out


# %% B4: 規模の層
def b4_size() -> dict:
    sev = _json(DOCV / "v55_b7_local.json")["severity"]
    rows = {k: {kk: sev[k].get(kk) for kk in ("odds_ratio", "ci95", "n_cells", "n_aoi", "n_units", "n_group_cells",
                                               "n_group_flooded", "n_group_units", "n_aoi_tercile")}
            for k in ("low", "middle", "high") if isinstance(sev.get(k), dict)}
    return {"tercile_cuts": sev.get("tercile_cuts"), "strata": rows,
            "rule": "B-13 B4: the size terciles of v55 (share of the AOI below the low-HAND threshold), table only"}


# %% B5: 文字モデルの地域構成を除いた基準
def b5_reference() -> dict:
    folds = v69._dump_folds()
    sse = dict.fromkeys(("M1c", "M3n", "M1c_no_nlearn", "M3n_no_nlearn_comp", "M0", "M1"), 0.0)
    n = 0
    for z in folds:
        comp = [str(c) for c in z["comp_cols"]]
        keep = [i for i, c in enumerate(comp) if c != "c_nlearn"]
        y = z["y_te"].astype(np.float64)
        xn_tr = np.column_stack([z["xc_tr"], z["xn_tr"]])
        xn_te = np.column_stack([z["xc_te"], z["xn_te"]])
        xk_tr = np.column_stack([z["xc_tr"][:, keep], z["xn_tr"]])
        xk_te = np.column_stack([z["xc_te"][:, keep], z["xn_te"]])
        preds = {"M1c": stage2_weighted(z["z_tr"], z["z_te"], z["xc_tr"], z["xc_te"], z["y_tr"])[0],
                 "M3n": stage2_weighted(z["z_tr"], z["z_te"], xn_tr, xn_te, z["y_tr"])[0],
                 "M1c_no_nlearn": stage2_weighted(z["z_tr"], z["z_te"], z["xc_tr"][:, keep], z["xc_te"][:, keep],
                                                  z["y_tr"])[0],
                 "M3n_no_nlearn_comp": stage2_weighted(z["z_tr"], z["z_te"], xk_tr, xk_te, z["y_tr"])[0],
                 "M0": z["p_M0_te"], "M1": z["p_M1_te"]}
        for k, p in preds.items():
            sse[k] += float(((p.astype(np.float64) - y) ** 2).sum())
        n += len(y)
    bs = {k: v / n for k, v in sse.items()}
    yard = bs["M0"] - bs["M1"]
    return {"brier": bs, "yardstick": yard, "R_standard": (bs["M1c"] - bs["M3n"]) / yard,
            "R_without_char_composition": (bs["M1c_no_nlearn"] - bs["M3n_no_nlearn_comp"]) / yard,
            "rule": "B-13 B5: second stage refitted on the dumps of B-8 C2 (corrected code, primary split); the "
                    "reference and the names model both without the AOI mean of the character model's score"}


# %% B7: 分母の代替と絶対の量
def b7_denominators(b1a: dict) -> dict:
    b17 = _json(DOCV / "v69_b12_local.json")["B17"]
    out: dict = {}
    for label in ("confirmatory", "corrected", "H4_nested", "H4_total"):
        r = b1a.get(label) or {}
        if "delta_bs" not in r:
            continue
        sp = "llmcut" if label.startswith("H4") else "primary"
        k = b17[sp]
        clim = k.get("brier_training_constant", k["brier_best_constant"])
        y3 = clim - k["brier_M1"]
        out[label] = {"delta_bs": r["delta_bs"], "R_Y0": r["R"], "R_Y3": r["delta_bs"] / y3, "Y3": y3,
                      "climatology_brier": clim,
                      "climatology": "training constant" if "brier_training_constant" in k else "test constant",
                      "bss_increment": r["delta_bs"] / clim}
    out["rule"] = ("B-13 B7: Y3 = Brier of the constant (training-side flooded share per fold; the test share for "
                   "the single-fold llmcut split) minus Brier of M1; Brier skill increment = delta BS / Brier of "
                   "that constant")
    return out


# %% B13: AOI で釣り合わせた指標
def _auc(y: np.ndarray, p: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(y, p))


def _prauc(y: np.ndarray, p: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score

    return float(average_precision_score(y, p))


def b13_aoi_balanced() -> dict:
    out: dict = {}
    for label, run in (("confirmatory", "world_v2"), ("corrected", "world_v2_lexfix")):
        d = v36.load(run, "primary")
        y = d["y"].astype(np.int64)
        aoi = np.asarray(d["aoi"]).astype(str)
        rows = []
        for a in np.unique(aoi):
            m = aoi == a
            ya = y[m]
            if ya.min() == ya.max():
                continue
            pc, pn = d["p_M1c"][m], d["p_M3n"][m]
            clim = float(ya.mean())
            bsc, bsn = float(np.mean((pc - ya) ** 2)), float(np.mean((pn - ya) ** 2))
            rows.append((_auc(ya, pn) - _auc(ya, pc), _prauc(ya, pn) - _prauc(ya, pc),
                         (bsc - bsn) / (clim * (1 - clim))))
        r = np.array(rows)
        out[label] = {"n_aoi": len(r), "median_delta_auc": float(np.median(r[:, 0])),
                      "median_delta_prauc": float(np.median(r[:, 1])),
                      "share_aoi_auc_up": float((r[:, 0] > 0).mean()),
                      "mean_bss_increment_aoi_weighted": float(r[:, 2].mean()),
                      "median_bss_increment": float(np.median(r[:, 2]))}
    out["rule"] = ("B-13 B13: per AOI with both classes, AUC and PR-AUC of M3n minus M1c, and the Brier skill "
                   "increment relative to the AOI's own flooded share (equal weight per AOI)")
    return out


# %% 実行
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    want = {x for x in args.only.split(",") if x}
    res = _json(OUT) if OUT.is_file() else {}
    res |= {"rule": "amendment B-13 (post hoc), local computations; confirmatory verdict unchanged", "seed": SEED}

    def run(key: str, fn: callable, *a: object) -> None:
        if want and key not in want:
            return
        t0 = time.time()
        res[key] = fn(*a)
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        print(f"{key} ({time.time() - t0:.0f} s) {json.dumps(res[key], default=float)[:600]}", flush=True)

    run("B1a", b1a_fixed)
    run("B4", b4_size)
    run("B5", b5_reference)
    run("B3", b3_coder2)
    run("B13", b13_aoi_balanced)
    if "B1a" in res:
        run("B7", b7_denominators, res["B1a"])
    run("B1op", b1_operating)
    run("B1bc", b1bc_refit)
    if "B1bc" in res:
        run("B1d", b1d_folds, res["B1bc"])
    run("B2", b2_bounds, res)


if __name__ == "__main__":
    main()
