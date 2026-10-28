# %% v59（事後・amendment B-8 の B1 と C1〜C5）: 手元で回した B-8 の実行の集計
"""確証判定は変えない。B-8 で計算の前に固定した量を、各実行の出力から計算する（追記 3 で実行は手元の Mac に移した）。

- B1: plasmode の上限・検出力・余白の換算を、無効果の基準 ΔBS_null を替えて読み直す。主の基準は B-7 の 1 で地域構成を
  作り直した置換 200 回の平均（v56）、並べて 0・AOI 内で中心化した置換の平均（v56）・B-5 の置換の平均（v47）。
  曲線は B-6 の基礎を固定した plasmode（world_v2_b6plasmode、形 A・B・CA・CB、回 0〜9）に、C1 の回 10〜49
  （追記 5: クラスタの world_v2_b8plas50、4 形）を足したもの。B-7 の当て直し（形 A、3 回、v56）も同じ基準で読む。
  モンテカルロの区間: 各 δ の中で回を復元抽出（2,000 回、種 20261028）し、平均の曲線を作り直して同じ読み方で読んだ
  オッズ比の 2.5〜97.5 百分位と、格子の中で届かない再標本の割合
- C1: 回 10〜49 の集計と、fold 0 の回 0〜9 を回し直した値と B-6 の値の差（クラスタと、手元の確認）
- C2（world_v2_b8bootref_local）: 二段目を bootstrap の中で当て直した R の 90%・95% 区間と、同じ重みでモデルを固定した区間
- C3（追記 5: クラスタの world_v2_b8margin、回 100〜139）: 余白（OR 2.336）を形 A に注入して全てを当て直した各回の判定（主比較と同じクラスタ
  bootstrap 4,000 回・判定表）と、90% 区間が ±0.05 に入った回の割合（Wilson の 95% 区間）
- C4（world_v2_b8gswcell_local）: セルに直接の GSW の対照 M1cW → M3oW と、B-6 の集落を通じた GSW の対照（v52）の読み
- C5（world_v2_b8lexclean_local）: 掃除した辞書での主比較 M1c → M3n（placebo M3np）

出力: stdout と docs/verification/v59_b8_aggregate.json（無い入力は "missing ..."）
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v59_b8_aggregate.json"
VER = ROOT / "docs" / "verification"
LEXFIX = RUNS_DIR / "world_v2_lexfix"
B6PLAS = RUNS_DIR / "world_v2_b6plasmode"
C1RUN = RUNS_DIR / "world_v2_b8plas50"  # 追記 5: クラスタ（4 形）
C1LOCAL = RUNS_DIR / "world_v2_b8plas50_local"  # 追記 3 の手元の確認（fold 0 の回 0〜9、形 B・CB）
C2DUMP = RUNS_DIR / "world_v2_b8stage2_local"
C2BOOT = RUNS_DIR / "world_v2_b8bootref_local"
C3RUN = RUNS_DIR / "world_v2_b8margin"  # 追記 5: クラスタ
C4RUN = RUNS_DIR / "world_v2_b8gswcell_local"
C5RUN = RUNS_DIR / "world_v2_b8lexclean_local"
SEED_MC = 20261028
B_MC = 2000
N_FOLDS = 10
C1_FIRSTS = (10, 20, 30, 40)
C1_FORMS = ("A", "B", "CA", "CB")
C2_REPS = 500
C3_FIRST, C3_REPS = 100, 40
Z80 = 1.959964 + 0.841621  # 検出力 0.8: (ΔBS_null + gain) / SE − 1.96 = Φ⁻¹(0.8)
REFS = ("b7_recomputed", "zero", "b7_centred", "b5")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v52 = _load("v52_b6_cluster")
v56 = _load("v56_b7_cluster")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _or_text(x: float | None, how: str, or_max: float) -> str:
    if x is not None:
        return f"{np.exp(x):.2f}"
    return f"> {or_max:.0f}" if how == "above grid" else "undefined"


# %% B1: 基準を替えた上限・検出力と、モンテカルロの区間
def _fold_gains(files: list[Path], forms: tuple[str, ...]) -> dict[str, np.ndarray]:
    """10 fold のファイル（同じ回の範囲）から、形ごとの回ごとの利得 ΔBS(δ, 回) − ΔBS(0)（[δ, 回]）。"""
    zs = [np.load(f) for f in files]
    have = [str(x) for x in zs[0]["forms"]]
    n = float(sum(int(z["n_te"]) for z in zs))
    sse_c, sse_n = sum(z["sse_c"] for z in zs), sum(z["sse_n"] for z in zs)
    d0 = float(sum(float(z["base_c"]) - float(z["base_n"]) for z in zs) / n)
    dbs = (sse_c - sse_n) / n  # [形, δ, 回]
    return {f: dbs[have.index(f)] - d0 for f in forms if f in have}


def plasmode_gains() -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    """形ごとの利得の行列 [δ, 回]: B-6 の回 0〜9 と C1 の回 10〜49（全ての fold で揃った範囲だけ）。"""
    d6 = B6PLAS / "predictions" / "primary"
    f6 = [d6 / f"fold{k}_plasmode.npz" for k in range(N_FOLDS)]
    deltas = np.asarray(np.load(f6[0])["deltas"], dtype=np.float64)
    g6 = _fold_gains(f6, ("A", "B", "CA", "CB"))
    gains: dict[str, np.ndarray] = {}
    info: dict = {"reps": {}}
    d1 = C1RUN / "predictions" / "primary"
    ranges, missing_ranges = [], []
    for first in C1_FIRSTS:
        fs = [d1 / f"fold{k}_plasmode_{first:03d}.npz" for k in range(N_FOLDS)]
        if all(f.is_file() for f in fs):
            ranges.append(_fold_gains(fs, C1_FORMS))
        else:
            missing_ranges.append(first)
    for f in C1_FORMS:
        gains[f] = np.concatenate([g6[f], *(r[f] for r in ranges)], axis=1)
        info["reps"][f] = int(gains[f].shape[1])
    info["c1_ranges_missing"] = missing_ranges
    return deltas, gains, info


def references() -> dict[str, float]:
    """無効果の基準 ΔBS_null（Brier の尺度）。"""
    p56 = _json(VER / "v56_b7_cluster.json")["permutations"]
    p47 = _json(VER / "v47_b5_cluster.json")["permutations"]
    y56 = float(p56["yardstick_lexfix"])
    return {"b7_recomputed": float(p56["recomputed_composition"]["R_perm_mean"]) * y56, "zero": 0.0,
            "b7_centred": float(p56["centred"]["R_perm_mean"]) * y56,
            "b5": float(p47["R_perm_mean"]) * float(p47["yardstick"])}


def _read(gmat: np.ndarray, grid: np.ndarray, kind: str, target: float) -> tuple[float | None, str]:
    """利得の行列 [δ, 回] の平均の曲線を読む。kind: bound（target = U90 − ΔBS_null を下から越える最大の δ）・
    reach（target に初めて届く δ）。"""
    g = np.concatenate([[0.0], gmat.mean(axis=1)])
    if kind == "bound":
        return v52._cross(g - target, grid)
    return v52._first_reach(g, target, grid)


def _mc(gmat: np.ndarray, grid: np.ndarray, kind: str, target: float, rng: np.random.Generator) -> dict:
    """各 δ の中で回を復元抽出して読み直したオッズ比の分布（届かない再標本は +∞）。"""
    nd, nr = gmat.shape
    vals = np.empty(B_MC)
    for b in range(B_MC):
        idx = rng.integers(0, nr, size=(nd, nr))
        x, how = _read(np.take_along_axis(gmat, idx, axis=1), grid, kind, target)
        vals[b] = np.exp(x) if x is not None else (np.inf if how == "above grid" else np.nan)
    ok = vals[~np.isnan(vals)]
    lo, hi = (float(np.quantile(ok, 0.025)), float(np.quantile(ok, 0.975))) if len(ok) else (np.nan, np.nan)
    return {"mc95": [lo, hi], "share_above_grid": float(np.mean(np.isinf(vals))),
            "share_undefined": float(np.mean(np.isnan(vals)))}


def _entry(gmat: np.ndarray, grid: np.ndarray, kind: str, target: float, or_max: float,
           rng: np.random.Generator) -> dict:
    x, how = _read(gmat, grid, kind, target)
    e = {"delta": x, "odds_ratio": float(np.exp(x)) if x is not None else None, "how": how,
         "odds_ratio_text": _or_text(x, how, or_max)} | _mc(gmat, grid, kind, target, rng)
    e["mc95_text"] = [f"{v:.2f}" if np.isfinite(v) else f"> {or_max:.0f}" for v in e["mc95"]]
    return e


def b1() -> dict | str:
    try:
        deltas, gains, info = plasmode_gains()
    except FileNotFoundError as e:
        return f"missing ({e})"
    lex = _json(LEXFIX / "metrics_primary.json")["primary"]
    yard, u90 = float(lex["delta_bs_M0_M1"]), float(lex["ci90_delta"][1])
    se = (float(lex["ci95_delta"][1]) - float(lex["ci95_delta"][0])) / 3.92
    margin = wa.SESOI * yard
    refs = references()
    grid = np.concatenate([[0.0], deltas])
    or_max = float(np.exp(deltas[-1]))
    rng = np.random.default_rng(SEED_MC)
    out: dict = {"deltas": deltas.tolist(), "yardstick_lexfix": yard, "U90": u90, "se_normal": se, "margin": margin,
                 "references": refs, "reference_primary": "b7_recomputed", "mc_resamples": B_MC, "seed": SEED_MC,
                 "or_grid_max": or_max} | info
    out["gain_mean"] = {f: g.mean(axis=1).tolist() for f, g in gains.items()}
    out["gain_mc_se"] = {f: (g.std(axis=1, ddof=1) / np.sqrt(g.shape[1])).tolist() for f, g in gains.items()}
    out["bounds"] = {ref: {f: _entry(g, grid, "bound", u90 - d, or_max, rng) for f, g in gains.items()}
                     for ref, d in refs.items()}
    d7 = refs["b7_recomputed"]
    out["mde80"] = {f: _entry(g, grid, "reach", Z80 * se - d7, or_max, rng) for f, g in gains.items()}
    out["margin_or"] = {f: _entry(g, grid, "reach", margin, or_max, rng) for f, g in gains.items()}
    out["bounds_10reps"] = {ref: {f: _entry(gains[f][:, :10], grid, "bound", u90 - d, or_max, rng)
                                  for f in C1_FORMS} for ref, d in refs.items()}
    pr = _json(VER / "v56_b7_cluster.json").get("plasmode_refit")
    if isinstance(pr, dict):  # B-7 の当て直し（形 A、3 回）を同じ基準で
        dr = np.asarray(pr["deltas"], dtype=np.float64)
        gr = np.asarray(pr["gain_refit_by_rep"], dtype=np.float64)
        grid_r = np.concatenate([[0.0], dr])
        omax_r = float(np.exp(dr[-1]))
        out["refit_A"] = {ref: _entry(gr, grid_r, "bound", u90 - d, omax_r, rng) for ref, d in refs.items()}
    return out


def c1_check(run: Path = C1RUN) -> dict | str:
    """fold 0 の回 0〜9 を回し直した値と B-6 の値の差、追記 4 の許容差（B-6 の 10 回の ΔBS の標準誤差の 0.1 倍）。"""
    f1 = run / "predictions" / "primary" / "fold0_plasmode.npz"
    f6 = B6PLAS / "predictions" / "primary" / "fold0_plasmode.npz"
    if not f1.is_file():
        return f"missing (fold 0 reps 0-9 rerun in {run.name})"
    a, b = np.load(f1), np.load(f6)
    fa, fb = [str(x) for x in a["forms"]], [str(x) for x in b["forms"]]
    out: dict = {"forms": fa}
    for key in ("sse_c", "sse_n"):
        x = np.stack([a[key][fa.index(f)] for f in fa])
        y = np.stack([b[key][fb.index(f)] for f in fa])
        out[f"max_abs_{key}"] = float(np.abs(x - y).max())
        out[f"max_rel_{key}"] = float((np.abs(x - y) / np.abs(y)).max())
    for key in ("base_c", "base_n"):
        out[f"abs_{key}"] = float(abs(float(a[key]) - float(b[key])))
    db = lambda z, f, fs: (z["sse_c"][fs.index(f)] - z["sse_n"][fs.index(f)]) / float(z["n_te"])
    out["max_abs_delta_bs_fold0"] = float(max(np.abs(db(a, f, fa) - db(b, f, fb)).max() for f in fa))
    worst = 0.0
    for f in fa:  # 追記 4: δ ごとの 10 回の平均の差 / B-6 の 10 回の標準誤差
        d1, d6 = db(a, f, fa), db(b, f, fb)
        se = d6.std(axis=1, ddof=1) / np.sqrt(d6.shape[1])
        worst = max(worst, float(np.max(np.abs(d1.mean(axis=1) - d6.mean(axis=1)) / np.where(se > 0, se, np.inf))))
    out["max_mean_diff_over_se"] = worst
    out["reproduced"] = bool(worst < 0.1)
    return out


# %% C2: 二段目を bootstrap の中で当て直した区間
def c2() -> dict | str:
    d = C2BOOT / "predictions" / "primary"
    files = sorted(d.glob("boot_*.npz")) if d.is_dir() else []
    if not files:
        return "missing (no bootstrap files)"
    acc = {k: np.concatenate([np.load(f)[k] for f in files]) for k in
           ("replicates", "refit_c", "refit_n", "fixed_c", "fixed_n", "m0", "m1", "w")}
    reps = acc["replicates"]
    if len(reps) != C2_REPS or len(np.unique(reps)) != C2_REPS:
        return f"incomplete ({len(np.unique(reps))} of {C2_REPS} replicates)"
    dd = C2DUMP / "predictions" / "primary"
    sse = {"M0": 0.0, "M1": 0.0, "M1c": 0.0, "M3n": 0.0}
    n = 0
    for k in range(N_FOLDS):
        z = np.load(dd / f"fold{k}_stage2.npz")
        y = z["y_te"].astype(np.float64)
        for a in sse:
            sse[a] += float(((z[f"p_{a}_te"].astype(np.float64) - y) ** 2).sum())
        n += len(y)
    yard = (sse["M0"] - sse["M1"]) / n
    r_hat = (sse["M1c"] - sse["M3n"]) / n / yard
    lex = _json(LEXFIX / "metrics_primary.json")
    out: dict = {"n_cells": n, "replicates": C2_REPS, "yardstick_local": yard, "R_local": r_hat,
                 "yardstick_lexfix": float(lex["primary"]["delta_bs_M0_M1"]), "R_lexfix": float(lex["primary"]["R"]),
                 "brier_local": {a: v / n for a, v in sse.items()},
                 "brier_lexfix": {a: float(lex["brier"][a]) for a in sse}}
    checks = [_json(f.with_suffix(".json")).get("check_max_abs_refit_vs_saved") for f in files]
    out["max_abs_unweighted_refit_vs_saved"] = float(np.max(np.asarray(checks, dtype=np.float64)))
    for name, c, nn in (("refit", "refit_c", "refit_n"), ("fixed", "fixed_c", "fixed_n")):
        r_b = (acc[c] - acc[nn]) / acc["w"] / yard
        out[name] = {"R": r_hat, "ci90_R": wa.ci(r_b, 0.90), "ci95_R": wa.ci(r_b, 0.95), "sd_R": float(r_b.std(ddof=1))}
    lo, hi = out["refit"]["ci90_R"]
    out["rule"] = ("equivalent-null even with the variability of refitting" if lo > -wa.SESOI and hi < wa.SESOI
                   else "the fixed-model interval understated the uncertainty")
    out["sd_ratio_refit_to_fixed"] = out["refit"]["sd_R"] / out["fixed"]["sd_R"]
    return {"corrected": out}


# %% C3: 余白での判定の操作特性
def _wilson(k: int, n: int) -> list[float]:
    z = stats.norm.ppf(0.975)
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(c - h), float(c + h)]


def _c3_folds() -> tuple[list[dict], list[int]]:
    """fold ごとに回の範囲のファイル（fold{k}_plasref_<最初の回>.npz）を集め、全ての fold で揃った回だけを返す
    （途中までの回を報告できるように、回を範囲に分けて回す）。"""
    d = C3RUN / "predictions" / "primary"
    per_fold: list[dict[int, tuple[np.lib.npyio.NpzFile, int]]] = []
    for k in range(N_FOLDS):
        reps: dict = {}
        for f in sorted(d.glob(f"fold{k}_plasref_*.npz")) if d.is_dir() else []:
            z = np.load(f)
            for jr, rep in enumerate(int(x) for x in z["rep_ids"]):
                reps[rep] = (z, jr)
        per_fold.append(reps)
    common = sorted(set.intersection(*(set(r) for r in per_fold))) if per_fold else []
    return per_fold, [r for r in common if C3_FIRST <= r < C3_FIRST + C3_REPS]


def c3() -> dict | str:
    per_fold, reps = _c3_folds()
    if not reps:
        return f"missing (replicates complete in all folds: 0; by fold {[len(r) for r in per_fold]})"
    z0 = per_fold[0][reps[0]][0]
    models = [str(x) for x in z0["models"]]
    n = np.concatenate([per_fold[k][reps[0]][0]["n_unit"] for k in range(N_FOLDS)]).astype(np.float64)
    idx = {a: j for j, a in enumerate(models)}
    rows = []
    for rep in reps:
        s = np.concatenate([per_fold[k][rep][0]["refit_unit"][0, per_fold[k][rep][1]] for k in range(N_FOLDS)],
                           axis=-1).T  # [単位, model]
        hat = s.sum(axis=0) / n.sum()
        delta = wa.make_delta(hat, wa.bootstrap_bs(s, n), idx)
        d01, b01 = delta("M0", "M1")
        v, _ = wa.verdict(delta, "M1c", "M3n", None, d01=d01, b01=b01)
        rows.append({"rep": rep, "R": v["R"], "ci90_R": v["ci90_R"], "decision": v["decision"],
                     "equivalent": bool(v["ci90_R"][0] > -wa.SESOI and v["ci90_R"][1] < wa.SESOI)})
    k = sum(r["equivalent"] for r in rows)
    wil = _wilson(k, len(rows))
    dec: dict[str, int] = {}
    for r in rows:
        dec[r["decision"]] = dec.get(r["decision"], 0) + 1
    return {"delta": float(z0["deltas"][0]), "odds_ratio": float(np.exp(z0["deltas"][0])), "rep_ids": reps,
            "reps_by_fold": [len(r) for r in per_fold],
            "reps": len(rows), "n_equivalent": k, "share_equivalent": k / len(rows), "wilson95": wil,
            "decisions": dec, "R_mean": float(np.mean([r["R"] for r in rows])), "by_rep": rows,
            "rule": ("the fixed-model interval is anti-conservative at the margin" if wil[0] > 0.05
                     else "no evidence that the fixed-model interval is anti-conservative at the margin")}


# %% C4・C5: 評価の実行
def _eval_run(run: Path, pairs: list) -> dict | str:
    if not v56._complete(run):
        return "missing"
    data = wa.load_split(run, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    out = {"run": run.name, "n_aoi": len(np.unique(data["aoi"]))} | v56._verdicts(y, preds, data["unit"], pairs)
    if v56._complete(LEXFIX):
        ref = wa.load_split(LEXFIX, "primary")
        same = len(ref["y"]) == len(y) and bool(np.array_equal(ref["y"], data["y"]))
        out["check_vs_lexfix"] = {"same_cells": same} | (
            {f"max_abs_diff_{a}": float(np.max(np.abs(ref[f"p_{a}"].astype(float) - preds[a].astype(float))))
             for a in ("M1c", "M3n", "M0", "M1") if a in preds and f"p_{a}" in ref} if same else {})
    return out


def c4() -> dict | str:
    out = _eval_run(C4RUN, [("gsw_cell", "M1cW", "M3oW", None), ("primary", "M1c", "M3n", None)])
    if not isinstance(out, dict):
        return out
    g6 = (_json(VER / "v52_b6_cluster.json").get("gsw_oracle") or {}).get("oracle_gsw") or {}
    r_cell, r_set = out["gsw_cell"]["R"], g6.get("R")
    out["gsw_settlement_b6"] = {"R": r_set, "ci90_R": g6.get("ci90_R"), "decision": g6.get("decision")}
    if r_set is None:
        out["rule"] = "missing (B-6 settlement GSW)"
    elif r_cell >= wa.SESOI > r_set:
        out["rule"] = "assignment is the bottleneck for GSW too"
    elif r_cell < wa.SESOI and r_set < wa.SESOI:
        out["rule"] = "recurrent water carries little about mapped extent once terrain is known, however assigned"
    else:
        out["rule"] = "neither reading of the rule applies"
    return out


def c5() -> dict | str:
    out = _eval_run(C5RUN, [("lexclean", "M1c", "M3n", "M3np")])
    if isinstance(out, dict):
        lex = _json(LEXFIX / "metrics_primary.json")["primary"]
        out["lexfix"] = {"R": float(lex["R"]), "ci90_R": lex.get("ci90_R")}
        rec = _json(C5RUN / "predictions" / "primary" / "fold0.json")
        out["lexicon"] = rec.get("lexicon")
    return out


def main() -> None:
    res: dict = {"rule": "amendment B-8 items B1, C1-C5 (addendum 3: local runs)"}
    res["B1"] = b1()
    res["C1_check"] = c1_check()
    res["C1_check_local"] = c1_check(C1LOCAL)
    res["C2"] = c2()
    res["C3"] = c3()
    res["C4"] = c4()
    res["C5"] = c5()
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    b = res["B1"]
    if isinstance(b, dict):
        print("B1 reps", b["reps"], "missing C1 ranges", b["c1_ranges_missing"])
        for ref in REFS:
            print(f"  bound [{ref:13s}]", {f: (e["odds_ratio_text"], e["mc95_text"], round(e["share_above_grid"], 3))
                                           for f, e in b["bounds"][ref].items()})
        print("  mde80", {f: (e["odds_ratio_text"], e["mc95_text"]) for f, e in b["mde80"].items()})
        print("  margin", {f: (e["odds_ratio_text"], e["mc95_text"]) for f, e in b["margin_or"].items()})
        if "refit_A" in b:
            print("  refit A", {r: (e["odds_ratio_text"], e["mc95_text"]) for r, e in b["refit_A"].items()})
    else:
        print("B1", b)
    print("C1 check", res["C1_check"])
    print("C1 check (local)", res["C1_check_local"])
    for key in ("C2", "C3"):
        v = res[key]
        print(key, json.dumps({k: x for k, x in v.items() if k != "by_rep"} if isinstance(v, dict) else v,
                              default=float)[:900])
    for key, name in (("C4", "gsw_cell"), ("C5", "lexclean")):
        v = res[key]
        if isinstance(v, dict):
            c = v[name]
            print(f"{key} {name} R {c['R']:+.4f} 90% [{c['ci90_R'][0]:+.4f}, {c['ci90_R'][1]:+.4f}] {c['decision']} "
                  f"rule: {v.get('rule', '')} check {v.get('check_vs_lexfix')}")
        else:
            print(key, v)


if __name__ == "__main__":
    main()
