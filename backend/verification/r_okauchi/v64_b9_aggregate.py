# %% v64（事後・amendment B-9 の B1・B2 と C1〜C7）: クラスタの実行の集計と、上限の読み直し
"""確証判定は変えない。B-9 で計算の前に固定した量を、各実行の出力から計算する。

- B1: 上限と検出力 80% の OR を、観測の増分の 90% 上限を (i) モデルを固定した区間（B-8 の B1 と同じ）、(ii) 二段目を
  bootstrap の中で当て直した区間（B-8 の C2）、(iii) 13 通りの fold の割り当てを合わせた区間（v47）に替えて読む。
  曲線は B-8 の B1 の 50 回（world_v2_b6plasmode + world_v2_b8plas50）に、形 A・CA は C1 の低い δ（world_v2_b9plaslow）を
  足した格子。無効果の基準は B-8 と同じ 4 つ。モンテカルロの区間は B-8 の B1 と同じ（種 20261029）
- B2: 物差しの 95% 区間の両端と点推定の 5% を余白として、形 A の曲線で余白に当たる OR
- C2: 浸水を減らす向きの注入（world_v2_b9plasneg、形 A・B）の下限（OR が x 未満の効果を除く）
- C3: 集落ごとに違う効果（world_v2_b9plashet、形 AH）の上限と余白に当たる OR
- C4: 全てのモデルを当て直した模擬の検出力（world_v2_b9power、OR 1.25・1.5・2.0、回 200〜239）
- C5〜C7: 共通の永久水域（world_v2_b9gswperm）、上限の抽出の種（world_v2_b9seed1〜3）、FABDEM（world_v2_b9fabdem）の
  主比較。C5・C7 は集落の近くの OR（B-8 の B6 の主の式、pairs bootstrap 2,000 回）、C7 は浸水の地図の技能も

出力: stdout と docs/verification/v64_b9_aggregate.json（無い入力は "missing ..."）
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v64_b9_aggregate.json"
VER = ROOT / "docs" / "verification"
LEXFIX = RUNS_DIR / "world_v2_lexfix"
LOW = RUNS_DIR / "world_v2_b9plaslow"
NEG = RUNS_DIR / "world_v2_b9plasneg"
HET = RUNS_DIR / "world_v2_b9plashet"
POWER = RUNS_DIR / "world_v2_b9power"
GSW = RUNS_DIR / "world_v2_b9gswperm"
SEEDS = [RUNS_DIR / f"world_v2_b9seed{s}" for s in (1, 2, 3)]
FAB = RUNS_DIR / "world_v2_b9fabdem"
SEED = 20261029
N_FOLDS = 10
FIRSTS = (0, 10, 20, 30, 40)
POWER_FIRST, POWER_REPS = 200, 40
POWER_SHARE = 0.8


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v59 = _load("v59_b8_aggregate")
v56 = v59.v56


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# %% 曲線
def run_gains(run: Path, forms: tuple[str, ...]) -> tuple[np.ndarray, dict[str, np.ndarray], list[int]] | None:
    """回の範囲ごとの 10 fold のファイルから形ごとの利得 [δ, 回]（全ての fold で揃った範囲だけ）。"""
    d = run / "predictions" / "primary"
    parts, done = [], []
    deltas = None
    for first in FIRSTS:
        stem = "plasmode" if first == 0 else f"plasmode_{first:03d}"
        fs = [d / f"fold{k}_{stem}.npz" for k in range(N_FOLDS)]
        if not all(f.is_file() for f in fs):
            continue
        z0 = np.load(fs[0])
        deltas = np.asarray(z0["deltas"], dtype=np.float64)
        parts.append(v59._fold_gains(fs, forms))  # noqa: SLF001
        done.append(first)
    if not parts:
        return None
    gains = {f: np.concatenate([p[f] for p in parts], axis=1) for f in forms if f in parts[0]}
    return deltas, gains, done


def combined_curves() -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    """B-8 の B1 の曲線（7 点 × 50 回）に、形 A・CA は C1 の低い δ（4 点 × 50 回）を足した格子。"""
    d6, g6, info = v59.plasmode_gains()
    low = run_gains(LOW, ("A", "CA"))
    out_g, out_d = {}, {}
    info = dict(info) | {"low_ranges": low[2] if low else []}
    for f, g in g6.items():
        if low and f in low[1] and low[1][f].shape[1] == g.shape[1]:
            d = np.concatenate([low[0], d6])
            gg = np.concatenate([low[1][f], g], axis=0)
            o = np.argsort(d)
            out_d[f], out_g[f] = d[o], gg[o]
        else:
            out_d[f], out_g[f] = d6, g
    info["grid"] = {f: out_d[f].tolist() for f in out_d}
    return out_d, out_g, info


# %% B1・B2
def intervals() -> dict[str, dict[str, float]]:
    """観測の増分の 90% 上限（Brier の尺度）と標準誤差: 固定・当て直し・fold の割り当て。"""
    lex = _json(LEXFIX / "metrics_primary.json")["primary"]
    yard = float(lex["delta_bs_M0_M1"])
    u90 = float(lex["ci90_delta"][1])
    se = (float(lex["ci95_delta"][1]) - float(lex["ci95_delta"][0])) / 3.92
    c2 = _json(VER / "v59_b8_aggregate.json")["C2"]["corrected"]
    y2 = float(c2["yardstick_local"])
    fs = _json(VER / "v47_b5_cluster.json")["fold_seeds"]["fixed_code"]
    return {"fixed": {"U90": u90, "se": se},
            "refit": {"U90": float(c2["refit"]["ci90_R"][1]) * y2, "se": float(c2["refit"]["sd_R"]) * y2},
            "fold": {"U90": float(fs["combined_ci90"][1]) * yard, "se": se}, "yardstick": {"point": yard}}


def b1_b2() -> dict:
    grids, gains, info = combined_curves()
    refs = v59.references()
    iv = intervals()
    rng = np.random.default_rng(SEED)
    out: dict = {"intervals": iv, "references": refs, "reference_primary": "b7_recomputed", "seed": SEED} | info
    # 図 6 の曲線（格子の点ごとの平均の利得とモンテカルロの標準誤差、形ごとの格子）
    out["gain_mean"] = {f: g.mean(axis=1).tolist() for f, g in gains.items()}
    out["gain_mc_se"] = {f: (g.std(axis=1, ddof=1) / np.sqrt(g.shape[1])).tolist() for f, g in gains.items()}
    out["bounds"] = {}
    for kind in ("fixed", "refit", "fold"):
        out["bounds"][kind] = {}
        for ref, dnull in refs.items():
            out["bounds"][kind][ref] = {}
            for f, g in gains.items():
                grid = np.concatenate([[0.0], grids[f]])
                or_max = float(np.exp(grids[f][-1]))
                out["bounds"][kind][ref][f] = v59._entry(g, grid, "bound", iv[kind]["U90"] - dnull, or_max, rng)  # noqa: SLF001
    d7 = refs["b7_recomputed"]
    out["mde80"] = {}
    for kind in ("fixed", "refit"):
        out["mde80"][kind] = {}
        for f, g in gains.items():
            grid = np.concatenate([[0.0], grids[f]])
            out["mde80"][kind][f] = v59._entry(g, grid, "reach", v59.Z80 * iv[kind]["se"] - d7,  # noqa: SLF001
                                               float(np.exp(grids[f][-1])), rng)
    # B2: 物差しの不確かさと余白の換算（形 A）
    lex = _json(LEXFIX / "metrics_primary.json")["primary"]
    lo, hi = (float(x) for x in lex["ci95_delta_M0_M1"])
    grid = np.concatenate([[0.0], grids["A"]])
    out["margin_by_yardstick"] = {
        q: v59._entry(gains["A"], grid, "reach", wa.SESOI * y, float(np.exp(grids["A"][-1])), rng)  # noqa: SLF001
        for q, y in (("q025", lo), ("point", iv["yardstick"]["point"]), ("q975", hi))}
    out["yardstick_ci95"] = [lo, hi]
    return out


# %% C2・C3
def c2_negative() -> dict | str:
    r = run_gains(NEG, ("A", "B"))
    if r is None:
        return "missing"
    deltas, gains, done = r
    refs, iv = v59.references(), intervals()
    rng = np.random.default_rng(SEED + 1)
    o = np.argsort(np.abs(deltas))
    ad = np.abs(deltas)[o]
    grid = np.concatenate([[0.0], ad])
    out: dict = {"deltas": deltas.tolist(), "ranges": done, "reps": {f: int(g.shape[1]) for f, g in gains.items()},
                 "gain_mean": {f: g[o].mean(axis=1).tolist() for f, g in gains.items()},
                 "gain_mc_se": {f: (g[o].std(axis=1, ddof=1) / np.sqrt(g.shape[1])).tolist() for f, g in gains.items()},
                 "abs_deltas": ad.tolist(), "lower_bounds": {}}
    for ref, dnull in refs.items():
        out["lower_bounds"][ref] = {}
        for f, g in gains.items():
            e = v59._entry(g[o], grid, "bound", iv["fixed"]["U90"] - dnull, float(np.exp(ad[-1])), rng)  # noqa: SLF001
            x = e["delta"]
            e["odds_ratio"] = float(np.exp(-x)) if x is not None else None
            e["odds_ratio_text"] = f"{np.exp(-x):.2f}" if x is not None else f"< {np.exp(-ad[-1]):.2f}"
            e["mc95_text"] = [f"{np.exp(-np.log(v)):.2f}" if np.isfinite(v) else f"< {np.exp(-ad[-1]):.2f}"
                              for v in e["mc95"][::-1]]
            out["lower_bounds"][ref][f] = e
    return out


def c3_het() -> dict | str:
    r = run_gains(HET, ("AH",))
    if r is None:
        return "missing"
    deltas, gains, done = r
    refs, iv = v59.references(), intervals()
    rng = np.random.default_rng(SEED + 2)
    grid = np.concatenate([[0.0], deltas])
    g = gains["AH"]
    om = float(np.exp(deltas[-1]))
    return {"deltas": deltas.tolist(), "ranges": done, "reps": int(g.shape[1]), "gain_mean": g.mean(axis=1).tolist(),
            "gain_mc_se": (g.std(axis=1, ddof=1) / np.sqrt(g.shape[1])).tolist(),
            "bound": {ref: v59._entry(g, grid, "bound", iv["fixed"]["U90"] - d, om, rng)  # noqa: SLF001
                      for ref, d in refs.items()},
            "margin_or": v59._entry(g, grid, "reach", wa.SESOI * iv["yardstick"]["point"], om, rng)}  # noqa: SLF001


# %% C4
def c4_power() -> dict | str:
    d = POWER / "predictions" / "primary"
    per_fold: list[dict] = []
    for k in range(N_FOLDS):
        reps: dict = {}
        for f in sorted(d.glob(f"fold{k}_plasref_*.npz")) if d.is_dir() else []:
            z = np.load(f)
            for jr, rep in enumerate(int(x) for x in z["rep_ids"]):
                reps[rep] = (z, jr)
        per_fold.append(reps)
    common = sorted(set.intersection(*(set(r) for r in per_fold))) if per_fold and all(per_fold) else []
    reps = [r for r in common if POWER_FIRST <= r < POWER_FIRST + POWER_REPS]
    if not reps:
        return f"missing (replicates complete in all folds: 0; by fold {[len(r) for r in per_fold]})"
    z0 = per_fold[0][reps[0]][0]
    models = [str(x) for x in z0["models"]]
    deltas = [float(x) for x in z0["deltas"]]
    idx = {a: j for j, a in enumerate(models)}
    n = np.concatenate([per_fold[k][reps[0]][0]["n_unit"] for k in range(N_FOLDS)]).astype(np.float64)
    out: dict = {"reps": len(reps), "rep_ids": reps, "by_delta": []}
    for j, delta in enumerate(deltas):
        rows = []
        for rep in reps:
            s = np.concatenate([per_fold[k][rep][0]["refit_unit"][j, per_fold[k][rep][1]] for k in range(N_FOLDS)],
                               axis=-1).T
            hat = s.sum(axis=0) / n.sum()
            delta_fn = wa.make_delta(hat, wa.bootstrap_bs(s, n), idx)
            d01, b01 = delta_fn("M0", "M1")
            v, _ = wa.verdict(delta_fn, "M1c", "M3n", None, d01=d01, b01=b01)
            rows.append({"rep": rep, "R": v["R"], "ci90_R": v["ci90_R"], "ci95_R": v.get("ci95_R"),
                         "decision": v["decision"], "positive": str(v["decision"]).startswith("positive"),
                         "equivalent": bool(v["ci90_R"][0] > -wa.SESOI and v["ci90_R"][1] < wa.SESOI)})
        k_pos = sum(r["positive"] for r in rows)
        k_eq = sum(r["equivalent"] for r in rows)
        out["by_delta"].append({"delta": delta, "odds_ratio": float(np.exp(delta)), "reps": len(rows),
                                "n_positive": k_pos, "share_positive": k_pos / len(rows),
                                "wilson95_positive": v59._wilson(k_pos, len(rows)),  # noqa: SLF001
                                "n_equivalent": k_eq, "share_equivalent": k_eq / len(rows),
                                "R_mean": float(np.mean([r["R"] for r in rows])), "rows": rows})
    hit = [b for b in out["by_delta"] if b["share_positive"] >= POWER_SHARE]
    out["simulated_detectable_or"] = min(b["odds_ratio"] for b in hit) if hit else None
    out["rule"] = ("smallest injected odds ratio detected in at least 80 % of replicates" if hit
                   else f"none of the injected odds ratios up to {max(np.exp(deltas)):.1f} was detected in 80 % "
                   "of replicates")
    return out


# %% C5〜C7: 評価の実行
def _eval(run: Path) -> dict | str:
    out = v59._eval_run(run, [("primary", "M1c", "M3n", "M3np")])  # noqa: SLF001
    if isinstance(out, dict):
        fold0 = run / "predictions" / "primary" / "fold0.json"
        if fold0.is_file():
            rec = _json(fold0)
            out["record"] = {k: rec.get(k) for k in ("b6", "argv", "repo_git_sha")}
    return out


def near_or(run: Path) -> dict | str:
    """B-8 の B6 の主の式（水系側 × 500 m、EMS）を、この実行の held-out の地形の予測で。"""
    if not v56._complete(run):  # noqa: SLF001
        return "missing"
    v55 = v59._load("v55_b7_local")  # noqa: SLF001
    v57 = v59._load("v57_b8_local")  # noqa: SLF001
    t = v55.table()
    data = wa.load_split(run, "primary")
    p = pd.DataFrame({"aoi": data["aoi"].astype(str), "lon": data["lon"].astype(np.float32),
                      "lat": data["lat"].astype(np.float32), "p_run_M1": data["p_M1"].astype(np.float64)})
    tt = t.assign(lon=t["lon"].astype(np.float32), lat=t["lat"].astype(np.float32))
    m = tt.merge(p, on=["aoi", "lon", "lat"], how="inner", validate="one_to_one")
    m["p_M1"] = m["p_run_M1"]
    o = v55.odds(m, cls="wetside")
    key = f"v64_{run.name}"
    b = v57.pairs_draws(m, "wetside", key)
    ok = np.isfinite(b)
    o["pairs_ci95"] = [float(np.exp(np.quantile(b[ok], 0.025))), float(np.exp(np.quantile(b[ok], 0.975)))]
    o["pairs_reps_ok"] = int(ok.sum())
    o["n_cells_matched"] = len(m)
    return o


def map_skill(run: Path) -> dict | str:
    """B の 5 の浸水の地図の技能（v63 と同じ関数）を、この実行の held-out の予測で（C7: FABDEM の地形のモデル）。

    HAND の行は v63 と同じ GLO-30 の HAND（参照）。
    """
    if not v56._complete(run):  # noqa: SLF001
        return "missing"
    v55 = v59._load("v55_b7_local")  # noqa: SLF001
    v63 = v59._load("v63_b9_local")  # noqa: SLF001
    t = v55.table()
    data = wa.load_split(run, "primary")
    p = pd.DataFrame({"aoi": data["aoi"].astype(str), "lon": data["lon"].astype(np.float32),
                      "lat": data["lat"].astype(np.float32)}
                     | {f"p_lf_{a}": data[f"p_{a}"].astype(np.float64) for a in ("M1", "M1c", "M3n")})
    tt = t.drop(columns=[c for c in t.columns if c.startswith("p_lf_")])
    tt = tt.assign(lon=tt["lon"].astype(np.float32), lat=tt["lat"].astype(np.float32))
    m = tt.merge(p, on=["aoi", "lon", "lat"], how="inner", validate="one_to_one")
    out = v63.b5_skill(m, np.random.default_rng(SEED + 7))
    out["n_cells_matched"] = len(m)
    return out


def main() -> None:
    res: dict = {"rule": "amendment B-9 items B1, B2 and C1-C7"}
    res["B1_B2"] = b1_b2()
    res["C2"] = c2_negative()
    res["C3"] = c3_het()
    res["C4"] = c4_power()
    res["C5"] = _eval(GSW)
    res["C6"] = {s.name: _eval(s) for s in SEEDS}
    res["C7"] = _eval(FAB)
    if "--near" in sys.argv:
        res["C5_near"] = near_or(GSW)
        res["C7_near"] = near_or(FAB)
        res["C7_skill"] = map_skill(FAB)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    b = res["B1_B2"]
    for kind in ("fixed", "refit", "fold"):
        e = b["bounds"][kind]["b7_recomputed"]
        print(f"B1 {kind}:", {f: (x["odds_ratio_text"], x["mc95_text"]) for f, x in e.items()}, flush=True)
    print("B1 mde80:", {k: {f: x["odds_ratio_text"] for f, x in v.items()} for k, v in b["mde80"].items()})
    print("B2:", {q: x["odds_ratio_text"] for q, x in b["margin_by_yardstick"].items()})
    for k in ("C2", "C3", "C4", "C5", "C6", "C7"):
        print(k, json.dumps(res[k], default=float)[:700], flush=True)


if __name__ == "__main__":
    main()
