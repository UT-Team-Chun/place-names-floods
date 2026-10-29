# %% v67（事後・amendment B-10 の C1〜C5）: クラスタの実行の集計
"""確証判定は変えない。B-10 で計算の前に固定した量を、各実行の出力から計算する（v64 と同じ道具）。

- C1: 名前の特徴 × 距離帯（M1cd → M3nd）、× 大陸（M1cr → M3nr）。world_v2_b10stk
- C2: 二段目の C を入れ子の交差検証で選んだ比較（M1ct → M3nt）と、fold ごとに選んだ C
- C3: 国の指標を入れた基準（M1kc → M3kn）と、国の指標の後の地域構成の増分（M1k → M1kc）
- C4: 形 AC（反転を集落の地面の中で相関させた注入、ρ = 0.5）の上限と余白に当たる OR。world_v2_b10plasac
- C5: JRC の洪水ハザード図を足した地形のモデル（物差し M0 → M1h、M1hc → M3hn）。world_v2_b10haz

判定は主比較と同じ判定表（確証の判定は置き換えない）。出力: stdout と docs/verification/v67_b10_aggregate.json
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
sys.path.insert(0, str(HERE))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v67_b10_aggregate.json"
STK = RUNS_DIR / "world_v2_b10stk"
HAZ = RUNS_DIR / "world_v2_b10haz"
PLASAC = RUNS_DIR / "world_v2_b10plasac"
SEED = 20261030
N_FOLDS = 10


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v64 = _load("v64_b9_aggregate")
v59, v56 = v64.v59, v64.v56


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _folds(run: Path) -> list[dict]:
    d = run / "predictions" / "primary"
    return [_json(d / f"fold{k}.json") for k in range(N_FOLDS) if (d / f"fold{k}.json").is_file()]


def c1_c3() -> dict | str:
    if not v56._complete(STK):
        return "missing"
    data = wa.load_split(STK, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    pairs = [("C1_dist", "M1cd", "M3nd", None), ("C1_cont", "M1cr", "M3nr", None),
             ("C2_tuned", "M1ct", "M3nt", None), ("C3_country", "M1kc", "M3kn", None),
             ("C3_comp", "M1k", "M1kc", None), ("primary_check", "M1c", "M3n", None)]
    out = {"run": STK.name, "n_aoi": len(np.unique(data["aoi"]))} | v56._verdicts(y, preds, data["unit"], pairs)
    recs = _folds(STK)
    out["C2_chosen_C"] = {str(r["fold"]): (r.get("b10_tuned") or {}).get("M3nt", {}).get("C") for r in recs}
    out["C2_chosen_C_reference"] = {str(r["fold"]): (r.get("b10_tuned") or {}).get("M1ct", {}).get("C") for r in recs}
    out["C3_countries_in_train"] = {str(r["fold"]): len(r.get("b10_countries_in_train") or []) for r in recs}
    if isinstance(out.get("C3_country"), dict) and isinstance(out.get("C3_comp"), dict):
        out["C3_country"]["composition_after_country"] = out["C3_comp"]
    lex = _json(RUNS_DIR / "world_v2_lexfix" / "metrics_primary.json")["primary"]
    pc = out.get("primary_check") or {}
    out["check_vs_lexfix"] = {"R_this_run": pc.get("R"), "R_lexfix": float(lex["R"])}
    return out


def c4_ac() -> dict | str:
    r = v64.run_gains(PLASAC, ("AC",))
    if r is None:
        return "missing"
    deltas, gains, done = r
    refs, iv = v59.references(), v64.intervals()
    rng = np.random.default_rng(SEED)
    grid = np.concatenate([[0.0], deltas])
    g = gains["AC"]
    om = float(np.exp(deltas[-1]))
    bound = {ref: v59._entry(g, grid, "bound", iv["fixed"]["U90"] - d, om, rng)
             for ref, d in refs.items()}
    d = PLASAC / "predictions" / "primary"
    recs = [_json(f) for f in sorted(d.glob("fold*_plasmode*.json"))]
    argv = (recs[0].get("argv") if recs else None) or []
    rho = float(argv[argv.index("--plasmode-rho") + 1]) if "--plasmode-rho" in argv else None
    return {"deltas": deltas.tolist(), "ranges": done, "reps": int(g.shape[1]),
            "gain_mean": g.mean(axis=1).tolist(), "gain_mc_se": (g.std(axis=1, ddof=1) / np.sqrt(g.shape[1])).tolist(),
            "bound": bound, "bound_text": bound["b7_recomputed"]["odds_ratio_text"],
            "margin_or": v59._entry(g, grid, "reach", wa.SESOI * iv["yardstick"]["point"], om, rng),
            "rho": rho}


def c5_hazard() -> dict | str:
    if not v56._complete(HAZ):
        return "missing"
    data = wa.load_split(HAZ, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    if "M1h" not in preds:
        return "missing (no M1h predictions)"
    preds["M1"] = preds.pop("M1h")  # 物差しは M0 → M1h（B-10 の C5）
    out = v56._verdicts(y, preds, data["unit"], [("hazard", "M1hc", "M3hn", None)])
    h = out.get("hazard")
    if not isinstance(h, dict):
        return "missing (comparison)"
    rec = (_folds(HAZ) or [{}])[0]
    lex = _json(RUNS_DIR / "world_v2_lexfix" / "metrics_primary.json")["primary"]
    return h | {"yardstick": out["yardstick"], "ci95_yardstick": out["ci95_yardstick"], "brier": out["brier"],
                "yardstick_without_hazard": float(lex["delta_bs_M0_M1"]),
                "hazard_record": (rec.get("b6") or {}).get("jrc_hazard"), "n_cells": out["n_cells"]}


def main() -> None:
    res: dict = {"rule": "amendment B-10 part C (exploratory; confirmatory verdict unchanged)"}
    c = c1_c3()
    if isinstance(c, dict):
        for k in ("C1_dist", "C1_cont", "C2_tuned", "C3_country"):
            res[k] = c.get(k)
        res["stack_run"] = {k: v for k, v in c.items() if k not in ("C1_dist", "C1_cont", "C2_tuned", "C3_country")}
    else:
        res["stack_run"] = c
    res["C4_ac"] = c4_ac()
    res["C5_hazard"] = c5_hazard()
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for k in ("C1_dist", "C1_cont", "C2_tuned", "C3_country", "C5_hazard"):
        v = res.get(k)
        if isinstance(v, dict):
            print(k, f"R={v['R']:.4f}", "ci90", [round(x, 4) for x in v["ci90_R"]], v["decision"], flush=True)
        else:
            print(k, v, flush=True)
    c4 = res["C4_ac"]
    print("C4", c4 if isinstance(c4, str) else (c4["bound_text"], c4["reps"]), flush=True)


if __name__ == "__main__":
    main()
