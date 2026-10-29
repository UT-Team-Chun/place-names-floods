# %% v74（事後・amendment B-14 の C1〜C3）: クラスタの実行の集計
"""確証判定は変えない。B-14 の追記 1（実装）と追記 3（集計）で計算の前に固定した量。

- C1: 観測できない地面（10 m の市街地と樹冠の割合の和が 0.5 以上）を除いた EMS のラベルでの主比較と判定表、
  その規則で絞った上限ありの表での集落の近くの OR
- C2: H3 の陽性対照（地点の HAND を集落を通じて、T0 と M1 に）
- C3: H4 の全体の増分の国の構成（言語モデルの切れ目の分割で M1s → M3rs と M1k → M3rk）

出力: stdout と docs/verification/v74_b14_cluster.json（部分ごとに足す。`--only C1` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v74_b14_cluster.json"
SEED = 20261103
B_FIXED = 4000
OBSERVABLE_MAX = 0.5
RUNS = {
    "C1_observable": ("world_v2_b14obs", "primary", [("names", "M1c", "M3n", "M3np")]),
    "C2_oracle": ("world_v2_b14oracle", "primary", [("site_hand_T0", "T0co", "T0o", None),
                                                    ("site_hand_M1", "M1co", "M3o", None)]),
    "C3_h4_country": ("world_v2_b14h4k", "llmcut", [("nread_total", "M1s", "M3rs", None),
                                                    ("nread_country", "M1k", "M3rk", None)]),
}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v71 = _load("v71_b13_cluster")
v56, wa = v71.v56, v71.wa
v70 = _load("v70_b13_local")
v69 = v70.v69


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _test_a(y: np.ndarray, preds: dict, unit: np.ndarray, base: str, alt: str, rng: np.random.Generator) -> dict:
    d = {"y": y, "unit": unit} | {f"p_{a}": preds[a] for a in ("M0", "M1", base, alt)}
    n, s, units = v69._group_sums(d, ("M0", "M1", base, alt))
    eg, dg = s[base] - s[alt], s["M0"] - s["M1"]
    big_n = float(n.sum())
    w = rng.multinomial(len(units), np.full(len(units), 1.0 / len(units)), size=B_FIXED).astype(np.float64)
    nw = w @ n
    return v70._summ(float(eg.sum() / big_n), float(dg.sum() / big_n), (w @ eg) / nw, (w @ dg) / nw)


def _eval(key: str) -> dict | str:
    tag, split, pairs = RUNS[key]
    run = RUNS_DIR / tag
    files = sorted((run / "predictions" / split).glob("fold*.npz"))
    if not files or (split == "primary" and len(files) < 10):
        return "missing"
    data = wa.load_split(run, split)
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    out = {"run": tag, "split": split, "n_aoi": len(np.unique(data["aoi"]))} | v56._verdicts(y, preds, data["unit"],
                                                                                                 pairs)
    rng = np.random.default_rng(SEED)
    for name, base, alt, _ in pairs:
        if base in preds and alt in preds:
            out[f"{name}_test_a"] = _test_a(y, preds, np.asarray(data["unit"]), base, alt, rng)
    rec = _json(run / "predictions" / split / "fold0.json")
    out["labels_record"] = rec.get("labels")
    out["lexicon"] = rec.get("lexicon")
    return out


def _lc_cover(t: pd.DataFrame) -> np.ndarray:
    """各セルの features の lc_built + lc_tree（欠けは 0、(aoi, ix, iy) で結ぶ）。"""
    out = np.full(len(t), np.nan)
    for a, ia in t.groupby("aoi").indices.items():
        f = pq.read_table(v69.LABELS / a / "features.parquet", columns=["lc_built", "lc_tree"]).to_pandas()
        c = pq.read_table(v69.LABELS / a / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        key = pd.MultiIndex.from_arrays([c["ix"].to_numpy(), c["iy"].to_numpy()])
        cov = np.nan_to_num(f["lc_built"].to_numpy(np.float64) + f["lc_tree"].to_numpy(np.float64), nan=0.0)
        s = pd.Series(cov, index=key)
        q = pd.MultiIndex.from_arrays([t["ix"].to_numpy()[ia].astype(np.int64),
                                       t["iy"].to_numpy()[ia].astype(np.int64)])
        out[ia] = s.reindex(q).to_numpy()
    return out


def c1_near() -> dict:
    v70b = _load("v70b_b13_part2")
    v55 = v70b.v69.v55
    gn = v70b._gn()
    t = pd.read_parquet(v70b.VER / "v37_table.parquet")
    t41 = pd.read_parquet(v70b.VER / "v41_table.parquet", columns=["aoi", "y", "ix", "iy"])
    if not (np.array_equal(t["aoi"].to_numpy(), t41["aoi"].to_numpy())
            and np.array_equal(t["y"].to_numpy(), t41["y"].to_numpy())):
        return {"note": "v37 and v41 tables are not aligned"}
    t["ix"], t["iy"] = t41["ix"].to_numpy(), t41["iy"].to_numpy()
    cov = _lc_cover(t)
    obs = np.nan_to_num(cov, nan=0.0) < OBSERVABLE_MAX
    g = t["P_gid1"].to_numpy(np.int64)
    fl = v71._corrected_flags(g, gn)
    tc = t.assign(lex_wetside_n1=fl.reindex(g).fillna(0.0).to_numpy())
    return {"share_cells_observable": float(obs.mean()), "n_cells_missing_cover": int(np.isnan(cov).sum()),
            "corrected_reading": {"all_cells": v55.odds(tc), "observable": v55.odds(tc, mask=obs)},
            "frozen_reading": {"all_cells": v55.odds(t), "observable": v55.odds(t, mask=obs)},
            "rule": "B-14 C1 (addendum 3): Sect. 5 estimator on the capped table of v37, restricted to cells with "
                    "lc_built + lc_tree < 0.5; corrected-lexicon reading as in B-13 C2, frozen reading alongside"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    want = {x for x in args.only.split(",") if x}

    def run(key: str, fn: callable) -> None:
        if want and key not in want:
            return
        t0 = time.time()
        val = fn()
        res = _json(OUT) if OUT.is_file() else {}
        res |= {"rule": "amendment B-14 (post hoc), cluster runs C1 to C3", "seed": SEED, key: val}
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        print(f"{key} ({time.time() - t0:.0f} s) {json.dumps(val, default=float)[:700]}", flush=True)

    for key in RUNS:
        run(key, lambda key=key: _eval(key))
    run("C1_near", c1_near)


if __name__ == "__main__":
    main()
