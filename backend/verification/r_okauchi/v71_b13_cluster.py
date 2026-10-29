# %% v71（事後・amendment B-13 の C1〜C4）: クラスタの実行の集計と、直した辞書での集落の近くのオッズ比
"""確証判定は変えない。B-13 の C で計算の前に固定した量。

- C1: 漏れのない対照（world_v2_b13gsw22・world_v2_b13recprior・world_v2_b13jrcset）。セルに直接と集落を通じた $R$ と
  判定表の判定、名の主比較も同じ表で
- C2: 誤りを直した辞書（world_v2_b13lexfix）の主比較と、同じ辞書で読んだ集落の近くのオッズ比（上限ありの EMS・計画した
  GFM・上限なしの EMS、手元で計算）
- C3: 場所グループで選んだ早期終了（world_v2_b13gres）の主比較、物差し、定数の予測との比較、選ばれた反復の回数
- C4: B-12 の C1（1 km 以内への注入、world_v2_b12plas1k、形 BK）の余白と上限

出力: docs/verification/v71_b13_cluster.json（部分ごとに足す。`--only C1,C3` で一部だけ）。
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

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import RUNS_DIR
from topsem.world.lexicon import WET_SIDE, Lexicon

OUT = ROOT / "docs" / "verification" / "v71_b13_cluster.json"
CORR_LEX = ROOT / "docs" / "verification" / "v70_lexicon_corrected.yaml"
SEED = 20261102
NAMES = ("names", "M1c", "M3n", "M3np")
RUNS = {
    "C1a_gsw22": ("world_v2_b13gsw22", [NAMES, ("gsw_cell", "M1cW", "M3oW", None), ("gsw_set", "M1cG", "M3oG", None),
                                        ("gsw_set_np", "M1cGn", "M3oGn", None)]),
    "C1b_recprior": ("world_v2_b13recprior", [NAMES, ("recur_cell", "M1cEc", "M3oEc", None),
                                              ("recur_set", "M1cE", "M3oE", None)]),
    "C1c_jrcset": ("world_v2_b13jrcset", [NAMES, ("jrc_cell", "M1cJc", "M3oJc", None), ("jrc_set", "M1cJ", "M3oJ", None)]),
    "C2_lexcorr": ("world_v2_b13lexfix", [NAMES]),
    "C3_grouped_es": ("world_v2_b13gres", [NAMES]),
}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v56 = _load("v56_b7_cluster")
wa = v56.wa


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _eval(key: str) -> dict | str:
    tag, pairs = RUNS[key]
    run = RUNS_DIR / tag
    if not v56._complete(run):
        return "missing"
    data = wa.load_split(run, "primary")
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    out = {"run": tag, "n_aoi": len(np.unique(data["aoi"]))} | v56._verdicts(y, preds, data["unit"], pairs)
    rec = _json(run / "predictions" / "primary" / "fold0.json")
    if rec.get("b6"):
        out["record"] = rec["b6"]
    out["lexicon"] = rec.get("lexicon")
    if key == "C3_grouped_es":
        # fold の番号は予測のファイルに無い。load_split と同じ順（名前の順）にファイルの長さで作る
        parts = sorted((run / "predictions" / "primary").glob("fold*.npz"))
        f = np.concatenate([np.full(len(np.load(p)["y"]), k) for k, p in enumerate(parts)])
        sse = sum(float(((y[f == k] - y[f != k].mean()) ** 2).sum()) for k in np.unique(f))
        out["brier_training_constant"] = sse / len(y)
        its = []
        for fp in sorted((run / "predictions" / "primary").glob("fold[0-9].json")):
            its += [r["best"] for r in (_json(fp).get("b13_grouped_es") or []) if "best" in r]
        out["grouped_es_iterations"] = {"n_fits": len(its), "median": float(np.median(its)) if its else None,
                                        "min": int(min(its)) if its else None, "max": int(max(its)) if its else None,
                                        "share_at_cap": float(np.mean(np.array(its) >= 200)) if its else None}
    return out


# %% C2: 直した辞書での集落の近くのオッズ比
def _corrected_flags(gids: np.ndarray, gn: pd.DataFrame) -> pd.Series:
    lex = Lexicon(CORR_LEX)
    u = np.unique(gids[gids >= 0])
    u = u[np.isin(u, gn.index.to_numpy())]
    flag = [float(any(lex.match(str(n), str(c))[k] > 0 for k in WET_SIDE))
            for n, c in zip(gn.loc[u, "asciiname"].to_numpy(), gn.loc[u, "cc"].to_numpy(), strict=True)]
    return pd.Series(flag, index=u)


def c2_near() -> dict:
    v70b = _load("v70b_b13_part2")
    v55, v68 = v70b.v69.v55, v70b.v68
    gn = v70b._gn()
    out: dict = {}
    # 上限ありの EMS（v37 の表: 最寄りの集落と M1 の held-out の予測）
    t = pd.read_parquet(v70b.VER / "v37_table.parquet")
    g = t["P_gid1"].to_numpy(np.int64)
    fl = _corrected_flags(g, gn)
    t = t.assign(lex_wetside_n1=fl.reindex(g).fillna(0.0).to_numpy())
    out["ems_capped"] = v55.odds(t, cls="wetside")
    out["share_cells_flagged_capped"] = float((t["lex_wetside_n1"].to_numpy() > 0).mean())
    # 計画した GFM の推定（v41 の表、クラスタ頑健）
    gt = v70b._planned_gfm()
    gg = gt["P_gid1"].to_numpy(np.int64)
    gt = gt.assign(lex_wetside_n1=_corrected_flags(gg, gn).reindex(gg).fillna(0.0).to_numpy())
    out["gfm_planned"] = v70b._odds_planned(gt)
    # 上限なしの EMS（B-11 の B2 と同じ表と pairs の再抽出）
    u = v70b._uncapped()
    ug = u["P_gid1"].to_numpy(np.int64)
    u = u.assign(lex_wetside_n1=_corrected_flags(ug, gn).reindex(ug).fillna(0.0).to_numpy())
    out["ems_uncapped"] = v68.pairs_or(u, "v71_c2_uncapped")
    out["rule"] = ("B-13 C2: water-related flags of the nearest settlement read with the corrected lexicon "
                   "(v70_lexicon_corrected.yaml); Sect. 5 estimator on the capped table (v37), the planned GFM estimate "
                   "(v41) and the uncapped table with the pairs bootstrap of B-11 B2")
    return out


# %% C4: B-12 の C1（1 km）
def c4_bk() -> dict | str:
    v64 = _load("v64_b9_aggregate")
    v59 = v64.v59
    r = v64.run_gains(RUNS_DIR / "world_v2_b12plas1k", ("BK",))
    if r is None:
        return "missing"
    deltas, gains, done = r
    refs, iv = v59.references(), v64.intervals()
    rng = np.random.default_rng(SEED + 4)
    grid = np.concatenate([[0.0], deltas])
    g = gains["BK"]
    om = float(np.exp(deltas[-1]))
    margin = v59._entry(g, grid, "reach", wa.SESOI * iv["yardstick"]["point"], om, rng)
    bounds = {kind: {ref: v59._entry(g, grid, "bound", iv[kind]["U90"] - d, om, rng)
                     for ref, d in refs.items()} for kind in ("fixed", "refit")}
    ort = margin.get("odds_ratio_text")
    text = (f"the margin corresponded to an injected odds ratio of {ort}" if margin.get("odds_ratio") is not None
            else f"the margin was not reached up to an injected odds ratio of {om:.0f}")
    return {"deltas": deltas.tolist(), "ranges": done, "reps": int(g.shape[1]), "gain_mean": g.mean(axis=1).tolist(),
            "margin_or": margin, "bounds": bounds, "text": text,
            "rule": "B-12 C1: shape BK (cells within 1,000 m of a water-related settlement), as v64 C3"}


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
        res |= {"rule": "amendment B-13 (post hoc), cluster runs and corrected-lexicon estimates", "seed": SEED,
                key: val}
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        print(f"{key} ({time.time() - t0:.0f} s) {json.dumps(val, default=float)[:700]}", flush=True)

    for key in RUNS:
        run(key, lambda key=key: _eval(key))
    run("C2_near", c2_near)
    run("C4_bk", c4_bk)


if __name__ == "__main__":
    main()
