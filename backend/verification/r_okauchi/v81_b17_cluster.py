# %% v81（事後・amendment B-17 の C1〜C6 と B3）: クラスタの実行の集計
"""確証判定は変えない。B-17 で計算の前に固定した量。

- C1・B3: 言語モデルの切れ目の分割で、プロトコルの placebo を 200 回引いた分布（辞書を使う三つの比較は凍結した読み、
  N-read だけの三つは二つの読みで同じ）。実の R − mean(R_p)、p = (1 + #{R_p ≥ R̂})/(1 + N)、判定の頻度
- C2・C3・B3: 置換（地域構成を作り直す、200 回）の分布に対する実の名の位置。主分割（直したコードの B-7、凍結した読み）、
  活性化を分けない分割（二つの読み）、GFM のラベル
- C4: 裸地で事象の前の情報だけの基準 M1b の上の名
- C5: ラベルの別定義で二段とも当て直した実行
- C6: 井戸と -ab を除いた辞書

出力: stdout と docs/verification/v81_b17_cluster.json（部分ごとに足す。`--only C1` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v81_b17_cluster.json"
SEED = 20261106


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


v78 = _load("v78_b16_cluster")
wa, v56, v77 = v78.wa, v78.v56, v78._load("v77_b16_local")
LEX_KEYS, H4_KEYS = v78.LEX_KEYS, v78.H4_KEYS


# %% C1・B3: placebo の 200 回
def _redraw(tag: str, pairs: list) -> dict | None:
    run = RUNS_DIR / tag
    if not v78._complete(run, "llmcut", 1):
        return None
    d = wa.load_split(run, "llmcut")
    return v56._verdicts(d["y"].astype(np.float64), {k[2:]: x for k, x in d.items() if k.startswith("p_")},
                         d["unit"], pairs)


def c1_placebo() -> dict:
    lex_pairs = [x for x in v78.LLMCUT_PAIRS if x[0] in LEX_KEYS]
    h4_pairs = [x for x in v78.LLMCUT_PAIRS if x[0] in H4_KEYS]
    conf = wa.load_split(RUNS_DIR / "world_v2", "llmcut")
    proto = v56._verdicts(conf["y"].astype(np.float64), {k[2:]: v for k, v in conf.items() if k.startswith("p_")},
                          conf["unit"], v78.LLMCUT_PAIRS)
    draws: dict[str, list] = {k: [] for k, *_ in v78.LLMCUT_PAIRS}
    used: dict[str, list] = {k: [] for k, *_ in v78.LLMCUT_PAIRS}
    # 辞書を使う三つ: B-16 の凍結した読みの 20 回（placf01〜20）と B-17 の 180 回（凍結した読み）
    tags_lex = [f"world_v2_b16placf{k:02d}" for k in range(1, 21)] + [f"world_v2_b17plac{k:03d}" for k in range(21, 201)]
    tags_h4 = [f"world_v2_b16plac{k:02d}" for k in range(1, 21)] + [f"world_v2_b17plac{k:03d}" for k in range(21, 201)]
    for tags, pairs in ((tags_lex, lex_pairs), (tags_h4, h4_pairs)):
        for tag in tags:
            v = _redraw(tag, pairs)
            if v is None:
                continue
            for name, *_ in pairs:
                e = v.get(name) or {}
                pl = e.get("placebo") or {}
                rp = pl.get("R_placebo")
                if rp is None:
                    continue
                draws[name].append({"R_p": float(rp), "decision": e.get("decision"), "R": float(e["R"])})
                used[name].append(tag)
    out: dict = {"protocol": {k: (proto.get(k) or {}).get("decision") for k in draws}}
    for name, rows in draws.items():
        if not rows:
            out[name] = {"n": 0}
            continue
        # 実の名の R はどの引き直しでも同じ（placebo の腕だけが変わる）。確証の実行に無い比較は引き直しの中央値
        r_real = float((proto.get(name) or {}).get("R", np.median([x["R"] for x in rows])))
        rps = np.array([x["R_p"] for x in rows])
        k = int((rps >= r_real).sum())
        out[name] = {"n": len(rows), "R": r_real, "mean_R_p": float(rps.mean()),
                     "R_minus_mean_R_p": r_real - float(rps.mean()),
                     "R_p_q025_q975": np.quantile(rps, [0.025, 0.975]).tolist(),
                     "p_one_sided": (1 + k) / (1 + len(rps)),
                     "verdicts": dict(Counter(str(x["decision"]) for x in rows)),
                     "max_abs_diff_R_vs_protocol": float(max(abs(x["R"] - r_real) for x in rows))}
    # 辞書を使う三つの、直したコードの 20 回（B-16 plac01〜20）の判定の内訳（記述）
    _, corr = v78._redraw_verdicts([f"world_v2_b16plac{k:02d}" for k in range(1, 21)], lex_pairs)
    out["lexicon_corrected_20"] = {k: dict(Counter(str(x) for x in v)) for k, v in corr.items()}
    out["rule"] = ("B-17 C1 and B3: language-model cut-off split; the comparisons that use the lexicon from the frozen-reading "
                   "redraws (B-16 placf01-20 and B-17 plac021-200), those that use only N-read from B-16 plac01-20 and "
                   "B-17 plac021-200 (the same under both readings); verdicts from v56._verdicts; p = (1 + #{R_p >= R})/(1 + N)")
    return out


# %% C2・C3・B3: 置換の分布
N_PERM_FILES = 40  # 10 fold × 50 回ずつ 4 本


def _perm(run: str) -> dict | None:
    if len(list((RUNS_DIR / run / "predictions" / "primary").glob("fold*_permsrc_*.npz"))) < N_PERM_FILES:
        return None
    v77.PERM_RUN = run
    return v77._perm_table()


def _yard_units(yard_run: str, check_run: str, units: np.ndarray, p: dict) -> dict | None:
    """物差し BS(M0) − BS(M1) の場所グループの和と、実の腕の一致の点検。

    M0・M1 は名を使わないので、同じ分割の実行ならどれでも同じ（world_v2・lexfix・b13lexfix、actgrp・b7actfix・
    b16strongact で差 0 を確かめた）。点検には読みとコードが置換の実行と同じ実行の M1c・M3n・M1s を使う。
    """
    out: dict = {}
    for role, run in (("yard", yard_run), ("check", check_run)):
        path = RUNS_DIR / run
        if not v78._complete(path, "primary", 10):
            return None
        d = wa.load_split(path, "primary")
        arms = tuple(a for a in ("M0", "M1", "M1c", "M3n", "M1s", "M3d") if f"p_{a}" in d)
        _n, s, u = v77.v69._group_sums(d, arms)
        pos = pd.Series(np.arange(len(u)), index=u.astype(str)).reindex(units).to_numpy()
        if np.isnan(pos).any() or len(u) != len(units):
            return None
        pos = pos.astype(int)
        if role == "yard":
            out["dg"] = s["M0"][pos] - s["M1"][pos]
        else:
            out["check_max_rel_diff"] = {a: float(np.max(np.abs(s[a][pos] - p[f"real_{a}"])
                                                       / np.maximum(p[f"real_{a}"], 1e-12)))
                                         for a in ("M1c", "M3n", "M1s", "M3d") if a in s}
    return out


# 置換の実行: (置換, 物差しの実行, 実の腕を点検する実行)
PERMS = {"primary_corrected": ("world_v2_b7perm", "world_v2", "world_v2_lexfix"),
         "primary_frozen": ("world_v2_b17perm_frozen", "world_v2", "world_v2"),
         "activation_whole_frozen": ("world_v2_b17permact_frozen", "world_v2_actgrp", "world_v2_actgrp"),
         "activation_whole_corrected": ("world_v2_b17permact", "world_v2_actgrp", "world_v2_b7actfix"),
         "gfm_corrected": ("world_v2_b17permgfm", "world_v2_gfm", "world_v2_gfm")}


def c2_perms() -> dict:
    out: dict = {}
    for key, (run, yard_run, check_run) in PERMS.items():
        p = _perm(run)
        if p is None:
            out[key] = {"note": f"{run} not complete locally"}
            continue
        units = p["units"]
        yard = _yard_units(yard_run, check_run, units, p)
        dg = yard["dg"] if yard else None
        big_n = float(p["n"].sum())
        res = {"run": run, "n_perms": int(p["n_perms"]), "n_units": len(units), "yardstick_run": yard_run,
               "check_run": check_run, "check_max_rel_diff": yard["check_max_rel_diff"] if yard else None,
               "yardstick": float(dg.sum() / big_n) if dg is not None else None}
        for comp, real_e, perm_e in (("nested", p["real_M1c"] - p["real_M3n"], p["ref"] - p["names"]),
                                     ("centred", p["real_M1s"] - p["real_M3d"], p["real_M1s"][None, :] - p["cent"])):
            d_hat = float(real_e.sum() / big_n)
            d_j = perm_e.sum(axis=1) / big_n
            k1 = int((d_j >= d_hat).sum())
            k2 = int((np.abs(d_j - d_j.mean()) >= abs(d_hat - d_j.mean())).sum())
            r = {"delta_real": d_hat, "delta_perm_mean": float(d_j.mean()),
                 "delta_perm_q025_q975": np.quantile(d_j, [0.025, 0.975]).tolist(),
                 "p_one_sided": (1 + k1) / (1 + len(d_j)), "p_two_sided": (1 + k2) / (1 + len(d_j))}
            if dg is not None:
                yard = float(dg.sum() / big_n)
                r |= {"R_real": d_hat / yard, "R_minus_mean_R_perm": (d_hat - float(d_j.mean())) / yard,
                      "R_perm_q025_q975": (np.quantile(d_j, [0.025, 0.975]) / yard).tolist()}
            res[comp] = r
        out[key] = res
    out["rule"] = ("B-17 C2, C3 and B3: permutations with composition rebuilt and the second stage fitted (200, seeds 1-200); "
                   "Delta_j summed over the run's location groups; one-sided p = (1 + #{Delta_j >= Delta})/(1 + N), "
                   "two-sided about the permutation mean; R with the yardstick of a run on the same cells")
    return out


# %% C4・C5・C6: 基準・ラベル・辞書の実行
def _run_verdicts(tag: str, pairs: list, *, yard: str | None = None) -> dict:
    run = RUNS_DIR / tag
    if not v78._complete(run, "primary", 10):
        return {"note": f"{tag} not complete locally"}
    d = wa.load_split(run, "primary")
    y = d["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in d.items() if k.startswith("p_")}
    own = v56._verdicts(y, preds | ({"M1": preds[yard]} if yard else {}), d["unit"], pairs)
    rec = _json(run / "predictions" / "primary" / "fold0.json")
    out = {"run": tag, "n_cells": len(y), "flooded_share": float(y.mean()), "own_yardstick": own,
           "record": rec.get("b6"), "label_variant": (rec.get("sensitivity") or {}).get("label_variant")}
    if yard:
        out["main_yardstick"] = v56._verdicts(y, preds, d["unit"], pairs)
        out["brier_M1_M1b"] = {"M1": float(np.mean((preds["M1"] - y) ** 2)), yard: float(np.mean((preds[yard] - y) ** 2))}
    return out


BASE_PAIRS = [("nested", "M1c", "M3n", "M3np"), ("centred", "M1s", "M3d", "M3dp")]
BARE_PAIRS = [("bare_nested", "M1bc", "M3bn", "M3bnp"), ("bare_centred", "M1bs", "M3bd", "M3bdp")]


def c4_bare() -> dict:
    return {key: _run_verdicts(tag, BARE_PAIRS, yard="M1b")
            for key, tag in (("primary", "world_v2_b17bare"), ("activation_whole", "world_v2_b17bareact"))} | {
        "rule": "B-17 C4: R with the bare-earth baseline's gain BS(M0) - BS(M1b) (M1 replaced by M1b in v56._verdicts) and "
                "with the main one; 4,000 draws of location groups, models fixed"}


def c5_labels() -> dict:
    return {v: _run_verdicts(f"world_v2_b17{v}", BASE_PAIRS) for v in ("t010", "t025", "landobs", "gswmask")} | {
        "rule": "B-17 C5: both stages refitted on each label definition (frozen reading); R with the definition's own "
                "yardstick; verdicts from v56._verdicts"}


def c6_nochah() -> dict:
    return {"nochah": _run_verdicts("world_v2_b17nochah", BASE_PAIRS),
            "rule": "B-17 C6: lexicon without chah, chahi, ab and aab in block fa (frozen reading otherwise)"}


PARTS = {"C1": c1_placebo, "C2": c2_perms, "C4": c4_bare, "C5": c5_labels, "C6": c6_nochah}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    want = [p for p in (a.only.split(",") if a.only else PARTS) if p]
    for k in want:
        t0 = time.time()
        r = PARTS[k]()
        r["seconds"] = round(time.time() - t0)
        res = _json(OUT) if OUT.is_file() else {}
        res[k] = r
        OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
        print(f"{k} ({r['seconds']} s) {json.dumps(r, ensure_ascii=False, default=float)[:2500]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
