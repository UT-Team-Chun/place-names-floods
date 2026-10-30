# %% v78（事後・amendment B-16 の C1・C3・C4）: クラスタの実行の集計
"""確証判定は変えない。B-16 の追記 1 で計算の前に固定した量。

- C1: 強い基準 M1p（FABDEM・河道からの高さ・上流の面積・JRC・GSW・WorldCover）の上の名の増分。主分割と活性化を
  分けない分割。R は M1p の物差し BS(M0) − BS(M1p)（全体の物差し BS(M0) − BS(M1) の R も並べる）。判定は v56 の
  _verdicts（場所グループ 4,000 回、モデル固定、プロトコルの規則）
- C3: 集落の地点の HAND の AOI 内の三分位の中で名を入れ替えた置換（200 回、地域構成と N-learn を作り直す）の、入れ子と
  中心化した比較の R の分布と片側の p
- C4: 言語モデルの切れ目の分割で、プロトコルの placebo の置換を 20 回引き直した判定

出力: stdout と docs/verification/v78_b16_cluster.json（部分ごとに足す。`--only C1` で一部だけ）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v78_b16_cluster.json"
STRONG = {"primary": "world_v2_b16strong", "activation_merged": "world_v2_b16strongact"}
# 追記 4: 中心化した積み重ね（M3pd・M3pdp）の回し直し。M0・M1・M1p・M1ps が C1 の実行と一致するときだけ使う
STRONG_D = {"primary": "world_v2_b16strongd", "activation_merged": "world_v2_b16strongactd"}
STRONG_D_CHECK = ("M0", "M1", "M1p", "M1ps")
STRAT = "world_v2_b16strat"
PLACEBO_SEEDS = tuple(20261105 + k for k in range(1, 21))
LLMCUT_PAIRS = [("primary", "M1c", "M3n", "M3np"), ("total", "M1s", "M3s", "M3ps"), ("flexible", "M1", "M3", "M3p"),
                ("H4", "M1rc", "M3rn", "M3rnp"), ("H4_total", "M1s", "M3rs", "M3rps"),
                ("H4_flexible", "M1", "M3r", "M3rp")]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v71 = _load("v71_b13_cluster")
v56, wa = v71.v56, v71.wa


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _complete(run: Path, split: str, n: int) -> bool:
    return len(list((run / "predictions" / split).glob("fold[0-9]*.npz"))) >= n


# %% C1: 強い基準
def c1_strong() -> dict:
    out: dict = {}
    for key, tag in STRONG.items():
        run = RUNS_DIR / tag
        if not _complete(run, "primary", 10):
            out[key] = "missing"
            continue
        data = wa.load_split(run, "primary")
        y = data["y"].astype(np.float64)
        preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
        check = _strong_d_check(key, data)
        if check.get("status") == "passed":
            dd = wa.load_split(RUNS_DIR / STRONG_D[key], "primary")
            preds |= {a: dd[f"p_{a}"] for a in ("M3pd", "M3pdp")}
        pairs = [("strong_nested", "M1pc", "M3pn", "M3pnp"), ("strong_centred", "M1ps", "M3pd", "M3pdp")]
        own = v56._verdicts(y, preds | {"M1": preds["M1p"]}, data["unit"], pairs)
        main = v56._verdicts(y, preds, data["unit"], pairs)
        rec = _json(run / "predictions" / "primary" / "fold0.json")
        out[key] = {"run": tag, "n_aoi": len(np.unique(data["aoi"])), "strong_yardstick": own,
                    "main_yardstick": {k: main[k] for k in ("yardstick", "ci95_yardstick", "strong_nested",
                                                            "strong_centred") if k in main},
                    "brier_M1_M1p": {"M1": float(np.mean((preds["M1"] - y) ** 2)) if "M1" in preds else None,
                                     "M1p": float(np.mean((preds["M1p"] - y) ** 2))},
                    "record": rec.get("b6"), "group_by_activation": rec.get("group_by_activation"),
                    "centred_rerun_check": check}
    out["rule"] = ("B-16 C1 (addendum 1, item 12; addendum 4): R with the strong baseline's yardstick BS(M0) - BS(M1p) (M1 "
                   "replaced by M1p in v56._verdicts) and with the main yardstick; 4,000 draws of location groups, models "
                   "fixed. The centred stacks M3pd and M3pdp come from the reruns of addendum 4, used only when M0, M1, M1p "
                   "and M1ps of the rerun equal those of the C1 run to 1e-9")
    return out


def _strong_d_check(key: str, data: dict) -> dict:
    """追記 4 の 2: 回し直した実行の M0・M1・M1p・M1ps が C1 の実行と一致するか（セルも同じ順か）。"""
    run = RUNS_DIR / STRONG_D[key]
    if not _complete(run, "primary", 10):
        return {"status": "missing", "run": STRONG_D[key]}
    dd = wa.load_split(run, "primary")
    same = len(dd["y"]) == len(data["y"]) and bool(np.array_equal(dd["y"], data["y"])) and bool(
        np.array_equal(dd["unit"], data["unit"]))
    diffs = ({a: float(np.max(np.abs(dd[f"p_{a}"].astype(np.float64) - data[f"p_{a}"].astype(np.float64))))
              for a in STRONG_D_CHECK} if same else {})
    ok = same and len(diffs) == len(STRONG_D_CHECK) and max(diffs.values()) <= FROZEN_TOL
    return {"status": "passed" if ok else "failed", "run": STRONG_D[key], "same_cells": same, "max_abs_diff": diffs,
            "tolerance": FROZEN_TOL}


# %% C3: 層別した置換
def c3_strata() -> dict:
    v77 = _load("v77_b16_local")
    v77.PERM_RUN = STRAT
    p = v77._perm_table()
    if p is None:
        return {"note": f"{STRAT} incomplete"}
    out: dict = {"n_perms": int(p["n_perms"]), "n_units": len(p["units"])}
    for key, real_e, perm_e in (("nested", p["real_M1c"] - p["real_M3n"], p["ref"] - p["names"]),
                                ("centred", p["real_M1s"] - p["real_M3d"], p["real_M1s"][None, :] - p["cent"])):
        d_hat, d_j = float(real_e.sum()), perm_e.sum(axis=1)
        k = int(np.sum(d_j >= d_hat))  # 同じ物差しなので R_j ≥ R̂ ⇔ Δ_j ≥ Δ̂
        out[key] = {"delta_sum_real": d_hat, "share_perm_ge_real": k / len(d_j),
                    "p_one_sided": (1 + k) / (1 + len(d_j)),
                    "perm_delta_sum_quantiles": np.quantile(d_j, [0.025, 0.5, 0.975]).tolist()}
    out["rule"] = ("B-16 C3 (addendum 1, item 13): names permuted within AOI terciles of the settlement-site HAND; "
                   "features and composition rebuilt; p = (1 + #{R_j >= R}) / (1 + 200)")
    return out


# %% C4: placebo の引き直し（追記 3 で直した集計）
LEX_KEYS = ("primary", "total", "flexible")  # 辞書を使う比較（凍結した読みで引き直す）
H4_KEYS = ("H4", "H4_total", "H4_flexible")  # N-read だけを使う比較（二つの読みで同じ）
FROZEN_ARMS = ("M0", "M1", "M3", "M3p", "M1c", "M3n", "M3np", "M1s", "M3s", "M3ps")
FROZEN_TOL = 1e-9


def _redraw_verdicts(tags: list[str], pairs: list) -> tuple[list, dict]:
    """引き直しの実行ごとの判定と、比較ごとの判定の並び。"""
    done, labels = [], {k: [] for k, *_ in pairs}
    for tag in tags:
        run = RUNS_DIR / tag
        if not _complete(run, "llmcut", 1):
            continue
        d = wa.load_split(run, "llmcut")
        v = v56._verdicts(d["y"].astype(np.float64), {k[2:]: x for k, x in d.items() if k.startswith("p_")},
                          d["unit"], pairs)
        rec = _json(run / "predictions" / "llmcut" / "fold0.json")
        done.append({"run": tag, "seed": rec.get("placebo_seed"),
                     "frozen_norway": bool((rec.get("b6") or {}).get("frozen_norway")),
                     "labels": {k: (v.get(k) or {}).get("decision") for k in labels}})
        for k, lst in labels.items():
            lst.append((v.get(k) or {}).get("decision"))
    return done, labels


def _changes(labels: dict, ref: dict) -> dict:
    """基準の判定（ref の比較ごとの decision）から変わった回数と、判定の内訳。"""
    return {k: {"n": len(v), "reference": (ref.get(k) or {}).get("decision"),
                "n_differs_from_protocol": sum(x != (ref.get(k) or {}).get("decision") for x in v),
                "labels": {lab: v.count(lab) for lab in sorted({str(x) for x in v})}}
            for k, v in labels.items()}


def _frozen_check() -> dict:
    """追記 3 の 2 (a): 凍結した読みの旗とプロトコルの種の実行が、確証の実行の予測と一致するか。"""
    run = RUNS_DIR / "world_v2_b16placf00"
    if not _complete(run, "llmcut", 1):
        return {"status": "missing"}
    d, conf = wa.load_split(run, "llmcut"), wa.load_split(RUNS_DIR / "world_v2", "llmcut")
    same = len(d["y"]) == len(conf["y"]) and bool(np.array_equal(d["y"], conf["y"]))
    diffs = ({a: float(np.max(np.abs(d[f"p_{a}"].astype(np.float64) - conf[f"p_{a}"].astype(np.float64))))
              for a in FROZEN_ARMS} if same else {})
    rec = _json(run / "predictions" / "llmcut" / "fold0.json")
    ok = same and len(diffs) == len(FROZEN_ARMS) and max(diffs.values()) <= FROZEN_TOL
    return {"status": "passed" if ok else "failed", "same_cells": same, "max_abs_diff": diffs,
            "frozen_norway": bool((rec.get("b6") or {}).get("frozen_norway")), "seed": rec.get("placebo_seed"),
            "tolerance": FROZEN_TOL}


def c4_placebo() -> dict:
    conf = wa.load_split(RUNS_DIR / "world_v2", "llmcut")
    lexfix = wa.load_split(RUNS_DIR / "world_v2_lexfix", "llmcut")
    out: dict = {
        "protocol": v56._verdicts(conf["y"].astype(np.float64), {k[2:]: v for k, v in conf.items() if k.startswith("p_")},
                                  conf["unit"], LLMCUT_PAIRS),
        "corrected_code_reference": v56._verdicts(lexfix["y"].astype(np.float64),
                                                  {k[2:]: v for k, v in lexfix.items() if k.startswith("p_")},
                                                  lexfix["unit"], LLMCUT_PAIRS)}
    # 既にある 20 回（直したコード、ノルウェーを読む）
    done, labels = _redraw_verdicts([f"world_v2_b16plac{s - 20261105:02d}" for s in PLACEBO_SEEDS], LLMCUT_PAIRS)
    out["redraws"] = done
    out["changes_vs_protocol_mixed"] = _changes(labels, out["protocol"])  # v78 の元の数（読みの直しが混ざる）
    out["changes_corrected_code"] = _changes(labels, out["corrected_code_reference"])
    # 凍結した読みの引き直し（追記 3 の 2）。(a) が通るまで (b) の判定を計算しない
    out["frozen_check"] = _frozen_check()
    lex_pairs = [x for x in LLMCUT_PAIRS if x[0] in LEX_KEYS]
    if out["frozen_check"]["status"] == "passed":
        fdone, flabels = _redraw_verdicts([f"world_v2_b16placf{s - 20261105:02d}" for s in PLACEBO_SEEDS], lex_pairs)
        out["frozen_redraws"] = fdone
        out["changes_frozen_reading"] = _changes(flabels, out["protocol"])
    # 本文の数（追記 3 の 3）: 辞書を使う比較は凍結した読みの引き直し、H4 は既にある引き直し（プロトコルと同じ基準）
    paper = {k: out["changes_vs_protocol_mixed"][k] | {"source": "existing redraws (N-read only)"} for k in H4_KEYS}
    for k in LEX_KEYS:
        if k in out.get("changes_frozen_reading", {}):
            paper[k] = out["changes_frozen_reading"][k] | {"source": "redraws with the frozen reading"}
    out["changes"] = paper
    out["rule"] = ("B-16 C4 (addendum 1, item 14; corrected by addendum 3): the protocol's placebo redrawn with seeds "
                   "20261105 + k (k = 1 to 20) on the language-model cut-off split; verdicts from the protocol's "
                   "aggregation (v56._verdicts = world_aggregate). The existing redraws ran on the corrected code (Norway "
                   "read), so their non-placebo arms equal world_v2_lexfix; the comparisons that use the lexicon are "
                   "redrawn with --frozen-norway (check run world_v2_b16placf00 must reproduce world_v2 to 1e-9)")
    return out


PARTS = {"C1": c1_strong, "C3": c3_strata, "C4": c4_placebo}


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
        print(f"{k} ({r['seconds']} s) {json.dumps(r, ensure_ascii=False, default=float)[:2000]}", flush=True)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
