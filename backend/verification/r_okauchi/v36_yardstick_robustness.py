# %% v36（事後・amendment B-2）: 物差しと推論の頑健性 — 名前の増分の固定マージン TOST を物差し 5 通り × bootstrap 2 通りで
"""確証判定（world_aggregate の label）は変えない。内部審査の指摘 1・2・5・6・12 に応じた記述。

物差し Y（余白 = 0.05 × Y、名前の増分 ΔBS の 90% 区間が (−余白, +余白) に入れば同等）:
- Y0: 登録どおり ΔBŜ(M0 → M1)（柔軟モデル同士の点推定）
- Y1: ΔBS(M0 → M1) の 95% 区間の下限（bootstrap と同じ種類の区間）
- Y2: ΔBŜ(M0s → M1s)（両方を二段目で再較正した後。M0s は world_v2_recal）
- Y3: ΔBŜ(定数 → M1s)。定数は fold ごとの訓練側の浸水率（10 fold の分割）、評価側の浸水率
  （1 fold の時間外し・llmcut。定数に有利で物差しを小さくする側 = 保守的）
- Y4: ΔBŜ(M0c → M1c)（二段目同士）
bootstrap は登録どおりのクラスタ（主・時間外し・llmcut は場所グループ、国外しは国。B = 4000、
world_aggregate と同じ重み）と、activation × 国（AOI の多数派の国）の二方向（B = 4000）。

ほかに: 地域構成の増分を再較正だけの基準から（M0s → M0c・T0s → T0c・M1s → M1c）、比の bootstrap 区間
（全ての R）、H5 の層（市街地・辞書の被覆・大陸）、H3・H4 の Holm 補正、決着しなかった比較の
「要る場所グループ数」の目安、fold の分け方・シード・グループの作り方を変えた実行の間の R の幅。

入力: data/runs/world_v2*/predictions/<split>/fold*.npz、metrics_<split>.json
出力: stdout の表と docs/verification/v36_yardstick_robustness.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.dataset import continent_map
from topsem.world.lexicon import Lexicon

OUT = ROOT / "docs" / "verification" / "v36_yardstick_robustness.json"
LANG: dict[str, str] = {}  # main で countryInfo から埋める
SPLITS = ("primary", "country", "temporal", "llmcut")
SINGLE_FOLD = ("temporal", "llmcut")
SESOI = wa.SESOI
B = wa.B
COMPARISONS = (
    ("primary", "M1c", "M3n"),
    ("H4", "M1rc", "M3rn"),
    ("total", "M1s", "M3s"),
    ("H4_total", "M1s", "M3rs"),
    ("flexible", "M1", "M3"),
    ("H4_flexible", "M1", "M3r"),
    # 判定の無い比較（記述）: 固定マージンの区間を判定と同じ形でそろえる（world_aggregate は比の bootstrap）
    ("placebo", "M1c", "M3np"),
    ("channel_lexicon", "M1c", "M3nx"),
    ("channel_learned", "M1c", "M3nl"),
    ("H1", "M0c", "M0n"),
    ("H1_flexible", "M0", "M0N"),
    ("H4_placebo", "M1rc", "M3rnp"),
    # B-2 追記 2: AOI 内で中心化した名前の特徴（地域構成の係数に依らない場所内の増分）
    ("demeaned", "M1s", "M3d"),
    ("demeaned_placebo", "M1s", "M3dp"),
)
COMPOSITION = (
    ("M1s", "M1rc"),  # N-read の地域構成（llmcut、B-2 追記 2）
    ("M0", "M0c"),
    ("M0s", "M0c"),
    ("T0", "T0c"),
    ("T0s", "T0c"),
    ("M1s", "M1c"),
)
SPREAD_RUNS = (
    "world_v2",
    "world_v2_fs1",
    "world_v2_fs2",
    "world_v2_actgrp",
    "world_v2_s5_k1",
    "world_v2_s6_buffer",
    "world_v2_s6_buffer_b1",
)


# %% 読み込み
def load(run: str, split: str) -> dict[str, np.ndarray] | None:
    files = sorted(
        (RUNS_DIR / run / "predictions" / split).glob("fold*.npz"),
        key=lambda p: int(p.stem[4:]),
    )
    if not files:
        return None
    acc: dict[str, list[np.ndarray]] = {}
    for f in files:
        with np.load(f, allow_pickle=False) as z:
            for k in z.files:
                acc.setdefault(k, []).append(z[k])
            acc.setdefault("fold", []).append(
                np.full(len(z["y"]), int(f.stem[4:]), dtype=np.int16)
            )
    return {k: np.concatenate(v) for k, v in acc.items()}


def add_recal(
    d: dict[str, np.ndarray], split: str, recal_run: str = "world_v2_recal"
) -> dict:
    """recal_run から物差しの腕（M0s・T0s、無ければ M0c・M0・M1 も）を足す。行の並びと、共有する腕の予測が
    同じかを確かめる。中心化した名前の特徴の実行（world_v2_dm*）は M0s・M0c を持たないので、ここで補う。"""
    if "p_M0s" in d:  # 自分の run に再較正の基準がある（activation で束ねた分割など）
        return {"available": True, "own": True}
    r = load(recal_run, split)
    if r is None:
        return {"available": False}
    if not (np.array_equal(r["y"], d["y"]) and np.array_equal(r["aoi"], d["aoi"])):
        return {"available": False, "reason": "rows differ from world_v2"}
    shared = sorted(k for k in r if k.startswith("p_") and k in d)
    diff = {
        k[2:]: float(np.max(np.abs(r[k].astype(np.float64) - d[k].astype(np.float64))))
        for k in shared
    }
    added = []
    for k in ("p_M0s", "p_T0s", "p_M0c", "p_M0", "p_M1"):
        if k in r and k not in d:
            d[k] = r[k]
            added.append(k[2:])
    return {"available": True, "max_abs_diff_vs_world_v2": diff, "arms_added": added}


def constant(d: dict[str, np.ndarray], split: str) -> np.ndarray:
    y = d["y"].astype(np.float64)
    if split in SINGLE_FOLD:
        return np.full(len(y), y.mean())
    p = np.empty(len(y))
    for k in np.unique(d["fold"]):
        m = d["fold"] == k
        p[m] = (y.sum() - y[m].sum()) / (len(y) - m.sum())
    return p


def aoi_country(d: dict[str, np.ndarray]) -> np.ndarray:
    """セルごとに、その AOI の多数派の国（最寄り集落の国コードの最頻値）。"""
    out = np.empty(len(d["aoi"]), dtype="<U2")
    for a in np.unique(d["aoi"]):
        m = d["aoi"] == a
        c = Counter(x for x in d["cc"][m].tolist() if x)
        out[m] = c.most_common(1)[0][0] if c else "??"
    return out


# %% bootstrap
def sums(
    y: np.ndarray, preds: dict[str, np.ndarray], key: np.ndarray
) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    uniq, codes = np.unique(key, return_inverse=True)
    arms = sorted(preds)
    s = np.column_stack(
        [
            np.bincount(
                codes,
                weights=(preds[a].astype(np.float64) - y) ** 2,
                minlength=len(uniq),
            )
            for a in arms
        ]
    )
    return s, np.bincount(codes, minlength=len(uniq)).astype(np.float64), arms, uniq


def two_way(
    y: np.ndarray,
    preds: dict[str, np.ndarray],
    act: np.ndarray,
    ctry: np.ndarray,
    b: int = B,
    seed: int = wa.SEED + 2,
) -> tuple[np.ndarray, dict]:
    """activation × 国の二方向 bootstrap（両軸の多項重みの積）の pooled BS [b, arms]。"""
    s, n, _, uniq = sums(y, preds, np.char.add(np.char.add(act, "|"), ctry))
    ua, ia = np.unique([u.split("|")[0] for u in uniq], return_inverse=True)
    uc, ic = np.unique([u.split("|")[1] for u in uniq], return_inverse=True)
    rng = np.random.default_rng(seed)
    out = np.full((b, s.shape[1]), np.nan)
    for r in range(b):
        w = (
            rng.multinomial(len(ua), np.full(len(ua), 1 / len(ua)))[ia]
            * rng.multinomial(len(uc), np.full(len(uc), 1 / len(uc)))[ic]
        ).astype(np.float64)
        if w @ n > 0:
            out[r] = (w @ s) / (w @ n)
    return out, {"n_activations": len(ua), "n_countries": len(uc)}


def q(x: np.ndarray, lo: float, hi: float) -> list[float]:
    x = x[np.isfinite(x)]
    return [float(np.quantile(x, lo)), float(np.quantile(x, hi))]


# %% 一つの分割
def analyse(
    split: str,
    lex: Lexicon,
    cont: dict[str, str],
    run: str = "world_v2",
    recal_run: str = "world_v2_recal",
) -> dict | None:
    d = load(run, split)
    if d is None:
        return None
    recal = add_recal(d, split, recal_run)
    y = d["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in d.items() if k.startswith("p_")}
    preds["const"] = constant(d, split)
    s, n, arms, _ = sums(y, preds, d["unit"])
    idx = {a: j for j, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    boots = {"cluster": wa.bootstrap_bs(s, n)}
    act = np.array([a.split("_")[0] for a in d["aoi"]])
    boots["two_way"], tw_info = two_way(y, preds, act, aoi_country(d))

    def delta(boot: np.ndarray, a: str, b: str) -> tuple[float, np.ndarray]:
        return float(hat[idx[a]] - hat[idx[b]]), boot[:, idx[a]] - boot[:, idx[b]]

    res: dict = {
        "n_cells": len(y),
        "n_units": len(n),
        "recal": recal,
        "two_way_axes": tw_info,
        "brier": {a: float(hat[idx[a]]) for a in arms},
    }
    for kind, boot in boots.items():
        y0, b01 = delta(boot, "M0", "M1")
        ys = {
            "Y0": y0,
            "Y1": q(b01, 0.025, 0.975)[0],
            "Y3": delta(boot, "const", "M1s")[0],
            "Y4": delta(boot, "M0c", "M1c")[0],
        }
        if "M0s" in idx:
            ys["Y2"] = delta(boot, "M0s", "M1s")[0]
        out = {"yardsticks": dict(sorted(ys.items()))}
        boots_by_key: dict[str, np.ndarray] = {}
        for key, a, b in COMPARISONS:
            if a not in idx or b not in idx:
                continue
            dn, bd = delta(boot, a, b)
            boots_by_key[key] = bd
            c90 = q(bd, 0.05, 0.95)
            row = {
                "delta_bs": dn,
                "ci90_delta": c90,
                "ci95_delta": q(bd, 0.025, 0.975),
                "smallest_yardstick_for_equivalence": max(abs(c90[0]), abs(c90[1]))
                / SESOI,
                "ci95_R_ratio_bootstrap": q(bd / b01, 0.025, 0.975),
            }
            for yk, yv in ys.items():
                ok = yv > 0
                row[yk] = {
                    "margin": SESOI * yv if ok else None,
                    "R": dn / yv if ok else None,
                    "ci90_R": [c / yv for c in c90] if ok else None,
                    "equivalent": bool(
                        ok and -SESOI * yv < c90[0] and c90[1] < SESOI * yv
                    ),
                }
            if kind == "cluster":
                m = SESOI * y0
                row["p_TOST_Y0"] = float(max(np.mean(bd >= m), np.mean(bd <= -m)))
                row["p_superiority"] = float(np.mean(bd <= 0))
            out[key] = row
        # B-2 追記 3: 凍結した判定表（§6、world_aggregate.decide）の label を、登録の物差し Y0 で。placebo 条件は
        # real − permuted の差の 95% 区間、分母の条件は ΔBS(M0 → M1) の 95% 区間
        den95 = q(b01, 0.025, 0.975)
        for key, pk in (
            ("primary", "placebo"),
            ("demeaned", "demeaned_placebo"),
            ("H4", "H4_placebo"),
        ):
            if key in boots_by_key and pk in boots_by_key and y0 > 0:
                bd, bp = boots_by_key[key], boots_by_key[pk]
                diff95 = q((bd - bp) / y0, 0.025, 0.975)
                dec = wa.decide(
                    out[key]["delta_bs"] / y0,
                    [c / y0 for c in out[key]["ci95_delta"]],
                    [c / y0 for c in out[key]["ci90_delta"]],
                    placebo_ci95=diff95,
                    den_ci95=den95,
                )
                out[key]["decision_Y0"] = dec["label"]
                out[key]["R_minus_Rp_Y0"] = (
                    out[key]["delta_bs"] - out[pk]["delta_bs"]
                ) / y0
                out[key]["ci95_R_minus_Rp_Y0"] = diff95
        comp = {}
        for a, b in COMPOSITION:
            if a in idx and b in idx:
                dc, bc = delta(boot, a, b)
                comp[f"{a}->{b}"] = {
                    "delta_bs": dc,
                    "ci95_delta": q(bc, 0.025, 0.975),
                    "relative_to_Y0": dc / y0,
                    "relative_to_Y2": dc / ys["Y2"] if ys.get("Y2", 0) > 0 else None,
                }
        out["composition"] = comp
        res[kind] = out
    res["strata"] = strata(d, y, preds, split, lex, cont)
    return res


def strata(
    d: dict[str, np.ndarray],
    y: np.ndarray,
    preds: dict[str, np.ndarray],
    split: str,
    lex: Lexicon,
    cont: dict[str, str],
) -> dict:
    """H5 の層（記述）: 層内の ΔBS(M1c → M3n)・90% 区間（層内のクラスタ bootstrap）・層内の物差しに対する比。"""
    layers = {"built_up_gt_10pct": d["s_built"] > wa.BUILT_UP}
    conts = np.array([cont.get(c, "??") for c in d["cc"]])
    langs = np.array([LANG.get(c, "??") for c in d["cc"]])
    out = {}
    for name, mask in layers.items():
        for label, m in ((name, mask), (f"not_{name}", ~mask)):
            out[label] = stratum(d, y, preds, m)
    for c in sorted(set(conts.tolist())):
        out[f"continent_{c}"] = stratum(d, y, preds, conts == c)
    for lang in sorted(set(langs.tolist())):  # 最寄り集落の国の第一言語（B-2 追記 2）
        out[f"language_{lang}"] = stratum(d, y, preds, langs == lang)
    return {k: v for k, v in out.items() if v is not None} | {"split": split}


def first_languages(path: Path) -> dict[str, str]:
    """GeoNames countryInfo.txt → {ISO2: 第一言語の ISO 639 コード}（例 en-AU → en）。"""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) > 15 and f[15]:
            out[f[0]] = f[15].split(",")[0].split("-")[0]
    return out


def stratum(
    d: dict[str, np.ndarray], y: np.ndarray, preds: dict[str, np.ndarray], m: np.ndarray
) -> dict | None:
    if m.sum() < wa.MIN_CONT_CELLS or y[m].sum() < wa.MIN_CONT_WET:
        return None
    sub = {a: preds[a][m] for a in ("M0", "M1", "M1c", "M3n")}
    s, n, arms, _ = sums(y[m], sub, d["unit"][m])
    j = {a: i for i, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    boot = wa.bootstrap_bs(s, n)
    dn = float(hat[j["M1c"]] - hat[j["M3n"]])
    y0 = float(hat[j["M0"]] - hat[j["M1"]])
    c90 = q(boot[:, j["M1c"]] - boot[:, j["M3n"]], 0.05, 0.95)
    return {
        "n_cells": int(m.sum()),
        "n_units": len(n),
        "wet_rate": float(y[m].mean()),
        "delta_bs_names": dn,
        "ci90_delta": c90,
        "yardstick_Y0": y0,
        "R": dn / y0 if y0 > 0 else None,
        "ci90_R": [c / y0 for c in c90] if y0 > 0 else None,
    }


# %% 実行間の幅・Holm・決着しなかった比較
def spread() -> dict:
    out = {}
    for run in SPREAD_RUNS:
        p = RUNS_DIR / run / "metrics_primary.json"
        if not p.is_file():
            continue
        met = json.loads(p.read_text())
        m, pl = met["primary"], met.get("placebo", {})
        out[run] = {
            k: m[k]
            for k in (
                "R",
                "ci90_R",
                "delta_bs_M0_M1",
                "margin_delta_bs",
                "decision",
                "equivalent_TOST",
            )
        }
        out[run] |= {
            "R_minus_Rp": pl.get("R_minus_Rp"),
            "ci95_R_minus_Rp": pl.get("ci95_R_minus_Rp"),
            "n_units": met["n_units"],
            "n_cells": met["n_cells"],
            "brier": {
                a: met["brier"][a]
                for a in ("M0", "M1", "M1s", "M1c", "M3n", "M3np")
                if a in met["brier"]
            },
        }
    return out


def country_leakage() -> dict:
    """国外し分割で、最寄りの集落が評価側の国にある訓練セル・訓練側の国にある評価セルの数と、
    複数の fold にまたがる activation の数（審査の指摘 6b、記述）。"""
    d = load("world_v2", "country")
    units = {
        int(k): set(np.unique(d["unit"][d["fold"] == k]).tolist())
        for k in np.unique(d["fold"])
    }
    all_units = set().union(*units.values())
    test_out = sum(
        int(np.sum(np.isin(d["cc"][d["fold"] == k], list(all_units - u))))
        for k, u in units.items()
    )
    train_in = sum(
        int(np.sum((d["fold"] != k) & np.isin(d["cc"], list(u))))
        for k, u in units.items()
    )
    act = np.array([a.split("_")[0] for a in d["aoi"]])
    span = sum(1 for a in np.unique(act) if len(np.unique(d["fold"][act == a])) > 1)
    return {
        "n_cells": len(d["y"]),
        "test_cells_nearest_settlement_in_training_country": test_out,
        "training_cells_nearest_settlement_in_held_out_country_summed_over_folds": train_in,
        "activations_spanning_folds": span,
        "n_activations": len(np.unique(act)),
    }


def holm(p: dict[str, float]) -> dict[str, float]:
    keys = sorted(p, key=p.get)
    adj, run = {}, 0.0
    for i, k in enumerate(keys):
        run = max(run, min(1.0, (len(keys) - i) * p[k]))
        adj[k] = run
    return adj


def inconclusive() -> dict:
    """決着しなかった比較: 達成した上限と、SESOI を除外するのに要る単位数の目安 n × (半幅 / (SESOI − |R̂|))²。"""
    out = {}
    for split in SPLITS:
        m = json.loads((RUNS_DIR / "world_v2" / f"metrics_{split}.json").read_text())
        for key in ("primary", "H4", "total", "H4_total", "flexible", "H4_flexible"):
            v = m.get(key)
            if not v or v["decision"] not in ("inconclusive", "positive-small"):
                continue
            hw = (v["ci90_R"][1] - v["ci90_R"][0]) / 2
            gap = SESOI - abs(v["R"])
            out[f"{split}/{key}"] = {
                "decision": v["decision"],
                "R": v["R"],
                "ci90_R": v["ci90_R"],
                "ci95_R": v["ci95_R"],
                "n_units": m["n_units"],
                "units_needed_to_exclude_SESOI": float(m["n_units"] * (hw / gap) ** 2)
                if gap > 0
                else None,
            }
    return out


# %% 表示
def main() -> None:
    lex = Lexicon(ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    LANG.update(first_languages(DATA_DIR / "global" / "geonames" / "countryInfo.txt"))
    res = {split: analyse(split, lex, cont) for split in SPLITS}
    # B-2 追記 2: activation で束ねた分割、AOI 内で中心化した名前の特徴（M3d）の実行
    for key, run, split, recal in (
        ("primary_actgrp", "world_v2_actgrp", "primary", "world_v2_actgrp"),
        *((f"{s}_dm", "world_v2_dm", s, "world_v2_recal") for s in SPLITS),
        ("primary_dm_actgrp", "world_v2_dm_actgrp", "primary", "world_v2_actgrp"),
    ):
        r = analyse(split, lex, cont, run=run, recal_run=recal)
        if r is not None:
            res[key] = r | {"run": run}
    for (
        key,
        r,
    ) in (
        res.items()
    ):  # 登録どおりの二方向 bootstrap（比・セルの国・1,000 回）と登録の判定
        run = r.get("run", "world_v2")
        split = key.split("_")[0]
        m = RUNS_DIR / run / f"metrics_{split}.json"
        if m.is_file():
            met = json.loads(m.read_text())
            r["prespecified_two_way"] = met.get("two_way")
            r["prespecified_decision"] = (met.get("primary") or {}).get("decision")
    for split, r in res.items():
        print(
            f"\n== {split}: {r['n_cells']:,} cells, {r['n_units']} units; recal {r['recal'].get('available')}"
            f"; two-way axes {r['two_way_axes']}"
        )
        if r["recal"].get("max_abs_diff_vs_world_v2"):
            worst = max(r["recal"]["max_abs_diff_vs_world_v2"].values())
            print(f"   recal vs world_v2: max |Δp| over shared arms = {worst:.2e}")
        for kind in ("cluster", "two_way"):
            ys = r[kind]["yardsticks"]
            print(
                f"  [{kind}] yardsticks "
                + "  ".join(f"{k}={v:.2e}" for k, v in ys.items())
            )
            for key, *_ in COMPARISONS:
                row = r[kind].get(key)
                if row is None:
                    continue
                cells = []
                for yk in ys:
                    e = row[yk]
                    cells.append(
                        f"{yk}:{'eq' if e['equivalent'] else '--'}"
                        + (
                            f"[{e['ci90_R'][0]:+.3f},{e['ci90_R'][1]:+.3f}]"
                            if e["ci90_R"]
                            else "[undef]"
                        )
                    )
                print(f"    {key:12s} ΔBS {row['delta_bs']:+.2e}  " + " ".join(cells))
            for k, v in r[kind]["composition"].items():
                rel2 = (
                    f"{v['relative_to_Y2']:+.3f}"
                    if v["relative_to_Y2"] is not None
                    else "  —"
                )
                print(
                    f"    comp {k:10s} ΔBS {v['delta_bs']:+.2e} [{v['ci95_delta'][0]:+.2e}, {v['ci95_delta'][1]:+.2e}]"
                    f"  /Y0 {v['relative_to_Y0']:+.3f}  /Y2 {rel2}"
                )
    llm = res["llmcut"]["cluster"].get("H4", {})
    lad = json.loads((RUNS_DIR / "world_v2" / "metrics_primary.json").read_text())[
        "ladder"
    ]
    p = {"H3_delta0_gt_delta4": lad["p_one_sided_delta0_gt_delta4"]}
    if llm:
        p["H4_TOST"] = llm["p_TOST_Y0"]
    extra = {
        "holm_H3_H4": {"raw": p, "adjusted": holm(p)},
        "inconclusive": inconclusive(),
        "spread_primary": spread(),
        "country_leakage": country_leakage(),
    }
    print("country split leakage:", extra["country_leakage"])
    print("\nHolm (H3, H4):", extra["holm_H3_H4"])
    print("inconclusive / positive-small:")
    for k, v in extra["inconclusive"].items():
        print(
            f"   {k:22s} {v['decision']:15s} R {v['R']:+.4f} 90% [{v['ci90_R'][0]:+.4f}, {v['ci90_R'][1]:+.4f}]"
            f"  units needed ≈ {v['units_needed_to_exclude_SESOI']}"
        )
    print("spread of the primary R across runs:")
    for k, v in extra["spread_primary"].items():
        print(
            f"   {k:24s} R {v['R']:+.4f} 90% [{v['ci90_R'][0]:+.4f}, {v['ci90_R'][1]:+.4f}]"
            f"  Y0 {v['delta_bs_M0_M1']:.2e}  {v['decision']}"
        )
    OUT.write_text(json.dumps({"splits": res, **extra}, indent=1), encoding="utf-8")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
