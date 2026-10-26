"""v2 全球検定の集計（preregistration_v2.md §6）と検出力シミュレーションの集計（§8）。

- ``--tag <run> --split primary``: predictions/<split>/fold*.npz をプールし、
  プール ΔBS・比 R・場所グループのクラスタ bootstrap・TOST・placebo 条件・梯子・H1・
  大陸別内訳・二方向頑健性を summary.json / metrics.json / diagnostics.png に出す
- ``--tag <run> --power``: sim/rep*.npz（合成ラベル）ごとに同じ判定を行い検出力を出す
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.evaluation.metrics import brier_score, log_loss, pr_auc, top_area_capture
from topsem.runs import start_run, write_json
from topsem.world.dataset import continent_map

SESOI = 0.05
B = 4000
SEED = 20261025
LADDER = [("T0", "T0N"), ("T1", "T1N"), ("T2", "T2N"), ("T3", "T3N"), ("M1", "M3")]  # 柔軟（副）
LADDER_S = [("T0c", "T0n"), ("T1c", "T1n"), ("T2c", "T2n"), ("T3c", "T3n"), ("M1c", "M3n")]
MIN_CONT_CELLS = 1000  # 大陸別内訳を出す最小セル数（記述のみ）
MIN_CONT_WET = 20
BUILT_UP = 0.1  # H5 の市街地層: CCI 建物群の被覆率 > 10%


def load_split(run_dir: Path, split: str) -> dict[str, np.ndarray]:
    parts = sorted((run_dir / "predictions" / split).glob("fold*.npz"))
    if not parts:
        msg = f"no predictions under {run_dir}/predictions/{split}"
        raise FileNotFoundError(msg)
    data: dict[str, list[np.ndarray]] = {}
    for p in parts:
        with np.load(p, allow_pickle=False) as z:
            for k in z.files:
                data.setdefault(k, []).append(z[k])
    return {k: np.concatenate(v) for k, v in data.items()}


def unit_sums(y: np.ndarray, preds: dict[str, np.ndarray], unit: np.ndarray,
              ) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """単位ごとの二乗誤差和（[G, arms]）とセル数（[G]）。"""
    codes, uniq = _codes(unit)
    arms = sorted(preds)
    s = np.zeros((len(uniq), len(arms)))
    for j, a in enumerate(arms):
        np.add.at(s[:, j], codes, (preds[a].astype(np.float64) - y) ** 2)
    n = np.bincount(codes, minlength=len(uniq)).astype(np.float64)
    return s, n, arms


def _codes(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uniq, codes = np.unique(x, return_inverse=True)
    return codes, uniq


def bootstrap_bs(s: np.ndarray, n: np.ndarray, b: int = B, seed: int = SEED) -> np.ndarray:
    """単位の多項重みで再標本化したプール BS（[b, arms]）。"""
    rng = np.random.default_rng(seed)
    g = len(n)
    w = rng.multinomial(g, np.full(g, 1.0 / g), size=b).astype(np.float64)
    return (w @ s) / (w @ n)[:, None]


def ci(x: np.ndarray, level: float) -> list[float]:
    a = (1 - level) / 2
    return [float(np.quantile(x, a)), float(np.quantile(x, 1 - a))]


def decide(r_hat: float, ci95: list[float], ci90: list[float], *,
           placebo_ci95: list[float] | None = None,
           den_ci95: list[float] | None = None) -> dict:
    """§6 の判定表（凍結）。上から順に最初に当てはまる label を返す。

    - undefined-denominator: ΔBS(M0→M1) の 95% CI 下限 ≤ 0（R が定義できない）
    - not-meaning: R の 95% CI 下限 > 0 だが placebo 条件（R − R_p の 95% CI 下限 > 0）を満たさない
    - positive-meaningful / positive-negligible / positive-small: R の 95% CI 下限 > 0 かつ placebo
      条件を満たす。R̂ ≥ SESOI なら meaningful、90% CI ⊂ (−SESOI, SESOI) なら negligible
    - negative: R の 95% CI 上限 < 0（名前で悪化）
    - equivalent-null: 上のいずれでもなく 90% CI ⊂ (−SESOI, SESOI)（TOST）
    - inconclusive: それ以外
    equivalent は TOST の成否（label と独立に報告）。
    """
    equivalent = bool(ci90[0] > -SESOI and ci90[1] < SESOI)
    placebo_ok = None if placebo_ci95 is None else bool(placebo_ci95[0] > 0)
    if den_ci95 is not None and den_ci95[0] <= 0:
        label = "undefined-denominator"
    elif ci95[0] > 0 and placebo_ok is False:
        label = "not-meaning"
    elif ci95[0] > 0:
        label = ("positive-meaningful" if r_hat >= SESOI else
                 "positive-negligible" if equivalent else "positive-small")
    elif ci95[1] < 0:
        label = "negative"
    elif equivalent:
        label = "equivalent-null"
    else:
        label = "inconclusive"
    return {"label": label, "positive": label.startswith("positive"),
            "equivalent": equivalent, "placebo_ok": placebo_ok}


def make_delta(hat: np.ndarray, boot: np.ndarray, idx: dict[str, int]):  # noqa: ANN201
    """腕の対 (base, alt) → (点推定 BS(base) − BS(alt), bootstrap 再標本の同差)。"""

    def delta(base: str, alt: str) -> tuple[float, np.ndarray]:
        return (float(hat[idx[base]] - hat[idx[alt]]), boot[:, idx[base]] - boot[:, idx[alt]])

    return delta


# 主（入れ子の増分モデル）と副（柔軟モデル同士）の比較の組: (key, 基準, 名前あり, placebo, placebo の key)
# 主 = 場所内の増分（名前の地域構成を基準側に入れた入れ子）、total = 地域構成込みの全体の増分、
# flexible = 柔軟モデル同士（副）
COMPARISONS = [
    ("primary", "M1c", "M3n", "M3np", "placebo"),
    ("H4", "M1rc", "M3rn", "M3rnp", "H4_placebo"),
    ("total", "M1s", "M3s", "M3ps", "total_placebo"),
    ("H4_total", "M1s", "M3rs", "M3rps", "H4_total_placebo"),
    ("flexible", "M1", "M3", "M3p", "flexible_placebo"),
    ("H4_flexible", "M1", "M3r", "M3rp", "H4_flexible_placebo"),
]
CHANNELS = {"channels": [("M1c", "M3nx"), ("M1c", "M3nl")],
            "channels_flexible": [("M1", "M3x"), ("M1", "M3l")]}
H1_PAIRS = {"H1": ("M0c", "M0n"), "H1_flexible": ("M0", "M0N")}


def name_verdicts(delta, idx: dict[str, int], d01: float,  # noqa: ANN001
                  b01: np.ndarray) -> dict:
    """名前の判定一式（H2・H4・副比較・チャネル分解・H1）。

    主は場所内の増分 M1c → M3n、H4 は M1rc → M3rn。分母は全て ΔBS(M0 → M1)（柔軟モデル）。
    """
    out: dict = {}
    for key, base, alt, plac_arm, pkey in COMPARISONS:
        if base in idx and alt in idx:  # H4 は llmcut 分割で読む
            out[key], plac = verdict(delta, base, alt, plac_arm if plac_arm in idx else None,
                                     d01=d01, b01=b01)
            if plac is not None:
                out[pkey] = plac
    for key, pairs in CHANNELS.items():
        for base, ch in pairs:
            if base in idx and ch in idx:
                dc, bc = delta(base, ch)
                out.setdefault(key, {})[ch] = {"delta_bs": dc, "R": dc / d01,
                                               "ci95_R": ci(bc / b01, 0.95)}
    for key, (base, alt) in H1_PAIRS.items():
        if base in idx and alt in idx:
            dn, bn = delta(base, alt)
            out[key] = {f"delta_bs_{base}_{alt}": dn, "ratio_to_terrain": dn / d01,
                        "ci95_ratio": ci(bn / b01, 0.95)}
    return out


def primary_pair(preds: dict[str, np.ndarray]) -> tuple[str, str]:
    """記述（大陸別・層別・二方向・診断図）に使う比較: 主（場所内の増分）を優先する。"""
    return next((pair for pair in (("M1c", "M3n"), ("M1s", "M3s"), ("M1", "M3"))
                 if all(a in preds for a in pair)), ("M1", "M3"))


def evaluate(d: dict[str, np.ndarray], cont: dict[str, str]) -> dict:
    y = d["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in d.items() if k.startswith("p_")}
    s, n, arms = unit_sums(y, preds, d["unit"])
    idx = {a: j for j, a in enumerate(arms)}
    bs_hat = s.sum(axis=0) / n.sum()
    delta = make_delta(bs_hat, bootstrap_bs(s, n), idx)
    out: dict = {"n_cells": len(y), "n_units": len(n), "wet_rate": float(y.mean()),
                 "brier": {a: float(bs_hat[idx[a]]) for a in arms}}
    d01, b01 = delta("M0", "M1")
    out |= name_verdicts(delta, idx, d01, b01)
    for key, rungs in (("ladder", LADDER_S), ("ladder_flexible", LADDER)):
        if all(a in idx and b in idx for a, b in rungs) and "T0" in idx:
            out[key] = ladder(delta, rungs)
    out["secondary_metrics"] = secondary_metrics(d["y"], preds)
    num = primary_pair(preds)
    if d.get("cc") is not None and all(a in preds for a in num):
        out["by_continent"] = by_continent(d, preds, cont, num)
    # 二方向頑健性（activation × 国、それぞれを独立に再標本化）
    act = np.array([a.split("_")[0] for a in d["aoi"]])
    out["two_way"] = two_way_ci(y, preds, act, d.get("cc", act), num=num)
    if all(a in preds for a in ("M0", "M1", *num)):
        out["H5_strata"] = strata(d, preds, num)
    return out


def verdict(delta, base: str, alt: str, placebo: str | None, *,  # noqa: ANN001, PLR0913
            d01: float, b01: np.ndarray) -> tuple[dict, dict | None]:
    """名前の増分 base → alt の比 R と判定（固定マージン法）。

    判定は ΔBS の尺度で、余白を SESOI × ΔBŜ(M0 → M1)（地形の増分の点推定）に固定して行う
    （非劣性試験の固定マージン法）。地形の増分の bootstrap 変動係数は 0.4–0.55 と大きく、比の
    bootstrap 区間は分母の揺らぎで膨らむため（pilot3）。比の bootstrap 区間は記述として併記。
    """
    d, b = delta(base, alt)
    c95, c90 = [x / d01 for x in ci(b, 0.95)], [x / d01 for x in ci(b, 0.90)]
    plac = None
    if placebo is not None:
        dp, bp = delta(base, placebo)
        plac = {f"delta_bs_{base}_{placebo}": dp, "R_placebo": dp / d01,
                "ci95_R_placebo": [x / d01 for x in ci(bp, 0.95)],
                "R_minus_Rp": (d - dp) / d01,
                "ci95_R_minus_Rp": [x / d01 for x in ci(b - bp, 0.95)]}
    den95 = ci(b01, 0.95)
    dec = decide(d / d01, c95, c90, placebo_ci95=plac["ci95_R_minus_Rp"] if plac else None,
                 den_ci95=den95)
    return ({"comparison": f"{base}->{alt}", f"delta_bs_{base}_{alt}": d,
             "ci95_delta": ci(b, 0.95), "ci90_delta": ci(b, 0.90), "delta_bs_M0_M1": d01,
             "ci95_delta_M0_M1": den95, "margin_delta_bs": SESOI * d01,
             "R": d / d01, "ci95_R": c95, "ci90_R": c90,
             "ci95_R_ratio_bootstrap": ci(b / b01, 0.95), "decision": dec["label"],
             "equivalent_TOST": dec["equivalent"], "placebo_ok": dec["placebo_ok"],
             "SESOI": SESOI}, plac)


def secondary_metrics(y: np.ndarray, preds: dict[str, np.ndarray]) -> dict:
    """副次評価量（v1 と同定義）。"""
    return {a: {"brier": brier_score(y, preds[a]), "log_loss": log_loss(y, preds[a]),
                "pr_auc": pr_auc(y, preds[a]), "capture10": top_area_capture(y, preds[a])}
            for a in ("M1c", "M3n", "M1s", "M3s", "M1", "M3") if a in preds}


def ladder(delta, rungs: list[tuple[str, str]]) -> dict:  # noqa: ANN001
    """H3: 梯子の各段での名前の増分と単調性・T0 での地形情報回収率（分母は柔軟 T0 → M1）。"""
    ds, rboot = [], []
    for a, b in rungs:
        dk, bk = delta(a, b)
        ds.append(dk)
        rboot.append(bk)
    rb = np.vstack(rboot)  # [5, B]
    d_t0_t4, b_t0_t4 = delta("T0", "M1")
    return {
        "rungs": [f"{a}->{b}" for a, b in rungs],
        "delta_names": ds,
        "ci95_delta_names": [ci(rb[k], 0.95) for k in range(len(rungs))],
        "delta0_minus_delta4": ds[0] - ds[-1],
        "ci95_delta0_minus_delta4": ci(rb[0] - rb[-1], 0.95),
        "p_one_sided_delta0_gt_delta4": float(np.mean(rb[0] - rb[-1] <= 0)),
        "spearman_k_delta": float(np.corrcoef(np.arange(len(ds)),
                                              np.argsort(np.argsort(ds)))[0, 1]),
        "terrain_info_recovered_at_T0": ds[0] / d_t0_t4,
        "ci95_recovered_at_T0": ci(rb[0] / b_t0_t4, 0.95),
    }


def by_continent(d: dict[str, np.ndarray], preds: dict[str, np.ndarray],
                 cont: dict[str, str], num: tuple[str, str]) -> dict:
    """大陸別の名前の増分と R（記述のみ・判定には使わない）。"""
    conts = np.array([cont.get(c, "??") for c in d["cc"]])
    out = {}
    for c in sorted(set(conts.tolist())):
        m = conts == c
        if m.sum() < MIN_CONT_CELLS or d["y"][m].sum() < MIN_CONT_WET:
            continue
        b1 = brier_score(d["y"][m], preds[num[0]][m])
        b3 = brier_score(d["y"][m], preds[num[1]][m])
        b0 = brier_score(d["y"][m], preds["M0"][m])
        bm1 = brier_score(d["y"][m], preds["M1"][m])
        out[c] = {"n": int(m.sum()), "delta_bs_names": b1 - b3,
                  "R": (b1 - b3) / (b0 - bm1) if b0 != bm1 else None}
    return out


def two_way_ci(y: np.ndarray, preds: dict[str, np.ndarray], e: np.ndarray,  # noqa: PLR0913
               w: np.ndarray, *, num: tuple[str, str] = ("M1s", "M3s"), b: int = 1000) -> dict:
    """activation × 国の二方向 bootstrap（両軸の多重度の積で重み付け）による R の CI。"""
    arms = ("M0", "M1", *num)
    if not all(a in preds for a in arms):
        return {}
    ce, ue = _codes(e)
    cw, uw = _codes(w)
    rng = np.random.default_rng(SEED + 1)
    err = {a: (preds[a].astype(np.float64) - y) ** 2 for a in arms}
    rs = []
    for _ in range(b):
        we = rng.multinomial(len(ue), np.full(len(ue), 1 / len(ue)))[ce]
        ww = rng.multinomial(len(uw), np.full(len(uw), 1 / len(uw)))[cw]
        wt = (we * ww).astype(np.float64)
        if wt.sum() == 0:
            continue
        bs = {a: float((wt * v).sum() / wt.sum()) for a, v in err.items()}
        den = bs["M0"] - bs["M1"]
        if den != 0:
            rs.append((bs[num[0]] - bs[num[1]]) / den)
    return {"comparison": f"{num[0]}->{num[1]}",
            "ci95_R": ci(np.array(rs), 0.95) if rs else None, "n_activations": len(ue),
            "n_countries": len(uw), "reps": len(rs)}


def strata(d: dict[str, np.ndarray], preds: dict[str, np.ndarray],
           num: tuple[str, str]) -> dict:
    """H5（探索的・記述のみ）: 層ごとの名前の増分 ΔBS と層内の地形の増分に対する比。"""
    if "s_lex_any" not in d:
        return {}
    y = d["y"].astype(np.float64)
    slope_med = float(np.nanmedian(d["s_slope"]))
    layers = {
        "nearest_name_in_lexicon": d["s_lex_any"] > 0,
        "nearest_name_wetside": d["s_wetside"] > 0,
        "low_relief": d["s_slope"] <= slope_med,
        "built_up_gt_10pct": d["s_built"] > BUILT_UP,
    }
    out = {"slope_median_deg": slope_med}
    for name, mask in layers.items():
        for label, m in ((name, mask), (f"not_{name}", ~mask)):
            if m.sum() < MIN_CONT_CELLS or y[m].sum() < MIN_CONT_WET:
                continue
            bs = {a: brier_score(y[m], preds[a][m]) for a in ("M0", "M1", *num)}
            den = bs["M0"] - bs["M1"]
            out[label] = {"n": int(m.sum()), "wet_rate": float(y[m].mean()),
                          "delta_bs_names": bs[num[0]] - bs[num[1]],
                          "R_within": (bs[num[0]] - bs[num[1]]) / den if den else None}
    return out


def diagnostics(d: dict[str, np.ndarray], preds: dict[str, np.ndarray], num: tuple[str, str],
                path: Path) -> None:
    """診断図: (a) 校正曲線 (b) 場所単位の名前の増分 (c) 場所単位の増分の地図。"""
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    y = d["y"].astype(np.float64)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    bins = np.quantile(preds[num[0]], np.linspace(0, 1, 21))
    for a in ("M0", num[0], num[1]):
        k = np.clip(np.digitize(preds[a], bins[1:-1]), 0, 19)
        mp = np.bincount(k, weights=preds[a], minlength=20) / np.maximum(np.bincount(k, minlength=20), 1)
        my = np.bincount(k, weights=y, minlength=20) / np.maximum(np.bincount(k, minlength=20), 1)
        ax[0].plot(mp, my, "o-", ms=3, label=a)
    ax[0].plot([0, 1], [0, 1], color="0.6", lw=0.8)
    ax[0].set(xlabel="predicted", ylabel="observed", title="(a) calibration (20 quantile bins)")
    ax[0].legend(fontsize=8)
    codes, uniq = _codes(d["unit"])
    n = np.bincount(codes).astype(float)
    e1 = np.bincount(codes, weights=(preds[num[0]] - y) ** 2) / n
    e3 = np.bincount(codes, weights=(preds[num[1]] - y) ** 2) / n
    gain = e1 - e3
    ax[1].scatter(n, gain, s=8, alpha=0.6)
    ax[1].axhline(0, color="0.5", lw=0.8)
    ax[1].set(xscale="log", xlabel="cells in place (unit)", ylabel=f"ΔBS {num[0]}→{num[1]}",
              title=f"(b) per-place increment ({len(uniq)} units)")
    if "lon" in d:
        lon = np.bincount(codes, weights=d["lon"]) / n
        lat = np.bincount(codes, weights=d["lat"]) / n
        lim = float(np.quantile(np.abs(gain), 0.95)) or 1e-6
        sc = ax[2].scatter(lon, lat, c=gain, cmap="RdBu", vmin=-lim, vmax=lim, s=12)
        fig.colorbar(sc, ax=ax[2], label="ΔBS")
        ax[2].set(xlabel="lon", ylabel="lat", title="(c) per-place increment, map")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


POWER_B = 1000  # §8: 検出力は縮約版（bootstrap 1000 回）


def power(run_dir: Path) -> dict:
    """合成ラベル各回に主解析と同じ判定（placebo 条件込み）を当て、判定の頻度を出す。"""
    reps = sorted((run_dir / "sim").glob("rep*.npz"))
    res = []
    for p in reps:
        with np.load(p, allow_pickle=False) as z:
            y = z["y"].astype(np.float64)
            preds = {a: z[a] for a in ("M0", "M1", "M1c", "M3n", "M3np", "M1s", "M3s", "M3ps")
                     if a in z.files}
            s, n, arms = unit_sums(y, preds, z["unit"])
        idx = {a: j for j, a in enumerate(arms)}
        delta = make_delta(s.sum(axis=0) / n.sum(), bootstrap_bs(s, n, b=POWER_B, seed=SEED), idx)
        d01, b01 = delta("M0", "M1")
        v, pl = verdict(delta, "M1c", "M3n", "M3np", d01=d01, b01=b01)
        vt, _ = verdict(delta, "M1s", "M3s", "M3ps", d01=d01, b01=b01)
        res.append({"rep": p.stem, "R": v["R"], "ci95": v["ci95_R"], "ci90": v["ci90_R"],
                    "decision": v["decision"], "equivalent": v["equivalent_TOST"],
                    "ci90_halfwidth": (v["ci90_R"][1] - v["ci90_R"][0]) / 2,
                    "R_minus_Rp": pl["R_minus_Rp"], "ci95_R_minus_Rp": pl["ci95_R_minus_Rp"],
                    "total": {"R": vt["R"], "ci95": vt["ci95_R"], "decision": vt["decision"]}})
    if not res:
        return {"n_reps": 0}
    dec = [x["decision"] for x in res]
    return {"n_reps": len(res), "bootstrap_B": POWER_B,
            "rate_positive": float(np.mean([x.startswith("positive") for x in dec])),
            "rate_equivalent_TOST": float(np.mean([x["equivalent"] for x in res])),
            "rate_equivalent_null_label": float(np.mean([x == "equivalent-null" for x in dec])),
            "decision_counts": {k: dec.count(k) for k in sorted(set(dec))},
            "mean_R": float(np.mean([x["R"] for x in res])),
            "median_ci90_halfwidth": float(np.median([x["ci90_halfwidth"] for x in res])),
            "reps": res}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--split", default="primary")
    ap.add_argument("--power", action="store_true")
    args = ap.parse_args()
    run_dir = RUNS_DIR / args.tag
    run = start_run(args.tag)
    if args.power:
        out = power(run_dir)
        write_json(run.path("power_summary.json"), {**run.provenance(), **out})
        print(json.dumps({k: v for k, v in out.items() if k != "reps"}, indent=1))
        return
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    d = load_split(run_dir, args.split)
    res = evaluate(d, cont)
    write_json(run.path(f"metrics_{args.split}.json"), res)
    preds = {k[2:]: v for k, v in d.items() if k.startswith("p_")}
    num = primary_pair(preds)
    if all(a in preds for a in ("M0", *num)):
        diagnostics(d, preds, num, run.path(f"diagnostics_{args.split}.png"))
        if args.split == "primary":  # CLAUDE.md の run 契約（diagnostics.png）
            diagnostics(d, preds, num, run.path("diagnostics.png"))
    write_json(run.path(f"summary_{args.split}.json"),
               {**run.provenance(), "split": args.split,
                **{k: res.get(k) for k in ("primary", "placebo", "H4", "H4_placebo", "flexible",
                                           "flexible_placebo", "ladder", "ladder_flexible",
                                           "H1", "H1_flexible")},
                "n_cells": res["n_cells"], "n_units": res["n_units"]})
    print(json.dumps({k: res.get(k) for k in ("primary", "placebo", "H4", "H4_placebo",
                                               "flexible")}, indent=1))


if __name__ == "__main__":
    main()
