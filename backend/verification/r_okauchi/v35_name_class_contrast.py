# %% v35（探索的・凍結後）: 水辺の名を持つ集落の近くは、予測より多く浸水しているか
"""確証解析（world_v2、主分割）の held-out 予測だけを使う記述解析。事前登録外（§6 の H5 の層別を
AOI 内の対比に直したもの）で、主判定は変えない。

各 AOI の中で、最寄りの集落名が辞書で水辺系のセル（s_wetside = 1）とそれ以外のセルを比べる:
- 観測の浸水率の差（生の対比）
- 観測 − 予測の差（予測で説明されない残差の対比）。予測は M0c（位置 + 名前の地域構成）と
  M1c（+ 地形・水文・土地被覆・集落配置）、どちらも名前の場所ごとの特徴は見ていない
両群に 20 セル以上ある AOI について、小さい群のセル数で重み付けして平均する（AOI 内の対比の
精度に比例）。区間は場所グループの cluster bootstrap（4,000 回、主解析と同じ単位）。
辞書に当たった名前（s_lex_any）についても同じ対比を出し、物差しとして地形の対比（AOI 内で HAND が
中央値以下のセルとそれ以外）も同じ手順で出す。

AOI 間の関係（地域の水辺名の多さと浸水の多さ）: AOI ごとの水辺系の割合と観測の浸水率（および
観測 − M1 の予測の平均）の Spearman 相関。区間は場所グループの cluster bootstrap。

入力: data/runs/world_v2/predictions/primary/fold*.npz
出力: stdout の表と docs/verification/v35_name_class_contrast.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
PRED = ROOT / "data" / "runs" / "world_v2" / "predictions" / "primary"
MIN_CELLS = 20
B = 4000
SEED = 20261026
ARMS = ("M0c", "M1c")


# %%
def load() -> dict[str, np.ndarray]:
    parts = [np.load(f) for f in sorted(PRED.glob("fold*.npz"))]
    keys = ("y", "aoi", "unit", "s_wetside", "s_lex_any", "s_hand", "p_M1",
            *(f"p_{a}" for a in ARMS))
    d = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    # 物差し: AOI 内で HAND が AOI の中央値以下のセル（低い側）
    d["low_hand"] = np.zeros(len(d["y"]), dtype=np.int8)
    for a in np.unique(d["aoi"]):
        m = d["aoi"] == a
        d["low_hand"][m] = (d["s_hand"][m] <= np.nanmedian(d["s_hand"][m])).astype(np.int8)
    return d


def per_aoi(d: dict[str, np.ndarray], flag: str) -> dict[str, np.ndarray]:
    """AOI ごとの対比（水辺系 − それ以外）と重み。"""
    rows = []
    for a in np.unique(d["aoi"]):
        m = d["aoi"] == a
        g = d[flag][m] > 0
        if g.sum() < MIN_CELLS or (~g).sum() < MIN_CELLS:
            continue
        y = d["y"][m].astype(np.float64)
        row = {"aoi": a, "unit": d["unit"][m][0], "w": float(min(g.sum(), (~g).sum())),
               "raw": y[g].mean() - y[~g].mean(), "rate_w": y[g].mean(), "rate_o": y[~g].mean()}
        for arm in ARMS:
            r = y - d[f"p_{arm}"][m].astype(np.float64)
            row[arm] = r[g].mean() - r[~g].mean()
        rows.append(row)
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def boot(t: dict[str, np.ndarray], key: str) -> tuple[float, list[float]]:
    units, inv = np.unique(t["unit"], return_inverse=True)
    rng = np.random.default_rng(SEED)
    est = float(np.average(t[key], weights=t["w"]))
    reps = []
    for _ in range(B):
        k = rng.multinomial(len(units), np.full(len(units), 1 / len(units)))[inv]
        w = t["w"] * k
        if w.sum() > 0:
            reps.append(float(np.average(t[key], weights=w)))
    return est, [float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))]


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman の順位相関（同順位は平均順位。水辺系 0 の AOI が多い）。"""
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def between_aoi(d: dict[str, np.ndarray]) -> dict:
    """AOI 間: 水辺系の割合 × 観測の浸水率 / 観測 − M1 の平均。場所グループの bootstrap。"""
    aois = np.unique(d["aoi"])
    share, obs, res, unit = [], [], [], []
    for a in aois:
        m = d["aoi"] == a
        share.append(d["s_wetside"][m].mean())
        obs.append(d["y"][m].mean())
        res.append(d["y"][m].mean() - d["p_M1"][m].astype(np.float64).mean())
        unit.append(d["unit"][m][0])
    share, obs, res, unit = map(np.array, (share, obs, res, unit))
    units, inv = np.unique(unit, return_inverse=True)
    rng = np.random.default_rng(SEED)
    out: dict = {"n_aoi": len(aois), "share": share.tolist(), "obs": obs.tolist(),
                 "resid_M1": res.tolist()}
    for key, v in (("rho_obs", obs), ("rho_resid_M1", res)):
        reps = []
        for _ in range(B):
            pick = np.concatenate([np.flatnonzero(inv == u) for u in
                                   rng.integers(0, len(units), len(units))])
            reps.append(spearman(share[pick], v[pick]))
        out[key] = {"estimate": spearman(share, v),
                    "ci95": [float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))]}
    q = np.quantile(share, [1 / 3, 2 / 3])
    out["tertiles"] = [{"share": float(share[s].mean()), "obs": float(obs[s].mean())}
                       for s in (share <= q[0], (share > q[0]) & (share <= q[1]), share > q[1])]
    return out


def main() -> None:
    d = load()
    out: dict = {"n_cells": len(d["y"]), "min_cells_per_group": MIN_CELLS, "bootstrap": B}
    for flag, name in (("s_wetside", "wetside"), ("s_lex_any", "lexicon_any"),
                       ("low_hand", "low_hand")):
        t = per_aoi(d, flag)
        res = {"n_aoi": len(t["aoi"]), "n_units": len(np.unique(t["unit"])),
               "mean_rate_flagged": float(np.average(t["rate_w"], weights=t["w"])),
               "mean_rate_other": float(np.average(t["rate_o"], weights=t["w"]))}
        for key in ("raw", *ARMS):
            est, c = boot(t, key)
            res[key] = {"estimate": est, "ci95": c}
        res["per_aoi"] = {k: (v.tolist() if v.dtype.kind == "f" else v.astype(str).tolist())
                          for k, v in t.items()}
        out[name] = res
        print(f"{name}: {res['n_aoi']} AOIs in {res['n_units']} places; flood rate "
              f"{res['mean_rate_flagged']:.4f} (flagged) vs {res['mean_rate_other']:.4f} (other)")
        for key in ("raw", *ARMS):
            e, c = res[key]["estimate"], res[key]["ci95"]
            lab = "observed" if key == "raw" else f"observed - {key}"
            print(f"  {lab:18s} difference {e:+.4f}  95% [{c[0]:+.4f}, {c[1]:+.4f}]")
    b = between_aoi(d)
    out["between_aoi"] = b
    for key in ("rho_obs", "rho_resid_M1"):
        e, c = b[key]["estimate"], b[key]["ci95"]
        print(f"between AOIs ({b['n_aoi']}): Spearman(water-side share, {key[4:]}) = {e:+.3f}  "
              f"95% [{c[0]:+.3f}, {c[1]:+.3f}]")
    print("  tertiles of water-side share:", [(round(r["share"], 3), round(r["obs"], 4))
                                              for r in b["tertiles"]])
    path = ROOT / "docs" / "verification" / "v35_name_class_contrast.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
