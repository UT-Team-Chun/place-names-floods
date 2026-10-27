# %% v48（事後・amendment B-6 の 1 と 7）: 既存の予測からの追加の要約
"""確証判定は変えない。B-6 で計算の前に固定した量を、既存の held-out 予測から計算する。

- 1: 局所の効果の天井 R_max(S) = Σ_{i∈S} (p_i − y_i)² / N / 物差し（p は M1c の held-out 予測）。S は最寄りの
  名前が水系側の (a) 500 m 以内・(b) 250 m 以内・(c) 1 km 以内・(d) 全セル。凍結したコード（world_v2）と直した
  コード（world_v2_lexfix）
- 7: AOI に等しい重みの R（AOI ごとの Brier の差の平均の比、場所グループのクラスタ bootstrap 4,000 回）。座標が
  角度の分の整数に丸められた集落に最寄りのセルを除いた、集落の近くのオッズ比（EMS の全 AOI、GFM で評価できる
  セルの GFM と EMS）。登録した被覆率（B-2 の 6・追記 3）と δ = 0 の判定の分布の再掲

出力: stdout と docs/verification/v48_b6_local.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR, RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v48_b6_local.json"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
V41T = DATA_DIR / "verification" / "v41_table.parquet"
B = 4000
SEED = 20261028  # B-6
SESOI = 0.05
MINUTE_TOL = 2e-3  # v37 と同じ（GeoNames の座標は小数 5 桁）


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v39 = _load("v39_within_aoi_odds")
v41 = _load("v41_gfm_independent_labels")


def table() -> pd.DataFrame:
    """v41 の表（凍結したコードの予測、fold 順）に、v37 の表から集落の座標と M1 の予測を足し、lexfix の予測も。"""
    t = pd.read_parquet(V41T)
    s = pd.read_parquet(V37T, columns=["aoi", "y", "P_d1", "P_gid1", "set_lat", "set_lon", "cc"])
    if not (np.array_equal(t["aoi"].to_numpy(), s["aoi"].to_numpy())
            and np.array_equal(t["y"].to_numpy(), s["y"].to_numpy())
            and np.allclose(t["P_d1"].to_numpy(), s["P_d1"].to_numpy())):
        msg = "v37 and v41 tables are not aligned"
        raise SystemExit(msg)
    for c in ("P_gid1", "set_lat", "set_lon", "cc"):
        t[c] = s[c].to_numpy()
    if not v41.attach_lexfix(t):
        msg = "world_v2_lexfix predictions missing"
        raise SystemExit(msg)
    # 水系側の印は、その実行の辞書の読み（予測のファイルの s_wetside = 最寄りの集落の lex_wetside_n1）
    for tag, col in (("world_v2", "wetside_frozen"), ("world_v2_lexfix", "wetside_lexfix")):
        pred = RUNS_DIR / tag / "predictions" / "primary"
        files = sorted(pred.glob("fold[0-9].npz"), key=lambda f: int(f.stem[4:]))
        t[col] = np.concatenate([np.load(f)["s_wetside"] for f in files]) > 0
    if not np.array_equal(t["wetside_frozen"].to_numpy(), t["lex_wetside_n1"].to_numpy() > 0):
        msg = "frozen wetside flags differ from the v41 table"
        raise SystemExit(msg)
    lat60, lon60 = t["set_lat"] * 60, t["set_lon"] * 60
    t["minute_rounded"] = ((lat60 - lat60.round()).abs() < MINUTE_TOL) & (
        (lon60 - lon60.round()).abs() < MINUTE_TOL)
    return t


# %% 1. 局所の効果の天井
def ceiling(t: pd.DataFrame, prefix: str, wetcol: str) -> dict:
    y = t["y"].to_numpy(np.float64)
    n = len(y)
    bs = {a: float(np.mean((t[f"{prefix}{a}"].to_numpy(np.float64) - y) ** 2)) for a in ("M0", "M1", "M1c")}
    yard = bs["M0"] - bs["M1"]
    err = (t[f"{prefix}M1c"].to_numpy(np.float64) - y) ** 2
    wet = t[wetcol].to_numpy(dtype=bool)
    d1 = t["P_d1"].to_numpy()
    out = {"yardstick": yard, "n_cells": n}
    for key, m in (("le500m", wet & (d1 <= 500)), ("le250m", wet & (d1 <= 250)), ("le1km", wet & (d1 <= 1000)),
                   ("all_wetside", wet)):
        out[key] = {"n_cells": int(m.sum()), "share_cells": float(m.mean()), "n_flooded": int(y[m].sum()),
                    "flooded_share": float(y[m].mean()) if m.any() else None,
                    "delta_bs_max": float(err[m].sum() / n), "R_max": float(err[m].sum() / n / yard),
                    "reaches_sesoi": bool(err[m].sum() / n / yard >= SESOI)}
    return out


# %% 7. AOI に等しい重みの R
def aoi_weighted(t: pd.DataFrame, prefix: str) -> dict:
    y = t["y"].to_numpy(np.float64)
    g = t.assign(**{a: (t[f"{prefix}{a}"].to_numpy(np.float64) - y) ** 2 for a in ("M0", "M1", "M1c", "M3n")})
    per = g.groupby("aoi").agg(M0=("M0", "mean"), M1=("M1", "mean"), M1c=("M1c", "mean"),
                               M3n=("M3n", "mean"), unit=("unit", "first"))
    num = (per["M1c"] - per["M3n"]).to_numpy()
    den = (per["M0"] - per["M1"]).to_numpy()
    r = float(num.mean() / den.mean())
    units, iu = np.unique(per["unit"].to_numpy(), return_inverse=True)
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(units), np.full(len(units), 1 / len(units)), size=B).astype(np.float64)
    wa = w[:, iu]  # 場所の重みを、その場所の AOI に
    rb = (wa @ num) / (wa @ den)
    q = lambda a, b: [float(np.quantile(rb, a)), float(np.quantile(rb, b))]
    return {"n_aoi": len(per), "n_units": len(units), "R": r, "ci95_R": q(0.025, 0.975), "ci90_R": q(0.05, 0.95),
            "equivalent_TOST": bool(-SESOI < np.quantile(rb, 0.05) and np.quantile(rb, 0.95) < SESOI),
            "delta_bs_names_mean": float(num.mean()), "yardstick_mean": float(den.mean())}


# %% 7. 丸めた座標を除いた集落の近くのオッズ比
def odds_wetside(t: pd.DataFrame, ycol: str) -> dict:
    y = t[ycol].to_numpy(dtype=np.float64)
    rate = t.assign(_y=y).groupby("aoi")["_y"].transform("mean").to_numpy()
    m = (rate > 0) & (rate < 1)
    s = t.loc[m].reset_index(drop=True)
    out = {}
    for key, banded, term in (("le500m", True, "wetside:le500m"), ("all", False, "wetside:all")):
        z, names = v39.design(s, ("wetside",), True, banded=banded)
        res = v39.fit(y[m], z, s["aoi"].to_numpy(), s["unit"].to_numpy())
        out[key] = v39.table(res, names)[term] | {"n_cells": int(m.sum()), "n_aoi": int(s["aoi"].nunique())}
    z, names = v39.design(s, tuple(v41.CATEGORIES), True, banded=True)  # 区分別（500 m 以内）
    tab = v39.table(v39.fit(y[m], z, s["aoi"].to_numpy(), s["unit"].to_numpy()), names)
    near = s["P_d1"].to_numpy() <= 500
    ym = y[m]
    out["categories_le500m"] = {}
    for c in v41.CATEGORIES:
        g = near & (s[f"lex_{c}_n1"].to_numpy() > 0)  # 群のセル数と浸水セル数（0 なら推定は退化）
        out["categories_le500m"][c] = tab[f"{c}:le500m"] | {"n_cells_group": int(g.sum()),
                                                            "n_flooded_group": int(ym[g].sum())}
    g = near & (s["lex_wetside_n1"].to_numpy() > 0)
    out["le500m"] |= {"n_cells_group": int(g.sum()), "n_flooded_group": int(ym[g].sum())}
    return out


def rounded(t: pd.DataFrame) -> dict:
    keep = ~t["minute_rounded"].to_numpy()
    out = {"share_cells_minute_rounded": float(1 - keep.mean())}
    e = t.assign(p_M1=t["p_M1"])
    out["ems_all"] = {"all_cells": odds_wetside(e, "y"), "not_rounded": odds_wetside(e.loc[keep], "y")}
    summ = json.loads((DATA_DIR / "global" / "gfm_summary.json").read_text(encoding="utf-8"))
    g = v41.attach_gfm(t, summ)
    kg = ~g["minute_rounded"].to_numpy()
    out["gfm_cells"] = {lab: {"all_cells": odds_wetside(g, col), "not_rounded": odds_wetside(g.loc[kg], col)}
                        for lab, col in (("gfm", "y_gfm"), ("ems", "y_ems"))}
    return out


def coverage() -> dict:
    v37 = json.loads((ROOT / "docs" / "verification" / "v37_positive_control.json").read_text(encoding="utf-8"))
    pw = json.loads((RUNS_DIR / "power_v3_null" / "power_summary.json").read_text(encoding="utf-8"))
    return {"null_coverage": v37.get("null_coverage"), "null_decision_counts": pw.get("decision_counts"),
            "null_rate_equivalent_TOST": pw.get("rate_equivalent_TOST"), "n_reps": pw.get("n_reps")}


def main() -> None:
    t = table()
    res = {"rule": "amendment B-6 items 1 and 7", "seed": SEED, "B": B,
           "ceiling": {"frozen": ceiling(t, "p_", "wetside_frozen"), "lexfix": ceiling(t, "p_lf_", "wetside_lexfix")},
           "aoi_weighted": {"frozen": aoi_weighted(t, "p_"), "lexfix": aoi_weighted(t, "p_lf_")},
           "rounded": rounded(t), "coverage": coverage()}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for run in ("frozen", "lexfix"):
        c = res["ceiling"][run]
        print(run, "ceiling", {k: (round(c[k]["R_max"], 4), c[k]["n_cells"]) for k in
                               ("le250m", "le500m", "le1km", "all_wetside")})
        a = res["aoi_weighted"][run]
        print(run, "AOI-weighted R", round(a["R"], 4), [round(x, 4) for x in a["ci90_R"]], a["equivalent_TOST"])
    r = res["rounded"]
    print("rounded share", round(r["share_cells_minute_rounded"], 3))
    for lab, blk in (("EMS all", r["ems_all"]), ("GFM", r["gfm_cells"]["gfm"]), ("EMS GFM cells", r["gfm_cells"]["ems"])):
        for k in ("all_cells", "not_rounded"):
            x = blk[k]["le500m"]
            print(lab, k, round(x["odds_ratio"], 3), [round(v, 3) for v in x["ci95_odds_ratio"]])


if __name__ == "__main__":
    main()
