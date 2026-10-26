# %% v42（事後・amendment B-5）: 初見の査読に応じた解析のうち、既存の成果物から計算するもの
"""確証判定は変えない。B-5 で計算の前に決めた量（4・6・8・9）を、world_v2 の held-out 予測と合成ラベルの記録から計算する。

- 4: 検定の反転によるオッズ比の上限（合成ラベルの各 δ の回の ΔBS の 5 パーセンタイルが観測の ΔBS に届く δ）。
  地域全体の効果（主分割の δ の格子）と、集落から 500 m 以内だけの効果（局所の格子）。経験分位と正規近似の両方
- 6: AOI ごとの AUC・PR-AUC の差（M3n − M1c）の分布、AOI ごとに切片を当て直した条件つき対数尤度の差（場所グループの
  bootstrap）、場所グループを 1 つずつ除いた R の範囲と寄与の大きい場所
- 8: 要素ごとの群（B-5 追記 1 で固定した形）の、集落から 500 m 以内と全セルのオッズ比（v39 と同じ回帰、地形で調整、
  Holm 補正）と、正の対照（HAND の平均・AOI 内の百分位、v37 の対比）
- 9: 南北アメリカとオセアニアの AOI を除いた主比較（v37 の band_increment）

入力: v37 と v41 の解析表（同じ行順、y と AOI の一致を確かめる）、GeoNames、凍結した辞書（B-4 の読み方）、
data/runs/power_v3*（合成ラベル）。区間は場所グループのクラスタ bootstrap（B = 4,000、シード 20261027）。
出力: stdout の表と docs/verification/v42_blind_review_analyses.json
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.dataset import continent_map
from topsem.world.lexicon import _SPLIT, Lexicon, _hit


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v37 = _load("v37_positive_control")
v39 = _load("v39_within_aoi_odds")
B = 4000
SEED = 20261027
v37.SEED, v37.B = SEED, B
OUT = ROOT / "docs" / "verification" / "v42_blind_review_analyses.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
GRID = {0.0: "power_v3_null", 0.2: "power_v3_d020", 0.4: "power_v3_d040", 0.6: "power_v3_d060",
        0.8: "power_v3_d080", 1.2: "power_v3_sesoi", 1.6: "power_v3_d160"}
GRID_LOCAL = {0.0: "power_v3_null", 1.0: "power_v3_loc_d100", 2.0: "power_v3_loc_d200",
              3.0: "power_v3_loc_d300"}
# B-5 追記 1 で固定した群（島・泉・堤は形で選び、湿地・渡り・低地は辞書の区分のまま）
ISLAND = set(["eiland", "djezira", "gezira", "jazira", "jazirat", "jezira", "chilumba", "chisumbu", "ostrov", "holm", "holmen", "insel", "werder", "wert", "werth", "woehrd", "woerth", "wohrd", "worth", "au", "aue", "auen", "nisaki", "nisi", "nisos", "ait", "eyot", "holme", "island", "isle", "isla", "islas", "saar", "saare", "luoto", "saari", "ile", "iles", "ilot", "sziget", "pulau", "isola", "isole", "isolotto", "kaoh", "koh", "aral", "araly", "sala", "salos", "salas", "nosy", "kyun", "waard", "waarden", "weerd", "oy", "oya", "ostrow", "ostrowek", "ostrowo", "kepa", "kepy", "ilha", "ilhas", "ilheu", "insua", "insula", "ostrovu", "ada", "ostrvo", "otok", "jasiirad", "ishull", "ishulli", "holmar", "kisiwa", "ko", "pulo", "adasi", "ostriv", "ostrovok", "con"])
SPRING = set(["fontein", "put", "ain", "aioun", "ayn", "ayoun", "abar", "bir", "biar", "hassi", "cheshma", "izvor", "izvori", "kladenci", "kladenets", "bunar", "pramen", "prameny", "studna", "studnice", "studnicka", "kilde", "born", "brunn", "brunnen", "quell", "quelle", "kefalovryso", "krini", "pigadi", "pigadia", "pigi", "vris", "vrises", "vrisi", "vrys", "vryses", "vrysi", "borne", "bourne", "spring", "springs", "well", "wells", "waterhole", "soak", "fuente", "fuentes", "manantial", "ojo", "pozo", "pozos", "allikas", "chah", "chahi", "chashma", "chashmeh", "lahde", "fontaine", "fontaines", "puits", "source", "sources", "forras", "kut", "sumur", "fontana", "fontane", "fonte", "pozzi", "pozzo", "sorgente", "bulak", "bulaq", "saltiniai", "saltinis", "avoti", "avots", "loharano", "ameyal", "studnia", "studzianka", "zrodla", "zrodlo", "cacimba", "fontes", "nascente", "poco", "pukyu", "puquio", "fantana", "fantanele", "izvoare", "izvorul", "studenac", "vrela", "vrelo", "vrutak", "izvir", "studenci", "studenec", "ceel", "burim", "burimi", "pus", "puse", "kalla", "kallan", "chemchemi", "kisima", "bukal", "kaynak", "kuyu", "pinar", "dzherelo", "kolodez", "kolodiazi", "krinitsa", "krynytsia", "gieng"])
DIKE = set(["wal", "hraz", "dige", "levee", "penger", "digue", "tanggul", "argine", "argini", "bonifica", "pylimas", "dambis", "dam", "dijk", "polder", "sluis", "dig", "nasip", "set"])
GROUPS = {"wetland": ("wet", None), "island": ("riverine", ISLAND), "spring_well": ("water", SPRING),
          "crossing": ("crossing", None), "dike_reclamation": ("riverine", DIKE), "lowland": ("low", None)}
AMERICAS_OCEANIA = ("NA", "SA", "OC")


# %% A. 解析表（v37 の表に v41 の M0・M3np を足す）
def table() -> pd.DataFrame:
    t = pd.read_parquet(DATA_DIR / "verification" / "v37_table.parquet")
    u = pd.read_parquet(DATA_DIR / "verification" / "v41_table.parquet", columns=["aoi", "y", "p_M0", "p_M3np"])
    if not (np.array_equal(t["y"].to_numpy(), u["y"].to_numpy())
            and np.array_equal(t["aoi"].to_numpy(dtype=str), u["aoi"].to_numpy(dtype=str))):
        msg = "v37 and v41 tables are not aligned"
        raise SystemExit(msg)
    t["p_M0"], t["p_M3np"] = u["p_M0"].to_numpy(), u["p_M3np"].to_numpy()
    return v37.derived(t)


# %% 4. 検定の反転
def inversion(grid: dict[float, str], observed: float) -> dict:
    rows = {}
    for d, run in grid.items():
        s = json.loads((RUNS_DIR / run / "power_summary.json").read_text(encoding="utf-8"))
        x = np.array([r["delta_bs"] for r in s["reps"]], dtype=np.float64)
        rows[d] = {"n": len(x), "q05": float(np.quantile(x, 0.05)), "mean": float(x.mean()),
                   "sd": float(x.std(ddof=1)), "q05_normal": float(x.mean() - 1.645 * x.std(ddof=1)),
                   "share_below_observed": float((x <= observed).mean())}

    def cross(key: str) -> float | None:
        ds = sorted(rows)
        for a, b in itertools.pairwise(ds):
            qa, qb = rows[a][key], rows[b][key]
            if qa <= observed < qb:
                return float(a + (observed - qa) / (qb - qa) * (b - a))
        return None if rows[ds[-1]][key] <= observed else float(ds[0])

    out = {"observed_delta_bs": observed, "grid": {str(k): v for k, v in rows.items()}}
    for key in ("q05", "q05_normal"):
        d = cross(key)
        out[f"delta_upper_{key}"] = d
        out[f"or_upper_{key}"] = float(np.exp(d)) if d is not None else None
    return out


# %% 6. 地域内の指標と場所ごとの影響
def _recal_ll(y: np.ndarray, p: np.ndarray, aoi: np.ndarray) -> np.ndarray:
    """AOI ごとに切片だけ当て直した対数尤度（セルごと）。"""
    z = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    out = np.zeros(len(y))
    for idx in pd.Series(np.arange(len(y))).groupby(aoi).indices.values():
        zi, yi = z[idx], y[idx]
        a = 0.0
        for _ in range(50):
            q = 1 / (1 + np.exp(-(zi + a)))
            g, h = (yi - q).sum(), (q * (1 - q)).sum()
            if h <= 0:
                break
            step = g / h
            a += step
            if abs(step) < 1e-10:
                break
        q = np.clip(1 / (1 + np.exp(-(zi + a))), 1e-12, 1 - 1e-12)
        out[idx] = yi * np.log(q) + (1 - yi) * np.log(1 - q)
    return out


def within_area(t: pd.DataFrame) -> dict:
    y, aoi = t["y"].to_numpy(), t["aoi"].to_numpy()
    rows = []
    for a, idx in pd.Series(np.arange(len(t))).groupby(aoi).indices.items():
        yi = y[idx]
        if yi.min() == yi.max():
            continue
        rows.append({"aoi": a, "unit": t["unit"].to_numpy()[idx[0]],
                     **{f"auc_{m}": roc_auc_score(yi, t[f"p_{m}"].to_numpy()[idx]) for m in ("M1c", "M3n")},
                     **{f"prauc_{m}": average_precision_score(yi, t[f"p_{m}"].to_numpy()[idx])
                        for m in ("M1c", "M3n")}})
    r = pd.DataFrame(rows)
    d_auc, d_pr = r["auc_M3n"] - r["auc_M1c"], r["prauc_M3n"] - r["prauc_M1c"]
    ll = {m: _recal_ll(y.astype(np.float64), t[f"p_{m}"].to_numpy(dtype=np.float64), aoi) for m in ("M1c", "M3n")}
    uniq, inv = np.unique(t["unit"].to_numpy(), return_inverse=True)
    per = np.bincount(inv, weights=ll["M3n"] - ll["M1c"])
    n = np.bincount(inv).astype(np.float64)
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(n), np.full(len(n), 1 / len(n)), size=B).astype(np.float64)
    bb = (w @ per) / (w @ n)
    # 場所グループを 1 つずつ除いた R（主比較と同じ ΔBS(M1c → M3n) / ΔBS(M0 → M1)）
    s = {m: np.bincount(inv, weights=(t[f"p_{m}"].to_numpy(dtype=np.float64) - y) ** 2) for m in ("M0", "M1", "M1c", "M3n")}
    tot = {m: v.sum() for m, v in s.items()}
    nn = n.sum()
    loo = []
    for g in range(len(uniq)):
        dn = ((tot["M1c"] - s["M1c"][g]) - (tot["M3n"] - s["M3n"][g])) / (nn - n[g])
        y0 = ((tot["M0"] - s["M0"][g]) - (tot["M1"] - s["M1"][g])) / (nn - n[g])
        loo.append(dn / y0)
    loo = np.array(loo)
    contrib = (s["M1c"] - s["M3n"]) / nn  # ΔBS への各場所の寄与
    top = np.argsort(-np.abs(contrib))[:5]
    return {"n_aoi_auc": len(r), "auc_diff_median": float(d_auc.median()),
            "auc_diff_iqr": [float(d_auc.quantile(0.25)), float(d_auc.quantile(0.75))],
            "auc_share_positive": float((d_auc > 0).mean()),
            "prauc_diff_median": float(d_pr.median()),
            "prauc_diff_iqr": [float(d_pr.quantile(0.25)), float(d_pr.quantile(0.75))],
            "cond_ll_per_cell": float(per.sum() / n.sum()),
            "cond_ll_ci95": [float(np.quantile(bb, 0.025)), float(np.quantile(bb, 0.975))],
            "loo_R_range": [float(loo.min()), float(loo.max())],
            "top_contributors": [{"unit": str(uniq[g]), "delta_bs_contribution": float(contrib[g]),
                                  "n_cells": int(n[g])} for g in top]}


# %% 8. 要素ごとの群
def group_flags(t: pd.DataFrame) -> dict[str, np.ndarray]:
    lex = Lexicon(LEX)
    gn = pd.read_parquet(DATA_DIR / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    gids = np.unique(t["P_gid1"].to_numpy())
    per: dict[str, set[int]] = {g: set() for g in GROUPS}
    for gid in gids:
        if gid < 0 or gid not in gn.index:
            continue
        name, cc = gn.at[gid, "asciiname"], gn.at[gid, "cc"]
        toks = [x for x in _SPLIT.split(str(name).lower()) if x]
        ents = lex.entries_for(str(cc))
        for g, (cat, forms) in GROUPS.items():
            if any(e.cat == cat and (forms is None or e.form in forms) and _hit(tok, e)
                   for e in ents for tok in toks):
                per[g].add(int(gid))
    p = t["P_gid1"].to_numpy()
    return {g: np.isin(p, list(s)) for g, s in per.items()}


def holm(ps: dict[str, float]) -> dict[str, float]:
    order = sorted(ps, key=ps.get)
    m, out, run = len(order), {}, 0.0
    for i, k in enumerate(order):
        run = max(run, min(1.0, (m - i) * ps[k]))
        out[k] = run
    return out


def groups(t: pd.DataFrame) -> dict:
    from scipy.stats import norm  # noqa: PLC0415

    fl = group_flags(t)
    y = t["y"].to_numpy(dtype=np.float64)
    rate = t.assign(_y=y).groupby("aoi")["_y"].transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = t.loc[keep].reset_index(drop=True)
    ys = y[keep]
    d1 = s["P_d1"].to_numpy()
    pm = np.clip(s["p_M1"].to_numpy(dtype=np.float64), 1e-5, 1 - 1e-5)
    lg = np.log(pm) - np.log(1 - pm)
    bands = [((d1 > 500) & (d1 <= 1000)).astype(float), (d1 > 1000).astype(float)]
    allm = np.ones(len(t), dtype=bool)
    near = t["P_d1"].to_numpy() <= 500
    out, pnear = {}, {}
    for g, f in fl.items():
        gi = f[keep]
        # 完全分離（群のセルに浸水が無い、または全て浸水）は推定できないので記す
        sep = {band: bool(m.sum() == 0 or ys[m].min() == ys[m].max())
               for band, m in (("le500m", gi & (d1 <= 500)), ("all", gi))}
        z_near = np.column_stack([*bands, gi & (d1 <= 500), gi & (d1 > 500) & (d1 <= 1000), gi & (d1 > 1000), lg])
        res = v39.fit(ys, z_near.astype(np.float64), s["aoi"].to_numpy(), s["unit"].to_numpy())
        tn = v39.table(res, ["band_500m_1km", "band_gt1km", "g:le500m", "g:500m_1km", "g:gt1km", "logit_M1"])
        res_all = v39.fit(ys, np.column_stack([gi, lg]).astype(np.float64), s["aoi"].to_numpy(), s["unit"].to_numpy())
        ta = v39.table(res_all, ["g:all", "logit_M1"])
        if not sep["le500m"]:
            b, se = tn["g:le500m"]["coef"], tn["g:le500m"]["se"]
            pnear[g] = float(2 * norm.sf(abs(b / se)))
        out[g] = {"n_cells": int(f.sum()), "n_cells_le500m": int((f & near).sum()),
                  "n_flooded_le500m": int(ys[gi & (d1 <= 500)].sum()), "separation": sep,
                  "n_settlements": len(np.unique(t["P_gid1"].to_numpy()[f])),
                  "or_le500m_adj": None if sep["le500m"] else tn["g:le500m"],
                  "or_all_adj": None if sep["all"] else ta["g:all"],
                  "pc_hand_mean_all": v37.contrast(t, "hand_mean", f, allm, "capped"),
                  "pc_hand_pct_all": v37.contrast(t, "hand_pct", f, allm, "capped"),
                  "pc_hand_mean_le500m": v37.contrast(t, "hand_mean", f, near, "capped")}
    for g, p in holm(pnear).items():
        out[g]["p_le500m"] = pnear[g]
        out[g]["p_le500m_holm"] = p
    return {"holm_n_groups": len(pnear), "by_group": out}


# %% 9. 入植地の国を除いた主比較
def settler(t: pd.DataFrame) -> dict:
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    cc = t.groupby("aoi")["cc"].agg(lambda s: s[s != ""].mode().iloc[0] if (s != "").any() else "??")
    aoi_cont = cc.map(lambda c: cont.get(c, "??"))
    drop = set(aoi_cont[aoi_cont.isin(AMERICAS_OCEANIA)].index)
    m = ~t["aoi"].isin(drop).to_numpy()
    return {"n_aoi_dropped": len(drop), "continents_dropped": list(AMERICAS_OCEANIA),
            **v37.band_increment(t, m)}


def main() -> None:
    t = table()
    obs = float(json.loads((RUNS_DIR / "world_v2" / "metrics_primary.json").read_text(encoding="utf-8"))
                ["primary"]["delta_bs_M1c_M3n"])
    res = {"seed": SEED, "B": B, "inversion_area_wide": inversion(GRID, obs),
           "inversion_local_500m": inversion(GRID_LOCAL, obs), "within_area": within_area(t),
           "groups": groups(t), "settler_excluded": settler(t)}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    show(res)


def show(res: dict) -> None:
    for k in ("inversion_area_wide", "inversion_local_500m"):
        r = res[k]
        print(f"{k}: OR upper (empirical q05) {r['or_upper_q05']}, (normal) {r['or_upper_q05_normal']}")
    w = res["within_area"]
    print(f"AUC diff median {w['auc_diff_median']:+.4f} IQR {w['auc_diff_iqr']}; cond LL/cell {w['cond_ll_per_cell']:+.2e} "
          f"{w['cond_ll_ci95']}; LOO R {w['loo_R_range']}")
    for g, v in res["groups"]["by_group"].items():
        n, a = v["or_le500m_adj"], v["or_all_adj"]
        ns = "not estimable" if n is None else f"{n['odds_ratio']:.2f} {[round(x, 2) for x in n['ci95_odds_ratio']]}"
        alls = "not estimable" if a is None else f"{a['odds_ratio']:.2f} {[round(x, 2) for x in a['ci95_odds_ratio']]}"
        print(f"{g:16s} cells {v['n_cells']:7d} near {v['n_cells_le500m']:6d} (flooded {v['n_flooded_le500m']})  "
              f"OR<=500m {ns} (Holm p {v.get('p_le500m_holm')})  OR all {alls}")
    s = res["settler_excluded"]
    print(f"settler excluded: dropped {s['n_aoi_dropped']} AOIs, R {s['R']} 90% {s['ci90_R']}")


if __name__ == "__main__":
    main()
