# %% v54（事後・amendment B-6 追記 1 の 9〜15）: 基準の参照点・残差の空間相関・読みの記述・AOI の採否・p 値
"""確証判定は変えない。B-6 追記 1 で計算の前に固定した記述を出す（8 は v52 で plasmode と一緒に）。

- 9: 基準の参照点（主分割、直した実行の held-out 予測）: AOI ごとに観測の浸水率と同じ割合のセルを浸水とする地図の
  CSI（HAND の低い順・M1・M3n の予測の高い順）と、定数の予測（fold ごとの訓練側の浸水率）の Brier
- 10: held-out の M1 の残差（AOI 内で中心化）の、AOI 内のセルの組の相関を距離帯別に
- 11: 読んでいない地域の言語（Wales・スペインの共同公用語の州・Brittany・Corsica、IE・NZ・AU）のセルの割合
- 12: 川沿いの構文（sur・upon・am・an der/dem/den・nad・nod・aan・on the）を名前に含むセルと集落の割合
- 13: 言語モデルの読みの区分の割合（unknown・解析の失敗を含む）と、大陸別の水系側の読みの割合（言語モデルと辞書）
- 14: AOI の採否の流れ（基準別の数、activation の数、除いた AOI の特徴）
- 15: プロトコルが求めた p 値（H2・H4 の TOST、H3 の 2 つの予測、H4 の改善の片側 p、H3・H4 の Holm 補正）

出力: stdout と docs/verification/v54_b6_descriptives.json
"""

from __future__ import annotations

import io
import itertools
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa
import world_eval as ev

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.dataset import continent_map
from topsem.world.lexicon import _SPLIT, WET_SIDE, Lexicon
from topsem.world.nread import NREAD_CATEGORIES

OUT = ROOT / "docs" / "verification" / "v54_b6_descriptives.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
V41T = DATA_DIR / "verification" / "v41_table.parquet"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
ALL = DATA_DIR / "global" / "geonames" / "allCountries.zip"
NREAD = DATA_DIR / "global" / "nread" / "nread.parquet"
SEED = 20261028
CELL_M = 250.0
MAX_CELLS_AOI = 1000
BINS_KM = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
NREAD_WETSIDE = ("water", "wetland", "riverine", "crossing", "lowland")
# 11: 読んでいない地域の言語（admin1 の符号は下の有名な町で確かめる）
REGIONS = {"GB": {"WLS": "Wales"},
           "ES": {"56": "Catalonia", "59": "Basque Country", "58": "Galicia", "60": "Valencia",
                  "07": "Balearic Islands", "32": "Navarre"},
           "FR": {"53": "Brittany", "94": "Corsica"}}
CHECK_TOWNS = {2653822: ("GB", "WLS"), 3128760: ("ES", "56"), 3128026: ("ES", "59"), 3109642: ("ES", "58"),
               2509954: ("ES", "60"), 2512989: ("ES", "07"), 3114472: ("ES", "32"), 2983990: ("FR", "53"),
               3038334: ("FR", "94")}
WHOLE_COUNTRIES = ("IE", "NZ", "AU")
RIVERSIDE_TOKENS = {"sur", "upon", "am", "nad", "nod", "aan"}
RIVERSIDE_PAIRS = {("an", "der"), ("an", "dem"), ("an", "den"), ("on", "the")}


def lexfix() -> dict[str, np.ndarray]:
    files = sorted((RUNS_DIR / "world_v2_lexfix" / "predictions" / "primary").glob("fold[0-9].npz"),
                   key=lambda f: int(f.stem[4:]))
    parts = [np.load(f) for f in files]
    out = {k: np.concatenate([q[k] for q in parts]) for k in ("y", "aoi", "unit", "p_M0", "p_M1", "p_M3n")}
    out["fold"] = np.concatenate([np.full(len(q["y"]), k) for k, q in enumerate(parts)])
    return out


# %% 9. 基準の参照点
def csi_top_share(score_low_is_wet: np.ndarray, y: np.ndarray, aoi: np.ndarray) -> dict:
    """AOI ごとに観測の浸水率と同じ数のセルを score の小さい順に浸水とする地図の CSI（全体と AOI の中央値）。"""
    flag = np.zeros(len(y), dtype=bool)
    per = []
    for idx in pd.Series(np.arange(len(y))).groupby(aoi).indices.values():
        k = int(y[idx].sum())
        if k == 0:
            continue
        order = idx[np.argsort(score_low_is_wet[idx], kind="stable")]
        flag[order[:k]] = True
        tp = int((flag[idx] & (y[idx] == 1)).sum())
        per.append(tp / (2 * k - tp))
    tp = int((flag & (y == 1)).sum())
    fp, fn = int((flag & (y == 0)).sum()), int((~flag & (y == 1)).sum())
    return {"csi": tp / (tp + fp + fn), "median_csi_aoi": float(np.median(per)), "n_aoi": len(per)}


def reference_points(d: dict, hand: np.ndarray) -> dict:
    y, aoi = d["y"].astype(np.int64), d["aoi"]
    out = {"hand_threshold_map": csi_top_share(hand, y, aoi),
           "M1_map": csi_top_share(-d["p_M1"].astype(np.float64), y, aoi),
           "M3n_map": csi_top_share(-d["p_M3n"].astype(np.float64), y, aoi)}
    yf = y.astype(np.float64)
    tr = np.array([yf[d["fold"] != k].mean() for k in range(int(d["fold"].max()) + 1)])
    out["brier"] = {"const_train_rate": float(np.mean((tr[d["fold"]] - yf) ** 2)),
                    "M0": float(np.mean((d["p_M0"] - yf) ** 2)), "M1": float(np.mean((d["p_M1"] - yf) ** 2))}
    return out


# %% 10. 残差の空間相関
def residual_correlogram(d: dict, ix: np.ndarray, iy: np.ndarray) -> dict:
    r = d["y"].astype(np.float64) - d["p_M1"].astype(np.float64)
    rng = np.random.default_rng(SEED)
    edges = np.array(BINS_KM) * 1000.0
    sxy = np.zeros(len(edges) - 1)
    sxx = np.zeros(len(edges) - 1)
    syy = np.zeros(len(edges) - 1)
    npairs = np.zeros(len(edges) - 1, dtype=np.int64)
    for idx in pd.Series(np.arange(len(r))).groupby(d["aoi"]).indices.values():
        if len(idx) > MAX_CELLS_AOI:
            idx = np.sort(rng.choice(idx, size=MAX_CELLS_AOI, replace=False))
        ri = r[idx] - r[idx].mean()  # AOI 内で中心化
        x, yy = ix[idx].astype(np.float64) * CELL_M, iy[idx].astype(np.float64) * CELL_M
        dist = np.hypot(x[:, None] - x[None, :], yy[:, None] - yy[None, :])
        iu = np.triu_indices(len(idx), k=1)
        dd, a, b = dist[iu], ri[iu[0]], ri[iu[1]]
        bi = np.digitize(dd, edges) - 1
        ok = (bi >= 0) & (bi < len(edges) - 1)
        np.add.at(sxy, bi[ok], (a * b)[ok])
        np.add.at(sxx, bi[ok], (a * a)[ok])
        np.add.at(syy, bi[ok], (b * b)[ok])
        np.add.at(npairs, bi[ok], 1)
    corr = sxy / np.sqrt(sxx * syy)
    return {"bins_km": [f"{lo:g}-{hi:g}" for lo, hi in itertools.pairwise(BINS_KM)],
            "corr": [float(c) for c in corr], "n_pairs": [int(n) for n in npairs],
            "max_cells_per_aoi": MAX_CELLS_AOI}


# %% 11. 読んでいない地域の言語
def admin1_for(gids: set[int]) -> dict[int, tuple[str, str]]:
    """GeoNames の全件（allCountries.zip）から、gid → (cc, admin1)。"""
    want = gids | set(CHECK_TOWNS)
    out = {}
    with zipfile.ZipFile(ALL) as z, z.open("allCountries.txt") as f:
        for line in io.TextIOWrapper(f, encoding="utf-8"):
            gid = int(line[:line.index("\t")])
            if gid in want:
                p = line.rstrip("\n").split("\t")
                out[gid] = (p[8], p[10])
    return out


def regional_languages(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    gids = {int(g) for g in t["P_gid1"].unique() if g >= 0}
    a1 = admin1_for(gids)
    check = {g: a1.get(g) == want for g, want in CHECK_TOWNS.items()}
    if not all(check.values()):
        msg = f"admin1 codes do not match the check towns: {check}"
        raise SystemExit(msg)
    cc = t["P_gid1"].map(lambda g: a1.get(int(g), ("", ""))[0])
    ad = t["P_gid1"].map(lambda g: a1.get(int(g), ("", ""))[1])
    n = len(t)
    out: dict = {"n_cells": n, "regions": {}, "countries": {}}
    for c, regs in REGIONS.items():
        for code, name in regs.items():
            m = (cc == c) & (ad == code)
            out["regions"][f"{c}-{name}"] = {"n_cells": int(m.sum()), "share": float(m.mean()),
                                             "n_settlements": int(t.loc[m, "P_gid1"].nunique())}
    for c in WHOLE_COUNTRIES:
        m = cc == c
        out["countries"][c] = {"n_cells": int(m.sum()), "share": float(m.mean()),
                               "n_settlements": int(t.loc[m, "P_gid1"].nunique())}
    tot = sum(v["n_cells"] for v in out["regions"].values()) + sum(v["n_cells"] for v in out["countries"].values())
    out["total_share"] = tot / n
    return out


# %% 12. 川沿いの構文
def riverside(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    names = gn.set_index("gid")["asciiname"]

    def has(name: str) -> bool:
        toks = [x for x in _SPLIT.split(str(name).lower()) if x]
        return bool(RIVERSIDE_TOKENS & set(toks)) or any(p in RIVERSIDE_PAIRS for p in itertools.pairwise(toks))

    g = pd.Series({int(x): has(names.get(int(x), "")) for x in t["P_gid1"].unique() if x >= 0})
    cell = t["P_gid1"].map(g).fillna(value=False).to_numpy(dtype=bool)
    ex = sorted(names.get(int(x), "") for x in g[g].index)[:15]
    return {"share_cells": float(cell.mean()), "n_settlements": int(g.sum()), "share_settlements": float(g.mean()),
            "examples": ex}


# %% 13. 言語モデルの読み
def nread_desc(t: pd.DataFrame, gn: pd.DataFrame) -> dict:
    nr = pd.read_parquet(NREAD).drop_duplicates("gid").set_index("gid")
    lex = Lexicon(LEX)
    cont = continent_map(DATA_DIR / "global" / "geonames" / "countryInfo.txt")
    g = gn.set_index("gid")
    gids = np.array([int(x) for x in t["P_gid1"].unique() if x >= 0])
    s = pd.DataFrame({"gid": gids})
    s["cat"] = [str(nr.at[x, "category"]) if x in nr.index else "missing" for x in gids]
    s["parsed"] = [bool(nr.at[x, "parsed_ok"]) if x in nr.index else False for x in gids]
    s["nr_wet"] = s["parsed"] & s["cat"].isin(NREAD_WETSIDE)
    s["lex_wet"] = [any(lex.match(g.at[x, "asciiname"], g.at[x, "cc"])[c] for c in WET_SIDE) for x in gids]
    s["continent"] = [cont.get(str(g.at[x, "cc"]), "") for x in gids]
    cells = t[["P_gid1"]].merge(s, left_on="P_gid1", right_on="gid", how="inner")
    out = {"n_settlements": len(s), "share_parse_failed_settlements": float((~s["parsed"]).mean()),
           "category_share_settlements": {c: float((s["parsed"] & (s["cat"] == c)).mean()) for c in NREAD_CATEGORIES},
           "category_share_cells": {c: float((cells["parsed"] & (cells["cat"] == c)).mean())
                                    for c in NREAD_CATEGORIES},
           "confidence_use": "feature nr_conf_n1 in the second stage; not used to threshold or filter readings",
           "wetside_by_continent_cells": {}}
    for k, d in cells.groupby("continent"):
        out["wetside_by_continent_cells"][k] = {"n_cells": len(d), "llm": float(d["nr_wet"].mean()),
                                                "lexicon": float(d["lex_wet"].mean()),
                                                "llm_unknown_or_failed": float(((d["cat"] == "unknown")
                                                                                | ~d["parsed"]).mean())}
    out["wetside_all_cells"] = {"llm": float(cells["nr_wet"].mean()), "lexicon": float(cells["lex_wet"].mean())}
    return out


# %% 14. AOI の採否の流れ
def inclusion_flow() -> dict:
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    dirs, reasons, infos = ev.select_aois(gn, ev.Sensitivity())
    cont = continent_map(ev.ROOT / "geonames" / "countryInfo.txt")
    iso = {}  # activation の国は名前で書かれている → ISO2（countryInfo の Country 列）
    for line in (ev.ROOT / "geonames" / "countryInfo.txt").read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            iso[f[4]] = f[0]
    keep = {d.name for d in dirs}
    rows = []
    for tag, info in infos.items():
        kind = ("riverine" if info.riverine_share >= 0.5 else
                "unspecified-as-riverine" if info.unspecified_as_riverine else
                "unspecified" if info.unspecified_share >= 0.5 else "other")
        rows.append({"tag": tag, "included": tag in keep, "activation": info.activation, "type": kind,
                     "n_eval": info.n_eval, "wet_frac": info.n_wet / info.n_eval if info.n_eval else None,
                     "continent": cont.get(iso.get(info.countries[0], "") or _cc_inside(gn, info.bounds), "?")
                     if info.countries else ""})
    d = pd.DataFrame(rows)
    inc, exc = d[d["included"]], d[~d["included"]]

    def profile(x: pd.DataFrame) -> dict:
        return {"n_aoi": len(x), "n_activations": int(x["activation"].nunique()),
                "type": x["type"].value_counts().to_dict(), "continent": x["continent"].value_counts().to_dict(),
                "median_n_eval": float(x["n_eval"].median()),
                "median_wet_frac": float(x["wet_frac"].dropna().median()) if x["wet_frac"].notna().any() else None}

    return {"n_aoi_with_labels": len(d), "reasons": reasons, "included": profile(inc), "excluded": profile(exc),
            "excluded_by_reason": {r: profile(exc[exc["tag"].isin([t for t in exc["tag"]
                                                                    if _why(infos[t], gn) == r])])
                                   for r in reasons if r != "ok"},
            "activations_only_in_excluded": len(set(exc["activation"]) - set(inc["activation"]))}


def _cc_inside(gn: pd.DataFrame, b: tuple[float, float, float, float]) -> str:
    """名前が countryInfo と一致しない国（例 Czech Republic）は、AOI の範囲内の集落の最も多い国にする。"""
    w, s_, e, n = b
    p = gn[(gn["fclass"] == "P") & gn["lon"].between(w, e) & gn["lat"].between(s_, n)]
    return str(p["cc"].mode().iloc[0]) if len(p) else ""


def _why(info, gn: pd.DataFrame) -> str:
    return ev.include(info, n_names_inside=ev._p_names_inside(gn, info.bounds))[1]


# %% 15. プロトコルが求めた p 値
def _boot(run: str, split: str) -> tuple:
    data = wa.load_split(RUNS_DIR / run, split)
    y = data["y"].astype(np.float64)
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    s, n, arms = wa.unit_sums(y, preds, data["unit"])
    idx = {a: j for j, a in enumerate(arms)}
    hat = s.sum(axis=0) / n.sum()
    return wa.make_delta(hat, wa.bootstrap_bs(s, n), idx), idx


def _tost_p(d: float, b: np.ndarray, margin: float) -> dict:
    p_up, p_lo = float(np.mean(b >= margin)), float(np.mean(b <= -margin))
    return {"delta_bs": d, "margin": margin, "p_upper": p_up, "p_lower": p_lo, "p_tost": max(p_up, p_lo),
            "resolution": 1 / len(b)}


def spearman_exact_p(rho: float, n: int = 5) -> float:
    """順位の全ての並べ替えでの P(ρ ≤ ρ_obs)（片側、ρ < 0 の予測）。"""
    k = np.arange(n)
    vals = [1 - 6 * float(((k - np.array(p)) ** 2).sum()) / (n * (n * n - 1)) for p in itertools.permutations(k)]
    return float(np.mean(np.array(vals) <= rho + 1e-12))


def p_values() -> dict:
    out = {}
    for run in ("world_v2", "world_v2_lexfix"):
        blk: dict = {}
        delta, _ = _boot(run, "primary")
        d01, _ = delta("M0", "M1")
        d, b = delta("M1c", "M3n")
        blk["H2_tost"] = _tost_p(d, b, wa.SESOI * d01)
        lad = json.loads((RUNS_DIR / run / "metrics_primary.json").read_text(encoding="utf-8")).get("ladder")
        if lad:
            blk["H3"] = {"p_delta0_gt_delta4": lad["p_one_sided_delta0_gt_delta4"],
                         "spearman": lad["spearman_k_delta"],
                         "p_spearman_exact": spearman_exact_p(lad["spearman_k_delta"])}
        if (RUNS_DIR / run / "predictions" / "llmcut").is_dir():
            delta4, idx4 = _boot(run, "llmcut")
            if "M1rc" in idx4 and "M3rn" in idx4:
                d01b, _ = delta4("M0", "M1")
                d4, b4 = delta4("M1rc", "M3rn")
                blk["H4_tost"] = _tost_p(d4, b4, wa.SESOI * d01b)
                blk["H4_improvement_p"] = float(np.mean(b4 <= 0))
        if "H3" in blk and "H4_improvement_p" in blk:
            ps = {"H3": blk["H3"]["p_delta0_gt_delta4"], "H4": blk["H4_improvement_p"]}
            order = sorted(ps, key=ps.get)
            adj, run_max = {}, 0.0
            for i, h in enumerate(order):
                run_max = max(run_max, min(1.0, (len(ps) - i) * ps[h]))
                adj[h] = run_max
            blk["holm"] = adj
        out[run] = blk
    return out


def main() -> None:
    d = lexfix()
    t = pd.read_parquet(V37T, columns=["aoi", "y", "P_gid1", "hand_mean"])
    t41 = pd.read_parquet(V41T, columns=["aoi", "y", "ix", "iy"])
    if not (np.array_equal(t["aoi"].to_numpy(), d["aoi"].astype(str)) and np.array_equal(t["y"].to_numpy(), d["y"])
            and np.array_equal(t41["aoi"].to_numpy(), t["aoi"].to_numpy())):
        msg = "tables and lexfix predictions are not aligned"
        raise SystemExit(msg)
    gn = pd.read_parquet(GN, columns=["gid", "asciiname", "cc"]).drop_duplicates("gid")
    res = {"rule": "amendment B-6 addendum 1, items 9-15", "seed": SEED,
           "reference_points": reference_points(d, t["hand_mean"].to_numpy(np.float64)),
           "residual_correlogram": residual_correlogram(d, t41["ix"].to_numpy(), t41["iy"].to_numpy()),
           "regional_languages": regional_languages(t, gn), "riverside": riverside(t, gn),
           "nread": nread_desc(t, gn), "inclusion_flow": inclusion_flow(), "p_values": p_values()}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: res[k] for k in ("reference_points", "residual_correlogram")}, indent=1)[:2500])
    print(json.dumps(res["regional_languages"], indent=1)[:1500])
    print(json.dumps(res["riverside"]))
    print(json.dumps(res["nread"]["wetside_all_cells"]), json.dumps(res["nread"]["category_share_cells"]))
    print(json.dumps({k: v for k, v in res["inclusion_flow"].items() if k != "excluded_by_reason"}, indent=1)[:2000])
    print(json.dumps(res["p_values"], indent=1))


if __name__ == "__main__":
    main()
