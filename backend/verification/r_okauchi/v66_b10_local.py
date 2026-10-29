# %% v66（事後・amendment B-10 の B1〜B15）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-10 で計算の前に固定した量を、既存の表と held-out の予測から計算する。

- B1: 「島・水辺の草地」の群（B-5 の 8 の island）の集落の場所（177 m 以内のセル）の HAND の AOI 内の差と、
  au・aue・auen を除いた群の全セル・場所の差と 500 m 以内の OR
- B2: 氾濫原の読み（河川・湿地）の OR を、国を一つずつ除いて・EMSR097 を除いて・形を除いて読み直して。形ごとの構成
- B3: 主の 500 m の推定の pairs bootstrap の SE から、両側 5%・検出力 80% で検出できる最小の OR
- B4: 主の 500 m の推定の場所グループの jackknife（CV3、t 分布）と wild cluster score bootstrap（Rademacher、帰無を
  課す、検定の反転）の 95% 区間
- B5: 誤分類の補正の格子（p1 × p0 の 16 点、OR* は点と pairs の両端）
- B6: 二人の符号化者の精度・再現率を、辞書の形ごと・国ごと（符号のある名が 10 以上）
- B7: EMS の画像（Sentinel-1・他のレーダー・光学・記録無し）と製品の世代ごとの R と 500 m 以内の OR
- B8: アフリカを除いた R
- B9: 活性化をまとめた分割（world_v2_b7actfix）の物差し 5 通り × bootstrap 2 通りの R と、分割ごとの Brier の表
- B10: GFM で当て直した実行（world_v2_gfm）の AOI の中だけの評価（B-9 の B4 と同じ）
- B11: 基礎を固定した注入（形 A・B、δ = ln 1.25・ln 1.5・ln 2・ln 3、各 20 回）のラベルでの第 5 節の推定
- B12: HAND だけの参照の二つの定義と、AOI の中央値・全てまとめた値の CSI・的中率・空振り率
- B13: ドイツ語・英語・フランス語の塊で読む国の AOI での、浸水率と名の構成の関係（activation の固定効果）
- B14: 他の言語で読んだセル（国の水準・行政区の水準）
- B15: N-read（Gemma）の区分の分布

出力: stdout と docs/verification/v66_b10_local.json（部分ごとに足す。`--only B1,B4` で一部だけ）。
pairs bootstrap の係数は data/verification/v66_pairs_<key>.npy に残して使い回す。
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR
from topsem.world.b6 import cell_keys, hashed_uniform, inject_labels, plasmode_salt
from topsem.world.b8 import correct_misclassification
from topsem.world.dataset import continent_map
from topsem.world.lexicon import _SPLIT, WET_SIDE, Lexicon, _hit

OUT = ROOT / "docs" / "verification" / "v66_b10_local.json"
VER = DATA_DIR / "verification"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
COUNTRY_INFO = DATA_DIR / "global" / "geonames" / "countryInfo.txt"
LABELS = DATA_DIR / "global" / "labels"
NE_ADMIN1 = DATA_DIR / "global" / "naturalearth" / "ne_10m_admin_1_states_provinces.zip"
NREAD = DATA_DIR / "global" / "nread" / "nread.parquet"
SEED = 20261030
SESOI_R = 0.05  # 余白（プロトコル §6、R の単位）
B_PAIRS = 2000
B_WILD = 2000
B_BOOT = 4000
NEAR_M = 500.0
TERM = "le500m"
N_WORKERS = int(os.environ.get("V66_WORKERS", "8"))
Z80 = 2.8  # 両側 5%・検出力 80%: 1.96 + 0.84
NOAU = {"au", "aue", "auen"}
# B2: 読み直しで除く形（河川の区分の、英語の bank、ドイツ語の au、堤と干拓地）
EXCL_FORMS = {"bank", "banks", "au", "aue", "auen", "levee", "dijk", "dam", "digue", "dige", "hraz", "argine",
              "argini", "dig", "nasip", "dambis", "pylimas", "tanggul", "penger", "polder", "bonifica"}
FLOODPLAIN = ("riverine", "wet")
P1_GRID = (0.6, 0.7, 0.8, 0.9)
P0_GRID = (0.05, 0.10, 0.20, 0.30)
MIN_CODED = 10
INJ_DELTAS = (math.log(1.25), math.log(1.5), math.log(2.0), math.log(3.0))
INJ_REPS = 20
# B14: 言語群の符号に含める言語（ISO 639）
LANG_GROUPS = {"sh": {"sh", "sr", "hr", "bs", "cnr"},
               "bantu": {"ny", "bem", "sn", "nd", "zu", "xh", "rw", "rn", "lg", "ln", "kg", "ts", "tn", "st", "ss",
                         "ve", "nr", "ki", "yao", "tum", "loz", "ng", "hz", "sw", "mg"}}
# B14: 地域の言語の行政区（Natural Earth の iso_a2 と、name・region・geonunit のどれかの値）
REGIONAL = {
    "Wales": ("GB", "geonunit", {"Wales"}),
    "Highland and Eilean Siar": ("GB", "name", {"Highland", "Eilean Siar"}),
    "Northern Ireland": ("GB", "geonunit", {"Northern Ireland"}),
    "Ireland": ("IE", None, None),
    "Catalan-speaking Spain": ("ES", "region", {"Cataluña", "Valenciana", "Islas Baleares"}),
    "Basque Country and Navarre": ("ES", "region", {"País Vasco", "Foral de Navarra"}),
    "Galicia": ("ES", "region", {"Galicia"}),
    "Brittany": ("FR", "region", {"Bretagne"}),
    "Alsace": ("FR", "name", {"Bas-Rhin", "Haut-Rhin"}),
    "Friesland": ("NL", "name", {"Friesland", "Fryslân"}),
    "Finnmark and Troms": ("NO", "name", {"Finnmark", "Troms", "Troms og Finnmark"}),
    "Norrbotten": ("SE", "name", {"Norrbotten"}),
    "Lapland (Finland)": ("FI", "name", {"Lapland", "Lappi"}),
}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v37 = _load("v37_positive_control")
v39 = _load("v39_within_aoi_odds")
v42 = _load("v42_blind_review_analyses")
v55 = _load("v55_b7_local")
v63 = _load("v63_b9_local")
v36 = _load("v36_yardstick_robustness")
v44 = _load("v44_wikipedia_etymology")
# 後から読んだスクリプトが同じ名のモジュールを読み直すので、子のプロセスへ渡す関数は登録済みのモジュールから取る
v37, v39, v55 = (sys.modules[n] for n in ("v37_positive_control", "v39_within_aoi_odds", "v55_b7_local"))
v37.SEED, v37.B = SEED, B_BOOT


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


# %% 共通: 名の読み（集落ごと）と pairs bootstrap
_LEX: Lexicon | None = None


def lexicon() -> Lexicon:
    global _LEX
    if _LEX is None:
        _LEX = Lexicon(LEX)
    return _LEX


def settlement_forms(gids: np.ndarray, cats: tuple[str, ...]) -> dict[int, set[tuple[str, str]]]:
    """集落ごとに、照合した (区分, 形) の集合（最寄りの名、n1 と同じ照合: 国の言語の項目だけ）。"""
    lex = lexicon()
    gn = pd.read_parquet(GN, columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    out: dict[int, set[tuple[str, str]]] = {}
    for gid in np.unique(gids):
        if gid < 0 or gid not in gn.index:
            continue
        name, cc = gn.at[gid, "asciiname"], gn.at[gid, "cc"]
        toks = [x for x in _SPLIT.split(str(name).lower()) if x]
        hits = {(e.cat, e.form) for e in lex.entries_for(str(cc)) if e.cat in cats
                for tok in toks if _hit(tok, e)}
        if hits:
            out[int(gid)] = hits
    return out


def flag_from_forms(gids: np.ndarray, forms: dict[int, set[tuple[str, str]]], keep) -> np.ndarray:
    ok = {g for g, h in forms.items() if any(keep(c, f) for c, f in h)}
    return np.isin(gids, list(ok)).astype(np.float64)


def _design_keep(t: pd.DataFrame, cls: str, *, banded: bool = True) -> tuple:
    y = t["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(t["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = t.loc[keep].reset_index(drop=True)
    z, names = v39.design(s, (cls,), True, banded=banded)
    return s, y[keep], z, names


def pairs_or(t: pd.DataFrame, cls: str, key: str, mask: np.ndarray | None = None) -> dict:
    """500 m 以内の区分 × 水系側のオッズ比（v55.odds と同じ式）と、場所の pairs cluster bootstrap（種 SEED）。"""
    s0 = t if mask is None else t.loc[mask].reset_index(drop=True)
    o = v55.odds(s0, cls=cls)
    cache = VER / f"v66_pairs_{key}.npy"
    if cache.is_file():
        b = np.load(cache)
    else:
        s, y, z, names = _design_keep(s0, cls)
        j = names.index(f"{cls}:{TERM}")
        aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
        unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
        g, na = len(unit_lab), int(aoi.max()) + 1
        rows = [np.flatnonzero(unit == u) for u in range(g)]
        rng = np.random.default_rng(SEED)
        sels = [rng.integers(0, g, size=g) for _ in range(B_PAIRS)]
        _threads1()
        t0 = time.time()
        with ProcessPoolExecutor(N_WORKERS, initializer=v55._pairs_init,
                                 initargs=(y, z, aoi, rows, na)) as ex:
            res = list(ex.map(v55._pairs_rep, sels, chunksize=8))
        b = np.array([r[j] if r is not None else np.nan for r in res])
        np.save(cache, b)
        print(f"  pairs {key}: {time.time() - t0:.0f} s", flush=True)
    ok = np.isfinite(b)
    o["pairs_ci95"] = [float(np.exp(np.quantile(b[ok], 0.025))), float(np.exp(np.quantile(b[ok], 0.975)))]
    o["pairs_log_or_sd"] = float(np.std(b[ok], ddof=1))
    o["pairs_reps_ok"] = int(ok.sum())
    return o


# %% B1: 島・水辺の草地
def b1_islands() -> dict:
    t = v42.table()
    gids = t["P_gid1"].to_numpy()
    forms = settlement_forms(gids, ("riverine",))
    isl = flag_from_forms(gids, forms, lambda c, f: c == "riverine" and f in v42.ISLAND) > 0
    noau = flag_from_forms(gids, forms, lambda c, f: c == "riverine" and f in v42.ISLAND - NOAU) > 0
    au_only = isl & ~noau
    allm = np.ones(len(t), dtype=bool)
    near = t["P_d1"].to_numpy() <= NEAR_M
    s = v37.settlement_table(t, {"island": isl.astype(np.float64), "island_noau": noau.astype(np.float64)})
    out: dict = {"n_cells": {"island": int(isl.sum()), "island_noau": int(noau.sum()), "au_only": int(au_only.sum())},
                 "n_settlements": {"island": len(np.unique(gids[isl])),
                                   "island_noau": len(np.unique(gids[noau]))}}
    for g, f in (("island", isl), ("island_noau", noau)):
        out[g] = {"site_hand_mean": v37.settlement_contrast(s, "hand_mean", g),
                  "site_hand_pct": v37.settlement_contrast(s, "hand_pct", g),
                  "all_hand_mean": v37.contrast(t, "hand_mean", f, allm, "capped"),
                  "le500m_hand_mean": v37.contrast(t, "hand_mean", f, near, "capped"),
                  "or_le500m": _group_or(t, f)}
    out["rule"] = "B-5 item 8 estimators; site = cells within 177 m of the settlement point (v37)"
    return out


def _group_or(t: pd.DataFrame, f: np.ndarray) -> dict | None:
    """v42 の群のオッズ比（距離帯 + 群 × 距離帯 + 凍結した M1 の logit、AOI の固定効果、場所でクラスタ頑健）。"""
    y = t["y"].to_numpy(dtype=np.float64)
    rate = t.assign(_y=y).groupby("aoi")["_y"].transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s, ys, gi = t.loc[keep].reset_index(drop=True), y[keep], f[keep]
    d1 = s["P_d1"].to_numpy()
    m = gi & (d1 <= NEAR_M)
    if m.sum() == 0 or ys[m].min() == ys[m].max():
        return None
    pm = np.clip(s["p_M1"].to_numpy(dtype=np.float64), 1e-5, 1 - 1e-5)
    bands = [((d1 > NEAR_M) & (d1 <= 1000)).astype(float), (d1 > 1000).astype(float)]
    z = np.column_stack([*bands, gi & (d1 <= NEAR_M), gi & (d1 > NEAR_M) & (d1 <= 1000), gi & (d1 > 1000),
                         np.log(pm) - np.log(1 - pm)]).astype(np.float64)
    res = v39.fit(ys, z, s["aoi"].to_numpy(), s["unit"].to_numpy())
    tab = v39.table(res, ["band_500m_1km", "band_gt1km", "g:le500m", "g:500m_1km", "g:gt1km", "logit_M1"])
    # 浸水した群のセル（500 m 以内）がどれだけの場所グループ・AOI・国にあるか（記述。B1 の追記、本文の注のため）
    fl = m & (ys == 1)
    where = {"n_units_flooded_le500m": int(s.loc[fl, "unit"].nunique()),
             "n_aoi_flooded_le500m": int(s.loc[fl, "aoi"].nunique())}
    if "cc" in s.columns:
        where["flooded_le500m_by_country"] = {str(k): int(v) for k, v in s.loc[fl, "cc"].value_counts().items()}
    return tab["g:le500m"] | {"n_cells_le500m": int(m.sum()), "n_flooded_le500m": int(ys[m].sum())} | where


# %% B2: 氾濫原の読み
def b2_floodplain(t: pd.DataFrame) -> dict:
    gids = t["P_gid1"].to_numpy()
    forms = settlement_forms(gids, FLOODPLAIN)
    t["lex_fpall_n1"] = flag_from_forms(gids, forms, lambda c, f: True)
    t["lex_fpx_n1"] = flag_from_forms(gids, forms, lambda c, f: not (c == "riverine" and f in EXCL_FORMS))
    fp = t["lex_floodplain_n1"].to_numpy() > 0
    out: dict = {"same_code_path_mismatch_share": float(((t["lex_fpall_n1"].to_numpy() > 0) != fp).mean()),
                 "main": pairs_or(t, "floodplain", "fp_main"),
                 "reread_all_forms": pairs_or(t, "fpall", "fp_all"),
                 "reread_without_forms": pairs_or(t, "fpx", "fp_excl") | {
                     "share_group_cells_dropped": float(((t["lex_fpx_n1"] == 0) & fp).sum() / max(fp.sum(), 1)),
                     "excluded_forms": sorted(EXCL_FORMS)}}
    y = t["y"].to_numpy()
    near = (t["P_d1"].to_numpy() <= NEAR_M) & fp
    cc = t["cc"].to_numpy(dtype=str)
    ctry = sorted(set(cc[near & (y == 1)]))
    out["leave_one_country_out"] = {}
    for c in ctry:
        out["leave_one_country_out"][c] = pairs_or(t, "floodplain", f"fp_no_{c}", cc != c) | {
            "n_flooded_near_in_country": int((near & (y == 1) & (cc == c)).sum())}
    act = t["aoi"].str.split("_").str[0].to_numpy()
    out["without_EMSR097"] = pairs_or(t, "floodplain", "fp_no_EMSR097", act != "EMSR097") | {
        "n_flooded_near_in_EMSR097": int((near & (y == 1) & (act == "EMSR097")).sum())}
    rows = []
    for (cat, form), idx in _form_index(gids, forms).items():
        m = near & np.isin(gids, idx)
        if m.any():
            rows.append({"category": cat, "form": form, "n_cells_le500m": int(m.sum()),
                         "n_flooded_le500m": int(y[m].sum()), "n_settlements": len(np.unique(gids[m])),
                         "excluded_in_reread": bool(cat == "riverine" and form in EXCL_FORMS)})
    out["composition_by_form"] = sorted(rows, key=lambda r: -r["n_flooded_le500m"])
    return out


def _form_index(gids: np.ndarray, forms: dict[int, set[tuple[str, str]]]) -> dict[tuple[str, str], list[int]]:
    idx: dict[tuple[str, str], list[int]] = {}
    for g, h in forms.items():
        for k in h:
            idx.setdefault(k, []).append(g)
    return idx


# %% B3・B4: 主の 500 m の推定の検出力と区間
def b3_mde(main: dict) -> dict:
    se_p = main["pairs_log_or_sd"]
    b, se_s = math.log(main["odds_ratio"]), (math.log(main["ci95"][1]) - math.log(main["ci95"][0])) / (2 * 1.96)
    return {"pairs": {"se_log_or": se_p, "mde_or_up": math.exp(Z80 * se_p), "mde_or_down": math.exp(-Z80 * se_p)},
            "sandwich": {"se_log_or": se_s, "mde_or_up": math.exp(Z80 * se_s), "mde_or_down": math.exp(-Z80 * se_s)},
            "point_log_or": b, "rule": "MDE = exp(±2.80 × SE): two-sided 5 %, 80 % power"}


_J: dict = {}


def _jk_init(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, unit: np.ndarray, j: int) -> None:
    _J.update(y=y, z=z, aoi=aoi, unit=unit, j=j)


def _jk_rep(g: int) -> float:
    m = _J["unit"] != g
    try:
        return float(v39.fit(_J["y"][m], _J["z"][m], _J["aoi"][m], _J["unit"][m])["beta"][_J["j"]])
    except np.linalg.LinAlgError:
        return float("nan")


def fit_offset(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, offset: np.ndarray,
               iters: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """AOI の固定効果 + z の最尤推定（offset つき、v39.fit と同じ IRLS）。戻り値: 係数と当てはめた確率。"""
    ua, ia = np.unique(aoi, return_inverse=True)
    n, k, a = len(y), z.shape[1], len(ua)
    d = sparse.csr_matrix((np.ones(n), (np.arange(n), ia)), shape=(n, a))
    rate = np.clip(np.bincount(ia, weights=y) / np.bincount(ia), 1e-6, 1 - 1e-6)
    alpha = np.log(rate / (1 - rate))
    beta = np.zeros(k)
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta + offset)))
        w = p * (1 - p)
        r = y - p
        haa = np.bincount(ia, weights=w)
        haz = d.T @ (z * w[:, None])
        hzz = z.T @ (z * w[:, None])
        ga, gz = np.bincount(ia, weights=r), z.T @ r
        s = hzz - haz.T @ (haz / haa[:, None])
        db = np.linalg.solve(s, gz - haz.T @ (ga / haa))
        da = (ga - haz @ db) / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max()) < 1e-10:
            break
    return beta, 1 / (1 + np.exp(-(alpha[ia] + z @ beta + offset)))


def score_stat(y: np.ndarray, z: np.ndarray, aoi: np.ndarray, unit_inv: np.ndarray, j: int,
               b0: float) -> tuple[float, np.ndarray]:
    """帰無 β_j = b0 を課した当てはめでの、β_j の効率的スコアのクラスタ和 s_g と統計量 (Σ s_g)² / Σ s_g²。"""
    zr = np.delete(z, j, axis=1)
    zj = z[:, j]
    _, p = fit_offset(y, zr, aoi, b0 * zj)
    w = p * (1 - p)
    r = y - p
    ia = pd.factorize(aoi)[0]
    haa = np.bincount(ia, weights=w)

    def dm(v: np.ndarray) -> np.ndarray:  # AOI の中で重みつきに中心化
        return v - (np.bincount(ia, weights=w * v) / haa)[ia]

    zj_d = dm(zj)
    zr_d = np.column_stack([dm(zr[:, c]) for c in range(zr.shape[1])])
    gam = np.linalg.solve(zr_d.T @ (zr_d * w[:, None]), zr_d.T @ (w * zj_d))
    zc = zj_d - zr_d @ gam
    sg = np.bincount(unit_inv, weights=zc * r)
    return float(sg.sum() ** 2 / (sg ** 2).sum()), sg


def b4_intervals(t: pd.DataFrame) -> dict:
    s, y, z, names = _design_keep(t, "wetside")
    j = names.index(f"wetside:{TERM}")
    aoi = s["aoi"].to_numpy()
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    g = len(unit_lab)
    full = v39.fit(y, z, aoi, s["unit"].to_numpy())
    b_hat = float(full["beta"][j])
    _threads1()
    t0 = time.time()
    with ProcessPoolExecutor(N_WORKERS, initializer=_jk_init, initargs=(y, z, aoi, unit, j)) as ex:
        jk = np.array(list(ex.map(_jk_rep, range(g), chunksize=4)))
    ok = np.isfinite(jk)
    se_jk = float(np.sqrt((g - 1) / g * ((jk[ok] - jk[ok].mean()) ** 2).sum()))
    tq = float(stats.t.ppf(0.975, g - 1))
    print(f"  jackknife: {time.time() - t0:.0f} s", flush=True)
    rng = np.random.default_rng(SEED)
    v = rng.choice([-1.0, 1.0], size=(B_WILD, g))
    grid = np.round(np.arange(-0.80, 0.8001, 0.02), 4)
    pv = {}
    t0 = time.time()
    for b0 in grid:
        w0, sg = score_stat(y, z, aoi, unit, j, float(b0))
        wb = (v @ sg) ** 2 / (sg ** 2).sum()
        pv[float(b0)] = float(np.mean(wb >= w0))
    print(f"  wild score bootstrap: {time.time() - t0:.0f} s", flush=True)
    acc = [b for b, p in pv.items() if p >= 0.05]
    contiguous = bool(acc) and len(acc) == round((max(acc) - min(acc)) / 0.02) + 1
    return {"odds_ratio": math.exp(b_hat),
            "jackknife_cv3": {"se_log_or": se_jk, "t_quantile": tq, "n_groups": g, "n_failed": int((~ok).sum()),
                              "ci95": [math.exp(b_hat - tq * se_jk), math.exp(b_hat + tq * se_jk)]},
            "wild_score_bootstrap": {"ci95": [math.exp(min(acc)), math.exp(max(acc))] if acc else None,
                                     "p_at_zero": pv.get(0.0), "grid_step_log_or": 0.02,
                                     "grid": [float(grid[0]), float(grid[-1])], "accepted_contiguous": contiguous,
                                     "draws": B_WILD, "weights": "Rademacher",
                                     "at_grid_edge": bool(acc) and (min(acc) <= grid[0] or max(acc) >= grid[-1])}}


# %% B5: 誤分類の補正の格子
def b5_grid(main: dict) -> dict:
    rows = []
    for label, o in (("point", main["odds_ratio"]), ("pairs_lo", main["pairs_ci95"][0]),
                     ("pairs_hi", main["pairs_ci95"][1])):
        for p1 in P1_GRID:
            for p0 in P0_GRID:
                c = float(correct_misclassification(o, p1, p0))
                rows.append({"or_obs": label, "or_obs_value": o, "p1": p1, "p0": p0,
                             "or_corrected": None if not np.isfinite(c) else c})
    return {"rows": rows, "formula": "OR = (OR*(1 - p0) - (1 - p1)) / (p1 - OR* p0); undefined when not positive"}


# %% B6: 符号化者ごとの精度・再現率（形ごと・国ごと）
def b6_validation() -> dict:
    samp = _json(v44.SAMPLE)
    sizes = samp["stratum_sizes"]
    ns = samp["sample_sizes"]
    meta = pd.DataFrame(samp["records"]).set_index("gid")
    key = {json.loads(x)["id"]: json.loads(x)["key"] for x in v44.TO_CODE.read_text(encoding="utf-8").splitlines()}
    c1 = pd.read_csv(v44.CODES)
    c1["cats"] = c1["cats"].fillna("").map(lambda s: {c for c in s.split(";") if c})
    c2 = pd.read_csv(ROOT / "docs" / "verification" / "v57b_codes_coder2.csv")
    c2["cats"] = c2["categories"].fillna("").map(lambda s: {c for c in str(s).split(";") if c})
    lex = lexicon()
    out: dict = {"min_coded": MIN_CODED, "weights": "stratum weights (population / sample) for recall"}
    for name, c in (("coder1", c1), ("coder2", c2)):
        d = c[c["status"] == "coded"].copy()
        d["gid"] = d["id"].map(key)
        d = d.join(meta, on="gid")
        d["truth"] = d["cats"].map(lambda s: bool(s & set(WET_SIDE)))
        d["w"] = d["stratum"].map({k: sizes[k] / ns[k] for k in sizes})
        d["forms"] = [sorted({e.form for e in lex.entries_for(str(r.cc)) if e.cat in WET_SIDE
                              for tok in _SPLIT.split(str(r.name).lower()) if tok and _hit(tok, e)})
                      for r in d.itertuples()]
        by_form = {}
        for f in sorted({f for fs in d["forms"] for f in fs}):
            sub = d[d["forms"].map(lambda fs, f=f: f in fs) & d["lex_wetside"]]
            if len(sub) >= MIN_CODED:
                by_form[f] = {"n_flagged": len(sub), "precision": float(sub["truth"].mean())}
        by_country = {}
        for cc, sub in d.groupby("cc"):
            if len(sub) < MIN_CODED:
                continue
            fl = sub["lex_wetside"].astype(bool)
            tr = sub["truth"].astype(bool)
            by_country[cc] = {"n_coded": len(sub), "n_flagged": int(fl.sum()),
                              "precision": float(tr[fl].mean()) if fl.any() else None,
                              "recall_weighted": float((sub["w"] * (tr & fl)).sum() / (sub["w"] * tr).sum())
                              if tr.any() else None,
                              "n_truth": int(tr.sum())}
        out[name] = {"n_coded": len(d), "by_form": by_form, "by_country": by_country}
    return out


# %% 共通: R と bootstrap（場所グループ）
def r_subset(t: pd.DataFrame, mask: np.ndarray, rng: np.random.Generator, prefix: str = "p_lf_") -> dict:
    s = t.loc[mask]
    y = s["y"].to_numpy(np.float64)
    lab, inv = np.unique(s["unit"].to_numpy(), return_inverse=True)
    g = len(lab)
    sse = {m: np.bincount(inv, weights=(s[f"{prefix}{m}"].to_numpy(np.float64) - y) ** 2, minlength=g)
           for m in ("M0", "M1", "M1c", "M3n")}
    n_u = np.bincount(inv, minlength=g).astype(np.float64)
    cnt = rng.multinomial(g, np.full(g, 1 / g), size=B_BOOT).astype(np.float64)
    den = cnt @ sse["M0"] - cnt @ sse["M1"]
    num = cnt @ sse["M1c"] - cnt @ sse["M3n"]
    rb = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    d0 = sse["M0"].sum() - sse["M1"].sum()
    return {"R": float((sse["M1c"].sum() - sse["M3n"].sum()) / d0) if d0 > 0 else None,
            "ci95": np.nanquantile(rb, [0.025, 0.975]).tolist(), "ci90": np.nanquantile(rb, [0.05, 0.95]).tolist(),
            "share_boot_yardstick_nonpositive": float(np.mean(den <= 0)),
            "yardstick": float(d0 / n_u.sum()), "delta_bs": float((sse["M1c"].sum() - sse["M3n"].sum()) / n_u.sum()),
            "n_cells": len(s), "n_aoi": int(s["aoi"].nunique()), "n_units": g, "replicates": B_BOOT}


# %% B7: センサーと世代
def b7_sensor(t: pd.DataFrame) -> dict:
    summ = _json(DATA_DIR / "global" / "gfm_summary.json")["aois"]
    sens = t["aoi"].map(lambda a: summ.get(a, {}).get("sensor", "none")).to_numpy(dtype=str)
    gen = t["aoi"].map({a: _json(LABELS / a / "meta.json").get("schema_generation") for a in t["aoi"].unique()})
    gen = gen.fillna(-1).astype(int).to_numpy()
    rng = np.random.default_rng(SEED)
    out: dict = {"sensor": {}, "generation": {}}
    for k in ("s1", "sar", "optical", "none"):
        m = sens == k
        if m.sum() == 0:
            continue
        out["sensor"][k] = {"R": r_subset(t, m, rng), "or_le500m": _or_if(t, m, f"sens_{k}")}
    for k in sorted(set(gen.tolist())):
        m = gen == k
        out["generation"][str(k)] = {"R": r_subset(t, m, rng), "or_le500m": _or_if(t, m, f"gen_{k}")}
    out["rule"] = "sensor from the EMS source table (gfm_summary sensor: s1, sar, optical, none); generation from meta"
    return out


def _or_if(t: pd.DataFrame, m: np.ndarray, key: str) -> dict | None:
    near = m & (t["P_d1"].to_numpy() <= NEAR_M) & (t["lex_wetside_n1"].to_numpy() > 0) & (t["y"].to_numpy() == 1)
    if near.sum() < 5:
        return {"n_flooded_near": int(near.sum()), "note": "too few flooded water-side cells within 500 m"}
    try:
        return pairs_or(t, "wetside", key, m) | {"n_flooded_near": int(near.sum())}
    except (np.linalg.LinAlgError, ValueError) as e:
        return {"error": f"{type(e).__name__}: {e}"}


# %% B8: アフリカを除く
def b8_no_africa(t: pd.DataFrame) -> dict:
    cont = continent_map(COUNTRY_INFO)
    aoi_cc = t.groupby("aoi")["cc"].agg(lambda s: s[s != ""].mode().iloc[0] if (s != "").any() else "")
    af = t["aoi"].map(aoi_cc).map(lambda c: cont.get(c, "")) == "AF"
    rng = np.random.default_rng(SEED)
    return {"n_aoi_africa": int(t.loc[af, "aoi"].nunique()), "without_africa": r_subset(t, ~af.to_numpy(), rng),
            "all": r_subset(t, np.ones(len(t), dtype=bool), rng),
            "rule": "AOI continent from the majority country of its cells' nearest settlements"}


# %% B9: 活性化をまとめた分割
def b9_activation() -> dict:
    lex = lexicon()
    cont = continent_map(COUNTRY_INFO)
    res = v36.analyse("primary", lex, cont, run="world_v2_b7actfix")
    rows = {}
    for label, run, split, code in (("main", "world_v2_lexfix", "primary", "corrected"),
                                    ("country", "world_v2", "country", "frozen"),
                                    ("time", "world_v2", "temporal", "frozen"),
                                    ("cutoff", "world_v2_lexfix", "llmcut", "corrected"),
                                    ("activation_merged", "world_v2_b7actfix", "primary", "corrected"),
                                    ("activation_merged_frozen", "world_v2_actgrp", "primary", "frozen")):
        d = v36.load(run, split)
        if d is None:
            continue
        y = d["y"].astype(np.float64)
        rows[label] = {"run": run, "split": split, "code": code, "n_cells": len(y),
                       "brier": {m: float(np.mean((d[f"p_{m}"].astype(np.float64) - y) ** 2))
                                 for m in ("M0", "M1", "M1s", "M1c", "M3n") if f"p_{m}" in d}}
    return {"actfix_yardsticks": res, "brier_by_split": rows}


# %% B10: GFM の AOI の中の評価
def b10_gfm_within() -> dict:
    d = v36.load("world_v2_gfm", "primary")
    t = pd.DataFrame({"y": d["y"], "aoi": d["aoi"], "unit": d["unit"],
                      **{f"p_lf_{m}": d[f"p_{m}"] for m in ("M0", "M1", "M1c", "M3n")}})
    rng = np.random.default_rng(SEED)
    inv, cnt = v63._boot_counts(t["unit"].to_numpy(), rng)
    out = v63.b4_within(t, inv, cnt) | {"run": "world_v2_gfm (labels: GFM, refitted)"}
    # 追記 2（B10 を見た後）: 判定表の placebo との差。同じ抽出で。EMS の直した実行（B-9 の B4、種 20261029）も並べる
    out["placebo"] = within_placebo("world_v2_gfm", SEED)
    out["placebo_ems_lexfix"] = within_placebo("world_v2_lexfix", 20261029)
    return out


def within_placebo(run: str, seed: int) -> dict:
    """AOI の中だけの評価（B-9 の B4 の採点）で、名前の R・placebo の R・差の区間と判定表の行。"""
    d = v36.load(run, "primary")
    y = np.asarray(d["y"], dtype=np.float64)
    aoi = pd.factorize(pd.Series(d["aoi"]))[0]
    n_a = np.bincount(aoi).astype(np.float64)
    rate = np.bincount(aoi, weights=y) / n_a
    inner = (rate > 0) & (rate < 1)
    inv, cnt = v63._boot_counts(np.asarray(d["unit"]), np.random.default_rng(seed))
    g = cnt.shape[1]
    sse = {}
    for m in ("M0", "M1", "M1c", "M3n", "M3np"):
        p = np.clip(np.asarray(d[f"p_{m}"], dtype=np.float64), 1e-7, 1 - 1e-7)
        q = v63._shift_to_rate(np.log(p / (1 - p)), aoi, rate, n_a)
        q = np.where(inner[aoi], q, rate[aoi])
        sse[m] = np.bincount(inv, weights=(q - y) ** 2, minlength=g)
    w = np.vstack([np.ones(g), cnt])
    den = w @ sse["M0"] - w @ sse["M1"]
    r = (w @ sse["M1c"] - w @ sse["M3n"]) / den
    rp = (w @ sse["M1c"] - w @ sse["M3np"]) / den
    dif = r - rp

    def ci(x: np.ndarray, q: tuple[float, float]) -> list[float]:
        return np.quantile(x[1:], q).tolist()
    c95, c90, d95 = ci(r, (0.025, 0.975)), ci(r, (0.05, 0.95)), ci(dif, (0.025, 0.975))
    # 判定表（プロトコル §6 の順。物差しは AOI の中の地形の利得で、その 95% 区間が 0 を含めば undefined）
    yard_ci = ci(den / (w @ np.bincount(inv, minlength=g).astype(np.float64)), (0.025, 0.975))
    if yard_ci[0] <= 0 <= yard_ci[1]:
        decision = "undefined-denominator"
    elif c95[0] > 0 and not d95[0] > 0:
        decision = "not-meaning"
    elif c95[0] > 0:
        decision = ("positive-meaningful" if r[0] >= SESOI_R else
                    "positive-negligible" if -SESOI_R <= c90[0] and c90[1] <= SESOI_R else "positive-small")
    elif c95[1] < 0:
        decision = "negative"
    elif -SESOI_R <= c90[0] and c90[1] <= SESOI_R:
        decision = "equivalent-null"
    else:
        decision = "inconclusive"
    return {"run": run, "seed": seed, "replicates": int(cnt.shape[0]), "R": float(r[0]), "ci90_R": c90, "ci95_R": c95,
            "R_placebo": float(rp[0]), "ci95_R_placebo": ci(rp, (0.025, 0.975)), "R_minus_placebo": float(dif[0]),
            "ci95_R_minus_placebo": d95, "ci95_yardstick_within": yard_ci, "decision": decision,
            "n_aoi_inner": int(inner.sum())}


# %% B11: 注入したラベルでの第 5 節の推定
_I: dict = {}


def _inj_init(t: pd.DataFrame) -> None:
    _I["t"] = t
    _I["keys"] = cell_keys(t["aoi"].to_numpy(dtype=str), t["ix"].to_numpy(), t["iy"].to_numpy())
    p = np.clip(t["p_lf_M1"].to_numpy(np.float64), 1e-6, 1 - 1e-6)
    _I["p"] = p
    wet = t["lex_wetside_n1"].to_numpy() > 0
    _I["target"] = {"A": wet, "B": wet & (t["P_d1"].to_numpy() <= NEAR_M)}


def _inj_rep(args: tuple[str, float, int]) -> dict:
    form, delta, rep = args
    t = _I["t"]
    u = hashed_uniform(_I["keys"], plasmode_salt(form, delta, rep))
    y2 = inject_labels(t["y"].to_numpy().astype(np.int8), _I["p"], _I["target"][form], delta, u)
    s = t.assign(y=y2.astype(np.float64))
    out = {"form": form, "delta": delta, "rep": rep, "flips": int((y2 != t["y"].to_numpy()).sum())}
    for banded, k in ((True, f"wetside:{TERM}"), (False, "wetside:all")):
        try:
            ss, yy, z, names = _design_keep(s, "wetside", banded=banded)
            r = v39.fit(yy, z, ss["aoi"].to_numpy(), ss["unit"].to_numpy())
            out["le500m" if banded else "all"] = float(r["beta"][names.index(k)])
        except np.linalg.LinAlgError:
            out["le500m" if banded else "all"] = float("nan")
    return out


def b11_injected(t: pd.DataFrame) -> dict:
    same = float(np.max(np.abs(t["p_M1"].to_numpy(np.float64) - t["p_lf_M1"].to_numpy(np.float64))))
    t = t[["aoi", "unit", "ix", "iy", "y", "P_d1", "lex_wetside_n1", "p_M1", "p_lf_M1"]].copy()
    jobs = [(f, d, r) for f in ("A", "B") for d in INJ_DELTAS for r in range(INJ_REPS)]
    _threads1()
    t0 = time.time()
    with ProcessPoolExecutor(N_WORKERS, initializer=_inj_init, initargs=(t,)) as ex:
        res = list(ex.map(_inj_rep, jobs, chunksize=2))
    print(f"  injected: {time.time() - t0:.0f} s", flush=True)
    base = {}
    for banded, k in ((True, f"wetside:{TERM}"), (False, "wetside:all")):
        ss, yy, z, names = _design_keep(t, "wetside", banded=banded)
        base["le500m" if banded else "all"] = float(v39.fit(yy, z, ss["aoi"].to_numpy(),
                                                          ss["unit"].to_numpy())["beta"][names.index(k)])
    rows = []
    for f in ("A", "B"):
        for d in INJ_DELTAS:
            rr = [r for r in res if r["form"] == f and r["delta"] == d]
            row = {"form": f, "nominal_or": math.exp(d), "reps": len(rr),
                   "mean_flips": float(np.mean([r["flips"] for r in rr]))}
            for k in ("le500m", "all"):
                b = np.array([r[k] for r in rr])
                row[f"mean_or_{k}"] = float(np.exp(np.nanmean(b)))
                row[f"mean_or_ratio_to_base_{k}"] = float(np.exp(np.nanmean(b) - base[k]))
                row[f"sd_log_or_{k}"] = float(np.nanstd(b, ddof=1))
            rows.append(row)
    return {"base_log_or": base, "base_or": {k: math.exp(v) for k, v in base.items()}, "rows": rows,
            "max_abs_diff_p_M1_frozen_vs_corrected": same,
            "rule": "injection on held-out cells with p = corrected M1 held-out prediction, same salts as B-6 forms A/B"}


# %% B12: HAND だけの参照
def _maps(y: np.ndarray, order: np.ndarray, k: int) -> tuple[int, int]:
    return int(y[order[:k]].sum()), k


def b12_hand(t: pd.DataFrame) -> dict:
    rng = np.random.default_rng(SEED)
    acc: dict[str, list[tuple[int, int]]] = {k: [] for k in ("HAND_min", "HAND_mean", "M1", "M1c", "M3n")}
    ties = []
    for _, x in t.groupby("aoi", sort=True):
        y = x["y"].to_numpy()
        k = int(y.sum())
        if k in (0, len(y)):
            continue
        u = rng.random(len(y))
        hmin, hmean = x["hand_min"].to_numpy(np.float64), x["hand_mean"].to_numpy(np.float64)
        o1 = np.lexsort((u, hmean, hmin))
        acc["HAND_min"].append(_maps(y, o1, k))
        acc["HAND_mean"].append(_maps(y, np.lexsort((u, hmean)), k))
        cut = hmin[o1[k - 1]]
        ties.append(float(np.mean(hmin == cut)))
        for m in ("M1", "M1c", "M3n"):
            acc[m].append(_maps(y, np.argsort(-x[f"p_lf_{m}"].to_numpy(np.float64), kind="stable"), k))
    out: dict = {}
    for m, r in acc.items():
        a = np.array(r, dtype=np.float64)
        tp, k = a[:, 0], a[:, 1]
        csi = tp / (2 * k - tp)
        out[m] = {"csi_median": float(np.median(csi)), "csi_q25_q75": np.quantile(csi, [0.25, 0.75]).tolist(),
                  "hit_median": float(np.median(tp / k)), "far_median": float(np.median((k - tp) / k)),
                  "csi_pooled": float(tp.sum() / (2 * k - tp).sum()), "hit_pooled": float(tp.sum() / k.sum()),
                  "far_pooled": float((k - tp).sum() / k.sum()), "n_aoi": len(a)}
    out["ties_at_cut_hand_min"] = {"median_share_cells_tied": float(np.median(ties)),
                                   "q75": float(np.quantile(ties, 0.75))}
    out["rule"] = ("maps of k = observed flooded cells per AOI; HAND_min: minimum HAND, ties by mean HAND then "
                   "random; HAND_mean: mean HAND, ties random; models by predicted probability")
    return out


# %% B13: 言語の地域の中の地域の間の関係
def _ols_fe_cluster(yv: np.ndarray, x: np.ndarray, fe: np.ndarray) -> dict:
    """fe の中で中心化した OLS と、fe でまとめたクラスタ頑健な標準誤差（CR1）。"""
    f = pd.factorize(fe)[0]
    cnt = np.bincount(f)
    keep = cnt[f] > 1
    yv, x, f = yv[keep], x[keep], pd.factorize(f[keep])[0]

    def dm(v: np.ndarray) -> np.ndarray:
        return v - (np.bincount(f, weights=v) / np.bincount(f))[f]

    yd = dm(yv)
    xd = np.column_stack([dm(x[:, c]) for c in range(x.shape[1])])
    xtx = np.linalg.inv(xd.T @ xd)
    b = xtx @ xd.T @ yd
    e = yd - xd @ b
    g = int(f.max()) + 1
    u = np.zeros((g, x.shape[1]))
    np.add.at(u, f, xd * e[:, None])
    n, k = len(yd), x.shape[1]
    c = g / (g - 1) * (n - 1) / (n - k - g)
    v = c * xtx @ (u.T @ u) @ xtx
    se = np.sqrt(np.diag(v))
    tq = stats.t.ppf(0.975, g - 1)
    return {"coef": b.tolist(), "se": se.tolist(), "ci95": [[float(b[i] - tq * se[i]), float(b[i] + tq * se[i])]
                                                           for i in range(k)],
            "n_aoi": n, "n_activations": g}


def b13_between(t: pd.DataFrame) -> dict:
    lex = lexicon()
    a = t.groupby("aoi").agg(n=("y", "size"), k=("y", "sum"), wet=("lex_wetside_n1", lambda s: float((s > 0).mean())),
                             cc=("cc", lambda s: s[s != ""].mode().iloc[0] if (s != "").any() else ""))
    pm = np.clip(t["p_lf_M1"].to_numpy(np.float64), 1e-6, 1 - 1e-6)
    a["m1_logit"] = pd.Series(np.log(pm / (1 - pm)), index=t.index).groupby(t["aoi"]).mean()
    a["elogit"] = np.log((a["k"] + 0.5) / (a["n"] - a["k"] + 0.5))
    a["act"] = a.index.str.split("_").str[0]
    a["langs"] = a["cc"].map(lambda c: set(lex.country_langs.get(str(c).upper(), ())))
    out: dict = {"model": "empirical logit of AOI flood rate ~ share water-side + mean logit M1 | activation FE",
                 "coef_names": ["share_water_side", "mean_logit_M1"]}
    for lang in ("de", "en", "fr"):
        s = a[a["langs"].map(lambda v, lang=lang: lang in v)]
        if len(s) < 10:
            continue
        out[lang] = _ols_fe_cluster(s["elogit"].to_numpy(), s[["wet", "m1_logit"]].to_numpy(), s["act"].to_numpy())
    s = a[a["langs"].map(lambda v: bool(v & {"de", "en", "fr"}))]
    out["pooled"] = _ols_fe_cluster(s["elogit"].to_numpy(), s[["wet", "m1_logit"]].to_numpy(), s["act"].to_numpy())
    return out


# %% B14: 他の言語で読んだセル
def b14_other_language(t: pd.DataFrame) -> dict:
    import pyogrio.raw
    import shapely

    lex = lexicon()
    first = v36.first_languages(COUNTRY_INFO)
    y = t["y"].to_numpy()
    wet = t["lex_wetside_n1"].to_numpy() > 0
    near = t["P_d1"].to_numpy() <= NEAR_M
    cc = t["cc"].to_numpy(dtype=str)
    rows = {}
    for c in sorted(set(cc) - {""}):
        langs = set(lex.country_langs.get(c, ()))
        fl = first.get(c, "")
        cover = set().union(*(LANG_GROUPS.get(x, {x}) for x in langs)) if langs else set()
        if fl and fl not in cover:
            m = cc == c
            rows[c] = {"first_language": fl, "lexicon_languages": sorted(langs), "n_cells": int(m.sum()),
                       "n_water_side": int((m & wet).sum()), "n_flooded_le500m": int((m & near & (y == 1)).sum())}
    path = f"/vsizip/{NE_ADMIN1}/ne_10m_admin_1_states_provinces.shp"
    _, _, wkb, fields = pyogrio.raw.read(path, columns=["iso_a2", "name", "region", "geonunit"])
    geoms = shapely.from_wkb(wkb)
    f = {k: np.asarray(v, dtype=object) for k, v in zip(("iso_a2", "name", "region", "geonunit"), fields, strict=True)}
    lon, lat = t["lon"].to_numpy(np.float64), t["lat"].to_numpy(np.float64)
    reg = {}
    for label, (iso, field, values) in REGIONAL.items():
        sel = f["iso_a2"] == iso
        if field is not None:
            sel &= np.isin(f[field], list(values))
        idx = np.flatnonzero(sel)
        cand = (cc == iso)
        inside = np.zeros(len(t), dtype=bool)
        for i in idx:
            g = geoms[i]
            mm = cand & ~inside
            if mm.any():
                inside[np.flatnonzero(mm)[shapely.contains_xy(g, lon[mm], lat[mm])]] = True
        reg[label] = {"matched_units": sorted({str(f["name"][i]) for i in idx}), "n_units": len(idx),
                      "n_cells": int(inside.sum()), "n_water_side": int((inside & wet).sum()),
                      "n_flooded_le500m": int((inside & near & (y == 1)).sum()),
                      "n_aoi": int(t.loc[inside, "aoi"].nunique())}
    return {"country_level": rows, "regional": reg, "n_cells_total": len(t),
            "rule": "cells assigned by their nearest settlement's country; regions by point-in-polygon of the cell"}


# %% B15: N-read の区分の分布
def b15_nread(t: pd.DataFrame) -> dict:
    nr = pd.read_parquet(NREAD, columns=["gid", "category"]).drop_duplicates("gid").set_index("gid")["category"]
    g = t["P_gid1"].to_numpy()
    cat = pd.Series(g).map(nr).fillna("(none)")
    by_cells = cat.value_counts()
    uniq = pd.Series(np.unique(g[g >= 0])).map(nr).fillna("(none)").value_counts()
    return {"by_names": {k: {"n": int(v), "share": float(v / uniq.sum())} for k, v in uniq.items()},
            "by_cells": {k: {"n": int(v), "share": float(v / by_cells.sum())} for k, v in by_cells.items()}}


# %% 実行
def main() -> None:
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else set()
    res = _json(OUT) if OUT.is_file() else {}
    res["rule"] = "amendment B-10 items B1-B15 (local; post hoc; confirmatory verdict unchanged)"
    res["seed"] = SEED

    def want(k: str) -> bool:
        return not only or k in only

    def save() -> None:
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")

    t = None
    if any(want(k) for k in ("B2", "B3", "B4", "B5", "B7", "B8", "B11", "B12", "B13", "B14", "B15")):
        t = v55.table()
    for k, fn in (("B1", lambda: b1_islands()), ("B2", lambda: b2_floodplain(t)),
                  ("B3", lambda: b3_mde(pairs_or(t, "wetside", "main_wetside"))),
                  ("B4", lambda: b4_intervals(t)), ("B5", lambda: b5_grid(pairs_or(t, "wetside", "main_wetside"))),
                  ("B6", b6_validation), ("B7", lambda: b7_sensor(t)), ("B8", lambda: b8_no_africa(t)),
                  ("B9", b9_activation), ("B10", b10_gfm_within), ("B11", lambda: b11_injected(t)),
                  ("B12", lambda: b12_hand(t)), ("B13", lambda: b13_between(t)),
                  ("B14", lambda: b14_other_language(t)), ("B15", lambda: b15_nread(t))):
        if not want(k):
            continue
        t0 = time.time()
        res[k] = fn()
        print(k, f"({time.time() - t0:.0f} s)", json.dumps(res[k], default=float)[:700], flush=True)
        save()


if __name__ == "__main__":
    main()
