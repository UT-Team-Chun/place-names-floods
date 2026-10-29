# %% v68（事後・amendment B-11 の B1〜B11）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-11 で計算の前に固定した量を、既存の表と held-out の予測から計算する。

- B1: 活性化をまとめた分割の、凍結 → 直したコードでの R の変化を場所グループに分ける（ノルウェーとそれ以外）
- B2: 集落の近くの水系側の OR を、上限をかけない全ての評価できるセルで（EMS、GFM、同じセルの EMS、Sentinel-1 以外）
- B3: 余白と上限を、第 5 節の推定量の OR の尺度に換算（B-10 の B11 の名目 → 実現の曲線）
- B4: GFM の評価域の検定（EMS で学習した予測を、GFM が評価できるセルの EMS ラベルで採点）と GFM の画像の数
- B5: 浸水の割合の分母（セルの数、観測された非恒久水域を分母にしたラベルでの OR と R）
- B6: fold の割り当て 13 通りの R からの Nadeau–Bengio の区間
- B7: 二段目の当て直し + 二方向の重み + 比の bootstrap（組み合わせた区間）
- B8: 探索的な推定の目録（v36〜v67 の JSON を機械的に走査）
- B9: 解析したセルの最寄りの集落と水系側の一致を GeoNames の地物の符号ごとに、丸めた座標の割合を国ごとに
- B10: パイロット（pilot3、δ = 0）の M1s → M3s と placebo の増分
- B11: 地名学者の照合のための、形 × 国で層別した標本（答えは別のファイル、公開しない）

出力: stdout と docs/verification/v68_b11_local.json（部分ごとに足す。`--only B1,B4` で一部だけ）。
重い計算: B2 と B5 の上限なしの pairs bootstrap（1,900 万セル）、B7（500 回 × 2）。作業プロセスの数は V68_WORKERS。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats
from scipy.spatial import cKDTree

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

from topsem.config import DATA_DIR, RUNS_DIR
from topsem.world.b8 import stage2_weighted
from topsem.world.lexicon import _SPLIT, WET_SIDE, Lexicon, _hit

OUT = ROOT / "docs" / "verification" / "v68_b11_local.json"
VER = DATA_DIR / "verification"
DOCV = ROOT / "docs" / "verification"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
LABELS = DATA_DIR / "global" / "labels"
GFM_SUMMARY = DATA_DIR / "global" / "gfm_summary.json"
SEED = 20261031
SESOI_R = 0.05
B_PAIRS = 2000
B_BOOT = 4000
B_COMB = 500
NEAR_M = 500.0
Z80 = 2.8  # 両側 5%・検出力 80%
N_WORKERS = int(os.environ.get("V68_WORKERS", "6"))
DENOM_MIN = 0.25  # B5: 分母（観測された非恒久水域）がこれ未満のセルは除く
AUDIT_PER_STRATUM = 3
AUDIT_MAX = 150


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v39 = _load("v39_within_aoi_odds")
v55 = _load("v55_b7_local")
v63 = _load("v63_b9_local")
v36 = _load("v36_yardstick_robustness")
v39, v55 = (sys.modules[n] for n in ("v39_within_aoi_odds", "v55_b7_local"))


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _threads1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        os.environ[k] = "1"


# %% 共通: 直したコードの読み（集落ごとの水系側）と pairs bootstrap
_LEX: Lexicon | None = None


def lexicon() -> Lexicon:
    global _LEX
    if _LEX is None:
        _LEX = Lexicon(LEX)
    return _LEX


def wetside_by_gid(gids: np.ndarray) -> pd.Series:
    """集落ごとの水系側（最寄りの名の照合、直したコード: Lexicon.match、WET_SIDE の区分のどれか）。"""
    lex = lexicon()
    gn = (
        pd.read_parquet(GN, columns=["gid", "asciiname", "cc"])
        .drop_duplicates("gid")
        .set_index("gid")
    )
    u = np.unique(gids[gids >= 0])
    u = u[np.isin(u, gn.index.to_numpy())]
    names, ccs = (
        gn.loc[u, "asciiname"].astype(str).to_numpy(),
        gn.loc[u, "cc"].astype(str).to_numpy(),
    )
    flag = np.array(
        [
            int(any(lex.match(n, c)[k] > 0 for k in WET_SIDE))
            for n, c in zip(names, ccs, strict=True)
        ],
        dtype=np.float64,
    )
    return pd.Series(flag, index=u)


_WW: dict = {}


def fit_w(
    y: np.ndarray, z: np.ndarray, aoi: np.ndarray, f: np.ndarray, iters: int = 50
) -> np.ndarray:
    """v39.fit の点推定を度数の重み f で解く。

    pairs bootstrap で k 回引いた場所グループの k 個の複写（v55._pairs_rep は複写ごとに AOI の切片を持つ）は、同じ
    データなので IRLS の各反復で同じ切片の更新を受ける。複写を並べる代わりに重み k を掛けると、ヘッセ行列・スコア・
    Schur 補行列・収束の判定がどれも同じになる（反復ごとに一致。違いは和の順序の丸めだけ）。上限なしの表（1,900 万
    セル）でメモリが足りるよう、AOI の指示行列を作らず、説明変数の列ごとに bincount と内積で集計する。
    """
    ua, ia = np.unique(aoi, return_inverse=True)
    a, k = len(ua), z.shape[1]
    rate = np.bincount(ia, weights=f * y, minlength=a) / np.bincount(
        ia, weights=f, minlength=a
    )
    alpha = np.log(rate / (1 - rate))
    beta = np.zeros(k)
    haz, hzz = np.empty((a, k)), np.empty((k, k))
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(alpha[ia] + z @ beta)))
        w = f * p * (1 - p)
        r = f * (y - p)
        del p
        haa = np.bincount(ia, weights=w, minlength=a)
        for jj in range(k):
            zw = z[:, jj] * w
            haz[:, jj] = np.bincount(ia, weights=zw, minlength=a)
            hzz[:, jj] = z.T @ zw
        del zw, w
        ga, gz = np.bincount(ia, weights=r, minlength=a), z.T @ r
        del r
        s_ = hzz - haz.T @ (haz / haa[:, None])
        db = np.linalg.solve(s_, gz - haz.T @ (ga / haa))
        da = (ga - haz @ db) / haa
        alpha += da
        beta += db
        if max(np.abs(da).max(), np.abs(db).max()) < 1e-10:
            break
    return beta


def _pairs_init_w(
    y: np.ndarray, z: np.ndarray, aoi: np.ndarray, unit: np.ndarray, g: int
) -> None:
    _WW.update(y=y, z=z, aoi=aoi, unit=unit, g=g)


def _pairs_rep_w(sel: np.ndarray) -> np.ndarray | None:
    cnt = np.bincount(sel, minlength=_WW["g"]).astype(np.float64)
    fc = cnt[_WW["unit"]]
    idx = np.flatnonzero(fc > 0)
    try:
        return fit_w(_WW["y"][idx], _WW["z"][idx], _WW["aoi"][idx], fc[idx])
    except np.linalg.LinAlgError:
        return None


def _pairs_arrays(s0: pd.DataFrame) -> tuple:
    y = s0["y"].to_numpy(dtype=np.float64)
    rate = pd.Series(y).groupby(s0["aoi"].to_numpy()).transform("mean").to_numpy()
    keep = (rate > 0) & (rate < 1)
    s = s0.loc[keep].reset_index(drop=True)
    y = y[keep]
    z, names = v39.design(s, ("wetside",), True, banded=True)
    aoi = pd.factorize(s["aoi"])[0].astype(np.int64)
    unit_lab, unit = np.unique(s["unit"].to_numpy(), return_inverse=True)
    rng = np.random.default_rng(SEED)
    sels = [rng.integers(0, len(unit_lab), size=len(unit_lab)) for _ in range(B_PAIRS)]
    return y, z, names, aoi, unit, len(unit_lab), sels


def pairs_or(t: pd.DataFrame, key: str, ycol: str = "y") -> dict:
    """500 m 以内の水系側のオッズ比（第 5 節の式、v55.odds）と、場所の pairs cluster bootstrap（種 SEED）。

    再標本は v55 と同じ（種 SEED から場所グループを B_PAIRS 回引く）。B2・B5 の上限なしの表（1,900 万セル）は複写を
    並べる形だとメモリが足りないので、度数の重みの形（fit_w）で解く（check_weighted で一致を確かめた）。
    """
    s0 = t if ycol == "y" else t.assign(y=t[ycol].to_numpy())
    o = v55.odds(s0, cls="wetside")
    cache = VER / f"v68_pairs_{key}.npy"
    if cache.is_file():
        b = np.load(cache)
    else:
        y, z, names, aoi, unit, g, sels = _pairs_arrays(s0)
        j = names.index("wetside:le500m")
        _threads1()
        t0 = time.time()
        with ProcessPoolExecutor(
            N_WORKERS, initializer=_pairs_init_w, initargs=(y, z, aoi, unit, g)
        ) as ex:
            res = list(ex.map(_pairs_rep_w, sels, chunksize=4))
        b = np.array([r[j] if r is not None else np.nan for r in res])
        np.save(cache, b)
        print(f"  pairs {key}: {time.time() - t0:.0f} s ({len(y):,} cells)", flush=True)
    ok = np.isfinite(b)
    sd = float(np.std(b[ok], ddof=1))
    o["pairs_ci95"] = [
        float(np.exp(np.quantile(b[ok], 0.025))),
        float(np.exp(np.quantile(b[ok], 0.975))),
    ]
    o["pairs_log_or_sd"] = sd
    o["pairs_reps_ok"] = int(ok.sum())
    o["mde80_up"], o["mde80_down"] = float(np.exp(Z80 * sd)), float(np.exp(-Z80 * sd))
    near = (s0["P_d1"].to_numpy() <= NEAR_M) & (s0["lex_wetside_n1"].to_numpy() > 0)
    fl = near & (s0["y"].to_numpy() == 1)
    o["n_flooded_near"] = int(fl.sum())
    o["n_units_flooded_near"] = int(s0.loc[fl, "unit"].nunique())
    return o


# %% 上限なしの表（B2・B5）
def uncapped_full() -> pd.DataFrame:
    """上限なしの感度の実行（world_v2_s2a_nocap）の held-out の予測に、セルの ix・iy・最寄りの集落・浸水の被覆・GFM を足す。"""
    cache = VER / "v68_uncapped.parquet"
    if cache.is_file():
        return pd.read_parquet(cache)
    u = pd.read_parquet(
        VER / "v57_uncapped.parquet"
    )  # aoi, unit, y, p_M1, lex_wetside_n1（凍結）, lon, lat, P_d1
    u = u.rename(columns={"lex_wetside_n1": "wetside_frozen"})
    cols_c = [
        "ix",
        "iy",
        "lon",
        "lat",
        "f_riverine",
        "f_other",
        "f_flood_any",
        "f_observed",
        "f_perm_water",
    ]
    add = {
        c: np.full(len(u), np.nan)
        for c in (
            "ix",
            "iy",
            "P_gid1",
            "f_riverine",
            "f_other",
            "f_flood_any",
            "f_observed",
            "f_perm_water",
            "gfm_eval",
            "gfm_y",
        )
    }
    t0 = time.time()
    for aoi, ia in u.groupby("aoi").indices.items():
        c = pq.read_table(LABELS / aoi / "cells.parquet", columns=cols_c).to_pandas()
        f = pq.read_table(
            LABELS / aoi / "features.parquet", columns=["P_gid1"]
        ).to_pandas()
        tree = cKDTree(c[["lon", "lat"]].to_numpy(np.float64))
        dist, j = tree.query(u.iloc[ia][["lon", "lat"]].to_numpy(np.float64))
        ok = dist <= 1e-5
        for k in (
            "ix",
            "iy",
            "f_riverine",
            "f_other",
            "f_flood_any",
            "f_observed",
            "f_perm_water",
        ):
            v = c[k].to_numpy(np.float64)[j]
            v[~ok] = np.nan
            add[k][ia] = v
        g = f["P_gid1"].to_numpy(np.float64)[j]
        g[~ok] = np.nan
        add["P_gid1"][ia] = g
        gp = LABELS / aoi / "gfm.parquet"
        if gp.is_file():
            gf = pq.read_table(
                gp, columns=["ix", "iy", "gfm_eval", "gfm_y"]
            ).to_pandas()
            key = pd.MultiIndex.from_arrays([gf["ix"].to_numpy(), gf["iy"].to_numpy()])
            gm = pd.DataFrame(
                {
                    "e": gf["gfm_eval"].astype(float).to_numpy(),
                    "y": gf["gfm_y"].astype(float).to_numpy(),
                },
                index=key,
            )
            q = pd.MultiIndex.from_arrays(
                [add["ix"][ia].astype(np.int64), add["iy"][ia].astype(np.int64)]
            )
            r = gm.reindex(q)
            add["gfm_eval"][ia] = r["e"].to_numpy()
            add["gfm_y"][ia] = r["y"].to_numpy()
    for k, v in add.items():
        u[k] = v
    print(
        f"uncapped table: {len(u):,} rows, unmatched {int(np.isnan(u['ix']).sum())} ({time.time() - t0:.0f} s)",
        flush=True,
    )
    u.to_parquet(cache, index=False)
    return u


def with_corrected_wetside(u: pd.DataFrame) -> pd.DataFrame:
    gids = u["P_gid1"].fillna(-1).to_numpy(np.int64)
    ws = wetside_by_gid(gids)
    u["lex_wetside_n1"] = ws.reindex(gids).fillna(0.0).to_numpy(np.float64)
    return u


# %% B1: 活性化をまとめた分割の変化
def b1_activation_change() -> dict:
    fr, co = (
        v36.load("world_v2_actgrp", "primary"),
        v36.load("world_v2_b7actfix", "primary"),
    )
    same = all(
        np.array_equal(np.asarray(fr[c]).astype(str), np.asarray(co[c]).astype(str))
        for c in ("aoi", "unit", "y")
    )
    out: dict = {"rows_aligned": bool(same)}
    if not same:
        return out | {"note": "rows of the two runs are not aligned"}
    y = fr["y"].astype(np.float64)
    out["max_abs_diff_p_M1"] = float(np.abs(fr["p_M1"] - co["p_M1"]).max())
    out["max_abs_diff_p_M0"] = float(np.abs(fr["p_M0"] - co["p_M0"]).max())
    units, ucode = np.unique(np.asarray(fr["unit"]).astype(str), return_inverse=True)
    nu = len(units)
    cc = np.asarray(fr.get("cc", np.array([""] * len(y)))).astype(str)
    nor_units = (
        set(units[np.unique(ucode[cc == "NO"])]) if (cc == "NO").any() else set()
    )

    def unit_d(d: dict) -> np.ndarray:
        e1 = (d["p_M1c"].astype(np.float64) - y) ** 2
        e3 = (d["p_M3n"].astype(np.float64) - y) ** 2
        return np.bincount(ucode, weights=e1 - e3, minlength=nu)

    yard = float(
        np.mean(
            (fr["p_M0"].astype(np.float64) - y) ** 2
            - (fr["p_M1"].astype(np.float64) - y) ** 2
        )
    )
    df, dc = unit_d(fr), unit_d(co)
    n = len(y)
    out |= {
        "R_frozen": float(df.sum() / n / yard),
        "R_corrected": float(dc.sum() / n / yard),
        "yardstick": yard,
        "change_R": float((dc - df).sum() / n / yard),
    }
    is_no = np.array([u in nor_units for u in units])
    out["change_R_norwegian_units"] = float((dc - df)[is_no].sum() / n / yard)
    out["change_R_other_units"] = float((dc - df)[~is_no].sum() / n / yard)
    out["n_norwegian_units"] = int(is_no.sum())
    order = np.argsort(-np.abs(dc - df))[:8]
    cell_cc = (
        pd.Series(cc)
        .groupby(ucode)
        .agg(lambda s: s.value_counts().index[0] if len(s) else "")
    )
    out["largest_unit_changes"] = [
        {
            "unit": str(units[i]),
            "change_R": float((dc - df)[i] / n / yard),
            "majority_cc": str(cell_cc.get(i, "")),
            "n_cells": int((ucode == i).sum()),
        }
        for i in order
    ]
    # 名前の特徴の変わったセルの割合（凍結と直したコードの読みが違うセル）
    if "s_wetside" in fr and "s_wetside" in co:
        out["share_cells_wetside_changed"] = float(
            (fr["s_wetside"] != co["s_wetside"]).mean()
        )
    return out


# %% B2: 上限なしの集落の近くの OR
KEEP = (
    "aoi",
    "unit",
    "wetside_frozen",
    "y",
    "P_d1",
    "P_gid1",
    "p_M1",
    "gfm_eval",
    "gfm_y",
    "f_riverine",
    "f_other",
    "f_flood_any",
    "f_observed",
    "f_perm_water",
)


def b2_uncapped() -> dict:
    u = uncapped_full()
    u = with_corrected_wetside(u[[c for c in KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    out: dict = {
        "n_cells": len(u),
        "share_wetside_changed_vs_frozen": float(
            ((u["lex_wetside_n1"] > 0) != (u["wetside_frozen"] > 0)).mean()
        ),
    }
    out["ems_all"] = pairs_or(u, "b2a_ems_uncapped")
    g = u.loc[u["gfm_eval"] == 1].reset_index(drop=True)
    out["n_cells_gfm"] = len(g)
    out["n_aoi_gfm"] = int(g["aoi"].nunique())
    out["gfm"] = pairs_or(g, "b2b_gfm_uncapped", ycol="gfm_y")
    out["ems_same_cells"] = pairs_or(g, "b2c_ems_gfmcells")
    sens = {a: s.get("sensor", "none") for a, s in _json(GFM_SUMMARY)["aois"].items()}
    nons1 = g["aoi"].map(lambda a: sens.get(a, "none") != "s1").to_numpy()
    h = g.loc[nons1].reset_index(drop=True)
    out["n_aoi_not_s1"] = int(h["aoi"].nunique())
    out["gfm_not_s1"] = pairs_or(h, "b2d_gfm_nots1", ycol="gfm_y")
    out["ems_not_s1"] = pairs_or(h, "b2d_ems_nots1")
    out["rule"] = (
        "B-11 B2: Sect. 5 estimator on all evaluable cells (uncapped run), corrected-code readings"
    )
    return out


# %% B3: 推定量の尺度への換算
def b3_estimator_scale() -> dict:
    v66 = _json(DOCV / "v66_b10_local.json")
    rows = v66["B11"]["rows"]
    curves = {}
    for form, col in (
        ("A", "mean_or_ratio_to_base_all"),
        ("B", "mean_or_ratio_to_base_le500m"),
    ):
        rr = sorted((r["nominal_or"], r[col]) for r in rows if r["form"] == form)
        x = np.log(np.array([1.0] + [a for a, _ in rr]))
        yv = np.log(np.array([1.0] + [b for _, b in rr]))
        curves[form] = {
            "nominal": np.exp(x).tolist(),
            "estimator": np.exp(yv).tolist(),
            "monotone": bool(np.all(np.diff(yv) > 0)),
        }

    def conv(form: str, nominal: float | None) -> dict | None:
        if nominal is None or not np.isfinite(nominal):
            return None
        c = curves[form]
        x, yv = np.log(c["nominal"]), np.log(c["estimator"])
        if np.log(nominal) > x[-1]:
            return {
                "nominal": nominal,
                "estimator": None,
                "text": f"> {c['estimator'][-1]:.2f}",
            }
        e = float(np.exp(np.interp(np.log(nominal), x, yv)))
        return {"nominal": nominal, "estimator": e, "text": f"{e:.2f}"}

    v59, v64 = (
        _json(DOCV / "v59_b8_aggregate.json"),
        _json(DOCV / "v64_b9_aggregate.json"),
    )
    out: dict = {
        "curves": curves,
        "rule": "monotone interpolation on log scale of B-10 B11 (nominal → estimator, "
        "ratio to the no-injection reference); beyond nominal 3 written as >",
    }
    mo = (v59.get("B1") or {}).get("margin_or") or {}
    out["margin"] = {
        f: conv(f, (mo.get(f) or {}).get("odds_ratio")) for f in ("A", "B")
    }
    bb = (v64.get("B1_B2") or {}).get("bounds") or {}
    out["bounds"] = {
        iv: {
            ref: {
                f: conv(
                    f, ((bb.get(iv) or {}).get(ref) or {}).get(f, {}).get("odds_ratio")
                )
                for f in ("A", "B")
            }
            for ref in (bb.get(iv) or {})
        }
        for iv in bb
    }
    # 仮定を緩めた上限（入っていれば同じ規則で）: B-10 の C4（反転の相関）、B-9 の C3（集落ごとに違う効果、平均の OR）
    v67 = (
        _json(DOCV / "v67_b10_aggregate.json")
        if (DOCV / "v67_b10_aggregate.json").is_file()
        else {}
    )
    c4 = v67.get("C4_ac") if isinstance(v67.get("C4_ac"), dict) else None
    out["bound_ac"] = conv(
        "A", ((c4 or {}).get("bound") or {}).get("b7_recomputed", {}).get("odds_ratio")
    )
    c3 = v64.get("C3") if isinstance(v64.get("C3"), dict) else None
    out["bound_het"] = conv(
        "A", ((c3 or {}).get("bound") or {}).get("b7_recomputed", {}).get("odds_ratio")
    )
    return out


# %% B4: GFM の評価域
def _gfm_join(d: dict) -> np.ndarray:
    """予測の行（aoi・lon・lat）に、GFM のラベルが定まるか（gfm_eval）を付ける。GFM の無い AOI は NaN。"""
    aoi = np.asarray(d["aoi"]).astype(str)
    # B-3 の GFM の解析に含めた AOI だけ（GFM のラベルが定まるセルの定義）
    inc = {a for a, s in _json(GFM_SUMMARY)["aois"].items() if s.get("include")}
    ev = np.full(len(aoi), np.nan)
    lon, lat = (
        np.asarray(d["lon"], dtype=np.float64),
        np.asarray(d["lat"], dtype=np.float64),
    )
    for a, ia in pd.Series(np.arange(len(aoi))).groupby(aoi).indices.items():
        gp = LABELS / a / "gfm.parquet"
        if a not in inc or not gp.is_file():
            continue
        c = pq.read_table(
            LABELS / a / "cells.parquet", columns=["ix", "iy", "lon", "lat"]
        ).to_pandas()
        gf = pq.read_table(gp, columns=["ix", "iy", "gfm_eval"]).to_pandas()
        tree = cKDTree(c[["lon", "lat"]].to_numpy(np.float64))
        dist, j = tree.query(np.column_stack([lon[ia], lat[ia]]))
        key = pd.MultiIndex.from_arrays([c["ix"].to_numpy()[j], c["iy"].to_numpy()[j]])
        m = pd.Series(
            gf["gfm_eval"].astype(float).to_numpy(),
            index=pd.MultiIndex.from_arrays([gf["ix"].to_numpy(), gf["iy"].to_numpy()]),
        )
        v = np.array(m.reindex(key).to_numpy(np.float64), copy=True)
        v[dist > 1e-5] = np.nan
        ev[ia] = v
    return ev


def _r_boot(
    y: np.ndarray,
    preds: dict[str, np.ndarray],
    unit: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    inv, cnt = v63._boot_counts(unit, rng)
    g = cnt.shape[1]
    sse = {
        k: np.bincount(inv, weights=(v.astype(np.float64) - y) ** 2, minlength=g)
        for k, v in preds.items()
    }
    nn = np.bincount(inv, minlength=g).astype(np.float64)
    w = np.vstack([np.ones(g), cnt])
    yard = (w @ sse["M0"] - w @ sse["M1"]) / (w @ nn)
    dbs = (w @ sse["M1c"] - w @ sse["M3n"]) / (w @ nn)
    r = dbs / yard
    q = lambda x, a, b: np.quantile(x[1:], [a, b]).tolist()
    return {
        "yardstick": float(yard[0]),
        "ci95_yardstick": q(yard, 0.025, 0.975),
        "delta_bs": float(dbs[0]),
        "ci95_delta_bs": q(dbs, 0.025, 0.975),
        "R": float(r[0]),
        "ci90_R": q(r, 0.05, 0.95),
        "ci95_R": q(r, 0.025, 0.975),
        "n_cells": len(y),
        "n_units": g,
        "replicates": int(cnt.shape[0]),
    }


def b4_gfm_domain() -> dict:
    d = v36.load("world_v2_lexfix", "primary")
    ev = _gfm_join(d)
    y = d["y"].astype(np.float64)
    unit = np.asarray(d["unit"]).astype(str)
    preds = {m: d[f"p_{m}"] for m in ("M0", "M1", "M1c", "M3n")}
    in_gfm = np.isfinite(ev)
    rng = np.random.default_rng(SEED)
    out: dict = {
        "gfm_aois_all_cells": _r_boot(
            y[in_gfm], {k: v[in_gfm] for k, v in preds.items()}, unit[in_gfm], rng
        )
    }
    m = ev == 1
    out["gfm_assessable_cells"] = _r_boot(
        y[m], {k: v[m] for k, v in preds.items()}, unit[m], rng
    )
    out["n_aoi_gfm"] = int(
        pd.Series(np.asarray(d["aoi"]).astype(str)[in_gfm]).nunique()
    )
    # GFM の画像の数（AOI ごと、要約の値か gfm.parquet の gfm_n_scenes の最大）
    n_sc = []
    for a in pd.unique(np.asarray(d["aoi"]).astype(str)[in_gfm]):
        gf = pq.read_table(
            LABELS / a / "gfm.parquet", columns=["gfm_n_scenes"]
        ).to_pandas()
        n_sc.append(int(gf["gfm_n_scenes"].max()))
    n_sc = np.array(n_sc)
    out["gfm_scenes_per_aoi"] = {
        "median": float(np.median(n_sc)),
        "q25": float(np.quantile(n_sc, 0.25)),
        "q75": float(np.quantile(n_sc, 0.75)),
        "min": int(n_sc.min()),
        "max": int(n_sc.max()),
        "share_more_than_one": float((n_sc > 1).mean()),
    }
    g3 = _json(DOCV / "v47_b5_cluster.json").get("gfm_c2") or {}
    out["gfm_refit_yardstick"] = g3.get("yardstick") if isinstance(g3, dict) else None
    return out


# %% B5: 浸水の割合の分母
def _relabel(
    df: pd.DataFrame, untyped: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    den = np.clip(
        df["f_observed"].to_numpy(np.float64) - df["f_perm_water"].to_numpy(np.float64),
        0.0,
        None,
    )
    ok = den >= DENOM_MIN
    dd = np.where(ok, den, 1.0)
    riv = df["f_riverine"].to_numpy(np.float64) + np.where(
        untyped, df["f_other"].to_numpy(np.float64), 0.0
    )
    y1 = riv / dd >= 0.5
    y0 = df["f_flood_any"].to_numpy(np.float64) / dd < 0.5
    return y1.astype(np.int8), (y1 | y0) & ok, den


def _capped_with_obs() -> pd.DataFrame:
    """上限の表（直したコードの読み、v55.table）に、セルの観測の割合と恒久水域の割合を足す。"""
    t = v55.table()  # f_riverine・f_other・f_flood_any・ix・iy
    obs = np.full(len(t), np.nan)
    pw = np.full(len(t), np.nan)
    for a, ia in t.groupby("aoi").indices.items():
        c = pq.read_table(
            LABELS / a / "cells.parquet",
            columns=["ix", "iy", "f_observed", "f_perm_water"],
        ).to_pandas()
        m = c.set_index(["ix", "iy"])
        r = m.reindex(
            pd.MultiIndex.from_arrays([t["ix"].to_numpy()[ia], t["iy"].to_numpy()[ia]])
        )
        obs[ia], pw[ia] = r["f_observed"].to_numpy(), r["f_perm_water"].to_numpy()
    t["f_observed"], t["f_perm_water"] = obs, pw
    return t


def b5_denominator() -> dict:
    _, infos, _ = v55._infos()
    t = _capped_with_obs()
    pw, obs = t["f_perm_water"].to_numpy(), t["f_observed"].to_numpy()
    near = t["P_d1"].to_numpy() <= NEAR_M
    ws = t["lex_wetside_n1"].to_numpy() > 0
    part = ((pw > 0) & (pw < 0.5)) | (obs < 1)
    out: dict = {
        "counts_capped": {
            "all": {
                "n": len(t),
                "partial": int(part.sum()),
                "share": float(part.mean()),
            },
            "le500m_wetside": {
                "n": int((near & ws).sum()),
                "partial": int((part & near & ws).sum()),
                "share": float(part[near & ws].mean()),
            },
            "le500m_other": {
                "n": int((near & ~ws).sum()),
                "partial": int((part & near & ~ws).sum()),
                "share": float(part[near & ~ws].mean()),
            },
        }
    }
    untyped = (
        t["aoi"]
        .map({a: bool(infos[a].unspecified_as_riverine) for a in t["aoi"].unique()})
        .to_numpy()
    )
    y_new, keep, _ = _relabel(t, untyped)
    out["relabel_capped"] = {
        "n_kept": int(keep.sum()),
        "n_dropped": int((~keep).sum()),
        "n_changed_0_to_1": int(((t["y"].to_numpy() == 0) & (y_new == 1) & keep).sum()),
        "n_changed_1_to_0": int(((t["y"].to_numpy() == 1) & (y_new == 0) & keep).sum()),
    }
    tk = t.loc[keep].reset_index(drop=True).assign(y_obsden=y_new[keep])
    out["or_capped"] = pairs_or(tk, "b5_capped_obsden", ycol="y_obsden")
    # 主比較の R（直した実行の既存の予測で採点し直す）
    d = v36.load("world_v2_lexfix", "primary")
    same = np.array_equal(
        np.asarray(d["y"]).astype(int), t["y"].to_numpy(int)
    ) and np.array_equal(np.asarray(d["aoi"]).astype(str), t["aoi"].to_numpy(str))
    out["predictions_aligned_with_table"] = bool(same)
    if same:
        preds = {m: np.asarray(d[f"p_{m}"])[keep] for m in ("M0", "M1", "M1c", "M3n")}
        out["R_rescored"] = _r_boot(
            y_new[keep].astype(np.float64),
            preds,
            t["unit"].to_numpy(str)[keep],
            np.random.default_rng(SEED),
        )
    return out


def b5_uncapped() -> dict:
    _, infos, _ = v55._infos()
    u = uncapped_full()
    u = with_corrected_wetside(u[[c for c in KEEP if c in u.columns]].copy())
    u = u.loc[
        u["P_gid1"].notna() & u["P_d1"].notna() & u["f_observed"].notna()
    ].reset_index(drop=True)
    near = u["P_d1"].to_numpy() <= NEAR_M
    ws = u["lex_wetside_n1"].to_numpy() > 0
    pw, obs = u["f_perm_water"].to_numpy(), u["f_observed"].to_numpy()
    part = ((pw > 0) & (pw < 0.5)) | (obs < 1)
    out: dict = {
        "counts_uncapped": {
            "all": {
                "n": len(u),
                "partial": int(part.sum()),
                "share": float(part.mean()),
            },
            "le500m_wetside": {
                "n": int((near & ws).sum()),
                "share": float(part[near & ws].mean()),
            },
            "le500m_other": {
                "n": int((near & ~ws).sum()),
                "share": float(part[near & ~ws].mean()),
            },
        }
    }
    untyped = (
        u["aoi"]
        .map(
            {
                a: bool(infos[a].unspecified_as_riverine) if a in infos else False
                for a in u["aoi"].unique()
            }
        )
        .to_numpy()
    )
    y_new, keep, _ = _relabel(u, untyped)
    uk = u.loc[keep].reset_index(drop=True).assign(y_obsden=y_new[keep])
    out["or_uncapped"] = pairs_or(uk, "b5_uncapped_obsden", ycol="y_obsden")
    return out


# %% A4: 二重字（ae・oe・ue）を母音として読むと失われる水系側の照合の内訳（記述、報告の訂正のため）
def a4_digraphs() -> dict:
    t = pd.read_parquet(VER / "v37_table.parquet", columns=["P_gid1"])
    gids = np.unique(t["P_gid1"].to_numpy(np.int64))
    gn = (
        pd.read_parquet(GN, columns=["gid", "asciiname", "cc"])
        .drop_duplicates("gid")
        .set_index("gid")
    )
    gn = gn.loc[np.intersect1d(gids, gn.index.to_numpy())]
    lex = lexicon()

    def wet_entries(name: str, cc: str) -> set:
        toks = [x for x in _SPLIT.split(name.lower()) if x]
        return {
            (e.cat, e.form)
            for e in lex.entries_for(cc)
            if e.cat in WET_SIDE and any(_hit(tk, e) for tk in toks)
        }

    lost: dict[str, int] = {}
    by_cc: dict[str, int] = {}
    n = 0
    for nm, cc in zip(gn["asciiname"].astype(str), gn["cc"].astype(str), strict=True):
        if not re.search(r"ae|oe|ue", nm.lower()):
            continue
        before = wet_entries(nm, cc)
        if not before:
            continue
        folded = re.sub(
            r"(?i)ue", "u", re.sub(r"(?i)oe", "o", re.sub(r"(?i)ae", "a", nm))
        )
        if wet_entries(folded, cc):
            continue
        n += 1
        by_cc[cc] = by_cc.get(cc, 0) + 1
        for cat, form in before:
            k = f"{cc}:{form}:{cat}"
            lost[k] = lost.get(k, 0) + 1
    return {
        "n_settlements_losing_wetside": n,
        "by_country": dict(sorted(by_cc.items(), key=lambda kv: -kv[1])),
        "by_form": dict(sorted(lost.items(), key=lambda kv: -kv[1])),
        "rule": "analysed settlements whose water-side match disappears when ae/oe/ue are read as a/o/u",
    }


# %% C1: 国外しの分割での、組み立て直した地域構成の置換（クラスタの world_v2_b11ctryperm を集める。v56 と同じ集計）
C1_RUN = "world_v2_b11ctryperm"
C1_PERMS = 100


def c1_country_perm() -> dict | str:
    d = RUNS_DIR / C1_RUN / "predictions" / "country"
    files = sorted(d.glob("fold*_permsrc_*.npz")) if d.is_dir() else []
    folds: dict[int, list] = {}
    for f in files:
        k = int(re.match(r"fold(\d+)_", f.name).group(1))
        folds.setdefault(k, []).append(f)
    if not folds:
        return "missing"
    met = _json(
        RUNS_DIR / "world_v2" / "metrics_country.json"
    )  # 基礎のモデルは辞書の直しで変わらない（B1）
    yard = float(met["primary"]["delta_bs_M0_M1"])
    real = {a: [] for a in ("M1c", "M3n")}
    sse = {k: [] for k in ("ref", "names")}
    n_unit = []
    for k in sorted(folds):
        parts = [np.load(f) for f in sorted(folds[k])]
        units = parts[0]["units"]
        if not all(np.array_equal(q["units"], units) for q in parts):
            return f"fold {k}: unit order differs between permutation ranges"
        n_unit.append(parts[0]["n_unit"].astype(np.float64))
        for a, acc in real.items():
            acc.append(parts[0][f"real_{a}"])
        for sk, acc in sse.items():
            acc.append(np.concatenate([q[f"sse_{sk}"] for q in parts], axis=0))
    n_perm = {sk: {v.shape[0] for v in acc} for sk, acc in sse.items()}
    if any(len(v) != 1 or next(iter(v)) != C1_PERMS for v in n_perm.values()):
        return f"incomplete (folds {sorted(folds)}, permutations {n_perm})"
    n = float(np.concatenate(n_unit).sum())
    tot = {a: float(np.concatenate(v).sum()) for a, v in real.items()}
    ssum = {sk: np.concatenate(v, axis=1).sum(axis=1) for sk, v in sse.items()}
    d_perm = (ssum["ref"] - ssum["names"]) / n
    d_real = (tot["M1c"] - tot["M3n"]) / n
    r = d_perm / yard
    return {
        "n_perms": C1_PERMS,
        "n_folds": len(folds),
        "n_cells": int(n),
        "yardstick_country_frozen_run": yard,
        "R_real": d_real / yard,
        "R_perm_mean": float(r.mean()),
        "R_perm_q025": float(np.quantile(r, 0.025)),
        "R_perm_q975": float(np.quantile(r, 0.975)),
        "percentile": float(np.mean(d_perm < d_real)),
        "p_one_sided": float(np.mean(d_perm >= d_real)),
        "rule": "B-11 C1: country split, reassignments 1-100 within AOIs with the regional composition recomputed "
        "(reference M1 + c(q), names M1 + c(q) + q); percentile of the real increment and share of reassignments "
        "at least as large",
    }


# %% 度数の重みの形の確かめ（B5 の上限の表の 2,000 回を、複写を並べる形の保存した値と比べる）
def check_weighted() -> dict:
    _, infos, _ = v55._infos()
    t = _capped_with_obs()
    untyped = (
        t["aoi"]
        .map({a: bool(infos[a].unspecified_as_riverine) for a in t["aoi"].unique()})
        .to_numpy()
    )
    y_new, keep, _ = _relabel(t, untyped)
    tk = t.loc[keep].reset_index(drop=True).assign(y=y_new[keep])
    ref = np.load(
        VER / "v68_pairs_b5_capped_obsden.npy"
    )  # v55._pairs_rep（複写を並べる形）
    y, z, names, aoi, unit, g, sels = _pairs_arrays(tk)
    j = names.index("wetside:le500m")
    _threads1()
    t0 = time.time()
    with ProcessPoolExecutor(
        N_WORKERS, initializer=_pairs_init_w, initargs=(y, z, aoi, unit, g)
    ) as ex:
        res = list(ex.map(_pairs_rep_w, sels, chunksize=4))
    b = np.array([r[j] if r is not None else np.nan for r in res])
    ok = np.isfinite(b) & np.isfinite(ref)
    out = {
        "n_reps": len(b),
        "n_both_finite": int(ok.sum()),
        "max_abs_diff_log_or": float(np.abs(b[ok] - ref[ok]).max()),
        "same_nan_pattern": bool(np.array_equal(np.isfinite(b), np.isfinite(ref))),
        "seconds": time.time() - t0,
        "n_cells": len(y),
    }
    return out


# %% B6: Nadeau–Bengio
def b6_nadeau_bengio() -> dict:
    fs = _json(DOCV / "v47_b5_cluster.json")["fold_seeds"]["fixed_code"]
    r = np.array(list(fs["R"].values()), dtype=np.float64)
    j = len(r)
    s2 = float(np.var(r, ddof=1))
    var = (1 / j + 1 / 9) * s2
    tq = float(stats.t.ppf(0.95, j - 1))
    m = float(r.mean())
    return {
        "J": j,
        "R_mean": m,
        "S": math.sqrt(s2),
        "se_corrected": math.sqrt(var),
        "ci90": [m - tq * math.sqrt(var), m + tq * math.sqrt(var)],
        "t_quantile": tq,
        "runs": fs["runs"],
        "rule": "Nadeau and Bengio (2003): (1/J + n_test/n_train) S^2, n_test/n_train = 1/9",
        "note": "captures variation from refitting under other fold assignments, not the sampling of test cells",
    }


# %% B7: 組み合わせた区間
B7_SRC = RUNS_DIR / "world_v2_b8stage2_local" / "predictions" / "primary"
B7_CODES = VER / "v68_b7_codes.npz"
B7_WORKERS = int(
    os.environ.get("V68_B7_WORKERS", "3")
)  # 1 プロセスが 10 fold の設計（約 5 GB）を持つ


def b7_combined() -> dict:
    tab = pd.read_parquet(
        VER / "v41_table.parquet", columns=["aoi", "unit", "fold", "y"]
    )
    lf = v36.load("world_v2_lexfix", "primary")
    if not (
        np.array_equal(np.asarray(lf["aoi"]).astype(str), tab["aoi"].to_numpy(str))
        and np.array_equal(np.asarray(lf["y"]).astype(int), tab["y"].to_numpy(int))
    ):
        return {"note": "lexfix predictions are not aligned with the v41 table"}
    tab["cc"] = np.asarray(lf["cc"]).astype(str)
    tab["act"] = tab["aoi"].str.split("_").str[0]
    acts, a_code = np.unique(tab["act"].to_numpy(), return_inverse=True)
    ctys, c_code = np.unique(tab["cc"].to_numpy(), return_inverse=True)
    units, u_code = np.unique(tab["unit"].to_numpy(str), return_inverse=True)
    tab["a_code"], tab["c_code"], tab["u_code"] = a_code, c_code, u_code
    order = tab.loc[
        tab.sort_values("aoi", kind="stable").index
    ]  # world_eval の表の順（学習側の並び）
    codes: dict[str, np.ndarray] = {}
    for k in range(10):
        z = np.load(B7_SRC / f"fold{k}_stage2.npz", allow_pickle=False)
        tr, te = order[order["fold"] != k], tab[tab["fold"] == k]
        ok = (
            np.array_equal(z["unit_tr"].astype(str), tr["unit"].to_numpy(str))
            and np.array_equal(z["y_tr"].astype(int), tr["y"].to_numpy(int))
            and np.array_equal(z["unit_te"].astype(str), te["unit"].to_numpy(str))
            and np.array_equal(z["y_te"].astype(int), te["y"].to_numpy(int))
        )
        if not ok:
            return {"note": f"fold {k}: stage-2 dump is not aligned with the table"}
        xn_tr, xn_te = (
            np.column_stack([z["xc_tr"], z["xn_tr"]]),
            np.column_stack([z["xc_te"], z["xn_te"]]),
        )
        _, (cc_, ic) = stage2_weighted(
            z["z_tr"], z["z_te"], z["xc_tr"], z["xc_te"], z["y_tr"]
        )
        _, (cn_, in_) = stage2_weighted(z["z_tr"], z["z_te"], xn_tr, xn_te, z["y_tr"])
        for nm, arr in (
            ("a_tr", tr["a_code"]),
            ("c_tr", tr["c_code"]),
            ("u_tr", tr["u_code"]),
            ("a_te", te["a_code"]),
            ("c_te", te["c_code"]),
            ("u_te", te["u_code"]),
        ):
            codes[f"{nm}{k}"] = arr.to_numpy(np.int32)
        codes[f"initc_coef{k}"], codes[f"initc_b{k}"] = cc_, np.array([ic])
        codes[f"initn_coef{k}"], codes[f"initn_b{k}"] = cn_, np.array([in_])
    np.savez(B7_CODES, na=len(acts), nc=len(ctys), nu=len(units), **codes)
    jobs = [(scheme, b) for scheme in ("two_way", "groups") for b in range(B_COMB)]
    _threads1()
    t0 = time.time()
    with ProcessPoolExecutor(B7_WORKERS, initializer=_b7_init) as ex:
        res = list(ex.map(_b7_rep, jobs, chunksize=4))
    print(f"  B7: {len(jobs)} replicates ({time.time() - t0:.0f} s)", flush=True)
    out: dict = {
        "n_activations": len(acts),
        "n_countries": len(ctys),
        "n_units": len(units),
        "replicates": B_COMB,
        "seed": SEED,
    }
    for scheme in ("two_way", "groups"):
        rr = np.array(
            [r for (s_, _), r in zip(jobs, res, strict=True) if s_ == scheme],
            dtype=np.float64,
        )
        r_ref, r_fix = rr[:, 0], rr[:, 1]
        ok = np.isfinite(r_ref)
        out[scheme] = {
            "ci90_R_refit_ratio": np.quantile(r_ref[ok], [0.05, 0.95]).tolist(),
            "ci95_R_refit_ratio": np.quantile(r_ref[ok], [0.025, 0.975]).tolist(),
            "ci90_R_fixed_ratio": np.quantile(r_fix[ok], [0.05, 0.95]).tolist(),
            "sd_refit": float(np.std(r_ref[ok], ddof=1)),
            "sd_fixed": float(np.std(r_fix[ok], ddof=1)),
            "reps_ok": int(ok.sum()),
        }
        lo, hi = out[scheme]["ci90_R_refit_ratio"]
        out[scheme]["within_margin"] = bool(-SESOI_R < lo and hi < SESOI_R)
    return out


_B7: dict = {}


def _b7_init() -> None:
    _threads1()
    c = np.load(B7_CODES)
    _B7.update(na=int(c["na"]), nc=int(c["nc"]), nu=int(c["nu"]), folds=[])
    for k in range(10):
        z = np.load(B7_SRC / f"fold{k}_stage2.npz", allow_pickle=False)
        d = {
            kk: z[kk]
            for kk in (
                "z_tr",
                "z_te",
                "xc_tr",
                "xc_te",
                "y_tr",
                "y_te",
                "p_M0_te",
                "p_M1_te",
                "p_M1c_te",
                "p_M3n_te",
            )
        }
        d["xn_tr"], d["xn_te"] = (
            np.column_stack([z["xc_tr"], z["xn_tr"]]),
            np.column_stack([z["xc_te"], z["xn_te"]]),
        )
        for nm in ("a_tr", "c_tr", "u_tr", "a_te", "c_te", "u_te"):
            d[nm] = c[f"{nm}{k}"]
        d["init_c"] = (c[f"initc_coef{k}"], float(c[f"initc_b{k}"][0]))
        d["init_n"] = (c[f"initn_coef{k}"], float(c[f"initn_b{k}"][0]))
        _B7["folds"].append(d)


def _b7_rep(job: tuple[str, int]) -> tuple[float, float]:
    scheme, b = job
    rng = np.random.default_rng(SEED + b + (0 if scheme == "two_way" else 100_000))
    if scheme == "two_way":
        wa = rng.multinomial(_B7["na"], np.full(_B7["na"], 1 / _B7["na"])).astype(
            np.float64
        )
        wc = rng.multinomial(_B7["nc"], np.full(_B7["nc"], 1 / _B7["nc"])).astype(
            np.float64
        )
    else:
        wu = rng.multinomial(_B7["nu"], np.full(_B7["nu"], 1 / _B7["nu"])).astype(
            np.float64
        )
    acc = np.zeros(6)  # refit_c, refit_n, fixed_c, fixed_n, m0, m1
    for f in _B7["folds"]:
        if scheme == "two_way":
            w_tr, w_te = wa[f["a_tr"]] * wc[f["c_tr"]], wa[f["a_te"]] * wc[f["c_te"]]
        else:
            w_tr, w_te = wu[f["u_tr"]], wu[f["u_te"]]
        if w_te.sum() == 0 or w_tr.sum() == 0:
            continue
        pc, _ = stage2_weighted(
            f["z_tr"],
            f["z_te"],
            f["xc_tr"],
            f["xc_te"],
            f["y_tr"],
            w_tr=w_tr,
            init=f["init_c"],
        )
        pn, _ = stage2_weighted(
            f["z_tr"],
            f["z_te"],
            f["xn_tr"],
            f["xn_te"],
            f["y_tr"],
            w_tr=w_tr,
            init=f["init_n"],
        )
        y = f["y_te"].astype(np.float64)
        for i, p in enumerate(
            (pc, pn, f["p_M1c_te"], f["p_M3n_te"], f["p_M0_te"], f["p_M1_te"])
        ):
            acc[i] += float((w_te * (p - y) ** 2).sum())
    den = acc[4] - acc[5]
    if den <= 0:
        return float("nan"), float("nan")
    return float((acc[0] - acc[1]) / den), float((acc[2] - acc[3]) / den)


# %% B8: 探索的な推定の目録
AMEND = [
    (36, 40, "2"),
    (41, 41, "3"),
    (42, 44, "5"),
    (45, 45, "3, add. 2"),
    (46, 46, "4"),
    (47, 47, "5"),
    (48, 54, "6"),
    (55, 56, "7"),
    (57, 62, "8"),
    (63, 65, "9"),
    (66, 67, "10"),
]
NOT_NAME = (
    "yardstick",
    "terrain",
    "oracle",
    "gsw_oracle",
    "recur",
    "hazard",
    "brier",
    "calib",
    "lab_oracle",
    "label_oracle",
    "positive_control_gsw",
    "jrc",
    "skill",
    "csi",
    "auc",
    "kappa",
    "agree",
    "mc95",
    "bound",
    "margin",
    "plasmode",
    "power",
    "gain",
    "coverage",
    "mde",
    "null_",
    "ceiling",
    "denominator",
)


def _amend(v: int) -> str:
    return next((a for lo, hi, a in AMEND if lo <= v <= hi), "?")


def _is_near(path: str) -> bool:
    """500 m 以内の推定か（道の要素が le500m・or_le500m…、または :le500m で終わる。500m_1km の帯や near_minus_far は外す）。"""
    parts = path.lower().split("/")
    return (
        any(
            q == "le500m"
            or q.endswith(":le500m")
            or q.startswith(("or_le500m", "le500m_"))
            or q == "near"
            for q in parts
        )
        and "near_minus_far" not in path.lower()
    )


def _walk(node, path: list[str], rows: list[dict], file: str, amend: str) -> None:
    if isinstance(node, dict):
        est = kind = ci = None
        if isinstance(node.get("odds_ratio"), (int, float)):
            for k in ("pairs_ci95", "ci95_odds_ratio", "ci95"):
                if isinstance(node.get(k), list) and len(node[k]) == 2:
                    est, kind, ci = float(node["odds_ratio"]), "odds_ratio", node[k]
                    break
        elif isinstance(node.get("R"), (int, float)) and isinstance(
            node.get("ci95_R") or node.get("ci95"), list
        ):
            est, kind, ci = (
                float(node["R"]),
                "R",
                node.get("ci95_R") or node.get("ci95"),
            )
        elif isinstance(node.get("estimate"), (int, float)) and isinstance(
            node.get("ci95"), list
        ):
            est, kind, ci = float(node["estimate"]), "difference", node["ci95"]
        if (
            est is not None
            and ci is not None
            and all(isinstance(x, (int, float)) for x in ci)
        ):
            p = "/".join(path)
            low = p.lower()
            name_effect = not any(tok in low for tok in NOT_NAME)
            null = 1.0 if kind == "odds_ratio" else 0.0
            excl = bool(ci[0] > null or ci[1] < null)
            pv = None
            if kind == "odds_ratio" and est > 0 and ci[0] > 0:
                se = (math.log(ci[1]) - math.log(ci[0])) / (2 * 1.96)
                pv = 2 * stats.norm.sf(abs(math.log(est)) / se) if se > 0 else None
            elif kind != "odds_ratio":
                se = (ci[1] - ci[0]) / (2 * 1.96)
                pv = 2 * stats.norm.sf(abs(est - null) / se) if se > 0 else None
            rows.append(
                {
                    "file": file,
                    "path": p,
                    "kind": kind,
                    "estimate": est,
                    "ci95_lo": float(ci[0]),
                    "ci95_hi": float(ci[1]),
                    "excludes_null": excl,
                    "p_approx": pv,
                    "amendment": amend,
                    "name_effect": name_effect,
                    "near_settlement": _is_near(p) and kind == "odds_ratio",
                }
            )
        for k, v in node.items():
            _walk(v, [*path, str(k)], rows, file, amend)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            if isinstance(v, (dict, list)):
                _walk(v, [*path, str(i)], rows, file, amend)


NUISANCE_LAST = {"logit_M1", "intercept", "const"}
SITE_LAST = {
    "low_hand",
    "low_half",
    "log_river",
    "hand_mean",
    "hand_pct",
    "occ",
    "occ_np",
    "seas",
    "seas_np",
    "elev_mean",
    "slope_mean",
    "hand_min",
}


def _family(path: str, kind: str) -> str:
    """推定の族: nuisance（距離帯・logit の係数）、site（ラベルを使わない地形・水の記述）、flood（浸水のラベル・予測）。"""
    parts = path.split("/")
    last = parts[-1]
    if last in NUISANCE_LAST or last.startswith("band_"):
        return "nuisance"
    if kind == "difference" and (
        last in SITE_LAST
        or last.startswith(("pc_hand", "site_hand", "all_hand", "le500m_hand"))
    ):
        return "site"
    if "gsw_positive" in path or ("positive_control/" in path and last in SITE_LAST):
        return "site"
    if "near_minus_far" in path:
        return "nuisance"
    return "flood"


def b8_census() -> dict:
    rows: list[dict] = []
    for f in sorted(DOCV.glob("v[3-6][0-9]*_*.json")):
        m = re.match(r"v(\d+)", f.name)
        v = int(m.group(1)) if m else 0
        if not 36 <= v <= 67:
            continue
        _walk(_json(f), [], rows, f.name, _amend(v))
    df = pd.DataFrame(rows)
    df["key"] = (
        df["kind"]
        + "|"
        + df["estimate"].round(6).astype(str)
        + "|"
        + df["ci95_lo"].round(6).astype(str)
        + "|"
        + df["ci95_hi"].round(6).astype(str)
    )
    df["duplicate"] = df.duplicated("key")
    # 分離（OR が 0 か無限大に向かう）: |log OR| > 10、または区間の比が 1e6 を超える
    sep = (df["kind"] == "odds_ratio") & (
        (np.abs(np.log(df["estimate"].clip(lower=1e-300))) > 10)
        | (df["ci95_hi"] / df["ci95_lo"].clip(lower=1e-300) > 1e6)
    )
    df["valid"] = (
        np.isfinite(df["ci95_lo"])
        & np.isfinite(df["ci95_hi"])
        & (df["ci95_lo"] <= df["ci95_hi"])
        & ~sep
    )
    df["family"] = [_family(p, k) for p, k in zip(df["path"], df["kind"], strict=True)]
    df.drop(columns="key").to_csv(DOCV / "v68_estimates.csv", index=False)
    base = df.loc[~df["duplicate"] & df["valid"] & df["name_effect"]]
    out: dict = {
        "n_rows": len(df),
        "n_unique": int((~df["duplicate"]).sum()),
        "n_invalid_interval": int((~df["valid"] & ~df["duplicate"]).sum()),
    }
    for fam in ("flood", "site", "nuisance"):
        u = base.loc[base["family"] == fam]
        pv = u["p_approx"].dropna().to_numpy()
        hist, edges = np.histogram(pv, bins=np.linspace(0, 1, 21))
        out[fam] = {
            "n": len(u),
            "n_excluding_null": int(u["excludes_null"].sum()),
            "share_excluding_null": float(u["excludes_null"].mean())
            if len(u)
            else None,
            "by_kind": {
                k: {"n": len(g), "excl": int(g["excludes_null"].sum())}
                for k, g in u.groupby("kind")
            },
            "by_amendment": {
                k: {"n": len(g), "excl": int(g["excludes_null"].sum())}
                for k, g in u.groupby("amendment")
            },
            "p_hist": {"edges": edges.tolist(), "counts": hist.tolist()},
        }
    out["rules"] = {
        "odds_ratio": "odds_ratio with pairs_ci95, else ci95_odds_ratio, else ci95",
        "R": "R with ci95_R or ci95",
        "difference": "estimate with ci95",
        "not_name_effect_tokens": list(NOT_NAME),
        "duplicates": "same kind, estimate and interval",
        "invalid": "interval not finite or inverted (e.g. an undefined ratio), or an odds ratio "
        "with |log OR| > 10 or an interval ratio above 1e6 (separation)",
        "near_settlement": "odds ratios whose path has le500m (not the 500 m to 1 km band) or near "
        "(not near_minus_far)",
        "families": "nuisance: distance-band and logit coefficients and near-minus-far contrasts; "
        "site: differences in terrain or water descriptors without flood labels; flood: the rest "
        "(odds ratios, R, flooded shares)",
        "refined_after_first_pass": "the first pass counted nuisance coefficients, inverted "
        "intervals, separated fits and the 500 m to 1 km band as name effects near settlements; "
        "these rules were added before any count was used in the paper",
    }
    fl = base.loc[
        (base["family"] == "flood") & base["near_settlement"] & base["p_approx"].notna()
    ]
    if len(fl):
        p = fl["p_approx"].to_numpy()
        o = np.argsort(p)
        q = np.empty_like(p)
        q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[
            ::-1
        ]
        out["near_settlement_or"] = {
            "n": len(p),
            "n_excluding_one": int(fl["excludes_null"].sum()),
            "min_bh_q": float(q.min()),
            "n_bh_q_below_0_05": int((q < 0.05).sum()),
            "n_bh_q_below_0_10": int((q < 0.10).sum()),
        }
    return out


# %% B9: 地物の符号と丸めた座標
def b9_gazetteer() -> dict:
    t = v55.table()
    gn = (
        pd.read_parquet(GN, columns=["gid", "fcode", "lat", "lon", "cc"])
        .drop_duplicates("gid")
        .set_index("gid")
    )
    gids = t["P_gid1"].to_numpy(np.int64)
    fc = gn["fcode"].reindex(gids).fillna("?").to_numpy(str)
    ws = t["lex_wetside_n1"].to_numpy() > 0
    out: dict = {"by_fcode_cells": {}, "by_fcode_settlements": {}}
    for code in sorted(set(fc)):
        m = fc == code
        out["by_fcode_cells"][code] = {
            "cells": int(m.sum()),
            "wetside_cells": int((m & ws).sum()),
        }
    us = pd.DataFrame({"gid": gids, "fc": fc, "ws": ws}).drop_duplicates("gid")
    for code, g in us.groupby("fc"):
        out["by_fcode_settlements"][code] = {
            "settlements": len(g),
            "wetside": int(g["ws"].sum()),
        }
    s = gn.loc[np.unique(gids[np.isin(gids, gn.index.to_numpy())])]
    # 整数の分（v37・v48 と同じ許容 2e-3 分。GeoNames の座標は度の小数 5 桁なので、51.98333 のような分の値は
    # 60 倍しても整数から 2e-4 ずれる。最初に書いた 1e-6 ではこれを数え落とした）
    tol = 2e-3
    rnd = (np.abs(s["lat"] * 60 - np.round(s["lat"] * 60)) < tol) & (
        np.abs(s["lon"] * 60 - np.round(s["lon"] * 60)) < tol
    )
    out["rounded_share_settlements"] = float(rnd.mean())
    out["rounded_by_country"] = {
        c: {"n": len(g), "share": float(rnd.loc[g.index].mean())}
        for c, g in s.groupby("cc")
        if len(g) >= 20
    }
    return out


# %% B10: パイロット
def b10_pilot() -> dict:
    out = {}
    for rep in sorted((RUNS_DIR / "pilot3_d000" / "sim").glob("rep*.npz")):
        z = np.load(rep, allow_pickle=False)
        y = z["y"].astype(np.float64)
        bs = {
            m: float(np.mean((z[m].astype(np.float64) - y) ** 2))
            for m in ("M0", "M1", "M1s", "M3s", "M3ps")
        }
        yard = bs["M0"] - bs["M1"]
        out[rep.stem] = {
            "brier": bs,
            "delta_bs_names": bs["M1s"] - bs["M3s"],
            "delta_bs_placebo": bs["M1s"] - bs["M3ps"],
            "R_names": (bs["M1s"] - bs["M3s"]) / yard,
            "R_placebo": (bs["M1s"] - bs["M3ps"]) / yard,
        }
    out["rule"] = (
        "pilot3 at delta = 0 (synthetic labels calibrated to real AOI flooded fractions), nested M1s -> M3s"
    )
    return out


# %% B11: 照合の標本
def b11_audit() -> dict:
    t = v55.table()
    lex = lexicon()
    gn = (
        pd.read_parquet(GN, columns=["gid", "name", "asciiname", "cc", "lat", "lon"])
        .drop_duplicates("gid")
        .set_index("gid")
    )
    gids = t["P_gid1"].to_numpy(np.int64)
    ws = t["lex_wetside_n1"].to_numpy() > 0
    fl = ws & (t["P_d1"].to_numpy() <= NEAR_M) & (t["y"].to_numpy() == 1)
    per = (
        pd.DataFrame({"gid": gids[ws], "fl": fl[ws]})
        .groupby("gid")
        .agg(cells=("fl", "size"), flooded=("fl", "sum"))
    )
    rows = []
    for gid, r in per.iterrows():
        if gid not in gn.index:
            continue
        name, cc = str(gn.at[gid, "asciiname"]), str(gn.at[gid, "cc"])
        toks = [x for x in _SPLIT.split(name.lower()) if x]
        forms = sorted(
            {
                e.form
                for e in lex.entries_for(cc)
                if e.cat in WET_SIDE
                for tok in toks
                if _hit(tok, e)
            }
        )
        for f in forms or ["?"]:
            rows.append(
                {
                    "gid": int(gid),
                    "cc": cc,
                    "form": f,
                    "cells": int(r["cells"]),
                    "flooded": int(r["flooded"]),
                }
            )
    d = pd.DataFrame(rows)
    st = d.groupby(["form", "cc"]).agg(
        cells=("cells", "sum"), flooded=("flooded", "sum"), n=("gid", "nunique")
    )
    st["w"] = st["cells"] / st["cells"].sum() + (
        st["flooded"] / st["flooded"].sum() if st["flooded"].sum() else 0
    )
    st = st.sort_values("w", ascending=False)
    rng = np.random.default_rng(SEED)
    picked, key = [], []
    for (form, cc), _ in st.iterrows():
        cand = d.loc[(d["form"] == form) & (d["cc"] == cc), "gid"].unique()
        cand = [g for g in cand if g not in {p["gid"] for p in picked}]
        take = (
            rng.choice(cand, size=min(AUDIT_PER_STRATUM, len(cand)), replace=False)
            if cand
            else []
        )
        for g in take:
            g = int(g)
            picked.append(
                {
                    "gid": g,
                    "name": str(gn.at[g, "name"]),
                    "asciiname": str(gn.at[g, "asciiname"]),
                    "cc": cc,
                    "lat": float(gn.at[g, "lat"]),
                    "lon": float(gn.at[g, "lon"]),
                }
            )
            key.append(
                {
                    "gid": g,
                    "stratum_form": form,
                    "stratum_cc": cc,
                    "stratum_weight": float(st.loc[(form, cc), "w"]),
                }
            )
            if len(picked) >= AUDIT_MAX:
                break
        if len(picked) >= AUDIT_MAX:
            break
    order = rng.permutation(len(picked))  # 答えの順が層の重みを示さないように並べ替える
    picked = [picked[i] for i in order]
    key = [key[i] for i in order]
    with (DOCV / "v68_audit_sample.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(picked[0]))
        w.writeheader()
        w.writerows(picked)
    with (VER / "v68_audit_key.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fh:  # 答え（公開しない）
        w = csv.DictWriter(fh, fieldnames=list(key[0]))
        w.writeheader()
        w.writerows(key)
    return {
        "n_names": len(picked),
        "n_strata_used": len({(k["stratum_form"], k["stratum_cc"]) for k in key}),
        "n_strata_total": len(st),
        "countries": sorted({p["cc"] for p in picked}),
        "share_flooded_le500m_covered": float(
            st.loc[
                list({(k["stratum_form"], k["stratum_cc"]) for k in key}), "flooded"
            ].sum()
            / max(st["flooded"].sum(), 1)
        ),
        "share_wetside_cells_covered": float(
            st.loc[
                list({(k["stratum_form"], k["stratum_cc"]) for k in key}), "cells"
            ].sum()
            / st["cells"].sum()
        ),
        "rule": "strata = matched form x country, weight = share of water-side cells + share of flooded cells "
        "within 500 m; up to 3 names per stratum at random (seed), up to 150 names",
    }


# %% 実行
def main() -> None:
    global N_WORKERS
    only = (
        set(sys.argv[sys.argv.index("--only") + 1].split(","))
        if "--only" in sys.argv
        else set()
    )
    if "--workers" in sys.argv:  # クラスタの script の段は環境変数を渡せないので引数で
        N_WORKERS = int(sys.argv[sys.argv.index("--workers") + 1])
    res = _json(OUT) if OUT.is_file() else {}
    res["rule"] = (
        "amendment B-11 items B1-B11 (local; post hoc; confirmatory verdict unchanged)"
    )
    res["seed"] = SEED

    def want(k: str) -> bool:
        return not only or k in only

    for k, fn in (
        ("B1", b1_activation_change),
        ("B3", b3_estimator_scale),
        ("B4", b4_gfm_domain),
        ("B6", b6_nadeau_bengio),
        ("B9", b9_gazetteer),
        ("B10", b10_pilot),
        ("B11", b11_audit),
        ("B8", b8_census),
        ("B5", b5_denominator),
        ("B5u", b5_uncapped),
        ("CHKW", check_weighted),
        ("A4dig", a4_digraphs),
        ("C1", c1_country_perm),
        ("B2", b2_uncapped),
        ("B7", b7_combined),
    ):
        if not want(k):
            continue
        t0 = time.time()
        res[k] = fn()
        print(
            k,
            f"({time.time() - t0:.0f} s)",
            json.dumps(res[k], default=float)[:600],
            flush=True,
        )
        cur = (
            _json(OUT) if OUT.is_file() else {}
        )  # 並行して回した別の部分を消さないよう、書く直前に読み直す
        cur.update({kk: vv for kk, vv in res.items() if kk in ("rule", "seed", k)})
        OUT.write_text(json.dumps(cur, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
