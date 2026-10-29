"""v75: 集落の近くの水系側のオッズ比から EMSR052・EMSR130 を除いた版（改訂 B-14 の追記 5）。

第 10 回の初見の批評（査読者 1 の (f)、査読者 3 の (d)）が求めた版。500 m 以内の水系側の浸水セルの半分ほどが
EMSR052（カンボジア）と EMSR130（ミャンマー）にある。

版（追記 5 の 1）:
  (a) 改訂 11 が主とした 2 版: 上限なしの EMS の評価できる全てのセル（B-11 の B2 の ems_all）と、GFM の評価できる
      全てのセル（B2 の gfm）。直したコードの読み（v68 の with_corrected_wetside）
  (b) 改訂 3 で計画した GFM の推定（v41 の表、上限あり、凍結した読み。B-14 の B13 (iv) と同じ当て直し）
それぞれで EMSR052 を除く、EMSR130 を除く、両方を除く。推定量は各版と同じ（AOI の固定効果、距離帯、M1 の logit）、
区間は場所グループでクラスタ頑健な 95% の Wald 区間（pairs bootstrap はしない）。除かない版も同じ区間で並べる。

出力: docs/verification/v75_b14_leaveout.json
"""

# %%
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

OUT = ROOT / "docs" / "verification" / "v75_b14_leaveout.json"
ACTS = ("EMSR052", "EMSR130")
NEAR_M = 500.0


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v70 = _load("v70_b13_local")
v69, v68 = v70.v69, v70.v68
v39, v55 = v69.v39, v68.v55


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _act(aoi: pd.Series) -> np.ndarray:
    return aoi.astype(str).str.split("_").str[0].to_numpy()


def _drops() -> dict[str, tuple[str, ...]]:
    return {"all": (), "without_EMSR052": ("EMSR052",), "without_EMSR130": ("EMSR130",), "without_both": ACTS}


# %% (a) 改訂 11 の主な 2 版（上限なし、直したコードの読み）
def _versions_a(t: pd.DataFrame, ycol: str) -> dict:
    """v68 の pairs_or と同じ表と推定量（v55.odds、場所でクラスタ頑健）。除く activation ごとに当て直す。"""
    s0 = t if ycol == "y" else t.assign(y=t[ycol].to_numpy())
    act = _act(s0["aoi"])
    out: dict = {}
    for key, drop in _drops().items():
        m = ~np.isin(act, drop)
        if drop and m.all():
            out[key] = {"not_applicable": True}
            continue
        t0 = time.time()
        o = v55.odds(s0, cls="wetside", mask=m)
        s = s0.loc[m]
        near = (s["P_d1"].to_numpy() <= NEAR_M) & (s["lex_wetside_n1"].to_numpy() > 0)
        o |= {"n_aoi_removed": int(s0.loc[~m, "aoi"].nunique()),
              "n_flooded_near": int((near & (s["y"].to_numpy() == 1)).sum())}
        out[key] = o
        print(f"    {key}: OR {o['odds_ratio']:.3f} ({o['ci95'][0]:.3f}-{o['ci95'][1]:.3f}), "
              f"{o['n_flooded_near']} flooded near, {time.time() - t0:.0f} s", flush=True)
    return out


def versions_a() -> dict:
    u = v68.uncapped_full()
    u = v68.with_corrected_wetside(u[[c for c in v68.KEEP if c in u.columns]].copy())
    u = u.loc[u["P_gid1"].notna() & u["P_d1"].notna()].reset_index(drop=True)
    print(f"  ems_all ({len(u):,} cells)", flush=True)
    out = {"ems_all": _versions_a(u, "y")}
    g = u.loc[u["gfm_eval"] == 1].reset_index(drop=True)
    del u
    print(f"  gfm ({len(g):,} cells)", flush=True)
    out["gfm"] = _versions_a(g, "gfm_y")
    return out


# %% (b) 改訂 3 で計画した GFM の推定（v41、上限あり、凍結した読み）
def versions_b() -> dict:
    """v73 の _gfm_jackknife と同じ表・設計・当て直し（v39.fit、場所でクラスタ頑健）。"""
    v41 = v70._load("v41_gfm_independent_labels")
    g = v41.attach_gfm(v41.build(), _json(v69.GFM_SUMMARY))
    y = g["y_gfm"].to_numpy(np.float64)
    rate = pd.Series(y).groupby(g["aoi"].to_numpy()).transform("mean").to_numpy()
    s = g.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)
    ys = s["y_gfm"].to_numpy(np.float64)
    z, names = v39.design(s, ("wetside",), True, banded=True)
    j = names.index("wetside:le500m")
    unit, aoi = s["unit"].to_numpy(str), s["aoi"].to_numpy()
    act = _act(s["aoi"])
    near = (s["P_d1"].to_numpy() <= NEAR_M) & (s["lex_wetside_n1"].to_numpy() > 0)
    out: dict = {"n_aoi_by_activation": {a: int(s.loc[act == a, "aoi"].nunique()) for a in ACTS}}
    for key, drop in _drops().items():
        m = ~np.isin(act, drop)
        if drop and m.all():
            out[key] = {"not_applicable": True}
            continue
        r = v39.fit(ys[m], z[m], aoi[m], unit[m])
        b, se = float(r["beta"][j]), float(r["se"][j])
        out[key] = {"odds_ratio": math.exp(b), "ci95": [math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)],
                    "n_cells": int(m.sum()), "n_aoi": int(pd.unique(aoi[m]).size), "n_units": int(r["n_groups"]),
                    "n_aoi_removed": int(pd.unique(aoi[~m]).size),
                    "n_flooded_near": int((near & m & (ys == 1)).sum())}
        print(f"    {key}: OR {out[key]['odds_ratio']:.3f}", flush=True)
    return out


# %%
def main() -> None:
    t0 = time.time()
    print("(b) planned GFM estimate (v41)", flush=True)
    res = {"planned_gfm": versions_b()}
    print("(a) main versions of amendment 11", flush=True)
    res |= versions_a()
    res["rule"] = ("B-14 addendum 5: water-side odds ratio within 500 m of settlements without the AOIs of EMSR052, of "
                   "EMSR130 or of both; estimator of each version (AOI fixed effects, distance bands, logit of M1); 95% "
                   "Wald intervals clustered by location group (no pairs bootstrap); 'all' is the version without "
                   "exclusions with the same interval")
    res["seconds"] = round(time.time() - t0)
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\n{'version':<14}{'removed':<18}{'OR':>7}  {'95% Wald':<16}{'flooded near':>13}")
    for v in ("planned_gfm", "ems_all", "gfm"):
        for k in _drops():
            o = res[v][k]
            if o.get("not_applicable"):
                print(f"{v:<14}{k:<18}{'n/a':>7}")
                continue
            print(f"{v:<14}{k:<18}{o['odds_ratio']:7.3f}  {o['ci95'][0]:.3f}-{o['ci95'][1]:.3f}    "
                  f"{o['n_flooded_near']:>10}")
    print(f"\n-> {OUT.relative_to(ROOT)} ({res['seconds']} s)")


if __name__ == "__main__":
    main()
