# %% v53（事後・amendment B-6 の 3(c)）: 繰り返す水の正の対照（浸水のラベルを使わない）
"""確証判定は変えない。B-6 の 3(c) で計算の前に固定した量を計算する。

集落の点から 500 m 以内の JRC Global Surface Water v1.4 の occurrence（水があった時間の割合、%）と seasonality
（水のあった月の数）の平均（v49、恒常水域を除いた変種も）を、B-2 の正の対照（集落の地点の HAND、v37）と同じ形で
比べる: AOI ごとに、名前が水系側（または区分・厳密な読み）の集落とそれ以外の集落の平均の差を、小さい群の集落数で
重み付けて平均し、場所グループのクラスタ bootstrap（4,000 回、シード 20261028）で区間を付ける。

集落の集合: 解析したセルの最寄りの集落（名前が特徴に入る集落）を AOI ごとに。B-2 と同じ集合（中心が集落から
177 m 以内のセルを持つ集落）でも。読みは B-4 の読み方（直したコード）。浸水のラベルは使わない。

出力: stdout と docs/verification/v53_gsw_positive_control.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.config import DATA_DIR
from topsem.world.b6 import STRICT_CATEGORIES, strict_lexicon
from topsem.world.lexicon import CATEGORIES, WET_SIDE, Lexicon

OUT = ROOT / "docs" / "verification" / "v53_gsw_positive_control.json"
LEX = ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
GSW = DATA_DIR / "global" / "gsw" / "settlement_gsw.parquet"
V37T = DATA_DIR / "verification" / "v37_table.parquet"
GN = DATA_DIR / "global" / "geonames" / "aoi_subset.parquet"
SEED = 20261028
VARS = ("occ", "seas", "occ_np", "seas_np")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v37 = _load("v37_positive_control")
v37.SEED = SEED  # 区間は B-6 のシードで（形は v37 の settlement_contrast のまま）


def settlements() -> tuple[pd.DataFrame, dict]:
    t = pd.read_parquet(V37T, columns=["aoi", "unit", "P_gid1", "P_d1"])
    t = t[t["P_gid1"] >= 0]
    near = t.loc[t["P_d1"] <= v37.POINT_M, ["aoi", "P_gid1"]].drop_duplicates()
    s = t[["aoi", "unit", "P_gid1"]].drop_duplicates().reset_index(drop=True)
    s["point_set"] = s.set_index(["aoi", "P_gid1"]).index.isin(near.set_index(["aoi", "P_gid1"]).index)
    gn = pd.read_parquet(GN, columns=["gid", "asciiname", "cc"]).drop_duplicates("gid").set_index("gid")
    lex = Lexicon(LEX)
    ls = strict_lexicon(lex)
    flags = []
    for g in s["P_gid1"].to_numpy():
        name, cc = str(gn.at[g, "asciiname"]), str(gn.at[g, "cc"])
        m, ms = lex.match(name, cc), ls.match(name, cc)
        flags.append({"wetside": float(any(m[c] for c in WET_SIDE)), **{c: float(m[c]) for c in CATEGORIES},
                      "strict": float(any(ms[c] for c in STRICT_CATEGORIES))})
    s = pd.concat([s, pd.DataFrame(flags)], axis=1)
    w = pd.read_parquet(GSW).set_index("gid")
    s = s.join(w, on="P_gid1")
    info = {"n_settlement_aoi_pairs": len(s), "n_settlements": int(s["P_gid1"].nunique()),
            "share_with_gsw": float(s["n_px"].fillna(0).gt(0).mean()),
            "share_point_set": float(s["point_set"].mean())}
    return s, info


def main() -> None:
    s, info = settlements()
    res: dict = {"rule": "amendment B-6 item 3(c)", "seed": SEED, "B": v37.B, "info": info,
                 "mean_by_group": {}, "all_named": {}, "point_set": {}}
    groups = ("wetside", "strict", *CATEGORIES)
    for g in groups:
        sel = s[g] > 0.5
        res["mean_by_group"][g] = {"n": int(sel.sum()),
                                   **{v: [float(s.loc[sel, v].mean()), float(s.loc[~sel, v].mean())] for v in VARS}}
    for key, sub in (("all_named", s), ("point_set", s[s["point_set"]])):
        for g in groups:
            res[key][g] = {v: v37.settlement_contrast(sub.dropna(subset=[v]), v, g) for v in VARS}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps(info))
    for key in ("all_named", "point_set"):
        for g in groups:
            row = []
            for v in ("occ", "seas", "occ_np"):
                c = res[key][g][v]
                row.append(f"{v} {c['estimate']:+.2f} [{c['ci95'][0]:+.2f}, {c['ci95'][1]:+.2f}]" if c else f"{v} -")
            print(f"{key:10s} {g:9s}", " | ".join(row))


if __name__ == "__main__":
    main()
