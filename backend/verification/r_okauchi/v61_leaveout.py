# %% v61（事後・amendment B-8 追記 6）: パキスタン（とバングラデシュ）の AOI を除いた集落の近くのオッズ比と、パキスタンの層の記述
"""確証判定は変えない。追記 6 で計算の前に固定した量（全て登録外の記述）:

- (a) パキスタンの AOI を除いた、水系側 × 500 m 以内のオッズ比（v57 の B6 の主と同じ式、クラスタ頑健な区間と
  pairs bootstrap 2,000 回）
- (b) パキスタンとバングラデシュの AOI を除いた同じもの
- (c) パキスタンの層の記述: activation ごとの評価セル・浸水セル・水系側の名の近くのセルと浸水セル

AOI の国は activation の記録の国（meta.json）。出力: stdout と docs/verification/v61_leaveout.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

from topsem.config import DATA_DIR

OUT = ROOT / "docs" / "verification" / "v61_leaveout.json"
LABELS = DATA_DIR / "global" / "labels"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    v57 = _load("v57_b8_local")
    t, _gn, _infos, _dirs = v57.analysis_table()
    countries = {}
    for a in t["aoi"].unique():
        act = json.loads((LABELS / a / "meta.json").read_text(encoding="utf-8")).get("activation") or {}
        countries[a] = {"countries": act.get("countries") or [], "name": act.get("name", ""),
                        "year": str(act.get("activation_time", ""))[:4]}
    in_c = lambda names: t["aoi"].map(lambda a: any(c in countries[a]["countries"] for c in names)).to_numpy()
    pk, pkbd = in_c(("Pakistan",)), in_c(("Pakistan", "Bangladesh"))
    out: dict = {"rule": "amendment B-8 addendum 6 (descriptive, not planned in the protocol)", "seed": v57.SEED,
                 "n_aoi_pakistan": int(t.loc[pk, "aoi"].nunique()),
                 "n_aoi_pakistan_bangladesh": int(t.loc[pkbd, "aoi"].nunique())}
    out["without_pakistan"] = v57._with_pairs(t, "wetside", "nopk", ~pk)
    print("without Pakistan", {k: out["without_pakistan"].get(k) for k in ("odds_ratio", "ci95", "pairs_ci95")}, flush=True)
    out["without_pakistan_bangladesh"] = v57._with_pairs(t, "wetside", "nopkbd", ~pkbd)
    print("without Pakistan and Bangladesh",
          {k: out["without_pakistan_bangladesh"].get(k) for k in ("odds_ratio", "ci95", "pairs_ci95")}, flush=True)
    s = t.loc[pk].copy()
    s["act"] = s["aoi"].str.split("_").str[0]
    s["near"] = (s["lex_wetside_n1"].to_numpy() > 0) & (s["P_d1"].to_numpy() <= 500)
    s["flooded_near"] = s["near"] & (s["y"].to_numpy() == 1)
    rows = []
    for act, g in s.groupby("act"):
        meta = countries[g["aoi"].iloc[0]]
        rows.append({"activation": act, "name": meta["name"], "year": meta["year"], "aois": int(g["aoi"].nunique()),
                     "cells": len(g), "flooded": int(g["y"].sum()), "cells_near": int(g["near"].sum()),
                     "flooded_near": int(g["flooded_near"].sum())})
    rows.sort(key=lambda r: -r["cells"])
    out["pakistan_by_activation"] = rows
    out["pakistan_flooded_share"] = float(s["y"].mean()) if len(s) else None
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    for r in rows:
        print(r)
    assert np.isfinite(out["without_pakistan"]["odds_ratio"])


if __name__ == "__main__":
    main()
