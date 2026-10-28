# %% v60（事後・記述、初見の査読 第 4 回の編集者の求め）: 水系側の名の集落の近くの浸水セルが、どの国・事象にあるか
"""確証判定は変えない。集落の近くのオッズ比（v55 の odds と同じ表・同じ AOI の採否）で、最寄りの集落の名が水系側で
500 m 以内のセルのうち浸水したもの（EMS）を、最寄りの集落の国と、事象（activation）ごとに数える。記述のみ。

出力: stdout と docs/verification/v60_where_near.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

from topsem.config import DATA_DIR

OUT = ROOT / "docs" / "verification" / "v60_where_near.json"
LABELS = DATA_DIR / "global" / "labels"
NEAR_M = 500.0
TOP = 15


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def summarise(t: pd.DataFrame, key: str) -> list[dict]:
    g = t.groupby(key).agg(cells_near=("near", "sum"), flooded_near=("flooded_near", "sum"),
                           aois=("aoi", "nunique")).reset_index()
    g = g.loc[g["cells_near"] > 0].sort_values(["flooded_near", "cells_near"], ascending=False)
    tot = float(g["flooded_near"].sum())
    return [{key: str(r[key]), "flooded_near": int(r["flooded_near"]), "cells_near": int(r["cells_near"]),
             "share_of_flooded_near": float(r["flooded_near"] / tot) if tot else 0.0, "aois": int(r["aois"])}
            for _, r in g.iterrows()]


def main() -> None:
    v55 = _load("v55_b7_local")
    t = v55.table()
    rate = t.groupby("aoi")["y"].transform("mean").to_numpy()
    t = t.loc[(rate > 0) & (rate < 1)].reset_index(drop=True)  # odds() と同じ AOI の採否
    t["near"] = (t["lex_wetside_n1"].to_numpy() > 0) & (t["P_d1"].to_numpy() <= NEAR_M)
    t["flooded_near"] = t["near"] & (t["y"].to_numpy() == 1)
    names = {}
    for a in t["aoi"].unique():
        act = (json.loads((LABELS / a / "meta.json").read_text(encoding="utf-8")).get("activation") or {})
        names[str(a).split("_")[0]] = {"name": act.get("name", ""), "year": str(act.get("activation_time", ""))[:4],
                                       "countries": act.get("countries") or []}
    by_cc = summarise(t, "cc")
    by_act = summarise(t, "activation")
    for r in by_act:
        r |= names.get(r["activation"], {})
    tot = int(t["flooded_near"].sum())
    top5 = sum(r["flooded_near"] for r in by_cc[:5])
    out = {"rule": "descriptive (editor, fourth round of LLM-agent critiques); same table and AOIs as the odds ratios",
           "near_m": NEAR_M, "flooded_near_total": tot, "cells_near_total": int(t["near"].sum()),
           "n_countries_with_flooded_near": sum(r["flooded_near"] > 0 for r in by_cc),
           "n_activations_with_flooded_near": sum(r["flooded_near"] > 0 for r in by_act),
           "share_top5_countries": top5 / tot if tot else 0.0,
           "by_country": by_cc, "by_activation": by_act[:TOP],
           "largest_activation_share": by_act[0]["share_of_flooded_near"] if by_act else 0.0}
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"flooded near water-side names: {tot} cells in {out['n_countries_with_flooded_near']} countries, "
          f"{out['n_activations_with_flooded_near']} activations; top 5 countries {100 * out['share_top5_countries']:.0f}%")
    for r in by_cc[:10]:
        print(f"  {r['cc']}: {r['flooded_near']} of {r['cells_near']} ({100 * r['share_of_flooded_near']:.0f}%)")
    for r in by_act[:8]:
        print(f"  {r['activation']} {r.get('name', '')[:50]} {r.get('year', '')}: {r['flooded_near']}")
    assert np.isfinite(out["share_top5_countries"])


if __name__ == "__main__":
    main()
