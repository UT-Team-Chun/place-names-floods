# %% v45（事後・amendment B-3 追記 2）: 仮説と逆向きの B は GFM の部分的な遮蔽で説明できるか
"""B の結果（水系側の群、500 m 以内、GFM で OR 0.55）を見た後に加えた確認。登録した B の解釈は変えない。

GFM で判定可能な割合 v（撮像の最大）が 0.9 以上のセル、1 のセルに限って、v41 の B と同じ回帰（AOI の切片、
距離帯、M1 のホールドアウト予測のロジット、場所グループでクラスタ）を GFM と EMS のラベルで当てる。

出力: stdout と docs/verification/v45_gfm_masking_check.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import v41_gfm_independent_labels as v41

OUT = v41.ROOT / "docs" / "verification" / "v45_gfm_masking_check.json"
CUTS = (0.9, 1.0)
KEYS = ("wetside_adjusted_M1", "wetside_raw", "wetside_allcells_adjusted_M1")


def attach_v(t: pd.DataFrame, summary: dict) -> np.ndarray:
    """v41.attach_gfm と同じ行順で、GFM で判定可能な割合 v を返す。"""
    keep = sorted(a for a, s in summary["aois"].items() if s.get("include"))
    g = pd.concat([pd.read_parquet(v41.LABELS / a / "gfm.parquet", columns=["ix", "iy", "gfm_v"])
                   .assign(aoi=a) for a in keep], ignore_index=True)
    m = t[["aoi", "ix", "iy"]].merge(g, on=["aoi", "ix", "iy"], how="left", validate="many_to_one")
    return m["gfm_v"].to_numpy(dtype=np.float64)


def main() -> None:
    summary = json.loads(v41.GFM_SUMMARY.read_text(encoding="utf-8"))
    t = v41.attach_gfm(v41.build(), summary)
    v = attach_v(t, summary)
    res: dict = {"rule": "amendment B-3 addendum 2 (post hoc, after seeing B)",
                 "n_cells_all": len(t), "share_v_ge_0.9": float((v >= 0.9).mean()),
                 "share_v_eq_1": float((v >= 1.0).mean())}
    for cut in CUTS:
        s = t.loc[v >= cut].reset_index(drop=True)
        blk = {"n_cells": len(s), "n_aoi": int(s["aoi"].nunique())}
        for lab, col in (("gfm", "y_gfm"), ("ems", "y_ems")):
            o = v41.odds(s, col)
            blk[lab] = {k: o[k]["terms"] for k in KEYS} | {"n_cells": o["n_cells"], "n_aoi": o["n_aoi"]}
        res[f"v_ge_{cut}"] = blk
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for cut in CUTS:
        blk = res[f"v_ge_{cut}"]
        for lab in ("gfm", "ems"):
            near = blk[lab]["wetside_adjusted_M1"]["wetside:le500m"]
            allc = blk[lab]["wetside_allcells_adjusted_M1"]["wetside:all"]
            print(f"v>={cut} {lab}: cells {blk[lab]['n_cells']:,}  <=500m OR {near['odds_ratio']:.2f} "
                  f"[{near['ci95_odds_ratio'][0]:.2f}, {near['ci95_odds_ratio'][1]:.2f}]  all OR "
                  f"{allc['odds_ratio']:.2f} [{allc['ci95_odds_ratio'][0]:.2f}, {allc['ci95_odds_ratio'][1]:.2f}]")


if __name__ == "__main__":
    main()
