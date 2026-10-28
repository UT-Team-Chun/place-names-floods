# %% v65（事後・記述）: 梯子の T0 の段の名前の増分の区間が、凍結したコードと直したコードで 4 割違う理由
"""初見の査読 第 5 回の編集者の指摘 7: 表 2（凍結したコード）の T0 の段の 95% 区間より、補足の表（直したコード）の
90% 区間が広い。直したコードはノルウェーの名に辞書を当てただけ（セルの 2.2%）なのに、なぜ区間がこれほど広がるか。

場所グループごとの寄与 Σ[(p_T0c − y)² − (p_T0n − y)²] を二つの実行で比べ、区間の幅（場所グループの間の散らばり）を
動かした場所グループとその国を挙げる。新しい推定はしない（記述）。

出力: docs/verification/v65_t0_rung.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _root() -> Path:
    p = Path(__file__).resolve()
    while not (p / ".git").exists():
        p = p.parent
    return p


ROOT = _root()
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa  # noqa: E402

from topsem.config import RUNS_DIR  # noqa: E402

OUT = ROOT / "docs" / "verification" / "v65_t0_rung.json"
RUNS = {"frozen": "world_v2", "corrected": "world_v2_lexfix"}
TOP = 5


# %% 場所グループごとの寄与
def unit_contrib(run: str) -> tuple[pd.DataFrame, int]:
    d = wa.load_split(RUNS_DIR / run, "primary")
    y = d["y"].astype(np.float64)
    c = (d["p_T0c"].astype(np.float64) - y) ** 2 - (d["p_T0n"].astype(np.float64) - y) ** 2
    df = pd.DataFrame({"unit": d["unit"].astype(str), "aoi": d["aoi"].astype(str), "c": c})
    u = df.groupby("unit").agg(c=("c", "sum"), n=("c", "size"), aoi=("aoi", "first"))
    return u, len(y)


def countries() -> dict[str, str]:
    t = pd.read_parquet(ROOT / "data" / "verification" / "v37_table.parquet", columns=["aoi", "cc"])
    return t.groupby("aoi")["cc"].agg(lambda s: s.value_counts().index[0]).to_dict()


def main() -> None:
    cc = countries()
    out: dict = {"rule": "descriptive (editor, fifth round of LLM-agent critiques, item 7)", "runs": RUNS}
    parts = {}
    for key, run in RUNS.items():
        u, n = unit_contrib(run)
        parts[key] = u
        lad = json.loads((RUNS_DIR / run / "metrics_primary.json").read_text(encoding="utf-8"))
        yard = float(lad["primary"]["delta_bs_M0_M1"])
        lo, hi = lad["ladder"]["ci95_delta_names"][0]
        out[key] = {"delta_bs_T0": float(u["c"].sum() / n), "ci95_delta_bs_T0": [lo, hi],
                    "ci95_R_T0": [lo / yard, hi / yard], "yardstick": yard,
                    "spread": float((u["c"] / n).std() * np.sqrt(len(u))), "n_units": len(u)}
    m = parts["frozen"].join(parts["corrected"], lsuffix="_frozen", rsuffix="_corrected")
    n_tot = float(m["n_frozen"].sum())
    m["change"] = (m["c_corrected"] - m["c_frozen"]) / n_tot
    top = m.reindex(m["change"].abs().sort_values(ascending=False).index).head(TOP)
    out["top_units"] = [{"unit": i, "aoi": r["aoi_frozen"], "country": cc.get(r["aoi_frozen"]),
                         "contrib_frozen": float(r["c_frozen"] / n_tot), "contrib_corrected": float(r["c_corrected"] / n_tot),
                         "change": float(r["change"])} for i, r in top.iterrows()]
    # 上位 2 つを除いた散らばり（区間の幅の比の目安）
    rest = m.drop(index=top.index[:2])
    out["spread_without_top2"] = {k: float((rest[f"c_{k}"] / n_tot).std() * np.sqrt(len(rest)))
                                  for k in ("frozen", "corrected")}
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
