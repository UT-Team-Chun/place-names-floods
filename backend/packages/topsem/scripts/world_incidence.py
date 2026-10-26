"""v2 の発生表（ラベル単独）: AOI ごとの評価セル・浸水セル・種別・海岸距離・採用理由・場所グループ。

protocol v2 の凍結前に許される集計だけを行う（特徴は読まない）。N_cap の決定（§2.4-4）と
検出力シミュレーション（§8）の入力、論文の Fig. 1 と採用フロー図の出所になる。

    uv run python packages/topsem/scripts/world_incidence.py
出力: data/global/incidence.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from topsem.config import DATA_DIR
from topsem.world.dataset import (
    GROUP_CAP,
    CoastIndex,
    aoi_info,
    group_quotas,
    include,
    location_groups,
)

ROOT = DATA_DIR / "global"
CAP_SHARE = 0.02  # 草案の規則（参考値として記録のみ。§2.4-4 は場所グループ上限に差し替え）


def n_cap(n_eval: list[int], share: float = CAP_SHARE) -> int:
    """草案の規則: 単一 AOI のセル比が share 以下に収まる最大の 1,000 刻み cap（参考値）。"""
    top = max(n_eval)
    best = 1000
    for cap in range(1000, top + 1000, 1000):
        capped = [min(n, cap) for n in n_eval]
        if max(capped) / sum(capped) <= share:
            best = cap
        else:
            break
    return min(best, top)


def main() -> None:
    gn = pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")
    p_names = gn[gn["fclass"] == "P"]
    coast = CoastIndex(ROOT / "coast" / "ne_10m_coastline.zip")
    rows = []
    for d in sorted(p for p in (ROOT / "labels").glob("EMSR*") if (p / "cells.parquet").is_file()):
        info = aoi_info(d, coast)
        w, s, e, n = info.bounds
        n_names = int((p_names["lon"].between(w, e) & p_names["lat"].between(s, n)).sum())
        ok, why = include(info, n_names_inside=n_names)
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        rows.append({
            "tag": info.tag, "activation": info.activation, "countries": list(info.countries),
            "event_time": info.event_time, "bounds": list(info.bounds), "n_eval": info.n_eval,
            "n_wet": info.n_wet, "riverine_share": info.riverine_share,
            "unspecified_share": info.unspecified_share, "coast_km": info.coast_km,
            "unspecified_as_riverine": info.unspecified_as_riverine,
            "footprint_implicit": info.footprint_implicit, "n_P_inside": n_names,
            "observed_event_layer": meta.get("observed_event_layer", True),
            "included": ok, "reason": why,
        })
    inc = [r for r in rows if r["included"]]
    grp = location_groups({r["tag"]: tuple(r["bounds"]) for r in inc})
    for r in inc:
        r["group"] = grp[r["tag"]]
    quota = group_quotas({r["tag"]: r["n_eval"] for r in inc}, grp) if inc else {}
    for r in inc:
        r["quota"] = quota[r["tag"]]
    summary = {
        "n_aoi_labelled": len(rows),
        "n_aoi_included": len(inc),
        "reasons": dict(Counter(r["reason"] for r in rows)),
        "n_no_event_products": sum(1 for r in rows if not r["observed_event_layer"]),
        "n_activations_included": len({r["activation"] for r in inc}),
        "n_location_groups": len(set(grp.values())),
        "n_unspecified_as_riverine": sum(1 for r in inc if r["unspecified_as_riverine"]),
        "group_cap": GROUP_CAP,
        "n_cells_after_cap": sum(quota.values()),
        "n_wet_after_cap_expected": sum(r["n_wet"] * r["quota"] / r["n_eval"] for r in inc),
        "draft_rule_N_cap": n_cap([r["n_eval"] for r in inc]) if inc else 0,
        "years": sorted({r["event_time"][:4] for r in inc}),
    }
    out = ROOT / "incidence.json"
    out.write_text(json.dumps({"summary": summary, "aois": rows}, ensure_ascii=False, indent=1),
                   encoding="utf-8")
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
