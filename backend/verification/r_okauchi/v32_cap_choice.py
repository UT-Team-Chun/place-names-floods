# %% v32: 1 AOI あたりの評価セル上限 N_cap の選び方（ラベル単独の発生表だけで比べる）
"""protocol v2 §2.4-4 の草案は「単一 AOI が全セルの 2% を超えない**最大の** 1,000 刻み」だった。
2026-10-25 の発生表（469 AOI・194 場所グループ）ではこれが 265,000 になり、ほぼ上限なし
（1,330 万セル）になる。上限は推定量の重み付け（どの場所が何割効くか）と計算量を同時に
決めるので、候補ごとに次を並べて凍結前に決め直す:

- 総セル数・期待浸水セル数（AOI ごとの浸水率 × 抽出数）
- 場所グループの実効数（Kish: (Σw)² / Σw²、w = グループのセル数）— 主推定量はセルでプール
  するので、重みの偏りが大きいほど実効的に少ないグループで決まる
- 最大の AOI・最大の場所グループが全セルに占める割合

使うのは発生表（AOI ごとの評価セル数・浸水セル数・場所グループ）だけで、特徴は見ない。
入力: data/global/incidence.json（world_incidence.py、NAS から fetch_nfs_path.sh で取得）
出力: stdout の表と docs/verification/v32_cap_choice.json
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

CAPS = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000]


# %%
def table(rows: list[dict], cap: int) -> dict:
    n = [min(r["n_eval"], cap) for r in rows]
    wet = sum(r["n_wet"] * min(1.0, cap / r["n_eval"]) for r in rows)
    grp: dict[int, float] = defaultdict(float)
    for r, k in zip(rows, n, strict=True):
        grp[r["group"]] += k
    w = list(grp.values())
    tot = sum(n)
    return {"cap": cap, "n_cells": tot, "n_wet_expected": round(wet),
            "wet_rate": wet / tot, "kish_groups": sum(w) ** 2 / sum(x * x for x in w),
            "n_groups": len(w), "max_aoi_share": max(n) / tot, "max_group_share": max(w) / tot,
            "n_aoi_capped": sum(r["n_eval"] > cap for r in rows)}


def table_group(rows: list[dict], gcap: int) -> dict:
    """場所グループ単位の上限: グループの評価セルを合わせて最大 gcap 個（AOI へは件数比で配分）。"""
    tot_g: dict[int, int] = defaultdict(int)
    wet_g: dict[int, float] = defaultdict(float)
    for r in rows:
        tot_g[r["group"]] += r["n_eval"]
        wet_g[r["group"]] += r["n_wet"]
    w = [min(v, gcap) for v in tot_g.values()]
    wet = sum(wet_g[g] * min(1.0, gcap / tot_g[g]) for g in tot_g)
    tot = sum(w)
    return {"group_cap": gcap, "n_cells": tot, "n_wet_expected": round(wet),
            "wet_rate": wet / tot, "kish_groups": tot ** 2 / sum(x * x for x in w),
            "n_groups": len(w), "max_group_share": max(w) / tot,
            "n_groups_capped": sum(v > gcap for v in tot_g.values())}


def main() -> None:
    inc = json.loads((ROOT / "data" / "global" / "incidence.json").read_text(encoding="utf-8"))
    rows = [r for r in inc["aois"] if r["included"]]
    sizes = sorted(r["n_eval"] for r in rows)
    q = {p: sizes[int(p / 100 * (len(sizes) - 1))] for p in (10, 25, 50, 75, 90, 99)}
    caps = [*CAPS, inc["summary"].get("draft_rule_N_cap", inc["summary"].get("N_cap"))]
    res = [table(rows, c) for c in caps]
    print(f"AOIs {len(rows)}, n_eval quantiles {q}, max {sizes[-1]:,}")
    print(f"{'cap':>8s} {'cells':>11s} {'wet':>9s} {'wet%':>6s} {'Kish':>7s} {'/grp':>5s} "
          f"{'maxAOI':>7s} {'maxGrp':>7s} {'capped':>6s}")
    for t in res:
        print(f"{t['cap']:8,d} {t['n_cells']:11,d} {t['n_wet_expected']:9,d} "
              f"{100 * t['wet_rate']:6.2f} {t['kish_groups']:7.1f} {t['n_groups']:5d} "
              f"{100 * t['max_aoi_share']:6.2f}% {100 * t['max_group_share']:6.2f}% "
              f"{t['n_aoi_capped']:6d}")
    gres = [table_group(rows, g) for g in (5_000, 10_000, 20_000, 50_000, 100_000)]
    print(f"\n{'gcap':>8s} {'cells':>11s} {'wet':>9s} {'wet%':>6s} {'Kish':>7s} {'/grp':>5s} "
          f"{'maxGrp':>7s} {'capped':>6s}")
    for t in gres:
        print(f"{t['group_cap']:8,d} {t['n_cells']:11,d} {t['n_wet_expected']:9,d} "
              f"{100 * t['wet_rate']:6.2f} {t['kish_groups']:7.1f} {t['n_groups']:5d} "
              f"{100 * t['max_group_share']:6.2f}% {t['n_groups_capped']:6d}")
    out = {"source": "data/global/incidence.json", "n_aoi": len(rows),
           "n_eval_quantiles": q, "n_eval_max": sizes[-1], "table": res, "group_cap": gres}
    dst = ROOT / "docs" / "verification" / "v32_cap_choice.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
