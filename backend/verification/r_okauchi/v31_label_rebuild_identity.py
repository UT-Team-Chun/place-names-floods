# %% v31: 第 0 世代対応後のラベル全再構築で、既存 AOI のラベルが変わっていないかの照合
"""ems.py に第 0 世代（2012–2013）の書式と属性名の大文字小文字の揺れへの対応を足したうえで、
全 AOI のラベルを作り直した（world_labels --force）。既存 957 AOI について、再構築の前
（イメージ a94ee9a、`_digest_20261025_before_gen0.json`）と後（meta.json の cells_sha256）で
cells.parquet の中身の指紋と格子を比べる。

- 指紋が同じ: ラベルは不変。特徴（features.parquet）もそのまま使える
- 指紋が違い格子が同じ: 浸水被覆などが変わった（属性名の揺れの是正など）。特徴は使えるが
  ラベル変更として記録する
- 格子が違う: AOI 層の選び方が変わった。特徴を作り直す必要がある

ラベル単独の照合で、特徴との関係は見ない（凍結前に許される工程）。クラスタ（NAS）で走らせる:
    cd infra/utens && ./world_submit.py --stage script \
        --extra "verification/r_okauchi/v31_label_rebuild_identity.py" --image ...
出力: stdout の表と data/global/labels/_v31_rebuild_identity.json
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "packages" / "topsem"))

from topsem.config import DATA_DIR

LABELS = DATA_DIR / "global" / "labels"
BEFORE = LABELS / "_digest_20261025_before_gen0.json"


# %%
def main() -> None:
    before = json.loads(BEFORE.read_text(encoding="utf-8"))["aois"]
    after: dict[str, dict] = {}
    for d in sorted(LABELS.glob("EMSR*")):
        if (d / "cells.parquet").is_file():
            m = json.loads((d / "meta.json").read_text(encoding="utf-8"))
            after[d.name] = {"cells_sha256": m.get("cells_sha256"), "grid": m["grid"],
                             "generation": m.get("schema_generation")}
    rows = {}
    for tag in sorted(set(before) | set(after)):
        b, a = before.get(tag), after.get(tag)
        if b is None:
            kind = "new"
        elif a is None:
            kind = "missing"
        elif a["grid"] != b["grid"]:
            kind = "grid_changed"
        elif a["cells_sha256"] != b["cells_sha256"]:
            kind = "cells_changed"
        else:
            kind = "identical"
        rows[tag] = {"kind": kind, "generation": (a or {}).get("generation")}
    count = Counter(r["kind"] for r in rows.values())
    by_gen = Counter((r["kind"], r["generation"]) for r in rows.values())
    out = {"made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "before": str(BEFORE.name), "counts": dict(count),
           "by_kind_generation": {f"{k}|gen{g}": v for (k, g), v in sorted(by_gen.items(),
                                                                         key=str)},
           "changed": {t: r for t, r in rows.items()
                       if r["kind"] in ("grid_changed", "cells_changed", "missing")}}
    (LABELS / "_v31_rebuild_identity.json").write_text(json.dumps(out, indent=1),
                                                       encoding="utf-8")
    print(f"{'kind':16s} {'n':>5s}")
    for k, v in sorted(count.items()):
        print(f"{k:16s} {v:5d}")
    for t, r in list(out["changed"].items())[:30]:
        print(f"  {t:24s} {r['kind']:14s} gen{r['generation']}")


if __name__ == "__main__":
    main()
