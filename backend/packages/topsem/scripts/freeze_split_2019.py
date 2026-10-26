"""2019 パイロットの LORO 分割を manifest として凍結する。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/freeze_split_2019.py

結果を data/splits/pilot-loro-2019/manifest.json に書く（以後 run はこれを読む）。
"""

from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.config import LABELS_DIR
from topsem.splits.manifest import make_watershed_holdout, save_split

EVENT = "2019_taihuu19gou"
EVENT_DATE = "2019-10-12"


def main() -> None:
    units = []
    for d in sorted((LABELS_DIR / EVENT / "derived").iterdir()):
        meta = json.loads((d / "metadata.json").read_text(encoding="utf-8"))
        units.append((EVENT, meta["river"]))
    spec = make_watershed_holdout(
        units,
        version="pilot-loro-2019",
        event_dates={EVENT: EVENT_DATE},
        notes="GSI 浸水推定段彩図 6 水系の leave-one-river-out（探査的 pilot 用）",
    )
    spec.created_utc = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    path = save_split(spec)
    print(f"{len(spec.folds)} folds -> {path}")


if __name__ == "__main__":
    main()
