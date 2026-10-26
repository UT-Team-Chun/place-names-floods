# %% paper_2 の派生データを Zenodo に預けるときの台帳（ファイル一覧・大きさ・sha256）を作る
"""手元で計算するだけで、どこにも送らない。出力: submission/zenodo_data_manifest.json

預ける対象（submission/README.md §3）:
- 確証解析と感度の run（集計 JSON・予測 npz・fold ごとの来歴 JSON）
- 検出力シミュレーションの集計（power_summary.json のみ。反復ごとの予測は大きいので除く）
- ラベル単独の発生表、N-read の入出力、EMS の製品台帳
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
PATTERNS = [
    "data/runs/world_v2/*.json",
    "data/runs/world_v2/predictions/*/*",
    "data/runs/world_v2_s*/*.json",
    "data/runs/world_v2_s*/predictions/*/*",
    "data/runs/power_v3_null/power_summary.json",
    "data/runs/power_v3_sesoi/power_summary.json",
    "data/global/incidence.json",
    "data/global/nread/*",
    "docs/data/ems_inventory_20261025.json",
    "docs/verification/v34_channel_agreement.json",
    "docs/verification/v35_name_class_contrast.json",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    files = sorted({p for pat in PATTERNS for p in ROOT.glob(pat) if p.is_file()})
    rows = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in files]
    total = sum(r["bytes"] for r in rows)
    out = Path(__file__).with_name("zenodo_data_manifest.json")
    out.write_text(json.dumps({"n_files": len(rows), "total_bytes": total, "files": rows},
                              indent=1), encoding="utf-8")
    print(f"{len(rows)} files, {total / 2**30:.2f} GiB -> {out.name}")


if __name__ == "__main__":
    main()
