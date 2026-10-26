# %% paper_2 の派生データを Zenodo に預けるときの台帳（ファイル一覧・大きさ・sha256）を作る
"""手元で計算するだけで、どこにも送らない。出力: submission/zenodo_data_manifest.json（公開リポジトリにも載せる）

預ける対象は make_zenodo_data_package.ARCHIVES と同じ（zip ごとの中身を一つの台帳にまとめる）:
セル表、確証・感度・事後（amendment B-2）の run、検出力の集計、記録類（発生表・N-read・GeoNames の部分集合・
EMS 台帳・検証 v31–v40 の出力）。
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from make_zenodo_data_package import ARCHIVES, ROOT


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    rows = []
    for archive, patterns in ARCHIVES.items():
        files = sorted({p for pat in patterns for p in ROOT.glob(pat) if p.is_file()})
        rows += [{"archive": archive, "path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size,
                  "sha256": sha256(p)} for p in files]
    total = sum(r["bytes"] for r in rows)
    out = Path(__file__).with_name("zenodo_data_manifest.json")
    out.write_text(json.dumps({"n_files": len(rows), "total_bytes": total, "files": rows},
                              indent=1), encoding="utf-8")
    print(f"{len(rows)} files, {total / 2**30:.2f} GiB -> {out.name}")


if __name__ == "__main__":
    main()
