# %% F10 対応: 全入力ファイルの sha256 manifest（査読提出パッケージ用）。
# ディレクトリー群（タイル・河川 line 群）は Merkle 集約（per-file sha → sorted 連結 → 1 ハッシュ）。
# 使い方（backend/ から）:: uv run python packages/topsem/scripts/make_input_manifest.py
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "data/labels").is_dir())
DATA = REPO / "data"
OUT = REPO / "docs/data/input_manifest.json"


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def merkle(files: list[Path]) -> tuple[str, int]:
    h = hashlib.sha256()
    for f in sorted(files):
        h.update(str(f.relative_to(DATA)).encode() + b"\0" + sha_file(f).encode() + b"\0")
    return h.hexdigest(), len(files)


manifest: dict[str, dict] = {}

# 1) population
for d in sorted((DATA / "population").glob("*")):
    if not d.is_dir():
        continue
    manifest[f"population/{d.name}"] = {
        f.relative_to(d).as_posix(): {"sha256": sha_file(f), "bytes": f.stat().st_size}
        for f in sorted(d.rglob("*")) if f.is_file()
    }

# 2) labels（derived 表は個別ハッシュ・raw/depth_sheets は集約）
for d in sorted((DATA / "labels").glob("*")):
    for sub in ("derived", "raw", "depth_sheets"):
        base = d / sub
        if not base.is_dir():
            continue
        rels = sorted(f for f in base.rglob("*") if f.is_file())
        dg, n = merkle(rels)
        manifest[f"labels/{d.name}/{sub}"] = {
            "merkle_sha256": dg, "n_files": n,
            "bytes": sum(f.stat().st_size for f in rels),
        }

# 3) splits（凍結 manifest は個別）
for m in sorted((DATA / "splits").glob("*/*/manifest.json")):
    manifest[f"splits/{m.parent.relative_to(DATA / 'splits')}"] = {
        "sha256": sha_file(m), "bytes": m.stat().st_size,
    }

# 4) covariates（大容量はディレクトリ単位の merkle）
for d in sorted((DATA / "covariates").glob("*")):
    if not d.is_dir():
        continue
    rels = [f for f in sorted(d.rglob("*")) if f.is_file()]
    dg, n = merkle(rels)
    manifest[f"covariates/{d.name}"] = {"merkle_sha256": dg, "n_files": n,
                                        "bytes": sum(f.stat().st_size for f in rels)}

# 5) corpus snapshot（app リポ側 manifest の写しがあれば同封）
snap = DATA / "corpus_snapshot.json"
if snap.is_file():
    manifest["corpus_snapshot"] = {"sha256": sha_file(snap), "bytes": snap.stat().st_size}

git_sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True, check=False).stdout.strip()
doc = {"repo_git_sha": git_sha,
       "generated_by": "packages/topsem/scripts/make_input_manifest.py",
       "entries": manifest}
OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True), encoding="utf-8")
def _bytes(e: dict) -> int:
    return e["bytes"] if "bytes" in e else sum(x["bytes"] for x in e.values())


print(f"datasets={len(manifest)} total_bytes={sum(_bytes(e) for e in manifest.values()):,}")
print("->", OUT)
