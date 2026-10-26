"""アプリ側リポジトリ（disaster-toponyms）のコーパスを git SHA で固定して読む。

研究リポジトリには地名データをコピーしない。代わりに app の HEAD SHA と主要
ファイルの sha256 を manifest に残し、どの状態のコーパスで回したかを run ごとに
特定できるようにする（CLAUDE.md「データの来歴」）。
"""

from __future__ import annotations

import hashlib
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from topsem.config import REPO_ROOT, app_repo_root

_PIN_PATH = REPO_ROOT / "corpus_pin.yml"

KEY_FILES: tuple[str, ...] = (
    "data/curated/toponyms/nationwide.yaml",
    "data/curated/elements.yaml",
    "data/curated/sources.yaml",
)
_TOPLYMNMS_REL = "data/curated/toponyms/nationwide.yaml"
_MSG_NO_APP_REPO = (
    "アプリ側リポジトリ（data/curated を含む）が見つからない。TOPSEM_APP_REPO を確認"
)


def git_sha(root: Path) -> str:
    """root が属する git リポジトリの HEAD SHA を返す。

    クラスタでは initContainer が /workspace/corpus をピン SHA に checkout 済み
    なので、まず root 自身の git を読む（r3-M3: TOPSEM_GIT_SHA は解析 repo の
    SHA であり app の来歴に使うのは誤りだった。局所実行で root が姉妹 repo の
    場合もここが正しい）。失敗時のみ TOPSEM_APP_SHA、それも無ければ unknown。
    """
    env = os.environ.get("TOPSEM_APP_SHA")
    if env and env != "unknown":
        return env
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"
    return proc.stdout.strip()


@dataclass(frozen=True)
class CorpusSnapshot:
    """スナップショット時点のコーパス（読み取り専用）。"""

    app_root: Path
    app_sha: str
    file_hashes: dict[str, str]
    toponyms: tuple[dict[str, Any], ...]
    elements: tuple[dict[str, Any], ...]
    pinned_sha: str = ""

    def manifest(self) -> dict[str, Any]:
        m: dict[str, Any] = {
            "app_repo": str(self.app_root),
            "app_git_sha": self.app_sha,
            "file_sha256": dict(self.file_hashes),
            "n_toponyms": len(self.toponyms),
            "n_elements": len(self.elements),
        }
        if self.pinned_sha:
            m["pinned_git_sha"] = self.pinned_sha
            m["pin_matches"] = self.pinned_sha == self.app_sha
        return m


def load_pin() -> dict[str, str]:
    """corpus_pin.yml（クラスタ実行が clone する版）を読む。無れば空。"""
    if _PIN_PATH.is_file():
        data = yaml.safe_load(_PIN_PATH.read_text(encoding="utf-8")) or {}
        return {str(k): str(v) for k, v in data.items()}
    return {}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_snapshot(app_root: Path | None = None) -> CorpusSnapshot:
    """キュレーション済み YAML を app リポジトリから読む（コピーは作らない）。"""
    root = app_root if app_root is not None else app_repo_root()
    if not (root / "data" / "curated").is_dir():
        msg = f"{root}: {_MSG_NO_APP_REPO}"
        raise FileNotFoundError(msg)
    hashes = {rel: _sha256(root / rel) for rel in KEY_FILES if (root / rel).is_file()}
    topos = yaml.safe_load((root / _TOPLYMNMS_REL).read_text(encoding="utf-8"))
    elems_path = root / "data" / "curated" / "elements.yaml"
    elems = yaml.safe_load(elems_path.read_text(encoding="utf-8")) if elems_path.is_file() else {}
    return CorpusSnapshot(
        app_root=root,
        app_sha=git_sha(root),
        file_hashes=hashes,
        toponyms=tuple((topos or {}).get("toponyms", [])),
        elements=tuple((elems or {}).get("elements", [])),
        pinned_sha=load_pin().get("app_git_sha", ""),
    )
