"""実行ディレクトリ契約 `data/runs/<run>/`（CLAUDE.md）。

すべての解析はここに少なくとも以下を出す:

- `summary.json` — 設定・来歴・主要指標（下の provenance() を必ず含む）
- `metrics.json` — fold / event 別の詳細
- `predictions/` — セル単位の生予測（再解析可能にする）
- `diagnostics.png` — 校正・イベント別内訳・特徴重要度
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from topsem.config import REPO_ROOT, RUNS_DIR
from topsem.corpus.snapshot import git_sha


@dataclass(frozen=True)
class RunContext:
    """1 実行 = 1 ディレクトリ。書き込み先名前をここから取る。"""

    tag: str
    run_dir: Path
    created_utc: str

    def path(self, name: str) -> Path:
        return self.run_dir / name

    def provenance(self) -> dict[str, Any]:
        import platform  # noqa: PLC0415

        dev = os.environ.get("TOPSEM_DEVICE")
        if not dev:
            dev = "cluster" if Path("/mnt/nas/toponym-semantics").is_dir() else "local"
        return {
            "run_tag": self.tag,
            "created_utc": self.created_utc,
            # 解析 repo 自身の来歴: クラスタでは .git を image に載せないため
            # 投入側が流す TOPSEM_GIT_SHA（解析 repo SHA）を第一出所にする。
            # snapshot 側の app 来歴とは混同しないこと（r3-M3 の教訓の対）
            "repo_git_sha": os.environ.get("TOPSEM_GIT_SHA") or git_sha(REPO_ROOT),
            "device": dev,
            "platform": platform.platform(),
            "argv": list(sys.argv),
            # A-4 等の発効状態を来歴に刻む（env が設計を左右するので必須）
            "env_flags": {k: v for k, v in sorted(os.environ.items())
                          if k.startswith("TOPSEM_")},
        }


def start_run(tag: str) -> RunContext:
    """`data/runs/<tag>/` を作り実行コンテキストを返す（既存なら再利用）。"""
    run_dir = RUNS_DIR / tag
    run_dir.mkdir(parents=True, exist_ok=True)
    return RunContext(
        tag=tag,
        run_dir=run_dir,
        created_utc=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
