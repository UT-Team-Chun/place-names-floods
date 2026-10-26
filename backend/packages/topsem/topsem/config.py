"""リポジトリ内の共通パス設定。

アプリ側リポジトリ（disaster-toponyms）は `TOPSEM_APP_REPO` 環境変数で指定向こう
（デフォルトは本リポジトリの姉妹ディレクトリ）。研究リポジトリ側にはコーパスを
コピーしない（CLAUDE.md「データの来歴」）。

クラスタ実行用オーバーライド:
- `TOPSEM_DATA_DIR`: 入力データ（labels/population/covariates/splits）の root。
  NAS を指してよい（例 /mnt/nas/toponym-semantics/data）
- `TOPSEM_RUNS_DIR`: run 出力先。既定は DATA_DIR/runs（クラスタではローカルへ
  出して job が NAS runs/ へ cp する）
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = Path(os.environ.get("TOPSEM_DATA_DIR") or REPO_ROOT / "data")
RUNS_DIR = Path(os.environ.get("TOPSEM_RUNS_DIR") or DATA_DIR / "runs")
CORPUS_DIR = DATA_DIR / "corpus"
LABELS_DIR = DATA_DIR / "labels"
SPLITS_DIR = DATA_DIR / "splits"
COVARIATES_DIR = DATA_DIR / "covariates"


def app_repo_root() -> Path:
    """アプリ側リポジトリのルート（存在チェックは呼び出し側の責務）。"""
    return Path(os.environ.get("TOPSEM_APP_REPO", str(REPO_ROOT.parent / "disaster-toponyms")))
