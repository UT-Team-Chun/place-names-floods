import os
from pathlib import Path

import pytest

from topsem.config import app_repo_root
from topsem.corpus.snapshot import CorpusSnapshot, load_snapshot
from topsem.runs import RunContext

needs_app_repo = pytest.mark.skipif(
    not (app_repo_root() / "data" / "curated").is_dir(),
    reason="disaster-toponyms が見つからない（TOPSEM_APP_REPO を設定）",
)


@needs_app_repo
def test_load_snapshot_reads_curated_yaml() -> None:
    snap = load_snapshot()
    assert isinstance(snap, CorpusSnapshot)
    assert len(snap.toponyms) > 0
    m = snap.manifest()
    assert len(m["app_git_sha"]) == 40
    assert Path(snap.app_root, "data/curated/toponyms/nationwide.yaml").is_file()


@needs_app_repo
def test_snapshot_is_deterministic() -> None:
    a = load_snapshot().manifest()
    b = load_snapshot().manifest()
    assert a == b


def test_missing_app_repo_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_snapshot(tmp_path)


def test_provenance_prefers_env_code_sha(monkeypatch: pytest.MonkeyPatch) -> None:
    """クラスタ（image に .git 無し）でも解析 repo SHA が来歴に残る回帰テスト。

    r3-M3 の対称形: app 来歴に TOPSEM_GIT_SHA を使うのは誤りだが、
    repo_git_sha（解析コードの来歴）では env が第一出所であるべき。
    """
    monkeypatch.setenv("TOPSEM_GIT_SHA", "a" * 40)
    ctx = RunContext(tag="t", run_dir=Path(), created_utc="now")
    assert ctx.provenance()["repo_git_sha"] == "a" * 40
    monkeypatch.delenv("TOPSEM_GIT_SHA", raising=False)
    sha = ctx.provenance()["repo_git_sha"]
    assert sha == "unknown" or (len(sha) == 40 and all(c in "0123456789abcdef" for c in sha))
    assert os.environ.get("TOPSEM_APP_SHA", "unset") == "unset"  # 混同防止の現認
