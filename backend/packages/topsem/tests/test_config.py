"""TOPSEM_DATA_DIR / TOPSEM_RUNS_DIR のオーバーライド挙動（クラスタ投入契約）。"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest


def test_env_overrides(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import topsem.config as cfg  # noqa: PLC0415

    monkeypatch.setenv("TOPSEM_DATA_DIR", str(tmp_path / "d"))
    monkeypatch.setenv("TOPSEM_RUNS_DIR", str(tmp_path / "r"))
    try:
        reloaded = importlib.reload(cfg)
        assert tmp_path / "d" == reloaded.DATA_DIR
        assert tmp_path / "r" == reloaded.RUNS_DIR
        assert tmp_path / "d" / "labels" == reloaded.LABELS_DIR
    finally:
        monkeypatch.delenv("TOPSEM_DATA_DIR")
        monkeypatch.delenv("TOPSEM_RUNS_DIR")
        importlib.reload(cfg)


def test_defaults_without_env() -> None:
    import topsem.config as cfg  # noqa: PLC0415

    assert cfg.DATA_DIR.name == "data"
    assert cfg.RUNS_DIR == cfg.DATA_DIR / "runs"
