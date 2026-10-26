"""Tests for Alembic migration generation and execution."""

import os
import subprocess
from pathlib import Path

import pytest
from sqlalchemy import inspect
from sqlmodel import create_engine

from gateways.rdb.config import DatabaseConfig


@pytest.fixture
def test_db_url() -> str:
    """Test用のデータベースURLを返す."""
    return "sqlite:///./test_migration.db"


@pytest.fixture
def test_engine(test_db_url: str):
    """Test用のエンジンを作成."""
    engine = create_engine(test_db_url, echo=False)
    yield engine
    engine.dispose()
    # テスト後にDBファイルを削除
    db_path = Path("./test_migration.db")
    if db_path.exists():
        db_path.unlink()


@pytest.fixture
def migration_dir() -> Path:
    """マイグレーションディレクトリのパスを返す."""
    return Path(__file__).parent.parent.parent / "gateways" / "rdb" / "migration"


@pytest.fixture
def alembic_ini_path() -> Path:
    """alembic.iniのパスを返す."""
    return Path(__file__).parent.parent.parent / "gateways" / "rdb" / "alembic.ini"


class TestMigration:
    """マイグレーション関連のテスト."""

    def test_migration_directory_exists(self, migration_dir: Path) -> None:
        """マイグレーションディレクトリが存在するか確認."""
        assert migration_dir.exists(), f"Migration directory does not exist: {migration_dir}"
        assert migration_dir.is_dir(), f"Migration path is not a directory: {migration_dir}"

    def test_alembic_ini_exists(self, alembic_ini_path: Path) -> None:
        """alembic.iniが存在するか確認."""
        assert alembic_ini_path.exists(), f"alembic.ini does not exist: {alembic_ini_path}"

    @pytest.mark.skip(
        reason="Versions directory may not exist until migration files are generated"
    )
    def test_versions_directory_exists(self, migration_dir: Path) -> None:
        """versions/ディレクトリが存在するか確認."""
        versions_dir = migration_dir / "versions"
        assert versions_dir.exists(), f"Versions directory does not exist: {versions_dir}"
        assert versions_dir.is_dir(), f"Versions path is not a directory: {versions_dir}"

    def test_env_py_exists(self, migration_dir: Path) -> None:
        """env.pyが存在するか確認."""
        env_py = migration_dir / "env.py"
        assert env_py.exists(), f"env.py does not exist: {env_py}"

    @pytest.mark.skip(reason="Migration files not yet generated")
    def test_migration_files_exist(self, migration_dir: Path) -> None:
        """マイグレーションファイルが少なくとも1つ存在するか確認."""
        versions_dir = migration_dir / "versions"
        migration_files = list(versions_dir.glob("*.py"))
        # __init__.pyを除外
        migration_files = [f for f in migration_files if f.name != "__init__.py"]
        assert len(migration_files) > 0, "No migration files found"

    @pytest.mark.skip(
        reason="Requires modification to env.py to support DATABASE_URL environment variable"
    )
    def test_alembic_current_command(self, alembic_ini_path: Path) -> None:
        """alembic currentコマンドが正常に実行できるか確認.

        Note: このテストを有効にするには、
        env.pyでDATABASE_URL環境変数を読み込むように修正が必要です。
        """
        # Test用のデータベースURLを環境変数に設定
        test_env = os.environ.copy()
        test_env["DATABASE_URL"] = "sqlite:///./test_migration_current.db"

        result = subprocess.run(  # noqa: S603
            ["alembic", "-c", str(alembic_ini_path), "current"],  # noqa: S607
            capture_output=True,
            text=True,
            cwd=alembic_ini_path.parent,
            env=test_env,
            check=True,
        )
        # エラーが発生していないか確認（存在しない場合は空の出力）
        assert result.returncode == 0, f"alembic current failed: {result.stderr}"

        # テスト後にDBファイルを削除
        db_path = Path(alembic_ini_path.parent) / "test_migration_current.db"
        if db_path.exists():
            db_path.unlink()

    @pytest.mark.skip(reason="Migration files not yet generated")
    def test_alembic_heads_command(self, alembic_ini_path: Path) -> None:
        """alembic headsコマンドが正常に実行できるか確認."""
        result = subprocess.run(  # noqa: S603
            ["alembic", "-c", str(alembic_ini_path), "heads"],  # noqa: S607
            capture_output=True,
            text=True,
            cwd=alembic_ini_path.parent,
            check=True,
        )
        assert result.returncode == 0, f"alembic heads failed: {result.stderr}"
        # 少なくとも1つのheadが存在することを確認
        assert len(result.stdout.strip()) > 0, "No migration heads found"

    @pytest.mark.skip(reason="Migration files not yet generated")
    def test_alembic_history_command(self, alembic_ini_path: Path) -> None:
        """alembic historyコマンドが正常に実行できるか確認."""
        result = subprocess.run(  # noqa: S603
            ["alembic", "-c", str(alembic_ini_path), "history"],  # noqa: S607
            capture_output=True,
            text=True,
            cwd=alembic_ini_path.parent,
            check=True,
        )
        assert result.returncode == 0, f"alembic history failed: {result.stderr}"
        # 少なくとも1つの履歴が存在することを確認
        assert len(result.stdout.strip()) > 0, "No migration history found"

    @pytest.mark.skip(
        reason="Requires modification to env.py to support DATABASE_URL environment variable"
    )
    def test_migration_upgrade_downgrade(self, alembic_ini_path: Path, tmp_path: Path) -> None:
        """マイグレーションのupgrade/downgradeが正常に動作するか確認.

        Note: このテストを有効にするには、
        env.pyでDATABASE_URL環境変数を読み込むように修正が必要です。
        """
        # Test用のデータベースファイルパス
        test_db_file = tmp_path / "test_migration_upgrade.db"
        test_db_url = f"sqlite:///{test_db_file}"

        # Test用のデータベースURLを環境変数に設定
        test_env = os.environ.copy()
        test_env["DATABASE_URL"] = test_db_url

        # データベースが空であることを確認
        test_engine = create_engine(test_db_url, echo=False)
        inspector = inspect(test_engine)
        initial_tables = inspector.get_table_names()
        assert len(initial_tables) == 0, "Database should be empty initially"

        # alembic upgrade head を実行
        upgrade_result = subprocess.run(  # noqa: S603
            ["alembic", "-c", str(alembic_ini_path), "upgrade", "head"],  # noqa: S607
            capture_output=True,
            text=True,
            cwd=alembic_ini_path.parent,
            env=test_env,
            check=True,
        )
        assert upgrade_result.returncode == 0, f"alembic upgrade failed: {upgrade_result.stderr}"

        # テーブルが作成されたことを確認
        inspector = inspect(test_engine)
        tables_after_upgrade = inspector.get_table_names()
        assert len(tables_after_upgrade) > 0, "No tables created after upgrade"
        # alembic_versionテーブルが存在することを確認
        assert "alembic_version" in tables_after_upgrade, "alembic_version table not found"

        # alembic downgrade base を実行
        downgrade_result = subprocess.run(  # noqa: S603
            ["alembic", "-c", str(alembic_ini_path), "downgrade", "base"],  # noqa: S607
            capture_output=True,
            text=True,
            cwd=alembic_ini_path.parent,
            env=test_env,
            check=True,
        )
        assert downgrade_result.returncode == 0, (
            f"alembic downgrade failed: {downgrade_result.stderr}"
        )

        # alembic_versionテーブル以外のテーブルが削除されたことを確認
        inspector = inspect(test_engine)
        tables_after_downgrade = inspector.get_table_names()
        # alembic_versionテーブルのみが残る
        assert set(tables_after_downgrade) <= {"alembic_version"}, (
            f"Tables not properly cleaned after downgrade: {tables_after_downgrade}"
        )

        # エンジンをクローズ
        test_engine.dispose()

    def test_database_config_get_url(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """DatabaseConfig.get_url()が正しいURL形式を返すか確認."""
        # テスト用の環境変数を設定
        monkeypatch.setenv("DB_USER", "test_user")
        monkeypatch.setenv("DB_PASS", "test_pass")
        monkeypatch.setenv("DB_HOST", "localhost")
        monkeypatch.setenv("DB_PORT", "5432")
        monkeypatch.setenv("DB_NAME", "test_db")

        config = DatabaseConfig()
        url = config.database_url
        assert isinstance(url, str), "URL should be a string"
        assert url.startswith(("postgresql://", "sqlite://")), (
            "URL should start with postgresql:// or sqlite://"
        )

    @pytest.mark.skip(reason="Migration files not yet generated")
    def test_migration_file_naming_convention(self, migration_dir: Path) -> None:
        """マイグレーションファイルの命名規則が正しいか確認."""
        versions_dir = migration_dir / "versions"
        migration_files = list(versions_dir.glob("*.py"))
        # __init__.pyを除外
        migration_files = [f for f in migration_files if f.name != "__init__.py"]

        for migration_file in migration_files:
            # ファイル名が {revision}_{slug}.py の形式であることを確認
            assert "_" in migration_file.stem, (
                f"Invalid migration file name: {migration_file.name}"
            )
            parts = migration_file.stem.split("_", 1)
            assert len(parts) == 2, f"Invalid migration file name format: {migration_file.name}"
            revision, slug = parts
            # revisionは英数字のみ
            assert revision.isalnum(), f"Invalid revision format: {revision}"
            # slugは空でない
            assert len(slug) > 0, f"Empty slug in migration file: {migration_file.name}"

    @pytest.mark.skip(reason="Migration files not yet generated")
    def test_migration_file_content_structure(self, migration_dir: Path) -> None:
        """マイグレーションファイルの内容構造が正しいか確認."""
        versions_dir = migration_dir / "versions"
        migration_files = list(versions_dir.glob("*.py"))
        # __init__.pyを除外
        migration_files = [f for f in migration_files if f.name != "__init__.py"]

        for migration_file in migration_files:
            content = migration_file.read_text()
            # 必須要素が含まれているか確認
            # revision = または revision: str = の形式を許容
            assert "revision = " in content or "revision: str =" in content, (
                f"Missing revision in {migration_file.name}"
            )
            # down_revision = または down_revision: Union[str, None] = の形式を許容
            assert (
                "down_revision = " in content or "down_revision: Union[str, None] =" in content
            ), f"Missing down_revision in {migration_file.name}"
            assert "def upgrade()" in content, (
                f"Missing upgrade() function in {migration_file.name}"
            )
            assert "def downgrade()" in content, (
                f"Missing downgrade() function in {migration_file.name}"
            )
