"""RDB configuration."""

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

# Define project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent.parent.parent
ENV_FILE = PROJECT_ROOT / ".envrc"


class DatabaseConfig(BaseSettings):
    """Integrated database configuration for PostgreSQL."""

    model_config = SettingsConfigDict(
        env_file=str(ENV_FILE) if ENV_FILE.exists() else None,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    db_user: str
    db_pass: str
    db_host: str
    db_port: int
    db_name: str

    # Connection pool settings
    pool_size: int = 10
    max_overflow: int = 2
    pool_recycle: int = 300

    # Operational settings
    database_echo: bool = False  # Set to True to see SQL queries in logs

    @property
    def database_url(self) -> str:
        return f"postgresql://{self.db_user}:{self.db_pass}@{self.db_host}:{self.db_port}/{self.db_name}"


# Backward compatibility aliases
RDBConfig = DatabaseConfig
PostgresConfig = DatabaseConfig
