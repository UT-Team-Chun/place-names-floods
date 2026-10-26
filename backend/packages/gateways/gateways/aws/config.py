"""AWS configuration from environment variables."""

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

# Define backend root directory (backend/)
BACKEND_ROOT = Path(__file__).parent.parent.parent.parent.parent
ENV_FILE = BACKEND_ROOT.parent / ".envrc"


class AWSConfig(BaseSettings):
    """AWS-related configuration."""

    model_config = SettingsConfigDict(
        env_file=str(ENV_FILE) if ENV_FILE.exists() else None,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    access_key_id: str | None = None
    secret_access_key: str | None = None
    region: str | None = "ap-northeast-1"
