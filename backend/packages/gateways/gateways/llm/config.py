"""LLM configuration from environment variables."""

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

# Define backend root directory (backend/)
BACKEND_ROOT = Path(__file__).parent.parent.parent.parent.parent
ENV_FILE = BACKEND_ROOT.parent / ".envrc"


class LLMConfig(BaseSettings):
    """LLM-related configuration."""

    model_config = SettingsConfigDict(
        env_file=str(ENV_FILE) if ENV_FILE.exists() else None,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    openai_api_key: str | None = None
    llm_provider: str = "openai"
    llm_model: str = "gpt-4"
    temperature: float = 0.7
    max_tokens: int = 1000
    embedding_model: str = "text-embedding-3-small"
