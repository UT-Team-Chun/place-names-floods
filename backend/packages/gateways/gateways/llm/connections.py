"""OpenAI connection client."""

from __future__ import annotations

from functools import lru_cache

from openai import OpenAI


@lru_cache(maxsize=1)
def get_openai_client(api_key: str | None = None) -> OpenAI:
    """Get singleton instance of OpenAI client.

    Args:
        api_key: OpenAI API key. If None, uses the environment variable.

    Returns:
        OpenAI: Cached OpenAI client instance.

    Raises:
        ValueError: If api_key is not provided and OPENAI_API_KEY env var is not set.

    """
    if api_key:
        return OpenAI(api_key=api_key)
    return OpenAI()  # Uses OPENAI_API_KEY from environment
