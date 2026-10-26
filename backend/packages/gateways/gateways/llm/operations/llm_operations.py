"""LLM operations for text generation and completions."""

from __future__ import annotations

from openai.types.chat import ChatCompletionMessageParam

from gateways.llm.connections import get_openai_client


class LLMOperations:
    """High-level LLM operations using OpenAI client."""

    def __init__(self, api_key: str, model: str = "gpt-4", temperature: float = 0.7) -> None:
        """Initialize LLM operations with OpenAI client.

        Args:
            api_key: OpenAI API key
            model: Default model name
            temperature: Default temperature

        """
        self._client = get_openai_client(api_key=api_key)
        self._default_model = model
        self._default_temperature = temperature

    # VLM版とかも追加したい。
    def generate(
        self,
        messages: list[ChatCompletionMessageParam],
        model: str | None = None,
        temperature: float | None = None,
    ) -> str:
        """Generate completion using OpenAI's chat API.

        Args:
            messages: List of chat messages for the completion.
            model: Model name to use. If None, uses the default model.
            temperature: Temperature to use. If None, uses the default temperature.

        Returns:
            Generated completion text.

        """
        completion = self._client.chat.completions.create(
            model=model or self._default_model,
            messages=messages,
            temperature=temperature or self._default_temperature,
        )
        return completion.choices[0].message.content or ""
