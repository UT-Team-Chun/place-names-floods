"""Tests for LLM operations module."""

from unittest.mock import MagicMock, patch

from gateways.llm.operations.llm_operations import LLMOperations


def test_llm_operations_initialization():
    """Test LLMOperations initialization."""
    with patch("gateways.llm.operations.llm_operations.get_openai_client") as mock_client:
        mock_client.return_value = MagicMock()

        operations = LLMOperations(
            api_key="test-key",
            model="gpt-4",
            temperature=0.7,
        )

        assert operations._default_model == "gpt-4"
        assert operations._default_temperature == 0.7
        mock_client.assert_called_once()


def test_llm_operations_generate():
    """Test LLMOperations.generate method."""
    with patch("gateways.llm.operations.llm_operations.get_openai_client") as mock_get_client:
        # Mock OpenAI client
        mock_client = MagicMock()
        mock_completion = MagicMock()
        mock_completion.choices = [MagicMock()]
        mock_completion.choices[0].message.content = "Test response"
        mock_client.chat.completions.create.return_value = mock_completion
        mock_get_client.return_value = mock_client

        operations = LLMOperations(api_key="test-key")

        messages = [{"role": "user", "content": "Hello"}]
        result = operations.generate(messages)  # type: ignore[arg-type]

        assert result == "Test response"
        mock_client.chat.completions.create.assert_called_once()


def test_llm_operations_generate_with_custom_params():
    """Test LLMOperations.generate with custom parameters."""
    with patch("gateways.llm.operations.llm_operations.get_openai_client") as mock_get_client:
        # Mock OpenAI client
        mock_client = MagicMock()
        mock_completion = MagicMock()
        mock_completion.choices = [MagicMock()]
        mock_completion.choices[0].message.content = "Custom response"
        mock_client.chat.completions.create.return_value = mock_completion
        mock_get_client.return_value = mock_client

        operations = LLMOperations(
            api_key="test-key",
            model="gpt-4",
            temperature=0.7,
        )

        messages = [{"role": "user", "content": "Hello"}]
        result = operations.generate(
            messages,  # type: ignore[arg-type]
            model="gpt-4-turbo",
            temperature=0.5,
        )

        assert result == "Custom response"

        # Verify custom parameters were used
        call_args = mock_client.chat.completions.create.call_args
        assert call_args.kwargs["model"] == "gpt-4-turbo"
        assert call_args.kwargs["temperature"] == 0.5
