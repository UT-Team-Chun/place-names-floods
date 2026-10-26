.PHONY: fmt
fmt:
	uv run ruff format .
	uv run ruff check .

.PHONY: fix
fix:
	uv run ruff format .
	uv run ruff check . --fix

.PHONY: lint
lint:
	uv run ruff check .
	uv run tach check
	uv run mypy .

.PHONY: test
test:
	uv run pytest . -vv

.PHONY: coverage
coverage:
	uv run coverage run -m pytest
	uv run coverage report
