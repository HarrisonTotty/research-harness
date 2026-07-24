# Task runner for the Python stack. Run `just` (or `just --list`) to see recipes.
# Every recipe drives the project's tooling through `uv` so the pinned
# environment is used without a manual activation step.

set shell := ["bash", "-euo", "pipefail", "-c"]

# List available recipes.
default:
    @just --list

# Create/refresh the virtual environment and install all dependency groups.
install:
    uv sync

# Install the git pre-commit hooks into this clone.
hooks:
    uv run pre-commit install --install-hooks

# Lint with ruff (no changes written).
lint:
    uv run ruff check .

# Format sources and apply lint autofixes.
format:
    uv run ruff check --fix .
    uv run ruff format .

# Verify formatting without writing changes (for CI).
format-check:
    uv run ruff format --check .

# Type-check with mypy.
typecheck:
    uv run mypy

# Run the test suite.
test *args:
    uv run pytest {{ args }}

# Run the test suite with coverage reporting.
coverage *args:
    uv run pytest --cov --cov-report=term-missing {{ args }}

# Run every pre-commit hook against all files.
pre-commit:
    uv run pre-commit run --all-files

# Full local gate: lint, format check, type check, and tests.
check: lint format-check typecheck test

# Run a parameterized experiment, e.g. `just experiment my_experiment --n 100`.
experiment name *args:
    uv run python -m experiments {{ name }} {{ args }}

# Remove caches and build artifacts.
clean:
    rm -rf .ruff_cache .mypy_cache .pytest_cache .coverage htmlcov dist build
    find . -type d -name __pycache__ -prune -exec rm -rf {} +

# --- Lean (theorem proving) ------------------------------------------------
# Lean recipes drive `lake`, which is provided by the toolchain pinned in
# `lean-toolchain` (install via `elan`). They are kept separate from the Python
# `check` gate because building Mathlib is slow.

# Resolve Lean dependencies and download the Mathlib build cache.
lean-update:
    lake update
    lake exe cache get

# Build the Lean library.
lean-build:
    lake build

# Run Mathlib's environment linters over the library.
lean-lint:
    lake exe runLinter Theorems

# Full Lean gate: build, then lint.
lean-check: lean-build lean-lint

# Remove Lean build outputs.
lean-clean:
    lake clean
