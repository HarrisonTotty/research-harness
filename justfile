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

# --- Documentation ---------------------------------------------------------
# Driven by `properdocs`, the maintained continuation of MkDocs 1.x, configured
# in `properdocs.yml`. The Python and Lean API references are generated during
# the build, so nothing under `docs/` needs regenerating by hand. Builds resolve
# cross-references into third-party inventories and so need network access.

# Build the documentation site into `site/`.
docs:
    uv run properdocs build

# Serve the documentation with live reload, e.g. `just docs-serve -a :9000`.
# `DOCS_SITE_URL` mounts the preview at `/` instead of under the deployed
# project-page path; `serve` replaces the host and port with the ones in use.
docs-serve *args:
    DOCS_SITE_URL=http://localhost/ uv run properdocs serve {{ args }}

# Build the documentation, failing on any warning (broken link, orphan page).
docs-check:
    uv run properdocs build --strict

# Remove caches and build artifacts.
clean:
    rm -rf .ruff_cache .mypy_cache .pytest_cache .coverage htmlcov dist build site
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
