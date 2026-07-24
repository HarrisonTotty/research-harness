# Computational Research Sandbox

This repository provides a suite of tools, workflows, and libaries to assist in
the exploration, validation, and publication of various research topics.

## Important Directories

* `data/results`: Raw results from experiments.
* `docs/agent/brainstorming`: General code feature brainstorming.
* `docs/agent/plans`: Phased action plans for implementing larger features.
* `docs/conj`: Proposed mathematical conjectures.
* `docs/ref`: Reference literature (papers).
* `docs/results`: Interpreted results from experiments.
* `docs/theorems`: Proven mathematical theorems (with proofs).
* `src/research`: Primary research library (Python) source code.
* `src/experiments`: Research experiments (Python) source code.
* `src/theorems`: Theorem proving (Lean) source code.
* `tests`: Unit tests.

## Common Commands

Both stacks are driven through `just` recipes. Python recipes invoke their tool
via `uv`; Lean recipes invoke `lake` (from the toolchain pinned in
`lean-toolchain`). Run `just` (or `just --list`) to see everything.

Python:

* `just check`: Full local gate — lint, format check, type check, and tests.
* `just coverage [args ...]`: Run the suite with coverage reporting.
* `just experiment <name> [args ...]`: Run a parameterized experiment.
* `just format-check`: Verify formatting without writing changes.
* `just format`: Apply `ruff` lint autofixes, then format sources.
* `just lint`: Lint with `ruff` (reports only, writes nothing).
* `just test [args ...]`: Run the `pytest` suite (extra args pass through).
* `just typecheck`: Type-check with `mypy`.

Lean (kept separate from `just check` because building Mathlib is slow):

* `just lean-update`: Resolve dependencies and fetch the Mathlib build cache.
* `just lean-build`: Build the Lean library.
* `just lean-lint`: Run Mathlib's environment linters over the library.
* `just lean-check`: Full Lean gate — build, then lint.
* `just lean-clean`: Remove Lean build outputs.
