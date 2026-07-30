# Computational Research Sandbox

This repository provides a suite of tools, workflows, and libaries to assist in
the exploration, validation, and publication of various research topics.

## Important Directories

* `data/results`: Raw results from experiments.
* `docs/agent/brainstorming`: General code feature brainstorming.
* `docs/agent/plans`: Phased action plans for implementing larger features.
* `docs/conj`: Proposed mathematical conjectures.
* `docs/experiments`: Experiment design docs with pre-registered predictions.
* `docs/fig`: Generated figures for placement in blog posts & papers.
* `docs/ref`: Reference literature (papers).
* `docs/results`: Interpreted results from experiments.
* `docs/theorems`: Proven mathematical theorems (with proofs).
* `scripts/docs`: Generators and hooks for the `properdocs` documentation build.
* `src/research`: Primary research library (Python) source code.
* `src/experiments`: Research experiments (Python) source code.
* `src/figures`: Paper & blog post figures and visualizations (Python) source code.
* `src/theorems`: Theorem proving (Lean) source code.
* `tests`: Unit tests.

## Common Commands

Both stacks are driven through `just` recipes. Python recipes invoke their tool
via `uv`; Lean recipes invoke `lake` (from the toolchain pinned in
`lean-toolchain`). Run `just` (or `just --list`) to see everything.

Python:

* `just check`: Full local gate — lint, format check, type check, and tests.
* `just coverage [args ...]`: Run the suite with coverage reporting.
* `just docs-check`: Build the documentation, failing on any warning.
* `just docs-serve [args ...]`: Serve the documentation with live reload.
* `just docs`: Build the documentation site into `site/`.
* `just experiment <name> [args ...]`: Run a parameterized experiment.
* `just figure <name> [args ...]`: Generate figures/visualizations.
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

The site is built with `properdocs` (the maintained continuation of MkDocs 1.x)
and configured in `properdocs.yml`. Both API references are generated during the
build — never write them into `docs/` by hand. Python modules are documented by
`mkdocstrings` from their sources; Lean modules by reading docstrings out of
`src/theorems` (no `lake` build involved), so a docstring is the only way a
declaration gets documented. The docs recipes are kept out of `just check`
because a build resolves cross-references over the network.
