# Harrison's Computational Research Repository & AI Harness

The following is a slightly cleaned-up and modernized version of the private
repository I use as the basis for the various research projects I've posted to
my personal blog.

It essentially provides me a sandbox for interactive investigation, validation,
authoring, etc.

One thing this repository cannot replicate on its own is an agent's interaction
with Logseq, which I use for notes.

## The Process

### Intake & Refinement

Relevant (or interesting) research papers are ingested by AI agents which record
concepts and their relationships to Logseq. These relationships are checked and
refined such that:

1. Existing information is not duplicated.
2. Non-standard author syntax is clarified (and sometimes converted).
3. Keywords and terms are sufficiently linked.

In addition, it is good to go back through and identify referenced terms,
keywords, and topics that don't yet have a dedicated page and write one. I also
save a copy of the paper PDF to `docs/ref/` incase I need to re-ingest.

### Development

The concepts and relationships formed by the intake process are converted into
computational structures in the form of Python and Lean modules. Extensive unit
tests are written to ensure accuracy.

### Experimentation

Parameterized experiments are designed by composing the computational models
into independent scripts with standardized CLI arguments and invoked via:

```
just experiment <name> [args ...]
```

The above will produce timestamped experimental results in the form of a
serialized Pandas `DataFrame` at `data/results/{experiment}.{timestamp}.json`,
along with a accompanying `{experiment}.{timestamp}.meta.json` recording the
input parameters and other relevant metadata.

### Exploration

Agents are instructed to explore the output of experiments and produce a
corresponding `docs/results/{experiment}.{timestamp}.md` providing an overview
of the experiment's findings - including anything unexpected in the data - and
suggestions for potential follow-up experiments.

### Feedback

Once a results document has been generated, I'll walk through it with an agent
tasked with critically analyzing it and ask questions, push back when necessary,
and refine interpretation. One critical task here is checking to see if any of
the data directly contradicts the literature (usually a sign of a bug in our
experiment).

Naturally, an "informal conjecture" will form from hypotheses refined by
successive testing via experiment. The key here is to conjure up experiments
that could disprove implications from previous results.

### Conjecture & Theorem Proving

Once a clear picture has formed, I'll work with an agent to develop a formal
conjecture explaining the results produced during the experimental phase and
save it in `docs/conj/` as well as Logseq.

An "attack plan" is produced for proving the conjecture mathematically from two
angles: _via chalkboard_ (i.e _traditionally_) and _computationally_ (via Lean).
Lean proofs are often _ugly for humans_ but are _safer_ and can provide the
scaffolding on which a more elegant, traditional proof may be built.

## Development

The Python stack uses [`uv`](https://docs.astral.sh/uv/) for environment and
dependency management, [`ruff`](https://docs.astral.sh/ruff/) for linting and
formatting, [`mypy`](https://mypy-lang.org/) for type checking, and
[`pytest`](https://docs.pytest.org/) for testing. The Lean stack uses
[`lake`](https://github.com/leanprover/lean4) with
[Mathlib](https://github.com/leanprover-community/mathlib4); the Lean version is
pinned in `lean-toolchain` and managed with
[`elan`](https://github.com/leanprover/elan). Everything is orchestrated through
[`just`](https://github.com/casey/just), so no manual environment activation is
needed. Run `just` (or `just --list`) to see all recipes.

### Python

| Command | Description |
| --- | --- |
| `just check` | Full local gate: lint, format check, type check, and tests. |
| `just clean` | Remove caches and build artifacts. |
| `just coverage [args ...]` | Run the suite with coverage reporting. |
| `just experiment <name> [args ...]` | Run a parameterized experiment. |
| `just format-check` | Verify formatting without writing changes. |
| `just format` | Apply `ruff` lint autofixes, then format sources. |
| `just hooks` | Install the git pre-commit hooks into this clone. |
| `just install` | Create/refresh the virtual environment and install all dependency groups. |
| `just lint` | Lint with `ruff` (reports only, writes nothing). |
| `just pre-commit` | Run every pre-commit hook against all files. |
| `just test [args ...]` | Run the `pytest` suite (extra args pass through). |
| `just typecheck` | Type-check with `mypy`. |

### Lean

Building [Mathlib](https://github.com/leanprover-community/mathlib4) is slow, so
these recipes are kept separate from `just check`. On a fresh clone, run
`just lean-update` once to resolve dependencies and download the prebuilt
Mathlib cache before the first `just lean-build`.

| Command | Description |
| --- | --- |
| `just lean-update` | Resolve dependencies and fetch the Mathlib build cache. |
| `just lean-build` | Build the Lean library. |
| `just lean-lint` | Run Mathlib's environment linters over the library. |
| `just lean-check` | Full Lean gate: build, then lint. |
| `just lean-clean` | Remove Lean build outputs. |
