# Harrison's Computational Research Repository & AI Harness

The following is a slightly cleaned-up and modernized version of the private
repository I use as the basis for the various research projects I've posted to
my personal blog.

It essentially provides me a sandbox for interactive investigation, validation,
authoring, etc.

One thing this repository cannot replicate on its own is an agent's interaction
with Logseq, which I use for notes.

## The Process

Most stages of the process are automated as Claude Code **skills** — checklisted
procedures in `.claude/skills/`, invoked as slash commands — which dispatch
specialized **subagents** defined in `.claude/agents/`. The subagents fall into
two roles:

* **Context isolators** (`source-reader`, `mathlib-scout`, `data-profiler`,
  `literature-checker`) read bulky material — papers, Mathlib source, raw
  result DataFrames — and return distilled, provenance-carrying summaries, so
  the raw material never enters the main conversation.
* **Fresh-eyes auditors** (`draft-auditor`, `test-auditor`,
  `statement-auditor`, `experiment-auditor`, `results-auditor`) independently
  re-check a deliverable against its specification and must return `CLEAN`
  before the skill may proceed. This guards against the failure mode where the
  work and its verification share a single misreading: an implementation and
  its tests written from the same misunderstanding of a page will pass a green
  gate together.

| Stage | Skill(s) | Isolators | Audit gate |
| --- | --- | --- | --- |
| Intake & Refinement | `/add-logseq-topic` | `source-reader` | `draft-auditor` |
| Development (Python) | `/add-python-topic` | — | `test-auditor` |
| Development (Lean) | `/add-lean-topic` | `mathlib-scout` | `statement-auditor` |
| Experimentation | `/design-experiment`, `/add-experiment` | — | `experiment-auditor` |
| Exploration | `/explore-results` | `data-profiler`, `literature-checker` | `results-auditor` |
| Feedback | `/critique-results` | `literature-checker` | — |

Each skill also maintains working files (fact maps, design specs, coverage
maps, drafts) in the session scratchpad as ground truth that survives context
compaction, and each ends by "closing the loop" — writing what was actually
built or found back to the relevant Logseq page.

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

Topic pages are written by `/add-logseq-topic <topic>`, which:

1. Checks the graph for existing coverage (extending rather than duplicating).
2. Researches from primary sources — the original paper or spec, standards,
   authoritative databases (OEIS/DLMF) — dispatching one `source-reader`
   subagent per source in parallel. Each returns distilled facts with
   provenance into a running fact → source map; raw source text stays out of
   the conversation.
3. Drafts the page in the house style: formal definition with all standard
   equivalent axiomatizations, attributed theorems written to double as
   property-test oracles, examples written as test fixtures, and
   implementation notes for the Python/Lean stages.
4. Gates publication on a `draft-auditor` pass — every claim in the draft must
   be supported by the fact map, with no strength drift or attribution
   mismatch — repeated until `CLEAN`.
5. Publishes to Logseq, verifies the published page against the draft, and
   reports the new red links it introduced (each a candidate for a future
   `/add-logseq-topic` run).

### Development

The concepts and relationships formed by the intake process are converted into
computational structures in the form of Python and Lean modules. Extensive unit
tests are written to ensure accuracy. In both tracks the Logseq page is the
specification — if it's missing or too thin, the skill stops and sends you back
to `/add-logseq-topic` rather than substituting its own research — and both
finish by rewriting the page's implementation notes from speculation into
verified fact.

`/add-python-topic <topic>` implements the object in `src/research`: the
page's definition fixes the representation, its derived vocabulary becomes
computed properties, its operations become transformation methods to/from
related objects, and every object is visualizable and serializable to a pandas
`DataFrame`. Tests are transcribed from the page — canonical examples as
fixtures, structural theorems as Hypothesis property tests — then audited by a
`test-auditor` subagent against a verbatim export of the page before
`just check` runs. The audit is the only independent read: the implementation
and suite were written from the same reading of the page, so a shared
misreading would otherwise pass green. A failing property test is treated as a
three-way question — implementation bug, mistranscribed test, or a wrong claim
on the page — and page errors are research findings, reported rather than
silently patched.

`/add-lean-topic <topic>` formalizes the object in `src/theorems`, with
Mathlib as the library of record. Parallel `mathlib-scout` subagents survey
the pinned Mathlib checkout per claim cluster, producing a coverage map that
assigns every page claim a decision: **reuse** (cite Mathlib), **extend** (a
repo declaration about Mathlib's object), **define** (Mathlib lacks it), or
**backlog**. Statements are written first with `sorry` bodies and gated on a
`statement-auditor` pass — a proved wrong statement is worth less than no
statement — before proof effort is spent. Proofs are checked for
counterexamples with `plausible` first, `just lean-check` must pass with a
clean axiom audit, and a `sorry` is never committed.

### Experimentation

Parameterized experiments are designed by composing the computational models
into independent scripts with standardized CLI arguments. This stage is split
across two skills so that design decisions are frozen before any code exists:

* `/design-experiment <hypothesis>` is interactive — every design element
  lands at an explicit checkpoint: falsifiable hypothesis, falsification
  criterion, parameter space, controls and baselines, replication and seed
  policy, and output schema (walked against the intended analysis, since an
  analysis the schema can't support is a design bug caught cheapest here).
  Predictions are **pre-registered** per region of the parameter space, citing
  the Logseq page or paper that motivates each, with no-prediction regions
  explicitly marked — that marking is what makes "unexpected" well-defined at
  exploration time. The result is a design doc in `docs/experiments/`, frozen
  once the first run exists.
* `/add-experiment <name>` implements the design doc as
  `src/experiments/<name>.py` using the `@experiment` harness (which supplies
  the standard CLI options and result/metadata plumbing), pilots it on a tiny
  grid into the session scratchpad — never `data/results/`, where a later
  exploration would resolve the pilot as a real run — and gates handoff on an
  `experiment-auditor` pass against the design doc: sweep ranges, seed policy,
  schema, and reproducibility metadata. A buggy experiment is exactly how
  "data that contradicts the literature" happens, so the audit sits between
  the pilot and the real sweep. The full sweep itself stays human-invoked:

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

This is the most autonomous stage, run as `/explore-results <name>[.<ts>]`
(given a second timestamp it also compares two runs of the same experiment),
under two hard rules: raw DataFrames never enter the conversation, and
contradiction with recorded literature is a bug until proven otherwise. The
workflow:

1. Resolve the run and extract the design doc's pre-registered predictions.
2. Profile the data through a `data-profiler` subagent — row count vs. the
   expected grid, NaN/inf, dtype surprises, coverage holes. An unexplained
   integrity finding stops the skill: a broken run gets reported, not
   interpreted.
3. Analyze via pandas scripts in the scratchpad that print aggregates, never
   raw frames — covering the design's analysis plan, then classifying every
   pre-registered prediction as **confirmed / contradicted / no prediction
   existed**.
4. Batch the contradictions and surprises into a `literature-checker`
   dispatch against the Logseq graph and `docs/ref/` papers, then triage
   bug-first: a `[contradicts]` disposition sends the agent hunting for the
   specific experiment bug that would produce it before any write-up.
5. Draft the results doc and gate publication on a `results-auditor` pass —
   every number recomputable from the data, no claim-strength drift, no
   surprise without its triage trail, no prediction left unaddressed.

### Feedback

Once a results document has been generated, I'll walk through it with an agent
tasked with critically analyzing it and ask questions, push back when necessary,
and refine interpretation. One critical task here is checking to see if any of
the data directly contradicts the literature (usually a sign of a bug in our
experiment).

Naturally, an "informal conjecture" will form from hypotheses refined by
successive testing via experiment. The key here is to conjure up experiments
that could disprove implications from previous results.

This walkthrough is `/critique-results <name>[.<ts>]` — a companion rather
than a pipeline, with the agent held in an explicitly skeptical stance: any
number I question is recomputed from the raw data with a fresh script (never
defended from the doc); alternative explanations I propose are steelmanned and
paired with the data that would discriminate between them; claims lacking a
literature disposition get a fresh `literature-checker` dispatch rather than
an assertion from memory; and each suggested follow-up is pressure-tested on
whether it could actually *falsify* the forming hypothesis, with
confirmation-only follow-ups replaced. Interpretation drifting into a region
the design marked "no prediction" is called out as exploration, not
confirmation. The session ends by editing the refined interpretation back into
the results doc (contested points left visible as open questions) and
recording the informal conjecture in Logseq — the seed the next phase
formalizes.

### Conjecture & Theorem Proving

Once a clear picture has formed, I'll work with an agent to develop a formal
conjecture explaining the results produced during the experimental phase and
save it in `docs/conj/` as well as Logseq.

An "attack plan" is produced for proving the conjecture mathematically from two
angles: _via chalkboard_ (i.e _traditionally_) and _computationally_ (via Lean).
Lean proofs are often _ugly for humans_ but are _safer_ and can provide the
scaffolding on which a more elegant, traditional proof may be built.

This phase has no dedicated skill yet — it starts from the informal-conjecture
notes that `/critique-results` accumulates in Logseq and is worked
interactively. The Lean side reuses the `/add-lean-topic` machinery where it
applies: `mathlib-scout` surveys for what the proof can lean on, and the
state-first / `statement-auditor` / prove discipline carries over.

### Figure Generation & Visualization

Figure generation for papers and blog posts is parameterized similarly to
experiments:

```
just figure <name> [args ...]
```

The above will produce one or more figures (or other visualization artifacts)
within `docs/fig/<name>/`. The format of files within each figure's directory
varies depending on the figure being generated.

Figure modules live in `src/figures` and are also written interactively rather
than through a dedicated skill. Scratch plots produced during exploration
never graduate directly into `docs/` — anything worth keeping is rebuilt as a
parameterized figure module here.

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

### Documentation

The site is built from `docs/` with [ProperDocs](https://properdocs.org/) — the
maintained continuation of MkDocs 1.x, by its original author — and configured
in `properdocs.yml`. It reads the same entry points as MkDocs, so the existing
theme and plugin ecosystem works unchanged.

Both API references are generated during the build rather than checked in — the
Python one from the sources under `src/` via
[mkdocstrings](https://mkdocstrings.github.io/), and the Lean one by reading
module and declaration docstrings out of `src/theorems/`. Neither needs
regenerating by hand, and the Lean reference needs no `lake` build. Cross-
references into third-party documentation are resolved from published
inventories, so a build needs network access.

| Command | Description |
| --- | --- |
| `just docs` | Build the documentation site into `site/`. |
| `just docs-serve [args ...]` | Serve the documentation with live reload. |
| `just docs-check` | Build the documentation, failing on any warning. |
