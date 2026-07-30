# Skills & agents for the Experimentation and Exploration phases

Brainstorming for tooling that covers the middle of the research process:
**Experimentation** (design and implement parameterized experiments),
**Exploration** (turn `data/results/*.json` into `docs/results/*.md`), and the
front half of **Feedback** (critical walkthrough with the user). Unlike intake
and development, experiment *design* is inherently interactive — the tooling
should structure the back-and-forth rather than try to eliminate it.
Exploration, by contrast, can be nearly autonomous.

## House patterns to reuse

The intake/development skills established four patterns worth carrying forward:

1. **Checklist skills with scratchpad working files** (`page.md`, `spec.md`)
   that survive context compaction and serve as ground truth for later steps.
2. **Fresh-eyes auditor agents** (`test-auditor`, `draft-auditor`,
   `statement-auditor`) that gate a step: the artifact and its source were
   produced by the same reading, so only an independent read breaks correlated
   errors. Machine-consumed findings format, `CLEAN` sentinel, loop until
   clean.
3. **Scout agents that keep bulk content out of the main context**
   (`source-reader`, `mathlib-scout`): return distilled, structured facts,
   never raw dumps.
4. **Logseq as specification**, with a close-the-loop step writing back what
   was actually built.

## Design principles specific to these phases

- **Pre-registration.** The design step records, *before the experiment runs*,
  what the literature/theory predicts for each region of the parameter space
  and what outcome would falsify the hypothesis. This makes the README's
  "anything unexpected in the data" well-defined at exploration time and
  guards against post-hoc rationalization.
- **Raw DataFrames never enter the main context.** All analysis happens
  through scripted pandas (via `uv run` in the scratchpad) or through a
  profiler subagent that returns summaries. Result files are large timestamped
  JSON blobs; reading them directly is both wasteful and error-prone.
- **Contradiction with literature is a bug until proven otherwise.** The
  README calls this out: data that contradicts known results usually means a
  broken experiment. Exploration must check this *before* the results doc
  presents a "surprising finding".
- **Falsification orientation.** Follow-up suggestions should preferentially
  be experiments that could *disprove* implications of the current results,
  not ones that would merely accumulate confirmation.

## Proposed artifact chain

```
hypothesis (Logseq / prior results doc)
  → docs/experiments/<name>.md          (design + pre-registration record)
  → src/experiments/<name>.py           (implementation, audited against design)
  → data/results/<name>.<ts>.json       (full sweeps; user-invoked via `just experiment`)
    + data/results/<name>.<ts>.meta.json
  → docs/results/<name>.<ts>.md         (exploration output, audited against data)
  → refined interpretation + follow-ups (feedback walkthrough)
```

New directory decision: `docs/experiments/` holds one checked-in design doc
per experiment, parallel to `docs/results/` and `docs/conj/`. The `.meta.json`
of each run records parameters; the design doc records *intent and
predictions*, which metadata cannot capture. Exploration reads both.

## Skills

### 1. `design-experiment` (interactive)

The structured back-and-forth phase. Input: a hypothesis, an open question
from a Logseq page, or a follow-up suggestion from a prior results doc.

Steps sketch:

1. Gather context — the relevant Logseq pages, prior `docs/results/` docs in
   the same line of inquiry, and what `src/research` structures exist to
   compose. (Prior-results survey can be a scout subagent.)
2. Draft the design with the user via explicit checkpoints: hypothesis;
   falsification criterion; parameter space and sweep granularity; controls /
   baselines; replication count and seed policy; expected DataFrame schema
   (one row per what? which columns answer the question?).
3. **Pre-register predictions**: for each sweep region, what does theory or
   literature predict, with the Logseq page or paper cited. Explicitly mark
   regions where no prediction exists — those are the interesting ones.
4. Write `docs/experiments/<name>.md` and stub the open question link back
   into Logseq.

The skill's value is a reference checklist of experiment-design failure modes
(confounded sweeps, seed reuse across arms, grids too coarse to see the
predicted effect, schema that can't support the planned analysis), applied
conversationally rather than as a gate.

### 2. `add-experiment`

Implements `src/experiments/<name>.py` from the design doc, following the
`@experiment` decorator conventions (`cli.py`): standard options, own click
parameters for the swept axes, `write_result` / `write_metadata`, seeds
threaded explicitly, metadata sufficient to reproduce the run (all parameters,
git commit, relevant library versions).

Steps sketch: read design doc → survey `src/experiments` conventions →
implement → smoke test → **pilot run** (tiny grid, verifies schema matches
the design and estimates runtime/variance so the user can size the real
sweep) → audit gate (`experiment-auditor`, below) → report, including the
suggested full-sweep invocation for the user to run.

The pilot run writes to the scratchpad via `--out`/`--meta-out` — never to
`data/results/`, where `explore-results` would later resolve it as a real
run. The full sweep itself stays user-invoked — runs cost real time and the
user decides when/what to run.

### 3. `explore-results` (mostly autonomous)

Input: a `data/results/<name>.<ts>.json` path (or just the experiment name,
resolving to the latest timestamp). Output: `docs/results/<name>.<ts>.md`.

Steps sketch:

1. Read the run's `.meta.json` and the experiment's design doc; extract the
   pre-registered predictions.
2. **Profile the data** via the `data-profiler` agent: integrity findings
   (row count vs. expected grid size, NaN/inf, dtype surprises, parameter
   coverage holes) plus a structured summary. Integrity failures stop the
   skill — a broken run gets reported, not interpreted.
3. Scripted analysis in the scratchpad (pandas/matplotlib via `uv run`):
   trends and fits per swept axis, invariant checks, outlier isolation,
   comparison of each pre-registered prediction against the data —
   classifying every prediction as *confirmed / contradicted / no prediction
   existed*.
4. **Literature check**: any *contradicted* finding goes to a
   `literature-checker` pass against Logseq + `docs/ref` before it may be
   written up as a finding rather than a suspected bug.
5. Draft the results doc: overview, per-prediction outcomes, unexpected
   observations, suspected-bug list, and falsification-oriented follow-up
   experiment suggestions. Steps 3–5 apply a shared "skeptic moves"
   reference (alternative explanations, bug-first triage, what-would-
   falsify-this) that `critique-results` reuses.
6. Audit gate (`results-auditor`, below); loop until `CLEAN`.
7. Update `docs/results/index.md` and close the loop in Logseq.

### 4. `critique-results` (interactive — the Feedback phase)

A walkthrough companion rather than a pipeline: load a results doc plus its
data and design doc, then work through it with the user in an explicitly
skeptical stance — recompute claims on demand, steelman alternative
explanations, hunt for literature contradictions the exploration pass missed,
and press on whether each follow-up suggestion could actually falsify the
forming hypothesis. Ends by recording the refined interpretation: edits to
the results doc and the emerging informal conjecture noted in Logseq (feeding
the later conjecture phase).

Could reuse the `literature-checker` agent and the "skeptic moves"
reference shared with `explore-results` steps 3–5.

## Agents

### `experiment-auditor` (gates `add-experiment`)

Fresh-eyes audit of the script against the design doc — the defense against
"contradicts literature because the experiment is buggy", applied at the
source. Given paths to `docs/experiments/<name>.md` and the script; checks:

- **[sweep]** the code actually varies what the design says it varies, over
  the stated ranges, and nothing else varies with it (confounds).
- **[seed]** randomness is seeded per the design's replication policy; no
  shared RNG state leaking across arms; seeds recorded in the output.
- **[schema]** the emitted DataFrame supports every analysis the design
  plans; units and column semantics match.
- **[meta]** metadata captures everything needed to reproduce the run.
- **[model]** the `src/research` calls match the semantics the design
  assumes (right operation, right convention/normalization).

Findings-list output, severity-ordered, `CLEAN` sentinel — same contract as
the existing auditors.

### `data-profiler` (scout; used by `explore-results`)

The `source-reader` analog for DataFrames: loads the result JSON with pandas
via Bash, returns a structured profile and integrity findings, never raw
rows beyond a handful of illustrative examples. Keeps hundreds of kilobytes
of JSON out of the orchestrating context.

### `results-auditor` (gates `explore-results`)

The `draft-auditor` analog for results docs. Given the results doc, the data
path, and the design doc; checks:

- **[num]** every numeric claim in the doc is recomputable from the
  DataFrame (recompute via pandas, not by eye).
- **[drift]** claim strength: "suggests" vs. "shows" vs. "proves"; a trend
  in a coarse sweep stated as a law.
- **[missing]** findings the data clearly supports that the doc omits —
  especially contradicted predictions quietly left out.
- **[prereg]** every pre-registered prediction from the design doc is
  addressed with an explicit outcome.
- **[bug-first]** surprising findings presented without the
  literature-check / suspected-bug triage.

### `literature-checker` (used by `explore-results` and `critique-results`)

Given a finding, search Logseq and `docs/ref` for claims it bears on; return
per-claim dispositions (consistent / contradicts / no coverage) with page and
paper citations. `mathlib-scout`-style: verified citations, no raw dumps.

## Build order

1. `explore-results` + `data-profiler` + `results-auditor` — highest leverage,
   most autonomous, and useful immediately for existing result files; the
   README's Exploration phase is currently entirely manual. Both the skill
   and the auditor degrade gracefully at this phase: with no design doc,
   `explore-results` skips prediction comparison and the auditor skips
   `[prereg]`; `[bug-first]` checks against the literature step inlined in
   the skill (see item 4) until the dedicated checker exists.
2. `design-experiment` + the `docs/experiments/` convention — establishes
   pre-registration, which `explore-results` reads in step 1 and compares
   against the data in step 3.
3. `add-experiment` + `experiment-auditor`.
4. `critique-results` + `literature-checker` (the checker can start as a step
   inside `explore-results` and be extracted once both skills need it).

## Open questions

- Should the design doc live in `docs/experiments/` (new taxonomy entry) or as a
  structured header inside the experiment module's docstring? A separate doc
  keeps predictions readable and diffable; the docstring keeps design and code
  from drifting apart. Leaning separate doc + an `experiment-auditor` check that
  they agree. - DECISION: I like the new directory.

- Multi-run exploration: results docs are per-timestamp, but interpretation
  often compares runs (before/after a parameter change). A `compare-runs` mode
  of `explore-results` (given two timestamps, diff parameters and outcome
  deltas) may be worth a step rather than a separate skill. - DECISION: I'll
  leave that up to you. → Resolved: built as a compare mode of
  `explore-results` — a second, baseline timestamp adds a parameter-diff and
  outcome-delta comparison section to the *newer* run's doc, keeping the
  per-timestamp naming convention.

- How much figure generation belongs in exploration? Scratch plots are analysis
  tools; anything worth keeping should graduate to `src/figures` via the
  existing figure process rather than being written ad hoc into `docs/`. -
  DECISION: I don't think explicit figure generation belongs in the exploration
  step.
