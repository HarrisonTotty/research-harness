---
name: design-experiment
description: Designs a parameterized experiment interactively with the user — hypothesis, falsification criterion, parameter space, controls, seed policy, output schema — and records it with pre-registered predictions as a design doc in docs/experiments/. Use when the user wants to design, plan, or spec an experiment from a hypothesis, an open question, or a follow-up suggested by a prior results doc.
argument-hint: [hypothesis or question]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Design a pre-registered experiment

Design an experiment answering **$ARGUMENTS** together with the user, and
record it as `docs/experiments/<name>.md`. Design is inherently
interactive: this skill structures the back-and-forth, it does not replace
it — every design element lands at an explicit checkpoint with the user
before the next one starts. The deliverable is a design doc complete
enough that `/add-experiment` can implement it without further design
decisions, with predictions recorded *before* the experiment ever runs.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Gather context
- [ ] Step 2: Draft the design with the user
- [ ] Step 3: Pre-register predictions
- [ ] Step 4: Write docs/experiments/<name>.md
- [ ] Step 5: Close the loop in Logseq and report
```

## Working files

Before gathering context, create two files in the session scratchpad
directory:

- `context.md` — the distilled inputs: the hypothesis and its provenance,
  the relevant Logseq claims, prior results in this line of inquiry, and
  the `src/research` API available to compose. Step 1 writes it; Steps
  2–3 read from it.
- `design.md` — the design as it accretes through checkpoints. Steps 2–3
  write it; Step 4 publishes it.

These files are the ground truth that survives context compaction.
Wherever conversation memory and file content disagree, trust the files.

## Step 1: Gather context

Collect and distill into `context.md`:

- **The Logseq pages behind the hypothesis** (`get_page`; fall back to
  `search_logseq` with singular variants): their claims, theorems, and
  open questions. These are the future sources for Step 3's predictions.
- **Prior `docs/results/` docs in the same line of inquiry** — especially
  their follow-up suggestions and suspected-bug lists, so a new
  experiment does not repeat a suspect run's mistake. If the survey is
  large, dispatch an `Explore` subagent to return a distilled summary
  rather than reading every doc into this conversation.
- **What `src/research` provides to compose.** Experiments compose
  existing library structures; if the hypothesis needs a structure the
  library lacks, stop and tell the user to run `/add-python-topic` first
  — do not plan an experiment around code that does not exist.

## Step 2: Draft the design with the user

Work through the design elements as explicit checkpoints: propose a
concrete choice, let the user push back, and record the agreement in
`design.md` before moving on. The elements, in order:

1. **Hypothesis** — one falsifiable statement.
2. **Falsification criterion** — the specific outcome that would disprove
   it. If no outcome could, return to 1.
3. **Parameter space** — axes, ranges, granularity. Granularity is chosen
   relative to the predicted effect, not to a round number.
4. **Controls and baselines** — the arm the effect is measured against.
5. **Replication and seeds** — replication count and whether arms share
   seeds (paired) or draw independently, as a decision rather than an
   accident.
6. **Output schema** — one row per what, and which columns answer the
   question. Walk the intended analysis against the columns: an analysis
   the schema cannot support is a design bug caught cheapest here.

Apply [reference/failure-modes.md](reference/failure-modes.md) as each
element lands — raise the failure mode as a question in the conversation,
conversationally rather than as a gate.

## Step 3: Pre-register predictions

For each region of the parameter space, record what theory or the
literature predicts, citing the Logseq page or paper from `context.md`
that says so — if a needed claim was not gathered, go check it now rather
than predicting from memory. Explicitly mark regions where **no
prediction exists**: those are the interesting ones, and marking them now
is what makes "unexpected" well-defined at exploration time. Predictions
are frozen once the first run exists; later insight goes in a dated
addendum, never as an edit to the originals.

## Step 4: Write `docs/experiments/<name>.md`

Publish `design.md` as `docs/experiments/<name>.md`, following the
template in [reference/design-doc.md](reference/design-doc.md). The doc
name must match the future module name (`<name>` hyphenated in the doc
and CLI, underscores in `src/experiments/<name>.py`) — check that neither
the doc nor the module already exists before writing.

## Step 5: Close the loop in Logseq and report

Add a pointer to the design doc on the Logseq page (or open-question
block) the hypothesis came from, respecting the MCP sharp edges in
[../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md).
Then report to the user: the doc path, the open risks the design accepts
(granularity limits, unpredicted regions), and the handoff — the next
step is `/add-experiment <name>`.
