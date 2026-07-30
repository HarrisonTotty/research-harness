---
name: explore-results
description: Explores an experiment run's raw results into an audited docs/results/<name>.<timestamp>.md — data profiled through a subagent, analysis scripted in the scratchpad, pre-registered predictions classified, surprises literature-checked bug-first, and the draft gated on a fresh-eyes audit. Given a second timestamp, also compares the two runs. Use when the user asks to explore, analyze, interpret, or write up experiment results.
argument-hint: [experiment[.timestamp]] [baseline-timestamp]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Explore an experiment run into a results doc

Turn the raw artifacts of run **$ARGUMENTS**
(`data/results/<name>.<ts>.json` + `.meta.json`) into
`docs/results/<name>.<ts>.md`. This phase is nearly autonomous — the
interactive scrutiny comes afterward via `/critique-results`. Two hard
rules throughout:

- **Raw DataFrames never enter this conversation.** Profiling goes
  through the `data-profiler` subagent; analysis happens in scripts that
  print summaries.
- **Contradiction with recorded literature is a bug until proven
  otherwise.** No "surprising finding" reaches the doc without its
  triage trail.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Resolve the run and its design doc
- [ ] Step 2: Profile the data
- [ ] Step 3: Scripted analysis
- [ ] Step 4: Literature-check the surprises
- [ ] Step 5: Draft the results doc
- [ ] Step 6: Audit the draft
- [ ] Step 7: Publish and close the loop
```

## Working files

Before resolving the run, create in the session scratchpad directory:

- `run.md` — resolved paths, metadata summary, and the design doc's
  extracted predictions. Step 1 writes it; every later step reads it.
- `profile.md` — the profiler's integrity findings and data profile.
  Step 2 writes it; Steps 3 and 5 read.
- `analysis/` — the pandas scripts and their printed outputs, plus
  `findings.md`, the accumulating structured findings. Steps 3–4 write;
  Step 5 reads.
- `draft.md` — the results doc draft. Step 5 writes it; Steps 6–7 read.

These files are the ground truth that survives context compaction.
Wherever conversation memory and file content disagree, trust the files.

## Step 1: Resolve the run and its design doc

Given only an experiment name, resolve the newest timestamp in
`data/results/` (the UTC stamp sorts lexicographically). Read the run's
`.meta.json` directly — it is small — and read the design doc
`docs/experiments/<name>.md`, extracting the pre-registered predictions
and analysis plan into `run.md`. If no design doc exists, note that in
`run.md` and degrade gracefully: Steps 3 and 5 skip prediction
bookkeeping, and the doc's prediction section states that the run was
not pre-registered.

**Compare mode:** a second timestamp names a baseline run of the same
experiment. Resolve both pairs of artifacts; the comparison rides along
as an extra analysis in Step 3 and an extra section in Step 5, in the
*newer* run's doc.

## Step 2: Profile the data

Dispatch a `data-profiler` subagent with the result and metadata paths
plus the expected grid, replication count, and schema from the design
doc; write its output to `profile.md`. (In compare mode, dispatch one
per run, in parallel; both outputs go in `profile.md`, each under a
heading naming its run, the baseline labeled as such.)

Any integrity finding the design doc does not explain **stops the
skill**: report the broken run to the user with the profiler's findings
and the suspected cause. A broken run gets reported, not interpreted.

## Step 3: Scripted analysis

Write pandas scripts under `analysis/` and run them from the repository
root with `uv run python`; scripts print aggregates and conclusions,
never raw frames. Cover, in order:

1. The design doc's analysis plan.
2. Trends and fits per swept axis; invariant checks the design or the
   model implies; outlier isolation (which cells, how far).
3. Every pre-registered prediction compared against the data and
   classified **confirmed / contradicted / no prediction existed**, the
   recomputed numbers recorded beside each in `findings.md`.
4. Compare mode: a parameter diff of the two `.meta.json` files
   (identical / changed / added / removed), then outcome deltas on
   shared measures over the shared grid.

Scratch plots are analysis tools: they live and die in `analysis/`, and
no figure is written into `docs/` — anything worth keeping graduates to
`src/figures` through the figure process.

## Step 4: Literature-check the surprises

Batch every contradicted prediction and unexpected observation into one
`literature-checker` dispatch. Then triage by disposition, recording
each in `findings.md`:

- **[contradicts]** — bug-first: re-examine the experiment for the
  specific error class that would produce the contradiction (see the
  bug-first moves in
  [reference/skeptic-moves.md](reference/skeptic-moves.md)). A bug found
  invalidates the run: stop and report it — the fix belongs in the
  experiment, not in prose. No bug found after honest triage: the
  finding may be written up, carrying the contradiction and the triage
  trail.
- **[consistent]** — the surprise was only a surprise relative to our
  predictions, not to the literature: record it with its supporting
  citation; no bug triage needed.
- **[no-coverage]** — stays a finding, flagged as uncharted and a prime
  candidate for follow-ups.

## Step 5: Draft the results doc

Write `draft.md` following
[reference/results-doc.md](reference/results-doc.md), applying
[reference/skeptic-moves.md](reference/skeptic-moves.md) to every claim
— alternative explanations considered, claim strength on the ladder,
follow-ups falsification-oriented. Every number is copied from an
`analysis/` output, never recalled from memory.

## Step 6: Audit the draft

Do not publish an unaudited draft: the draft and the analysis came from
the same reading of the data. Dispatch a `results-auditor` subagent with
the paths to `draft.md`, the result and metadata files (in compare mode,
both runs' pairs, labeling which is the baseline), and the design doc
(when it exists); it returns a severity-ordered findings list, or
`CLEAN`. A `[num]` finding usually means the draft misquotes the
analysis — recompute, then fix the draft. Re-dispatch until `CLEAN` —
only a clean audit unlocks Step 7.

## Step 7: Publish and close the loop

Install `draft.md` as `docs/results/<name>.<ts>.md` (the site nav picks
the page up automatically). Close the loop in Logseq: on the page the
experiment's hypothesis came from, note the run, the doc, and the
headline outcome, respecting the MCP sharp edges in
[../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md).
Then report to the user: the per-prediction outcomes, the surprises and
their triage, any suspected bugs, the follow-up suggestions, and the
handoff — the next step is `/critique-results <name>.<ts>`.
