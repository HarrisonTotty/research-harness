---
name: add-experiment
description: Implements a designed experiment as src/experiments/<name>.py following the @experiment harness conventions, pilots it on a tiny grid into the scratchpad, and gates it on a fresh-eyes audit against the design doc. Use when the user asks to implement, code up, or build an experiment that has a design doc in docs/experiments/.
argument-hint: [experiment]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash
---

# Implement a designed experiment

Implement **$ARGUMENTS** as `src/experiments/<name>.py`. The design doc
`docs/experiments/<name>.md` is the specification: its parameter space
fixes the CLI options, its schema fixes the emitted DataFrame, its seed
policy fixes the randomness plumbing, and its analysis plan is what the
output must be able to answer. If the doc is missing, or lacks a
falsification criterion, parameter space, output schema, or predictions,
stop and tell the user to run `/design-experiment` first — do not design
on the fly; the doc is the specification.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Read the design doc
- [ ] Step 2: Survey the harness and research API
- [ ] Step 3: Implement the module
- [ ] Step 4: Pilot run into the scratchpad
- [ ] Step 5: Audit against the design
- [ ] Step 6: Run the full gate
- [ ] Step 7: Update the design doc and report
```

## Working files

Create a `pilot/` directory in the session scratchpad before Step 4.
Every pilot artifact goes there — **never into `data/results/`**, where
`/explore-results` would later resolve a pilot as the newest real run.

## Step 1: Read the design doc

Read `docs/experiments/<name>.md` in full. Note every analysis in the
analysis plan — the schema check in Step 4 and the audit in Step 5 are
made against what those analyses need, not against what looks reasonable.

## Step 2: Survey the harness and research API

- The harness contract in `src/experiments/cli.py` and `io.py`: the
  `@experiment` decorator supplies `--log-level`, `--log-file`, `--out`,
  and `--meta-out` and hands the body an `ExperimentContext`; the module
  must define exactly one command (the dispatcher in `__main__.py`
  rejects zero or several); the CLI name hyphenates what the module name
  spells with underscores.
- Existing `src/experiments/` modules, as the style to extend.
- The `src/research` structures the design composes — read the actual
  signatures and docstrings of every call you plan to make. The design's
  semantic assumptions about them are audited in Step 5; implement from
  source, not from memory of the design's phrasing.

## Step 3: Implement the module

Write `src/experiments/<name>.py`: one `@experiment`-decorated function,
with its own click options for the swept axes **defaulting to the
design's full ranges**, so the bare `just experiment <name>` runs the
designed sweep and the pilot overrides downward. House Python rules
apply; the module docstring cites the design doc. Skill-specific
obligations:

- Seeds threaded explicitly per the design's policy (a base seed option,
  per-arm derivation as designed), and every seed emitted into the
  result or metadata.
- `write_metadata` records everything needed to reproduce the run: every
  parameter (including defaults), all seeds, the git commit, and
  relevant library versions.
- `write_result` emits one row per the design's stated unit, with the
  design's exact column names and units.

## Step 4: Pilot run into the scratchpad

Run a tiny grid (a few cells, minimal replication) with `--out` and
`--meta-out` pointed into the scratchpad `pilot/` directory. Verify
against the design, reading the pilot files with a short
`uv run python` script:

- schema exact: columns, dtypes, one row per unit, row count = cells ×
  replications;
- metadata complete: parameters, seeds, commit;
- runtime per cell, extrapolated to the full sweep, and observed
  variance across replications — the user sizes the real sweep from
  these.

A discrepancy that traces to the design itself (a schema that cannot
support a planned analysis, a range that explodes runtime) goes back to
the user — do not silently redesign.

## Step 5: Audit against the design

Do not hand over an unaudited experiment: the script and the pilot
verification came from the same reading of the design, and a buggy
experiment is exactly how "data that contradicts the literature"
happens. Dispatch an `experiment-auditor` subagent with the paths to the
design doc and the script; it returns a severity-ordered findings list,
or `CLEAN`. A finding usually means editing the script to match the
design; a finding that exposes a design defect is reported to the user,
not silently patched into either file. Re-dispatch until `CLEAN` — only
a clean audit unlocks Step 6.

## Step 6: Run the full gate

Run `just check` and fix findings until it passes clean.

## Step 7: Update the design doc and report

Fill the design doc's **Implementation notes** section: module path,
pilot observations (runtime per cell, observed variance, any granularity
or replication concern the pilot raised), and the exact suggested
full-sweep invocation. Then report to the user: what was built, the
pilot numbers, and the invocation — `just experiment <name> [args ...]`.
The full sweep stays user-invoked: runs cost real time, so report the
command and stop rather than running it.
