---
name: critique-results
description: Walks through a results document with the user in an explicitly skeptical stance — recomputing claims on demand, steelmanning alternative explanations, hunting missed literature contradictions, and pressure-testing follow-ups — then records the refined interpretation in the doc and Logseq. Use when the user wants to review, critique, discuss, or push back on a results doc (the Feedback phase).
argument-hint: [experiment[.timestamp]]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Critique a results doc with the user

Walk through the results doc for **$ARGUMENTS** with the user. This is a
companion, not a pipeline: the user asks questions and pushes back; you
are the skeptic-in-residence. The deliverable is *refined
interpretation* — an edited results doc and an informal-conjecture note
in Logseq — not a report.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Setup: load the run's artifacts
- [ ] Walkthrough: the skeptical pass with the user
- [ ] Record: doc edits and the Logseq note
```

## Working files

Create a `recompute/` directory in the session scratchpad before the
walkthrough: every on-demand analysis script and its printed output goes
there. Raw DataFrames stay out of this conversation here too —
recomputation happens in scripts (run from the repository root via
`uv run python`) that print what they conclude.

## Setup

Read the results doc `docs/results/<name>.<ts>.md`, the design doc
`docs/experiments/<name>.md` (when it exists), the run's `.meta.json`
(small), and any prior results docs the doc links. Resolve — but do not
read — the result file `data/results/<name>.<ts>.json`: its path feeds
the recompute scripts.

Open by offering an agenda — the doc's strongest claim, its
contradicted or unexpected findings, and its follow-ups — but follow the
user's lead.

## Walkthrough

The stance, held throughout
([../explore-results/reference/skeptic-moves.md](../explore-results/reference/skeptic-moves.md)
is the repertoire):

- **Recompute on demand.** Any number the user questions is recomputed
  from the raw data with a fresh script — never defended from the doc.
- **Steelman alternatives.** When the user proposes another explanation,
  argue it as well as it can be argued, then identify what existing or
  new data would discriminate between it and the doc's interpretation.
- **Hunt missed contradictions.** For claims whose doc trail lacks a
  literature disposition — or when the user raises a connection the
  exploration pass never checked — dispatch a `literature-checker`
  subagent rather than asserting from memory. A `[contradicts]` that
  survives scoping reopens bug-first triage, even this late.
- **Pressure-test follow-ups.** For each suggested follow-up, press on
  whether it could actually falsify the forming hypothesis; rank by
  discriminating power, and prefer replacing a confirmation-only
  suggestion over keeping it.
- **Honor pre-registration.** When interpretation drifts toward a region
  the design marked "no prediction", say so — that is exploration, not
  confirmation, and the refined text should label it.

## Record

Close the loop before ending — agreement that never lands in a file is
lost:

- **Edit the results doc** with everything agreed during the walkthrough:
  claim-strength adjustments, added caveats and alternatives, revised or
  replaced follow-ups. Leave contested points visible as open questions
  in the doc rather than resolving them by fiat.
- **Note the informal conjecture in Logseq** — the hypothesis as refined
  by this session, on the page the experiment's line of inquiry lives
  on, linking the results doc(s) feeding it. Respect the MCP sharp edges
  in
  [../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md).
  This note is the seed the conjecture phase later formalizes into
  `docs/conj/`.
- **Report**: what changed in the doc, the conjecture note as written,
  and the disagreements left open.
