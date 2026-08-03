---
name: add-conjecture
description: Formalizes an informal conjecture from Logseq into an audited docs/conj/<name>.md plus a committed Lean Prop def — literature-checked first, statement drafted with the user at explicit checkpoints, stress-tested by counterexample search beyond the experimental ranges, and gated on fresh-eyes audits. Use when the user wants to state, register, or formalize a conjecture from an informal note or a line of inquiry.
argument-hint: [conjecture or line of inquiry]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Formalize a conjecture

Turn the informal conjecture behind **$ARGUMENTS** into an audited
`docs/conj/<name>.md` and a committed Lean statement. Choosing the right
statement is a judgment call, so this skill is interactive — it
structures the session and holds the discipline. Two hard rules
throughout:

- **Falsify before you formalize.** The statement gets hammered
  computationally before any Lean or doc effort is spent on it; a
  counterexample here is cheap and is itself a finding.
- **A conjecture is a `Prop`, not a `sorry`.** The statement is
  committed as a `def ... : Prop` — it elaborates and type-checks with
  no proof obligation, and the theorem phase later proves the literal
  def, so the proposition audited today is the term proven then.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Gather
- [ ] Step 2: Check the literature and Mathlib
- [ ] Step 3: Draft the formal statement with the user
- [ ] Step 4: Stress-test the statement
- [ ] Step 5: State it in Lean
- [ ] Step 6: Write docs/conj/<name>.md
- [ ] Step 7: Close the loop in Logseq and report
```

## Working files

Before gathering, create in the session scratchpad directory:

- `conjecture.md` — the skill's single handoff artifact, accreting one
  section per step: candidate statements and the evidence map (Step 1);
  literature and Mathlib dispositions (Step 2); the agreed statement
  verbatim and its vocabulary map (Step 3); stress-test coverage and
  the revision log (Step 4). Every auditor dispatch reads this file —
  nothing is handed to a subagent from conversational memory.
- `hunt/` — the counterexample-search harness. It persists across the
  refute→refine loop's dispatches; each hunter is told what is already
  there.

These files are the ground truth that survives context compaction.
Wherever conversation memory and file content disagree, trust the files.

## Step 1: Gather

Collect and distill into `conjecture.md`:

- **The informal-conjecture note** (`get_page`; fall back to
  `search_logseq`) — usually written by `/critique-results` on a
  line-of-inquiry page. Quote the claim verbatim. If no note exists and
  the user is stating the conjecture directly, record their wording and
  confirm it back before proceeding.
- **The results docs the note links**, and their design docs: the
  evidence map is built here — which experiments, over which ranges
  *actually swept*, found what. Quote ranges from the docs, never from
  memory; the covered ranges also parameterize Step 4's search.
- **The relevant topic pages** — the vocabulary the statement will be
  expressed in, and any theorems that already bear on it.

## Step 2: Check the literature and Mathlib

An informal conjecture that is already settled is not a conjecture.
Dispatch in parallel, and append both outputs to `conjecture.md`:

- `literature-checker` with the claim (and any distinct candidate
  statements) — is it known, contradicted, or uncovered?
- `mathlib-scout` — does Mathlib already state or prove it, or its
  natural generalization?

Either disposition reroutes the session, and the check runs first
because it does: **known-true** is intake work — tell the user to
`/add-logseq-topic` the result and `/add-lean-topic` it; **known-false**
is a bug hunt in the experiments that suggested it — report it, since
the fix belongs in the experiment pipeline, not here. Both reroutes stop
the skill.

## Step 3: Draft the formal statement with the user

Work through the statement at explicit checkpoints, applying
[reference/statement-craft.md](reference/statement-craft.md)
conversationally as each element lands:

1. **Quantifiers and side conditions** — every one explicit, including
   the ones the experiments' generators imposed silently.
2. **The generality dial** — the *minimal* statement the evidence
   supports vs. the *natural* general form. When they differ, state
   both: the general form as the conjecture, the minimal one as its
   implied special case.
3. **Vocabulary** — which `src/research` and Mathlib names the statement
   is expressed in, taken from the Step 2 scout dispositions.

Record the agreed statement verbatim and its vocabulary map in
`conjecture.md` before moving on — Step 5's audit is meaningless
without them.

## Step 4: Stress-test the statement

Dispatch `counterexample-hunter` with the statement, its vocabulary, the
covered ranges from the evidence map, and what is already in `hunt/` —
instructed to search *beyond* the covered ranges (the experiments
motivated the conjecture; they cannot also be its stress test).

- **REFUTED** — the heart of the skill: refine the statement with the
  user, log the revision (counterexample and the change it forced) in
  `conjecture.md`, and re-dispatch. A refinement that changes the
  statement's substance re-runs Step 2 first — the original dispatches
  answered a different claim.
- **SURVIVED** — the coverage report goes in `conjecture.md`'s evidence
  section. A "survived" without a coverage report is no answer;
  re-dispatch.

## Step 5: State it in Lean

Write `src/theorems/Theorems/Conjectures/<Name>.lean`: copyright header,
module docstring, minimal imports, then the `Prop` def — UpperCamelCase
(it is a Prop-valued definition), named for the statement, in the
namespace of its principal subject, with a docstring citing the evidence
trail. Add theorems for special cases that are already provable (each
one both strengthens the evidence and seeds the attack; no `sorry` — an
unprovable special case is just not stated). Add the module to
`src/theorems/Theorems.lean` and confirm it elaborates with
`lake build Theorems.Conjectures.<Name>`; check the unfolded statement
with `plausible` in a scratch `example` (not committed).

Gate on `statement-auditor`, with `conjecture.md` playing the role
`coverage.md` plays in `/add-lean-topic` — it must already carry the
agreed statement verbatim and the vocabulary map, or the auditor's
`[reuse]` check has nothing to bite on. Tell the dispatch the artifact
under audit is a `def ... : Prop` in the Mathlib open-problem style, so
the shape is not itself read as a finding. Fix every finding and
re-dispatch until `CLEAN`, then run `just lean-check`.

## Step 6: Write `docs/conj/<name>.md`

Write the doc from `conjecture.md` following
[reference/conjecture-doc.md](reference/conjecture-doc.md): formal
statement (prose plus the Lean name), evidence map, stress-test
coverage, known special cases, an initial attack-plan skeleton
(candidate strategies, known obstructions), and the revision log.

Gate on `conjecture-auditor` with the paths to the doc, `conjecture.md`,
the Lean module, and the cited results docs. Fix every finding and
re-dispatch until `CLEAN`.

## Step 7: Close the loop in Logseq and report

Add the conjecture's line to `docs/conj/index.md`. Create the
conjecture's Logseq page and link it from the line-of-inquiry note that
seeded it and from the relevant topic pages, respecting the MCP sharp
edges in
[../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md).

Then report to the user: the statement as registered (prose and Lean
name), the stress-test coverage and where it ends, special cases already
proven, the revision history (what counterexamples forced), and the
handoff — the next step is `/attack-conjecture <name>`.
