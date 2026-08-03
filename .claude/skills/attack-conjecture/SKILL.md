---
name: attack-conjecture
description: Works a proof session against a docs/conj/ conjecture — the lemma DAG planned with the user, proved Lean-first with audited statements and no committed sorry, every status change written through to the doc as it happens, and a closed proof promoted to docs/theorems/ behind a fresh-eyes proof audit. Use when the user wants to attack, prove, or continue proving a registered conjecture.
argument-hint: [conjecture]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Attack a conjecture

Work a proof session against `docs/conj/`'s page for **$ARGUMENTS**. A
session rarely ends with a full proof — the deliverable is *recorded
progress*: lemmas closed, obstructions hit, the plan updated so the next
session starts where this one stopped. Proof work is a dialogue; the
skill holds the discipline:

- **Lean scaffolds the chalkboard.** Each lemma runs Lean-first; the
  prose proof at promotion is written *from* the closed Lean proof —
  two proofs written in parallel share nothing, and the prose one
  silently proves a nicer statement.
- **No `sorry` is ever left behind.** An unproven lemma's statement
  returns to the plan file at session end, not the Lean tree.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Load the conjecture and prior attack state
- [ ] Step 2: Plan the lemma DAG with the user
- [ ] Step 3: Prove, Lean-first per lemma
- [ ] Step 4: Record the session
- [ ] Step 5: Promote on closure
```

## Working file

Before planning, create `attack.md` in the session scratchpad directory:
the lemma DAG with per-lemma entries in the format of
[reference/attack-plan.md](reference/attack-plan.md), rebuilt from the
doc's attack-plan section — the scratchpad is session-scoped; the doc is
the durable copy. Wherever conversation memory and the file disagree,
trust the file; wherever the file and the doc disagree, trust the doc.

## Step 1: Load the conjecture and prior attack state

Read `docs/conj/<name>.md`, its Lean module
(`src/theorems/Theorems/Conjectures/<Name>.lean`), and the doc's
attack-plan section into `attack.md`. Confirm the pinned Mathlib
checkout exists at `.lake/packages/mathlib/` (`just lean-update` if
not). Two special states short-circuit:

- The page carries a **Refuted** banner: there is nothing to attack —
  point the user at the successor statement the banner names.
- The main theorem is already proven but unpromoted (a prior session
  closed it at the end of its budget): jump straight to Step 5, so
  promotion's writing and audit run on a fresh context instead of the
  tail of an exhausted one.

## Step 2: Plan the lemma DAG with the user

Decompose toward the goal with the user: candidate lemma chains,
recorded as DAG entries per
[reference/attack-plan.md](reference/attack-plan.md). Chalkboard
sketches live here as the plan's prose — the README's "two angles" are
two views of one DAG, not two plans. Per planning round:

- **One `mathlib-scout` dispatch covering all candidate lemmas** — its
  contract is per-claim dispositions, so a batch is its native input;
  per-lemma dispatches just multiply round-trips. Each lemma's entry
  records the vocabulary its dispatch found.
- `literature-checker` when an obstruction smells known — a published
  counterexample or a known hardness result is cheaper to read than to
  rediscover.
- **Risk-first ordering**: attack the load-bearing uncertain lemma
  first — if it falls, the plan falls cheaply.

Write the agreed plan through to the doc's attack-plan section before
proving starts.

## Step 3: Prove, Lean-first per lemma

Per lemma, in plan order:

1. **State** it with `sorry` in the working tree, placed by the
   lemma-placement test in
   [reference/attack-plan.md](reference/attack-plan.md) — topic module
   if it would survive the conjecture being forgotten, conjecture
   module otherwise.
2. **Audit** new statements with `statement-auditor`, batched per
   planning round rather than one dispatch per lemma — here `attack.md`
   plays the coverage-map role: each entry carries the lemma's informal
   statement, its role in the DAG, and its scout vocabulary, which is
   already the multi-claim shape that auditor expects. `CLEAN` before
   proof effort.
3. **Falsify before hard effort**: `plausible` against the statement; a
   `counterexample-hunter` dispatch for anything with sweepable
   structure.
4. **Prove**, under the per-lemma effort cap (two serious tactic
   strategies, then record as obstructed and move on — sessions end
   with breadth, not one heroic stuck proof). The `/add-lean-topic`
   stuck-proof triage applies: a stuck proof is a wrong statement, a
   genuinely hard proof, or a false claim — diagnose before pushing.

A refutation of the *conjecture itself* mid-attack ends the session into
`/add-conjecture`'s refinement loop — and since an attack has now
accumulated, the refinement is a new page, the old one marked refuted →
superseded. At session end every remaining `sorry` is deleted: the
statement returns to the plan file as its entry's Lean text, never a
committed `sorry`.

## Step 4: Record the session

Not a sync at session end — **write each lemma's status change through
to the doc's attack-plan section as it happens**: `attack.md` dies with
the session, and an attack session is exactly the kind that ends
abruptly, context exhausted mid-proof. This step is the closing pass:

- Obstruction records completed — what was tried, why it failed, what
  the next session should try first.
- `just lean-check` on whatever landed; every warning is a defect.
- A one-line session note on the conjecture's Logseq page (lemmas
  closed, obstructions), respecting the MCP sharp edges in
  [../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md).

Then report: lemmas closed this session, obstructions with their why,
and the recommended first move for the next session.

## Step 5: Promote on closure

When the main theorem proves the `Prop` def — `theorem` in the topic
module `src/theorems/Theorems/<Topic>.lean`, stating the def from
`Conjectures/<Name>.lean` literally (the def stays put; its location
encodes its history):

1. **Axiom audit** per
   [../add-lean-topic/reference/proving.md](../add-lean-topic/reference/proving.md):
   nothing beyond `propext`, `Classical.choice`, `Quot.sound`.
2. **Write `docs/theorems/<name>.md`** per
   [reference/theorem-doc.md](reference/theorem-doc.md) — the human
   proof written *from* the Lean proof, at prose altitude, naming the
   Lean declaration.
3. **Gate on `proof-auditor`** with the page and module paths; fix
   every finding and re-dispatch until `CLEAN`.
4. **Move the page**: delete `docs/conj/<name>.md`, move its index line
   from `docs/conj/index.md` to `docs/theorems/index.md`, and update
   the conjecture's Logseq page (status → proven, linking the theorem
   page).

Then report: the theorem as promoted, its axiom audit, where the def and
theorem live, and any lemmas that graduated to topic modules along the
way.
