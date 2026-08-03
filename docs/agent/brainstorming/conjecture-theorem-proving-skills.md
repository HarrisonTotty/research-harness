# Skills & agents for the Conjecture & Theorem Proving phase

Brainstorming for tooling that covers the back end of the research process:
turning the informal-conjecture notes that `/critique-results` accumulates in
Logseq into formal conjectures in `docs/conj/`, attacking them from the two
angles the README names — *via chalkboard* and *computationally via Lean* —
and promoting proven ones to `docs/theorems/` + `src/theorems/`. Like
experiment design, both halves are inherently interactive: choosing the right
statement is a judgment call, and proof work is a dialogue. The tooling should
structure the sessions, hold the discipline, and keep bulk output (Mathlib
source, counterexample sweeps) out of the main context.

## House patterns to reuse

All four established patterns carry forward directly:

1. **Checklist skills with scratchpad working files** — here a statement/
   evidence map and a lemma-DAG attack file.
2. **Fresh-eyes auditor agents** with the findings-list / `CLEAN` contract.
   The correlated-error argument is at its sharpest in this phase: a proved
   wrong statement is the most expensive artifact the repository can produce.
3. **Scout/isolator agents** — `mathlib-scout` reused as-is; one new isolator
   for counterexample search.
4. **Close the loop in Logseq** — the informal note that seeded the conjecture
   gets linked to the formal artifact, and status flows back on promotion or
   refutation.

Plus two disciplines from `/add-lean-topic`: **state before proving**
(statements audited by `statement-auditor` before any proof effort — the
discipline the README explicitly says carries over to this phase) and **no
committed `sorry`** (that skill's own rule, adopted here too).

## Design principles specific to this phase

- **A conjecture is a `Prop`, not a `sorry`.** The no-committed-`sorry` rule
  seems to conflict with formalizing an unproven statement — resolved the way
  Mathlib states open problems: the conjecture is committed as a `def ... :
  Prop`, which elaborates and type-checks with no proof obligation. The
  theorem phase later proves `theorem <name> : <conjecture_def>`. This buys
  the key property of the whole phase: **the proposition audited at
  conjecture time is the literal term proven at theorem time** — zero drift
  between what the evidence supported and what the proof establishes.
- **Falsify before you formalize, and again before you prove.** Proof effort
  is the most expensive resource in the pipeline. Before it is spent, the
  statement gets hammered computationally: a Python counterexample search
  over ranges *beyond* what the experiments covered (the experiments
  motivated the conjecture; they cannot also be its stress test), plus
  `plausible` on the Lean statement. A counterexample here is cheap and is
  itself a finding — it forces a refinement that makes the conjecture better.
- **Check the literature before claiming a conjecture.** An informal
  conjecture that is already a known theorem is not a conjecture — it is
  intake work (`/add-logseq-topic` the result, then `/add-lean-topic` it).
  Already-refuted is a bug hunt in the experiments that suggested it. Either
  disposition reroutes the session, so the check comes first.
- **Evidence is quoted, not remembered.** The conjecture doc's evidence
  section cites specific results docs with the ranges actually swept. The
  gap between tested range and claimed scope is the point of a conjecture —
  but it must be *visible*, never blurred ("verified for n ≤ 12" stated as
  "holds computationally").
- **Lean scaffolds the chalkboard.** Per the README: Lean proofs are ugly but
  safe. The attack runs Lean-first per lemma; the prose proof is written from
  the closed Lean proof, not in parallel from scratch — otherwise the two
  proofs share nothing and the prose one silently proves a nicer statement.

## Proposed artifact chain

```
informal conjecture note (Logseq, written by /critique-results)
  → docs/conj/<name>.md                      (formal statement + evidence map
                                              + attack plan; statement audited)
  → src/theorems/Theorems/Conjectures/<Name>.lean
                                             (the Prop def + any already-proven
                                              special cases)
  → attack sessions                          (lemma DAG worked Lean-first;
                                              progress/obstructions recorded
                                              back into the attack plan)
  → on closure:
      src/theorems/Theorems/<Topic>.lean     (theorem proving the Prop def)
      docs/theorems/<name>.md                (human proof, audited against Lean)
      docs/conj/<name>.md removed            (page "moves", per the indexes)
  → on refutation:
      the counterexample recorded; statement refined (revision logged
      in place pre-attack; a new page once an attack has accumulated) or
      the page kept in docs/conj/ with a Refuted banner
```

Directory decision implied: `src/theorems/Theorems/Conjectures/` holds the
Prop defs (one module per conjecture), imported like any other module so
`just lean-check` gates them. On promotion the def *stays* — the proving
theorem references it, preserving the statement's audit trail — and the
module stays put too: `Conjectures/` is the permanent home of
statements-as-defs, with the proving theorem living in its topic module
(decided below).

## Skills

### 1. `add-conjecture` (interactive)

From an informal Logseq note (or a line of inquiry the user names) to an
audited `docs/conj/<name>.md` + committed Lean statement.

Steps sketch:

1. **Gather** — the informal-conjecture note, the results docs it links,
   their design docs, and the relevant topic pages. Working file
   `conjecture.md` in the scratchpad — the skill's single handoff artifact,
   accreting one section per step: candidate statements and the evidence map
   (this step); literature/Mathlib dispositions (step 2); the agreed
   statement verbatim and its vocabulary map (step 3); stress-test coverage
   and the revision log (step 4). Every auditor dispatch reads this file —
   nothing is handed to a subagent from conversational memory.
2. **Literature check first** — `literature-checker` on the informal claim,
   and `mathlib-scout` for whether Mathlib already states or proves it.
   Known-true or known-false reroutes the session (see principles).
3. **Draft the formal statement with the user** at explicit checkpoints:
   every quantifier and side condition; the generality dial (the *minimal*
   statement the evidence supports vs. the *natural* general form — state
   both when they differ, as a conjecture and its implied special case);
   which `src/research` / Mathlib vocabulary the statement is expressed in.
4. **Stress test** — dispatch `counterexample-hunter` with the statement and
   the ranges the experiments covered, instructed to search *beyond* them.
   Refuted → refine with the user and repeat (this loop is the heart of the
   skill); survived → the coverage report goes in the evidence section.
   A refinement that changes the statement's substance re-runs the step-2
   known-result check — the original dispatches answered a different claim.
5. **State in Lean** — the `Prop` def in `Theorems/Conjectures/<Name>.lean`,
   plus theorems for special cases that are already provable (each one both
   strengthens the evidence and seeds the attack). `plausible` against the
   unfolded statement. Gate on `statement-auditor`, with the working file
   `conjecture.md` playing the role `coverage.md` plays in `/add-lean-topic` —
   which means the working file must carry what that contract expects: the
   agreed statement verbatim and the vocabulary map to `src/research` /
   Mathlib declarations from the step-2 scout dispatches (its `[reuse]`
   check is meaningless without them). The dispatch must also say the
   artifact under audit is a `def ... : Prop` — the auditor's usual input is
   sorry-bodied theorems, and the open-problem shape shouldn't read as a
   finding.
6. **Write the doc** — `docs/conj/<name>.md` per `reference/conjecture-doc.md`:
   formal statement (prose + the Lean name), evidence map, stress-test
   coverage, known special cases, and an initial attack-plan skeleton
   (candidate strategies, known obstructions). Gate on `conjecture-auditor`.
7. **Close the loop** — Logseq page for the conjecture linked from the
   line-of-inquiry note and topic pages; update `docs/conj/index.md`.

### 2. `attack-conjecture` (interactive companion)

Works a proof session against a `docs/conj/` page. A session rarely ends
with a full proof — the deliverable is *recorded progress*: lemmas closed,
obstructions hit, the plan updated. Nothing is lost when the session ends
mid-attack, and no `sorry` is ever left behind (an unproven lemma's statement
returns to the plan file, not the Lean tree).

Steps sketch:

1. **Load** — the conjecture doc, its Lean module, prior attack state.
   Working file `attack.md` in the scratchpad: the lemma DAG with per-lemma
   status (open / stated / proven / obstructed), rebuilt from the doc's
   attack-plan section (the scratchpad is session-scoped; the doc is the
   durable copy). If the load finds the main theorem already proven but
   unpromoted — a prior session closed it at the end of its budget — jump
   straight to step 5, so promotion's writing and audit run on a fresh
   context instead of the tail of an exhausted one.
2. **Plan** — decompose toward the goal with the user: candidate lemma
   chains, one `mathlib-scout` dispatch per planning round covering *all*
   candidate lemmas (its contract is per-claim dispositions, so a batch is
   its native input — per-lemma dispatches just multiply round-trips), risk
   ordering (attack the load-bearing uncertain lemma first — if it falls,
   the plan falls cheaply). Chalkboard sketches live here as the plan's
   prose; the README's "two angles" are two views of one DAG, not two plans.
3. **Prove, Lean-first per lemma** — state with `sorry`, `statement-auditor`
   gate on new statements, batched per planning round rather than one
   dispatch per lemma (here `attack.md` plays the coverage-map role: each
   lemma's entry carries its informal statement, its role in the DAG, and the
   Mathlib vocabulary its scout dispatch found — a multi-claim coverage map
   is already the shape that auditor expects), `plausible` before hard
   effort, then prove.
   The `/add-lean-topic` Step-6 triage applies: a stuck proof is a wrong
   statement, a hard proof, or a false claim — diagnose before pushing.
   A refutation of the *conjecture itself* mid-attack ends the session into
   `add-conjecture`'s refinement loop.
4. **Record** — write each lemma's status change through to the doc's
   attack-plan section *as it happens*, not in one sync at session end:
   `attack.md` dies with the session, and an attack session is exactly the
   kind that ends abruptly (context exhausted mid-proof). This step is then
   a closing pass — obstructions with their why, what the next session
   should try first — plus `just lean-check` on whatever landed.
5. **Promote on closure** — when the main theorem proves the Prop def:
   axiom audit; write `docs/theorems/<name>.md` per
   `reference/theorem-doc.md` — the human proof written *from* the Lean
   proof, at prose altitude, naming the Lean declaration; gate on
   `proof-auditor`; move the page out of `docs/conj/`; update both indexes
   and the Logseq conjecture page (status → proven). Promotion is a step
   here, not a third skill — it happens exactly once, at the moment the
   context is already loaded.

## Agents

### `counterexample-hunter` (isolator)

Used by `add-conjecture`, and by `attack-conjecture` for lemma-level checks.

The `data-profiler` analog for conjecture stress tests. Given a formal
statement, the `src/research` vocabulary it is expressed in, and the ranges
already covered by experiments; writes and runs search scripts in the
scratchpad (via `uv run` from the repo root) and returns either:

- **refuted**, with the *minimal* counterexample found (shrunk, reproducible:
  the object, the parameters, a one-liner to reconfirm), or
- **survived**, with an honest coverage report — ranges searched, search
  strategy (exhaustive / random / structured), case count, and the regions
  *not* reachable at this budget.

Raw sweep output never enters the conversation; "survived" without a coverage
report is treated as no answer. Search scripts go in a stable `hunt/`
subdirectory of the scratchpad, and each dispatch is told what's already
there: the refute→refine loop re-dispatches the hunter against near-identical
statements, and a hunter that rebuilds its harness from scratch each round
spends its budget on plumbing instead of search. Tools: `Read`, `Grep`,
`Glob`, `Bash`.

### `conjecture-auditor` (gates `add-conjecture`)

The `draft-auditor` analog for conjecture docs. Given the doc, the working
file, the Lean module, and the paths of cited results docs; checks:

- **[stmt]** the prose statement and the Lean `Prop` def say the same thing —
  quantifiers, side conditions, and vocabulary (this is *in addition to* the
  `statement-auditor` pass: that audits Lean against the claim source; this
  audits the published doc against the Lean).
- **[evidence]** every cited experiment, range, and finding matches what the
  results docs actually say; no evidence invented, none quietly omitted.
- **[scope]** tested-range vs. claimed-scope drift — the doc must state
  where the evidence ends; "verified" language never extends past the
  stress-test coverage report.
- **[known]** the statement is not already proven or refuted per the
  literature-checker and mathlib-scout dispatches (their outputs are inputs
  to the audit) — and those dispatches answered the statement *as published*,
  not a pre-refinement ancestor of it.
- **[trail]** the stress test is reported with its coverage, and any
  refinement forced by a counterexample appears in the revision log.

Tools: `Read`, `Grep`, `Bash` (Bash only to confirm the Lean name
elaborates — no proof work). The dispatch outputs it audits against —
literature dispositions, scout results, the hunter's coverage report — are
read from `conjecture.md`, never relayed from the conversation.

### `proof-auditor` (gates promotion in `attack-conjecture`)

Fresh-eyes audit of the human-facing theorem page against the machine-checked
proof. The failure mode is the prose telling a nicer story than what was
proven. Given `docs/theorems/<name>.md` and the Lean module; checks:

- **[stmt]** the page's statement is the Lean declaration's statement —
  hypotheses included, none dropped for elegance.
- **[gap]** every prose proof step corresponds to something the Lean proof
  actually does; a prose step with no Lean counterpart is flagged (it may be
  fine bridging text, but the auditor names it and the author decides).
- **[dep]** the stated dependencies match what the Lean proof uses.
- **[axiom]** the axiom audit is clean (`propext`, `Classical.choice`,
  `Quot.sound` only) and the page's Lean declaration name resolves.

Tools: `Read`, `Grep`, `Bash`. The `[axiom]` check re-runs `#print axioms`
itself via `lake env` rather than trusting the promoting session's report —
cheap, because step 4's `just lean-check` already built everything.

### Reused as-is

`mathlib-scout` (proof-ingredient surveys), `statement-auditor` (every new
Lean statement, conjecture or lemma), `literature-checker` (the known-result
check, and obstruction hunting during attack planning).

## Reference files

- `add-conjecture/reference/conjecture-doc.md` — the doc template matching
  `docs/conj/index.md`'s contract, plus the revision-log convention and the
  **Refuted** banner (refuted pages stay in `docs/conj/`, carrying the
  counterexample and a pointer to any successor statement).
- `add-conjecture/reference/statement-craft.md` — failure modes of statement
  formulation: hidden quantifier order, boundary/degenerate cases (n = 0,
  empty structures — the prose analog of Lean junk values), the generality
  dial, vocabulary mismatch between prose and `src/research` semantics.
- `attack-conjecture/reference/attack-plan.md` — the lemma-DAG format,
  status vocabulary, risk-first ordering, obstruction records, the
  per-lemma effort cap (two serious tactic strategies, then record as
  obstructed and move on), and the lemma-placement test (topic module if it
  would survive the conjecture being forgotten, conjecture module
  otherwise).
- `attack-conjecture/reference/theorem-doc.md` — the promotion template
  matching `docs/theorems/index.md`'s contract.

## Build order

1. `add-conjecture` + `counterexample-hunter` + `conjecture-auditor` — the
   phase's entry point, and valuable alone: a precisely stated, stress-tested,
   Lean-committed conjecture is a publishable artifact even before any attack.
2. `attack-conjecture` + `proof-auditor` — needs conjecture docs to exist
   first; the promotion path (and `theorem-doc.md`) rides along.

## Open questions

- **Refuted conjectures: keep or remove?** A refutation is a finding — the
  counterexample often matters more than the conjecture did. Leaning: keep
  the page in `docs/conj/` with a **Refuted** banner, the counterexample, and
  a pointer to whatever refined statement replaced it, rather than deleting.
  (The index's "moves to Theorems" covers the proven path only.) -
  DECISION: Keep refuted conjectures, per the lean.

- **Statement revision policy.** When the stress test or attack refutes a
  statement, is the refinement the *same* conjecture (revision log in place)
  or a *new* one (new page, old one marked refuted → superseded)? Leaning:
  in-place with a revision log while the conjecture is young (pre-attack),
  new page once an attack plan has accumulated — lemmas proven against the
  old statement are real dependencies. - DECISION: Agreed with the lean:
  in-place revision pre-attack, new page once an attack has accumulated.

- **Where does the Prop def go on promotion?** Options: stays in
  `Theorems/Conjectures/<Name>.lean` forever with the theorem elsewhere
  importing it, or the whole module migrates into the topic module. Leaning:
  keep `Conjectures/` as the permanent home of *statements-as-defs* and let
  the proving theorem live with its topic — the def's location then encodes
  its history, and nothing needs renaming at promotion time. - DECISION:
  Agreed — `Conjectures/` is the permanent home of statements-as-defs.

- **Lemma placement during attack.** Proven lemmas that are genuinely about a
  topic belong in that topic's module (they are `/add-lean-topic` "extend"
  material); lemmas that only serve this proof could live in the conjecture's
  module. Leaning: topic module when it would survive the conjecture being
  forgotten, conjecture module otherwise — same test Mathlib uses for
  private lemmas. - DECISION: Agreed with the survive-being-forgotten test.

- **Does `attack-conjecture` need a budget convention?** Proof search can eat
  a session. The plan's risk-first ordering helps, but an explicit per-lemma
  effort cap ("two serious tactic strategies, then record as obstructed and
  move on") may be worth encoding in `attack-plan.md` so sessions end with
  breadth rather than one heroic stuck proof. - DECISION: Yes — encode the
  per-lemma effort cap in `attack-plan.md`.
