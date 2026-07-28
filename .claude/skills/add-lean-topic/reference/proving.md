# Statement and proof contract for topic formalizations

## Statements: transcribe, then defend against Lean

Transcribe each statement from the verbatim claim text in `coverage.md`, then
re-read it as mathematics with Lean's totality in mind — this is where the
statement auditor will look hardest:

- Lean's total functions return junk outside their mathematical domain
  (division by zero is zero, `Nat` subtraction truncates, suprema of empty
  sets are arbitrary). Add the domain hypotheses that keep the statement out
  of the junk regime; prefer reformulating away `Nat` subtraction entirely.
- Preserve claim strength exactly: an `↔` only where the page states an
  equivalence; the page's quantifiers and hypothesis order; no generality
  the page does not license (generalize only when the proof is unchanged and
  note it in the docstring).
- Be explicit about numeric types wherever literals, division, or coercions
  appear, and keep coercions at the boundary of the statement.
- Follow the house Lean rules for naming (conclusion first, `_of_` for
  hypotheses) and layout; they apply unchanged.

## The sorry lifecycle

`sorry` exists in exactly one window: after Step 4 states a declaration and
before Step 6 finishes it. Inside that window it is scaffolding that lets
statements elaborate and be audited before proof effort is spent. At the end
of Step 6 the count is zero — an unfinished proof means the statement is
deleted and its coverage-map row moved to backlog. A committed `sorry` is
never an acceptable intermediate state.

## Canonical examples

- Each canonical example becomes a `def` with the page's data, plus theorems
  proving exactly what the page says the example certifies — named for the
  certification, since that is the example's reason to exist.
- Counterexamples prove the failing direction explicitly: the property
  violated *and* the property retained, so the one-sidedness the page
  records is enforced.
- Small finite certificates may close with `decide` — it is kernel-checked.
  `native_decide` is not: it trades the kernel for trusted evaluation and
  must never close a result the repository presents as proved.

## When a proof gets stuck

A stuck proof is a three-way question; diagnose before pushing harder:

1. **The statement is subtly wrong.** Check it with `plausible` and against
   the canonical examples first — minutes of counterexample hunting are
   cheaper than hours of proving. A counterexample to your transcription
   means fixing the statement (and re-auditing, Step 5).
2. **The proof is genuinely hard.** Delete the statement, move the row to
   backlog, report it. Do not commit `sorry`; do not paper over a failing
   step with heavier automation without knowing why it fails.
3. **The page's claim is false.** A `plausible` counterexample that survives
   against the *faithful* transcription is a research finding about the
   graph: stop and report it (SKILL.md, Step 6); never quietly weaken the
   statement to something provable the page did not claim.

## Axiom audit

After the gate passes, print axiom dependencies for each top-level theorem —
e.g. `#print axioms Matroid.myTheorem` in a scratch file, or temporarily at
the end of the module — and confirm nothing beyond `propext`,
`Classical.choice`, and `Quot.sound`:

- `sorryAx` — unfinished work escaped Step 6; the gate's warnings should
  have caught it, but the axiom check is what makes sure.
- `Lean.ofReduceBool` / `Lean.trustCompiler` — native evaluation was trusted
  rather than kernel-checked somewhere in the proof; replace it.

Do not commit the `#print axioms` lines; they are a check, not content.
