# Theorem-page template for `docs/theorems/`

The promotion target: one page per proven conjecture, named
`docs/theorems/<name>.md` (keeping the conjecture page's name). Written
only after the axiom audit is clean, and written *from* the closed Lean
proof — the page carries what `docs/theorems/index.md` promises: the
statement, a human proof, the Lean declaration name, and the
dependency context.

Section skeleton:

```markdown
# <Theorem title>

One-paragraph summary: what is now proven, and the line of inquiry it
closes — linking the evidence that motivated it
(`../results/<name>.<ts>.md`).

## Statement

The statement in prose, in the same form as the Lean declaration —
every hypothesis included, none dropped for elegance.

**Lean**: `<Namespace.theorem_name>`
(`src/theorems/Theorems/<Topic>.lean`), proving
`<Namespace.DefName>`
(`src/theorems/Theorems/Conjectures/<Name>.lean`).

## Proof

Written for a human reader, at prose altitude — the argument, not a
transcription of the tactic script. Every prose step corresponds to
something the Lean proof actually does; bridging text that has no Lean
counterpart is fine, but it reads as exposition, never as a proof step.
Named lemmas from the attack are cited by their Lean declarations.

## Dependencies

What the proof uses — the Mathlib results and repo lemmas it actually
invokes — and what this theorem is used to prove, once anything is.

## History

The conjecture's trail: when it was registered, the revisions its
stress tests forced, and the obstructions the attack recorded on the
way. The page replaces `docs/conj/<name>.md`, so the trail worth
keeping is summarized here.
```

Rules:

- **The Lean proof is the source.** The prose proof is written from the
  closed, axiom-audited Lean proof — never in parallel from scratch.
  Where the Lean proof is ugly, the prose may reorganize, but only over
  steps the Lean actually takes; the `proof-auditor` names every prose
  step with no Lean counterpart.
- The statement section restates the *proven* form. If the attack
  narrowed or reshaped the statement from the conjecture's original,
  the History section says so — the page never quietly presents the
  nicer unproven version.
- The `Prop` def stays in `Conjectures/` forever; the theorem cites it
  rather than restating it, so the audited proposition and the proven
  one are the same term.
- Links are relative site links (`../results/...`, `../conj/...`) so
  the page builds cleanly into the docs site.
