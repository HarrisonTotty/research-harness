# Attack-plan contract

The format of the conjecture doc's **Attack plan** section and of the
session working file `attack.md` — one DAG, two copies: the doc is
durable, the scratchpad copy dies with the session. Every status change
is written through to the doc *as it happens*, never batched to session
end; attack sessions are exactly the kind that end abruptly.

## The lemma DAG

A short strategy paragraph first — the chalkboard view: the shape of the
argument, in prose. The README's "two angles" (chalkboard and Lean) are
two views of this one DAG, not two plans. Then one block per lemma:

```markdown
- **L3** <informal statement, one sentence> — `proven`
  - needs: L1, L2 (or: feeds the goal directly)
  - vocabulary: `<Mathlib decls from its scout dispatch>`
  - placement: topic (`Theorems/<Topic>.lean`) | conjecture module
  - lean: `<decl name>` once stated; the full statement text while
    unproven, so no `sorry` needs to survive in the tree
  - notes: <obstruction record, proof idea, session breadcrumbs>
```

## Status vocabulary

- `open` — an idea in the DAG; no Lean statement yet.
- `stated` — Lean statement written and audited. Transient *within* a
  session: at session end a stated-but-unproven lemma's statement text
  returns to its DAG entry and the `sorry` is deleted.
- `proven` — closed in Lean and placed; the entry names the
  declaration.
- `obstructed` — the effort cap was hit or a real obstruction found;
  the notes say what was tried, why it failed, and what to try next.
- `refuted` — falsified by `plausible` or the hunter; kept in the DAG
  as a record, since a refuted lemma reshapes the plan around it.

## Risk-first ordering

Attack the load-bearing uncertain lemma first: the lemma whose failure
invalidates the most downstream plan. If it falls, the plan falls
cheaply — before sessions were spent proving lemmas that fed it.
Easy-but-peripheral lemmas are budget filler, not progress.

## The per-lemma effort cap

**Two serious tactic strategies, then record as obstructed and move
on.** A "serious strategy" is a distinct proof approach (a different
induction, a different intermediate object, a different automation
family) — not two variations of one stuck tactic block. The cap exists
so sessions end with breadth: three lemmas obstructed-with-notes beat
one heroic stuck proof, because the notes are what the next session
(or the user at the chalkboard) needs. Before either strategy, the
stuck-proof triage applies — wrong statement, hard proof, or false
claim — and an obstruction record that skips the triage is incomplete.

## The lemma-placement test

Would the lemma survive the conjecture being forgotten? If a topic
reader would want it regardless — it is about the topic's objects, in
their vocabulary — it belongs in the topic module
(`src/theorems/Theorems/<Topic>.lean`); it is `/add-lean-topic`
"extend" material that happens to be discovered during an attack. If it
only serves this proof — stated in terms of the conjecture's specific
setup — it lives in the conjecture's module. Same test Mathlib applies
to private lemmas. Record the placement in the entry when the lemma is
stated, not when it is filed.

## Obstruction records

An obstruction record carries: the strategies tried (both of them),
where each one broke (the goal state or missing lemma, in a phrase),
the triage verdict, and the recommended next move. "Couldn't prove it"
is not a record; "induction on circuits stalls at the contraction step
— needs a rank identity Mathlib lacks (nearest: `Matroid.eRk_contract`);
try the closure characterization next" is.
