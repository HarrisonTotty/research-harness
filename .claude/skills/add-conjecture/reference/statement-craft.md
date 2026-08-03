# Failure modes of conjecture formulation

Raise these conversationally as the statement takes shape — each is a
question to settle with the user, not a gate. They are ordered by how
expensive they are to discover late.

## Hidden quantifier order

"For every X there is a Y" and "there is a Y that works for every X"
read almost identically in casual prose and are different conjectures.
Informal notes almost never mark the order; the experiments usually
tested the weaker one (a Y was found per X — nothing says it was the
same Y). Make every quantifier explicit and confirm which order the
evidence actually supports.

## Boundary and degenerate cases

n = 0 and n = 1, the empty structure, the single-element structure, the
disconnected case — the prose analogs of Lean junk values. Decide
explicitly: does the claim hold there, is it excluded by a side
condition, or is it vacuous? A conjecture refuted at n = 0 by
convention rather than substance wastes a stress-test round and a
revision; a side condition added silently to dodge one weakens the
claim without anyone deciding to.

## Silent generator conditions

The experiments' data came from a generator, and generators impose
conditions nobody wrote down — connectedness, simplicity, a fixed
ground-set size, nondegenerate parameters. For each such condition,
decide: is it a *hypothesis* of the conjecture, or an accident of the
sampling? An accident left out of the statement means the evidence
supports a narrower claim than the page states; a hypothesis quietly
inherited means the conjecture is weaker than the interesting one.

## The generality dial

The *minimal* statement the evidence supports and the *natural* general
form are usually different statements. Neither is automatically right:
the minimal one is better evidenced, the general one is the one worth
proving. When they differ, state both — the general form as the
conjecture, the minimal one as its implied special case — so the
stress test can attack the general form where the evidence is thinnest.

## Claim-strength choices

Equality vs. inequality, exact vs. asymptotic, "for all n" vs. "for n
sufficiently large", monotone vs. eventually monotone. Experiments over
finite ranges cannot distinguish these; the choice is a judgment about
mechanism, so it must be the user's, made explicitly. Prefer the form
that is falsifiable at the ranges the hunter can reach.

## Vocabulary mismatch

The statement must mean the same thing in prose, in `src/research`, and
in Lean. The usual traps: normalization conventions (is the invariant
scaled?), indexing (0- vs. 1-based), strict vs. non-strict comparisons,
and near-miss concepts (rank vs. corank, circuit vs. cocircuit). The
vocabulary map — each prose term paired with the `src/research` symbol
and Mathlib declaration it denotes — is written before the Lean is, and
the scout dispositions are its source of truth.
