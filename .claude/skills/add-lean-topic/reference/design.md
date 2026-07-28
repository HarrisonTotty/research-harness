# Design contract for topic formalizations

## The reuse ladder

For every claim on the page, take the highest rung that applies:

1. **Reuse.** Mathlib states it (possibly in a more general form — record
   what specializes it). The coverage map cites the declaration; nothing is
   built. Never restate or redefine a covered item "for readability": a
   duplicate connects to none of Mathlib's downstream API and forks the
   library.
2. **Extend.** Mathlib has the object but not this definition or theorem.
   Build it here as a declaration *about Mathlib's object*, in the object's
   namespace, following Mathlib conventions — these are the upstream
   candidates.
3. **Define.** Mathlib lacks the object itself. Define it fresh, still in
   Mathlib's conventions (the house rules already require this so material
   can migrate upstream without a rewrite), in the namespace of its
   principal subject.
4. **Backlog.** No tractable plan (the proof is a research project, or the
   definition needs machinery that is itself missing). Recorded in the
   coverage map and reported in Step 8 — never stated with a committed
   `sorry`, never stubbed.

## Working with Mathlib's axiomatization

Mathlib chose one formulation as the stored form; the page may lead with
another. Always work with Mathlib's: the page's other axiomatizations become
constructors and characterization lemmas *about* Mathlib's object (many
already exist — the scout reports them), not a parallel development. A
missing cryptomorphism from the page's table is an **extend** item: state it
as a theorem or constructor against Mathlib's definition.

## The coverage map

`coverage.md` must contain one row per page block across the Definition,
derived-vocabulary, operations, structural-theorems, and canonical-examples
sections:

- **Claim** — the page block, quoted verbatim. The statement auditor and
  Step 4 both work from this text, not from memory.
- **Disposition** — the scout's verdict (covered / partial / missing) with
  the locally verified declaration names and module paths.
- **Decision** — reuse / extend / define / backlog, with a one-line reason.
- **Plan** — for extend/define rows: the intended declaration name, a
  statement sketch, and the expected proof approach.

Open questions from the page get no rows — they are candidate future work to
report, not to formalize.

## Module shape

- One module per topic: `src/theorems/Theorems/<Topic>.lean`, imported from
  the root `src/theorems/Theorems.lean` so `lake build` covers it.
- Import the specific Mathlib modules the coverage map names, plus
  `Mathlib.Tactic` for proofs; do not import all of Mathlib.
- Declarations live in the namespace of their principal subject — for
  extensions of a Mathlib object, that is the object's own Mathlib
  namespace — so dot notation works and upstreaming is a copy, not a
  rewrite.
- The module docstring's implementation notes must say which claims are
  covered by Mathlib (and where) versus formalized here — the file should
  make sense to a reader who has not seen the coverage map.
- Every declaration's docstring carries the attribution the page records
  (name, year, source) for the claim it formalizes, per the house rules on
  citing research sources.
