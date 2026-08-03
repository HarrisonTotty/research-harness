# Conjecture-page template for `docs/conj/`

One page per conjecture, named `docs/conj/<name>.md` (hyphenated), with
the Lean `Prop` def in `src/theorems/Theorems/Conjectures/<Name>.lean`
(UpperCamelCase, matching the def). The page carries everything
`docs/conj/index.md` promises: the formal statement, the evidence, known
special cases, and what a proof would need.

Section skeleton:

```markdown
# <Conjecture title>

One-paragraph summary: the claim informally, the line of inquiry it came
from, and its current status.

## Statement

The formal statement in prose — every quantifier and side condition
explicit, in the same form as the Lean def.

**Lean**: `<Namespace.DefName>`
(`src/theorems/Theorems/Conjectures/<Name>.lean`)

Where the natural general form and the minimal evidenced statement
differ, the general form is the conjecture and the minimal form is
stated beneath it as the implied special case.

## Evidence

One entry per supporting experiment: the results doc
(`../results/<name>.<ts>.md`), the ranges *actually swept*, and what was
found. Quoted from the docs, never paraphrased from memory — the gap
between tested range and claimed scope is the point of a conjecture, and
it must stay visible.

## Stress-test coverage

The counterexample search that the statement survived: per range, the
strategy (exhaustive / structured / random) and case count, plus the
regions not reached at the search budget. "Verified" language never
extends past this section.

## Known special cases

Proven instances, each with its Lean declaration and the module it lives
in.

## Attack plan

Maintained by /attack-conjecture; seeded at creation with candidate
strategies and known obstructions. Per-lemma status is written through
as it changes — this section is the durable copy of the attack state.

## Revision log

One dated entry per statement revision: the counterexample (with its
self-contained reconfirmation command) and the change it forced.
```

Rules:

- **Attack-plan format.** The section follows the lemma-DAG contract in
  [attack-plan.md](../../attack-conjecture/reference/attack-plan.md) —
  linked here, not from the page: published pages never link into
  `.claude/`, and this template's skeleton must instantiate with
  site-relative links only.
- **Reconfirmation commands are self-contained.** The revision log and
  the Refuted banner quote the hunter's reconfirmation command, which
  imports from `src/research` only — never a path into a session
  scratchpad, which dies with the session.
- **Revision policy.** Pre-attack, a refuted statement is revised in
  place with a revision-log entry — the conjecture is still young.
  Once an attack plan has accumulated, a refutation gets a *new* page:
  lemmas proven against the old statement are real dependencies, so the
  old page is kept and marked refuted → superseded.
- **Refuted banner.** A refuted page stays in `docs/conj/` — the
  counterexample often matters more than the conjecture did. Its
  summary paragraph is replaced by a `!!! failure "Refuted"` admonition
  carrying the minimal counterexample, the reconfirmation one-liner,
  and a link to the successor statement (or "no successor" explicitly).
  The rest of the page is left as the historical record.
- **Promotion.** When the proof closes, the page's content moves to
  `docs/theorems/<name>.md` (per
  [theorem-doc.md](../../attack-conjecture/reference/theorem-doc.md))
  and this page is deleted; the index line moves with it. The Lean def
  never moves.
- Each page gets a line in `docs/conj/index.md` — name, one-line
  statement, status (open / refuted → successor) — added at creation.
- Links are relative site links (`../results/...`, `../theorems/...`)
  so the page builds cleanly into the docs site.
