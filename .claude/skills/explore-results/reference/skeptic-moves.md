# Skeptic moves

The shared skeptical repertoire: `/explore-results` applies it to every
claim before drafting (Steps 3–5), and `/critique-results` applies it
again, conversationally, to every claim the user walks through.

## Alternative explanations

Before a claim reads "the effect is real", run the alternatives:

- **Confound** — does anything co-vary with the swept axis that could
  produce the pattern on its own?
- **Grid artifact** — would a finer or shifted grid plausibly dissolve
  it? A trend supported by three points is a hypothesis, not a trend.
- **Seed artifact** — does it survive across replications, and is the
  replication count large enough that it should?
- **Numerics** — dtype coercion, accumulation error, a tolerance doing
  silent work in a comparison.
- **Selection** — was this region singled out after seeing the data?
  Check against the pre-registered regions.
- **Implementation** — does the `src/research` call mean what the
  analysis assumes it means (convention, normalization, defaults)?

Name the alternatives considered in the write-up; an unconsidered
alternative is a `[drift]` audit finding waiting to happen.

## Bug-first triage

A claim that contradicts recorded literature (a `literature-checker`
`[contradicts]` disposition) is a bug until proven otherwise. In order:

1. **Scope the contradiction.** Quote the contradicted claim and check
   the scopes actually clash: finite-`n` observation vs. asymptotic law,
   differing conventions or normalizations — scope mismatches resolve
   many apparent contradictions outright.
2. **Hunt the specific bug class** that would produce exactly this
   contradiction: the sweep's confounds, seed leakage, a model call with
   the wrong semantics, a schema column meaning something else.
3. **Reproduce minimally** — if cheap, re-run the smallest slice that
   exhibits the contradiction (into the scratchpad, never
   `data/results/`) and hand-check one instance.

Only after all three may the claim appear as a finding — and the doc
must carry the trail: what was contradicted, what was checked, what was
ruled out.

## What would falsify this

Every claim in a results doc names the observation that would kill it.
Follow-up suggestions are then chosen for falsifying power: prefer the
experiment that could produce the killing observation over the one that
can only add confirmation. The sharpest follow-ups test *implications*
of the current interpretation in a regime not yet swept — if the
interpretation is right, it constrains what must happen there.

## Claim-strength ladder

- **consistent with** — the data does not contradict it; nothing more.
- **suggests** — a visible trend; alternatives not yet excluded.
- **shows** — reserved for effects that span the swept range, survive
  replication, and have had the alternative explanations above
  addressed.
- **proves** — never available from a finite sweep.

Every claim is scoped to the swept ranges; anything beyond them is
labeled extrapolation.
