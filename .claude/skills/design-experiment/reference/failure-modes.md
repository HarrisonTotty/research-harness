# Experiment-design failure modes

The checklist behind Step 2's checkpoints. For each mode: what goes
wrong, and the question to raise with the user while the corresponding
design element is on the table. These are applied conversationally — the
point is to ask the question before the sweep is expensive, not to gate.

## Confounded sweep

Two things change together, so the effect cannot be attributed. Often
subtle: an axis that also changes problem size or density, a derived
quantity computed once and reused across arms, shared caches or
warm-started state. **Ask:** for each swept axis, what else changes when
it moves — and is that intended?

## Accidental seed policy

Same seeds reused across arms (arms correlated), or fresh seeds per arm
(arm confounded with draw) — either can be right, but only as a choice.
Paired seeds reduce variance for comparisons; independent seeds keep
arms exchangeable. **Ask:** should the arms see the same random worlds
or different ones, and where does the answer get recorded?

## Grid too coarse to see the effect

The predicted effect lives at a scale smaller than the step size, so
even a real effect produces a flat-looking sweep. **Ask:** if the effect
is exactly where theory says, how many grid points land inside it?

## Schema starves the analysis

The planned analysis needs per-instance values but the schema stores
aggregates, or needs a column nobody emitted. Cheapest bug to catch at
design time, most expensive after a full sweep. **Ask:** walk each item
of the analysis plan against the schema — which columns does it read?

## Missing baseline

An effect claim with nothing to measure against. **Ask:** compared to
what — and is that arm in the sweep?

## Ceiling and floor saturation

The measure saturates over part of the range, flattening real
differences into ties. **Ask:** near the range's ends, can the measure
still move in both directions?

## Combinatorial runtime blowup

Cells × replications × cost per cell, estimated *before* committing to
granularity — the pilot run refines the estimate but should not be the
first time anyone multiplies the numbers. **Ask:** what is the cell
count, and what does a cell cost?

## Ambiguous measurement

Units, normalization, or convention left implicit — which definition of
the quantity the code computes, normalized by what. Becomes a `[model]`
or `[schema]` audit finding later if unstated now. **Ask:** for each
measured column, which exact definition and units?

## Underpowered replication

Effects claimed later that are smaller than run-to-run noise. Variance
is usually unknown at design time — that is what the pilot estimates —
but the replication count should state what noise level it assumes.
**Ask:** how big is the predicted effect relative to expected noise, and
what replication does that need?

## Post-hoc region selection

Regions of interest chosen after seeing the data, then presented as if
targeted. Pre-registration is the defense: the interesting regions are
named now. **Ask:** which regions do we expect to matter — and are the
rest labeled "no prediction"?
