# Design-doc template for `docs/experiments/`

One design doc per experiment, named `docs/experiments/<name>.md` with
`<name>` matching the experiment's CLI name (hyphenated; the module in
`src/experiments/` uses underscores). The doc records *intent and
predictions* — everything a run's `.meta.json` cannot capture. It is
written before the experiment runs and read again by `/explore-results`,
which compares its predictions against the data.

Section skeleton:

```markdown
# <Experiment title>

One-paragraph summary: what is varied, what is measured, and why now.

## Hypothesis

The single falsifiable statement under test.

## Provenance

Where the hypothesis came from: Logseq page links, prior results docs
(`../results/<name>.<ts>.md`), papers. Follow-ups name the results doc
that suggested them.

## Falsification criterion

The specific observable outcome that would disprove the hypothesis.

## Parameter space

| Axis | Range | Granularity | Rationale |
| ---- | ----- | ----------- | --------- |

Rationale says why this range and step size can see the predicted effect.

## Controls and baselines

The reference arm(s) the effect is measured against, and what each
controls for.

## Replication and seeds

Replication count per grid cell; whether arms share seeds (paired) or
draw independently, and why; how seeds are derived and where they are
recorded.

## Output schema

One row per <unit>.

| Column | Dtype | Meaning / units |
| ------ | ----- | --------------- |

## Pre-registered predictions

| Region | Prediction | Source |
| ------ | ---------- | ------ |

One row per region of the parameter space. Source cites the Logseq page
or paper making the prediction. Regions with no prediction get an
explicit "no prediction" row — they are the interesting ones. **Frozen
once the first run exists**; later insight goes in a dated addendum
subsection, never as an edit to these rows.

## Analysis plan

The analyses the exploration step will run — each must be supported by
the schema above.

## Implementation notes

Written by /add-experiment, not during design: module path, pilot-run
observations (runtime per cell, observed variance), and the suggested
full-sweep invocation.
```

Rules:

- Every prediction carries a source; a prediction nobody can cite is a
  guess and belongs in a "no prediction" row's notes instead.
- Ranges and granularity are numbers, not adjectives.
- Links are relative site links (`../results/...`, `../conj/...`) so the
  page builds cleanly into the docs site.
