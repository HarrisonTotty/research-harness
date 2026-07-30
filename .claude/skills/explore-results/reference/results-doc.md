# Results-doc template for `docs/results/`

One doc per explored run, named `docs/results/<name>.<ts>.md`. The doc
must let a reader reproduce the run, see what the data showed (including
where it is inconclusive), and know which conjecture it bears on — the
contract stated in `docs/results/index.md`. Numbers are copied from
analysis outputs, never recalled.

Section skeleton:

```markdown
# <Experiment title>: <ts> run

Two-to-three sentence overview — the headline outcomes, phrased on the
claim-strength ladder.

## Run

- Invocation: the exact reproducing command
  (`just experiment <name> [args ...]`).
- Design doc: link to `../experiments/<name>.md` (or an explicit note
  that the run was not pre-registered).
- Commit and parameters: from the `.meta.json`, as a table — the range
  over which every conclusion below is claimed.

## Data integrity

The profiler's verdict. `SOUND`, or the findings and why each is benign
(explained by the design). Caveats stated here must survive into the
findings they touch.

## Prediction outcomes

| Region | Predicted (source) | Observed | Outcome |
| ------ | ------------------ | -------- | ------- |

One row per pre-registered prediction; Outcome is confirmed /
contradicted / no prediction existed. For a non-pre-registered run,
replace the table with an explicit note that the run was not
pre-registered.

## Findings

Per finding: the claim (scoped to the swept range), the supporting
numbers, and the alternatives considered (the skeptic moves applied). A
finding built on a contradicted prediction carries its literature
disposition and triage trail inline.

## Unexpected observations

Observations outside every prediction, each with its literature-check
disposition ([consistent] / [contradicts] + triage / [no-coverage]).

## Suspected bugs

Anything triaged as a probable experiment defect, with the evidence.
State "None." explicitly rather than omitting the section.

## Comparison with <baseline ts>          (compare mode only)

Parameter diff of the two runs (changed / added / removed; identical
parameters summarized in one line), then outcome deltas on the shared
measures over the shared grid.

## Follow-up experiments

Falsification-oriented: for each suggestion, the implication it tests
and the observation that would kill the current interpretation.

## Conjecture links

Which `../conj/` pages this run supports or refutes — or an explicit
"none yet", making the run a candidate seed for one.
```

Rules:

- Claim strength on the ladder in
  [skeptic-moves.md](skeptic-moves.md); "proves" never appears.
- Every conclusion is scoped to the parameters actually swept;
  extrapolation is labeled as such.
- Links are relative site links (`../experiments/...`, `../conj/...`) so
  the page builds cleanly into the docs site.
