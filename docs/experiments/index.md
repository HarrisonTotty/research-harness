# Experiments

Design and pre-registration records for the parameterized experiments in
`src/experiments/`. A page here is written *before* its experiment runs: it
fixes what the experiment varies, what the literature predicts for each region
of the parameter space, and what outcome would falsify the hypothesis — so that
"unexpected" is well-defined when the [results](../results/index.md) are
explored, and interpretation cannot quietly become post-hoc.

A page in this section should carry:

- the hypothesis, and the specific outcome that would falsify it;
- the parameter space — axes, ranges, granularity — with its controls,
  replication count, and seed policy;
- the output schema the planned analysis depends on;
- pre-registered predictions per parameter region, cited to the literature,
  with regions where no prediction exists explicitly marked.

Implementations live in `src/experiments/`, documented under the
[Python API](../reference/python/index.md); each design doc is named after its
experiment's CLI name (the module name, with hyphens for underscores). Run one
with:

```console
$ just experiment <name> [args ...]
```

Each run writes `data/results/{experiment}.{timestamp}.json` and a
`.meta.json` sidecar recording its parameters; the interpreted write-ups land
in [Results](../results/index.md).
