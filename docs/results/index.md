# Results

Interpreted output of experiment runs — what a run showed and what it means.
The raw artifacts stay in `data/results/`, written by the experiment itself; the
pages here explain them.

A page in this section should carry:

- which experiment produced it, and the exact invocation, so the run can be
  reproduced;
- the parameters covered, and the range over which the conclusion is claimed;
- what the data shows, including where it is inconclusive;
- which [conjecture](../conj/index.md) it supports or refutes.

Experiments are written against the harness in `src/experiments/` and the
library in `src/research/`, both documented under the
[Python API](../reference/python/index.md). Run one with:

```console
$ just experiment <name> [args ...]
```
