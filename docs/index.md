# Research Harness

A computational research sandbox for the exploration, validation, and
publication of research topics. Work moves through the repository in a fixed
order — conjecture, computational evidence, machine-checked proof — and each
stage has a home here.

## Sections

<div class="grid cards" markdown>

- **[Conjectures](conj/index.md)** — proposed statements, not yet proven.
- **[Theorems](theorems/index.md)** — proven results, with their proofs.
- **[Results](results/index.md)** — interpreted output of experiment runs.
- **[Literature](ref/index.md)** — source papers held in the repository.
- **[Python API](reference/python/index.md)** — the `research` library and the
  `experiments` harness.
- **[Lean API](reference/lean/index.md)** — the `Theorems` Lean library.

</div>

## How the reference is built

Neither API reference is written by hand, and neither is checked into
`docs/`; both are generated on every build.

The site itself is built with [ProperDocs](https://properdocs.org/), the
maintained continuation of MkDocs 1.x.

The **Python** reference renders one page per module under `src/`, taking
signatures, type annotations, and docstrings directly from the sources via
[mkdocstrings](https://mkdocstrings.github.io/). Docstrings are Google-style,
and the Sphinx roles they use for cross-references — `` :func:`name` `` and
friends — are rewritten into working links.

The **Lean** reference is read out of the source text of `src/theorems`:
module docstrings, declaration docstrings, and the signature each one
introduces. It deliberately does not run `lake`, so it costs nothing to build
and needs no Mathlib cache — but it also reports signatures as written rather
than as Lean elaborates them.

## Building the docs

```console
$ just docs-serve     # live-reloading preview at http://127.0.0.1:8000
$ just docs           # build the static site into `site/`
$ just docs-check     # build with `--strict`, failing on any warning
```

The build resolves cross-references into the Python, pandas, matplotlib, and
Click documentation, so it needs network access on first run.
