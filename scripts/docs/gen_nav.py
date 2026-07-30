"""Build the site's top-level navigation from whatever ``docs/`` contains.

The prose sections of this repository — conjectures, theorems, interpreted
results, reference literature — grow a file at a time, and an explicit ``nav``
in ``properdocs.yml`` would need editing on every addition. This emits the root
``SUMMARY.md`` for ``mkdocs-literate-nav`` instead, listing only the sections
that currently hold something. Each section's own contents are then discovered
by ``literate-nav``, from a nested ``SUMMARY.md`` where one is generated and by
directory listing otherwise.

The reference-literature section is a directory of PDFs rather than pages, so it
also gets a generated index linking to them.
"""

from pathlib import Path

import mkdocs_gen_files

_ROOT = Path(__file__).resolve().parents[2]
_DOCS_ROOT = _ROOT / "docs"

_PROSE_SECTIONS = (
    ("conj", "Conjectures"),
    ("theorems", "Theorems"),
    ("experiments", "Experiments"),
    ("results", "Results"),
)
"""Hand-written Markdown sections, in nav order."""

_LITERATURE = "ref"
"""Directory of reference papers, indexed rather than listed as pages."""

_LITERATURE_TITLE = "Literature"

_API_SECTIONS = (
    ("reference/python", "Python API"),
    ("reference/lean", "Lean API"),
)
"""Generated API references, in nav order."""


def _has_pages(directory: Path) -> bool:
    """Report whether ``directory`` holds at least one Markdown page."""
    return directory.is_dir() and any(directory.rglob("*.md"))


def _write_literature_index(directory: Path) -> bool:
    """Generate an index of the PDFs in ``directory``; report whether it has any.

    A hand-written ``index.md`` wins: the generated listing is a fallback for
    when nobody has written a proper bibliography yet.
    """
    if not directory.is_dir():
        return False
    papers = sorted(directory.glob("*.pdf"))
    if not papers:
        return False
    if (directory / "index.md").exists():
        return True

    lines = [
        f"# {_LITERATURE_TITLE}",
        "",
        "Source papers held in this repository, served alongside the site.",
        "",
    ]
    lines += [f"- [{paper.stem}]({paper.name})" for paper in papers]

    with mkdocs_gen_files.open(Path(_LITERATURE, "index.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return True


entries = ["* [Home](index.md)"]

entries += [
    f"* [{title}]({directory}/)"
    for directory, title in _PROSE_SECTIONS
    if _has_pages(_DOCS_ROOT / directory)
]

if _write_literature_index(_DOCS_ROOT / _LITERATURE):
    entries.append(f"* [{_LITERATURE_TITLE}]({_LITERATURE}/)")

entries += [f"* [{title}]({directory}/)" for directory, title in _API_SECTIONS]

with mkdocs_gen_files.open("SUMMARY.md", "w") as handle:
    handle.write("\n".join(entries) + "\n")
