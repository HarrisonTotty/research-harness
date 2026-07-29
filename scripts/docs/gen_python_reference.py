"""Generate one mkdocstrings page per Python module under ``src``.

Run by ``mkdocs-gen-files`` during the build. Pages are virtual — nothing is
written into ``docs/`` on disk — so the reference cannot drift from the sources
and never needs regenerating by hand.
"""

from pathlib import Path

import mkdocs_gen_files
from mkdocs_gen_files.nav import Nav

_ROOT = Path(__file__).resolve().parents[2]
_SOURCE_ROOT = _ROOT / "src"
_OUTPUT_ROOT = Path("reference", "python")
_PACKAGES = ("research", "experiments", "figures")
"""Import packages to document, in nav order."""

_KEPT_DUNDERS = frozenset({"__init__", "__main__"})
"""Underscore-prefixed modules that are still part of the public surface."""

# Upstream leaves `Nav.__init__` unannotated, so a strict call is rejected.
navigation = Nav()  # type: ignore[no-untyped-call]

for package in _PACKAGES:
    for source in sorted((_SOURCE_ROOT / package).rglob("*.py")):
        relative = source.relative_to(_SOURCE_ROOT)
        parts = list(relative.with_suffix("").parts)
        page = relative.with_suffix(".md")

        if parts[-1] == "__init__":
            parts.pop()
            page = page.with_name("index.md")
        elif parts[-1].startswith("_") and parts[-1] not in _KEPT_DUNDERS:
            continue

        if not parts:
            continue

        identifier = ".".join(parts)
        navigation[tuple(parts)] = page.as_posix()

        with mkdocs_gen_files.open(_OUTPUT_ROOT / page, "w") as handle:
            handle.write(f"::: {identifier}\n")

        mkdocs_gen_files.set_edit_path(
            _OUTPUT_ROOT / page, source.relative_to(_ROOT).as_posix()
        )

_OVERVIEW = """# Python API

Generated from the sources under `src/` on every documentation build, so this
never falls out of step with the code. Signatures, type annotations, and
docstrings are read directly from the modules; the Sphinx roles the docstrings
use for cross-references are rewritten into working links.

{packages}
"""

with mkdocs_gen_files.open(_OUTPUT_ROOT / "index.md", "w") as handle:
    handle.write(
        _OVERVIEW.format(
            packages="\n".join(
                f"- [`{package}`]({package}/index.md)" for package in _PACKAGES
            )
        )
    )

with mkdocs_gen_files.open(_OUTPUT_ROOT / "SUMMARY.md", "w") as handle:
    # `index.md` first so `mkdocs-section-index` makes it the section's landing
    # page rather than a sibling entry.
    handle.write("* [Python API](index.md)\n")
    handle.writelines(navigation.build_literate_nav())
