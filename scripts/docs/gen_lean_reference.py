"""Generate one reference page per Lean module under ``src/theorems``.

Run by ``mkdocs-gen-files`` during the build. Pages are rendered from the source
text by :mod:`lean_parser`, so the reference stays in step with the library
without a ``lake`` build — see that module for what a purely syntactic reading
can and cannot report.
"""

import re
import sys
from pathlib import Path

import mkdocs_gen_files
from mkdocs_gen_files.nav import Nav

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lean_parser import LeanDeclaration, LeanProse, parse_lean_source

_ROOT = Path(__file__).resolve().parents[2]
_SOURCE_ROOT = _ROOT / "src" / "theorems"
_OUTPUT_ROOT = Path("reference", "lean")

_HEADING = re.compile(r"^(#{1,6})(\s)", re.MULTILINE)
_MAX_HEADING_LEVEL = 6

_MODULE_DOC_DEMOTION = 1
"""Module docstrings sit under the page's ``# <module name>`` heading."""

_DECLARATION_DOC_DEMOTION = 3
"""Declaration docstrings sit under their own ``### <name>`` heading."""


def demote_headings(markdown: str, levels: int) -> str:
    """Push every ATX heading in ``markdown`` down by ``levels``.

    Lean module docstrings are written as standalone documents starting at ``#``.
    Nesting them beneath the page heading keeps one ``<h1>`` per page and gives
    the table of contents a correct shape.
    """

    def deepen(match: re.Match[str]) -> str:
        level = min(len(match.group(1)) + levels, _MAX_HEADING_LEVEL)
        return "#" * level + match.group(2)

    return _HEADING.sub(deepen, markdown)


def module_name(source: Path) -> str:
    """Return the Lean module name for a file under ``src/theorems``."""
    return ".".join(source.relative_to(_SOURCE_ROOT).with_suffix("").parts)


def page_path(source: Path) -> Path:
    """Return the output page for ``source``, relative to the Lean reference root.

    A file whose name matches a sibling directory — ``Theorems.lean`` alongside
    ``Theorems/`` — is that directory's root module, and becomes its index page.
    """
    relative = source.relative_to(_SOURCE_ROOT)
    if source.with_suffix("").is_dir():
        return relative.with_suffix("") / "index.md"
    return relative.with_suffix(".md")


def render(source: Path) -> str:
    """Render the Markdown page for one Lean source file."""
    name = module_name(source)
    module = parse_lean_source(name, source.read_text(encoding="utf-8"))

    lines = [f"# {name}", ""]
    if module.docstring:
        lines += [demote_headings(module.docstring, _MODULE_DOC_DEMOTION), ""]

    for entry in module.entries:
        if isinstance(entry, LeanProse):
            lines += [demote_headings(entry.body, _MODULE_DOC_DEMOTION), ""]
            continue
        lines += _render_declaration(entry)

    return "\n".join(lines).rstrip() + "\n"


def _render_declaration(declaration: LeanDeclaration) -> list[str]:
    """Render one declaration as a heading, a signature block, and its prose."""
    lines = [
        f"### {declaration.name} {{ #{declaration.name} }}",
        "",
        "```lean4",
        declaration.signature,
        "```",
        "",
    ]
    if declaration.docstring:
        lines += [
            demote_headings(declaration.docstring, _DECLARATION_DOC_DEMOTION),
            "",
        ]
    return lines


_OVERVIEW = """# Lean API

Generated from the sources under `src/theorems/` on every documentation build.

Unlike the Python reference, this is read out of the *source text* rather than
from an elaborated library: module docstrings, declaration docstrings, and the
signature each one introduces. That keeps the docs buildable without a `lake`
build or a Mathlib cache, at the cost of showing signatures as they were written
rather than as Lean elaborates them. For fully resolved types, hover the
declaration in an editor or build the library with `just lean-build`.
"""

# Upstream leaves `Nav.__init__` unannotated, so a strict call is rejected.
navigation = Nav()  # type: ignore[no-untyped-call]

if _SOURCE_ROOT.is_dir():
    for lean_source in sorted(_SOURCE_ROOT.rglob("*.lean")):
        page = page_path(lean_source)
        navigation[tuple(module_name(lean_source).split("."))] = page.as_posix()

        with mkdocs_gen_files.open(_OUTPUT_ROOT / page, "w") as handle:
            handle.write(render(lean_source))

        mkdocs_gen_files.set_edit_path(
            _OUTPUT_ROOT / page, lean_source.relative_to(_ROOT).as_posix()
        )

with mkdocs_gen_files.open(_OUTPUT_ROOT / "index.md", "w") as handle:
    handle.write(_OVERVIEW)

with mkdocs_gen_files.open(_OUTPUT_ROOT / "SUMMARY.md", "w") as handle:
    # `index.md` first so `mkdocs-section-index` makes it the section's landing
    # page rather than a sibling entry.
    handle.write("* [Lean API](index.md)\n")
    handle.writelines(navigation.build_literate_nav())
