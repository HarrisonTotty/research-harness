"""Extract module and declaration documentation from Lean 4 source text.

Lean's own documentation generator (``doc-gen4``) elaborates the library *and*
every module it imports, which for a Mathlib-based project means an hours-long
build producing gigabytes of HTML. This module takes the cheap route instead: it
reads the source text and pulls out module docstrings (``/-! ... -/``),
declaration docstrings (``/-- ... -/``), and the signature of the declaration
each one precedes. That is enough to render a browsable reference from the
repository's Lean conventions without invoking ``lake``.

The parser is deliberately syntactic. It tracks ``namespace``/``end`` nesting to
qualify names, and it stops a signature where the definition's body begins, but
it never elaborates a term. Signatures are therefore reported exactly as they
were written rather than as Lean would normalize them, and elaborated
information — inferred instance arguments, resolved notation, types of
``example`` blocks — is simply absent.
"""

import re
import textwrap
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

__all__ = [
    "LeanDeclaration",
    "LeanEntry",
    "LeanModule",
    "LeanProse",
    "parse_lean_source",
]

_COMMENT_OPEN = "/-"
_COMMENT_CLOSE = "-/"
_MODULE_DOC_MARKER = "!"
_DECLARATION_DOC_MARKER = "-"

_DECLARATION_KEYWORDS = frozenset(
    {
        "abbrev",
        "axiom",
        "class",
        "def",
        "inductive",
        "instance",
        "lemma",
        "opaque",
        "structure",
        "theorem",
    }
)
"""Keywords introducing a named declaration worth documenting.

``example`` is intentionally absent: it is anonymous, so there is nothing to
link to or index.
"""

_MODIFIER_KEYWORDS = frozenset(
    {
        "noncomputable",
        "nonrec",
        "partial",
        "private",
        "protected",
        "scoped",
        "unsafe",
    }
)
"""Words that may precede a declaration keyword without changing its shape."""

_BODY_MARKER = ":="
"""Token separating a declaration's signature from its body."""

_WHERE_MARKER = re.compile(r"(?:^|\s)where\b")
"""``where`` opens a structure body or a set of equations, ending the signature."""

_NAME = re.compile(r"^[^\s(){}\[\]:=,⦃⦄⟨⟩]+$")
"""A bare declaration name, once punctuation has been split away."""

_NAME_SEPARATOR = re.compile(r"[\s(\[{:]")
"""Characters that may abut a declaration name without whitespace."""

_SIGNATURE_INTERRUPTERS = ("|", "@[", "--", "attribute ", "end ", "namespace ")
"""Line prefixes that cannot continue a signature already in progress."""


@dataclass(frozen=True, slots=True)
class LeanDeclaration:
    """A single documented declaration."""

    kind: str
    """The introducing keyword, e.g. ``def`` or ``theorem``."""

    name: str
    """The declaration's name, qualified by any enclosing namespaces."""

    signature: str
    """The declaration as written, up to but excluding its body."""

    docstring: str | None
    """The ``/-- ... -/`` docstring attached to it, if any."""


@dataclass(frozen=True, slots=True)
class LeanProse:
    """A ``/-! ... -/`` block appearing after the module docstring.

    Mathlib convention uses these to introduce sections partway down a file, so
    they are kept in source order alongside the declarations they precede.
    """

    body: str
    """The block's Markdown content, dedented."""


type LeanEntry = LeanDeclaration | LeanProse
"""One item of a module's body, in source order."""


@dataclass(frozen=True, slots=True)
class LeanModule:
    """Everything worth documenting about one Lean source file."""

    name: str
    """Fully-qualified module name, e.g. ``Theorems.Basic``."""

    docstring: str | None
    """The file's leading ``/-! ... -/`` module docstring, if any."""

    entries: tuple[LeanEntry, ...]
    """Documented declarations and section prose, in source order."""


@dataclass(frozen=True, slots=True)
class _Block:
    """A lexed chunk of source: either a comment block or a line of code."""

    kind: str
    """One of ``module_doc``, ``declaration_doc``, ``comment``, or ``code``."""

    text: str
    """Comment body, or the code line itself."""


def _scan_comment(source: str, start: int) -> tuple[str, int]:
    """Return the body of the comment opening at ``start``, and the index after it.

    Lean block comments nest, so this tracks depth rather than searching for the
    first ``-/``. An unterminated comment yields everything to end of file.
    """
    depth = 0
    index = start
    body_start = start
    while index < len(source):
        if source.startswith(_COMMENT_OPEN, index):
            depth += 1
            index += len(_COMMENT_OPEN)
            if depth == 1:
                # Skip the `!` or `-` that distinguishes a docstring from a
                # plain comment, so it does not leak into the body.
                if index < len(source) and source[index] in (
                    _MODULE_DOC_MARKER,
                    _DECLARATION_DOC_MARKER,
                ):
                    index += 1
                body_start = index
        elif source.startswith(_COMMENT_CLOSE, index):
            depth -= 1
            if depth == 0:
                return source[body_start:index], index + len(_COMMENT_CLOSE)
            index += len(_COMMENT_CLOSE)
        else:
            index += 1
    return source[body_start:], len(source)


def _lex(source: str) -> Iterator[_Block]:
    """Split ``source`` into comment blocks and the code lines between them."""
    pending: list[str] = []
    index = 0
    while index < len(source):
        if source.startswith(_COMMENT_OPEN, index):
            marker = source[index + len(_COMMENT_OPEN) : index + len(_COMMENT_OPEN) + 1]
            kind = {
                _MODULE_DOC_MARKER: "module_doc",
                _DECLARATION_DOC_MARKER: "declaration_doc",
            }.get(marker, "comment")
            for line in "".join(pending).splitlines():
                yield _Block("code", line)
            pending = []
            body, index = _scan_comment(source, index)
            yield _Block(kind, body)
        else:
            pending.append(source[index])
            index += 1
    for line in "".join(pending).splitlines():
        yield _Block("code", line)


def dedent_block(body: str) -> str:
    """Normalize a comment body into Markdown.

    Removing the common indent has to account for the first line, which may sit
    on the same line as the opening delimiter and so carry no indentation of its
    own. Including it in the common prefix would then dedent everything to zero
    and flatten any indented block — a nested list or code sample — inside the
    docstring, so in that case it is stripped on its own and the common prefix is
    taken over the remaining lines.
    """
    if not body.strip():
        return ""

    if body.lstrip(" \t").startswith("\n"):
        # Nothing shares the delimiter's line, so every line carries its own
        # indentation and the common prefix is already meaningful.
        return textwrap.dedent(body.strip("\n")).rstrip()

    lines = body.strip("\n").splitlines()
    remainder = textwrap.dedent("\n".join(lines[1:])) if lines[1:] else ""
    return "\n".join([lines[0].strip(), *remainder.splitlines()]).rstrip()


def _split_declaration(line: str) -> tuple[str, str] | None:
    """Return the ``(kind, name)`` a code line declares, or ``None``.

    Anonymous declarations — ``instance : Foo Bar``, ``example ...`` — have
    nothing to name or link, and are reported as ``None``.
    """
    tokens = line.split()
    position = 0
    while position < len(tokens) and tokens[position] in _MODIFIER_KEYWORDS:
        position += 1
    if position + 1 >= len(tokens) or tokens[position] not in _DECLARATION_KEYWORDS:
        return None
    name = _NAME_SEPARATOR.split(tokens[position + 1])[0]
    if not _NAME.match(name):
        return None
    return tokens[position], name


def _truncate_signature(line: str) -> tuple[str, bool]:
    """Cut ``line`` where the declaration body starts.

    Returns the retained text and whether the signature ended on this line.
    """
    cut = len(line)
    terminated = False
    body = line.find(_BODY_MARKER)
    if body != -1:
        cut, terminated = body, True
    where = _WHERE_MARKER.search(line)
    if where is not None and where.start() < cut:
        cut, terminated = where.start(), True
    return line[:cut].rstrip(), terminated


def _read_signature(lines: Sequence[str], start: int) -> tuple[str, int]:
    """Collect the signature beginning at ``lines[start]``.

    Signatures wrap freely across lines, so this consumes lines until the body
    marker appears or something that cannot be part of a signature does.
    """
    collected: list[str] = []
    index = start
    while index < len(lines):
        line = lines[index].strip()
        if index > start and (
            not line
            or line.startswith(_SIGNATURE_INTERRUPTERS)
            or _split_declaration(line) is not None
        ):
            break
        retained, terminated = _truncate_signature(line)
        if retained:
            collected.append(retained)
        index += 1
        if terminated:
            break
    return " ".join(collected), index


def parse_lean_source(name: str, source: str) -> LeanModule:
    """Parse one Lean source file into its documentable parts.

    Args:
        name: Fully-qualified module name, e.g. ``Theorems.Basic``.
        source: The file's full text.

    Returns:
        The module docstring plus every named declaration and section prose
        block, in source order. Declarations without a docstring are still
        included, so an undocumented definition is visible as a gap rather than
        silently missing.
    """
    blocks = list(_lex(source))
    entries: list[LeanEntry] = []
    module_docstring: str | None = None
    namespaces: list[str] = []
    pending_docstring: str | None = None

    position = 0
    while position < len(blocks):
        block = blocks[position]
        position += 1

        if block.kind == "module_doc":
            body = dedent_block(block.text)
            if module_docstring is None and not entries:
                module_docstring = body
            elif body:
                entries.append(LeanProse(body=body))
            continue
        if block.kind == "declaration_doc":
            pending_docstring = dedent_block(block.text)
            continue
        if block.kind == "comment":
            continue

        line = block.text.strip()
        if not line or line.startswith(("@[", "--")):
            # Attributes and blank lines may sit between a docstring and the
            # declaration it documents, so `pending_docstring` survives them.
            continue
        if line.startswith("namespace "):
            namespaces.append(line.split()[1])
            pending_docstring = None
            continue
        if line.startswith("end ") or line == "end":
            closing = line.removeprefix("end").strip()
            if namespaces and closing == namespaces[-1]:
                namespaces.pop()
            pending_docstring = None
            continue

        declaration = _split_declaration(line)
        if declaration is None:
            pending_docstring = None
            continue

        kind, short_name = declaration
        # `_read_signature` needs the surrounding code lines; the lexer already
        # guarantees consecutive code blocks are consecutive source lines.
        following: list[str] = []
        scan = position - 1
        while scan < len(blocks) and blocks[scan].kind == "code":
            following.append(blocks[scan].text)
            scan += 1
        signature, consumed = _read_signature(following, 0)
        position = position - 1 + consumed

        entries.append(
            LeanDeclaration(
                kind=kind,
                name=".".join([*namespaces, short_name]),
                signature=signature,
                docstring=pending_docstring or None,
            )
        )
        pending_docstring = None

    return LeanModule(
        name=name,
        docstring=module_docstring,
        entries=tuple(entries),
    )
