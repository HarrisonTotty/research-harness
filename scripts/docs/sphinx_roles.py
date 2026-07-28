"""Griffe extension rewriting Sphinx docstring roles into cross-references.

The Python sources document themselves with Google-style sections but reach for
reStructuredText inline roles — ``:func:`experiment``` , ``:class:`pandas.DataFrame```
— to cross-reference other objects. Markdown has no notion of those, so without
this they render as literal text.

Each role is rewritten into the ``mkdocs-autorefs`` link syntax, which
mkdocstrings resolves against the surrounding scope for project symbols (via the
``scoped_crossrefs`` option) and against the configured inventories for
third-party ones. Roles pointing at private names become plain code spans
instead: private members are filtered out of the generated pages, so a link to
one would dangle.
"""

import re
from collections.abc import Iterable
from typing import override

from griffe import Extension, Object

__all__ = ["SphinxRoles"]

_ROLE = re.compile(
    r":(?:py:)?(?:mod|module|func|function|class|exc|meth|method"
    r"|attr|attribute|data|const|obj|type):`(?P<body>[^`]+)`"
)
"""A Sphinx cross-reference role and its target."""

_EXPLICIT_TITLE = re.compile(r"^(?P<title>.+?)\s*<(?P<target>[^<>]+)>$")
"""The ``:role:`Text <target>``` form, which names its own link text."""

_ABBREVIATE = "~"
"""Sphinx's prefix for "show only the last component of the target"."""


def _is_private(target: str) -> bool:
    """Report whether any component of ``target`` is a private name.

    Dunders are public API by convention, so ``__init__`` does not count.
    """
    return any(
        part.startswith("_") and not (part.startswith("__") and part.endswith("__"))
        for part in target.split(".")
    )


class SphinxRoles(Extension):
    """Rewrite Sphinx roles in every docstring Griffe loads.

    Args:
        literal_targets: Targets to render as code spans rather than links.
            Some published inventories omit entries a docstring nonetheless
            refers to — Click, for instance, documents ``click.Command`` but
            registers no ``click`` module — and a link to one of those can never
            resolve.
    """

    def __init__(self, literal_targets: Iterable[str] = ()) -> None:
        """Store the targets that must not be turned into links."""
        super().__init__()
        self._literal_targets = frozenset(literal_targets)

    def _rewrite(self, match: re.Match[str]) -> str:
        """Turn a single matched role into Markdown."""
        body = match.group("body").strip()

        explicit = _EXPLICIT_TITLE.match(body)
        if explicit is not None:
            title, target = explicit.group("title"), explicit.group("target")
        elif body.startswith(_ABBREVIATE):
            target = body.removeprefix(_ABBREVIATE)
            title = target.rsplit(".", maxsplit=1)[-1]
        else:
            title = target = body

        if target in self._literal_targets or _is_private(target):
            return f"`{title}`"
        return f"[{title}][{target}]"

    @override
    def on_object(self, *, obj: Object, **kwargs: object) -> None:
        """Rewrite the roles in ``obj``'s docstring, if it has one."""
        docstring = obj.docstring
        if docstring is not None:
            docstring.value = _ROLE.sub(self._rewrite, docstring.value)
