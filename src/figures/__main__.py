"""Command-line entry point dispatching named figure generators.

Invoked as ``python -m figures <name> [args ...]`` (or via ``just figure``).
The launcher maps ``<name>`` to the module ``figures.<name>`` — hyphens in the
name become underscores — imports it, and invokes the single :mod:`click`
command it defines (the one produced by decorating a function with
:func:`figures.figure`). Arguments after the name are forwarded to that
command unmodified.
"""

import importlib
import pkgutil
import sys
from collections.abc import Sequence

import click

_PACKAGE: str = __package__ or "figures"
"""Import path of this package; ``__package__`` is optional-typed but always set."""

_LIBRARY_MODULES: frozenset[str] = frozenset({"cli", "style"})
"""Modules in this package that are the shared harness, not figures."""


class _FigureError(Exception):
    """A named figure could not be resolved to a single command."""


def _available_figures() -> list[str]:
    """Return the sorted names of importable figure modules."""
    package = importlib.import_module(_PACKAGE)
    return sorted(
        info.name
        for info in pkgutil.iter_modules(package.__path__)
        if info.name not in _LIBRARY_MODULES and not info.name.startswith("_")
    )


def _resolve_command(name: str) -> click.Command:
    """Resolve a figure name to the command its module defines.

    Args:
        name: Figure name as given on the command line; hyphens are treated
            as underscores when locating the module.

    Returns:
        The single :class:`click.Command` defined by ``figures.<name>``.

    Raises:
        _FigureError: If no such module exists, or it does not define exactly
            one command.
    """
    module_name = f"{_PACKAGE}.{name.replace('-', '_')}"
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name != module_name:
            raise  # a dependency of the figure is missing — a real error
        available = _available_figures()
        listing = "\n".join(f"  - {n}" for n in available) or "  (none found)"
        msg = f"no figure named {name!r}\n\navailable figures:\n{listing}"
        raise _FigureError(msg) from exc

    commands: list[click.Command] = []
    seen: set[int] = set()
    for value in vars(module).values():
        if isinstance(value, click.Command) and id(value) not in seen:
            seen.add(id(value))
            commands.append(value)
    if not commands:
        msg = (
            f"figure module {module_name!r} defines no click command; "
            "decorate a function with @figure"
        )
        raise _FigureError(msg)
    if len(commands) > 1:
        msg = (
            f"figure module {module_name!r} defines multiple click commands; "
            "expected exactly one"
        )
        raise _FigureError(msg)
    return commands[0]


def main(argv: Sequence[str] | None = None) -> int:
    """Dispatch the named figure generator and return a process exit code.

    Reads ``argv`` when provided, otherwise the process arguments. The first
    positional value selects the figure; any following values are forwarded
    to its command unmodified. The command runs in :mod:`click`'s standalone
    mode, so its own exit code (including usage errors and ``--help``) is
    surfaced as this function's return value.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("usage: figure <name> [args ...]", file=sys.stderr)
        return 2
    name, *rest = args
    try:
        command = _resolve_command(name)
    except _FigureError as exc:
        print(exc, file=sys.stderr)
        return 2
    # ``Command.main`` runs in standalone mode and always terminates via
    # ``SystemExit`` — on success, usage error, or ``--help`` alike — so its exit
    # code is recovered here rather than returned.
    try:
        command.main(args=rest, prog_name=f"figure {name}")
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        print(code, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
