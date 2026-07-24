"""Command-line entry point dispatching named experiments.

Invoked as ``python -m experiments <name> [args ...]`` (or via ``just
experiment``). Each experiment is a self-contained module; this launcher only
resolves the name and forwards the remaining arguments.
"""

import sys
from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    """Dispatch the named experiment and return a process exit code.

    Reads ``argv`` when provided, otherwise the process arguments. The first
    positional value selects the experiment; any following values are passed
    through to it unmodified.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("usage: experiment <name> [args ...]", file=sys.stderr)
        return 2
    name, *rest = args
    # Experiment resolution and invocation are implemented as experiments are
    # added under this package.
    raise NotImplementedError(f"no experiment registered under {name!r} (args: {rest})")


if __name__ == "__main__":
    raise SystemExit(main())
