"""Strand dynamics on the canonical graphs of every fixed-point-free cell.

Implements the design doc ``docs/experiments/strand-dynamics.md``: for every
fixed-point-free decorated permutation of ``[n]``, ``n`` over the swept
range, build the bridge graph and the Le-graph of its positroid cell, run
the strand-gradient automaton of :mod:`research.strand_dynamics` from each
graph's own coloring, and record the transient and period of the trajectory
— one row per (cell, graph). The bridge arm replicates "Finding 020" of the
predecessor repository; the Le arm asks what survives a change of
representative. The run is deterministic, so there are no seeds.
"""

import enum
import importlib.metadata
import itertools
import platform
import shutil
import subprocess
import time
from collections.abc import Hashable, Iterator, Sequence
from pathlib import Path

import click
import pandas as pd

from experiments.cli import ExperimentContext, experiment
from research.decorated_permutation import DecoratedPermutation
from research.le_diagram import LeDiagram
from research.plabic_graph import PlabicGraph
from research.strand_dynamics import StrandAutomaton

_NULLABLE_COLUMNS = ("transient", "period", "permutation_period", "cycle_permutations")
"""Orbit columns, null on a row whose trajectory did not close within budget."""


class GraphKind(enum.StrEnum):
    """The canonical graph chosen to represent a positroid cell."""

    BRIDGE = "bridge"
    LE = "le"


def _derangements(n: int) -> Iterator[DecoratedPermutation]:
    """Yield the fixed-point-free decorated permutations of ``[n]``, lazily.

    In lexicographic order of the one-line images; factorial in ``n`` and
    single-use.
    """
    for targets in itertools.permutations(range(1, n + 1)):
        if all(target != i for i, target in enumerate(targets, start=1)):
            yield DecoratedPermutation(targets)


def _canonical_graph(
    decorated: DecoratedPermutation, kind: GraphKind
) -> PlabicGraph[Hashable]:
    """Return the graph of ``kind`` whose trip permutation is ``decorated``.

    The Le-graph of a Le-diagram has the *inverse* of the diagram's
    decorated permutation as its trip permutation (Postnikov Corollary 20.1,
    as recorded on :meth:`research.le_diagram.LeDiagram.to_plabic_graph`),
    so the diagram is built from ``decorated.inverse()``.
    """
    if kind is GraphKind.BRIDGE:
        return PlabicGraph.from_bridge_decomposition(decorated)
    return LeDiagram.from_decorated_permutation(decorated.inverse()).to_plabic_graph()


def _is_noncrossing_matching(decorated: DecoratedPermutation) -> bool:
    """Return whether ``decorated`` is an involution with pairwise non-crossing arcs."""
    targets = decorated.targets
    arcs = [(i, t) for i, t in enumerate(targets, start=1) if i < t]
    if any(targets[t - 1] != i for i, t in arcs) or 2 * len(arcs) != len(targets):
        return False
    return not any(a < c < b < d for (a, b), (c, d) in itertools.permutations(arcs, 2))


def _row(
    decorated: DecoratedPermutation, kind: GraphKind, max_steps: int
) -> dict[str, object]:
    """Run one (cell, graph) unit and return its result row."""
    graph = _canonical_graph(decorated, kind)
    if graph.trip_permutation != decorated.targets:
        msg = (
            f"the {kind.value} graph of {decorated.targets!r} has trip "
            f"permutation {graph.trip_permutation!r}; the cell was not built"
        )
        raise RuntimeError(msg)
    orbit = StrandAutomaton(graph).orbit(max_steps=max_steps)
    return {
        "n": decorated.size,
        "k": decorated.anti_exceedance_count,
        "permutation": ",".join(map(str, decorated.targets)),
        "graph": kind.value,
        "dimension": decorated.dimension,
        "internal_vertices": len(graph.internal_vertices),
        "edges": len(graph.edges),
        "noncrossing_matching": _is_noncrossing_matching(decorated),
        "converged": orbit is not None,
        "transient": None if orbit is None else orbit.transient,
        "period": None if orbit is None else orbit.period,
        "permutation_period": None if orbit is None else orbit.permutation_period,
        "cycle_permutations": None if orbit is None else orbit.cycle_permutations,
    }


def _git_commit() -> str | None:
    """Return the checked-out commit, suffixed ``-dirty`` on a modified tree.

    Untracked files count as modifications, so this must be read before the
    run writes its own artifacts — otherwise every run marks itself dirty.
    """
    git = shutil.which("git")
    if git is None:
        return None
    here = Path(__file__).parent
    try:
        commit = subprocess.run(  # noqa: S603 - fixed argv, resolved executable
            [git, "rev-parse", "HEAD"],
            cwd=here,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
        status = subprocess.run(  # noqa: S603 - fixed argv, resolved executable
            [git, "status", "--porcelain"],
            cwd=here,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
    except subprocess.SubprocessError:
        return None
    return f"{commit}-dirty" if status else commit


@experiment()
@click.option(
    "--n-min",
    type=click.IntRange(min=2),
    default=2,
    show_default=True,
    help="Smallest number of boundary vertices swept.",
)
@click.option(
    "--n-max",
    type=click.IntRange(min=2),
    default=8,
    show_default=True,
    help="Largest number of boundary vertices swept (9 takes tens of minutes).",
)
@click.option(
    "--graph",
    "graphs",
    type=click.Choice([kind.value for kind in GraphKind]),
    multiple=True,
    default=[kind.value for kind in GraphKind],
    show_default=True,
    help="Canonical graph(s) to run each cell on; repeat the option for several.",
)
@click.option(
    "--max-steps",
    type=click.IntRange(min=1),
    default=100_000,
    show_default=True,
    help="Rule applications allowed per trajectory before it counts as unconverged.",
)
def strand_dynamics(
    ctx: ExperimentContext,
    n_min: int,
    n_max: int,
    graphs: Sequence[str],
    max_steps: int,
) -> None:
    """Measure strand-gradient transients and periods over all derangement cells."""
    if n_min > n_max:
        msg = f"--n-min ({n_min}) must not exceed --n-max ({n_max})"
        raise click.UsageError(msg)
    kinds = [kind for kind in GraphKind if kind.value in graphs]
    commit = _git_commit()
    started = time.monotonic()
    rows: list[dict[str, object]] = []
    for n in range(n_min, n_max + 1):
        before = len(rows)
        for decorated in _derangements(n):
            rows.extend(_row(decorated, kind, max_steps) for kind in kinds)
        ctx.logger.info(
            "n=%d: %d rows (%.1fs elapsed)",
            n,
            len(rows) - before,
            time.monotonic() - started,
        )
    frame = pd.DataFrame(rows)
    for column in _NULLABLE_COLUMNS:
        frame[column] = frame[column].astype("Int64")
    unconverged = int((~frame["converged"]).sum())
    if unconverged:
        ctx.logger.warning(
            "%d rows did not converge in %d steps", unconverged, max_steps
        )
    result_path = ctx.write_result(frame)
    metadata_path = ctx.write_metadata(
        {
            "experiment": ctx.name,
            "design_doc": "docs/experiments/strand-dynamics.md",
            "timestamp": ctx.timestamp,
            "parameters": {
                "n_min": n_min,
                "n_max": n_max,
                "graphs": [kind.value for kind in kinds],
                "max_steps": max_steps,
            },
            "seeds": [],
            "deterministic": True,
            "git_commit": commit,
            "versions": {
                "python": platform.python_version(),
                "pandas": importlib.metadata.version("pandas"),
                "research-harness": importlib.metadata.version("research-harness"),
            },
            "rows": len(frame),
            "unconverged_rows": unconverged,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "result_path": result_path,
        }
    )
    ctx.logger.info("wrote %d rows to %s (%s)", len(frame), result_path, metadata_path)
