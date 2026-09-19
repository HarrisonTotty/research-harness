r"""Strand dynamics: a cellular automaton on the colors of a plabic graph.

The topology of a :class:`research.plabic_graph.PlabicGraph` — its rotation
system — is held fixed and only the black/white coloring of the internal
vertices evolves. The update rule, *strand gradient*, is driven by the
strands (trips) of the current coloring, so it couples every vertex to the
global routing: a flip changes the rules of the road at that vertex,
reroutes every strand through it, and thereby changes the quantity that
decides the next flips. The scheme is original to this project, not taken
from the literature; it was first explored as "Finding 020" of the
predecessor repository, whose ``strand_gradient`` rule this module ports
literally.

One synchronous step on a coloring ``c``:

1. **Load.** Trace the ``n`` one-way trips of ``c`` (right at black, left at
   white; Postnikov section 13). A trip's *length* is the number of darts it
   traverses. The load of an internal vertex is the sum of the lengths of
   the trips that visit it, a trip counting once however often it returns.
   Roundtrips carry no load.
2. **Gradient.** An internal vertex flips color iff its load is positive and
   strictly exceeds the mean load of its internal neighbors, one term per
   connecting half-edge. A vertex with no internal neighbor (a lollipop)
   never flips. The comparison is done in integers, ``load * degree >
   sum``, so no floating-point tie can arise.

The state space is finite and the rule deterministic, so every trajectory is
eventually periodic; :meth:`StrandAutomaton.orbit` measures its transient
and period, together with the period of the trip permutation along the
cycle. A recoloring generally leaves the positroid cell of the initial
graph: the trip permutation is an observable of the dynamics, not an
invariant.
"""

import dataclasses
import functools
from collections.abc import Hashable, Sequence
from dataclasses import dataclass

from research.plabic_graph import BLACK, WHITE, PlabicGraph

__all__ = ["Orbit", "StrandAutomaton"]

type State = tuple[int, ...]
"""Colors (1 black, -1 white) aligned with ``PlabicGraph.internal_vertices``."""


@dataclass(frozen=True, slots=True)
class Orbit:
    """The eventually periodic shape of one trajectory.

    ``transient`` is the number of steps before the trajectory first enters
    its cycle and ``period`` the length of that cycle, both exact (the full
    coloring is compared, so a repeat guarantees identical futures).
    ``permutation_period`` is the minimal period of the trip-permutation
    sequence around the cycle — a divisor of ``period`` — and
    ``cycle_permutations`` the number of distinct trip permutations met on
    the cycle.
    """

    transient: int
    period: int
    permutation_period: int
    cycle_permutations: int


@dataclass(frozen=True, slots=True)
class _Tables:
    """Dart-indexed routing tables of a fixed rotation system.

    For a dart ``d``, ``head[d]`` is the index of the internal vertex it
    arrives at, or minus the boundary position (``1..n``) it ends on;
    ``right[d]`` and ``left[d]`` are the darts a trip continues on from
    there at a black, respectively white, vertex. ``starts`` holds the dart
    leaving each boundary vertex and ``neighbors`` the internal neighbors of
    each internal vertex, one entry per connecting half-edge.
    """

    head: tuple[int, ...]
    right: tuple[int, ...]
    left: tuple[int, ...]
    starts: tuple[int, ...]
    neighbors: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class StrandAutomaton[V: Hashable]:
    """The strand-gradient automaton on the colorings of ``graph``.

    States are tuples of colors aligned with ``graph.internal_vertices``;
    the automaton itself is immutable and every method is a pure function
    of the state passed in, defaulting to the graph's own coloring.
    """

    graph: PlabicGraph[V]

    @functools.cached_property
    def _tables(self) -> _Tables:
        """Build the routing tables (edge ``e`` has darts ``2e``, ``2e + 1``)."""
        index = {v: i for i, v in enumerate(self.graph.internal_vertices)}
        position = {b: -p for p, b in enumerate(self.graph.boundary, start=1)}
        size = 2 * len(self.graph.edges)
        head, right, left = [0] * size, [0] * size, [0] * size
        for v, darts in self.graph.rotations:
            for slot, dart in enumerate(darts):
                arriving = dart ^ 1
                head[arriving] = index[v] if v in index else position[v]
                right[arriving] = darts[(slot + 1) % len(darts)]
                left[arriving] = darts[slot - 1]
        neighbors = tuple(
            tuple(head[dart] for dart in darts if head[dart] >= 0)
            for v, darts in self.graph.rotations
            if v in index
        )
        starts = tuple(darts[0] for v, darts in self.graph.rotations if v in position)
        return _Tables(tuple(head), tuple(right), tuple(left), starts, neighbors)

    @property
    def initial_state(self) -> State:
        """The coloring carried by ``graph`` itself."""
        return tuple(color for _, color in self.graph.colors)

    def _resolve(self, state: Sequence[int] | None) -> State:
        """Default ``state`` to the graph's coloring and validate it."""
        if state is None:
            return self.initial_state
        resolved = tuple(state)
        expected = len(self.graph.colors)
        if len(resolved) != expected or any(c not in (BLACK, WHITE) for c in resolved):
            msg = (
                f"a state assigns 1 (black) or -1 (white) to each of the "
                f"{expected} internal vertices; got {resolved!r}"
            )
            raise ValueError(msg)
        return resolved

    def _sweep(self, state: State) -> tuple[list[int], tuple[int, ...]]:
        """Trace the one-way trips of ``state``: vertex loads and trip permutation."""
        tables = self._tables
        head, right, left = tables.head, tables.right, tables.left
        load = [0] * len(state)
        permutation: list[int] = []
        for start in tables.starts:
            dart, length = start, 1
            visited: set[int] = set()
            while (vertex := head[dart]) >= 0:
                visited.add(vertex)
                dart = right[dart] if state[vertex] == BLACK else left[dart]
                length += 1
            permutation.append(-vertex)
            for seen in visited:
                load[seen] += length
        return load, tuple(permutation)

    def _advance(self, state: State) -> tuple[State, tuple[int, ...]]:
        """Return the successor of ``state`` and the trip permutation of ``state``."""
        load, permutation = self._sweep(state)
        successor = tuple(
            -color
            if load[v] > 0 and load[v] * len(around) > sum(load[u] for u in around)
            else color
            for v, (color, around) in enumerate(
                zip(state, self._tables.neighbors, strict=True)
            )
        )
        return successor, permutation

    def loads(self, state: Sequence[int] | None = None) -> tuple[int, ...]:
        """Return the trip load of every internal vertex under ``state``.

        Raises:
            ValueError: If ``state`` is not a coloring of the internal vertices.
        """
        return tuple(self._sweep(self._resolve(state))[0])

    def trip_permutation(self, state: Sequence[int] | None = None) -> tuple[int, ...]:
        """Return the trip permutation of ``state`` on boundary positions ``1..n``.

        Postnikov's direction, as
        :attr:`research.plabic_graph.PlabicGraph.trip_permutation`.

        Raises:
            ValueError: If ``state`` is not a coloring of the internal vertices.
        """
        return self._sweep(self._resolve(state))[1]

    def step(self, state: Sequence[int] | None = None) -> State:
        """Return the coloring one synchronous strand-gradient step after ``state``.

        Raises:
            ValueError: If ``state`` is not a coloring of the internal vertices.
        """
        return self._advance(self._resolve(state))[0]

    def recolored(self, state: Sequence[int]) -> PlabicGraph[V]:
        """Return ``graph`` with its internal vertices recolored by ``state``.

        Raises:
            ValueError: If ``state`` is not a coloring of the internal vertices.
        """
        colors = tuple(
            (v, color)
            for (v, _), color in zip(
                self.graph.colors, self._resolve(state), strict=True
            )
        )
        return dataclasses.replace(self.graph, colors=colors)

    def orbit(
        self, state: Sequence[int] | None = None, *, max_steps: int
    ) -> Orbit | None:
        """Follow the trajectory of ``state`` until a coloring repeats.

        Args:
            state: The starting coloring; the graph's own by default.
            max_steps: The most rule applications to spend.

        Returns:
            The orbit, or ``None`` if no coloring has repeated after
            ``max_steps`` steps. Every visited coloring is kept in memory.

        Raises:
            ValueError: If ``state`` is not a coloring of the internal
                vertices or ``max_steps`` is negative.
        """
        if max_steps < 0:
            msg = f"max_steps must be non-negative; got {max_steps}"
            raise ValueError(msg)
        current = self._resolve(state)
        first_seen: dict[State, int] = {current: 0}
        permutations: list[tuple[int, ...]] = []
        for time in range(1, max_steps + 1):
            current, permutation = self._advance(current)
            permutations.append(permutation)
            if current in first_seen:
                transient = first_seen[current]
                cycle = permutations[transient:]
                return Orbit(
                    transient=transient,
                    period=time - transient,
                    permutation_period=_minimal_period(cycle),
                    cycle_permutations=len(set(cycle)),
                )
            first_seen[current] = time
        return None


def _minimal_period[T](cycle: Sequence[T]) -> int:
    """Return the least divisor ``p`` of ``len(cycle)`` that is a period of it."""
    size = len(cycle)
    return next(
        p
        for p in range(1, size + 1)
        if size % p == 0 and all(cycle[i] == cycle[i - p] for i in range(p, size))
    )
