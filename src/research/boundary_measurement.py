"""The boundary measurement map: planar networks in a disk to positroid cells.

A planar directed network with positive edge weights determines, for each
boundary source ``b_i`` and sink ``b_j``, a boundary measurement ``M_ij``,
and those numbers are ratios of Plucker coordinates of a point of the
totally nonnegative Grassmannian (Postnikov, *Total positivity,
Grassmannians, and networks*, arXiv:math/0609764, 2006, Definitions
4.4-4.6). This module implements the two finite formulations the Boundary
Measurement Map page singles out: Talaska's flow formula on perfectly
oriented networks (:class:`PlanarNetwork`; Talaska, arXiv:0801.4822, 2008,
Theorem 3.2) and Lam's dimer partition function on planar bipartite
networks (:class:`PlanarBipartiteNetwork`; Lam, arXiv:1506.00603, 2015,
section 4.1), bridged by Lam's Propositions 5.1 and 5.3. Networks carry an
explicit exact-rational embedding, because planarity is essential to the
flow formula (Talaska, section 5) and the embedding also yields faces, face
weights, and trips. Matrix-level companions — the signed boundary
measurement matrix, Chevalley generators and bridge removal, the twist, and
positroid cell dimension — operate on the same label-to-column mappings
that :meth:`research.positroid.Positroid.from_matrix` consumes.

All arithmetic is exact (`Fraction`), the enumeration helpers are
exponential and guarded by an explicit size bound, and every formula
docstring cites the source the page attributes it to.
"""

import functools
import itertools
from collections.abc import Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Literal, override

import pandas as pd

from research._linalg import det_q
from research._plot import ensure_axes, scatter_labeled
from research.grassmann_necklace import GrassmannNecklace
from research.positroid import Positroid

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "PlanarBipartiteNetwork",
    "PlanarNetwork",
    "acyclic_baseline_network",
    "apply_chevalley",
    "bridge_parameter",
    "cell_dimension",
    "geometric_series_network",
    "has_bridge",
    "left_twist",
    "lollipop_network",
    "measurement_matrix",
    "remove_bridge",
    "right_twist",
    "square_network",
]

type Point = tuple[Fraction, Fraction]
type Columns[T: Hashable] = dict[T, tuple[Fraction, ...]]

_MAX_ENUMERATION_EDGES = 18
"""Flow and matching enumeration bound (page cost note: explicit guards)."""

_SPHERE_EULER_CHARACTERISTIC = 2
"""Euler characteristic certifying a genus-zero (disk) embedding."""

_SQUARE_SIDES = 4
"""The quadrilateral face size the square move operates on."""


# --------------------------------------------------------------------------- #
# Exact planar geometry
# --------------------------------------------------------------------------- #
def _cross(o: Point, a: Point, b: Point) -> Fraction:
    """Return the signed area cross product ``(a - o) x (b - o)``."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _on_segment(p: Point, a: Point, b: Point) -> bool:
    """Return whether ``p`` lies on the closed segment ``a``-``b``."""
    if _cross(a, b, p) != 0:
        return False
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[
        1
    ] <= max(a[1], b[1])


def _segments_conflict(a: Point, b: Point, c: Point, d: Point) -> bool:
    """Return whether segments ``a-b`` and ``c-d`` share a non-endpoint point.

    Segments sharing an endpoint conflict only if they overlap beyond it.
    Exact rational orientation tests; used to enforce that the network is
    drawn inside the disk without crossings (Postnikov Definition 4.1;
    planarity is essential per Talaska section 5).
    """
    shared = {a, b} & {c, d}
    if {a, b} == {c, d}:
        return True
    if shared:
        hinge = shared.pop()
        (p,) = {a, b} - {hinge}
        (q,) = {c, d} - {hinge}
        return (
            _cross(hinge, p, q) == 0
            and (p[0] - hinge[0]) * (q[0] - hinge[0])
            + (p[1] - hinge[1]) * (q[1] - hinge[1])
            > 0
        )
    d1 = _cross(c, d, a)
    d2 = _cross(c, d, b)
    d3 = _cross(a, b, c)
    d4 = _cross(a, b, d)
    if ((d1 > 0) != (d2 > 0) and d1 != 0 and d2 != 0) and (
        (d3 > 0) != (d4 > 0) and d3 != 0 and d4 != 0
    ):
        return True
    return (
        _on_segment(a, c, d)
        or _on_segment(b, c, d)
        or _on_segment(c, a, b)
        or _on_segment(d, a, b)
    )


def _direction_key(d: Point) -> tuple[int, Point]:
    """Return a sort key placing directions in counterclockwise order.

    Order starts at the positive x-axis; within a half-plane two directions
    compare by their cross product, exactly.
    """
    half = 0 if d[1] > 0 or (d[1] == 0 and d[0] > 0) else 1
    return (half, d)


def _sort_ccw[K](darts: list[tuple[K, Point]]) -> list[K]:
    """Sort dart keys by direction, counterclockwise from the x-axis."""

    def compare(u: tuple[K, Point], v: tuple[K, Point]) -> int:
        hu, du = _direction_key(u[1])
        hv, dv = _direction_key(v[1])
        if hu != hv:
            return -1 if hu < hv else 1
        c = du[0] * dv[1] - du[1] * dv[0]
        if c > 0:
            return -1
        if c < 0:
            return 1
        return 0

    return [key for key, _ in sorted(darts, key=functools.cmp_to_key(compare))]


def _solve_q(
    rows: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    """Solve a square rational linear system exactly by Gaussian elimination.

    Raises:
        ValueError: If the system is singular.
    """
    d = len(rows)
    m = [[*row, rhs[r]] for r, row in enumerate(rows)]
    for col in range(d):
        pivot_row = next((r for r in range(col, d) if m[r][col]), None)
        if pivot_row is None:
            msg = "singular linear system in twist computation"
            raise ValueError(msg)
        m[col], m[pivot_row] = m[pivot_row], m[col]
        pivot = m[col][col]
        m[col] = [a / pivot for a in m[col]]
        for r in range(d):
            if r != col and m[r][col]:
                factor = m[r][col]
                m[r] = [a - factor * b for a, b in zip(m[r], m[col], strict=True)]
    return tuple(m[r][d] for r in range(d))


# --------------------------------------------------------------------------- #
# The boundary measurement matrix and matrix-level operations
# --------------------------------------------------------------------------- #
def measurement_matrix(
    source_set: Iterable[int],
    n: int,
    measurements: Mapping[tuple[int, int], Fraction],
) -> tuple[tuple[Fraction, ...], ...]:
    """Return Postnikov's boundary measurement matrix ``A(N)``, as columns.

    Postnikov Definition 4.6: for source set ``I = {i_1 < ... < i_k}`` the
    matrix has the identity in the source columns and entry
    ``(-1)^s M_{i_r, j}`` in row ``r``, sink column ``j``, where ``s``
    counts the elements of ``I`` strictly between ``i_r`` and ``j``. That
    sign choice makes ``Delta_{(I - {i}) + {j}}(A(N)) = M_ij`` and
    ``Delta_I(A(N)) = 1`` (Postnikov, remark after Definition 4.6).

    Args:
        source_set: The source positions ``I``, a subset of ``1..n``.
        n: The number of boundary vertices.
        measurements: ``M_ij`` keyed by ``(i, j)`` for ``i`` a source and
            ``j`` a sink position; missing sink pairs default to zero.

    Returns:
        The ``n`` columns of ``A(N)``, index ``j - 1`` holding column ``j``.

    Raises:
        ValueError: If ``source_set`` is not a subset of ``1..n``.
    """
    sources = sorted(set(source_set))
    if any(i < 1 or i > n for i in sources):
        msg = f"source set {sources} must be a subset of 1..{n}"
        raise ValueError(msg)
    k = len(sources)
    row_of = {i: r for r, i in enumerate(sources)}
    columns: list[tuple[Fraction, ...]] = []
    for j in range(1, n + 1):
        if j in row_of:
            columns.append(
                tuple(Fraction(1) if r == row_of[j] else Fraction(0) for r in range(k))
            )
            continue
        entries: list[Fraction] = []
        for i in sources:
            lo, hi = min(i, j), max(i, j)
            s = sum(1 for other in sources if lo < other < hi)
            value = measurements.get((i, j), Fraction(0))
            entries.append(-value if s % 2 else value)
        columns.append(tuple(entries))
    return tuple(columns)


def apply_chevalley[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]],
    i: int,
    a: Fraction,
    *,
    variant: Literal["x", "y"] = "x",
) -> Columns[T]:
    """Return the point ``X . x_i(a)`` (or ``X . y_i(a)``), as columns.

    Lam Lemma 7.6: adding a bridge of weight ``a`` from ``i`` to ``i + 1``
    acts on the point by the Chevalley generator ``x_i(a)`` (white at
    ``i``), i.e. column ``i + 1`` gains ``a`` times column ``i``; the
    opposite colouring gives ``y_i(a)``, column ``i`` gaining ``a`` times
    column ``i + 1``. Columns are keyed by ground-set label in cyclic
    order; ``i`` is a position in ``1..n``.

    Raises:
        ValueError: If ``i`` is not in ``1..n-1`` (the wrap-around bridge
            at ``i = n`` is not implemented).
    """
    labels = list(columns)
    if not 1 <= i < len(labels):
        msg = f"bridge position {i} must satisfy 1 <= i < n = {len(labels)}"
        raise ValueError(msg)
    source, target = (i, i + 1) if variant == "x" else (i + 1, i)
    src, dst = labels[source - 1], labels[target - 1]
    updated = {
        label: tuple(Fraction(entry) for entry in column)
        for label, column in columns.items()
    }
    updated[dst] = tuple(
        t + a * s for t, s in zip(updated[dst], updated[src], strict=True)
    )
    return updated


def _bounded_affine_length(f: Sequence[int]) -> int:
    """Return the Coxeter length of a bounded affine permutation.

    The inversion count ``#{(i, j) : i < j, f(i) > f(j)}`` with ``i`` in
    ``1..n`` and ``j`` ranging over the integers, finite because
    ``j <= f(j)`` (Lam section 7 conventions).
    """
    n = len(f)

    def value(j: int) -> int:
        return f[(j - 1) % n] + n * ((j - 1) // n)

    return sum(
        1 for i in range(1, n + 1) for j in range(i + 1, i + n) if value(i) > value(j)
    )


def cell_dimension[T: Hashable](positroid: Positroid[T]) -> int:
    """Return the dimension of the positroid's totally nonnegative cell.

    Lam Theorem 7.12(4): ``dim = k(n - k) - l(f_M)`` where ``f_M`` is the
    bounded affine permutation of the positroid, reached here through the
    Grassmann necklace pipeline (Lam Theorem 7.12(2)).
    """
    n = positroid.size
    k = positroid.rank()
    necklace = GrassmannNecklace.from_matroid(positroid)
    return k * (n - k) - _bounded_affine_length(
        necklace.to_bounded_affine_permutation()
    )


def has_bridge[T: Hashable](positroid: Positroid[T], i: int) -> bool:
    """Return whether the cell's points have a bridge at position ``i``.

    Lam section 7: ``X`` has a bridge at ``i`` when its bounded affine
    permutation satisfies ``i < i + 1 <= f(i) < f(i + 1) <= i + n``.

    Raises:
        ValueError: If ``i`` is not in ``1..n-1`` (the wrap-around case is
            not implemented).
    """
    n = positroid.size
    if not 1 <= i < n:
        msg = f"bridge position {i} must satisfy 1 <= i < n = {n}"
        raise ValueError(msg)
    f = GrassmannNecklace.from_matroid(positroid).to_bounded_affine_permutation()
    return i + 1 <= f[i - 1] < f[i] <= i + n


def bridge_parameter[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]], i: int
) -> Fraction:
    """Return Lam's bridge parameter ``a`` at position ``i`` of the point.

    Lam Proposition 7.10: for totally nonnegative ``X`` with a bridge at
    ``i``, ``a = Delta_{I_{i+1}}(X) / Delta_{I_{i+1} + {i} - {i+1}}(X)`` is
    positive and well defined, where ``I_{i+1}`` is the Grassmann necklace
    entry at ``i + 1``.

    Raises:
        ValueError: If the point is not totally nonnegative or has no
            bridge at ``i``.
    """
    labels = list(columns)
    exact = {
        label: tuple(Fraction(entry) for entry in column)
        for label, column in columns.items()
    }
    positroid = Positroid.from_matrix(exact)
    if not has_bridge(positroid, i):
        msg = f"the point has no bridge at position {i} (Lam section 7)"
        raise ValueError(msg)
    necklace = positroid.grassmann_necklace[i % len(labels)]
    positions = sorted(labels.index(label) + 1 for label in necklace)
    swapped = sorted((set(positions) | {i}) - {i + 1})

    def minor(position_set: list[int]) -> Fraction:
        rows = [
            [exact[labels[p - 1]][r] for p in position_set]
            for r in range(len(position_set))
        ]
        return det_q(rows)

    return minor(positions) / minor(swapped)


def remove_bridge[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]], i: int
) -> Columns[T]:
    """Return ``X . x_i(-a)``, the bridge-removed point.

    Lam Proposition 7.10: removing the bridge keeps the point totally
    nonnegative and lands in a strictly smaller positroid cell, with
    ``f_{X'} = f_X s_i``.
    """
    return apply_chevalley(columns, i, -bridge_parameter(columns, i))


# --------------------------------------------------------------------------- #
# The twist (Muller-Speyer, zero-column-free form)
# --------------------------------------------------------------------------- #
def _twist[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]], *, forward: bool
) -> Columns[T]:
    """Compute a twist by exact linear solves; see the public wrappers."""
    exact: Columns[T] = {
        label: tuple(Fraction(entry) for entry in column)
        for label, column in columns.items()
    }
    labels = list(exact)
    n = len(labels)
    k = len(next(iter(exact.values()), ()))
    for label in labels:
        if all(entry == 0 for entry in exact[label]):
            msg = (
                f"column {label!r} is zero; the twist is implemented only "
                f"under Muller-Speyer's section 1.8 simplifying assumption "
                f"that the matrix has no zero columns"
            )
            raise ValueError(msg)
    result: Columns[T] = {}
    for start, label in enumerate(labels):
        step = 1 if forward else -1
        basis = [exact[label]]
        constraints: list[tuple[Fraction, ...]] = []
        offset = step
        while len(basis) < k:
            candidate = exact[labels[(start + offset) % n]]
            if _rank([*basis, candidate]) > len(basis):
                basis.append(candidate)
                constraints.append(candidate)
            offset += step
        rows = [exact[label], *constraints]
        rhs = [Fraction(1)] + [Fraction(0)] * (k - 1)
        result[label] = _solve_q(rows, rhs)
    return result


def _rank(vectors: Sequence[tuple[Fraction, ...]]) -> int:
    """Return the rank of the vectors by exact elimination."""
    m = [list(v) for v in vectors]
    rank = 0
    width = len(m[0]) if m else 0
    for col in range(width):
        pivot_row = next((r for r in range(rank, len(m)) if m[r][col]), None)
        if pivot_row is None:
            continue
        m[rank], m[pivot_row] = m[pivot_row], m[rank]
        pivot = m[rank]
        for r in range(len(m)):
            if r != rank and m[r][col]:
                factor = m[r][col] / pivot[col]
                m[r] = [a - factor * b for a, b in zip(m[r], pivot, strict=True)]
        rank += 1
    return rank


def right_twist[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]],
) -> Columns[T]:
    """Return the right twist of a point, column by column.

    Muller-Speyer section 1.8, under their stated simplifying assumption of
    no zero columns: column ``i`` of the right twist pairs to ``1`` with
    ``A_i`` and to ``0`` with each later ``A_j`` (indices cyclic) that
    leaves the span of ``A_i, ..., A_{j-1}``. On the ``Gr(2, n)`` uniform
    positroid the twisted Pluckers are cyclic shifts up to monomials
    (Muller-Speyer Appendix B).

    Raises:
        ValueError: If some column is zero.
    """
    return _twist(columns, forward=True)


def left_twist[T: Hashable](
    columns: Mapping[T, Sequence[Fraction | int]],
) -> Columns[T]:
    """Return the left twist of a point, column by column.

    Muller-Speyer section 1.8: as :func:`right_twist` with the reversed
    spans. The two twists are mutually inverse on the open positroid
    variety (Muller-Speyer Theorem 6.7).

    Raises:
        ValueError: If some column is zero.
    """
    return _twist(columns, forward=False)


# --------------------------------------------------------------------------- #
# Planar directed networks (Postnikov Definition 4.1, Talaska Definition 2.1)
# --------------------------------------------------------------------------- #
type _Dart = tuple[str, int, int]
"""A dart: kind ("e" real / "f" frame), edge or arc index, +-1 direction."""


@dataclass(frozen=True)
class PlanarNetwork[V: Hashable]:
    """A planar directed network in a disk with positive rational weights.

    Postnikov Definition 4.1 restricted to Talaska's Definition 2.1
    variant: boundary vertices are labelled clockwise and each is incident
    to at most one edge; isolated boundary vertices are declared sources
    or sinks explicitly.

    **What the stored embedding is and is not.** The exact straight-line
    coordinates serve three purposes: no two edges may cross, which is
    what makes the flow formula applicable at all (Talaska section 5);
    the angular order of the edges at each vertex gives the rotation
    system, from which faces, face weights, and trips are read; and the
    frame of boundary arcs closes that rotation system into a sphere,
    whose Euler characteristic is checked, so the clockwise labelling has
    to be one the drawing can realize. It is **not** required to be a
    literal picture of the disk: the boundary vertices need not enclose
    the internal ones, and in the Le-graphs they do not, since the wires
    of a Gamma-graph reach the boundary only on the right and below. The
    disk is therefore combinatorial — the genus-zero embedding with the
    frame — and every derived quantity depends only on that.

    Representation restrictions, all narrowings of Postnikov's class:
    no loops or repeated parallel edges (subdivide with a weight-1 vertex
    instead — path and flow weights are unchanged), and no components
    disconnected from the boundary circle.

    Build through :meth:`from_edges` (which validates) or the named example
    constructors; calling the dataclass constructor directly skips
    validation. Boundary measurements are computed only by the finite
    routes — the acyclic path sum (Postnikov section 4) and the flow
    formula (Talaska Theorem 3.2) — never the signed infinite series.
    """

    boundary: tuple[V, ...]
    edges: tuple[tuple[V, V, Fraction], ...]
    positions: tuple[tuple[V, Point], ...]
    isolated_source_set: frozenset[V] = frozenset()

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        k, n = self.network_type
        return (
            f"PlanarNetwork(n={n}, k={k}, edges={len(self.edges)}, "
            f"internal={len(self.internal_vertices)})"
        )

    # ---------------------------------------------------------- construction
    @classmethod
    def from_edges(
        cls,
        boundary: Sequence[V],
        positions: Mapping[V, tuple[Fraction | int, Fraction | int]],
        edges: Iterable[tuple[V, V, Fraction | int]],
        *,
        isolated_sources: Iterable[V] = (),
    ) -> PlanarNetwork[V]:
        """Build and validate a planar directed network.

        Args:
            boundary: Boundary vertex labels in clockwise order (Postnikov
                Definition 4.1).
            positions: Exact coordinates for every vertex; the embedding.
            edges: ``(tail, head, weight)`` triples with positive weights.
            isolated_sources: Isolated boundary vertices declared sources
                (Postnikov Definition 4.1 declares every boundary vertex,
                including isolated ones); undeclared isolated vertices are
                sinks.

        Returns:
            The validated network.

        Raises:
            ValueError: If a weight is nonpositive (Definition 4.1), a
                boundary vertex meets more than one edge (Talaska
                Definition 2.1), edges repeat or loop, segments cross or
                pass through vertices (the network must be drawn in the
                disk without crossings), the boundary order is not
                realized by the embedding, or the graph is disconnected
                from the boundary circle.
        """
        rational_edges = tuple(
            (tail, head, Fraction(weight)) for tail, head, weight in edges
        )
        placed = {v: (Fraction(p[0]), Fraction(p[1])) for v, p in positions.items()}
        ordered = [*boundary, *(v for v in placed if v not in set(boundary))]
        network = cls(
            tuple(boundary),
            rational_edges,
            tuple((v, placed[v]) for v in ordered),
            frozenset(isolated_sources),
        )
        network._validate()
        return network

    def _validate(self) -> None:
        """Check the definition; raise ``ValueError`` naming the condition."""
        self._validate_edges()
        self._validate_boundary()
        self._validate_embedding()

    def _validate_edges(self) -> None:
        """Check weights, loops, repeats, and placement (Definition 4.1)."""
        pos = self._pos
        if len(set(self.boundary)) != len(self.boundary):
            msg = f"boundary labels must be distinct, got {self.boundary!r}"
            raise ValueError(msg)
        if len({p for _, p in self.positions}) != len(self.positions):
            msg = "vertex positions must be distinct points"
            raise ValueError(msg)
        seen: set[tuple[V, V]] = set()
        for tail, head, weight in self.edges:
            if weight <= 0:
                msg = (
                    f"edge ({tail!r}, {head!r}) has weight {weight}; a "
                    f"network requires strictly positive weights "
                    f"(Postnikov Definition 4.1)"
                )
                raise ValueError(msg)
            if tail == head:
                msg = (
                    f"loop at {tail!r}: loops are not representable with "
                    f"straight-line edges; subdivide with a weight-1 vertex"
                )
                raise ValueError(msg)
            if tail not in pos or head not in pos:
                msg = f"edge ({tail!r}, {head!r}) has an unplaced endpoint"
                raise ValueError(msg)
            if (tail, head) in seen or (head, tail) in seen:
                msg = (
                    f"repeated edge between {tail!r} and {head!r}: parallel "
                    f"edges are not representable with straight-line edges; "
                    f"subdivide with a weight-1 vertex"
                )
                raise ValueError(msg)
            seen.add((tail, head))

    def _validate_boundary(self) -> None:
        """Check the boundary vertex conditions (Talaska Definition 2.1)."""
        pos = self._pos
        for b in self.boundary:
            if b not in pos:
                msg = f"boundary vertex {b!r} has no position"
                raise ValueError(msg)
            degree = sum(1 for t, h, _ in self.edges if b in (t, h))
            if degree > 1:
                msg = (
                    f"boundary vertex {b!r} meets {degree} edges; each "
                    f"boundary vertex is incident to at most one edge "
                    f"(Talaska Definition 2.1)"
                )
                raise ValueError(msg)
        stray = self.isolated_source_set - set(self.boundary)
        if stray:
            msg = (
                f"declared isolated sources {sorted(map(repr, stray))} "
                f"are not boundary vertices"
            )
            raise ValueError(msg)
        for b in self.isolated_source_set:
            if any(b in (t, h) for t, h, _ in self.edges):
                msg = (
                    f"boundary vertex {b!r} is declared an isolated source "
                    f"but is incident to an edge (Postnikov Definition 4.1 "
                    f"declares only isolated vertices this way)"
                )
                raise ValueError(msg)

    def _validate_embedding(self) -> None:
        """Check planarity, connectivity, and the disk boundary order.

        Note what is *not* checked, per the class docstring: the boundary
        vertices are not required to be in convex position, nor to enclose
        the internal ones. Only the genus-zero condition on the
        frame-augmented rotation system constrains the clockwise
        labelling, and that is the condition every derived quantity
        actually rests on.
        """
        pos = self._pos
        segments = [(pos[t], pos[h]) for t, h, _ in self.edges]
        for (a, b), (c, d) in itertools.combinations(segments, 2):
            if _segments_conflict(a, b, c, d):
                msg = (
                    "edges cross in the drawing; the network must be drawn "
                    "inside the disk without crossings (Postnikov "
                    "Definition 4.1; planarity is essential, Talaska "
                    "section 5)"
                )
                raise ValueError(msg)
        for v, p in self.positions:
            for (a, b), (t, h, _) in zip(segments, self.edges, strict=True):
                if v not in (t, h) and _on_segment(p, a, b):
                    msg = f"vertex {v!r} lies on the edge ({t!r}, {h!r})"
                    raise ValueError(msg)
        reached = self._frame_reachable()
        missing = {v for v, _ in self.positions} - reached
        if missing:
            msg = (
                f"vertices {sorted(map(repr, missing))} are disconnected "
                f"from the boundary circle; components floating in the "
                f"disk are not supported (representation restriction)"
            )
            raise ValueError(msg)
        faces = self._faces
        vertex_count = len(self.positions)
        edge_count = len(self.edges) + len(self.boundary)
        euler = vertex_count - edge_count + len(faces)
        if euler != _SPHERE_EULER_CHARACTERISTIC:
            msg = (
                "the embedding does not realize the boundary circle in the "
                "given clockwise order (Postnikov Definition 4.1)"
            )
            raise ValueError(msg)

    def _frame_reachable(self) -> set[V]:
        """Return vertices reachable from the boundary via edges or arcs."""
        adjacency: dict[V, set[V]] = {v: set() for v, _ in self.positions}
        for tail, head, _ in self.edges:
            adjacency[tail].add(head)
            adjacency[head].add(tail)
        for a, b in zip(
            self.boundary, self.boundary[1:] + self.boundary[:1], strict=True
        ):
            adjacency[a].add(b)
            adjacency[b].add(a)
        frontier = list(self.boundary)
        reached = set(frontier)
        while frontier:
            for neighbor in adjacency[frontier.pop()]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    frontier.append(neighbor)
        return reached

    # ------------------------------------------------------------- structure
    @functools.cached_property
    def _pos(self) -> dict[V, Point]:
        """The embedding as a mapping."""
        return dict(self.positions)

    @functools.cached_property
    def internal_vertices(self) -> frozenset[V]:
        """The vertices strictly inside the disk."""
        return frozenset(v for v, _ in self.positions) - frozenset(self.boundary)

    @functools.cached_property
    def _out_edges(self) -> dict[V, tuple[int, ...]]:
        """Outgoing edge indices per vertex."""
        table: dict[V, list[int]] = {v: [] for v, _ in self.positions}
        for index, (tail, _, _) in enumerate(self.edges):
            table[tail].append(index)
        return {v: tuple(indices) for v, indices in table.items()}

    @functools.cached_property
    def _in_edges(self) -> dict[V, tuple[int, ...]]:
        """Incoming edge indices per vertex."""
        table: dict[V, list[int]] = {v: [] for v, _ in self.positions}
        for index, (_, head, _) in enumerate(self.edges):
            table[head].append(index)
        return {v: tuple(indices) for v, indices in table.items()}

    @functools.cached_property
    def source_set(self) -> frozenset[V]:
        """The boundary sources: edge directed away, or declared isolated.

        Postnikov Definition 4.1; the source set ``I`` of the page.
        """
        sources = set(self.isolated_source_set)
        for b in self.boundary:
            if self._out_edges[b]:
                sources.add(b)
        return frozenset(sources)

    @property
    def sink_set(self) -> frozenset[V]:
        """The boundary sinks — the complement of the source set."""
        return frozenset(self.boundary) - self.source_set

    @property
    def network_type(self) -> tuple[int, int]:
        """The type ``(k, n)``: source count and boundary size."""
        return len(self.source_set), len(self.boundary)

    @functools.cached_property
    def is_acyclic(self) -> bool:
        """Whether the network has no closed directed path (Postnikov §4)."""
        state: dict[V, int] = {}

        def visit(v: V) -> bool:
            state[v] = 1
            for index in self._out_edges[v]:
                head = self.edges[index][1]
                mark = state.get(head, 0)
                if mark == 1 or (mark == 0 and not visit(head)):
                    return False
            state[v] = 2
            return True

        return all(visit(v) for v, _ in self.positions if state.get(v, 0) == 0)

    @functools.cached_property
    def is_perfectly_oriented(self) -> bool:
        """Whether every internal vertex has out-degree 1 or in-degree 1.

        Postnikov Definition 9.2 / the orientation class of Talaska's flow
        formula (Theorem 3.2) and of perfect orientations of plabic graphs
        (Postnikov Definition 11.5).
        """
        return all(
            len(self._out_edges[v]) == 1 or len(self._in_edges[v]) == 1
            for v in self.internal_vertices
        )

    # ------------------------------------------------------------ embedding
    @functools.cached_property
    def _rotations(self) -> dict[V, tuple[_Dart, ...]]:
        """Counterclockwise dart order at each vertex, frame arcs included.

        Internal rotations come from the exact embedding; at a boundary
        vertex the disk's boundary arcs enclose the at-most-one real edge,
        in the order (arc to next, arc to previous, real edge) fixed by the
        clockwise labelling (interior to the walker's right).
        """
        n = len(self.boundary)
        boundary_index = {b: t for t, b in enumerate(self.boundary)}
        rotations: dict[V, tuple[_Dart, ...]] = {}
        for v, origin in self.positions:
            if v in boundary_index:
                t = boundary_index[v]
                darts: list[_Dart] = [("f", t, 1), ("f", (t - 1) % n, -1)]
                darts.extend(
                    ("e", index, 1 if self.edges[index][0] == v else -1)
                    for index in self._out_edges[v] + self._in_edges[v]
                )
                rotations[v] = tuple(darts)
                continue
            with_directions: list[tuple[_Dart, Point]] = []
            for index in self._out_edges[v]:
                target = self._pos[self.edges[index][1]]
                delta = (target[0] - origin[0], target[1] - origin[1])
                with_directions.append((("e", index, 1), delta))
            for index in self._in_edges[v]:
                target = self._pos[self.edges[index][0]]
                delta = (target[0] - origin[0], target[1] - origin[1])
                with_directions.append((("e", index, -1), delta))
            rotations[v] = tuple(_sort_ccw(with_directions))
        return rotations

    def _dart_tail(self, dart: _Dart) -> V:
        """Return the vertex a dart leaves."""
        kind, index, sign = dart
        if kind == "f":
            offset = index if sign == 1 else (index + 1) % len(self.boundary)
            return self.boundary[offset]
        tail, head, _ = self.edges[index]
        return tail if sign == 1 else head

    def _dart_head(self, dart: _Dart) -> V:
        """Return the vertex a dart points to."""
        kind, index, sign = dart
        return self._dart_tail((kind, index, -sign))

    @functools.cached_property
    def _faces(self) -> tuple[tuple[_Dart, ...], ...]:
        """The faces of the frame-augmented embedding, as dart orbits.

        Faces are orbits of the next-dart permutation derived from the
        rotation system; the count enters the Euler check in validation and
        :attr:`face_count`.
        """
        n = len(self.boundary)
        all_darts: list[_Dart] = [
            ("e", index, sign) for index in range(len(self.edges)) for sign in (1, -1)
        ] + [("f", t, sign) for t in range(n) for sign in (1, -1)]
        remaining = set(all_darts)
        faces: list[tuple[_Dart, ...]] = []
        while remaining:
            start = remaining.pop()
            orbit = [start]
            current = start
            while True:
                kind, index, sign = current
                reverse: _Dart = (kind, index, -sign)
                rotation = self._rotations[self._dart_head(current)]
                position = rotation.index(reverse)
                successor = rotation[position - 1]
                if successor == start:
                    break
                orbit.append(successor)
                remaining.discard(successor)
                current = successor
            faces.append(tuple(orbit))
        return tuple(faces)

    @functools.cached_property
    def face_count(self) -> int:
        """The number of faces of the network in the disk, ``|F(G)|``.

        The outer face of the frame-augmented embedding (the region outside
        the disk, bounded by boundary arcs alone) is excluded, so this is
        the ``|F(G)|`` of Postnikov Theorem 12.7.
        """
        return len(self._faces) - 1

    def faces(self) -> tuple[tuple[tuple[int, int], ...], ...]:
        """Return the faces of the network in the disk, as edge darts.

        One entry per face (the region outside the disk omitted), aligned
        index-by-index with :meth:`face_weights`; each face is the
        ``(edge index, direction)`` darts of the real edges on its
        boundary, ``+1`` when the traversal follows the edge's direction.
        Boundary arcs of the disk are omitted (they carry no weight,
        Postnikov section 11).
        """
        collected: list[tuple[tuple[int, int], ...]] = []
        for orbit in self._faces:
            if all(kind == "f" for kind, _, _ in orbit):
                continue
            collected.append(
                tuple((index, sign) for kind, index, sign in orbit if kind == "e")
            )
        return tuple(collected)

    def face_weights(self) -> tuple[Fraction, ...]:
        """Return the face weight ``y_f`` of each face of the disk.

        Postnikov section 11: with the exterior boundary of the face
        oriented clockwise, ``y_f`` multiplies the weights of agreeing
        edges and divides by the disagreeing ones. The orientation is
        pinned by Lemma 11.4 — the product of the face weights to the
        right of a path equals the product of its edge weights — which the
        test suite checks on the square fixture. Aligned with
        :meth:`faces`.
        """
        weights: list[Fraction] = []
        for face in self.faces():
            y = Fraction(1)
            for index, sign in face:
                y *= 1 / self.edges[index][2] if sign == 1 else self.edges[index][2]
            weights.append(y)
        return tuple(weights)

    # ----------------------------------------------------------- measurement
    @functools.cached_property
    def _boundary_measurements(self) -> dict[tuple[V, V], Fraction]:
        """All measurements ``M_ij``, by the acyclic or the flow route."""
        result: dict[tuple[V, V], Fraction] = {}
        sources = sorted(self.source_set, key=list(self.boundary).index)
        sinks = [b for b in self.boundary if b in self.sink_set]
        if self.is_acyclic:
            for i in sources:
                paths = self._path_sums(i)
                for j in sinks:
                    result[(i, j)] = paths.get(j, Fraction(0))
            return result
        if not self.is_perfectly_oriented:
            msg = (
                "boundary measurements of a cyclic network are computable "
                "only through the flow formula, which requires a perfectly "
                "oriented network (Talaska Theorem 3.2); the signed series "
                "of Postnikov eq. (4.1) is not implemented"
            )
            raise ValueError(msg)
        denominator = self._flow_sum(self.source_set)
        source_list = list(self.source_set)
        for i in sources:
            for j in sinks:
                swapped = frozenset(set(source_list) - {i} | {j})
                result[(i, j)] = self._flow_sum(swapped) / denominator
        return result

    def _path_sums(self, source: V) -> dict[V, Fraction]:
        """Return the weighted path sums from ``source`` to every vertex.

        The acyclic boundary measurement (Postnikov section 4): a finite
        sum over directed paths, computed by memoized traversal.
        """
        memo: dict[V, dict[V, Fraction]] = {}

        def sums_from(v: V) -> dict[V, Fraction]:
            if v in memo:
                return memo[v]
            totals: dict[V, Fraction] = {v: Fraction(1)}
            for index in self._out_edges[v]:
                _, head, weight = self.edges[index]
                for target, value in sums_from(head).items():
                    totals[target] = totals.get(target, Fraction(0)) + weight * value
            memo[v] = totals
            return totals

        return sums_from(source)

    def boundary_measurement(self, i: V, j: V) -> Fraction:
        """Return the boundary measurement ``M_ij``.

        Postnikov Definition 4.4, evaluated by the acyclic path sum
        (Postnikov section 4) or, on a cyclic perfectly oriented network,
        by the flow formula applied to ``Delta_{(I - {i}) + {j}} / Delta_I``
        (Talaska Theorem 3.2 with Postnikov Definition 4.6).

        Raises:
            ValueError: If ``i`` is not a source or ``j`` not a sink, or
                the network is cyclic but not perfectly oriented.
        """
        if i not in self.source_set or j not in self.sink_set:
            msg = (
                f"measurements are indexed by a source and a sink "
                f"(Postnikov section 4); got ({i!r}, {j!r})"
            )
            raise ValueError(msg)
        return self._boundary_measurements[(i, j)]

    def boundary_measurements(self) -> dict[tuple[V, V], Fraction]:
        """Return all boundary measurements, keyed by (source, sink)."""
        return dict(self._boundary_measurements)

    def to_matrix(self) -> Columns[V]:
        """Return the boundary measurement matrix as label-to-column mapping.

        Postnikov Definition 4.6 via :func:`measurement_matrix`; the
        mapping order is the clockwise boundary order, which is exactly the
        input contract of :meth:`research.positroid.Positroid.from_matrix`.
        """
        index = {b: t + 1 for t, b in enumerate(self.boundary)}
        measurements = {
            (index[i], index[j]): value
            for (i, j), value in self._boundary_measurements.items()
        }
        columns = measurement_matrix(
            [index[i] for i in self.source_set], len(self.boundary), measurements
        )
        return {b: columns[index[b] - 1] for b in self.boundary}

    def plucker(self, subset: Iterable[V]) -> Fraction:
        """Return the Plucker coordinate ``Delta_J`` of the network's point.

        The minor of the boundary measurement matrix in the columns ``J``
        (Postnikov Definition 4.6), computed exactly.

        Raises:
            ValueError: If ``J`` is not a ``k``-subset of the boundary.
        """
        labels = set(subset)
        k = len(self.source_set)
        if len(labels) != k or not labels <= set(self.boundary):
            msg = (
                f"Plucker coordinates are indexed by {k}-subsets of the "
                f"boundary, got {sorted(map(repr, labels))}"
            )
            raise ValueError(msg)
        matrix = self.to_matrix()
        chosen = [b for b in self.boundary if b in labels]
        rows = [[matrix[b][r] for b in chosen] for r in range(k)]
        return det_q(rows)

    # ----------------------------------------------------------------- flows
    def _interior_edge_indices(self) -> list[int]:
        """Return indices of edges not incident to the boundary."""
        on_boundary = set(self.boundary)
        return [
            index
            for index, (tail, head, _) in enumerate(self.edges)
            if tail not in on_boundary and head not in on_boundary
        ]

    def _balanced(self, chosen: set[int]) -> bool:
        """Return whether in-degree equals out-degree at every interior vertex."""
        return all(
            sum(1 for e in self._in_edges[v] if e in chosen)
            == sum(1 for e in self._out_edges[v] if e in chosen)
            for v in self.internal_vertices
        )

    def flows(self, subset: Iterable[V]) -> tuple[frozenset[int], ...]:
        """Enumerate the flows from the source set to ``J``, as edge indices.

        Talaska Definition 3.1: edge sets with in-degree equal to
        out-degree at every interior vertex, using the boundary edges at
        ``I - J`` (sources) and ``J - I`` (sinks) only. Exponential in the
        number of interior edges; guarded.

        Raises:
            ValueError: If ``J`` is not a ``k``-subset of the boundary, or
                the interior edge count exceeds the enumeration guard.
        """
        labels = set(subset)
        k = len(self.source_set)
        if len(labels) != k or not labels <= set(self.boundary):
            msg = (
                f"flows are indexed by {k}-subsets of the boundary, got "
                f"{sorted(map(repr, labels))}"
            )
            raise ValueError(msg)
        forced: set[int] = set()
        for i in self.source_set - labels:
            if not self._out_edges[i]:
                return ()
            forced.add(self._out_edges[i][0])
        for j in labels - self.source_set:
            if not self._in_edges[j]:
                return ()
            forced.add(self._in_edges[j][0])
        interior = self._interior_edge_indices()
        if len(interior) > _MAX_ENUMERATION_EDGES:
            msg = (
                f"{len(interior)} interior edges exceed the enumeration "
                f"guard of {_MAX_ENUMERATION_EDGES} (page cost note: flow "
                f"enumeration is exponential)"
            )
            raise ValueError(msg)
        found: list[frozenset[int]] = []
        for r in range(len(interior) + 1):
            for combo in itertools.combinations(interior, r):
                chosen = forced | set(combo)
                if self._balanced(chosen):
                    found.append(frozenset(chosen))
        return tuple(found)

    def conservative_flows(self) -> tuple[frozenset[int], ...]:
        """Enumerate the conservative flows — no edge touches the boundary.

        Talaska Definition 3.1; these are the collections of pairwise
        vertex-disjoint oriented cycles, and the flow formula's
        denominator.
        """
        return self.flows(self.source_set)

    def _flow_sum(self, subset: Iterable[V]) -> Fraction:
        """Return the flow generating function ``sum_F wt(F)`` for ``J``."""
        return sum(
            (
                functools.reduce(
                    lambda acc, index: acc * self.edges[index][2],
                    flow,
                    Fraction(1),
                )
                for flow in self.flows(subset)
            ),
            Fraction(0),
        )

    def plucker_via_flows(self, subset: Iterable[V]) -> Fraction:
        """Return ``Delta_J`` by Talaska's flow formula.

        Talaska Theorem 3.2: the ratio of the flow generating function of
        ``J`` to the conservative-flow generating function, unsigned and
        subtraction-free. Valid for perfectly oriented planar networks
        only.

        Raises:
            ValueError: If the network is not perfectly oriented (the
                theorem's hypothesis) or ``J`` is not a ``k``-subset.
        """
        if not self.is_perfectly_oriented:
            msg = (
                "the flow formula requires a perfectly oriented network "
                "(Talaska Theorem 3.2)"
            )
            raise ValueError(msg)
        return self._flow_sum(subset) / self._flow_sum(self.source_set)

    # ------------------------------------------------------- transformations
    def gauge_transform(self, scalars: Mapping[V, Fraction]) -> PlanarNetwork[V]:
        """Return the gauge-transformed network, same boundary measurements.

        Postnikov eq. (4.2): positive ``t_v`` at internal vertices (boundary
        fixed at 1) replace each edge weight ``x_e`` by ``x_e t_u / t_v``
        for ``e = (u, v)``. Vertices missing from ``scalars`` keep 1.

        Raises:
            ValueError: If a scalar is nonpositive or attached to a
                boundary vertex.
        """
        for v, t in scalars.items():
            if v not in self.internal_vertices:
                msg = (
                    f"gauge scalars attach to internal vertices only "
                    f"(Postnikov eq. (4.2)); {v!r} is not internal"
                )
                raise ValueError(msg)
            if t <= 0:
                msg = f"gauge scalar at {v!r} must be positive, got {t}"
                raise ValueError(msg)

        def scale(v: V) -> Fraction:
            return scalars.get(v, Fraction(1))

        rescaled = tuple(
            (tail, head, weight * scale(tail) / scale(head))
            for tail, head, weight in self.edges
        )
        return PlanarNetwork(
            self.boundary, rescaled, self.positions, self.isolated_source_set
        )

    def to_positroid(self) -> Positroid[V]:
        """Return the positroid cell of the network's point.

        Postnikov Theorem 4.10 (a fixed graph maps into a single cell) and
        Corollary 5.4 (the point is totally nonnegative), realized through
        :meth:`research.positroid.Positroid.from_matrix`.
        """
        return Positroid.from_matrix(self.to_matrix())

    def to_bipartite(self) -> PlanarBipartiteNetwork[V]:
        """Return the plabic (bipartite) network of a perfect orientation.

        Postnikov's perfect-orientation correspondence: an internal vertex
        with a unique outgoing edge is black, with a unique incoming edge
        white (Definition 11.5); vertices with one of each are colored by
        propagation. Weights are inverted on black-to-white edges (Lam
        Proposition 5.3), so the bipartite point equals this network's
        point projectively.

        Raises:
            ValueError: If the network is not perfectly oriented, a
                boundary vertex is isolated (Lam's class has degree-one
                boundary vertices), or no consistent bipartition exists
                (general plabic graphs are not supported).
        """
        if not self.is_perfectly_oriented:
            msg = (
                "only a perfectly oriented network corresponds to a plabic "
                "network (Postnikov Definition 11.5)"
            )
            raise ValueError(msg)
        colors: dict[V, bool] = {}
        for v in self.internal_vertices:
            out, into = len(self._out_edges[v]), len(self._in_edges[v])
            if out == 1 and into != 1:
                colors[v] = True
            elif into == 1 and out != 1:
                colors[v] = False
        neighbors: dict[V, set[V]] = {v: set() for v, _ in self.positions}
        for tail, head, _ in self.edges:
            neighbors[tail].add(head)
            neighbors[head].add(tail)
        for b in self.boundary:
            if not neighbors[b]:
                msg = (
                    f"boundary vertex {b!r} is isolated; Lam's bipartite "
                    f"networks have degree-one boundary vertices "
                    f"(Lam section 4.1)"
                )
                raise ValueError(msg)
        every_vertex = [v for v, _ in self.positions]
        frontier = list(colors)
        while frontier or len(colors) < len(every_vertex):
            if not frontier:
                # A component with only pass-through vertices is unforced;
                # either coloring is a valid plabic structure, so seed one.
                seed = next(v for v in every_vertex if v not in colors)
                colors[seed] = True
                frontier.append(seed)
            v = frontier.pop()
            for w in neighbors[v]:
                expected = not colors[v]
                if w not in colors:
                    colors[w] = expected
                    frontier.append(w)
                elif colors[w] != expected:
                    msg = (
                        f"vertices {v!r} and {w!r} are adjacent with equal "
                        f"colors; the orientation has no bipartite plabic "
                        f"form (general plabic graphs are not supported)"
                    )
                    raise ValueError(msg)
        undirected = tuple(
            (tail, head, 1 / weight if colors[tail] else weight)
            for tail, head, weight in self.edges
        )
        return PlanarBipartiteNetwork.from_edges(
            self.boundary,
            dict(self.positions),
            undirected,
            black_vertices=[v for v, black in colors.items() if black],
        )

    @classmethod
    def from_le_diagram(
        cls,
        filling: Sequence[Sequence[int]],
        n: int,
        weights: Mapping[tuple[int, int], Fraction | int] | None = None,
    ) -> PlanarNetwork[int | tuple[int, int]]:
        """Build the Le-graph network of a Le-diagram (Postnikov §6).

        The Gamma-graph of Definition 6.3: internal vertices at the 1-boxes
        and at line crossings, leftward horizontal and downward vertical
        edges, each row line running from its boundary source to the row's
        leftmost 1-box and each column line from its topmost 1-box to the
        column's boundary sink. Rows and columns are labelled by the border
        path of the shape inside the ``k x (n - k)`` rectangle, so the
        source set is the lexicographically minimal basis ``I(lambda)``.
        The Gamma-tableau weight ``T(i, j)`` sits on the horizontal edge
        entering box ``(i, j)`` from the right and all other edges have
        weight 1 — a gauge representative (Postnikov §6: a unique gauge
        makes all vertical weights 1); the image cell does not depend on
        the choice (Theorem 4.10).

        Args:
            filling: Rows of 0/1 values, weakly decreasing lengths — the
                Le-diagram ``D`` in English notation.
            n: The number of boundary vertices; the shape must fit in the
                ``k x (n - k)`` rectangle, ``k = len(filling)``.
            weights: Optional Gamma-tableau, keyed by 1-indexed ``(i, j)``
                boxes; positive exactly on the 1-boxes (Postnikov §6).
                Defaults to all 1.

        Returns:
            The acyclic network ``N_T`` with boundary ``1..n``, whose image
            is the cell of ``D`` (Theorem 6.5).

        Raises:
            ValueError: If the shape or filling is invalid, the
                Le-condition of Definition 6.1 fails, or the tableau is
                not positive exactly on the 1-boxes.
        """
        shape = [len(row) for row in filling]
        k = len(shape)
        width = n - k
        if any(a < b for a, b in itertools.pairwise(shape)):
            msg = f"row lengths {shape} must be weakly decreasing (a shape)"
            raise ValueError(msg)
        if shape and (shape[0] > width or k > n):
            msg = (
                f"shape {shape} must fit in the {k} x {width} rectangle "
                f"(Postnikov Theorem 6.5)"
            )
            raise ValueError(msg)
        boxes = {
            (i, j): value
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
        }
        if any(value not in (0, 1) for value in boxes.values()):
            msg = "a Le-diagram is filled with 0s and 1s (Definition 6.1)"
            raise ValueError(msg)
        for (i2, j1), a in boxes.items():
            for (i1, j2), c in boxes.items():
                if (
                    i1 < i2
                    and j1 < j2
                    and a
                    and c
                    and (i2, j2) in boxes
                    and not boxes[(i2, j2)]
                ):
                    msg = (
                        f"Le-condition violated (Postnikov Definition 6.1): "
                        f"boxes ({i2}, {j1}) and ({i1}, {j2}) hold 1 but "
                        f"({i2}, {j2}) holds 0"
                    )
                    raise ValueError(msg)
        tableau = {box: Fraction(1) for box, value in boxes.items() if value}
        for box, weight in (weights or {}).items():
            value = Fraction(weight)
            if box not in boxes or (value > 0) != bool(boxes[box]):
                msg = (
                    f"a Gamma-tableau is positive exactly on the 1-boxes "
                    f"(Postnikov section 6); got T{box!r} = {weight}"
                )
                raise ValueError(msg)
            if value > 0:
                tableau[box] = value
        return PlanarNetwork._le_graph(shape, boxes, tableau, n)

    @staticmethod
    def _le_graph(
        shape: list[int],
        boxes: Mapping[tuple[int, int], int],
        tableau: Mapping[tuple[int, int], Fraction],
        n: int,
    ) -> PlanarNetwork[int | tuple[int, int]]:
        """Assemble the Gamma-graph; see :meth:`from_le_diagram`."""
        k = len(shape)
        width = n - k
        half = Fraction(1, 2)
        row_label = {i: i + width - shape[i - 1] for i in range(1, k + 1)}
        left_labels = sorted(set(range(1, n + 1)) - set(row_label.values()))
        column_label = {
            width - offset: label for offset, label in enumerate(left_labels)
        }
        ones = [box for box, value in boxes.items() if value]
        leftmost = {
            i: min(j for r, j in ones if r == i)
            for i in range(1, k + 1)
            if any(r == i for r, _ in ones)
        }
        topmost = {
            j: min(r for r, c in ones if c == j)
            for j in range(1, width + 1)
            if any(c == j for _, c in ones)
        }
        bottom_row = {
            j: max((i for i in range(1, k + 1) if shape[i - 1] >= j), default=0)
            for j in range(1, width + 1)
        }
        grid: set[tuple[int, int]] = set()
        for i in range(1, k + 1):
            for j in range(1, shape[i - 1] + 1):
                if (
                    i in leftmost
                    and leftmost[i] <= j
                    and j in topmost
                    and topmost[j] <= i
                ):
                    grid.add((i, j))
        positions: dict[int | tuple[int, int], tuple[Fraction, Fraction]] = {
            box: (Fraction(box[1]), Fraction(-box[0])) for box in grid
        }
        for i in range(1, k + 1):
            positions[row_label[i]] = (shape[i - 1] + half, Fraction(-i))
        for j in range(1, width + 1):
            positions[column_label[j]] = (Fraction(j), -bottom_row[j] - half)
        edges: list[tuple[int | tuple[int, int], int | tuple[int, int], Fraction]] = []
        for i in range(1, k + 1):
            row_vertices = sorted((b for b in grid if b[0] == i), reverse=True)
            chain: list[int | tuple[int, int]] = [row_label[i], *row_vertices]
            for tail, head in itertools.pairwise(chain):
                weight = (
                    tableau[head]
                    if isinstance(head, tuple) and boxes[head]
                    else Fraction(1)
                )
                edges.append((tail, head, weight))
        for j in range(1, width + 1):
            column_vertices = sorted(b for b in grid if b[1] == j)
            if not column_vertices:
                continue
            chain = [*column_vertices, column_label[j]]
            for tail, head in itertools.pairwise(chain):
                edges.append((tail, head, Fraction(1)))
        isolated: list[int | tuple[int, int]] = [
            row_label[i] for i in range(1, k + 1) if i not in leftmost
        ]
        boundary: tuple[int | tuple[int, int], ...] = tuple(range(1, n + 1))
        return PlanarNetwork.from_edges(
            boundary,
            positions,
            edges,
            isolated_sources=isolated,
        )

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per vertex or edge.

        Columns: ``kind`` (``vertex``/``edge``); vertex rows carry
        ``vertex`` (repr for non-string labels stays as-is for hashables
        serializable by pandas), ``role`` (``source``/``sink``/
        ``internal``) and exact coordinates ``x``, ``y`` as fraction
        strings; edge rows carry ``tail``, ``head`` and ``weight`` as a
        fraction string. Vertex rows list the boundary first, in clockwise
        order, so the boundary order survives a records-oriented JSON round
        trip through ``experiments.io.write_result``.

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        ordered = [
            *self.boundary,
            *(v for v, _ in self.positions if v in self.internal_vertices),
        ]
        rows: list[dict[str, object]] = []
        for v in ordered:
            if v in self.internal_vertices:
                role = "internal"
            else:
                role = "source" if v in self.source_set else "sink"
            x, y = self._pos[v]
            rows.append(
                {
                    "kind": "vertex",
                    "vertex": v,
                    "role": role,
                    "x": str(x),
                    "y": str(y),
                    "tail": None,
                    "head": None,
                    "weight": None,
                }
            )
        rows.extend(
            {
                "kind": "edge",
                "vertex": None,
                "role": None,
                "x": None,
                "y": None,
                "tail": tail,
                "head": head,
                "weight": str(weight),
            }
            for tail, head, weight in self.edges
        )
        return pd.DataFrame(rows)

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> PlanarNetwork[Hashable]:
        """Rebuild a network from a frame produced by :meth:`to_dataframe`.

        Re-validates through :meth:`from_edges`.

        Args:
            df: Frame with the ``kind``/``vertex``/``role``/``x``/``y``/
                ``tail``/``head``/``weight`` columns.

        Returns:
            The decoded network.

        Raises:
            ValueError: If required columns are missing or the decoded
                data fails network validation.
        """
        required = {"kind", "vertex", "role", "x", "y", "tail", "head", "weight"}
        missing = required - set(df.columns)
        if missing:
            msg = f"dataframe is missing required columns {sorted(missing)}"
            raise ValueError(msg)
        vertex_rows = df[df["kind"] == "vertex"]
        edge_rows = df[df["kind"] == "edge"]
        boundary: list[Hashable] = []
        positions: dict[Hashable, tuple[Fraction, Fraction]] = {}
        isolated_sources: list[Hashable] = []
        incident: set[Hashable] = set()
        for tail, head in zip(edge_rows["tail"], edge_rows["head"], strict=True):
            incident.add(tail)
            incident.add(head)
        for vertex, role, x, y in zip(
            vertex_rows["vertex"],
            vertex_rows["role"],
            vertex_rows["x"],
            vertex_rows["y"],
            strict=True,
        ):
            positions[vertex] = (Fraction(str(x)), Fraction(str(y)))
            if role != "internal":
                boundary.append(vertex)
            if role == "source" and vertex not in incident:
                isolated_sources.append(vertex)
        edges = [
            (tail, head, Fraction(str(weight)))
            for tail, head, weight in zip(
                edge_rows["tail"], edge_rows["head"], edge_rows["weight"], strict=True
            )
        ]
        return PlanarNetwork.from_edges(
            boundary, positions, edges, isolated_sources=isolated_sources
        )

    # ------------------------------------------------------- visualization
    def plot_network(self, ax: Axes | None = None) -> Axes:
        """Draw the network at its stored embedding onto the axes.

        Boundary vertices are labeled with their clockwise labels, edges
        drawn as arrows annotated with their weights. Draws onto ``ax`` or
        a fresh figure; never calls ``show``.
        """
        ax = ensure_axes(ax)
        points = [
            (float(self._pos[v][0]), float(self._pos[v][1])) for v, _ in self.positions
        ]
        scatter_labeled(ax, points, [str(v) for v, _ in self.positions])
        for tail, head, weight in self.edges:
            (x0, y0), (x1, y1) = self._pos[tail], self._pos[head]
            ax.annotate(
                "",
                xy=(float(x1), float(y1)),
                xytext=(float(x0), float(y0)),
                arrowprops={"arrowstyle": "-|>", "color": "tab:gray"},
            )
            ax.annotate(
                str(weight),
                ((float(x0) + float(x1)) / 2, (float(y0) + float(y1)) / 2),
                fontsize=7,
                ha="center",
            )
        ax.set_aspect("equal")
        ax.set_axis_off()
        return ax


# --------------------------------------------------------------------------- #
# Planar bipartite (dimer) networks — Lam's matching formulation
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PlanarBipartiteNetwork[V: Hashable]:
    """A weighted planar bipartite network in a disk (Lam section 4.1).

    Boundary vertices sit on the disk's boundary circle labelled clockwise,
    each of degree one; every vertex is black or white and every edge joins
    the two colors; interior vertices are never isolated (Lam's standing
    assumption, without which the matching formulation is empty). The
    boundary measurement of a k-subset is the dimer partition function over
    almost perfect matchings with that boundary subset.

    Build through :meth:`from_edges` (which validates) or the named example
    constructors; calling the dataclass constructor directly skips
    validation.
    """

    boundary: tuple[V, ...]
    edges: tuple[tuple[V, V, Fraction], ...]
    positions: tuple[tuple[V, Point], ...]
    black_vertices: frozenset[V]

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        return (
            f"PlanarBipartiteNetwork(n={len(self.boundary)}, "
            f"edges={len(self.edges)}, "
            f"internal={len(self.internal_vertices)})"
        )

    # ---------------------------------------------------------- construction
    @classmethod
    def from_edges(
        cls,
        boundary: Sequence[V],
        positions: Mapping[V, tuple[Fraction | int, Fraction | int]],
        edges: Iterable[tuple[V, V, Fraction | int]],
        *,
        black_vertices: Iterable[V],
    ) -> PlanarBipartiteNetwork[V]:
        """Build and validate a planar bipartite network.

        Args:
            boundary: Boundary vertex labels in clockwise order.
            positions: Exact coordinates for every vertex; the embedding.
            edges: ``(u, v, weight)`` triples with positive weights;
                undirected.
            black_vertices: The black color class; all other vertices are
                white.

        Returns:
            The validated network.

        Raises:
            ValueError: If a weight is nonpositive, an edge joins two
                vertices of one color (the network must be bipartite, Lam
                section 4.1), a boundary vertex does not have degree one, an
                interior vertex is isolated (Lam's standing assumption), or
                the embedding checks shared with :class:`PlanarNetwork`
                fail.
        """
        placed = {v: (Fraction(p[0]), Fraction(p[1])) for v, p in positions.items()}
        ordered = [*boundary, *(v for v in placed if v not in set(boundary))]
        network = cls(
            tuple(boundary),
            tuple((u, v, Fraction(weight)) for u, v, weight in edges),
            tuple((v, placed[v]) for v in ordered),
            frozenset(black_vertices),
        )
        network._validate()
        return network

    def _validate(self) -> None:
        """Check the definition; raise ``ValueError`` naming the condition."""
        black = self.black_vertices
        for u, v, _ in self.edges:
            if (u in black) == (v in black):
                color = "black" if u in black else "white"
                msg = (
                    f"edge ({u!r}, {v!r}) joins two {color} vertices; the "
                    f"network must be bipartite (Lam section 4.1)"
                )
                raise ValueError(msg)
        degree = {v: 0 for v, _ in self.positions}
        for u, v, _ in self.edges:
            degree[u] += 1
            degree[v] += 1
        for b in self.boundary:
            if degree[b] != 1:
                msg = (
                    f"boundary vertex {b!r} has degree {degree[b]}; Lam's "
                    f"networks have degree-one boundary vertices "
                    f"(Lam section 4.1)"
                )
                raise ValueError(msg)
        for v in self.internal_vertices:
            if degree[v] == 0:
                msg = (
                    f"interior vertex {v!r} is isolated; almost perfect "
                    f"matchings could not exist (Lam's standing assumption, "
                    f"section 4.1)"
                )
                raise ValueError(msg)
        # The geometric and disk checks coincide with the directed class;
        # an arbitrary orientation of the edges does not affect them.
        directed = PlanarNetwork(self.boundary, self.edges, self.positions)
        directed._validate()

    # ------------------------------------------------------------- structure
    @functools.cached_property
    def _pos(self) -> dict[V, Point]:
        """The embedding as a mapping."""
        return dict(self.positions)

    @functools.cached_property
    def internal_vertices(self) -> frozenset[V]:
        """The vertices strictly inside the disk."""
        return frozenset(v for v, _ in self.positions) - frozenset(self.boundary)

    @functools.cached_property
    def _incident(self) -> dict[V, tuple[int, ...]]:
        """Incident edge indices per vertex."""
        table: dict[V, list[int]] = {v: [] for v, _ in self.positions}
        for index, (u, v, _) in enumerate(self.edges):
            table[u].append(index)
            table[v].append(index)
        return {v: tuple(indices) for v, indices in table.items()}

    def _other(self, index: int, v: V) -> V:
        """Return the endpoint of edge ``index`` other than ``v``."""
        u, w, _ = self.edges[index]
        return w if u == v else u

    # -------------------------------------------------------------- matchings
    def almost_perfect_matchings(self) -> tuple[frozenset[int], ...]:
        """Enumerate the almost perfect matchings, as edge-index sets.

        Lam section 4.1: edge subsets using each interior vertex exactly
        once; boundary vertices may or may not be used. Exponential;
        guarded.

        Raises:
            ValueError: If the edge count exceeds the enumeration guard.
        """
        if len(self.edges) > _MAX_ENUMERATION_EDGES:
            msg = (
                f"{len(self.edges)} edges exceed the enumeration guard of "
                f"{_MAX_ENUMERATION_EDGES} (page cost note: matching "
                f"enumeration is exponential)"
            )
            raise ValueError(msg)
        interior = sorted(self.internal_vertices, key=str)
        found: list[frozenset[int]] = []

        def extend(position: int, used: set[V], chosen: set[int]) -> None:
            if position == len(interior):
                found.append(frozenset(chosen))
                return
            v = interior[position]
            if v in used:
                extend(position + 1, used, chosen)
                return
            for index in self._incident[v]:
                other = self._other(index, v)
                if other in used:
                    continue
                extend(position + 1, used | {v, other}, chosen | {index})

        extend(0, set(), set())
        return tuple(found)

    def boundary_subset(self, matching: Iterable[int]) -> frozenset[V]:
        """Return ``I(Pi)`` — the boundary subset of a matching.

        Lam section 4.1: the black boundary vertices used by the matching
        together with the white boundary vertices *not* used by it (the
        "not used" clause is the page's, not a typo).
        """
        used: set[V] = set()
        for index in matching:
            u, v, _ = self.edges[index]
            used.add(u)
            used.add(v)
        return frozenset(
            b for b in self.boundary if (b in self.black_vertices) == (b in used)
        )

    @functools.cached_property
    def _pluckers(self) -> dict[frozenset[V], Fraction]:
        """The dimer partition functions, keyed by boundary subset."""
        totals: dict[frozenset[V], Fraction] = {}
        for matching in self.almost_perfect_matchings():
            weight = Fraction(1)
            for index in matching:
                weight *= self.edges[index][2]
            subset = self.boundary_subset(matching)
            totals[subset] = totals.get(subset, Fraction(0)) + weight
        if not totals:
            msg = (
                "the network has no almost perfect matching; Lam's standing "
                "assumption (section 4.1) fails and the boundary "
                "measurements are undefined"
            )
            raise ValueError(msg)
        return totals

    def pluckers(self) -> dict[frozenset[V], Fraction]:
        """Return all nonzero Plucker coordinates of the network's point.

        Lam section 4.1: ``Delta_I(N) = sum over matchings with
        I(Pi) = I of wt(Pi)`` — for Lam this is the definition of the
        boundary measurement, and Theorem 4.1 identifies the numbers as
        Plucker coordinates of a point of ``Gr(k, n)``.
        """
        return dict(self._pluckers)

    def plucker(self, subset: Iterable[V]) -> Fraction:
        """Return the dimer partition function ``Delta_I`` (Lam §4.1)."""
        return self._pluckers.get(frozenset(subset), Fraction(0))

    def to_positroid(self) -> Positroid[V]:
        """Return the positroid of the network's point.

        The bases are the boundary subsets with a positive partition
        function; the point is totally nonnegative by Lam Theorem 5.2
        (positive weights, perfectly orientable via Proposition 5.1), and
        :meth:`research.positroid.Positroid.from_bases` validates the
        positroid property.
        """
        return Positroid.from_bases(self.boundary, self._pluckers.keys())

    # ------------------------------------------------------- transformations
    def to_perfect_orientation(
        self, matching: Iterable[int] | None = None
    ) -> PlanarNetwork[V]:
        """Return the perfectly oriented directed network of a matching.

        Lam Proposition 5.1: orient each edge white-to-black when it is not
        in the matching and black-to-white when it is; the result is a
        perfect orientation whose source set is ``I(Pi)``. Weights are
        inverted on the black-to-white edges (Lam Proposition 5.3), so the
        directed network represents the same point: its Plucker coordinates
        are this network's divided by ``Delta_{I(Pi)}``.

        Args:
            matching: Edge indices of an almost perfect matching; defaults
                to the first enumerated one.

        Returns:
            The perfectly oriented planar directed network.
        """
        if matching is None:
            chosen = set(self.almost_perfect_matchings()[0])
        else:
            chosen = set(matching)
        oriented: list[tuple[V, V, Fraction]] = []
        for index, (u, v, weight) in enumerate(self.edges):
            white, black = (v, u) if u in self.black_vertices else (u, v)
            if index in chosen:
                oriented.append((black, white, 1 / weight))
            else:
                oriented.append((white, black, weight))
        return PlanarNetwork.from_edges(self.boundary, dict(self.positions), oriented)

    # ------------------------------------------------------------------ trips
    @functools.cached_property
    def _ccw_edges(self) -> dict[V, tuple[int, ...]]:
        """Incident edge indices in counterclockwise order at each vertex."""
        rotations: dict[V, tuple[int, ...]] = {}
        for v, origin in self.positions:
            with_directions: list[tuple[int, Point]] = []
            for index in self._incident[v]:
                target = self._pos[self._other(index, v)]
                delta = (target[0] - origin[0], target[1] - origin[1])
                with_directions.append((index, delta))
            rotations[v] = tuple(_sort_ccw(with_directions))
        return rotations

    def _next_leg(self, v: V, arrived_by: int) -> int:
        """Return the edge continuing a trip that reached ``v``.

        Lam section 7.1: turn maximally right at a black vertex and
        maximally left at a white one — the counterclockwise successor
        (respectively predecessor) of the arrival edge in the vertex
        rotation, as fixed by the clockwise boundary convention.
        """
        rotation = self._ccw_edges[v]
        position = rotation.index(arrived_by)
        step = 1 if v in self.black_vertices else -1
        return rotation[(position + step) % len(rotation)]

    def trips(self) -> tuple[tuple[tuple[V, int], ...], ...]:
        """Return the trip decomposition, boundary trips first.

        Lam section 7.1: from each boundary vertex, follow the unique edge
        and apply the turning rules until the trip exits at the boundary;
        any darts not covered belong to closed cyclic trips, traced last.
        Each trip is a sequence of ``(tail vertex, edge index)`` darts;
        together they cover each edge once in each direction.
        """
        on_boundary = set(self.boundary)
        trips: list[tuple[tuple[V, int], ...]] = []
        covered: set[tuple[V, int]] = set()

        def trace(start: V, first_edge: int) -> tuple[tuple[V, int], ...]:
            darts: list[tuple[V, int]] = []
            v, edge = start, first_edge
            while True:
                darts.append((v, edge))
                covered.add((v, edge))
                v = self._other(edge, v)
                if v in on_boundary:
                    return tuple(darts)
                edge = self._next_leg(v, edge)
                if (v, edge) == darts[0]:
                    return tuple(darts)

        trips.extend(trace(b, self._incident[b][0]) for b in self.boundary)
        trips.extend(
            trace(v, edge)
            for v, _ in self.positions
            if v not in on_boundary
            for edge in self._incident[v]
            if (v, edge) not in covered
        )
        return tuple(trips)

    def trip_permutation(self) -> tuple[int, ...]:
        """Return the trip permutation, on boundary positions ``1..n``.

        Lam section 7.1: ``pi_G(i) = j`` when the trip entering at boundary
        vertex ``i`` leaves at ``j``. Matches the decorated permutation of
        the image cell for reduced graphs (Lam Theorem 7.12(3); the page's
        square-graph fixture certifies ``(3, 4, 1, 2)``).
        """
        position = {b: t for t, b in enumerate(self.boundary, start=1)}
        all_trips = self.trips()
        targets: list[int] = []
        for b in self.boundary:
            darts = all_trips[position[b] - 1]
            tail, edge = darts[-1]
            targets.append(position[self._other(edge, tail)])
        return tuple(targets)

    def is_reduced(self) -> bool:
        """Return whether the graph is reduced, by Lam's trip conditions.

        Lam section 7.1, stated for leafless planar bipartite graphs: no
        trip is a closed cycle, no trip uses an edge twice (except at a
        boundary leaf — a lollipop head), and no two trips share two edges
        appearing in the same order in both.

        The middle clause does fire, and does so on its own. A trip
        cannot repeat a *dart* (the trip map is a bijection on darts), but
        it can leave along an edge, travel the rest of a cycle, and come
        back out along that same edge in the other direction — not the
        immediate reversal :meth:`_next_leg` performs at a leaf, so the
        boundary-leaf exception does not cover it. The square 4-cycle with
        its two boundary legs on adjacent corners is the smallest
        instance; two stacked squares on four legs is the smallest on
        which this clause alone decides, with the other two silent.

        Raises:
            ValueError: If the graph has an interior leaf hanging off an
                internal vertex — outside the scope Lam states the
                condition for.
        """
        for v in self.internal_vertices:
            incident = self._incident[v]
            if len(incident) == 1 and self._other(incident[0], v) not in set(
                self.boundary
            ):
                msg = (
                    f"interior leaf {v!r} hangs off an internal vertex; "
                    f"Lam's reducedness condition (section 7.1) is stated "
                    f"for leafless planar bipartite graphs"
                )
                raise ValueError(msg)
        all_trips = self.trips()
        boundary_trips = all_trips[: len(self.boundary)]
        if len(all_trips) > len(self.boundary):
            return False
        for darts in boundary_trips:
            edge_sequence = [edge for _, edge in darts]
            for edge in set(edge_sequence):
                occurrences = [
                    index for index, other in enumerate(edge_sequence) if other == edge
                ]
                if len(occurrences) == 1:
                    continue
                hinge = self._other(edge, darts[occurrences[0]][0])
                first_pair = occurrences[:2]
                bounce = (
                    occurrences == [first_pair[0], first_pair[0] + 1]
                    and len(self._incident[hinge]) == 1
                )
                if not bounce:
                    return False
        for first, second in itertools.combinations(boundary_trips, 2):
            first_edges = [edge for _, edge in first]
            second_edges = [edge for _, edge in second]
            shared = set(first_edges) & set(second_edges)
            for e, f in itertools.combinations(shared, 2):
                if (first_edges.index(e) < first_edges.index(f)) == (
                    second_edges.index(e) < second_edges.index(f)
                ):
                    return False
        return True

    # ------------------------------------------------------------ the moves
    def square_move(
        self, corners: tuple[V, V, V, V]
    ) -> PlanarBipartiteNetwork[V | tuple[V, str]]:
        """Apply the square move (urban renewal) at a quadrilateral face.

        Lam section 4.5, eq. (19), the gauge-fixed spider move: the four
        corner weights ``a, b, c, d`` (cyclic, so ``a`` and ``c`` are
        opposite) become ``a/(ac + bd)``, ..., with every other edge at the
        spider gauge-fixed to weight 1. Following urban renewal, the
        corners are recolored, the moved weight lands on the opposite side
        of the square — the labeling of the primed weights is fixed by the
        invariance oracle (Lam Proposition 4.8) on Example 4.3 — and each
        corner's outside edges transfer to a new same-colored vertex joined
        to it by a weight-1 edge.

        Args:
            corners: The four internal corner vertices, in cyclic order
                around the face.

        Returns:
            The moved network, representing the same point (Lam
            Proposition 4.8 / Postnikov Lemma 12.2).

        Raises:
            ValueError: If the corners are not an alternately colored
                quadrilateral face of internal vertices, or some outside
                edge at a corner does not have weight 1.
        """
        distinct = set(corners)
        if len(distinct) != _SQUARE_SIDES or not distinct <= self.internal_vertices:
            msg = f"corners {corners!r} must be four distinct internal vertices"
            raise ValueError(msg)
        colors = [v in self.black_vertices for v in corners]
        if colors[0] == colors[1] or colors != [colors[0], colors[1]] * 2:
            msg = (
                f"corners {corners!r} must alternate colors around the "
                f"square (Lam section 4.5)"
            )
            raise ValueError(msg)
        square_edges, outside = self._spider_edges(corners)
        weights = [square_edges[offset][1] for offset in range(_SQUARE_SIDES)]
        denominator = weights[0] * weights[2] + weights[1] * weights[3]
        center_x = sum((self._pos[v][0] for v in corners), Fraction(0)) / 4
        center_y = sum((self._pos[v][1] for v in corners), Fraction(0)) / 4
        moved: dict[V | tuple[V, str], Point] = dict(self.positions)
        new_edges: list[tuple[V | tuple[V, str], V | tuple[V, str], Fraction]] = []
        square_indices = {index for index, _ in square_edges.values()}
        carried = {index for indices in outside.values() for index in indices}
        for index, (u, v, weight) in enumerate(self.edges):
            if index in square_indices or index in carried:
                continue
            new_edges.append((u, v, weight))
        for offset in range(_SQUARE_SIDES):
            u, v = corners[offset], corners[(offset + 1) % _SQUARE_SIDES]
            side = weights[(offset + 2) % _SQUARE_SIDES]
            new_edges.append((u, v, side / denominator))
        black: set[V | tuple[V, str]] = set(self.black_vertices)
        for v in corners:
            carrier: tuple[V, str] = (v, "leg")
            x, y = self._pos[v]
            moved[carrier] = (x + (x - center_x) / 2, y + (y - center_y) / 2)
            new_edges.append((v, carrier, Fraction(1)))
            for index in outside[v]:
                a, b_vertex, weight = self.edges[index]
                other = b_vertex if a == v else a
                new_edges.append((other, carrier, weight))
            if v in self.black_vertices:
                black.discard(v)
                black.add(carrier)
            else:
                black.add(v)
        return PlanarBipartiteNetwork.from_edges(
            self.boundary, moved, new_edges, black_vertices=black
        )

    def _spider_edges(
        self, corners: tuple[V, V, V, V]
    ) -> tuple[dict[int, tuple[int, Fraction]], dict[V, list[int]]]:
        """Split the edges at the spider into square sides and outside legs.

        Enforces the gauge-fixed hypothesis of Lam eq. (19): every
        non-square edge at a corner has weight 1.

        Raises:
            ValueError: If a side edge is missing or an outside edge is
                not gauge-fixed.
        """
        distinct = set(corners)
        pair_index: dict[frozenset[V], int] = {}
        for offset in range(_SQUARE_SIDES):
            u, v = corners[offset], corners[(offset + 1) % _SQUARE_SIDES]
            pair_index[frozenset((u, v))] = offset
        square_edges: dict[int, tuple[int, Fraction]] = {}
        outside: dict[V, list[int]] = {v: [] for v in corners}
        for index, (u, v, weight) in enumerate(self.edges):
            key = frozenset((u, v))
            if key in pair_index:
                square_edges[pair_index[key]] = (index, weight)
                continue
            for endpoint in (u, v):
                if endpoint in distinct:
                    if weight != 1:
                        msg = (
                            f"outside edge ({u!r}, {v!r}) has weight "
                            f"{weight}; the gauge-fixed square move needs "
                            f"the spider's other edges at weight 1 "
                            f"(Lam section 4.5, eq. (19))"
                        )
                        raise ValueError(msg)
                    outside[endpoint].append(index)
        if len(square_edges) != _SQUARE_SIDES:
            msg = (
                f"corners {corners!r} must bound a quadrilateral face with "
                f"all four side edges present"
            )
            raise ValueError(msg)
        return square_edges, outside

    def add_bridge(
        self,
        i: int,
        weight: Fraction | int,
        *,
        variant: Literal["x", "y"] = "x",
    ) -> PlanarBipartiteNetwork[V | tuple[V, str]]:
        """Add a bridge between boundary positions ``i`` and ``i + 1``.

        Lam Lemma 7.6: a bridge of weight ``a``, white at ``i`` and black
        at ``i + 1``, multiplies the network's point by the Chevalley
        generator ``x_i(a)``; the opposite colouring (``variant="y"``)
        gives ``y_i(a)``. The bridge hangs between the two boundary legs,
        with weight-1 buffer vertices inserted where the new endpoints
        would clash with the leg colors.

        Args:
            i: Boundary position in ``1..n-1`` (the wrap-around bridge is
                not implemented).
            weight: The bridge weight ``a``, positive.
            variant: ``"x"`` for white at ``i``, ``"y"`` for black at
                ``i``.

        Returns:
            The bridged network.

        Raises:
            ValueError: If ``i`` is out of range.
        """
        n = len(self.boundary)
        if not 1 <= i < n:
            msg = f"bridge position {i} must satisfy 1 <= i < n = {n}"
            raise ValueError(msg)
        a = Fraction(weight)
        moved: dict[V | tuple[V, str], Point] = dict(self.positions)
        new_edges: list[tuple[V | tuple[V, str], V | tuple[V, str], Fraction]] = []
        black: set[V | tuple[V, str]] = set(self.black_vertices)
        posts: list[tuple[V, str]] = []
        for position, wants_black in ((i, variant == "y"), (i + 1, variant == "x")):
            b = self.boundary[position - 1]
            leg = self._incident[b][0]
            neighbor = self._other(leg, b)
            leg_weight = self.edges[leg][2]
            post: tuple[V, str] = (b, "bridge")
            posts.append(post)
            bx, by = self._pos[b]
            nx, ny = self._pos[neighbor]
            moved[post] = (bx + (nx - bx) / 4, by + (ny - by) / 4)
            if wants_black:
                black.add(post)
            chain: list[V | tuple[V, str]] = [b, post]
            if (b in black) == wants_black:
                buffer_near: tuple[V, str] = (b, "buffer near")
                moved[buffer_near] = (bx + (nx - bx) / 8, by + (ny - by) / 8)
                if not wants_black:
                    black.add(buffer_near)
                chain = [b, buffer_near, post]
            if (neighbor in black) == wants_black:
                buffer_far: tuple[V, str] = (b, "buffer far")
                moved[buffer_far] = (bx + (nx - bx) / 2, by + (ny - by) / 2)
                if not wants_black:
                    black.add(buffer_far)
                chain.append(buffer_far)
            chain.append(neighbor)
            weights = [leg_weight] + [Fraction(1)] * (len(chain) - 2)
            for (tail, head), w in zip(itertools.pairwise(chain), weights, strict=True):
                new_edges.append((tail, head, w))
        skipped = {
            self._incident[self.boundary[i - 1]][0],
            self._incident[self.boundary[i]][0],
        }
        new_edges.extend(
            (u, v, w)
            for index, (u, v, w) in enumerate(self.edges)
            if index not in skipped
        )
        new_edges.append((posts[0], posts[1], a))
        return PlanarBipartiteNetwork.from_edges(
            self.boundary, moved, new_edges, black_vertices=black
        )

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per vertex or edge.

        The encoding of :meth:`PlanarNetwork.to_dataframe` with ``role``
        taking ``boundary``/``internal`` and an additional ``color`` column
        on vertex rows; boundary rows come first, in clockwise order.

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        ordered = [
            *self.boundary,
            *(v for v, _ in self.positions if v in self.internal_vertices),
        ]
        rows: list[dict[str, object]] = []
        for v in ordered:
            x, y = self._pos[v]
            rows.append(
                {
                    "kind": "vertex",
                    "vertex": v,
                    "role": "internal" if v in self.internal_vertices else "boundary",
                    "color": "black" if v in self.black_vertices else "white",
                    "x": str(x),
                    "y": str(y),
                    "tail": None,
                    "head": None,
                    "weight": None,
                }
            )
        rows.extend(
            {
                "kind": "edge",
                "vertex": None,
                "role": None,
                "color": None,
                "x": None,
                "y": None,
                "tail": u,
                "head": v,
                "weight": str(weight),
            }
            for u, v, weight in self.edges
        )
        return pd.DataFrame(rows)

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> PlanarBipartiteNetwork[Hashable]:
        """Rebuild a network from a frame produced by :meth:`to_dataframe`.

        Re-validates through :meth:`from_edges`.

        Args:
            df: Frame with the columns written by :meth:`to_dataframe`.

        Returns:
            The decoded network.

        Raises:
            ValueError: If required columns are missing or the decoded
                data fails network validation.
        """
        required = {
            "kind",
            "vertex",
            "role",
            "color",
            "x",
            "y",
            "tail",
            "head",
            "weight",
        }
        missing = required - set(df.columns)
        if missing:
            msg = f"dataframe is missing required columns {sorted(missing)}"
            raise ValueError(msg)
        vertex_rows = df[df["kind"] == "vertex"]
        edge_rows = df[df["kind"] == "edge"]
        boundary: list[Hashable] = []
        positions: dict[Hashable, tuple[Fraction, Fraction]] = {}
        black: list[Hashable] = []
        for vertex, role, color, x, y in zip(
            vertex_rows["vertex"],
            vertex_rows["role"],
            vertex_rows["color"],
            vertex_rows["x"],
            vertex_rows["y"],
            strict=True,
        ):
            positions[vertex] = (Fraction(str(x)), Fraction(str(y)))
            if role == "boundary":
                boundary.append(vertex)
            if color == "black":
                black.append(vertex)
        edges = [
            (tail, head, Fraction(str(weight)))
            for tail, head, weight in zip(
                edge_rows["tail"], edge_rows["head"], edge_rows["weight"], strict=True
            )
        ]
        return PlanarBipartiteNetwork.from_edges(
            boundary, positions, edges, black_vertices=black
        )

    # ------------------------------------------------------- visualization
    def plot_network(self, ax: Axes | None = None) -> Axes:
        """Draw the network with filled black and open white vertices.

        Edges are annotated with their weights; boundary vertices carry
        their labels. Draws onto ``ax`` or a fresh figure; never calls
        ``show``.
        """
        ax = ensure_axes(ax)
        for u, v, weight in self.edges:
            (x0, y0), (x1, y1) = self._pos[u], self._pos[v]
            ax.plot(
                [float(x0), float(x1)],
                [float(y0), float(y1)],
                color="tab:gray",
                zorder=1,
            )
            ax.annotate(
                str(weight),
                ((float(x0) + float(x1)) / 2, (float(y0) + float(y1)) / 2),
                fontsize=7,
                ha="center",
            )
        for v, point in self.positions:
            filled = v in self.black_vertices
            ax.scatter(
                [float(point[0])],
                [float(point[1])],
                facecolors="black" if filled else "white",
                edgecolors="black",
                zorder=2,
            )
            ax.annotate(
                str(v),
                (float(point[0]), float(point[1])),
                textcoords="offset points",
                xytext=(0.0, 6.0),
                ha="center",
                fontsize=8,
            )
        ax.set_aspect("equal")
        ax.set_axis_off()
        return ax


# --------------------------------------------------------------------------- #
# Canonical example constructors — the page's test fixtures
# --------------------------------------------------------------------------- #
def square_network(
    a: Fraction | int = 1,
    b: Fraction | int = 1,
    c: Fraction | int = 1,
    d: Fraction | int = 1,
) -> PlanarBipartiteNetwork[int | str]:
    """Return the ``Gr(2, 4)`` square graph — the workhorse fixture.

    Lam Example 4.3: boundary vertices 1 (top), 2 (right), 3 (bottom),
    4 (left); an interior square with white vertices left and right, black
    top and bottom; edge weights ``a`` (left-top), ``b`` (top-right),
    ``c`` (right-bottom), ``d`` (bottom-left), legs of weight 1. Its
    Plucker coordinates are ``Delta_12 = a``, ``Delta_13 = ac + bd``,
    ``Delta_14 = b``, ``Delta_23 = d``, ``Delta_24 = 1``,
    ``Delta_34 = c``, its image the top cell of ``Gr(2, 4)``, and its trip
    permutation ``(3, 4, 1, 2)`` (Lam section 7.1).
    """
    return PlanarBipartiteNetwork.from_edges(
        (1, 2, 3, 4),
        {
            1: (0, 2),
            2: (2, 0),
            3: (0, -2),
            4: (-2, 0),
            "T": (0, 1),
            "R": (1, 0),
            "B": (0, -1),
            "L": (-1, 0),
        },
        [
            (1, "T", 1),
            (2, "R", 1),
            (3, "B", 1),
            (4, "L", 1),
            ("L", "T", a),
            ("T", "R", b),
            ("R", "B", c),
            ("B", "L", d),
        ],
        black_vertices=["T", "B", 2, 4],
    )


def lollipop_network() -> PlanarBipartiteNetwork[int | str]:
    """Return the lollipop graph on four boundary vertices.

    Lam Example 4.2: the smallest possible graph, all boundary vertices of
    degree one, each attached to its own interior lollipop head. It has a
    single almost perfect matching (all four edges) with
    ``I(Pi) = {3, 4}``, so the point is the torus-fixed point
    ``span(e_3, e_4)`` of ``Gr(2, 4)`` — certifying zero-dimensional cells
    and the "not used" clause in ``I(Pi)``.
    """
    return PlanarBipartiteNetwork.from_edges(
        (1, 2, 3, 4),
        {
            1: (0, 2),
            2: (2, 0),
            3: (0, -2),
            4: (-2, 0),
            "p1": (0, 1),
            "p2": (1, 0),
            "p3": (0, -1),
            "p4": (-1, 0),
        },
        [(1, "p1", 1), (2, "p2", 1), (3, "p3", 1), (4, "p4", 1)],
        black_vertices=[3, 4, "p1", "p2"],
    )


def geometric_series_network(
    x: Fraction | int = 1,
    y: Fraction | int = 1,
    z: Fraction | int = 1,
    t: Fraction | int = 1,
) -> PlanarNetwork[int | str]:
    """Return Postnikov's geometric-series network.

    Postnikov Example 4.5: two boundary vertices and a cycle, with
    ``M_12 = xyt / (1 + yz)`` — the smallest witness that a cyclic network
    needs the subtraction-free rational form; at all weights 1 the value
    is 1/2. The return edge of the cycle is subdivided by a weight-1
    vertex (a representation restriction that changes no path or flow
    weight), and the network is perfectly oriented, so the flow formula
    computes the measurement.
    """
    return PlanarNetwork.from_edges(
        (1, 2),
        {1: (-2, 0), 2: (2, 0), "u": (-1, 0), "v": (1, 0), "w": (0, -1)},
        [
            (1, "u", x),
            ("u", "v", y),
            ("v", 2, t),
            ("v", "w", z),
            ("w", "u", 1),
        ],
    )


def acyclic_baseline_network(
    a: Fraction | int = 1,
    b: Fraction | int = 1,
    c: Fraction | int = 1,
) -> PlanarNetwork[int | str]:
    """Return Lam's acyclic three-wire baseline network.

    Lam section 2.3: sources 1, 2, 3 and sinks 1', 2', 3' joined by three
    rightward wires, weights ``a``, ``b``, ``c`` on three diagonal edges
    and 1 elsewhere, with path matrix ``M(N) = [[1 + ac, a, 0], [c, 1, 0],
    [bc, b, 1]]`` (rows sources, columns sinks) — the Lindstrom side of
    the story, where no denominator appears at all.
    """
    return PlanarNetwork.from_edges(
        (1, "1'", "2'", "3'", 3, 2),
        {
            1: (0, 0),
            2: (0, -2),
            3: (0, -4),
            "1'": (6, 0),
            "2'": (6, -2),
            "3'": (6, -4),
            "A1": (1, 0),
            "C1": (3, 0),
            "A2": (1, -2),
            "B2": (2, -2),
            "C2": (3, -2),
            "B3": (2, -4),
        },
        [
            (1, "A1", 1),
            ("A1", "C1", 1),
            ("C1", "1'", 1),
            (2, "A2", 1),
            ("A2", "B2", 1),
            ("B2", "C2", 1),
            ("C2", "2'", 1),
            (3, "B3", 1),
            ("B3", "3'", 1),
            ("A1", "A2", a),
            ("B3", "B2", b),
            ("C2", "C1", c),
        ],
    )
