r"""Plabic graphs: planar bicolored graphs in a disk, as rotation systems.

A plabic graph is a planar graph drawn in a disk with ``n`` boundary
vertices labelled clockwise and every internal vertex colored black or white
(Postnikov, *Total positivity, Grassmannians, and networks*,
arXiv:math/0609764, 2006, Definition 11.5). Reduced plabic graphs modulo
local moves index the cells of the totally nonnegative Grassmannian, and
their face labels are the maximal weakly separated collections
(Oh-Postnikov-Speyer, *Weak separation and plabic graphs*, arXiv:1109.4434,
2015). :class:`PlabicGraph` stores the graph the way the Plabic Graph page
recommends: a *rotation system* — the counterclockwise cyclic order of
half-edges (darts) at every vertex — plus the clockwise boundary cycle and
the color map. No coordinates are involved, so loops and parallel edges are
first-class (the hollow digon is a fixture), and trips, faces, the moves
(M1)-(M3), the reductions (R1)-(R2), and the reducedness criterion are all
local dart manipulations.

The module is the unweighted, combinatorial companion of
:mod:`research.boundary_measurement`, whose coordinate-embedded
:class:`~research.boundary_measurement.PlanarNetwork` and
:class:`~research.boundary_measurement.PlanarBipartiteNetwork` carry the
weights and the boundary measurement map; both convert into a
:class:`PlabicGraph`. The trip permutation leaves through
:class:`research.decorated_permutation.DecoratedPermutation` to
:class:`research.positroid.Positroid` and
:class:`research.grassmann_necklace.GrassmannNecklace`.

Two conventions are fixed throughout (the page's convention warning):
boundary vertices are labelled **clockwise**, and a trip turns sharpest
**right at black**, sharpest **left at white** — the choice shared by
Postnikov and Fomin-Williams-Zelevinsky. The trip permutation is returned in
Postnikov's direction, literally (``pi(i) = j`` when the trip from ``b_i``
ends at ``b_j``); the positroid and necklace classes exchange permutations
in the inverse (Ardila-Rincon-Williams) direction, and the conversions here
do that inversion in one documented place.

Enumeration helpers (perfect orientations, matchings, square-move classes,
weakly separated cliques) are exponential and sit behind explicit size
guards, following the page's cost note.
"""

import functools
import itertools
import math
from collections.abc import Collection, Hashable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Literal, override

import pandas as pd

from research._plot import ensure_axes, unit_circle
from research.boundary_measurement import PlanarBipartiteNetwork, PlanarNetwork
from research.decorated_permutation import DecoratedPermutation
from research.grassmann_necklace import GrassmannNecklace
from research.positroid import Positroid

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "PlabicGraph",
    "gr24_square_pair",
    "hollow_digon",
    "is_weakly_separated",
    "is_weakly_separated_collection",
    "lollipop_graph",
    "maximal_weakly_separated_collections",
    "nonreduced_wiring_example",
    "top_cell_graph",
]

BLACK = 1
"""Postnikov's ``col(v) = 1`` (Definition 11.5): turn right, one outgoing edge."""

WHITE = -1
"""Postnikov's ``col(v) = -1`` (Definition 11.5): turn left, one incoming edge."""

_COLUMNS = ("vertex", "kind", "slot", "edge", "end")
_KIND = {BLACK: "black", WHITE: "white"}
_COLOR = {"black": BLACK, "white": WHITE}

_MAX_ENUMERATION_EDGES = 40
"""Perfect-orientation and matching enumeration bound (page cost note)."""

_MAX_CLIQUE_CANDIDATES = 64
"""Bron-Kerbosch candidate bound for weakly separated collections."""

_SPHERE_EULER_CHARACTERISTIC = 2
"""Euler characteristic certifying a genus-zero (disk) embedding."""

_SQUARE_SIDES = 4
"""The quadrilateral face size the square move (M1) operates on."""

_DIGON_SIDES = 2
"""The face size of the hollow digon removed by (R1)."""

_TRIVALENT = 3
"""The vertex degree the square move and (R1) require."""

_MIDDLE_DEGREE = 2
"""The degree of a vertex removable by (M3)."""

_WEAK_SEPARATION_MAX_CHANGES = 2
"""Cyclic sign changes allowed between two weakly separated sets."""

_RELAXATION_SWEEPS = 300
"""Gauss-Seidel sweeps of the barycentric layout used for plotting."""

type _Point = tuple[Fraction, Fraction]


def _mate(dart: int) -> int:
    """Return the opposite dart of the same edge (frame darts are negative)."""
    if dart >= 0:
        return dart ^ 1
    return dart - 1 if dart % 2 else dart + 1


def _half_plane(direction: _Point) -> int:
    """Return 0 for directions with angle in ``[0, pi)``, else 1."""
    dx, dy = direction
    return 0 if dy > 0 or (dy == 0 and dx > 0) else 1


def _sorted_ccw(origin: _Point, darts: Sequence[tuple[int, _Point]]) -> list[int]:
    """Sort darts counterclockwise by the exact direction to their far end."""

    def compare(a: tuple[int, _Point], b: tuple[int, _Point]) -> int:
        da = (a[1][0] - origin[0], a[1][1] - origin[1])
        db = (b[1][0] - origin[0], b[1][1] - origin[1])
        if _half_plane(da) != _half_plane(db):
            return _half_plane(da) - _half_plane(db)
        cross = da[0] * db[1] - da[1] * db[0]
        return -1 if cross > 0 else (1 if cross < 0 else 0)

    return [dart for dart, _ in sorted(darts, key=functools.cmp_to_key(compare))]


def _hashable(label: object) -> Hashable:
    """Undo the tuple-to-list decay of a records-oriented JSON round trip."""
    if isinstance(label, list):
        return tuple(_hashable(item) for item in label)
    try:
        hash(label)
    except TypeError:
        msg = f"vertex label {label!r} is not hashable"
        raise ValueError(msg) from None
    return label


# --------------------------------------------------------------------------- #
# Mutable scratch copy used by the moves
# --------------------------------------------------------------------------- #
class _Draft:
    """A mutable half-edge copy of a plabic graph, for composing moves.

    Darts keep their ids while a draft is edited (new edges take fresh even
    ids), and :meth:`finalize` renumbers the surviving edges compactly and
    validates the result through :meth:`PlabicGraph._build`.
    """

    def __init__[W: Hashable](
        self,
        boundary: Sequence[W],
        rotations: Mapping[W, Sequence[int]],
        colors: Mapping[W, int],
    ) -> None:
        """Copy the rotation system into mutable containers."""
        self.boundary = tuple(boundary)
        self.rot: dict[Hashable, list[int]] = {
            v: list(darts) for v, darts in rotations.items()
        }
        self.colors: dict[Hashable, int] = {v: c for v, c in colors.items()}  # noqa: C416
        self.owner: dict[int, Hashable] = {
            dart: v for v, darts in self.rot.items() for dart in darts
        }
        self._next_edge = max((dart >> 1 for dart in self.owner), default=-1) + 1
        self._fresh = itertools.count()

    def new_edge(self) -> tuple[int, int]:
        """Allocate a fresh edge and return its two darts."""
        dart = 2 * self._next_edge
        self._next_edge += 1
        return dart, dart + 1

    def fresh_label(self) -> Hashable:
        """Return a vertex label ``("sq", k)`` not present in the draft."""
        while True:
            label = ("sq", next(self._fresh))
            if label not in self.rot:
                return label

    def add_vertex(self, label: Hashable, darts: Sequence[int], color: int) -> None:
        """Add an internal vertex owning ``darts`` in counterclockwise order."""
        if label in self.rot:
            msg = f"vertex label {label!r} is already in use"
            raise ValueError(msg)
        self.rot[label] = list(darts)
        self.colors[label] = color
        for dart in darts:
            self.owner[dart] = label

    def drop_vertex(self, label: Hashable) -> None:
        """Remove a vertex together with the darts it still owns."""
        for dart in self.rot.pop(label):
            if self.owner.get(dart) == label:
                del self.owner[dart]
        self.colors.pop(label, None)

    def replace(self, old: int, new: int) -> None:
        """Substitute dart ``new`` for ``old`` in place in its rotation."""
        v = self.owner.pop(old)
        self.rot[v][self.rot[v].index(old)] = new
        self.owner[new] = v

    def contract(self, dart: int) -> None:
        """Merge the far endpoint of ``dart``'s edge into its near endpoint.

        The far vertex's darts, read counterclockwise from just after the
        contracted edge, take the contracted dart's slot — the rotation of
        the merged vertex in the drawing.
        """
        far_dart = _mate(dart)
        near, far = self.owner[dart], self.owner[far_dart]
        far_rot = self.rot[far]
        cut = far_rot.index(far_dart)
        moved = far_rot[cut + 1 :] + far_rot[:cut]
        near_rot = self.rot[near]
        slot = near_rot.index(dart)
        self.rot[near] = near_rot[:slot] + moved + near_rot[slot + 1 :]
        for carried in moved:
            self.owner[carried] = near
        del self.owner[dart], self.owner[far_dart], self.rot[far], self.colors[far]

    def uncontract(
        self, label: Hashable, vertex: Hashable, block: Sequence[int]
    ) -> None:
        """Split ``block`` (consecutive darts of ``vertex``) off to ``label``."""
        near, far = self.new_edge()
        rotation = self.rot[vertex]
        start = rotation.index(block[0])
        rotated = rotation[start:] + rotation[:start]
        self.rot[vertex] = [near, *rotated[len(block) :]]
        self.owner[near] = vertex
        self.add_vertex(label, [far, *block], self.colors[vertex])

    def splice_out(self, vertex: Hashable) -> None:
        """Remove a degree-two vertex, gluing its two edges into one (M3)."""
        first, second = self.rot[vertex]
        near, far = self.new_edge()
        self.replace(_mate(first), near)
        self.replace(_mate(second), far)
        self.drop_vertex(vertex)

    def _contractible(self) -> int | None:
        """Return a dart of some unicolored internal non-loop edge, if any."""
        for dart, v in self.owner.items():
            w = self.owner[_mate(dart)]
            shade = self.colors.get(v)
            if v != w and shade is not None and shade == self.colors.get(w):
                return dart
        return None

    def _removable(self) -> Hashable | None:
        """Return some degree-two internal vertex that is not a bare loop."""
        for v, darts in self.rot.items():
            middle = v in self.colors and len(darts) == _MIDDLE_DEGREE
            if middle and _mate(darts[0]) != darts[1]:
                return v
        return None

    def normalize(self) -> None:
        """Exhaust (M2) contractions and (M3) middle-vertex removals."""
        while True:
            dart = self._contractible()
            if dart is not None:
                self.contract(dart)
                continue
            vertex = self._removable()
            if vertex is None:
                return
            self.splice_out(vertex)

    def finalize(self) -> PlabicGraph[Hashable]:
        """Renumber the edges compactly and build the validated graph."""
        live = sorted({dart >> 1 for dart in self.owner})
        renumber = {old: new for new, old in enumerate(live)}
        edges: list[tuple[Hashable, Hashable]] = []
        for old in live:
            if 2 * old not in self.owner or 2 * old + 1 not in self.owner:
                msg = f"edge {old} lost one of its endpoints during a move"
                raise ValueError(msg)
            edges.append((self.owner[2 * old], self.owner[2 * old + 1]))
        rotations = {
            v: [2 * renumber[dart >> 1] + (dart & 1) for dart in darts]
            for v, darts in self.rot.items()
        }
        return PlabicGraph._build(self.boundary, edges, rotations, self.colors)


# --------------------------------------------------------------------------- #
# The structure
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PlabicGraph[V: Hashable]:
    """A plabic graph, stored as a rotation system in a disk.

    Postnikov Definition 11.5, with the connectivity clause of
    Fomin-Williams-Zelevinsky Definition 7.1.1. Edge ``e`` of ``edges`` has
    two darts: ``2e`` at ``edges[e][0]`` and ``2e + 1`` at ``edges[e][1]``.
    ``rotations`` lists, for every vertex (boundary vertices first, in
    clockwise boundary order), its darts in **counterclockwise** order,
    rotated to start at the smallest dart. ``colors`` gives ``col(v)`` in
    ``{1, -1}`` for the internal vertices: 1 black, -1 white.

    The definition is checked as the page numbers it:

    * (P1) ``n`` distinct boundary vertices labelled clockwise, all other
      vertices internal; loops and multiple edges allowed. Being *drawn in
      the disk* is checked combinatorially: closing the rotation system
      with the frame of boundary arcs must give a sphere (Euler
      characteristic 2).
    * (P2) each boundary vertex is incident to a single edge.
    * (P3) each internal vertex carries a color in ``{1, -1}``.
    * FWZ Definition 7.1.1: every internal vertex is joined by a path to
      the boundary. This is a **representation restriction** as much as a
      convention: a bare rotation system cannot say which face a floating
      component sits in, so Postnikov's isolated components — and with
      them the dipole reduction (R3) — are not representable.

    Build through :meth:`from_rotation_system` or the other ``from_*``
    constructors, which validate; calling the dataclass constructor
    directly skips validation.
    """

    boundary: tuple[V, ...]
    edges: tuple[tuple[V, V], ...]
    rotations: tuple[tuple[V, tuple[int, ...]], ...]
    colors: tuple[tuple[V, int], ...]

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        black = sum(1 for _, color in self.colors if color == BLACK)
        return (
            f"PlabicGraph(n={len(self.boundary)}, edges={len(self.edges)}, "
            f"black={black}, white={len(self.colors) - black}, "
            f"trips={self.trip_permutation!r})"
        )

    # ---------------------------------------------------------- construction
    @classmethod
    def _build[W: Hashable](
        cls,
        boundary: Sequence[W],
        edges: Sequence[tuple[W, W]],
        rotations: Mapping[W, Sequence[int]],
        colors: Mapping[W, int],
    ) -> PlabicGraph[W]:
        """Canonicalize a dart-level rotation system, then validate it."""
        on_boundary = set(boundary)
        ordered = [*boundary, *(v for v in rotations if v not in on_boundary)]
        canonical: list[tuple[W, tuple[int, ...]]] = []
        for v in ordered:
            darts = list(rotations.get(v, ()))
            if darts:
                start = darts.index(min(darts))
                darts = darts[start:] + darts[:start]
            canonical.append((v, tuple(darts)))
        internal = [v for v in ordered if v not in on_boundary]
        extra = [v for v in colors if v not in rotations or v in on_boundary]
        graph = PlabicGraph(
            tuple(boundary),
            tuple((u, w) for u, w in edges),
            tuple(canonical),
            tuple((v, colors[v]) for v in [*internal, *extra] if v in colors),
        )
        graph._validate()
        return graph

    @classmethod
    def from_rotation_system(
        cls,
        boundary: Sequence[V],
        edges: Sequence[tuple[V, V]],
        rotations: Mapping[V, Sequence[int]],
        colors: Mapping[V, int],
    ) -> PlabicGraph[V]:
        """Build and validate a plabic graph from its rotation system.

        The primary formulation (Postnikov Definition 11.5), in the
        combinatorial-embedding form of the page's representation note.

        Args:
            boundary: Boundary vertex labels in clockwise order (P1).
            edges: Endpoint pairs; loops and repeated pairs are allowed
                (P1).
            rotations: For each internal vertex, the indices of its
                incident edges in **counterclockwise** order; a loop's
                index appears twice. Boundary vertices need no entry —
                each has the single edge (P2) found in ``edges``.
            colors: ``col(v)`` for each internal vertex: 1 black, -1 white
                (P3).

        Returns:
            The validated graph.

        Raises:
            ValueError: Naming the violated axiom — (P1), (P2), (P3), or
                the connectivity clause of FWZ Definition 7.1.1.
        """
        on_boundary = set(boundary)
        darts: dict[V, list[int]] = {b: [] for b in boundary}
        for index, (u, w) in enumerate(edges):
            if u in on_boundary:
                darts[u].append(2 * index)
            if w in on_boundary:
                darts[w].append(2 * index + 1)
        for v, incident in rotations.items():
            if v in on_boundary:
                continue
            seen_loops: set[int] = set()
            darts[v] = []
            for index in incident:
                if not 0 <= index < len(edges) or v not in edges[index]:
                    msg = (
                        f"(P1) violated: the rotation at {v!r} lists edge "
                        f"{index}, which is not an edge at that vertex"
                    )
                    raise ValueError(msg)
                u, w = edges[index]
                second = (u == w and index in seen_loops) or u != v
                seen_loops.add(index)
                darts[v].append(2 * index + int(second))
        return cls._build(boundary, edges, darts, colors)

    @classmethod
    def from_planar_bipartite_network[W: Hashable](
        cls, network: PlanarBipartiteNetwork[W]
    ) -> PlabicGraph[W]:
        """Forget the weights and coordinates of a planar bipartite network.

        The rotation at each vertex is the exact counterclockwise angular
        order of its straight-line edges; internal colors are kept and the
        boundary colors (a device of Lam's matching formulation) dropped,
        since plabic boundary vertices are uncolored (Postnikov Definition
        11.5, (P3)).
        """
        edges = [(u, w) for u, w, _ in network.edges]
        internal = network.internal_vertices
        colors = {
            v: BLACK if v in network.black_vertices else WHITE
            for v, _ in network.positions
            if v in internal
        }
        rotations = _geometric_rotations(dict(network.positions), edges)
        return PlabicGraph._build(network.boundary, edges, rotations, colors)

    @classmethod
    def from_planar_network[W: Hashable](
        cls, network: PlanarNetwork[W]
    ) -> PlabicGraph[Hashable]:
        """Forget a perfect orientation, keeping the colors it determines.

        Postnikov Definition 11.5 (perfect orientation): a vertex with
        exactly one outgoing edge is black and one with exactly one
        incoming edge is white; a pass-through vertex (one in, one out)
        may take either color by (M3) and is made white. An isolated
        boundary source becomes a white lollipop and an isolated sink a
        black one, labelled ``(b, "lollipop")``, so that the boundary
        vertex has its single edge (P2) and the orientation's source set
        is kept.

        Raises:
            ValueError: If some internal vertex has neither a unique
                outgoing nor a unique incoming edge — the network is not
                perfectly oriented.
        """
        return _from_oriented(network, split=False)

    @classmethod
    def from_le_diagram(
        cls, filling: Sequence[Sequence[int]], n: int
    ) -> PlabicGraph[Hashable]:
        """Return the Le-graph ``G_D`` of a Le-diagram — a reduced graph.

        Postnikov section 20 (the page's "canonical graph per cell"): the
        Gamma-graph of the diagram, built by
        :meth:`research.boundary_measurement.PlanarNetwork.from_le_diagram`,
        with each 4-valent vertex split into a pair of trivalent vertices.
        The split used here is the one compatible with the Gamma-graph's
        orientation and with planarity: the two incoming edges (from the
        right and from above) stay on a black vertex, the two outgoing
        edges (left and down) move to a white vertex ``(v, "white")``, and
        the pair is joined by an edge. The page records that Postnikov
        prescribes a coloring without reproducing his Figure 20.1; this
        choice is pinned by the oracle of his Corollary 20.1 (the
        decorated trip permutation is the cell's), checked exhaustively in
        the tests. Empty rows give white lollipops and empty columns black
        ones.

        Args:
            filling: Rows of 0/1 values with weakly decreasing lengths.
            n: The number of boundary vertices; the shape must fit in the
                ``k x (n - k)`` rectangle, ``k = len(filling)``.

        Raises:
            ValueError: If the filling is not a Le-diagram of that type.
        """
        return _from_oriented(PlanarNetwork.from_le_diagram(filling, n), split=True)

    @classmethod
    def from_wiring_diagram(
        cls, word: Sequence[int], n: int
    ) -> PlabicGraph[int | tuple[int, int]]:
        """Return the plabic graph of the wiring diagram of a word in ``S_n``.

        Fomin-Williams-Zelevinsky Examples 7.3.4 and 7.3.11, Remark 7.3.8:
        ``n`` horizontal lines at heights ``1..n`` (bottom to top), and for
        the ``j``-th letter ``s_i`` a vertical edge between lines ``i`` and
        ``i + 1`` whose upper endpoint ``(j, 1)`` is black and lower
        endpoint ``(j, -1)`` white ("a vertical black-over-white pair").
        The ``2n`` boundary vertices are labelled clockwise: ``h`` is the
        left end of line ``h`` and ``2n + 1 - h`` its right end. With the
        rules of the road, the trips from the left follow the wires and
        the trips from the right run straight. The word is a reduced
        expression iff the graph is reduced (FWZ Exercise 7.3.9, Remark
        7.3.10).

        Args:
            word: Letters ``i`` of the generators ``s_i``, ``1 <= i < n``.
            n: The number of wires.

        Raises:
            ValueError: If ``n < 1`` or a letter is out of range.
        """
        if n < 1 or any(not 1 <= letter < n for letter in word):
            msg = f"a word in S_{n} uses letters 1..{n - 1}; got {list(word)!r}"
            raise ValueError(msg)
        type Label = int | tuple[int, int]
        lines: dict[int, list[Label]] = {h: [h] for h in range(1, n + 1)}
        for j, letter in enumerate(word):
            lines[letter].append((j, WHITE))
            lines[letter + 1].append((j, BLACK))
        edges: list[tuple[Label, Label]] = []
        east: dict[Label, int] = {}
        west: dict[Label, int] = {}
        for h, chain in lines.items():
            for a, b in itertools.pairwise([*chain, 2 * n + 1 - h]):
                east[a] = west[b] = len(edges)
                edges.append((a, b))
        rotations: dict[Label, list[int]] = {}
        for j in range(len(word)):
            bridge = len(edges)
            edges.append(((j, BLACK), (j, WHITE)))
            rotations[(j, WHITE)] = [east[(j, WHITE)], bridge, west[(j, WHITE)]]
            rotations[(j, BLACK)] = [east[(j, BLACK)], west[(j, BLACK)], bridge]
        colors: dict[Label, int] = {v: v[1] for v in rotations if isinstance(v, tuple)}
        return PlabicGraph.from_rotation_system(
            tuple(range(1, 2 * n + 1)), edges, rotations, colors
        )

    # ------------------------------------------------------------ validation
    def _validate(self) -> None:
        """Check the definition; raise ``ValueError`` naming the axiom."""
        self._validate_vertices()
        self._validate_darts()
        self._validate_disk()

    def _validate_vertices(self) -> None:
        """Check the boundary labels and the color map: (P1), (P3)."""
        if len(set(self.boundary)) != len(self.boundary):
            msg = (
                f"(P1) violated: the boundary vertices b_1..b_n must be "
                f"distinct, got {self.boundary!r}"
            )
            raise ValueError(msg)
        colored = dict(self.colors)
        vertices = [v for v, _ in self.rotations]
        for v in vertices:
            if v not in self._on_boundary and v not in colored:
                msg = f"(P3) violated: internal vertex {v!r} carries no color"
                raise ValueError(msg)
        for v, color in self.colors:
            if v in self._on_boundary or v not in set(vertices):
                msg = (
                    f"(P3) violated: only internal vertices are colored, but "
                    f"{v!r} is not an internal vertex"
                )
                raise ValueError(msg)
            if color not in (BLACK, WHITE):
                msg = (
                    f"(P3) violated: col({v!r}) = {color!r} is not in "
                    f"{{1, -1}} (black / white)"
                )
                raise ValueError(msg)

    def _validate_darts(self) -> None:
        """Check the rotation system and boundary degrees: (P1), (P2)."""
        expected: dict[int, V] = {}
        for index, (u, w) in enumerate(self.edges):
            expected[2 * index] = u
            expected[2 * index + 1] = w
        listed = [(dart, v) for v, darts in self.rotations for dart in darts]
        for b in self.boundary:
            degree = sum(1 for _, v in listed if v == b)
            if degree != 1:
                msg = (
                    f"(P2) violated: boundary vertex {b!r} is incident to "
                    f"{degree} edges; each boundary vertex is incident to a "
                    f"single edge"
                )
                raise ValueError(msg)
        if sorted(dart for dart, _ in listed) != sorted(expected) or any(
            expected[dart] != v for dart, v in listed
        ):
            msg = (
                "(P1) violated: the rotations must list every half-edge "
                "exactly once, at the endpoint its edge names"
            )
            raise ValueError(msg)

    def _validate_disk(self) -> None:
        """Check FWZ connectivity and that the drawing is planar: (P1)."""
        reached = set(self.boundary)
        frontier = list(self.boundary)
        while frontier:
            v = frontier.pop()
            for dart in self._rotation[v]:
                w = self._at[_mate(dart)]
                if w not in reached:
                    reached.add(w)
                    frontier.append(w)
        stranded = [v for v, _ in self.rotations if v not in reached]
        if stranded:
            msg = (
                f"FWZ Definition 7.1.1 violated: internal vertices "
                f"{stranded!r} are not connected by a path to the boundary "
                f"(isolated components are not representable)"
            )
            raise ValueError(msg)
        if not self.boundary:
            return
        orbits = len(self._orbits)
        euler = len(self.rotations) - (len(self.edges) + len(self.boundary)) + orbits
        if euler != _SPHERE_EULER_CHARACTERISTIC:
            msg = (
                f"(P1) violated: the rotation system is not a planar drawing "
                f"in the disk with the boundary labelled clockwise (Euler "
                f"characteristic {euler}, expected 2)"
            )
            raise ValueError(msg)

    # ------------------------------------------------------------- structure
    @functools.cached_property
    def _on_boundary(self) -> frozenset[V]:
        """The boundary vertices as a set."""
        return frozenset(self.boundary)

    @functools.cached_property
    def _rotation(self) -> dict[V, tuple[int, ...]]:
        """The counterclockwise darts at each vertex."""
        return dict(self.rotations)

    @functools.cached_property
    def _color(self) -> dict[V, int]:
        """The color map on internal vertices."""
        return dict(self.colors)

    @functools.cached_property
    def _at(self) -> dict[int, V]:
        """The vertex each dart leaves from."""
        return {dart: v for v, darts in self.rotations for dart in darts}

    @functools.cached_property
    def _turn(self) -> dict[int, int]:
        """The dart a trip continues on after arriving by each dart.

        Rules of the road (Postnikov section 13; FWZ Definition 7.1.8): at
        a black vertex the sharpest right turn is the counterclockwise
        successor of the arrival dart, at a white vertex the sharpest left
        turn is its predecessor; at a degree-one vertex both are a U-turn.
        """
        turn: dict[int, int] = {}
        for v, darts in self.rotations:
            if v in self._on_boundary:
                continue
            step = 1 if self._color[v] == BLACK else -1
            for slot, dart in enumerate(darts):
                turn[dart] = darts[(slot + step) % len(darts)]
        return turn

    @property
    def size(self) -> int:
        """The number ``n`` of boundary vertices."""
        return len(self.boundary)

    @functools.cached_property
    def internal_vertices(self) -> tuple[V, ...]:
        """The internal vertices, in stored order."""
        return tuple(v for v, _ in self.colors)

    def degree(self, vertex: V) -> int:
        """Return the number of half-edges at ``vertex`` (a loop counts twice)."""
        return len(self._rotation[vertex])

    def color(self, vertex: V) -> int:
        """Return ``col(v)``: 1 black, -1 white (Postnikov Definition 11.5).

        Raises:
            ValueError: If ``vertex`` is a boundary vertex, which carries
                no color.
        """
        if vertex not in self._color:
            msg = f"{vertex!r} is not an internal vertex and carries no color"
            raise ValueError(msg)
        return self._color[vertex]

    def _neighbor(self, vertex: V) -> V:
        """Return the far end of a degree-one vertex's only edge."""
        return self._at[_mate(self._rotation[vertex][0])]

    @functools.cached_property
    def lollipops(self) -> dict[V, int]:
        """The boundary vertices carrying a lollipop, with its color.

        FWZ Definition 7.1.1: a lollipop is a degree-one internal vertex
        attached to a boundary vertex. A black lollipop makes its label a
        loop of the positroid and a white one a coloop (Postnikov, proof
        of Proposition 16.4).
        """
        found: dict[V, int] = {}
        for b in self.boundary:
            head = self._neighbor(b)
            if head in self._color and self.degree(head) == 1:
                found[b] = self._color[head]
        return found

    @functools.cached_property
    def leaves(self) -> tuple[V, ...]:
        """The degree-one internal vertices that are not lollipops.

        FWZ Definition 7.1.1: "any other degree-1 internal vertex is a
        leaf".
        """
        return tuple(
            v
            for v in self.internal_vertices
            if self.degree(v) == 1 and self._neighbor(v) not in self._on_boundary
        )

    @property
    def is_leafless(self) -> bool:
        """Whether the graph has no leaves (FWZ Definition 7.1.1)."""
        return not self.leaves

    # ------------------------------------------------------------------ trips
    @functools.cached_property
    def _trips(self) -> tuple[tuple[int, ...], ...]:
        """All trips as dart sequences: one-way trips first, then roundtrips."""
        trips: list[tuple[int, ...]] = []
        covered: set[int] = set()
        for b in self.boundary:
            dart = self._rotation[b][0]
            walk = [dart]
            while _mate(dart) in self._turn:
                dart = self._turn[_mate(dart)]
                walk.append(dart)
            covered.update(walk)
            trips.append(tuple(walk))
        for start in sorted(self._turn):
            if start in covered:
                continue
            walk = [start]
            dart = self._turn[_mate(start)]
            while dart != start:
                walk.append(dart)
                dart = self._turn[_mate(dart)]
            covered.update(walk)
            trips.append(tuple(walk))
        return tuple(trips)

    def trips(self) -> tuple[tuple[int, ...], ...]:
        """Return the trips as sequences of darts.

        Postnikov section 13; FWZ Definition 7.1.8. The first ``n`` entries
        are the one-way trips, entry ``i - 1`` starting at ``b_i``; any
        further entries are the roundtrips. Dart ``d`` traverses edge
        ``d // 2`` from its endpoint ``d % 2`` to the other. The trips
        together use every edge once in each direction.
        """
        return self._trips

    @property
    def roundtrips(self) -> tuple[tuple[int, ...], ...]:
        """The closed trips avoiding the boundary (Postnikov section 13)."""
        return self._trips[len(self.boundary) :]

    @functools.cached_property
    def trip_permutation(self) -> tuple[int, ...]:
        """The trip permutation ``pi_G`` on boundary positions ``1..n``.

        Postnikov section 13: ``pi_G(i) = j`` when the trip starting at
        ``b_i`` ends at ``b_j``. Invariant under the moves (M1)-(M3)
        (Postnikov Lemma 13.1).
        """
        position = {b: index for index, b in enumerate(self.boundary, start=1)}
        return tuple(
            position[self._at[_mate(trip[-1])]]
            for trip in self._trips[: len(self.boundary)]
        )

    def _essential_edges(self) -> list[int]:
        """Return the edges joining internal vertices of different colors."""
        return [
            index
            for index, (u, w) in enumerate(self.edges)
            if u in self._color
            and w in self._color
            and self._color[u] != self._color[w]
        ]

    def essential_self_intersections(self) -> tuple[int, ...]:
        """Return the edges at which some trip essentially meets itself.

        Postnikov section 13: a trip has an essential self-intersection at
        an edge it passes twice, in opposite directions, when the edge's
        endpoints have different colors. Boundary vertices are uncolored,
        so a lollipop's U-turn is not one.
        """
        trip_of = {dart: t for t, trip in enumerate(self._trips) for dart in trip}
        return tuple(
            index
            for index in self._essential_edges()
            if trip_of[2 * index] == trip_of[2 * index + 1]
        )

    def bad_double_crossings(self) -> tuple[tuple[int, int], ...]:
        """Return the edge pairs ``(e1, e2)`` witnessing a bad double crossing.

        Postnikov section 13: two distinct trips with essential
        intersections at ``e1`` and ``e2`` (each passes both edges, in
        opposite directions, and both edges join different colors), both
        trips directed from ``e1`` to ``e2``. Only one-way trips are
        compared — along a roundtrip "from ``e1`` to ``e2``" has no
        meaning, and a roundtrip already fails the reducedness criterion.
        """
        n = len(self.boundary)
        place = {
            dart: (t, step)
            for t, trip in enumerate(self._trips[:n])
            for step, dart in enumerate(trip)
        }
        shared: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
        for index in self._essential_edges():
            if 2 * index not in place or 2 * index + 1 not in place:
                continue
            (s, s_step), (t, t_step) = sorted((place[2 * index], place[2 * index + 1]))
            if s != t:
                shared.setdefault((s, t), []).append((index, s_step, t_step))
        return tuple(
            (first[0], second[0]) if first[1] < second[1] else (second[0], first[0])
            for crossings in shared.values()
            for first, second in itertools.combinations(crossings, 2)
            if (first[1] < second[1]) == (first[2] < second[2])
        )

    def _require_leafless(self, statement: str) -> None:
        """Raise unless the graph is leafless, the hypothesis of ``statement``."""
        if self.leaves:
            msg = (
                f"{statement} is stated for leafless plabic graphs, but "
                f"{list(self.leaves)!r} are leaves"
            )
            raise ValueError(msg)

    def is_reduced(self) -> bool:
        """Return whether the graph is reduced, by the trip criterion.

        Postnikov Theorem 13.2: a leafless plabic graph without isolated
        components is reduced iff (1) it has no roundtrips, (2) no trip
        has an essential self-intersection, (3) no pair of trips has a bad
        double crossing, and (4) every fixed point of ``pi_G`` comes from
        a boundary lollipop. (Isolated components cannot occur in this
        representation.)

        Raises:
            ValueError: If the graph has leaves — outside the theorem's
                hypothesis.
        """
        self._require_leafless("the reducedness criterion (Postnikov Theorem 13.2)")
        fixed = {
            b
            for index, b in enumerate(self.boundary, start=1)
            if self.trip_permutation[index - 1] == index
        }
        return (
            not self.roundtrips
            and not self.essential_self_intersections()
            and not self.bad_double_crossings()
            and fixed <= set(self.lollipops)
        )

    def has_resonance_property(self) -> bool:
        """Return whether the graph has the resonance property.

        Kodama-Williams, via FWZ Theorem 7.11.5 (a leafless graph is
        reduced iff it has the property): label each edge by the two
        one-way trips through it, a trip being named by the boundary
        position where it **ends**; then around every internal
        non-lollipop vertex the labels must read ``{i1 i2}, {i2 i3}, ...,
        {i(m-1) im}, {i1 im}`` clockwise for some ``i1 < ... < im``.

        Raises:
            ValueError: If the graph has leaves — outside the theorem's
                hypothesis.
        """
        self._require_leafless("the resonance criterion (FWZ Theorem 7.11.5)")
        end = {
            dart: self.trip_permutation[t]
            for t, trip in enumerate(self._trips[: len(self.boundary)])
            for dart in trip
        }
        for v in self.internal_vertices:
            darts = self._rotation[v]
            if len(darts) == 1:
                continue
            if any(d not in end or _mate(d) not in end for d in darts):
                return False
            clockwise = [frozenset({end[d], end[_mate(d)]}) for d in reversed(darts)]
            indices = sorted(set().union(*clockwise))
            if len(indices) != len(darts):
                return False
            wanted = [frozenset(pair) for pair in itertools.pairwise(indices)]
            wanted.append(frozenset({indices[0], indices[-1]}))
            if not any(
                clockwise[shift:] + clockwise[:shift] == wanted
                for shift in range(len(darts))
            ):
                return False
        return True

    # ------------------------------------------------------------------ faces
    @functools.cached_property
    def _orbits(self) -> tuple[tuple[int, ...], ...]:
        """The face orbits of the framed rotation system, outer face last.

        The frame closes the disk: arc ``i`` runs clockwise from boundary
        position ``i`` to ``i + 1`` with darts ``-(2i + 1)`` (clockwise)
        and ``-(2i + 2)`` (counterclockwise), and the counterclockwise
        rotation at a boundary vertex reads (arc to the next vertex, arc
        to the previous one, its edge). Following "arrive, then take the
        counterclockwise successor" keeps a face on the right, so every
        orbit is the set of darts with that face to their right, and the
        counterclockwise arcs bound the outside of the disk.
        """
        n = len(self.boundary)
        successor: dict[int, int] = {}
        for v, darts in self.rotations:
            ring = list(darts)
            if v in self._on_boundary:
                i = self.boundary.index(v)
                ring = [-(2 * i + 1), -(2 * ((i - 1) % n) + 2), *darts]
            for slot, dart in enumerate(ring):
                successor[dart] = ring[(slot + 1) % len(ring)]
        starts = [-(2 * i + 1) for i in range(n)]
        starts += sorted(d for d in successor if d >= 0)
        starts += [-(2 * i + 2) for i in range(n)]
        seen: set[int] = set()
        orbits: list[tuple[int, ...]] = []
        for start in starts:
            if start in seen:
                continue
            orbit = [start]
            dart = successor[_mate(start)]
            while dart != start:
                orbit.append(dart)
                dart = successor[_mate(dart)]
            seen.update(orbit)
            orbits.append(tuple(orbit))
        return tuple(orbits)

    @functools.cached_property
    def _faces(self) -> tuple[tuple[int, ...], ...]:
        """The faces inside the disk, each as its orbit with frame darts."""
        if not self.boundary:
            return ((),)
        return self._orbits[:-1]

    def faces(self) -> tuple[tuple[int, ...], ...]:
        """Return the faces of the graph inside the disk.

        Each face is the cyclic sequence of darts that have it on their
        right; dart ``d`` runs along edge ``d // 2``. Faces touching the
        boundary come first, ordered by the first boundary arc they meet
        (the arc from ``b_1`` to ``b_2`` first). The count needs no
        drawing (the page's representation note).
        """
        return tuple(tuple(d for d in face if d >= 0) for face in self._faces)

    @property
    def face_count(self) -> int:
        """The number of faces ``|F(G)|`` inside the disk.

        For a reduced graph this is the cell dimension plus one (Postnikov
        Theorem 12.7), that is ``l(pi) + 1`` (Oh-Postnikov-Speyer Theorem
        6.8).
        """
        return len(self._faces)

    @functools.cached_property
    def _face_of(self) -> dict[int, int]:
        """The index of the face to the right of each real dart."""
        return {
            dart: index
            for index, face in enumerate(self._faces)
            for dart in face
            if dart >= 0
        }

    def _left_faces(self, trip: Sequence[int]) -> frozenset[int]:
        """Return the faces to the left of a one-way trip.

        The trip is a curve from boundary to boundary; two faces across an
        edge lie on opposite sides exactly when the trip runs along that
        edge once. A lollipop's U-turn separates nothing, and its side is
        fixed by the lollipop block of the page: a white lollipop is a
        coloop (in every label), a black one a loop (in none).
        """
        runs: dict[int, int] = {}
        for dart in trip:
            runs[dart >> 1] = runs.get(dart >> 1, 0) + 1
        start = self._at[trip[0]]
        if start in self.lollipops and len(trip) == _DIGON_SIDES:
            everything = frozenset(range(len(self._faces)))
            return everything if self.lollipops[start] == WHITE else frozenset()
        side = {self._face_of[_mate(trip[0])]: True}
        frontier = list(side)
        while frontier:
            face = frontier.pop()
            for dart in self._faces[face]:
                if dart < 0:
                    continue
                across = self._face_of[_mate(dart)]
                flipped = side[face] ^ (runs.get(dart >> 1, 0) == 1)
                if across not in side:
                    side[across] = flipped
                    frontier.append(across)
                elif side[across] != flipped:
                    msg = (
                        "a trip crosses itself, so its left side is not "
                        "defined; face labels need a reduced graph"
                    )
                    raise ValueError(msg)
        return frozenset(face for face, left in side.items() if left)

    def _face_labels(self, labeling: str) -> tuple[frozenset[V], ...]:
        """Compute the face labels without re-checking reducedness."""
        members: list[set[V]] = [set() for _ in self._faces]
        for t, trip in enumerate(self._trips[: len(self.boundary)]):
            name = (
                self.boundary[self.trip_permutation[t] - 1]
                if labeling == "target"
                else self.boundary[t]
            )
            for face in self._left_faces(trip):
                members[face].add(name)
        return tuple(frozenset(label) for label in members)

    def face_labels(
        self, labeling: Literal["target", "source"] = "target"
    ) -> tuple[frozenset[V], ...]:
        """Return the label of every face, aligned with :meth:`faces`.

        FWZ Definition 7.12.2 (the construction first published by Scott):
        the *target* label of a face is the set of ``i`` such that the
        face lies to the left of the trip **ending** at ``b_i``; the
        *source* labeling uses the trip starting at ``b_i``. Every face of
        a reduced graph of type ``(k, n)`` gets a ``k``-subset (FWZ
        Theorem 7.12.4). Labels are sets of boundary labels.

        Raises:
            ValueError: If the graph is not reduced (or has leaves) — the
                page defines face labels for reduced graphs only.
        """
        if not self.is_reduced():
            msg = "face labels are defined for reduced plabic graphs (FWZ 7.12.2)"
            raise ValueError(msg)
        return self._face_labels(labeling)

    def boundary_face_labels(self) -> tuple[frozenset[V], ...]:
        """Return the target labels of the boundary faces, in necklace order.

        Entry ``i - 1`` is the label of the face touching the boundary arc
        that ends at ``b_i`` (the arc between ``b_(i-1)`` and ``b_i``). In
        this order the boundary faces read off the Grassmann necklace
        ``(I_1, ..., I_n)`` of the positroid (Oh-Postnikov-Speyer section
        6).

        Raises:
            ValueError: If the graph is not reduced.
        """
        labels = self.face_labels("target")
        n = len(self.boundary)
        arc_face = {
            -dart // 2: index
            for index, face in enumerate(self._faces)
            for dart in face
            if dart < 0
        }
        return tuple(labels[arc_face[(i - 1) % n]] for i in range(n))

    # ---------------------------------------------- orientations and matchings
    def _guard_enumeration(self, what: str) -> None:
        """Raise when the graph is too large to enumerate ``what``."""
        if len(self.edges) > _MAX_ENUMERATION_EDGES:
            msg = (
                f"{len(self.edges)} edges exceed the enumeration guard of "
                f"{_MAX_ENUMERATION_EDGES} for {what} (page cost note)"
            )
            raise ValueError(msg)

    def _orientations(self) -> Iterator[tuple[int, ...]]:
        """Yield the perfect orientations lazily, as tuples of tail darts."""
        special = dict.fromkeys(self._color, 0)
        remaining = {v: self.degree(v) for v in self._color}
        chosen: list[int] = []

        def counts(dart: int, tail: bool) -> bool:
            v = self._at[dart]
            return v in self._color and (self._color[v] == BLACK) == tail

        def feasible(v: V) -> bool:
            if v not in self._color:
                return True
            return special[v] <= 1 and special[v] + remaining[v] >= 1

        def extend(index: int) -> Iterator[tuple[int, ...]]:
            if index == len(self.edges):
                yield tuple(chosen)
                return
            for tail in (2 * index, 2 * index + 1):
                head = _mate(tail)
                touched = [self._at[tail], self._at[head]]
                bumps = [(tail, counts(tail, True)), (head, counts(head, False))]
                for dart, bump in bumps:
                    if self._at[dart] in self._color:
                        remaining[self._at[dart]] -= 1
                        special[self._at[dart]] += int(bump)
                if all(feasible(v) for v in touched):
                    chosen.append(tail)
                    yield from extend(index + 1)
                    chosen.pop()
                for dart, bump in bumps:
                    if self._at[dart] in self._color:
                        remaining[self._at[dart]] += 1
                        special[self._at[dart]] -= int(bump)

        return extend(0)

    def perfect_orientations(self) -> tuple[tuple[int, ...], ...]:
        """Enumerate the perfect orientations of the graph.

        Postnikov Definition 11.5: every black vertex has exactly one
        outgoing edge and every white vertex exactly one incoming edge.
        Each orientation is a tuple with one entry per edge: the dart at
        the edge's tail (``2e`` when edge ``e`` points from ``edges[e][0]``
        to ``edges[e][1]``, else ``2e + 1``). Exponential; guarded.

        Raises:
            ValueError: If the edge count exceeds the enumeration guard.
        """
        self._guard_enumeration("perfect orientations")
        return tuple(self._orientations())

    @functools.cached_property
    def is_perfectly_orientable(self) -> bool:
        """Whether a perfect orientation exists (Postnikov Definition 11.5).

        Every reduced graph is perfectly orientable (Postnikov Theorem
        12.7). Stops at the first orientation found; guarded like
        :meth:`perfect_orientations`.
        """
        self._guard_enumeration("perfect orientations")
        return next(self._orientations(), None) is not None

    def source_set(self, orientation: Sequence[int]) -> frozenset[V]:
        """Return the boundary sources ``I_O`` of a perfect orientation.

        Postnikov Definition 11.5; ``orientation`` is an entry of
        :meth:`perfect_orientations`.
        """
        return frozenset(
            self._at[tail]
            for tail in orientation
            if self._at[tail] in self._on_boundary
        )

    @functools.cached_property
    def graph_type(self) -> tuple[int, int]:
        r"""The type ``(k, n)`` of a perfectly orientable graph.

        Postnikov Lemma 9.4, as used in Definition 11.5: every perfect
        orientation has ``k`` sources, where
        ``k - (n - k) = \sum_v col(v) (deg(v) - 2)`` over internal
        vertices. (The page notes the arXiv text of Definition 11.5 prints
        ``k + (n - k)``, an evident typo.)

        Raises:
            ValueError: If the graph is not perfectly orientable.
        """
        if not self.is_perfectly_orientable:
            msg = "the type (k, n) is defined for perfectly orientable graphs"
            raise ValueError(msg)
        n = len(self.boundary)
        excess = sum(color * (self.degree(v) - 2) for v, color in self.colors)
        return (n + excess) // 2, n

    def matroid(self) -> Positroid[V]:
        """Return the matroid ``M_G`` of the graph, a positroid.

        Postnikov section 11: ``M_G`` is the set of source sets ``I_O``
        over all perfect orientations ``O``; it is a positroid whenever
        ``G`` is perfectly orientable (Proposition 11.7), which
        :meth:`research.positroid.Positroid.from_bases` re-verifies through
        Oh's theorem. Exponential; guarded.

        Raises:
            ValueError: If the graph is not perfectly orientable or is too
                large to enumerate.
        """
        bases = {self.source_set(o) for o in self.perfect_orientations()}
        if not bases:
            msg = "the graph is not perfectly orientable, so M_G is undefined"
            raise ValueError(msg)
        return Positroid.from_bases(self.boundary, bases)

    @functools.cached_property
    def is_bipartite(self) -> bool:
        """Whether every edge between internal vertices joins both colors."""
        return all(
            self._color[u] != self._color[w]
            for u, w in self.edges
            if u in self._color and w in self._color
        )

    def _require_matchable(self) -> None:
        """Raise unless almost perfect matchings and their boundary make sense."""
        if not self.is_bipartite:
            msg = "almost perfect matchings need a bipartite graph (Lemma 11.10)"
            raise ValueError(msg)
        if any(
            u in self._on_boundary and w in self._on_boundary for u, w in self.edges
        ):
            msg = (
                "an edge joins two boundary vertices, so the boundary of a "
                "matching is not defined"
            )
            raise ValueError(msg)

    def almost_perfect_matchings(self) -> tuple[frozenset[int], ...]:
        """Enumerate the almost perfect matchings, as sets of edge indices.

        Postnikov Lemma 11.10, for bipartite graphs: edge subsets covering
        every internal vertex exactly once; boundary vertices may or may
        not be covered. Exponential; guarded.

        Raises:
            ValueError: If the graph is not bipartite, has an edge between
                two boundary vertices, or is too large to enumerate.
        """
        self._require_matchable()
        self._guard_enumeration("almost perfect matchings")
        interior = list(self.internal_vertices)
        found: list[frozenset[int]] = []

        def extend(position: int, used: frozenset[V], chosen: frozenset[int]) -> None:
            if position == len(interior):
                found.append(chosen)
                return
            v = interior[position]
            if v in used:
                extend(position + 1, used, chosen)
                return
            for dart in self._rotation[v]:
                other = self._at[_mate(dart)]
                if other not in used:
                    extend(position + 1, used | {v, other}, chosen | {dart >> 1})

        extend(0, frozenset(), frozenset())
        return tuple(found)

    def matching_boundary(self, matching: Collection[int]) -> frozenset[V]:
        """Return the boundary ``dM`` of an almost perfect matching.

        Postnikov Lemma 11.10: for bipartite ``G`` the bases of ``M_G`` are
        the boundaries of the almost perfect matchings. Under the
        correspondence with perfect orientations (the matched edge at a
        black vertex is its one outgoing edge, at a white vertex its one
        incoming edge), ``b_i`` is a source exactly when it is matched to
        a white vertex or left unmatched next to a black one — the rule
        Lam states as ``I(Pi)`` (arXiv:1506.00603, section 4.1), with the
        boundary vertex colored opposite to its neighbor.

        Raises:
            ValueError: If the graph is not bipartite or has an edge
                between two boundary vertices.
        """
        self._require_matchable()
        used = {v for index in matching for v in self.edges[index]}
        return frozenset(
            b
            for b in self.boundary
            if (self._color[self._neighbor(b)] == WHITE) == (b in used)
        )

    # -------------------------------------------------------- transformations
    def to_decorated_permutation(self) -> DecoratedPermutation:
        """Return the decorated trip permutation ``pi_G^:``.

        Postnikov section 13; FWZ Definition 7.1.8: the trip permutation
        in Postnikov's direction, each fixed point decorated by the color
        of its lollipop. A white lollipop is Postnikov's ``col = -1``, the
        counted color this library calls clockwise (a coloop).

        Raises:
            ValueError: If some fixed point does not come from a lollipop
                — then the decoration is undefined (and the graph is not
                reduced, Postnikov Theorem 13.2 (4)).
        """
        white: set[int] = set()
        for index, b in enumerate(self.boundary, start=1):
            if self.trip_permutation[index - 1] != index:
                continue
            if b not in self.lollipops:
                msg = (
                    f"the fixed point {b!r} of the trip permutation does not "
                    f"come from a lollipop, so it has no decoration"
                )
                raise ValueError(msg)
            if self.lollipops[b] == WHITE:
                white.add(index)
        return DecoratedPermutation(self.trip_permutation, frozenset(white))

    def _require_reduced(self, statement: str) -> None:
        """Raise unless the graph is reduced, the hypothesis of ``statement``."""
        if not self.is_reduced():
            msg = f"{statement} holds for reduced plabic graphs only"
            raise ValueError(msg)

    def to_positroid(self) -> Positroid[V]:
        """Return the positroid of a reduced graph, read off its trips.

        Postnikov Proposition 16.4: for reduced ``G`` the Grassmann
        necklace of ``M_G`` is the necklace of ``pi_G^:``. The positroid
        class takes decorated permutations in the inverse
        (Ardila-Rincon-Williams) direction, so ``pi_G^:`` is inverted
        here. Agrees with :meth:`matroid` without enumerating
        orientations.

        Raises:
            ValueError: If the graph is not reduced.
        """
        self._require_reduced("reading the cell off the trip permutation")
        decorated = self.to_decorated_permutation().inverse()
        return Positroid.from_decorated_permutation(self.boundary, decorated)

    def to_grassmann_necklace(self) -> GrassmannNecklace[V]:
        """Return the Grassmann necklace of a reduced graph.

        Postnikov Proposition 16.4, through the same inversion as
        :meth:`to_positroid`; equal to :meth:`boundary_face_labels`
        (Oh-Postnikov-Speyer section 6).

        Raises:
            ValueError: If the graph is not reduced.
        """
        self._require_reduced("reading the necklace off the trip permutation")
        decorated = self.to_decorated_permutation().inverse()
        return GrassmannNecklace.from_decorated_permutation(self.boundary, decorated)

    def cyclic_shift(self, steps: int = 1) -> PlabicGraph[V]:
        """Return the same graph with the boundary labelling rotated.

        The new ``b_1`` is the old ``b_(1 + steps)``: the vertices keep
        their labels and the clockwise order is unchanged, only its
        starting point moves. The trip permutation on positions is
        conjugated accordingly.
        """
        n = len(self.boundary)
        if n == 0:
            return self
        cut = steps % n
        boundary = self.boundary[cut:] + self.boundary[:cut]
        return PlabicGraph._build(boundary, self.edges, self._rotation, self._color)

    # -------------------------------------------------------------- the moves
    def _draft(self) -> _Draft:
        """Return a mutable copy for composing moves."""
        return _Draft(self.boundary, self._rotation, self._color)

    def _corners(self, face: int) -> list[V]:
        """Return the vertices around an interior face, or raise."""
        if not 0 <= face < len(self._faces):
            msg = f"face index {face} is out of range 0..{len(self._faces) - 1}"
            raise ValueError(msg)
        if any(dart < 0 for dart in self._faces[face]):
            msg = (
                f"face {face} touches the boundary of the disk; a move needs "
                f"an internal face"
            )
            raise ValueError(msg)
        return [self._at[dart] for dart in self._faces[face]]

    def _is_alternating_square(self, corners: Sequence[V]) -> bool:
        """Whether ``corners`` are four distinct internal alternating vertices."""
        if len(corners) != _SQUARE_SIDES or len(set(corners)) != _SQUARE_SIDES:
            return False
        if any(v not in self._color for v in corners):
            return False
        shades = [self._color[v] for v in corners]
        return shades[0] != shades[1] and shades == [shades[0], shades[1]] * 2

    def square_move(self, face: int) -> PlabicGraph[V]:
        """Apply the square move (M1) at a face.

        Postnikov section 12, (M1): four trivalent vertices of alternating
        colors forming a square — switch all four colors. The rotation
        system is untouched. Preserves the trip permutation (Postnikov
        Lemma 13.1) and mutates the target face labels by exchanging the
        square's label ``Sac`` for ``Sbd`` (Oh-Postnikov-Speyer).

        Args:
            face: Index into :meth:`faces` of the square.

        Raises:
            ValueError: If the face is not an internal square on four
                distinct trivalent vertices of alternating colors.
        """
        corners = self._corners(face)
        if not self._is_alternating_square(corners) or any(
            self.degree(v) != _TRIVALENT for v in corners
        ):
            msg = (
                f"(M1) needs four trivalent vertices of alternating colors "
                f"forming a square; face {face} has corners {corners!r}"
            )
            raise ValueError(msg)
        flipped = {v: -c if v in corners else c for v, c in self.colors}
        return PlabicGraph._build(self.boundary, self.edges, self._rotation, flipped)

    def _internal_edge(self, edge: int) -> tuple[V, V]:
        """Return the endpoints of an edge between internal vertices, or raise."""
        if not 0 <= edge < len(self.edges):
            msg = f"edge index {edge} is out of range 0..{len(self.edges) - 1}"
            raise ValueError(msg)
        return self.edges[edge]

    def contract_edge(self, edge: int) -> PlabicGraph[V]:
        """Contract a unicolored edge (M2).

        Postnikov section 12, (M2) (FWZ number it (M3)): contract an edge
        joining two internal vertices of the same color. The merged vertex
        keeps the label of ``edges[edge][0]``.

        Raises:
            ValueError: If the edge is a loop or its endpoints are not two
                internal vertices of the same color.
        """
        u, w = self._internal_edge(edge)
        if (
            u == w
            or self._color.get(u) is None
            or self._color.get(u) != self._color.get(w)
        ):
            msg = (
                f"(M2) contracts an edge joining two distinct internal "
                f"vertices of the same color; edge {edge} joins {u!r} and {w!r}"
            )
            raise ValueError(msg)
        draft = self._draft()
        draft.contract(2 * edge)
        return _narrow(draft.finalize(), self)

    def uncontract_vertex(
        self, vertex: V, edges: Collection[int], label: V
    ) -> PlabicGraph[V]:
        """Split a vertex in two along a new unicolored edge (M2, reversed).

        Postnikov section 12, (M2): the inverse of contraction. The listed
        edges — which must be consecutive around ``vertex`` — move to a new
        vertex ``label`` of the same color, joined to ``vertex`` by a new
        edge.

        Args:
            vertex: The internal vertex to split.
            edges: Indices of the incident non-loop edges to move.
            label: The label of the new vertex.

        Raises:
            ValueError: If ``vertex`` is not internal, an index is not a
                non-loop edge at it, the edges are not consecutive in its
                rotation, or ``label`` is taken.
        """
        if vertex not in self._color:
            msg = f"(M2) splits an internal vertex; {vertex!r} is not one"
            raise ValueError(msg)
        rotation = self._rotation[vertex]
        moving = [d for d in rotation if d >> 1 in edges]
        loops = [d for d in moving if self._at[_mate(d)] == vertex]
        if len(moving) != len(set(edges)) or loops:
            msg = (
                f"(M2) moves non-loop edges incident to {vertex!r}; got "
                f"{sorted(edges)!r}"
            )
            raise ValueError(msg)
        block = _cyclic_block(rotation, moving)
        if block is None:
            msg = (
                f"(M2) needs the moved edges {sorted(edges)!r} to be "
                f"consecutive around {vertex!r}"
            )
            raise ValueError(msg)
        draft = self._draft()
        draft.uncontract(label, vertex, block)
        return _narrow(draft.finalize(), self)

    def remove_middle_vertex(self, vertex: V) -> PlabicGraph[V]:
        """Remove a degree-two internal vertex, gluing its edges (M3).

        Postnikov section 12, (M3) (FWZ number it (M2)).

        Raises:
            ValueError: If ``vertex`` is not an internal vertex of degree
                two on two distinct edges.
        """
        darts = self._rotation.get(vertex, ())
        if (
            vertex not in self._color
            or len(darts) != _MIDDLE_DEGREE
            or _mate(darts[0]) == darts[1]
        ):
            msg = (
                f"(M3) removes an internal vertex of degree two; {vertex!r} is not one"
            )
            raise ValueError(msg)
        draft = self._draft()
        draft.splice_out(vertex)
        return _narrow(draft.finalize(), self)

    def insert_middle_vertex(self, edge: int, color: int, label: V) -> PlabicGraph[V]:
        """Insert a degree-two vertex of either color into an edge (M3).

        Postnikov section 12, (M3), reversed.

        Raises:
            ValueError: If the edge index is out of range, the color is
                not 1 or -1, or ``label`` is taken.
        """
        self._internal_edge(edge)
        if color not in (BLACK, WHITE):
            msg = f"(M3) inserts a black (1) or white (-1) vertex; got {color!r}"
            raise ValueError(msg)
        draft = self._draft()
        first = draft.new_edge()
        second = draft.new_edge()
        draft.replace(2 * edge, first[0])
        draft.replace(2 * edge + 1, second[0])
        draft.add_vertex(label, [first[1], second[1]], color)
        return _narrow(draft.finalize(), self)

    def parallel_edge_reduction(self, face: int) -> PlabicGraph[V]:
        """Remove a hollow digon (R1).

        Postnikov section 12, (R1): two parallel edges joining two
        trivalent vertices of different colors bound a face; remove both
        vertices and glue the two hanging edges into one. Changes the trip
        permutation (Postnikov Lemma 13.1).

        Args:
            face: Index into :meth:`faces` of the digon.

        Raises:
            ValueError: If the face is not a digon on two trivalent
                vertices of different colors.
        """
        corners = self._corners(face)
        ok = (
            len(corners) == _DIGON_SIDES
            and all(v in self._color and self.degree(v) == _TRIVALENT for v in corners)
            and self._color[corners[0]] != self._color[corners[1]]
        )
        if not ok:
            msg = (
                f"(R1) removes two parallel edges joining trivalent vertices "
                f"of different colors; face {face} has corners {corners!r}"
            )
            raise ValueError(msg)
        digon = {dart >> 1 for dart in self._faces[face]}
        draft = self._draft()
        near, far = draft.new_edge()
        for v, fresh in zip(corners, (near, far), strict=True):
            (hanging,) = [d for d in self._rotation[v] if d >> 1 not in digon]
            draft.replace(_mate(hanging), fresh)
        for v in corners:
            draft.drop_vertex(v)
        return _narrow(draft.finalize(), self)

    def leaf_reduction(self, leaf: V, labels: Sequence[V]) -> PlabicGraph[V]:
        """Remove a leaf together with its neighbor (R2).

        Postnikov section 12, (R2): remove an internal leaf ``u`` whose
        neighbor ``v`` has the opposite color and degree at least 3; the
        other edges of ``v`` are disconnected into new endpoints of color
        ``col(u)``. Changes the trip permutation (Postnikov Lemma 13.1).

        Args:
            leaf: The degree-one internal vertex ``u``.
            labels: Labels for the new endpoints, one per other half-edge
                of ``v``, in counterclockwise order starting after ``u``.

        Raises:
            ValueError: If the hypotheses of (R2) fail, the labels do not
                match, or the result leaves the representable class (a
                piece cut off from the boundary — FWZ Definition 7.1.1).
        """
        if leaf not in self._color or self.degree(leaf) != 1:
            msg = f"(R2) removes a degree-one internal vertex; {leaf!r} is not one"
            raise ValueError(msg)
        v = self._neighbor(leaf)
        if (
            v not in self._color
            or self._color[v] == self._color[leaf]
            or self.degree(v) < _TRIVALENT
        ):
            msg = (
                f"(R2) needs the neighbor of {leaf!r} to be internal, of the "
                f"opposite color, and of degree at least 3; got {v!r}"
            )
            raise ValueError(msg)
        rotation = self._rotation[v]
        cut = rotation.index(_mate(self._rotation[leaf][0]))
        others = rotation[cut + 1 :] + rotation[:cut]
        if len(labels) != len(others):
            msg = (
                f"(R2) at {v!r} creates {len(others)} endpoints; got "
                f"{len(labels)} labels"
            )
            raise ValueError(msg)
        draft = self._draft()
        color = self._color[leaf]
        draft.rot[v] = [rotation[cut]]
        draft.drop_vertex(v)
        draft.drop_vertex(leaf)
        for label, dart in zip(labels, others, strict=True):
            draft.add_vertex(label, [dart], color)
        return _narrow(draft.finalize(), self)

    def normalized(self) -> PlabicGraph[V]:
        """Return the graph with all (M2) contractions and (M3) removals done.

        Contracts every unicolored edge between distinct internal vertices
        and removes every degree-two internal vertex until none is left;
        the result is bipartite apart from loops. Only moves (M2) and (M3)
        are used, so the trip permutation (Postnikov Lemma 13.1) and the
        face labels (Oh-Postnikov-Speyer) are unchanged.
        """
        draft = self._draft()
        draft.normalize()
        return _narrow(draft.finalize(), self)

    def _mutable_squares(self) -> list[int]:
        """Return the internal alternating square faces of a normalized graph."""
        return [
            index
            for index, face in enumerate(self._faces)
            if all(dart >= 0 for dart in face)
            and self._is_alternating_square([self._at[dart] for dart in face])
        ]

    def _mutated(self, face: int) -> PlabicGraph[Hashable]:
        """Square-move a face of a normalized graph, then renormalize.

        Corners of degree above three are first made trivalent by (M2)
        uncontraction of their outside edges; these are the "auxiliary
        moves" that accompany a square move.
        """
        draft = self._draft()
        for leaving in self._faces[face]:
            v = self._at[leaving]
            rotation = self._rotation[v]
            slot = rotation.index(leaving)
            outside = [
                rotation[(slot + offset) % len(rotation)]
                for offset in range(1, len(rotation) - 1)
            ]
            if len(outside) > 1:
                draft.uncontract(draft.fresh_label(), v, outside)
            draft.colors[v] = -draft.colors[v]
        draft.normalize()
        return draft.finalize()

    def square_move_class(
        self, *, limit: int = 2000
    ) -> tuple[PlabicGraph[Hashable], ...]:
        """Explore the reduced graphs reachable by square moves.

        Breadth-first search from the normalized graph, applying the
        square move (M1) with its auxiliary (M2)/(M3) moves at every
        internal square, and keeping one normalized graph per distinct set
        of target face labels. By Oh-Postnikov-Speyer Theorems 1.4 and 1.5
        the collections found are exactly the maximal weakly separated
        collections of the graph's positroid. The start graph comes first.

        Args:
            limit: Size guard on the number of graphs kept (page cost
                note: the classes grow fast).

        Raises:
            ValueError: If the graph is not reduced, or the class exceeds
                ``limit``.
        """
        self._require_reduced("exploring a move-equivalence class by face labels")
        start: PlabicGraph[Hashable] = _widen(self.normalized())
        found = {frozenset(start._face_labels("target")): start}
        frontier = [start]
        while frontier:
            graph = frontier.pop(0)
            for face in graph._mutable_squares():
                moved = graph._mutated(face)
                key = frozenset(moved._face_labels("target"))
                if key in found:
                    continue
                if len(found) >= limit:
                    msg = (
                        f"the square-move class exceeds the guard of {limit} "
                        f"graphs (page cost note); raise ``limit`` explicitly"
                    )
                    raise ValueError(msg)
                found[key] = moved
                frontier.append(moved)
        return tuple(found.values())

    # ---------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per half-edge.

        Columns: ``vertex`` (the label), ``kind`` (``boundary``, ``black``
        or ``white``), ``slot`` (the half-edge's index in the vertex's
        counterclockwise rotation), ``edge`` (the edge index) and ``end``
        (0 or 1: which end of the edge this is). Rows list the boundary
        vertices first, in clockwise order, then the internal vertices in
        stored order, so the boundary order survives a records-oriented
        JSON round trip through ``experiments.io.write_result``.

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        rows = [
            {
                "vertex": v,
                "kind": _KIND[self._color[v]] if v in self._color else "boundary",
                "slot": slot,
                "edge": dart >> 1,
                "end": dart & 1,
            }
            for v, darts in self.rotations
            for slot, dart in enumerate(darts)
        ]
        return pd.DataFrame(rows, columns=list(_COLUMNS))

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> PlabicGraph[Hashable]:
        """Rebuild a graph from a frame produced by :meth:`to_dataframe`.

        Re-validates the definition. Labels that a JSON round trip turned
        from tuples into lists are turned back, and a frame with no rows
        and no columns decodes to the empty graph.

        Raises:
            ValueError: If required columns are missing, a kind is
                unknown, or the decoded data fails validation.
        """
        if df.empty and len(df.columns) == 0:
            # Records-oriented JSON of the empty graph has no columns.
            return PlabicGraph._build((), (), {}, {})
        missing = set(_COLUMNS) - set(df.columns)
        if missing:
            msg = f"dataframe is missing required columns {sorted(missing)}"
            raise ValueError(msg)
        boundary: list[Hashable] = []
        colors: dict[Hashable, int] = {}
        slots: dict[Hashable, list[tuple[int, int]]] = {}
        ends: dict[int, Hashable] = {}
        for raw, kind, slot, edge, end in zip(
            df["vertex"], df["kind"], df["slot"], df["edge"], df["end"], strict=True
        ):
            v = _hashable(raw)
            if kind == "boundary":
                if v not in slots:
                    boundary.append(v)
            elif kind in _COLOR:
                colors[v] = _COLOR[kind]
            else:
                msg = f"unknown vertex kind {kind!r}; expected boundary/black/white"
                raise ValueError(msg)
            dart = 2 * int(edge) + int(end)
            slots.setdefault(v, []).append((int(slot), dart))
            ends[dart] = v
        count = len(ends) // 2
        if sorted(ends) != list(range(2 * count)):
            msg = "(P1) violated: every edge needs exactly its two half-edge rows"
            raise ValueError(msg)
        edges = [(ends[2 * e], ends[2 * e + 1]) for e in range(count)]
        rotations = {v: [d for _, d in sorted(pairs)] for v, pairs in slots.items()}
        return PlabicGraph._build(boundary, edges, rotations, colors)

    # ---------------------------------------------------------- visualization
    def _layout(self) -> dict[V, tuple[float, float]]:
        """Place the boundary on a circle and relax the rest to barycenters.

        A Tutte-style embedding from the rotation system alone (the page's
        visualization note): boundary vertices clockwise from the top,
        each internal vertex at the mean of its neighbors, and lollipop
        heads and leaves pulled in from the vertex they hang off.
        """
        circle = unit_circle(len(self.boundary), phase=math.pi / 2, clockwise=True)
        place: dict[V, tuple[float, float]] = dict(
            zip(self.boundary, circle, strict=True)
        )
        for v in self.internal_vertices:
            place[v] = (0.0, 0.0)
        for _ in range(_RELAXATION_SWEEPS):
            for v in self.internal_vertices:
                around = [place[self._at[_mate(d)]] for d in self._rotation[v]]
                place[v] = (
                    sum(x for x, _ in around) / len(around),
                    sum(y for _, y in around) / len(around),
                )
        for v in self.internal_vertices:
            if self.degree(v) == 1:
                x, y = place[self._neighbor(v)]
                place[v] = (0.78 * x + 0.04, 0.78 * y + 0.04)
        return place

    def _edge_curves(
        self, place: Mapping[V, tuple[float, float]]
    ) -> list[list[tuple[float, float]]]:
        """Return a polyline per edge, bending parallel edges and loops apart."""
        groups: dict[frozenset[V], list[int]] = {}
        for index, (u, w) in enumerate(self.edges):
            groups.setdefault(frozenset({u, w}), []).append(index)
        curves: list[list[tuple[float, float]]] = [[] for _ in self.edges]
        samples = [step / 16 for step in range(17)]
        for members in groups.values():
            for rank, index in enumerate(members):
                u, w = self.edges[index]
                (x0, y0), (x1, y1) = place[u], place[w]
                if u == w:
                    radius = 0.08 * (rank + 1)
                    curves[index] = [
                        (
                            x0 + radius * (1 - math.cos(2 * math.pi * t)),
                            y0 + radius * math.sin(2 * math.pi * t),
                        )
                        for t in samples
                    ]
                    continue
                bend = 0.35 * (rank - (len(members) - 1) / 2)
                if index != members[0] and self.edges[members[0]][0] != u:
                    bend = -bend
                cx = (x0 + x1) / 2 - bend * (y1 - y0)
                cy = (y0 + y1) / 2 + bend * (x1 - x0)
                curves[index] = [
                    (
                        (1 - t) ** 2 * x0 + 2 * t * (1 - t) * cx + t**2 * x1,
                        (1 - t) ** 2 * y0 + 2 * t * (1 - t) * cy + t**2 * y1,
                    )
                    for t in samples
                ]
        return curves

    def plot_graph(self, ax: Axes | None = None) -> Axes:
        """Draw the graph in the disk: filled black and open white vertices.

        The layout is computed from the rotation system (boundary on a
        circle, clockwise from the top; internal vertices at barycenters),
        so it is a schematic rather than a guaranteed crossing-free
        drawing. Draws onto ``ax`` or a fresh figure; never calls ``show``.
        """
        ax = ensure_axes(ax)
        place = self._layout()
        ring = [
            (math.cos(2 * math.pi * t / 120), math.sin(2 * math.pi * t / 120))
            for t in range(121)
        ]
        ax.plot([x for x, _ in ring], [y for _, y in ring], color="tab:gray", lw=0.8)
        for curve in self._edge_curves(place):
            ax.plot([x for x, _ in curve], [y for _, y in curve], color="black", lw=1.0)
        for v in self.internal_vertices:
            ax.scatter(
                [place[v][0]],
                [place[v][1]],
                facecolors="black" if self._color[v] == BLACK else "white",
                edgecolors="black",
                zorder=3,
            )
        for index, b in enumerate(self.boundary, start=1):
            x, y = place[b]
            ax.scatter([x], [y], color="tab:gray", s=12, zorder=3)
            ax.annotate(
                str(b) if b == index else f"{index}:{b}",
                (1.12 * x, 1.12 * y),
                ha="center",
                va="center",
                fontsize=8,
            )
        ax.set_aspect("equal")
        ax.set_axis_off()
        return ax

    def plot_trips(self, ax: Axes | None = None) -> Axes:
        """Draw the graph with its one-way trips overlaid as strands.

        Each trip is drawn through the midpoints of the edges it follows,
        with an arrow into the boundary vertex where it ends — the
        alternating strand diagram dual to the graph (Postnikov section
        14; Scott's "Postnikov arrangement"). Draws onto ``ax`` or a fresh
        figure; never calls ``show``.
        """
        ax = self.plot_graph(ax)
        place = self._layout()
        curves = self._edge_curves(place)
        for t, trip in enumerate(self._trips[: len(self.boundary)]):
            start = place[self.boundary[t]]
            finish = place[self.boundary[self.trip_permutation[t] - 1]]
            middle = [curves[d >> 1][len(curves[d >> 1]) // 2] for d in trip]
            path = [start, *middle, finish]
            color = f"C{t % 10}"
            ax.plot(
                [x for x, _ in path[:-1]],
                [y for _, y in path[:-1]],
                color=color,
                lw=1.4,
                alpha=0.8,
                zorder=2,
            )
            ax.annotate(
                "",
                xy=path[-1],
                xytext=path[-2],
                arrowprops={"arrowstyle": "-|>", "color": color, "alpha": 0.8},
            )
        return ax


# --------------------------------------------------------------------------- #
# Construction helpers
# --------------------------------------------------------------------------- #
def _narrow[V: Hashable](
    graph: PlabicGraph[Hashable], like: PlabicGraph[V]
) -> PlabicGraph[V]:
    """Re-type a draft's result whose labels all have the type of ``like``'s."""
    del like
    return graph  # type: ignore[return-value]  # labels were supplied as V by the caller


def _widen[V: Hashable](graph: PlabicGraph[V]) -> PlabicGraph[Hashable]:
    """Forget the label type (the dataclass is invariant in ``V``)."""
    return PlabicGraph(
        tuple(graph.boundary),
        tuple((u, w) for u, w in graph.edges),
        tuple((v, darts) for v, darts in graph.rotations),
        tuple((v, color) for v, color in graph.colors),
    )


def _cyclic_block(rotation: Sequence[int], moving: Sequence[int]) -> list[int] | None:
    """Return ``moving`` as a consecutive cyclic run of ``rotation``, if it is one."""
    if not moving or len(moving) >= len(rotation):
        return None
    wanted = set(moving)
    for start, dart in enumerate(rotation):
        if dart in wanted and rotation[start - 1] not in wanted:
            run = [rotation[(start + k) % len(rotation)] for k in range(len(moving))]
            return run if set(run) == wanted else None
    return None


def _geometric_rotations[V: Hashable](
    positions: Mapping[V, _Point], edges: Sequence[tuple[V, V]]
) -> dict[V, list[int]]:
    """Read the counterclockwise dart order at each vertex off a drawing."""
    incident: dict[V, list[tuple[int, _Point]]] = {v: [] for v in positions}
    for index, (u, w) in enumerate(edges):
        incident[u].append((2 * index, positions[w]))
        incident[w].append((2 * index + 1, positions[u]))
    return {v: _sorted_ccw(positions[v], darts) for v, darts in incident.items()}


def _split_crossing(draft: _Draft, v: Hashable) -> None:
    """Split a 2-in/2-out vertex into a black and a white trivalent vertex.

    Darts at even ids leave along the orientation and odd ones arrive. The
    two arriving darts must be neighbors in the rotation; they stay on
    ``v`` (black: its one outgoing edge is the new one) and the two
    leaving darts move to ``(v, "white")``.
    """
    rotation = draft.rot[v]
    for start in range(len(rotation)):
        ring = rotation[start:] + rotation[:start]
        if ring[0] & 1 and ring[1] & 1:
            draft.colors[v] = BLACK
            draft.uncontract((v, "white"), v, ring[2:])
            draft.colors[(v, "white")] = WHITE
            return
    msg = f"the two incoming edges at {v!r} are not adjacent; cannot split it"
    raise ValueError(msg)


def _from_oriented[W: Hashable](
    network: PlanarNetwork[W], *, split: bool
) -> PlabicGraph[Hashable]:
    """Build the plabic graph of a (Gamma- or perfectly) oriented network."""
    edges = [(tail, head) for tail, head, _ in network.edges]
    rotations = _geometric_rotations(dict(network.positions), edges)
    draft = _Draft(network.boundary, rotations, {})
    for v in network.internal_vertices:
        arriving = sum(dart & 1 for dart in rotations[v])
        leaving = len(rotations[v]) - arriving
        if arriving == 1:
            draft.colors[v] = WHITE
        elif leaving == 1:
            draft.colors[v] = BLACK
        elif split and arriving == leaving == _DIGON_SIDES:
            _split_crossing(draft, v)
        else:
            msg = (
                f"internal vertex {v!r} has {leaving} outgoing and {arriving} "
                f"incoming edges; a perfect orientation needs exactly one "
                f"outgoing (black) or one incoming (white) (Postnikov "
                f"Definition 11.5)"
            )
            raise ValueError(msg)
    for b in network.boundary:
        if not rotations[b]:
            near, far = draft.new_edge()
            draft.rot[b] = [near]
            draft.owner[near] = b
            color = WHITE if b in network.source_set else BLACK
            draft.add_vertex((b, "lollipop"), [far], color)
    return draft.finalize()


# --------------------------------------------------------------------------- #
# Weak separation
# --------------------------------------------------------------------------- #
def is_weakly_separated(
    first: Collection[int], second: Collection[int], n: int
) -> bool:
    r"""Return whether two equal-size subsets of ``[n]`` are weakly separated.

    Oh-Postnikov-Speyer Definition 3.1 (following Scott 2005): ``I`` and
    ``J`` are weakly separated if there are no cyclically ordered
    ``a, b, a', b'`` with ``a, a'`` in ``I \ J`` and ``b, b'`` in
    ``J \ I`` — equivalently a chord separates ``I \ J`` from ``J \ I``.
    Read around the circle, the elements of the symmetric difference then
    switch sides at most twice.

    Raises:
        ValueError: If the sets are not equal-size subsets of ``[n]`` (the
            unequal-size notion of Leclerc-Zelevinsky is a different
            definition).
    """
    a, b = frozenset(first), frozenset(second)
    ground = frozenset(range(1, n + 1))
    if len(a) != len(b) or not a <= ground or not b <= ground:
        msg = (
            f"weak separation (Oh-Postnikov-Speyer Definition 3.1) compares "
            f"equal-size subsets of [{n}]; got {sorted(a)!r} and {sorted(b)!r}"
        )
        raise ValueError(msg)
    sides = [i in a for i in range(1, n + 1) if (i in a) != (i in b)]
    changes = sum(1 for k, side in enumerate(sides) if side != sides[k - 1])
    return changes <= _WEAK_SEPARATION_MAX_CHANGES


def is_weakly_separated_collection(
    collection: Iterable[Collection[int]], n: int
) -> bool:
    """Return whether a family of subsets is pairwise weakly separated.

    The page's definition of a weakly separated collection
    (Oh-Postnikov-Speyer Definition 3.1).
    """
    return all(
        is_weakly_separated(a, b, n) for a, b in itertools.combinations(collection, 2)
    )


def maximal_weakly_separated_collections(
    n: int,
    candidates: Iterable[Collection[int]],
    *,
    containing: Iterable[Collection[int]] = (),
) -> tuple[frozenset[frozenset[int]], ...]:
    """Enumerate the maximal weakly separated collections among ``candidates``.

    Bron-Kerbosch clique enumeration (with pivoting) over the
    weak-separation graph — the shortcut the page used to verify its
    counting fixtures. With ``candidates`` the bases of a positroid and
    ``containing`` its Grassmann necklace, these are the collections of
    Oh-Postnikov-Speyer Theorem 1.5; by purity (their Theorem 1.3) all have
    the same size.

    Args:
        n: The ground set is ``[n]``.
        candidates: The subsets a collection may use.
        containing: Subsets every collection must contain; candidates not
            weakly separated from all of them are discarded.

    Raises:
        ValueError: If there are too many candidates to enumerate (page
            cost note), or ``containing`` is not itself weakly separated.
    """
    forced = frozenset(frozenset(entry) for entry in containing)
    if not is_weakly_separated_collection(forced, n):
        msg = "the required subsets are not pairwise weakly separated"
        raise ValueError(msg)
    pool = sorted({frozenset(c) for c in candidates} - forced, key=sorted)
    pool = [c for c in pool if all(is_weakly_separated(c, f, n) for f in forced)]
    if len(pool) > _MAX_CLIQUE_CANDIDATES:
        msg = (
            f"{len(pool)} candidate subsets exceed the enumeration guard of "
            f"{_MAX_CLIQUE_CANDIDATES} (page cost note)"
        )
        raise ValueError(msg)
    friends = {
        c: {d for d in pool if d != c and is_weakly_separated(c, d, n)} for c in pool
    }
    found: list[frozenset[frozenset[int]]] = []

    def expand(
        clique: frozenset[frozenset[int]],
        open_: set[frozenset[int]],
        closed: set[frozenset[int]],
    ) -> None:
        if not open_ and not closed:
            found.append(clique | forced)
            return
        pivot = max(open_ | closed, key=lambda c: len(friends[c] & open_))
        for c in sorted(open_ - friends[pivot], key=sorted):
            expand(clique | {c}, open_ & friends[c], closed & friends[c])
            open_ = open_ - {c}
            closed = closed | {c}

    expand(frozenset(), set(pool), set())
    return tuple(found)


# --------------------------------------------------------------------------- #
# Canonical example constructors — the page's test fixtures
# --------------------------------------------------------------------------- #
def lollipop_graph(colors: Sequence[int]) -> PlabicGraph[int | tuple[int, str]]:
    """Return the graph consisting only of ``n`` lollipops.

    The page's lollipop fixture: reduced, with the identity trip
    permutation decorated by the lollipop colors, a single face
    (dimension ``F - 1 = 0``), each black lollipop a loop and each white
    one a coloop of the positroid. ``colors[i - 1]`` is the color of the
    lollipop at ``b_i``: 1 black, -1 white.
    """
    n = len(colors)
    heads: list[int | tuple[int, str]] = [(i, "lollipop") for i in range(1, n + 1)]
    return PlabicGraph.from_rotation_system(
        tuple(range(1, n + 1)),
        [(i, heads[i - 1]) for i in range(1, n + 1)],
        {head: [index] for index, head in enumerate(heads)},
        dict(zip(heads, colors, strict=True)),
    )


def gr24_square_pair() -> tuple[PlabicGraph[int | str], PlabicGraph[int | str]]:
    """Return the two reduced graphs of the top cell of ``Gr(2, 4)``.

    The page's square-move pair: a square ``T, R, B, L`` with a leg to
    each of ``b_1..b_4`` (top, right, bottom, left). The first graph has
    target face labels ``{12, 23, 34, 14, 13}``, the second — its square
    move — ``{12, 23, 34, 14, 24}``: the mutation ``13 <-> 24`` realizing
    ``p13 p24 = p12 p34 + p14 p23``.
    """
    edges: list[tuple[int | str, int | str]] = [
        (1, "T"),
        (2, "R"),
        (3, "B"),
        (4, "L"),
        ("L", "T"),
        ("T", "R"),
        ("R", "B"),
        ("B", "L"),
    ]
    rotations: dict[int | str, list[int]] = {
        "T": [0, 4, 5],
        "R": [1, 5, 6],
        "B": [6, 7, 2],
        "L": [4, 3, 7],
    }
    boundary: tuple[int | str, ...] = (1, 2, 3, 4)
    colors: dict[int | str, int] = {"T": BLACK, "R": WHITE, "B": BLACK, "L": WHITE}
    first = PlabicGraph.from_rotation_system(boundary, edges, rotations, colors)
    return first, first.square_move(first.face_count - 1)


def top_cell_graph(k: int, n: int) -> PlabicGraph[Hashable]:
    """Return the Le-graph of the top cell of ``Gr(k, n)``.

    The Le-diagram of the full ``k x (n - k)`` rectangle filled with 1's
    (Postnikov section 20): a reduced graph with trip permutation
    ``i -> i + k`` and ``k(n - k) + 1`` faces (FWZ Corollary 7.10.7).

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    if not 0 <= k <= n:
        msg = f"the top cell of Gr(k, n) needs 0 <= k <= n; got k={k}, n={n}"
        raise ValueError(msg)
    return PlabicGraph.from_le_diagram([[1] * (n - k) for _ in range(k)], n)


def nonreduced_wiring_example() -> PlabicGraph[int | tuple[int, int]]:
    """Return the plabic graph of the non-reduced word ``s_1 s_1`` in ``S_2``.

    The page's negative fixture for the reducedness criterion: the graph
    satisfies every plabic-graph axiom but is not reduced, because the two
    wires cross twice (FWZ Exercise 7.3.9, Remarks 7.3.10 and 7.8.6).
    """
    return PlabicGraph.from_wiring_diagram([1, 1], 2)


def hollow_digon() -> PlabicGraph[int | str]:
    """Return the hollow digon — the smallest non-reduced fixture.

    Two trivalent internal vertices of opposite colors (``u`` black,
    ``w`` white) joined by two parallel edges, with legs to ``b_1`` and
    ``b_2``: a legitimate plabic graph to which (R1) applies as it stands.
    """
    return PlabicGraph.from_rotation_system(
        (1, 2),
        [(1, "u"), ("u", "w"), ("u", "w"), ("w", 2)],
        {"u": [1, 0, 2], "w": [3, 1, 2]},
        {"u": BLACK, "w": WHITE},
    )
