"""Tests for research.plabic_graph, from the Logseq Plabic Graph page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Postnikov Lemma 13.1,
Theorems 13.2, 13.4, 12.7 and 17.1, Corollary 14.7, Propositions 11.7 and 16.4,
Lemma 11.10, Corollary 20.1; the matching-sum formula of Postnikov, refined
by Talaska and Speyer (Williams ICM Theorem 2.17); Oh-Postnikov-Speyer
Theorems 1.3, 1.4, 1.5 and 6.8; Fomin-Williams-Zelevinsky Theorems 7.11.5
and 7.13.2); round-trip laws come from the API contract; and every numbered
axiom has a rejection test naming it.

Direction discipline: the page states the trip permutation and its
statistics in Postnikov's direction, and ``to_decorated_permutation`` returns
it literally. Every necklace read off it here uses
``NecklaceConvention.POSTNIKOV`` on that literal permutation, so the tests
do not depend on the inversion ``to_positroid`` performs internally.

Three sources of reduced graphs drive the property tests: Le-graphs of
Le-diagrams (closed up from random fillings), arbitrary members of the
square-move class of a large cell, and wiring-diagram graphs of reduced
words, all then scrambled by random moves (M1)-(M3) and rotations of the
boundary labelling.

Move-equivalence is decided through the page's Oh-Postnikov-Speyer block:
graphs with equal target face-label sets are strongly equivalent, and the
square-move class of a graph lists every label set reachable from it.

Non-reduced graphs come from non-reduced words and from splicing a hollow
digon, a hollow monogon, or a triple edge into a reduced graph (the digon
ones scrambled again by moves) — (R1) read backwards — so that
Definition 12.5 / FWZ Definition 7.1.6 is an oracle independent of the trip
criterion being tested.

Not transcribed, because the structure they quantify over is not implemented
(spec, transformation backlog): Lemma 13.6 in general (no reduction search —
certified on the two reduction fixtures), Theorem 12.1 and the injectivity
half of Theorem 12.7 (no weighted moves; the image half is tested through
weighted Le-networks of research.boundary_measurement, and the weighted
matching formula on its square network and on weighted Le-graphs, both with
multi-term sums), Corollary 16.5 (no weighted non-reduced source), statement (1) of the
cluster-structure theorem and the finite-type counts 16/42/128 (no
cluster-algebra module), the weight formulas of (M1), (R3), and T-duality.
"""

import functools
import itertools
import math
from collections.abc import Hashable, Iterable, Mapping, Sequence
from fractions import Fraction
from typing import cast

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from experiments import io
from research import boundary_measurement as bm
from research import decorated_permutation as dp
from research import plabic_graph as pg
from research import positroid as ps

matplotlib.use("Agg")

POSTNIKOV = dp.NecklaceConvention.POSTNIKOV
BLACK, WHITE = pg.BLACK, pg.WHITE
Graph = pg.PlabicGraph
type AnyGraph = pg.PlabicGraph[Hashable]


# --------------------------------------------------------------------------- #
# Helpers and strategies
# --------------------------------------------------------------------------- #
def _subsets(*labels: str) -> frozenset[frozenset[int]]:
    """Parse labels written as on the page: ``"13"`` is ``{1, 3}``."""
    return frozenset(frozenset(int(c) for c in label) for label in labels)


def _le_close(filling: list[list[int]]) -> list[list[int]]:
    """Set to 1 every 0 with a 1 above it and a 1 to its left, until stable."""
    rows = [list(row) for row in filling]
    changed = True
    while changed:
        changed = False
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                above = any(len(rows[r]) > j and rows[r][j] for r in range(i))
                left = any(row[c] for c in range(j))
                if not value and above and left:
                    row[j] = 1
                    changed = True
    return rows


def _all_le_diagrams(k: int, n: int) -> list[list[list[int]]]:
    """Enumerate every Le-diagram of type ``(k, n)`` by brute force."""
    found = []
    for shape in itertools.product(range(n - k, -1, -1), repeat=k):
        if any(a < b for a, b in itertools.pairwise(shape)):
            continue
        cells = [(i, j) for i, width in enumerate(shape) for j in range(width)]
        for bits in itertools.product((0, 1), repeat=len(cells)):
            filling = [[0] * width for width in shape]
            for (i, j), bit in zip(cells, bits, strict=True):
                filling[i][j] = bit
            if _le_close(filling) == filling:
                found.append(filling)
    return found


def _any[V: Hashable](graph: pg.PlabicGraph[V]) -> AnyGraph:
    """Forget the label type (the class is invariant in it)."""
    return cast("AnyGraph", graph)


def _graph(
    boundary: Sequence[Hashable],
    edges: Sequence[tuple[Hashable, Hashable]],
    rotations: Mapping[Hashable, Sequence[int]],
    colors: Mapping[Hashable, int],
) -> AnyGraph:
    """Build through ``from_rotation_system`` with mixed label types."""
    return Graph.from_rotation_system(boundary, edges, rotations, colors)


def _le_graphs(n: int) -> list[AnyGraph]:
    """Return the Le-graph of every Le-diagram on ``n`` boundary vertices."""
    return [
        Graph.from_le_diagram(filling, n)
        for k in range(n + 1)
        for filling in _all_le_diagrams(k, n)
    ]


@st.composite
def le_fillings(draw: st.DrawFn, max_n: int = 6) -> tuple[list[list[int]], int]:
    n = draw(st.integers(1, max_n))
    k = draw(st.integers(0, n))
    shape = sorted(
        draw(st.lists(st.integers(0, n - k), min_size=k, max_size=k)), reverse=True
    )
    raw = [draw(st.lists(st.integers(0, 1), min_size=w, max_size=w)) for w in shape]
    return _le_close(raw), n


@st.composite
def words(
    draw: st.DrawFn, max_n: int = 4, max_length: int = 5
) -> tuple[list[int], int]:
    n = draw(st.integers(2, max_n))
    return draw(st.lists(st.integers(1, n - 1), max_size=max_length)), n


def _permutation_of_word(word: Sequence[int], n: int) -> list[int]:
    wires = list(range(n))
    for letter in word:
        wires[letter - 1], wires[letter] = wires[letter], wires[letter - 1]
    return wires


def _is_reduced_word(word: Sequence[int], n: int) -> bool:
    """A word is reduced iff its length is the inversion number of its product."""
    wires = _permutation_of_word(word, n)
    inversions = sum(1 for a, b in itertools.combinations(wires, 2) if a > b)
    return inversions == len(word)


def _fresh(graph: AnyGraph, stem: str) -> tuple[str, int]:
    """Return a vertex label not yet used by ``graph``."""
    taken = {v for v, _ in graph.rotations}
    return next((stem, k) for k in itertools.count() if (stem, k) not in taken)


def _mutable_squares(graph: AnyGraph) -> list[int]:
    """Return the internal faces on four distinct, alternately colored vertices."""
    internal = set(graph.internal_vertices)
    found = []
    for face, darts in enumerate(graph.faces()):
        corners = [graph.edges[d >> 1][d & 1] for d in darts]
        if len(corners) != 4 or len(set(corners)) != 4 or not set(corners) <= internal:
            continue
        if any((v, v) in graph.edges for v in corners):
            continue  # a loop at a corner cannot be moved by ``uncontract_vertex``
        shades = [graph.color(v) for v in corners]
        if shades[0] != shades[1] and shades == [shades[0], shades[1]] * 2:
            found.append(face)
    return found


def _face_on(graph: AnyGraph, corners: Sequence[Hashable]) -> int:
    """Return the index of the four-sided face whose corners are ``corners``."""
    (face,) = [
        f
        for f, darts in enumerate(graph.faces())
        if len(darts) == 4
        and {graph.edges[d >> 1][d & 1] for d in darts} == set(corners)
    ]
    return face


def _expose(graph: AnyGraph, face: int) -> tuple[AnyGraph, int]:
    """Make the corners of a square trivalent by (M2); return graph and face."""
    corners = [graph.edges[d >> 1][d & 1] for d in graph.faces()[face]]
    for v in corners:
        sides = {d >> 1 for d in graph.faces()[_face_on(graph, corners)]}
        outside = [d >> 1 for d in dict(graph.rotations)[v] if d >> 1 not in sides]
        if len(outside) > 1:
            graph = graph.uncontract_vertex(v, outside, _fresh(graph, "leg"))
    return graph, _face_on(graph, corners)


def _scramble(draw: st.DrawFn, graph: AnyGraph, steps: int = 3) -> AnyGraph:
    """Apply a random sequence of moves (M1)-(M3) and a boundary rotation.

    Middle vertices are never inserted next to a degree-one vertex: that
    move is legal, but it turns a lollipop into a leaf and so leaves the
    leafless class the criteria are stated for.
    """
    for _ in range(draw(st.integers(0, steps))):
        kind = draw(st.sampled_from(["square", "insert", "split", "normalize"]))
        if kind == "square":
            squares = _mutable_squares(graph)
            if squares:
                exposed, face = _expose(graph, draw(st.sampled_from(squares)))
                if _square_move_applies(exposed, face):
                    graph = exposed.square_move(face)
        elif kind == "insert":
            edges = [
                index
                for index, edge in enumerate(graph.edges)
                if all(graph.degree(v) > 1 or v in graph.boundary for v in edge)
            ]
            if edges:
                color = draw(st.sampled_from([BLACK, WHITE]))
                graph = graph.insert_middle_vertex(
                    draw(st.sampled_from(edges)), color, _fresh(graph, "new")
                )
        elif kind == "split":
            big = [v for v in graph.internal_vertices if graph.degree(v) >= 3]
            if big:
                graph = _split_somewhere(draw, graph, draw(st.sampled_from(big)))
        elif kind == "normalize":
            graph = graph.normalized()
    if graph.size:
        graph = graph.cyclic_shift(draw(st.integers(0, graph.size - 1)))
    return graph


def _square_move_applies(graph: AnyGraph, face: int) -> bool:
    try:
        graph.square_move(face)
    except ValueError:
        return False
    return True


def _split_somewhere(draw: st.DrawFn, graph: AnyGraph, vertex: Hashable) -> AnyGraph:
    darts = dict(graph.rotations)[vertex]
    if len({d >> 1 for d in darts}) < len(darts):
        return graph  # a loop sits at the vertex; leave it alone
    start = draw(st.integers(0, len(darts) - 1))
    size = draw(st.integers(1, len(darts) - 1))
    block = [darts[(start + offset) % len(darts)] >> 1 for offset in range(size)]
    return graph.uncontract_vertex(vertex, block, _fresh(graph, "split"))


def _edge_rotations(graph: AnyGraph) -> dict[Hashable, list[int]]:
    """Return the internal rotations as edge indices, for rebuilding a graph."""
    return {
        v: [d >> 1 for d in darts]
        for v, darts in graph.rotations
        if v not in graph.boundary
    }


def _digon_edges(graph: AnyGraph) -> list[int]:
    """Return the non-loop edges away from every degree-one internal vertex."""
    return [
        index
        for index, (x, y) in enumerate(graph.edges)
        if x != y and all(graph.degree(v) > 1 or v in graph.boundary for v in (x, y))
    ]


def _insert_digon(graph: AnyGraph, edge: int, *, flip: bool = False) -> AnyGraph:
    """Splice a hollow digon ``du = dw`` into an edge — (R1) read backwards."""
    x, y = graph.edges[edge]
    count = len(graph.edges)
    edges = list(graph.edges)
    edges[edge] = (x, "du")
    edges += [("du", "dw"), ("du", "dw"), ("dw", y)]
    rotations = _edge_rotations(graph)
    if y not in graph.boundary:
        rotations[y] = [
            count + 2 if d == 2 * edge + 1 else d >> 1 for d in dict(graph.rotations)[y]
        ]
    rotations["du"] = [count, edge, count + 1]
    rotations["dw"] = [count + 2, count, count + 1]
    colors = dict(graph.colors)
    colors["du"], colors["dw"] = (WHITE, BLACK) if flip else (BLACK, WHITE)
    return Graph.from_rotation_system(graph.boundary, edges, rotations, colors)


def _insert_monogon(graph: AnyGraph, edge: int, color: int) -> AnyGraph:
    """Splice a vertex carrying a loop — a hollow monogon — into an edge."""
    x, y = graph.edges[edge]
    count = len(graph.edges)
    edges = list(graph.edges)
    edges[edge] = (x, "dm")
    edges += [("dm", "dm"), ("dm", y)]
    rotations = _edge_rotations(graph)
    if y not in graph.boundary:
        rotations[y] = [
            count + 1 if d == 2 * edge + 1 else d >> 1 for d in dict(graph.rotations)[y]
        ]
    rotations["dm"] = [edge, count, count, count + 1]
    colors = dict(graph.colors)
    colors["dm"] = color
    return Graph.from_rotation_system(graph.boundary, edges, rotations, colors)


def _insert_triple(graph: AnyGraph, edge: int, *, flip: bool = False) -> AnyGraph:
    """Splice in a vertex ``ta`` joined to a trivalent ``tb`` by three edges.

    Uncontracting two of the three parallel edges off ``ta`` exposes a
    hollow digon, so the result is non-reduced by Definition 12.5; among
    the trip conditions it is condition (2) that sees it.
    """
    x, y = graph.edges[edge]
    count = len(graph.edges)
    edges = list(graph.edges)
    edges[edge] = (x, "ta")
    edges += [("ta", y), ("ta", "tb"), ("ta", "tb"), ("ta", "tb")]
    rotations = _edge_rotations(graph)
    if y not in graph.boundary:
        rotations[y] = [
            count if d == 2 * edge + 1 else d >> 1 for d in dict(graph.rotations)[y]
        ]
    rotations["ta"] = [edge, count, count + 1, count + 2, count + 3]
    rotations["tb"] = [count + 1, count + 3, count + 2]
    colors = dict(graph.colors)
    colors["ta"], colors["tb"] = (WHITE, BLACK) if flip else (BLACK, WHITE)
    return Graph.from_rotation_system(graph.boundary, edges, rotations, colors)


def _triple_edge_fixture() -> AnyGraph:
    """A black vertex on both legs, joined to a white vertex by three edges."""
    return _graph(
        (1, 2),
        [(1, "a"), (2, "a"), ("a", "b"), ("a", "b"), ("a", "b")],
        {"a": [0, 1, 2, 3, 4], "b": [2, 4, 3]},
        {"a": BLACK, "b": WHITE},
    )


def _face_on_digon(graph: AnyGraph) -> int:
    """Return the face bounded by the two parallel edges of a spliced digon."""
    (face,) = [
        f
        for f, darts in enumerate(graph.faces())
        if {graph.edges[d >> 1][d & 1] for d in darts} == {"du", "dw"}
    ]
    return face


def _attach_leaf(graph: AnyGraph, vertex: Hashable, slot: int) -> AnyGraph:
    """Hang a leaf of the opposite color off ``vertex`` — (R2) read backwards."""
    count = len(graph.edges)
    rotations = _edge_rotations(graph)
    rotations[vertex] = [*rotations[vertex][:slot], count, *rotations[vertex][slot:]]
    rotations["leaf"] = [count]
    colors = dict(graph.colors)
    colors["leaf"] = -graph.color(vertex)
    return Graph.from_rotation_system(
        graph.boundary, [*graph.edges, (vertex, "leaf")], rotations, colors
    )


def _neighbor_rotations(graph: AnyGraph) -> dict[Hashable, tuple[str, ...]]:
    """Return each vertex's cyclic neighbor sequence — blind to edge numbering."""
    found: dict[Hashable, tuple[str, ...]] = {}
    for v, darts in graph.rotations:
        around = [repr(graph.edges[d >> 1][1 - (d & 1)]) for d in darts]
        turns = [tuple(around[k:] + around[:k]) for k in range(len(around))]
        found[v] = min(turns)
    return found


@functools.cache
def _top_cell_class(k: int, n: int) -> tuple[AnyGraph, ...]:
    return pg.top_cell_graph(k, n).square_move_class()


@st.composite
def rich_fillings(draw: st.DrawFn) -> tuple[list[list[int]], int]:
    """Le-diagrams of large cells: near-rectangular shapes, mostly ones."""
    n = draw(st.integers(4, 6))
    k = draw(st.integers(2, n - 2))
    trims = sorted(draw(st.lists(st.integers(0, 1), min_size=k, max_size=k)))
    bit = st.just(1) if draw(st.booleans()) else st.sampled_from([1, 1, 1, 0])
    raw = [
        draw(st.lists(bit, min_size=n - k - trim, max_size=n - k - trim))
        for trim in trims
    ]
    return _le_close(raw), n


@st.composite
def reduced_graphs(draw: st.DrawFn, max_n: int = 6) -> AnyGraph:
    """Reduced graphs from three sources, scrambled by moves and rotations.

    Le-graphs of arbitrary cells, arbitrary members of the square-move class
    of a large cell (these carry squares, so (M1) is exercised), and wiring
    graphs of reduced words.
    """
    source = draw(st.sampled_from(["le", "class", "word"]))
    if source == "le":
        filling, n = draw(le_fillings(max_n=max_n))
        graph = Graph.from_le_diagram(filling, n)
    elif source == "class":
        filling, n = draw(rich_fillings())
        members = Graph.from_le_diagram(filling, n).square_move_class()
        graph = draw(st.sampled_from(members))
    else:
        word, n = draw(words(max_n=max(2, max_n // 2)))
        assume(_is_reduced_word(word, n))
        graph = _any(Graph.from_wiring_diagram(word, n))
    return _scramble(draw, graph)


@st.composite
def digon_graphs(draw: st.DrawFn) -> AnyGraph:
    """Non-reduced by Definition 12.5: a reduced graph with a digon spliced in.

    The graph is then scrambled by moves, so it is only *move-equivalent* to
    a graph to which (R1) applies.
    """
    graph = draw(reduced_graphs(max_n=5))
    edges = _digon_edges(graph)
    assume(edges)
    spliced = _insert_digon(
        graph, draw(st.sampled_from(edges)), flip=draw(st.booleans())
    )
    return _scramble(draw, spliced)


@st.composite
def monogon_graphs(draw: st.DrawFn) -> AnyGraph:
    """Non-reduced by FWZ Definition 7.1.6: a hollow monogon spliced in."""
    graph = draw(reduced_graphs(max_n=5))
    edges = _digon_edges(graph)
    assume(edges)
    color = draw(st.sampled_from([BLACK, WHITE]))
    return _insert_monogon(graph, draw(st.sampled_from(edges)), color)


@st.composite
def triple_edge_graphs(draw: st.DrawFn) -> AnyGraph:
    """Non-reduced by Definition 12.5, with an essential self-intersection."""
    graph = draw(reduced_graphs(max_n=5))
    edges = _digon_edges(graph)
    assume(edges)
    return _insert_triple(graph, draw(st.sampled_from(edges)), flip=draw(st.booleans()))


@st.composite
def leafless_graphs(draw: st.DrawFn) -> AnyGraph:
    """Leafless graphs, reduced or not."""
    source = draw(st.sampled_from(["reduced", "digon", "monogon", "triple", "word"]))
    if source == "reduced":
        return draw(reduced_graphs())
    if source == "digon":
        return draw(digon_graphs())
    if source == "monogon":
        return draw(monogon_graphs())
    if source == "triple":
        return draw(triple_edge_graphs())
    word, n = draw(words())
    return _scramble(draw, _any(Graph.from_wiring_diagram(word, n)))


@st.composite
def squared_graphs(draw: st.DrawFn, *, reduced: bool = False) -> tuple[AnyGraph, int]:
    """A graph with a trivalent alternating square, and that square's face.

    Members of the square-move class of a large cell (any cell, not only the
    top one), optionally made non-reduced by a spliced digon.
    """
    filling, n = draw(rich_fillings())
    members = Graph.from_le_diagram(filling, n).square_move_class()
    graph = draw(st.sampled_from(members))
    graph = graph.cyclic_shift(draw(st.integers(0, n - 1)))
    assume(_mutable_squares(graph))
    if not reduced and draw(st.booleans()):
        graph = _insert_digon(graph, draw(st.sampled_from(_digon_edges(graph))))
    squares = _mutable_squares(graph)
    assume(squares)
    return _expose(graph, draw(st.sampled_from(squares)))


@functools.cache
def _dense_le_cases() -> tuple[tuple[tuple[tuple[int, ...], ...], int], ...]:
    """Return the small Le-diagrams with at least four ones (crossings appear)."""
    return tuple(
        (tuple(tuple(row) for row in filling), n)
        for k, n in [(2, 4), (2, 5), (3, 5)]
        for filling in _all_le_diagrams(k, n)
        if sum(map(sum, filling)) >= 4
    )


def _weighted_bipartite_le_graph(
    filling: Sequence[Sequence[int]],
    n: int,
    weights: Mapping[tuple[int, int], Fraction],
) -> tuple[
    bm.PlanarNetwork[int | tuple[int, int]],
    AnyGraph,
    dict[frozenset[Hashable], Fraction],
]:
    """Return the weighted Gamma-network, a bipartite Le-graph, and its weights.

    The Le-graph keeps the network's edges (same indices, tail first) and
    adds weight-1 connector and lollipop edges. An edge leaving a black
    vertex carries the inverse weight (Lam Proposition 5.3, as in
    ``PlanarNetwork.to_bipartite``; a boundary vertex counts as colored
    opposite to its neighbor), and every unicolored edge is subdivided by
    (M3), the new half getting weight 1.
    """
    oriented = bm.PlanarNetwork.from_le_diagram(filling, n, weights)
    graph = Graph.from_le_diagram(filling, n)
    plan = []
    for index, (tail, head) in enumerate(graph.edges):
        weight = (
            oriented.edges[index][2] if index < len(oriented.edges) else Fraction(1)
        )
        if tail in graph.boundary:
            black_tail = graph.color(head) == WHITE
        else:
            black_tail = graph.color(tail) == BLACK
        plan.append((tail, head, 1 / weight if black_tail else weight))
    weight_of: dict[frozenset[Hashable], Fraction] = {}
    for tail, head, weight in plan:
        internal = tail not in graph.boundary and head not in graph.boundary
        if internal and graph.color(tail) == graph.color(head):
            middle = ("mid", tail, head)
            graph = graph.insert_middle_vertex(
                graph.edges.index((tail, head)), -graph.color(tail), middle
            )
            weight_of[frozenset({tail, middle})] = weight
            weight_of[frozenset({middle, head})] = Fraction(1)
        else:
            weight_of[frozenset({tail, head})] = weight
    return oriented, graph, weight_of


def _postnikov_necklace(graph: AnyGraph) -> tuple[frozenset[Hashable], ...]:
    """Return the Postnikov necklace of ``pi_G^:`` on the boundary labels."""
    decorated = graph.to_decorated_permutation()
    return tuple(
        frozenset(graph.boundary[i - 1] for i in entry)
        for entry in decorated.grassmann_necklace(POSTNIKOV)
    )


def _positions(
    graph: AnyGraph, labels: Iterable[Iterable[Hashable]]
) -> list[frozenset[int]]:
    """Translate sets of boundary labels into sets of positions ``1..n``."""
    position = {b: i for i, b in enumerate(graph.boundary, start=1)}
    return [frozenset(position[b] for b in label) for label in labels]


def _leaf_fixture() -> AnyGraph:
    """A black trivalent vertex on two legs carrying a white leaf ``u``."""
    return _graph(
        (1, 2),
        [(1, "v"), ("v", 2), ("v", "u")],
        {"v": [1, 0, 2], "u": [2]},
        {"v": BLACK, "u": WHITE},
    )


def _top_cell_permutation(k: int, n: int) -> dp.DecoratedPermutation:
    if k == 0:
        return dp.DecoratedPermutation(tuple(range(1, n + 1)))
    if k == n:
        return dp.DecoratedPermutation(
            tuple(range(1, n + 1)), frozenset(range(1, n + 1))
        )
    return dp.DecoratedPermutation(tuple((i + k - 1) % n + 1 for i in range(1, n + 1)))


# --------------------------------------------------------------------------- #
# Canonical examples — what each block of the page certifies
# --------------------------------------------------------------------------- #
class TestLollipopGraphs:
    COLORS = (BLACK, WHITE, WHITE, BLACK, WHITE)

    def test_lollipop_graph_is_reduced_with_identity_trip_permutation(self):
        graph = pg.lollipop_graph(self.COLORS)
        assert graph.is_reduced() is True
        assert graph.trip_permutation == (1, 2, 3, 4, 5)

    def test_fixed_points_are_decorated_by_the_lollipop_colors(self):
        decorated = pg.lollipop_graph(self.COLORS).to_decorated_permutation()
        assert decorated == dp.DecoratedPermutation.from_colors(
            (1, 2, 3, 4, 5), dict(enumerate(self.COLORS, start=1))
        )

    def test_single_face_matches_dimension_zero(self):
        graph = pg.lollipop_graph(self.COLORS)
        assert graph.face_count == 1
        assert bm.cell_dimension(graph.matroid()) == graph.face_count - 1 == 0

    def test_black_lollipops_are_loops_and_white_lollipops_coloops(self):
        positroid = pg.lollipop_graph(self.COLORS).matroid()
        assert positroid.loops == frozenset({1, 4})
        assert positroid.coloops == frozenset({2, 3, 5})

    def test_face_count_formula_floor(self):
        graph = pg.lollipop_graph(self.COLORS)
        decorated = graph.to_decorated_permutation()
        k, n = graph.graph_type
        assert graph.face_count == k * (n - k) - decorated.alignment_number + 1 == 1


class TestGr24SquareMovePair:
    FIRST = _subsets("12", "23", "34", "14", "13")
    SECOND = _subsets("12", "23", "34", "14", "24")

    def test_face_label_collections_are_the_pages(self):
        first, second = pg.gr24_square_pair()
        assert frozenset(first.face_labels()) == self.FIRST
        assert frozenset(second.face_labels()) == self.SECOND

    def test_both_graphs_are_reduced_graphs_of_the_top_cell(self):
        for graph in pg.gr24_square_pair():
            assert graph.is_reduced() is True
            assert graph.to_positroid().bases == {
                frozenset(pair) for pair in itertools.combinations((1, 2, 3, 4), 2)
            }

    def test_top_cell_has_exactly_these_two_graphs_up_to_strong_equivalence(self):
        first, _ = pg.gr24_square_pair()
        found = {frozenset(g.face_labels()) for g in first.square_move_class()}
        assert found == {self.FIRST, self.SECOND}

    def test_both_collections_are_maximal_weakly_separated_of_size_five(self):
        candidates = list(itertools.combinations(range(1, 5), 2))
        maximal = set(pg.maximal_weakly_separated_collections(4, candidates))
        assert maximal == {self.FIRST, self.SECOND}
        assert {len(c) for c in maximal} == {2 * 2 + 1}

    def test_collections_differ_by_the_single_mutation_13_24(self):
        assert _subsets("13") == self.FIRST - self.SECOND
        assert _subsets("24") == self.SECOND - self.FIRST

    def test_the_square_move_is_the_mutation(self):
        first, second = pg.gr24_square_pair()
        (square,) = [f for f in range(first.face_count) if len(first.faces()[f]) == 4]
        assert first.face_labels()[square] == frozenset({1, 3})
        assert first.square_move(square) == second

    def test_square_move_realizes_the_three_term_plucker_relation(self):
        """The exchanged labels are one side of the relation, the kept ones the rest."""
        network = bm.square_network(2, 3, 5, 7)
        first, second = pg.gr24_square_pair()
        assert Graph.from_planar_bipartite_network(network) == first
        before, after = frozenset(first.face_labels()), frozenset(second.face_labels())
        (old,) = before - after
        (new,) = after - before
        kept = before & after
        pairs = {
            frozenset({a, b}) for a, b in itertools.combinations(kept, 2) if not a & b
        }
        assert (old, new) == (frozenset({1, 3}), frozenset({2, 4}))
        assert pairs == {_subsets("12", "34"), _subsets("14", "23")}
        assert network.plucker(old) * network.plucker(new) == sum(
            math.prod((network.plucker(label) for label in pair), start=Fraction(1))
            for pair in pairs
        )


class TestGr2nTriangulations:
    @pytest.mark.parametrize(("n", "catalan"), [(4, 2), (5, 5), (6, 14)])
    def test_reduced_graphs_of_the_top_cell_are_counted_by_catalan(self, n, catalan):
        assert catalan == math.comb(2 * (n - 2), n - 2) // (n - 1)
        classes = pg.top_cell_graph(2, n).square_move_class()
        assert len(classes) == catalan

    @pytest.mark.parametrize(("n", "catalan"), [(4, 2), (5, 5), (6, 14)])
    def test_maximal_weakly_separated_collections_are_counted_by_catalan(
        self, n, catalan
    ):
        candidates = list(itertools.combinations(range(1, n + 1), 2))
        maximal = pg.maximal_weakly_separated_collections(n, candidates)
        assert len(maximal) == catalan

    @pytest.mark.parametrize("n", [4, 5, 6])
    def test_face_labels_are_the_sides_and_diagonals_of_a_triangulation(self, n):
        sides = {frozenset({i, i % n + 1}) for i in range(1, n + 1)}
        for graph in pg.top_cell_graph(2, n).square_move_class():
            labels = _positions(graph, graph.face_labels())
            diagonals = [sorted(label) for label in set(labels) - sides]
            assert sides <= set(labels)
            assert len(diagonals) == n - 3
            for (a, c), (b, d) in itertools.combinations(diagonals, 2):
                assert not a < b < c < d
                assert not b < a < d < c


class TestCountsForKThree:
    @pytest.mark.parametrize(("n", "count"), [(6, 34), (7, 259)])
    def test_counts_from_clique_enumeration(self, n, count):
        candidates = list(itertools.combinations(range(1, n + 1), 3))
        maximal = pg.maximal_weakly_separated_collections(n, candidates)
        assert len(maximal) == count
        assert {len(c) for c in maximal} == {3 * (n - 3) + 1}

    @pytest.mark.parametrize(("n", "count"), [(6, 34), (7, 259)])
    def test_counts_from_actual_plabic_graphs(self, n, count):
        classes = pg.top_cell_graph(3, n).square_move_class()
        assert len(classes) == count
        assert {graph.face_count for graph in classes} == {3 * (n - 3) + 1}


class TestNonReducedWiringDiagram:
    def test_the_graph_satisfies_every_axiom(self):
        graph = pg.nonreduced_wiring_example()
        assert Graph.from_dataframe(graph.to_dataframe()) == graph
        assert graph.size == 4

    def test_the_graph_is_not_reduced(self):
        assert _is_reduced_word([1, 1], 2) is False
        assert pg.nonreduced_wiring_example().is_reduced() is False

    def test_the_two_lines_meeting_twice_are_a_bad_double_crossing(self):
        graph = pg.nonreduced_wiring_example()
        assert len(graph.bad_double_crossings()) == 1
        assert graph.roundtrips == ()
        assert graph.essential_self_intersections() == ()

    @pytest.mark.parametrize(
        ("word", "n"), [([], 2), ([1], 2), ([1, 2, 1], 3), ([2, 1, 3, 2], 4)]
    )
    def test_graphs_of_reduced_words_are_reduced(self, word, n):
        assert _is_reduced_word(word, n) is True
        assert Graph.from_wiring_diagram(word, n).is_reduced() is True

    @given(words())
    @settings(max_examples=100, deadline=None)
    def test_word_is_reduced_iff_its_graph_is_reduced(self, word_n):
        """FWZ Exercise 7.3.9, Remark 7.3.10."""
        word, n = word_n
        graph = Graph.from_wiring_diagram(word, n)
        assert graph.is_reduced() == _is_reduced_word(word, n)

    @given(words())
    @settings(max_examples=100, deadline=None)
    def test_bad_features_mean_two_lines_meet_more_than_once(self, word_n):
        """FWZ Remark 7.8.6, on wiring diagrams."""
        word, n = word_n
        graph = Graph.from_wiring_diagram(word, n)
        bad = bool(
            graph.roundtrips
            or graph.essential_self_intersections()
            or graph.bad_double_crossings()
        )
        wires = list(range(n))
        meetings: dict[frozenset[int], int] = {}
        for letter in word:
            pair = frozenset({wires[letter - 1], wires[letter]})
            meetings[pair] = meetings.get(pair, 0) + 1
            wires[letter - 1], wires[letter] = wires[letter], wires[letter - 1]
        assert bad == any(count > 1 for count in meetings.values())

    @given(words())
    @settings(max_examples=50, deadline=None)
    def test_trips_from_the_left_follow_the_wires(self, word_n):
        word, n = word_n
        graph = Graph.from_wiring_diagram(word, n)
        wires = _permutation_of_word(word, n)
        for height, wire in enumerate(wires, start=1):
            assert graph.trip_permutation[wire] == 2 * n + 1 - height


class TestHollowDigon:
    def test_two_trivalent_vertices_of_opposite_colors_joined_twice(self):
        graph = pg.hollow_digon()
        assert graph.edges.count(("u", "w")) == 2
        assert graph.degree("u") == graph.degree("w") == 3
        assert graph.color("u") == -graph.color("w")

    def test_it_is_a_legitimate_plabic_graph_that_is_not_reduced(self):
        graph = pg.hollow_digon()
        assert Graph.from_dataframe(graph.to_dataframe()) == graph
        assert graph.is_reduced() is False

    def test_parallel_edge_reduction_applies_as_is(self):
        graph = pg.hollow_digon()
        (digon,) = [f for f in range(graph.face_count) if len(graph.faces()[f]) == 2]
        reduced = graph.parallel_edge_reduction(digon)
        assert reduced.edges == ((1, 2),)
        assert reduced.internal_vertices == ()


# --------------------------------------------------------------------------- #
# Structural theorems
# --------------------------------------------------------------------------- #
class TestMoveInvariance:
    """Postnikov Lemma 13.1."""

    @given(squared_graphs())
    @settings(max_examples=100, deadline=None)
    def test_square_move_preserves_the_trip_permutation(self, graph_face):
        graph, face = graph_face
        assert graph.square_move(face).trip_permutation == graph.trip_permutation

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_contraction_preserves_the_trip_permutation(self, graph, data):
        unicolored = [
            index
            for index, (u, w) in enumerate(graph.edges)
            if u != w
            and u in graph.internal_vertices
            and w in graph.internal_vertices
            and graph.color(u) == graph.color(w)
        ]
        assume(unicolored)
        moved = graph.contract_edge(data.draw(st.sampled_from(unicolored)))
        assert moved.trip_permutation == graph.trip_permutation

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_uncontraction_preserves_the_trip_permutation(self, graph, data):
        big = [v for v in graph.internal_vertices if graph.degree(v) >= 3]
        assume(big)
        moved = _split_somewhere(data.draw, graph, data.draw(st.sampled_from(big)))
        assert moved.trip_permutation == graph.trip_permutation

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_middle_vertex_insertion_preserves_the_trip_permutation(self, graph, data):
        assume(graph.edges)
        edge = data.draw(st.integers(0, len(graph.edges) - 1))
        color = data.draw(st.sampled_from([BLACK, WHITE]))
        moved = graph.insert_middle_vertex(edge, color, "middle")
        assert moved.trip_permutation == graph.trip_permutation

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_middle_vertex_removal_preserves_the_trip_permutation(self, graph, data):
        assume(graph.edges)
        edge = data.draw(st.integers(0, len(graph.edges) - 1))
        graph = graph.insert_middle_vertex(edge, WHITE, "middle")
        middles = [
            v
            for v in graph.internal_vertices
            if graph.degree(v) == 2
            and len({d >> 1 for d in dict(graph.rotations)[v]}) == 2
        ]
        moved = graph.remove_middle_vertex(data.draw(st.sampled_from(middles)))
        assert moved.trip_permutation == graph.trip_permutation

    def test_parallel_edge_reduction_changes_the_trip_permutation(self):
        graph = pg.hollow_digon()
        assert graph.trip_permutation == (1, 2)
        assert graph.parallel_edge_reduction(2).trip_permutation == (2, 1)

    @given(reduced_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_parallel_edge_reduction_changes_the_trip_permutation_in_general(
        self, graph, data
    ):
        edges = _digon_edges(graph)
        assume(edges)
        spliced = _insert_digon(
            graph, data.draw(st.sampled_from(edges)), flip=data.draw(st.booleans())
        )
        reduced = spliced.parallel_edge_reduction(_face_on_digon(spliced))
        assert reduced.trip_permutation == graph.trip_permutation
        assert reduced.trip_permutation != spliced.trip_permutation

    def test_parallel_edge_reduction_cannot_change_a_permutation_of_one_letter(self):
        """The page's "reductions change it" is generic, not universal.

        A digon spliced into the edge of the only lollipop: (R1) applies,
        and the trip permutation is the identity of ``S_1`` before and
        after. What changes is that the fixed point gains a decoration.
        """
        spliced = _insert_digon(_any(pg.lollipop_graph([WHITE])), 0)
        reduced = spliced.parallel_edge_reduction(_face_on_digon(spliced))
        assert spliced.trip_permutation == reduced.trip_permutation == (1,)
        assert spliced.lollipops == {}
        assert reduced.lollipops == {1: WHITE}

    def test_parallel_edge_reduction_keeps_the_permutation_when_one_trip_uses_the_digon(
        self,
    ):
        """The same exception on a leafless graph with two boundary vertices.

        (R1) changes the trip permutation exactly when the two trips through
        the digon are distinct; here a single trip runs through it both ways.
        """
        exposed = _triple_edge_fixture().uncontract_vertex("a", [2, 3], "a2")
        (digon,) = [
            f
            for f, darts in enumerate(exposed.faces())
            if {exposed.edges[d >> 1][d & 1] for d in darts} == {"a2", "b"}
        ]
        reduced = exposed.parallel_edge_reduction(digon)
        assert exposed.trip_permutation == reduced.trip_permutation == (2, 1)

    def test_leaf_reduction_changes_the_trip_permutation(self):
        graph = _leaf_fixture()
        assert graph.trip_permutation == (2, 1)
        assert graph.leaf_reduction("u", ["x", "y"]).trip_permutation == (1, 2)

    @given(reduced_graphs(max_n=5), st.data())
    @settings(max_examples=100, deadline=None)
    def test_leaf_reduction_changes_the_trip_permutation_in_general(self, graph, data):
        plain = [
            v
            for v in graph.internal_vertices
            if graph.degree(v) >= 2
            and len(set(_edge_rotations(graph)[v])) == graph.degree(v)
        ]
        assume(plain)
        v = data.draw(st.sampled_from(plain))
        leafy = _attach_leaf(graph, v, data.draw(st.integers(0, graph.degree(v) - 1)))
        labels = [("cut", k) for k in range(graph.degree(v))]
        reduced = leafy.leaf_reduction("leaf", labels)
        assert reduced.trip_permutation != leafy.trip_permutation


class TestReducednessCriterion:
    """Postnikov Theorem 13.2, against independent sources of (non-)reduced graphs."""

    @given(le_fillings())
    @settings(max_examples=100, deadline=None)
    def test_le_graphs_are_reduced(self, filling_n):
        """Postnikov section 20: ``G_D`` is a reduced plabic graph."""
        filling, n = filling_n
        assert Graph.from_le_diagram(filling, n).is_reduced() is True

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_moves_do_not_change_reducedness(self, graph, data):
        """Definition 12.5 is a property of the move-equivalence class."""
        assert _scramble(data.draw, graph).is_reduced() == graph.is_reduced()

    @given(digon_graphs())
    @settings(max_examples=100, deadline=None)
    def test_graph_move_equivalent_to_one_with_a_hollow_digon_is_not_reduced(
        self, graph
    ):
        """Definition 12.5: some move-equivalent graph admits (R1)."""
        assert graph.is_reduced() is False

    @given(monogon_graphs())
    @settings(max_examples=100, deadline=None)
    def test_graph_with_a_hollow_monogon_is_not_reduced(self, graph):
        """FWZ Definition 7.1.6."""
        assert graph.is_reduced() is False

    @given(reduced_graphs(max_n=5))
    @settings(max_examples=50, deadline=None)
    def test_no_graph_in_the_class_of_a_reduced_graph_has_a_monogon_or_digon(
        self, graph
    ):
        """FWZ Definition 7.1.6, searched over the whole square-move class."""
        assert graph.is_reduced() is True
        for other in graph.square_move_class():
            internal = [
                darts
                for darts in other.faces()
                if all(other.edges[d >> 1][d & 1] not in other.boundary for d in darts)
                and all(
                    other.edges[d >> 1][1 - (d & 1)] not in other.boundary
                    for d in darts
                )
            ]
            assert all(len(darts) > 2 for darts in internal)

    def test_condition_1_alone_fails_on_a_hollow_monogon(self):
        graph = _graph(
            (1, 2), [(1, "v"), ("v", "v"), ("v", 2)], {"v": [0, 1, 1, 2]}, {"v": BLACK}
        )
        assert len(graph.roundtrips) == 1
        assert graph.essential_self_intersections() == ()
        assert graph.bad_double_crossings() == ()
        assert graph.trip_permutation == (2, 1)
        assert graph.is_reduced() is False
        assert graph.has_resonance_property() is False

    def test_condition_2_alone_fails_on_a_triple_edge(self):
        graph = _triple_edge_fixture()
        assert graph.roundtrips == ()
        assert graph.essential_self_intersections() == (2, 3, 4)
        assert graph.bad_double_crossings() == ()
        assert graph.trip_permutation == (2, 1)
        assert graph.is_reduced() is False
        assert graph.has_resonance_property() is False

    def test_the_triple_edge_is_move_equivalent_to_a_graph_with_a_hollow_digon(self):
        """Definition 12.5 sees the triple edge without the trip criterion."""
        exposed = _triple_edge_fixture().uncontract_vertex("a", [2, 3], "a2")
        (digon,) = [
            f
            for f, darts in enumerate(exposed.faces())
            if {exposed.edges[d >> 1][d & 1] for d in darts} == {"a2", "b"}
        ]
        assert exposed.parallel_edge_reduction(digon).edges != ()

    @given(triple_edge_graphs())
    @settings(max_examples=100, deadline=None)
    def test_graph_with_a_spliced_triple_edge_is_not_reduced(self, graph):
        """Definition 12.5, on graphs whose only tell is condition (2) or worse."""
        assert graph.essential_self_intersections() != ()
        assert graph.is_reduced() is False

    def test_condition_3_alone_fails_on_the_non_reduced_wiring_diagram(self):
        graph = pg.nonreduced_wiring_example()
        assert graph.roundtrips == ()
        assert graph.essential_self_intersections() == ()
        assert len(graph.bad_double_crossings()) == 1
        assert all(graph.trip_permutation[i] != i + 1 for i in range(graph.size))
        assert graph.is_reduced() is False
        assert graph.has_resonance_property() is False

    def test_conditions_3_and_4_fail_on_the_hollow_digon(self):
        graph = pg.hollow_digon()
        assert graph.roundtrips == ()
        assert graph.essential_self_intersections() == ()
        assert len(graph.bad_double_crossings()) == 1
        assert graph.trip_permutation == (1, 2)
        assert graph.lollipops == {}
        assert graph.has_resonance_property() is False

    def test_conditions_1_and_2_fail_on_a_graph_with_a_loop_behind_a_digon(self):
        graph = _graph(
            (1,),
            [(1, "a"), ("a", "b"), ("a", "b"), ("b", "c"), ("c", "c")],
            {"a": [0, 1, 2], "b": [1, 3, 2], "c": [3, 4, 4]},
            {"a": BLACK, "b": WHITE, "c": BLACK},
        )
        assert graph.roundtrips != ()
        assert graph.essential_self_intersections() != ()
        assert graph.is_reduced() is False
        assert graph.has_resonance_property() is False

    def test_criterion_is_not_stated_for_graphs_with_leaves(self):
        with pytest.raises(ValueError, match=r"Theorem 13\.2.*leafless"):
            _leaf_fixture().is_reduced()


class TestFundamentalTheorem:
    """Postnikov Theorem 13.4 and its corollaries."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_same_decorated_trip_permutation_iff_move_equivalent(self, n):
        pool = [g.cyclic_shift(s) for g in _le_graphs(n) for s in range(min(n, 2))]
        classes = [
            {
                frozenset(_positions(h, h.face_labels()))
                for h in graph.square_move_class()
            }
            for graph in pool
        ]
        for (first, reachable), second in itertools.product(
            zip(pool, classes, strict=True), pool
        ):
            same = first.to_decorated_permutation() == second.to_decorated_permutation()
            labels = frozenset(_positions(second, second.face_labels()))
            assert (labels in reachable) == same

    @given(reduced_graphs(max_n=5))
    @settings(max_examples=50, deadline=None)
    def test_move_equivalent_graphs_share_the_decorated_trip_permutation(self, graph):
        decorated = graph.to_decorated_permutation()
        for other in graph.square_move_class():
            assert other.is_reduced() is True
            assert other.to_decorated_permutation() == decorated

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
    def test_every_decorated_permutation_comes_from_a_reduced_graph(self, n):
        """Postnikov Corollary 14.7; FWZ Corollary 7.10.4."""
        graphs = _le_graphs(n)
        assert all(graph.is_reduced() for graph in graphs)
        realized = {graph.to_decorated_permutation() for graph in graphs}
        assert realized == set(dp.enumerate_decorated_permutations(n))

    @pytest.mark.parametrize(("n", "count"), [(1, 2), (2, 5), (3, 16), (4, 65)])
    def test_move_classes_are_counted_by_decorated_permutations(self, n, count):
        """FWZ Exercise 7.1.19: ``n! sum 1/k!`` (OEIS A000522)."""
        assert count == sum(
            math.factorial(n) // math.factorial(k) for k in range(n + 1)
        )
        realized = {graph.to_decorated_permutation() for graph in _le_graphs(n)}
        assert len(realized) == count

    def test_reductions_reach_a_reduced_graph_on_the_fixtures(self):
        """Postnikov Lemma 13.6, certified on the two reduction fixtures."""
        assert pg.hollow_digon().parallel_edge_reduction(2).is_reduced() is True
        assert _leaf_fixture().leaf_reduction("u", ["x", "y"]).is_reduced() is True


class TestTripPermutationReadsOffTheCell:
    """Postnikov Proposition 16.4, Theorem 17.1."""

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_necklace_of_the_matroid_is_the_necklace_of_the_permutation(self, graph):
        assert graph.matroid().grassmann_necklace == _postnikov_necklace(graph)

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_type_is_the_number_of_anti_exceedances(self, graph):
        decorated = graph.to_decorated_permutation()
        k, n = graph.graph_type
        assert (decorated.anti_exceedance_count, decorated.size) == (k, n)
        assert graph.matroid().rank() == k

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_conversions_agree_with_the_matroid_of_orientations(self, graph):
        assert graph.to_positroid() == graph.matroid()
        assert graph.to_grassmann_necklace().entries == _postnikov_necklace(graph)

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_cells_biject_with_decorated_permutations_of_each_type(self, n):
        for k in range(n + 1):
            graphs = [Graph.from_le_diagram(d, n) for d in _all_le_diagrams(k, n)]
            cells = {graph.matroid() for graph in graphs}
            permutations = {graph.to_decorated_permutation() for graph in graphs}
            assert len(cells) == len(permutations) == len(graphs)
            assert {p.anti_exceedance_count for p in permutations} <= {k}


class TestParameterizationOfCells:
    """Postnikov Theorem 12.7, the combinatorial part."""

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_reduced_graphs_are_perfectly_orientable(self, graph):
        assert graph.is_perfectly_orientable is True

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_cell_dimension_is_face_count_minus_one(self, graph):
        assert bm.cell_dimension(graph.matroid()) == graph.face_count - 1

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_every_cell_arises_from_a_reduced_graph(self, n):
        cells = {graph.matroid() for graph in _le_graphs(n)}
        assert cells == set(ps.enumerate_positroids(n))


class TestMatroidViaOrientationsAndMatchings:
    """Postnikov Proposition 11.7, Lemma 11.10; Williams ICM Theorem 2.17.

    The last is the quantitative form (Postnikov, refined by Talaska 2008
    and Speyer 2016): Plucker coordinates are weighted matching sums.
    """

    @given(leafless_graphs())
    @settings(max_examples=100, deadline=None)
    def test_matroid_of_a_perfectly_orientable_graph_is_a_positroid(self, graph):
        assume(graph.is_perfectly_orientable)
        matroid = graph.matroid()
        assert isinstance(matroid, ps.Positroid)
        assert ps.is_positroid(matroid.to_matroid()) is True
        assert matroid.bases == {
            graph.source_set(o) for o in graph.perfect_orientations()
        }

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_all_positroids_arise_from_plabic_graphs(self, n):
        assert {g.matroid() for g in _le_graphs(n)} == set(ps.enumerate_positroids(n))

    @given(leafless_graphs())
    @settings(max_examples=100, deadline=None)
    def test_bases_of_a_bipartite_graph_are_boundaries_of_matchings(self, graph):
        bipartite = graph.normalized()
        assume(bipartite.is_bipartite and bipartite.is_perfectly_orientable)
        assume(
            not any(set(edge) <= set(bipartite.boundary) for edge in bipartite.edges)
        )
        boundaries = {
            bipartite.matching_boundary(m) for m in bipartite.almost_perfect_matchings()
        }
        assert boundaries == bipartite.matroid().bases

    @given(*[st.fractions(min_value=Fraction(1, 9), max_value=9)] * 4)
    @settings(max_examples=50, deadline=None)
    def test_plucker_coordinates_are_weighted_matching_sums_on_the_square(
        self, a, b, c, d
    ):
        """Postnikov/Talaska/Speyer via Williams ICM Theorem 2.17, on Lam's square.

        The left side is the boundary measurement point itself — the minors
        of a perfect orientation's measurement matrix — so the comparison
        is projective.
        """
        network = bm.square_network(a, b, c, d)
        graph = Graph.from_planar_bipartite_network(network)
        sums: dict[frozenset[int | str], Fraction] = {}
        for matching in graph.almost_perfect_matchings():
            weight = math.prod(
                (network.edges[index][2] for index in matching), start=Fraction(1)
            )
            key = graph.matching_boundary(matching)
            sums[key] = sums.get(key, Fraction(0)) + weight
        oriented = network.to_perfect_orientation()
        scale = sums[frozenset(oriented.source_set)]
        for pair in itertools.combinations((1, 2, 3, 4), 2):
            assert oriented.plucker(pair) * scale == sums[frozenset(pair)]

    @given(st.data())
    @settings(max_examples=60, deadline=None)
    def test_plucker_coordinates_are_weighted_matching_sums_on_le_graphs(self, data):
        """Postnikov/Talaska/Speyer via Williams ICM Theorem 2.17, on Le-graphs.

        ``p_I(Meas(N))`` is the sum of ``w(M)`` over matchings with boundary
        ``I``; on these graphs the sums have several terms.
        """
        filling, n = data.draw(st.sampled_from(_dense_le_cases()))
        weight = st.fractions(min_value=Fraction(1, 9), max_value=9)
        weights = {
            (i, j): data.draw(weight)
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
            if value
        }
        oriented, graph, weight_of = _weighted_bipartite_le_graph(filling, n, weights)
        sums: dict[frozenset[Hashable], Fraction] = {}
        for matching in graph.almost_perfect_matchings():
            term = math.prod(
                (weight_of[frozenset(graph.edges[index])] for index in matching),
                start=Fraction(1),
            )
            key = graph.matching_boundary(matching)
            sums[key] = sums.get(key, Fraction(0)) + term
        scale = sums[frozenset(oriented.source_set)]
        for subset in itertools.combinations(range(1, n + 1), len(oriented.source_set)):
            expected = sums.get(frozenset(subset), Fraction(0))
            assert oriented.plucker(subset) * scale == expected

    def test_matching_sums_on_the_full_rectangle_have_several_terms(self):
        """Guard against the additive half of the formula going untested."""
        _, graph, _ = _weighted_bipartite_le_graph([[1, 1], [1, 1]], 4, {})
        boundaries = [
            graph.matching_boundary(m) for m in graph.almost_perfect_matchings()
        ]
        assert len(boundaries) > len(set(boundaries))


class TestFaceCount:
    """Oh-Postnikov-Speyer Theorem 6.8; FWZ Corollaries 7.10.5, 7.10.7."""

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_face_count_is_length_plus_one(self, graph):
        decorated = graph.to_decorated_permutation()
        k, n = graph.graph_type
        assert graph.face_count == k * (n - k) - decorated.alignment_number + 1

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_at_most_k_n_minus_k_plus_one_faces_with_equality_iff_top_cell(self, graph):
        k, n = graph.graph_type
        assert graph.face_count <= k * (n - k) + 1
        is_top = graph.to_decorated_permutation() == _top_cell_permutation(k, n)
        assert (graph.face_count == k * (n - k) + 1) == is_top

    @pytest.mark.parametrize(
        ("k", "n"), [(0, 3), (1, 3), (2, 4), (2, 5), (3, 6), (4, 4)]
    )
    def test_top_cell_graph_attains_the_bound(self, k, n):
        graph = pg.top_cell_graph(k, n)
        assert graph.to_decorated_permutation() == _top_cell_permutation(k, n)
        assert graph.face_count == k * (n - k) + 1


class TestFaceLabels:
    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_every_face_label_is_a_k_subset(self, graph):
        """FWZ Theorem 7.12.4."""
        k, _ = graph.graph_type
        assert {len(label) for label in graph.face_labels("target")} == {k}
        assert {len(label) for label in graph.face_labels("source")} == {k}

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_boundary_faces_read_off_the_grassmann_necklace_in_order(self, graph):
        """Oh-Postnikov-Speyer section 6."""
        assert graph.boundary_face_labels() == graph.matroid().grassmann_necklace

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_face_labels_are_weakly_separated(self, graph):
        """Danilov-Karzanov-Koshevoy; Oh-Postnikov-Speyer; FWZ Theorem 7.13.2."""
        for labeling in ("target", "source"):
            labels = _positions(graph, graph.face_labels(labeling))
            assert pg.is_weakly_separated_collection(labels, graph.size) is True

    @given(reduced_graphs())
    @settings(max_examples=100, deadline=None)
    def test_faces_have_distinct_labels(self, graph):
        assert len(set(graph.face_labels())) == graph.face_count

    def test_face_labels_need_a_reduced_graph(self):
        with pytest.raises(ValueError, match="reduced"):
            pg.hollow_digon().face_labels()


class TestPurityAndMaximalCollections:
    """Oh-Postnikov-Speyer Theorems 1.3, 1.4, 1.5, 6.6; FWZ Theorem 7.13.16."""

    @given(st.one_of(le_fillings(max_n=6), rich_fillings()))
    @settings(max_examples=60, deadline=None)
    def test_maximal_collections_in_a_positroid_all_have_length_plus_one(
        self, filling_n
    ):
        filling, n = filling_n
        graph = Graph.from_le_diagram(filling, n)
        positroid = graph.matroid()
        k, _ = graph.graph_type
        length = k * (n - k) - graph.to_decorated_permutation().alignment_number
        maximal = pg.maximal_weakly_separated_collections(
            n,
            _positions(graph, positroid.bases),
            containing=_positions(graph, positroid.grassmann_necklace),
        )
        assert {len(collection) for collection in maximal} == {length + 1}

    @given(st.one_of(le_fillings(max_n=6), rich_fillings()))
    @settings(max_examples=60, deadline=None)
    def test_reduced_graphs_biject_with_maximal_collections(self, filling_n):
        filling, n = filling_n
        graph = Graph.from_le_diagram(filling, n)
        positroid = graph.matroid()
        from_graphs = [
            frozenset(_positions(g, g.face_labels())) for g in graph.square_move_class()
        ]
        maximal = pg.maximal_weakly_separated_collections(
            n,
            _positions(graph, positroid.bases),
            containing=_positions(graph, positroid.grassmann_necklace),
        )
        assert len(set(from_graphs)) == len(from_graphs)
        assert set(from_graphs) == set(maximal)

    @given(reduced_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_moves_m2_and_m3_fix_the_face_labels(self, graph, data):
        labels = frozenset(graph.face_labels())
        assert frozenset(graph.normalized().face_labels()) == labels
        # Subdividing a lollipop's edge makes a leaf, where labels are not computed.
        edges = [
            index
            for index, edge in enumerate(graph.edges)
            if all(graph.degree(v) > 1 or v in graph.boundary for v in edge)
        ]
        assume(edges)
        inserted = graph.insert_middle_vertex(
            data.draw(st.sampled_from(edges)), WHITE, "middle"
        )
        assert frozenset(inserted.face_labels()) == labels

    @given(squared_graphs(reduced=True))
    @settings(max_examples=100, deadline=None)
    def test_square_move_exchanges_sac_for_sbd(self, graph_face):
        graph, face = graph_face
        before = _positions(graph, graph.face_labels())
        moved = graph.square_move(face)
        after = _positions(moved, moved.face_labels())
        (old,) = set(before) - set(after)
        (new,) = set(after) - set(before)
        assert old == before[face]
        assert new == after[face]
        (a, c), (b, d) = sorted(old - new), sorted(new - old)
        assert old & new == old - {a, c} == new - {b, d}
        assert a < b < c < d or b < a < d < c

    @pytest.mark.parametrize(("k", "n"), [(2, 5), (3, 6)])
    def test_boundary_faces_of_the_top_cell_are_the_cyclic_intervals(self, k, n):
        """The frozen variables of the cluster structure (Scott; OPS Theorem 7.1)."""
        intervals = tuple(
            frozenset((i + offset - 1) % n + 1 for offset in range(k))
            for i in range(1, n + 1)
        )
        for graph in pg.top_cell_graph(k, n).square_move_class():
            assert graph.boundary_face_labels() == intervals


class TestResonanceCriterion:
    """Kodama-Williams, via FWZ Theorem 7.11.5."""

    @given(leafless_graphs())
    @settings(max_examples=150, deadline=None)
    def test_leafless_graph_is_reduced_iff_it_has_the_resonance_property(self, graph):
        assert graph.has_resonance_property() == graph.is_reduced()

    def test_resonance_is_not_stated_for_graphs_with_leaves(self):
        with pytest.raises(ValueError, match=r"7\.11\.5.*leafless"):
            _leaf_fixture().has_resonance_property()


class TestLeGraph:
    """Postnikov section 20, Corollary 20.1."""

    @given(le_fillings())
    @settings(max_examples=100, deadline=None)
    def test_anti_exceedance_set_is_the_source_set_of_the_shape(self, filling_n):
        filling, n = filling_n
        network = bm.PlanarNetwork.from_le_diagram(filling, n)
        decorated = Graph.from_le_diagram(filling, n).to_decorated_permutation()
        assert decorated.anti_exceedances == network.source_set

    @given(le_fillings())
    @settings(max_examples=100, deadline=None)
    def test_le_graph_has_the_cell_of_its_gamma_network(self, filling_n):
        filling, n = filling_n
        network = bm.PlanarNetwork.from_le_diagram(filling, n)
        graph = Graph.from_le_diagram(filling, n)
        assert graph.to_positroid() == network.to_positroid()
        # Not this page's claim: dim = |D| is Postnikov Theorem 6.5 (Le-Diagram
        # page), combined with dim = |F| - 1 (Theorem 12.7).
        assert graph.face_count == sum(map(sum, filling)) + 1

    @given(le_fillings(max_n=5), st.data())
    @settings(max_examples=60, deadline=None)
    def test_positively_weighted_le_network_lands_in_the_cell_of_the_graph(
        self, filling_n, data
    ):
        """Postnikov Theorem 12.7, image half: Meas of a reduced graph hits its cell.

        Corollary 16.5 (the same for non-reduced graphs) stays untested:
        the weighted networks available here are all Le-networks.
        """
        filling, n = filling_n
        weight = st.fractions(min_value=Fraction(1, 9), max_value=9)
        weights = {
            (i, j): data.draw(weight)
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
            if value
        }
        network = bm.PlanarNetwork.from_le_diagram(filling, n, weights)
        assert network.to_positroid() == Graph.from_le_diagram(filling, n).matroid()

    @pytest.mark.parametrize(("k", "n"), [(1, 3), (2, 4), (2, 5)])
    def test_diagrams_of_a_type_biject_with_permutations_of_that_type(self, k, n):
        diagrams = _all_le_diagrams(k, n)
        realized = {
            Graph.from_le_diagram(d, n).to_decorated_permutation() for d in diagrams
        }
        expected = {
            p
            for p in dp.enumerate_decorated_permutations(n)
            if p.anti_exceedance_count == k
        }
        assert len(diagrams) == len(realized)
        assert realized == expected


# --------------------------------------------------------------------------- #
# Derived vocabulary
# --------------------------------------------------------------------------- #
class TestDerivedVocabulary:
    def test_lollipops_and_leaves_are_told_apart(self):
        graph = _leaf_fixture()
        assert graph.leaves == ("u",)
        assert graph.is_leafless is False
        assert graph.lollipops == {}
        reduced = graph.leaf_reduction("u", ["x", "y"])
        assert reduced.lollipops == {1: WHITE, 2: WHITE}
        assert reduced.is_leafless is True

    def test_subdividing_a_lollipop_edge_turns_the_lollipop_into_a_leaf(self):
        graph = _any(pg.lollipop_graph([WHITE])).insert_middle_vertex(
            0, BLACK, "middle"
        )
        assert graph.lollipops == {}
        assert graph.leaves == ((1, "lollipop"),)
        assert graph.trip_permutation == (1,)

    def test_trips_use_every_edge_once_in_each_direction(self):
        graph = pg.top_cell_graph(2, 5)
        darts = sorted(d for trip in graph.trips() for d in trip)
        assert darts == list(range(2 * len(graph.edges)))

    def test_essential_self_intersection_needs_different_internal_colors(self):
        graph = pg.lollipop_graph([BLACK])
        assert graph.trips() == ((0, 1),)
        assert graph.essential_self_intersections() == ()

    def test_type_follows_postnikov_lemma_9_4(self):
        first, _ = pg.gr24_square_pair()
        assert sum(c * (first.degree(v) - 2) for v, c in first.colors) == 0
        assert first.graph_type == (2, 4)
        assert {len(first.source_set(o)) for o in first.perfect_orientations()} == {2}

    def test_perfect_orientation_has_one_out_at_black_and_one_in_at_white(self):
        graph = pg.top_cell_graph(2, 4)
        for orientation in graph.perfect_orientations():
            tails = set(orientation)
            for v, darts in graph.rotations:
                if v in graph.boundary:
                    continue
                outgoing = sum(1 for d in darts if d in tails)
                wanted = outgoing if graph.color(v) == BLACK else len(darts) - outgoing
                assert wanted == 1

    def test_graph_with_white_leaves_on_a_black_vertex_is_not_perfectly_orientable(
        self,
    ):
        graph = _graph(
            (1,),
            [(1, "v"), ("v", "x"), ("v", "y")],
            {"v": [0, 1, 2], "x": [1], "y": [2]},
            {"v": BLACK, "x": WHITE, "y": WHITE},
        )
        assert graph.is_perfectly_orientable is False
        with pytest.raises(ValueError, match="perfectly orientable"):
            _ = graph.graph_type
        with pytest.raises(ValueError, match="not perfectly orientable"):
            graph.matroid()

    def test_fixed_point_without_lollipop_has_no_decoration(self):
        with pytest.raises(ValueError, match="lollipop"):
            pg.hollow_digon().to_decorated_permutation()

    def test_color_of_a_boundary_vertex_is_undefined(self):
        with pytest.raises(ValueError, match="carries no color"):
            pg.hollow_digon().color(1)


class TestWeakSeparation:
    """Oh-Postnikov-Speyer Definition 3.1."""

    @pytest.mark.parametrize(
        ("first", "second", "expected"),
        [
            ({1, 3}, {2, 4}, False),
            ({1, 2}, {3, 4}, True),
            ({1, 3}, {1, 2}, True),
            ({1, 3, 5}, {2, 4, 6}, False),
            ({1, 2, 5}, {3, 4, 5}, True),
            ({1, 4}, {1, 4}, True),
        ],
    )
    def test_chord_separates_the_differences(self, first, second, expected):
        assert pg.is_weakly_separated(first, second, 6) is expected

    @given(st.data())
    def test_definition_by_forbidden_cyclic_pattern(self, data):
        n = data.draw(st.integers(2, 7))
        k = data.draw(st.integers(1, n))
        subsets = st.sets(st.integers(1, n), min_size=k, max_size=k)
        first, second = data.draw(subsets), data.draw(subsets)
        only_first, only_second = first - second, second - first
        pattern = any(
            (a < b < a2 < b2)
            or (b < a2 < b2 < a)
            or (a2 < b2 < a < b)
            or (b2 < a < b < a2)
            for a, a2 in itertools.permutations(only_first, 2)
            for b, b2 in itertools.permutations(only_second, 2)
        )
        assert pg.is_weakly_separated(first, second, n) == (not pattern)

    def test_unequal_sizes_are_a_different_definition(self):
        with pytest.raises(ValueError, match="equal-size"):
            pg.is_weakly_separated({1}, {1, 2}, 3)

    def test_required_subsets_must_be_weakly_separated(self):
        with pytest.raises(ValueError, match="not pairwise weakly separated"):
            pg.maximal_weakly_separated_collections(4, [], containing=[{1, 3}, {2, 4}])

    def test_clique_enumeration_is_guarded(self):
        candidates = list(itertools.combinations(range(1, 9), 4))
        with pytest.raises(ValueError, match="enumeration guard"):
            pg.maximal_weakly_separated_collections(8, candidates)


# --------------------------------------------------------------------------- #
# Round-trip laws
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(leafless_graphs())
    @settings(max_examples=100, deadline=None)
    def test_dataframe_round_trip(self, graph):
        frame = graph.to_dataframe()
        assert list(frame.columns) == ["vertex", "kind", "slot", "edge", "end"]
        assert len(frame) == 2 * len(graph.edges)
        assert Graph.from_dataframe(frame) == graph

    @pytest.mark.parametrize(
        "example",
        [
            pg.hollow_digon(),
            pg.gr24_square_pair()[0],
            pg.lollipop_graph([BLACK, WHITE]),
            pg.nonreduced_wiring_example(),
            _graph((), [], {}, {}),
        ],
    )
    def test_experiment_io_round_trip(self, example, tmp_path):
        path = io.write_result(example.to_dataframe(), tmp_path / "plabic.json")
        decoded = Graph.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == example

    def test_constructors_agree_on_the_square(self):
        network = bm.square_network()
        first, _ = pg.gr24_square_pair()
        from_bipartite = Graph.from_planar_bipartite_network(network)
        from_orientation = Graph.from_planar_network(network.to_perfect_orientation())
        assert from_bipartite == first
        assert from_orientation.trip_permutation == first.trip_permutation
        assert frozenset(from_orientation.face_labels()) == frozenset(
            first.face_labels()
        )

    def test_bipartite_network_trips_agree(self):
        for network in (bm.square_network(), bm.lollipop_network()):
            graph = Graph.from_planar_bipartite_network(network)
            assert graph.trip_permutation == network.trip_permutation()
            assert graph.to_positroid() == network.to_positroid()

    @given(reduced_graphs(), st.integers(-7, 7))
    @settings(max_examples=100, deadline=None)
    def test_cyclic_shift_is_invertible_and_conjugates_the_permutation(
        self, graph, steps
    ):
        shifted = graph.cyclic_shift(steps)
        assert shifted.cyclic_shift(-steps) == graph
        assert shifted.to_decorated_permutation() == (
            graph.to_decorated_permutation().cyclic_shift(-steps)
        )

    @given(squared_graphs())
    @settings(max_examples=100, deadline=None)
    def test_square_move_is_an_involution(self, graph_face):
        graph, face = graph_face
        assert graph.square_move(face).square_move(face) == graph

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_middle_vertex_insertion_then_removal_restores_the_graph(self, graph, data):
        assume(graph.edges)
        edge = data.draw(st.integers(0, len(graph.edges) - 1))
        back = graph.insert_middle_vertex(edge, BLACK, "middle").remove_middle_vertex(
            "middle"
        )
        assert back.boundary == graph.boundary
        assert back.colors == graph.colors
        assert _neighbor_rotations(back) == _neighbor_rotations(graph)
        assert back.trip_permutation == graph.trip_permutation

    @given(leafless_graphs(), st.data())
    @settings(max_examples=100, deadline=None)
    def test_uncontraction_then_contraction_restores_the_graph(self, graph, data):
        big = [v for v in graph.internal_vertices if graph.degree(v) >= 3]
        assume(big)
        split = _split_somewhere(data.draw, graph, data.draw(st.sampled_from(big)))
        assume(split != graph)
        back = split.contract_edge(len(split.edges) - 1)
        assert back.boundary == graph.boundary
        assert dict(back.colors) == dict(graph.colors)
        assert _neighbor_rotations(back) == _neighbor_rotations(graph)
        assert back.trip_permutation == graph.trip_permutation

    @given(leafless_graphs())
    @settings(max_examples=100, deadline=None)
    def test_normalization_is_idempotent_and_leaves_nothing_to_contract(self, graph):
        normal = graph.normalized()
        assert normal.normalized() == normal
        assert all(
            normal.degree(v) != 2
            or len({d >> 1 for d in dict(normal.rotations)[v]}) == 1
            for v in normal.internal_vertices
        )
        assert all(
            u == w or normal.color(u) != normal.color(w)
            for u, w in normal.edges
            if u in normal.internal_vertices and w in normal.internal_vertices
        )

    def test_repr_is_compact(self):
        assert repr(pg.hollow_digon()) == (
            "PlabicGraph(n=2, edges=4, black=1, white=1, trips=(1, 2))"
        )


# --------------------------------------------------------------------------- #
# Axiom violations — each numbered axiom rejected by name
# --------------------------------------------------------------------------- #
class TestAxiomViolations:
    def test_p1_boundary_vertices_must_be_distinct(self):
        with pytest.raises(ValueError, match=r"\(P1\).*distinct"):
            _graph((1, 1), [(1, "v")], {"v": [0]}, {"v": BLACK})

    def test_p1_rotation_must_list_edges_at_the_vertex(self):
        with pytest.raises(ValueError, match=r"\(P1\).*not an edge at that vertex"):
            _graph(
                (1, 2),
                [(1, "v"), (2, "w")],
                {"v": [0, 1], "w": [1]},
                {"v": BLACK, "w": WHITE},
            )

    def test_p1_rotation_must_list_every_half_edge(self):
        with pytest.raises(ValueError, match=r"\(P1\).*every half-edge"):
            _graph((1, 2), [(1, "v"), ("v", 2)], {"v": [0]}, {"v": BLACK})

    def test_p1_graph_must_be_drawn_in_the_disk(self):
        first, _ = pg.gr24_square_pair()
        twisted: dict[Hashable, list[int]] = {
            "T": [0, 5, 4],
            "R": [1, 5, 6],
            "B": [6, 7, 2],
            "L": [4, 3, 7],
        }
        with pytest.raises(ValueError, match=r"\(P1\).*planar drawing in the disk"):
            _graph(first.boundary, first.edges, twisted, dict(first.colors))

    def test_p1_boundary_must_be_labelled_clockwise(self):
        first, _ = pg.gr24_square_pair()
        rotations: dict[Hashable, list[int]] = {
            "T": [0, 4, 5],
            "R": [1, 5, 6],
            "B": [6, 7, 2],
            "L": [4, 3, 7],
        }
        with pytest.raises(ValueError, match=r"\(P1\).*labelled clockwise"):
            _graph((1, 3, 2, 4), first.edges, rotations, dict(first.colors))

    @pytest.mark.parametrize(
        "edges", [[(1, "v"), (1, "v"), ("v", 2)], [(2, "v")]], ids=["two", "none"]
    )
    def test_p2_boundary_vertex_is_incident_to_a_single_edge(self, edges):
        rotations: dict[Hashable, list[int]] = {"v": list(range(len(edges)))}
        with pytest.raises(ValueError, match=r"\(P2\).*single edge"):
            _graph((1, 2), edges, rotations, {"v": BLACK})

    def test_p3_color_must_be_black_or_white(self):
        with pytest.raises(ValueError, match=r"\(P3\).*not in \{1, -1\}"):
            _graph((1,), [(1, "v")], {"v": [0]}, {"v": 0})

    def test_p3_every_internal_vertex_carries_a_color(self):
        with pytest.raises(ValueError, match=r"\(P3\).*carries no color"):
            _graph((1,), [(1, "v")], {"v": [0]}, {})

    def test_p3_boundary_vertices_are_not_colored(self):
        with pytest.raises(ValueError, match=r"\(P3\).*only internal vertices"):
            _graph((1,), [(1, "v")], {"v": [0]}, {"v": BLACK, 1: WHITE})

    def test_fwz_every_internal_vertex_is_connected_to_the_boundary(self):
        with pytest.raises(ValueError, match=r"FWZ Definition 7\.1\.1"):
            _graph(
                (1,),
                [(1, "v"), ("x", "y")],
                {"v": [0], "x": [1], "y": [1]},
                {"v": BLACK, "x": BLACK, "y": WHITE},
            )


class TestOperationPreconditions:
    def test_square_move_needs_trivalent_alternating_square(self):
        graph = pg.hollow_digon()
        with pytest.raises(ValueError, match=r"\(M1\)"):
            graph.square_move(2)
        with pytest.raises(ValueError, match="touches the boundary"):
            graph.square_move(0)
        with pytest.raises(ValueError, match="out of range"):
            graph.square_move(9)

    def test_contraction_needs_a_unicolored_internal_edge(self):
        with pytest.raises(ValueError, match=r"\(M2\)"):
            pg.hollow_digon().contract_edge(1)
        with pytest.raises(ValueError, match=r"\(M2\)"):
            pg.hollow_digon().contract_edge(0)

    def test_uncontraction_needs_consecutive_edges_and_a_fresh_label(self):
        graph = pg.top_cell_graph(2, 5).normalized()
        big = max(graph.internal_vertices, key=graph.degree)
        darts = dict(graph.rotations)[big]
        assert len(darts) >= 4
        with pytest.raises(ValueError, match=r"\(M2\).*consecutive"):
            graph.uncontract_vertex(big, [darts[0] >> 1, darts[2] >> 1], "fresh")
        with pytest.raises(ValueError, match="already in use"):
            graph.uncontract_vertex(big, [darts[0] >> 1], big)
        with pytest.raises(ValueError, match=r"\(M2\).*internal vertex"):
            graph.uncontract_vertex(1, [0], "fresh")

    def test_middle_vertex_moves_check_their_hypotheses(self):
        graph = pg.hollow_digon()
        with pytest.raises(ValueError, match=r"\(M3\).*degree two"):
            graph.remove_middle_vertex("u")
        with pytest.raises(ValueError, match=r"\(M3\).*black \(1\) or white \(-1\)"):
            graph.insert_middle_vertex(0, 2, "m")

    def test_parallel_edge_reduction_needs_a_digon(self):
        first, _ = pg.gr24_square_pair()
        with pytest.raises(ValueError, match=r"\(R1\)"):
            first.parallel_edge_reduction(first.face_count - 1)

    def test_leaf_reduction_checks_its_hypotheses(self):
        graph = _leaf_fixture()
        with pytest.raises(ValueError, match=r"\(R2\).*degree-one"):
            graph.leaf_reduction("v", [])
        with pytest.raises(ValueError, match=r"\(R2\).*2 endpoints"):
            graph.leaf_reduction("u", ["x"])
        with pytest.raises(ValueError, match=r"\(R2\).*opposite color"):
            pg.lollipop_graph([BLACK]).leaf_reduction((1, "lollipop"), [])

    def test_square_move_class_is_guarded_and_needs_a_reduced_graph(self):
        with pytest.raises(ValueError, match="exceeds the guard of 1"):
            pg.top_cell_graph(2, 4).square_move_class(limit=1)
        with pytest.raises(ValueError, match="reduced"):
            pg.hollow_digon().square_move_class()

    def test_matchings_need_a_bipartite_graph(self):
        graph = pg.nonreduced_wiring_example()
        assert graph.is_bipartite is False
        with pytest.raises(ValueError, match="bipartite"):
            graph.almost_perfect_matchings()
        with pytest.raises(ValueError, match="two boundary vertices"):
            Graph.from_wiring_diagram([], 2).almost_perfect_matchings()

    def test_constructors_reject_bad_input(self):
        with pytest.raises(ValueError, match="letters"):
            Graph.from_wiring_diagram([2], 2)
        with pytest.raises(ValueError, match="Le-condition"):
            Graph.from_le_diagram([[0, 1], [1, 0]], 4)
        with pytest.raises(ValueError, match="0 <= k <= n"):
            pg.top_cell_graph(3, 2)
        with pytest.raises(ValueError, match="perfect orientation"):
            Graph.from_planar_network(
                bm.PlanarNetwork.from_le_diagram([[1, 1], [1, 1]], 4)
            )

    def test_from_dataframe_rejects_malformed_frames(self):
        frame = pg.hollow_digon().to_dataframe()
        with pytest.raises(ValueError, match="missing required columns"):
            Graph.from_dataframe(frame.drop(columns=["slot"]))
        with pytest.raises(ValueError, match="unknown vertex kind"):
            Graph.from_dataframe(frame.assign(kind="grey"))
        with pytest.raises(ValueError, match="two half-edge rows"):
            Graph.from_dataframe(frame.iloc[1:])


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestPlots:
    @pytest.mark.parametrize("method", ["plot_graph", "plot_trips"])
    @pytest.mark.parametrize(
        "graph",
        [pg.hollow_digon(), pg.top_cell_graph(2, 5), pg.lollipop_graph([BLACK, WHITE])],
    )
    def test_plot_draws_onto_the_given_axes(self, method, graph):
        figure, ax = plt.subplots()
        try:
            assert getattr(graph, method)(ax) is ax
            assert ax.lines
        finally:
            plt.close(figure)

    def test_plot_creates_axes_when_none_is_given(self):
        ax = pg.gr24_square_pair()[0].plot_graph()
        try:
            assert ax.get_aspect() == 1.0
        finally:
            plt.close("all")
