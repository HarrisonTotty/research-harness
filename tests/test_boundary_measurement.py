"""Tests for research.boundary_measurement, from the Boundary Measurement Map page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Postnikov's
nonnegativity, fixed-graph, and Le-diagram theorems, Talaska's flow
formula, Lam's matching/orientation bridges, moves, bridges, and the
Muller-Speyer twist); round-trip laws come from the API contract; and every
named condition has a rejection test naming it. Hypothesis strategies draw
random positive rational weights on the page's fixture graph families and
random Le-diagrams — generating arbitrary embedded planar networks is out
of scope, so property coverage is over these families.
"""

import itertools
from collections.abc import Hashable, Iterable, Mapping
from fractions import Fraction

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from research import boundary_measurement as bmm
from research import positroid as ps
from research._linalg import det_q
from research.grassmann_necklace import GrassmannNecklace

matplotlib.use("Agg")


# --------------------------------------------------------------------------- #
# Strategies: random positive weights on the page's fixture families
# --------------------------------------------------------------------------- #
positive_fractions = st.builds(Fraction, st.integers(1, 9), st.integers(1, 9))

square_weights = st.tuples(
    positive_fractions, positive_fractions, positive_fractions, positive_fractions
)


def _satisfies_le(filling: list[list[int]]) -> bool:
    boxes = {
        (i, j): value
        for i, row in enumerate(filling, start=1)
        for j, value in enumerate(row, start=1)
    }
    return not any(
        a and c and boxes.get((i2, j2)) == 0
        for (i2, j1), a in boxes.items()
        for (i1, j2), c in boxes.items()
        if i1 < i2 and j1 < j2
    )


@st.composite
def le_diagrams(draw: st.DrawFn) -> tuple[list[list[int]], int]:
    """Draw a valid Le-diagram filling with its rectangle size ``n``."""
    k = draw(st.integers(1, 2))
    width = draw(st.integers(1, 3))
    shape = sorted((draw(st.integers(0, width)) for _ in range(k)), reverse=True)
    filling = [[draw(st.integers(0, 1)) for _ in range(row)] for row in shape]
    assume(_satisfies_le(filling))
    return filling, k + width


def _pluck[T: Hashable](
    columns: Mapping[T, tuple[Fraction, ...]],
    subset: Iterable[T],
) -> Fraction:
    """Return the sorted-column minor of a label-to-column mapping."""
    chosen = [label for label in columns if label in set(subset)]
    height = len(next(iter(columns.values())))
    return det_q([[columns[label][r] for label in chosen] for r in range(height)])


def _projectively_equal[K](
    first: Mapping[K, Fraction],
    second: Mapping[K, Fraction],
    keys: Iterable[K],
) -> bool:
    """Whether two Plucker vectors agree up to one global scalar."""
    ratio = None
    for key in keys:
        a, b = first.get(key, Fraction(0)), second.get(key, Fraction(0))
        if (a == 0) != (b == 0):
            return False
        if a != 0:
            if ratio is None:
                ratio = a / b
            elif a / b != ratio:
                return False
    return True


def _plucker_vector(
    network: bmm.PlanarBipartiteNetwork[int | str]
    | bmm.PlanarBipartiteNetwork[int | str | tuple[int | str, str]],
    keys: Iterable[frozenset[int]],
) -> dict[frozenset[int], Fraction]:
    """Restrict a network's Plucker mapping to the given boundary subsets."""
    values = network.pluckers()
    return {key: values.get(key, Fraction(0)) for key in keys}


def _oriented_square(
    a: Fraction | int, b: Fraction | int, c: Fraction | int, d: Fraction | int
) -> tuple[bmm.PlanarBipartiteNetwork[int | str], bmm.PlanarNetwork[int | str]]:
    """The square network oriented from its all-legs matching (I = {2, 4})."""
    square = bmm.square_network(a, b, c, d)
    matching = next(
        m
        for m in square.almost_perfect_matchings()
        if square.boundary_subset(m) == frozenset({2, 4})
    )
    return square, square.to_perfect_orientation(matching)


PAIRS = [frozenset(pair) for pair in itertools.combinations([1, 2, 3, 4], 2)]


_SQUARE_SIDES = [("L", "T"), ("T", "R"), ("R", "B"), ("B", "L")]


def _square_minus_side(
    side: tuple[str, str],
) -> bmm.PlanarBipartiteNetwork[int | str]:
    """The square fixture with one square-side edge deleted.

    Still reduced, but its trip permutation is no longer an involution —
    the square's own ``(3, 4, 1, 2)`` is, which makes it blind to a
    reversal of Lam's turning rule.
    """
    every: list[tuple[int | str, int | str, Fraction | int]] = [
        (1, "T", 1),
        (2, "R", 1),
        (3, "B", 1),
        (4, "L", 1),
        ("L", "T", 2),
        ("T", "R", 3),
        ("R", "B", 5),
        ("B", "L", 7),
    ]
    return bmm.PlanarBipartiteNetwork.from_edges(
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
        [edge for edge in every if {edge[0], edge[1]} != set(side)],
        black_vertices=["T", "B", 2, 4],
    )


def _inverse_permutation(targets: tuple[int, ...]) -> tuple[int, ...]:
    """Return the inverse of a permutation given in one-line notation."""
    out = [0] * len(targets)
    for position, value in enumerate(targets, start=1):
        out[value - 1] = position
    return tuple(out)


def _stacked_squares() -> bmm.PlanarBipartiteNetwork[str]:
    """Two stacked squares on four legs — the edge-repeat reducedness witness.

    The only fixture on which Lam's middle §7.1 clause decides alone: no
    trip is closed and no two trips share two edges in the same order, so
    the verdict rests entirely on the trips that use an edge twice away
    from a boundary leaf.
    """
    return bmm.PlanarBipartiteNetwork.from_edges(
        ("bB", "bC", "bD", "bA"),
        {
            "p00": (0, 0),
            "p01": (0, 3),
            "p02": (0, 6),
            "p10": (3, 0),
            "p11": (3, 3),
            "p12": (3, 6),
            "bA": (-4, 0),
            "bB": (-4, 6),
            "bC": (0, 10),
            "bD": (7, 3),
        },
        [
            ("p00", "p10", 1),
            ("p01", "p11", 1),
            ("p02", "p12", 1),
            ("p00", "p01", 1),
            ("p01", "p02", 1),
            ("p10", "p11", 1),
            ("p11", "p12", 1),
            ("bA", "p00", 1),
            ("bB", "p02", 1),
            ("bC", "p02", 1),
            ("bD", "p11", 1),
        ],
        black_vertices=["p00", "p02", "p11"],
    )


def _trips_share_two_edges_in_order(
    network: bmm.PlanarBipartiteNetwork[str],
) -> bool:
    """Whether two boundary trips share two edges in the same order."""
    boundary_trips = network.trips()[: len(network.boundary)]
    for index, first in enumerate(boundary_trips):
        for second in boundary_trips[index + 1 :]:
            one = [edge for _, edge in first]
            two = [edge for _, edge in second]
            common = sorted(set(one) & set(two))
            for a in common:
                for b in common:
                    if a < b and (one.index(a) < one.index(b)) == (
                        two.index(a) < two.index(b)
                    ):
                        return True
    return False


def _two_boundary_square() -> bmm.PlanarBipartiteNetwork[int | str]:
    """The square with only two boundary legs — Lam-non-reduced fixture."""
    return bmm.PlanarBipartiteNetwork.from_edges(
        (1, 3),
        {1: (0, 2), 3: (0, -2), "T": (0, 1), "R": (1, 0), "B": (0, -1), "L": (-1, 0)},
        [
            (1, "T", 1),
            (3, "B", 1),
            ("L", "T", 2),
            ("T", "R", 3),
            ("R", "B", 5),
            ("B", "L", 7),
        ],
        black_vertices=["T", "B"],
    )


# --------------------------------------------------------------------------- #
# Canonical examples: each asserts what the page says the example certifies
# --------------------------------------------------------------------------- #
class TestSquareGraphFixture:
    """Lam Example 4.3 — the workhorse fixture."""

    def test_the_six_plucker_coordinates(self):
        square = bmm.square_network(2, 3, 5, 7)
        a, b, c, d = Fraction(2), Fraction(3), Fraction(5), Fraction(7)
        assert square.pluckers() == {
            frozenset({1, 2}): a,
            frozenset({1, 3}): a * c + b * d,
            frozenset({1, 4}): b,
            frozenset({2, 3}): d,
            frozenset({2, 4}): Fraction(1),
            frozenset({3, 4}): c,
        }

    @given(square_weights)
    @settings(max_examples=50)
    def test_certifies_the_plucker_relation(self, weights):
        square = bmm.square_network(*weights)
        p = square.pluckers()

        def delta(i: int, j: int) -> Fraction:
            return p.get(frozenset({i, j}), Fraction(0))

        assert (
            delta(1, 2) * delta(3, 4)
            - delta(1, 3) * delta(2, 4)
            + delta(1, 4) * delta(2, 3)
            == 0
        )

    @given(square_weights)
    @settings(max_examples=25)
    def test_realizing_matrix_has_identity_in_source_columns(self, weights):
        a, b, c, d = weights
        _, net = _oriented_square(a, b, c, d)
        assert net.to_matrix() == {
            1: (b, -a),
            2: (Fraction(1), Fraction(0)),
            3: (c, d),
            4: (Fraction(0), Fraction(1)),
        }

    @pytest.mark.parametrize(
        "weights",
        [
            (1, 1, 1, 1),
            (2, 3, 5, 7),
            (Fraction(1, 3), Fraction(2, 5), Fraction(7, 2), Fraction(9, 4)),
        ],
    )
    def test_lands_in_the_top_cell_with_trip_permutation(self, weights):
        square, net = _oriented_square(*weights)
        positroid = net.to_positroid()
        assert positroid.bases == set(PAIRS)
        decorated = positroid.to_decorated_permutation()
        assert decorated.targets == (3, 4, 1, 2)
        assert decorated == ps.uniform_positroid(2, 4).to_decorated_permutation()
        assert square.trip_permutation() == (3, 4, 1, 2)

    def test_degenerations_are_outside_the_map_domain(self):
        # The page's caveat: zero weights are not values of Meas — the
        # constructor rejects them, and the degenerate positroids are
        # reached only through the matrix directly.
        with pytest.raises(ValueError, match="strictly positive weights"):
            bmm.square_network(0, 1, 1, 1)
        b, c, d = Fraction(1), Fraction(1), Fraction(1)
        five = ps.Positroid.from_matrix({1: (b, 0), 2: (1, 0), 3: (c, d), 4: (0, 1)})
        assert five.bases == {
            frozenset(s) for s in [{1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}]
        }
        a = Fraction(1)
        four = ps.Positroid.from_matrix({1: (0, -a), 2: (1, 0), 3: (1, 0), 4: (0, 1)})
        assert four.bases == {frozenset(s) for s in [{1, 2}, {1, 3}, {2, 4}, {3, 4}]}


class TestOtherCanonicalExamples:
    """The page's remaining fixtures, each with its source attribution."""

    def test_lollipop_certifies_zero_dimensional_cells_and_not_used_clause(self):
        """Lam Example 4.2 (2015): the smallest graph, a torus-fixed point."""
        lollipop = bmm.lollipop_network()
        matchings = lollipop.almost_perfect_matchings()
        assert len(matchings) == 1
        assert len(matchings[0]) == 4
        assert lollipop.boundary_subset(matchings[0]) == frozenset({3, 4})
        positroid = lollipop.to_positroid()
        assert positroid.bases == {frozenset({3, 4})}
        assert bmm.cell_dimension(positroid) == 0

    def test_geometric_series_needs_the_subtraction_free_form(self):
        """Postnikov Example 4.5 (2006): M_12 = 1/2 at all weights one.

        The smallest witness that a cyclic network needs the
        subtraction-free rational form — the signed series is infinite and
        alternating, the answer a finite positive rational.
        """
        network = bmm.geometric_series_network()
        assert not network.is_acyclic
        assert network.boundary_measurement(1, 2) == Fraction(1, 2)

    @given(
        positive_fractions, positive_fractions, positive_fractions, positive_fractions
    )
    @settings(max_examples=50)
    def test_geometric_series_measurement_is_xyt_over_one_plus_yz(self, x, y, z, t):
        """Postnikov Example 4.5 (2006): M_12 = xyt / (1 + yz)."""
        network = bmm.geometric_series_network(x, y, z, t)
        assert network.boundary_measurement(1, 2) == x * y * t / (1 + y * z)

    def test_geometric_series_conservative_flows_are_empty_and_cycle(self):
        """Postnikov Example 4.5 (2006): the denominator is 1 + yz.

        Its two conservative flows (Talaska Definition 3.1) are the empty
        one and the single three-edge cycle.
        """
        network = bmm.geometric_series_network()
        flows = network.conservative_flows()
        assert sorted(len(flow) for flow in flows) == [0, 3]

    def test_acyclic_baseline_has_the_lindstrom_path_matrix(self):
        """Lam section 2.3 (2015): the acyclic baseline, no denominator."""
        a, b, c = Fraction(2), Fraction(3), Fraction(5)
        network = bmm.acyclic_baseline_network(a, b, c)
        assert network.is_acyclic
        measurements = network.boundary_measurements()
        assert measurements == {
            (1, "1'"): 1 + a * c,
            (1, "2'"): a,
            (1, "3'"): Fraction(0),
            (2, "1'"): c,
            (2, "2'"): Fraction(1),
            (2, "3'"): Fraction(0),
            (3, "1'"): b * c,
            (3, "2'"): b,
            (3, "3'"): Fraction(1),
        }

    def test_planarity_is_necessary_so_crossing_chords_are_rejected(self):
        # Talaska Example 5.4's lesson: the flow formula lives on planar
        # networks only, so the embedding validation refuses crossings.
        with pytest.raises(ValueError, match="planarity is essential"):
            bmm.PlanarNetwork.from_edges(
                (1, 2, 3, 4),
                {1: (0, 2), 2: (2, 0), 3: (0, -2), 4: (-2, 0)},
                [(1, 3, 1), (2, 4, 1)],
            )

    @given(
        positive_fractions, positive_fractions, positive_fractions, positive_fractions
    )
    @settings(max_examples=50)
    def test_talaska_conservative_flow_factorization(self, w, y, z, t):
        # Talaska Examples 2.5-3.3, as far as the page records them: the
        # ten pairwise-disjoint sub-collections sum to the stated
        # factorization, and only the (1 + Z) factor cancels.
        terms = [
            Fraction(1),
            w,
            y,
            z,
            t,
            w * z,
            w * y,
            y * z,
            w * y * z,
            z * t,
        ]
        factored = (1 + z) * ((1 + w) * (1 + y) + t)
        assert sum(terms) == factored
        # The (1 + Z) cancellation and the surviving denominator are read
        # off this factorization; her network itself is not on the page.


class TestMeasurementMatrix:
    """Postnikov Definition 4.6 and the two sign fixtures."""

    def test_postnikov_example_4_7_negates_exactly_one_entry(self):
        m12, m14, m32, m34 = (
            Fraction(2),
            Fraction(3),
            Fraction(5),
            Fraction(7),
        )
        columns = bmm.measurement_matrix(
            [1, 3],
            4,
            {(1, 2): m12, (1, 4): m14, (3, 2): m32, (3, 4): m34},
        )
        assert columns == (
            (Fraction(1), Fraction(0)),
            (m12, m32),
            (Fraction(0), Fraction(1)),
            (-m14, m34),
        )

    def test_talaska_example_2_7_negates_exactly_one_entry_at_n_5(self):
        values = {
            (1, 2): Fraction(2),
            (1, 3): Fraction(3),
            (1, 5): Fraction(5),
            (4, 2): Fraction(7),
            (4, 3): Fraction(11),
            (4, 5): Fraction(13),
        }
        columns = bmm.measurement_matrix([1, 4], 5, values)
        assert columns == (
            (Fraction(1), Fraction(0)),
            (values[(1, 2)], values[(4, 2)]),
            (values[(1, 3)], values[(4, 3)]),
            (Fraction(0), Fraction(1)),
            (-values[(1, 5)], values[(4, 5)]),
        )

    @given(
        positive_fractions, positive_fractions, positive_fractions, positive_fractions
    )
    @settings(max_examples=50)
    def test_companion_minor_is_postnikov_proposition_5_2(self, m12, m14, m32, m34):
        columns = bmm.measurement_matrix(
            [1, 3],
            4,
            {(1, 2): m12, (1, 4): m14, (3, 2): m32, (3, 4): m34},
        )
        minor = det_q(
            [
                [columns[1][0], columns[3][0]],
                [columns[1][1], columns[3][1]],
            ]
        )
        assert minor == m12 * m34 + m14 * m32


# --------------------------------------------------------------------------- #
# Property tests transcribed from the structural theorems
# --------------------------------------------------------------------------- #
class TestStructuralTheorems:
    @given(le_diagrams(), square_weights)
    @settings(max_examples=30, deadline=None)
    def test_nonnegativity_postnikov_corollary_5_4(self, diagram, weights):
        """Postnikov Corollary 5.4 (2006): the image is in Gr_kn^tnn."""
        filling, n = diagram
        le_network = bmm.PlanarNetwork.from_le_diagram(filling, n)
        assert le_network.to_positroid().rank() == len(filling)
        _, net = _oriented_square(*weights)
        assert net.to_positroid().rank() == 2

    @given(le_diagrams(), st.data())
    @settings(max_examples=25, deadline=None)
    def test_fixed_graph_hits_one_cell_postnikov_theorem_4_10(self, diagram, data):
        """Postnikov Theorem 4.10 (2006): Meas_G maps into a single cell."""
        filling, n = diagram
        ones = [
            (i, j)
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
            if value
        ]
        first = data.draw(st.tuples(*(positive_fractions for _ in ones)), label="first")
        second = data.draw(
            st.tuples(*(positive_fractions for _ in ones)), label="second"
        )
        cells = [
            bmm.PlanarNetwork.from_le_diagram(
                filling, n, dict(zip(ones, draw, strict=True))
            ).to_positroid()
            for draw in (first, second)
        ]
        assert cells[0] == cells[1]

    @given(le_diagrams())
    @settings(max_examples=30, deadline=None)
    def test_le_diagram_cell_dimension_postnikov_theorem_6_5(self, diagram):
        """Postnikov Theorem 6.5 (2006): dim S = |D|, the number of 1s."""
        filling, n = diagram
        network = bmm.PlanarNetwork.from_le_diagram(filling, n)
        ones = sum(sum(row) for row in filling)
        assert bmm.cell_dimension(network.to_positroid()) == ones

    def test_every_cell_is_hit_postnikov_theorems_4_8_and_6_5(self):
        """Postnikov Theorems 4.8 and 6.5 (2006), exhaustively at n = 4.

        The image of the boundary measurement map is the whole totally
        nonnegative Grassmannian: sweeping every Le-diagram in every
        ``k x (n - k)`` rectangle, the Le-graph networks realize each of
        the 65 positroids on ``[4]`` exactly once — the surjectivity of
        Theorem 4.8 through the Le-diagram/cell bijection of Theorem 6.5.
        """
        n = 4
        realized = []
        for k in range(n + 1):
            width = n - k
            shapes = [
                shape
                for shape in itertools.product(range(width + 1), repeat=k)
                if all(x >= y for x, y in itertools.pairwise(shape))
            ]
            for shape in shapes:
                for bits in itertools.product((0, 1), repeat=sum(shape)):
                    box = iter(bits)
                    filling = [[next(box) for _ in range(row)] for row in shape]
                    if not _satisfies_le(filling):
                        continue
                    realized.append(
                        bmm.PlanarNetwork.from_le_diagram(filling, n).to_positroid()
                    )
        assert len(set(realized)) == len(realized)
        assert set(realized) == set(ps.enumerate_positroids(n))

    @given(le_diagrams(), st.data())
    @settings(max_examples=25, deadline=None)
    def test_acyclic_le_graph_is_i_polynomial_postnikov_theorem_4_11(
        self, diagram, data
    ):
        """Postnikov Theorem 4.11 (2006), the acyclic clause at i = 1.

        The theorem asserts, for each ``i``, an acyclic graph with source
        set ``I_i`` on which the map is an ``I_i``-polynomial
        parameterization. The Le-graph realizes the ``i = 1`` instance:
        it is acyclic and its source set is the lexicographically minimal
        basis ``I_1`` (the first Grassmann necklace entry). Two
        discriminating consequences are checked — the parameterization is
        injective, so distinct tableaux give distinct points (the
        "isomorphism" of Definition 4.9), and raising one tableau entry
        cannot decrease any ``Delta_J``, as nonnegative coefficients
        require. (``Delta_{I_1} = 1`` and integrality at integer weights
        are *not* asserted: the first holds for every network by the
        identity block of Definition 4.6, the second for every acyclic
        one, so neither discriminates.) Source sets ``I_i`` for ``i > 1``
        need the cyclically rotated diagram, i.e. the
        positroid-to-Le-diagram direction (backlog).
        """
        filling, n = diagram
        ones = [
            (i, j)
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
            if value
        ]
        assume(ones)
        weights = {box: data.draw(st.integers(1, 6), label=str(box)) for box in ones}
        network = bmm.PlanarNetwork.from_le_diagram(filling, n, weights)
        positroid = network.to_positroid()
        first_entry = GrassmannNecklace.from_matroid(positroid).entries[0]
        assert network.is_acyclic
        assert frozenset(network.source_set) == first_entry

        def point(
            net: bmm.PlanarNetwork[int | tuple[int, int]],
        ) -> dict[frozenset[int | tuple[int, int]], Fraction]:
            base = min(positroid.bases, key=lambda s: sorted(map(str, s)))
            scale = net.plucker(base)
            return {basis: net.plucker(basis) / scale for basis in positroid.bases}

        raised = dict(weights)
        raised[ones[0]] += 1
        grown = bmm.PlanarNetwork.from_le_diagram(filling, n, raised)
        assert point(grown) != point(network)
        for basis in positroid.bases:
            assert grown.plucker(basis) >= network.plucker(basis)

    @pytest.mark.parametrize(
        ("network", "expected"),
        [
            (bmm.square_network(2, 3, 5, 7), set(PAIRS)),
            (bmm.lollipop_network(), {frozenset({3, 4})}),
        ],
    )
    def test_graph_matroid_is_the_image_cell_postnikov_prop_11_7(
        self, network, expected
    ):
        """Postnikov Proposition 11.7 and Lemma 11.10 (2006).

        ``M_G``, the set of source sets of the perfect orientations, is a
        totally nonnegative matroid, and the image of the map lies in its
        cell. The discriminating comparison is against the *directed*
        network's positroid, which reaches its bases through
        ``Positroid.from_matrix`` on the signed measurement matrix — a
        route independent of ``M_G``. (Comparing ``M_G`` against the
        bipartite ``to_positroid()`` would be definitional: that method is
        ``from_bases`` over the dimer partition function's support, which
        is ``{I(Pi)}`` by construction.) The lollipop is the non-uniform
        instance: its ``M_G`` is the single basis ``{3, 4}``, so which
        subsets *fail* to be bases is exercised — on the top-cell square
        alone every 2-subset is a basis and nothing discriminates.
        """
        matchings = network.almost_perfect_matchings()
        graph_matroid = {
            frozenset(network.to_perfect_orientation(m).source_set) for m in matchings
        }
        assert graph_matroid == expected
        for matching in matchings:
            directed = network.to_perfect_orientation(matching)
            assert directed.to_positroid().bases == graph_matroid

    def test_cell_closures_have_euler_characteristic_one_psw_corollary_6_3(self):
        """Postnikov-Speyer-Williams Corollary 6.3 (2009), at n = 4.

        Each cell closure of the totally nonnegative Grassmannian has
        Euler characteristic 1 — the computable consequence of
        Postnikov-Speyer-Williams Theorem 6.2 (the decomposition is a
        finite CW complex, itself a topological statement no finite test
        can carry). The closure of a cell is the set of cells of the same
        rank whose necklace is above it in the order
        :class:`research.grassmann_necklace.GrassmannNecklace` implements,
        which that module attributes to Lam's *Totally nonnegative
        Grassmannian and Grassmann polytopes* Theorem 6.2 — a different
        paper and result from the PSW Theorem 6.2 above, which shares its
        number.
        """
        n = 4
        cells = ps.enumerate_positroids(n)
        necklaces = {cell: GrassmannNecklace.from_matroid(cell) for cell in cells}
        dimensions = {cell: bmm.cell_dimension(cell) for cell in cells}
        for cell in cells:
            closure = [
                other
                for other in cells
                if other.rank() == cell.rank() and necklaces[cell] <= necklaces[other]
            ]
            assert sum((-1) ** dimensions[other] for other in closure) == 1

    @pytest.mark.parametrize(
        "network",
        [
            bmm.square_network(2, 3, 5, 7),
            bmm.lollipop_network(),
            _two_boundary_square(),
        ],
    )
    def test_type_from_colours_and_degrees_lam_eq_20(self, network):
        """Lam eq. (20) (2015): k = (n + sum_b(deg-2) + sum_w(2-deg)) / 2.

        The page records Postnikov's Definition 11.5 as printing this
        type-(k, n) condition with ``k + (n - k)`` where his own Lemma 9.4
        has ``k - (n - k)``, and settles the erratum in favour of Lam's
        equation — so Lam's is what is checked here, against the rank of
        the cell the network maps into.
        """
        degree = dict.fromkeys((v for v, _ in network.positions), 0)
        for u, v, _ in network.edges:
            degree[u] += 1
            degree[v] += 1
        signed = sum(
            degree[v] - 2 if v in network.black_vertices else 2 - degree[v]
            for v in network.internal_vertices
        )
        k, remainder = divmod(len(network.boundary) + signed, 2)
        assert remainder == 0
        assert k == network.to_positroid().rank()

    @pytest.mark.parametrize(
        "network",
        [
            _oriented_square(2, 3, 5, 7)[1],
            bmm.geometric_series_network(2, 3, 5, 7),
            bmm.acyclic_baseline_network(2, 3, 5),
            bmm.PlanarNetwork.from_le_diagram([[1, 1], [1, 1]], 4),
            bmm.PlanarNetwork.from_le_diagram([[1, 0], [1]], 4),
        ],
    )
    def test_gauge_quotient_dimension_postnikov_lemma_11_1(self, network):
        """Postnikov Lemma 11.1 (2006): the gauge quotient is R_{>0}^{|F|-1}.

        The parameter count of the networks modulo gauge equals both the
        edge count less the *internal* vertex count and ``|F(G)| - 1``.
        The page prints the first as ``|E| - |V|`` over all vertices,
        which is false on its own fixtures — the oriented square has
        ``|E| = |V| = 8`` against ``|F| - 1 = 4``, and the lollipop gives
        ``-4`` — so the internal-vertex reading is the one transcribed
        here; see the discrepancy note reported with this module.
        """
        internal = len(network.internal_vertices)
        assert len(network.edges) - internal == network.face_count - 1

    @given(square_weights, st.data())
    @settings(max_examples=20)
    def test_face_weights_identify_gauge_classes_postnikov_lemma_11_2(
        self, weights, data
    ):
        """Postnikov Lemma 11.2 / Lam Lemma 4.7: gauge classes = face weights.

        The identification is a bijection, so both halves must hold: a
        gauge transformation leaves the face weights alone (Lam Lemma
        4.6), and distinct gauge classes have distinct face weights.
        Scaling a single edge weight leaves the gauge class — no choice of
        vertex scalars reproduces it, since gauge acts trivially on the
        face weights — so it must move them.
        """
        _, net = _oriented_square(*weights)
        scalars = {
            v: data.draw(positive_fractions, label=str(v))
            for v in sorted(net.internal_vertices, key=str)
        }
        assert net.gauge_transform(scalars).face_weights() == net.face_weights()
        index = data.draw(st.integers(0, len(net.edges) - 1), label="edge")
        factor = data.draw(positive_fractions.filter(lambda f: f != 1), label="factor")
        edges = list(net.edges)
        tail, head, weight = edges[index]
        edges[index] = (tail, head, weight * factor)
        scaled = bmm.PlanarNetwork.from_edges(net.boundary, dict(net.positions), edges)
        # Compared as vectors, not multisets: an edge borders exactly two
        # faces, whose weights pick up `factor` and `1 / factor`, so the
        # two can swap places and leave the multiset alone.
        assert scaled.face_weights() != net.face_weights()

    @pytest.mark.parametrize(
        "network",
        [
            bmm.square_network(2, 3, 5, 7),
            bmm.lollipop_network(),
            *(_square_minus_side(side) for side in _SQUARE_SIDES),
        ],
    )
    def test_trips_double_cover_every_edge_lam_section_7_1(self, network):
        """Lam section 7.1 (2015): the trips cover each edge twice.

        Following a directed edge and turning maximally right at black and
        left at white decomposes the graph into paths and cycles covering
        each edge twice, once in each direction.
        """
        traversals: dict[int, list[object]] = {}
        for trip in network.trips():
            for tail, edge in trip:
                traversals.setdefault(edge, []).append(tail)
        assert set(traversals) == set(range(len(network.edges)))
        for tails in traversals.values():
            assert len(tails) == 2
            assert tails[0] != tails[1]

    @pytest.mark.parametrize(
        "network",
        [
            bmm.square_network(2, 3, 5, 7),
            *(_square_minus_side(side) for side in _SQUARE_SIDES),
        ],
    )
    def test_trip_permutation_indexes_the_cell_lam_corollary_7_14(self, network):
        """Lam Corollary 7.14 and Theorem 7.12(3) (2015), reduced graphs only.

        For a reduced planar bipartite graph the trip permutation indexes
        the cell the network maps onto. It matches this repository's
        decorated permutation after inversion, the two conventions running
        opposite ways. The four squares with one side deleted carry the
        weight here: reversing Lam's turning rule replaces every trip
        permutation by its inverse, and the square's own value —
        ``(3, 4, 1, 2)``, the one the page prints, checked against Lam
        directly in the companion test — is an involution, so it fixes the
        values but not the handedness. Theirs are not involutions, so they
        do.
        """
        assert network.is_reduced() is True
        decorated = network.to_positroid().to_decorated_permutation().targets
        assert network.trip_permutation() == _inverse_permutation(decorated)

    def test_trip_permutation_of_the_square_is_lams_printed_value(self):
        """Lam section 7.1 (2015): pi_G = (3, 4, 1, 2) for the square."""
        assert bmm.square_network(2, 3, 5, 7).trip_permutation() == (3, 4, 1, 2)

    @given(square_weights)
    @settings(max_examples=25)
    def test_reduced_graph_dimension_postnikov_theorem_12_7(self, weights):
        """Postnikov Theorem 12.7 (2006): dim = |F(G)| - 1 for reduced G."""
        square, net = _oriented_square(*weights)
        assert square.is_reduced() is True
        assert bmm.cell_dimension(net.to_positroid()) == net.face_count - 1

    def test_non_reduced_graph_breaks_the_dimension_count(self):
        # Postnikov Remark 12.8: without reducedness the parameterization
        # collapses; the two-boundary square hits a 0-cell from 2 faces.
        network = _two_boundary_square()
        assert network.is_reduced() is False
        directed = network.to_perfect_orientation()
        assert bmm.cell_dimension(directed.to_positroid()) == 0
        assert directed.face_count - 1 == 2

    @given(square_weights)
    @settings(max_examples=25, deadline=None)
    def test_non_reduced_graphs_still_cover_their_cell_postnikov_cor_16_5(
        self, weights
    ):
        """Postnikov Corollary 16.5 (2006), on zero-dimensional cells.

        A perfectly orientable graph, reduced or not, maps onto its whole
        cell. Surjectivity onto a positive-dimensional cell is not
        certifiable by sampling, so the instances are cells that are a
        single point: every weight choice must land on it, covering the
        cell entirely. Two are needed to cover the claim — the two-boundary
        square is the *non-reduced* instance but sits in ``Gr(2, 2)``,
        which is itself a point, so the weighted lollipop supplies a
        0-cell inside ``Gr(2, 4)``, where five other bases were available
        and none is taken.
        """
        a, b, c, d = weights
        squashed = bmm.PlanarBipartiteNetwork.from_edges(
            (1, 3),
            {
                1: (0, 2),
                3: (0, -2),
                "T": (0, 1),
                "R": (1, 0),
                "B": (0, -1),
                "L": (-1, 0),
            },
            [
                (1, "T", 1),
                (3, "B", 1),
                ("L", "T", a),
                ("T", "R", b),
                ("R", "B", c),
                ("B", "L", d),
            ],
            black_vertices=["T", "B"],
        )
        assert squashed.is_reduced() is False
        assert squashed.to_positroid().bases == {frozenset({1, 3})}
        lollipop = bmm.PlanarBipartiteNetwork.from_edges(
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
            [(1, "p1", a), (2, "p2", b), (3, "p3", c), (4, "p4", d)],
            black_vertices=[3, 4, "p1", "p2"],
        )
        assert lollipop.to_positroid().bases == {frozenset({3, 4})}
        assert bmm.cell_dimension(lollipop.to_positroid()) == 0
        # A genuinely non-reduced graph in a cell with alternatives: the
        # bridged square lands on the same whole top cell as the reduced
        # square (Corollary 16.5), and its parameter count exceeds that
        # cell's dimension, so it cannot be injective (Remark 12.8).
        reduced = bmm.square_network(a, b, c, d)
        bridged = reduced.add_bridge(1, a)
        assert bridged.is_reduced() is False
        assert bridged.to_positroid().bases == reduced.to_positroid().bases
        parameters = bridged.to_perfect_orientation().face_count - 1
        assert parameters > bmm.cell_dimension(bridged.to_positroid())

    @given(square_weights, square_weights)
    @settings(max_examples=20, deadline=None)
    def test_injectivity_on_the_gauge_quotient_muller_speyer_prop_7_6(
        self, first, second
    ):
        """Muller-Speyer Propositions 5.14 and 7.6 (2017): open immersion.

        For a reduced graph the map embeds the gauge quotient into the
        open positroid variety — in particular it is injective there. Face
        weights coordinatize the gauge quotient (Postnikov Lemma 11.2;
        Muller-Speyer Proposition 5.5), so networks with distinct face
        weights must map to distinct points.
        """
        _, one = _oriented_square(*first)
        _, other = _oriented_square(*second)
        assume(one.face_weights() != other.face_weights())
        assert not _projectively_equal(
            {subset: one.plucker(subset) for subset in PAIRS},
            {subset: other.plucker(subset) for subset in PAIRS},
            PAIRS,
        )

    @given(square_weights)
    @settings(max_examples=25)
    def test_flow_formula_talaska_theorem_3_2(self, weights):
        """Talaska Theorem 3.2 (2008): Delta_J = flows(J) / conservative."""
        _, net = _oriented_square(*weights)
        for subset in PAIRS:
            assert net.plucker_via_flows(subset) == net.plucker(subset)

    @given(
        positive_fractions, positive_fractions, positive_fractions, positive_fractions
    )
    @settings(max_examples=25)
    def test_flow_formula_on_the_cyclic_fixture(self, x, y, z, t):
        """Talaska Theorem 3.2 (2008), against the page's summed series.

        On a cyclic network the matrix route is itself flow-based, so the
        independent oracle is Postnikov Example 4.5's closed form
        ``M_12 = xyt / (1 + yz)`` — the value the signed series sums to.
        """
        network = bmm.geometric_series_network(x, y, z, t)
        assert network.plucker_via_flows(frozenset({2})) == x * y * t / (1 + y * z)

    @given(le_diagrams())
    @settings(max_examples=25, deadline=None)
    def test_acyclic_collapse_makes_the_denominator_one(self, diagram):
        """Talaska's acyclic remark: no cycles, so conservative sum is 1."""
        filling, n = diagram
        network = bmm.PlanarNetwork.from_le_diagram(filling, n)
        assert network.is_acyclic
        flows = network.conservative_flows()
        assert list(flows) == [frozenset()]

    @given(square_weights, st.data())
    @settings(max_examples=25)
    def test_gauge_invariance_postnikov_eq_4_2(self, weights, data):
        """Postnikov eq. (4.2): gauge transformations preserve every M_ij.

        A single vertex scalar other than 1 is drawn as well, so the
        transformation is asserted to actually move the edge weights — the
        invariance laws are all satisfied by a no-op otherwise.
        """
        _, net = _oriented_square(*weights)
        scalars = {
            v: data.draw(positive_fractions, label=str(v))
            for v in sorted(net.internal_vertices, key=str)
        }
        rescaled = net.gauge_transform(scalars)
        assert rescaled.boundary_measurements() == net.boundary_measurements()
        vertex = data.draw(
            st.sampled_from(sorted(net.internal_vertices, key=str)), label="vertex"
        )
        factor = data.draw(positive_fractions.filter(lambda f: f != 1), label="factor")
        moved = net.gauge_transform({vertex: factor})
        assert moved.boundary_measurements() == net.boundary_measurements()
        assert moved.edges != net.edges

    def test_gauge_scales_out_edges_up_and_in_edges_down_postnikov_eq_4_2(self):
        """Postnikov eq. (4.2): x'_e = x_e t_u / t_v for e = (u, v).

        The direction, which the invariance laws cannot see — inverting
        every scalar is again a gauge transformation, so they hold either
        way. Scaling ``t_R = 2`` on the oriented square must multiply the
        two edges leaving ``R`` and divide the one entering it.
        """
        _, net = _oriented_square(2, 3, 5, 7)
        moved = net.gauge_transform({"R": Fraction(2)})
        touching = {
            (tail, head): weight
            for tail, head, weight in moved.edges
            if "R" in (tail, head)
        }
        assert touching == {
            (2, "R"): Fraction(1, 2),
            ("R", "T"): Fraction(6),
            ("R", "B"): Fraction(10),
        }

    @given(square_weights, st.data())
    @settings(max_examples=25)
    def test_face_weights_are_gauge_invariant_lam_lemma_4_6(self, weights, data):
        """Lam Lemma 4.6 (2015): face weights survive gauge equivalence."""
        _, net = _oriented_square(*weights)
        scalars = {
            v: data.draw(positive_fractions, label=str(v))
            for v in sorted(net.internal_vertices, key=str)
        }
        assert net.gauge_transform(scalars).face_weights() == net.face_weights()

    @given(square_weights)
    @settings(max_examples=25)
    def test_path_weight_is_the_face_product_postnikov_lemma_11_4(self, weights):
        """Postnikov Lemma 11.4 (2006): wt(P, y) = prod of P's edge weights.

        The product of the face weights to the right of a path equals the
        product of its edge weights — the statement that pins the
        exterior-clockwise orientation convention of ``face_weights``.
        Checked on the oriented square for the paths 2 -> R -> T -> 1
        (one face to the right, weight ``b``), 4 -> L -> B -> 3 (weight
        ``d``), and 2 -> R -> B -> 3 (four faces to the right, product
        ``c``).
        """
        a, b, c, d = weights
        _, net = _oriented_square(a, b, c, d)

        def face_weight_of(edge_labels: set[frozenset[int | str]]) -> Fraction:
            for face, weight in zip(net.faces(), net.face_weights(), strict=True):
                found = {
                    frozenset((net.edges[index][0], net.edges[index][1]))
                    for index, _ in face
                }
                if found == edge_labels:
                    return weight
            msg = f"no face bounded by {edge_labels}"
            raise AssertionError(msg)

        corner_12 = face_weight_of(
            {frozenset((1, "T")), frozenset(("R", "T")), frozenset((2, "R"))}
        )
        corner_34 = face_weight_of(
            {frozenset((3, "B")), frozenset(("B", "L")), frozenset((4, "L"))}
        )
        corner_23 = face_weight_of(
            {frozenset((2, "R")), frozenset(("R", "B")), frozenset((3, "B"))}
        )
        assert corner_12 == b
        assert corner_34 == d
        product = Fraction(1)
        for value in net.face_weights():
            product *= value
        assert product / corner_23 == c

    @given(square_weights)
    @settings(max_examples=25)
    def test_face_weights_multiply_to_one_postnikov_section_11(self, weights):
        """Postnikov section 11: the single relation prod_f y_f = 1."""
        _, net = _oriented_square(*weights)
        product = Fraction(1)
        for value in net.face_weights():
            product *= value
        assert product == 1

    @given(square_weights)
    @settings(max_examples=15, deadline=None)
    def test_matchings_biject_with_perfect_orientations_lam_prop_5_1(self, weights):
        """Lam Proposition 5.1 (2015): matchings <-> perfect orientations.

        Both directions: each matching orients the graph perfectly with
        source set ``I(Pi)``, distinct matchings give distinct
        orientations, and an independent brute-force enumeration of all
        perfect orientations finds exactly as many.
        """
        square = bmm.square_network(*weights)
        matchings = square.almost_perfect_matchings()
        oriented = [square.to_perfect_orientation(m) for m in matchings]
        for matching, net in zip(matchings, oriented, strict=True):
            assert net.is_perfectly_oriented
            assert net.source_set == square.boundary_subset(matching)
        edge_sets = {
            tuple(sorted((str(t), str(h)) for t, h, _ in net.edges)) for net in oriented
        }
        assert len(edge_sets) == len(matchings)
        perfect = 0
        for signs in itertools.product((False, True), repeat=len(square.edges)):
            out_count = dict.fromkeys(square.internal_vertices, 0)
            in_count = dict.fromkeys(square.internal_vertices, 0)
            for (u, v, _), flipped in zip(square.edges, signs, strict=True):
                tail, head = (v, u) if flipped else (u, v)
                if tail in out_count:
                    out_count[tail] += 1
                if head in in_count:
                    in_count[head] += 1
            perfect += all(
                out_count[v] == 1 if v in square.black_vertices else in_count[v] == 1
                for v in square.internal_vertices
            )
        assert perfect == len(matchings)

    @given(square_weights)
    @settings(max_examples=15, deadline=None)
    def test_flow_and_matching_pluckers_agree_lam_prop_5_3(self, weights):
        """Lam Proposition 5.3 (2015): wt(F) = wt(Pi) / wt(Pi_O)."""
        square = bmm.square_network(*weights)
        for matching in square.almost_perfect_matchings():
            net = square.to_perfect_orientation(matching)
            base = square.plucker(net.source_set)
            for subset in PAIRS:
                assert net.plucker(subset) == square.plucker(subset) / base

    @given(square_weights)
    @settings(max_examples=15, deadline=None)
    def test_orientation_invariance_postnikov_theorem_10_1(self, weights):
        """Postnikov Theorem 10.1 (2006): the point ignores the orientation."""
        square = bmm.square_network(*weights)
        matchings = square.almost_perfect_matchings()
        reference = {
            subset: square.to_perfect_orientation(matchings[0]).plucker(subset)
            for subset in PAIRS
        }
        for matching in matchings[1:]:
            other = {
                subset: square.to_perfect_orientation(matching).plucker(subset)
                for subset in PAIRS
            }
            assert _projectively_equal(reference, other, PAIRS)

    @given(square_weights)
    @settings(max_examples=20, deadline=None)
    def test_square_move_preserves_the_point_lam_prop_4_8(self, weights):
        """Lam Proposition 4.8 / Postnikov Lemma 12.2: moves fix the point."""
        square = bmm.square_network(*weights)
        moved = square.square_move(("T", "R", "B", "L"))
        assert _projectively_equal(
            _plucker_vector(square, PAIRS), _plucker_vector(moved, PAIRS), PAIRS
        )

    @given(square_weights)
    @settings(max_examples=20, deadline=None)
    def test_square_move_rescales_weights_and_recolors_lam_eq_19(self, weights):
        """Lam section 4.5, eq. (19): the moved weights are w / (ac + bd).

        Every square-side weight of the moved network is an original
        weight divided by ``ac + bd`` and the four corners are recolored;
        the placement of the primed weights (each landing on the opposite
        side) is the implementation's documented convention — the page
        records eq. (19)'s values but not the figure fixing their
        placement.
        """
        a, b, c, d = (Fraction(w) for w in weights)
        denominator = a * c + b * d
        square = bmm.square_network(a, b, c, d)
        moved = square.square_move(("T", "R", "B", "L"))
        new_weights = {frozenset((u, v)): w for u, v, w in moved.edges}
        assert new_weights[frozenset(("T", "R"))] == d / denominator
        assert new_weights[frozenset(("R", "B"))] == a / denominator
        assert new_weights[frozenset(("B", "L"))] == b / denominator
        assert new_weights[frozenset(("L", "T"))] == c / denominator
        for corner in ("T", "R", "B", "L"):
            assert (corner in moved.black_vertices) != (corner in square.black_vertices)

    @given(
        square_weights,
        st.integers(1, 3),
        positive_fractions,
        st.sampled_from(["x", "y"]),
    )
    @settings(max_examples=20, deadline=None)
    def test_bridges_are_chevalley_generators_lam_lemma_7_6(
        self, weights, i, a, variant
    ):
        """Lam Lemma 7.6 (2015): a bridge acts by x_i(a) or y_i(a)."""
        square, net = _oriented_square(*weights)
        bridged = square.add_bridge(i, a, variant=variant)
        acted = bmm.apply_chevalley(net.to_matrix(), i, a, variant=variant)
        expected = {subset: _pluck(acted, subset) for subset in PAIRS}
        assert _projectively_equal(_plucker_vector(bridged, PAIRS), expected, PAIRS)

    @given(square_weights, st.integers(1, 3))
    @settings(max_examples=20, deadline=None)
    def test_bridge_removal_lam_proposition_7_10(self, weights, i):
        """Lam Proposition 7.10 (2015): the reverse step shrinks the cell."""
        _, net = _oriented_square(*weights)
        columns = net.to_matrix()
        positroid = net.to_positroid()
        assert bmm.has_bridge(positroid, i) is True
        parameter = bmm.bridge_parameter(columns, i)
        assert parameter > 0
        removed = bmm.remove_bridge(columns, i)
        smaller = ps.Positroid.from_matrix(removed)
        assert smaller.bases < positroid.bases
        assert bmm.cell_dimension(smaller) == bmm.cell_dimension(positroid) - 1
        f = GrassmannNecklace.from_matroid(positroid).to_bounded_affine_permutation()
        swapped = list(f)
        swapped[i - 1], swapped[i] = swapped[i], swapped[i - 1]
        f_removed = GrassmannNecklace.from_matroid(
            smaller
        ).to_bounded_affine_permutation()
        assert f_removed == tuple(swapped)

    @pytest.mark.parametrize(
        ("bases", "window", "clauses"),
        [
            ([{3, 4}], (1, 2, 7, 8), (False, True, True)),
            ([{1}], (5, 2, 3, 4), (True, False, True)),
            ([{1, 2}], (5, 6, 3, 4), (True, True, False)),
        ],
    )
    def test_bridge_condition_needs_all_three_clauses_lam_section_7(
        self, bases, window, clauses
    ):
        """Lam section 7 (2015): i < i+1 <= f(i) < f(i+1) <= i+n.

        Each of the three inequalities is the sole reason some cell has no
        bridge at ``i = 1``, so none of them is redundant: the lollipop's
        cell fails the lower bound, the rank-one cell on ``{1}`` fails only
        ``f(i) < f(i+1)``, and the cell on ``{1, 2}`` fails only the upper
        bound. The satisfied direction is the top cell of ``Gr(2, 4)``,
        exercised by the bridge-removal test above.
        """
        i, n = 1, 4
        positroid = ps.Positroid.from_bases(range(1, n + 1), bases)
        f = GrassmannNecklace.from_matroid(positroid).to_bounded_affine_permutation()
        assert f == window
        assert (i + 1 <= f[i - 1], f[i - 1] < f[i], f[i] <= i + n) == clauses
        assert bmm.has_bridge(positroid, i) is False

    @pytest.mark.parametrize(
        ("bases", "window", "tight"),
        [
            ([{2}, {3}, {4}], (1, 3, 4, 6), "lower"),
            ([{2, 3}, {2, 4}, {3, 4}], (1, 4, 6, 7), "upper"),
        ],
    )
    def test_bridge_condition_bounds_are_inclusive_lam_section_7(
        self, bases, window, tight
    ):
        """Lam section 7 (2015): both bounds are ``<=``, not ``<``.

        Each cell here meets one of the two bounds with equality and still
        has a bridge, so neither may be tightened. (The middle ``<`` is
        not a third case: ``f`` is injective, so there ``<`` and ``<=``
        agree.)
        """
        i, n = 2, 4
        positroid = ps.Positroid.from_bases(range(1, n + 1), bases)
        f = GrassmannNecklace.from_matroid(positroid).to_bounded_affine_permutation()
        assert f == window
        assert (f[i - 1] == i + 1) is (tight == "lower")
        assert (f[i] == i + n) is (tight == "upper")
        assert bmm.has_bridge(positroid, i) is True

    @given(le_diagrams())
    @settings(max_examples=25, deadline=None)
    def test_necklace_and_affine_permutation_round_trip_lam_7_12_2(self, diagram):
        """Lam Theorem 7.12(2) (2015): positroids <-> necklaces <-> f_M."""
        filling, n = diagram
        positroid = bmm.PlanarNetwork.from_le_diagram(filling, n).to_positroid()
        necklace = GrassmannNecklace.from_matroid(positroid)
        window = necklace.to_bounded_affine_permutation()
        rebuilt = GrassmannNecklace.from_bounded_affine_permutation(
            positroid.elements, window
        )
        assert rebuilt == necklace


class TestTwist:
    """Muller-Speyer, in the zero-column-free form the page records."""

    @pytest.mark.parametrize("n", [4, 5])
    def test_gr2n_twist_is_a_cyclic_shift_muller_speyer_appendix_b(self, n):
        columns = {j: (Fraction(1), Fraction(j)) for j in range(1, n + 1)}
        twisted = bmm.right_twist(columns)

        def cyclic_minor(i: int, j: int) -> Fraction:
            a, b = sorted(((i - 1) % n + 1, (j - 1) % n + 1))
            minor = _pluck(columns, {a, b})
            return minor if (i - 1) % n + 1 == a else -minor

        for i, j in itertools.combinations(range(1, n + 1), 2):
            expected = cyclic_minor(i + 1, j + 1) / (
                cyclic_minor(i, i + 1) * cyclic_minor(j, j + 1)
            )
            assert _pluck(twisted, {i, j}) == expected

    def test_twists_are_mutually_inverse_muller_speyer_theorem_6_7(self):
        columns = {j: (Fraction(1), Fraction(j)) for j in range(1, 6)}
        assert bmm.right_twist(bmm.left_twist(columns)) == columns
        assert bmm.left_twist(bmm.right_twist(columns)) == columns

    def test_boundary_face_pluckers_invert_muller_speyer_corollary_6_8(self):
        n = 5
        columns = {j: (Fraction(1), Fraction(j)) for j in range(1, n + 1)}
        twisted = bmm.right_twist(columns)
        for a in range(1, n + 1):
            necklace_entry = {(a - 1) % n + 1, a % n + 1}
            assert _pluck(twisted, necklace_entry) == 1 / _pluck(
                columns, necklace_entry
            )

    def test_twist_of_the_square_point_round_trips(self):
        _, net = _oriented_square(2, 3, 5, 7)
        columns = net.to_matrix()
        assert bmm.left_twist(bmm.right_twist(columns)) == columns

    def test_twist_skips_columns_inside_the_running_span_muller_speyer_1_8(self):
        """Muller-Speyer section 1.8 (2017): the span condition on A_j.

        Column ``i`` of the right twist is orthogonal to ``A_j`` only
        ``when A_j is not in span{A_i, ..., A_{j-1}}``. Every moment-curve
        fixture has each consecutive column pair independent, so that
        clause is never reached there; this point has ``A_2 = 2 A_1``, so
        column 2 must be skipped when the constraints for column 1 are
        assembled. It is still totally nonnegative — its positroid is the
        page's five-basis degeneration — and the twists remain mutually
        inverse on it (Theorem 6.7).
        """
        columns = {
            1: (Fraction(1), Fraction(0)),
            2: (Fraction(2), Fraction(0)),
            3: (Fraction(1), Fraction(1)),
            4: (Fraction(0), Fraction(1)),
        }
        assert ps.Positroid.from_matrix(columns).bases == {
            frozenset(pair) for pair in [{1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}]
        }
        twisted = bmm.right_twist(columns)
        assert _pluck(twisted, {1, 2}) == 0
        assert bmm.left_twist(twisted) == columns
        assert bmm.right_twist(bmm.left_twist(columns)) == columns


# --------------------------------------------------------------------------- #
# Round-trip laws from the API contract
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @pytest.mark.parametrize(
        "network",
        [
            bmm.geometric_series_network(2, 3, 5, 7),
            bmm.acyclic_baseline_network(2, 3, 5),
            bmm.PlanarNetwork.from_le_diagram([[1, 0], [1]], 4),
        ],
    )
    def test_planar_network_dataframe_round_trip(self, network):
        assert bmm.PlanarNetwork.from_dataframe(network.to_dataframe()) == network

    @pytest.mark.parametrize(
        "network",
        [bmm.square_network(2, 3, 5, 7), bmm.lollipop_network()],
    )
    def test_bipartite_network_dataframe_round_trip(self, network):
        rebuilt = bmm.PlanarBipartiteNetwork.from_dataframe(network.to_dataframe())
        assert rebuilt == network

    def test_orientation_and_bipartite_form_invert_on_the_square(self):
        square = bmm.square_network(2, 3, 5, 7)
        matching = next(
            m
            for m in square.almost_perfect_matchings()
            if square.boundary_subset(m) == frozenset({2, 4})
        )
        back = square.to_perfect_orientation(matching).to_bipartite()
        assert back.black_vertices == square.black_vertices
        assert {(frozenset((u, v)), w) for u, v, w in back.edges} == {
            (frozenset((u, v)), w) for u, v, w in square.edges
        }
        assert all(back.plucker(subset) == square.plucker(subset) for subset in PAIRS)

    def test_le_graph_source_set_is_the_lexicographic_minimum(self):
        network = bmm.PlanarNetwork.from_le_diagram([[1, 1], [1, 1]], 4)
        assert network.source_set == {1, 2}
        assert bmm.cell_dimension(network.to_positroid()) == 4


# --------------------------------------------------------------------------- #
# Rejection tests: each named condition refuses a violating input
# --------------------------------------------------------------------------- #
class TestValidation:
    def test_nonpositive_weight_is_rejected(self):
        with pytest.raises(ValueError, match="strictly positive weights"):
            bmm.geometric_series_network(-1, 1, 1, 1)

    def test_boundary_vertex_with_two_edges_is_rejected(self):
        with pytest.raises(ValueError, match=r"Talaska Definition 2\.1"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1), "v": (1, 1)},
                [(1, "u", 1), ("v", 1, 1)],
            )

    def test_loops_are_rejected(self):
        with pytest.raises(ValueError, match="loop at"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1)},
                [(1, "u", 1), ("u", "u", 1)],
            )

    def test_repeated_edges_are_rejected(self):
        with pytest.raises(ValueError, match="repeated edge"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1), "v": (0, 2)},
                [(1, "u", 1), ("u", "v", 1), ("v", "u", 1)],
            )

    def test_vertex_on_an_edge_is_rejected(self):
        with pytest.raises(ValueError, match="lies on the edge"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 0)},
                [(1, 2, 1)],
            )

    def test_disconnected_pieces_are_rejected(self):
        with pytest.raises(ValueError, match="disconnected from the boundary"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1), "v": (0, 2)},
                [("u", "v", 1)],
            )

    def test_scrambled_boundary_order_is_rejected(self):
        _, net = _oriented_square(1, 1, 1, 1)
        with pytest.raises(ValueError, match="clockwise order"):
            bmm.PlanarNetwork.from_edges(
                (1, 3, 2, 4),
                dict(net.positions),
                list(net.edges),
            )

    def test_declared_isolated_source_with_an_edge_is_rejected(self):
        with pytest.raises(ValueError, match="isolated source"):
            bmm.PlanarNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0)},
                [(1, 2, 1)],
                isolated_sources=[1],
            )

    def test_monochromatic_edge_is_rejected(self):
        with pytest.raises(ValueError, match="bipartite"):
            bmm.PlanarBipartiteNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1), "v": (0, 2)},
                [(1, "u", 1), ("u", "v", 1), ("v", 2, 1)],
                black_vertices=["u", "v"],
            )

    def test_bipartite_boundary_degree_must_be_one(self):
        with pytest.raises(ValueError, match="degree-one boundary"):
            bmm.PlanarBipartiteNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1)},
                [(1, "u", 1)],
                black_vertices=["u"],
            )

    def test_isolated_interior_vertex_is_rejected(self):
        with pytest.raises(ValueError, match="standing assumption"):
            bmm.PlanarBipartiteNetwork.from_edges(
                (1, 2),
                {1: (-2, 0), 2: (2, 0), "u": (0, 1), "w": (0, -1)},
                [(1, "u", 1), ("u", 2, 1)],
                black_vertices=["u"],
            )

    def test_le_condition_violation_is_rejected(self):
        with pytest.raises(ValueError, match="Le-condition"):
            bmm.PlanarNetwork.from_le_diagram([[0, 1], [1, 0]], 4)

    def test_le_shape_must_be_weakly_decreasing(self):
        with pytest.raises(ValueError, match="weakly decreasing"):
            bmm.PlanarNetwork.from_le_diagram([[1], [1, 1]], 4)

    def test_le_shape_must_fit_the_rectangle(self):
        with pytest.raises(ValueError, match="rectangle"):
            bmm.PlanarNetwork.from_le_diagram([[1, 1, 1]], 3)

    def test_le_filling_must_be_zero_one(self):
        with pytest.raises(ValueError, match="0s and 1s"):
            bmm.PlanarNetwork.from_le_diagram([[2]], 2)

    @pytest.mark.parametrize(
        ("filling", "n", "weights"),
        [
            # Postnikov section 6 states T(i, j) > 0 *iff* the box holds a
            # 1, so both directions reject: a 1-box left at zero, and a
            # positive weight on a 0-box or off the shape entirely.
            ([[1]], 2, {(1, 1): Fraction(0)}),
            ([[1, 0], [1]], 4, {(1, 2): Fraction(3)}),
            ([[1, 0], [1]], 4, {(2, 2): Fraction(3)}),
        ],
    )
    def test_gamma_tableau_must_be_positive_exactly_on_ones(self, filling, n, weights):
        with pytest.raises(ValueError, match="Gamma-tableau"):
            bmm.PlanarNetwork.from_le_diagram(filling, n, weights)

    def test_flow_formula_requires_perfect_orientation(self):
        crossing = bmm.PlanarNetwork.from_edges(
            (1, 2, 3, 4),
            {1: (0, 2), 2: (2, 0), 3: (0, -2), 4: (-2, 0), "X": (0, 0)},
            [(1, "X", 2), (4, "X", 3), ("X", 2, 5), ("X", 3, 7)],
        )
        assert crossing.is_acyclic
        with pytest.raises(ValueError, match="perfectly oriented"):
            crossing.plucker_via_flows(frozenset({1, 4}))

    def test_cyclic_non_perfect_measurements_are_refused(self):
        network = bmm.PlanarNetwork.from_edges(
            (1, 2, 3),
            {
                1: (-2, 0),
                2: (2, 0),
                3: (0, -2),
                "u": (-1, 0),
                "v": (1, 0),
                "w": (0, 1),
            },
            [
                (1, "u", 2),
                ("u", "v", 3),
                ("v", "w", 5),
                ("w", "u", 7),
                ("v", 2, 11),
                ("u", 3, 13),
            ],
        )
        assert not network.is_acyclic
        assert not network.is_perfectly_oriented
        with pytest.raises(ValueError, match="signed series"):
            network.boundary_measurement(1, 2)

    def test_measurement_indices_must_be_source_and_sink(self):
        network = bmm.geometric_series_network()
        with pytest.raises(ValueError, match="source and a sink"):
            network.boundary_measurement(2, 1)

    def test_plucker_subsets_must_have_size_k(self):
        network = bmm.geometric_series_network()
        with pytest.raises(ValueError, match="1-subsets"):
            network.plucker({1, 2})

    def test_matching_enumeration_guard_trips(self):
        chain = 21
        positions = {j: (j, 0) for j in range(chain)}
        positions[0] = (0, 0)
        edges = [(j, j + 1, 1) for j in range(chain - 1)]
        network = bmm.PlanarBipartiteNetwork.from_edges(
            (0, chain - 1),
            positions,
            edges,
            black_vertices=[j for j in range(chain) if j % 2],
        )
        with pytest.raises(ValueError, match="enumeration guard"):
            network.almost_perfect_matchings()

    def test_gauge_scalars_must_be_internal_and_positive(self):
        network = bmm.geometric_series_network()
        with pytest.raises(ValueError, match="internal vertices only"):
            network.gauge_transform({1: Fraction(2)})
        with pytest.raises(ValueError, match="must be positive"):
            network.gauge_transform({"u": Fraction(0)})

    def test_chevalley_position_is_range_checked(self):
        columns = {j: (Fraction(1), Fraction(j)) for j in range(1, 5)}
        with pytest.raises(ValueError, match="bridge position"):
            bmm.apply_chevalley(columns, 4, Fraction(1))

    def test_bridge_parameter_requires_a_bridge(self):
        lollipop = bmm.lollipop_network()
        columns = {
            1: (Fraction(0), Fraction(0)),
            2: (Fraction(0), Fraction(0)),
            3: (Fraction(1), Fraction(0)),
            4: (Fraction(0), Fraction(1)),
        }
        assert lollipop.to_positroid().bases == {frozenset({3, 4})}
        with pytest.raises(ValueError, match="no bridge"):
            bmm.bridge_parameter(columns, 1)

    def test_twist_requires_no_zero_columns(self):
        columns = {
            1: (Fraction(0), Fraction(0)),
            2: (Fraction(1), Fraction(0)),
            3: (Fraction(0), Fraction(1)),
        }
        with pytest.raises(ValueError, match="zero"):
            bmm.right_twist(columns)

    def test_square_move_requires_alternating_colors(self):
        square = bmm.square_network(2, 3, 5, 7)
        with pytest.raises(ValueError, match="alternate colors"):
            square.square_move(("T", "B", "R", "L"))

    def test_square_move_requires_gauge_fixed_legs(self):
        square = bmm.square_network(2, 3, 5, 7)
        relegged = bmm.PlanarBipartiteNetwork.from_edges(
            square.boundary,
            dict(square.positions),
            [(u, v, Fraction(2) if 1 in (u, v) else w) for u, v, w in square.edges],
            black_vertices=square.black_vertices,
        )
        with pytest.raises(ValueError, match="weight 1"):
            relegged.square_move(("T", "R", "B", "L"))

    def test_reducedness_rejects_an_edge_used_twice_at_a_non_leaf_hinge(self):
        """Lam section 7.1 (2015): the trip-repeats-an-edge condition.

        A reduced graph has no trip using an edge twice *except* at a
        boundary leaf. This fixture isolates that clause: the other two
        are checked here to be silent — no trip is closed, and no two
        trips share two edges in the same order — so the rejection can
        only come from the edge repeats, which occur at interior hinges
        and non-consecutively, away from the leaf bounce the clause
        excepts (the lollipop, still reduced, is that excepted case).
        Postnikov Theorem 12.7 corroborates from outside the trip
        machinery: a reduced graph would have ``dim = |F(G)| - 1``.
        """
        network = _stacked_squares()
        trips = network.trips()
        assert len(trips) == len(network.boundary)  # the closed-trip clause is silent
        assert _trips_share_two_edges_in_order(network) is False  # and so is the third
        sequence = [edge for _, edge in trips[0]]
        repeated = sorted(e for e in set(sequence) if sequence.count(e) > 1)
        assert repeated
        for edge in repeated:
            places = [i for i, other in enumerate(sequence) if other == edge]
            assert places[1] > places[0] + 1  # not an immediate bounce
            tail = trips[0][places[0]][0]
            (hinge,) = (v for v in network.edges[edge][:2] if v != tail)
            degree = sum(1 for u, v, _ in network.edges if hinge in (u, v))
            assert degree > 1  # not a boundary leaf
        assert network.is_reduced() is False
        directed = network.to_perfect_orientation()
        assert bmm.cell_dimension(network.to_positroid()) != directed.face_count - 1
        assert bmm.lollipop_network().is_reduced() is True

    def test_reducedness_scope_excludes_interior_leaves(self):
        network = bmm.PlanarBipartiteNetwork.from_edges(
            (1, 2),
            {1: (-2, 0), 2: (2, 0), "u": (0, 0), "leaf": (0, 1)},
            [(1, "u", 1), ("u", 2, 1), ("u", "leaf", 1)],
            black_vertices=["u"],
        )
        with pytest.raises(ValueError, match="leafless"):
            network.is_reduced()

    def test_dataframe_columns_are_required(self):
        with pytest.raises(ValueError, match="missing required columns"):
            bmm.PlanarNetwork.from_dataframe(pd.DataFrame({"kind": []}))
        with pytest.raises(ValueError, match="missing required columns"):
            bmm.PlanarBipartiteNetwork.from_dataframe(pd.DataFrame({"kind": []}))


# --------------------------------------------------------------------------- #
# Visualization smoke tests
# --------------------------------------------------------------------------- #
class TestPlotting:
    def test_planar_network_plot_draws_on_given_axes(self):
        _, ax = plt.subplots()
        try:
            returned = bmm.geometric_series_network().plot_network(ax)
            assert returned is ax
        finally:
            plt.close("all")

    def test_bipartite_plot_creates_axes_when_none_given(self):
        try:
            ax = bmm.square_network().plot_network()
            assert ax is not None
        finally:
            plt.close("all")
