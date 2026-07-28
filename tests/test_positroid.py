"""Tests for research.positroid, transcribed from the Logseq Positroid page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Oh's theorem, the
Ardila-Rincon-Williams closure results, the non-crossing decomposition, the
polytope characterization, and the enumeration bijections); round-trip laws
come from the API contract; and every named condition has a rejection test
naming it. The Hypothesis generator walks decorated permutations, which by
the page's bijection reach every positroid on ``[n]``.
"""

import functools
import itertools
from fractions import Fraction

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from experiments import io
from research import matroid as mt
from research import positroid as ps

matplotlib.use("Agg")


# --------------------------------------------------------------------------- #
# Strategies: positroids through the decorated-permutation bijection
# --------------------------------------------------------------------------- #
def _draw_decorated_permutation(draw: st.DrawFn, n: int) -> ps.DecoratedPermutation:
    labels = tuple(range(1, n + 1))
    targets = tuple(draw(st.permutations(labels))) if n else ()
    fixed = [i for i in labels if targets[i - 1] == i]
    clockwise = (
        frozenset(draw(st.sets(st.sampled_from(fixed)))) if fixed else frozenset()
    )
    return ps.DecoratedPermutation(targets, clockwise)


def _positroid_on(draw: st.DrawFn, elements: tuple[int, ...]) -> ps.Positroid[int]:
    decorated = _draw_decorated_permutation(draw, len(elements))
    return ps.Positroid.from_decorated_permutation(elements, decorated)


@st.composite
def positroids(draw: st.DrawFn, max_n: int = 6) -> ps.Positroid[int]:
    n = draw(st.integers(0, max_n))
    return _positroid_on(draw, tuple(range(1, n + 1)))


@st.composite
def positroid_pairs(
    draw: st.DrawFn,
) -> tuple[ps.Positroid[int], ps.Positroid[int]]:
    n1 = draw(st.integers(0, 4))
    n2 = draw(st.integers(0, 4))
    first = _positroid_on(draw, tuple(range(1, n1 + 1)))
    second = _positroid_on(draw, tuple(range(n1 + 1, n1 + n2 + 1)))
    return first, second


def _position(p: ps.Positroid[int], element: int) -> int:
    return p.elements.index(element) + 1


def _interval_labels(p: ps.Positroid[int], i: int, j: int) -> list[int]:
    n = len(p.elements)
    length = (j - i) % n + 1
    return [p.elements[(i - 1 + t) % n] for t in range(length)]


def _blocks_cross(first: set[int], second: set[int]) -> bool:
    def pattern(x: set[int], y: set[int]) -> bool:
        return any(
            a < b < c < d
            for a, c in itertools.combinations(sorted(x), 2)
            for b, d in itertools.combinations(sorted(y), 2)
        )

    return pattern(first, second) or pattern(second, first)


_CROSSING_BASES = [[1, 2], [1, 4], [2, 3], [3, 4]]
"""``U_{1,2} + U_{1,2}`` on the crossing blocks ``{1,3}`` and ``{2,4}``."""

_NONCROSSING_BASES = [[1, 3], [1, 4], [2, 3], [2, 4]]
"""The same matroid on the non-crossing blocks ``{1,2}`` and ``{3,4}``."""


# --------------------------------------------------------------------------- #
# Canonical examples: each asserts what the page says the example certifies
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_uniform_positroid_is_the_matroid_of_the_top_cell(self):
        p = ps.uniform_positroid(2, 4)
        assert p.rank() == 2
        assert p.bases == {
            frozenset(pair) for pair in itertools.combinations(range(4), 2)
        }
        assert ps.is_positroid(p) is True

    def test_totally_positive_matrix_realizes_the_uniform_positroid(self):
        # A Vandermonde matrix with increasing positive nodes is totally
        # positive, so every maximal minor is positive: the top cell.
        d, n = 2, 4
        columns = {j: tuple(j**r for r in range(d)) for j in range(1, n + 1)}
        p = ps.Positroid.from_matrix(columns)
        assert p.bases == {
            frozenset(pair) for pair in itertools.combinations(range(1, n + 1), 2)
        }

    @pytest.mark.parametrize("position", [1, 2, 3, 4])
    def test_shifted_schubert_matroids_are_positroids(self, position):
        p = ps.shifted_schubert_positroid(range(1, 5), [1, 2], position)
        assert isinstance(p, ps.Positroid)
        assert ps.is_positroid(p) is True

    def test_crossing_direct_sum_is_not_a_positroid(self):
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_bases((1, 2, 3, 4), _CROSSING_BASES)

    def test_noncrossing_direct_sum_is_a_positroid(self):
        p = ps.Positroid.from_bases((1, 2, 3, 4), _NONCROSSING_BASES)
        assert p.connected_components() == {
            frozenset({1, 2}),
            frozenset({3, 4}),
        }

    def test_reordering_to_crossing_blocks_destroys_the_positroid(self):
        p = ps.Positroid.from_bases((1, 2, 3, 4), _NONCROSSING_BASES)
        with pytest.raises(ValueError, match="Oh's theorem"):
            p.with_cyclic_order((1, 3, 2, 4))

    def test_empty_positroid_round_trips_through_every_view(self):
        p = ps.uniform_positroid(0, 0)
        assert p.grassmann_necklace == ()
        assert p.to_decorated_permutation() == ps.DecoratedPermutation(())
        assert ps.Positroid.from_dataframe(p.to_dataframe()) == p


# --------------------------------------------------------------------------- #
# Constructor validation: each named condition rejected with its name
# --------------------------------------------------------------------------- #
class TestConstructorValidation:
    def test_negative_maximal_minor_is_rejected_by_name(self):
        with pytest.raises(ValueError, match="maximal minor"):
            ps.Positroid.from_matrix({1: (1, 0), 2: (0, 1), 3: (1, -1)})

    def test_rank_deficient_matrix_is_rejected_by_name(self):
        with pytest.raises(ValueError, match="full rank"):
            ps.Positroid.from_matrix({1: (1, 0), 2: (2, 0)})

    def test_inconsistent_column_dimensions_are_rejected(self):
        with pytest.raises(ValueError, match="dimension"):
            ps.Positroid.from_matrix({1: (1, 0), 2: (1,)})

    def test_fractional_matrix_entries_are_exact(self):
        p = ps.Positroid.from_matrix(
            {1: (Fraction(1, 3),), 2: (Fraction(2, 7),), 3: (0,)}
        )
        assert p.loops == {3}
        assert p.rank() == 1

    def test_necklace_with_wrong_length_is_rejected(self):
        with pytest.raises(ValueError, match="one entry per"):
            ps.Positroid.from_grassmann_necklace((1, 2), [{1}])

    def test_necklace_violating_the_membership_transition_is_rejected(self):
        # 1 lies in I_1, so I_2 must be (I_1 - {1}) + {j}; {1, 3} is not.
        with pytest.raises(ValueError, match="Postnikov section 16"):
            ps.Positroid.from_grassmann_necklace((1, 2, 3), [{1, 2}, {1, 3}, {1, 3}])

    def test_necklace_violating_the_fixed_transition_is_rejected(self):
        # 1 is not in I_1, so I_2 must equal I_1.
        with pytest.raises(ValueError, match="must equal it"):
            ps.Positroid.from_grassmann_necklace((1, 2), [{2}, {1}])

    def test_non_bijective_decorated_permutation_is_rejected(self):
        with pytest.raises(ValueError, match="bijection"):
            ps.DecoratedPermutation((1, 1))

    def test_decorating_a_non_fixed_point_is_rejected(self):
        with pytest.raises(ValueError, match="not fixed"):
            ps.DecoratedPermutation((2, 1), frozenset({1}))

    def test_decorated_permutation_size_mismatch_is_rejected(self):
        decorated = ps.DecoratedPermutation((2, 1))
        with pytest.raises(ValueError, match="ground set has"):
            ps.Positroid.from_decorated_permutation((1, 2, 3), decorated)

    def test_from_vectors_gate_rejects_a_crossing_realization(self):
        vectors = {1: (1, 0), 2: (0, 1), 3: (1, 0), 4: (0, 1)}
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_vectors(vectors)

    def test_from_independent_sets_gate_rejects_the_crossing_sum(self):
        independent = [[], [1], [2], [3], [4], *_CROSSING_BASES]
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_independent_sets((1, 2, 3, 4), independent)

    def test_from_circuits_gate_rejects_the_crossing_sum(self):
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_circuits((1, 2, 3, 4), [[1, 3], [2, 4]])

    def test_from_graph_edges_gate_rejects_interleaved_parallel_edges(self):
        edges = {1: ("u", "v"), 2: ("x", "y"), 3: ("u", "v"), 4: ("x", "y")}
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_graph_edges(edges)

    def test_from_matroid_rejects_and_accepts_by_order(self):
        crossing = mt.Matroid.from_bases((1, 2, 3, 4), _CROSSING_BASES)
        with pytest.raises(ValueError, match="Oh's theorem"):
            ps.Positroid.from_matroid(crossing)
        accepted = ps.Positroid.from_matroid(
            mt.Matroid.from_bases((1, 2, 3, 4), _NONCROSSING_BASES)
        )
        assert isinstance(accepted, ps.Positroid)


# --------------------------------------------------------------------------- #
# Grassmann necklaces and decorated permutations
# --------------------------------------------------------------------------- #
class TestCryptomorphicIndexings:
    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_oh_theorem_necklace_round_trip(self, p):
        # Oh 2011: M(I(M)) = M for every positroid M.
        rebuilt = ps.Positroid.from_grassmann_necklace(p.elements, p.grassmann_necklace)
        assert rebuilt == p

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_necklace_satisfies_the_postnikov_conditions(self, p):
        # Postnikov section 16: the two transition conditions, indices mod n.
        necklace = p.grassmann_necklace
        n = len(p.elements)
        for i in range(n):
            element = p.elements[i]
            current, following = necklace[i], necklace[(i + 1) % n]
            if element in current:
                kept = current - {element}
                assert kept <= following
                assert len(following - kept) == 1
            else:
                assert following == current

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_necklace_entries_are_bases(self, p):
        # Each entry is the Gale-minimal basis of a shifted order, hence a
        # basis (Positroid page, "From a matroid").
        assert all(entry in p.bases for entry in p.grassmann_necklace)

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_positroid_is_intersection_of_shifted_schubert_matroids(self, p):
        # Oh 2011: each condition B >=_j I_j cuts out one cyclically shifted
        # Schubert matroid, and the positroid is their intersection.
        assume(p.elements)
        necklace = p.grassmann_necklace
        blocks = [
            ps.shifted_schubert_positroid(p.elements, entry, j + 1).bases
            for j, entry in enumerate(necklace)
        ]
        assert functools.reduce(frozenset.intersection, blocks) == p.bases

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_decorated_permutation_round_trip(self, p):
        decorated = p.to_decorated_permutation()
        assert ps.Positroid.from_decorated_permutation(p.elements, decorated) == p

    @settings(max_examples=50, deadline=None)
    @given(st.data())
    def test_decorated_permutation_bijection_is_injective(self, data):
        n = data.draw(st.integers(0, 5))
        decorated = _draw_decorated_permutation(data.draw, n)
        p = ps.Positroid.from_decorated_permutation(tuple(range(1, n + 1)), decorated)
        assert p.to_decorated_permutation() == decorated

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_weak_excedance_count_is_the_rank(self, p):
        # Positroid page: rank-d positroids correspond to decorated
        # permutations with exactly d weak excedances.
        assert p.to_decorated_permutation().weak_excedance_count == p.rank()

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_loops_and_coloops_are_the_colored_fixed_points(self, p):
        # Positroid page: loops and coloops are exactly the fixed points of
        # the two colors; this library maps coloops to clockwise.
        decorated = p.to_decorated_permutation()
        clockwise = {p.elements[i - 1] for i in decorated.clockwise_fixed}
        counterclockwise = {p.elements[i - 1] for i in decorated.counterclockwise_fixed}
        assert p.coloops == clockwise
        assert p.loops == counterclockwise


# --------------------------------------------------------------------------- #
# Closure properties (Ardila-Rincon-Williams)
# --------------------------------------------------------------------------- #
class TestClosureProperties:
    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_closure_under_duality(self, p):
        # ARW Prop. 3.5, plus the duality involution from the Matroid page.
        dual = p.dual()
        assert isinstance(dual, ps.Positroid)
        assert ps.is_positroid(dual) is True
        assert dual.dual() == p

    @settings(max_examples=50, deadline=None)
    @given(positroids(), st.data())
    def test_closure_under_restriction_and_contraction(self, p, data):
        # ARW Prop. 3.5, with the inherited cyclic order.
        subset = (
            data.draw(st.sets(st.sampled_from(p.elements))) if p.elements else set()
        )
        restricted = p.restrict(subset)
        contracted = p.contract(subset)
        assert isinstance(restricted, ps.Positroid)
        assert isinstance(contracted, ps.Positroid)
        assert ps.is_positroid(restricted) is True
        assert ps.is_positroid(contracted) is True

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_closure_under_cyclic_shift(self, p):
        # ARW Lemma 3.3: every rotation of the cyclic order stays a positroid.
        for steps in range(len(p.elements)):
            shifted = p.cyclic_shift(steps)
            assert ps.is_positroid(shifted) is True
            assert shifted == p

    @settings(max_examples=50, deadline=None)
    @given(positroid_pairs())
    def test_closure_under_direct_sum(self, pair):
        # ARW Prop. 3.4: concatenation places the summands on cyclic
        # intervals, so the direct sum of positroids is a positroid.
        first, second = pair
        combined = first.direct_sum(second)
        assert isinstance(combined, ps.Positroid)
        assert ps.is_positroid(combined) is True

    def test_direct_sum_with_a_plain_matroid_stays_plain(self):
        p = ps.uniform_positroid(1, 2)
        plain = mt.Matroid.from_bases(("a", "b"), [["a"], ["b"]])
        assert not isinstance(p.direct_sum(plain), ps.Positroid)

    def test_minor_delete_and_simplification_return_positroids(self):
        p = ps.uniform_positroid(2, 4)
        assert isinstance(p.minor(deletions=[0], contractions=[1]), ps.Positroid)
        assert isinstance(p.delete([2]), ps.Positroid)
        assert isinstance(p.simplification(), ps.Positroid)

    def test_to_matroid_forgets_the_positroid_reading(self):
        p = ps.uniform_positroid(2, 4)
        plain = p.to_matroid()
        assert type(plain) is mt.Matroid
        assert plain == p


# --------------------------------------------------------------------------- #
# Non-crossing decomposition (ARW Thm. 7.6)
# --------------------------------------------------------------------------- #
class TestNonCrossingDecomposition:
    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_components_form_a_noncrossing_partition(self, p):
        blocks = [
            {_position(p, e) for e in block} for block in p.connected_components()
        ]
        assert all(
            not _blocks_cross(a, b) for a, b in itertools.combinations(blocks, 2)
        )

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_components_restrict_to_connected_positroids(self, p):
        for block in p.connected_components():
            restricted = p.restrict(block)
            assert restricted.is_connected is True
            assert ps.is_positroid(restricted) is True

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_direct_sum_of_component_restrictions_rebuilds_the_positroid(self, p):
        assume(p.elements)
        ordered = sorted(p.connected_components(), key=min)
        pieces = [p.restrict(block) for block in ordered]
        assert functools.reduce(lambda a, b: a.direct_sum(b), pieces) == p


# --------------------------------------------------------------------------- #
# Polytope characterization (ARW Prop. 5.6, vertex form, via Oh)
# --------------------------------------------------------------------------- #
class TestCyclicIntervalCharacterization:
    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_bases_are_cut_out_by_cyclic_interval_ranks(self, p):
        # ARW Prop. 5.6 restricted to 0/1 vertices with the tight bounds
        # a_ij = r([i, j]): the d-subsets satisfying every cyclic-interval
        # inequality are exactly the bases.
        n = len(p.elements)
        d = p.rank()
        bounds = p.cyclic_rank_bounds()
        for combo in itertools.combinations(p.elements, d):
            chosen = set(combo)
            satisfies = all(
                len(chosen & set(_interval_labels(p, i, j))) <= bounds[(i, j)]
                for i in range(1, n + 1)
                for j in range(1, n + 1)
            )
            assert satisfies == (frozenset(combo) in p.bases)

    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
    def test_recognition_matches_the_interval_criterion_on_all_matroids(self, n):
        # Both directions of the characterization, over every matroid on n
        # elements: the Oh membership test agrees with the cyclic-interval
        # rank criterion.
        for m in mt.enumerate_matroids(n):
            d = m.rank()
            intervals = [
                [m.elements[(i + t) % n] for t in range((j - i) % n + 1)]
                for i in range(n)
                for j in range(n)
            ]
            criterion = all(
                (frozenset(combo) in m.bases)
                == all(len(set(combo) & set(iv)) <= m.rank(iv) for iv in intervals)
                for combo in itertools.combinations(m.elements, d)
            )
            assert ps.is_positroid(m) == criterion


# --------------------------------------------------------------------------- #
# Enumeration (OEIS A000522 and A075834)
# --------------------------------------------------------------------------- #
class TestEnumeration:
    @pytest.mark.parametrize(
        ("n", "expected"), [(0, 1), (1, 2), (2, 5), (3, 16), (4, 65), (5, 326)]
    )
    def test_counts_match_oeis_a000522(self, n, expected):
        found = ps.enumerate_positroids(n)
        assert len(found) == expected
        assert len(set(found)) == expected

    def test_connected_counts_match_ardila_rincon_williams(self):
        # ARW Cor. 7.11 states the SIF bijection for n >= 2 only; the paper
        # reports pc(1) = 2 (both single-element positroids are connected)
        # and notes its sequence equals OEIS A075834 "except for the first
        # term". The Positroid page omits that caveat — reported upstream.
        counts = [
            sum(1 for p in ps.enumerate_positroids(n) if p.is_connected)
            for n in range(1, 5)
        ]
        assert counts == [2, 1, 2, 7]

    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_connected_iff_decorated_permutation_is_sif(self, p):
        # ARW Cor. 7.11 (via Cor. 7.9): a positroid is connected iff its
        # decorated permutation is stabilized-interval-free.
        decorated = p.to_decorated_permutation()
        assert p.is_connected == decorated.is_stabilized_interval_free

    def test_sif_permutations_of_three_elements_are_the_two_cycles(self):
        sif = [
            targets
            for targets in itertools.permutations((1, 2, 3))
            if ps.DecoratedPermutation(targets).is_stabilized_interval_free
        ]
        assert sif == [(2, 3, 1), (3, 1, 2)]


# --------------------------------------------------------------------------- #
# DataFrame serialization and experiment IO
# --------------------------------------------------------------------------- #
class TestDataFrames:
    @settings(max_examples=50, deadline=None)
    @given(positroids())
    def test_dataframe_round_trip(self, p):
        decoded = ps.Positroid.from_dataframe(p.to_dataframe())
        assert isinstance(decoded, ps.Positroid)
        assert decoded == p
        assert decoded.elements == p.elements

    def test_experiment_io_round_trip_preserves_the_cyclic_order(self, tmp_path):
        p = ps.Positroid.from_bases((1, 2, 3, 4), _NONCROSSING_BASES)
        path = io.write_result(p.to_dataframe(), tmp_path / "positroid.json")
        decoded = ps.Positroid.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == p
        assert decoded.elements == p.elements


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestVisualization:
    def test_plot_decorated_permutation_draws_on_provided_axes(self):
        _, ax = plt.subplots()
        try:
            returned = ps.uniform_positroid(2, 4).plot_decorated_permutation(ax)
            assert returned is ax
            # Four element labels plus four arrow annotations: pi = (3,4,1,2)
            # has no fixed points, so every position carries one arrow.
            assert len(ax.texts) == 8
        finally:
            plt.close("all")

    def test_plot_decorated_permutation_creates_axes_when_omitted(self):
        try:
            ax = ps.uniform_positroid(1, 3).plot_decorated_permutation()
            assert ax.get_title() == "Decorated permutation"
        finally:
            plt.close("all")

    def test_plot_connected_components_draws_noncrossing_loops(self):
        _, ax = plt.subplots()
        try:
            p = ps.Positroid.from_bases((1, 2, 3, 4), _NONCROSSING_BASES)
            returned = p.plot_connected_components(ax)
            assert returned is ax
            assert len(ax.texts) == 4
            assert len(ax.lines) == 2
        finally:
            plt.close("all")

    def test_plot_connected_components_creates_axes_when_omitted(self):
        try:
            ax = ps.uniform_positroid(2, 4).plot_connected_components()
            assert ax.get_title() == "Connected components (non-crossing partition)"
        finally:
            plt.close("all")
