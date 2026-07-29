"""Tests for research.grassmann_necklace, from the Logseq Grassmann Necklace page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Postnikov Lemmas 16.2,
16.3 and 17.6, Corollary 17.7, Oh's theorem and his Theorems 19 and Lemma
21, the Ardila-Rincon-Williams envelope, Lam Theorem 6.2, the
Knutson-Lam-Speyer juggling condition, and the A046802 enumeration);
round-trip laws come from the API contract; and every named condition has a
rejection test naming it. The Hypothesis generator walks decorated
permutations, which by Lemma 16.2 reach every necklace on ``[n]``.
"""

import itertools

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from experiments import io
from research import grassmann_necklace as gn
from research import matroid as mt
from research import positroid as ps

matplotlib.use("Agg")


# --------------------------------------------------------------------------- #
# Strategies: necklaces through the decorated-permutation bijection
# --------------------------------------------------------------------------- #
def _draw_decorated_permutation(draw: st.DrawFn, n: int) -> ps.DecoratedPermutation:
    labels = tuple(range(1, n + 1))
    targets = tuple(draw(st.permutations(labels))) if n else ()
    fixed = [i for i in labels if targets[i - 1] == i]
    clockwise = (
        frozenset(draw(st.sets(st.sampled_from(fixed)))) if fixed else frozenset()
    )
    return ps.DecoratedPermutation(targets, clockwise)


@st.composite
def necklaces(draw: st.DrawFn, max_n: int = 6) -> gn.GrassmannNecklace[int]:
    n = draw(st.integers(0, max_n))
    decorated = _draw_decorated_permutation(draw, n)
    return gn.GrassmannNecklace.from_decorated_permutation(
        tuple(range(1, n + 1)), decorated
    )


def _gale_leq(
    candidate: frozenset[int], bound: frozenset[int], start: int, n: int
) -> bool:
    """Return ``candidate <=_start bound`` in the shifted Gale order."""
    shifted_candidate = sorted((x - start) % n for x in candidate)
    shifted_bound = sorted((x - start) % n for x in bound)
    return all(c <= b for c, b in zip(shifted_candidate, shifted_bound, strict=True))


def _all_necklaces(n: int) -> list[gn.GrassmannNecklace[int]]:
    return [
        necklace
        for rank in range(n + 1)
        for necklace in gn.enumerate_necklaces(rank, n)
    ]


# --------------------------------------------------------------------------- #
# Canonical examples: each asserts what the page says the example certifies
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_uniform_necklace_lists_the_cyclic_intervals(self):
        necklace = gn.uniform_necklace(2, 4)
        assert necklace.entries == (
            frozenset({1, 2}),
            frozenset({2, 3}),
            frozenset({3, 4}),
            frozenset({4, 1}),
        )
        assert necklace.necklace_type == (2, 4)

    def test_uniform_necklace_yields_the_full_uniform_positroid(self):
        # OPS section 4: M(I) is all of C([n], k) for the interval necklace.
        necklace = gn.uniform_necklace(2, 4)
        assert necklace.to_positroid().bases == {
            frozenset(pair) for pair in itertools.combinations(range(1, 5), 2)
        }

    @pytest.mark.parametrize(("rank", "n"), [(1, 3), (2, 4)])
    def test_uniform_necklace_is_the_unique_circular_bruhat_top(self, rank, n):
        # Postnikov Lemma 17.6 via the Corollary 17.7 comparison test.
        top = gn.uniform_necklace(rank, n)
        population = gn.enumerate_necklaces(rank, n)
        assert all(necklace.circular_bruhat_leq(top) for necklace in population)
        tops = [
            candidate
            for candidate in population
            if all(other.circular_bruhat_leq(candidate) for other in population)
        ]
        assert tops == [top]

    @pytest.mark.parametrize(("rank", "n"), [(1, 3), (2, 4)])
    def test_constant_necklaces_are_the_circular_bruhat_minimal_elements(self, rank, n):
        # Postnikov Lemma 17.6: the C(n, k) minimal elements are constant.
        population = gn.enumerate_necklaces(rank, n)
        minimal = {
            candidate
            for candidate in population
            if not any(
                other != candidate and other.circular_bruhat_leq(candidate)
                for other in population
            )
        }
        expected = {
            gn.constant_necklace(n, combo)
            for combo in itertools.combinations(range(1, n + 1), rank)
        }
        assert minimal == expected

    def test_constant_necklace_collapses_to_the_single_basis_positroid(self):
        # Page-verified fixture: B(({1,3}, ..., {1,3})) = {{1,3}} on [4].
        necklace = gn.constant_necklace(4, {1, 3})
        assert necklace.is_constant is True
        assert necklace.to_positroid().bases == {frozenset({1, 3})}

    def test_postnikov_figure_16_1_certifies_the_inverse_map(self):
        # Both fixed-point colors: 6 (white, coloop) in every entry, 4
        # (black, loop) in none.
        necklace = gn.postnikov_figure_16_1()
        assert necklace.entries == (
            frozenset({1, 2, 6}),
            frozenset({2, 3, 6}),
            frozenset({1, 3, 6}),
            frozenset({1, 5, 6}),
            frozenset({1, 5, 6}),
            frozenset({1, 2, 6}),
        )
        assert necklace.coloops == {6}
        assert necklace.loops == {4}

    def test_oh_worked_example_certifies_the_basis_construction(self):
        # A nontrivial necklace with a repeated entry; the page lists the
        # six bases and notes the necklace round-trips.
        necklace = gn.oh_worked_example()
        assert necklace.entries[1] == necklace.entries[3]
        expected = {
            frozenset(basis)
            for basis in [
                {1, 2, 4},
                {1, 2, 5},
                {1, 3, 4},
                {1, 3, 5},
                {2, 4, 5},
                {3, 4, 5},
            ]
        }
        assert necklace.to_positroid().bases == expected
        assert gn.GrassmannNecklace.from_matroid(necklace.to_positroid()) == necklace

    def test_lam_example_6_1_certifies_the_window_map(self):
        necklace = gn.lam_example_6_1()
        assert necklace.entries == (
            frozenset({1, 3}),
            frozenset({2, 3}),
            frozenset({3, 4}),
            frozenset({4, 6}),
            frozenset({5, 6}),
            frozenset({1, 6}),
        )
        assert necklace.to_bounded_affine_permutation() == (2, 4, 6, 5, 7, 9)

    def test_kls_example_3_14_certifies_the_juggling_dictionary(self):
        # Ball number 2, states recomputed from the window (2, 3, 5, 8),
        # siteswap 4112 up to rotation.
        necklace = gn.kls_example_3_14()
        assert necklace.rank == 2
        assert necklace.juggling_states == (
            frozenset({1, 4}),
            frozenset({1, 3}),
            frozenset({1, 2}),
            frozenset({1, 2}),
        )
        siteswap = necklace.siteswap
        rotations = {
            siteswap[shift:] + siteswap[:shift] for shift in range(len(siteswap))
        }
        assert (4, 1, 1, 2) in rotations

    def test_the_page_non_example_is_rejected_naming_n1(self):
        # (N1) genuinely constrains beyond equal cardinality: every entry
        # has two elements, yet the transition at i = 1 is illegal.
        entries = [{1, 2}, {1, 3}, {1, 3}]
        assert {len(entry) for entry in entries} == {2}
        with pytest.raises(ValueError, match=r"\(N1\)"):
            gn.GrassmannNecklace.from_entries((1, 2, 3), entries)

    @pytest.mark.parametrize("n", [0, 1, 3])
    def test_edge_types_have_a_single_necklace(self, n):
        # Page edge-case block: k = 0 gives (0, ..., 0) and k = n gives
        # ([n], ..., [n]).
        assert gn.enumerate_necklaces(0, n) == [gn.constant_necklace(n, ())]
        assert gn.enumerate_necklaces(n, n) == [
            gn.constant_necklace(n, range(1, n + 1))
        ]


# --------------------------------------------------------------------------- #
# Constructor validation: each named condition rejected with its name
# --------------------------------------------------------------------------- #
class TestConstructorValidation:
    def test_violating_the_fixed_transition_is_rejected_naming_n2(self):
        # 1 is not in I_1, so I_2 must equal I_1.
        with pytest.raises(ValueError, match=r"\(N2\)"):
            gn.GrassmannNecklace.from_entries((1, 2), [{2}, {1}])

    def test_wrong_entry_count_is_rejected(self):
        with pytest.raises(ValueError, match="one entry per"):
            gn.GrassmannNecklace.from_entries((1, 2), [{1}])

    def test_unknown_entry_label_is_rejected(self):
        with pytest.raises(ValueError, match="not in the ground set"):
            gn.GrassmannNecklace.from_entries((1, 2), [{3}, {3}])

    def test_decorated_permutation_size_mismatch_is_rejected(self):
        decorated = ps.DecoratedPermutation((2, 1))
        with pytest.raises(ValueError, match="ground set has"):
            gn.GrassmannNecklace.from_decorated_permutation((1, 2, 3), decorated)

    def test_window_length_mismatch_is_rejected(self):
        with pytest.raises(ValueError, match="one value per"):
            gn.GrassmannNecklace.from_bounded_affine_permutation((1, 2), (2, 3, 4))

    def test_window_out_of_bounds_is_rejected(self):
        with pytest.raises(ValueError, match=r"i <= f\(i\) <= i \+ n"):
            gn.GrassmannNecklace.from_bounded_affine_permutation((1, 2, 3), (1, 2, 7))

    def test_window_with_repeated_residues_is_rejected(self):
        with pytest.raises(ValueError, match="bijection"):
            gn.GrassmannNecklace.from_bounded_affine_permutation((1, 2, 3), (2, 2, 6))

    def test_interval_rank_position_out_of_range_is_rejected(self):
        with pytest.raises(ValueError, match="positions must be in"):
            gn.uniform_necklace(1, 3).interval_rank(0, 2)

    def test_circular_bruhat_comparison_across_types_is_rejected(self):
        with pytest.raises(ValueError, match="one type"):
            gn.uniform_necklace(1, 4).circular_bruhat_leq(gn.uniform_necklace(2, 4))

    def test_bad_uniform_and_enumeration_ranks_are_rejected(self):
        with pytest.raises(ValueError, match="0 <= k <= n"):
            gn.uniform_necklace(3, 2)
        with pytest.raises(ValueError, match="0 <= k <= n"):
            gn.enumerate_necklaces(-1, 2)

    def test_constant_necklace_with_unknown_label_is_rejected(self):
        with pytest.raises(ValueError, match="not in the ground set"):
            gn.constant_necklace(3, {5})


# --------------------------------------------------------------------------- #
# Structural theorems
# --------------------------------------------------------------------------- #
class TestStructuralTheorems:
    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
    def test_every_matroid_has_a_necklace_postnikov_lemma_16_3(self, n):
        # Postnikov Lemma 16.3: the Gale-minimal shifted bases of *any*
        # matroid satisfy (N1)-(N2), with type (rank, n).
        for matroid in mt.enumerate_matroids(n):
            necklace = gn.GrassmannNecklace.from_matroid(matroid)
            revalidated = gn.GrassmannNecklace.from_entries(
                matroid.elements, necklace.entries
            )
            assert revalidated == necklace
            assert necklace.necklace_type == (matroid.rank(), n)

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_oh_theorem_round_trip_on_necklaces(self, necklace):
        # Oh 2011: I(M(I)) = I for every Grassmann necklace I.
        assert gn.GrassmannNecklace.from_matroid(necklace.to_positroid()) == necklace

    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
    def test_positroid_envelope_contains_the_matroid_arw_prop_4_4(self, n):
        # ARW Prop. 4.4: every basis of M is a basis of M(I(M)), with
        # equality exactly for positroids.
        for matroid in mt.enumerate_matroids(n):
            envelope = gn.GrassmannNecklace.from_matroid(matroid).to_positroid()
            assert matroid.bases <= envelope.bases
            assert (matroid.bases == envelope.bases) == ps.is_positroid(matroid)

    @pytest.mark.parametrize("n", [0, 1, 2, 3])
    def test_positroid_envelope_is_smallest_arw_prop_4_4(self, n):
        # ARW Prop. 4.4, minimality: any positroid containing all bases of
        # M contains all bases of the envelope. Positroids on the matroid's
        # own ground set are reached through decorated permutations.
        elements = tuple(range(n))
        candidates = []
        for targets in itertools.permutations(range(1, n + 1)):
            fixed = [i for i in range(1, n + 1) if targets[i - 1] == i]
            for count in range(len(fixed) + 1):
                for clockwise in itertools.combinations(fixed, count):
                    decorated = ps.DecoratedPermutation(targets, frozenset(clockwise))
                    candidates.append(
                        ps.Positroid.from_decorated_permutation(
                            elements, decorated, validate=False
                        )
                    )
        for matroid in mt.enumerate_matroids(n):
            envelope = gn.GrassmannNecklace.from_matroid(matroid).to_positroid()
            for candidate in candidates:
                if matroid.bases <= candidate.bases:
                    assert envelope.bases <= candidate.bases

    @settings(max_examples=50, deadline=None)
    @given(st.data())
    def test_lemma_16_2_round_trip_from_decorated_permutations(self, data):
        # Postnikov Lemma 16.2: pi -> I(pi) -> pi is the identity.
        n = data.draw(st.integers(0, 5))
        decorated = _draw_decorated_permutation(data.draw, n)
        necklace = gn.GrassmannNecklace.from_decorated_permutation(
            tuple(range(1, n + 1)), decorated
        )
        assert necklace.to_decorated_permutation() == decorated

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_lemma_16_2_round_trip_from_necklaces(self, necklace):
        # Postnikov Lemma 16.2: I -> pi(I) -> I is the identity.
        rebuilt = gn.GrassmannNecklace.from_decorated_permutation(
            necklace.elements, necklace.to_decorated_permutation()
        )
        assert rebuilt == necklace

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_lam_theorem_6_2_window_round_trip(self, necklace):
        # Lam Thm 6.2: the window map and its inverse compose to the
        # identity, and windows satisfy the KLS section 3.2 bounds.
        window = necklace.to_bounded_affine_permutation()
        n = len(necklace.elements)
        assert all(a <= f <= a + n for a, f in enumerate(window, start=1))
        rebuilt = gn.GrassmannNecklace.from_bounded_affine_permutation(
            necklace.elements, window
        )
        assert rebuilt == necklace

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_kls_window_reduction_recovers_the_decorated_permutation(self, necklace):
        # KLS section 3.2: reduce f mod n and color fixed points by
        # f(i) = i (loop) vs f(i) = i + n (coloop). The reduction is in
        # Postnikov's direction — the inverse of the stored targets.
        n = len(necklace.elements)
        window = necklace.to_bounded_affine_permutation()
        decorated = necklace.to_decorated_permutation()
        for a, landing in enumerate(window, start=1):
            if landing == a:
                assert a in decorated.counterclockwise_fixed
            elif landing == a + n:
                assert a in decorated.clockwise_fixed
            else:
                assert decorated.targets[(landing - 1) % n] == a

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_juggling_states_satisfy_the_kls_condition(self, necklace):
        # KLS section 3 overview: J_{r+1} contains (J_r \ {1}) - 1, the
        # cyclic subtraction acting on landing times.
        states = necklace.juggling_states
        n = len(necklace.elements)
        for r, state in enumerate(states):
            following = states[(r + 1) % n]
            decremented = {(s - 2) % n + 1 for s in state if s != 1}
            assert decremented <= following

    @settings(max_examples=25, deadline=None)
    @given(necklaces(max_n=5))
    def test_oh_theorem_19_upper_necklace_cuts_the_same_positroid(self, necklace):
        # Oh Thm 19: the intersection of the cyclically shifted *dual*
        # Schubert matroids of the upper necklace is the positroid.
        n = len(necklace.elements)
        k = necklace.rank
        upper = necklace.upper_necklace
        cut = {
            frozenset(combo)
            for combo in itertools.combinations(range(1, n + 1), k)
            if all(_gale_leq(frozenset(combo), upper[i], i + 1, n) for i in range(n))
        }
        assert cut == necklace.to_positroid().bases

    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_loops_and_coloops_match_the_positroid(self, necklace):
        # Correspondences block: elements in no entry are the loops of the
        # positroid, elements in every entry its coloops.
        positroid = necklace.to_positroid()
        assert necklace.loops == positroid.loops
        assert necklace.coloops == positroid.coloops

    @settings(max_examples=25, deadline=None)
    @given(necklaces(max_n=5), st.integers(0, 5))
    def test_cyclic_shift_stays_a_necklace_arw_lemma_3_3(self, necklace, steps):
        # ARW Lemma 3.3 oracle: the necklace of the shifted positroid must
        # still satisfy (N1)-(N2), on the rotated ground set.
        shifted = necklace.cyclic_shift(steps)
        revalidated = gn.GrassmannNecklace.from_entries(
            shifted.elements, shifted.entries
        )
        assert revalidated == shifted
        assert shifted.necklace_type == necklace.necklace_type

    def test_lam_order_is_a_partial_order(self):
        # Lam section 6.3 order on all type-(2,4) necklaces: reflexive,
        # antisymmetric, and transitive.
        population = gn.enumerate_necklaces(2, 4)
        for necklace in population:
            # An equal-value copy, so reflexivity is not a self-comparison.
            copy = gn.GrassmannNecklace(necklace.elements, necklace.entry_masks)
            assert copy <= necklace
        for first, second in itertools.permutations(population, 2):
            if first <= second <= first:
                pytest.fail(f"antisymmetry fails on {first!r}, {second!r}")
        for first, second, third in itertools.product(population, repeat=3):
            if first <= second <= third:
                assert first <= third

    def test_lattice_path_matroids_are_positroids_oh_lemma_21(self):
        # Oh Lemma 21: {H : I <= H <= J} in the Gale order is a positroid
        # for every comparable pair, here exhaustively for k = 2, n = 4.
        n, k = 4, 2
        subsets = [frozenset(c) for c in itertools.combinations(range(1, n + 1), k)]
        pairs = [
            (low, high)
            for low in subsets
            for high in subsets
            if _gale_leq(low, high, 1, n)
        ]
        assert pairs
        for low, high in pairs:
            bases = [
                combo
                for combo in subsets
                if _gale_leq(low, combo, 1, n) and _gale_leq(combo, high, 1, n)
            ]
            lattice_path = ps.Positroid.from_bases(tuple(range(1, n + 1)), bases)
            assert ps.is_positroid(lattice_path) is True


# --------------------------------------------------------------------------- #
# Enumeration (OEIS A046802)
# --------------------------------------------------------------------------- #
class TestEnumeration:
    @pytest.mark.parametrize(
        ("n", "row"),
        [
            (0, [1]),
            (1, [1, 1]),
            (2, [1, 3, 1]),
            (3, [1, 7, 7, 1]),
            (4, [1, 15, 33, 15, 1]),
            (5, [1, 31, 131, 131, 31, 1]),
        ],
    )
    def test_counts_match_oeis_a046802(self, n, row):
        counts = []
        for rank in range(n + 1):
            found = gn.enumerate_necklaces(rank, n)
            assert len(set(found)) == len(found)
            counts.append(len(found))
        assert counts == row

    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
    def test_necklaces_biject_with_positroids(self, n):
        # Oh's theorem / Postnikov Thm 17.1: necklaces of type (k, n) are
        # equinumerous with rank-k positroids on [n], via I(M).
        positroids = ps.enumerate_positroids(n)
        from_positroids = {
            gn.GrassmannNecklace.from_matroid(positroid) for positroid in positroids
        }
        assert len(from_positroids) == len(positroids)
        assert from_positroids == set(_all_necklaces(n))


# --------------------------------------------------------------------------- #
# DataFrame serialization and experiment IO
# --------------------------------------------------------------------------- #
class TestDataFrames:
    @settings(max_examples=50, deadline=None)
    @given(necklaces())
    def test_dataframe_round_trip(self, necklace):
        decoded = gn.GrassmannNecklace.from_dataframe(necklace.to_dataframe())
        assert decoded == necklace
        assert decoded.elements == necklace.elements

    def test_experiment_io_round_trip_preserves_the_cyclic_order(self, tmp_path):
        necklace = gn.oh_worked_example()
        path = io.write_result(necklace.to_dataframe(), tmp_path / "necklace.json")
        decoded = gn.GrassmannNecklace.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == necklace
        assert decoded.elements == necklace.elements

    def test_missing_columns_are_rejected(self):
        with pytest.raises(ValueError, match="missing required columns"):
            gn.GrassmannNecklace.from_dataframe(pd.DataFrame({"element": [1]}))

    def test_empty_frame_decodes_to_the_empty_necklace(self):
        decoded = gn.GrassmannNecklace.from_dataframe(pd.DataFrame())
        assert decoded == gn.GrassmannNecklace((), ())


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestVisualization:
    def test_plot_membership_grid_draws_on_provided_axes(self):
        _, ax = plt.subplots()
        try:
            returned = gn.uniform_necklace(2, 4).plot_membership_grid(ax)
            assert returned is ax
            assert len(ax.images) == 1
        finally:
            plt.close("all")

    def test_plot_membership_grid_creates_axes_when_omitted(self):
        try:
            ax = gn.oh_worked_example().plot_membership_grid()
            assert ax.get_title() == "Grassmann necklace membership"
        finally:
            plt.close("all")

    def test_plot_juggling_pattern_draws_one_arc_per_throw(self):
        _, ax = plt.subplots()
        try:
            # The KLS example has no empty seconds: four throws, four arcs.
            returned = gn.kls_example_3_14().plot_juggling_pattern(ax)
            assert returned is ax
            assert len(ax.patches) == 4
        finally:
            plt.close("all")

    def test_plot_juggling_pattern_titles_the_siteswap(self):
        try:
            ax = gn.kls_example_3_14().plot_juggling_pattern()
            assert ax.get_title() == "Juggling pattern (siteswap 1,1,2,4)"
        finally:
            plt.close("all")
