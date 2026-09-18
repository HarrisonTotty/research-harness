"""Tests for research.decorated_permutation, from the Logseq Decorated Permutation page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Postnikov Lemma 16.2,
Theorem 17.1, Proposition 17.10 and Lemma 17.6, the affinization bijection
and length of Fomin-Williams-Zelevinsky Lemma 7.9.6 and Definition 7.9.10,
the Ardila-Rincon-Williams connectivity corollaries, and the enumeration
block); round-trip laws come from the API contract; and every part of the
definition has a rejection test naming it.

Direction discipline: the page states Postnikov's statistics for Postnikov's
permutation, while :mod:`research.positroid` and
:mod:`research.grassmann_necklace` exchange decorated permutations in the
Ardila-Rincon-Williams direction — the inverse. Every cross-check against
those modules therefore passes ``x.inverse()``.

Not transcribed in general, because the structure they quantify over cannot
be enumerated yet (spec, transformation backlog): the fundamental theorem of
reduced plabic graphs and the existence corollary (certified only in the
easy direction, on one square move), the general face count (certified only
on the two reduced graphs the library ships), the circular Bruhat order's
covers, its Bruhat interval, and its order reversal
(the order itself). The asymptotic statements are checked at finite ``n``
with explicit tolerances: ``e`` against its exact error bound, ``1/e`` and
``1/e^2`` along Callan's recurrence, where the observed error shrinks like
``1/n``.
"""

import itertools
import math
from collections.abc import Collection, Iterator, Sequence
from fractions import Fraction

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from experiments import io
from research import boundary_measurement as bm
from research import decorated_permutation as dp
from research import grassmann_necklace as gn
from research import matroid as mt
from research import positroid as ps

matplotlib.use("Agg")

POSTNIKOV = dp.NecklaceConvention.POSTNIKOV
ARW = dp.NecklaceConvention.ARW

# OEIS A046802, rows n = 0..6, as printed on the page.
A046802_ROWS: tuple[tuple[int, ...], ...] = (
    (1,),
    (1, 1),
    (1, 3, 1),
    (1, 7, 7, 1),
    (1, 15, 33, 15, 1),
    (1, 31, 131, 131, 31, 1),
    (1, 63, 473, 883, 473, 63, 1),
)
# OEIS A000522 as printed on the page.
A000522 = (1, 2, 5, 16, 65, 326, 1957, 13700)
# OEIS A075834 (offset 0) as printed on the page.
A075834 = (1, 1, 1, 2, 7, 34, 206, 1476)

type _Series = list[list[Fraction]]

TYPES = [(k, n) for n in range(1, 7) for k in range(1, n + 1)]


# --------------------------------------------------------------------------- #
# Strategies and helpers
# --------------------------------------------------------------------------- #
@st.composite
def decorated_permutations(
    draw: st.DrawFn, min_n: int = 0, max_n: int = 7
) -> dp.DecoratedPermutation:
    n = draw(st.integers(min_n, max_n))
    labels = tuple(range(1, n + 1))
    targets = tuple(draw(st.permutations(labels))) if n else ()
    fixed = [i for i in labels if targets[i - 1] == i]
    clockwise = (
        frozenset(draw(st.sets(st.sampled_from(fixed)))) if fixed else frozenset()
    )
    return dp.DecoratedPermutation(targets, clockwise)


def _all(n: int) -> list[dp.DecoratedPermutation]:
    return list(dp.enumerate_decorated_permutations(n))


def _of_type(k: int, n: int) -> list[dp.DecoratedPermutation]:
    return [x for x in _all(n) if x.permutation_type == (k, n)]


def _positroid(x: dp.DecoratedPermutation) -> ps.Positroid[int]:
    """Return the positroid Postnikov's ``x`` indexes (the library takes ARW)."""
    return ps.Positroid.from_decorated_permutation(
        tuple(range(1, x.size + 1)), x.inverse()
    )


def _dimension_counts(items: list[dp.DecoratedPermutation]) -> tuple[int, ...]:
    dimensions = [x.dimension for x in items]
    return tuple(dimensions.count(d) for d in range(max(dimensions) + 1))


def _evaluate(coefficients: tuple[int, ...], q: int) -> int:
    return sum(c * q**power for power, c in enumerate(coefficients))


def _hat(k: int, n: int) -> tuple[int, ...]:
    """Return ``q^{k-n} E_{k,n}(q)``; the division must be exact."""
    coefficients = dp.q_eulerian_polynomial(k, n)
    assert not any(coefficients[: n - k])
    return coefficients[n - k :]


def _classical_eulerian(k: int, n: int) -> int:
    return sum((-1) ** i * math.comb(n + 1, i) * (k - i) ** n for i in range(k + 1))


def _insert_fixed_point(
    x: dp.DecoratedPermutation, position: int, *, clockwise: bool
) -> dp.DecoratedPermutation:
    def lift(value: int) -> int:
        return value + 1 if value >= position else value

    targets = [position] * (x.size + 1)
    for i, target in enumerate(x.targets, start=1):
        targets[lift(i) - 1] = lift(target)
    colored = {lift(i) for i in x.clockwise_fixed} | (
        {position} if clockwise else set()
    )
    return dp.DecoratedPermutation(tuple(targets), frozenset(colored))


def _cyclic_interval(start: int, end: int, n: int) -> list[int]:
    """Walk the closed cyclic interval ``[start, end]`` clockwise, label by label."""
    walk = [start]
    while walk[-1] != end:
        walk.append(walk[-1] % n + 1)
    return walk


def _positroid_bases_of_rank(k: int, n: int) -> set[frozenset[frozenset[int]]]:
    """Enumerate rank-``k`` positroids on ``[n]`` with no decorated permutation.

    Every family of ``k``-subsets is tried as a basis system and kept when it
    satisfies the basis axioms and Oh's positroid test.
    """
    elements = tuple(range(1, n + 1))
    subsets = [frozenset(c) for c in itertools.combinations(elements, k)]
    found: set[frozenset[frozenset[int]]] = set()
    for size in range(1, len(subsets) + 1):
        for family in itertools.combinations(subsets, size):
            matroid = _matroid_or_none(elements, family)
            if matroid is not None and ps.is_positroid(matroid):
                found.add(frozenset(family))
    return found


def _matroid_or_none(
    elements: tuple[int, ...], family: Sequence[frozenset[int]]
) -> mt.Matroid[int] | None:
    try:
        return mt.Matroid.from_bases(elements, family)
    except ValueError:
        return None


def _bounded_windows(n: int) -> list[tuple[int, ...]]:
    """Enumerate bounded affine windows from the definition alone."""
    return [
        window
        for window in itertools.product(*(range(i, i + n + 1) for i in range(1, n + 1)))
        if sorted((value - 1) % n for value in window) == list(range(n))
    ]


def _sif_count(n: int) -> int:
    return sum(
        1
        for targets in itertools.permutations(range(1, n + 1))
        if dp.DecoratedPermutation(targets).is_stabilized_interval_free
    )


def _callan_counts(size: int) -> list[int]:
    """Extend 1, 1, 1 by Callan's recurrence up to index ``size``."""
    counts = [1, 1, 1]
    for n in range(3, size + 1):
        counts.append(
            (n - 1) * counts[n - 1]
            + sum((j - 1) * counts[j] * counts[n - j] for j in range(2, n - 1))
        )
    return counts


def _set_partitions(items: list[int]) -> Iterator[list[list[int]]]:
    if not items:
        yield []
        return
    head, rest = items[0], items[1:]
    for partition in _set_partitions(rest):
        yield [[head], *partition]
        for index in range(len(partition)):
            yield [
                *partition[:index],
                [head, *partition[index]],
                *partition[index + 1 :],
            ]


def _is_noncrossing(partition: Sequence[Collection[int]]) -> bool:
    for first, second in itertools.combinations(partition, 2):
        for a, c in itertools.combinations(sorted(first), 2):
            if any(a < b < c for b in second) and any(d < a or d > c for d in second):
                return False
    return True


# --------------------------------------------------------------------------- #
# Canonical examples: each asserts what the page says the example certifies
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_postnikov_section_16_example_is_the_page_data(self):
        example = dp.postnikov_section_16_example()
        assert example.targets == (3, 1, 5, 4, 2, 6)
        assert example.colors == {4: 1, 6: -1}

    def test_postnikov_section_16_necklace_begins_and_ends_as_stated(self):
        necklace = dp.postnikov_section_16_example().grassmann_necklace(POSTNIKOV)
        assert necklace[0] == frozenset({1, 2, 6})
        assert necklace[5] == frozenset({6, 1, 2})

    def test_postnikov_section_16_white_point_in_every_entry_black_in_none(self):
        necklace = dp.postnikov_section_16_example().grassmann_necklace(POSTNIKOV)
        assert all(6 in entry for entry in necklace)
        assert all(4 not in entry for entry in necklace)

    def test_postnikov_section_16_necklace_maps_back_to_the_example(self):
        example = dp.postnikov_section_16_example()
        rebuilt = dp.DecoratedPermutation.from_grassmann_necklace(
            example.grassmann_necklace(POSTNIKOV), convention=POSTNIKOV
        )
        assert rebuilt == example

    def test_fwz_example_7_9_2_has_anti_exceedances_1_4_and_overlined_3(self):
        example = dp.fwz_example_7_9_2()
        assert example.targets == (5, 2, 3, 6, 4, 1)
        assert example.anti_exceedance_count == 3
        assert example.anti_exceedances == frozenset({1, 4, 3})

    def test_fwz_example_7_9_2_has_one_fixed_point_of_each_color(self):
        example = dp.fwz_example_7_9_2()
        assert example.clockwise_fixed == frozenset({3})
        assert example.counterclockwise_fixed == frozenset({2})

    def test_fwz_definition_7_4_20_example_has_a_single_overlined_fixed_point(self):
        example = dp.fwz_definition_7_4_20_example()
        assert example.targets == (3, 4, 5, 1, 2, 6)
        assert example.size == 6
        assert example.fixed_points == frozenset({6})
        assert example.clockwise_fixed == frozenset({6})

    def test_arw_section_4_3_example_is_the_page_data(self):
        example = dp.arw_section_4_3_example()
        assert example.targets == (1, 7, 9, 3, 2, 6, 5, 10, 4, 8)
        assert example.counterclockwise_fixed == frozenset({1})
        assert example.clockwise_fixed == frozenset({6})

    def test_arw_section_4_3_example_has_d_4_weak_excedances_on_n_10(self):
        example = dp.arw_section_4_3_example()
        assert example.size == 10
        assert example.weak_excedance_count == 4
        assert {len(entry) for entry in example.grassmann_necklace(ARW)} == {4}

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(6) for k in range(n + 1)]
    )
    def test_top_cell_is_the_unique_alignment_free_one_of_its_type(self, k, n):
        alignment_free = [x for x in _of_type(k, n) if x.alignment_number == 0]
        assert alignment_free == [dp.top_cell(k, n)]

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(7) for k in range(n + 1)]
    )
    def test_top_cell_maps_i_to_i_plus_k_and_has_dimension_k_n_minus_k(self, k, n):
        top = dp.top_cell(k, n)
        assert top.targets == tuple((i + k - 1) % n + 1 for i in range(1, n + 1))
        assert top.permutation_type == (k, n)
        assert top.dimension == k * (n - k)

    @pytest.mark.parametrize("n", [1, 2, 5])
    def test_top_cell_fixed_points_are_black_for_k_0_and_white_for_k_n(self, n):
        everything = frozenset(range(1, n + 1))
        assert dp.top_cell(0, n).counterclockwise_fixed == everything
        assert dp.top_cell(n, n).clockwise_fixed == everything

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(6) for k in range(n + 1)]
    )
    def test_zero_dimensional_cells_are_the_binom_n_k_minimal_elements(self, k, n):
        expected = {
            dp.zero_dimensional_cell(n, white)
            for white in itertools.combinations(range(1, n + 1), k)
        }
        minimal = {x for x in _of_type(k, n) if x.dimension == 0}
        assert minimal == expected
        assert len(minimal) == math.comb(n, k)
        assert {x.alignment_number for x in minimal} == {k * (n - k)}

    @given(
        st.integers(0, 6).flatmap(
            lambda n: st.tuples(
                st.just(n), st.sets(st.integers(1, n)) if n else st.just(set())
            )
        )
    )
    def test_zero_dimensional_cell_has_the_constant_necklace(self, case):
        n, white = case
        cell = dp.zero_dimensional_cell(n, white)
        assert cell.grassmann_necklace(POSTNIKOV) == (frozenset(white),) * n
        assert _positroid(cell).bases == frozenset({frozenset(white)})

    def test_gr_2_4_has_33_decorated_permutations(self):
        assert len(_of_type(2, 4)) == 33

    def test_gr_2_4_has_10_decorated_permutations_up_to_cyclic_rotation(self):
        orbits = {
            frozenset(x.cyclic_shift(steps) for steps in range(4))
            for x in _of_type(2, 4)
        }
        assert len(orbits) == 10

    def test_gr_2_4_is_graded_by_dimension_as_williams_table_1(self):
        # A_{2,4}(q) = q^4 + 4q^3 + 10q^2 + 12q + 6.
        assert _dimension_counts(_of_type(2, 4)) == (6, 12, 10, 4, 1)
        assert dp.dimension_generating_polynomial(2, 4) == (6, 12, 10, 4, 1)

    def test_counting_fixture_n_6_total_and_split_by_type(self):
        everything = _all(6)
        by_type = [x.anti_exceedance_count for x in everything]
        assert len(everything) == 1957
        assert tuple(by_type.count(k) for k in range(7)) == (
            1,
            63,
            473,
            883,
            473,
            63,
            1,
        )

    def test_non_example_decorating_a_non_fixed_point_is_rejected(self):
        with pytest.raises(ValueError, match="only fixed points carry a decoration"):
            dp.DecoratedPermutation((2, 1), frozenset({1}))


# --------------------------------------------------------------------------- #
# Derived vocabulary
# --------------------------------------------------------------------------- #
class TestDerivedVocabulary:
    @given(decorated_permutations())
    def test_anti_exceedance_count_is_the_fwz_position_count(self, x):
        # FWZ Definition 7.9.1: #{i : pi(i) < i or pi(i) is overlined}.
        positions = {i for i, t in enumerate(x.targets, start=1) if t < i}
        assert x.anti_exceedance_count == len(positions | x.clockwise_fixed)

    @given(decorated_permutations())
    def test_anti_exceedances_and_weak_excedances_are_exchanged_by_inversion(self, x):
        assert x.inverse().anti_exceedances == x.weak_excedances
        assert x.inverse().weak_excedances == x.anti_exceedances

    def test_anti_exceedance_and_weak_excedance_counts_generally_differ(self):
        three_cycle = dp.DecoratedPermutation((2, 3, 1))
        assert three_cycle.anti_exceedance_count == 1
        assert three_cycle.weak_excedance_count == 2

    @given(decorated_permutations())
    def test_number_of_weak_k_excedances_is_independent_of_k(self, x):
        # ARW Definition 4.5; the entries I_k are the weak k-excedance sets.
        sizes = {len(entry) for entry in x.grassmann_necklace(ARW)}
        assert sizes <= {x.weak_excedance_count}

    def test_a_two_cycle_is_a_crossing_by_the_closed_interval_definition(self):
        # i = 1, j = 2: pi(j) = 1 lies in [1, 2] and j = 2 lies in [2, 1].
        assert dp.DecoratedPermutation((2, 1)).crossings == ((1, 2),)

    def test_disjoint_two_cycles_cross_only_themselves(self):
        assert dp.DecoratedPermutation((2, 1, 4, 3)).crossings == ((1, 2), (3, 4))

    def test_every_pair_of_chords_of_the_gr_2_4_top_cell_crosses(self):
        assert dp.top_cell(2, 4).crossings == (
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4),
        )

    @given(decorated_permutations(min_n=2))
    def test_crossings_are_exactly_the_pairs_meeting_postnikovs_condition(self, x):
        n = x.size
        expected = set()
        for i, j in itertools.permutations(range(1, n + 1), 2):
            image_i, image_j = x.targets[i - 1], x.targets[j - 1]
            if i == image_i or j == image_j:
                continue
            if image_j in _cyclic_interval(i, image_i, n) and j in _cyclic_interval(
                image_i, i, n
            ):
                expected.add((min(i, j), max(i, j)))
        assert set(x.crossings) == expected

    @given(decorated_permutations(min_n=2))
    def test_alignments_are_exactly_the_pairs_meeting_postnikovs_condition(self, x):
        n = x.size
        expected = set()
        for i, j in itertools.permutations(range(1, n + 1), 2):
            image_i, image_j = x.targets[i - 1], x.targets[j - 1]
            black_or_moved = i != image_i or x.colors[i] == 1
            white_or_moved = j != image_j or x.colors[j] == -1
            if (
                black_or_moved
                and white_or_moved
                and image_i in _cyclic_interval(i, image_j, n)
                and j in _cyclic_interval(image_j, i, n)
            ):
                expected.add((i, j))
        assert set(x.alignments) == expected
        assert len(x.alignments) == len(expected)

    @given(decorated_permutations())
    def test_no_pair_of_chords_is_both_crossing_and_aligned(self, x):
        # Derived consistency check, not a page statement: the cover relation
        # trades a crossing for an alignment, so no pair can be both.
        aligned = {(min(i, j), max(i, j)) for i, j in x.alignments}
        assert not aligned & set(x.crossings)

    @given(decorated_permutations())
    def test_a_loop_never_participates_in_a_crossing(self, x):
        touched = {i for pair in x.crossings for i in pair}
        assert not touched & x.fixed_points

    @given(decorated_permutations())
    def test_every_counterclockwise_loop_aligns_with_every_clockwise_loop(self, x):
        # Postnikov section 17; his counterclockwise loop is a black fixed
        # point and his clockwise loop a white one.
        expected = set(itertools.product(x.counterclockwise_fixed, x.clockwise_fixed))
        assert expected <= set(x.alignments)

    @given(decorated_permutations())
    def test_fixed_points_enter_alignments_only_in_their_permitted_role(self, x):
        for i, j in x.alignments:
            assert i not in x.clockwise_fixed
            assert j not in x.counterclockwise_fixed

    @given(decorated_permutations())
    def test_alignment_number_equals_the_length_of_the_affinization(self, x):
        # FWZ Definition 7.9.10: inversion classes correspond to alignments.
        assert x.affine_length == x.alignment_number

    @given(decorated_permutations(max_n=6), st.data())
    def test_adding_an_uncounted_fixed_point_adds_exactly_k_alignments(self, x, data):
        # Williams section 5: her "clockwise" fixed point is the uncounted
        # color — counterclockwise in this library.
        position = data.draw(st.integers(1, x.size + 1))
        grown = _insert_fixed_point(x, position, clockwise=False)
        assert grown.alignment_number == x.alignment_number + x.anti_exceedance_count

    @given(decorated_permutations())
    def test_regular_means_every_fixed_point_is_counted(self, x):
        assert x.is_regular == (not x.counterclockwise_fixed)

    @given(decorated_permutations())
    def test_colors_send_clockwise_to_minus_one_and_the_rest_to_one(self, x):
        assert {i for i, c in x.colors.items() if c == -1} == x.clockwise_fixed
        assert {i for i, c in x.colors.items() if c == 1} == x.counterclockwise_fixed

    def test_stabilized_interval_free_matches_the_definition_on_examples(self):
        assert dp.DecoratedPermutation((2, 3, 1)).is_stabilized_interval_free
        # (2, 1, 3) stabilizes the proper intervals {1, 2} and {3}.
        assert not dp.DecoratedPermutation((2, 1, 3)).is_stabilized_interval_free


# --------------------------------------------------------------------------- #
# Structural theorems
# --------------------------------------------------------------------------- #
class TestNecklaceBijection:
    """Postnikov Lemma 16.2 and Ardila-Rincon-Williams Proposition 4.6."""

    @given(decorated_permutations(), st.sampled_from([POSTNIKOV, ARW]))
    def test_necklace_maps_are_inverse_on_decorated_permutations(self, x, convention):
        necklace = x.grassmann_necklace(convention)
        rebuilt = dp.DecoratedPermutation.from_grassmann_necklace(
            necklace, convention=convention
        )
        assert rebuilt == x

    @pytest.mark.parametrize("n", range(6))
    @pytest.mark.parametrize("convention", [POSTNIKOV, ARW])
    def test_necklace_map_is_a_bijection_onto_all_necklaces_of_size_n(
        self, n, convention
    ):
        images = [x.grassmann_necklace(convention) for x in _all(n)]
        every_necklace = {
            necklace.entries
            for rank in range(n + 1)
            for necklace in gn.enumerate_necklaces(rank, n)
        }
        assert len(set(images)) == len(images)
        assert set(images) == every_necklace

    @given(decorated_permutations())
    def test_type_k_n_corresponds_to_necklaces_of_type_k_n(self, x):
        k, n = x.permutation_type
        assert all(len(entry) == k for entry in x.grassmann_necklace(POSTNIKOV))
        assert len(x.grassmann_necklace(POSTNIKOV)) == n

    @given(decorated_permutations())
    def test_d_weak_excedances_correspond_to_necklaces_of_type_d_n(self, x):
        d = x.weak_excedance_count
        assert all(len(entry) == d for entry in x.grassmann_necklace(ARW))

    @given(decorated_permutations())
    def test_postnikov_shifted_anti_exceedance_rule_gives_the_first_entry(self, x):
        assert x.size == 0 or x.grassmann_necklace(POSTNIKOV)[0] == x.anti_exceedances

    @given(decorated_permutations())
    def test_arw_map_uses_the_inverse_permutation_of_postnikovs(self, x):
        assert x.grassmann_necklace(ARW) == x.inverse().grassmann_necklace(POSTNIKOV)

    @given(decorated_permutations())
    def test_black_fixed_points_in_no_entry_and_white_in_all(self, x):
        for entry in x.grassmann_necklace(POSTNIKOV):
            assert x.clockwise_fixed <= entry
            assert not x.counterclockwise_fixed & entry

    @given(decorated_permutations())
    def test_arw_necklace_agrees_with_the_grassmann_necklace_module(self, x):
        necklace = gn.GrassmannNecklace.from_decorated_permutation(
            tuple(range(1, x.size + 1)), x
        )
        assert x.grassmann_necklace(ARW) == necklace.entries


class TestCellIndexing:
    """Postnikov Theorem 17.1 and Proposition 17.10."""

    @pytest.mark.parametrize("n", range(5))
    def test_distinct_decorated_permutations_index_distinct_cells(self, n):
        cells = [_positroid(x).bases for x in _all(n)]
        assert len(set(cells)) == len(cells)

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 5) for k in range(n + 1)]
    )
    def test_type_k_n_decorated_permutations_biject_with_rank_k_positroids(self, k, n):
        cells = [_positroid(x).bases for x in _of_type(k, n)]
        assert len(set(cells)) == len(cells)
        assert set(cells) == _positroid_bases_of_rank(k, n)

    @given(decorated_permutations(max_n=6))
    def test_type_k_n_indexes_a_cell_of_the_rank_k_grassmannian(self, x):
        k, n = x.permutation_type
        positroid = _positroid(x)
        assert positroid.rank() == k
        assert positroid.size == n

    @given(decorated_permutations(min_n=1, max_n=6))
    def test_dimension_is_k_n_minus_k_minus_alignments(self, x):
        k, n = x.permutation_type
        assert x.dimension == k * (n - k) - x.alignment_number
        assert x.dimension == bm.cell_dimension(_positroid(x))

    @given(decorated_permutations())
    def test_dimension_lies_between_zero_and_k_n_minus_k(self, x):
        k, n = x.permutation_type
        assert 0 <= x.dimension <= k * (n - k)

    @pytest.mark.parametrize(("k", "n"), [(1, 3), (2, 4), (2, 5)])
    def test_top_cell_indexes_the_uniform_positroid(self, k, n):
        # Derived cross-check, not a page statement: the unique maximum of
        # the closure order is the dense cell, whose positroid is uniform.
        # uniform_positroid lives on 0..n-1; relabel onto 1..n.
        uniform = {
            frozenset(e + 1 for e in basis)
            for basis in ps.uniform_positroid(k, n).bases
        }
        assert _positroid(dp.top_cell(k, n)).bases == uniform


class TestAffinization:
    """FWZ Definition 7.9.3 and Lemma 7.9.6; Knutson-Lam-Speyer section 3.2."""

    def test_affinization_of_fwz_example_7_9_2_follows_the_four_cases(self):
        # (5, _2, ^3, 6, 4, 1): 5 > 1; underlined 2; overlined 3 -> 3 + 6;
        # 6 > 4; 4 < 5 -> 4 + 6; 1 < 6 -> 1 + 6.
        window = dp.fwz_example_7_9_2().to_bounded_affine_permutation()
        assert window == (5, 2, 9, 6, 10, 7)

    @given(decorated_permutations())
    def test_affinization_is_a_k_n_bounded_affine_permutation(self, x):
        k, n = x.permutation_type
        window = x.to_bounded_affine_permutation()
        assert all(i <= value <= i + n for i, value in enumerate(window, start=1))
        assert sorted((value - 1) % n + 1 for value in window) == list(range(1, n + 1))
        assert sum(value - i for i, value in enumerate(window, start=1)) == k * n

    @given(decorated_permutations())
    def test_reduction_modulo_n_inverts_the_affinization(self, x):
        window = x.to_bounded_affine_permutation()
        assert dp.DecoratedPermutation.from_bounded_affine_permutation(window) == x

    @pytest.mark.parametrize("n", range(1, 6))
    def test_affinization_inverts_reduction_on_every_bounded_window(self, n):
        for window in _bounded_windows(n):
            rebuilt = dp.DecoratedPermutation.from_bounded_affine_permutation(window)
            assert rebuilt.to_bounded_affine_permutation() == window

    @pytest.mark.parametrize("n", range(1, 6))
    def test_affinization_is_a_bijection_onto_the_k_n_bounded_windows(self, n):
        images = [x.to_bounded_affine_permutation() for x in _all(n)]
        assert len(set(images)) == len(images)
        assert set(images) == set(_bounded_windows(n))
        for k in range(n + 1):
            of_average_k = {
                w for w in _bounded_windows(n) if sum(w) - n * (n + 1) // 2 == k * n
            }
            assert of_average_k == {
                x.to_bounded_affine_permutation() for x in _of_type(k, n)
            }

    def test_kls_coloring_counts_f_i_equal_to_i_plus_n(self):
        # f(1) = 1 is the uncounted color, f(2) = 2 + 2 the counted one.
        rebuilt = dp.DecoratedPermutation.from_bounded_affine_permutation((1, 4))
        assert rebuilt == dp.DecoratedPermutation((1, 2), frozenset({2}))

    @given(decorated_permutations())
    def test_affinization_agrees_with_the_grassmann_necklace_module(self, x):
        necklace = gn.GrassmannNecklace.from_decorated_permutation(
            tuple(range(1, x.size + 1)), x.inverse()
        )
        assert x.to_bounded_affine_permutation() == (
            necklace.to_bounded_affine_permutation()
        )

    @pytest.mark.parametrize("n", range(1, 6))
    def test_affine_length_vanishes_exactly_for_the_top_permutation(self, n):
        # Derived, not a transcription of the face count: the length equals
        # the alignment number (FWZ Definition 7.9.10) and the top element
        # i -> i + k is the unique alignment-free one (Postnikov Lemma 17.6).
        for x in _all(n):
            assert (x.affine_length == 0) == (x == dp.top_cell(*x.permutation_type))


class TestFaceCount:
    """FWZ Corollaries 7.10.5 and 7.10.7, certified on the library's reduced graphs.

    The general statement quantifies over reduced plabic graphs, which the
    library cannot yet enumerate per decorated permutation; these are
    instance checks.
    """

    def test_square_graph_attains_the_face_bound_at_the_top_permutation(self):
        network = bm.square_network()
        decorated = dp.DecoratedPermutation(network.trip_permutation())
        a, b = decorated.permutation_type
        assert network.is_reduced()
        assert decorated == dp.top_cell(2, 4)
        assert decorated.affine_length == 0
        assert network.to_perfect_orientation().face_count == a * (b - a) + 1

    def test_lollipop_graph_has_a_b_minus_a_minus_length_plus_one_faces(self):
        # Lam Example 4.2: lollipops with I(Pi) = {3, 4}, so 3 and 4 are the
        # counted (white-lollipop, coloop) fixed points.
        network = bm.lollipop_network()
        decorated = dp.zero_dimensional_cell(4, {3, 4})
        a, b = decorated.permutation_type
        assert network.is_reduced()
        assert network.trip_permutation() == decorated.targets
        # Tie the decoration to the graph: its positroid's decorated
        # permutation (library direction, hence the inverse).
        assert network.to_positroid().to_decorated_permutation().inverse() == decorated
        assert network.to_perfect_orientation().face_count == (
            a * (b - a) - decorated.affine_length + 1
        )
        assert network.to_perfect_orientation().face_count < a * (b - a) + 1


class TestFundamentalTheorem:
    """Postnikov Theorem 13.4, easy direction, certified on one square move.

    The converse and the existence statement (Corollary 14.7) quantify over
    all reduced plabic graphs and are backlog.
    """

    def test_a_square_move_keeps_the_decorated_trip_permutation(self):
        before = bm.square_network()
        after = before.square_move(("T", "R", "B", "L"))
        assert before.is_reduced()
        assert after.is_reduced()
        # Fixed-point free, so the trip permutation needs no decoration.
        assert after.trip_permutation() == before.trip_permutation() == (3, 4, 1, 2)
        assert (
            after.to_positroid().to_decorated_permutation()
            == before.to_positroid().to_decorated_permutation()
        )


class TestConnectivity:
    """Ardila-Rincon-Williams Corollaries 7.9 and 7.11."""

    @given(decorated_permutations(min_n=1))
    def test_no_stabilized_cyclic_interval_iff_stabilized_interval_free(self, x):
        assert x.stabilizes_proper_cyclic_interval == (
            not x.is_stabilized_interval_free
        )

    @given(decorated_permutations(min_n=1, max_n=6))
    def test_positroid_is_connected_iff_no_proper_cyclic_interval_is_stabilized(
        self, x
    ):
        connected = len(_positroid(x).connected_components()) == 1
        assert connected == (not x.stabilizes_proper_cyclic_interval)

    @given(decorated_permutations(min_n=1, max_n=6))
    def test_noncrossing_partition_is_the_partition_into_components(self, x):
        assert x.noncrossing_partition == _positroid(x).connected_components()

    @given(decorated_permutations())
    def test_noncrossing_partition_is_noncrossing_and_keeps_i_with_its_image(self, x):
        partition = x.noncrossing_partition
        assert _is_noncrossing([sorted(block) for block in partition])
        assert sorted(i for block in partition for i in block) == list(
            range(1, x.size + 1)
        )
        for i, target in enumerate(x.targets, start=1):
            assert any({i, target} <= block for block in partition)

    @pytest.mark.parametrize("n", range(1, 6))
    def test_noncrossing_partition_is_the_finest_such_partition(self, n):
        candidates = [
            [frozenset(block) for block in partition]
            for partition in _set_partitions(list(range(1, n + 1)))
            if _is_noncrossing(partition)
        ]
        for targets in itertools.permutations(range(1, n + 1)):
            x = dp.DecoratedPermutation(targets)
            closed = [
                partition
                for partition in candidates
                if all(
                    any({i, t} <= block for block in partition)
                    for i, t in enumerate(targets, start=1)
                )
            ]
            for partition in closed:
                assert all(
                    any(fine <= coarse for coarse in partition)
                    for fine in x.noncrossing_partition
                )

    @pytest.mark.parametrize(
        ("n", "expected"), [(2, 1), (3, 2), (4, 7), (5, 34), (6, 206)]
    )
    def test_connected_positroids_are_equinumerous_with_sif_permutations(
        self, n, expected
    ):
        connected = [x for x in _all(n) if not x.stabilizes_proper_cyclic_interval]
        sif = [
            targets
            for targets in itertools.permutations(range(1, n + 1))
            if dp.DecoratedPermutation(targets).is_stabilized_interval_free
        ]
        assert len(connected) == len(sif) == expected == A075834[n]

    def test_both_single_element_positroids_are_connected(self):
        # The page's caveat: the ARW sequence starts 2, departing from
        # A075834 only at n = 1.
        connected = [x for x in _all(1) if not x.stabilizes_proper_cyclic_interval]
        assert len(connected) == 2
        assert A075834[1] == 1

    def test_brute_force_sif_counts_are_a075834_and_satisfy_callans_recurrence(self):
        counts = [_sif_count(n) for n in range(8)]
        assert tuple(counts) == A075834
        assert counts == _callan_counts(7)

    def test_connected_proportions_approach_one_over_e_and_one_over_e_squared(self):
        # ARW Theorems 10.6-10.7. Along Callan's recurrence the error shrinks
        # like 1/n (about 0.0019 and 0.0007 at n = 200), hence the tolerances.
        counts = _callan_counts(200)
        errors_e, errors_e2 = [], []
        for n in (50, 100, 200):
            errors_e.append(abs(Fraction(counts[n], math.factorial(n)) - 1 / math.e))
            errors_e2.append(
                abs(
                    Fraction(counts[n], dp.count_decorated_permutations(n))
                    - math.exp(-2)
                )
            )
        assert errors_e[2] < errors_e[1] < errors_e[0]
        assert errors_e2[2] < errors_e2[1] < errors_e2[0]
        assert errors_e[2] < 2e-3
        assert errors_e2[2] < 1e-3


class TestOperations:
    @given(decorated_permutations())
    def test_inverse_is_an_involution_keeping_the_decoration(self, x):
        assert x.inverse().inverse() == x
        assert x.inverse().clockwise_fixed == x.clockwise_fixed

    @given(decorated_permutations())
    def test_inverse_composes_with_the_permutation_to_the_identity(self, x):
        inverse = x.inverse()
        assert all(
            inverse.targets[t - 1] == i for i, t in enumerate(x.targets, start=1)
        )

    @given(decorated_permutations(max_n=4), decorated_permutations(max_n=4))
    def test_direct_sum_restricts_to_each_summand(self, a, b):
        total = a.direct_sum(b)
        n = a.size
        assert total.targets[:n] == a.targets
        assert tuple(t - n for t in total.targets[n:]) == b.targets
        assert total.clockwise_fixed == a.clockwise_fixed | {
            i + n for i in b.clockwise_fixed
        }

    @given(
        decorated_permutations(min_n=1, max_n=3),
        decorated_permutations(min_n=1, max_n=3),
    )
    def test_direct_sum_indexes_the_direct_sum_of_positroids(self, a, b):
        # ARW Proposition 7.8, for the non-crossing partition {[n], [n+1, n+m]};
        # library direction throughout (the direct sum commutes with inversion).
        first = ps.Positroid.from_decorated_permutation(tuple(range(1, a.size + 1)), a)
        second = ps.Positroid.from_decorated_permutation(
            tuple(range(a.size + 1, a.size + b.size + 1)), b
        )
        total = ps.Positroid.from_matroid(first.direct_sum(second))
        assert total.to_decorated_permutation() == a.direct_sum(b)

    @given(decorated_permutations(min_n=1), st.integers(-8, 8))
    def test_cyclic_shift_rotates_every_chord(self, x, steps):
        n = x.size
        shifted = x.cyclic_shift(steps)
        for i, target in enumerate(x.targets, start=1):
            assert shifted.targets[(i - 1 + steps) % n] == (target - 1 + steps) % n + 1
        assert shifted.clockwise_fixed == {
            (i - 1 + steps) % n + 1 for i in x.clockwise_fixed
        }

    @given(decorated_permutations())
    def test_cyclic_shift_by_n_is_the_identity(self, x):
        assert x.cyclic_shift(x.size) == x
        assert x.cyclic_shift(0) == x

    @given(decorated_permutations(), st.integers(0, 7))
    def test_cyclic_shift_preserves_type_and_alignment_number(self, x, steps):
        shifted = x.cyclic_shift(steps)
        assert shifted.permutation_type == x.permutation_type
        assert shifted.alignment_number == x.alignment_number

    @given(decorated_permutations(), st.integers(0, 7))
    def test_cyclic_shift_rotates_the_necklace(self, x, steps):
        n = x.size
        necklace = x.grassmann_necklace(POSTNIKOV)
        shifted = x.cyclic_shift(steps).grassmann_necklace(POSTNIKOV)
        for r in range(n):
            assert shifted[(r + steps) % n] == {
                (i - 1 + steps) % n + 1 for i in necklace[r]
            }


class TestEnumeration:
    @pytest.mark.parametrize("n", range(7))
    def test_enumeration_yields_each_decorated_permutation_once(self, n):
        everything = _all(n)
        assert len(set(everything)) == len(everything) == A000522[n]

    @pytest.mark.parametrize("n", range(8))
    def test_count_matches_a000522(self, n):
        assert dp.count_decorated_permutations(n) == A000522[n]

    @pytest.mark.parametrize("n", range(1, 12))
    def test_count_satisfies_postnikovs_recurrence(self, n):
        assert dp.count_decorated_permutations(0) == 1
        assert dp.count_decorated_permutations(n) == (
            n * dp.count_decorated_permutations(n - 1) + 1
        )

    def test_count_over_n_factorial_tends_to_e(self):
        # e - N_n / n! = sum_{k > n} 1/k! < 2 / (n + 1)!.
        n = 15
        ratio = dp.count_decorated_permutations(n) / math.factorial(n)
        assert ratio == pytest.approx(
            math.e, rel=0, abs=2 / math.factorial(n + 1) + 1e-15
        )

    @pytest.mark.parametrize("n", range(8))
    def test_exponential_generating_function_is_e_to_the_x_over_one_minus_x(self, n):
        # [x^n] e^x / (1 - x) is the Cauchy product sum_{j <= n} 1/j!.
        coefficient = sum(Fraction(1, math.factorial(j)) for j in range(n + 1))
        assert coefficient * math.factorial(n) == A000522[n]

    @pytest.mark.parametrize("n", range(8))
    def test_ordinary_generating_function_is_arws_p_of_x(self, n):
        # [x^n] x^k / (1 - x)^{k+1} = binom(n, k).
        assert (
            sum(math.factorial(k) * math.comb(n, k) for k in range(n + 1))
            == (A000522[n])
        )

    def test_bivariate_generating_function_expands_to_the_a046802_rows(self):
        # e^{xy} (x - 1) / (x - e^{y(x-1)}) = e^{xy} / (1 - S) with
        # S = sum_{m >= 1} y^m (x - 1)^{m-1} / m!, expanded as a power series
        # in y whose coefficients are polynomials in x (index = power of x).
        order = 6

        def multiply(a: _Series, b: _Series) -> _Series:
            product = [[Fraction(0)] * (order + 1) for _ in range(order + 1)]
            for p, row_a in enumerate(a):
                for q, row_b in enumerate(b[: order + 1 - p]):
                    for s, u in enumerate(row_a):
                        for t, v in enumerate(row_b[: order + 1 - s]):
                            product[p + q][s + t] += u * v
            return product

        def zero() -> _Series:
            return [[Fraction(0)] * (order + 1) for _ in range(order + 1)]

        s_series = zero()
        for m in range(1, order + 1):
            for t in range(m):
                s_series[m][t] = Fraction(
                    math.comb(m - 1, t) * (-1) ** (m - 1 - t), math.factorial(m)
                )
        exponential = zero()
        for m in range(order + 1):
            exponential[m][m] = Fraction(1, math.factorial(m))
        geometric, power = zero(), zero()
        geometric[0][0] = power[0][0] = Fraction(1)
        for _ in range(order):
            power = multiply(power, s_series)
            geometric = [
                [g + p for g, p in zip(row_g, row_p, strict=True)]
                for row_g, row_p in zip(geometric, power, strict=True)
            ]
        series = multiply(exponential, geometric)
        for n in range(order + 1):
            row = tuple(series[n][k] * math.factorial(n) for k in range(n + 1))
            assert row == tuple(Fraction(c) for c in A046802_ROWS[n])
            assert not any(series[n][n + 1 :])

    @pytest.mark.parametrize("n", range(7))
    def test_counts_by_type_are_the_a046802_row(self, n):
        counts = [x.anti_exceedance_count for x in _all(n)]
        assert tuple(counts.count(k) for k in range(n + 1)) == A046802_ROWS[n]

    @pytest.mark.parametrize("n", range(7))
    def test_counts_by_weak_excedances_are_the_same_row(self, n):
        # ARW Proposition 4.6 counts by weak excedances; inversion matches them.
        counts = [x.weak_excedance_count for x in _all(n)]
        assert tuple(counts.count(k) for k in range(n + 1)) == A046802_ROWS[n]

    @pytest.mark.parametrize(("k", "n"), TYPES)
    def test_williams_closed_form_matches_the_triangle(self, k, n):
        assert dp.count_decorated_permutations_of_type(k, n) == A046802_ROWS[n][k]

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(7) for k in range(n + 1)]
    )
    def test_counts_by_type_are_binomial_sums_of_eulerian_numbers(self, k, n):
        total = sum(
            math.comb(n, r) * _classical_eulerian(k, n - r) for r in range(n + 1)
        )
        assert total == A046802_ROWS[n][k]

    @pytest.mark.parametrize(("k", "n"), TYPES)
    def test_dimension_generating_polynomial_counts_cells_by_dimension(self, k, n):
        assert dp.dimension_generating_polynomial(k, n) == _dimension_counts(
            _of_type(k, n)
        )

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(2, 9) for k in range(1, n)]
    )
    def test_dimension_generating_polynomial_is_dual_symmetric(self, k, n):
        assert dp.dimension_generating_polynomial(k, n) == (
            dp.dimension_generating_polynomial(n - k, n)
        )

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 9) for k in range(1, n + 1)]
    )
    def test_euler_characteristic_is_one(self, k, n):
        assert _evaluate(dp.dimension_generating_polynomial(k, n), -1) == 1

    @pytest.mark.parametrize(("k", "n"), TYPES)
    def test_q_eulerian_polynomial_is_the_regular_part(self, k, n):
        regular = [x for x in _of_type(k, n) if x.is_regular]
        assert dp.q_eulerian_polynomial(k, n) == _dimension_counts(regular)

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 9) for k in range(1, n + 1)]
    )
    def test_dimension_polynomial_is_a_binomial_sum_of_q_eulerian_polynomials(
        self, k, n
    ):
        total = [0] * (k * (n - k) + 1)
        for i in range(n - k + 1):
            for power, c in enumerate(dp.q_eulerian_polynomial(k, n - i)):
                total[power] += math.comb(n, i) * c
        assert tuple(total) == dp.dimension_generating_polynomial(k, n)

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 9) for k in range(1, n + 1)]
    )
    def test_normalized_q_eulerian_interpolates_eulerian_narayana_binomial(self, k, n):
        hat = _hat(k, n)
        assert _evaluate(hat, 1) == _classical_eulerian(k, n)
        assert _evaluate(hat, 0) == math.comb(n, k) * math.comb(n, k - 1) // n
        assert abs(_evaluate(hat, -1)) == math.comb(n - 1, k - 1)

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 7) for k in range(1, n + 1)]
    )
    def test_classical_eulerian_formula_counts_weak_excedances(self, k, n):
        count = sum(
            1
            for targets in itertools.permutations(range(1, n + 1))
            if sum(t >= i for i, t in enumerate(targets, start=1)) == k
        )
        assert _classical_eulerian(k, n) == count

    @pytest.mark.parametrize(
        ("k", "n"), [(k, n) for n in range(1, 9) for k in range(1, n + 1)]
    )
    def test_normalized_q_eulerian_is_symmetric_in_k(self, k, n):
        assert _hat(k, n) == _hat(n + 1 - k, n)

    @pytest.mark.parametrize("n", range(1, 7))
    def test_permanent_coefficients_count_by_type(self, n):
        # Williams section 7: 1 + x on the diagonal, x above it, 1 below.
        coefficients = [0] * (n + 1)
        for targets in itertools.permutations(range(n)):
            term = [1]
            for i, j in enumerate(targets):
                factor = [1, 1] if i == j else [0, 1] if j > i else [1]
                term = [
                    sum(
                        term[a] * factor[p - a]
                        for a in range(len(term))
                        if 0 <= p - a < len(factor)
                    )
                    for p in range(len(term) + len(factor) - 1)
                ]
            for power, c in enumerate(term):
                coefficients[power] += c
        assert tuple(coefficients) == A046802_ROWS[n]
        for k in range(1, n + 1):
            assert coefficients[k] == _evaluate(
                dp.dimension_generating_polynomial(k, n), 1
            )

    def test_moments_and_free_cumulants_of_one_plus_exp_one(self):
        # ARW Theorem 11.1: moments = positroid counts, free cumulants =
        # connected positroid counts, tied by the non-crossing
        # moment-cumulant recursion m_n = sum_s k_s sum m_{i_1} ... m_{i_s}.
        size = 6
        cumulants = [0] + [
            sum(1 for x in _all(n) if not x.stabilizes_proper_cyclic_interval)
            for n in range(1, size + 1)
        ]
        moments = [1]
        for n in range(1, size + 1):
            total = 0
            for s in range(1, n + 1):
                for parts in itertools.product(range(n - s + 1), repeat=s):
                    if sum(parts) == n - s:
                        total += cumulants[s] * math.prod(moments[i] for i in parts)
            moments.append(total)
        assert moments == [dp.count_decorated_permutations(n) for n in range(size + 1)]
        # Moments of 1 + Exp(1): E[(1 + X)^n] = sum binom(n, j) j!.
        assert moments == [
            sum(math.comb(n, j) * math.factorial(j) for j in range(n + 1))
            for n in range(size + 1)
        ]


# --------------------------------------------------------------------------- #
# Round-trip laws
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(decorated_permutations())
    def test_dataframe_round_trip(self, x):
        assert dp.DecoratedPermutation.from_dataframe(x.to_dataframe()) == x

    @given(decorated_permutations(min_n=1), st.randoms(use_true_random=False))
    def test_dataframe_rows_may_come_in_any_order(self, x, random):
        frame = x.to_dataframe()
        order = list(frame.index)
        random.shuffle(order)
        assert dp.DecoratedPermutation.from_dataframe(frame.loc[order]) == x

    def test_dataframe_has_one_row_per_position(self):
        frame = dp.fwz_example_7_9_2().to_dataframe()
        assert list(frame.columns) == ["position", "target", "decoration"]
        assert frame["position"].tolist() == [1, 2, 3, 4, 5, 6]
        assert frame["target"].tolist() == [5, 2, 3, 6, 4, 1]
        assert frame["decoration"].tolist() == [
            None,
            "counterclockwise",
            "clockwise",
            None,
            None,
            None,
        ]

    @pytest.mark.parametrize(
        "example",
        [dp.fwz_example_7_9_2(), dp.top_cell(2, 4), dp.DecoratedPermutation(())],
    )
    def test_experiment_io_round_trip(self, example, tmp_path):
        path = io.write_result(example.to_dataframe(), tmp_path / "decorated.json")
        decoded = dp.DecoratedPermutation.from_dataframe(
            pd.read_json(path, dtype=False)
        )
        assert decoded == example

    @given(decorated_permutations())
    def test_constructors_agree(self, x):
        by_colors = dp.DecoratedPermutation.from_colors(x.targets, x.colors)
        by_lines = dp.DecoratedPermutation.from_underline_overline(
            x.targets,
            overlined=x.clockwise_fixed,
            underlined=x.counterclockwise_fixed,
        )
        by_window = dp.DecoratedPermutation.from_bounded_affine_permutation(
            x.to_bounded_affine_permutation()
        )
        assert by_colors == by_lines == by_window == x

    @given(decorated_permutations())
    def test_value_semantics(self, x):
        twin = dp.DecoratedPermutation(x.targets, x.clockwise_fixed)
        assert twin == x
        assert hash(twin) == hash(x)


# --------------------------------------------------------------------------- #
# Rejections: every part of the definition, named
# --------------------------------------------------------------------------- #
class TestRejections:
    @pytest.mark.parametrize("targets", [(1, 1), (0, 1), (2, 3), (1, 2, 4)])
    def test_non_bijection_is_rejected(self, targets):
        with pytest.raises(ValueError, match=r"must be a bijection of \[n\]"):
            dp.DecoratedPermutation(targets)

    def test_decoration_outside_the_ground_set_is_rejected(self):
        with pytest.raises(ValueError, match="only fixed points carry a decoration"):
            dp.DecoratedPermutation((1, 2), frozenset({3}))

    def test_color_outside_plus_minus_one_is_rejected(self):
        with pytest.raises(ValueError, match=r"values in \{1, -1\}"):
            dp.DecoratedPermutation.from_colors((1, 2), {1: 1, 2: 0})

    def test_uncolored_fixed_point_is_rejected(self):
        with pytest.raises(ValueError, match="defined on exactly the fixed points"):
            dp.DecoratedPermutation.from_colors((1, 2), {1: 1})

    def test_colored_non_fixed_point_is_rejected(self):
        with pytest.raises(ValueError, match="defined on exactly the fixed points"):
            dp.DecoratedPermutation.from_colors((2, 1, 3), {1: 1, 3: -1})

    def test_fixed_point_with_both_lines_is_rejected(self):
        with pytest.raises(
            ValueError, match="exactly one of an overline or an underline"
        ):
            dp.DecoratedPermutation.from_underline_overline(
                (1, 2), overlined={1, 2}, underlined={2}
            )

    def test_fixed_point_with_no_line_is_rejected(self):
        with pytest.raises(
            ValueError, match="exactly one of an overline or an underline"
        ):
            dp.DecoratedPermutation.from_underline_overline((1, 2), overlined={1})

    @pytest.mark.parametrize("window", [(0, 2), (1, 5), (4, 2)])
    def test_unbounded_window_is_rejected(self, window):
        with pytest.raises(ValueError, match=r"i <= f\(i\) <= i \+ n"):
            dp.DecoratedPermutation.from_bounded_affine_permutation(window)

    def test_window_with_repeated_residues_is_rejected(self):
        with pytest.raises(ValueError, match="bijection of \\[n\\] modulo n"):
            dp.DecoratedPermutation.from_bounded_affine_permutation((2, 4, 5))

    def test_necklace_entry_outside_the_ground_set_is_rejected(self):
        with pytest.raises(ValueError, match=r"subsets of \[n\]"):
            dp.DecoratedPermutation.from_grassmann_necklace(
                [{3}, {3}], convention=POSTNIKOV
            )

    def test_necklace_violating_the_exchange_condition_is_rejected(self):
        # 1 is in I_1 = {1, 2}, so I_2 must be ({2}) + {j}; {3} is not.
        with pytest.raises(ValueError, match=r"\(N1\)"):
            dp.DecoratedPermutation.from_grassmann_necklace(
                [{1, 2}, {3}, {3}], convention=POSTNIKOV
            )

    def test_necklace_changing_without_its_index_is_rejected(self):
        # 1 is not in I_1 = {2}, so I_2 must equal I_1.
        with pytest.raises(ValueError, match="must equal it"):
            dp.DecoratedPermutation.from_grassmann_necklace([{2}, {1}], convention=ARW)

    def test_frame_missing_a_column_is_rejected(self):
        with pytest.raises(ValueError, match="missing columns"):
            dp.DecoratedPermutation.from_dataframe(pd.DataFrame({"position": [1]}))

    def test_frame_with_a_gap_in_positions_is_rejected(self):
        frame = pd.DataFrame(
            {"position": [1, 3], "target": [3, 1], "decoration": [None, None]}
        )
        with pytest.raises(ValueError, match=r"each position 1\.\.n exactly once"):
            dp.DecoratedPermutation.from_dataframe(frame)

    def test_frame_with_an_unknown_decoration_is_rejected(self):
        frame = pd.DataFrame({"position": [1], "target": [1], "decoration": ["white"]})
        with pytest.raises(ValueError, match="decoration must be"):
            dp.DecoratedPermutation.from_dataframe(frame)

    def test_frame_with_an_undecorated_fixed_point_is_rejected(self):
        frame = pd.DataFrame({"position": [1], "target": [1], "decoration": [None]})
        with pytest.raises(ValueError, match="has no decoration"):
            dp.DecoratedPermutation.from_dataframe(frame)

    def test_frame_decorating_a_non_fixed_point_is_rejected(self):
        frame = pd.DataFrame(
            {
                "position": [1, 2],
                "target": [2, 1],
                "decoration": ["clockwise", None],
            }
        )
        with pytest.raises(ValueError, match="only fixed points carry a decoration"):
            dp.DecoratedPermutation.from_dataframe(frame)

    @pytest.mark.parametrize(("k", "n"), [(-1, 3), (4, 3)])
    def test_top_cell_outside_the_type_range_is_rejected(self, k, n):
        with pytest.raises(ValueError, match="0 <= k <= n"):
            dp.top_cell(k, n)

    def test_zero_dimensional_cell_of_negative_size_is_rejected(self):
        with pytest.raises(ValueError, match="non-negative"):
            dp.zero_dimensional_cell(-1, ())

    def test_negative_size_is_rejected_by_the_enumeration(self):
        with pytest.raises(ValueError, match="non-negative"):
            next(dp.enumerate_decorated_permutations(-1))

    def test_negative_size_is_rejected_by_the_count(self):
        with pytest.raises(ValueError, match="non-negative"):
            dp.count_decorated_permutations(-1)

    @pytest.mark.parametrize(
        "formula",
        [
            dp.count_decorated_permutations_of_type,
            dp.dimension_generating_polynomial,
            dp.q_eulerian_polynomial,
        ],
    )
    @pytest.mark.parametrize(("k", "n"), [(0, 3), (4, 3)])
    def test_williams_closed_forms_reject_types_outside_their_sum(self, formula, k, n):
        with pytest.raises(ValueError, match="1 <= k <= n"):
            formula(k, n)


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestVisualization:
    def test_plot_chord_diagram_draws_on_provided_axes(self):
        _, ax = plt.subplots()
        try:
            returned = dp.postnikov_section_16_example().plot_chord_diagram(ax)
            assert returned is ax
            # The boundary circle plus one loop per fixed point (4 and 6).
            assert len(ax.lines) == 3
        finally:
            plt.close("all")

    def test_plot_chord_diagram_draws_one_arrow_per_chord_and_loop(self):
        _, ax = plt.subplots()
        try:
            dp.postnikov_section_16_example().plot_chord_diagram(ax)
            arrows = [text for text in ax.texts if text.get_text() == ""]
            assert len(arrows) == 6
        finally:
            plt.close("all")

    def test_plot_chord_diagram_distinguishes_the_two_loop_senses(self):
        _, ax = plt.subplots()
        try:
            dp.postnikov_section_16_example().plot_chord_diagram(ax)
            styles = sorted(line.get_linestyle() for line in ax.lines[1:])
            assert styles == ["-", "--"]
        finally:
            plt.close("all")

    def test_plot_chord_diagram_creates_axes_when_omitted(self):
        try:
            ax = dp.top_cell(2, 4).plot_chord_diagram()
            assert ax.get_title() == "Chord diagram"
        finally:
            plt.close("all")

    def test_plot_chord_diagram_handles_the_empty_permutation(self):
        try:
            ax = dp.DecoratedPermutation(()).plot_chord_diagram()
            assert len(ax.lines) == 1
        finally:
            plt.close("all")
