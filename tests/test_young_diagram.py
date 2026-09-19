"""Tests for research.young_diagram, from the Logseq Young Diagram page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Stanley EC1 sections
1.7-1.8, 3.4, 3.21 and the cited exercises; Frame-Robinson-Thrall 1954;
Adin-Roichman; Postnikov sections 2.1, 2.3, 17, 19 and Theorem 6.5; the
dominance-order block); round-trip laws come from the API contract; and
every clause of each formulation has a rejection test naming it.

Not transcribed (spec, "Not transcribed"): Schubert cells, the closure order
and Proposition 17.4 (no Schubert-cell or cell-order structure; Schubert Cell
is a red link), Grassmannian permutations proper (red link; only the
inversion twin is tested), Schensted's theorem and the Robinson-Schensted
bijection itself (red link; tested as the two identities it implies), the
shifted analogue and disconnected diagrams (Shifted/Skew Young Diagram are
red links), q-quotients (Frame-Robinson-Thrall Lemma 2, Theorems 3-4: no
quotient structure), defect zero (a statement about representations; only
the fixture's hook-graph side is tested), the dominance lattice's
join-irreducibles (the page does not define "slippery step") and its Möbius
values through the atoms (Blass-Sagan Corollary 5.2, paraphrased only),
Early's H-steps and V-steps (his conventions are only partly quoted), the
non-uniqueness and product
statements about differential posets (not about Young diagrams alone), and
the generating-function identities as formal power series — those are tested
through the diagram statements that prove them.
"""

import collections
import functools
import itertools
import math
from collections.abc import Iterable
from typing import cast

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st
from matplotlib.lines import Line2D

from experiments import io
from research import le_diagram as ld
from research import young_diagram as yd

matplotlib.use("Agg")

Y = yd.YoungDiagram

# OEIS A000041, n = 0..12, as printed on the page.
A000041 = (1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77)
# OEIS A000085, n = 0..10, as printed on the page.
A000085 = (1, 1, 2, 4, 10, 26, 76, 232, 764, 2620, 9496)
# OEIS A000700, n = 0..16, as printed on the page.
A000700 = (1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5)
# Stanley EC1 Exercise 1.71, n = 0..7, as printed on the page.
COVER_COUNTS = (1, 2, 4, 7, 12, 19, 30, 45)
# OEIS A000108, m = 1..6, as printed on the page.
CATALAN = (1, 2, 5, 14, 42, 132)

BOXES = [(k, n) for n in range(9) for k in range(n + 1)]
SMALL_BOXES = [(k, n) for k, n in BOXES if n <= 7]


# --------------------------------------------------------------------------- #
# Strategies and helpers
# --------------------------------------------------------------------------- #
@functools.cache
def _partitions(n: int) -> tuple[yd.YoungDiagram, ...]:
    return tuple(yd.partitions(n))


def _up_to(n: int) -> list[yd.YoungDiagram]:
    return [d for m in range(n + 1) for d in _partitions(m)]


@st.composite
def diagrams(draw: st.DrawFn, max_rows: int = 6, max_width: int = 6) -> yd.YoungDiagram:
    lengths = draw(st.lists(st.integers(0, max_width), max_size=max_rows))
    return Y.from_rows(sorted(lengths, reverse=True))


@st.composite
def nested_pairs(draw: st.DrawFn) -> tuple[yd.YoungDiagram, yd.YoungDiagram]:
    """A pair ``(nu, lambda)`` with ``nu`` contained in ``lambda``."""
    upper = draw(diagrams())
    return upper.meet(draw(diagrams())), upper


@st.composite
def same_size_pairs(draw: st.DrawFn) -> tuple[yd.YoungDiagram, yd.YoungDiagram]:
    n = draw(st.integers(0, 9))
    options = st.sampled_from(_partitions(n))
    return draw(options), draw(options)


@st.composite
def diagrams_with_cell(draw: st.DrawFn) -> tuple[yd.YoungDiagram, int, int]:
    diagram = draw(diagrams().filter(lambda d: d.size > 0))
    i, j = draw(st.sampled_from(sorted(diagram.cells)))
    return diagram, i, j


@st.composite
def boxed_diagrams(draw: st.DrawFn) -> tuple[yd.YoungDiagram, int, int]:
    n = draw(st.integers(0, 8))
    k = draw(st.integers(0, n))
    lengths = draw(st.lists(st.integers(0, n - k), min_size=k, max_size=k))
    return Y.from_rows(sorted(lengths, reverse=True)), k, n


def _ydata(line: Line2D) -> list[float]:
    return [float(y) for y in cast("Iterable[float]", line.get_ydata())]


def _is_order_convex(cells: frozenset[tuple[int, int]]) -> bool:
    return all(
        (i, j) in cells
        for a, b in itertools.product(cells, repeat=2)
        for i in range(a[0], b[0] + 1)
        for j in range(a[1], b[1] + 1)
    )


def _components(
    cells: frozenset[tuple[int, int]], *, by_order: bool
) -> set[frozenset[tuple[int, int]]]:
    """Components under edge-adjacency, or under componentwise comparability."""

    def linked(a: tuple[int, int], b: tuple[int, int]) -> bool:
        if by_order:
            return (a[0] <= b[0] and a[1] <= b[1]) or (b[0] <= a[0] and b[1] <= a[1])
        return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1

    remaining, found = set(cells), set()
    while remaining:
        frontier = [remaining.pop()]
        component = set(frontier)
        while frontier:
            current = frontier.pop()
            near = {cell for cell in remaining if linked(current, cell)}
            remaining -= near
            component |= near
            frontier.extend(near)
        found.add(frozenset(component))
    return found


def _is_line_convex(cells: frozenset[tuple[int, int]]) -> bool:
    lines = [sorted(j for i, j in cells if i == row) for row in {i for i, _ in cells}]
    lines += [sorted(i for i, j in cells if j == col) for col in {j for _, j in cells}]
    return all(line == list(range(line[0], line[-1] + 1)) for line in lines)


def _is_down_closed(cells: frozenset[tuple[int, int]]) -> bool:
    return all(
        (i, j) in cells
        for a, b in cells
        for i in range(1, a + 1)
        for j in range(1, b + 1)
    )


def _major_index(tableau: tuple[tuple[int, ...], ...]) -> int:
    row_of = {entry: i for i, row in enumerate(tableau) for entry in row}
    return sum(i for i in range(1, len(row_of)) if row_of[i + 1] > row_of[i])


def _q_binomial(n: int, k: int) -> tuple[int, ...]:
    """Coefficients by the q-Pascal recurrence [n,k] = [n-1,k-1] + q^k [n-1,k]."""
    if k < 0 or k > n:
        return ()
    if n == 0:
        return (1,)
    total = [0] * (k * (n - k) + 1)
    for power, c in enumerate(_q_binomial(n - 1, k - 1)):
        total[power] += c
    for power, c in enumerate(_q_binomial(n - 1, k)):
        total[power + k] += c
    return tuple(total)


def _partitions_with_parts_at_most(n: int, k: int) -> int:
    """The coefficient of q^n in prod_{i<=k} 1/(1 - q^i), by expansion."""
    series = [1] + [0] * n
    for part in range(1, k + 1):
        for total in range(part, n + 1):
            series[total] += series[total - part]
    return series[n]


def _interval(lower: yd.YoungDiagram, upper: yd.YoungDiagram) -> list[yd.YoungDiagram]:
    return [d for d in _up_to(upper.size) if upper.contains(d) and d.contains(lower)]


@functools.cache
def _mobius_by_recursion(lower: yd.YoungDiagram, upper: yd.YoungDiagram) -> int:
    if lower == upper:
        return 1
    return -sum(
        _mobius_by_recursion(lower, mid)
        for mid in _interval(lower, upper)
        if mid != upper
    )


def _chain_lengths(n: int) -> tuple[int, int]:
    """Steps in the shortest and longest maximal chains of the dominance order."""

    @functools.cache
    def below(diagram: yd.YoungDiagram) -> tuple[int, int]:
        covered = diagram.dominance_lower_covers()
        if not covered:
            return (0, 0)
        lengths = [below(d) for d in covered]
        return (1 + min(s for s, _ in lengths), 1 + max(g for _, g in lengths))

    return below(Y((n,)))


NESTED_SMALL = [(nu, lam) for lam in _up_to(6) for nu in _interval(Y(), lam)]


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_empty_diagram_certifies_the_zero_case_conventions(self):
        empty = yd.empty_diagram()
        assert empty.cells == frozenset()
        assert empty.length == 0
        assert empty.count_standard_tableaux() == 1
        assert list(yd.partitions(0)) == [empty]

    def test_empty_diagram_is_the_bottom_covered_by_exactly_the_single_box(self):
        empty = yd.empty_diagram()
        assert all(d.contains(empty) for d in _up_to(6))
        assert empty.upper_covers() == (Y((1,)),)
        assert len(empty.addable_cells) == 1
        assert empty.corners == ()

    def test_all_in_one_conjugate_hooks_and_tableau_count(self):
        example = yd.all_in_one_diagram()
        assert example.rows == (4, 3, 1)
        assert example.conjugate() == Y((3, 2, 2, 1))
        assert example.hook_lengths == ((6, 4, 3, 1), (4, 2, 1), (1,))
        assert example.hook_product == 576
        assert example.count_standard_tableaux() == math.factorial(8) // 576 == 70

    def test_all_in_one_durfee_side_and_frobenius_coordinates(self):
        example = yd.all_in_one_diagram()
        assert example.durfee_side == 2
        assert example.frobenius_coordinates == ((3, 1), (2, 0))

    def test_all_in_one_boundary_word_in_the_3_by_5_box(self):
        example = yd.all_in_one_diagram()
        assert example.boundary_word(3, 8) == (1, 2, 1, 2, 1, 1, 2, 1)
        assert example.to_subset(3, 8) == frozenset({2, 4, 7})

    def test_all_in_one_example_tableau_is_standard_of_that_shape(self):
        tableau = ((1, 2, 5, 8), (3, 4, 6), (7,))
        assert tableau in set(yd.all_in_one_diagram().standard_tableaux())

    def test_frt_example_1_hook_graph_and_degree(self):
        example = yd.frt_example_1_diagram()
        assert example.rows == (6, 4, 2)
        assert example.hook_lengths == ((8, 7, 5, 4, 2, 1), (5, 4, 2, 1), (2, 1))
        assert example.count_standard_tableaux() == 2673
        assert math.factorial(12) // example.hook_product == 2673

    @pytest.mark.parametrize("p", [3, 11])
    def test_frt_example_1_has_no_hook_length_divisible_by(self, p):
        hooks = [h for row in yd.frt_example_1_diagram().hook_lengths for h in row]
        assert [h for h in hooks if h % p == 0] == []

    def test_frt_example_2_removing_the_6_hook_at_2_3(self):
        example = yd.frt_example_2_diagram()
        assert example.rows == (7, 6, 5, 3)
        assert example.hook_lengths == (
            (10, 9, 8, 6, 5, 3, 1),
            (8, 7, 6, 4, 3, 1),
            (6, 5, 4, 2, 1),
            (3, 2, 1),
        )
        assert example.hook_length(2, 3) == 6
        assert len(example.rim_hook(2, 3)) == 6
        smaller = example.remove_rim_hook(2, 3)
        assert smaller == Y((7, 4, 2, 2))
        assert smaller.cells == example.cells - example.rim_hook(2, 3)
        assert smaller.hook_lengths == (
            (10, 9, 6, 5, 3, 2, 1),
            (6, 5, 2, 1),
            (3, 2),
            (2, 1),
        )

    def test_frt_example_4_three_core(self):
        example = yd.frt_example_4_diagram()
        assert example.rows == (9, 7, 4, 3, 2)
        assert example.core(3) == Y((3, 1))

    def test_box_fixture_certifies_the_vertical_step_labeling(self):
        example = yd.postnikov_figure_2_1_diagram()
        assert example.rows == (4, 4, 2, 1)
        assert example.fits_in_box(4, 10)
        assert example.to_subset(4, 10) == frozenset({3, 4, 7, 9})
        assert example.to_le_diagram(4, 10).vertical_steps == frozenset({3, 4, 7, 9})
        assert example.size == 11

    def test_box_fixture_certifies_the_inverse_formula(self):
        subset = [3, 4, 7, 9]
        literal = tuple(
            len(set(range(i, 11)) - set(subset)) for i in subset
        )  # lambda_s = |[i_s, n] \ I|
        assert literal == (4, 4, 2, 1)
        assert Y.from_subset(subset, 10) == yd.postnikov_figure_2_1_diagram()

    def test_all_diagrams_in_the_2_by_3_box(self):
        box = list(yd.diagrams_in_box(2, 5))
        assert len(box) == math.comb(5, 2) == 10
        assert sorted(d.size for d in box) == [0, 1, 2, 2, 3, 3, 4, 4, 5, 6]
        assert yd.box_size_distribution(2, 5) == (1, 1, 2, 2, 2, 1, 1)

    def test_all_diagrams_in_the_2_by_2_box(self):
        assert len(list(yd.diagrams_in_box(2, 4))) == 6
        assert yd.box_size_distribution(2, 4) == (1, 1, 2, 1, 1)

    def test_self_conjugate_fixture_certifies_the_diagonal_hook_bijection(self):
        example = yd.stanley_figure_1_16_diagram()
        assert example.rows == (5, 4, 4, 3, 1)
        assert example.is_self_conjugate
        assert example.diagonal_hook_lengths == (9, 5, 3)
        assert sum(example.diagonal_hook_lengths) == example.size == 17

    def test_durfee_dissection_fixture(self):
        example = yd.stanley_figure_1_18_diagram()
        assert example.rows == (7, 5, 3, 3, 2)
        assert example.durfee_dissection() == (3, Y((4, 2)), Y((3, 2)))
        assert example.size == 20 == 9 + 6 + 5

    def test_frobenius_notation_fixture(self):
        example = yd.stanley_figure_1_33_diagram()
        assert example.rows == (7, 7, 4, 2, 1)
        assert example.frobenius_coordinates == ((6, 5, 1), (4, 2, 0))
        assert example.size == 21 == 3 + (6 + 5 + 1) + (4 + 2 + 0)

    def test_boundary_sequence_fixture(self):
        example = yd.adin_roichman_example_9_7_diagram()
        assert example.rows == (3, 1)
        assert example.boundary_sequence == (1, 0, 1, 1, 0)

    @pytest.mark.parametrize(("m", "catalan"), list(enumerate(CATALAN, start=1)))
    def test_two_equal_rows_are_counted_by_catalan_numbers(self, m, catalan):
        assert yd.rectangle(2, m).count_standard_tableaux() == catalan

    def test_square_and_staircase_tableau_counts(self):
        assert yd.rectangle(3, 3).count_standard_tableaux() == 42
        assert yd.staircase(3) == Y((3, 2, 1))
        assert yd.staircase(3).is_self_conjugate
        assert yd.staircase(3).count_standard_tableaux() == 16

    def test_first_incomparable_dominance_pairs_are_at_n_6(self):
        incomparable = {
            frozenset({a, b})
            for a, b in itertools.combinations(_partitions(6), 2)
            if not a.dominates(b) and not b.dominates(a)
        }
        assert incomparable == {
            frozenset({Y((3, 3)), Y((4, 1, 1))}),
            frozenset({Y((2, 2, 2)), Y((3, 1, 1, 1))}),
        }

    def test_incomparable_pairs_at_n_6_are_conjugates_of_one_another(self):
        first = {Y((3, 3)), Y((4, 1, 1))}
        assert {d.conjugate() for d in first} == {Y((2, 2, 2)), Y((3, 1, 1, 1))}

    def test_non_down_closed_cells_are_a_finite_cell_set_but_not_a_young_diagram(self):
        cells = yd.non_down_closed_cells()
        assert cells == frozenset({(1, 1), (2, 2)})
        assert not _is_down_closed(cells)
        with pytest.raises(ValueError, match="down-closure violated"):
            Y.from_cells(cells)

    def test_skew_non_example_is_order_convex_but_not_an_order_ideal(self):
        cells = yd.skew_non_example_cells()
        assert cells == Y((6, 4, 3, 1)).cells - Y((4, 2, 1)).cells
        assert _is_order_convex(cells)
        assert not _is_down_closed(cells)
        with pytest.raises(ValueError, match="down-closure violated"):
            Y.from_cells(cells)

    def test_non_decreasing_rows_violate_weak_decrease(self):
        assert yd.non_decreasing_rows() == (2, 3)
        with pytest.raises(ValueError, match="weak decrease violated"):
            Y.from_rows(yd.non_decreasing_rows())

    def test_special_shapes(self):
        assert yd.hook_shape(9, 3) == Y((6, 1, 1, 1))
        assert yd.staircase(4) == Y((4, 3, 2, 1))
        assert yd.rectangle(4, 6) == Y((6, 6, 6, 6))
        assert yd.rectangle(0, 3) == yd.rectangle(3, 0) == yd.empty_diagram()


# --------------------------------------------------------------------------- #
# Definition: equivalence of the formulations
# --------------------------------------------------------------------------- #
class TestFormulations:
    @given(diagrams())
    def test_cells_are_the_set_of_formulation_b(self, diagram):
        """Adin-Roichman Definition 2.10."""
        expected = {
            (i, j)
            for i in range(1, diagram.length + 1)
            for j in range(1, diagram.rows[i - 1] + 1)
        }
        assert diagram.cells == expected
        assert len(diagram.cells) == diagram.size

    @given(diagrams())
    def test_cells_are_down_closed(self, diagram):
        """Stanley EC1 section 1.5, p. 50."""
        assert _is_down_closed(diagram.cells)

    def test_every_down_closed_subset_of_a_small_grid_is_a_diagram(self):
        """Stanley EC1 section 1.5, p. 50."""
        grid = [(i, j) for i in range(1, 4) for j in range(1, 4)]
        down_closed = [
            frozenset(subset)
            for size in range(len(grid) + 1)
            for subset in itertools.combinations(grid, size)
            if _is_down_closed(frozenset(subset))
        ]
        assert {Y.from_cells(cells).cells for cells in down_closed} == set(down_closed)
        assert len(down_closed) == math.comb(6, 3)

    def test_trailing_zeros_are_ignored(self):
        """Stanley EC1 section 1.7, p. 65."""
        assert Y.from_rows((3, 3, 2, 1, 0, 0)) == Y.from_rows((3, 3, 2, 1))

    @given(diagrams())
    def test_boundary_sequence_determines_the_diagram(self, diagram):
        """Adin-Roichman Definition 9.6; the page's equivalence block."""
        sequence = diagram.boundary_sequence
        assert Y.from_boundary_sequence(sequence) == diagram
        assert sequence.count(1) == diagram.largest_part
        assert sequence.count(0) == diagram.length

    @pytest.mark.parametrize("n", range(11))
    def test_boundary_sequence_is_injective_on_partitions_of(self, n):
        sequences = {d.boundary_sequence for d in _partitions(n)}
        assert len(sequences) == len(_partitions(n))

    @given(diagrams())
    def test_frobenius_coordinates_are_strict_and_sum_to_the_size(self, diagram):
        """Stanley EC1 Exercise 1.70(a)."""
        arms, legs = diagram.frobenius_coordinates
        rank = diagram.durfee_side
        assert len(arms) == len(legs) == rank
        assert all(a > b for a, b in itertools.pairwise(arms))
        assert all(a > b for a, b in itertools.pairwise(legs))
        assert all(x >= 0 for x in (*arms, *legs))
        assert rank + sum(arms) + sum(legs) == diagram.size

    @given(diagrams())
    def test_frobenius_coordinates_round_trip(self, diagram):
        """Stanley EC1 Exercise 1.70(a); Frobenius 1900."""
        assert Y.from_frobenius(*diagram.frobenius_coordinates) == diagram

    @given(diagrams())
    def test_conjugation_swaps_the_frobenius_rows(self, diagram):
        arms, legs = diagram.frobenius_coordinates
        assert diagram.conjugate().frobenius_coordinates == (legs, arms)

    @pytest.mark.parametrize("n", range(9))
    def test_frobenius_coordinates_biject_with_all_strict_arrays(self, n):
        """Stanley EC1 Exercise 1.70(a); Frobenius 1900."""
        arrays = set()
        for rank in range(n + 1):
            pool = itertools.combinations(range(n, -1, -1), rank)
            for arms, legs in itertools.product(list(pool), repeat=2):
                if rank + sum(arms) + sum(legs) == n:
                    arrays.add((arms, legs))
        assert {d.frobenius_coordinates for d in _partitions(n)} == arrays
        assert len(arrays) == len(_partitions(n))


# --------------------------------------------------------------------------- #
# Conjugation
# --------------------------------------------------------------------------- #
class TestConjugation:
    @given(diagrams())
    def test_conjugation_is_a_size_preserving_involution(self, diagram):
        """Stanley EC1 section 1.8, p. 68 (unattributed there)."""
        assert diagram.conjugate().conjugate() == diagram
        assert diagram.conjugate().size == diagram.size

    @given(diagrams())
    def test_conjugation_transposes_the_cell_set(self, diagram):
        """Stanley EC1 section 1.8."""
        assert diagram.conjugate().cells == {(j, i) for i, j in diagram.cells}

    @given(diagrams())
    def test_conjugate_columns_count_the_rows_at_least_j(self, diagram):
        """Adin-Roichman section 2.4; Grinberg-Reiner Definition 2.2.8."""
        conjugate = diagram.conjugate()
        for j in range(1, diagram.largest_part + 1):
            assert conjugate.rows[j - 1] == sum(1 for part in diagram.rows if part >= j)

    @given(diagrams(), st.integers(0, 7))
    def test_at_most_k_parts_iff_conjugate_has_largest_part_at_most_k(self, diagram, k):
        """Stanley EC1 section 1.8, p. 68; the Ferrers-Sylvester (1853) argument."""
        assert (diagram.length <= k) == (diagram.conjugate().largest_part <= k)
        assert diagram.conjugate().length == diagram.largest_part

    @given(diagrams())
    def test_conjugate_multiplicities_are_the_row_differences(self, diagram):
        """Stanley EC1 section 1.8: m_i(lambda') = lambda_i - lambda_{i+1}."""
        padded = (*diagram.rows, 0)
        expected = {
            i: padded[i - 1] - padded[i]
            for i in range(1, diagram.length + 1)
            if padded[i - 1] > padded[i]
        }
        assert diagram.conjugate().multiplicities == expected

    def test_ferrers_proof_example_reads_the_columns_as_lines(self):
        """Sylvester 1853 (Ferrers's proof), as quoted in the page's overview."""
        assert Y((3, 3, 2, 1)).conjugate() == Y((4, 3, 2))

    def test_conjugation_shortcut_example(self):
        """Stanley EC1 section 1.8."""
        assert Y((4, 3, 1, 1, 1)).conjugate() == Y((5, 2, 2, 1))

    @given(diagrams())
    def test_row_column_double_count(self, diagram):
        """Stanley EC1 Proposition 1.8.3."""
        by_rows = sum((i - 1) * part for i, part in enumerate(diagram.rows, start=1))
        by_columns = sum(math.comb(part, 2) for part in diagram.conjugate().rows)
        assert by_rows == by_columns

    @given(diagrams())
    def test_self_conjugate_means_equal_to_the_conjugate(self, diagram):
        assert diagram.is_self_conjugate == (diagram.conjugate() == diagram)


# --------------------------------------------------------------------------- #
# Partition counts and identities proved on the diagram
# --------------------------------------------------------------------------- #
class TestPartitionCounts:
    @pytest.mark.parametrize(("n", "count"), list(enumerate(A000041)))
    def test_number_of_partitions(self, n, count):
        """Stanley EC1 eq. 1.77; OEIS A000041."""
        listed = _partitions(n)
        assert len(listed) == len(set(listed)) == count
        assert all(d.size == n for d in listed)

    @pytest.mark.parametrize("n", range(13))
    def test_at_most_k_parts_has_generating_function_prod_up_to_k(self, n):
        """Stanley EC1 eq. 1.76."""
        for k in range(n + 2):
            at_most_k_parts = sum(1 for d in _partitions(n) if d.length <= k)
            assert at_most_k_parts == _partitions_with_parts_at_most(n, k)

    @pytest.mark.parametrize("n", range(1, 14))
    def test_recurrence_for_exactly_k_parts(self, n):
        """Stanley EC1 section 1.7."""

        def p(k: int, m: int) -> int:
            return sum(1 for d in _partitions(m) if d.length == k) if m >= 0 else 0

        for k in range(1, n + 1):
            assert p(k, n) == p(k - 1, n - 1) + p(k, n - k)

    @pytest.mark.parametrize(("n", "count"), list(enumerate(A000700)))
    def test_number_of_self_conjugate_partitions(self, n, count):
        """Stanley EC1 Proposition 1.8.4; OEIS A000700."""
        assert sum(1 for d in _partitions(n) if d.is_self_conjugate) == count

    @pytest.mark.parametrize("n", range(17))
    def test_diagonal_hooks_biject_self_conjugate_with_distinct_odd_parts(self, n):
        """Stanley EC1 Proposition 1.8.4."""
        self_conjugate = [d for d in _partitions(n) if d.is_self_conjugate]
        distinct_odd = {
            d.rows
            for d in _partitions(n)
            if all(part % 2 for part in d.rows) and len(set(d.rows)) == d.length
        }
        images = [d.diagonal_hook_lengths for d in self_conjugate]
        assert len(set(images)) == len(images)
        assert set(images) == distinct_odd

    @pytest.mark.parametrize("n", range(17))
    def test_self_conjugate_partitions_by_durfee_side(self, n):
        """Stanley EC1 Exercise 1.76."""
        for k in range(n + 1):
            counted = sum(
                1 for d in _partitions(n) if d.is_self_conjugate and d.durfee_side == k
            )
            rest = n - k * k
            expected = (
                _partitions_with_parts_at_most(rest // 2, k)
                if rest >= 0 and rest % 2 == 0
                else 0
            )
            assert counted == expected

    @pytest.mark.parametrize("n", range(17))
    def test_self_conjugate_dissection_leaves_a_partition_and_its_conjugate(self, n):
        """Stanley EC1 Exercise 1.76."""
        self_conjugate = [d for d in _partitions(n) if d.is_self_conjugate]
        for diagram in self_conjugate:
            side, right, below = diagram.durfee_dissection()
            assert below == right.conjugate()
            assert below.largest_part <= side

    @given(st.integers(0, 5), diagrams(max_rows=5))
    def test_a_partition_and_its_conjugate_around_a_square_are_self_conjugate(
        self, side, right
    ):
        """Stanley EC1 Exercise 1.76."""
        right = right.meet(yd.rectangle(side, 6))
        diagram = Y.from_durfee_dissection(side, right, right.conjugate())
        assert diagram.is_self_conjugate
        assert diagram.durfee_side == side

    @pytest.mark.parametrize("n", range(15))
    def test_distinct_parts_equinumerous_with_odd_parts(self, n):
        """Euler's partition theorem; Stanley EC1 Proposition 1.8.5."""
        distinct = sum(1 for d in _partitions(n) if len(set(d.rows)) == d.length)
        odd = sum(1 for d in _partitions(n) if all(part % 2 for part in d.rows))
        assert distinct == odd

    @pytest.mark.parametrize("n", range(13))
    def test_companion_identity_a_by_number_of_parts(self, n):
        """Stanley EC1 Proposition 1.8.6(a), eqs. 1.82 and 1.84."""
        for k in range(n + 1):
            exactly_k = sum(1 for d in _partitions(n) if d.length == k)
            assert exactly_k == _partitions_with_parts_at_most(n - k, k)

    @pytest.mark.parametrize("n", range(13))
    def test_companion_identity_c_by_number_of_distinct_parts(self, n):
        """Stanley EC1 Proposition 1.8.6(c), eq. 1.83."""
        for k in range(n + 1):
            distinct_k = sum(
                1 for d in _partitions(n) if d.length == k and len(set(d.rows)) == k
            )
            rest = n - math.comb(k + 1, 2)
            expected = _partitions_with_parts_at_most(rest, k) if rest >= 0 else 0
            assert distinct_k == expected

    @pytest.mark.parametrize("n", range(31))
    def test_pentagonal_number_formula(self, n):
        """Stanley EC1 Proposition 1.8.7."""
        signed = sum(
            (-1) ** d.length for d in _partitions(n) if len(set(d.rows)) == d.length
        )
        pentagonal = {m * (3 * m - 1) // 2: (-1) ** m for m in range(-6, 7)}
        assert signed == pentagonal.get(n, 0)

    def test_pentagonal_expansion_as_printed(self):
        """Stanley EC1 Proposition 1.8.7."""
        nonzero = {
            n: sum(
                (-1) ** d.length for d in _partitions(n) if len(set(d.rows)) == d.length
            )
            for n in range(16)
        }
        assert {n: c for n, c in nonzero.items() if c} == {
            0: 1,
            1: -1,
            2: -1,
            5: 1,
            7: 1,
            12: -1,
            15: -1,
        }


# --------------------------------------------------------------------------- #
# Durfee dissection
# --------------------------------------------------------------------------- #
class TestDurfee:
    @given(diagrams())
    def test_durfee_side_is_the_largest_i_with_row_i_at_least_i(self, diagram):
        """Stanley EC1 section 1.8, pp. 71-72."""
        side = diagram.durfee_side
        assert all(diagram.rows[i - 1] >= i for i in range(1, side + 1))
        assert diagram.length <= side or diagram.rows[side] < side + 1
        assert side == sum(1 for i, j in diagram.cells if i == j)
        assert yd.rectangle(side, side).cells <= diagram.cells
        assert not yd.rectangle(side + 1, side + 1).cells <= diagram.cells

    @given(diagrams())
    def test_durfee_dissection_constraints_and_size(self, diagram):
        """Stanley EC1, proof of Proposition 1.8.6(b) (Durfee-square identity)."""
        side, right, below = diagram.durfee_dissection()
        assert right.length <= side
        assert below.largest_part <= side
        assert diagram.size == side * side + right.size + below.size
        assert diagram.length == side + below.length

    @given(diagrams())
    def test_durfee_dissection_round_trip(self, diagram):
        """Stanley EC1, proof of Proposition 1.8.6(b): obtained uniquely."""
        assert Y.from_durfee_dissection(*diagram.durfee_dissection()) == diagram

    @given(st.integers(0, 4), diagrams(max_rows=4), diagrams(max_width=4))
    def test_every_admissible_triple_gives_a_diagram_of_that_rank(
        self, side, right, below
    ):
        """Stanley EC1, proof of Proposition 1.8.6(b): obtained uniquely."""
        right = right.meet(yd.rectangle(side, 6))
        below = below.meet(yd.rectangle(6, side))
        diagram = Y.from_durfee_dissection(side, right, below)
        assert diagram.durfee_dissection() == (side, right, below)


# --------------------------------------------------------------------------- #
# Young's lattice
# --------------------------------------------------------------------------- #
class TestYoungsLattice:
    @given(diagrams(), diagrams())
    def test_containment_is_cell_containment_and_rowwise(self, a, b):
        """Stanley EC1 Example 3.4.4(b)."""
        rowwise = all(
            x <= y for x, y in itertools.zip_longest(a.rows, b.rows, fillvalue=0)
        )
        assert b.contains(a) == (a.cells <= b.cells) == rowwise

    @given(diagrams(), diagrams())
    def test_meet_is_intersection_and_join_is_union(self, a, b):
        """Stanley EC1 section 3.4."""
        assert a.meet(b).cells == a.cells & b.cells
        assert a.join(b).cells == a.cells | b.cells

    @given(diagrams(), diagrams())
    def test_meet_and_join_in_row_lengths(self, a, b):
        """Stanley EC1 section 3.4."""
        pairs = list(itertools.zip_longest(a.rows, b.rows, fillvalue=0))
        assert a.join(b) == Y.from_rows(max(x, y) for x, y in pairs)
        assert a.meet(b) == Y.from_rows(min(x, y) for x, y in pairs)

    @given(diagrams(), diagrams(), diagrams())
    def test_distributivity(self, a, b, c):
        """Stanley EC1 Example 3.4.4(b), Proposition 3.4.3."""
        assert a.meet(b.join(c)) == a.meet(b).join(a.meet(c))
        assert a.join(b.meet(c)) == a.join(b).meet(a.join(c))

    @given(diagrams(), diagrams())
    def test_absorption(self, a, b):
        """Stanley EC1 Example 3.4.4(b): the page's oracle."""
        assert a.meet(a.join(b)) == a
        assert a.join(a.meet(b)) == a

    @pytest.mark.parametrize("n", range(6))
    def test_every_cover_changes_the_size_by_exactly_one(self, n):
        """The page's oracle for the rank function |lambda| (the page's inference).

        Covers are recomputed from containment alone, among all diagrams with
        up to three more cells.
        """
        for lower in _partitions(n):
            above = [d for d in _up_to(n + 3) if d != lower and d.contains(lower)]
            covers = {
                d
                for d in above
                if not any(mid != d and d.contains(mid) for mid in above)
            }
            assert {d.size for d in covers} == {n + 1}
            assert covers == set(lower.upper_covers())
            assert all(lower in d.lower_covers() for d in covers)

    @given(diagrams())
    def test_covers_add_an_addable_cell_or_delete_a_corner(self, diagram):
        """Stanley EC1 section 3.4, applied to Y (the page's inference)."""
        assert [d.cells - diagram.cells for d in diagram.upper_covers()] == [
            {cell} for cell in diagram.addable_cells
        ]
        assert [diagram.cells - d.cells for d in diagram.lower_covers()] == [
            {cell} for cell in diagram.corners
        ]

    @given(diagrams())
    def test_corners_are_the_cells_of_hook_length_one(self, diagram):
        """Adin-Roichman section 6.1."""
        expected = {
            (i, j)
            for i, row in enumerate(diagram.hook_lengths, start=1)
            for j, h in enumerate(row, start=1)
            if h == 1
        }
        assert set(diagram.corners) == expected

    @given(diagrams())
    def test_addable_cells_and_corners_alternate_along_the_boundary(self, diagram):
        """Stanley EC1 Example 3.21.2(2)."""
        marked = [(cell, "x") for cell in diagram.addable_cells]
        marked += [(cell, "o") for cell in diagram.corners]
        # Boundary order from the NE end: down the rows, right to left in a row.
        marked.sort(key=lambda item: (item[0][0], -item[0][1]))
        kinds = "".join(kind for _, kind in marked)
        assert kinds == "xo" * len(diagram.corners) + "x"

    @given(diagrams())
    def test_youngs_lattice_is_one_differential(self, diagram):
        """Stanley EC1 Example 3.21.2(2); differential posets: Fomin, Stanley."""
        assert len(diagram.addable_cells) - len(diagram.corners) == 1
        assert len(diagram.upper_covers()) == len(diagram.lower_covers()) + 1

    @given(diagrams())
    def test_down_up_minus_up_down_is_the_identity(self, diagram):
        """Stanley EC1 Proposition 3.21.3."""
        down_up = collections.Counter(
            low for up in diagram.upper_covers() for low in up.lower_covers()
        )
        up_down = collections.Counter(
            up for low in diagram.lower_covers() for up in low.upper_covers()
        )
        down_up.subtract(up_down)
        assert {d: c for d, c in down_up.items() if c} == {diagram: 1}

    @pytest.mark.parametrize(("n", "count"), list(enumerate(COVER_COUNTS)))
    def test_number_of_covers(self, n, count):
        """Stanley EC1 Exercise 1.71."""
        pairs = sum(len(d.upper_covers()) for d in _partitions(n))
        assert pairs == sum(len(_partitions(m)) for m in range(n + 1)) == count

    @pytest.mark.parametrize(("lower", "upper"), NESTED_SMALL)
    def test_mobius_function_of_intervals(self, lower, upper):
        """Stanley EC1 Example 3.9.6, specialized to Y."""
        assert yd.mobius(lower, upper) == _mobius_by_recursion(lower, upper)

    @pytest.mark.parametrize(("lower", "upper"), NESTED_SMALL)
    def test_multichain_determinant_counts_z2_and_z3(self, lower, upper):
        """Stanley EC1 Exercise 3.149; Kreweras 1965."""
        interval = _interval(lower, upper)
        pairs = sum(1 for a in interval for b in interval if b.contains(a))
        assert yd.count_multichains(lower, upper, 1) == len(interval)
        assert yd.count_multichains(lower, upper, 2) == pairs

    @pytest.mark.parametrize("n", range(11))
    def test_chain_counts_square_to_n_factorial_and_sum_to_involutions(self, n):
        """Stanley EC1 Theorems 3.21.8 and 3.21.10."""
        chains = [
            d.count_standard_tableaux(yd.TableauFormula.CHAINS) for d in _partitions(n)
        ]
        assert sum(e * e for e in chains) == math.factorial(n)
        assert sum(chains) == A000085[n]

    @pytest.mark.parametrize("n", range(8))
    def test_chain_count_is_the_number_of_linear_extensions(self, n):
        """Stanley EC1 section 3.5, eq. 3.11; Adin-Roichman section 2.5.1.

        Brute force: fillings of the cells by 1..n increasing along rows and
        down columns, independent of the corner-removal recursion.
        """
        for diagram in _partitions(n):
            cells = sorted(diagram.cells)
            extensions = set()
            for entries in itertools.permutations(range(1, n + 1)):
                filling = dict(zip(cells, entries, strict=True))
                increasing = all(
                    filling[i, j] > filling.get((i - 1, j), 0)
                    and filling[i, j] > filling.get((i, j - 1), 0)
                    for i, j in cells
                )
                if increasing:
                    extensions.add(
                        tuple(
                            tuple(filling[i, j] for j in range(1, length + 1))
                            for i, length in enumerate(diagram.rows, start=1)
                        )
                    )
            chains = diagram.count_standard_tableaux(yd.TableauFormula.CHAINS)
            assert len(extensions) == chains
            assert set(diagram.standard_tableaux()) == extensions

    @pytest.mark.parametrize("n", range(6))
    def test_closed_hasse_walks_from_the_bottom(self, n):
        """Stanley EC1 Theorem 3.21.7."""
        walks = collections.Counter({yd.empty_diagram(): 1})
        for _ in range(2 * n):
            step: collections.Counter[yd.YoungDiagram] = collections.Counter()
            for diagram, count in walks.items():
                for neighbor in (*diagram.upper_covers(), *diagram.lower_covers()):
                    step[neighbor] += count
            walks = step
        double_factorial = math.prod(range(2 * n - 1, 0, -2))
        assert walks[yd.empty_diagram()] == double_factorial

    @given(nested_pairs())
    def test_skew_cells_are_the_set_difference(self, pair):
        """Adin-Roichman Definition 2.12."""
        inner, outer = pair
        expected = {
            (i, j)
            for i, j in outer.cells
            if inner.padded_rows(outer.length)[i - 1] + 1 <= j
        }
        assert outer.skew_cells(inner) == outer.cells - inner.cells == expected

    def test_skew_shapes_are_exactly_the_order_convex_diagrams(self):
        """Adin-Roichman Observation 2.13."""
        grid = [(i, j) for i in range(1, 4) for j in range(1, 4)]
        box = list(yd.diagrams_in_box(3, 6))
        skew = {
            outer.skew_cells(inner)
            for outer in box
            for inner in box
            if outer.contains(inner)
        }
        for size in range(len(grid) + 1):
            for subset in itertools.combinations(grid, size):
                cells = frozenset(subset)
                assert (cells in skew) == _is_order_convex(cells)

    @given(diagrams().filter(lambda d: d.size > 0))
    def test_an_ordinary_diagram_is_path_connected_and_line_convex(self, diagram):
        """Adin-Roichman Observation 2.13 (a necessary condition only)."""
        assert _components(diagram.cells, by_order=False) == {diagram.cells}
        assert _is_line_convex(diagram.cells)

    def test_path_connected_and_line_convex_does_not_characterize_ordinary(self):
        """The page's gloss on Adin-Roichman Observation 2.13; the suite's witness."""
        cells = Y((2, 2)).skew_cells(Y((1,)))
        assert cells == {(1, 2), (2, 1), (2, 2)}
        assert _components(cells, by_order=False) == {cells}
        assert _is_line_convex(cells)
        assert not _is_down_closed(cells)

    @given(nested_pairs())
    def test_a_skew_diagram_is_line_convex(self, pair):
        """Adin-Roichman, as quoted in the page's Skew Young Diagram block."""
        inner, outer = pair
        assert _is_line_convex(outer.skew_cells(inner))

    @given(nested_pairs())
    def test_skew_path_components_are_the_order_components(self, pair):
        """Adin-Roichman, as quoted in the page's Skew Young Diagram block."""
        inner, outer = pair
        cells = outer.skew_cells(inner)
        assert _components(cells, by_order=False) == _components(cells, by_order=True)

    def test_a_skew_diagram_need_not_be_path_connected(self):
        """Adin-Roichman: "not necessarily path-connected" (the page's skew block)."""
        cells = yd.skew_non_example_cells()
        assert len(_components(cells, by_order=False)) > 1


# --------------------------------------------------------------------------- #
# Diagrams in a box
# --------------------------------------------------------------------------- #
class TestBox:
    @pytest.mark.parametrize(("k", "n"), BOXES)
    def test_box_count_and_size_distribution(self, k, n):
        """Stanley EC1 Proposition 1.7.3."""
        box = list(yd.diagrams_in_box(k, n))
        assert len(box) == len(set(box)) == math.comb(n, k)
        assert all(d.fits_in_box(k, n) for d in box)
        assert yd.box_size_distribution(k, n) == _q_binomial(n, k)

    @pytest.mark.parametrize(("k", "n"), SMALL_BOXES)
    def test_box_diagrams_are_the_partitions_with_bounded_parts(self, k, n):
        """Stanley EC1 Proposition 1.7.3."""
        expected = {
            d for d in _up_to(k * (n - k)) if d.length <= k and d.largest_part <= n - k
        }
        assert set(yd.diagrams_in_box(k, n)) == expected

    @pytest.mark.parametrize(("k", "n"), SMALL_BOXES)
    def test_bijection_with_k_subsets(self, k, n):
        """Postnikov section 2.1."""
        subsets = [d.to_subset(k, n) for d in yd.diagrams_in_box(k, n)]
        expected = {frozenset(c) for c in itertools.combinations(range(1, n + 1), k)}
        assert len(set(subsets)) == len(subsets)
        assert set(subsets) == expected

    @given(boxed_diagrams())
    def test_subset_round_trip_through_the_inverse_formula(self, boxed):
        """Postnikov section 2.1."""
        diagram, k, n = boxed
        subset = diagram.to_subset(k, n)
        literal = [len(set(range(i, n + 1)) - subset) for i in sorted(subset)]
        assert Y.from_rows(literal) == diagram
        assert Y.from_subset(subset, n) == diagram

    @given(boxed_diagrams())
    def test_subset_is_the_set_of_vertical_step_labels(self, boxed):
        """Postnikov section 2.1; Stanley EC1 p. 67."""
        diagram, k, n = boxed
        # Walk the border from the upper-right to the lower-left corner.
        labels, label, column = set(), 0, n - k
        for length in diagram.padded_rows(k):
            label += column - length + 1
            column = length
            labels.add(label)
        assert diagram.to_subset(k, n) == labels
        assert sorted(labels) == [
            (n - k) - length + s
            for s, length in enumerate(diagram.padded_rows(k), start=1)
        ]

    @given(boxed_diagrams())
    def test_subset_agrees_with_the_le_diagram_module(self, boxed):
        diagram, k, n = boxed
        assert diagram.to_le_diagram(k, n).vertical_steps == diagram.to_subset(k, n)

    @pytest.mark.parametrize(("k", "n"), SMALL_BOXES)
    def test_containment_reverses_the_componentwise_order_on_subsets(self, k, n):
        """Postnikov, proof of Corollary 17.7."""
        box = list(yd.diagrams_in_box(k, n))
        for small, large in itertools.product(box, repeat=2):
            i_small, i_large = small.to_subset(k, n), large.to_subset(k, n)
            prefix = all(
                len(i_small & set(range(1, b + 1)))
                <= len(i_large & set(range(1, b + 1)))
                for b in range(1, n + 1)
            )
            gale = all(
                x <= y for x, y in zip(sorted(i_large), sorted(i_small), strict=True)
            )
            assert large.contains(small) == prefix == gale

    @pytest.mark.parametrize(("k", "n"), SMALL_BOXES)
    def test_extreme_subsets(self, k, n):
        """Postnikov section 2.1; the page's Corollary 17.7 block."""
        assert yd.rectangle(k, n - k).to_subset(k, n) == frozenset(range(1, k + 1))
        assert yd.empty_diagram().to_subset(k, n) == frozenset(range(n - k + 1, n + 1))

    @given(boxed_diagrams())
    def test_boundary_word_has_size_many_inversions(self, boxed):
        """Stanley EC1 section 1.7, p. 67; Postnikov section 19."""
        diagram, k, n = boxed
        word = diagram.boundary_word(k, n)
        inversions = sum(1 for a, b in itertools.combinations(word, 2) if a > b)
        assert word.count(2) == k
        assert word.count(1) == n - k
        assert {p for p, step in enumerate(word, start=1) if step == 2} == (
            diagram.to_subset(k, n)
        )
        assert inversions == diagram.size

    @given(boxed_diagrams())
    def test_le_diagram_shape_round_trip(self, boxed):
        diagram, k, n = boxed
        le = diagram.to_le_diagram(k, n)
        assert le.diagram_type == (k, n)
        assert le.rank == 0
        assert Y.from_le_diagram(le) == diagram

    @pytest.mark.parametrize(("k", "n"), [(k, n) for k, n in BOXES if n <= 5])
    def test_every_le_diagram_of_a_shape_has_that_shape(self, k, n):
        """API consistency with ``enumerate_le_diagrams_of_shape``; no theorem."""
        for diagram in yd.diagrams_in_box(k, n):
            fillings = ld.enumerate_le_diagrams_of_shape(diagram.padded_rows(k), n)
            assert {Y.from_le_diagram(le) for le in fillings} == {diagram}

    @pytest.mark.parametrize(("k", "n"), [(k, n) for k, n in BOXES if 1 <= n <= 5])
    def test_subset_is_the_lexicographically_minimal_base_of_the_cells(self, k, n):
        """Postnikov section 2.3 and Lemma 17.3 (one-sided), via Theorem 6.5."""
        for diagram in yd.diagrams_in_box(k, n):
            expected = sorted(diagram.to_subset(k, n))
            for le in ld.enumerate_le_diagrams_of_shape(diagram.padded_rows(k), n):
                bases = le.to_positroid().bases
                assert min(sorted(base) for base in bases) == expected
                assert all(
                    Y.from_subset(base, n).cells <= diagram.cells for base in bases
                )

    def test_top_cell_filling_is_a_le_diagram_on_the_rectangle(self):
        shape = yd.rectangle(2, 3)
        assert shape.to_le_diagram(2, 5, shape.cells) == ld.top_cell_diagram(2, 5)


# --------------------------------------------------------------------------- #
# Hooks and tableaux
# --------------------------------------------------------------------------- #
class TestHooks:
    @given(diagrams_with_cell())
    def test_hook_is_the_cell_its_arm_and_its_leg(self, case):
        """Frame-Robinson-Thrall 1954, eq. 1.1; Adin-Roichman Definition 5.2."""
        diagram, i, j = case
        expected = {
            (s, t)
            for s, t in diagram.cells
            if (s == i and t >= j) or (t == j and s >= i)
        }
        assert diagram.hook(i, j) == expected
        assert diagram.arm_length(i, j) == diagram.rows[i - 1] - j
        assert diagram.leg_length(i, j) == diagram.conjugate().rows[j - 1] - i
        assert (
            diagram.hook_length(i, j)
            == len(expected)
            == diagram.rows[i - 1] + diagram.conjugate().rows[j - 1] - i - j + 1
            == diagram.hook_lengths[i - 1][j - 1]
        )

    @given(diagrams())
    def test_hook_product_is_the_product_of_the_hook_graph(self, diagram):
        assert diagram.hook_product == math.prod(
            diagram.hook_length(i, j) for i, j in diagram.cells
        )

    @pytest.mark.parametrize("n", range(11))
    def test_hook_length_formula_matches_corner_removal(self, n):
        """Frame-Robinson-Thrall 1954, Theorem 1; Adin-Roichman Theorem 5.3."""
        for diagram in _partitions(n):
            assert math.factorial(n) % diagram.hook_product == 0
            assert diagram.count_standard_tableaux(
                yd.TableauFormula.HOOK
            ) == diagram.count_standard_tableaux(yd.TableauFormula.CHAINS)

    @pytest.mark.parametrize("n", range(8))
    def test_tableau_count_is_the_number_of_standard_tableaux(self, n):
        """Adin-Roichman Definition 2.3 and section 2.5.1."""
        for diagram in _partitions(n):
            tableaux = list(diagram.standard_tableaux())
            assert len(set(tableaux)) == len(tableaux)
            assert len(tableaux) == diagram.count_standard_tableaux()

    @pytest.mark.parametrize("n", range(7))
    def test_standard_tableaux_are_order_preserving_bijections(self, n):
        """Adin-Roichman Definition 2.3."""
        for diagram in _partitions(n):
            for tableau in diagram.standard_tableaux():
                assert tuple(len(row) for row in tableau) == diagram.rows
                assert sorted(e for row in tableau for e in row) == list(
                    range(1, n + 1)
                )
                assert all(a < b for row in tableau for a, b in itertools.pairwise(row))
                assert all(
                    upper[j] < entry
                    for upper, lower in itertools.pairwise(tableau)
                    for j, entry in enumerate(lower)
                )

    @given(diagrams(max_rows=5, max_width=5))
    def test_the_closed_formulas_and_the_chain_recursion_agree(self, diagram):
        """Adin-Roichman Theorems 5.1, 5.3, 5.4 (Claim 5.5) and section 2.5.1."""
        counts = {diagram.count_standard_tableaux(f) for f in yd.TableauFormula}
        assert len(counts) == 1

    @given(diagrams())
    def test_tableau_count_is_transpose_invariant(self, diagram):
        """Adin-Roichman Observation 2.8."""
        assert (
            diagram.count_standard_tableaux()
            == diagram.conjugate().count_standard_tableaux()
        )

    @given(diagrams_with_cell())
    def test_hook_permutation_lemma(self, case):
        """Frame-Robinson-Thrall 1954, Lemma 1."""
        diagram, i, j = case
        h = diagram.hook_length(i, j)
        along_arm = [
            diagram.hook_length(i, t) for t in range(j, diagram.rows[i - 1] + 1)
        ]
        along_leg = [
            h - diagram.hook_length(s, j)
            for s in range(i + 1, diagram.conjugate().rows[j - 1] + 1)
        ]
        assert sorted(along_arm + along_leg) == list(range(1, h + 1))

    @given(diagrams_with_cell())
    def test_leg_length_from_the_hook_graph(self, case):
        """Frame-Robinson-Thrall 1954, Theorem 7."""
        diagram, i, j = case
        h = diagram.hook_length(i, j)
        to_the_right = set(diagram.hook_lengths[i - 1][j:])
        missing = [m for m in range(1, h) if m not in to_the_right]
        assert diagram.leg_length(i, j) == len(missing)

    @pytest.mark.parametrize("n", range(8))
    def test_major_index_generating_function(self, n):
        """Adin-Roichman Theorem 10.26 (Stanley EC2 Corollary 7.21.5)."""
        for diagram in _partitions(n):
            counts = collections.Counter(
                _major_index(t) for t in diagram.standard_tableaux()
            )
            polynomial = diagram.major_index_polynomial()
            assert {p: c for p, c in enumerate(polynomial) if c} == dict(counts)

    @pytest.mark.parametrize("n", range(1, 10))
    def test_hook_shape_closed_form(self, n):
        """Adin-Roichman Observation 3.1."""
        for k in range(n):
            assert yd.hook_shape(n, k).count_standard_tableaux() == math.comb(n - 1, k)

    @pytest.mark.parametrize("n", range(1, 12))
    def test_two_rowed_closed_form(self, n):
        """Adin-Roichman Proposition 3.3."""
        for k in range(n // 2 + 1):
            expected = math.comb(n, k) - (math.comb(n, k - 1) if k else 0)
            assert Y.from_rows((n - k, k)).count_standard_tableaux() == expected

    @pytest.mark.parametrize(("m", "catalan"), list(enumerate(CATALAN, start=1)))
    def test_two_rowed_catalan_cases(self, m, catalan):
        """Adin-Roichman Proposition 3.3; OEIS A000108."""
        assert Y((m, m)).count_standard_tableaux() == catalan
        assert Y.from_rows((m, m - 1)).count_standard_tableaux() == catalan

    @pytest.mark.parametrize("n", range(11))
    def test_robinson_schensted_identities(self, n):
        """Adin-Roichman Theorem 4.9, Corollary 4.12; OEIS A000085."""
        counts = [d.count_standard_tableaux() for d in _partitions(n)]
        assert sum(f * f for f in counts) == math.factorial(n)
        assert sum(counts) == A000085[n]


class TestRimHooks:
    @given(diagrams())
    def test_rim_nodes_have_no_cell_to_their_south_east(self, diagram):
        """Frame-Robinson-Thrall 1954, section 2."""
        assert diagram.rim == {
            (s, t) for s, t in diagram.cells if (s + 1, t + 1) not in diagram.cells
        }

    @given(diagrams_with_cell())
    def test_rim_hook_has_hook_length_many_rim_nodes_from_head_to_foot(self, case):
        """Frame-Robinson-Thrall 1954, section 2."""
        diagram, i, j = case
        rim_hook = diagram.rim_hook(i, j)
        head = (i, diagram.rows[i - 1])
        foot = (diagram.conjugate().rows[j - 1], j)
        assert len(rim_hook) == diagram.hook_length(i, j)
        assert rim_hook <= diagram.rim
        assert head in rim_hook
        assert foot in rim_hook
        assert all(
            head[0] <= s <= foot[0] and foot[1] <= t <= head[1] for s, t in rim_hook
        )

    @given(diagrams_with_cell())
    def test_removing_a_rim_hook_leaves_a_diagram_with_hook_length_fewer_cells(
        self, case
    ):
        """Frame-Robinson-Thrall 1954, section 3 (the page's rim-hook removal block)."""
        diagram, i, j = case
        smaller = diagram.remove_rim_hook(i, j)
        assert smaller.size == diagram.size - diagram.hook_length(i, j)
        assert smaller.cells == diagram.cells - diagram.rim_hook(i, j)

    @given(diagrams_with_cell())
    def test_removing_the_rim_hook_matches_removing_the_right_hook(self, case):
        """Frame-Robinson-Thrall section 3 treat the two removals as the same diagram.

        The page does not define right-hook removal; "delete the hook and
        close up along the diagonal" is this suite's reading of it.
        """
        diagram, i, j = case
        # Delete the hook, then close up: cells south-east of it move up-left.
        hook = diagram.hook(i, j)
        closed = {
            (s - 1, t - 1) if s > i and t > j else (s, t)
            for s, t in diagram.cells - hook
        }
        assert diagram.remove_rim_hook(i, j).cells == closed

    @given(diagrams_with_cell().filter(lambda case: case[2] > 1))
    def test_q_hook_off_the_first_column_diminishes_one_first_column_hook(self, case):
        """Frame-Robinson-Thrall 1954, proof of Theorem 6."""
        diagram, i, j = case
        q = diagram.hook_length(i, j)
        before = [row[0] for row in diagram.hook_lengths]
        after = [row[0] for row in diagram.remove_rim_hook(i, j).hook_lengths]
        candidates = [
            sorted([*before[:r], before[r] - q, *before[r + 1 :]])
            for r in range(len(before))
        ]
        assert sorted(after) in candidates

    @given(diagrams(), st.integers(1, 5), st.randoms(use_true_random=False))
    def test_core_is_what_remains_after_all_q_hooks_in_any_order(
        self, diagram, q, random
    ):
        """Frame-Robinson-Thrall 1954, section 6."""
        current = diagram
        while True:
            q_hooks = [
                (i, j)
                for i, row in enumerate(current.hook_lengths, start=1)
                for j, h in enumerate(row, start=1)
                if h == q
            ]
            if not q_hooks:
                break
            current = current.remove_rim_hook(*random.choice(q_hooks))
        assert current == diagram.core(q)
        assert (diagram.size - current.size) % q == 0

    @given(diagrams(), st.integers(1, 5))
    def test_core_from_the_first_column_hook_lengths(self, diagram, q):
        """Frame-Robinson-Thrall 1954, Theorem 6."""
        rows = diagram.length
        alpha = diagram.core(q).padded_rows(rows)
        reduced = [
            row[0] - q * sum(1 for h in row if h % q == 0)
            for row in diagram.hook_lengths
        ]
        assert len(set(reduced)) == rows
        assert all(value >= 0 for value in reduced)
        assert sorted(reduced) == sorted(
            alpha[i - 1] - i + rows for i in range(1, rows + 1)
        )


# --------------------------------------------------------------------------- #
# Dominance order
# --------------------------------------------------------------------------- #
class TestDominance:
    @given(same_size_pairs())
    def test_dominance_is_the_partial_sum_order(self, pair):
        """Grinberg-Reiner Definition 2.2.7."""
        a, b = pair
        width = max(a.length, b.length)
        expected = all(
            sum(a.padded_rows(width)[:k]) >= sum(b.padded_rows(width)[:k])
            for k in range(1, width + 1)
        )
        assert a.dominates(b) == expected

    @pytest.mark.parametrize("n", range(9))
    def test_dominance_is_a_partial_order(self, n):
        listed = _partitions(n)
        assert all(a.dominates(a) for a in listed)
        for a, b in itertools.combinations(listed, 2):
            assert not (a.dominates(b) and b.dominates(a))
        for a, b, c in itertools.product(listed, repeat=3):
            if a.dominates(b) and b.dominates(c):
                assert a.dominates(c)

    @pytest.mark.parametrize("n", range(9))
    def test_meet_and_join_are_greatest_lower_and_least_upper_bounds(self, n):
        """Brylawski 1973 (Proposition 2.2 per Behrisch et al.)."""
        listed = _partitions(n)
        for a, b in itertools.product(listed, repeat=2):
            meet, join = a.dominance_meet(b), a.dominance_join(b)
            lower = [c for c in listed if a.dominates(c) and b.dominates(c)]
            upper = [c for c in listed if c.dominates(a) and c.dominates(b)]
            assert meet in lower
            assert all(meet.dominates(c) for c in lower)
            assert join in upper
            assert all(c.dominates(join) for c in upper)

    @given(same_size_pairs())
    def test_meet_has_the_minima_of_the_partial_sums(self, pair):
        """Brylawski 1973, per Latapy-Phan, proof of Proposition 1."""
        a, b = pair
        width = max(a.length, b.length)
        meet = a.dominance_meet(b).padded_rows(width)
        for k in range(1, width + 1):
            assert sum(meet[:k]) == min(
                sum(a.padded_rows(width)[:k]), sum(b.padded_rows(width)[:k])
            )

    @given(same_size_pairs())
    def test_join_is_the_conjugate_of_the_meet_of_the_conjugates(self, pair):
        """Wikipedia, 'Dominance order'; verified on the page for n <= 9."""
        a, b = pair
        assert a.dominance_join(b) == (
            a.conjugate().dominance_meet(b.conjugate()).conjugate()
        )

    def test_join_example_where_partial_sum_maxima_fail(self):
        """Wikipedia, 'Dominance order'."""
        a, b = Y((3, 1, 1, 1)), Y((2, 2, 2))
        assert a.dominance_join(b) == Y((3, 2, 1))
        maxima = [
            max(x, y)
            for x, y in zip(
                itertools.accumulate(a.padded_rows(4)),
                itertools.accumulate(b.padded_rows(4)),
                strict=True,
            )
        ]
        parts = [y - x for x, y in itertools.pairwise([0, *maxima])]
        assert parts == [3, 1, 2, 0]

    @pytest.mark.parametrize("n", range(11))
    def test_conjugation_is_an_anti_automorphism(self, n):
        """Grinberg-Reiner Exercise 2.2.9; Brylawski 1973, abstract."""
        for a, b in itertools.product(_partitions(n), repeat=2):
            assert a.dominates(b) == b.conjugate().dominates(a.conjugate())

    @pytest.mark.parametrize("n", range(10))
    def test_covers_are_single_box_moves(self, n):
        """Brylawski 1973, as restated in Behrisch et al. Theorem 3."""
        listed = _partitions(n)
        for a in listed:
            below = [b for b in listed if b != a and a.dominates(b)]
            covered = {
                b for b in below if not any(c != b and c.dominates(b) for c in below)
            }
            assert set(a.dominance_lower_covers()) == covered

    @pytest.mark.parametrize("n", range(10))
    def test_blass_sagan_description_of_the_covers(self, n):
        """Blass-Sagan section 5: adjacent rows not in a wall, or the ends of a wall."""
        listed = _partitions(n)
        for beta in listed:
            rows = beta.rows
            walls = [
                (i, j)
                for i, j in itertools.combinations(range(len(rows)), 2)
                if rows[i] == rows[j]
                and (i == 0 or rows[i - 1] > rows[i])
                and (j == len(rows) - 1 or rows[j + 1] < rows[j])
            ]
            in_wall = {r for i, j in walls for r in range(i, j + 1)}
            raised = [(i, j) for i, j in walls]
            raised += [
                (i, i + 1)
                for i in range(len(rows) - 1)
                if i not in in_wall and i + 1 not in in_wall
            ]
            described = {
                Y.from_rows(part + (r == i) - (r == j) for r, part in enumerate(rows))
                for i, j in raised
            }
            covering = {a for a in listed if beta in a.dominance_lower_covers()}
            assert described == covering

    @given(diagrams())
    def test_partial_sum_tuples_are_nondecreasing_and_concave(self, diagram):
        """Wikipedia's associated (n+1)-tuple; Brylawski's, the page presumes."""
        n = diagram.size
        sums = [0, *itertools.accumulate(diagram.padded_rows(n))]
        assert len(sums) == n + 1
        assert all(a <= b for a, b in itertools.pairwise(sums))
        assert all(
            2 * b >= a + c for a, b, c in zip(sums, sums[1:], sums[2:], strict=False)
        )

    @pytest.mark.parametrize("n", range(1, 9))
    def test_dominance_is_the_closure_of_the_transition_rule_from_n(self, n):
        """Latapy-Phan section 1."""

        @functools.cache
        def reachable(start: yd.YoungDiagram) -> frozenset[yd.YoungDiagram]:
            below = {start}
            for covered in start.dominance_lower_covers():
                below |= reachable(covered)
            return frozenset(below)

        assert reachable(Y((n,))) == set(_partitions(n))
        for a in _partitions(n):
            assert reachable(a) == {b for b in _partitions(n) if a.dominates(b)}

    @pytest.mark.parametrize("n", range(10))
    def test_dominance_is_a_chain_iff_n_at_most_5(self, n):
        """The page's small-n thresholds (Wikipedia; computed on the page)."""
        total = all(
            a.dominates(b) or b.dominates(a)
            for a, b in itertools.combinations(_partitions(n), 2)
        )
        assert total == (n <= 5)

    @pytest.mark.parametrize("n", range(1, 11))
    def test_dominance_is_graded_iff_n_at_most_6(self, n):
        """The page's small-n thresholds (Wikipedia; Early; computed for n <= 10)."""
        shortest, longest = _chain_lengths(n)
        assert (shortest == longest) == (n <= 6)

    @pytest.mark.parametrize("n", range(9))
    def test_dominance_is_distributive_and_modular_iff_n_at_most_6(self, n):
        """The page's small-n thresholds (modularity is the page's own computation)."""
        listed = _partitions(n)
        distributive = all(
            a.dominance_meet(b.dominance_join(c))
            == a.dominance_meet(b).dominance_join(a.dominance_meet(c))
            for a, b, c in itertools.product(listed, repeat=3)
        )
        modular = all(
            a.dominance_join(b.dominance_meet(c))
            == a.dominance_join(b).dominance_meet(c)
            for a, b, c in itertools.product(listed, repeat=3)
            if c.dominates(a)
        )
        assert distributive == (n <= 6)
        assert modular == (n <= 6)

    @pytest.mark.parametrize(
        ("n", "longest"), [(6, 8), (7, 11), (8, 14), (9, 17), (10, 20)]
    )
    def test_maximal_chain_lengths(self, n, longest):
        """Greene-Kleitman 1986 as quoted by Early; 2n - 4 computed on the page."""
        a = max(a for a in range(1, n + 1) if math.comb(a + 1, 2) <= n)
        c = n - math.comb(a + 1, 2)
        assert 0 <= c <= a
        assert _chain_lengths(n) == (2 * n - 4, longest)
        assert longest == (a**3 - a) // 3 + c * a

    @pytest.mark.parametrize("n", range(9))
    def test_dominance_mobius_function_takes_values_0_1_minus_1(self, n):
        """Brylawski 1973, abstract."""
        listed = _partitions(n)

        @functools.cache
        def mobius(low: yd.YoungDiagram, high: yd.YoungDiagram) -> int:
            if low == high:
                return 1
            return -sum(
                mobius(low, mid)
                for mid in listed
                if mid != high and high.dominates(mid) and mid.dominates(low)
            )

        values = {
            mobius(low, high)
            for low, high in itertools.product(listed, repeat=2)
            if high.dominates(low)
        }
        assert values <= {-1, 0, 1}

    def test_dominance_between_different_sizes_is_rejected(self):
        with pytest.raises(ValueError, match="same n"):
            Y((2, 1)).dominates(Y((2,)))
        with pytest.raises(ValueError, match="same n"):
            Y((2, 1)).dominance_meet(Y((2,)))
        with pytest.raises(ValueError, match="same n"):
            Y((2, 1)).dominance_join(Y((2,)))


# --------------------------------------------------------------------------- #
# Round-trip laws
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(diagrams())
    def test_dataframe_round_trip(self, diagram):
        assert Y.from_dataframe(diagram.to_dataframe()) == diagram

    @given(diagrams(), st.randoms(use_true_random=False))
    def test_dataframe_rows_may_come_in_any_order(self, diagram, random):
        frame = diagram.to_dataframe()
        order = list(range(len(frame)))
        random.shuffle(order)
        assert Y.from_dataframe(frame.iloc[order]) == diagram

    @given(diagrams())
    def test_dataframe_is_tidy_one_row_per_cell(self, diagram):
        frame = diagram.to_dataframe()
        assert list(frame.columns) == ["row", "column"]
        assert len(frame) == diagram.size
        assert set(zip(frame["row"], frame["column"], strict=True)) == diagram.cells

    @pytest.mark.parametrize(
        "example",
        [
            yd.all_in_one_diagram(),
            yd.postnikov_figure_2_1_diagram(),
            yd.empty_diagram(),
        ],
    )
    def test_experiment_io_round_trip(self, example, tmp_path):
        path = io.write_result(example.to_dataframe(), tmp_path / "young.json")
        assert Y.from_dataframe(pd.read_json(path, dtype=False)) == example

    @given(diagrams())
    def test_constructors_agree(self, diagram):
        k = diagram.length
        n = k + diagram.largest_part
        assert Y(diagram.rows) == diagram
        assert Y.from_rows((*diagram.rows, 0, 0)) == diagram
        assert Y.from_cells(diagram.cells) == diagram
        assert Y.from_boundary_sequence(diagram.boundary_sequence) == diagram
        assert Y.from_frobenius(*diagram.frobenius_coordinates) == diagram
        assert Y.from_subset(diagram.to_subset(k, n), n) == diagram
        assert Y.from_durfee_dissection(*diagram.durfee_dissection()) == diagram
        assert Y.from_le_diagram(diagram.to_le_diagram(k, n)) == diagram

    @given(st.lists(st.integers(0, 1), max_size=12))
    def test_untrimmed_boundary_sequences_count_the_ones_before_each_zero(self, steps):
        expected = [sum(steps[:p]) for p, step in enumerate(steps) if step == 0]
        assert Y.from_boundary_sequence(steps) == Y.from_rows(reversed(expected))

    @given(diagrams())
    def test_add_then_remove_a_cell(self, diagram):
        for cell in diagram.addable_cells:
            assert diagram.add_cell(cell).remove_cell(cell) == diagram
        for cell in diagram.corners:
            assert diagram.remove_cell(cell).add_cell(cell) == diagram

    @given(diagrams())
    def test_diagrams_are_hashable_values(self, diagram):
        assert hash(Y.from_rows((*diagram.rows, 0))) == hash(diagram)
        assert repr(diagram) == f"YoungDiagram(rows={diagram.rows!r})"


# --------------------------------------------------------------------------- #
# Rejections
# --------------------------------------------------------------------------- #
class TestRejections:
    def test_rows_that_increase_are_rejected(self):
        with pytest.raises(ValueError, match="weak decrease violated"):
            Y((2, 3))
        with pytest.raises(ValueError, match="weak decrease violated"):
            Y.from_rows((1, 0, 1))

    def test_negative_parts_are_rejected(self):
        with pytest.raises(ValueError, match="positivity violated"):
            Y.from_rows((2, -1))

    def test_stored_zero_rows_are_rejected(self):
        with pytest.raises(ValueError, match="positivity violated"):
            Y((2, 0))

    def test_cells_outside_the_quarter_plane_are_rejected(self):
        with pytest.raises(ValueError, match="quarter-plane violated"):
            Y.from_cells({(0, 1)})
        with pytest.raises(ValueError, match="quarter-plane violated"):
            Y.from_cells({(1, 1), (1, 0)})

    @pytest.mark.parametrize(
        "cells", [{(1, 2)}, {(2, 1)}, {(1, 1), (2, 2)}, {(1, 1), (1, 3)}]
    )
    def test_cells_that_are_not_down_closed_are_rejected(self, cells):
        with pytest.raises(ValueError, match="down-closure violated"):
            Y.from_cells(cells)

    def test_boundary_sequence_entries_must_be_0_or_1(self):
        with pytest.raises(ValueError, match="boundary sequence violated"):
            Y.from_boundary_sequence((1, 2, 0))

    @pytest.mark.parametrize(
        ("arms", "legs"),
        [((2, 1), (0,)), ((1, 1), (1, 0)), ((1, 0), (0, 0)), ((0,), (-1,))],
    )
    def test_invalid_frobenius_coordinates_are_rejected(self, arms, legs):
        with pytest.raises(ValueError, match="Frobenius coordinates violated"):
            Y.from_frobenius(arms, legs)

    def test_subsets_outside_the_ground_set_are_rejected(self):
        with pytest.raises(ValueError, match="subset violated"):
            Y.from_subset({0, 2}, 4)
        with pytest.raises(ValueError, match="subset violated"):
            Y.from_subset({5}, 4)

    def test_diagrams_that_do_not_fit_the_box_are_rejected(self):
        with pytest.raises(ValueError, match="box violated"):
            Y((3, 1)).to_subset(2, 4)
        with pytest.raises(ValueError, match="box violated"):
            Y((1, 1, 1)).boundary_word(2, 5)
        with pytest.raises(ValueError, match="box violated"):
            Y((3, 1)).to_le_diagram(1, 5)
        with pytest.raises(ValueError, match="box violated"):
            Y((1,)).fits_in_box(3, 2)
        with pytest.raises(ValueError, match="box violated"):
            list(yd.diagrams_in_box(3, 2))
        with pytest.raises(ValueError, match="box violated"):
            Y((1, 1)).padded_rows(1)

    def test_fillings_that_are_not_le_diagrams_are_rejected(self):
        with pytest.raises(ValueError, match="Le-property violated"):
            Y((2, 2)).to_le_diagram(2, 4, {(1, 2), (2, 1)})

    def test_invalid_durfee_dissections_are_rejected(self):
        with pytest.raises(ValueError, match="Durfee dissection violated"):
            Y.from_durfee_dissection(1, Y((1, 1)), Y())
        with pytest.raises(ValueError, match="Durfee dissection violated"):
            Y.from_durfee_dissection(1, Y(), Y((2,)))
        with pytest.raises(ValueError, match="Durfee dissection violated"):
            Y.from_durfee_dissection(-1, Y(), Y())

    def test_queries_about_cells_outside_the_diagram_are_rejected(self):
        example = yd.all_in_one_diagram()
        for method in (
            example.arm_length,
            example.leg_length,
            example.hook_length,
            example.hook,
            example.rim_hook,
            example.remove_rim_hook,
        ):
            with pytest.raises(ValueError, match="is not a cell"):
                method(3, 2)

    def test_adding_or_removing_the_wrong_cell_is_rejected(self):
        example = yd.all_in_one_diagram()
        with pytest.raises(ValueError, match="not an addable position"):
            example.add_cell((3, 3))
        with pytest.raises(ValueError, match="not a corner"):
            example.remove_cell((1, 1))

    def test_skew_and_interval_functions_need_containment(self):
        with pytest.raises(ValueError, match="mu inside lambda"):
            Y((2,)).skew_cells(Y((1, 1)))
        with pytest.raises(ValueError, match="mu inside lambda"):
            yd.mobius(Y((1, 1)), Y((2,)))
        with pytest.raises(ValueError, match="nu inside lambda"):
            yd.count_multichains(Y((1, 1)), Y((2,)), 1)
        with pytest.raises(ValueError, match="n >= 0"):
            yd.count_multichains(Y(), Y((2,)), -1)

    def test_invalid_sizes_are_rejected(self):
        with pytest.raises(ValueError, match="non-negative"):
            yd.partitions(-1)
        with pytest.raises(ValueError, match="n >= 0"):
            yd.staircase(-1)
        with pytest.raises(ValueError, match="non-negative sides"):
            yd.rectangle(-1, 2)
        with pytest.raises(ValueError, match="0 <= k <= n - 1"):
            yd.hook_shape(3, 3)
        with pytest.raises(ValueError, match="positive hook length"):
            yd.all_in_one_diagram().core(0)

    def test_dataframe_missing_column_is_rejected(self):
        frame = yd.all_in_one_diagram().to_dataframe().drop(columns=["column"])
        with pytest.raises(ValueError, match="missing columns"):
            Y.from_dataframe(frame)

    def test_dataframe_with_a_repeated_cell_is_rejected(self):
        frame = yd.all_in_one_diagram().to_dataframe()
        with pytest.raises(ValueError, match="repeats a cell"):
            Y.from_dataframe(pd.concat([frame, frame.iloc[:1]]))

    def test_dataframe_decoding_revalidates_down_closure(self):
        frame = yd.all_in_one_diagram().to_dataframe()
        with pytest.raises(ValueError, match="down-closure violated"):
            Y.from_dataframe(frame.iloc[1:])


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestPlots:
    @pytest.mark.parametrize("convention", list(yd.DrawingConvention))
    @pytest.mark.parametrize("diagram", [yd.all_in_one_diagram(), yd.empty_diagram()])
    def test_plot_diagram_draws_one_outline_per_cell(self, convention, diagram):
        figure, ax = plt.subplots()
        try:
            assert diagram.plot_diagram(ax, convention=convention) is ax
            assert len(ax.lines) == diagram.size
            assert convention.value in ax.get_title()
        finally:
            plt.close(figure)

    def test_plot_hook_lengths_writes_every_hook_length(self):
        example = yd.all_in_one_diagram()
        figure, ax = plt.subplots()
        try:
            assert example.plot_hook_lengths(ax) is ax
            assert sorted(int(text.get_text()) for text in ax.texts) == sorted(
                h for row in example.hook_lengths for h in row
            )
        finally:
            plt.close(figure)

    @pytest.mark.parametrize("method", ["plot_diagram", "plot_hook_lengths"])
    def test_plots_create_axes_when_none_are_given(self, method):
        ax = getattr(yd.all_in_one_diagram(), method)()
        try:
            assert ax.figure is not None
        finally:
            plt.close("all")

    def test_english_rows_go_down_and_french_rows_go_up(self):
        example = yd.adin_roichman_example_9_7_diagram()
        heights = {}
        for convention in (yd.DrawingConvention.ENGLISH, yd.DrawingConvention.FRENCH):
            figure, ax = plt.subplots()
            try:
                example.plot_diagram(ax, convention=convention)
                # Cells are drawn in reading order: the last outline is (2, 1).
                ys = _ydata(ax.lines[-1])
                heights[convention] = (min(ys), max(ys))
            finally:
                plt.close(figure)
        assert heights[yd.DrawingConvention.ENGLISH] == (-2, -1)
        assert heights[yd.DrawingConvention.FRENCH] == (1, 2)

    def test_russian_corner_cell_sits_at_the_bottom_of_the_v(self):
        figure, ax = plt.subplots()
        try:
            yd.all_in_one_diagram().plot_diagram(
                ax, convention=yd.DrawingConvention.RUSSIAN
            )
            lowest = min(min(_ydata(line)) for line in ax.lines)
            first = min(_ydata(ax.lines[0]))
        finally:
            plt.close(figure)
        assert lowest == first == 0

    def test_str_is_the_english_picture(self):
        assert str(yd.all_in_one_diagram()) == "[][][][]\n[][][]\n[]"
        assert str(yd.empty_diagram()) == "()"
