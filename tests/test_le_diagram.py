"""Tests for research.le_diagram, from the Logseq Le-Diagram page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Postnikov Theorem 6.5 and
Corollary 20.1, Ardila-Rincon-Williams Lemma 4.8, Steingrimsson-Williams
Theorem 11, Corteel-Williams Theorem 3.4, Williams's enumeration block and
transpose duality); round-trip laws come from the API contract;
and every part of the definition has a rejection test naming it.

Direction discipline: ``to_decorated_permutation`` is the
Ardila-Rincon-Williams pipe dream, the direction the positroid and necklace
modules exchange. Postnikov's statements (anti-exceedances, the Le-graph's
trip permutation, the cell dimension through alignments) are about the
inverse permutation, so those tests pass ``x.inverse()``.

Not transcribed (spec, "Not transcribed"): the Bruhat-interval bijection
(Theorem 19.1 — the page does not define ``w_lambda`` and Bruhat Order is a
red link), the zig-zag shortcut (the page's quotation elides part of the
rule), the PASEP stationary distribution ("unrestricted rows" is undefined
on the page), the cominuscule generalization with its pattern conditions and
enumeration (no type-B/D structure), the q-Eulerian interpolation values
(the page does not define the renormalization; tested with the Decorated
Permutation suite), and Postnikov's Figure 6.2 example (the page records its
shape and rank but not its filling). Weakened by necessity: the page asserts
the Steingrimsson-Williams map Phi is a bijection (their Theorems 7 and 11,
Lemma 5) without defining Phi, so those are tested as the counts and
distributions they imply.
"""

import itertools
import math
import random
from collections.abc import Sized
from fractions import Fraction
from typing import cast

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import assume, find, given, settings
from hypothesis import strategies as st

from experiments import io
from research import boundary_measurement as bm
from research import decorated_permutation as dp
from research import grassmann_necklace as gn
from research import le_diagram as ld
from research import positroid as ps

matplotlib.use("Agg")

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
# OEIS A000522, the row sums.
A000522 = (1, 2, 5, 16, 65, 326, 1957)

ALL_TYPES = [(k, n) for n in range(7) for k in range(n + 1)]
POSITIVE_TYPES = [(k, n) for k, n in ALL_TYPES if k >= 1]
SMALL_TYPES = [(k, n) for k, n in ALL_TYPES if n <= 5]


# --------------------------------------------------------------------------- #
# Strategies and helpers
# --------------------------------------------------------------------------- #
@st.composite
def shapes(
    draw: st.DrawFn, min_n: int = 0, max_n: int = 7
) -> tuple[tuple[int, ...], int]:
    """Shapes in a box; half the draws are roomy (two or more rows of width >= 2).

    Only a type with ``2 <= k <= n - 2`` has room for a 2 x 2 block, so a
    uniform choice of type leaves the Le-property almost never binding.
    """
    roomy = max_n >= 4 and draw(st.booleans())
    n = draw(st.integers(max(min_n, 4) if roomy else min_n, max_n))
    k = draw(st.integers(2, n - 2) if roomy else st.integers(0, n))
    low = 2 if roomy else 0
    lengths = draw(st.lists(st.integers(low, n - k), min_size=k, max_size=k))
    return tuple(sorted(lengths, reverse=True)), n


@st.composite
def raw_fillings(draw: st.DrawFn, max_n: int = 7) -> tuple[list[list[int]], int]:
    """Arbitrary 0/1 fillings of a shape, Le-diagrams or not."""
    shape, n = draw(shapes(max_n=max_n))
    rows = [draw(st.lists(st.integers(0, 1), min_size=w, max_size=w)) for w in shape]
    return rows, n


@st.composite
def non_le_fillings(draw: st.DrawFn) -> tuple[list[list[int]], int]:
    """Arbitrary fillings with the forbidden pattern planted at a random box."""
    n = draw(st.integers(4, 7))
    k = draw(st.integers(2, n - 2))
    wide = st.integers(2, n - k)
    narrow = st.lists(st.integers(0, n - k), min_size=k - 2, max_size=k - 2)
    shape = sorted([draw(wide), draw(wide), *draw(narrow)], reverse=True)
    rows = [draw(st.lists(st.integers(0, 1), min_size=w, max_size=w)) for w in shape]
    i2 = draw(st.integers(1, sum(1 for row in rows if len(row) > 1) - 1))
    j2 = draw(st.integers(1, len(rows[i2]) - 1))
    i1 = draw(st.integers(0, i2 - 1))
    j1 = draw(st.integers(0, j2 - 1))
    rows[i1][j2], rows[i2][j1], rows[i2][j2] = 1, 1, 0
    return rows, n


@st.composite
def le_diagrams(draw: st.DrawFn, min_n: int = 0, max_n: int = 7) -> ld.LeDiagram:
    """Random Le-diagrams through ``from_filling``.

    Boxes are drawn in reading order; a drawn 0 is kept unless a 1 already
    sits above it and to its left, where the Le-property forces a 1.
    """
    shape, n = draw(shapes(min_n=min_n, max_n=max_n))
    # Mostly 1's in half the draws, so 2 x 2 blocks of 1's actually occur.
    entry = draw(st.sampled_from([st.integers(0, 1), st.sampled_from([1, 1, 1, 0])]))
    rows: list[list[int]] = []
    for width in shape:
        row: list[int] = []
        for j in range(width):
            forced = any(row) and any(above[j] for above in rows)
            row.append(1 if forced else draw(entry))
        rows.append(row)
    return ld.LeDiagram.from_filling(rows, n)


@st.composite
def decorated_permutations(draw: st.DrawFn, max_n: int = 7) -> dp.DecoratedPermutation:
    n = draw(st.integers(0, max_n))
    targets = tuple(draw(st.permutations(range(1, n + 1))))
    fixed = [i for i, target in enumerate(targets, start=1) if target == i]
    clockwise = frozenset(i for i in fixed if draw(st.booleans()))
    return dp.DecoratedPermutation(targets, clockwise)


@st.composite
def le_tableaux(
    draw: st.DrawFn,
) -> tuple[ld.LeDiagram, dict[tuple[int, int], Fraction]]:
    """A Le-diagram with a Le-tableau: positive rationals on its 1-boxes."""
    diagram = draw(le_diagrams(min_n=1, max_n=5))
    positive = st.fractions(min_value=Fraction(1, 9), max_value=9)
    weights = {box: draw(positive) for box in sorted(diagram.ones)}
    return diagram, weights


def _all(k: int, n: int) -> list[ld.LeDiagram]:
    return list(ld.enumerate_le_diagrams(k, n))


def _brute_force(k: int, n: int) -> set[ld.LeDiagram]:
    """Every Le-diagram of type (k, n), by filtering all 0/1 fillings."""
    found = set()
    for shape in ld.enumerate_shapes(k, n):
        for bits in itertools.product((0, 1), repeat=sum(shape)):
            entries = iter(bits)
            rows = [[next(entries) for _ in range(width)] for width in shape]
            if ld.satisfies_le_condition(rows, ld.LeCondition.FORBIDDEN_PATTERN):
                found.add(ld.LeDiagram.from_filling(rows, n))
    return found


def _evaluate(coefficients: tuple[int, ...], q: int) -> int:
    return sum(c * q**power for power, c in enumerate(coefficients))


def _weak_excedance_count(perm: tuple[int, ...]) -> int:
    """Steingrimsson-Williams Theorem 11: positions with ``pi(i) >= i``."""
    return sum(1 for i, target in enumerate(perm, start=1) if target >= i)


def _crossing_count(perm: tuple[int, ...]) -> int:
    """Corteel-Williams: ``j < i <= pi(j) < pi(i)`` or ``pi(i) < pi(j) < i < j``."""
    n = len(perm)
    return sum(
        1
        for i, j in itertools.permutations(range(1, n + 1), 2)
        if j < i <= perm[j - 1] < perm[i - 1] or perm[i - 1] < perm[j - 1] < i < j
    )


def _q_eulerian(k: int, n: int) -> tuple[int, ...]:
    """Williams's E_{k,n}(q), from the page's definition rather than a formula.

    It "q-counts permutations with k weak excedances by alignments": a
    permutation with ``l`` alignments contributes ``q^{k(n-k)-l}``, the
    grading of the alignment block. Fixed points of a permutation are weak
    excedances, and Williams's alignments are Postnikov's for the inverse.
    """
    counts = [0] * (k * (n - k) + 1)
    for targets in itertools.permutations(range(1, n + 1)):
        fixed = frozenset(i for i, t in enumerate(targets, start=1) if t == i)
        x = dp.DecoratedPermutation(targets, fixed)
        if x.weak_excedance_count == k:
            counts[k * (n - k) - x.inverse().alignment_number] += 1
    return tuple(counts)


def _same_row(
    diagram: ld.LeDiagram, tail: int | tuple[int, int], head: int | tuple[int, int]
) -> bool:
    """Whether a Gamma-network edge is horizontal: it stays in one row."""
    rows = dict(zip(diagram.row_labels, range(1, diagram.k + 1), strict=True))

    def row_of(vertex: int | tuple[int, int]) -> int:
        return vertex[0] if isinstance(vertex, tuple) else rows.get(vertex, 0)

    return row_of(tail) == row_of(head)


def _padded_sum(terms: list[tuple[int, tuple[int, ...]]]) -> tuple[int, ...]:
    """Sum ``coefficient * polynomial`` and strip trailing zeros."""
    length = max(len(polynomial) for _, polynomial in terms)
    total = [0] * length
    for coefficient, polynomial in terms:
        for power, value in enumerate(polynomial):
            total[power] += coefficient * value
    while len(total) > 1 and total[-1] == 0:
        total.pop()
    return tuple(total)


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_empty_fillings_are_the_binom_n_k_rank_0_diagrams(self, k, n):
        empties = [ld.empty_filling(shape, n) for shape in ld.enumerate_shapes(k, n)]
        assert len(set(empties)) == math.comb(n, k)
        assert {d.rank for d in empties} == {0}
        assert {d for d in _all(k, n) if d.rank == 0} == set(empties)

    @pytest.mark.parametrize(("k", "n"), SMALL_TYPES)
    def test_empty_fillings_index_the_single_basis_positroids(self, k, n):
        """Rank 0 <-> point cells: one basis each, all binom(n, k) of them."""
        bases = {
            ld.empty_filling(shape, n).to_positroid().bases
            for shape in ld.enumerate_shapes(k, n)
        }
        assert all(len(family) == 1 for family in bases)
        assert len(bases) == math.comb(n, k)

    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_all_ones_rectangle_has_the_maximal_rank_k_n_minus_k(self, k, n):
        top = ld.top_cell_diagram(k, n)
        assert top.shape == (n - k,) * k
        assert top.rank == k * (n - k)
        assert max(d.rank for d in _all(k, n)) == k * (n - k)
        assert [d for d in _all(k, n) if d.rank == k * (n - k)] == [top]

    @pytest.mark.parametrize(("k", "n"), SMALL_TYPES)
    def test_all_ones_rectangle_indexes_the_uniform_matroid(self, k, n):
        """U_{k,n}: every k-subset of the boundary labels is a basis."""
        subsets = itertools.combinations(range(1, n + 1), k)
        expected = frozenset(map(frozenset, subsets))
        assert ld.top_cell_diagram(k, n).to_positroid().bases == expected

    def test_arw_example_is_the_page_data(self):
        example = ld.arw_section_4_3_diagram()
        assert example.shape == (5, 5, 3, 2)
        assert example.diagram_type == (4, 10)
        assert example.to_plus_filling() == ("0+0+0", "+++++", "000", "++")

    def test_arw_example_maps_to_the_stated_decorated_permutation(self):
        image = ld.arw_section_4_3_diagram().to_decorated_permutation()
        assert image.targets == (1, 7, 9, 3, 2, 6, 5, 10, 4, 8)
        assert image == dp.arw_section_4_3_example()

    def test_arw_example_has_both_fixed_point_decorations(self):
        """All-zero row 3 gives overlined 6; label 1 lies outside the shape."""
        example = ld.arw_section_4_3_diagram()
        image = example.to_decorated_permutation()
        assert example.filling[2] == (0, 0, 0)
        assert example.row_labels[2] == 6
        assert example.pipes()[1] == ()
        assert image == dp.DecoratedPermutation.from_underline_overline(
            image.targets, overlined={6}, underlined={1}
        )

    @pytest.mark.parametrize("rows", ld.forbidden_pattern_fillings())
    def test_forbidden_pattern_in_a_2_by_2_square_is_rejected(self, rows):
        assert rows[0][1] == 1
        assert rows[1] == (1, 0)
        with pytest.raises(ValueError, match="Le-property"):
            ld.LeDiagram.from_filling(rows, 4)

    def test_forbidden_pattern_covers_both_values_at_the_corner(self):
        assert {rows[0][0] for rows in ld.forbidden_pattern_fillings()} == {0, 1}

    def test_a_2_4_is_williams_table_1(self):
        """q^4 + 4q^3 + 10q^2 + 12q + 6, 33 cells in total."""
        assert ld.rank_generating_polynomial(2, 4) == (6, 12, 10, 4, 1)
        assert sum(ld.rank_generating_polynomial(2, 4)) == 33

    def test_a_2_5_is_williams_table_1(self):
        """q^6 + 5q^5 + 15q^4 + 30q^3 + 40q^2 + 30q + 10."""
        assert ld.rank_generating_polynomial(2, 5) == (10, 30, 40, 30, 15, 5, 1)


# --------------------------------------------------------------------------- #
# Definition: the four statements of the Le-property
# --------------------------------------------------------------------------- #
class TestDefinition:
    @given(non_le_fillings())
    def test_every_statement_rejects_a_planted_forbidden_pattern(self, case):
        rows, _ = case
        verdicts = [ld.satisfies_le_condition(rows, c) for c in ld.LeCondition]
        assert verdicts == [False] * 4

    @given(raw_fillings())
    def test_the_four_statements_of_the_le_property_agree(self, case):
        rows, _ = case
        verdicts = {
            condition: ld.satisfies_le_condition(rows, condition)
            for condition in ld.LeCondition
        }
        assert len(set(verdicts.values())) == 1, verdicts

    @given(raw_fillings())
    def test_constructor_accepts_every_filling_meeting_the_definition(self, case):
        rows, n = case
        assume(ld.satisfies_le_condition(rows, ld.LeCondition.DEFINITION_6_1))
        assert ld.LeDiagram.from_filling(rows, n).filling == tuple(map(tuple, rows))

    @given(non_le_fillings())
    def test_constructor_rejects_every_filling_failing_the_definition(self, case):
        rows, n = case
        with pytest.raises(ValueError, match="Le-property"):
            ld.LeDiagram.from_filling(rows, n)

    @given(le_diagrams())
    def test_every_diagram_satisfies_every_statement(self, diagram):
        assert all(diagram.satisfies(condition) for condition in ld.LeCondition)

    def test_strategy_reaches_diagrams_with_a_blocked_zero(self):
        """Guard against a vacuous strategy: the Le-property must bite."""
        blocked = find(
            le_diagrams(),
            lambda d: (
                not d.satisfies(ld.LeCondition.BLOCKED_ZERO)
                or any(
                    (i, j) not in d.ones and any((a, j) in d.ones for a in range(1, i))
                    for i, width in enumerate(d.shape, start=1)
                    for j in range(1, width + 1)
                )
            ),
        )
        assert blocked.satisfies(ld.LeCondition.BLOCKED_ZERO) is True

    @pytest.mark.parametrize(("k", "n"), SMALL_TYPES)
    def test_enumeration_is_exactly_the_brute_force_filter(self, k, n):
        listed = _all(k, n)
        assert len(listed) == len(set(listed))
        assert set(listed) == _brute_force(k, n)

    def test_french_notation_is_the_horizontal_reflection(self):
        """Lam-Williams: "no 0 which has a + below it and a + to its left"."""
        french = [[1, 1], [0, 1, 0]]
        diagram = ld.LeDiagram.from_french_filling(french, 5)
        assert diagram.filling == ((0, 1, 0), (1, 1))
        assert diagram.to_french_filling() == ((1, 1), (0, 1, 0))

    def test_french_forbidden_pattern_is_rejected(self):
        with pytest.raises(ValueError, match="Le-property"):
            ld.LeDiagram.from_french_filling([[1, 0], [0, 1]], 4)


class TestRejections:
    def test_shape_must_be_weakly_decreasing(self):
        with pytest.raises(ValueError, match="shape violated"):
            ld.LeDiagram((1, 2), 5)

    def test_shape_must_be_non_negative(self):
        with pytest.raises(ValueError, match="shape violated"):
            ld.LeDiagram((1, -1), 5)

    def test_shape_must_fit_the_k_by_n_minus_k_rectangle(self):
        with pytest.raises(ValueError, match="type violated"):
            ld.LeDiagram((3, 1), 4)

    def test_more_rows_than_n_is_rejected(self):
        with pytest.raises(ValueError, match="type violated"):
            ld.LeDiagram((0, 0, 0), 2)

    def test_ones_must_lie_inside_the_shape(self):
        with pytest.raises(ValueError, match="filling violated"):
            ld.LeDiagram((2, 1), 4, frozenset({(2, 2)}))

    def test_entries_must_be_zero_or_one(self):
        with pytest.raises(ValueError, match="filling violated"):
            ld.LeDiagram.from_filling([[0, 2]], 3)

    def test_plus_rows_must_use_zero_and_plus(self):
        with pytest.raises(ValueError, match="filling violated"):
            ld.LeDiagram.from_plus_filling(["0-"], 3)

    def test_le_property_violation_names_the_three_boxes(self):
        with pytest.raises(
            ValueError, match=r"Le-property.*\(2, 1\).*\(1, 2\).*\(2, 2\)"
        ):
            ld.LeDiagram((2, 2), 4, frozenset({(1, 2), (2, 1)}))

    def test_statements_need_a_shape(self):
        with pytest.raises(ValueError, match="shape violated"):
            ld.satisfies_le_condition([[1], [1, 0]], ld.LeCondition.HOOK)

    @pytest.mark.parametrize(
        "build",
        [
            lambda: ld.top_cell_diagram(3, 2),
            lambda: ld.rank_generating_polynomial(3, 2),
            lambda: list(ld.enumerate_shapes(3, 2)),
            lambda: list(ld.enumerate_le_diagrams(3, 2)),
            lambda: list(ld.enumerate_permutation_tableaux(-1, 2)),
        ],
    )
    def test_types_need_k_between_0_and_n(self, build):
        with pytest.raises(ValueError, match="0 <= k <= n"):
            build()

    def test_weight_is_defined_for_permutation_tableaux_only(self):
        with pytest.raises(ValueError, match="permutation tableaux only"):
            _ = ld.LeDiagram.from_plus_filling(["+0"], 3).weight


# --------------------------------------------------------------------------- #
# Derived vocabulary
# --------------------------------------------------------------------------- #
class TestDerivedVocabulary:
    @given(le_diagrams())
    def test_rank_counts_the_ones(self, diagram):
        assert diagram.rank == sum(map(sum, diagram.filling))

    @given(le_diagrams())
    def test_boundary_path_has_k_vertical_steps_among_n(self, diagram):
        path = diagram.boundary_path
        assert len(path) == diagram.n
        assert path.count("vertical") == diagram.k
        assert diagram.vertical_steps == {
            label for label, step in enumerate(path, start=1) if step == "vertical"
        }

    @given(le_diagrams())
    def test_shape_is_recovered_from_the_vertical_steps(self, diagram):
        """Postnikov 2.1: row i ends where the i-th vertical step sits."""
        width = diagram.n - diagram.k
        labels = sorted(diagram.vertical_steps)
        recovered = tuple(i + width - label for i, label in enumerate(labels, start=1))
        assert recovered == diagram.shape

    def test_boundary_path_of_the_arw_example(self):
        """Shape 5532 in 4 x 6, read from the north-east corner."""
        example = ld.arw_section_4_3_diagram()
        assert example.vertical_steps == {2, 3, 6, 8}
        assert example.column_labels == (10, 9, 7, 5, 4, 1)

    @given(le_diagrams())
    def test_permutation_tableau_means_every_column_has_a_one(self, diagram):
        columns = range(1, diagram.n - diagram.k + 1)
        expected = all(any(j == column for _, j in diagram.ones) for column in columns)
        assert diagram.is_permutation_tableau == expected

    @pytest.mark.parametrize(("k", "n"), POSITIVE_TYPES)
    def test_permutation_tableau_has_exactly_n_minus_k_columns(self, k, n):
        """SW section 1: the column condition forces n - k columns."""
        columns = range(n - k)
        tableaux = [
            d
            for d in _all(k, n)
            if all(any(row[j : j + 1] == (1,) for row in d.filling) for j in columns)
        ]
        assert tableaux
        assert {d.shape[0] for d in tableaux} == {n - k}
        assert all(d.is_permutation_tableau for d in tableaux)

    def test_permutation_tableau_may_have_empty_rows(self):
        tableau = ld.LeDiagram.from_plus_filling(["+", ""], 3)
        assert tableau.is_permutation_tableau is True
        assert tableau.weight == 0


# --------------------------------------------------------------------------- #
# Structural theorems
# --------------------------------------------------------------------------- #
class TestCellIndexing:
    """Postnikov Theorem 6.5."""

    @given(le_tableaux())
    @settings(max_examples=60, deadline=None)
    def test_every_le_tableau_lands_in_the_cell_of_its_diagram(self, case):
        """Sampling positive Le-tableaux gives a constant positroid."""
        diagram, weights = case
        network = diagram.to_planar_network(weights)
        assert network.to_positroid() == diagram.to_positroid()

    @pytest.mark.parametrize("n", range(1, 6))
    def test_every_diagram_keeps_its_cell_under_random_le_tableaux(self, n):
        """The same oracle, exhaustively over diagrams with seeded tableaux."""
        rng = random.Random(n)  # noqa: S311 - seeded test data, not cryptography.
        for diagram in (d for k in range(n + 1) for d in _all(k, n)):
            weights = {
                box: Fraction(rng.randint(1, 9), rng.randint(1, 9))
                for box in sorted(diagram.ones)
            }
            network = diagram.to_planar_network(weights)
            assert network.to_positroid() == diagram.to_positroid(), diagram

    @given(le_tableaux())
    @settings(max_examples=60, deadline=None)
    def test_network_sources_are_the_vertical_steps(self, case):
        diagram, weights = case
        assert diagram.to_planar_network(weights).source_set == diagram.vertical_steps

    @given(le_diagrams(max_n=6))
    @settings(deadline=None)
    def test_cell_dimension_is_the_rank(self, diagram):
        assert bm.cell_dimension(diagram.to_positroid()) == diagram.rank

    @pytest.mark.parametrize("n", range(1, 6))
    def test_boundary_measurements_send_le_diagrams_bijectively_onto_the_cells(self, n):
        """Distinct diagrams parameterize distinct cells, and they exhaust them.

        The cells are read off the Gamma-networks' boundary measurements,
        independently of the pipe dream; there are A000522(n) cells in all.
        """
        diagrams = [d for k in range(n + 1) for d in _all(k, n)]
        cells = [d.to_planar_network().to_positroid() for d in diagrams]
        assert len(set(cells)) == len(cells) == A000522[n]
        assert cells == [d.to_positroid() for d in diagrams]

    @given(le_diagrams(min_n=1, max_n=6))
    @settings(deadline=None)
    def test_shape_lambda_means_the_cell_lies_in_the_schubert_cell(self, diagram):
        """Theorem 6.5: shape lambda iff the cell lies in Omega_lambda.

        Read as: ``I(lambda)`` is the lexicographically minimal basis of the
        cell's positroid, and the first entry of its Grassmann necklace.
        """
        bases = diagram.to_positroid().bases
        assert min(sorted(basis) for basis in bases) == sorted(diagram.vertical_steps)
        necklace = diagram.to_grassmann_necklace()
        assert necklace.entries[0] == diagram.vertical_steps

    @given(le_diagrams(max_n=6))
    @settings(deadline=None)
    def test_type_k_n_indexes_a_rank_k_positroid_on_n_elements(self, diagram):
        positroid = diagram.to_positroid()
        assert (positroid.rank(), positroid.size) == diagram.diagram_type

    @given(le_diagrams(max_n=6))
    @settings(deadline=None)
    def test_positroid_round_trip(self, diagram):
        assert ld.LeDiagram.from_positroid(diagram.to_positroid()) == diagram

    @given(le_diagrams(max_n=6))
    @settings(deadline=None)
    def test_necklace_round_trip(self, diagram):
        necklace = diagram.to_grassmann_necklace()
        assert ld.LeDiagram.from_grassmann_necklace(necklace) == diagram


class TestDecoratedPermutationBijection:
    """Postnikov Corollary 20.1 / Theorem 20.3 and ARW Lemma 4.8."""

    @given(le_diagrams(min_n=1, max_n=6))
    @settings(max_examples=60, deadline=None)
    def test_pipe_dream_agrees_with_the_le_graph_trip_permutation(self, diagram):
        """Postnikov's direction is the inverse of the ARW pipe dream."""
        trip = diagram.to_plabic_graph().to_decorated_permutation()
        assert trip == diagram.to_decorated_permutation().inverse()

    @pytest.mark.parametrize("n", range(1, 6))
    def test_every_le_graph_trip_permutation_is_the_inverse_pipe_dream(self, n):
        """Corollary 20.1 exhaustively, 4-valent splits included."""
        for diagram in (d for k in range(n + 1) for d in _all(k, n)):
            trip = diagram.to_plabic_graph().to_decorated_permutation()
            assert trip == diagram.to_decorated_permutation().inverse(), diagram

    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_every_anti_exceedance_set_is_i_lambda(self, k, n):
        """Corollary 20.1's anti-exceedance clause, exhaustively."""
        for diagram in _all(k, n):
            image = diagram.to_decorated_permutation()
            assert image.inverse().anti_exceedances == diagram.vertical_steps
            assert image.weak_excedances == diagram.vertical_steps

    @given(le_diagrams())
    def test_anti_exceedance_set_is_i_lambda(self, diagram):
        postnikov = diagram.to_decorated_permutation().inverse()
        assert postnikov.anti_exceedances == diagram.vertical_steps
        assert postnikov.permutation_type == diagram.diagram_type

    @given(le_diagrams())
    def test_type_d_n_gives_d_weak_excedances(self, diagram):
        image = diagram.to_decorated_permutation()
        assert (image.weak_excedance_count, image.size) == diagram.diagram_type

    @given(le_diagrams())
    def test_weak_excedances_are_the_vertical_edge_labels(self, diagram):
        """The page's I(lambda) block, read in the ARW direction.

        ``I(lambda)`` is the anti-exceedance set of the cell's decorated
        permutation, so the weak excedances of its inverse, the pipe dream.
        (SW Lemma 5 states the analogue for their undefined map Phi.)
        """
        image = diagram.to_decorated_permutation()
        assert image.weak_excedances == diagram.vertical_steps
        others = frozenset(range(1, diagram.n + 1)) - image.weak_excedances
        assert others == frozenset(diagram.column_labels)

    @given(le_diagrams())
    def test_pipe_dream_round_trip_from_diagrams(self, diagram):
        image = diagram.to_decorated_permutation()
        assert ld.LeDiagram.from_decorated_permutation(image) == diagram

    @given(decorated_permutations())
    def test_pipe_dream_round_trip_from_decorated_permutations(self, decorated):
        diagram = ld.LeDiagram.from_decorated_permutation(decorated)
        assert diagram.to_decorated_permutation() == decorated
        assert diagram.diagram_type == (decorated.weak_excedance_count, decorated.size)

    @pytest.mark.parametrize("n", range(7))
    def test_pipe_dream_is_a_bijection_onto_decorated_permutations(self, n):
        for k in range(n + 1):
            images = [d.to_decorated_permutation() for d in _all(k, n)]
            expected = {
                x
                for x in dp.enumerate_decorated_permutations(n)
                if x.weak_excedance_count == k
            }
            assert len(set(images)) == len(images)
            assert set(images) == expected

    @given(le_diagrams(max_n=6))
    @settings(deadline=None)
    def test_round_trip_against_positroid_from_decorated_permutation(self, diagram):
        image = diagram.to_decorated_permutation()
        positroid = ps.Positroid.from_decorated_permutation(
            range(1, diagram.n + 1), image
        )
        assert positroid.to_decorated_permutation() == image
        assert ld.LeDiagram.from_positroid(positroid) == diagram

    @given(le_diagrams())
    def test_empty_rows_and_columns_are_the_two_kinds_of_fixed_point(self, diagram):
        """Postnikov section 20: empty rows white (counted), empty columns black."""
        image = diagram.to_decorated_permutation()
        zero_rows = {
            label
            for label, row in zip(diagram.row_labels, diagram.filling, strict=True)
            if not any(row)
        }
        zero_columns = {
            label
            for j, label in enumerate(diagram.column_labels, start=1)
            if all(column != j for _, column in diagram.ones)
        }
        assert image.clockwise_fixed == zero_rows
        assert image.counterclockwise_fixed == zero_columns

    @given(le_diagrams())
    def test_pipes_run_south_east_through_the_shape(self, diagram):
        for path in diagram.pipes().values():
            for (i1, j1), (i2, j2) in itertools.pairwise(path):
                assert (i2 - i1, j2 - j1) in {(0, 1), (1, 0)}


class TestEnumeration:
    @pytest.mark.parametrize(("k", "n"), POSITIVE_TYPES)
    def test_graded_counts_match_williams_theorem_4_1(self, k, n):
        closed_form = dp.dimension_generating_polynomial(k, n)
        assert ld.rank_generating_polynomial(k, n) == closed_form

    @pytest.mark.parametrize("n", range(7))
    def test_total_counts_are_the_a046802_row(self, n):
        counts = tuple(sum(ld.rank_generating_polynomial(k, n)) for k in range(n + 1))
        assert counts == A046802_ROWS[n]

    @pytest.mark.parametrize("n", range(7))
    def test_row_sums_are_a000522(self, n):
        total = sum(1 for k in range(n + 1) for _ in ld.enumerate_le_diagrams(k, n))
        assert total == A000522[n] == dp.count_decorated_permutations(n)

    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_euler_characteristic_is_one(self, k, n):
        """Williams Corollary 4.10: A_{k,n}(-1) = 1."""
        assert _evaluate(ld.rank_generating_polynomial(k, n), -1) == 1

    @pytest.mark.parametrize(("k", "n"), POSITIVE_TYPES)
    def test_alignment_grading(self, k, n):
        """Williams Corollary 5.2: coefficient of q^{k(n-k)-l} counts l alignments.

        ``CB_kn`` is cut out by Williams's ``K(pi) = k`` — the weak
        excedance count — and her alignments are Postnikov's for the
        inverse permutation.
        """
        coefficients = ld.rank_generating_polynomial(k, n)
        by_alignments = [0] * len(coefficients)
        for x in dp.enumerate_decorated_permutations(n):
            if x.weak_excedance_count == k:
                by_alignments[x.inverse().alignment_number] += 1
        assert tuple(reversed(by_alignments)) == coefficients

    @given(le_diagrams())
    def test_rank_is_k_n_minus_k_minus_alignments(self, diagram):
        """The pointwise law behind the alignment grading.

        Theorem 6.5 (``dim = |D|``) combined with Postnikov Proposition
        17.10 (``dim = k(n - k) - alignments``, on the Decorated Permutation
        page); the Le-Diagram page itself states only the graded count.
        """
        postnikov = diagram.to_decorated_permutation().inverse()
        k, n = diagram.diagram_type
        assert diagram.rank == k * (n - k) - postnikov.alignment_number
        assert diagram.rank == postnikov.dimension

    @pytest.mark.parametrize(("k", "n"), POSITIVE_TYPES)
    def test_q_eulerian_sum_williams_lemma_5_1(self, k, n):
        """A_{k,n}(q) = sum_i binom(n, i) E_{k,n-i}(q)."""
        terms = [(math.comb(n, i), _q_eulerian(k, n - i)) for i in range(n - k + 1)]
        assert _padded_sum(terms) == ld.rank_generating_polynomial(k, n)

    @pytest.mark.parametrize(("k", "n"), POSITIVE_TYPES)
    def test_q_eulerian_inversion_williams_corollary_5_2(self, k, n):
        """E_{k,n}(q) = sum_i (-1)^i binom(n, i) A_{k,n-i}(q)."""
        terms = [
            ((-1) ** i * math.comb(n, i), ld.rank_generating_polynomial(k, n - i))
            for i in range(n - k + 1)
        ]
        assert _padded_sum(terms) == _padded_sum([(1, _q_eulerian(k, n))])


class TestPermutationStatistics:
    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_total_counts_are_binomial_sums_of_eulerian_numbers(self, k, n):
        """Postnikov Proposition 23.1: N_kn = sum_r binom(n, r) A_{k,n-r}."""
        eulerian = [
            sum(
                1
                for perm in itertools.permutations(range(1, m + 1))
                if _weak_excedance_count(perm) == k
            )
            for m in range(n + 1)
        ]
        expected = sum(math.comb(n, r) * eulerian[n - r] for r in range(n + 1))
        assert sum(1 for _ in ld.enumerate_le_diagrams(k, n)) == expected

    @pytest.mark.parametrize("n", range(1, 7))
    def test_extremal_alignments_are_counted_by_narayana_and_catalan(self, n):
        """Williams Proposition 6.1 and Corollary 6.2.

        A permutation is a decorated permutation whose fixed points all
        count as weak excedances; Williams's alignments are Postnikov's for
        the inverse.
        """
        alignments: dict[int, list[int]] = {k: [] for k in range(1, n + 1)}
        for targets in itertools.permutations(range(1, n + 1)):
            fixed = frozenset(i for i, t in enumerate(targets, start=1) if t == i)
            x = dp.DecoratedPermutation(targets, fixed)
            alignments[x.weak_excedance_count].append(x.inverse().alignment_number)
        maxima = {k: max(values) for k, values in alignments.items()}
        extremal = {k: values.count(maxima[k]) for k, values in alignments.items()}
        narayana = {
            k: math.comb(n, k) * math.comb(n, k - 1) // n for k in range(1, n + 1)
        }
        assert maxima == {k: (k - 1) * (n - k) for k in range(1, n + 1)}
        assert extremal == narayana
        assert sum(extremal.values()) == math.comb(2 * n, n) // (n + 1)


class TestTransposeDuality:
    """Williams section 4."""

    @given(le_diagrams())
    def test_transpose_is_a_rank_preserving_le_diagram_of_type_n_minus_k(self, diagram):
        k, n = diagram.diagram_type
        mirrored = diagram.transpose()
        assert mirrored.diagram_type == (n - k, n)
        assert mirrored.rank == diagram.rank
        assert mirrored.ones == {(j, i) for i, j in diagram.ones}
        assert all(mirrored.satisfies(condition) for condition in ld.LeCondition)

    @given(le_diagrams())
    def test_transpose_is_an_involution(self, diagram):
        assert diagram.transpose().transpose() == diagram

    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_graded_counts_are_symmetric_in_k_and_n_minus_k(self, k, n):
        mirrored = ld.rank_generating_polynomial(n - k, n)
        assert ld.rank_generating_polynomial(k, n) == mirrored

    @pytest.mark.parametrize(("k", "n"), SMALL_TYPES)
    def test_transposition_maps_le_kn_onto_le_n_minus_k_n(self, k, n):
        assert {d.transpose() for d in _all(k, n)} == set(_all(n - k, n))


class TestPermutationTableaux:
    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_tableaux_count_permutations_with_k_weak_excedances(self, k, n):
        """Steingrimsson-Williams Theorem 11, as the count it implies.

        The page asserts a bijection Phi but does not define it, so only
        equinumerosity can be transcribed.
        """
        tableaux = sum(1 for _ in ld.enumerate_permutation_tableaux(k, n))
        permutations = sum(
            1
            for perm in itertools.permutations(range(1, n + 1))
            if _weak_excedance_count(perm) == k
        )
        assert tableaux == permutations

    @pytest.mark.parametrize("n", range(7))
    def test_tableaux_of_expanse_n_are_counted_by_n_factorial(self, n):
        total = sum(
            1 for k in range(n + 1) for _ in ld.enumerate_permutation_tableaux(k, n)
        )
        assert total == math.factorial(n)

    @pytest.mark.parametrize(("k", "n"), ALL_TYPES)
    def test_weight_is_equidistributed_with_crossings(self, k, n):
        """Corteel-Williams Theorem 3.4 = SW Theorem 7, as the counts it implies.

        The page's Phi is undefined, so the restricted bijection is tested
        only as equidistribution of weight and crossings.
        """
        by_weight: dict[int, int] = {}
        for tableau in ld.enumerate_permutation_tableaux(k, n):
            by_weight[tableau.weight] = by_weight.get(tableau.weight, 0) + 1
        by_crossings: dict[int, int] = {}
        for perm in itertools.permutations(range(1, n + 1)):
            if _weak_excedance_count(perm) == k:
                c = _crossing_count(perm)
                by_crossings[c] = by_crossings.get(c, 0) + 1
        assert by_weight == by_crossings

    @pytest.mark.parametrize(("k", "n"), SMALL_TYPES)
    def test_enumerated_tableaux_are_the_diagrams_meeting_the_column_condition(
        self, k, n
    ):
        expected = [d for d in _all(k, n) if d.is_permutation_tableau]
        assert list(ld.enumerate_permutation_tableaux(k, n)) == expected


class TestRestrictionByShape:
    @given(shapes(max_n=6))
    def test_shape_restriction_lists_exactly_the_diagrams_of_that_shape(self, case):
        shape, n = case
        listed = list(ld.enumerate_le_diagrams_of_shape(shape, n))
        assert len(set(listed)) == len(listed)
        assert set(listed) == {d for d in _all(len(shape), n) if d.shape == shape}


# --------------------------------------------------------------------------- #
# Round-trip laws
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(le_diagrams())
    def test_dataframe_round_trip(self, diagram):
        assert ld.LeDiagram.from_dataframe(diagram.to_dataframe()) == diagram

    @given(le_diagrams(), st.randoms(use_true_random=False))
    def test_dataframe_rows_may_come_in_any_order(self, diagram, random):
        frame = diagram.to_dataframe()
        order = list(range(len(frame)))
        random.shuffle(order)
        assert ld.LeDiagram.from_dataframe(frame.iloc[order]) == diagram

    @given(le_diagrams(min_n=1))
    def test_dataframe_is_tidy_one_row_per_step_and_per_box(self, diagram):
        frame = diagram.to_dataframe()
        assert list(frame.columns) == [
            "kind",
            "label",
            "direction",
            "row",
            "column",
            "value",
        ]
        assert (frame["kind"] == "step").sum() == diagram.n
        assert (frame["kind"] == "box").sum() == sum(diagram.shape)

    @pytest.mark.parametrize(
        "example",
        [
            ld.arw_section_4_3_diagram(),
            ld.top_cell_diagram(2, 4),
            ld.empty_filling((0, 0), 2),
            ld.LeDiagram((), 0),
        ],
    )
    def test_experiment_io_round_trip(self, example, tmp_path):
        path = io.write_result(example.to_dataframe(), tmp_path / "le.json")
        decoded = ld.LeDiagram.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == example

    @given(le_diagrams())
    def test_constructors_agree(self, diagram):
        k, n = diagram.diagram_type
        assert ld.LeDiagram.from_filling(diagram.filling, n) == diagram
        assert ld.LeDiagram.from_plus_filling(diagram.to_plus_filling(), n) == diagram
        assert (
            ld.LeDiagram.from_french_filling(diagram.to_french_filling(), n) == diagram
        )
        assert ld.LeDiagram(diagram.shape, n, diagram.ones) == diagram
        assert len(diagram.shape) == k

    def test_dataframe_missing_column_is_rejected(self):
        frame = ld.top_cell_diagram(1, 2).to_dataframe().drop(columns=["value"])
        with pytest.raises(ValueError, match="missing columns"):
            ld.LeDiagram.from_dataframe(frame)

    def test_dataframe_with_a_missing_step_is_rejected(self):
        frame = ld.top_cell_diagram(1, 3).to_dataframe()
        with pytest.raises(ValueError, match="step label"):
            ld.LeDiagram.from_dataframe(frame[frame["label"] != 1])

    def test_dataframe_with_a_missing_box_is_rejected(self):
        frame = ld.top_cell_diagram(1, 3).to_dataframe()
        with pytest.raises(ValueError, match="one box row per box"):
            ld.LeDiagram.from_dataframe(frame[frame["column"] != 2])

    def test_dataframe_with_an_unknown_direction_is_rejected(self):
        frame = ld.top_cell_diagram(1, 2).to_dataframe()
        frame.loc[0, "direction"] = "diagonal"
        with pytest.raises(ValueError, match="direction must be"):
            ld.LeDiagram.from_dataframe(frame)

    def test_dataframe_with_an_unknown_kind_is_rejected(self):
        frame = ld.top_cell_diagram(1, 2).to_dataframe()
        frame.loc[0, "kind"] = "pipe"
        with pytest.raises(ValueError, match="kind must be"):
            ld.LeDiagram.from_dataframe(frame)

    def test_dataframe_decoding_revalidates_the_le_property(self):
        frame = ld.top_cell_diagram(2, 4).to_dataframe()
        corner = (frame["row"] == 2) & (frame["column"] == 2)
        frame.loc[corner, "value"] = 0
        with pytest.raises(ValueError, match="Le-property"):
            ld.LeDiagram.from_dataframe(frame)


# --------------------------------------------------------------------------- #
# Representation and visualization
# --------------------------------------------------------------------------- #
class TestPresentation:
    def test_repr_shows_the_type_and_the_rows(self):
        text = repr(ld.arw_section_4_3_diagram())
        assert text == "LeDiagram(type=(4, 10), rows=('0+0+0', '+++++', '000', '++'))"

    def test_str_draws_the_rows_and_marks_empty_ones(self):
        assert str(ld.LeDiagram.from_plus_filling(["+0", ""], 4)) == "+0\n."

    def test_diagrams_are_hashable_values(self):
        assert len({ld.top_cell_diagram(2, 4), ld.top_cell_diagram(2, 4)}) == 1

    def test_networks_and_graphs_need_a_boundary(self):
        with pytest.raises(ValueError, match="boundary"):
            ld.LeDiagram((), 0).to_planar_network()

    @pytest.mark.parametrize(
        "weights",
        [{(1, 2): Fraction(1)}, {(1, 1): Fraction(0)}, {(1, 1): Fraction(-2)}],
        ids=["positive-on-a-0-box", "zero-on-a-1-box", "negative-on-a-1-box"],
    )
    def test_le_tableau_must_be_positive_exactly_on_the_ones(self, weights):
        diagram = ld.LeDiagram.from_plus_filling(["+0"], 3)
        with pytest.raises(ValueError, match="positive exactly on the 1-boxes"):
            diagram.to_planar_network(weights)

    def test_le_tableau_may_list_its_zeros_on_the_zero_boxes(self):
        """The page's literal form: T is defined on every box of the shape."""
        diagram = ld.LeDiagram.from_plus_filling(["+0"], 3)
        total = diagram.to_planar_network({(1, 1): 2, (1, 2): 0})
        assert total == diagram.to_planar_network({(1, 1): 2})

    def test_gamma_network_weights_the_horizontal_edge_into_each_one_box(self):
        """Definition 6.3: T(i, j) on the edge into (i, j), vertical edges 1.

        Rows ``++ / 0+`` at n = 4 with T(1,1) = 2, T(1,2) = 3, T(2,2) = 5:
        source 1 reaches sink 3 by the path of weight 3 and sink 4 by the
        path of weight 3 * 2; source 2 reaches sink 3 with weight 5 and never
        reaches sink 4.
        """
        diagram = ld.LeDiagram.from_plus_filling(["++", "0+"], 4)
        network = diagram.to_planar_network({(1, 1): 2, (1, 2): 3, (2, 2): 5})
        assert network.boundary_measurements() == {
            (1, 3): Fraction(3),
            (1, 4): Fraction(6),
            (2, 3): Fraction(5),
            (2, 4): Fraction(0),
        }

    @given(le_tableaux())
    @settings(max_examples=40, deadline=None)
    def test_le_tableaux_parameterize_by_one_positive_real_per_one_box(self, case):
        """Section 6: the tableaux of D form R^{|D|}_{>0}; each weight is used."""
        diagram, weights = case
        network = diagram.to_planar_network(weights)
        horizontal = {
            head: weight
            for tail, head, weight in network.edges
            if head in diagram.ones and _same_row(diagram, tail, head)
        }
        vertical = [
            weight
            for tail, head, weight in network.edges
            if not _same_row(diagram, tail, head)
        ]
        assert len(weights) == diagram.rank
        assert horizontal == weights
        assert set(vertical) <= {Fraction(1)}

    @pytest.mark.parametrize(
        "method", ["plot_diagram", "plot_hook_diagram", "plot_pipe_dream"]
    )
    @pytest.mark.parametrize(
        "diagram",
        [ld.arw_section_4_3_diagram(), ld.empty_filling((0,), 2), ld.LeDiagram((), 0)],
    )
    def test_plots_draw_on_given_axes_and_return_them(self, method, diagram):
        figure, ax = plt.subplots()
        try:
            assert getattr(diagram, method)(ax) is ax
            assert ax.get_title() != ""
        finally:
            plt.close(figure)

    @pytest.mark.parametrize(
        "method", ["plot_diagram", "plot_hook_diagram", "plot_pipe_dream"]
    )
    def test_plots_create_axes_when_none_are_given(self, method):
        ax = getattr(ld.top_cell_diagram(2, 4), method)()
        try:
            assert ax.figure is not None
        finally:
            plt.close("all")

    def test_hook_diagram_draws_one_dot_per_one(self):
        example = ld.arw_section_4_3_diagram()
        ax = example.plot_hook_diagram()
        try:
            (dots,) = ax.collections
            assert len(cast("Sized", dots.get_offsets())) == example.rank
        finally:
            plt.close("all")

    def test_necklace_of_the_top_cell_is_the_uniform_necklace(self):
        assert ld.top_cell_diagram(2, 5).to_grassmann_necklace() == gn.uniform_necklace(
            2, 5
        )
