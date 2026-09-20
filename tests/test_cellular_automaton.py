"""Tests for research.cellular_automaton, from the Logseq Cellular Automaton page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural theorems (Curtis-Hedlund-Lyndon,
reversible = injective, Moore-Myhill, Kari 2005 Theorem 7, the balance
theorem, the pair-graph decision procedures, Ito-Osato-Nasu, the finite
phase-space bound, Martin-Odlyzko-Wolfram section 3, Wolfram's single-seed
formula) and the operations block (composition, shifts, symmetry
conjugation, the Fredkin construction, the Laurent-polynomial calculus);
round-trip laws come from the API contract; and every clause of each
formulation has a rejection test naming it.

Not transcribed (spec, "Not transcribed"): the Garden of Eden theorem over
amenable groups, surjunctivity and the entropy statements (the group is fixed
to Z^d); the undecidability and universality theorems; the structure theorem
for reversible CA; Kurka's classes and the transitivity, sensitivity and
expansivity clauses of the linear-dynamics theorem (Kari 2005, Theorem 18 —
only its equicontinuity clause has a finite oracle, through "equicontinuous
iff G^n = G^m"); Moore's growth bound; the other dynamical implications; the
conservation-law algorithm (cited, not stated); the inverse-neighborhood
bound; the parts of Theorem 7 about G_F being surjective (they quantify over
all finite configurations) and "G_P injective implies G_P surjective" (on one
ring that is finiteness, not the theorem);
the two-dimensional tiling fixture of Theorem 7 and SNAKE-XOR (not specified
on the page); von Neumann's 29-state automaton ("too large for a unit-test
fixture"); Life's 88-cell orphan; the R-pentomino's stabilization beyond a
constant population over a finite window after generation 1103 (its gliders
never stop moving); the rule 90 fractal dimension; rule 150's
positive expansivity; and the Strand Dynamics non-example, which is not a
cellular automaton and so has nothing to construct.
"""

import functools
import itertools
import math
from collections import Counter
from collections.abc import Iterable, Sequence

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from experiments import io
from research import cellular_automaton as ca

matplotlib.use("Agg")

CA = ca.CellularAutomaton

# The line-level decisions build graphs on k^(m-1) vertices; one example may
# exceed Hypothesis's default 200 ms deadline without being wrong.
slow = settings(deadline=None)

# The 30 surjective elementary rules, as printed on the page.
SURJECTIVE = (
    15, 30, 45, 51, 60, 75, 85, 86, 89, 90, 101, 102, 105, 106, 120, 135, 149,
    150, 153, 154, 165, 166, 169, 170, 180, 195, 204, 210, 225, 240,
)  # fmt: skip
# The 6 injective elementary rules, as printed on the page.
INJECTIVE = (15, 51, 85, 170, 204, 240)
# The rules whose shortest orphan has length 9, as printed on the page.
LONGEST_ORPHANS = (37, 91, 164, 218)
# The number-conserving elementary rules, as printed on the page.
NUMBER_CONSERVING = (170, 184, 204, 226, 240)
# The linear elementary rules by coefficient vector 000, 001, ..., 111.
LINEAR_RULES = (0, 170, 204, 102, 240, 90, 60, 150)
# Pi_N for N = 3..30 (MOW Table 1 = OEIS A085587), as printed on the page.
PI_N = (
    1, 1, 3, 2, 7, 1, 7, 6, 31, 4, 63, 14, 15, 1, 15, 14, 511, 12, 63, 62,
    2047, 8, 1023, 126, 511, 28, 16383, 30,
)  # fmt: skip
# OEIS A051023 (rule 30 center column), the digits printed on the page.
RULE_30_COLUMN = "1101110011000101100100111010111001110101"

TRAFFIC_LIGHTS = (
    "..###..",
    ".......",
    "#.....#",
    "#.....#",
    "#.....#",
    ".......",
    "..###..",
)
HONEY_FARM = (
    "......#......",
    ".....#.#.....",
    ".....#.#.....",
    "......#......",
    ".............",
    ".##.......##.",
    "#..#.....#..#",
    ".##.......##.",
    ".............",
    "......#......",
    ".....#.#.....",
    ".....#.#.....",
    "......#......",
)


# --------------------------------------------------------------------------- #
# Strategies and helpers
# --------------------------------------------------------------------------- #
@functools.cache
def _elementary(number: int) -> ca.CellularAutomaton:
    return CA.from_wolfram_number(number)


def _all_elementary() -> list[ca.CellularAutomaton]:
    return [_elementary(number) for number in range(256)]


def _seed(size: int) -> tuple[int, ...]:
    return (1,) + (0,) * (size - 1)


def _normalize(cells: Iterable[tuple[int, ...]]) -> frozenset[tuple[int, ...]]:
    points = list(cells)
    top = min(r for r, _ in points)
    left = min(c for _, c in points)
    return frozenset((r - top, c - left) for r, c in points)


def _picture(rows: Sequence[str]) -> frozenset[tuple[int, ...]]:
    return frozenset(
        (r, c)
        for r, row in enumerate(rows)
        for c, mark in enumerate(row)
        if mark == "#"
    )


def _prime_factors(n: int) -> list[int]:
    return [
        p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, p))
    ]


@functools.cache
def _images(rule: ca.CellularAutomaton, length: int) -> frozenset[tuple[int, ...]]:
    """Brute force: every image word of the given length."""
    width = rule.to_contiguous().neighborhood_size
    return frozenset(
        rule.apply_to_word(word)
        for word in itertools.product(range(rule.states), repeat=length + width - 1)
    )


def _all_words(rule: ca.CellularAutomaton, length: int) -> set[tuple[int, ...]]:
    return set(itertools.product(range(rule.states), repeat=length))


@st.composite
def line_automata(
    draw: st.DrawFn,
    *,
    states: int | None = None,
    max_states: int = 3,
    max_neighbors: int = 3,
    reach: int = 1,
) -> ca.CellularAutomaton:
    """Random one-dimensional automata; ``states`` fixes the alphabet size."""
    k = draw(st.integers(1, max_states)) if states is None else states
    positions = draw(
        st.lists(
            st.integers(-reach, reach), min_size=0, max_size=max_neighbors, unique=True
        )
    )
    table = draw(
        st.lists(
            st.integers(0, k - 1),
            min_size=k ** len(positions),
            max_size=k ** len(positions),
        )
    )
    return CA(1, k, tuple((x,) for x in positions), tuple(table))


@st.composite
def plane_automata(draw: st.DrawFn) -> ca.CellularAutomaton:
    offsets = draw(
        st.lists(
            st.sampled_from(ca.moore_neighborhood(2)),
            min_size=1,
            max_size=4,
            unique=True,
        )
    )
    table = draw(
        st.lists(
            st.integers(0, 1), min_size=2 ** len(offsets), max_size=2 ** len(offsets)
        )
    )
    return CA(2, 2, tuple(offsets), tuple(table))


@st.composite
def ring_configurations(
    draw: st.DrawFn,
    automaton: ca.CellularAutomaton,
    min_size: int = 1,
    max_size: int = 8,
) -> tuple[int, ...]:
    size = draw(st.integers(min_size, max_size))
    return tuple(
        draw(
            st.lists(st.integers(0, automaton.states - 1), min_size=size, max_size=size)
        )
    )


@st.composite
def right_permutive_automata(draw: st.DrawFn) -> ca.CellularAutomaton:
    """Rules f(a_1..a_n) = pi_{a_1..a_(n-1)}(a_n), one permutation per prefix."""
    states = draw(st.integers(2, 3))
    width = draw(st.integers(1, 3))
    table: list[int] = []
    for _ in range(states ** (width - 1)):
        table.extend(draw(st.permutations(range(states))))
    left = draw(st.integers(-1, 0))
    return CA(1, states, tuple((left + i,) for i in range(width)), tuple(table))


@st.composite
def linear_automata(draw: st.DrawFn, max_modulus: int = 6) -> ca.CellularAutomaton:
    modulus = draw(st.integers(2, max_modulus))
    positions = draw(st.lists(st.integers(-1, 1), min_size=1, max_size=3, unique=True))
    coefficients: dict[Sequence[int], int] = {
        (x,): draw(st.integers(0, modulus - 1)) for x in positions
    }
    return CA.from_linear(1, modulus, coefficients)


def injective_pool() -> list[ca.CellularAutomaton]:
    patt = ca.patt_rule()
    return [
        *ca.reversible_elementary_rules(),
        *patt.symmetry_class(),
        ca.shift_automaton(1, 3, 0),
        _elementary(30).second_order(),
        _elementary(110).second_order(),
    ]


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_identity_complement_and_shifts_are_the_trivial_reversible_class(self):
        rules = ca.reversible_elementary_rules()
        assert [rule.wolfram_number for rule in rules] == [204, 51, 170, 240, 85, 15]
        assert all(rule.is_reversible() for rule in rules)
        assert all(rule.minimal().neighborhood_size == 1 for rule in rules)

    def test_rule_170_is_the_left_shift_and_rule_240_reads_the_left_cell(self):
        assert _elementary(170) == CA.from_block_map(2, 1, 1, lambda n: n[2])
        assert _elementary(240) == CA.from_block_map(2, 1, 1, lambda n: n[0])
        assert _elementary(170).same_global_map(ca.shift_automaton(1, 2, 0))

    def test_every_other_elementary_rule_fails_injectivity(self):
        injective = [r.wolfram_number for r in _all_elementary() if r.is_injective()]
        assert injective == sorted(
            rule.wolfram_number for rule in ca.reversible_elementary_rules()
        )

    def test_rule_90_is_additive_surjective_and_not_injective(self):
        rule = ca.rule_90()
        assert rule == CA.from_block_map(2, 1, 1, lambda n: (n[0] + n[2]) % 2)
        assert rule.is_additive
        assert rule.is_surjective()
        assert not rule.is_injective()

    @pytest.mark.parametrize("size", range(1, 15))
    def test_rule_90_line_versus_ring_discrepancy(self, size):
        space = ca.rule_90().phase_space(size)
        fraction = (1, 2) if size % 2 else (3, 4)
        assert len(space.gardens_of_eden) * fraction[1] == 2**size * fraction[0]
        assert not space.is_bijective

    @pytest.mark.parametrize("size", range(1, 15))
    def test_rule_150_is_bijective_on_a_ring_exactly_when_3_does_not_divide_it(
        self, size
    ):
        rule = ca.rule_150()
        assert rule == CA.from_block_map(2, 1, 1, lambda n: sum(n) % 2)
        assert rule.is_additive
        assert rule.phase_space(size).is_bijective == (size % 3 != 0)

    def test_rule_102_is_pre_injective_but_not_injective(self):
        rule = ca.rule_102()
        assert rule == CA.from_block_map(2, 1, 1, lambda n: (n[1] + n[2]) % 2)
        assert rule.is_pre_injective()
        assert not rule.is_injective()

    @pytest.mark.parametrize("size", range(1, 8))
    def test_rule_102_maps_the_two_constant_configurations_to_one_image(self, size):
        rule = ca.rule_102()
        assert rule.step((0,) * size) == rule.step((1,) * size)

    @pytest.mark.parametrize(
        ("constructor", "number"),
        [
            (ca.rule_30, 30),
            (ca.rule_90, 90),
            (ca.rule_102, 102),
            (ca.rule_110, 110),
            (ca.rule_150, 150),
            (ca.rule_184, 184),
            (ca.rule_226, 226),
            (ca.rule_232, 232),
        ],
    )
    def test_named_elementary_rules_carry_their_wolfram_numbers(
        self, constructor, number
    ):
        assert constructor() == CA.from_wolfram_number(number)
        assert constructor().wolfram_number == number

    def test_rule_30_is_surjective_and_not_injective(self):
        assert ca.rule_30().is_surjective()
        assert not ca.rule_30().is_injective()

    def test_rule_30_symmetry_class(self):
        numbers = {rule.wolfram_number for rule in ca.rule_30().symmetry_class()}
        assert numbers == {30, 86, 135, 149}

    def test_rule_30_center_column_from_a_single_one(self):
        rule = ca.rule_30()
        rows = rule.evolve_finite({(0,): 1}, len(RULE_30_COLUMN) - 1)
        column = "".join(str(row.get((0,), 0)) for row in rows)
        assert column == RULE_30_COLUMN

    def test_rule_110_local_rule_is_cooks(self):
        rule = ca.rule_110()
        ones = {(0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)}
        for neighbors in itertools.product(range(2), repeat=3):
            assert rule.local_rule(*neighbors) == int(neighbors in ones)

    def test_rule_110_is_unbalanced_and_not_surjective(self):
        rule = ca.rule_110()
        assert Counter(rule.table) == {0: 3, 1: 5}
        assert not rule.is_surjective()

    def test_rule_110_has_the_unique_shortest_orphan_01010(self):
        assert ca.rule_110().shortest_orphans() == ((0, 1, 0, 1, 0),)

    def test_rule_110_diamond_words_share_the_image_111110(self):
        rule = ca.rule_110()
        first, second = (0, 0, 1, 1, 0, 1, 0, 0), (0, 0, 1, 0, 1, 1, 0, 0)
        assert rule.apply_to_word(first) == (1, 1, 1, 1, 1, 0)
        assert rule.apply_to_word(second) == (1, 1, 1, 1, 1, 0)
        assert first[:2] == second[:2]
        assert first[-2:] == second[-2:]
        assert not rule.is_pre_injective()

    def test_rule_110_symmetry_class(self):
        numbers = {rule.wolfram_number for rule in ca.rule_110().symmetry_class()}
        assert numbers == {110, 124, 137, 193}

    @given(st.lists(st.integers(0, 1), min_size=2, max_size=10))
    def test_rule_226_replaces_pattern_01_by_pattern_10(self, cells):
        size = len(cells)
        expected = list(cells)
        for i in range(size):
            if (cells[i], cells[(i + 1) % size]) == (0, 1):
                expected[i], expected[(i + 1) % size] = 1, 0
        assert ca.rule_226().step(cells) == tuple(expected)

    def test_rule_184_is_the_mirror_image_of_rule_226(self):
        assert ca.rule_226().reflect() == ca.rule_184()

    def test_number_conserving_elementary_rules_on_rings_up_to_9_cells(self):
        conserving = [
            rule.wolfram_number
            for rule in _all_elementary()
            if all(rule.conserves_number_on_ring(size) for size in range(1, 10))
        ]
        assert conserving == sorted(NUMBER_CONSERVING)

    def test_rule_184_is_not_surjective_with_unique_shortest_orphan_1100(self):
        assert not ca.rule_184().is_surjective()
        assert ca.rule_184().shortest_orphans() == ((1, 1, 0, 0),)

    def test_rule_232_is_balanced_but_neither_surjective_nor_pre_injective(self):
        rule = ca.rule_232()
        assert rule == CA.from_block_map(2, 1, 1, lambda n: int(sum(n) >= 2))
        assert rule.is_balanced
        assert not rule.is_surjective()
        assert not rule.is_pre_injective()

    def test_rule_232_orphans(self):
        rule = ca.rule_232()
        assert rule.is_orphan((0, 1, 0, 0, 1))
        assert rule.shortest_orphans() == (
            (0, 1, 0, 0, 1),
            (0, 1, 1, 0, 1),
            (1, 0, 0, 1, 0),
            (1, 0, 1, 1, 0),
        )

    def test_rule_232_words_00000_and_00100_are_mutually_erasable(self):
        rule = ca.rule_232()
        assert rule.apply_to_word((0, 0, 0, 0, 0)) == rule.apply_to_word(
            (0, 0, 1, 0, 0)
        )

    def test_patt_rule_flips_the_second_cell_exactly_in_the_landscape_0_1_0(self):
        rule = ca.patt_rule()
        assert rule.offsets == ((-1,), (0,), (1,), (2,))
        for a, b, c, d in itertools.product(range(2), repeat=4):
            flipped = (a, c, d) == (0, 1, 0)
            assert rule.local_rule(a, b, c, d) == (1 - b if flipped else b)

    def test_patt_rule_is_a_nontrivial_reversible_rule(self):
        rule = ca.patt_rule()
        assert rule.is_injective()
        assert rule.minimal().neighborhood_size == 4

    @pytest.mark.parametrize("size", range(5, 15))
    def test_patt_rule_is_an_involution_on_rings(self, size):
        rule = ca.patt_rule()
        space = rule.phase_space(size)
        assert all(
            space.successor[image] == i for i, image in enumerate(space.successor)
        )

    def test_exactly_eight_nontrivial_injective_maps_on_four_contiguous_cells(self):
        offsets = ((-1,), (0,), (1,), (2,))
        injective = []
        # Injective implies surjective implies balanced (Corollary 14), and
        # bijective on rings (Theorem 7): both prune before the pair graph.
        for ones in itertools.combinations(range(16), 8):
            rule = CA(1, 2, offsets, tuple(int(i in ones) for i in range(16)))
            if rule.phase_space(3).is_bijective and rule.is_injective():
                injective.append(rule)
        nontrivial = [r for r in injective if r.minimal().neighborhood_size > 1]
        assert len(nontrivial) == 8
        assert sum(0 in rule.quiescent_states for rule in nontrivial) == 4
        assert ca.patt_rule() in nontrivial
        # "variants (obtained by reflection or complementation) of a single
        # cellular automaton": the conjugation class of Patt's rule, and the
        # same four tables with the output complemented.
        variants = {rule.table for rule in ca.patt_rule().symmetry_class()}
        variants |= {tuple(1 - value for value in table) for table in variants}
        assert {rule.table for rule in nontrivial} == variants

    def test_patt_table_on_a_gapped_neighborhood_is_not_injective(self):
        gapped = ca.patt_table_on_gapped_neighborhood()
        assert gapped.offsets == ((-4,), (-2,), (-1,), (0,))
        assert gapped.table == ca.patt_rule().table
        assert not gapped.is_injective()

    @pytest.mark.parametrize(
        ("pattern", "picture"),
        [
            (ca.life_block, ("##", "##")),
            (ca.life_beehive, (".##.", "#..#", ".##.")),
            (ca.life_blinker, ("###",)),
            (ca.life_glider, (".#.", "..#", "###")),
            (ca.life_t_tetromino, ("###", ".#.")),
            (ca.life_r_pentomino, (".##", "##.", ".#.")),
            (lambda: ca.life_row(7), ("#######",)),
            (lambda: ca.life_row(10), ("##########",)),
        ],
    )
    def test_life_patterns_are_the_named_cell_sets(self, pattern, picture):
        cells = pattern()
        assert set(cells.values()) == {1}
        assert _normalize(cells) == _picture(picture)

    def test_life_rule_is_gardners(self):
        life = ca.game_of_life()
        assert (life.dimension, life.states) == (2, 2)
        assert life.offsets == ca.moore_neighborhood(2)
        for cells in itertools.product(range(2), repeat=9):
            around = sum(cells) - cells[4]
            alive = around in (2, 3) if cells[4] else around == 3
            assert life.local_rule(*cells) == int(alive)

    def test_life_empty_and_single_cell_have_the_same_image(self):
        life = ca.game_of_life()
        assert life.step_finite({}) == {}
        assert life.step_finite({(0, 0): 1}) == {}

    @pytest.mark.parametrize("pattern", [ca.life_block, ca.life_beehive])
    def test_life_still_lifes(self, pattern):
        assert ca.game_of_life().step_finite(pattern()) == pattern()

    def test_life_blinker_has_period_2(self):
        rows = ca.game_of_life().evolve_finite(ca.life_blinker(), 2)
        assert rows[1] != rows[0]
        assert rows[2] == rows[0]

    def test_life_glider_has_period_4_and_moves_one_cell_diagonally(self):
        rows = ca.game_of_life().evolve_finite(ca.life_glider(), 4)
        assert set(rows[4]) == {(r + 1, c + 1) for r, c in rows[0]}
        assert all(_normalize(rows[t]) != _normalize(rows[0]) for t in (1, 2, 3))

    def test_life_t_tetromino_becomes_traffic_lights_at_generation_9(self):
        rows = ca.game_of_life().evolve_finite(ca.life_t_tetromino(), 11)
        assert _normalize(rows[9]) == _picture(TRAFFIC_LIGHTS)
        assert _normalize(rows[8]) != _picture(TRAFFIC_LIGHTS)
        assert rows[11] == rows[9]

    def test_life_row_of_7_becomes_a_honey_farm_at_generation_14(self):
        rows = ca.game_of_life().evolve_finite(ca.life_row(7), 15)
        assert _normalize(rows[14]) == _picture(HONEY_FARM)
        assert _normalize(rows[13]) != _picture(HONEY_FARM)
        assert rows[15] == rows[14]

    def test_life_row_of_10_becomes_the_pentadecathlon_of_period_15(self):
        rows = ca.game_of_life().evolve_finite(ca.life_row(10), 40)
        assert rows[25] == rows[10]
        assert all(rows[10 + p] != rows[10] for p in range(1, 15))

    def test_life_r_pentomino_stabilizes_at_generation_1103_with_population_116(self):
        life = ca.game_of_life()
        cells = ca.life_r_pentomino()
        populations = [len(cells)]
        for _ in range(1130):
            cells = life.step_finite(cells)
            populations.append(len(cells))
        assert populations[1102] != 116
        assert set(populations[1103:]) == {116}

    def test_bridge_example_nand_on_the_directed_3_cycle(self):
        rule = ca.one_way_nand()
        assert rule.offsets == ((-1,), (0,))
        space = rule.phase_space(3)
        gardens = {space.configuration(i) for i in space.gardens_of_eden}
        assert gardens == {(0, 0, 1), (0, 1, 0), (1, 0, 0)}
        assert space.cycle_lengths == (2, 3)
        assert rule.step((1, 1, 1)) == (0, 0, 0)
        assert rule.step((0, 0, 0)) == (1, 1, 1)
        three_cycle = {(1, 1, 0), (1, 0, 1), (0, 1, 1)}
        assert {rule.step(c) for c in three_cycle} == three_cycle
        assert space.fixed_points == ()

    def test_bridge_example_is_the_mirror_image_of_the_radius_half_nand(self):
        radius_half = CA.from_local_rule(
            1, 2, ca.one_way_neighborhood(1), lambda n: 1 - n[0] * n[1]
        )
        assert ca.one_way_nand().reflect().same_global_map(radius_half)


# --------------------------------------------------------------------------- #
# Census of the elementary rules
# --------------------------------------------------------------------------- #
class TestElementaryCensus:
    def test_there_are_256_rules_with_distinct_tables(self):
        assert len({rule.table for rule in _all_elementary()}) == 256

    def test_88_symmetry_classes_of_sizes_4_2_1(self):
        classes = {
            frozenset(r.wolfram_number for r in rule.symmetry_class())
            for rule in _all_elementary()
        }
        assert len(classes) == 88
        assert Counter(len(c) for c in classes) == {4: 44, 2: 36, 1: 8}

    def test_32_legal_rules(self):
        assert sum(rule.is_legal for rule in _all_elementary()) == 32

    def test_30_surjective_rules(self):
        surjective = [r.wolfram_number for r in _all_elementary() if r.is_surjective()]
        assert surjective == list(SURJECTIVE)

    def test_6_injective_rules(self):
        injective = [r.wolfram_number for r in _all_elementary() if r.is_injective()]
        assert injective == list(INJECTIVE)

    def test_70_balanced_tables_of_which_30_are_surjective(self):
        balanced = [rule for rule in _all_elementary() if rule.is_balanced]
        assert len(balanced) == 70
        assert sum(rule.is_surjective() for rule in balanced) == 30

    def test_additive_legal_rules(self):
        additive_legal = [
            rule.wolfram_number
            for rule in _all_elementary()
            if rule.is_additive and rule.is_legal
        ]
        assert additive_legal == [0, 90, 150, 204]


# --------------------------------------------------------------------------- #
# Formulations
# --------------------------------------------------------------------------- #
class TestFormulations:
    @pytest.mark.parametrize("number", range(256))
    def test_wolfram_number_reads_the_bits_f111_to_f000(self, number):
        rule = _elementary(number)
        bits = "".join(
            str(rule.local_rule(a, b, c))
            for a, b, c in itertools.product((1, 0), repeat=3)
        )
        assert int(bits, 2) == number
        assert rule.wolfram_number == number
        assert number == sum(
            rule.local_rule(a, b, c) << (4 * a + 2 * b + c)
            for a, b, c in itertools.product(range(2), repeat=3)
        )

    @given(line_automata())
    def test_d1_and_d2_agree_when_the_memory_set_is_the_entries_of_n(self, rule):
        from_memory = CA.from_memory_set(
            1,
            rule.states,
            rule.offsets,
            lambda pattern: rule.local_rule(*(pattern[o] for o in rule.offsets)),
        )
        assert from_memory.same_global_map(rule)
        assert from_memory.offsets == tuple(sorted(rule.offsets))

    @given(line_automata())
    def test_local_rule_constructor_tabulates_f(self, rule):
        rebuilt = CA.from_local_rule(
            1, rule.states, rule.offsets, lambda n: rule.local_rule(*n)
        )
        assert rebuilt == rule

    def test_d4_is_d1_with_neighborhood_minus_m_to_n(self):
        block = CA.from_block_map(3, 2, 1, lambda n: (n[0] + 2 * n[3]) % 3)
        assert block.dimension == 1
        assert block.offsets == ((-2,), (-1,), (0,), (1,))
        assert block == CA.from_local_rule(
            1, 3, [(-2,), (-1,), (0,), (1,)], lambda n: (n[0] + 2 * n[3]) % 3
        )

    def test_d5_is_d1_on_the_moore_neighborhood_with_a_quiescent_state(self):
        rule = CA.from_tessellation(1, 3, 2, max)
        assert rule.offsets == ((-1,), (0,), (1,))
        assert 2 in rule.quiescent_states

    @given(line_automata())
    def test_quiescent_states_satisfy_f_q_q_equals_q(self, rule):
        n = rule.neighborhood_size
        assert rule.quiescent_states == tuple(
            q for q in range(rule.states) if rule.local_rule(*[q] * n) == q
        )

    @pytest.mark.parametrize("dimension", [1, 2, 3])
    def test_standard_neighborhood_sizes(self, dimension):
        assert len(ca.von_neumann_neighborhood(dimension)) == 2 * dimension + 1
        assert len(ca.moore_neighborhood(dimension)) == 3**dimension
        assert len(ca.radius_neighborhood(dimension, 2)) == 5**dimension
        assert len(ca.one_way_neighborhood(dimension)) == 2**dimension

    @pytest.mark.parametrize("dimension", [1, 2, 3])
    def test_standard_neighborhoods_are_norm_balls(self, dimension):
        box = list(itertools.product(range(-3, 4), repeat=dimension))
        assert set(ca.von_neumann_neighborhood(dimension)) == {
            y for y in box if sum(abs(c) for c in y) <= 1
        }
        assert set(ca.moore_neighborhood(dimension)) == {
            y for y in box if max(abs(c) for c in y) <= 1
        }
        for radius in range(3):
            ball = ca.radius_neighborhood(dimension, radius)
            assert len(ball) == len(set(ball))
            assert set(ball) == {y for y in box if max(abs(c) for c in y) <= radius}
        assert set(ca.one_way_neighborhood(dimension)) == {
            y for y in box if all(c in (0, 1) for c in y)
        }

    def test_one_way_is_the_one_dimensional_radius_half_neighborhood(self):
        assert ca.one_way_neighborhood(1) == ((0,), (1,))

    @given(line_automata())
    def test_any_finite_superset_of_a_memory_set_is_a_memory_set(self, rule):
        extended = rule.with_neighborhood([*rule.offsets, (5,), (-7,)])
        assert extended.neighborhood_size == rule.neighborhood_size + 2
        assert extended.same_global_map(rule)

    @given(line_automata())
    def test_minimal_memory_set_is_contained_in_memory_sets_and_essential(self, rule):
        minimal = rule.minimal()
        assert set(minimal.offsets) <= set(rule.offsets)
        assert minimal.same_global_map(rule)
        assert rule.with_neighborhood([*rule.offsets, (9,)]).minimal() == minimal
        n, k = minimal.neighborhood_size, minimal.states
        for position in range(n):
            assert any(
                minimal.local_rule(*cells)
                != minimal.local_rule(*cells[:position], other, *cells[position + 1 :])
                for cells in itertools.product(range(k), repeat=n)
                for other in range(k)
            )

    @given(line_automata(reach=3), st.data())
    def test_contiguous_form_has_the_same_global_map(self, rule, data):
        config = data.draw(ring_configurations(rule))
        contiguous = rule.to_contiguous()
        positions = [x for (x,) in contiguous.offsets]
        assert positions == list(range(positions[0], positions[-1] + 1))
        assert contiguous.step(config) == rule.step(config)


# --------------------------------------------------------------------------- #
# Structural theorems
# --------------------------------------------------------------------------- #
class TestCurtisHedlundLyndon:
    @given(line_automata(reach=3), st.data())
    def test_step_commutes_with_every_shift_of_a_ring(self, rule, data):
        config = data.draw(ring_configurations(rule))
        by = data.draw(st.integers(-10, 10))
        shape = (len(config),)
        assert rule.step(ca.translate(config, shape, (by,))) == ca.translate(
            rule.step(config), shape, (by,)
        )

    @given(plane_automata(), st.data())
    def test_step_commutes_with_every_translation_of_a_torus(self, rule, data):
        shape = data.draw(st.tuples(st.integers(1, 4), st.integers(1, 4)))
        cells = math.prod(shape)
        config = tuple(
            data.draw(st.lists(st.integers(0, 1), min_size=cells, max_size=cells))
        )
        by = data.draw(st.tuples(st.integers(-5, 5), st.integers(-5, 5)))
        assert rule.step(ca.translate(config, shape, by), shape) == ca.translate(
            rule.step(config, shape), shape, by
        )

    @given(line_automata())
    def test_every_automaton_commutes_with_the_shift_automaton(self, rule):
        shift = ca.shift_automaton(1, rule.states, 0)
        assert rule.compose(shift).same_global_map(shift.compose(rule))

    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_shift_automaton_has_neighborhood_e_i_and_identity_rule(self, axis):
        shift = ca.shift_automaton(3, 4, axis)
        unit = tuple(int(i == axis) for i in range(3))
        assert shift == CA(3, 4, (unit,), (0, 1, 2, 3))

    @given(
        st.integers(0, 1),
        st.lists(st.integers(0, 2), min_size=6, max_size=6),
    )
    def test_shift_automaton_translates_a_torus_along_its_axis(self, axis, cells):
        # sigma_i(c)(x) = c(x + e_i), written out in row-major index arithmetic.
        shift = ca.shift_automaton(2, 3, axis)
        dr, dc = (1, 0) if axis == 0 else (0, 1)
        expected = tuple(
            cells[((r + dr) % 2) * 3 + (c + dc) % 3] for r in range(2) for c in range(3)
        )
        assert shift.step(cells, (2, 3)) == expected
        assert ca.translate(cells, (2, 3), (dr, dc)) == expected

    @given(
        plane_automata(),
        st.sets(st.tuples(st.integers(3, 6), st.integers(3, 6)), max_size=6),
    )
    def test_finite_step_in_the_plane_agrees_with_a_large_torus(self, rule, support):
        if 0 not in rule.quiescent_states:
            return
        torus = [0] * 100
        for r, c in support:
            torus[r * 10 + c] = 1
        image = rule.step(torus, (10, 10))
        expected = {
            (r, c): 1 for r in range(10) for c in range(10) if image[r * 10 + c]
        }
        assert rule.step_finite(dict.fromkeys(support, 1)) == expected

    @given(plane_automata(), st.lists(st.integers(0, 1), min_size=12, max_size=12))
    def test_torus_step_is_the_d1_global_function(self, rule, cells):
        # G(c)(x) = f(c(x + x_1), ..., c(x + x_n)) on the 3 x 4 torus.
        expected = tuple(
            rule.local_rule(
                *(cells[((r + dr) % 3) * 4 + (c + dc) % 4] for dr, dc in rule.offsets)
            )
            for r in range(3)
            for c in range(4)
        )
        assert rule.step(cells, (3, 4)) == expected

    @given(st.lists(st.integers(0, 2), min_size=1, max_size=8))
    def test_shift_automaton_is_the_unit_translation(self, cells):
        shift = ca.shift_automaton(1, 3, 0)
        assert shift.step(cells) == ca.translate(cells, (len(cells),), (1,))


class TestReversibility:
    def test_injective_elementary_rules_are_those_bijective_on_every_small_ring(self):
        # Also the bounded form of Kari 2005, Theorem 7, dimension one:
        # G_P injective implies G injective (a non-injective elementary rule
        # already collides on some ring of at most 12 cells).
        bijective = [
            rule.wolfram_number
            for rule in _all_elementary()
            if all(rule.phase_space(size).is_bijective for size in range(1, 13))
        ]
        assert bijective == list(INJECTIVE)

    def test_no_nontrivial_injective_binary_map_on_two_or_three_cells(self):
        for width in (2, 3):
            offsets = tuple((x,) for x in range(width))
            for table in itertools.product(range(2), repeat=2**width):
                rule = CA(1, 2, offsets, table)
                if rule.is_injective():
                    assert rule.minimal().neighborhood_size == 1

    def test_each_reversible_elementary_rule_has_its_inverse_among_the_six(self):
        rules = ca.reversible_elementary_rules()
        identity = _elementary(204)
        for rule in rules:
            assert any(rule.compose(other).same_global_map(identity) for other in rules)

    def test_compositions_of_reversible_rules_are_reversible(self):
        rules = ca.reversible_elementary_rules()
        for first, second in itertools.product(rules, repeat=2):
            assert first.compose(second).is_reversible()

    @given(st.sampled_from(injective_pool()), st.integers(1, 5))
    def test_injective_implies_bijective_on_rings(self, rule, size):
        assert rule.is_injective()
        assert rule.phase_space(size).is_bijective

    @given(st.sampled_from(injective_pool()))
    def test_injective_implies_pre_injective(self, rule):
        assert rule.is_pre_injective()

    def test_reversibility_is_injectivity_on_the_line(self):
        for rule in _all_elementary():
            assert rule.is_reversible() == rule.is_injective()


class TestGardenOfEden:
    def test_surjective_iff_pre_injective_on_all_elementary_rules(self):
        for rule in _all_elementary():
            assert rule.is_surjective() == rule.is_pre_injective()

    @slow
    @given(line_automata())
    def test_surjective_iff_pre_injective(self, rule):
        assert rule.is_surjective() == rule.is_pre_injective()

    @slow
    @given(right_permutive_automata())
    def test_surjective_iff_pre_injective_on_surjective_rules(self, rule):
        assert rule.is_surjective()
        assert rule.is_pre_injective()

    @slow
    @given(line_automata())
    def test_every_injective_automaton_is_surjective(self, rule):
        assert not rule.is_injective() or rule.is_surjective()

    @given(st.sampled_from(injective_pool()))
    def test_every_injective_automaton_in_the_pool_is_surjective(self, rule):
        assert rule.is_surjective()

    @pytest.mark.parametrize("number", SURJECTIVE)
    def test_surjective_rules_have_no_diamond_up_to_segment_length_8(self, number):
        rule = _elementary(number)
        assert all(rule.diamond(length) is None for length in range(1, 9))

    @pytest.mark.parametrize("number", sorted(set(range(256)) - set(SURJECTIVE)))
    def test_non_surjective_rules_have_a_short_diamond_and_a_short_orphan(self, number):
        rule = _elementary(number)
        diamond = rule.diamond(2)
        assert diamond is not None
        first, second = diamond
        assert first != second
        assert (first[:2], first[-2:]) == (second[:2], second[-2:])
        assert rule.apply_to_word(first) == rule.apply_to_word(second)
        orphans = rule.shortest_orphans()
        assert 1 <= len(orphans[0]) <= 9
        assert (len(orphans[0]) == 9) == (number in LONGEST_ORPHANS)

    @pytest.mark.parametrize("number", range(256))
    def test_surjective_iff_no_orphan_by_brute_force(self, number):
        # An extension of an orphan is an orphan, and the page bounds the
        # shortest elementary orphans by 9: length 9 decides it.
        rule = _elementary(number)
        assert rule.is_surjective() == (_images(rule, 9) == _all_words(rule, 9))
        assert rule.is_surjective() == (rule.shortest_orphans() == ())

    @slow
    @given(line_automata(max_neighbors=2), st.integers(1, 5))
    def test_surjective_rules_have_no_orphan(self, rule, length):
        if rule.is_surjective():
            assert _images(rule, length) == _all_words(rule, length)
            assert not any(rule.is_orphan(w) for w in _all_words(rule, length))

    @pytest.mark.parametrize("number", sorted(set(range(256)) - set(SURJECTIVE)))
    def test_shortest_orphans_match_brute_force(self, number):
        rule = _elementary(number)
        orphans = rule.shortest_orphans()
        length = len(orphans[0])
        assert _images(rule, length - 1) == _all_words(rule, length - 1)
        assert set(orphans) == _all_words(rule, length) - _images(rule, length)
        assert all(rule.is_orphan(orphan) for orphan in orphans)

    @given(line_automata(reach=2), st.data())
    def test_finite_configurations_are_preserved(self, rule, data):
        quiescent = rule.quiescent_states
        if not quiescent:
            return
        q = data.draw(st.sampled_from(quiescent))
        cells = data.draw(
            st.dictionaries(
                st.tuples(st.integers(-6, 6)),
                st.integers(0, rule.states - 1),
                max_size=6,
            )
        )
        image = rule.step_finite(cells, q)
        assert all(state != q for state in image.values())
        ring = [q] * 40
        for (x,), state in cells.items():
            ring[x % 40] = state
        expected = rule.step(ring)
        assert {(x,): s for x, s in enumerate(expected) if s != q} == {
            (x % 40,): s for (x,), s in image.items()
        }


class TestTheorem7:
    @pytest.mark.parametrize("number", SURJECTIVE)
    def test_surjective_line_rule_is_surjective_on_periodic_configurations(
        self, number
    ):
        # Dimension one: G surjective implies G_P surjective. A periodic
        # preimage may need a longer period (rule 90 on odd rings has none of
        # the same period), so periods N, 2N, ... up to 12 cells are searched.
        rule = _elementary(number)
        images = {
            size: frozenset(rule.phase_space(size).successor) for size in range(1, 13)
        }
        for size in range(1, 5):
            for config in itertools.product(range(2), repeat=size):
                assert any(
                    rule.phase_space(size * j).index(config * j) in images[size * j]
                    for j in range(1, 12 // size + 1)
                )

    def test_rule_90_seed_on_3_cells_first_has_a_preimage_at_period_12(self):
        rule = ca.rule_90()
        found = []
        for j in range(1, 5):
            space = rule.phase_space(3 * j)
            found.append(space.index(_seed(3) * j) in set(space.successor))
        assert found == [False, False, False, True]

    @pytest.mark.parametrize("number", sorted(set(range(256)) - set(SURJECTIVE)))
    def test_non_surjective_line_rule_is_not_surjective_on_periodic_configurations(
        self, number
    ):
        # Contrapositive of "G_P surjective implies G surjective", bounded:
        # the periodic configuration spelled by a shortest orphan has no
        # preimage of period L, 2L, ... up to 10 cells.
        rule = _elementary(number)
        orphan = rule.shortest_orphans()[0]
        for j in range(1, max(1, 10 // len(orphan)) + 1):
            space = rule.phase_space(len(orphan) * j)
            assert space.index(orphan * j) not in set(space.successor)

    def test_radius_half_xor_is_the_sum_modulo_2_on_the_radius_half_neighborhood(self):
        rule = ca.radius_half_xor()
        assert rule.offsets == ca.one_way_neighborhood(1) == ((0,), (1,))
        assert rule == CA.from_local_rule(
            1, 2, ca.one_way_neighborhood(1), lambda n: (n[0] + n[1]) % 2
        )

    def test_radius_half_xor_is_injective_on_finite_configurations(self):
        assert ca.radius_half_xor().is_pre_injective()

    def test_radius_half_xor_single_seed_has_no_finite_preimage_in_a_window(self):
        rule = ca.radius_half_xor()
        for support in itertools.product(range(2), repeat=12):
            cells: dict[tuple[int, ...], int] = {
                (x,): 1 for x, state in enumerate(support, start=-6) if state
            }
            assert rule.step_finite(cells) != {(0,): 1}

    @given(st.sets(st.tuples(st.integers(-20, 20)), max_size=8))
    def test_radius_half_xor_images_of_finite_configurations_have_even_weight(
        self, support
    ):
        image = ca.radius_half_xor().step_finite(dict.fromkeys(support, 1))
        assert len(image) % 2 == 0


class TestBalanceTheorem:
    @slow
    @given(right_permutive_automata(), st.data())
    def test_surjective_rule_has_uniform_preimage_counts(self, rule, data):
        word = data.draw(st.lists(st.integers(0, rule.states - 1), max_size=5))
        width = rule.to_contiguous().neighborhood_size
        assert rule.is_surjective()
        assert rule.preimage_count(word) == rule.states ** (width - 1)

    @slow
    @given(line_automata())
    def test_surjective_implies_balanced_table(self, rule):
        assert not rule.is_surjective() or rule.minimal().is_balanced

    def test_balanced_table_does_not_imply_surjective(self):
        assert ca.rule_232().is_balanced
        assert not ca.rule_232().is_surjective()
        assert _elementary(116).is_balanced
        assert not _elementary(116).is_surjective()

    def test_uniform_preimage_count_4_characterizes_the_surjective_rules(self):
        words = list(itertools.product(range(2), repeat=8))
        uniform = [
            rule.wolfram_number
            for rule in _all_elementary()
            if all(rule.preimage_count(word) == 4 for word in words)
        ]
        assert uniform == list(SURJECTIVE)

    @given(line_automata(states=2), st.lists(st.integers(0, 1), max_size=5))
    def test_preimage_count_counts_preimage_words(self, rule, word):
        width = rule.to_contiguous().neighborhood_size
        preimages = [
            w
            for w in itertools.product(range(2), repeat=len(word) + width - 1)
            if rule.apply_to_word(w) == tuple(word)
        ]
        assert rule.preimage_count(word) == len(preimages)

    @pytest.mark.parametrize("number", SURJECTIVE)
    @pytest.mark.parametrize("size", [3, 6, 8])
    def test_surjective_radius_r_rule_has_at_most_s_to_the_2r_preimages(
        self, number, size
    ):
        space = _elementary(number).phase_space(size)
        assert max(space.in_degrees) <= 2 ** (2 * 1)


class TestDecidability:
    @slow
    @given(line_automata())
    def test_pair_graph_injectivity_implies_ring_injectivity(self, rule):
        if rule.is_injective():
            assert all(rule.phase_space(size).is_bijective for size in range(1, 5))

    @slow
    @given(line_automata(states=2), st.integers(1, 6))
    def test_ring_collision_refutes_injectivity(self, rule, size):
        if not rule.phase_space(size).is_bijective:
            assert not rule.is_injective()

    @pytest.mark.parametrize(
        "method",
        ["is_surjective", "is_injective", "is_pre_injective", "is_reversible"],
    )
    def test_line_decisions_are_refused_in_dimension_two(self, method):
        with pytest.raises(ValueError, match="one-dimensional automata only"):
            getattr(ca.game_of_life(), method)()

    @pytest.mark.parametrize("table", [(0, 0), (1, 1), (0, 1), (1, 0)])
    def test_width_one_rules_are_injective_exactly_when_permutations(self, table):
        # Regression: a single de Bruijn vertex hid every difference from the
        # pair graph, so constant maps were reported (pre-)injective.
        rule = CA(1, 2, ((0,),), table)
        permutation = sorted(table) == [0, 1]
        assert rule.is_injective() == permutation
        assert rule.is_pre_injective() == permutation
        assert rule.is_surjective() == permutation

    def test_constant_map_with_empty_neighborhood_is_not_pre_injective(self):
        constant = CA(1, 2, (), (0,))
        assert not constant.is_pre_injective()
        assert not constant.is_injective()
        assert not constant.is_surjective()

    def test_de_bruijn_graph_of_width_m(self):
        edges = ca.rule_110().de_bruijn_graph()
        assert len(edges) == 8
        assert {edge[0] for edge in edges} == set(itertools.product(range(2), repeat=2))
        for source, target, label in edges:
            assert source[1:] == target[:-1]
            assert label == ca.rule_110().local_rule(*source, target[-1])


class TestLinear:
    @settings(max_examples=60, deadline=None)
    @given(linear_automata())
    def test_ito_osato_nasu_surjectivity(self, rule):
        coefficients = rule.linear_coefficients.values()
        assert rule.is_surjective() == (math.gcd(rule.states, *coefficients) == 1)

    @settings(max_examples=60, deadline=None)
    @given(linear_automata())
    def test_ito_osato_nasu_injectivity(self, rule):
        coefficients = list(rule.linear_coefficients.values())
        expected = all(
            sum(c % p != 0 for c in coefficients) == 1
            for p in _prime_factors(rule.states)
        )
        assert rule.is_injective() == expected

    def test_the_eight_linear_elementary_rules(self):
        for number, vector in zip(
            LINEAR_RULES, itertools.product(range(2), repeat=3), strict=True
        ):
            rule = _elementary(number)
            assert rule == CA.from_linear(
                1, 2, dict(zip([(-1,), (0,), (1,)], vector, strict=True))
            )
            assert rule.is_surjective() == (sum(vector) >= 1)
            assert rule.is_injective() == (sum(vector) == 1)

    def test_linear_elementary_rules_are_exactly_the_additive_ones(self):
        additive = [r.wolfram_number for r in _all_elementary() if r.is_additive]
        assert additive == sorted(LINEAR_RULES)

    @given(linear_automata(), st.data())
    def test_additive_rules_obey_superposition(self, rule, data):
        first = data.draw(ring_configurations(rule, max_size=6))
        second = tuple(
            data.draw(
                st.lists(
                    st.integers(0, rule.states - 1),
                    min_size=len(first),
                    max_size=len(first),
                )
            )
        )
        total = tuple((a + b) % rule.states for a, b in zip(first, second, strict=True))
        assert rule.step(total) == tuple(
            (a + b) % rule.states
            for a, b in zip(rule.step(first), rule.step(second), strict=True)
        )

    @slow
    @given(st.integers(2, 4), st.data())
    def test_products_of_polynomials_correspond_to_composition(self, modulus, data):
        coefficient_maps = st.dictionaries(
            st.tuples(st.integers(-1, 1)),
            st.integers(0, modulus - 1),
            min_size=1,
            max_size=3,
        )
        p, q = data.draw(coefficient_maps), data.draw(coefficient_maps)
        product: dict[tuple[int, ...], int] = {}
        for (x,), a in p.items():
            for (y,), b in q.items():
                product[(x + y,)] = (product.get((x + y,), 0) + a * b) % modulus
        composed = CA.from_linear(1, modulus, p).compose(CA.from_linear(1, modulus, q))
        assert composed.linear_coefficients == product

    @given(linear_automata(), st.data())
    def test_on_a_ring_the_polynomial_acts_modulo_x_to_the_n_minus_1(self, rule, data):
        # p(Z) = sum c_i Z^(-x_i) times s(Z) = sum c(x) Z^x represents G(c);
        # on a ring of N cells exponents are read modulo N (MOW, section 2).
        config = data.draw(ring_configurations(rule))
        size, m = len(config), rule.states
        product = [0] * size
        for (offset,), c in rule.linear_coefficients.items():
            for exponent, a in enumerate(config):
                product[(exponent - offset) % size] += c * a
        assert rule.step(config) == tuple(value % m for value in product)

    def test_rule_90_is_the_dipolynomial_x_plus_x_inverse(self):
        assert ca.rule_90().linear_coefficients == {(-1,): 1, (0,): 0, (1,): 1}

    @pytest.mark.parametrize("modulus", [2, 3, 4])
    def test_linear_rule_has_equal_powers_iff_primes_divide_off_center_coefficients(
        self, modulus
    ):
        # Kari 2005, Theorem 18 (equicontinuous iff every prime factor of m
        # divides the off-center coefficients) with "equicontinuous iff
        # G^n = G^m for some n != m". For m <= 4 equal powers, when they
        # exist, already occur among G^0..G^3; their absence there is the
        # bounded half of the check.
        for vector in itertools.product(range(modulus), repeat=3):
            rule = CA.from_linear(
                1, modulus, dict(zip([(-1,), (0,), (1,)], vector, strict=True))
            )
            powers = [rule.power(exponent) for exponent in range(4)]
            repeats = any(
                a.same_global_map(b) for a, b in itertools.combinations(powers, 2)
            )
            expected = all(
                vector[0] % p == 0 and vector[2] % p == 0
                for p in _prime_factors(modulus)
            )
            assert repeats == expected


class TestComposition:
    @given(line_automata(states=2), line_automata(states=2), st.data())
    def test_composition_applies_the_right_factor_first(self, outer, inner, data):
        config = data.draw(ring_configurations(outer))
        composed = outer.compose(inner)
        assert composed.step(config) == outer.step(inner.step(config))

    @given(line_automata(states=3), line_automata(states=3))
    def test_composition_neighborhood_is_the_set_of_sums(self, outer, inner):
        sums = {(x + y,) for (x,) in outer.offsets for (y,) in inner.offsets}
        assert set(outer.compose(inner).offsets) == sums

    @given(
        line_automata(states=2, max_neighbors=2),
        line_automata(states=2, max_neighbors=2),
        line_automata(states=2, max_neighbors=2),
    )
    def test_composition_is_associative_with_identity(self, a, b, c):
        assert a.compose(b).compose(c).same_global_map(a.compose(b.compose(c)))
        assert a.compose(a.power(0)).same_global_map(a)
        assert a.power(0).compose(a).same_global_map(a)

    @given(line_automata(states=2, max_neighbors=2), st.data())
    def test_power_iterates_the_global_map(self, rule, data):
        config = data.draw(ring_configurations(rule))
        assert rule.power(3).step(config) == rule.evolve(config, 3)[-1]

    @given(line_automata())
    def test_equality_of_global_maps_is_decided_up_to_neighborhood(self, rule):
        padded = rule.with_neighborhood([(11,), *rule.offsets])
        assert padded != rule
        assert padded.same_global_map(rule)

    def test_distinct_global_maps_are_told_apart(self):
        assert not ca.rule_90().same_global_map(ca.rule_150())
        assert not _elementary(170).same_global_map(_elementary(240))


class TestSymmetryConjugation:
    @given(st.integers(0, 255))
    def test_reflection_and_complement_formulas(self, number):
        rule = _elementary(number)
        mirrored, complemented = rule.reflect(), rule.complement()
        for a, b, c in itertools.product(range(2), repeat=3):
            assert mirrored.local_rule(a, b, c) == rule.local_rule(c, b, a)
            assert complemented.local_rule(a, b, c) == 1 - rule.local_rule(
                1 - a, 1 - b, 1 - c
            )

    @given(line_automata(states=2, reach=2))
    def test_reflection_and_complement_are_commuting_involutions(self, rule):
        assert rule.reflect().reflect().same_global_map(rule)
        assert rule.complement().complement() == rule
        assert rule.reflect().complement() == rule.complement().reflect()

    @given(line_automata(reach=2), st.data())
    def test_reflection_conjugates_by_the_mirror(self, rule, data):
        config = data.draw(ring_configurations(rule))
        assert rule.reflect().step(config[::-1]) == rule.step(config)[::-1]

    @given(line_automata(states=2, reach=2), st.data())
    def test_complement_conjugates_by_the_state_flip(self, rule, data):
        config = data.draw(ring_configurations(rule))
        flipped = tuple(1 - s for s in config)
        assert rule.complement().step(flipped) == tuple(
            1 - s for s in rule.step(config)
        )


class TestSecondOrder:
    @given(st.integers(0, 255))
    def test_second_order_elementary_rules_are_invertible(self, number):
        assert _elementary(number).second_order().is_injective()

    @given(line_automata(max_neighbors=2), st.data())
    def test_second_order_rule_is_tau_minus_the_past(self, rule, data):
        current = data.draw(ring_configurations(rule, max_size=6))
        size, k = len(current), rule.states
        past = tuple(
            data.draw(st.lists(st.integers(0, k - 1), min_size=size, max_size=size))
        )
        pairs = tuple(k * q + p for q, p in zip(current, past, strict=True))
        image = rule.second_order().step(pairs)
        future = tuple(
            (t - p) % k for t, p in zip(rule.step(current), past, strict=True)
        )
        assert image == tuple(k * f + q for f, q in zip(future, current, strict=True))

    @given(line_automata(max_neighbors=2), st.data())
    def test_second_order_rule_is_undone_by_running_it_on_the_swapped_pair(
        self, rule, data
    ):
        # q^{t-1} = tau q^t - q^{t+1}: the same rule applied to the pair
        # (q^t, q^{t+1}) recovers (q^{t-1}, q^t).
        k = rule.states
        second = rule.second_order()
        size = data.draw(st.integers(1, 6))
        pairs = tuple(
            data.draw(st.lists(st.integers(0, k * k - 1), min_size=size, max_size=size))
        )

        def swap(cells: tuple[int, ...]) -> tuple[int, ...]:
            return tuple(k * (cell % k) + cell // k for cell in cells)

        assert second.step(swap(second.step(pairs))) == swap(pairs)

    def test_second_order_automaton_lives_on_the_alphabet_a_times_a(self):
        assert _elementary(30).second_order().states == 4
        assert ca.one_way_nand().second_order().states == 4


class TestFinitePhaseSpace:
    @given(line_automata(), st.data())
    def test_transient_plus_period_is_at_most_k_to_the_n(self, rule, data):
        config = data.draw(ring_configurations(rule))
        orbit = rule.orbit(config)
        assert orbit.period >= 1
        assert orbit.transient + orbit.period <= rule.states ** len(config)

    @given(line_automata(), st.data())
    def test_orbit_measures_the_first_entry_into_the_cycle(self, rule, data):
        config = data.draw(ring_configurations(rule, max_size=6))
        orbit = rule.orbit(config)
        rows = rule.evolve(config, orbit.transient + 2 * orbit.period)
        assert rows[orbit.transient] == rows[orbit.transient + orbit.period]
        assert len(set(rows[: orbit.transient + orbit.period])) == (
            orbit.transient + orbit.period
        )

    @given(line_automata(), st.integers(1, 5))
    def test_phase_space_is_cycles_with_trees_attached(self, rule, size):
        space = rule.phase_space(size)
        total = rule.states**size
        assert len(space.successor) == total
        assert sum(space.in_degrees) == total
        assert sum(space.cycle_lengths) == len(space.periodic_states)
        assert set(space.fixed_points) <= set(space.periodic_states)
        for index in range(total):
            config = space.configuration(index)
            assert space.index(config) == index
            assert space.successor[index] == space.index(rule.step(config))
            orbit = rule.orbit(config)
            assert (orbit.transient == 0) == (index in space.periodic_states)
            assert orbit.transient <= space.max_transient

    @given(plane_automata())
    def test_torus_phase_space_matches_step(self, rule):
        space = rule.phase_space((2, 2))
        for index in range(16):
            config = space.configuration(index)
            assert space.successor[index] == space.index(rule.step(config, (2, 2)))

    @given(line_automata(), st.data())
    def test_spatially_periodic_configurations_stay_periodic(self, rule, data):
        config = data.draw(ring_configurations(rule, max_size=4))
        assert rule.step(config * 3) == rule.step(config) * 3


class TestRule90OnARing:
    @pytest.mark.parametrize("size", range(3, 13))
    def test_theorem_3_1_fraction_of_gardens_of_eden(self, size):
        space = ca.rule_90().phase_space(size)
        numerator, denominator = (1, 2) if size % 2 else (3, 4)
        assert len(space.gardens_of_eden) * denominator == 2**size * numerator

    @pytest.mark.parametrize("size", range(3, 13))
    def test_theorem_3_2_predecessor_counts(self, size):
        space = ca.rule_90().phase_space(size)
        expected = 2 if size % 2 else 4
        assert set(space.in_degrees) == {0, expected}

    @pytest.mark.parametrize("size", range(3, 13))
    def test_lemma_3_1_odd_weight_configurations_are_never_generated(self, size):
        space = ca.rule_90().phase_space(size)
        gardens = set(space.gardens_of_eden)
        for index in range(2**size):
            if sum(space.configuration(index)) % 2:
                assert index in gardens

    @pytest.mark.parametrize("size", range(3, 13))
    def test_theorems_3_3_and_3_4_trees_and_cycle_fraction(self, size):
        rule = ca.rule_90()
        space = rule.phase_space(size)
        d2 = size & -size
        height = 1 if size % 2 else d2 // 2
        # Balanced: every leaf (Garden of Eden) sits at the full height.
        for index in space.gardens_of_eden:
            assert rule.orbit(space.configuration(index)).transient == height
        assert space.max_transient == height
        if size % 2:
            assert len(space.periodic_states) * 2 == 2**size
            assert all(space.in_degrees[i] == 2 for i in space.periodic_states)
        else:
            assert len(space.periodic_states) * 2**d2 == 2**size
            assert all(space.in_degrees[i] == 4 for i in space.periodic_states)

    @pytest.mark.parametrize("size", range(3, 13))
    def test_lemma_3_4_every_cycle_length_divides_pi_n(self, size):
        pi_n = ca.rule_90().orbit(_seed(size)).period
        space = ca.rule_90().phase_space(size)
        assert all(pi_n % length == 0 for length in space.cycle_lengths)

    @pytest.mark.parametrize("size", [4, 8, 16, 32])
    def test_lemma_3_5_pi_n_is_1_for_powers_of_two(self, size):
        assert ca.rule_90().orbit(_seed(size)).period == 1

    @pytest.mark.parametrize("size", [6, 10, 12, 14, 18, 20, 22, 24, 26, 28, 30])
    def test_lemma_3_6_pi_n_doubles_from_half_the_ring(self, size):
        rule = ca.rule_90()
        assert rule.orbit(_seed(size)).period == 2 * rule.orbit(_seed(size // 2)).period

    @pytest.mark.parametrize("size", range(3, 39, 2))
    def test_theorem_3_5_pi_n_divides_2_to_the_sord_minus_1(self, size):
        sord = next(j for j in range(1, size) if pow(2, j, size) in (1, size - 1))
        pi_n = ca.rule_90().orbit(_seed(size)).period
        assert (2**sord - 1) % pi_n == 0
        assert pi_n <= 2 ** ((size - 1) // 2) - 1
        assert (pi_n == 2**sord - 1) == (size < 37)

    @pytest.mark.parametrize("size", range(3, 16))
    def test_lemma_3_7_fixed_points(self, size):
        space = ca.rule_90().phase_space(size)
        if size % 3 == 0:
            assert len(space.fixed_points) == 4
        else:
            assert space.fixed_points == (0,)

    def test_pi_n_table(self):
        rule = ca.rule_90()
        periods = tuple(rule.orbit(_seed(size)).period for size in range(3, 31))
        assert periods == PI_N

    def test_cycle_census_for_10_cells(self):
        space = ca.rule_90().phase_space(10)
        assert Counter(space.cycle_lengths) == {6: 40, 3: 5, 1: 1}

    def test_pi_13_differs_between_rules_90_and_150(self):
        assert ca.rule_90().orbit(_seed(13)).period == 63
        assert ca.rule_150().orbit(_seed(13)).period == 21


class TestRule90FromASingleSeed:
    def test_rows_are_pascals_triangle_modulo_2(self):
        rows = ca.rule_90().evolve_finite({(0,): 1}, 63)
        for t, row in enumerate(rows):
            expected = {(2 * j - t,): 1 for j in range(t + 1) if math.comb(t, j) % 2}
            assert row == expected

    def test_number_of_ones_is_2_to_the_binary_weight_of_t(self):
        rows = ca.rule_90().evolve_finite({(0,): 1}, 63)
        for t, row in enumerate(rows):
            assert len(row) == 2 ** t.bit_count()


class TestWolframRuleClasses:
    def test_legal_means_null_fixed_and_reflection_symmetric(self):
        for rule in _all_elementary():
            symmetric = all(
                rule.local_rule(a, b, c) == rule.local_rule(c, b, a)
                for a, b, c in itertools.product(range(2), repeat=3)
            )
            assert rule.is_legal == (rule.local_rule(0, 0, 0) == 0 and symmetric)

    def test_totalistic_rules_depend_only_on_the_sum(self):
        assert ca.rule_232().is_totalistic
        assert ca.rule_150().is_totalistic
        assert not ca.rule_110().is_totalistic
        assert sum(rule.is_totalistic for rule in _all_elementary()) == 2**4

    def test_peripheral_rules_ignore_the_cells_own_value(self):
        assert ca.rule_90().is_peripheral
        assert not ca.rule_150().is_peripheral
        assert ca.shift_automaton(1, 2, 0).is_peripheral
        assert sum(rule.is_peripheral for rule in _all_elementary()) == 2**4


# --------------------------------------------------------------------------- #
# Round trips
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(line_automata())
    def test_dataframe_round_trip(self, rule):
        assert CA.from_dataframe(rule.to_dataframe()) == rule

    @given(plane_automata())
    def test_dataframe_round_trip_in_the_plane(self, rule):
        assert CA.from_dataframe(rule.to_dataframe()) == rule

    @given(line_automata(), st.integers(0, 2**16))
    def test_dataframe_rows_may_come_in_any_order(self, rule, seed):
        frame = rule.to_dataframe().sample(frac=1, random_state=seed)
        assert CA.from_dataframe(frame) == rule

    def test_dataframe_is_tidy(self):
        frame = ca.rule_110().to_dataframe()
        assert list(frame.columns) == ["kind", "index", "axis", "value"]
        assert Counter(frame["kind"]) == {
            "dimension": 1, "states": 1, "offset": 3, "rule": 8
        }  # fmt: skip
        rules = frame[frame["kind"] == "rule"]
        assert list(rules["value"]) == list(ca.rule_110().table)

    @pytest.mark.parametrize(
        "example",
        [ca.rule_110, ca.patt_rule, ca.game_of_life, ca.one_way_nand],
    )
    def test_dataframe_survives_write_result(self, example, tmp_path):
        path = io.write_result(example().to_dataframe(), tmp_path / "automaton.json")
        assert CA.from_dataframe(pd.read_json(path, dtype=False)) == example()

    def test_constant_map_with_empty_neighborhood_round_trips(self):
        constant = CA(2, 3, (), (2,))
        assert CA.from_dataframe(constant.to_dataframe()) == constant
        assert constant.step((0, 1, 2, 0), (2, 2)) == (2, 2, 2, 2)

    @given(line_automata())
    def test_minimal_form_is_idempotent(self, rule):
        assert rule.minimal().minimal() == rule.minimal()

    def test_repr_is_compact(self):
        assert repr(ca.rule_110()) == "CellularAutomaton(elementary rule 110)"
        assert repr(ca.patt_rule()) == (
            "CellularAutomaton(d=1, k=2, N=(-1, 0, 1, 2), f=0010110100001111)"
        )
        assert len(repr(ca.game_of_life())) < 200


# --------------------------------------------------------------------------- #
# Rejections
# --------------------------------------------------------------------------- #
class TestRejections:
    def test_d1_dimension_must_be_positive(self):
        with pytest.raises(ValueError, match=r"\(D1\) dimension violated"):
            CA(0, 2, (), (0,))

    def test_d1_state_set_must_be_non_empty(self):
        with pytest.raises(ValueError, match=r"\(D1\) state set violated"):
            CA(1, 0, (), (0,))
        with pytest.raises(ValueError, match=r"\(D1\) state set violated"):
            CA.from_local_rule(1, 0, [(0,)], lambda n: 0)

    def test_d1_neighbors_must_lie_in_z_d(self):
        with pytest.raises(
            ValueError, match=r"\(D1\) neighborhood vector violated.*Z\^2"
        ):
            CA(2, 2, ((0,),), (0, 1))

    def test_d1_neighbors_must_be_distinct(self):
        with pytest.raises(ValueError, match=r"\(D1\) neighborhood vector.*distinct"):
            CA(1, 2, ((0,), (0,)), (0, 1, 1, 0))

    def test_d1_local_rule_must_be_total(self):
        with pytest.raises(ValueError, match=r"\(D1\) local rule violated.*8 table"):
            CA(1, 2, ((-1,), (0,), (1,)), (0, 1, 1, 0))

    def test_d1_local_rule_must_take_values_in_s(self):
        with pytest.raises(ValueError, match=r"\(D1\) local rule violated.*values"):
            CA(1, 2, ((0,),), (0, 2))

    def test_d2_memory_set_must_be_a_set_in_z_d(self):
        with pytest.raises(ValueError, match=r"\(D2\) memory set violated"):
            CA.from_memory_set(1, 2, [(0,), (0,)], lambda pattern: 0)
        with pytest.raises(ValueError, match=r"\(D2\) memory set violated"):
            CA.from_memory_set(1, 2, [(0, 0)], lambda pattern: 0)

    def test_d4_memory_and_anticipation_are_non_negative(self):
        with pytest.raises(ValueError, match=r"\(D4\) violated"):
            CA.from_block_map(2, -1, 1, lambda n: 0)

    def test_d5_quiescent_state_must_be_fixed(self):
        with pytest.raises(ValueError, match=r"\(D5\) quiescent state violated"):
            CA.from_tessellation(1, 2, 0, lambda n: 1)

    def test_finite_dynamics_need_a_quiescent_state(self):
        with pytest.raises(ValueError, match=r"\(D5\) quiescent state violated"):
            ca.one_way_nand().step_finite({(0,): 1})
        with pytest.raises(ValueError, match=r"\(D5\) quiescent state violated"):
            ca.rule_90().evolve_finite({}, 3, quiescent=1)

    def test_wolfram_number_out_of_range_is_rejected(self):
        with pytest.raises(ValueError, match=r"0\.\.255"):
            CA.from_wolfram_number(256)

    def test_wolfram_number_is_defined_for_elementary_automata_only(self):
        with pytest.raises(ValueError, match="elementary automata"):
            _ = ca.patt_rule().wolfram_number

    def test_linear_coefficients_need_an_additive_rule(self):
        with pytest.raises(ValueError, match="additive rule"):
            _ = ca.rule_110().linear_coefficients

    def test_complement_is_defined_for_two_states(self):
        with pytest.raises(ValueError, match="two states"):
            ca.shift_automaton(1, 3, 0).complement()

    def test_reflection_is_defined_on_the_line(self):
        with pytest.raises(ValueError, match="one-dimensional"):
            ca.game_of_life().reflect()

    def test_composition_needs_matching_dimension_and_states(self):
        with pytest.raises(ValueError, match="one dimension and one state set"):
            ca.rule_90().compose(ca.game_of_life())

    def test_negative_power_is_rejected(self):
        with pytest.raises(ValueError, match="non-negative"):
            ca.rule_90().power(-1)

    def test_neighborhood_extension_must_be_a_superset(self):
        with pytest.raises(ValueError, match="must contain the old one"):
            ca.rule_90().with_neighborhood([(0,), (1,)])

    def test_bad_configurations_are_rejected(self):
        with pytest.raises(ValueError, match="a configuration assigns"):
            ca.rule_90().step((0, 2, 0))
        with pytest.raises(ValueError, match="explicit shape"):
            ca.game_of_life().step((0, 0, 0, 0))
        with pytest.raises(ValueError, match="positive periods"):
            ca.game_of_life().step((0, 0, 0, 0), (4,))
        with pytest.raises(ValueError, match="a finite configuration maps"):
            ca.rule_90().step_finite({(0, 0): 1})

    def test_negative_steps_are_rejected(self):
        with pytest.raises(ValueError, match="steps must be non-negative"):
            ca.rule_90().evolve((0, 1), -1)
        with pytest.raises(ValueError, match="steps must be non-negative"):
            ca.rule_90().evolve_finite({}, -1)

    def test_bad_words_are_rejected(self):
        with pytest.raises(ValueError, match="a word lists states"):
            ca.rule_90().apply_to_word((0, 3))
        with pytest.raises(ValueError, match="at least m - 1 = 2 cells"):
            ca.rule_90().apply_to_word((0,))
        with pytest.raises(ValueError, match="non-empty"):
            ca.rule_90().diamond(0)

    def test_bad_phase_space_lookups_are_rejected(self):
        space = ca.rule_90().phase_space(3)
        with pytest.raises(ValueError, match="configuration numbers run over"):
            space.configuration(8)
        with pytest.raises(ValueError, match="a configuration assigns"):
            space.index((0, 1))

    def test_bad_neighborhood_requests_are_rejected(self):
        with pytest.raises(ValueError, match=r"\(D1\) dimension violated"):
            ca.moore_neighborhood(0)
        with pytest.raises(ValueError, match="radius is non-negative"):
            ca.radius_neighborhood(1, -1)
        with pytest.raises(ValueError, match="axes"):
            ca.shift_automaton(2, 2, 2)
        with pytest.raises(ValueError, match="lengths disagree"):
            ca.translate((0, 1), (2,), (1, 1))
        with pytest.raises(ValueError, match="non-negative length"):
            ca.life_row(-1)

    def test_dataframe_missing_column_is_rejected(self):
        frame = ca.rule_110().to_dataframe().drop(columns=["axis"])
        with pytest.raises(ValueError, match="missing columns"):
            CA.from_dataframe(frame)

    def test_dataframe_with_a_gap_is_rejected(self):
        frame = ca.rule_110().to_dataframe()
        with pytest.raises(ValueError, match="gap"):
            CA.from_dataframe(
                frame[~((frame["kind"] == "rule") & (frame["index"] == 5))]
            )

    def test_dataframe_without_scalars_is_rejected(self):
        frame = ca.rule_110().to_dataframe()
        with pytest.raises(ValueError, match="one dimension row"):
            CA.from_dataframe(frame[frame["kind"] != "states"])

    def test_dataframe_with_an_unknown_kind_is_rejected(self):
        frame = ca.rule_110().to_dataframe()
        frame.loc[0, "kind"] = "mystery"
        with pytest.raises(ValueError, match="unknown row kind"):
            CA.from_dataframe(frame)

    def test_dataframe_decoding_revalidates_d1(self):
        frame = ca.rule_110().to_dataframe()
        frame.loc[frame["kind"] == "rule", "value"] = 7
        with pytest.raises(ValueError, match=r"\(D1\) local rule violated"):
            CA.from_dataframe(frame)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
class TestPlots:
    def test_space_time_diagram_has_one_row_per_time_step(self):
        figure, ax = plt.subplots()
        try:
            assert ca.rule_90().plot_space_time(_seed(11), 7, ax) is ax
            (image,) = ax.images
            array = image.get_array()
            assert array is not None
            assert array.shape == (8, 11)
            assert array[1].tolist() == list(ca.rule_90().step(_seed(11)))
        finally:
            plt.close(figure)

    def test_configuration_picture_has_the_torus_shape(self):
        figure, ax = plt.subplots()
        try:
            life = ca.game_of_life()
            assert life.plot_configuration((0, 1, 0, 0, 1, 0), (2, 3), ax) is ax
            (image,) = ax.images
            array = image.get_array()
            assert array is not None
            assert array.tolist() == [[0, 1, 0], [0, 1, 0]]
        finally:
            plt.close(figure)

    def test_de_bruijn_picture_labels_every_edge_and_vertex(self):
        figure, ax = plt.subplots()
        try:
            rule = ca.rule_110()
            assert rule.plot_de_bruijn(ax) is ax
            texts = Counter(text.get_text() for text in ax.texts)
            for a, b, c in itertools.product(range(2), repeat=3):
                assert texts[f"{a}{b}{c} -> {rule.local_rule(a, b, c)}"] == 1
            for vertex in ("00", "01", "10", "11"):
                assert texts[vertex] == 1
        finally:
            plt.close(figure)

    @pytest.mark.parametrize(
        ("rule", "method", "args"),
        [
            (ca.rule_30, "plot_space_time", ((0, 1, 0, 0), 3)),
            (ca.game_of_life, "plot_configuration", ((0, 1, 1, 0), (2, 2))),
            (ca.rule_30, "plot_de_bruijn", ()),
        ],
    )
    def test_plots_create_axes_when_none_are_given(self, rule, method, args):
        ax = getattr(rule(), method)(*args)
        try:
            assert ax.figure is not None
        finally:
            plt.close(ax.figure)

    def test_plots_refuse_the_wrong_dimension(self):
        with pytest.raises(ValueError, match="one-dimensional"):
            ca.game_of_life().plot_space_time((0, 1), 1)
        with pytest.raises(ValueError, match="two-dimensional automaton"):
            ca.rule_90().plot_configuration((0, 1), (2,))
