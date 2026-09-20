"""Tests for research.elementary_cellular_automaton, from its Logseq page.

Fixtures assert exactly what the page's canonical-example blocks certify;
the census, symmetry, permutivity, additive-ring (Martin-Odlyzko-Wolfram),
single-seed, rule 30, rule 110, rule 184/232 and class-table tests transcribe
the structural theorems and their *Oracle* blocks; round-trip laws come from
the API contract; and every formulation (E1)-(E7) has a rejection test naming
it. The family is finite, so most laws are checked over all 256 rules rather
than sampled; Hypothesis covers the ring configurations.

One test is marked PAGE DISCREPANCY: the page states the rule 30 ring
criterion (Wolfram 1986, section 9) with the opposite polarity to what
exhaustive computation gives (spec, "Page discrepancies").

Not transcribed (spec, "Not transcribed" and its amendments): universality,
P-completeness, the small universal Turing machines and the undecidability
consequences; fractal dimensions, the density formula (3.3), propagation
speeds, Jen's two-column theorem and equal block frequencies; Wolfram's 24
nested / 10 random counts (the sets are not listed) and his simple/complex
split with the Fig. 16 list (no decidable definition); the Li-Packard
transition probabilities; "simulates no other rule with blocks up to eight";
rule 18's kink dynamics, growth constant and density, and rule 22's density;
the 16-state bound on the regular-language automaton; conserved quantities on
blocks longer than 1; the null-boundary embedding (MOW 4F); the chaos,
transitivity and dense-periodic-point clauses of Cattaneo-Finelli-Margara
Corollary 3.3 and Schule-Stoop Propositions 16 and 18, and sensitivity,
equicontinuity and expansivity as dynamics (only the published tables are
carried); Cattaneo-Finelli-Margara Theorems 5.1 and 5.2 (outside the family);
rule 110's glider periods; the blocking words, limit sets and attractors of
rules 73, 128, 136 and 232; the dynamics of rule 180 composed with the shift;
the fit log2 Pi_N ~ 0.61 (N + 1); and the undecidability of the Culik-Yu
classes.
"""

import itertools
import math
from collections import Counter, defaultdict
from collections.abc import Callable

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from experiments import io
from research import cellular_automaton as ca
from research import elementary_cellular_automaton as eca

matplotlib.use("Agg")

E = eca.ElementaryCellularAutomaton
ALL = eca.all_rules()

type Bits = tuple[int, ...]
type LocalRule = Callable[[int, int, int], int]
type Tree = tuple[Tree, ...]

# --- oracles, as printed on the page --------------------------------------- #
REPRESENTATIVES = (
    *range(16), 18, 19, *range(22, 31), *range(32, 39), *range(40, 47), 50, 51,
    54, 56, 57, 58, 60, 62, 72, 73, 74, 76, 77, 78, 90, 94, 104, 105, 106, 108,
    110, 122, 126, 128, 130, 132, 134, 136, 138, 140, 142, 146, 150, 152, 154,
    156, 160, 162, 164, 168, 170, 172, 178, 184, 200, 204, 232,
)  # fmt: skip
LEGAL = (
    0, 4, 18, 22, 32, 36, 50, 54, 72, 76, 90, 94, 104, 108, 122, 126, 128, 132,
    146, 150, 160, 164, 178, 182, 200, 204, 218, 222, 232, 236, 250, 254,
)  # fmt: skip
TOTALISTIC = (0, 1, 22, 23, 104, 105, 126, 127, 128, 129, 150, 151, 232, 233, 254, 255)
ADDITIVE = (0, 60, 90, 102, 150, 170, 204, 240)
AFFINE_COMPLEMENTS = (255, 195, 165, 153, 105, 85, 51, 15)
LEFT_PERMUTIVE = (
    15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240,
)  # fmt: skip
RIGHT_PERMUTIVE = (
    85, 86, 89, 90, 101, 102, 105, 106, 149, 150, 153, 154, 165, 166, 169, 170,
)  # fmt: skip
BIPERMUTIVE = (90, 105, 150, 165)
SURJECTIVE = tuple(sorted({*LEFT_PERMUTIVE, *RIGHT_PERMUTIVE, 51, 204}))
SURJECTIVE_CLASSES = (15, 30, 45, 51, 60, 90, 105, 106, 150, 154, 170, 204)
REVERSIBLE = (15, 51, 85, 170, 204, 240)
NUMBER_CONSERVING = (170, 184, 204, 226, 240)
SELF_EQUIVALENT = (23, 51, 77, 105, 150, 178, 204, 232)

KURKA = {
    eca.KurkaClass.EQUICONTINUOUS: (
        0, 1, 4, 5, 8, 12, 19, 29, 36, 51, 72, 76, 108, 200, 204,
    ),
    eca.KurkaClass.ALMOST_EQUICONTINUOUS: (
        13, 23, 28, 32, 33, 40, 44, 50, 73, 77, 78, 94, 104, 128, 132, 136, 140,
        156, 160, 164, 168, 172, 178, 232,
    ),
    eca.KurkaClass.SENSITIVE: (
        2, 3, 6, 7, 9, 10, 11, 14, 15, 18, 22, 24, 25, 26, 27, 30, 34, 35, 37, 38,
        41, 42, 43, 45, 46, 54, 56, 57, 58, 60, 62, 74, 106, 110, 122, 126, 130,
        134, 138, 142, 146, 152, 154, 162, 170, 184,
    ),
    eca.KurkaClass.POSITIVELY_EXPANSIVE: (90, 105, 150),
}  # fmt: skip
WOLFRAM = {
    eca.WolframClass.I: (0, 8, 32, 40, 128, 136, 160, 168),
    eca.WolframClass.III: (18, 22, 30, 45, 60, 90, 105, 122, 126, 146, 150),
    eca.WolframClass.IV: (41, 54, 106, 110),
}
LI_PACKARD = {
    eca.LiPackardClass.NULL: (0, 8, 32, 40, 128, 136, 160, 168),
    eca.LiPackardClass.FIXED_POINT: (
        2, 4, 10, 12, 13, 24, 34, 36, 42, 44, 46, 56, 57, 58, 72, 76, 77, 78, 104,
        130, 132, 138, 140, 152, 162, 164, 170, 172, 184, 200, 204, 232,
    ),
    eca.LiPackardClass.PERIODIC: (
        1, 3, 5, 6, 7, 9, 11, 14, 15, 19, 23, 25, 27, 28, 29, 33, 35, 37, 38, 41,
        43, 50, 51, 62, 74, 94, 108, 134, 142, 156, 178,
    ),
    eca.LiPackardClass.LOCALLY_CHAOTIC: (26, 73, 154),
    eca.LiPackardClass.CHAOTIC: (
        18, 22, 30, 45, 54, 60, 90, 105, 106, 110, 122, 126, 146, 150,
    ),
}  # fmt: skip

# Pi_N for N = 3..30 (MOW Table 1): rule 90 = OEIS A085587, rule 150 = A085588.
PI_90 = (
    1, 1, 3, 2, 7, 1, 7, 6, 31, 4, 63, 14, 15, 1, 15, 14, 511, 12, 63, 62, 2047,
    8, 1023, 126, 511, 28, 16383, 30,
)  # fmt: skip
PI_150 = (
    1, 2, 3, 1, 7, 4, 7, 6, 31, 2, 21, 14, 15, 8, 15, 14, 511, 12, 63, 62, 2047,
    4, 1023, 42, 511, 28, 16383, 30,
)  # fmt: skip
RULE_90_MAX_TRANSIENT = (1, 2, 1, 1, 1, 4, 1, 1, 1, 2)  # N = 3..12
RULE_30_MAX_CYCLE = (
    1, 1, 1, 8, 5, 1, 63, 40, 171, 15, 154, 102, 832, 1428, 1455, 6016,
)  # fmt: skip
RULE_110_MAX_CYCLE = (1, 1, 1, 2, 1, 9, 14, 16, 7, 25, 110, 18, 351, 91, 295, 32)
RULE_30_COLUMN = "1101110011000101100100111010111001110101"  # OEIS A051023
RULE_30_COUNTS = (1, 3, 3, 6, 4, 9, 5, 12, 7, 12, 11, 14, 12, 19, 13, 22)  # A070952
RULE_30_ROWS = (1, 7, 25, 111, 401, 1783, 6409, 28479)  # A110240
RULE_110_COUNTS = (1, 2, 3, 3, 5, 3, 5, 6, 8, 5, 6, 8, 8, 8, 11, 11)  # A071049
RULE_110_ROWS = (1, 3, 7, 13, 31, 49, 115, 215)  # A006978
RULE_150_COUNTS = (1, 3, 3, 5, 3, 9, 5, 11, 3, 9, 9, 15, 5, 15, 11, 21)  # A071053
RULE_22_COUNTS = (1, 3, 2, 6, 2, 6, 4, 12)  # A071044
GOULD = (1, 2, 2, 4, 2, 4, 4, 8)  # A001316
RULE_18_REACHABLE = (1, 1, 4, 7, 11, 19, 36, 67, 121, 216)  # N = 1..10

NORMAL_FORMS = {
    30: ("p", "q", "r", "qr"),
    45: ("1", "p", "r", "qr"),
    54: ("p", "q", "r", "pr"),
    60: ("p", "q"),
    90: ("p", "r"),
    102: ("q", "r"),
    105: ("1", "p", "q", "r"),
    110: ("q", "r", "qr", "pqr"),
    150: ("p", "q", "r"),
    184: ("p", "pq", "qr"),
    232: ("pq", "pr", "qr"),
    22: ("p", "q", "r", "pqr"),
    126: ("p", "q", "r", "pq", "pr", "qr"),
    18: ("p", "r", "pq", "qr"),
}
LOGIC_FORMS = {
    30: lambda p, q, r: p ^ (q | r),
    45: lambda p, q, r: p ^ (q | (1 - r)),
    54: lambda p, q, r: q ^ (p | r),
    110: lambda p, q, r: (q | r) & (1 - (p & q & r)),
    184: lambda p, q, r: p ^ (p & q) ^ (q & r),
    232: lambda p, q, r: (p & q) | ((p | q) & r),
    18: lambda p, q, r: (1 - q) & (p ^ r),
}
EQUIVALENCE_FIXTURES = (
    {30, 86, 135, 149},
    {110, 124, 137, 193},
    {45, 75, 89, 101},
    {60, 102, 153, 195},
    {90, 165},
    {184, 226},
    {54, 147},
    {22, 151},
    {18, 183},
    {73, 109},
    {126, 129},
    {170, 240},
    {15, 85},
    {2, 16, 191, 247},
)
DIPOLYNOMIALS = {
    0: (),
    204: (0,),
    240: (1,),
    170: (-1,),
    60: (1, 0),
    102: (0, -1),
    90: (1, -1),
    150: (1, 0, -1),
}


def _left_shift() -> ca.CellularAutomaton:
    return E(170).to_cellular_automaton()


def _right_shift() -> ca.CellularAutomaton:
    return E(240).to_cellular_automaton()


def _left_shift_squared() -> ca.CellularAutomaton:
    return E(170).compose(E(170))


def _right_shift_squared() -> ca.CellularAutomaton:
    return E(240).compose(E(240))


# Powers sigma^k of the shift, as general automata (rule 170 is the left shift).
SHIFT_POWERS = {
    _left_shift: "sigma",
    _right_shift: "sigma^-1",
    _left_shift_squared: "sigma^2",
    _right_shift_squared: "sigma^-2",
}


# --- families of local rules, by shape ---------------------------------------- #
def _linear(c: Bits) -> LocalRule:
    return lambda p, q, r: (c[0] * p + c[1] * q + c[2] * r) % 2


def _of_the_sum(g: Bits) -> LocalRule:
    return lambda p, q, r: g[p + q + r]


def _p_xor(g: Bits) -> LocalRule:
    return lambda p, q, r: p ^ g[2 * q + r]


def _xor_r(g: Bits) -> LocalRule:
    return lambda p, q, r: g[2 * p + q] ^ r


def _of_the_center(g: Bits) -> LocalRule:
    return lambda p, q, r: g[q]


# --- independent reference implementations --------------------------------- #
def word(text: str) -> Bits:
    return tuple(int(digit) for digit in text)


def kari_rule(number: int, p: int, q: int, r: int) -> int:
    """f(p, q, r) read off Kari's bit sequence f(111) f(110) ... f(000)."""
    return int(f"{number:08b}"[7 - (4 * p + 2 * q + r)])


def reference_step(number: int, cells: Bits) -> Bits:
    """(E7) cell by cell: G(x)_i = f(x_(i-1), x_i, x_(i+1)), indices modulo N."""
    n = len(cells)
    return tuple(
        kari_rule(number, cells[(i - 1) % n], cells[i], cells[(i + 1) % n])
        for i in range(n)
    )


def rotate(config: int, cells: int, k: int) -> int:
    """Move every cell k places to the left, cyclically."""
    k %= cells
    full = (1 << cells) - 1
    return (config << k | config >> (cells - k)) & full


def mirror(config: int, cells: int) -> int:
    return int(format(config, f"0{cells}b")[::-1], 2)


def cyclic_runs_of_ones(cells: Bits) -> list[int]:
    """The lengths of the runs of 1s between cyclically consecutive 0s."""
    n = len(cells)
    zeros = [i for i, value in enumerate(cells) if value == 0]
    return [b - a - 1 for a, b in zip(zeros, [*zeros[1:], zeros[0] + n], strict=True)]


def poly_mod(a: int, b: int) -> int:
    """Remainder over GF(2); polynomials are ints, bit i the coefficient of x^i."""
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


def lambda_1(rule: eca.ElementaryCellularAutomaton, cells: int) -> int:
    """MOW Lemma 4.4: gcd(x^N - 1, T(x)), with T cleared of x^-1 by the unit x."""
    numerator = sum(1 << (exponent + 1) for exponent in rule.dipolynomial)
    return poly_gcd((1 << cells) | 1, numerator)


def polynomial(config: int, cells: int) -> int:
    """A(x) = sum a_i x^i for the ring configuration (cell i is digit i)."""
    return mirror(config, cells)


def transient_trees(space: ca.PhaseSpace) -> dict[int, Tree]:
    """The canonical (sorted nested tuple) form of the tree on every cycle node."""
    periodic = set(space.periodic_states)
    children = defaultdict(list)
    for state, image in enumerate(space.successor):
        if state not in periodic:
            children[image].append(state)

    def canonical(node: int) -> Tree:
        return tuple(sorted(canonical(child) for child in children[node]))

    return {node: canonical(node) for node in periodic}


def unreachable_under_rule_18(cells: Bits) -> bool:
    """MOW Lemma 5.1, clauses (a), (b), (c), on a ring."""
    n = len(cells)
    blocks = [i for i in range(n) if cells[i] == cells[(i + 1) % n] == 1]
    if any(cells[(i + 2) % n] == 1 for i in blocks):
        return True
    if not blocks:
        return sum(cells) % 2 == 1
    for i, j in itertools.product(blocks, repeat=2):
        length = n - 2 if i == j else (j - i - 2) % n
        between = [cells[(i + 2 + k) % n] for k in range(length)]
        if sum(between) % 2 == 1:
            return True
    return False


# --- strategies -------------------------------------------------------------- #
rules = st.integers(0, 255).map(E)
additive_rules = st.sampled_from(ADDITIVE).map(E)


@st.composite
def ring_configs(draw, min_cells=1, max_cells=16):
    cells = draw(st.integers(min_cells, max_cells))
    return draw(st.integers(0, 2**cells - 1)), cells


# --------------------------------------------------------------------------- #
# Definition
# --------------------------------------------------------------------------- #
class TestDefinition:
    def test_rule_110_fixture(self):
        rule = eca.rule_110()
        assert rule.bits == "01101110"
        ones = {"110", "101", "011", "010", "001"}
        zeros = {"111", "100", "000"}
        assert {n for n in ones | zeros if rule.local_rule(*word(n))} == ones

    def test_wolfram_number_closed_form(self):
        assert all(
            rule.rule
            == sum(
                rule.local_rule(p, q, r) << (4 * p + 2 * q + r)
                for p, q, r in itertools.product(range(2), repeat=3)
            )
            for rule in ALL
        )

    def test_local_rule_is_kari_bit_sequence(self):
        assert all(
            rule.local_rule(p, q, r) == kari_rule(rule.rule, p, q, r)
            for rule in ALL
            for p, q, r in itertools.product(range(2), repeat=3)
        )

    def test_there_are_256_rules_and_the_numbering_is_a_bijection(self):
        assert len({rule.table for rule in ALL}) == 256
        assert [rule.rule for rule in ALL] == list(range(256))

    def test_table_and_bits_are_the_same_eight_bits(self):
        assert all(
            rule.table[index] == int(rule.bits[7 - index])
            for rule in ALL
            for index in range(8)
        )

    def test_the_256_normal_forms_are_pairwise_distinct(self):
        assert len({rule.anf for rule in ALL}) == 256

    @pytest.mark.parametrize(("number", "monomials"), NORMAL_FORMS.items())
    def test_normal_form_fixture(self, number, monomials):
        assert E(number).anf == monomials

    @pytest.mark.parametrize(("number", "form"), LOGIC_FORMS.items())
    def test_logic_form_fixture(self, number, form):
        assert E.from_local_rule(form) == E(number)

    def test_normal_form_evaluates_to_the_rule(self):
        values: dict[str, LocalRule] = {
            "1": lambda p, q, r: 1,
            "p": lambda p, q, r: p,
            "q": lambda p, q, r: q,
            "r": lambda p, q, r: r,
            "pq": lambda p, q, r: p & q,
            "pr": lambda p, q, r: p & r,
            "qr": lambda p, q, r: q & r,
            "pqr": lambda p, q, r: p & q & r,
        }
        assert all(
            rule.local_rule(p, q, r) == sum(values[m](p, q, r) for m in rule.anf) % 2
            for rule in ALL
            for p, q, r in itertools.product(range(2), repeat=3)
        )

    def test_nks_logic_form_for_rule_110_is_a_misprint_for_rule_124(self):
        printed = E.from_local_rule(lambda p, q, r: (p | q) ^ (p & q & r))
        assert printed == E(124)
        assert printed == eca.rule_110().reflect()

    def test_reprinted_logic_form_for_rule_18_is_a_misprint_for_rule_222(self):
        assert E.from_local_rule(lambda p, q, r: q | (p ^ r)) == E(222)

    def test_rules_are_adjacent_in_the_cube_at_hamming_distance_one(self):
        assert all(
            [rule.hamming_distance(other) for other in rule.neighbors()] == [1] * 8
            and len(set(rule.neighbors())) == 8
            for rule in ALL
        )

    def test_hamming_distance_counts_differing_table_entries(self):
        assert E(110).hamming_distance(E(102)) == 1
        assert E(0).hamming_distance(E(255)) == 8

    def test_de_bruijn_graph_has_the_eight_edges_pq_to_qr(self):
        rule = eca.rule_110()
        assert set(rule.de_bruijn_graph()) == {
            ((p, q), (q, r), rule.local_rule(p, q, r))
            for p, q, r in itertools.product(range(2), repeat=3)
        }

    @pytest.mark.parametrize("shift", SHIFT_POWERS, ids=SHIFT_POWERS.get)
    @pytest.mark.parametrize("number", [30, 110])
    def test_de_bruijn_graph_forgets_powers_of_the_shift(self, number, shift):
        rule = E(number)
        shifted = rule.to_cellular_automaton().compose(shift())
        assert shifted.minimal().de_bruijn_graph() == rule.de_bruijn_graph()

    def test_minimal_memory_set_fixtures(self):
        assert E(170).minimal_memory_set == {1}
        assert E(204).minimal_memory_set == {0}
        assert E(0).minimal_memory_set == frozenset()
        assert E(110).minimal_memory_set == {-1, 0, 1}

    @given(rules, rules)
    def test_distinct_tables_give_distinct_global_maps(self, first, second):
        same = first.to_cellular_automaton().same_global_map(
            second.to_cellular_automaton()
        )
        assert same == (first == second)

    def test_elementary_automata_are_those_with_memory_set_in_the_neighborhood(self):
        assert E.from_cellular_automaton(ca.radius_half_xor()) == E(102)
        assert E.from_cellular_automaton(ca.shift_automaton(1, 2, 0)) == E(170)

    @given(rules, ring_configs(max_cells=12))
    def test_ring_step_is_the_local_rule_with_indices_modulo_n(self, rule, ring):
        config, cells = ring
        expected = reference_step(rule.rule, eca.cells_from_config(config, cells))
        assert eca.cells_from_config(rule.step(config, cells), cells) == expected

    @given(rules, ring_configs(max_cells=8))
    def test_ring_is_the_restriction_to_spatially_periodic_configurations(
        self, rule, ring
    ):
        config, cells = ring
        doubled = config << cells | config
        image = rule.step(config, cells)
        assert rule.step(doubled, 2 * cells) == image << cells | image

    @pytest.mark.parametrize("cells", [3, 4, 5])
    def test_distinct_rules_are_distinct_maps_on_rings_of_three_or_more(self, cells):
        assert len({rule.phase_space(cells).successor for rule in ALL}) == 256

    @pytest.mark.parametrize(("cells", "occurring"), [(1, (0, 7)), (2, (0, 2, 5, 7))])
    def test_rules_coincide_on_tiny_rings_when_the_occurring_bits_agree(
        self, cells, occurring
    ):
        maps = defaultdict(set)
        for rule in ALL:
            key = tuple(rule.table[index] for index in occurring)
            maps[key].add(rule.phase_space(cells).successor)
        assert all(len(successors) == 1 for successors in maps.values())
        assert len(set().union(*maps.values())) == 2 ** len(occurring)


# --------------------------------------------------------------------------- #
# Derived vocabulary
# --------------------------------------------------------------------------- #
class TestVocabulary:
    def test_legal_rules_are_the_32_of_wolfram_fig_3(self):
        assert tuple(rule.rule for rule in ALL if rule.is_legal) == LEGAL

    def test_legal_means_quiescent_zero_and_reflection_symmetric(self):
        assert all(
            rule.is_legal == (rule.local_rule(0, 0, 0) == 0 and rule.reflect() == rule)
            for rule in ALL
        )

    def test_additive_rules_are_the_eight_of_nks_p_952(self):
        assert tuple(rule.rule for rule in ALL if rule.is_additive) == ADDITIVE

    def test_additive_means_linear_modulo_two(self):
        linear = {
            E.from_local_rule(_linear(c)) for c in itertools.product(range(2), repeat=3)
        }
        assert linear == {rule for rule in ALL if rule.is_additive}

    def test_only_0_90_150_204_are_additive_among_legal_rules(self):
        assert [r.rule for r in ALL if r.is_legal and r.is_additive] == [
            0, 90, 150, 204,
        ]  # fmt: skip

    def test_affine_rules_are_the_additive_rules_and_their_complements(self):
        assert {rule.rule for rule in ALL if rule.is_affine} == {
            *ADDITIVE,
            *AFFINE_COMPLEMENTS,
        }
        assert sorted(255 - number for number in ADDITIVE) == sorted(AFFINE_COMPLEMENTS)

    def test_rule_105_is_affine_not_additive(self):
        rule = E(105)
        assert rule.anf == ("1", "p", "q", "r")
        assert rule.is_affine
        assert not rule.is_additive

    def test_totalistic_rules_are_the_16_on_the_page(self):
        assert tuple(rule.rule for rule in ALL if rule.is_totalistic) == TOTALISTIC

    def test_totalistic_means_depending_on_the_sum_only(self):
        by_sum = {
            E.from_local_rule(_of_the_sum(g))
            for g in itertools.product(range(2), repeat=4)
        }
        assert by_sum == {rule for rule in ALL if rule.is_totalistic}

    def test_legal_peripheral_rules_are_0_90_160_250(self):
        assert [r.rule for r in ALL if r.is_legal and r.is_peripheral] == [
            0, 90, 160, 250,
        ]  # fmt: skip

    def test_permutive_rules_are_the_lists_on_the_page(self):
        assert tuple(r.rule for r in ALL if r.is_left_permutive) == LEFT_PERMUTIVE
        assert tuple(r.rule for r in ALL if r.is_right_permutive) == RIGHT_PERMUTIVE
        assert tuple(r.rule for r in ALL if r.is_bipermutive) == BIPERMUTIVE

    def test_left_permutive_means_p_xor_g(self):
        shaped = {
            E.from_local_rule(_p_xor(g)) for g in itertools.product(range(2), repeat=4)
        }
        assert shaped == {rule for rule in ALL if rule.is_left_permutive}

    def test_right_permutive_means_g_xor_r(self):
        shaped = {
            E.from_local_rule(_xor_r(g)) for g in itertools.product(range(2), repeat=4)
        }
        assert shaped == {rule for rule in ALL if rule.is_right_permutive}

    def test_rule_102_is_not_bipermutive_on_the_three_cell_table(self):
        rule = eca.rule_102()
        assert rule.is_right_permutive
        assert not rule.is_left_permutive
        assert not rule.is_bipermutive

    def test_trivial_rules_depend_on_the_center_cell_only(self):
        center_only = {
            E.from_local_rule(_of_the_center(g))
            for g in itertools.product(range(2), repeat=2)
        }
        assert center_only == {rule for rule in ALL if rule.is_trivial}

    def test_balanced_means_four_ones_and_density_one_half(self):
        assert all(
            rule.is_balanced == (sum(rule.table) == 4) == (2 * rule.density == 1)
            for rule in ALL
        )
        assert all(8 * rule.density == sum(rule.table) for rule in ALL)

    def test_hot_bits_are_the_images_of_000_and_111(self):
        assert all(
            rule.hot_bits == (rule.local_rule(0, 0, 0), rule.local_rule(1, 1, 1))
            for rule in ALL
        )

    def test_there_are_64_mean_field_clusters_and_36_up_to_conjugation(self):
        clusters = {rule.mean_field_cluster for rule in ALL}
        paired = {
            frozenset({rule.mean_field_cluster, rule.conjugate().mean_field_cluster})
            for rule in ALL
        }
        assert len(clusters) == 64
        assert len(paired) == 36

    def test_mean_field_cluster_counts_blocks_by_number_of_ones(self):
        assert eca.rule_110().mean_field_cluster == (0, 2, 3, 0)
        assert eca.rule_232().mean_field_cluster == (0, 0, 3, 1)

    def test_two_adic_part(self):
        assert [eca.two_adic_part(n) for n in (1, 6, 8, 12, 40)] == [1, 2, 8, 4, 8]

    def test_multiplicative_suborder(self):
        assert [eca.multiplicative_suborder(n) for n in (3, 5, 7, 9, 37)] == [
            1, 2, 3, 3, 18,
        ]  # fmt: skip
        assert all(
            pow(2, eca.multiplicative_suborder(n), n) in {1 % n, n - 1}
            for n in range(1, 60, 2)
        )

    def test_space_time_rows_are_consecutive_configurations(self):
        rule = eca.rule_110()
        rows = rule.evolve(0b0001000, 7, 5)
        assert rows[0] == 0b0001000
        assert all(rows[t + 1] == rule.step(rows[t], 7) for t in range(5))

    def test_center_column_is_cell_zero_of_the_single_seed_pattern(self):
        rows = eca.rule_30().evolve_from_single_cell(20)
        assert eca.rule_30().center_column(20) == tuple(
            row[t] for t, row in enumerate(rows)
        )


# --------------------------------------------------------------------------- #
# The symmetry group
# --------------------------------------------------------------------------- #
class TestSymmetry:
    def test_reflection_and_conjugation_by_definition(self):
        assert all(
            rule.reflect().local_rule(p, q, r) == rule.local_rule(r, q, p)
            and rule.conjugate().local_rule(p, q, r)
            == 1 - rule.local_rule(1 - p, 1 - q, 1 - r)
            for rule in ALL
            for p, q, r in itertools.product(range(2), repeat=3)
        )

    def test_reflection_and_conjugation_are_commuting_involutions(self):
        assert all(
            rule.reflect().reflect() == rule
            and rule.conjugate().conjugate() == rule
            and rule.reflect().conjugate() == rule.conjugate().reflect()
            for rule in ALL
        )

    def test_table_permutations_of_li_packard(self):
        def t(rule: eca.ElementaryCellularAutomaton, index: int) -> int:
            return rule.table[index]

        def printed(rule: eca.ElementaryCellularAutomaton) -> Bits:
            return tuple(int(digit) for digit in rule.bits)

        assert all(
            printed(f.reflect()) == tuple(t(f, i) for i in (7, 3, 5, 1, 6, 2, 4, 0))
            and printed(f.conjugate()) == tuple(1 - t(f, i) for i in range(8))
            and printed(f.reflect().conjugate())
            == tuple(1 - t(f, i) for i in (0, 4, 2, 6, 1, 5, 3, 7))
            for f in ALL
        )

    def test_there_are_88_classes_of_sizes_8_36_44(self):
        """Li-Packard 1990, Appendix."""
        classes = {rule.equivalence_class() for rule in ALL}
        assert len(classes) == 88
        assert Counter(len(members) for members in classes) == {1: 8, 2: 36, 4: 44}
        assert 8 + (8 + 8 + 56) // 2 + 176 // 4 == 88

    def test_burnside_fixed_point_counts(self):
        """Burnside's lemma, computed on the page beside Li-Packard 1990, Appendix."""
        images: tuple[Callable[[eca.ElementaryCellularAutomaton], object], ...] = (
            lambda f: f,
            lambda f: f.reflect(),
            lambda f: f.conjugate(),
            lambda f: f.reflect().conjugate(),
        )
        fixed = [sum(image(rule) == rule for rule in ALL) for image in images]
        assert fixed == [256, 64, 16, 16]
        assert sum(fixed) // 4 == 88

    def test_the_eight_self_equivalent_rules(self):
        assert (
            tuple(rule.rule for rule in ALL if len(rule.equivalence_class()) == 1)
            == SELF_EQUIVALENT
        )

    def test_minimal_representatives_are_the_88_on_the_page(self):
        assert tuple(rule.rule for rule in eca.representatives()) == REPRESENTATIVES
        assert all(
            rule.canonical().rule == min(m.rule for m in rule.equivalence_class())
            for rule in ALL
        )

    def test_the_two_representative_conventions_differ_on_five_classes(self):
        differing = {
            rule.rule: rule.li_packard_representative().rule
            for rule in eca.representatives()
            if rule.li_packard_representative() != rule
        }
        assert differing == {62: 131, 94: 133, 110: 137, 122: 161, 126: 129}

    @pytest.mark.parametrize("members", EQUIVALENCE_FIXTURES, ids=str)
    def test_equivalence_class_fixture(self, members):
        assert {rule.rule for rule in E(min(members)).equivalence_class()} == members

    def test_named_images(self):
        assert eca.rule_30().reflect() == E(86)
        assert eca.rule_30().conjugate() == E(135)
        assert eca.rule_184().reflect() == E(226)
        assert eca.rule_184().conjugate() == E(226)

    def test_class_of_rule_138_contains_244_not_224(self):
        assert {rule.rule for rule in E(138).equivalence_class()} == {
            138, 174, 208, 244,
        }  # fmt: skip

    @given(rules, ring_configs())
    def test_reflected_rule_draws_the_mirror_image(self, rule, ring):
        config, cells = ring
        assert rule.reflect().step(mirror(config, cells), cells) == mirror(
            rule.step(config, cells), cells
        )

    @given(rules, ring_configs())
    def test_conjugate_rule_draws_the_black_white_exchange(self, rule, ring):
        config, cells = ring
        full = (1 << cells) - 1
        assert rule.conjugate().step(config ^ full, cells) == (
            rule.step(config, cells) ^ full
        )


# --------------------------------------------------------------------------- #
# Census of the 256 rules
# --------------------------------------------------------------------------- #
class TestCensus:
    def test_census(self):
        counts = {
            "rules": len(ALL),
            "classes": len(eca.representatives()),
            "legal": sum(rule.is_legal for rule in ALL),
            "totalistic": sum(rule.is_totalistic for rule in ALL),
            "additive": sum(rule.is_additive for rule in ALL),
            "surjective": sum(rule.is_surjective for rule in ALL),
            "reversible": sum(rule.is_reversible for rule in ALL),
            "number-conserving": sum(rule.is_number_conserving for rule in ALL),
        }
        assert counts == {
            "rules": 256,
            "classes": 88,
            "legal": 32,
            "totalistic": 16,
            "additive": 8,
            "surjective": 30,
            "reversible": 6,
            "number-conserving": 5,
        }


# --------------------------------------------------------------------------- #
# Superposition and shift commutation
# --------------------------------------------------------------------------- #
def _superposes(
    rule: eca.ElementaryCellularAutomaton,
    cells: int,
    combine: Callable[[int, int], int],
) -> bool:
    return all(
        rule.step(combine(x, y), cells)
        == combine(rule.step(x, cells), rule.step(y, cells))
        for x in range(1 << cells)
        for y in range(x)
    )


class TestSuperposition:
    def test_xor_superposition_holds_for_exactly_the_additive_rules(self):
        superposing = tuple(
            rule.rule
            for rule in ALL
            if all(_superposes(rule, cells, int.__xor__) for cells in (3, 4, 5))
        )
        assert superposing == ADDITIVE

    @given(additive_rules, ring_configs(), st.data())
    def test_xor_superposition_on_every_ring(self, rule, ring, data):
        x, cells = ring
        y = data.draw(st.integers(0, 2**cells - 1))
        assert rule.step(x ^ y, cells) == rule.step(x, cells) ^ rule.step(y, cells)

    def test_or_superposition_among_legal_rules_on_a_ring_of_five(self):
        assert [
            number for number in LEGAL if _superposes(E(number), 5, int.__or__)
        ] == [0, 204, 250, 254]

    def test_and_superposition_among_legal_rules_on_a_ring_of_five(self):
        assert [
            number for number in LEGAL if _superposes(E(number), 5, int.__and__)
        ] == [0, 128, 160, 204]


class TestShiftCommutation:
    @given(rules, ring_configs(), st.integers(-40, 40))
    def test_curtis_hedlund_lyndon_on_rings(self, rule, ring, k):
        config, cells = ring
        assert rule.step(rotate(config, cells, k), cells) == rotate(
            rule.step(config, cells), cells, k
        )


# --------------------------------------------------------------------------- #
# Surjectivity = permutivity
# --------------------------------------------------------------------------- #
class TestSurjectivity:
    def test_surjective_rules_are_the_30_on_the_page(self):
        assert tuple(rule.rule for rule in ALL if rule.is_surjective) == SURJECTIVE
        assert len(SURJECTIVE) == 30

    def test_corollary_3_3_permutive_iff_surjective_and_non_trivial(self):
        """Cattaneo-Finelli-Margara 2000, Corollary 3.3, (1) iff (4)."""
        assert all(
            rule.is_permutive == (rule.is_surjective and not rule.is_trivial)
            for rule in ALL
        )

    def test_rules_51_and_204_are_the_only_surjective_non_permutive_rules(self):
        assert [
            rule.rule for rule in ALL if rule.is_surjective and not rule.is_permutive
        ] == [51, 204]

    def test_surjective_rules_form_the_12_classes_of_schule_stoop(self):
        assert (
            tuple(sorted({E(number).canonical().rule for number in SURJECTIVE}))
            == SURJECTIVE_CLASSES
        )

    def test_devaney_chaotic_classes_are_the_ten_other_than_51_and_204(self):
        """Schule-Stoop 2012, Corollary 19."""
        chaotic = {rule.canonical().rule for rule in ALL if rule.is_devaney_chaotic}
        assert chaotic == set(SURJECTIVE_CLASSES) - {51, 204}

    def test_uniform_four_to_one_preimages_characterize_the_30(self):
        """Cattaneo-Finelli-Margara 2000, Theorem 3.4 (Hedlund 1969)."""
        words = list(itertools.product(range(2), repeat=8))
        uniform = tuple(
            rule.rule
            for rule in ALL
            if all(rule.preimage_count(image) == 4 for image in words)
        )
        assert uniform == SURJECTIVE

    @given(st.sampled_from(SURJECTIVE).map(E), st.lists(st.integers(0, 1), min_size=1))
    def test_hedlund_every_word_of_every_length_has_four_preimages(self, rule, image):
        """Cattaneo-Finelli-Margara 2000, Theorem 3.4 (Hedlund 1969)."""
        assert rule.preimage_count(image) == 4

    def test_subset_construction_agrees_with_brute_force_on_words_of_length_10(self):
        def covers_every_word(rule: eca.ElementaryCellularAutomaton) -> bool:
            images = {rule.step(w, 12) >> 1 & 0x3FF for w in range(1 << 12)}
            return len(images) == 1 << 10

        assert all(rule.is_surjective == covers_every_word(rule) for rule in ALL)

    def test_brute_force_word_images_match_apply_to_word(self):
        rule = eca.rule_110()
        assert all(
            eca.cells_from_config(rule.step(w, 12) >> 1 & 0x3FF, 10)
            == rule.apply_to_word(eca.cells_from_config(w, 12))
            for w in range(0, 1 << 12, 37)
        )

    def test_garden_of_eden_theorem_surjective_iff_no_diamond(self):
        """Moore 1962 and Myhill 1963 (Kari 2005, Theorem 6).

        Pre-injectivity is decided exactly on the pair graph; the diamond
        search is bounded to differing segments of length 1..6, which already
        finds a diamond in every non-surjective elementary rule.
        """
        assert all(
            rule.is_surjective
            == rule.to_cellular_automaton().is_pre_injective()
            == all(rule.diamond(k) is None for k in range(1, 7))
            for rule in ALL
        )

    def test_nks_surjective_iff_additive_in_the_first_or_last_dependent_position(self):
        """Wolfram 2002, pp. 959-960: "precisely for those 30 rules"."""

        def flips_output(rule: eca.ElementaryCellularAutomaton, offset: int) -> bool:
            bit = {-1: 4, 0: 2, 1: 1}[offset]
            return all(
                rule.table[index] != rule.table[index | bit]
                for index in range(8)
                if not index & bit
            )

        def additive_at_an_end(rule: eca.ElementaryCellularAutomaton) -> bool:
            depends = sorted(rule.minimal_memory_set)
            return any(
                flips_output(rule, offset) for offset in depends[:1] + depends[-1:]
            )

        assert (
            tuple(rule.rule for rule in ALL if additive_at_an_end(rule)) == SURJECTIVE
        )

    def test_rule_232_is_balanced_yet_00_has_six_preimages(self):
        """Kari's lecture notes, Example 6."""
        rule = eca.rule_232()
        assert rule.is_balanced
        assert rule.preimage_count(word("00")) == 6
        assert not rule.is_surjective

    def test_rule_184_is_balanced_yet_1100_is_an_orphan(self):
        rule = eca.rule_184()
        assert rule.is_balanced
        assert rule.is_orphan(word("1100"))
        assert not rule.is_surjective

    def test_the_40_balanced_non_trivial_non_permutive_rules_fail_by_length_six(self):
        suspects = [
            rule
            for rule in ALL
            if rule.is_balanced and not rule.is_trivial and not rule.is_permutive
        ]
        short = [
            image
            for length in range(1, 7)
            for image in itertools.product(range(2), repeat=length)
        ]
        assert len(suspects) == 40
        assert all(
            any(rule.preimage_count(image) != 4 for image in short) for rule in suspects
        )

    @pytest.mark.parametrize("number", [43, 113, 142, 212])
    def test_word_010_has_three_preimages(self, number):
        """Cattaneo-Finelli-Margara 2000, proof of Theorem 3.5."""
        assert E(number).preimage_count(word("010")) == 3

    @pytest.mark.parametrize("shift", SHIFT_POWERS, ids=SHIFT_POWERS.get)
    def test_composition_with_a_power_of_the_shift_keeps_surjectivity_and_injectivity(
        self, shift
    ):
        assert all(
            (shifted.is_surjective(), shifted.is_injective())
            == (rule.is_surjective, rule.is_injective)
            for rule in ALL
            for shifted in [rule.to_cellular_automaton().compose(shift())]
        )


# --------------------------------------------------------------------------- #
# Reversible rules and finite rings
# --------------------------------------------------------------------------- #
class TestReversible:
    def test_reversible_rules_are_the_six_of_nks_p_436(self):
        assert tuple(rule.rule for rule in ALL if rule.is_reversible) == REVERSIBLE
        assert tuple(rule.rule for rule in ALL if rule.is_injective) == REVERSIBLE

    def test_the_six_are_the_identity_the_shifts_and_their_complements(self):
        named = {
            15: lambda p, q, r: 1 - p,
            51: lambda p, q, r: 1 - q,
            85: lambda p, q, r: 1 - r,
            170: lambda p, q, r: r,
            204: lambda p, q, r: q,
            240: lambda p, q, r: p,
        }
        assert {n: E.from_local_rule(f).rule for n, f in named.items()} == {
            n: n for n in REVERSIBLE
        }

    def test_injective_rules_are_surjective(self):
        """Kari 2005, Corollary 3."""
        assert all(rule.is_surjective for rule in ALL if rule.is_injective)

    def test_rules_bijective_on_every_ring_up_to_12_are_exactly_the_six(self):
        """Kari 2005, Theorem 7: rings are a sound test of injectivity."""
        bijective = tuple(
            rule.rule
            for rule in ALL
            if all(rule.is_bijective_on_ring(cells) for cells in range(1, 13))
        )
        assert bijective == REVERSIBLE

    @given(rules, st.integers(1, 9))
    def test_is_bijective_on_ring_agrees_with_the_phase_space(self, rule, cells):
        assert rule.is_bijective_on_ring(cells) == rule.phase_space(cells).is_bijective

    def test_shifts_shift(self):
        assert eca.cells_from_config(E(170).step(0b00100, 5), 5) == (0, 1, 0, 0, 0)
        assert eca.cells_from_config(E(240).step(0b00100, 5), 5) == (0, 0, 0, 1, 0)
        assert E(204).step(0b00100, 5) == 0b00100
        assert E(51).step(0b00100, 5) == 0b11011


class TestFiniteRings:
    @given(rules, ring_configs(max_cells=12))
    def test_trajectories_are_eventually_periodic_within_2_to_the_n(self, rule, ring):
        """Wolfram 1983, section 4."""
        config, cells = ring
        orbit = rule.orbit(config, cells)
        assert orbit.period >= 1
        assert orbit.transient + orbit.period <= 2**cells

    @given(rules, ring_configs(max_cells=10))
    def test_orbit_reports_the_first_repeat(self, rule, ring):
        config, cells = ring
        orbit = rule.orbit(config, cells)
        rows = rule.evolve(config, cells, orbit.transient + orbit.period)
        assert rows[-1] == rows[orbit.transient]
        assert len(set(rows[:-1])) == len(rows) - 1

    @given(rules, st.integers(1, 9))
    def test_phase_space_numbers_configurations_as_the_general_class(self, rule, cells):
        space = rule.phase_space(cells)
        general = rule.to_cellular_automaton().phase_space(cells)
        assert space == general


# --------------------------------------------------------------------------- #
# Rule 90 on a ring (MOW 1984, section 3)
# --------------------------------------------------------------------------- #
RING_SIZES = range(3, 13)


class TestRule90OnRings:
    @pytest.mark.parametrize("cells", RING_SIZES)
    def test_lemma_3_1_odd_configurations_are_never_generated(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Lemma 3.1."""
        space = eca.rule_90().phase_space(cells)
        assert all(image.bit_count() % 2 == 0 for image in space.successor)

    @pytest.mark.parametrize("cells", RING_SIZES)
    def test_theorem_3_1_fraction_of_gardens_of_eden(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.1."""
        space = eca.rule_90().phase_space(cells)
        unreachable = {1: 2, 0: 4}[cells % 2]
        assert len(space.gardens_of_eden) * unreachable == (
            (unreachable - 1) * 2**cells
        )

    @pytest.mark.parametrize("cells", RING_SIZES)
    def test_theorem_3_2_reachable_configurations_have_two_or_four_predecessors(
        self, cells
    ):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.2."""
        space = eca.rule_90().phase_space(cells)
        assert set(space.in_degrees) == {0, {1: 2, 0: 4}[cells % 2]}

    @pytest.mark.parametrize("cells", [n for n in RING_SIZES if n % 2])
    def test_theorem_3_3_odd_rings_hang_a_single_arc_on_every_cycle_node(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.3 and its corollary."""
        space = eca.rule_90().phase_space(cells)
        assert space.max_transient == 1
        assert all(space.in_degrees[state] == 2 for state in space.periodic_states)
        assert 2 * len(space.periodic_states) == 2**cells

    @pytest.mark.parametrize("cells", [n for n in RING_SIZES if n % 2 == 0])
    def test_theorem_3_4_even_rings_have_trees_of_height_half_d2(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.4 and its corollary."""
        space = eca.rule_90().phase_space(cells)
        d2 = eca.two_adic_part(cells)
        assert 2 * space.max_transient == d2
        assert all(space.in_degrees[state] == 4 for state in space.periodic_states)
        assert len(space.periodic_states) * 2**d2 == 2**cells

    def test_maximal_transients(self):
        assert (
            tuple(eca.rule_90().phase_space(n).max_transient for n in RING_SIZES)
            == RULE_90_MAX_TRANSIENT
        )

    @pytest.mark.parametrize("number", ADDITIVE)
    @pytest.mark.parametrize("cells", range(3, 11))
    def test_lemma_3_3_trees_at_all_cycle_nodes_are_identical(self, number, cells):
        """MOW 1984, Lemma 3.3 (and Theorem 4.3), for any additive rule."""
        trees = transient_trees(E(number).phase_space(cells))
        assert len(set(trees.values())) == 1

    @pytest.mark.parametrize("cells", RING_SIZES)
    def test_lemma_3_4_every_cycle_length_divides_pi_n(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Lemma 3.4."""
        rule = eca.rule_90()
        pi = rule.cycle_length_from_seed(cells)
        assert all(pi % length == 0 for length in rule.phase_space(cells).cycle_lengths)

    @pytest.mark.parametrize("cells", [2, 4, 8, 16, 32, 64])
    def test_lemma_3_5_pi_is_one_on_powers_of_two(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Lemma 3.5."""
        assert eca.rule_90().cycle_length_from_seed(cells) == 1

    @pytest.mark.parametrize("cells", [n for n in range(6, 31, 2) if n & (n - 1)])
    def test_lemma_3_6_pi_doubles_from_half_the_ring(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Lemma 3.6."""
        rule = eca.rule_90()
        assert rule.cycle_length_from_seed(cells) == 2 * rule.cycle_length_from_seed(
            cells // 2
        )

    @pytest.mark.parametrize("cells", range(3, 36, 2))
    def test_theorem_3_5_pi_divides_two_to_the_suborder_minus_one(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.5."""
        pi = eca.rule_90().cycle_length_from_seed(cells)
        assert (2 ** eca.multiplicative_suborder(cells) - 1) % pi == 0
        assert pi <= 2 ** ((cells - 1) // 2) - 1

    def test_theorem_3_5_first_proper_divisor_is_at_37_with_quotient_3(self):
        """Martin-Odlyzko-Wolfram 1984, section 3, Theorem 3.5."""
        rule = eca.rule_90()
        quotients = {
            n: (2 ** eca.multiplicative_suborder(n) - 1)
            // rule.cycle_length_from_seed(n)
            for n in range(3, 38, 2)
        }
        assert {n for n, quotient in quotients.items() if quotient != 1} == {37}
        assert quotients[37] == 3

    @pytest.mark.parametrize("cells", RING_SIZES)
    def test_lemma_3_7_fixed_points(self, cells):
        """Martin-Odlyzko-Wolfram 1984, section 3, Lemma 3.7."""
        blocks = ("000", "011", "101", "110")
        expected = (
            {int(block * (cells // 3), 2) for block in blocks}
            if cells % 3 == 0
            else {0}
        )
        assert set(eca.rule_90().phase_space(cells).fixed_points) == expected

    def test_table_1_pi_n(self):
        """Martin-Odlyzko-Wolfram 1984, Table 1 (OEIS A085587)."""
        rule = eca.rule_90()
        assert tuple(rule.cycle_length_from_seed(n) for n in range(3, 31)) == PI_90

    def test_table_2_cycles_on_a_ring_of_ten(self):
        """Martin-Odlyzko-Wolfram 1984, Table 2."""
        assert Counter(eca.rule_90().phase_space(10).cycle_lengths) == {
            1: 1, 3: 5, 6: 40,
        }  # fmt: skip


# --------------------------------------------------------------------------- #
# Rule 150 and general additive rules on a ring (MOW 1984, section 4)
# --------------------------------------------------------------------------- #
class TestAdditiveRulesOnRings:
    @pytest.mark.parametrize(("number", "exponents"), DIPOLYNOMIALS.items())
    def test_dipolynomial_dictionary(self, number, exponents):
        assert E(number).dipolynomial == exponents

    @pytest.mark.parametrize(("number", "exponents"), DIPOLYNOMIALS.items())
    @pytest.mark.parametrize("cells", [5, 6, 7])
    def test_a_step_is_multiplication_by_t_modulo_x_n_minus_1(
        self, number, exponents, cells
    ):
        rule = E(number)

        def product(a: Bits) -> Bits:
            return tuple(
                sum(a[(i - exponent) % cells] for exponent in exponents) % 2
                for i in range(cells)
            )

        assert all(
            eca.cells_from_config(rule.step(config, cells), cells)
            == product(eca.cells_from_config(config, cells))
            for config in range(1 << cells)
        )

    @pytest.mark.parametrize("number", ADDITIVE)
    @pytest.mark.parametrize("cells", range(3, 10))
    def test_lemma_4_4_reachable_iff_divisible_by_lambda_1(self, number, cells):
        """Martin-Odlyzko-Wolfram 1984, section 4, Lemma 4.4."""
        rule = E(number)
        divisor = lambda_1(rule, cells)
        reachable = set(rule.phase_space(cells).successor)
        assert reachable == {
            config
            for config in range(1 << cells)
            if poly_mod(polynomial(config, cells), divisor) == 0
        }

    @pytest.mark.parametrize("number", ADDITIVE)
    @pytest.mark.parametrize("cells", range(3, 10))
    def test_theorems_4_2_and_4_3_reachable_fraction_and_predecessors(
        self, number, cells
    ):
        """Martin-Odlyzko-Wolfram 1984, section 4, Theorems 4.2 and 4.3."""
        rule = E(number)
        degree = lambda_1(rule, cells).bit_length() - 1
        space = rule.phase_space(cells)
        assert len(set(space.successor)) * 2**degree == 2**cells
        assert set(space.in_degrees) <= {0, 2**degree}

    @pytest.mark.parametrize("number", ADDITIVE)
    @pytest.mark.parametrize("cells", range(3, 11))
    def test_lemma_4_2_every_cycle_length_divides_pi_n(self, number, cells):
        """Martin-Odlyzko-Wolfram 1984, section 4, Lemma 4.2."""
        rule = E(number)
        pi = rule.cycle_length_from_seed(cells)
        assert all(pi % length == 0 for length in rule.phase_space(cells).cycle_lengths)

    @pytest.mark.parametrize("number", [90, 150])
    @pytest.mark.parametrize("cells", range(3, 30, 2))
    def test_theorem_4_1_symmetric_t_on_odd_rings(self, number, cells):
        """Martin-Odlyzko-Wolfram 1984, section 4, Theorem 4.1."""
        pi = E(number).cycle_length_from_seed(cells)
        assert (2 ** eca.multiplicative_suborder(cells) - 1) % pi == 0

    def test_rule_150_table_1_pi_n(self):
        """Martin-Odlyzko-Wolfram 1984, Table 1 (OEIS A085588)."""
        rule = eca.rule_150()
        assert tuple(rule.cycle_length_from_seed(n) for n in range(3, 31)) == PI_150

    @pytest.mark.parametrize("cells", range(1, 13))
    def test_rule_150_transients(self, cells):
        expected = 0
        if cells % 3 == 0:
            expected = 1 if cells % 2 else eca.two_adic_part(cells)
        assert eca.rule_150().transient_from_seed(cells) == expected

    def test_pi_n_differs_between_rules_90_and_150(self):
        assert [
            n for n, a, b in zip(range(3, 31), PI_90, PI_150, strict=True) if a != b
        ] == [4, 6, 8, 12, 13, 16, 24, 26]
        assert (PI_90[13 - 3], PI_150[13 - 3]) == (63, 21)

    def test_pi_n_of_150_divides_that_of_90_when_n_is_coprime_to_6(self):
        assert all(
            eca.rule_90().cycle_length_from_seed(n)
            % eca.rule_150().cycle_length_from_seed(n)
            == 0
            for n in range(3, 31)
            if math.gcd(n, 6) == 1
        )

    @pytest.mark.parametrize(
        ("number", "bijective_when"),
        [
            (150, lambda n: n % 3 != 0),
            (105, lambda n: n % 3 != 0),
            (45, lambda n: n % 2 == 1),
            (154, lambda n: n % 2 == 1),
            (30, lambda n: False),
            (60, lambda n: False),
            (90, lambda n: False),
            (106, lambda n: False),
        ],
    )
    def test_ring_bijectivity_depends_on_n(self, number, bijective_when):
        rule = E(number)
        assert [rule.is_bijective_on_ring(n) for n in range(1, 17)] == [
            bijective_when(n) for n in range(1, 17)
        ]

    @pytest.mark.parametrize("cells", [3, 6, 9, 12])
    def test_rule_150_reaches_a_quarter_with_four_predecessors_when_3_divides_n(
        self, cells
    ):
        space = eca.rule_150().phase_space(cells)
        assert set(space.in_degrees) == {0, 4}
        assert 4 * len(set(space.successor)) == 2**cells


class TestAdditiveRulesOnTheLine:
    def test_rule_0_is_the_only_non_surjective_additive_rule(self):
        """Kari 2005, Corollary 4 (Ito, Osato and Nasu 1983), with m = 2."""
        assert [n for n in ADDITIVE if not E(n).is_surjective] == [0]

    def test_rules_170_204_240_are_the_injective_additive_rules(self):
        """Kari 2005, Corollary 4 (Ito, Osato and Nasu 1983), with m = 2."""
        assert [n for n in ADDITIVE if E(n).is_injective] == [170, 204, 240]

    def test_equicontinuity_of_the_additive_rules(self):
        """Kari 2005, Theorem 18 (Cattaneo-Formenti-Manzini-Margara 2000;
        Manzini-Margara 1999), with m = 2, read off the Schule-Stoop table."""
        classes = {n: E(n).kurka_class for n in ADDITIVE}
        assert [
            n for n, c in classes.items() if c is eca.KurkaClass.EQUICONTINUOUS
        ] == [0, 204]
        assert [
            n for n, c in classes.items() if c is eca.KurkaClass.POSITIVELY_EXPANSIVE
        ] == [90, 150]
        assert [n for n, c in classes.items() if c is eca.KurkaClass.SENSITIVE] == [
            60, 102, 170, 240,
        ]  # fmt: skip


# --------------------------------------------------------------------------- #
# Growth from a single 1
# --------------------------------------------------------------------------- #
class TestSingleSeed:
    def test_rule_90_draws_pascals_triangle_modulo_two(self):
        rows = eca.rule_90().evolve_from_single_cell(127)
        assert all(
            rows[t][t + j]
            == (math.comb(t, (t + j) // 2) % 2 if (t + j) % 2 == 0 else 0)
            for t in range(128)
            for j in range(-t, t + 1)
        )

    def test_rule_90_row_counts_are_two_to_the_binary_weight(self):
        rows = eca.rule_90().evolve_from_single_cell(127)
        assert [sum(row) for row in rows] == [2 ** t.bit_count() for t in range(128)]
        assert tuple(sum(row) for row in rows[:8]) == GOULD

    def test_rule_150_rows_are_trinomial_coefficients_modulo_two(self):
        rows = eca.rule_150().evolve_from_single_cell(63)
        coefficients = [1]
        for t in range(64):
            assert rows[t] == tuple(c % 2 for c in coefficients)
            padded = [0, 0, *coefficients, 0, 0]
            coefficients = [sum(padded[i : i + 3]) for i in range(len(padded) - 2)]

    def test_rule_150_row_counts(self):
        rows = eca.rule_150().evolve_from_single_cell(15)
        assert tuple(sum(row) for row in rows) == RULE_150_COUNTS

    def test_rule_18_agrees_with_rule_90_from_a_single_1_up_to_t_64(self):
        assert eca.rule_18().evolve_from_single_cell(64) == (
            eca.rule_90().evolve_from_single_cell(64)
        )

    def test_rule_22_row_counts(self):
        rows = eca.rule_22().evolve_from_single_cell(7)
        assert tuple(sum(row) for row in rows) == RULE_22_COUNTS

    def test_rule_128_erases_a_single_1_but_fixes_all_ones(self):
        rule = eca.rule_128()
        assert rule.evolve_from_single_cell(3)[1:] == ((0,) * 3, (0,) * 5, (0,) * 7)
        assert rule.step(0b1111111, 7) == 0b1111111

    @given(rules, st.integers(0, 12))
    def test_single_seed_rows_agree_with_a_wide_ring(self, rule, steps):
        width = 4 * steps + 3
        center = 2 * steps + 1
        wide = rule.evolve(1 << (width - 1 - center), width, steps)
        assert rule.evolve_from_single_cell(steps) == tuple(
            eca.cells_from_config(row, width)[center - t : center + t + 1]
            for t, row in enumerate(wide)
        )


# --------------------------------------------------------------------------- #
# Rule 30
# --------------------------------------------------------------------------- #
class TestRule30:
    def test_left_permutive_hence_surjective_with_four_predecessors(self):
        rule = eca.rule_30()
        assert rule.is_left_permutive
        assert not rule.is_additive
        assert rule.is_surjective
        assert all(
            rule.preimage_count(image) == 4
            for image in itertools.product(range(2), repeat=6)
        )

    def test_class_and_wolfram_class(self):
        rule = eca.rule_30()
        assert {member.rule for member in rule.equivalence_class()} == {
            30, 86, 135, 149,
        }  # fmt: skip
        assert rule.wolfram_class is eca.WolframClass.III

    def test_bijective_on_no_ring_from_4_to_16(self):
        assert not any(eca.rule_30().is_bijective_on_ring(n) for n in range(4, 17))

    @pytest.mark.parametrize("cells", range(4, 14))
    def test_ring_in_degrees(self, cells):
        degrees = Counter(eca.rule_30().phase_space(cells).in_degrees)
        assert set(degrees) == ({0, 1, 2, 3} if cells % 3 == 0 else {0, 1, 2})
        assert degrees[3] == (1 if cells % 3 == 0 else 0)
        assert degrees[0] - degrees[2] == (2 if cells % 3 == 0 else 0)

    def test_ring_of_six_has_12_gardens_and_10_doubly_reached(self):
        degrees = Counter(eca.rule_30().phase_space(6).in_degrees)
        assert (degrees[0], degrees[2]) == (12, 10)

    def test_all_ones_is_the_configuration_with_three_predecessors(self):
        space = eca.rule_30().phase_space(9)
        assert [c for c, degree in enumerate(space.in_degrees) if degree == 3] == [
            2**9 - 1
        ]

    # PAGE DISCREPANCY: the page says a configuration "has a unique
    # predecessor unless it contains two 0s separated by 3n+1 consecutive
    # 1s". Exhaustively, the polarity is the other way round; this asserts
    # what the computation gives and is reported, not transcribed.
    @pytest.mark.parametrize("cells", range(4, 14))
    def test_unique_predecessor_iff_a_run_of_3n_plus_1_ones_PAGE_DISCREPANCY(  # noqa: N802
        self, cells
    ):
        space = eca.rule_30().phase_space(cells)
        unique = {c for c, degree in enumerate(space.in_degrees) if degree == 1}
        assert unique == {
            config
            for config in range(2**cells - 1)
            if any(
                run % 3 == 1
                for run in cyclic_runs_of_ones(eca.cells_from_config(config, cells))
            )
        }

    def test_maximal_cycle_lengths_on_rings(self):
        rule = eca.rule_30()
        assert (
            tuple(max(rule.phase_space(n).cycle_lengths) for n in range(1, 17))
            == RULE_30_MAX_CYCLE
        )

    def test_maximal_cycle_length_on_a_ring_of_17_is_a334497_not_table_9_2(self):
        """Wolfram 1986: 10846 in Table 9.1, 10845 in Table 9.2; A334497 has 10846."""
        assert max(eca.rule_30().phase_space(17).cycle_lengths) == 10846

    def test_center_column_is_a051023(self):
        column = eca.rule_30().center_column(len(RULE_30_COLUMN) - 1)
        assert "".join(map(str, column)) == RULE_30_COLUMN

    def test_row_counts_are_a070952(self):
        rows = eca.rule_30().evolve_from_single_cell(15)
        assert tuple(sum(row) for row in rows) == RULE_30_COUNTS

    def test_rows_as_binary_integers_are_a110240(self):
        rows = eca.rule_30().evolve_from_single_cell(7)
        assert tuple(int("".join(map(str, row)), 2) for row in rows) == RULE_30_ROWS

    def test_solving_for_the_left_cell_is_wolfram_1986_eq_3_3(self):
        rule = eca.rule_30()
        assert all(
            rule.solve_left(value, q, r) == value ^ (q | r)
            for value, q, r in itertools.product(range(2), repeat=3)
        )

    def test_two_adjacent_columns_determine_the_pattern_to_their_left(self):
        rule, steps = eca.rule_30(), 24
        rows = rule.evolve_from_single_cell(steps)

        def cell(t: int, j: int) -> int:
            return rows[t][t + j] if abs(j) <= t else 0

        assert all(
            rule.solve_left(cell(t + 1, j), cell(t, j), cell(t, j + 1))
            == cell(t, j - 1)
            for t in range(steps)
            for j in range(-t - 1, t + 2)
        )

    def test_rule_30_is_rule_150_plus_the_term_qr(self):
        assert set(eca.rule_30().anf) - set(eca.rule_150().anf) == {"qr"}
        assert set(eca.rule_150().anf) <= set(eca.rule_30().anf)


# --------------------------------------------------------------------------- #
# Rule 110
# --------------------------------------------------------------------------- #
class TestRule110:
    def test_cooks_description(self):
        described = E.from_local_rule(lambda p, q, r: r if q == 0 else 1 - (p & r))
        assert described == eca.rule_110()

    def test_class(self):
        assert {member.rule for member in eca.rule_110().equivalence_class()} == {
            110, 124, 137, 193,
        }  # fmt: skip

    def test_not_surjective_three_neighborhoods_to_0_and_five_to_1(self):
        rule = eca.rule_110()
        assert Counter(rule.table) == {0: 3, 1: 5}
        assert not rule.is_surjective

    def test_shortest_orphan_is_01010(self):
        assert eca.rule_110().shortest_orphans() == (word("01010"),)

    def test_diamond_of_karis_example(self):
        rule = eca.rule_110()
        first, second = word("00110100"), word("00101100")
        assert rule.apply_to_word(first) == word("111110")
        assert rule.apply_to_word(second) == word("111110")
        assert (first[:2], first[-2:]) == (second[:2], second[-2:])

    def test_diamond_search_returns_a_diamond(self):
        rule = eca.rule_110()
        found = rule.diamond(4)
        assert found is not None
        first, second = found
        assert first != second
        assert (len(first), len(second)) == (8, 8)
        assert (first[:2], first[-2:]) == (second[:2], second[-2:])
        assert rule.apply_to_word(first) == rule.apply_to_word(second)

    def test_subset_construction_reaches_8_of_the_16_subsets(self):
        """Kari's lecture notes, Example 14."""
        rule = eca.rule_110()
        edges = rule.de_bruijn_graph()
        reached = {frozenset(source for source, _, _ in edges)}
        frontier = list(reached)
        while frontier:
            subset = frontier.pop()
            for label in (0, 1):
                following = frozenset(
                    target
                    for source, target, value in edges
                    if source in subset and value == label
                )
                if following not in reached:
                    reached.add(following)
                    frontier.append(following)
        assert len(reached) == 8

    def test_differs_from_rule_102_in_exactly_one_case(self):
        assert eca.rule_110().hamming_distance(eca.rule_102()) == 1

    def test_sensitive_class_k3_and_wolfram_class_iv(self):
        rule = eca.rule_110()
        assert rule.kurka_class is eca.KurkaClass.SENSITIVE
        assert rule.wolfram_class is eca.WolframClass.IV

    def test_row_counts_are_a071049(self):
        rows = eca.rule_110().evolve_from_single_cell(15)
        assert tuple(sum(row) for row in rows) == RULE_110_COUNTS

    def test_rows_as_binary_integers_are_a006978(self):
        rows = eca.rule_110().evolve_from_single_cell(7)
        assert all(not any(row[t + 1 :]) for t, row in enumerate(rows))
        assert (
            tuple(int("".join(map(str, row[: t + 1])), 2) for t, row in enumerate(rows))
            == RULE_110_ROWS
        )

    def test_ether_block_returns_after_exactly_seven_steps(self):
        ether = eca.config_from_cells(eca.rule_110_ether())
        rows = eca.rule_110().evolve(ether, 14, 7)
        assert rows[7] == ether
        assert ether not in rows[1:7]

    def test_ether_cell_formula(self):
        block = eca.rule_110_ether()
        rows = eca.rule_110().evolve(eca.config_from_cells(block), 14, 7)
        assert "".join(map(str, block)) == "10011011111000"
        assert all(
            eca.cells_from_config(rows[t], 14)[x] == block[(x + 4 * t) % 14]
            for t in range(8)
            for x in range(14)
        )

    def test_maximal_cycle_lengths_on_rings(self):
        rule = eca.rule_110()
        assert (
            tuple(max(rule.phase_space(n).cycle_lengths) for n in range(1, 17))
            == RULE_110_MAX_CYCLE
        )


# --------------------------------------------------------------------------- #
# Rule 184, rule 232 and number conservation
# --------------------------------------------------------------------------- #
def _all_00_blocks_gone(config: int, cells: int) -> bool:
    digits = eca.cells_from_config(config, cells)
    return not any(digits[i] == digits[(i + 1) % cells] == 0 for i in range(cells))


class TestTrafficAndDensity:
    @given(ring_configs(min_cells=2))
    def test_rule_184_moves_every_1_with_a_0_to_its_right(self, ring):
        config, cells = ring
        c = eca.cells_from_config(config, cells)
        moved = tuple(
            int(
                (c[i] == 1 and c[(i + 1) % cells] == 1)
                or (c[i] == 0 and c[(i - 1) % cells] == 1)
            )
            for i in range(cells)
        )
        assert eca.cells_from_config(eca.rule_184().step(config, cells), cells) == moved

    @given(ring_configs(min_cells=2))
    def test_rule_226_replaces_01_by_10(self, ring):
        config, cells = ring
        c = eca.cells_from_config(config, cells)
        swapped = list(c)
        for i in range(cells):
            if c[i] == 0 and c[(i + 1) % cells] == 1:
                swapped[i], swapped[(i + 1) % cells] = 1, 0
        assert eca.cells_from_config(eca.rule_226().step(config, cells), cells) == (
            tuple(swapped)
        )

    def test_rule_184_is_self_dual(self):
        rule = eca.rule_184()
        assert rule.reflect().conjugate() == rule
        assert eca.rule_226() == rule.reflect() == rule.conjugate()

    @given(ring_configs())
    def test_rule_184_conserves_the_number_of_cars(self, ring):
        config, cells = ring
        assert eca.rule_184().step(config, cells).bit_count() == config.bit_count()

    @pytest.mark.parametrize("cells", range(2, 13))
    def test_fuks_proposition_2_00_blocks_disappear(self, cells):
        """Fuks 1997, Proposition 2."""
        rule = eca.rule_184()
        assert all(
            _all_00_blocks_gone(rule.evolve(config, cells, (cells - 2) // 2)[-1], cells)
            for config in range(1 << cells)
            if 2 * config.bit_count() > cells
        )

    @pytest.mark.parametrize("cells", range(3, 17))
    def test_fuks_proposition_4_exhaustively(self, cells):
        """Fuks 1997, Proposition 4."""
        full = (1 << cells) - 1
        alternating = {int("01" * cells, 2) & full, int("10" * cells, 2) & full}
        outcomes = {
            -1: {0},
            0: alternating,
            1: {full},
        }
        assert all(
            eca.density_classifier(config, cells)
            in outcomes[
                (2 * config.bit_count() > cells) - (2 * config.bit_count() < cells)
            ]
            for config in range(1 << cells)
        )

    @given(ring_configs(min_cells=17, max_cells=40))
    def test_fuks_proposition_4_on_larger_rings(self, ring):
        """Fuks 1997, Proposition 4."""
        config, cells = ring
        full = (1 << cells) - 1
        ones = config.bit_count()
        expected = {
            -1: {0},
            0: {int("01" * cells, 2) & full, int("10" * cells, 2) & full},
            1: {full},
        }[(2 * ones > cells) - (2 * ones < cells)]
        assert eca.density_classifier(config, cells) in expected

    def test_density_classifier_is_the_two_rule_schedule(self):
        config, cells = 0b0110100, 7
        after_traffic = eca.rule_184().evolve(config, cells, 2)[-1]
        assert (
            eca.density_classifier(config, cells)
            == (eca.rule_232().evolve(after_traffic, cells, 3)[-1])
        )


class TestNumberConservation:
    def test_number_conserving_rules_are_the_five_of_nks_p_1022(self):
        assert (
            tuple(rule.rule for rule in ALL if rule.is_number_conserving)
            == NUMBER_CONSERVING
        )

    def test_theorem_2_1_identity_holds_iff_the_rule_conserves_number(self):
        """Boccara-Fuks, Theorem 2.1, for three inputs."""

        def identity(f: LocalRule) -> bool:
            return all(
                f(x1, x2, x3)
                == x1 + (f(0, x2, x3) - f(0, x1, x2)) + (f(0, 0, x2) - f(0, 0, x1))
                for x1, x2, x3 in itertools.product(range(2), repeat=3)
            )

        assert all(
            identity(rule.local_rule) == rule.conserves_number_on_ring(5)
            for rule in ALL
        )

    def test_rings_of_five_cells_decide_conservation(self):
        """Boccara-Fuks, Remark 2.1."""
        assert (
            tuple(rule.rule for rule in ALL if rule.conserves_number_on_ring(5))
            == NUMBER_CONSERVING
        )

    @pytest.mark.parametrize("number", NUMBER_CONSERVING)
    def test_the_five_conserve_on_rings_of_up_to_nine_cells(self, number):
        assert all(E(number).conserves_number_on_ring(n) for n in range(3, 10))

    @pytest.mark.parametrize("number", NUMBER_CONSERVING)
    def test_corollaries_2_1_and_2_2(self, number):
        """Boccara-Fuks, Corollaries 2.1 and 2.2 (necessary conditions)."""
        rule = E(number)
        assert (rule.local_rule(0, 0, 0), rule.local_rule(1, 1, 1)) == (0, 1)
        assert rule.is_balanced

    def test_boccara_fuks_closed_form_for_rule_184(self):
        assert eca.rule_184() == E.from_local_rule(
            lambda x1, x2, x3: x2 + min(x1, 1 - x2) - min(x2, 1 - x3)
        )


# --------------------------------------------------------------------------- #
# Published class tables
# --------------------------------------------------------------------------- #
class TestClassTables:
    @pytest.mark.parametrize(("label", "members"), KURKA.items())
    def test_kurka_classes_of_schule_stoop(self, label, members):
        assert (
            tuple(r.rule for r in eca.representatives() if r.kurka_class is label)
            == members
        )

    def test_kurka_class_sizes(self):
        assert Counter(r.kurka_class.value for r in eca.representatives()) == {
            "K1": 15, "K2": 24, "K3": 46, "K4": 3,
        }  # fmt: skip
        assert Counter(rule.kurka_class.value for rule in ALL) == {
            "K1": 32, "K2": 62, "K3": 158, "K4": 4,
        }  # fmt: skip

    def test_proposition_13_positively_expansive_iff_bipermutive(self):
        """Schule-Stoop 2012, Proposition 13."""
        assert all(
            (rule.kurka_class is eca.KurkaClass.POSITIVELY_EXPANSIVE)
            == rule.is_bipermutive
            == rule.is_positively_expansive
            for rule in ALL
        )

    def test_corollary_17_permutive_iff_surjective_and_sensitive(self):
        """Schule-Stoop 2012, Proposition 16 and Corollary 17."""
        sensitive = {
            eca.KurkaClass.SENSITIVE,
            eca.KurkaClass.POSITIVELY_EXPANSIVE,
        }
        assert all(
            rule.is_permutive == (rule.is_surjective and rule.kurka_class in sensitive)
            for rule in ALL
        )

    @pytest.mark.parametrize(("label", "members"), WOLFRAM.items())
    def test_wolfram_classes_of_martinez(self, label, members):
        assert (
            tuple(r.rule for r in eca.representatives() if r.wolfram_class is label)
            == members
        )

    def test_wolfram_class_ii_is_the_remaining_65(self):
        assert Counter(r.wolfram_class.value for r in eca.representatives()) == {
            "W1": 8, "W2": 65, "W3": 11, "W4": 4,
        }  # fmt: skip

    @pytest.mark.parametrize(("label", "members"), LI_PACKARD.items())
    def test_li_packard_classes(self, label, members):
        assert (
            tuple(r.rule for r in eca.representatives() if r.li_packard_class is label)
            == members
        )

    def test_li_packard_class_sizes_over_all_256_rules(self):
        assert Counter(rule.li_packard_class.value for rule in ALL) == {
            "null": 24,
            "fixed point": 97,
            "periodic": 89,
            "locally chaotic": 10,
            "chaotic": 36,
        }

    def test_null_class_is_wolfram_class_i(self):
        assert all(
            (rule.li_packard_class is eca.LiPackardClass.NULL)
            == (rule.wolfram_class is eca.WolframClass.I)
            for rule in ALL
        )

    def test_classes_are_constant_on_equivalence_classes(self):
        assert all(
            (member.kurka_class, member.wolfram_class, member.li_packard_class)
            == (rule.kurka_class, rule.wolfram_class, rule.li_packard_class)
            for rule in ALL
            for member in rule.equivalence_class()
        )

    def test_the_sources_disagree_on_class_iv(self):
        assert all(
            E(n).li_packard_class is eca.LiPackardClass.CHAOTIC for n in (54, 110, 137)
        )
        assert all(E(n).wolfram_class is eca.WolframClass.IV for n in (41, 54, 106))


# --------------------------------------------------------------------------- #
# Orphans and rule 18's Garden of Eden
# --------------------------------------------------------------------------- #
class TestOrphans:
    @pytest.mark.parametrize(
        ("number", "orphans"),
        [
            (110, ("01010",)),
            (184, ("1100",)),
            (232, ("01001", "01101", "10010", "10110")),
            (54, ("01101", "10101", "10110")),
            (22, ("10010101", "10101001")),
        ],
    )
    def test_shortest_orphans(self, number, orphans):
        assert E(number).shortest_orphans() == tuple(word(text) for text in orphans)

    @given(rules)
    def test_surjective_iff_no_orphan(self, rule):
        assert rule.is_surjective == (rule.shortest_orphans() == ())

    @given(rules, st.lists(st.integers(0, 1), min_size=1, max_size=8))
    def test_orphans_are_the_words_without_preimage(self, rule, image):
        assert rule.is_orphan(image) == (rule.preimage_count(image) == 0)


class TestRule18:
    @pytest.mark.parametrize("cells", range(3, 13))
    def test_lemma_5_1_characterizes_the_unreachable_configurations(self, cells):
        """Martin-Odlyzko-Wolfram 1984, Lemma 5.1."""
        gardens = set(eca.rule_18().phase_space(cells).gardens_of_eden)
        assert gardens == {
            config
            for config in range(1 << cells)
            if unreachable_under_rule_18(eca.cells_from_config(config, cells))
        }

    def test_number_of_reachable_configurations(self):
        rule = eca.rule_18()
        assert (
            tuple(len(set(rule.phase_space(n).successor)) for n in range(1, 11))
            == RULE_18_REACHABLE
        )

    def test_legal_and_non_additive(self):
        rule = eca.rule_18()
        assert rule == E.from_local_rule(lambda p, q, r: (1 - q) & (p ^ r))
        assert rule.is_legal
        assert not rule.is_additive


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_reversible_rules_fixture(self):
        rules_ = eca.reversible_rules()
        assert [rule.rule for rule in rules_] == [204, 51, 170, 240, 85, 15]
        assert all(rule.is_reversible for rule in rules_)

    def test_reversibility_says_nothing_about_equicontinuity(self):
        assert E(204).kurka_class is eca.KurkaClass.EQUICONTINUOUS
        assert E(170).kurka_class is eca.KurkaClass.SENSITIVE
        assert E(170).is_devaney_chaotic

    def test_rule_90_is_additive_bipermutive_and_positively_expansive(self):
        rule = eca.rule_90()
        assert rule.anf == ("p", "r")
        assert rule.is_additive
        assert rule.is_bipermutive
        assert rule.is_positively_expansive
        assert rule.kurka_class is eca.KurkaClass.POSITIVELY_EXPANSIVE

    def test_rule_90_is_surjective_on_the_line_but_not_on_rings(self):
        rule = eca.rule_90()
        assert rule.is_surjective
        assert all(
            rule.preimage_count(image) == 4
            for image in itertools.product(range(2), repeat=6)
        )
        assert [
            len(rule.phase_space(n).gardens_of_eden) * 4 // 2**n for n in (5, 6, 7, 8)
        ] == [2, 3, 2, 3]

    def test_rule_150_is_additive_totalistic_and_self_equivalent(self):
        rule = eca.rule_150()
        assert rule.anf == ("p", "q", "r")
        assert rule.is_additive
        assert rule.is_totalistic
        assert rule.reflect() == rule == rule.conjugate()
        assert rule.is_positively_expansive

    def test_rule_150_ring_reversibility_depends_on_n(self):
        rule = eca.rule_150()
        assert [n for n in range(1, 17) if not rule.is_bijective_on_ring(n)] == [
            3, 6, 9, 12, 15,
        ]  # fmt: skip

    def test_rules_60_and_102_are_one_sided(self):
        left, right = eca.rule_60(), eca.rule_102()
        assert (left.anf, right.anf) == (("p", "q"), ("q", "r"))
        assert left == right.reflect()
        assert (left.is_left_permutive, left.is_right_permutive) == (True, False)
        assert (right.is_left_permutive, right.is_right_permutive) == (False, True)

    @pytest.mark.parametrize("example", [eca.rule_60, eca.rule_102])
    def test_one_sided_permutive_is_chaotic_but_only_k3(self, example):
        rule = example()
        assert rule.is_devaney_chaotic
        assert rule.kurka_class is eca.KurkaClass.SENSITIVE
        assert not rule.is_positively_expansive

    def test_rule_102_is_karis_radius_half_xor(self):
        assert E.from_cellular_automaton(ca.radius_half_xor()) == eca.rule_102()

    def test_rule_45_is_bijective_on_rings_exactly_for_odd_n(self):
        rule = eca.rule_45()
        assert rule == E.from_local_rule(lambda p, q, r: p ^ (q | (1 - r)))
        assert [n for n in range(1, 17) if rule.is_bijective_on_ring(n)] == list(
            range(1, 17, 2)
        )

    def test_rule_54_fixture(self):
        rule = eca.rule_54()
        assert rule.is_legal
        assert {member.rule for member in rule.equivalence_class()} == {54, 147}

    def test_rule_184_balanced_number_conserving_not_surjective(self):
        rule = eca.rule_184()
        assert rule.is_number_conserving
        assert not rule.is_trivial
        assert rule.is_balanced
        assert not rule.is_surjective
        assert rule.kurka_class is eca.KurkaClass.SENSITIVE
        assert rule.li_packard_class is eca.LiPackardClass.FIXED_POINT
        assert rule.wolfram_class is eca.WolframClass.II

    def test_rule_232_is_the_totalistic_self_equivalent_majority_rule(self):
        rule = eca.rule_232()
        assert rule == E.from_local_rule(lambda p, q, r: int(p + q + r >= 2))
        assert rule.is_totalistic
        assert len(rule.equivalence_class()) == 1
        assert rule.is_balanced
        assert not rule.is_surjective
        assert rule.kurka_class is eca.KurkaClass.ALMOST_EQUICONTINUOUS

    def test_rules_128_and_136_are_conjunctions(self):
        assert eca.rule_128() == E.from_local_rule(lambda p, q, r: p & q & r)
        assert eca.rule_136() == E.from_local_rule(lambda p, q, r: q & r)

    def test_rule_22_is_totalistic_with_orphans_of_length_eight(self):
        rule = eca.rule_22()
        assert rule.is_totalistic
        assert {len(orphan) for orphan in rule.shortest_orphans()} == {8}

    def test_rule_126_on_a_ring_of_eight(self):
        space = eca.rule_126().phase_space(8)
        assert space.fixed_points == (0,)
        assert Counter(space.cycle_lengths) == {1: 1, 6: 4, 2: 2}
        assert len(space.periodic_states) == 29
        assert len(space.gardens_of_eden) == 190

    def test_rule_73_shows_the_class_tables_are_not_refinements(self):
        rule = eca.rule_73()
        assert {member.rule for member in rule.equivalence_class()} == {73, 109}
        assert rule.li_packard_class is eca.LiPackardClass.LOCALLY_CHAOTIC
        assert rule.wolfram_class is eca.WolframClass.II
        assert rule.kurka_class is eca.KurkaClass.ALMOST_EQUICONTINUOUS

    def test_rule_180_is_left_permutive_hence_devaney_chaotic(self):
        rule = eca.rule_180()
        assert rule.is_left_permutive
        assert rule.is_devaney_chaotic
        assert rule.kurka_class is eca.KurkaClass.SENSITIVE
        assert rule.compose(E(170)).is_surjective()


# --------------------------------------------------------------------------- #
# Operations
# --------------------------------------------------------------------------- #
class TestOperations:
    @given(rules, ring_configs(min_cells=1, max_cells=10))
    def test_composition_with_itself_gives_alternate_time_steps(self, rule, ring):
        config, cells = ring
        squared = rule.compose(rule)
        assert squared.step(eca.cells_from_config(config, cells)) == (
            eca.cells_from_config(rule.evolve(config, cells, 2)[-1], cells)
        )

    @given(rules, rules, ring_configs(max_cells=10))
    def test_compose_applies_the_other_rule_first(self, first, second, ring):
        config, cells = ring
        composed = first.compose(second)
        assert composed.step(eca.cells_from_config(config, cells)) == (
            eca.cells_from_config(first.step(second.step(config, cells), cells), cells)
        )

    def test_composition_generally_has_range_two(self):
        composed = eca.rule_30().compose(eca.rule_110()).minimal()
        assert composed.offsets == ((-2,), (-1,), (0,), (1,), (2,))

    def test_second_order_variant_is_invertible_for_every_rule(self):
        """Wolfram 1983, section 4, after Fredkin and Margolus."""
        assert all(rule.second_order().is_reversible() for rule in ALL)

    @given(rules, ring_configs(min_cells=3, max_cells=8), st.data())
    def test_second_order_variant_is_f_of_s_xor_the_past(self, rule, ring, data):
        present, cells = ring
        past = data.draw(st.integers(0, 2**cells - 1))
        pairs = tuple(
            2 * now + before
            for now, before in zip(
                eca.cells_from_config(present, cells),
                eca.cells_from_config(past, cells),
                strict=True,
            )
        )
        future = rule.step(present, cells) ^ past
        assert rule.second_order().step(pairs) == tuple(
            2 * new + now
            for new, now in zip(
                eca.cells_from_config(future, cells),
                eca.cells_from_config(present, cells),
                strict=True,
            )
        )

    @given(st.sampled_from(LEFT_PERMUTIVE).map(E))
    def test_solve_left_inverts_the_rule_in_p(self, rule):
        assert all(
            rule.local_rule(rule.solve_left(value, q, r), q, r) == value
            for value, q, r in itertools.product(range(2), repeat=3)
        )

    @given(st.sampled_from(RIGHT_PERMUTIVE).map(E))
    def test_solve_right_inverts_the_rule_in_r(self, rule):
        assert all(
            rule.local_rule(p, q, rule.solve_right(p, q, value)) == value
            for p, q, value in itertools.product(range(2), repeat=3)
        )

    def test_rule_22_simulates_rule_146_by_blocks_at_half_speed(self):
        def encode(config: int, cells: int) -> int:
            doubled = []
            for value in eca.cells_from_config(config, cells):
                doubled += [0, value]
            return eca.config_from_cells(doubled)

        cells = 7
        assert all(
            E(22).evolve(encode(config, cells), 2 * cells, 2)[-1]
            == encode(E(146).step(config, cells), cells)
            for config in range(1 << cells)
        )


# --------------------------------------------------------------------------- #
# Round trips
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @given(rules)
    def test_dataframe_round_trip(self, rule):
        assert E.from_dataframe(rule.to_dataframe()) == rule

    @given(rules, st.integers(0, 2**16))
    def test_dataframe_rows_may_come_in_any_order(self, rule, seed):
        frame = rule.to_dataframe().sample(frac=1, random_state=seed)
        assert E.from_dataframe(frame) == rule

    def test_dataframe_is_the_rule_table(self):
        frame = eca.rule_110().to_dataframe()
        assert list(frame.columns) == ["p", "q", "r", "value"]
        assert frame.to_numpy().tolist() == [
            [0, 0, 0, 0],
            [0, 0, 1, 1],
            [0, 1, 0, 1],
            [0, 1, 1, 1],
            [1, 0, 0, 0],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ]

    @pytest.mark.parametrize("example", [eca.rule_30, eca.rule_110, eca.rule_184])
    def test_dataframe_survives_write_result(self, example, tmp_path):
        path = io.write_result(example().to_dataframe(), tmp_path / "rule.json")
        assert E.from_dataframe(pd.read_json(path, dtype=False)) == example()

    def test_constructors_agree_on_every_rule(self):
        assert all(
            rule
            == E.from_wolfram_number(rule.rule)
            == E.from_table(rule.table)
            == E.from_bits(rule.bits)
            == E.from_bits([int(digit) for digit in rule.bits])
            == E.from_anf(rule.anf)
            == E.from_local_rule(rule.local_rule)
            == E.from_de_bruijn(rule.de_bruijn_graph())
            == E.from_cellular_automaton(rule.to_cellular_automaton())
            for rule in ALL
        )

    def test_conversion_to_the_general_class_keeps_the_wolfram_number(self):
        assert all(
            rule.to_cellular_automaton() == ca.CellularAutomaton.from_wolfram_number(n)
            and rule.to_cellular_automaton().wolfram_number == n
            for n, rule in enumerate(ALL)
        )

    @given(rules)
    def test_conversion_accepts_the_minimal_neighborhood(self, rule):
        assert E.from_cellular_automaton(rule.to_cellular_automaton().minimal()) == rule

    @given(ring_configs())
    def test_cells_round_trip(self, ring):
        config, cells = ring
        digits = eca.cells_from_config(config, cells)
        assert digits == word(format(config, f"0{cells}b"))
        assert eca.config_from_cells(digits) == config

    def test_repr_is_compact(self):
        assert repr(eca.rule_110()) == "ElementaryCellularAutomaton(rule=110)"


# --------------------------------------------------------------------------- #
# Rejections
# --------------------------------------------------------------------------- #
class TestRejections:
    @pytest.mark.parametrize("number", [-1, 256])
    def test_e2_rejects_numbers_outside_0_255(self, number):
        with pytest.raises(ValueError, match=r"\(E2\)"):
            E(number)
        with pytest.raises(ValueError, match=r"\(E2\)"):
            E.from_wolfram_number(number)

    @pytest.mark.parametrize(
        "table", [(0, 1) * 3, (0, 1) * 5, (0, 1, 2, 0, 0, 0, 0, 0)]
    )
    def test_e1_rejects_a_table_that_is_not_eight_bits(self, table):
        with pytest.raises(ValueError, match=r"\(E1\)"):
            E.from_table(table)

    def test_e1_rejects_a_local_rule_with_values_outside_0_1(self):
        with pytest.raises(ValueError, match=r"\(E1\)"):
            E.from_local_rule(lambda p, q, r: p + q + r)

    def test_e1_rejects_arguments_that_are_not_bits(self):
        with pytest.raises(ValueError, match=r"\(E1\)"):
            eca.rule_30().local_rule(0, 2, 0)

    @pytest.mark.parametrize("monomials", [("p", "s"), ("p", "p"), ("qp",)])
    def test_e3_rejects_unknown_or_repeated_monomials(self, monomials):
        with pytest.raises(ValueError, match=r"\(E3\)"):
            E.from_anf(monomials)

    @pytest.mark.parametrize("bits", ["0110111", "011011100", "0110111x", [0] * 7])
    def test_e4_rejects_a_sequence_that_is_not_eight_binary_digits(self, bits):
        with pytest.raises(ValueError, match=r"\(E4\)"):
            E.from_bits(bits)

    def test_e5_rejects_a_missing_edge(self):
        edges = eca.rule_110().de_bruijn_graph()[:-1]
        with pytest.raises(ValueError, match=r"\(E5\)"):
            E.from_de_bruijn(edges)

    def test_e5_rejects_a_repeated_edge(self):
        edges = eca.rule_110().de_bruijn_graph()
        with pytest.raises(ValueError, match=r"\(E5\)"):
            E.from_de_bruijn([*edges, edges[0]])

    def test_e5_rejects_an_edge_that_is_not_pq_to_qr(self):
        edges = list(eca.rule_110().de_bruijn_graph())
        edges[0] = ((0, 0), (1, 0), 0)
        with pytest.raises(ValueError, match=r"\(E5\)"):
            E.from_de_bruijn(edges)

    def test_e5_rejects_a_label_outside_0_1(self):
        edges = list(eca.rule_110().de_bruijn_graph())
        edges[0] = (edges[0][0], edges[0][1], 2)
        with pytest.raises(ValueError, match=r"\(E5\)"):
            E.from_de_bruijn(edges)

    @pytest.mark.parametrize(
        "automaton",
        [
            ca.game_of_life,
            lambda: ca.CellularAutomaton.from_linear(1, 3, {(-1,): 1, (1,): 1}),
            lambda: eca.rule_30().compose(eca.rule_30()),
            lambda: ca.CellularAutomaton(1, 2, ((2,),), (0, 1)),
        ],
    )
    def test_e6_rejects_automata_without_the_elementary_memory_set(self, automaton):
        with pytest.raises(ValueError, match=r"\(E6\)"):
            E.from_cellular_automaton(automaton())

    @pytest.mark.parametrize(("config", "cells"), [(0, 0), (8, 3), (-1, 3)])
    def test_e7_rejects_a_configuration_off_the_ring(self, config, cells):
        with pytest.raises(ValueError, match=r"\(E7\)"):
            eca.rule_30().step(config, cells)

    def test_e7_rejects_cells_that_are_not_bits(self):
        with pytest.raises(ValueError, match=r"\(E7\)"):
            eca.config_from_cells((0, 2, 1))

    @pytest.mark.parametrize(
        "method", ["phase_space", "is_bijective_on_ring", "conserves_number_on_ring"]
    )
    def test_e7_rejects_an_empty_ring(self, method):
        with pytest.raises(ValueError, match=r"\(E7\)"):
            getattr(eca.rule_30(), method)(0)

    def test_negative_steps_are_rejected(self):
        with pytest.raises(ValueError, match="non-negative"):
            eca.rule_30().evolve(0, 3, -1)
        with pytest.raises(ValueError, match="non-negative"):
            eca.rule_30().evolve_from_single_cell(-1)
        with pytest.raises(ValueError, match="non-negative"):
            eca.rule_30().center_column(-1)

    def test_pi_and_upsilon_are_defined_for_additive_rules_only(self):
        with pytest.raises(ValueError, match="additive"):
            eca.rule_30().cycle_length_from_seed(5)
        with pytest.raises(ValueError, match="additive"):
            eca.rule_30().transient_from_seed(5)

    def test_dipolynomial_is_defined_for_additive_rules_only(self):
        with pytest.raises(ValueError, match="additive"):
            _ = eca.rule_110().dipolynomial

    def test_sideways_evolution_needs_permutivity_on_that_side(self):
        with pytest.raises(ValueError, match="leftmost permutive"):
            eca.rule_102().solve_left(0, 0, 0)
        with pytest.raises(ValueError, match="rightmost permutive"):
            eca.rule_30().solve_right(0, 0, 0)

    def test_density_classifier_needs_two_cells(self):
        with pytest.raises(ValueError, match="at least two cells"):
            eca.density_classifier(0, 1)

    def test_suborder_and_two_adic_part_reject_bad_arguments(self):
        with pytest.raises(ValueError, match="odd positive"):
            eca.multiplicative_suborder(6)
        with pytest.raises(ValueError, match="positive"):
            eca.two_adic_part(0)

    def test_dataframe_missing_column_is_rejected(self):
        frame = eca.rule_30().to_dataframe().drop(columns="value")
        with pytest.raises(ValueError, match="missing columns"):
            E.from_dataframe(frame)

    def test_dataframe_with_a_missing_neighborhood_is_rejected(self):
        frame = eca.rule_30().to_dataframe().iloc[:-1]
        with pytest.raises(ValueError, match=r"\(E1\)"):
            E.from_dataframe(frame)

    def test_dataframe_with_a_repeated_neighborhood_is_rejected(self):
        frame = eca.rule_30().to_dataframe()
        with pytest.raises(ValueError, match="repeats"):
            E.from_dataframe(pd.concat([frame, frame.iloc[:1]]))

    def test_dataframe_with_a_non_bit_value_is_rejected(self):
        frame = eca.rule_30().to_dataframe()
        frame.loc[0, "value"] = 2
        with pytest.raises(ValueError, match=r"\(E1\)"):
            E.from_dataframe(frame)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
class TestPlots:
    def test_space_time_diagram_has_one_row_per_time_step(self):
        figure, ax = plt.subplots()
        try:
            rule = eca.rule_90()
            assert rule.plot_space_time(0b00000100000, 11, 7, ax) is ax
            (image,) = ax.images
            array = image.get_array()
            assert array is not None
            assert array.shape == (8, 11)
            assert array[1].tolist() == list(
                eca.cells_from_config(rule.step(0b00000100000, 11), 11)
            )
        finally:
            plt.close(figure)

    def test_single_cell_picture_is_the_widening_triangle(self):
        figure, ax = plt.subplots()
        try:
            rule = eca.rule_30()
            assert rule.plot_single_cell(6, ax) is ax
            (image,) = ax.images
            array = image.get_array()
            assert array is not None
            assert array.shape == (7, 13)
            assert [int(array[t][6]) for t in range(7)] == list(rule.center_column(6))
        finally:
            plt.close(figure)

    def test_rule_icon_reads_as_the_bit_string(self):
        figure, ax = plt.subplots()
        try:
            assert eca.rule_110().plot_rule_icon(ax) is ax
            black = [
                patch for patch in ax.patches if patch.get_facecolor()[:3] == (0, 0, 0)
            ]
            assert len(ax.patches) == 32
            assert len(black) == 12 + 5
        finally:
            plt.close(figure)

    def test_de_bruijn_picture_labels_every_edge(self):
        figure, ax = plt.subplots()
        try:
            rule = eca.rule_110()
            assert rule.plot_de_bruijn(ax) is ax
            texts = Counter(text.get_text() for text in ax.texts)
            assert all(
                texts[f"{p}{q}{r} -> {rule.local_rule(p, q, r)}"] == 1
                for p, q, r in itertools.product(range(2), repeat=3)
            )
        finally:
            plt.close(figure)

    @pytest.mark.parametrize(
        ("method", "args"),
        [
            ("plot_space_time", (0b0100, 4, 3)),
            ("plot_single_cell", (4,)),
            ("plot_rule_icon", ()),
            ("plot_de_bruijn", ()),
        ],
    )
    def test_plots_create_axes_when_none_are_given(self, method, args):
        ax = getattr(eca.rule_30(), method)(*args)
        try:
            assert ax.figure is not None
        finally:
            plt.close(ax.figure)
