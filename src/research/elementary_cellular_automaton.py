r"""Elementary cellular automata: the 256 rules, their symmetries and ring dynamics.

An elementary cellular automaton is stored as its Wolfram number — definition
(E2) of the Elementary Cellular Automaton page (Wolfram 1983, Fig. 1) — since
every other formulation on the page is the same eight bits,
``t_(4p+2q+r) = f(p, q, r)``: the rule table (E1), the algebraic normal form
(E3), the vertex of the 8-cube (E4) and the labeled de Bruijn graph (E5) are
constructors onto the number and accessors off it. The topological form (E6)
is the conversion to and from :class:`research.cellular_automaton.CellularAutomaton`,
which also carries every line-level decision — surjectivity, injectivity,
preimage counts, orphans, diamonds. They are delegated, never re-derived on
a ring: rules 30 and 90 are surjective on the line and bijective on no small
ring, and rule 150 is bijective only when 3 does not divide ``N``.

What this module adds to the general class is what is special to the
family: the Klein four-group of reflection and conjugation acting on rule
numbers, permutivity and the classifications that hang on it (Cattaneo-
Finelli-Margara 2000; Schule-Stoop 2012), the published class tables, and
bitmask dynamics on rings (E7). A ring configuration is an ``int`` below
``2^N`` read as ``N`` binary digits with the leftmost cell most significant —
the numbering of :class:`research.cellular_automaton.PhaseSpace` — so a step
is two rotations and one AND/XOR per monomial of the normal form, which is
what makes the Martin-Odlyzko-Wolfram tables reachable.
"""

import enum
import functools
import itertools
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING

import pandas as pd

from research._plot import ensure_axes
from research.cellular_automaton import (
    CellularAutomaton,
    Edge,
    Orbit,
    PhaseSpace,
    Word,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "ElementaryCellularAutomaton",
    "KurkaClass",
    "LiPackardClass",
    "WolframClass",
    "all_rules",
    "cells_from_config",
    "config_from_cells",
    "density_classifier",
    "multiplicative_suborder",
    "representatives",
    "reversible_rules",
    "rule_18",
    "rule_22",
    "rule_30",
    "rule_45",
    "rule_54",
    "rule_60",
    "rule_73",
    "rule_90",
    "rule_102",
    "rule_110",
    "rule_110_ether",
    "rule_126",
    "rule_128",
    "rule_136",
    "rule_150",
    "rule_180",
    "rule_184",
    "rule_226",
    "rule_232",
    "two_adic_part",
]

_RULES = 256
_NEIGHBORHOODS = 8
_FULL_TABLE = _RULES - 1
_LEFT, _CENTER, _RIGHT = 4, 2, 1
_OFFSETS = ((-1,), (0,), (1,))
_BITS = (0, 1)
_BINARY = 2

_MONOMIALS = ("1", "p", "q", "r", "pq", "pr", "qr", "pqr")
"""The monomials of the normal form (E3), in the page's order."""

_MONOMIAL_MASKS = (0, _LEFT, _CENTER, _RIGHT, 6, 5, 3, 7)
"""The variables of each monomial as a mask ``4[p] + 2[q] + [r]``."""

_COLUMNS = ("p", "q", "r", "value")


class KurkaClass(enum.Enum):
    """Kurka's equicontinuity classes (K1)-(K4) (Schule-Stoop 2012, section III)."""

    EQUICONTINUOUS = "K1"
    ALMOST_EQUICONTINUOUS = "K2"
    SENSITIVE = "K3"
    POSITIVELY_EXPANSIVE = "K4"


class WolframClass(enum.Enum):
    """Wolfram's empirical classes (W1)-(W4) (Wolfram 1984; Kari 2005, section 2.5)."""

    I = "W1"  # noqa: E741 - the class is written with a Roman numeral
    II = "W2"
    III = "W3"
    IV = "W4"


class LiPackardClass(enum.Enum):
    """The five classes of Li and Packard (1990, section 3)."""

    NULL = "null"
    FIXED_POINT = "fixed point"
    PERIODIC = "periodic"
    LOCALLY_CHAOTIC = "locally chaotic"
    CHAOTIC = "chaotic"


# The three tables are keyed by the smallest rule number of each of the 88
# equivalence classes, as printed on the page.
_KURKA: dict[KurkaClass, tuple[int, ...]] = {
    KurkaClass.EQUICONTINUOUS: (
        0, 1, 4, 5, 8, 12, 19, 29, 36, 51, 72, 76, 108, 200, 204,
    ),
    KurkaClass.ALMOST_EQUICONTINUOUS: (
        13, 23, 28, 32, 33, 40, 44, 50, 73, 77, 78, 94, 104, 128, 132, 136, 140,
        156, 160, 164, 168, 172, 178, 232,
    ),
    KurkaClass.SENSITIVE: (
        2, 3, 6, 7, 9, 10, 11, 14, 15, 18, 22, 24, 25, 26, 27, 30, 34, 35, 37, 38,
        41, 42, 43, 45, 46, 54, 56, 57, 58, 60, 62, 74, 106, 110, 122, 126, 130,
        134, 138, 142, 146, 152, 154, 162, 170, 184,
    ),
    KurkaClass.POSITIVELY_EXPANSIVE: (90, 105, 150),
}  # fmt: skip

_WOLFRAM: dict[WolframClass, tuple[int, ...]] = {
    WolframClass.I: (0, 8, 32, 40, 128, 136, 160, 168),
    WolframClass.III: (18, 22, 30, 45, 60, 90, 105, 122, 126, 146, 150),
    WolframClass.IV: (41, 54, 106, 110),
}

_LI_PACKARD: dict[LiPackardClass, tuple[int, ...]] = {
    LiPackardClass.NULL: (0, 8, 32, 40, 128, 136, 160, 168),
    LiPackardClass.FIXED_POINT: (
        2, 4, 10, 12, 13, 24, 34, 36, 42, 44, 46, 56, 57, 58, 72, 76, 77, 78, 104,
        130, 132, 138, 140, 152, 162, 164, 170, 172, 184, 200, 204, 232,
    ),
    LiPackardClass.PERIODIC: (
        1, 3, 5, 6, 7, 9, 11, 14, 15, 19, 23, 25, 27, 28, 29, 33, 35, 37, 38, 41,
        43, 50, 51, 62, 74, 94, 108, 134, 142, 156, 178,
    ),
    LiPackardClass.LOCALLY_CHAOTIC: (26, 73, 154),
    LiPackardClass.CHAOTIC: (
        18, 22, 30, 45, 54, 60, 90, 105, 106, 110, 122, 126, 146, 150,
    ),
}  # fmt: skip


def _by_representative[C: enum.Enum](table: dict[C, tuple[int, ...]]) -> dict[int, C]:
    return {rule: label for label, rules in table.items() for rule in rules}


_KURKA_OF = _by_representative(_KURKA)
_WOLFRAM_OF = _by_representative(_WOLFRAM)
_LI_PACKARD_OF = _by_representative(_LI_PACKARD)


# --------------------------------------------------------------------------- #
# Ring configurations
# --------------------------------------------------------------------------- #
def _require_ring(cells: int) -> None:
    if cells < 1:
        msg = f"(E7) ring violated: a ring has N >= 1 cells; got {cells}"
        raise ValueError(msg)


def _require_config(config: int, cells: int) -> None:
    _require_ring(cells)
    if not 0 <= config < 1 << cells:
        msg = (
            f"(E7) ring violated: a configuration of {cells} cells is an integer "
            f"in 0..2^{cells} - 1; got {config}"
        )
        raise ValueError(msg)


def _require_steps(steps: int) -> None:
    if steps < 0:
        msg = f"steps must be non-negative; got {steps}"
        raise ValueError(msg)


def config_from_cells(cells: Sequence[int]) -> int:
    """Return the ring configuration with these cells, leftmost most significant.

    The inverse of :func:`cells_from_config`; the result is also the number
    of the configuration in :class:`research.cellular_automaton.PhaseSpace`.

    Raises:
        ValueError: If a cell is not 0 or 1.
    """
    if any(value not in _BITS for value in cells):
        msg = f"(E7) ring violated: cells take the values 0 and 1; got {cells!r}"
        raise ValueError(msg)
    config = 0
    for value in cells:
        config = config << 1 | value
    return config


def cells_from_config(config: int, cells: int) -> tuple[int, ...]:
    """Return the ``cells`` cells of a ring configuration, leftmost first.

    Raises:
        ValueError: If ``config`` is not a configuration of that ring.
    """
    _require_config(config, cells)
    return tuple(config >> shift & 1 for shift in range(cells - 1, -1, -1))


def two_adic_part(n: int) -> int:
    """Return ``D_2(n)``, the largest power of 2 dividing ``n``.

    Martin-Odlyzko-Wolfram 1984, section 3.

    Raises:
        ValueError: If ``n`` is not positive.
    """
    if n < 1:
        msg = f"D_2(n) is defined for positive n; got {n}"
        raise ValueError(msg)
    return n & -n


def multiplicative_suborder(n: int) -> int:
    r"""Return ``sord_n(2)``, the least ``j >= 1`` with ``2^j = \pm 1 (mod n)``.

    Martin-Odlyzko-Wolfram 1984, section 3 (the multiplicative suborder).

    Raises:
        ValueError: If ``n`` is not an odd positive integer; 2 is then not
            a unit modulo ``n`` and no such ``j`` exists.
    """
    if n < 1 or n % _BINARY == 0:
        msg = f"sord_n(2) is defined for odd positive n; got {n}"
        raise ValueError(msg)
    power, j = _BINARY % n, 1
    while power not in {1 % n, (n - 1) % n}:
        power, j = power * _BINARY % n, j + 1
    return j


# --------------------------------------------------------------------------- #
# The automaton
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ElementaryCellularAutomaton:
    r"""An elementary cellular automaton, stored as its Wolfram number.

    Definition (E2) of the Elementary Cellular Automaton page (Wolfram 1983,
    Fig. 1; Kari 2005, section 2.5): ``rule`` is
    ``R = sum f(p, q, r) 2^(4p + 2q + r)`` for the local rule
    ``f : {0,1}^3 -> {0,1}`` on (left neighbor, cell, right neighbor), with
    global map ``G_f(x)_i = f(x_(i-1), x_i, x_(i+1))``. The map ``f -> R``
    is a bijection onto ``0..255``.

    The definition is one range check, so ``__post_init__`` makes it and
    direct construction is validated too; the ``from_*`` classmethods
    convert the page's other formulations (E1), (E3)-(E6). Equality is
    equality of rule numbers, which is equality of global maps (page, (E6)).

    Raises:
        ValueError: Naming (E2) when ``rule`` is outside ``0..255``.
    """

    rule: int

    def __post_init__(self) -> None:
        """Validate (E2)."""
        if not 0 <= self.rule < _RULES:
            msg = (
                f"(E2) Wolfram number violated: rule numbers run over 0..255; "
                f"got {self.rule}"
            )
            raise ValueError(msg)

    # -------------------------------------------------------- constructors
    @classmethod
    def from_wolfram_number(cls, rule: int) -> ElementaryCellularAutomaton:
        """Build the rule with this Wolfram number, formulation (E2).

        Raises:
            ValueError: Naming (E2) when ``rule`` is outside ``0..255``.
        """
        return cls(rule)

    @classmethod
    def from_table(cls, table: Sequence[int]) -> ElementaryCellularAutomaton:
        """Build the rule with outputs ``table[4p + 2q + r] = f(p, q, r)`` (E1).

        The order ``(t_0, ..., t_7)`` is that of
        :attr:`research.cellular_automaton.CellularAutomaton.table`; use
        :meth:`from_bits` for the printed order ``t_7 ... t_0``.

        Raises:
            ValueError: Naming (E1) unless there are eight outputs in ``{0, 1}``.
        """
        outputs = tuple(table)
        if len(outputs) != _NEIGHBORHOODS or any(v not in _BITS for v in outputs):
            msg = (
                f"(E1) rule table violated: f : {{0,1}}^3 -> {{0,1}} has eight "
                f"outputs in {{0, 1}}; got {outputs!r}"
            )
            raise ValueError(msg)
        return cls(sum(value << index for index, value in enumerate(outputs)))

    @classmethod
    def from_local_rule(
        cls, rule: Callable[[int, int, int], int]
    ) -> ElementaryCellularAutomaton:
        """Build the rule by tabulating ``f(p, q, r)``, formulation (E1).

        Raises:
            ValueError: Naming (E1) if ``rule`` returns anything but 0 or 1.
        """
        return cls.from_table(
            [rule(p, q, r) for p, q, r in itertools.product(_BITS, repeat=3)]
        )

    @classmethod
    def from_bits(cls, bits: str | Sequence[int]) -> ElementaryCellularAutomaton:
        """Build the rule from its table in printed order ``t_7 ... t_0`` (E4).

        Li-Packard 1990, sections 1-2: the vertex of the 8-cube; also the
        eight-bit sequence ``f(111) f(110) ... f(000)`` of Kari 2005,
        section 2.5. ``"01101110"`` is rule 110.

        Raises:
            ValueError: Naming (E4) unless there are eight binary digits.
        """
        digits = tuple(bits)
        if len(digits) != _NEIGHBORHOODS or any(
            digit not in {0, 1, "0", "1"} for digit in digits
        ):
            msg = (
                f"(E4) vertex of the 8-cube violated: a rule table is a binary "
                f"sequence t_7 ... t_0 of length eight; got {bits!r}"
            )
            raise ValueError(msg)
        return cls.from_table([int(str(digit)) for digit in reversed(digits)])

    @classmethod
    def from_anf(cls, monomials: Iterable[str]) -> ElementaryCellularAutomaton:
        """Build the rule with this algebraic normal form over ``Z_2`` (E3).

        ``monomials`` lists the monomials with coefficient 1, each one of
        ``"1", "p", "q", "r", "pq", "pr", "qr", "pqr"``; the empty list is
        rule 0. The truth table is the Mobius transform of the coefficients
        over the subsets of ``{p, q, r}`` (page, "How the formulations
        connect").

        Raises:
            ValueError: Naming (E3) on an unknown or a repeated monomial.
        """
        chosen = tuple(monomials)
        unknown = [monomial for monomial in chosen if monomial not in _MONOMIALS]
        if unknown or len(set(chosen)) != len(chosen):
            msg = (
                f"(E3) Boolean polynomial violated: the monomials are distinct "
                f"members of {_MONOMIALS!r}; got {chosen!r}"
            )
            raise ValueError(msg)
        masks = [_MONOMIAL_MASKS[_MONOMIALS.index(monomial)] for monomial in chosen]
        return cls.from_table(
            [
                sum(mask & ~index == 0 for mask in masks) % _BINARY
                for index in range(_NEIGHBORHOODS)
            ]
        )

    @classmethod
    def from_de_bruijn(cls, edges: Iterable[Edge]) -> ElementaryCellularAutomaton:
        """Build the rule from a 0/1 labeling of the eight de Bruijn edges (E5).

        Kari's lecture notes, section 2.6, width 3: the edge ``pqr`` goes
        from ``pq`` to ``qr`` and carries ``f(p, q, r)``. The graph
        determines ``G_f`` only up to a power of the shift; the rule
        returned is the one with neighborhood ``(-1, 0, 1)``.

        Raises:
            ValueError: Naming (E5) unless ``edges`` are exactly the eight
                edges ``pq -> qr``, each labeled 0 or 1.
        """
        listed = [(tuple(s), tuple(t), label) for s, t, label in edges]
        expected: dict[tuple[Word, Word], int] = {
            ((p, q), (q, r)): _LEFT * p + _CENTER * q + r
            for p, q, r in itertools.product(_BITS, repeat=3)
        }
        labels = {(source, target): label for source, target, label in listed}
        if (
            len(listed) == len(expected)
            and set(labels) == set(expected)
            and all(label in _BITS for label in labels.values())
        ):
            return cls(sum(labels[edge] << index for edge, index in expected.items()))
        msg = (
            f"(E5) labeled de Bruijn graph violated: a rule labels each of the "
            f"eight edges pqr : pq -> qr with 0 or 1; got {listed!r}"
        )
        raise ValueError(msg)

    @classmethod
    def from_cellular_automaton(
        cls, automaton: CellularAutomaton
    ) -> ElementaryCellularAutomaton:
        """Build the rule with the same global map as ``automaton`` (E6).

        Page, (E6): the elementary automata are exactly the cellular
        automata on ``{0,1}^Z`` admitting ``{-1, 0, 1}`` as a memory set, so
        any neighborhood is accepted as long as the minimal memory set fits.

        Raises:
            ValueError: Naming (E6) if the automaton is not one-dimensional
                with two states, or depends on a cell outside ``{-1, 0, 1}``.
        """
        if automaton.dimension == 1 and automaton.states == _BINARY:
            minimal = automaton.minimal()
            if set(minimal.offsets) <= set(_OFFSETS):
                return cls(minimal.with_neighborhood(_OFFSETS).wolfram_number)
        msg = (
            f"(E6) topological form violated: an elementary automaton acts on "
            f"{{0,1}}^Z with memory set {{-1, 0, 1}}; got {automaton!r}"
        )
        raise ValueError(msg)

    # -------------------------------------------------------- the local rule
    @property
    def table(self) -> tuple[int, ...]:
        """The outputs ``(t_0, ..., t_7)`` with ``t_(4p+2q+r) = f(p, q, r)`` (E1)."""
        return tuple(self.rule >> index & 1 for index in range(_NEIGHBORHOODS))

    @property
    def bits(self) -> str:
        """The table in printed order ``t_7 ... t_0``: the binary expansion (E4).

        Wolfram 2002, p. 865: ``IntegerDigits[n, 2, 8]``.
        """
        return f"{self.rule:08b}"

    def local_rule(self, p: int, q: int, r: int) -> int:
        """Return ``f(p, q, r) = floor(R / 2^(4p + 2q + r)) mod 2`` (E2).

        Raises:
            ValueError: If an argument is not 0 or 1.
        """
        if any(value not in _BITS for value in (p, q, r)):
            msg = f"(E1) rule table violated: f takes bits; got {(p, q, r)!r}"
            raise ValueError(msg)
        return self.rule >> (_LEFT * p + _CENTER * q + r) & 1

    @functools.cached_property
    def _anf_masks(self) -> tuple[int, ...]:
        """The variable masks of the monomials with coefficient 1, page order."""
        return tuple(
            mask
            for mask in _MONOMIAL_MASKS
            if sum(self.rule >> sub & 1 for sub in range(mask + 1) if sub & ~mask == 0)
            % _BINARY
        )

    @property
    def anf(self) -> tuple[str, ...]:
        """The monomials of the algebraic normal form over ``Z_2`` (E3).

        In the order ``1, p, q, r, pq, pr, qr, pqr``; rule 30 gives
        ``("p", "q", "r", "qr")`` (Wolfram 2002, p. 869). The coefficient
        of a monomial is the sum over its sub-assignments of the table — the
        Mobius transform the page cites under "How the formulations connect".
        """
        return tuple(
            _MONOMIALS[_MONOMIAL_MASKS.index(mask)] for mask in self._anf_masks
        )

    @property
    def density(self) -> Fraction:
        """Langton's ``lambda``: "the density of 1's in the rule table".

        Li-Packard 1990, after Langton.
        """
        return Fraction(self.rule.bit_count(), _NEIGHBORHOODS)

    @property
    def is_balanced(self) -> bool:
        """Whether the table has four 1s, ``lambda = 1/2`` (page, "Balanced rule").

        Necessary but not sufficient for surjectivity: rules 184 and 232.
        """
        return self.rule.bit_count() * _BINARY == _NEIGHBORHOODS

    @property
    def hot_bits(self) -> tuple[int, int]:
        """The hot bits ``(t_0, t_7)``, the images of 000 and 111.

        Li-Packard 1990, section 5.
        """
        return self.rule & 1, self.rule >> (_NEIGHBORHOODS - 1)

    @property
    def mean_field_cluster(self) -> tuple[int, int, int, int]:
        """The mean-field cluster ``[n_0 n_1 n_2 n_3]`` (Gutowitz's notation).

        Li-Packard 1990, section 5: ``n_i`` is the number of blocks with
        ``i`` ones that the rule maps to 1.
        """
        counts = [0, 0, 0, 0]
        for index in range(_NEIGHBORHOODS):
            counts[index.bit_count()] += self.rule >> index & 1
        return counts[0], counts[1], counts[2], counts[3]

    # -------------------------------------------------------- rule predicates
    def _depends_on(self, variable: int) -> bool:
        return any(
            (self.rule >> index ^ self.rule >> (index | variable)) & 1
            for index in range(_NEIGHBORHOODS)
            if not index & variable
        )

    def _permutive_in(self, variable: int) -> bool:
        return all(
            (self.rule >> index ^ self.rule >> (index | variable)) & 1
            for index in range(_NEIGHBORHOODS)
            if not index & variable
        )

    @property
    def is_legal(self) -> bool:
        """Whether ``f(0,0,0) = 0`` and the rule equals its reflection.

        Wolfram 1983, section 2 (Wolfram 1984, conditions (2.4), (2.5));
        there are 32 legal rules.
        """
        return self.rule & 1 == 0 and self.reflect() == self

    @property
    def is_additive(self) -> bool:
        """Whether ``f = a p + b q + c r mod 2``, the superposition principle.

        Wolfram 1983, section 2; Martin-Odlyzko-Wolfram 1984. The eight
        additive rules are 0, 60, 90, 102, 150, 170, 204, 240 (Wolfram 2002,
        p. 952).
        """
        return all(mask in {_LEFT, _CENTER, _RIGHT} for mask in self._anf_masks)

    @property
    def is_affine(self) -> bool:
        """Whether the rule is additive or the complement of an additive rule.

        The page's term (after Cattaneo-Finelli-Margara 2000, Table 4):
        rule 105 is affine, not additive.
        """
        return all(mask in {0, _LEFT, _CENTER, _RIGHT} for mask in self._anf_masks)

    @property
    def dipolynomial(self) -> tuple[int, ...]:
        """The exponents of ``T(x)``, descending, for an additive rule.

        Martin-Odlyzko-Wolfram 1984, section 2: with a ring configuration
        written ``A(x) = sum a_i x^i``, one step is
        ``A -> T(x) A(x) mod (x^N - 1)`` and multiplication by ``x`` shifts
        one site to the right, so ``p`` contributes ``x``, ``q``
        contributes 1 and ``r`` contributes ``x^-1``. Rule 90 gives
        ``(1, -1)`` and rule 150 ``(1, 0, -1)``.

        Raises:
            ValueError: If the rule is not additive.
        """
        if not self.is_additive:
            msg = f"only an additive rule has a dipolynomial T(x); got {self!r}"
            raise ValueError(msg)
        exponent = {_LEFT: 1, _CENTER: 0, _RIGHT: -1}
        return tuple(exponent[mask] for mask in self._anf_masks)

    @property
    def is_totalistic(self) -> bool:
        """Whether ``f`` depends only on ``p + q + r`` (Wolfram 1984, eq. (2.3))."""
        return all(
            len({self.rule >> index & 1 for index in indices}) == 1
            for indices in ((1, 2, 4), (3, 5, 6))
        )

    @property
    def is_peripheral(self) -> bool:
        """Whether the new value ignores the cell's own old value (Wolfram 1983)."""
        return not self._depends_on(_CENTER)

    @property
    def is_trivial(self) -> bool:
        """Whether ``f(p, q, r) = g(q)`` depends on the center cell only.

        Cattaneo-Finelli-Margara 2000: rules 0, 51, 204, 255.
        """
        return not self._depends_on(_LEFT) and not self._depends_on(_RIGHT)

    @property
    def is_left_permutive(self) -> bool:
        """Whether ``f(0, q, r) != f(1, q, r)`` for all ``q, r``.

        Cattaneo-Finelli-Margara 2000, Definition 2.2 (after Hedlund), on
        the three-cell table; equivalently ``f = p xor g(q, r)``.
        """
        return self._permutive_in(_LEFT)

    @property
    def is_right_permutive(self) -> bool:
        """Whether ``f(p, q, 0) != f(p, q, 1)`` for all ``p, q``.

        Cattaneo-Finelli-Margara 2000, Definition 2.2.
        """
        return self._permutive_in(_RIGHT)

    @property
    def is_permutive(self) -> bool:
        """Whether the rule is leftmost or rightmost permutive."""
        return self.is_left_permutive or self.is_right_permutive

    @property
    def is_bipermutive(self) -> bool:
        """Whether the rule is leftmost and rightmost permutive: 90, 105, 150, 165."""
        return self.is_left_permutive and self.is_right_permutive

    @property
    def is_number_conserving(self) -> bool:
        """Whether the number of 1s is conserved on every ring.

        Boccara-Fuks, Theorem 2.1 for three inputs: ``f`` is
        number-conserving iff for all ``(x1, x2, x3)``
        ``f(x1,x2,x3) = x1 + [f(0,x2,x3) - f(0,x1,x2)] + [f(0,0,x2) - f(0,0,x1)]``.
        """
        f = self.local_rule
        return all(
            f(x1, x2, x3)
            == x1 + f(0, x2, x3) - f(0, x1, x2) + f(0, 0, x2) - f(0, 0, x1)
            for x1, x2, x3 in itertools.product(_BITS, repeat=3)
        )

    @property
    def is_devaney_chaotic(self) -> bool:
        """Whether ``G_f`` is Devaney-chaotic: iff the rule is permutive.

        Cattaneo-Finelli-Margara 2000, Corollary 3.3 — special to the
        elementary family (their Theorems 5.1, 5.2 give non-permutive
        chaotic rules outside it).
        """
        return self.is_permutive

    @property
    def is_positively_expansive(self) -> bool:
        """Whether ``G_f`` is positively expansive: iff the rule is bipermutive.

        Schule-Stoop 2012, Proposition 13 — not true of one-dimensional
        cellular automata in general.
        """
        return self.is_bipermutive

    # -------------------------------------------------------- published classes
    @property
    def kurka_class(self) -> KurkaClass:
        """The equicontinuity class, from the table of Schule-Stoop 2012.

        Propositions 8, 9, 11, 12, transcribed by smallest representative;
        the class is invariant under reflection and conjugation because the
        global maps are conjugate by an isometry (Cattaneo-Finelli-Margara
        2000, Proposition 2.1). A lookup, not a decision procedure.
        """
        return _KURKA_OF[self.canonical().rule]

    @property
    def wolfram_class(self) -> WolframClass:
        """The empirical Wolfram class as tabulated by Martinez 2013, Table 2.

        An observation about random initial conditions — "It is not
        intended as a rigorous mathematical treatment" (Wolfram 1984) — and
        the sources disagree on class IV; this is Martinez's table.
        """
        return _WOLFRAM_OF.get(self.canonical().rule, WolframClass.II)

    @property
    def li_packard_class(self) -> LiPackardClass:
        """The class of Li-Packard 1990, Table 2.

        Converted on the page to smallest-number representatives.
        """
        return _LI_PACKARD_OF[self.canonical().rule]

    # -------------------------------------------------------- the symmetry group
    def reflect(self) -> ElementaryCellularAutomaton:
        """Return the reflected rule ``rho f(p, q, r) = f(r, q, p)``.

        Li-Packard 1990, section 2: on the table,
        ``(t_7 t_3 t_5 t_1 t_6 t_2 t_4 t_0)``.
        """
        return ElementaryCellularAutomaton.from_local_rule(
            lambda p, q, r: self.local_rule(r, q, p)
        )

    def conjugate(self) -> ElementaryCellularAutomaton:
        """Return the conjugate rule ``gamma f = 1 - f(1-p, 1-q, 1-r)``.

        Li-Packard 1990, section 2; Wolfram 2002, p. 883, as
        ``1 - Reverse[list]``.
        """
        reversed_bits = int(self.bits[::-1], _BINARY)
        return ElementaryCellularAutomaton(_FULL_TABLE ^ reversed_bits)

    def equivalence_class(self) -> frozenset[ElementaryCellularAutomaton]:
        """Return the orbit under the Klein four-group of ``rho`` and ``gamma``.

        One, two or four rules; there are 88 classes (Li-Packard 1990,
        Appendix).
        """
        mirrored = self.reflect()
        return frozenset({self, mirrored, self.conjugate(), mirrored.conjugate()})

    def canonical(self) -> ElementaryCellularAutomaton:
        """Return the class member with the smallest rule number.

        Convention (A) of the page: Wolfram 2002, Martinez, Schule-Stoop.
        """
        return min(self.equivalence_class(), key=lambda member: member.rule)

    def li_packard_representative(self) -> ElementaryCellularAutomaton:
        """Return Li-Packard's representative, convention (B) of the page.

        The smallest number among the members with the smaller ``lambda``;
        it differs from :meth:`canonical` for exactly five classes.
        """
        return min(
            self.equivalence_class(), key=lambda member: (member.density, member.rule)
        )

    def hamming_distance(self, other: ElementaryCellularAutomaton) -> int:
        """Return the number of neighborhoods on which the two tables differ (E4)."""
        return (self.rule ^ other.rule).bit_count()

    def neighbors(self) -> tuple[ElementaryCellularAutomaton, ...]:
        """Return the eight rules at Hamming distance 1, by flipped bit ``t_0..t_7``.

        Li-Packard 1990, sections 1-2: the neighbors in the 8-cube, one
        "mutation" away.
        """
        return tuple(
            ElementaryCellularAutomaton(self.rule ^ 1 << index)
            for index in range(_NEIGHBORHOODS)
        )

    # -------------------------------------------------------- the general class
    @functools.cached_property
    def _automaton(self) -> CellularAutomaton:
        return CellularAutomaton.from_wolfram_number(self.rule)

    def to_cellular_automaton(self) -> CellularAutomaton:
        """Return the quadruple ``(1, {0,1}, (-1, 0, 1), f)``.

        The case of definition (D1) of the Cellular Automaton page that the
        elementary automata specialize.
        """
        return self._automaton

    def compose(self, other: ElementaryCellularAutomaton) -> CellularAutomaton:
        """Return ``G_self o G_other``, applying ``other`` first.

        Wolfram 1984, section 2: composing two elementary rules generally
        gives an automaton of range 2, so the result is a general
        :class:`~research.cellular_automaton.CellularAutomaton` on
        ``(-2, ..., 2)``.
        """
        return self._automaton.compose(other._automaton)

    def second_order(self) -> CellularAutomaton:
        """Return the second-order variant ``S_n = F[S_(n-1)] xor S_(n-2)``.

        Wolfram 1983, section 4 (after Fredkin and Margolus): invertible for
        every elementary rule, and an ordinary automaton on the alphabet
        ``A x A`` rather than an elementary one.
        """
        return self._automaton.second_order()

    @property
    def minimal_memory_set(self) -> frozenset[int]:
        """The cells ``f`` really depends on, a subset of ``{-1, 0, 1}`` (E6).

        ``{1}`` for rule 170, ``{0}`` for rule 204, empty for rule 0.
        """
        return frozenset(
            offset
            for offset, variable in zip(
                (-1, 0, 1), (_LEFT, _CENTER, _RIGHT), strict=True
            )
            if self._depends_on(variable)
        )

    def de_bruijn_graph(self) -> tuple[Edge, ...]:
        """Return the eight labeled edges ``(pq, qr, f(p, q, r))`` (E5).

        Kari's lecture notes, section 2.6.
        """
        return self._automaton.de_bruijn_graph()

    @property
    def is_surjective(self) -> bool:
        """Whether ``G_f`` is surjective on ``{0,1}^Z``.

        Decided on the de Bruijn graph (Sutner; Kari 2005, Theorem 9), never
        on a ring. For an elementary rule this holds iff the rule is
        permutive or is 51 or 204 (Cattaneo-Finelli-Margara 2000,
        Corollary 3.3).
        """
        return self._automaton.is_surjective()

    @property
    def is_injective(self) -> bool:
        """Whether ``G_f`` is injective, by the pair graph.

        Kari's lecture notes, Proposition 26. True of rules 15, 51, 85,
        170, 204, 240 only (Wolfram 2002, p. 436).
        """
        return self._automaton.is_injective()

    @property
    def is_reversible(self) -> bool:
        """Whether ``G_f`` is reversible; equivalent to injective.

        Kari 2005, Corollary 3.
        """
        return self._automaton.is_reversible()

    def apply_to_word(self, word: Sequence[int]) -> Word:
        """Return the image of a word, two symbols shorter.

        Raises:
            ValueError: If ``word`` has fewer than two symbols or a symbol
                that is not a bit.
        """
        return self._automaton.apply_to_word(word)

    def preimage_count(self, word: Sequence[int]) -> int:
        """Return the number of words of length ``len(word) + 2`` mapping to it.

        Exactly 4 for every word iff the rule is surjective
        (Cattaneo-Finelli-Margara 2000, Theorem 3.4, after Hedlund).

        Raises:
            ValueError: If a symbol is not a bit.
        """
        return self._automaton.preimage_count(word)

    def is_orphan(self, word: Sequence[int]) -> bool:
        """Return whether ``word`` has no preimage word (page, "orphan").

        Raises:
            ValueError: If a symbol is not a bit.
        """
        return self._automaton.is_orphan(word)

    def shortest_orphans(self) -> tuple[Word, ...]:
        """Return every orphan of minimal length, sorted; empty when surjective."""
        return self._automaton.shortest_orphans()

    def diamond(self, length: int) -> tuple[Word, Word] | None:
        """Return two words differing in ``length`` interior cells with one image.

        Page, "diamond": distinct words with equal image that begin and end
        alike, here on two border cells each side.

        Raises:
            ValueError: If ``length`` is not positive.
        """
        return self._automaton.diamond(length)

    def solve_left(self, value: int, q: int, r: int) -> int:
        """Return the ``p`` with ``f(p, q, r) = value``, running the rule sideways.

        Wolfram 1986, eq. (3.3): for rule 30, ``p = value xor (q or r)``, so
        the time sequences of two adjacent cells determine the pattern to
        their left.

        Raises:
            ValueError: If the rule is not leftmost permutive, or an
                argument is not a bit.
        """
        if not self.is_left_permutive:
            msg = f"only a leftmost permutive rule can be solved for p; got {self!r}"
            raise ValueError(msg)
        return value ^ self.local_rule(0, q, r)

    def solve_right(self, p: int, q: int, value: int) -> int:
        """Return the ``r`` with ``f(p, q, r) = value``; the mirror of `solve_left`.

        Raises:
            ValueError: If the rule is not rightmost permutive, or an
                argument is not a bit.
        """
        if not self.is_right_permutive:
            msg = f"only a rightmost permutive rule can be solved for r; got {self!r}"
            raise ValueError(msg)
        return value ^ self.local_rule(p, q, 0)

    # -------------------------------------------------------- rings (E7)
    def _step(self, config: int, cells: int) -> int:
        """Apply the rule to a validated ring configuration through its normal form."""
        full = (1 << cells) - 1
        left = config >> 1 | (config & 1) << cells - 1
        right = config << 1 & full | config >> cells - 1
        result = 0
        for mask in self._anf_masks:
            term = full
            if mask & _LEFT:
                term &= left
            if mask & _CENTER:
                term &= config
            if mask & _RIGHT:
                term &= right
            result ^= term
        return result

    def step(self, config: int, cells: int) -> int:
        """Return ``G_f(config)`` on the ring of ``cells`` cells (E7).

        Wolfram 1983, section 4: periodic boundary conditions, "as if the
        sites lay on a circle of circumference N" — the restriction of
        ``G_f`` to configurations of spatial period ``N``. The configuration
        is an integer read as ``N`` binary digits, leftmost cell most
        significant; rings of one and two cells are allowed.

        Raises:
            ValueError: Naming (E7) if ``config`` is not a configuration of
                the ring.
        """
        _require_config(config, cells)
        return self._step(config, cells)

    def evolve(self, config: int, cells: int, steps: int) -> tuple[int, ...]:
        """Return the ``steps + 1`` rows ``c, G(c), ..., G^steps(c)`` on a ring.

        The rows of a space-time diagram, top row first (Kari 2005,
        section 2.5).

        Raises:
            ValueError: If ``steps`` is negative or ``config`` is not a
                configuration of the ring.
        """
        _require_config(config, cells)
        _require_steps(steps)
        rows = [config]
        for _ in range(steps):
            rows.append(self._step(rows[-1], cells))
        return tuple(rows)

    def orbit(self, config: int, cells: int) -> Orbit:
        """Follow a ring trajectory until a configuration repeats.

        Wolfram 1983, section 4: there are only ``2^N`` configurations, so
        ``transient + period <= 2^N``. Every visited configuration is kept
        in memory.

        Raises:
            ValueError: If ``config`` is not a configuration of the ring.
        """
        _require_config(config, cells)
        first_seen: dict[int, int] = {}
        current, time = config, 0
        while current not in first_seen:
            first_seen[current] = time
            current = self._step(current, cells)
            time += 1
        transient = first_seen[current]
        return Orbit(transient=transient, period=time - transient)

    def phase_space(self, cells: int) -> PhaseSpace:
        """Return the phase space on the ring of ``cells`` cells.

        ``2^N`` steps and as much memory; configuration numbers are the
        ring configurations of this module.

        Raises:
            ValueError: If ``cells`` is not positive.
        """
        _require_ring(cells)
        successor = tuple(self._step(config, cells) for config in range(1 << cells))
        return PhaseSpace(_BINARY, (cells,), successor)

    def is_bijective_on_ring(self, cells: int) -> bool:
        """Return whether the rule permutes the configurations of this ring.

        Says nothing about surjectivity on the line (rules 30 and 90), but
        in dimension one a rule bijective on every ring is injective (Kari
        2005, Theorem 7). Stops at the first collision; at most ``2^N``
        steps.

        Raises:
            ValueError: If ``cells`` is not positive.
        """
        _require_ring(cells)
        seen: set[int] = set()
        for config in range(1 << cells):
            image = self._step(config, cells)
            if image in seen:
                return False
            seen.add(image)
        return True

    def conserves_number_on_ring(self, cells: int) -> bool:
        """Return whether every configuration of this ring keeps its number of 1s.

        Boccara-Fuks, Definition 2.1, one ring at a time; rings of 5 cells
        decide it (their Remark 2.1).

        Raises:
            ValueError: If ``cells`` is not positive.
        """
        _require_ring(cells)
        return all(
            self._step(config, cells).bit_count() == config.bit_count()
            for config in range(1 << cells)
        )

    def _seed_orbit(self, cells: int) -> Orbit:
        if not self.is_additive:
            msg = (
                f"Pi_N and Upsilon_N are defined for additive rules "
                f"(Martin-Odlyzko-Wolfram 1984); got {self!r}"
            )
            raise ValueError(msg)
        return self.orbit(1, cells)

    def cycle_length_from_seed(self, cells: int) -> int:
        """Return ``Pi_N``, the length of the cycle reached from a single 1.

        Martin-Odlyzko-Wolfram 1984, section 3: for an additive rule every
        cycle length on the ring divides it (Lemmas 3.4, 4.2).

        Raises:
            ValueError: If the rule is not additive or ``cells`` is not
                positive.
        """
        return self._seed_orbit(cells).period

    def transient_from_seed(self, cells: int) -> int:
        """Return ``Upsilon_N``, the length of the transient from a single 1.

        Wolfram 1983, section 4; Martin-Odlyzko-Wolfram 1984.

        Raises:
            ValueError: If the rule is not additive or ``cells`` is not
                positive.
        """
        return self._seed_orbit(cells).transient

    # -------------------------------------------------------- a single 1 on the line
    def evolve_from_single_cell(self, steps: int) -> tuple[tuple[int, ...], ...]:
        """Return the evolution of a single 1 on the line, as widening rows.

        Row ``t`` lists the ``2t + 1`` cells at positions ``-t, ..., t``, so
        the cell at position ``j`` is ``rows[t][t + j]`` and cells outside
        the row hold ``f^t(0, 0, 0)``. Computed on a ring of ``2 steps + 1``
        cells, which the light cone of the window never wraps around; the
        cost is ``O(steps^2)`` cell updates.

        Raises:
            ValueError: If ``steps`` is negative.
        """
        _require_steps(steps)
        width = 2 * steps + 1
        rows = self.evolve(1 << steps, width, steps)
        return tuple(
            cells_from_config(row, width)[steps - t : steps + t + 1]
            for t, row in enumerate(rows)
        )

    def center_column(self, steps: int) -> tuple[int, ...]:
        """Return the values of cell 0 at times ``0..steps`` from a single 1.

        Page, "Center column"; for rule 30 this is OEIS A051023. The direct
        method, ``O(steps^2)`` cell updates (Wolfram 2019).

        Raises:
            ValueError: If ``steps`` is negative.
        """
        return tuple(
            row[t] for t, row in enumerate(self.evolve_from_single_cell(steps))
        )

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to the rule table: one row per neighborhood.

        Columns ``p``, ``q``, ``r``, ``value`` (all integers), rows in the
        order 000, 001, ..., 111 so that row ``4p + 2q + r`` holds
        ``f(p, q, r)``. The frame always has eight rows and survives
        ``experiments.io.write_result`` (records-oriented JSON).

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        records = [
            (p, q, r, self.local_rule(p, q, r))
            for p, q, r in itertools.product(_BITS, repeat=3)
        ]
        return pd.DataFrame.from_records(records, columns=list(_COLUMNS))

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> ElementaryCellularAutomaton:
        """Rebuild a rule from :meth:`to_dataframe` output; rows in any order.

        Raises:
            ValueError: If a column is missing, a neighborhood is missing or
                repeated, or an entry is not a bit (naming (E1)).
        """
        missing = [column for column in _COLUMNS if column not in df.columns]
        if missing:
            msg = f"elementary-rule frame is missing columns {missing!r}"
            raise ValueError(msg)
        outputs: dict[tuple[int, ...], int] = {}
        for p, q, r, value in zip(df["p"], df["q"], df["r"], df["value"], strict=True):
            key = (int(p), int(q), int(r))
            if key in outputs:
                msg = f"elementary-rule frame repeats the neighborhood {key!r}"
                raise ValueError(msg)
            outputs[key] = int(value)
        neighborhoods = list(itertools.product(_BITS, repeat=3))
        if set(outputs) != set(neighborhoods):
            msg = (
                f"(E1) rule table violated: the frame must list the eight "
                f"neighborhoods in {{0,1}}^3; got {sorted(outputs)!r}"
            )
            raise ValueError(msg)
        return ElementaryCellularAutomaton.from_table(
            [outputs[key] for key in neighborhoods]
        )

    # -------------------------------------------------------- visualization
    def plot_space_time(
        self, config: int, cells: int, steps: int, ax: Axes | None = None
    ) -> Axes:
        """Draw the space-time diagram of a ring onto ``ax``, time downward.

        "Horizontal rows of a space-time diagram are consecutive
        configurations. The top row is the initial configuration." (Kari
        2005, section 2.5.)

        Args:
            config: The initial ring configuration.
            cells: The number of cells of the ring.
            steps: The number of rule applications to draw.
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.

        Raises:
            ValueError: If ``steps`` is negative or ``config`` is not a
                configuration of the ring.
        """
        rows = self.evolve(config, cells, steps)
        ax = ensure_axes(ax)
        ax.imshow(
            [list(cells_from_config(row, cells)) for row in rows],
            cmap="binary",
            interpolation="nearest",
            vmin=0,
            vmax=1,
            aspect="equal",
        )
        ax.set_xlabel("cell")
        ax.set_ylabel("time")
        ax.set_title(f"Rule {self.rule}")
        return ax

    def plot_single_cell(self, steps: int, ax: Axes | None = None) -> Axes:
        """Draw the pattern grown from a single 1 onto ``ax``, time downward.

        The picture by which Wolfram 1983 (section 2) sorts the legal rules
        into simple and complex; position 0 is the center column.

        Args:
            steps: The number of rule applications to draw.
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.

        Raises:
            ValueError: If ``steps`` is negative.
        """
        _require_steps(steps)
        width = 2 * steps + 1
        rows = self.evolve(1 << steps, width, steps)
        ax = ensure_axes(ax)
        ax.imshow(
            [list(cells_from_config(row, width)) for row in rows],
            cmap="binary",
            interpolation="nearest",
            vmin=0,
            vmax=1,
            aspect="equal",
            extent=(-steps - 0.5, steps + 0.5, steps + 0.5, -0.5),
        )
        ax.set_xlabel("position")
        ax.set_ylabel("time")
        ax.set_title(f"Rule {self.rule} from a single 1")
        return ax

    def plot_rule_icon(self, ax: Axes | None = None) -> Axes:
        """Draw the rule icon onto ``ax``: each neighborhood above its output.

        Wolfram 1983, Fig. 1: the eight neighborhoods 111, 110, ..., 000
        from left to right, so the bottom row reads as :attr:`bits`.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        shades = ("white", "black")

        def square(x: float, y: float, value: int) -> None:
            ax.fill(
                [x, x + 1, x + 1, x],
                [y, y, y + 1, y + 1],
                facecolor=shades[value],
                edgecolor="0.5",
            )

        for column, index in enumerate(range(_NEIGHBORHOODS - 1, -1, -1)):
            origin = 4 * column
            for offset, variable in enumerate((_LEFT, _CENTER, _RIGHT)):
                square(origin + offset, 1, int(bool(index & variable)))
            square(origin + 1, 0, self.rule >> index & 1)
        ax.set_xlim(-0.5, 4 * _NEIGHBORHOODS - 0.5)
        ax.set_ylim(-0.5, 2.5)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(f"Rule {self.rule}")
        return ax

    def plot_de_bruijn(self, ax: Axes | None = None) -> Axes:
        """Draw the labeled de Bruijn graph (E5) onto ``ax``.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        return self._automaton.plot_de_bruijn(ax)


# --------------------------------------------------------------------------- #
# The rule space
# --------------------------------------------------------------------------- #
def all_rules() -> tuple[ElementaryCellularAutomaton, ...]:
    """Return the 256 elementary rules in order of Wolfram number."""
    return tuple(ElementaryCellularAutomaton(rule) for rule in range(_RULES))


def representatives() -> tuple[ElementaryCellularAutomaton, ...]:
    """Return the 88 rules that are the smallest member of their class, ascending.

    Li-Packard 1990, Appendix: 8 classes of one rule, 36 of two, 44 of four.
    """
    return tuple(rule for rule in all_rules() if rule.canonical() == rule)


# --------------------------------------------------------------------------- #
# Two-rule schedules
# --------------------------------------------------------------------------- #
def density_classifier(config: int, cells: int) -> int:
    """Return ``G_232^m (G_184^n (config))`` on a ring: Fuks's density classifier.

    Fuks 1997, Proposition 4: with ``n = floor((L - 2) / 2)`` and
    ``m = floor((L - 1) / 2)`` the result consists of only 0s if the density
    of 1s is below 1/2, of only 1s if it is above, and alternates
    ``...0101...`` at exactly 1/2. No single two-state rule does this (Land
    and Belew 1995, per Fuks); the schedule evades that by using two.

    Raises:
        ValueError: If the ring has fewer than two cells (``n`` would be
            negative) or ``config`` is not one of its configurations.
    """
    _require_config(config, cells)
    if cells < _BINARY:
        msg = f"Fuks's schedule needs a ring of at least two cells; got {cells}"
        raise ValueError(msg)
    traffic, majority = rule_184(), rule_232()
    config = traffic.evolve(config, cells, (cells - 2) // 2)[-1]
    return majority.evolve(config, cells, (cells - 1) // 2)[-1]


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
def reversible_rules() -> tuple[ElementaryCellularAutomaton, ...]:
    """Return rules 204, 51, 170, 240, 85, 15: the complete reversible list.

    Identity, complement, left shift, right shift and the complemented
    shifts — ``q``, ``not q``, ``r``, ``p``, ``not r``, ``not p`` (Wolfram
    2002, pp. 436, 883). Rules 204 and 51 are the only surjective rules that
    are not permutive; 204 is equicontinuous while 170 is sensitive.
    """
    return tuple(
        ElementaryCellularAutomaton(rule) for rule in (204, 51, 170, 240, 85, 15)
    )


def rule_18() -> ElementaryCellularAutomaton:
    """Return rule 18, ``not q and (p xor r)``: complex, legal, non-additive.

    Its ring Gardens of Eden are characterized exactly
    (Martin-Odlyzko-Wolfram 1984, Lemma 5.1).
    """
    return ElementaryCellularAutomaton(18)


def rule_22() -> ElementaryCellularAutomaton:
    """Return rule 22, totalistic: shortest orphans of length 8 in a small rule."""
    return ElementaryCellularAutomaton(22)


def rule_30() -> ElementaryCellularAutomaton:
    """Return rule 30, ``p xor (q or r)`` (Wolfram 1986).

    Leftmost permutive and non-additive: surjective on the line but
    bijective on no ring with ``4 <= N <= 16``; the center-column fixture.
    """
    return ElementaryCellularAutomaton(30)


def rule_45() -> ElementaryCellularAutomaton:
    """Return rule 45, ``p xor (q or not r)``, the companion of rule 30.

    Wolfram 1986; bijective on rings exactly for ``N`` odd (page, computed
    for ``N <= 16``).
    """
    return ElementaryCellularAutomaton(45)


def rule_54() -> ElementaryCellularAutomaton:
    """Return rule 54, ``q xor (p or r)``: legal, a candidate for class 4.

    Whether it is universal is open (Kari 2005, Open problem 1); the page
    records that it certifies nothing rigorous yet.
    """
    return ElementaryCellularAutomaton(54)


def rule_60() -> ElementaryCellularAutomaton:
    """Return rule 60, ``p xor q``: additive, permutive on the left only.

    Devaney-chaotic but in class (K3): differences propagate only to the
    right (page, derived from Kari's notes, Example 45).
    """
    return ElementaryCellularAutomaton(60)


def rule_73() -> ElementaryCellularAutomaton:
    """Return rule 73: the class tables are not refinements of one another.

    Locally chaotic for Li-Packard, Wolfram class II in Martinez's table,
    almost equicontinuous (Schule-Stoop 2012, Corollaries 4 and 5).
    """
    return ElementaryCellularAutomaton(73)


def rule_90() -> ElementaryCellularAutomaton:
    """Return rule 90, ``p xor r``: additive, bipermutive, positively expansive.

    From a single 1 it draws Pascal's triangle modulo two (Wolfram 1983);
    surjective on the line while half or three quarters of all ring
    configurations are unreachable (Martin-Odlyzko-Wolfram 1984, section 3).
    """
    return ElementaryCellularAutomaton(90)


def rule_102() -> ElementaryCellularAutomaton:
    """Return rule 102, ``q xor r``: Kari's radius-1/2 xor automaton.

    Kari's lecture notes, Examples 2 and 45: differences only propagate to
    the left.
    """
    return ElementaryCellularAutomaton(102)


def rule_110() -> ElementaryCellularAutomaton:
    """Return rule 110, ``(q or r) and not (p and q and r)``: universal.

    Cook 2004. Not surjective, with shortest orphan 01010 (Kari's lecture
    notes, Example 14); class ``{110, 124, 137, 193}``.
    """
    return ElementaryCellularAutomaton(110)


def rule_110_ether() -> tuple[int, ...]:
    """Return the 14-cell ether block ``10011011111000`` of rule 110.

    Cook 2004, section 3.1; Wolfram 2002, pp. 290, 964: "blocks of 14 cells
    that repeat every 7 steps". With the block written ``b_0 ... b_13`` the
    cell at position ``x`` on step ``t`` is ``b_((x + 4t) mod 14)``.
    """
    return tuple(int(digit) for digit in "10011011111000")


def rule_126() -> ElementaryCellularAutomaton:
    """Return rule 126, Wolfram's typical class 3 rule (2002, p. 250).

    On a ring of 8 cells: the null fixed point, four 6-cycles and two
    2-cycles, with 190 of 256 configurations unreachable (Wolfram 1983,
    section 4).
    """
    return ElementaryCellularAutomaton(126)


def rule_128() -> ElementaryCellularAutomaton:
    """Return rule 128, ``p and q and r``: class 1 is not nilpotent.

    Every finite configuration dies, yet the all-ones configuration is
    fixed (Kari 2005, section 2.5).
    """
    return ElementaryCellularAutomaton(128)


def rule_136() -> ElementaryCellularAutomaton:
    """Return rule 136, ``q and r`` (Kari's lecture notes, Example 41)."""
    return ElementaryCellularAutomaton(136)


def rule_150() -> ElementaryCellularAutomaton:
    """Return rule 150, ``p xor q xor r``: additive, totalistic, self-equivalent.

    Reversibility on a ring can depend on ``N``: bijective iff 3 does not
    divide ``N`` (page, from Martin-Odlyzko-Wolfram 1984, Theorem 4.2).
    """
    return ElementaryCellularAutomaton(150)


def rule_180() -> ElementaryCellularAutomaton:
    """Return rule 180, leftmost permutive (the Coven-Hedlund example).

    Kari's lecture notes, Example 48: composed with the left shift it
    becomes non-sensitive, so the elementary theorems are about the
    neighborhood ``(-1, 0, 1)`` exactly.
    """
    return ElementaryCellularAutomaton(180)


def rule_184() -> ElementaryCellularAutomaton:
    """Return rule 184, the traffic rule: number-conserving, not surjective.

    Fuks 1997; Boccara-Fuks eq. (3). Balanced with the orphan 1100, so a
    balanced table does not imply surjectivity.
    """
    return ElementaryCellularAutomaton(184)


def rule_226() -> ElementaryCellularAutomaton:
    """Return rule 226, the reflection and the conjugate of rule 184.

    It "replaces pattern 01 by pattern 10" (Kari's lecture notes,
    Example 26).
    """
    return ElementaryCellularAutomaton(226)


def rule_232() -> ElementaryCellularAutomaton:
    """Return rule 232, the majority rule: balanced, not surjective, (K2).

    Kari's lecture notes, Examples 6, 42, 44; with rule 184 it solves
    density classification (Fuks 1997).
    """
    return ElementaryCellularAutomaton(232)
