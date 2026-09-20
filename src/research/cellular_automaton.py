r"""Cellular automata on ``Z^d``: local rules, finite dynamics, line-level decisions.

A cellular automaton is stored exactly as Kari's quadruple ``(d, S, N, f)``
— definition (D1) of the Cellular Automaton page (Kari 1994, section 2.1) —
with ``S = {0, ..., k - 1}``: a dimension, a number of states, a neighborhood
vector of distinct elements of ``Z^d`` and the local rule as a lookup table.
The other formulations on the page are constructors onto that data: the
memory-set form (D2) over the group ``Z^d``, the sliding-block form (D4),
Moore's tessellation structures with a quiescent state (D5), Wolfram numbers
for the elementary rules, and linear rules over ``Z_m``. The topological
definition (D3) carries no finite data; it survives as the shift-commutation
law the tests check (Curtis-Hedlund-Lyndon).

Three kinds of configuration are supported, and the page's warning that they
behave differently is the reason they are kept apart:

* *periodic* configurations, as flat row-major tuples on a ring or torus
  (:meth:`CellularAutomaton.step`, :meth:`CellularAutomaton.orbit`,
  :meth:`CellularAutomaton.phase_space`) — Kari's ``G_P`` on a fixed period;
* *finite* configurations on the unbounded lattice, as sparse mappings of
  the non-quiescent cells (:meth:`CellularAutomaton.step_finite`) — ``G_F``;
* *all* configurations of the line, reached only through the de Bruijn and
  pair graphs (:meth:`CellularAutomaton.is_surjective`,
  :meth:`CellularAutomaton.is_injective`,
  :meth:`CellularAutomaton.is_pre_injective`). Surjectivity and injectivity
  of the global map are never decided by exhausting a ring: rules 30 and 90
  are surjective on the line and bijective on no small ring.

Transformations stay inside the type — composition, powers, neighborhood
extension and minimization, the reflection/complement symmetries, the
second-order (Fredkin) construction — because every related structure on the
page (graph dynamical systems, block automata, lattice gases) is still a red
link.
"""

import dataclasses
import functools
import itertools
import math
from collections import Counter, deque
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, override

import pandas as pd

from research._bitmask import bits
from research._plot import ensure_axes, unit_circle

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "CellularAutomaton",
    "Orbit",
    "PhaseSpace",
    "game_of_life",
    "life_beehive",
    "life_blinker",
    "life_block",
    "life_glider",
    "life_r_pentomino",
    "life_row",
    "life_t_tetromino",
    "moore_neighborhood",
    "one_way_nand",
    "one_way_neighborhood",
    "patt_rule",
    "patt_table_on_gapped_neighborhood",
    "radius_half_xor",
    "radius_neighborhood",
    "reversible_elementary_rules",
    "rule_30",
    "rule_90",
    "rule_102",
    "rule_110",
    "rule_150",
    "rule_184",
    "rule_226",
    "rule_232",
    "shift_automaton",
    "translate",
    "von_neumann_neighborhood",
]

type Offset = tuple[int, ...]
"""An element of ``Z^d``: a cell, or a neighborhood offset."""

type Configuration = tuple[int, ...]
"""A periodic configuration: the cells of a ring or torus in row-major order."""

type Word = tuple[int, ...]
"""A pattern on an interval of ``Z``."""

type Edge = tuple[Word, Word, int]
"""A labeled de Bruijn edge ``(source, target, label)``."""

type LocalRule = Callable[[tuple[int, ...]], int]
"""A local rule ``f : S^n -> S`` as a callable on the tuple ``(a_1, ..., a_n)``."""

_BINARY = 2
_PLANE = 2
_ELEMENTARY_OFFSETS: tuple[Offset, ...] = ((-1,), (0,), (1,))
_WOLFRAM_NUMBERS = 256
_REPR_TABLE_LIMIT = 32
_DIGIT_STATES = 10

_DIMENSION = "dimension"
_STATES = "states"
_OFFSET = "offset"
_RULE = "rule"
_COLUMNS = ("kind", "index", "axis", "value")


# --------------------------------------------------------------------------- #
# Indexing helpers
# --------------------------------------------------------------------------- #
def _index(values: Iterable[int], states: int) -> int:
    """Return the base-``states`` number with these digits, first most significant."""
    index = 0
    for value in values:
        index = index * states + value
    return index


def _digits(index: int, states: int, count: int) -> tuple[int, ...]:
    """Return the ``count`` base-``states`` digits of ``index``, inverse of `_index`."""
    digits = [0] * count
    for position in range(count - 1, -1, -1):
        index, digits[position] = divmod(index, states)
    return tuple(digits)


@functools.cache
def _torus_neighbors(
    offsets: tuple[Offset, ...], shape: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Return, for every cell of the torus, the flat indices of its neighbors."""
    strides = [math.prod(shape[axis + 1 :]) for axis in range(len(shape))]
    return tuple(
        tuple(
            sum(
                ((c + o) % size) * stride
                for c, o, size, stride in zip(cell, offset, shape, strides, strict=True)
            )
            for offset in offsets
        )
        for cell in itertools.product(*(range(size) for size in shape))
    )


# --------------------------------------------------------------------------- #
# Standard neighborhoods
# --------------------------------------------------------------------------- #
def _require_dimension(dimension: int) -> None:
    if dimension < 1:
        msg = f"(D1) dimension violated: d must be a positive integer; got {dimension}"
        raise ValueError(msg)


def radius_neighborhood(dimension: int, radius: int) -> tuple[Offset, ...]:
    r"""Return the radius-``r`` neighborhood ``{y : \|y\|_inf <= r}``, sorted.

    Kari 2005, section 2.2.

    Raises:
        ValueError: If ``dimension`` is not positive or ``radius`` is negative.
    """
    _require_dimension(dimension)
    if radius < 0:
        msg = f"a neighborhood radius is non-negative; got {radius}"
        raise ValueError(msg)
    return tuple(itertools.product(range(-radius, radius + 1), repeat=dimension))


def moore_neighborhood(dimension: int) -> tuple[Offset, ...]:
    r"""Return the Moore neighborhood ``{y : \|y\|_inf <= 1}``: ``3^d`` cells, sorted.

    Kari 2005, section 2.2; Moore 1962 took these nine cells in the plane.

    Raises:
        ValueError: If ``dimension`` is not positive.
    """
    return radius_neighborhood(dimension, 1)


def von_neumann_neighborhood(dimension: int) -> tuple[Offset, ...]:
    r"""Return the von Neumann neighborhood ``{y : \|y\|_1 <= 1}``: ``2d + 1`` cells.

    Kari 2005, section 2.2; von Neumann (section 2.1.2) chose the four
    nearest neighbors over the eight.

    Raises:
        ValueError: If ``dimension`` is not positive.
    """
    return tuple(
        offset
        for offset in moore_neighborhood(dimension)
        if sum(abs(coordinate) for coordinate in offset) <= 1
    )


def one_way_neighborhood(dimension: int) -> tuple[Offset, ...]:
    """Return the radius-1/2 neighborhood, every coordinate in ``{0, 1}``, sorted.

    Kari 2005, section 2.2: a one-dimensional radius-1/2 automaton is called
    one-way.

    Raises:
        ValueError: If ``dimension`` is not positive.
    """
    _require_dimension(dimension)
    return tuple(itertools.product(range(2), repeat=dimension))


def translate(
    config: Sequence[int], shape: Sequence[int], by: Sequence[int]
) -> Configuration:
    """Return the torus configuration ``c'(x) = c(x + by)``.

    For ``by = e_i`` this is Kari's shift ``sigma_i`` (neighborhood
    ``(e_i)``, identity local rule) restricted to the torus.

    Raises:
        ValueError: If ``shape`` and ``by`` disagree in length, or ``config``
            does not fill ``shape``.
    """
    sizes = tuple(shape)
    if len(by) != len(sizes) or len(config) != math.prod(sizes):
        msg = (
            f"cannot translate {len(config)} cells of shape {sizes!r} by "
            f"{tuple(by)!r}: the lengths disagree"
        )
        raise ValueError(msg)
    neighbors = _torus_neighbors((tuple(by),), sizes)
    return tuple(config[source] for (source,) in neighbors)


# --------------------------------------------------------------------------- #
# Finite phase spaces
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, slots=True)
class Orbit:
    """The eventually periodic shape of one trajectory on a ring or torus.

    ``transient`` is the number of steps before the trajectory first enters
    its cycle and ``period`` the length of that cycle (Cellular Automaton
    page, "Phase space"; the usage of ``research.strand_dynamics.Orbit``).
    """

    transient: int
    period: int


@dataclass(frozen=True)
class PhaseSpace:
    """The phase space of an automaton on one ring or torus.

    The directed graph on all ``k^N`` configurations with an edge
    ``c -> G(c)`` (Kari 2005, section 2.4). A configuration is numbered by
    reading its cells in row-major order as base-``k`` digits, first cell
    most significant; ``successor[i]`` is the number of the image of
    configuration ``i``. It "consists of disjoint cycles (periodic states),
    with directed trees (transient states) attached" (Macauley-Mortveit,
    section 2).
    """

    states: int
    shape: tuple[int, ...]
    successor: tuple[int, ...]

    @override
    def __repr__(self) -> str:
        """Return the sizes only; the successor table is ``k^N`` entries."""
        return (
            f"PhaseSpace(states={self.states}, shape={self.shape!r}, "
            f"configurations={len(self.successor)})"
        )

    @property
    def cells(self) -> int:
        """The number ``N`` of cells of the ring or torus."""
        return math.prod(self.shape)

    def configuration(self, index: int) -> Configuration:
        """Return the configuration numbered ``index``.

        Raises:
            ValueError: If ``index`` is not a configuration number.
        """
        if not 0 <= index < len(self.successor):
            last = len(self.successor) - 1
            msg = f"configuration numbers run over 0..{last}; got {index}"
            raise ValueError(msg)
        return _digits(index, self.states, self.cells)

    def index(self, config: Sequence[int]) -> int:
        """Return the number of ``config``, the inverse of :meth:`configuration`.

        Raises:
            ValueError: If ``config`` is not a configuration of this space.
        """
        if len(config) != self.cells or any(
            not 0 <= value < self.states for value in config
        ):
            msg = (
                f"a configuration assigns one of {self.states} states to each of "
                f"{self.cells} cells; got {tuple(config)!r}"
            )
            raise ValueError(msg)
        return _index(config, self.states)

    @functools.cached_property
    def in_degrees(self) -> tuple[int, ...]:
        """The number of predecessors of every configuration."""
        counts = Counter(self.successor)
        return tuple(counts[i] for i in range(len(self.successor)))

    @functools.cached_property
    def gardens_of_eden(self) -> tuple[int, ...]:
        """The configurations with no predecessor on this ring, ascending."""
        return tuple(i for i, degree in enumerate(self.in_degrees) if degree == 0)

    @property
    def is_bijective(self) -> bool:
        """Whether the automaton permutes the configurations of this ring."""
        return not self.gardens_of_eden

    @functools.cached_property
    def _peeled(self) -> tuple[tuple[int, ...], int]:
        """Peel the trees leaf-first: the periodic states and the tallest tree."""
        remaining = list(self.in_degrees)
        height = [0] * len(self.successor)
        queue = deque(self.gardens_of_eden)
        removed = [False] * len(self.successor)
        while queue:
            node = queue.popleft()
            removed[node] = True
            image = self.successor[node]
            height[image] = max(height[image], height[node] + 1)
            remaining[image] -= 1
            if remaining[image] == 0:
                queue.append(image)
        periodic = tuple(i for i, gone in enumerate(removed) if not gone)
        return periodic, max((height[i] for i in periodic), default=0)

    @property
    def periodic_states(self) -> tuple[int, ...]:
        """The temporally periodic configurations (those on cycles), ascending."""
        return self._peeled[0]

    @property
    def max_transient(self) -> int:
        """The longest transient: the height of the tallest tree on a cycle node."""
        return self._peeled[1]

    @functools.cached_property
    def cycle_lengths(self) -> tuple[int, ...]:
        """The length of every cycle, one entry per cycle, ascending."""
        seen: set[int] = set()
        lengths: list[int] = []
        for start in self.periodic_states:
            if start in seen:
                continue
            length, node = 0, start
            while node not in seen:
                seen.add(node)
                node = self.successor[node]
                length += 1
            lengths.append(length)
        return tuple(sorted(lengths))

    @property
    def fixed_points(self) -> tuple[int, ...]:
        """The configurations with ``G(c) = c``, ascending."""
        return tuple(i for i, image in enumerate(self.successor) if i == image)


# --------------------------------------------------------------------------- #
# The automaton
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CellularAutomaton:
    r"""A cellular automaton ``(d, S, N, f)`` on ``Z^d`` with ``S = {0..k-1}``.

    Definition (D1) of the Cellular Automaton page (Kari 1994, section 2.1):
    ``offsets`` is the neighborhood vector ``N = (x_1, ..., x_n)`` of
    distinct elements of ``Z^d`` and ``table`` the local rule
    ``f : S^n -> S``, with ``f(a_1, ..., a_n)`` stored at index
    ``sum a_i k^(n-i)`` (first neighbor most significant). The global
    function is ``G(c)(x) = f(c(x + x_1), ..., c(x + x_n))``. For the
    elementary neighborhood ``(-1, 0, 1)`` the index is ``4a + 2b + c``, so
    the table read as bits is the Wolfram number. The empty neighborhood is
    allowed: it carries the constant maps, whose minimal memory set is empty.

    The definition is one linear scan, so ``__post_init__`` checks it and
    direct construction is validated too; the ``from_*`` classmethods convert
    the page's other formulations. Equality is equality of the quadruple;
    use :meth:`same_global_map` to compare global functions.

    Raises:
        ValueError: Naming the violated clause of (D1).
    """

    dimension: int
    states: int
    offsets: tuple[Offset, ...]
    table: tuple[int, ...]

    def __post_init__(self) -> None:
        """Validate the four clauses of (D1)."""
        _require_dimension(self.dimension)
        if self.states < 1:
            msg = (
                f"(D1) state set violated: S must be finite and non-empty; "
                f"got {self.states} states"
            )
            raise ValueError(msg)
        for offset in self.offsets:
            if len(offset) != self.dimension:
                msg = (
                    f"(D1) neighborhood vector violated: {offset!r} is not an "
                    f"element of Z^{self.dimension}"
                )
                raise ValueError(msg)
        if len(set(self.offsets)) != len(self.offsets):
            msg = (
                f"(D1) neighborhood vector violated: its elements must be "
                f"distinct; got {self.offsets!r}"
            )
            raise ValueError(msg)
        expected = self.states ** len(self.offsets)
        if len(self.table) != expected:
            msg = (
                f"(D1) local rule violated: f : S^{len(self.offsets)} -> S needs "
                f"{expected} table entries; got {len(self.table)}"
            )
            raise ValueError(msg)
        if any(not 0 <= value < self.states for value in self.table):
            msg = (
                f"(D1) local rule violated: f must take values in "
                f"S = 0..{self.states - 1}; got {self.table!r}"
            )
            raise ValueError(msg)

    @override
    def __repr__(self) -> str:
        """Return a compact form: the Wolfram number when elementary."""
        if self._is_elementary:
            return f"CellularAutomaton(elementary rule {self.wolfram_number})"
        offsets = ", ".join(
            str(o[0]) if self.dimension == 1 else repr(o) for o in self.offsets
        )
        shown = self.table[:_REPR_TABLE_LIMIT]
        body = (
            "".join(map(str, shown))
            if self.states <= _DIGIT_STATES
            else ",".join(map(str, shown))
        )
        hidden = len(self.table) - len(shown)
        suffix = f"...(+{hidden})" if hidden else ""
        return (
            f"CellularAutomaton(d={self.dimension}, k={self.states}, "
            f"N=({offsets}), f={body}{suffix})"
        )

    # -------------------------------------------------------- constructors
    @classmethod
    def from_local_rule(
        cls,
        dimension: int,
        states: int,
        offsets: Iterable[Sequence[int]],
        rule: LocalRule,
    ) -> CellularAutomaton:
        """Build from (D1) with the local rule given as a callable.

        Kari 1994, section 2.1. ``rule`` is called as ``rule((a_1, ..., a_n))``
        on every tuple of ``S^n``, in the order of ``offsets``.

        Raises:
            ValueError: Naming the violated clause of (D1).
        """
        neighborhood = tuple(tuple(offset) for offset in offsets)
        if states < 1:
            msg = (
                f"(D1) state set violated: S must be finite and non-empty; "
                f"got {states} states"
            )
            raise ValueError(msg)
        table = tuple(
            rule(neighbors)
            for neighbors in itertools.product(range(states), repeat=len(neighborhood))
        )
        return cls(dimension, states, neighborhood, table)

    @classmethod
    def from_memory_set(
        cls,
        dimension: int,
        states: int,
        memory_set: Iterable[Sequence[int]],
        mu: Callable[[Mapping[Offset, int]], int],
    ) -> CellularAutomaton:
        r"""Build from (D2) over the group ``G = Z^d``.

        Ceccherini-Silberstein-Coornaert 2017, Definition 3.1:
        ``tau(x)(g) = mu((g^{-1} x)|_S)``. In additive notation
        ``(g^{-1} x)(s) = x(g + s)``, so the memory set is the set of
        entries of ``N`` and ``mu = f`` (page, "How the formulations
        connect"). ``mu`` receives the pattern ``{s: x(g + s)}``; the
        neighborhood vector is the memory set in sorted order.

        Raises:
            ValueError: If the memory set repeats an element or leaves
                ``Z^d``, naming (D2).
        """
        elements = [tuple(element) for element in memory_set]
        if len(set(elements)) != len(elements) or any(
            len(element) != dimension for element in elements
        ):
            msg = (
                f"(D2) memory set violated: S must be a finite set of elements "
                f"of Z^{dimension}; got {elements!r}"
            )
            raise ValueError(msg)
        offsets = tuple(sorted(elements))

        def rule(neighbors: tuple[int, ...]) -> int:
            return mu(dict(zip(offsets, neighbors, strict=True)))

        return cls.from_local_rule(dimension, states, offsets, rule)

    @classmethod
    def from_block_map(
        cls,
        states: int,
        memory: int,
        anticipation: int,
        rule: LocalRule,
    ) -> CellularAutomaton:
        """Build from (D4), a sliding block map with memory and anticipation.

        Beal-Berstel-Eilers-Perrin, section 2.2:
        ``y_i = f(x_{i-m} ... x_i ... x_{i+n})``, which is (D1) with
        ``d = 1`` and ``N = (-m, ..., n)``. Wolfram's ``(k, r)`` rules
        (1984, eq. 2.1) are the case ``m = n = r``.

        Raises:
            ValueError: If ``memory`` or ``anticipation`` is negative, naming
                (D4), or a clause of (D1) fails.
        """
        if memory < 0 or anticipation < 0:
            msg = (
                f"(D4) violated: memory and anticipation are non-negative; got "
                f"m = {memory}, n = {anticipation}"
            )
            raise ValueError(msg)
        offsets = [(x,) for x in range(-memory, anticipation + 1)]
        return cls.from_local_rule(1, states, offsets, rule)

    @classmethod
    def from_tessellation(
        cls,
        dimension: int,
        states: int,
        quiescent: int,
        rule: LocalRule,
    ) -> CellularAutomaton:
        """Build from (D5), a tessellation structure with a quiescent state.

        Moore 1962, p. 21: the rule reads the ``3^N`` cells of the Moore
        neighborhood (in sorted order) and a cell whose neighbors are all
        quiescent stays quiescent, ``f(q0, ..., q0) = q0``.

        Raises:
            ValueError: If ``quiescent`` is not a quiescent state, naming
                (D5), or a clause of (D1) fails.
        """
        automaton = cls.from_local_rule(
            dimension, states, moore_neighborhood(dimension), rule
        )
        automaton._require_quiescent(quiescent)
        return automaton

    @classmethod
    def from_wolfram_number(cls, number: int) -> CellularAutomaton:
        """Build the elementary automaton with the given rule number.

        Wolfram 1983, Fig. 1 (Kari 2005, section 2.5): ``d = 1``,
        ``S = {0, 1}``, ``N = (-1, 0, 1)`` and
        ``R = sum f(a, b, c) 2^(4a + 2b + c)``.

        Raises:
            ValueError: If ``number`` is outside ``0..255``.
        """
        if not 0 <= number < _WOLFRAM_NUMBERS:
            msg = f"elementary rule numbers run over 0..255; got {number}"
            raise ValueError(msg)
        table = tuple(number >> index & 1 for index in range(_BINARY**3))
        return cls(1, _BINARY, _ELEMENTARY_OFFSETS, table)

    @classmethod
    def from_linear(
        cls,
        dimension: int,
        modulus: int,
        coefficients: Mapping[Sequence[int], int],
    ) -> CellularAutomaton:
        """Build the linear automaton ``f(a_1..a_n) = sum c_i a_i`` over ``Z_m``.

        Kari 2005, section 9. ``coefficients`` maps each offset to its
        coefficient (reduced modulo ``m``); as a Laurent polynomial it is
        ``p(Z) = sum c_i Z^(-x_i)``. Offsets are taken in sorted order.

        Raises:
            ValueError: Naming the violated clause of (D1).
        """
        ordered = sorted(
            (tuple(offset), c % modulus) for offset, c in coefficients.items()
        )
        weights = [c for _, c in ordered]

        def rule(neighbors: tuple[int, ...]) -> int:
            return sum(c * a for c, a in zip(weights, neighbors, strict=True)) % modulus

        return cls.from_local_rule(
            dimension, modulus, [offset for offset, _ in ordered], rule
        )

    # -------------------------------------------------------- the local rule
    @property
    def neighborhood_size(self) -> int:
        """The number ``n`` of neighbors the local rule reads."""
        return len(self.offsets)

    @property
    def _origin(self) -> Offset:
        return (0,) * self.dimension

    @property
    def _is_elementary(self) -> bool:
        return self.states == _BINARY and self.offsets == _ELEMENTARY_OFFSETS

    def local_rule(self, *neighbors: int) -> int:
        """Return ``f(a_1, ..., a_n)``.

        Raises:
            ValueError: If the arguments are not ``n`` states.
        """
        if len(neighbors) != len(self.offsets) or any(
            not 0 <= value < self.states for value in neighbors
        ):
            msg = (
                f"the local rule takes {len(self.offsets)} states in "
                f"0..{self.states - 1}; got {neighbors!r}"
            )
            raise ValueError(msg)
        return self.table[_index(neighbors, self.states)]

    @property
    def wolfram_number(self) -> int:
        """The Wolfram number ``sum f(a, b, c) 2^(4a + 2b + c)`` (Wolfram 1983, Fig. 1).

        Raises:
            ValueError: If the automaton is not elementary (``d = 1``,
                ``S = {0, 1}``, ``N = (-1, 0, 1)``); the page defines the
                number for those only.
        """
        if not self._is_elementary:
            msg = (
                f"the Wolfram number is defined for elementary automata "
                f"(d = 1, S = {{0, 1}}, N = (-1, 0, 1)); got {self!r}"
            )
            raise ValueError(msg)
        return sum(bit << index for index, bit in enumerate(self.table))

    @property
    def quiescent_states(self) -> tuple[int, ...]:
        """The states ``q`` with ``f(q, ..., q) = q`` (Kari 2005, section 2.3)."""
        n = len(self.offsets)
        return tuple(
            q
            for q in range(self.states)
            if self.table[_index([q] * n, self.states)] == q
        )

    def _require_quiescent(self, quiescent: int) -> None:
        if quiescent not in self.quiescent_states:
            msg = (
                f"(D5) quiescent state violated: f(q0, ..., q0) = q0 fails for "
                f"q0 = {quiescent}; the quiescent states are {self.quiescent_states!r}"
            )
            raise ValueError(msg)

    @property
    def is_balanced(self) -> bool:
        """Whether every state has exactly ``|S|^(n-1)`` preimage tuples.

        Page, "Balanced". Necessary for surjectivity (Kari's notes,
        Corollary 14) but not sufficient: rule 232 is balanced.
        """
        counts = Counter(self.table)
        return all(
            counts[state] * self.states == len(self.table)
            for state in range(self.states)
        )

    def _depends_on(self, position: int) -> bool:
        """Return whether ``f`` depends on its argument at ``position``."""
        weight = self.states ** (len(self.offsets) - 1 - position)
        return any(
            self.table[index] != self.table[index - digit * weight]
            for index in range(len(self.table))
            if (digit := index // weight % self.states)
        )

    @property
    def is_totalistic(self) -> bool:
        """Whether ``f`` depends only on the neighborhood sum (Wolfram 1984)."""
        by_sum: dict[int, int] = {}
        n = len(self.offsets)
        return all(
            by_sum.setdefault(sum(_digits(index, self.states, n)), value) == value
            for index, value in enumerate(self.table)
        )

    @property
    def is_additive(self) -> bool:
        """Whether ``f`` is linear modulo ``k``, ``f = sum c_i a_i`` (Wolfram 1983).

        Such rules obey a superposition principle; over ``Z_m`` they are the
        linear automata of Kari 2005, section 9.
        """
        n, k = len(self.offsets), self.states
        if k == 1:
            return True
        weights = [self.table[k ** (n - 1 - position)] for position in range(n)]
        return all(
            value
            == sum(c * a for c, a in zip(weights, _digits(index, k, n), strict=True))
            % k
            for index, value in enumerate(self.table)
        )

    @property
    def is_peripheral(self) -> bool:
        """Whether the new value ignores the cell's own old value (Wolfram 1983)."""
        if self._origin not in self.offsets:
            return True
        return not self._depends_on(self.offsets.index(self._origin))

    @property
    def linear_coefficients(self) -> dict[Offset, int]:
        """The coefficients ``c_i`` by offset: the polynomial ``sum c_i Z^(-x_i)``.

        Kari 2005, section 9 (Laurent-polynomial calculus): products of
        these polynomials correspond to composition.

        Raises:
            ValueError: If the rule is not additive.
        """
        if not self.is_additive:
            msg = f"only an additive rule has linear coefficients; got {self!r}"
            raise ValueError(msg)
        n, k = len(self.offsets), self.states
        if k == 1:
            return dict.fromkeys(self.offsets, 0)
        return {
            offset: self.table[k ** (n - 1 - position)]
            for position, offset in enumerate(self.offsets)
        }

    # -------------------------------------------------------- neighborhoods
    def _retabulate(self, offsets: tuple[Offset, ...]) -> CellularAutomaton:
        """Re-express ``f`` on ``offsets``; neighbors absent there read state 0."""
        position = {offset: i for i, offset in enumerate(offsets)}
        sources = [position.get(offset) for offset in self.offsets]
        table = tuple(
            self.table[
                _index(
                    (0 if source is None else neighbors[source] for source in sources),
                    self.states,
                )
            ]
            for neighbors in itertools.product(range(self.states), repeat=len(offsets))
        )
        return CellularAutomaton(self.dimension, self.states, offsets, table)

    def with_neighborhood(self, offsets: Iterable[Sequence[int]]) -> CellularAutomaton:
        """Return the same global map on a larger neighborhood vector.

        Ceccherini-Silberstein-Coornaert 2017, section 3.1: "Any finite
        superset of a memory set is again one."

        Raises:
            ValueError: If ``offsets`` omits a current neighbor, or a clause
                of (D1) fails.
        """
        neighborhood = tuple(tuple(offset) for offset in offsets)
        missing = [offset for offset in self.offsets if offset not in neighborhood]
        if missing:
            msg = f"the new neighborhood must contain the old one; it omits {missing!r}"
            raise ValueError(msg)
        return self._retabulate(neighborhood)

    def minimal(self) -> CellularAutomaton:
        """Return the same global map on its minimal memory set, sorted.

        Ceccherini-Silberstein-Coornaert 2017, section 3.1: there is a
        unique memory set of minimal cardinality. It consists of the
        neighbors ``f`` depends on, so this form is canonical.
        """
        needed = sorted(
            offset for i, offset in enumerate(self.offsets) if self._depends_on(i)
        )
        return self._retabulate(tuple(needed))

    def same_global_map(self, other: CellularAutomaton) -> bool:
        """Return whether both automata define the same global function.

        Kari 2005, section 2.4: equality of two cellular automata is
        decidable — here by comparing minimal forms.
        """
        return (
            self.dimension == other.dimension
            and self.states == other.states
            and self.minimal() == other.minimal()
        )

    # -------------------------------------------------------- periodic dynamics
    def _resolve_shape(
        self, cells: int, shape: Sequence[int] | None
    ) -> tuple[int, ...]:
        if shape is None:
            if self.dimension != 1:
                msg = (
                    f"a torus configuration of a {self.dimension}-dimensional "
                    f"automaton needs an explicit shape"
                )
                raise ValueError(msg)
            return (cells,)
        return tuple(shape)

    def _check_shape(self, shape: tuple[int, ...]) -> None:
        if len(shape) != self.dimension or any(size < 1 for size in shape):
            msg = (
                f"a torus for a {self.dimension}-dimensional automaton has "
                f"{self.dimension} positive periods; got {shape!r}"
            )
            raise ValueError(msg)

    def _check_config(self, config: Configuration, shape: tuple[int, ...]) -> None:
        self._check_shape(shape)
        if len(config) != math.prod(shape) or any(
            not 0 <= value < self.states for value in config
        ):
            msg = (
                f"a configuration assigns a state in 0..{self.states - 1} to each of "
                f"the {math.prod(shape)} cells of shape {shape!r}; got {config!r}"
            )
            raise ValueError(msg)

    def _step(self, config: Configuration, shape: tuple[int, ...]) -> Configuration:
        table, k = self.table, self.states
        result: list[int] = []
        for neighbors in _torus_neighbors(self.offsets, shape):
            index = 0
            for source in neighbors:
                index = index * k + config[source]
            result.append(table[index])
        return tuple(result)

    def step(
        self, config: Sequence[int], shape: Sequence[int] | None = None
    ) -> Configuration:
        """Return ``G(c)`` for a periodic configuration on a ring or torus.

        Page, "Periodic configuration": a CA on a ring ``Z/N`` (or a torus)
        is ``G_P`` on a fixed period. Offsets wrap around, so rings shorter
        than the neighborhood are allowed.

        Args:
            config: The cells in row-major order.
            shape: The periods of the torus; ``(len(config),)`` by default,
                which requires ``d = 1``.

        Raises:
            ValueError: If ``config`` is not a configuration of that torus.
        """
        cells = tuple(config)
        periods = self._resolve_shape(len(cells), shape)
        self._check_config(cells, periods)
        return self._step(cells, periods)

    def evolve(
        self, config: Sequence[int], steps: int, shape: Sequence[int] | None = None
    ) -> tuple[Configuration, ...]:
        """Return the ``steps + 1`` configurations ``c, G(c), ..., G^steps(c)``.

        Raises:
            ValueError: If ``steps`` is negative or ``config`` is not a
                configuration of the torus.
        """
        if steps < 0:
            msg = f"steps must be non-negative; got {steps}"
            raise ValueError(msg)
        cells = tuple(config)
        periods = self._resolve_shape(len(cells), shape)
        self._check_config(cells, periods)
        rows = [cells]
        for _ in range(steps):
            rows.append(self._step(rows[-1], periods))
        return tuple(rows)

    def orbit(self, config: Sequence[int], shape: Sequence[int] | None = None) -> Orbit:
        """Follow the trajectory of ``config`` until a configuration repeats.

        The state space is finite, so this always terminates, with
        ``transient + period <= k^N`` (pigeonhole; Wolfram 1983, section 4).
        Every visited configuration is kept in memory.

        Raises:
            ValueError: If ``config`` is not a configuration of the torus.
        """
        current = tuple(config)
        periods = self._resolve_shape(len(current), shape)
        self._check_config(current, periods)
        first_seen: dict[Configuration, int] = {}
        time = 0
        while current not in first_seen:
            first_seen[current] = time
            current = self._step(current, periods)
            time += 1
        transient = first_seen[current]
        return Orbit(transient=transient, period=time - transient)

    def phase_space(self, shape: Sequence[int] | int) -> PhaseSpace:
        """Return the phase space on the torus of the given shape.

        Costs ``k^N`` rule sweeps of ``N`` cells and as much memory; an
        ``int`` is read as a ring size.

        Raises:
            ValueError: If ``shape`` is not a torus for this dimension.
        """
        periods = (shape,) if isinstance(shape, int) else tuple(shape)
        self._check_shape(periods)
        cells = math.prod(periods)
        successor = tuple(
            _index(self._step(config, periods), self.states)
            for config in itertools.product(range(self.states), repeat=cells)
        )
        return PhaseSpace(self.states, periods, successor)

    def conserves_number_on_ring(self, shape: Sequence[int] | int) -> bool:
        """Return whether the sum of the cell values is preserved on this torus.

        Number conservation (Boccara-Fuks 2002) in the page's computational
        sense, one ring at a time; exhaustive over its ``k^N`` configurations.

        Raises:
            ValueError: If ``shape`` is not a torus for this dimension.
        """
        periods = (shape,) if isinstance(shape, int) else tuple(shape)
        self._check_shape(periods)
        return all(
            sum(self._step(config, periods)) == sum(config)
            for config in itertools.product(
                range(self.states), repeat=math.prod(periods)
            )
        )

    # -------------------------------------------------------- finite dynamics
    def _step_finite(
        self, cells: Mapping[Offset, int], quiescent: int
    ) -> dict[Offset, int]:
        n, k = len(self.offsets), self.states
        base = _index([quiescent] * n, k)
        weights = [k ** (n - 1 - position) for position in range(n)]
        shifts: dict[Offset, int] = {}
        for cell, state in cells.items():
            if state == quiescent:
                continue
            for offset, weight in zip(self.offsets, weights, strict=True):
                target = tuple(c - o for c, o in zip(cell, offset, strict=True))
                shifts[target] = shifts.get(target, 0) + (state - quiescent) * weight
        image = {cell: self.table[base + shift] for cell, shift in shifts.items()}
        return {cell: state for cell, state in image.items() if state != quiescent}

    def _check_finite(self, cells: Mapping[Offset, int], quiescent: int) -> None:
        self._require_quiescent(quiescent)
        for cell, state in cells.items():
            if len(cell) != self.dimension or not 0 <= state < self.states:
                msg = (
                    f"a finite configuration maps cells of Z^{self.dimension} to "
                    f"states in 0..{self.states - 1}; got {cell!r}: {state!r}"
                )
                raise ValueError(msg)

    def step_finite(
        self, cells: Mapping[Offset, int], quiescent: int = 0
    ) -> dict[Offset, int]:
        """Return ``G_F(c)`` for a finite configuration on the unbounded lattice.

        Kari 2005, section 2.3: a configuration is finite if its support
        ``{x : c(x) != q}`` is finite, and ``G`` preserves finiteness.

        Args:
            cells: The non-quiescent cells (quiescent entries are ignored).
            quiescent: The quiescent state ``q``.

        Returns:
            The non-quiescent cells of the image.

        Raises:
            ValueError: If ``quiescent`` is not a quiescent state, naming
                (D5), or ``cells`` is not a configuration.
        """
        self._check_finite(cells, quiescent)
        return self._step_finite(cells, quiescent)

    def evolve_finite(
        self, cells: Mapping[Offset, int], steps: int, quiescent: int = 0
    ) -> tuple[dict[Offset, int], ...]:
        """Return the ``steps + 1`` finite configurations ``c, ..., G_F^steps(c)``.

        Raises:
            ValueError: If ``steps`` is negative, ``quiescent`` is not a
                quiescent state, or ``cells`` is not a configuration.
        """
        if steps < 0:
            msg = f"steps must be non-negative; got {steps}"
            raise ValueError(msg)
        self._check_finite(cells, quiescent)
        rows = [{cell: s for cell, s in cells.items() if s != quiescent}]
        for _ in range(steps):
            rows.append(self._step_finite(rows[-1], quiescent))
        return tuple(rows)

    # -------------------------------------------------------- the line (d = 1)
    def _require_line(self, what: str) -> None:
        if self.dimension != 1:
            msg = (
                f"{what} is implemented for one-dimensional automata only "
                f"(surjectivity and reversibility are undecidable from "
                f"dimension two on, Kari 1994); got dimension {self.dimension}"
            )
            raise ValueError(msg)

    def to_contiguous(self) -> CellularAutomaton:
        """Return the same global map on the contiguous neighborhood ``(min..max)``.

        Amoroso-Patt 1972, section VI: any neighborhood can be replaced by
        a contiguous one without loss of generality. The empty neighborhood
        becomes ``(0)``.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("the contiguous form")
        positions = [offset[0] for offset in self.offsets] or [0]
        return self._retabulate(
            tuple((x,) for x in range(min(positions), max(positions) + 1))
        )

    @functools.cached_property
    def _line(self) -> CellularAutomaton:
        return self.to_contiguous()

    @property
    def _width(self) -> int:
        """The width ``m`` of the contiguous neighborhood."""
        return len(self._line.offsets)

    @property
    def _nodes(self) -> int:
        """The number ``k^(m-1)`` of de Bruijn vertices."""
        return int(self.states ** (self._width - 1))

    def de_bruijn_graph(self) -> tuple[Edge, ...]:
        """Return the labeled de Bruijn graph of width ``m``, one edge per word.

        Kari's lecture notes, section 2.6: vertices ``S^(m-1)``, the edge
        ``s_1..s_m`` going from ``s_1..s_(m-1)`` to ``s_2..s_m`` and labeled
        ``f(s_1, ..., s_m)``, for the contiguous form of the neighborhood.
        The graph forgets where the neighborhood sits.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("the de Bruijn graph")
        return tuple(
            (word[:-1], word[1:], self._line.table[index])
            for index, word in enumerate(
                itertools.product(range(self.states), repeat=self._width)
            )
        )

    @functools.cached_property
    def _label_successors(self) -> tuple[tuple[int, ...], ...]:
        """``[label][node]``: bitmask of nodes reached by an edge with that label."""
        k, nodes = self.states, self._nodes
        masks = [[0] * nodes for _ in range(k)]
        for index, label in enumerate(self._line.table):
            masks[label][index // k] |= 1 << (index % nodes)
        return tuple(tuple(row) for row in masks)

    def _advance_subset(self, subset: int, label: int) -> int:
        successors = self._label_successors[label]
        image = 0
        for node in bits(subset):
            image |= successors[node]
        return image

    def is_surjective(self) -> bool:
        """Decide surjectivity of the global map on ``S^Z``.

        Subset construction on the de Bruijn graph (Sutner 1991, as
        described in Kari's notes, section 2.6): a word is an orphan iff it
        drives the full vertex set to the empty set, and ``G`` is
        surjective iff it has no orphan (Kari's notes, Proposition 11).
        Worst case exponential in the ``k^(m-1)`` vertices.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("surjectivity")
        full = (1 << self._nodes) - 1
        seen = {full}
        stack = [full]
        while stack:
            subset = stack.pop()
            for label in range(self.states):
                image = self._advance_subset(subset, label)
                if image == 0:
                    return False
                if image not in seen:
                    seen.add(image)
                    stack.append(image)
        return True

    @functools.cached_property
    def _pair_graph(self) -> tuple[int, tuple[tuple[int, ...], ...]]:
        """The vertex count ``V`` and successor lists; ``(u1, u2)`` is ``u1 * V + u2``.

        Built on a contiguous form of width at least two: with a single
        de Bruijn vertex two distinct configurations would trace the same
        vertex path, and the off-diagonal criteria could never fire.
        """
        line = self._line
        if len(line.offsets) == 1:
            ((left,),) = line.offsets
            line = line.with_neighborhood([(left,), (left + 1,)])
        k, table = self.states, line.table
        nodes = len(table) // k
        by_label: list[list[list[int]]] = [[[] for _ in range(k)] for _ in range(nodes)]
        for index, label in enumerate(table):
            by_label[index // k][label].append(index % nodes)
        graph = tuple(
            tuple(
                v1 * nodes + v2
                for label in range(k)
                for v1 in by_label[u1][label]
                for v2 in by_label[u2][label]
            )
            for u1 in range(nodes)
            for u2 in range(nodes)
        )
        return nodes, graph

    def _reachable(
        self, graph: Sequence[Sequence[int]], sources: Iterable[int]
    ) -> set[int]:
        seen = set(sources)
        stack = list(seen)
        while stack:
            for image in graph[stack.pop()]:
                if image not in seen:
                    seen.add(image)
                    stack.append(image)
        return seen

    def is_pre_injective(self) -> bool:
        """Decide pre-injectivity: no two distinct almost equal configurations collide.

        A diamond is a path of the pair graph that leaves the diagonal and
        returns to it (Kari's notes, section 2.6), so ``G`` is pre-injective
        iff no off-diagonal vertex is both reachable from and co-reachable
        to the diagonal. By the Garden of Eden theorem this agrees with
        :meth:`is_surjective`, which is computed independently.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("pre-injectivity")
        nodes, graph = self._pair_graph
        diagonal = [u * nodes + u for u in range(nodes)]
        reverse: list[list[int]] = [[] for _ in graph]
        for source, images in enumerate(graph):
            for image in images:
                reverse[image].append(source)
        between = self._reachable(graph, diagonal) & self._reachable(reverse, diagonal)
        return all(node // nodes == node % nodes for node in between)

    def is_injective(self) -> bool:
        """Decide injectivity of the global map on ``S^Z``.

        Pair graph (Kari's notes, Proposition 26): ``G`` is non-injective
        iff the pair graph has a cycle through an off-diagonal vertex. Two
        distinct configurations with one image are exactly a bi-infinite
        path through an off-diagonal vertex, so the vertices that survive
        trimming of dead ends in both directions are inspected.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("injectivity")
        nodes, graph = self._pair_graph
        reverse: list[list[int]] = [[] for _ in graph]
        for source, images in enumerate(graph):
            for image in images:
                reverse[image].append(source)
        out_degree = [len(images) for images in graph]
        in_degree = [len(sources) for sources in reverse]
        alive = [True] * len(graph)
        queue = deque(
            node
            for node in range(len(graph))
            if not (out_degree[node] and in_degree[node])
        )
        while queue:
            node = queue.popleft()
            if not alive[node]:
                continue
            alive[node] = False
            for image in graph[node]:
                in_degree[image] -= 1
                if alive[image] and in_degree[image] == 0:
                    queue.append(image)
            for source in reverse[node]:
                out_degree[source] -= 1
                if alive[source] and out_degree[source] == 0:
                    queue.append(source)
        return all(
            node // nodes == node % nodes for node in range(len(graph)) if alive[node]
        )

    def is_reversible(self) -> bool:
        """Decide reversibility, which on ``Z^d`` is injectivity.

        Hedlund 1969 and Richardson 1972, independently (per Kari 2005):
        an injective automaton on ``Z^d`` is bijective and its inverse is a
        cellular automaton.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        return self.is_injective()

    def _check_word(self, word: Sequence[int]) -> Word:
        result = tuple(word)
        if any(not 0 <= value < self.states for value in result):
            msg = f"a word lists states in 0..{self.states - 1}; got {result!r}"
            raise ValueError(msg)
        return result

    def apply_to_word(self, word: Sequence[int]) -> Word:
        """Return the image of a word, shorter by ``m - 1`` symbols.

        ``m`` is the width of the contiguous neighborhood; the image of
        ``s_1..s_L`` is ``f(s_1..s_m) f(s_2..s_(m+1)) ...``.

        Raises:
            ValueError: If the automaton is not one-dimensional or the word
                is shorter than ``m - 1``.
        """
        self._require_line("the image of a word")
        symbols = self._check_word(word)
        width = self._width
        if len(symbols) < width - 1:
            msg = (
                f"a word must cover at least m - 1 = {width - 1} cells; got {symbols!r}"
            )
            raise ValueError(msg)
        table = self._line.table
        return tuple(
            table[_index(symbols[start : start + width], self.states)]
            for start in range(len(symbols) - width + 1)
        )

    def preimage_count(self, word: Sequence[int]) -> int:
        """Return the number of words of length ``len(word) + m - 1`` mapping to it.

        The pattern count of the balance theorem (Kari's notes,
        Proposition 13): for a surjective rule it is ``|S|^(m-1)`` for
        every word.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("preimage counting")
        k, nodes, table = self.states, self._nodes, self._line.table
        counts = [1] * nodes
        for label in self._check_word(word):
            following = [0] * nodes
            for index, value in enumerate(table):
                if value == label:
                    following[index % nodes] += counts[index // k]
            counts = following
        return sum(counts)

    def is_orphan(self, word: Sequence[int]) -> bool:
        """Return whether ``word`` occurs in no image (a Garden of Eden pattern).

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        return self.preimage_count(word) == 0

    def shortest_orphans(self) -> tuple[Word, ...]:
        """Return every orphan of minimal length, sorted; empty when surjective.

        Breadth-first over the subset construction of
        :meth:`is_surjective`, keeping every word of the current length, so
        the cost is ``k^L`` for shortest orphans of length ``L``.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        if not self.is_surjective():
            frontier: dict[Word, int] = {(): (1 << self._nodes) - 1}
            while True:
                following = {
                    (*word, label): self._advance_subset(subset, label)
                    for word, subset in frontier.items()
                    for label in range(self.states)
                }
                orphans = sorted(
                    word for word, subset in following.items() if not subset
                )
                if orphans:
                    return tuple(orphans)
                frontier = following
        return ()

    def diamond(self, length: int) -> tuple[Word, Word] | None:
        """Return a diamond whose differing segment has the given length, if any.

        A pair of words ``u w v`` and ``u w' v`` with ``|u| = |v| = m - 1``,
        ``|w| = |w'| = length`` and ``w != w'`` that have the same image —
        two mutually erasable patterns (Ceccherini-Silberstein-Coornaert
        2017, Proposition 3.18). Exhaustive over ``k^(2m - 2 + length)``
        words; the first collision in lexicographic order is returned.

        Raises:
            ValueError: If the automaton is not one-dimensional or
                ``length`` is not positive.
        """
        self._require_line("the diamond search")
        if length < 1:
            msg = (
                f"the differing segment of a diamond is non-empty; got length {length}"
            )
            raise ValueError(msg)
        border = self._width - 1
        seen: dict[tuple[Word, Word, Word], Word] = {}
        for word in itertools.product(range(self.states), repeat=2 * border + length):
            key = (word[:border], word[len(word) - border :], self.apply_to_word(word))
            if key in seen:
                return seen[key], word
            seen[key] = word
        return None

    @property
    def is_legal(self) -> bool:
        """Whether the rule is legal in Wolfram's sense (1983).

        The null configuration is fixed, ``f(0, ..., 0) = 0``, and the rule
        is reflection-symmetric.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("legality")
        return self.table[0] == 0 and self.reflect().same_global_map(self)

    # -------------------------------------------------------- transformations
    def _require_compatible(self, other: CellularAutomaton) -> None:
        if (self.dimension, self.states) != (other.dimension, other.states):
            msg = (
                f"composition needs one dimension and one state set; got "
                f"{self!r} and {other!r}"
            )
            raise ValueError(msg)

    def compose(self, other: CellularAutomaton) -> CellularAutomaton:
        """Return ``self o other``, the automaton applying ``other`` first.

        Kari 2005, section 2.4: the composition is a cellular automaton
        with neighborhood ``{x + y : x in N_1, y in N_2}`` (sorted here).

        Raises:
            ValueError: If the dimensions or state sets differ.
        """
        self._require_compatible(other)
        sums = sorted(
            {
                tuple(a + b for a, b in zip(x, y, strict=True))
                for x in self.offsets
                for y in other.offsets
            }
        )
        if not self.offsets:
            return dataclasses.replace(self)
        position = {offset: i for i, offset in enumerate(sums)}
        reads = [
            [
                position[tuple(a + b for a, b in zip(x, y, strict=True))]
                for y in other.offsets
            ]
            for x in self.offsets
        ]
        k = self.states
        table = tuple(
            self.table[
                _index(
                    (
                        other.table[_index((cells[i] for i in read), k)]
                        for read in reads
                    ),
                    k,
                )
            ]
            for cells in itertools.product(range(k), repeat=len(sums))
        )
        return CellularAutomaton(self.dimension, k, tuple(sums), table)

    def power(self, exponent: int) -> CellularAutomaton:
        """Return ``G^exponent``; the identity automaton for exponent 0.

        Raises:
            ValueError: If ``exponent`` is negative.
        """
        if exponent < 0:
            msg = f"a power of a cellular automaton is non-negative; got {exponent}"
            raise ValueError(msg)
        result = CellularAutomaton(
            self.dimension, self.states, (self._origin,), tuple(range(self.states))
        )
        for _ in range(exponent):
            result = self.compose(result)
        return result

    def reflect(self) -> CellularAutomaton:
        """Return the left-right mirror image, ``f(a, b, c) -> f(c, b, a)``.

        Kari 2005 (symmetry conjugation): the conjugate of ``G`` by the
        reflection ``x -> -x``; the neighborhood vector is negated and
        reversed, so ``(-1, 0, 1)`` is kept.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        self._require_line("reflection")
        negated = tuple((-offset[0],) for offset in self.offsets)
        mirrored = CellularAutomaton(self.dimension, self.states, negated, self.table)
        return mirrored._retabulate(negated[::-1])

    def complement(self) -> CellularAutomaton:
        """Return the state complement, ``f -> 1 - f(1 - a, 1 - b, 1 - c)``.

        Kari 2005 (symmetry conjugation), defined for two states.

        Raises:
            ValueError: If the automaton does not have exactly two states.
        """
        if self.states != _BINARY:
            msg = f"state complementation is defined for two states; got {self.states}"
            raise ValueError(msg)
        return dataclasses.replace(
            self, table=tuple(1 - value for value in reversed(self.table))
        )

    def symmetry_class(self) -> frozenset[CellularAutomaton]:
        """Return the orbit under reflection and complementation.

        The two involutions generate a four-element group; on the
        elementary rules it has 88 orbits (Kari 2005).

        Raises:
            ValueError: If the automaton is not one-dimensional with two
                states.
        """
        mirrored = self.reflect()
        return frozenset({self, mirrored, self.complement(), mirrored.complement()})

    def second_order(self) -> CellularAutomaton:
        r"""Return the second-order (Fredkin) automaton on the alphabet ``A x A``.

        Toffoli-Margolus 1990, section 5.4: from any automaton ``tau`` over
        states ``{0..r-1}``, ``q^{t+1} = tau q^t - q^{t-1} (mod r)`` is
        invertible, since ``q^{t-1} = tau q^t - q^{t+1}``. A cell of the
        result holds the pair ``(q^t, q^{t-1})`` encoded as
        ``r * q^t + q^{t-1}``; the origin joins the neighborhood if absent.
        """
        k = self.states
        offsets = self.offsets
        if self._origin not in offsets:
            offsets = (*offsets, self._origin)
        center = offsets.index(self._origin)
        n = len(self.offsets)

        def rule(pairs: tuple[int, ...]) -> int:
            current = [pair // k for pair in pairs]
            updated = (self.table[_index(current[:n], k)] - pairs[center] % k) % k
            return updated * k + current[center]

        return CellularAutomaton.from_local_rule(self.dimension, k * k, offsets, rule)

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame of the quadruple ``(d, S, N, f)``.

        Columns: ``kind``, ``index``, ``axis``, ``value`` (the last three
        nullable integers). One ``"dimension"`` row (``value = d``); one
        ``"states"`` row (``value = k``); one ``"offset"`` row per neighbor
        and axis (``index`` the 0-based position in ``N``, ``value`` the
        coordinate); one ``"rule"`` row per table entry (``index`` the table
        index, first neighbor most significant, ``value`` the output). The
        frame is never empty and survives ``experiments.io.write_result``
        (records-oriented JSON).

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        records: list[tuple[str, int | None, int | None, int]] = [
            (_DIMENSION, None, None, self.dimension),
            (_STATES, None, None, self.states),
        ]
        records.extend(
            (_OFFSET, position, axis, coordinate)
            for position, offset in enumerate(self.offsets)
            for axis, coordinate in enumerate(offset)
        )
        records.extend(
            (_RULE, index, None, value) for index, value in enumerate(self.table)
        )
        frame = pd.DataFrame.from_records(records, columns=list(_COLUMNS))
        for column in _COLUMNS[1:]:
            frame[column] = frame[column].astype("Int64")
        return frame

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> CellularAutomaton:
        """Rebuild an automaton from :meth:`to_dataframe` output.

        Rows may come in any order.

        Raises:
            ValueError: If a column is missing, the scalar rows are not
                unique, offsets or table indices have gaps, or (D1) fails.
        """
        missing = [column for column in _COLUMNS if column not in df.columns]
        if missing:
            msg = f"cellular-automaton frame is missing columns {missing!r}"
            raise ValueError(msg)
        scalars: dict[str, list[int]] = {_DIMENSION: [], _STATES: []}
        coordinates: dict[tuple[int, int], int] = {}
        outputs: dict[int, int] = {}
        for kind, index, axis, value in zip(
            df["kind"], df["index"], df["axis"], df["value"], strict=True
        ):
            if kind in scalars:
                scalars[kind].append(int(value))
            elif kind == _OFFSET:
                coordinates[int(index), int(axis)] = int(value)
            elif kind == _RULE:
                outputs[int(index)] = int(value)
            else:
                msg = f"cellular-automaton frame has an unknown row kind {kind!r}"
                raise ValueError(msg)
        if any(len(values) != 1 for values in scalars.values()):
            msg = "cellular-automaton frame needs one dimension row and one states row"
            raise ValueError(msg)
        dimension, states = scalars[_DIMENSION][0], scalars[_STATES][0]
        size = len({position for position, _ in coordinates})
        try:
            offsets = tuple(
                tuple(coordinates[position, axis] for axis in range(dimension))
                for position in range(size)
            )
            table = tuple(outputs[index] for index in range(len(outputs)))
        except KeyError as error:
            msg = f"cellular-automaton frame has a gap at {error.args[0]!r}"
            raise ValueError(msg) from error
        if len(coordinates) != size * dimension:
            msg = f"cellular-automaton frame has offsets outside Z^{dimension}"
            raise ValueError(msg)
        return CellularAutomaton(dimension, states, offsets, table)

    # -------------------------------------------------------- visualization
    def plot_space_time(
        self, config: Sequence[int], steps: int, ax: Axes | None = None
    ) -> Axes:
        """Draw the space-time diagram of a ring onto ``ax``, time downward.

        The picture of Wolfram 1983, section 2: row ``t`` is ``G^t(c)``.

        Args:
            config: The initial ring configuration.
            steps: The number of rule applications to draw.
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.

        Raises:
            ValueError: If the automaton is not one-dimensional, ``steps``
                is negative, or ``config`` is not a ring configuration.
        """
        self._require_line("the space-time diagram")
        rows = self.evolve(config, steps)
        ax = ensure_axes(ax)
        ax.imshow(
            [list(row) for row in rows],
            cmap="binary",
            interpolation="nearest",
            vmin=0,
            vmax=max(self.states - 1, 1),
            aspect="equal",
        )
        ax.set_xlabel("cell")
        ax.set_ylabel("time")
        ax.set_title("Space-time diagram")
        return ax

    def plot_configuration(
        self, config: Sequence[int], shape: Sequence[int], ax: Axes | None = None
    ) -> Axes:
        """Draw a two-dimensional torus configuration onto ``ax``.

        Args:
            config: The cells in row-major order.
            shape: The two periods ``(rows, columns)``.
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.

        Raises:
            ValueError: If the automaton is not two-dimensional or
                ``config`` is not a configuration of the torus.
        """
        if self.dimension != _PLANE:
            msg = (
                f"a configuration picture needs a two-dimensional automaton; "
                f"got dimension {self.dimension}"
            )
            raise ValueError(msg)
        cells, periods = tuple(config), tuple(shape)
        self._check_config(cells, periods)
        columns = periods[1]
        ax = ensure_axes(ax)
        ax.imshow(
            [
                list(cells[start : start + columns])
                for start in range(0, len(cells), columns)
            ],
            cmap="binary",
            interpolation="nearest",
            vmin=0,
            vmax=max(self.states - 1, 1),
            aspect="equal",
        )
        ax.set_axis_off()
        ax.set_title("Configuration")
        return ax

    def plot_de_bruijn(self, ax: Axes | None = None) -> Axes:
        """Draw the labeled de Bruijn graph onto ``ax``, vertices on a circle.

        Kari's lecture notes, section 2.6. Every edge is annotated
        ``word -> label``; loops are written beside their vertex.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.

        Raises:
            ValueError: If the automaton is not one-dimensional.
        """
        edges = self.de_bruijn_graph()
        vertices = sorted({edge[0] for edge in edges})
        layout = dict(zip(vertices, unit_circle(len(vertices)), strict=True))
        ax = ensure_axes(ax)
        for source, target, label in edges:
            (x0, y0), (x1, y1) = layout[source], layout[target]
            text = f"{''.join(map(str, (*source, *target[-1:])))} -> {label}"
            if source == target:
                ax.annotate(text, (1.25 * x0, 1.25 * y0), ha="center", fontsize=7)
                continue
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops={
                    "arrowstyle": "->",
                    "connectionstyle": "arc3,rad=0.15",
                    "color": "0.4",
                },
            )
            ax.annotate(
                text,
                (0.6 * x0 + 0.4 * x1, 0.6 * y0 + 0.4 * y1),
                ha="center",
                fontsize=7,
            )
        for vertex, (x, y) in layout.items():
            ax.annotate(
                "".join(map(str, vertex)) or "()",
                (x, y),
                ha="center",
                va="center",
                fontsize=9,
                bbox={"boxstyle": "round", "facecolor": "white"},
            )
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title("De Bruijn graph")
        return ax


# --------------------------------------------------------------------------- #
# Shifts
# --------------------------------------------------------------------------- #
def shift_automaton(dimension: int, states: int, axis: int) -> CellularAutomaton:
    """Return the shift ``sigma_i``: neighborhood ``(e_i)``, identity local rule.

    Kari 2005, section 2.4: translations are compositions of shifts, and
    every cellular automaton commutes with them.

    Raises:
        ValueError: If ``axis`` is not an axis of ``Z^d``, or (D1) fails.
    """
    _require_dimension(dimension)
    if not 0 <= axis < dimension:
        msg = f"Z^{dimension} has axes 0..{dimension - 1}; got {axis}"
        raise ValueError(msg)
    unit = tuple(int(i == axis) for i in range(dimension))
    return CellularAutomaton(dimension, states, (unit,), tuple(range(states)))


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
def reversible_elementary_rules() -> tuple[CellularAutomaton, ...]:
    """Return rules 204, 51, 170, 240, 85, 15: the only reversible elementary CA.

    Cellular Automaton page, canonical examples: identity, complement, the
    two shifts (170 is Kari's left shift, 240 is ``f = a``) and their
    complements — the trivial reversible class.
    """
    return tuple(
        CellularAutomaton.from_wolfram_number(number)
        for number in (204, 51, 170, 240, 85, 15)
    )


def rule_30() -> CellularAutomaton:
    """Return rule 30: surjective, not injective (page, canonical examples)."""
    return CellularAutomaton.from_wolfram_number(30)


def rule_90() -> CellularAutomaton:
    """Return rule 90, ``f(a, b, c) = a + c mod 2``.

    Additive and surjective on the line, yet with Gardens of Eden on every
    ring (Martin-Odlyzko-Wolfram 1984, section 3): the page's
    line-versus-ring fixture.
    """
    return CellularAutomaton.from_wolfram_number(90)


def rule_102() -> CellularAutomaton:
    """Return rule 102, ``f = b + c mod 2``, the XOR automaton.

    Pre-injective but not injective (Ceccherini-Silberstein-Coornaert 2017,
    Examples 3.2, 3.15).
    """
    return CellularAutomaton.from_wolfram_number(102)


def rule_110() -> CellularAutomaton:
    """Return rule 110 (Cook 2004, section 1): universal, not surjective."""
    return CellularAutomaton.from_wolfram_number(110)


def rule_150() -> CellularAutomaton:
    """Return rule 150, ``f = a + b + c mod 2``.

    Bijective on the ring ``Z/N`` exactly when 3 does not divide ``N``
    (page, canonical examples; consistent with Martin-Odlyzko-Wolfram
    Theorem 4.2).
    """
    return CellularAutomaton.from_wolfram_number(150)


def rule_184() -> CellularAutomaton:
    """Return rule 184, the traffic rule, mirror image of rule 226."""
    return CellularAutomaton.from_wolfram_number(184)


def rule_226() -> CellularAutomaton:
    """Return rule 226, which "replaces pattern 01 by pattern 10".

    Kari's lecture notes, Example 26; number-conserving.
    """
    return CellularAutomaton.from_wolfram_number(226)


def rule_232() -> CellularAutomaton:
    """Return rule 232, majority vote: balanced table, not surjective.

    Ceccherini-Silberstein-Coornaert 2017, Examples 3.12, 3.17.
    """
    return CellularAutomaton.from_wolfram_number(232)


def radius_half_xor() -> CellularAutomaton:
    """Return the radius-1/2 XOR, ``f(x, y) = x + y mod 2`` on ``N = (0, 1)``.

    Kari 2005, section 3: ``G_F`` is injective but not surjective.
    """
    return CellularAutomaton.from_local_rule(
        1, _BINARY, [(0,), (1,)], lambda cells: sum(cells) % _BINARY
    )


def _patt_flip(cells: tuple[int, ...]) -> int:
    """Flip the second cell exactly in the landscape ``0 . 1 0``."""
    a, b, c, d = cells
    return b ^ int((a, c, d) == (0, 1, 0))


def patt_rule() -> CellularAutomaton:
    """Return Patt's rule on the four contiguous cells ``(-1, 0, 1, 2)``.

    Patt 1971, via Amoroso-Patt 1972, section III and Appendix A: binary,
    the second cell flips exactly in the landscape ``0 . 1 0`` — the
    smallest contiguous neighborhood carrying a nontrivial reversible
    binary one-dimensional automaton.
    """
    return CellularAutomaton.from_block_map(_BINARY, 1, 2, _patt_flip)


def patt_table_on_gapped_neighborhood() -> CellularAutomaton:
    """Return Patt's table read on the neighborhood ``(-4, -2, -1, 0)``.

    Amoroso-Patt 1972, section VI: the same table is not injective there,
    so injectivity depends on the shape of the neighborhood.
    """
    return CellularAutomaton.from_local_rule(
        1, _BINARY, [(-4,), (-2,), (-1,), (0,)], _patt_flip
    )


def one_way_nand() -> CellularAutomaton:
    """Return the one-way NAND automaton on ``N = (-1, 0)``.

    The page's bridge example: on a ring of 3 cells it is the NAND network
    on the directed 3-cycle of Aledo-Martinez-Valverde 2015, Example 7.
    """
    return CellularAutomaton.from_local_rule(
        1, _BINARY, [(-1,), (0,)], lambda cells: 1 - cells[0] * cells[1]
    )


_LIFE_CENTER = 4
_LIFE_SURVIVE = (2, 3)
_LIFE_BIRTH = 3


def _life_rule(cells: tuple[int, ...]) -> int:
    neighbors = sum(cells) - cells[_LIFE_CENTER]
    if cells[_LIFE_CENTER]:
        return int(neighbors in _LIFE_SURVIVE)
    return int(neighbors == _LIFE_BIRTH)


def game_of_life() -> CellularAutomaton:
    """Return Conway's Game of Life as a tessellation structure with quiescent 0.

    Gardner 1970: two states, Moore neighborhood; a live cell survives with
    2 or 3 live neighbors and a dead cell is born with exactly 3.
    """
    return CellularAutomaton.from_tessellation(_PLANE, _BINARY, 0, _life_rule)


def _alive(cells: Iterable[Offset]) -> dict[Offset, int]:
    return dict.fromkeys(cells, 1)


def life_block() -> dict[Offset, int]:
    """Return the block, a still life (cells as ``(row, column)``)."""
    return _alive([(0, 0), (0, 1), (1, 0), (1, 1)])


def life_beehive() -> dict[Offset, int]:
    """Return the beehive, a still life."""
    return _alive([(0, 1), (0, 2), (1, 0), (1, 3), (2, 1), (2, 2)])


def life_blinker() -> dict[Offset, int]:
    """Return the blinker, of period 2."""
    return life_row(3)


def life_glider() -> dict[Offset, int]:
    """Return the glider: period 4, one cell diagonally per period."""
    return _alive([(0, 1), (1, 2), (2, 0), (2, 1), (2, 2)])


def life_t_tetromino() -> dict[Offset, int]:
    """Return the T-tetromino, which becomes traffic lights at generation 9."""
    return _alive([(0, 0), (0, 1), (0, 2), (1, 1)])


def life_r_pentomino() -> dict[Offset, int]:
    """Return the R-pentomino: generation 1103 has population 116 (Life Lexicon)."""
    return _alive([(0, 1), (0, 2), (1, 0), (1, 1), (2, 1)])


def life_row(length: int) -> dict[Offset, int]:
    """Return a horizontal row of ``length`` live cells.

    Page, canonical examples: a row of 7 becomes a honey farm, a row of 10
    the pentadecathlon.

    Raises:
        ValueError: If ``length`` is negative.
    """
    if length < 0:
        msg = f"a row has a non-negative length; got {length}"
        raise ValueError(msg)
    return _alive((0, column) for column in range(length))
