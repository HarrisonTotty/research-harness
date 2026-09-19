r"""Young diagrams: finite order ideals of the quarter-plane grid.

A Young diagram is the left-justified array of boxes of an integer partition
``lambda_1 >= lambda_2 >= ...``, with ``lambda_i`` boxes in row ``i``
(Stanley, *Enumerative Combinatorics* vol. 1, 2nd ed., section 1.7, p. 65);
equivalently a finite down-closed set of cells of ``P x P`` (Stanley EC1
section 1.5, p. 50). :class:`YoungDiagram` stores the positive row lengths
and derives the rest of the Young Diagram page as views — the cell set, the
conjugate, hooks and rims, Frobenius coordinates, the boundary path and
Postnikov's box embedding ``I(lambda)``, the lattice operations of Young's
lattice and of the dominance order, and the tableau counts — next to the
page's enumerators and canonical examples.

Conventions, all following the page:

* **Coordinates.** Cells are matrix-style ``(row, column)``, 1-indexed, in
  English notation (Stanley, Frame-Robinson-Thrall, Adin-Roichman,
  Postnikov, and :mod:`research.le_diagram`). Mathlib's ``YoungDiagram`` is
  0-indexed; it is not followed here.
* **Zeros.** Trailing zero parts are ignored, ``(3, 3, 2, 1, 0, 0) =
  (3, 3, 2, 1)`` (Stanley EC1 p. 65), so only positive rows are stored and
  equality is canonical. The box-aware methods take ``k`` and ``n`` as
  arguments and pad with empty rows on demand (Postnikov section 2.1 writes
  ``lambda`` with exactly ``k`` parts).
* **Frobenius coordinates** are arms first, legs second (Stanley EC1
  Exercise 1.70(a)); Brunat-Nath write the reverse order.
* **Rank.** ``durfee_side`` is the rank of the partition (the Durfee side);
  ``size`` is its poset rank in Young's lattice.
"""

import enum
import functools
import itertools
import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, override

import pandas as pd

from research._linalg import det_q
from research._plot import ensure_axes
from research.le_diagram import LeDiagram

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "DrawingConvention",
    "TableauFormula",
    "YoungDiagram",
    "adin_roichman_example_9_7_diagram",
    "all_in_one_diagram",
    "box_size_distribution",
    "count_multichains",
    "diagrams_in_box",
    "empty_diagram",
    "frt_example_1_diagram",
    "frt_example_2_diagram",
    "frt_example_4_diagram",
    "hook_shape",
    "mobius",
    "non_decreasing_rows",
    "non_down_closed_cells",
    "partitions",
    "postnikov_figure_2_1_diagram",
    "rectangle",
    "skew_non_example_cells",
    "staircase",
    "stanley_figure_1_16_diagram",
    "stanley_figure_1_18_diagram",
    "stanley_figure_1_33_diagram",
]

type Cell = tuple[int, int]
type Tableau = tuple[tuple[int, ...], ...]

_COLUMNS = ("row", "column")


class DrawingConvention(enum.Enum):
    """The three ways of drawing a diagram (Adin-Roichman section 2.4).

    * ``ENGLISH`` — "row indices increase from top to bottom and column
      indices increase from left to right."
    * ``FRENCH`` — "row indices increase from bottom to top (and column
      indices increase from left to right)."
    * ``RUSSIAN`` — "rotated 45 degrees": the French picture turned so that
      the corner cell sits at the bottom of a V.
    """

    ENGLISH = "english"
    FRENCH = "french"
    RUSSIAN = "russian"


class TableauFormula(enum.Enum):
    """The page's equivalent ways of counting standard Young tableaux.

    * ``HOOK`` — the hook length formula ``|lambda|! / prod h_c``
      (Frame-Robinson-Thrall 1954, Theorem 1; Adin-Roichman Theorem 5.3).
    * ``PRODUCT`` — the Frobenius-Young product formula over the
      first-column hook lengths (Adin-Roichman Theorem 5.1).
    * ``DETERMINANT`` — ``|lambda|! det[1 / (lambda_i - i + j)!]``
      (Adin-Roichman Theorem 5.4).
    * ``CHAINS`` — the corner-removal recursion ``e(A) = sum e(A')`` over
      the diagrams ``A'`` covered by ``A`` (Stanley EC1 section 3.5,
      eq. 3.11), with ``e(lambda) = f^lambda`` by Adin-Roichman
      section 2.5.1.
    """

    HOOK = "hook"
    PRODUCT = "product"
    DETERMINANT = "determinant"
    CHAINS = "chains"


def _conjugate_rows(rows: Sequence[int]) -> tuple[int, ...]:
    """Return ``lambda'_j = |{i : lambda_i >= j}|`` (Adin-Roichman section 2.4)."""
    width = rows[0] if rows else 0
    return tuple(sum(1 for length in rows if length >= j) for j in range(1, width + 1))


@functools.cache
def _count_chains(rows: tuple[int, ...]) -> int:
    """Return ``e(lambda)`` by corner removal (Stanley EC1 eq. 3.11)."""
    if not rows:
        return 1
    total = 0
    for index, length in enumerate(rows):
        below = rows[index + 1] if index + 1 < len(rows) else 0
        if length > below:
            smaller = (*rows[:index], length - 1, *rows[index + 1 :])
            total += _count_chains(tuple(r for r in smaller if r))
    return total


def _binomial(a: int, b: int) -> int:
    """Return ``binom(a, b)``, zero when ``a < 0`` or ``b`` is out of range."""
    if a < 0 or b < 0 or b > a:
        return 0
    return math.comb(a, b)


def _poly_mul(p: Sequence[int], q: Sequence[int]) -> list[int]:
    """Return the product of two coefficient lists (index = power)."""
    product = [0] * (len(p) + len(q) - 1)
    for a, x in enumerate(p):
        for b, y in enumerate(q):
            product[a + b] += x * y
    return product


def _poly_div_exact(numerator: Sequence[int], denominator: Sequence[int]) -> list[int]:
    """Return the quotient of an exact division by a monic polynomial.

    Raises:
        ArithmeticError: If the division leaves a remainder.
    """
    remainder = list(numerator)
    degree = len(denominator) - 1
    quotient = [0] * (len(remainder) - degree)
    for power in range(len(quotient) - 1, -1, -1):
        factor = remainder[power + degree]
        quotient[power] = factor
        for offset, coefficient in enumerate(denominator):
            remainder[power + offset] -= factor * coefficient
    if any(remainder):
        msg = f"{list(numerator)!r} is not divisible by {list(denominator)!r}"
        raise ArithmeticError(msg)
    return quotient


def _check_box(k: int, n: int) -> None:
    """Reject boxes outside ``0 <= k <= n``."""
    if not 0 <= k <= n:
        msg = f"box violated: a k x (n - k) box needs 0 <= k <= n; got k = {k}, n = {n}"
        raise ValueError(msg)


# --------------------------------------------------------------------------- #
# The structure
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class YoungDiagram:
    """A Young diagram, stored as its positive, weakly decreasing row lengths.

    ``rows`` lists ``lambda_1 >= ... >= lambda_t > 0`` (Stanley EC1
    section 1.7, p. 65; Adin-Roichman Definitions 2.9-2.10). The empty
    tuple is the empty diagram, the partition of 0.

    Like :class:`~research.le_diagram.LeDiagram`, the definition is one
    linear scan, so ``__post_init__`` checks it and direct construction is
    validated too. Direct construction rejects zero rows so that equal
    diagrams are equal values; :meth:`from_rows` strips them instead, and
    the other ``from_*`` classmethods convert the page's equivalent
    formulations.

    Raises:
        ValueError: Naming the violated clause — weak decrease or
            positivity.
    """

    rows: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Validate weak decrease and positivity of the row lengths."""
        rows = self.rows
        if any(a < b for a, b in itertools.pairwise(rows)):
            msg = (
                f"weak decrease violated: no row may be longer than the row "
                f"above it (Stanley EC1 section 1.7); got {rows!r}"
            )
            raise ValueError(msg)
        if any(length <= 0 for length in rows):
            msg = (
                f"positivity violated: stored rows are the positive parts "
                f"(use from_rows to strip zeros); got {rows!r}"
            )
            raise ValueError(msg)

    @override
    def __str__(self) -> str:
        """Return the English picture, ``[]`` per cell; ``()`` when empty."""
        return "\n".join("[]" * length for length in self.rows) or "()"

    # -------------------------------------------------------- constructors
    @classmethod
    def from_rows(cls, rows: Iterable[int]) -> YoungDiagram:
        """Build from the parts of a partition, zeros allowed and ignored.

        Stanley EC1 section 1.7, p. 65: a partition is a sequence with
        ``lambda_1 >= lambda_2 >= ... >= 0``, and trailing zeros are
        ignored — formulations (A) and (B) of the page.

        Raises:
            ValueError: If the parts are not weakly decreasing, or one is
                negative.
        """
        parts = tuple(rows)
        if any(part < 0 for part in parts):
            msg = (
                f"positivity violated: the parts of a partition are "
                f"non-negative (Stanley EC1 section 1.7); got {parts!r}"
            )
            raise ValueError(msg)
        if any(a < b for a, b in itertools.pairwise(parts)):
            msg = (
                f"weak decrease violated: no row may be longer than the row "
                f"above it (Stanley EC1 section 1.7); got {parts!r}"
            )
            raise ValueError(msg)
        return cls(tuple(part for part in parts if part))

    @classmethod
    def from_cells(cls, cells: Iterable[Cell]) -> YoungDiagram:
        """Build from a finite down-closed set of cells of ``P x P``.

        Formulation (C) of the page (Stanley EC1 section 1.5, p. 50): if
        ``(i, j)`` is a cell and ``i' <= i``, ``j' <= j`` then ``(i', j')``
        is a cell. It suffices that the cell above and the cell to the left
        of every cell are present.

        Raises:
            ValueError: If a cell lies outside the positive quarter-plane,
                or the set is not down-closed.
        """
        cell_set = frozenset(cells)
        outside = sorted(cell for cell in cell_set if cell[0] < 1 or cell[1] < 1)
        if outside:
            msg = (
                f"quarter-plane violated: cells are 1-indexed pairs of "
                f"positive integers; got {outside!r}"
            )
            raise ValueError(msg)
        for i, j in sorted(cell_set):
            for neighbor in ((i - 1, j), (i, j - 1)):
                if min(neighbor) >= 1 and neighbor not in cell_set:
                    msg = (
                        f"down-closure violated: {(i, j)!r} is a cell but "
                        f"{neighbor!r} is not (Stanley EC1 section 1.5)"
                    )
                    raise ValueError(msg)
        height = max((i for i, _ in cell_set), default=0)
        return cls(
            tuple(
                sum(1 for row, _ in cell_set if row == i) for i in range(1, height + 1)
            )
        )

    @classmethod
    def from_boundary_sequence(cls, sequence: Iterable[int]) -> YoungDiagram:
        """Build from a 0/1 boundary sequence read from the SW to the NE corner.

        Adin-Roichman Definition 9.6: each east step is a 1 and each north
        step a 0. The row lengths, bottom to top, are the numbers of 1's
        before each 0 (the page's "Equivalence of the formulations"); 0's
        before the first 1 are empty rows, and 1's after the last 0 bound
        no row, so neither changes the diagram.

        Raises:
            ValueError: If an entry is not 0 or 1.
        """
        steps = tuple(sequence)
        bad = sorted(set(steps) - {0, 1})
        if bad:
            msg = (
                f"boundary sequence violated: steps are 0 (north) or 1 "
                f"(east) (Adin-Roichman Definition 9.6); got {bad!r}"
            )
            raise ValueError(msg)
        lengths: list[int] = []
        ones = 0
        for step in steps:
            if step:
                ones += 1
            else:
                lengths.append(ones)
        return cls.from_rows(reversed(lengths))

    @classmethod
    def from_frobenius(cls, arms: Iterable[int], legs: Iterable[int]) -> YoungDiagram:
        """Build from diagonal (Frobenius) coordinates ``(a_1 ... a_r | b_1 ... b_r)``.

        Formulation (E) of the page (Stanley EC1 Exercise 1.70(a), which
        attributes the bijection to Frobenius 1900): ``a_i = lambda_i - i``
        and ``b_i = lambda'_i - i``, arms first.

        Raises:
            ValueError: Unless both rows have the same length and are
                strictly decreasing and non-negative.
        """
        a, b = tuple(arms), tuple(legs)
        strict = all(
            x > y for row in (a, b) for x, y in itertools.pairwise(row)
        ) and all(x >= 0 for x in (*a, *b))
        if len(a) != len(b) or not strict:
            msg = (
                f"Frobenius coordinates violated: need a_1 > ... > a_r >= 0 "
                f"and b_1 > ... > b_r >= 0 of equal length (Stanley EC1 "
                f"Exercise 1.70(a)); got ({a!r} | {b!r})"
            )
            raise ValueError(msg)
        rank = len(a)
        upper = [arm + i for i, arm in enumerate(a, start=1)]
        depth = b[0] + 1 if b else 0
        lower = [
            sum(1 for j, leg in enumerate(b, start=1) if leg + j >= i)
            for i in range(rank + 1, depth + 1)
        ]
        return cls(tuple(upper + lower))

    @classmethod
    def from_subset(cls, subset: Iterable[int], n: int) -> YoungDiagram:
        r"""Build the diagram ``lambda(I)`` of a subset ``I`` of ``[n]``.

        Postnikov (arXiv:math/0609764) section 2.1, the inverse of
        :meth:`to_subset`: ``lambda_s = |[i_s, n] \ I|`` for the ``s``-th
        smallest element ``i_s``. The diagram fits the ``k x (n - k)`` box
        with ``k = |I|``.

        Raises:
            ValueError: If ``I`` is not a subset of ``{1, ..., n}``.
        """
        members = frozenset(subset)
        if not members <= frozenset(range(1, n + 1)):
            msg = (
                f"subset violated: I must be a subset of [n] with n = {n} "
                f"(Postnikov section 2.1); got {sorted(members)!r}"
            )
            raise ValueError(msg)
        return cls.from_rows(
            sum(1 for label in range(i, n + 1) if label not in members)
            for i in sorted(members)
        )

    @classmethod
    def from_durfee_dissection(
        cls, side: int, right: YoungDiagram, below: YoungDiagram
    ) -> YoungDiagram:
        """Build the diagram with Durfee side ``k``, ``mu`` to its right, ``nu`` below.

        Stanley EC1, proof of Proposition 1.8.6(b): "every lambda of rank k
        is obtained uniquely from such mu and nu", where ``l(mu) <= k`` and
        ``nu_1 <= k``; then ``|lambda| = k^2 + |mu| + |nu|``.

        Raises:
            ValueError: If ``k < 0``, ``mu`` has more than ``k`` parts, or
                ``nu`` has a part larger than ``k``.
        """
        if side < 0 or right.length > side or below.largest_part > side:
            msg = (
                f"Durfee dissection violated: need k >= 0, l(mu) <= k and "
                f"nu_1 <= k (Stanley EC1 Proposition 1.8.6(b)); got k = "
                f"{side}, mu = {right.rows!r}, nu = {below.rows!r}"
            )
            raise ValueError(msg)
        top = tuple(side + part for part in right.padded_rows(side))
        return cls(top + below.rows)

    @classmethod
    def from_le_diagram(cls, diagram: LeDiagram) -> YoungDiagram:
        """Return the shape of a Le-diagram, its empty rows dropped.

        Postnikov Theorem 6.5: "the Le-diagram D has shape lambda if and
        only if" its cell lies in the Schubert cell ``Omega_lambda``.
        """
        return cls.from_rows(diagram.shape)

    # ------------------------------------------------- computed properties
    @property
    def size(self) -> int:
        """The number ``|lambda|`` of cells (Stanley EC1 p. 65)."""
        return sum(self.rows)

    @property
    def length(self) -> int:
        """The number ``l(lambda)`` of rows, i.e. of parts (Stanley EC1 p. 65)."""
        return len(self.rows)

    @property
    def largest_part(self) -> int:
        """The largest part ``lambda_1`` — the width; 0 for the empty diagram."""
        return self.rows[0] if self.rows else 0

    @property
    def multiplicities(self) -> dict[int, int]:
        """The numbers ``m_i`` of parts equal to ``i``, for the ``i`` that occur.

        Stanley EC1 eq. 1.74: ``lambda = <1^{m_1}, 2^{m_2}, ...>``.
        """
        counts: dict[int, int] = {}
        for part in self.rows:
            counts[part] = counts.get(part, 0) + 1
        return counts

    @functools.cached_property
    def cells(self) -> frozenset[Cell]:
        """The cell set ``[lambda] = {(i, j) : 1 <= j <= lambda_i}``.

        Adin-Roichman Definition 2.10; Postnikov section 2.1.
        """
        return frozenset(
            (i, j)
            for i, length in enumerate(self.rows, start=1)
            for j in range(1, length + 1)
        )

    @functools.cached_property
    def _columns(self) -> tuple[int, ...]:
        """The column lengths ``lambda'``."""
        return _conjugate_rows(self.rows)

    @property
    def is_self_conjugate(self) -> bool:
        """Whether ``lambda = lambda'`` (Stanley EC1 section 1.8)."""
        return self._columns == self.rows

    def _require_cell(self, row: int, column: int) -> None:
        """Reject a cell that is not in the diagram."""
        if not (1 <= row <= len(self.rows) and 1 <= column <= self.rows[row - 1]):
            msg = f"{(row, column)!r} is not a cell of the diagram {self.rows!r}"
            raise ValueError(msg)

    def arm_length(self, row: int, column: int) -> int:
        """Return the number ``lambda_i - j`` of cells to the right of ``(i, j)``.

        Frame-Robinson-Thrall (*The hook graphs of the symmetric group*,
        1954), section 1.

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        self._require_cell(row, column)
        return self.rows[row - 1] - column

    def leg_length(self, row: int, column: int) -> int:
        """Return the number ``lambda'_j - i`` of cells below ``(i, j)``.

        Frame-Robinson-Thrall 1954, section 1.

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        self._require_cell(row, column)
        return self._columns[column - 1] - row

    def hook_length(self, row: int, column: int) -> int:
        """Return ``h_ij = 1 + (lambda_i - j) + (lambda'_j - i)``.

        Frame-Robinson-Thrall 1954, eq. 1.1; Adin-Roichman Definition 5.2.

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        return 1 + self.arm_length(row, column) + self.leg_length(row, column)

    def hook(self, row: int, column: int) -> frozenset[Cell]:
        """Return the hook of a cell: it, its arm and its leg.

        Frame-Robinson-Thrall 1954, section 1: "this node and all nodes to
        the right of it or below it".

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        arm = self.arm_length(row, column)
        leg = self.leg_length(row, column)
        return frozenset(
            [(row, column + step) for step in range(arm + 1)]
            + [(row + step, column) for step in range(1, leg + 1)]
        )

    @functools.cached_property
    def hook_lengths(self) -> tuple[tuple[int, ...], ...]:
        """The hook graph ``H[lambda]``: the hook lengths, row by row.

        Frame-Robinson-Thrall 1954, section 1. Costs ``O(|lambda|)``.
        """
        columns = self._columns
        return tuple(
            tuple(1 + (length - j) + (columns[j - 1] - i) for j in range(1, length + 1))
            for i, length in enumerate(self.rows, start=1)
        )

    @property
    def hook_product(self) -> int:
        """The hook product ``H_lambda`` (Frame-Robinson-Thrall 1954, Theorem 1)."""
        return math.prod(h for row in self.hook_lengths for h in row)

    @property
    def diagonal_hook_lengths(self) -> tuple[int, ...]:
        """The hook lengths ``h_ii`` of the main-diagonal cells, top to bottom.

        Stanley EC1 Proposition 1.8.4: on a self-conjugate diagram they are
        distinct and odd, which gives the bijection onto partitions into
        distinct odd parts.
        """
        return tuple(self.hook_lengths[i][i] for i in range(self.durfee_side))

    @property
    def corners(self) -> tuple[Cell, ...]:
        """The removable cells, top to bottom.

        Adin-Roichman section 6.1: "a cell which is last in its row and in
        its column (equivalently, has hook length 1)".
        """
        rows = self.rows
        return tuple(
            (i, length)
            for i, length in enumerate(rows, start=1)
            if i == len(rows) or length > rows[i]
        )

    @property
    def addable_cells(self) -> tuple[Cell, ...]:
        """The positions where one box can be added, top to bottom.

        Stanley EC1 section 3.4 applied to Young's lattice (the page's
        "Adding or removing a box"): the minimal elements of the
        complement. They alternate with the corners along the boundary
        (Stanley EC1 Example 3.21.2(2)).
        """
        rows = self.rows
        inner = [
            (i, length + 1)
            for i, length in enumerate(rows, start=1)
            if i == 1 or rows[i - 2] > length
        ]
        return (*inner, (len(rows) + 1, 1))

    @property
    def rim(self) -> frozenset[Cell]:
        """The rim nodes: cells ``(s, t)`` with ``(s + 1, t + 1)`` not a cell.

        Frame-Robinson-Thrall 1954, section 2.
        """
        return frozenset(
            (s, t) for s, t in self.cells if (s + 1, t + 1) not in self.cells
        )

    def rim_hook(self, row: int, column: int) -> frozenset[Cell]:
        """Return the rim hook of a cell: ``h_ij`` rim nodes from head to foot.

        Frame-Robinson-Thrall 1954, section 2: the rim nodes running from
        the head ``(i, lambda_i)`` of the ``(i, j)``-hook to its foot
        ``(lambda'_j, j)`` — the rim nodes weakly below row ``i`` and
        weakly right of column ``j``.

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        self._require_cell(row, column)
        return frozenset((s, t) for s, t in self.rim if s >= row and t >= column)

    @property
    def durfee_side(self) -> int:
        """The rank of the partition: the side of its Durfee square.

        Stanley EC1 section 1.8, pp. 71-72: "the largest i for which
        lambda_i >= i". Unrelated to the poset rank :attr:`size`.
        """
        return sum(1 for i, length in enumerate(self.rows, start=1) if length >= i)

    @property
    def frobenius_coordinates(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        """The diagonal coordinates ``(a_1 ... a_r | b_1 ... b_r)``, arms first.

        Formulation (E) of the page (Stanley EC1 Exercise 1.70(a)):
        ``a_i = lambda_i - i`` and ``b_i = lambda'_i - i`` for ``i`` up to
        the Durfee side ``r``; ``r + sum(a_i + b_i) = |lambda|``.
        """
        rank = self.durfee_side
        arms = tuple(self.rows[i] - (i + 1) for i in range(rank))
        legs = tuple(self._columns[i] - (i + 1) for i in range(rank))
        return arms, legs

    @property
    def boundary_sequence(self) -> tuple[int, ...]:
        """The boundary sequence from the SW to the NE corner of the diagram.

        Adin-Roichman Definition 9.6: each east step a 1, each north step a
        0. It starts with a 1 and ends with a 0 unless the diagram is
        empty, whose sequence is empty.
        """
        steps: list[int] = []
        below = 0
        for length in reversed(self.rows):
            steps.extend([1] * (length - below))
            steps.append(0)
            below = length
        return tuple(steps)

    def fits_in_box(self, k: int, n: int) -> bool:
        """Return whether the diagram fits the ``k x (n - k)`` rectangle.

        Postnikov section 2.1: ``lambda`` is contained in ``(n - k)^k``,
        i.e. has at most ``k`` parts and largest part at most ``n - k``.

        Raises:
            ValueError: Unless ``0 <= k <= n``.
        """
        _check_box(k, n)
        return self.length <= k and self.largest_part <= n - k

    def padded_rows(self, k: int) -> tuple[int, ...]:
        """Return exactly ``k`` row lengths, padding with empty rows.

        Postnikov section 2.1 writes ``lambda`` with exactly ``k`` parts;
        this is the shape convention of :mod:`research.le_diagram`.

        Raises:
            ValueError: If the diagram has more than ``k`` rows.
        """
        if self.length > k:
            msg = f"box violated: {self.rows!r} has more than k = {k} rows"
            raise ValueError(msg)
        return self.rows + (0,) * (k - self.length)

    def _require_box(self, k: int, n: int) -> None:
        """Reject a box the diagram does not fit."""
        if not self.fits_in_box(k, n):
            msg = (
                f"box violated: {self.rows!r} does not fit the k x (n - k) "
                f"rectangle with k = {k}, n = {n} (Postnikov section 2.1)"
            )
            raise ValueError(msg)

    def to_subset(self, k: int, n: int) -> frozenset[int]:
        """Return Postnikov's ``I(lambda)``, a ``k``-subset of ``[n]``.

        Postnikov section 2.1: label the ``n`` steps of the boundary path
        from the upper-right to the lower-left corner of the rectangle by
        ``1..n``; ``I(lambda)`` is the set of labels of the ``k`` vertical
        steps. In closed form ``i_s = (n - k) - lambda_s + s`` (Stanley EC1
        section 1.7, p. 67: "positions ``j - lambda_i + i``").

        Raises:
            ValueError: If the diagram does not fit the box.
        """
        self._require_box(k, n)
        return frozenset(
            (n - k) - length + s
            for s, length in enumerate(self.padded_rows(k), start=1)
        )

    def boundary_word(self, k: int, n: int) -> tuple[int, ...]:
        """Return Stanley's boundary word of the diagram in the box.

        Stanley EC1 section 1.7, p. 67: walk "from the upper right-hand
        corner to the lower left-hand corner of the rectangle", writing 1
        for a horizontal and 2 for a vertical step. The word has ``|lambda|``
        inversions.

        Raises:
            ValueError: If the diagram does not fit the box.
        """
        vertical = self.to_subset(k, n)
        return tuple(2 if label in vertical else 1 for label in range(1, n + 1))

    def contains(self, other: YoungDiagram) -> bool:
        """Return whether ``other`` is contained in this diagram.

        The order of Young's lattice (Stanley EC1 Example 3.4.4(b)):
        ``mu_i <= lambda_i`` for all ``i``.
        """
        return other.length <= self.length and all(
            small <= large for small, large in zip(other.rows, self.rows, strict=False)
        )

    def _require_same_size(self, other: YoungDiagram) -> None:
        """Reject a dominance comparison between partitions of different sizes."""
        if self.size != other.size:
            msg = (
                f"the dominance order compares partitions of the same n "
                f"(Grinberg-Reiner Definition 2.2.7); got {self.rows!r} of "
                f"size {self.size} and {other.rows!r} of size {other.size}"
            )
            raise ValueError(msg)

    def dominates(self, other: YoungDiagram) -> bool:
        """Return whether this partition dominates ``other``.

        Grinberg-Reiner (*Hopf algebras in combinatorics*) Definition
        2.2.7: ``lambda_1 + ... + lambda_k >= mu_1 + ... + mu_k`` for all
        ``k``.

        Raises:
            ValueError: If the sizes differ — the order is defined on the
                partitions of one ``n``.
        """
        self._require_same_size(other)
        pairs = itertools.zip_longest(
            itertools.accumulate(self.rows),
            itertools.accumulate(other.rows),
            fillvalue=self.size,
        )
        return all(mine >= theirs for mine, theirs in pairs)

    def count_standard_tableaux(
        self, formula: TableauFormula = TableauFormula.HOOK
    ) -> int:
        """Return the number ``f^lambda`` of standard Young tableaux.

        The four formulas of :class:`TableauFormula` agree (Adin-Roichman
        Claim 5.5 for the three closed forms; Adin-Roichman section 2.5.1
        with Stanley EC1 eq. 3.11 for the recursion); each is evaluated
        literally and exactly, so their agreement is a checkable fact.
        ``CHAINS`` is memoized across calls and costs one step per
        subdiagram.
        """
        n, t = self.size, self.length
        match formula:
            case TableauFormula.HOOK:
                return math.factorial(n) // self.hook_product
            case TableauFormula.PRODUCT:
                first = [length + t - i for i, length in enumerate(self.rows, start=1)]
                differences = math.prod(
                    a - b for a, b in itertools.combinations(first, 2)
                )
                factorials = math.prod(math.factorial(h) for h in first)
                return math.factorial(n) * differences // factorials
            case TableauFormula.DETERMINANT:
                matrix = [
                    [
                        Fraction(1, math.factorial(length - i + j))
                        if length - i + j >= 0
                        else Fraction(0)
                        for j in range(1, t + 1)
                    ]
                    for i, length in enumerate(self.rows, start=1)
                ]
                return int(math.factorial(n) * det_q(matrix))
            case TableauFormula.CHAINS:
                return _count_chains(self.rows)

    def standard_tableaux(self) -> Iterator[Tableau]:
        """Yield every standard Young tableau of this shape, lazily.

        Adin-Roichman Definition 2.3: an order-preserving bijection from
        the cells to ``1..n``; equivalently a maximal chain from the empty
        diagram in Young's lattice (section 2.5.1), which is how they are
        generated — the largest entry sits in a corner. Each tableau is a
        tuple of rows; the iterator is single-use and has ``f^lambda``
        items.
        """
        if not self.rows:
            yield ()
            return
        n = self.size
        for i, _ in self.corners:
            for smaller in self.remove_cell((i, self.rows[i - 1])).standard_tableaux():
                padded = [*smaller, *([()] * (self.length - len(smaller)))]
                padded[i - 1] = (*padded[i - 1], n)
                yield tuple(padded)

    def major_index_polynomial(self) -> tuple[int, ...]:
        r"""Return the coefficients of ``sum_T q^maj(T)``, index = power.

        Adin-Roichman Theorem 10.26 (citing Stanley EC2 Corollary 7.21.5):
        the sum over standard tableaux equals
        ``q^{sum_i binom(lambda'_i, 2)} [n]_q! / prod_c [h_c]_q``. Evaluated
        by exact polynomial division.
        """
        numerator = [1]
        for m in range(1, self.size + 1):
            numerator = _poly_mul(numerator, [1] * m)
        denominator = [1]
        for row in self.hook_lengths:
            for h in row:
                denominator = _poly_mul(denominator, [1] * h)
        shift = sum(math.comb(column, 2) for column in self._columns)
        return (*([0] * shift), *_poly_div_exact(numerator, denominator))

    # ------------------------------------------------------ transformations
    def conjugate(self) -> YoungDiagram:
        """Return the conjugate (transpose) ``lambda'``.

        Stanley EC1 section 1.8, p. 68: the diagram obtained "by
        interchanging rows and columns"; ``lambda'_j = |{i : lambda_i >=
        j}|`` (Adin-Roichman section 2.4). "A fundamental involution on the
        set of partitions of n".
        """
        return YoungDiagram(self._columns)

    def meet(self, other: YoungDiagram) -> YoungDiagram:
        """Return the intersection, the meet in Young's lattice.

        Stanley EC1 section 3.4: "the lattice operations on order ideals
        are just ordinary intersection and union"; in row lengths,
        ``min(lambda_i, mu_i)``.
        """
        return YoungDiagram(
            tuple(min(a, b) for a, b in zip(self.rows, other.rows, strict=False))
        )

    def join(self, other: YoungDiagram) -> YoungDiagram:
        """Return the union, the join in Young's lattice.

        Stanley EC1 section 3.4; in row lengths, ``max(lambda_i, mu_i)``.
        """
        pairs = itertools.zip_longest(self.rows, other.rows, fillvalue=0)
        return YoungDiagram(tuple(max(a, b) for a, b in pairs))

    def add_cell(self, cell: Cell) -> YoungDiagram:
        """Return the diagram with one box added at an addable position.

        The covers of ``lambda`` in Young's lattice (Stanley EC1
        section 3.4; the page's "Adding or removing a box").

        Raises:
            ValueError: If ``cell`` is not an addable position.
        """
        if cell not in self.addable_cells:
            msg = (
                f"{cell!r} is not an addable position of {self.rows!r}; "
                f"those are {self.addable_cells!r}"
            )
            raise ValueError(msg)
        row = cell[0]
        return YoungDiagram((*self.rows[: row - 1], cell[1], *self.rows[row:]))

    def remove_cell(self, cell: Cell) -> YoungDiagram:
        """Return the diagram with one corner deleted.

        The elements covered by ``lambda`` in Young's lattice (the page's
        "Adding or removing a box").

        Raises:
            ValueError: If ``cell`` is not a corner.
        """
        if cell not in self.corners:
            msg = (
                f"{cell!r} is not a corner of {self.rows!r}; those are {self.corners!r}"
            )
            raise ValueError(msg)
        row = cell[0]
        shorter = (*self.rows[: row - 1], cell[1] - 1, *self.rows[row:])
        return YoungDiagram(tuple(length for length in shorter if length))

    def upper_covers(self) -> tuple[YoungDiagram, ...]:
        """Return the diagrams covering this one in Young's lattice."""
        return tuple(self.add_cell(cell) for cell in self.addable_cells)

    def lower_covers(self) -> tuple[YoungDiagram, ...]:
        """Return the diagrams this one covers in Young's lattice."""
        return tuple(self.remove_cell(cell) for cell in self.corners)

    def skew_cells(self, inner: YoungDiagram) -> frozenset[Cell]:
        r"""Return the skew diagram ``[lambda / mu] = [lambda] \ [mu]``.

        Adin-Roichman Definition 2.12, for ``mu`` contained in ``lambda``.

        Raises:
            ValueError: If ``inner`` is not contained in this diagram.
        """
        if not self.contains(inner):
            msg = (
                f"a skew diagram lambda/mu needs mu inside lambda "
                f"(Adin-Roichman Definition 2.12); got {inner.rows!r} and "
                f"{self.rows!r}"
            )
            raise ValueError(msg)
        return self.cells - inner.cells

    def durfee_dissection(self) -> tuple[int, YoungDiagram, YoungDiagram]:
        """Return ``(k, mu, nu)``: the Durfee side, what lies right of it, and below.

        Stanley EC1, proof of Proposition 1.8.6(b): ``l(mu) <= k``,
        ``nu_1 <= k`` and ``|lambda| = k^2 + |mu| + |nu|``. Inverted by
        :meth:`from_durfee_dissection`.
        """
        side = self.durfee_side
        right = YoungDiagram.from_rows(length - side for length in self.rows[:side])
        return side, right, YoungDiagram(self.rows[side:])

    def remove_rim_hook(self, row: int, column: int) -> YoungDiagram:
        """Return the diagram with the rim hook of a cell removed.

        Frame-Robinson-Thrall 1954, section 3: the result is a diagram
        with ``h_ij`` fewer cells — the same one obtained by removing the
        right hook and closing up.

        Raises:
            ValueError: If ``(row, column)`` is not a cell.
        """
        return YoungDiagram.from_cells(self.cells - self.rim_hook(row, column))

    def core(self, q: int) -> YoungDiagram:
        """Return the ``q``-core: what remains after removing all ``q``-hooks.

        Frame-Robinson-Thrall 1954, section 6: "The q-core of [lambda] is
        the diagram [alpha] that remains after all q-hooks have been
        removed". Here ``q`` is a positive integer, a hook length.

        Raises:
            ValueError: If ``q < 1``.
        """
        if q < 1:
            msg = f"a q-core needs a positive hook length q; got {q}"
            raise ValueError(msg)
        diagram = self
        while True:
            cell = next(
                (
                    (i, j)
                    for i, row in enumerate(diagram.hook_lengths, start=1)
                    for j, h in enumerate(row, start=1)
                    if h == q
                ),
                None,
            )
            if cell is None:
                return diagram
            diagram = diagram.remove_rim_hook(*cell)

    def dominance_meet(self, other: YoungDiagram) -> YoungDiagram:
        """Return the meet in the dominance lattice.

        Brylawski (*The lattice of integer partitions*, 1973), as
        attributed by Latapy-Phan (proof of their Proposition 1): the
        partial sums of the meet are the minima of the partial sums.

        Raises:
            ValueError: If the sizes differ.
        """
        self._require_same_size(other)
        sums = [
            min(mine, theirs)
            for mine, theirs in itertools.zip_longest(
                itertools.accumulate(self.rows),
                itertools.accumulate(other.rows),
                fillvalue=self.size,
            )
        ]
        return YoungDiagram.from_rows(b - a for a, b in itertools.pairwise([0, *sums]))

    def dominance_join(self, other: YoungDiagram) -> YoungDiagram:
        """Return the join in the dominance lattice, ``(lambda' ^ mu')'``.

        The page's "Join" block (Wikipedia, "Dominance order", verified
        there by computation for ``n <= 9``); it rests on conjugation
        being a lattice anti-automorphism (Brylawski 1973, abstract). The
        componentwise maximum of partial sums fails.

        Raises:
            ValueError: If the sizes differ.
        """
        self._require_same_size(other)
        return self.conjugate().dominance_meet(other.conjugate()).conjugate()

    def dominance_lower_covers(self) -> tuple[YoungDiagram, ...]:
        """Return the partitions this one covers in the dominance order.

        Brylawski 1973, as restated in Behrisch et al.
        (arXiv:2012.09926) Theorem 3: ``alpha`` covers ``beta`` iff
        ``alpha_j = beta_j + 1``, ``alpha_k = beta_k - 1`` for some
        ``j < k``, all other parts agree, and either ``k = j + 1`` or
        ``beta_j = beta_k``.
        """
        alpha = (*self.rows, 0)
        covered: list[YoungDiagram] = []
        for j, k in itertools.combinations(range(len(alpha)), 2):
            beta = list(alpha)
            beta[j] -= 1
            beta[k] += 1
            is_partition = all(a >= b for a, b in itertools.pairwise(beta))
            if is_partition and (k == j + 1 or beta[j] == beta[k]):
                covered.append(YoungDiagram.from_rows(beta))
        return tuple(covered)

    def to_le_diagram(self, k: int, n: int, ones: Iterable[Cell] = ()) -> LeDiagram:
        """Return the Le-diagram of type ``(k, n)`` on this shape.

        Postnikov 2006 fills the diagrams inside the ``k x (n - k)``
        rectangle with 0's and 1's to obtain Le-diagrams (the page's "Role
        in this graph"; Theorem 6.5). Inverted on shapes by
        :meth:`from_le_diagram`.

        Args:
            k: The number of rows of the box; the shape is padded to it.
            n: The ambient size.
            ones: The boxes filled with 1; all 0's when omitted.

        Raises:
            ValueError: If the diagram does not fit the box, or the
                filling is not a Le-diagram.
        """
        self._require_box(k, n)
        return LeDiagram(self.padded_rows(k), n, frozenset(ones))

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame: one row per cell, in reading order.

        Columns: ``row`` and ``column`` (1-indexed, English notation). The
        empty diagram is the frame with no rows. The encoding survives
        ``experiments.io.write_result`` (records-oriented JSON), where the
        empty diagram loses its columns.

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        cells = sorted(self.cells)
        return pd.DataFrame(
            {
                "row": pd.array([i for i, _ in cells], dtype="Int64"),
                "column": pd.array([j for _, j in cells], dtype="Int64"),
            }
        )

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> YoungDiagram:
        """Rebuild a diagram from :meth:`to_dataframe` output.

        Rows may come in any order.

        Raises:
            ValueError: If a column is missing, a cell is repeated, or the
                cells are not a Young diagram.
        """
        if df.empty and len(df.columns) == 0:
            # Records-oriented JSON of the empty diagram has no columns.
            return YoungDiagram()
        missing = [column for column in _COLUMNS if column not in df.columns]
        if missing:
            msg = f"Young-diagram frame is missing columns {missing!r}"
            raise ValueError(msg)
        cells = [(int(i), int(j)) for i, j in zip(df["row"], df["column"], strict=True)]
        if len(set(cells)) != len(cells):
            msg = f"Young-diagram frame repeats a cell; got {sorted(cells)!r}"
            raise ValueError(msg)
        return YoungDiagram.from_cells(cells)

    # -------------------------------------------------------- visualization
    def _cell_outline(
        self, cell: Cell, convention: DrawingConvention
    ) -> list[tuple[float, float]]:
        """Return the closed outline of a cell in plot coordinates."""
        i, j = cell
        french = [(j - 1, i - 1), (j, i - 1), (j, i), (j - 1, i), (j - 1, i - 1)]
        match convention:
            case DrawingConvention.ENGLISH:
                return [(x, -y) for x, y in french]
            case DrawingConvention.FRENCH:
                return [(float(x), float(y)) for x, y in french]
            case DrawingConvention.RUSSIAN:
                scale = math.sqrt(2) / 2
                return [((x - y) * scale, (x + y) * scale) for x, y in french]

    def _draw_cells(self, ax: Axes, convention: DrawingConvention, title: str) -> None:
        """Draw every cell's outline and frame the axes around them."""
        for cell in sorted(self.cells):
            outline = self._cell_outline(cell, convention)
            ax.plot(
                [x for x, _ in outline],
                [y for _, y in outline],
                color="0.3",
                linewidth=1,
            )
        ax.set_aspect("equal")
        ax.margins(0.1)
        ax.set_axis_off()
        ax.set_title(title)

    def plot_diagram(
        self,
        ax: Axes | None = None,
        *,
        convention: DrawingConvention = DrawingConvention.ENGLISH,
    ) -> Axes:
        """Draw the diagram of boxes onto ``ax``.

        Stanley EC1 section 1.7: the Young diagram replaces the dots of the
        Ferrers diagram "by juxtaposed squares"; the three drawing
        conventions are those of Adin-Roichman section 2.4.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.
            convention: English (default), French, or Russian.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        self._draw_cells(ax, convention, f"Young diagram ({convention.value})")
        return ax

    def plot_hook_lengths(self, ax: Axes | None = None) -> Axes:
        """Draw the hook graph ``H[lambda]`` onto ``ax``, in English notation.

        Frame-Robinson-Thrall 1954, section 1: every cell carries its hook
        length.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        self._draw_cells(ax, DrawingConvention.ENGLISH, "Hook graph")
        for i, row in enumerate(self.hook_lengths, start=1):
            for j, h in enumerate(row, start=1):
                ax.annotate(
                    str(h), (j - 0.5, 0.5 - i), ha="center", va="center", fontsize=10
                )
        return ax


# --------------------------------------------------------------------------- #
# Intervals of Young's lattice
# --------------------------------------------------------------------------- #
def mobius(lower: YoungDiagram, upper: YoungDiagram) -> int:
    r"""Return the Möbius function ``mu_Y(nu, lambda)`` of Young's lattice.

    Stanley EC1 Example 3.9.6, specialized to ``Y`` (the page's inference):
    ``(-1)^{|lambda / nu|}`` if no two cells of ``lambda / nu`` are
    componentwise comparable, and 0 otherwise.

    Raises:
        ValueError: If ``nu`` is not contained in ``lambda``.
    """
    skew = upper.skew_cells(lower)
    comparable = any(
        (a[0] <= b[0] and a[1] <= b[1]) or (b[0] <= a[0] and b[1] <= a[1])
        for a, b in itertools.combinations(skew, 2)
    )
    return 0 if comparable else (-1) ** len(skew)


def count_multichains(lower: YoungDiagram, upper: YoungDiagram, n: int) -> int:
    r"""Return ``Z(n + 1)``, the zeta polynomial of the interval ``[nu, lambda]``.

    Stanley EC1 Exercise 3.149 (Kreweras 1965; MacMahon for ``nu`` empty):
    ``Z(n + 1) = det[binom(lambda_i - nu_j + n, i - j + n)]`` over
    ``1 <= i, j <= m`` with ``lambda_{m+1} = 0``, a binomial with negative
    top being 0. ``Z(2)`` counts the diagrams in the interval and ``Z(3)``
    the pairs ``nu <= alpha <= beta <= lambda``.

    Raises:
        ValueError: If ``nu`` is not contained in ``lambda``, or ``n < 0``.
    """
    if n < 0 or not upper.contains(lower):
        msg = (
            f"a multichain count needs nu inside lambda and n >= 0 (Stanley "
            f"EC1 Exercise 3.149); got {lower.rows!r}, {upper.rows!r}, n = {n}"
        )
        raise ValueError(msg)
    m = upper.length
    inner = lower.padded_rows(m)
    matrix = [
        [
            Fraction(_binomial(length - inner[j - 1] + n, i - j + n))
            for j in range(1, m + 1)
        ]
        for i, length in enumerate(upper.rows, start=1)
    ]
    return int(det_q(matrix))


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
def empty_diagram() -> YoungDiagram:
    """Return the empty diagram, the partition of 0.

    Young Diagram page, canonical examples: no cells, ``f = 1``, the bottom
    of Young's lattice, covered by exactly one diagram ``(1)``.
    """
    return YoungDiagram()


def all_in_one_diagram() -> YoungDiagram:
    """Return ``(4, 3, 1)``, the page's all-in-one fixture.

    Adin-Roichman section 2.4 and the picture after Definition 5.2 (hook
    lengths ``6 4 3 1 / 4 2 1 / 1``); Stanley EC1 Figure 1.15 (boundary
    word ``12121121`` in a ``3 x 5`` box).
    """
    return YoungDiagram((4, 3, 1))


def frt_example_1_diagram() -> YoungDiagram:
    """Return ``(6, 4, 2)``, the original hook-graph example.

    Frame-Robinson-Thrall 1954, Example 1: ``f = 12!/H = 2673``.
    """
    return YoungDiagram((6, 4, 2))


def frt_example_2_diagram() -> YoungDiagram:
    """Return ``(7, 6, 5, 3)``, the rim-hook removal example.

    Frame-Robinson-Thrall 1954, Example 2: removing the 6-hook at cell
    ``(2, 3)`` leaves ``(7, 4, 2, 2)``.
    """
    return YoungDiagram((7, 6, 5, 3))


def frt_example_4_diagram() -> YoungDiagram:
    """Return ``(9, 7, 4, 3, 2)``, whose 3-core is ``(3, 1)``.

    Frame-Robinson-Thrall 1954, Example 4.
    """
    return YoungDiagram((9, 7, 4, 3, 2))


def postnikov_figure_2_1_diagram() -> YoungDiagram:
    """Return ``(4, 4, 2, 1)``, the box fixture inside ``6^4``.

    Postnikov Figure 2.1 and section 2.1: ``I(lambda) = {3, 4, 7, 9}`` in
    ``[10]``.
    """
    return YoungDiagram((4, 4, 2, 1))


def stanley_figure_1_16_diagram() -> YoungDiagram:
    """Return the self-conjugate ``(5, 4, 4, 3, 1)``.

    Stanley EC1 Figure 1.16: diagonal hooks of sizes ``9, 5, 3``.
    """
    return YoungDiagram((5, 4, 4, 3, 1))


def stanley_figure_1_18_diagram() -> YoungDiagram:
    """Return ``(7, 5, 3, 3, 2)``, the Durfee-dissection example.

    Stanley EC1 Figure 1.18: Durfee side 3, ``mu = (4, 2)``,
    ``nu = (3, 2)``.
    """
    return YoungDiagram((7, 5, 3, 3, 2))


def stanley_figure_1_33_diagram() -> YoungDiagram:
    """Return ``(7, 7, 4, 2, 1)``, the Frobenius-notation example.

    Stanley EC1 Solution 1.70(a), Figure 1.33: ``(6 5 1 | 4 2 0)``.
    """
    return YoungDiagram((7, 7, 4, 2, 1))


def adin_roichman_example_9_7_diagram() -> YoungDiagram:
    """Return ``(3, 1)``, the boundary-sequence example.

    Adin-Roichman Example 9.7: boundary sequence ``(1, 0, 1, 1, 0)``.
    """
    return YoungDiagram((3, 1))


def staircase(n: int) -> YoungDiagram:
    """Return the staircase ``delta_n = (n, n - 1, ..., 1)``.

    Adin-Roichman Definition 4.5.

    Raises:
        ValueError: If ``n < 0``.
    """
    if n < 0:
        msg = f"a staircase needs n >= 0; got {n}"
        raise ValueError(msg)
    return YoungDiagram(tuple(range(n, 0, -1)))


def rectangle(k: int, width: int) -> YoungDiagram:
    """Return the rectangle ``width^k`` with ``k`` rows.

    Postnikov section 2.1: ``(n - k)^k`` is the full ``k x (n - k)`` box.

    Raises:
        ValueError: If a side is negative.
    """
    if k < 0 or width < 0:
        msg = f"a rectangle needs non-negative sides; got {k} x {width}"
        raise ValueError(msg)
    return YoungDiagram.from_rows((width,) * k)


def hook_shape(n: int, k: int) -> YoungDiagram:
    """Return the hook shape ``(n - k, 1^k)`` with ``n`` cells.

    Adin-Roichman section 3.1: "the union of one row and one column".

    Raises:
        ValueError: Unless ``0 <= k <= n - 1``.
    """
    if not 0 <= k <= n - 1:
        msg = f"a hook shape (n - k, 1^k) needs 0 <= k <= n - 1; got n = {n}, k = {k}"
        raise ValueError(msg)
    return YoungDiagram((n - k, *([1] * k)))


def non_down_closed_cells() -> frozenset[Cell]:
    """Return ``{(1, 1), (2, 2)}``, a diagram that is not a Young diagram.

    Young Diagram page, non-examples: a diagram in Adin-Roichman's sense
    (any finite subset of ``Z^2``) that violates down-closure.
    :meth:`YoungDiagram.from_cells` rejects it.
    """
    return frozenset({(1, 1), (2, 2)})


def skew_non_example_cells() -> frozenset[Cell]:
    """Return the skew shape ``(6, 4, 3, 1) / (4, 2, 1)``.

    Young Diagram page, non-examples (Adin-Roichman, example for
    Definition 2.12): order-convex but not an order ideal, so *skew* is
    strictly weaker than *ordinary*.
    """
    return YoungDiagram((6, 4, 3, 1)).skew_cells(YoungDiagram((4, 2, 1)))


def non_decreasing_rows() -> tuple[int, ...]:
    """Return the row sequence ``(2, 3)``, which violates weak decrease.

    Young Diagram page, non-examples. :meth:`YoungDiagram.from_rows`
    rejects it.
    """
    return (2, 3)


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #
def _partition_rows(n: int, largest: int) -> Iterator[tuple[int, ...]]:
    """Yield the partitions of ``n`` with parts at most ``largest``."""
    if n == 0:
        yield ()
        return
    for first in range(min(n, largest), 0, -1):
        for rest in _partition_rows(n - first, first):
            yield (first, *rest)


def partitions(n: int) -> Iterator[YoungDiagram]:
    """Yield the ``p(n)`` diagrams with ``n`` cells, lazily.

    In reverse lexicographic order, from ``(n)`` to ``(1^n)``; ``p(0) = 1``
    (Stanley EC1 p. 65). The iterator is single-use.

    Raises:
        ValueError: If ``n < 0``.
    """
    if n < 0:
        msg = f"partitions are of a non-negative integer; got {n}"
        raise ValueError(msg)
    return (YoungDiagram(rows) for rows in _partition_rows(n, n))


def diagrams_in_box(k: int, n: int) -> Iterator[YoungDiagram]:
    """Yield the ``binom(n, k)`` diagrams inside the ``k x (n - k)`` box, lazily.

    Stanley EC1 Proposition 1.7.3 at ``q = 1``; Postnikov section 2.1. The
    iterator is single-use.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    _check_box(k, n)
    lengths = itertools.combinations_with_replacement(range(n - k, -1, -1), k)
    return (YoungDiagram.from_rows(rows) for rows in lengths)


def box_size_distribution(k: int, n: int) -> tuple[int, ...]:
    """Return the numbers ``p(n - k, k, m)`` of box diagrams of each size ``m``.

    Stanley EC1 Proposition 1.7.3: these are the coefficients of the
    Gaussian binomial coefficient ``binom(n, k)_q``, index = power. Counted
    here by enumeration.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    counts = [0] * (k * (n - k) + 1) if 0 <= k <= n else []
    for diagram in diagrams_in_box(k, n):
        counts[diagram.size] += 1
    return tuple(counts)
