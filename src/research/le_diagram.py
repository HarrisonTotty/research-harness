r"""Le-diagrams: 0/1 fillings of Young diagrams avoiding the mirrored-L pattern.

A Le-diagram is a filling of a Young diagram with 0's and 1's in which no 0
has simultaneously a 1 above it in its column and a 1 to its left in its row
(Postnikov, *Total positivity, Grassmannians, and networks*,
arXiv:math/0609764, 2006, Definition 6.1). The Le-diagrams whose shape fits in
the ``k x (n - k)`` rectangle index the cells of the totally nonnegative
Grassmannian, the cell's dimension being the number of 1's (Postnikov
Theorem 6.5). :class:`LeDiagram` stores the shape, the ambient ``n``, and the
set of 1-boxes, and derives the rest of the Le-Diagram page as views — the
boundary path and ``I(lambda)``, the pipe dream to a decorated permutation
and its inverse, the transpose, the Gamma-network and the Le-graph — next to
the page's four equivalent statements of the Le-property, its enumerations,
and its canonical examples.

Conventions, all following the page:

* **Coordinates.** Boxes are matrix-style ``(row, column)``, 1-indexed, in
  English notation (Postnikov, Williams, Ardila-Rincon-Williams).
  Lam-Williams draw the same diagrams in French notation;
  :meth:`LeDiagram.from_french_filling` reflects them.
* **Type.** The shape has exactly ``k`` rows, empty rows included, and sits
  in the ``k x (n - k)`` rectangle, so empty rows and empty columns — the
  two kinds of fixed point — are both recorded.
* **Boundary labels.** The ``n`` steps of the south-east border path, from
  the north-east to the south-west corner of the rectangle, carry ``1..n``;
  the vertical steps label the rows and the horizontal steps the columns.
* **Direction and color.** :meth:`LeDiagram.to_decorated_permutation` is the
  Ardila-Rincon-Williams pipe dream, which is the direction
  :mod:`research.positroid` and :mod:`research.grassmann_necklace` exchange;
  Postnikov's decorated trip permutation of the Le-graph is its
  :meth:`~research.decorated_permutation.DecoratedPermutation.inverse`. An
  all-zero row gives the fixed point Ardila-Rincon-Williams overline, this
  library's clockwise (counted) color.
"""

import enum
import functools
import itertools
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Literal, override

import pandas as pd

from research._plot import ensure_axes
from research.boundary_measurement import PlanarNetwork
from research.decorated_permutation import DecoratedPermutation
from research.grassmann_necklace import GrassmannNecklace
from research.plabic_graph import PlabicGraph
from research.positroid import Positroid

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "LeCondition",
    "LeDiagram",
    "arw_section_4_3_diagram",
    "empty_filling",
    "enumerate_le_diagrams",
    "enumerate_le_diagrams_of_shape",
    "enumerate_permutation_tableaux",
    "enumerate_shapes",
    "forbidden_pattern_fillings",
    "rank_generating_polynomial",
    "satisfies_le_condition",
    "top_cell_diagram",
]

type Box = tuple[int, int]
type Step = Literal["vertical", "horizontal"]

_VERTICAL: Step = "vertical"
_HORIZONTAL: Step = "horizontal"
_STEP = "step"
_BOX = "box"
_COLUMNS = ("kind", "label", "direction", "row", "column", "value")


class LeCondition(enum.Enum):
    """The page's four equivalent statements of the Le-property.

    * ``DEFINITION_6_1`` — Postnikov Definition 6.1: for boxes ``(i', j)``,
      ``(i', j')``, ``(i, j')`` with ``i < i'`` and ``j < j'`` filled with
      ``a``, ``b``, ``c``, if ``a, c != 0`` then ``b != 0``.
    * ``FORBIDDEN_PATTERN`` — Williams section 2; Ardila-Rincon-Williams
      Definition 4.7: no 0 has a 1 above it in its column and a 1 to its
      left in its row.
    * ``HOOK`` — Postnikov section 6: draw a hook (right and down) from
      every 1-box; each box where a horizontal and a vertical hook line
      intersect holds a 1.
    * ``BLOCKED_ZERO`` — Postnikov section 6: every entry to the left of a
      blocked 0 (one with a 1 above it in its column) is 0.
    """

    DEFINITION_6_1 = "definition-6.1"
    FORBIDDEN_PATTERN = "forbidden-pattern"
    HOOK = "hook"
    BLOCKED_ZERO = "blocked-zero"


def _definition_6_1_violation(
    rows: Sequence[Sequence[int]],
) -> tuple[Box, Box, Box] | None:
    """Return boxes ``(a, b, c)`` violating Postnikov Definition 6.1, if any."""
    for i2, row in enumerate(rows, start=1):
        for j2, b in enumerate(row, start=1):
            if b:
                continue
            for i1, j1 in itertools.product(range(1, i2), range(1, j2)):
                if row[j1 - 1] and rows[i1 - 1][j2 - 1]:
                    return ((i2, j1), (i2, j2), (i1, j2))
    return None


def _has_forbidden_pattern(rows: Sequence[Sequence[int]]) -> bool:
    """Return whether some 0 sees a 1 above it and a 1 to its left."""
    return any(
        not value
        and any(row[:j])
        and any(above[j] for above in rows[:i] if len(above) > j)
        for i, row in enumerate(rows)
        for j, value in enumerate(row)
    )


def _violates_hook_condition(rows: Sequence[Sequence[int]]) -> bool:
    """Return whether two hook lines intersect at a box without a dot."""
    horizontal: set[Box] = set()
    vertical: set[Box] = set()
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            if not value:
                continue
            horizontal.update((i, right) for right in range(j + 1, len(row)))
            vertical.update(
                (below, j) for below in range(i + 1, len(rows)) if len(rows[below]) > j
            )
    return any(not rows[i][j] for i, j in horizontal & vertical)


def _violates_blocked_zero_condition(rows: Sequence[Sequence[int]]) -> bool:
    """Return whether some blocked 0 has a nonzero entry to its left."""
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            blocked = not value and any(
                above[j] for above in rows[:i] if len(above) > j
            )
            if blocked and any(row[:j]):
                return True
    return False


def satisfies_le_condition(
    rows: Sequence[Sequence[int]], condition: LeCondition
) -> bool:
    """Return whether a 0/1 filling satisfies one statement of the Le-property.

    The four statements (:class:`LeCondition`) are equivalent on every
    filling of a Young diagram (Le-Diagram page, Definition); each is
    evaluated here literally, so their agreement is a checkable fact rather
    than an identity of code.

    Args:
        rows: The filling, rows top to bottom in English notation, of weakly
            decreasing lengths; nonzero entries count as 1.
        condition: Which statement to evaluate.

    Raises:
        ValueError: If the row lengths are not weakly decreasing.
    """
    lengths = [len(row) for row in rows]
    if any(a < b for a, b in itertools.pairwise(lengths)):
        msg = f"shape violated: row lengths {lengths} must be weakly decreasing"
        raise ValueError(msg)
    match condition:
        case LeCondition.DEFINITION_6_1:
            return _definition_6_1_violation(rows) is None
        case LeCondition.FORBIDDEN_PATTERN:
            return not _has_forbidden_pattern(rows)
        case LeCondition.HOOK:
            return not _violates_hook_condition(rows)
        case LeCondition.BLOCKED_ZERO:
            return not _violates_blocked_zero_condition(rows)


def _border_labels(
    shape: Sequence[int], n: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return the border-path labels of the rows and of the columns.

    Postnikov section 2.1: the ``n`` steps of the boundary path from the
    upper-right to the lower-left corner are labelled ``1..n``. Row ``i``
    owns the vertical step ``i + (n - k) - shape[i - 1]``, and the remaining
    labels go to the columns from right to left.
    """
    k = len(shape)
    width = n - k
    rows = tuple(i + width - length for i, length in enumerate(shape, start=1))
    rest = sorted(set(range(1, n + 1)) - set(rows))
    return rows, tuple(reversed(rest))


# --------------------------------------------------------------------------- #
# The structure
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class LeDiagram:
    """A Le-diagram of type ``(k, n)``: a shape in a box with its 1-boxes.

    ``shape`` lists the ``k`` row lengths, weakly decreasing with empty rows
    allowed, inside the ``k x (n - k)`` rectangle; ``ones`` holds the boxes
    filled with 1, as 1-indexed ``(row, column)`` pairs, and every other box
    of the shape holds 0 (Postnikov Definition 6.1, with the type ``(d, n)``
    of Ardila-Rincon-Williams Definition 4.7).

    Like :class:`~research.decorated_permutation.DecoratedPermutation`, the
    definition is cheap enough to be checked by ``__post_init__``, so direct
    construction is validated too; the ``from_*`` classmethods convert the
    page's other presentations.

    Raises:
        ValueError: Naming the violated part of the definition — the shape,
            the type, the filling, or the Le-property.
    """

    shape: tuple[int, ...]
    n: int
    ones: frozenset[Box] = frozenset()

    def __post_init__(self) -> None:
        """Validate the shape, the type, the filling, and the Le-property."""
        shape, n = self.shape, self.n
        decreasing = all(a >= b for a, b in itertools.pairwise(shape))
        if not decreasing or any(length < 0 for length in shape):
            msg = (
                f"shape violated: row lengths must be weakly decreasing and "
                f"non-negative (a partition); got {shape!r}"
            )
            raise ValueError(msg)
        k = len(shape)
        if k > n or (shape and shape[0] > n - k):
            msg = (
                f"type violated: shape {shape!r} must fit in the k x (n - k) "
                f"rectangle with k = {k} rows and n = {n} "
                f"(Postnikov Definition 6.1)"
            )
            raise ValueError(msg)
        stray = sorted(
            box
            for box in self.ones
            if not (1 <= box[0] <= k and 1 <= box[1] <= shape[box[0] - 1])
        )
        if stray:
            msg = (
                f"filling violated: the 1-boxes {stray!r} lie outside the "
                f"shape {shape!r}"
            )
            raise ValueError(msg)
        violation = _definition_6_1_violation(self.filling)
        if violation is not None:
            a, b, c = violation
            msg = (
                f"Le-property violated (Postnikov Definition 6.1): boxes {a!r} "
                f"and {c!r} hold 1 but {b!r} holds 0"
            )
            raise ValueError(msg)

    @override
    def __repr__(self) -> str:
        """Return a compact form: the type and the rows in 0/+ notation."""
        return f"LeDiagram(type={self.diagram_type!r}, rows={self.to_plus_filling()!r})"

    @override
    def __str__(self) -> str:
        """Return the rows in 0/+ notation, one per line; ``.`` is an empty row."""
        return "\n".join(row or "." for row in self.to_plus_filling())

    # -------------------------------------------------------- constructors
    @classmethod
    def from_filling(cls, rows: Iterable[Iterable[int]], n: int) -> LeDiagram:
        """Build from the rows of a 0/1 filling, top to bottom.

        Postnikov Definition 6.1, English notation. The number of rows is
        ``k``, so empty rows must be given as empty sequences.

        Args:
            rows: The entries of each row, left to right.
            n: The ambient size; the shape must fit in ``k x (n - k)``.

        Raises:
            ValueError: If an entry is not 0 or 1, or the definition fails.
        """
        filling = [list(row) for row in rows]
        bad = sorted({value for row in filling for value in row} - {0, 1})
        if bad:
            msg = (
                f"filling violated: a Le-diagram is filled with 0's and 1's "
                f"(Postnikov Definition 6.1); got entries {bad!r}"
            )
            raise ValueError(msg)
        ones = frozenset(
            (i, j)
            for i, row in enumerate(filling, start=1)
            for j, value in enumerate(row, start=1)
            if value
        )
        return cls(tuple(len(row) for row in filling), n, ones)

    @classmethod
    def from_plus_filling(cls, rows: Iterable[str], n: int) -> LeDiagram:
        """Build from rows written with ``0`` and ``+``.

        Ardila-Rincon-Williams Definition 4.7 fill with 0's and +'s instead
        of 0's and 1's — "pure notation, same condition" (Le-Diagram page).

        Args:
            rows: One string per row, top to bottom, e.g. ``"0+0+0"``.
            n: The ambient size; the shape must fit in ``k x (n - k)``.

        Raises:
            ValueError: If a character is not ``0`` or ``+``, or the
                definition fails.
        """
        strings = list(rows)
        bad = sorted(set("".join(strings)) - {"0", "+"})
        if bad:
            msg = (
                f"filling violated: rows are written with '0' and '+' "
                f"(Ardila-Rincon-Williams Definition 4.7); got {bad!r}"
            )
            raise ValueError(msg)
        return cls.from_filling(([int(c == "+") for c in row] for row in strings), n)

    @classmethod
    def from_french_filling(cls, rows: Iterable[Iterable[int]], n: int) -> LeDiagram:
        """Build from a filling drawn in French notation.

        Lam-Williams (*Total positivity for cominuscule Grassmannians*,
        2008, Theorem 5.1) draw the type-A diagrams in French notation,
        where the condition reads "no 0 which has a + below it and a + to
        its left": "Postnikov's Le-diagrams are obtained from ours by
        reflecting in a horizontal axis."

        Args:
            rows: The rows as drawn, top to bottom — shortest first.
            n: The ambient size.

        Raises:
            ValueError: If the reflected filling is not a Le-diagram.
        """
        return cls.from_filling(reversed([list(row) for row in rows]), n)

    @classmethod
    def from_decorated_permutation(cls, decorated: DecoratedPermutation) -> LeDiagram:
        """Build the Le-diagram whose pipe dream is ``decorated``.

        The inverse of :meth:`to_decorated_permutation`, which
        Ardila-Rincon-Williams Lemma 4.8 states is a bijection between
        Le-diagrams of type ``(d, n)`` and decorated permutations on ``n``
        letters with ``d`` weak excedances. The shape is read off the weak
        excedances, which label the vertical border steps: ``I(lambda)`` is
        the anti-exceedance set of Postnikov's inverse permutation
        (Corollary 20.1). The boxes are then filled in
        reading order, following the two pipes that meet in each box: the
        box is a crossing (0) exactly when the pipe arriving from the west
        has the smaller target. The lemma supplies the bijection, not this
        filling rule, so the result is confirmed by running the pipe dream
        forwards.

        Args:
            decorated: A decorated permutation in the library's stored
                (Ardila-Rincon-Williams) direction.

        Raises:
            ValueError: If no Le-diagram is recovered — impossible for a
                valid decorated permutation, by the lemma.
        """
        n = decorated.size
        vertical = sorted(decorated.weak_excedances)
        width = n - len(vertical)
        shape = tuple(i + width - label for i, label in enumerate(vertical, start=1))
        row_labels, column_labels = _border_labels(shape, n)
        north = list(column_labels)
        ones: set[Box] = set()
        for i, length in enumerate(shape, start=1):
            west = row_labels[i - 1]
            for j in range(1, length + 1):
                if decorated.targets[west - 1] > decorated.targets[north[j - 1] - 1]:
                    ones.add((i, j))
                    west, north[j - 1] = north[j - 1], west
        diagram = cls(shape, n, frozenset(ones))
        if diagram.to_decorated_permutation() != decorated:
            msg = (
                f"no Le-diagram was recovered for {decorated!r} "
                f"(Ardila-Rincon-Williams Lemma 4.8)"
            )
            raise ValueError(msg)
        return diagram

    @classmethod
    def from_grassmann_necklace[T: Hashable](
        cls, necklace: GrassmannNecklace[T]
    ) -> LeDiagram:
        """Build the Le-diagram of a Grassmann necklace.

        Through the necklace's decorated permutation (Postnikov Lemma 16.2)
        and :meth:`from_decorated_permutation`; boundary label ``i`` is the
        necklace's ``i``-th ground-set element.
        """
        return cls.from_decorated_permutation(necklace.to_decorated_permutation())

    @classmethod
    def from_positroid[T: Hashable](cls, positroid: Positroid[T]) -> LeDiagram:
        """Build the Le-diagram indexing a positroid's cell.

        Postnikov Theorem 6.5 (Le-diagrams biject with the cells), through
        the positroid's decorated permutation; boundary label ``i`` is the
        ``i``-th ground-set element in the cyclic order.
        """
        return cls.from_decorated_permutation(positroid.to_decorated_permutation())

    # ------------------------------------------------- computed properties
    @property
    def k(self) -> int:
        """The number of rows of the ambient rectangle, empty rows included."""
        return len(self.shape)

    @property
    def diagram_type(self) -> tuple[int, int]:
        """The type ``(k, n)``: the shape lies in the ``k x (n - k)`` rectangle.

        Ardila-Rincon-Williams Definition 4.7 ("type ``(d, n)``"); the set
        of these diagrams is Postnikov's ``Le_kn``.
        """
        return (len(self.shape), self.n)

    @property
    def rank(self) -> int:
        """The number ``|D|`` of 1's — the dimension of the indexed cell.

        Williams (*Enumeration of totally positive Grassmann cells*, 2005):
        "the rank of ``(lambda, D)_{k,n}``"; Postnikov Theorem 6.5:
        ``dim S_M^tnn = |D|``.
        """
        return len(self.ones)

    @functools.cached_property
    def filling(self) -> tuple[tuple[int, ...], ...]:
        """The rows of the filling as 0/1 tuples, top to bottom."""
        return tuple(
            tuple(int((i, j) in self.ones) for j in range(1, length + 1))
            for i, length in enumerate(self.shape, start=1)
        )

    @property
    def row_labels(self) -> tuple[int, ...]:
        """The border labels of rows ``1..k`` — the vertical steps, ascending."""
        return _border_labels(self.shape, self.n)[0]

    @property
    def column_labels(self) -> tuple[int, ...]:
        """The border labels of columns ``1..n-k`` — descending left to right."""
        return _border_labels(self.shape, self.n)[1]

    @property
    def boundary_path(self) -> tuple[Step, ...]:
        """The boundary lattice path, the step labelled ``i`` at index ``i - 1``.

        Postnikov section 2.1: the path from the upper-right to the
        lower-left corner of the rectangle along the border of the shape.
        The shape is recovered from it.
        """
        vertical = set(self.row_labels)
        return tuple(
            _VERTICAL if label in vertical else _HORIZONTAL
            for label in range(1, self.n + 1)
        )

    @property
    def vertical_steps(self) -> frozenset[int]:
        """The set ``I(lambda)`` of labels of the ``k`` vertical steps.

        Postnikov section 2.1. It is the source set of the Gamma-network,
        the anti-exceedance set of Postnikov's decorated permutation (his
        Corollary 20.1), and hence the set of weak excedances of its
        inverse, :meth:`to_decorated_permutation`.
        """
        return frozenset(self.row_labels)

    @property
    def is_permutation_tableau(self) -> bool:
        """Whether every column of the rectangle contains at least one 1.

        Steingrimsson-Williams (*Permutation tableaux and permutation
        patterns*, 2007, section 1): a permutation tableau is a Le-diagram
        in a ``k x (n - k)`` rectangle with that extra condition, which
        forces the shape to have exactly ``n - k`` columns while rows may
        be empty.
        """
        return {j for _, j in self.ones} == set(range(1, self.n - self.k + 1))

    @property
    def weight(self) -> int:
        """The weight of a permutation tableau: its 1's minus its columns.

        Corteel-Williams (*Tableaux combinatorics for the asymmetric
        exclusion process*, 2007, Theorems 3.1 and 3.4), where it matches
        the number of crossings of the corresponding permutation.

        Raises:
            ValueError: If the diagram is not a permutation tableau — the
                page defines the weight only there.
        """
        if not self.is_permutation_tableau:
            msg = (
                f"the weight is defined for permutation tableaux only "
                f"(Corteel-Williams Theorem 3.1); {self!r} has a column "
                f"without a 1"
            )
            raise ValueError(msg)
        return self.rank - (self.n - self.k)

    def satisfies(self, condition: LeCondition) -> bool:
        """Return whether the filling meets one statement of the Le-property.

        Always true for a validated diagram — the four statements are
        equivalent (Le-Diagram page, Definition).
        """
        return satisfies_le_condition(self.filling, condition)

    # ------------------------------------------------------ transformations
    def transpose(self) -> LeDiagram:
        """Return the reflection over the main diagonal, of type ``(n - k, n)``.

        Williams 2005, section 4: reflecting ``(lambda, D)_{k,n}`` gives a
        Le-diagram ``(lambda', D')_{n-k,n}`` of the same rank, whence
        ``A_{k,n}(q) = A_{n-k,n}(q)``. An involution.
        """
        conjugate = tuple(
            sum(1 for length in self.shape if length >= j)
            for j in range(1, self.n - self.k + 1)
        )
        return LeDiagram(conjugate, self.n, frozenset((j, i) for i, j in self.ones))

    def to_plus_filling(self) -> tuple[str, ...]:
        """Return the rows in the 0/+ notation of Ardila-Rincon-Williams."""
        return tuple("".join("+" if v else "0" for v in row) for row in self.filling)

    def to_french_filling(self) -> tuple[tuple[int, ...], ...]:
        """Return the rows as Lam-Williams draw them: reflected top to bottom."""
        return tuple(reversed(self.filling))

    @functools.cached_property
    def _pipe_dream(self) -> tuple[dict[int, tuple[Box, ...]], dict[int, int]]:
        """Trace every pipe: the boxes it visits and the label where it exits."""
        north = list(self.column_labels)
        paths: dict[int, list[Box]] = {label: [] for label in range(1, self.n + 1)}
        exits: dict[int, int] = {}
        for i, row in enumerate(self.filling, start=1):
            west = self.row_labels[i - 1]
            for j, value in enumerate(row, start=1):
                paths[west].append((i, j))
                paths[north[j - 1]].append((i, j))
                if value:
                    west, north[j - 1] = north[j - 1], west
            exits[west] = self.row_labels[i - 1]
        exits.update(zip(north, self.column_labels, strict=True))
        return {label: tuple(path) for label, path in paths.items()}, exits

    def pipes(self) -> dict[int, tuple[Box, ...]]:
        """Return the boxes each pipe of the pipe dream passes through.

        Ardila-Rincon-Williams Lemma 4.8, steps (1)-(4): every 1 becomes an
        elbow joint and every 0 a crossing; the border labels are copied to
        the opposite (north and west) border, and the pipe entering at
        label ``i`` travels south-east. Traveling that way, an elbow sends
        the pipe from the west down and the pipe from the north right.

        Returns:
            For each entry label, the boxes visited in order — empty for a
            pipe that never enters the shape.
        """
        return dict(self._pipe_dream[0])

    def to_decorated_permutation(self) -> DecoratedPermutation:
        """Return the decorated permutation of the pipe dream.

        Ardila-Rincon-Williams Lemma 4.8: the pipe entering at ``i`` and
        exiting at ``j`` defines ``pi(i) = j``; a fixed point on two
        horizontal edges is underlined and one on two vertical edges
        overlined. The overline is this library's clockwise (counted)
        color. The map is a bijection from Le-diagrams of type ``(d, n)``
        to decorated permutations with ``d`` weak excedances, in the
        direction :mod:`research.positroid` uses; Postnikov's
        ``pi^:(G_D)`` (Corollary 20.1) is the inverse permutation.
        """
        _, exits = self._pipe_dream
        targets = tuple(exits[label] for label in range(1, self.n + 1))
        overlined = frozenset(
            label for label in self.row_labels if exits[label] == label
        )
        return DecoratedPermutation(targets, overlined)

    def to_grassmann_necklace(self) -> GrassmannNecklace[int]:
        """Return the Grassmann necklace on ``1..n`` of the indexed cell.

        Through :meth:`to_decorated_permutation` (Postnikov Lemma 16.2);
        its first entry is ``I(lambda)``.
        """
        return GrassmannNecklace.from_decorated_permutation(
            range(1, self.n + 1), self.to_decorated_permutation()
        )

    def to_positroid(self) -> Positroid[int]:
        """Return the positroid on ``1..n`` of the indexed cell.

        Postnikov Theorem 6.5 and Corollary 20.1, through
        :meth:`to_decorated_permutation`; agrees with the matroid of every
        matrix that :meth:`to_planar_network` produces from a Le-tableau.
        """
        return Positroid.from_decorated_permutation(
            range(1, self.n + 1), self.to_decorated_permutation()
        )

    def to_planar_network(
        self, weights: Mapping[Box, Fraction | int] | None = None
    ) -> PlanarNetwork[int | Box]:
        """Return the Gamma-network of a Le-tableau on this diagram.

        Postnikov Definition 6.3 and section 6: one internal vertex per
        1-box with its hook, horizontal edges oriented left and vertical
        edges down, the vertical border steps ``I(lambda)`` as sources; the
        horizontal edge into box ``(i, j)`` carries ``T(i, j)`` and vertical
        edges weight 1. Its boundary measurements parameterize the cell
        (Theorem 6.5). Weights are exact rationals, never floats.

        Args:
            weights: The Le-tableau ``T``, positive exactly on the 1-boxes;
                boxes left out default to 1.

        Raises:
            ValueError: If the tableau is not positive exactly on the
                1-boxes, or ``n = 0`` (no boundary to embed).
        """
        return PlanarNetwork.from_le_diagram(self.filling, self.n, weights)

    def to_plabic_graph(self) -> PlabicGraph[Hashable]:
        """Return the Le-graph ``G_D``, a reduced plabic graph.

        Postnikov section 20: the Gamma-graph with every 4-valent vertex
        split into two trivalent ones. Its decorated trip permutation is
        the inverse of :meth:`to_decorated_permutation` (Corollary 20.1).

        Raises:
            ValueError: If ``n = 0`` (no boundary to embed).
        """
        return PlabicGraph.from_le_diagram(self.filling, self.n)

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame: one row per border step, then per box.

        Columns: ``kind`` (``"step"``/``"box"``); step rows carry ``label``
        (``1..n``, ascending) and ``direction`` (``"vertical"`` or
        ``"horizontal"``), which together encode ``n``, ``k`` and the shape
        even when it has no boxes; box rows carry ``row`` and ``column``
        (1-indexed, English notation, reading order) and ``value`` (0 or
        1). Unused cells are null. The encoding survives
        ``experiments.io.write_result`` (records-oriented JSON); the type
        ``(0, 0)`` is the frame with no rows and no columns.

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        records: list[dict[str, object]] = [
            {
                "kind": _STEP,
                "label": label,
                "direction": direction,
                "row": None,
                "column": None,
                "value": None,
            }
            for label, direction in enumerate(self.boundary_path, start=1)
        ]
        records.extend(
            {
                "kind": _BOX,
                "label": None,
                "direction": None,
                "row": i,
                "column": j,
                "value": value,
            }
            for i, row in enumerate(self.filling, start=1)
            for j, value in enumerate(row, start=1)
        )
        if not records:
            return pd.DataFrame()
        frame = pd.DataFrame.from_records(records, columns=list(_COLUMNS))
        for column in ("label", "row", "column", "value"):
            frame[column] = frame[column].astype("Int64")
        return frame

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> LeDiagram:
        """Rebuild a Le-diagram from :meth:`to_dataframe` output.

        Rows may come in any order.

        Raises:
            ValueError: If a column is missing, the step labels are not
                exactly ``1..n``, a direction or kind is unknown, the boxes
                are not exactly those of the shape the steps encode, or the
                decoded filling is not a Le-diagram.
        """
        if df.empty and len(df.columns) == 0:
            # Records-oriented JSON of the type-(0, 0) diagram has no columns.
            return LeDiagram((), 0)
        missing = [column for column in _COLUMNS if column not in df.columns]
        if missing:
            msg = f"Le-diagram frame is missing columns {missing!r}"
            raise ValueError(msg)
        kinds = set(df["kind"])
        if not kinds <= {_STEP, _BOX}:
            msg = f"kind must be {_STEP!r} or {_BOX!r}; got {sorted(kinds)!r}"
            raise ValueError(msg)
        steps = df[df["kind"] == _STEP]
        directions = {
            int(label): str(direction)
            for label, direction in zip(steps["label"], steps["direction"], strict=True)
        }
        n = len(steps)
        if sorted(directions) != list(range(1, n + 1)):
            msg = (
                f"Le-diagram frame needs each step label 1..n exactly once; "
                f"got {sorted(int(label) for label in steps['label'])!r}"
            )
            raise ValueError(msg)
        unknown = sorted(set(directions.values()) - {_VERTICAL, _HORIZONTAL})
        if unknown:
            msg = f"direction must be {_VERTICAL!r} or {_HORIZONTAL!r}; got {unknown!r}"
            raise ValueError(msg)
        vertical = [
            label for label in range(1, n + 1) if directions[label] == _VERTICAL
        ]
        width = n - len(vertical)
        shape = tuple(i + width - label for i, label in enumerate(vertical, start=1))
        boxes = df[df["kind"] == _BOX]
        values = {
            (int(i), int(j)): int(value)
            for i, j, value in zip(
                boxes["row"], boxes["column"], boxes["value"], strict=True
            )
        }
        expected = {
            (i, j)
            for i, length in enumerate(shape, start=1)
            for j in range(1, length + 1)
        }
        if set(values) != expected or len(boxes) != len(expected):
            msg = (
                f"Le-diagram frame needs exactly one box row per box of the "
                f"shape {shape!r}; got {sorted(values)!r}"
            )
            raise ValueError(msg)
        return LeDiagram.from_filling(
            (
                [values[i, j] for j in range(1, length + 1)]
                for i, length in enumerate(shape, start=1)
            ),
            n,
        )

    # -------------------------------------------------------- visualization
    def _draw_frame(self, ax: Axes, title: str) -> None:
        """Draw the boxes, the ambient rectangle, and the border labels."""
        k, width = self.k, self.n - self.k
        ax.plot(
            [0, width, width, 0, 0],
            [0, 0, -k, -k, 0],
            color="0.8",
            linewidth=1,
            linestyle=":",
            zorder=0,
        )
        for i, length in enumerate(self.shape, start=1):
            for j in range(1, length + 1):
                ax.plot(
                    [j - 1, j, j, j - 1, j - 1],
                    [1 - i, 1 - i, -i, -i, 1 - i],
                    color="0.3",
                    linewidth=1,
                    zorder=1,
                )
        for i, label in enumerate(self.row_labels, start=1):
            ax.annotate(
                str(label),
                (self.shape[i - 1] + 0.2, 0.5 - i),
                ha="left",
                va="center",
                fontsize=8,
            )
        for j, label in enumerate(self.column_labels, start=1):
            depth = sum(1 for length in self.shape if length >= j)
            ax.annotate(
                str(label), (j - 0.5, -depth - 0.2), ha="center", va="top", fontsize=8
            )
        ax.set_xlim(-0.6, width + 0.8)
        ax.set_ylim(-k - 0.8, 0.6)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(title)

    def plot_diagram(self, ax: Axes | None = None) -> Axes:
        """Draw the filled Young diagram onto ``ax``.

        English notation with 0's and 1's in the boxes (Postnikov
        Definition 6.1), the ambient ``k x (n - k)`` rectangle dotted, and
        the border-path labels ``1..n`` on the south-east border
        (Postnikov section 2.1).

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        self._draw_frame(ax, "Le-diagram")
        for i, row in enumerate(self.filling, start=1):
            for j, value in enumerate(row, start=1):
                ax.annotate(
                    str(value),
                    (j - 0.5, 0.5 - i),
                    ha="center",
                    va="center",
                    fontsize=10,
                    color="black" if value else "0.6",
                )
        return ax

    def plot_hook_diagram(self, ax: Axes | None = None) -> Axes:
        """Draw the hook diagram onto ``ax``: a dot and a hook per 1-box.

        Postnikov section 6 (the hook formulation): from every dot two
        lines run right and down to the border of the shape; in a
        Le-diagram every intersection of two hook lines carries a dot.
        This is the Gamma-graph of Definition 6.3 drawn on the diagram.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        self._draw_frame(ax, "Hook diagram")
        for i, j in sorted(self.ones):
            depth = sum(1 for length in self.shape if length >= j)
            ax.plot(
                [j - 0.5, self.shape[i - 1]],
                [0.5 - i, 0.5 - i],
                color="0.3",
                linewidth=1.5,
                zorder=2,
            )
            ax.plot(
                [j - 0.5, j - 0.5],
                [0.5 - i, -depth],
                color="0.3",
                linewidth=1.5,
                zorder=2,
            )
        if self.ones:
            ax.scatter(
                [j - 0.5 for _, j in self.ones],
                [0.5 - i for i, _ in self.ones],
                color="black",
                zorder=3,
            )
        return ax

    def plot_pipe_dream(self, ax: Axes | None = None) -> Axes:
        """Draw the pipe dream onto ``ax``: elbows at 1's, crossings at 0's.

        Ardila-Rincon-Williams Lemma 4.8, steps (1)-(3): each 1 is replaced
        by an elbow joint and each 0 by a crossing, and the border labels
        are copied to the opposite (north and west) border, where the
        pipes enter.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        self._draw_frame(ax, "Pipe dream")
        for i, row in enumerate(self.filling, start=1):
            for j, value in enumerate(row, start=1):
                x, y = j - 0.5, 0.5 - i
                west, east = (x - 0.5, y), (x + 0.5, y)
                north, south = (x, y + 0.5), (x, y - 0.5)
                pairs = (
                    [(west, south), (north, east)]
                    if value
                    else [(west, east), (north, south)]
                )
                for start, end in pairs:
                    ax.plot(
                        [start[0], end[0]],
                        [start[1], end[1]],
                        color="0.3",
                        linewidth=1.5,
                        zorder=2,
                    )
        for i, label in enumerate(self.row_labels, start=1):
            ax.annotate(
                str(label), (-0.2, 0.5 - i), ha="right", va="center", fontsize=8
            )
        for j, label in enumerate(self.column_labels, start=1):
            ax.annotate(
                str(label), (j - 0.5, 0.2), ha="center", va="bottom", fontsize=8
            )
        return ax


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
def _check_type(k: int, n: int) -> None:
    """Reject types outside ``0 <= k <= n``."""
    if not 0 <= k <= n:
        msg = f"a type (k, n) needs 0 <= k <= n; got k = {k}, n = {n}"
        raise ValueError(msg)


def empty_filling(shape: Iterable[int], n: int) -> LeDiagram:
    """Return the shape filled with all 0's — a 0-dimensional cell.

    Le-Diagram page, canonical examples: every shape in the box filled with
    0's is a Le-diagram of rank 0, and there are exactly ``binom(n, k)`` of
    them, one per partition in the box.

    Raises:
        ValueError: If the shape does not fit the type ``(len(shape), n)``.
    """
    return LeDiagram(tuple(shape), n)


def top_cell_diagram(k: int, n: int) -> LeDiagram:
    """Return the full ``k x (n - k)`` rectangle filled with 1's — the top cell.

    Le-Diagram page, canonical examples: the Le-property holds vacuously
    and ``|D| = k(n - k)``; the positroid is the uniform matroid
    ``U_{k,n}``.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    _check_type(k, n)
    width = n - k
    ones = frozenset(itertools.product(range(1, k + 1), range(1, width + 1)))
    return LeDiagram((width,) * k, n, ones)


def arw_section_4_3_diagram() -> LeDiagram:
    """Return the Ardila-Rincon-Williams section 4.3 worked example.

    Shape 5532 in the ``4 x 6`` rectangle (``d = 4``, ``n = 10``) with rows
    ``0+0+0 / +++++ / 000 / ++`` (their Figures 1-2); its pipe dream is the
    decorated permutation ``(_1, 7, 9, 3, 2, ^6, 5, 10, 4, 8)`` of
    :func:`research.decorated_permutation.arw_section_4_3_example`.
    """
    return LeDiagram.from_plus_filling(["0+0+0", "+++++", "000", "++"], 10)


def forbidden_pattern_fillings() -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return the page's non-example: the forbidden pattern in a 2 x 2 square.

    Le-Diagram page, canonical examples: ``(1, 2) = 1``, ``(2, 1) = 1``,
    ``(2, 2) = 0`` with either value at ``(1, 1)`` — the 0 at ``(2, 2)``
    has a 1 above it and a 1 to its left. These are raw fillings, not
    Le-diagrams: :meth:`LeDiagram.from_filling` rejects both.
    """
    return (((0, 1), (1, 0)), ((1, 1), (1, 0)))


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #
def enumerate_shapes(k: int, n: int) -> Iterator[tuple[int, ...]]:
    """Yield the ``binom(n, k)`` shapes with ``k`` rows in the ``k x (n - k)`` box.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    _check_type(k, n)
    for lengths in itertools.combinations_with_replacement(range(n - k, -1, -1), k):
        yield tuple(lengths)


def enumerate_le_diagrams_of_shape(shape: Sequence[int], n: int) -> Iterator[LeDiagram]:
    """Yield every Le-diagram of one shape and type ``(len(shape), n)``, lazily.

    Restriction by shape (Postnikov Theorem 6.5): these index exactly the
    cells contained in the Schubert cell ``Omega_lambda``. Boxes are filled
    in reading order, a 0 being allowed unless a 1 already sits above it
    and to its left (the forbidden-pattern statement), so only Le-diagrams
    are ever built; the iterator is single-use.

    Raises:
        ValueError: If the shape does not fit the type.
    """
    base = LeDiagram(tuple(shape), n)
    boxes = [
        (i, j)
        for i, length in enumerate(base.shape, start=1)
        for j in range(1, length + 1)
    ]

    def extend(index: int, ones: frozenset[Box]) -> Iterator[frozenset[Box]]:
        if index == len(boxes):
            yield ones
            return
        i, j = boxes[index]
        yield from extend(index + 1, ones | {(i, j)})
        seen_left = any(row == i for row, _ in ones)
        seen_above = any(column == j for _, column in ones)
        if not (seen_left and seen_above):
            yield from extend(index + 1, ones)

    for ones in extend(0, frozenset()):
        yield LeDiagram(base.shape, n, ones)


def enumerate_le_diagrams(k: int, n: int) -> Iterator[LeDiagram]:
    """Yield every Le-diagram of type ``(k, n)`` — Postnikov's ``Le_kn``, lazily.

    Their number is the entry ``(n, k)`` of OEIS A046802 (Postnikov
    Proposition 23.1), which grows quickly; the iterator is single-use.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    for shape in enumerate_shapes(k, n):
        yield from enumerate_le_diagrams_of_shape(shape, n)


def enumerate_permutation_tableaux(k: int, n: int) -> Iterator[LeDiagram]:
    """Yield the permutation tableaux among the Le-diagrams of type ``(k, n)``.

    Column restriction (Steingrimsson-Williams Theorem 11): they biject
    with the permutations in ``S_n`` with ``k`` weak excedances. Lazy and
    single-use.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    return (d for d in enumerate_le_diagrams(k, n) if d.is_permutation_tableau)


def rank_generating_polynomial(k: int, n: int) -> tuple[int, ...]:
    """Return the brute-force coefficients of ``A_{k,n}(q)``, index = power.

    Williams 2005: ``A_{k,n}(q)`` is the generating function of the
    Le-diagrams in the ``k x (n - k)`` box graded by rank (cells by
    dimension). This counts them by enumeration; the closed form of her
    Theorem 4.1 is
    :func:`research.decorated_permutation.dimension_generating_polynomial`.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    _check_type(k, n)
    counts = [0] * (k * (n - k) + 1)
    for diagram in enumerate_le_diagrams(k, n):
        counts[diagram.rank] += 1
    return tuple(counts)
