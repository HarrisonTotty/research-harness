r"""Decorated permutations: permutations with two-colored fixed points.

A decorated permutation is a permutation of ``[n]`` together with a coloring
of each fixed point by one of two colors (Postnikov, *Total positivity,
Grassmannians, and networks*, arXiv:math/0609764, 2006, Definition 13.3). It
is the permutation-shaped index of the cells of the totally nonnegative
Grassmannian: positroids, Grassmann necklaces, and bounded affine
permutations carry the same information. :class:`DecoratedPermutation`
stores the pair literally and derives the rest of the Decorated Permutation
page as views — anti-exceedances and type, crossings and alignments, the
affinization and its length, the cell dimension, and the two Grassmann
necklace maps — alongside the page's enumeration formulas and canonical
examples.

Two conventions matter everywhere, and the sources disagree on both (the
page's convention warnings):

* **Color.** This library calls the *counted* color — the fixed points that
  count toward the type and are coloops of the positroid — **clockwise**.
  That is Postnikov's white (``col = -1``, a clockwise loop), the
  Fomin-Williams-Zelevinsky overline, the Knutson-Lam-Speyer color ``+1``
  (``f(i) = i + n``), and what Williams and Ardila-Rincon-Williams call
  "counterclockwise". The other color (loops of the positroid) is
  counterclockwise here: Postnikov's black (``col = 1``), the underline.
* **Direction.** Postnikov's statistics (anti-exceedances, alignments, the
  affinization, his necklace map) and Ardila-Rincon-Williams's (weak
  excedances, their necklace map) describe the same positroid through
  *inverse* permutations. Every method here applies its source's formula to
  the stored permutation literally, so nothing is inverted behind the
  caller's back; :class:`NecklaceConvention` selects the necklace map and
  :meth:`DecoratedPermutation.inverse` moves between the two directions.
  :mod:`research.positroid` and :mod:`research.grassmann_necklace` exchange
  decorated permutations in the Ardila-Rincon-Williams direction, so for
  such an object ``sigma`` the positroid's rank is
  ``sigma.weak_excedance_count`` and its cell dimension is
  ``sigma.inverse().dimension``.

Conversions to the :class:`research.grassmann_necklace.GrassmannNecklace`
and :class:`research.positroid.Positroid` classes live on those classes,
which import this module; here necklaces are plain tuples of frozensets and
bounded affine permutations plain window tuples.
"""

import enum
import functools
import itertools
import math
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

from research._cyclic import check_necklace_conditions
from research._plot import ensure_axes, scatter_labeled, unit_circle

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "DecoratedPermutation",
    "NecklaceConvention",
    "arw_section_4_3_example",
    "count_decorated_permutations",
    "count_decorated_permutations_of_type",
    "dimension_generating_polynomial",
    "enumerate_decorated_permutations",
    "fwz_definition_7_4_20_example",
    "fwz_example_7_9_2",
    "postnikov_section_16_example",
    "q_eulerian_polynomial",
    "top_cell",
    "zero_dimensional_cell",
]

_CLOCKWISE = "clockwise"
_COUNTERCLOCKWISE = "counterclockwise"
_COLUMNS = ("position", "target", "decoration")


class NecklaceConvention(enum.Enum):
    r"""Which way a Grassmann necklace transition is read as a permutation.

    For a transition ``I_{i+1} = (I_i \ {i}) + {j}`` with ``j != i``:

    * ``POSTNIKOV`` sets ``pi(i) = j`` (Postnikov section 16, Lemma 16.2),
      and ``I_r`` is the shifted anti-exceedance set.
    * ``ARW`` sets ``pi(j) = i`` — the inverse permutation — and ``I_k`` is
      the set of weak ``k``-excedances (Ardila-Rincon-Williams section 4.2,
      Proposition 4.6). This is the direction the rest of the library uses.
    """

    POSTNIKOV = "postnikov"
    ARW = "arw"


def _cyclic_between(x: int, start: int, end: int, n: int) -> bool:
    """Return whether ``x`` lies in the closed cyclic interval ``[start, end]``.

    The interval runs clockwise (increasing labels, wrapping at ``n``) from
    ``start`` to ``end`` inclusive; ``[a, a]`` is the single point ``a``.
    """
    return (x - start) % n <= (end - start) % n


def _shifted_less(a: int, b: int, start: int, n: int) -> bool:
    """Return whether ``a <_start b`` in the cyclic shift of the usual order.

    ``<_r`` orders ``[n]`` as ``r < r + 1 < ... < n < 1 < ... < r - 1``
    (Postnikov section 16).
    """
    return (a - start) % n < (b - start) % n


# --------------------------------------------------------------------------- #
# The structure
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DecoratedPermutation:
    """A bijection of ``[n]`` whose fixed points are each colored.

    ``targets[i - 1]`` is the image of position ``i`` in ``1..n``; fixed
    points in ``clockwise_fixed`` are colored clockwise and the remaining
    fixed points counterclockwise (Postnikov Definition 13.3, storing the
    coloring as the set of fixed points of one color). Clockwise is the
    *counted* color — Postnikov's white ``col = -1``, the overline of
    Fomin-Williams-Zelevinsky, a coloop of the positroid; see the module
    docstring for the full dictionary and for the direction convention.

    Unlike the larger structures of this library the definition is checked
    by ``__post_init__`` (two linear scans), so direct construction is
    validated too; the ``from_*`` classmethods convert the page's other
    formulations.

    Raises:
        ValueError: If ``targets`` is not a bijection of ``[n]`` or a
            decorated position is not a fixed point.
    """

    targets: tuple[int, ...]
    clockwise_fixed: frozenset[int] = frozenset()

    def __post_init__(self) -> None:
        """Validate the definition (both checks are linear-time)."""
        n = len(self.targets)
        if sorted(self.targets) != list(range(1, n + 1)):
            msg = (
                f"a decorated permutation must be a bijection of [n] "
                f"(Postnikov Definition 13.3); got targets {self.targets!r}"
            )
            raise ValueError(msg)
        stray = self.clockwise_fixed - self.fixed_points
        if stray:
            msg = (
                f"only fixed points carry a decoration (Postnikov Definition "
                f"13.3), but positions {sorted(stray)!r} are not fixed by "
                f"{self.targets!r}"
            )
            raise ValueError(msg)

    # -------------------------------------------------------- constructors
    @classmethod
    def from_colors(
        cls, targets: Iterable[int], colors: Mapping[int, int]
    ) -> DecoratedPermutation:
        """Build from Postnikov's coloring function ``col``.

        Postnikov Definition 13.3: ``col`` maps the fixed points to
        ``{1, -1}``, with ``1`` black and ``-1`` white. White fixed points
        are the counted ones (Postnikov section 16), so they become this
        library's clockwise color.

        Args:
            targets: The one-line images ``pi(1), ..., pi(n)``.
            colors: The value of ``col`` at every fixed point.

        Raises:
            ValueError: If ``targets`` is not a bijection, ``colors`` is not
                defined on exactly the fixed points, or a color is not
                ``1`` or ``-1``.
        """
        images = tuple(targets)
        bad = {i: c for i, c in colors.items() if c not in (1, -1)}
        if bad:
            msg = (
                f"the coloring function takes values in {{1, -1}} "
                f"(Postnikov Definition 13.3); got {bad!r}"
            )
            raise ValueError(msg)
        plain = cls(images)
        if frozenset(colors) != plain.fixed_points:
            msg = (
                f"the coloring function is defined on exactly the fixed "
                f"points {sorted(plain.fixed_points)!r} (Postnikov Definition "
                f"13.3); got colors for {sorted(colors)!r}"
            )
            raise ValueError(msg)
        return cls(images, frozenset(i for i, c in colors.items() if c == -1))

    @classmethod
    def from_underline_overline(
        cls,
        targets: Iterable[int],
        *,
        overlined: Collection[int] = (),
        underlined: Collection[int] = (),
    ) -> DecoratedPermutation:
        """Build from the overline/underline notation.

        Fomin-Williams-Zelevinsky Definition 7.4.20 (and
        Ardila-Rincon-Williams Definition 4.5): every fixed point is
        decorated by an overline or an underline. The overlined fixed points
        are the counted ones (FWZ Definition 7.9.1), so they become this
        library's clockwise color.

        Args:
            targets: The one-line images with the decorations stripped.
            overlined: The fixed points written with an overline.
            underlined: The fixed points written with an underline.

        Raises:
            ValueError: If ``targets`` is not a bijection or the two sets do
                not partition the fixed points.
        """
        images = tuple(targets)
        over, under = frozenset(overlined), frozenset(underlined)
        plain = cls(images)
        if over & under or over | under != plain.fixed_points:
            msg = (
                f"every fixed point carries exactly one of an overline or an "
                f"underline (Fomin-Williams-Zelevinsky Definition 7.4.20); "
                f"fixed points are {sorted(plain.fixed_points)!r}, got "
                f"overlined {sorted(over)!r} and underlined {sorted(under)!r}"
            )
            raise ValueError(msg)
        return cls(images, over)

    @classmethod
    def from_bounded_affine_permutation(
        cls, window: Iterable[int]
    ) -> DecoratedPermutation:
        """Build from the window ``(f(1), ..., f(n))`` of a bounded affine permutation.

        Knutson-Lam-Speyer section 3.2: reduce ``f`` modulo ``n`` and color
        the fixed points by whether ``f(i) = i`` or ``f(i) = i + n``; the
        ``f(i) = i + n`` fixed points are the counted ones, this library's
        clockwise color. Inverse to :meth:`to_bounded_affine_permutation`
        (Fomin-Williams-Zelevinsky Lemma 7.9.6).

        Raises:
            ValueError: If some ``f(i)`` is outside ``[i, i + n]`` or the
                residues modulo ``n`` are not a bijection of ``[n]``.
        """
        values = tuple(window)
        n = len(values)
        for i, value in enumerate(values, start=1):
            if not i <= value <= i + n:
                msg = (
                    f"a bounded affine permutation needs i <= f(i) <= i + n "
                    f"(Knutson-Lam-Speyer section 3.2); got f({i}) = {value} "
                    f"for n = {n}"
                )
                raise ValueError(msg)
        targets = tuple((value - 1) % n + 1 for value in values)
        if sorted(targets) != list(range(1, n + 1)):
            msg = (
                f"a bounded affine permutation must reduce to a bijection of "
                f"[n] modulo n (Knutson-Lam-Speyer section 3.2); got window "
                f"{values!r}"
            )
            raise ValueError(msg)
        clockwise = frozenset(
            i for i, value in enumerate(values, start=1) if value == i + n
        )
        return cls(targets, clockwise)

    @classmethod
    def from_grassmann_necklace(
        cls,
        entries: Sequence[Collection[int]],
        *,
        convention: NecklaceConvention,
    ) -> DecoratedPermutation:
        r"""Build the decorated permutation of a Grassmann necklace on ``[n]``.

        Postnikov section 16: if ``I_{i+1} = (I_i \ {i}) + {j}`` with
        ``j != i`` then ``pi(i) = j``; if ``I_{i+1} = I_i`` and ``i`` is not
        in ``I_i`` then ``i`` is a black (here counterclockwise) fixed
        point; if ``I_{i+1} = I_i`` and ``i`` is in ``I_i`` then ``i`` is a
        white (clockwise) fixed point. Under ``NecklaceConvention.ARW`` the
        transition instead sets ``pi(j) = i`` (Ardila-Rincon-Williams
        section 4.2). Inverse to :meth:`grassmann_necklace` (Postnikov
        Lemma 16.2; Ardila-Rincon-Williams Proposition 4.6).

        Args:
            entries: The necklace ``(I_1, ..., I_n)`` as subsets of ``[n]``.
            convention: Which way the transitions are read.

        Raises:
            ValueError: If an entry is not a subset of ``[n]`` or the
                necklace conditions of Postnikov section 16 fail.
        """
        n = len(entries)
        sets = [frozenset(entry) for entry in entries]
        ground = frozenset(range(1, n + 1))
        for r, entry in enumerate(sets, start=1):
            if not entry <= ground:
                msg = (
                    f"necklace entries are subsets of [n] = [{n}]; I_{r} "
                    f"contains {sorted(entry - ground)!r}"
                )
                raise ValueError(msg)
        masks = [sum(1 << (e - 1) for e in entry) for entry in sets]
        check_necklace_conditions(tuple(range(1, n + 1)), masks)
        targets = list(range(1, n + 1))
        clockwise: set[int] = set()
        for i in range(1, n + 1):
            current = sets[i - 1]
            if i not in current:
                continue
            (j,) = sets[i % n] - (current - {i})
            if j == i:
                clockwise.add(i)
            elif convention is NecklaceConvention.POSTNIKOV:
                targets[i - 1] = j
            else:
                targets[j - 1] = i
        return cls(tuple(targets), frozenset(clockwise))

    # ------------------------------------------------- computed properties
    @property
    def size(self) -> int:
        """The number ``n`` of letters the permutation acts on."""
        return len(self.targets)

    @functools.cached_property
    def fixed_points(self) -> frozenset[int]:
        """The positions ``i`` with ``pi(i) = i``."""
        return frozenset(
            i for i, target in enumerate(self.targets, start=1) if target == i
        )

    @property
    def counterclockwise_fixed(self) -> frozenset[int]:
        """The fixed points not colored clockwise (loops, here)."""
        return self.fixed_points - self.clockwise_fixed

    @property
    def colors(self) -> dict[int, int]:
        """Postnikov's coloring function ``col`` on the fixed points.

        Postnikov Definition 13.3 with his section 16 drawing convention:
        clockwise (white, counted) fixed points get ``-1`` and
        counterclockwise (black) ones ``1``.
        """
        return {i: -1 if i in self.clockwise_fixed else 1 for i in self.fixed_points}

    @functools.cached_property
    def _sources(self) -> tuple[int, ...]:
        """The one-line images of the inverse permutation."""
        sources = [0] * len(self.targets)
        for i, target in enumerate(self.targets, start=1):
            sources[target - 1] = i
        return tuple(sources)

    @property
    def anti_exceedances(self) -> frozenset[int]:
        """The ``i`` with ``pi^{-1}(i) > i``, plus the clockwise fixed points.

        Postnikov section 16: "``i`` is an anti-exceedance of ``pi`` if
        ``pi^{-1}(i) > i`` or ``pi(i) = i`` and ``col(i) = -1``" — white
        fixed points, this library's clockwise ones, count.
        """
        strict = frozenset(
            i for i, source in enumerate(self._sources, start=1) if source > i
        )
        return strict | self.clockwise_fixed

    @property
    def anti_exceedance_count(self) -> int:
        """The number of anti-exceedances.

        Equals ``#{i : pi(i) < i or i is overlined}``
        (Fomin-Williams-Zelevinsky Definition 7.9.1).
        """
        return len(self.anti_exceedances)

    @property
    def permutation_type(self) -> tuple[int, int]:
        """The type ``(k, n)``: size ``n`` with exactly ``k`` anti-exceedances.

        Postnikov section 17. For a decorated permutation in the
        Ardila-Rincon-Williams direction the positroid's rank is
        :attr:`weak_excedance_count` instead — the type of :meth:`inverse`.
        """
        return (self.anti_exceedance_count, len(self.targets))

    @property
    def weak_excedances(self) -> frozenset[int]:
        """Positions with ``pi(i) > i``, plus the clockwise fixed points.

        Williams (*Enumeration of totally positive Grassmann cells*, 2005,
        section 5): a weak excedence at ``i`` means ``pi(i) > i``, or
        ``pi(i) = i`` with the counted decoration — clockwise under this
        library's convention. Rank-``d`` positroids correspond to decorated
        permutations with exactly ``d`` weak excedances under the
        Ardila-Rincon-Williams necklace map (their Proposition 4.6).
        """
        strict = frozenset(
            i for i, target in enumerate(self.targets, start=1) if target > i
        )
        return strict | self.clockwise_fixed

    @property
    def weak_excedance_count(self) -> int:
        """The number of weak excedances — the rank in the ARW direction."""
        return len(self.weak_excedances)

    @functools.cached_property
    def crossings(self) -> tuple[tuple[int, int], ...]:
        """The pairs ``i < j`` whose chords cross, in sorted order.

        Postnikov section 17: chords ``(b_i, b_pi(i))`` and
        ``(b_j, b_pi(j))`` form a crossing when ``pi(j)`` lies in the cyclic
        interval ``[i, pi(i)]`` and ``j`` in ``[pi(i), i]`` — for either
        labelling of the pair. Endpoints may coincide (a 2-cycle crosses
        itself), but "a loop can never participate in a crossing".
        Quadratic in ``n``.
        """
        n = len(self.targets)
        moved = [i for i in range(1, n + 1) if i not in self.fixed_points]
        found: set[tuple[int, int]] = set()
        for i, j in itertools.permutations(moved, 2):
            image_i, image_j = self.targets[i - 1], self.targets[j - 1]
            if _cyclic_between(image_j, i, image_i, n) and _cyclic_between(
                j, image_i, i, n
            ):
                found.add((min(i, j), max(i, j)))
        return tuple(sorted(found))

    @functools.cached_property
    def alignments(self) -> tuple[tuple[int, int], ...]:
        """The ordered pairs ``(i, j)`` forming an alignment, in sorted order.

        Postnikov section 17: ``(i, j)`` is an alignment when ``pi(i)`` lies
        in the cyclic interval ``[i, pi(j)]`` and ``j`` in ``[pi(j), i]``,
        where a fixed ``i`` must be black (a counterclockwise loop) and a
        fixed ``j`` white (a clockwise loop). In particular "any
        counterclockwise loop forms an alignment with any clockwise loop".
        Quadratic in ``n``.
        """
        n = len(self.targets)
        found: list[tuple[int, int]] = []
        for i, j in itertools.permutations(range(1, n + 1), 2):
            if i in self.clockwise_fixed:
                continue
            if j in self.fixed_points and j not in self.clockwise_fixed:
                continue
            image_i, image_j = self.targets[i - 1], self.targets[j - 1]
            if _cyclic_between(image_i, i, image_j, n) and _cyclic_between(
                j, image_j, i, n
            ):
                found.append((i, j))
        return tuple(found)

    @property
    def alignment_number(self) -> int:
        """The alignment number ``A``: the total number of aligned pairs.

        Postnikov section 17. Equals :attr:`affine_length`
        (Fomin-Williams-Zelevinsky Definition 7.9.10).
        """
        return len(self.alignments)

    @property
    def affine_length(self) -> int:
        """The length of the affinization, counted from its inversions.

        Fomin-Williams-Zelevinsky Definition 7.9.10: an inversion is a pair
        ``i < j`` with ``f(i) > f(j)``, two inversions are equivalent when
        they differ by a common multiple of ``n``, and the length is the
        number of classes — one representative per ``i`` in ``1..n``. Only
        ``j < i + n`` can be inverted because ``j <= f(j)`` and
        ``f(i) <= i + n``. Quadratic in ``n``.
        """
        window = self.to_bounded_affine_permutation()
        n = len(window)

        def value(j: int) -> int:
            return window[(j - 1) % n] + n * ((j - 1) // n)

        return sum(
            1
            for i in range(1, n + 1)
            for j in range(i + 1, i + n)
            if window[i - 1] > value(j)
        )

    @property
    def dimension(self) -> int:
        """The dimension ``k(n - k) - A`` of the cell this permutation indexes.

        Postnikov Proposition 17.10, with ``(k, n)`` the
        :attr:`permutation_type`; also the rank in the circular Bruhat
        order. This is the Postnikov direction: for a decorated permutation
        ``sigma`` obtained from :mod:`research.positroid` or
        :mod:`research.grassmann_necklace`, the cell dimension is
        ``sigma.inverse().dimension``.
        """
        k, n = self.permutation_type
        return k * (n - k) - self.alignment_number

    @property
    def is_regular(self) -> bool:
        """Whether every fixed point carries the counted color.

        Williams section 5 calls a decorated permutation regular when all
        its fixed points are "counterclockwise" in her labels — the counted
        color, clockwise here — so it is a plain permutation with the
        weak-excedance convention ``pi(i) >= i``.
        """
        return self.clockwise_fixed == self.fixed_points

    @property
    def is_stabilized_interval_free(self) -> bool:
        """Whether ``pi(I) != I`` for every proper interval ``I`` of ``[n]``.

        Ardila-Rincon-Williams Definition 7.10, after Callan. Connected
        positroids correspond to stabilized-interval-free permutations
        (their Corollary 7.11; OEIS A075834). Quadratic in ``n``.
        """
        n = len(self.targets)
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                if b - a + 1 == n:
                    continue
                if set(self.targets[a - 1 : b]) == set(range(a, b + 1)):
                    return False
        return True

    @property
    def stabilizes_proper_cyclic_interval(self) -> bool:
        """Whether ``pi(I) = I`` for some proper cyclic interval ``I`` of ``[n]``.

        Ardila-Rincon-Williams Corollary 7.11: a positroid is connected iff
        its decorated permutation stabilizes no proper cyclic interval;
        since a cyclic interval or its complement is a genuine interval,
        that is the negation of :attr:`is_stabilized_interval_free`.
        Intervals are nonempty. Quadratic in ``n``.
        """
        n = len(self.targets)
        for start in range(1, n + 1):
            interval: set[int] = set()
            images: set[int] = set()
            for offset in range(n - 1):
                position = (start - 1 + offset) % n + 1
                interval.add(position)
                images.add(self.targets[position - 1])
                if interval == images:
                    return True
        return False

    @functools.cached_property
    def noncrossing_partition(self) -> frozenset[frozenset[int]]:
        """The finest non-crossing partition with ``i`` and ``pi(i)`` together.

        Ardila-Rincon-Williams Corollary 7.9: this is the partition of the
        positroid into connected components — the connected components of
        the chord diagram. Computed by merging crossing cycles until no two
        blocks cross; blocks ``B`` and ``C`` cross when some ``a < b < c < d``
        has ``a, c`` in one and ``b, d`` in the other. Polynomial in ``n``.
        """
        n = len(self.targets)
        blocks: list[set[int]] = []
        seen: set[int] = set()
        for start in range(1, n + 1):
            if start in seen:
                continue
            cycle = {start}
            current = self.targets[start - 1]
            while current != start:
                cycle.add(current)
                current = self.targets[current - 1]
            seen |= cycle
            blocks.append(cycle)
        merged = True
        while merged:
            merged = False
            for first, second in itertools.combinations(blocks, 2):
                if _blocks_cross(first, second):
                    first.update(second)
                    blocks.remove(second)
                    merged = True
                    break
        return frozenset(frozenset(block) for block in blocks)

    # ------------------------------------------------------ transformations
    def inverse(self) -> DecoratedPermutation:
        """Return the inverse permutation with the same decoration.

        Anti-exceedances and weak excedances are exchanged by inversion
        (``pi^{-1}(i) > i`` says ``pi^{-1}`` has a strict excedance at
        ``i``, and fixed points are shared), and the
        Ardila-Rincon-Williams necklace map uses the inverse permutation of
        Postnikov's (their section 4.2) — so this is the change of
        direction convention. An involution.
        """
        return DecoratedPermutation(self._sources, self.clockwise_fixed)

    def direct_sum(self, other: DecoratedPermutation) -> DecoratedPermutation:
        """Return the direct sum, with ``other`` placed on ``n+1..n+m``.

        Ardila-Rincon-Williams Definition 7.7: the direct sum of decorated
        permutations of disjoint sets is the decorated permutation of the
        union restricting to each. Here the two sets are the consecutive
        intervals ``[1, n]`` and ``[n + 1, n + m]``, the simplest
        non-crossing case of their Proposition 7.8.
        """
        n = len(self.targets)
        targets = self.targets + tuple(t + n for t in other.targets)
        clockwise = self.clockwise_fixed | {i + n for i in other.clockwise_fixed}
        return DecoratedPermutation(targets, frozenset(clockwise))

    def cyclic_shift(self, steps: int = 1) -> DecoratedPermutation:
        """Return the rotation by ``steps``: conjugation by ``i -> i + steps``.

        Cyclic rotation of the ground set acts on decorated permutations
        (Knutson-Lam-Speyer Example 3.7): the chord ``i -> pi(i)`` becomes
        ``i + steps -> pi(i) + steps`` modulo ``n``, decorations carried
        along. ``cyclic_shift(n)`` is the identity.
        """
        n = len(self.targets)
        if n == 0:
            return self
        targets = [0] * n
        for i, target in enumerate(self.targets, start=1):
            targets[(i - 1 + steps) % n] = (target - 1 + steps) % n + 1
        clockwise = frozenset((i - 1 + steps) % n + 1 for i in self.clockwise_fixed)
        return DecoratedPermutation(tuple(targets), clockwise)

    def to_bounded_affine_permutation(self) -> tuple[int, ...]:
        """Return the window ``(f(1), ..., f(n))`` of the affinization.

        Fomin-Williams-Zelevinsky Definition 7.9.3: ``f(i) = pi(i)`` if
        ``pi(i) > i``; ``f(i) = i`` at an underlined (counterclockwise)
        fixed point; ``f(i) = pi(i) + n`` if ``pi(i) < i``; ``f(i) = i + n``
        at an overlined (clockwise) fixed point; extended by
        ``f(i + n) = f(i) + n``. A bijection onto the ``(a, n)``-bounded
        affine permutations, ``a`` the anti-exceedance count (their Lemma
        7.9.6, after Knutson-Lam-Speyer).
        """
        n = len(self.targets)
        window: list[int] = []
        for i, target in enumerate(self.targets, start=1):
            if target > i:
                window.append(target)
            elif target < i or i in self.clockwise_fixed:
                window.append(target + n)
            else:
                window.append(i)
        return tuple(window)

    def grassmann_necklace(
        self, convention: NecklaceConvention
    ) -> tuple[frozenset[int], ...]:
        """Return the Grassmann necklace ``(I_1, ..., I_n)`` on ``[n]``.

        Under ``NecklaceConvention.POSTNIKOV``, ``I_r`` is the shifted
        anti-exceedance set ``{i : i <_r pi^{-1}(i), or i is a white fixed
        point}`` (Postnikov section 16); black (counterclockwise) fixed
        points belong to no entry and white (clockwise) ones to all. Under
        ``NecklaceConvention.ARW``, ``I_k`` is the set of weak
        ``k``-excedances ``{i : i <_k pi(i), or i is a counted fixed point}``
        (Ardila-Rincon-Williams section 4.2) — the Postnikov necklace of
        :meth:`inverse`. Inverse to :meth:`from_grassmann_necklace`
        (Postnikov Lemma 16.2). Quadratic in ``n``.
        """
        n = len(self.targets)
        partners = (
            self._sources
            if convention is NecklaceConvention.POSTNIKOV
            else self.targets
        )
        return tuple(
            frozenset(
                i
                for i, partner in enumerate(partners, start=1)
                if _shifted_less(i, partner, r, n)
            )
            | self.clockwise_fixed
            for r in range(1, n + 1)
        )

    # -------------------------------------------------------- serialization
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per position.

        Columns are ``position`` (``1..n``, ascending), ``target``
        (``pi(position)``), and ``decoration`` — ``"clockwise"`` or
        ``"counterclockwise"`` at a fixed point and null elsewhere. The
        encoding survives ``experiments.io.write_result`` (records-oriented
        JSON).

        Returns:
            The tidy frame; invert with :meth:`from_dataframe`.
        """
        decorations: list[str | None] = [
            None
            if target != i
            else _CLOCKWISE
            if i in self.clockwise_fixed
            else _COUNTERCLOCKWISE
            for i, target in enumerate(self.targets, start=1)
        ]
        return pd.DataFrame(
            {
                "position": list(range(1, len(self.targets) + 1)),
                "target": list(self.targets),
                "decoration": pd.Series(decorations, dtype="object"),
            }
        )

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> DecoratedPermutation:
        """Rebuild a decorated permutation from :meth:`to_dataframe` output.

        Rows may come in any order, and a frame with no rows and no columns
        decodes to the empty decorated permutation.

        Raises:
            ValueError: If a column is missing, the positions are not
                exactly ``1..n``, a decoration is not one of the two names,
                a fixed point is undecorated, or the decoded data is not a
                decorated permutation.
        """
        if df.empty and len(df.columns) == 0:
            # Records-oriented JSON of the empty permutation has no columns.
            return DecoratedPermutation(())
        missing = [column for column in _COLUMNS if column not in df.columns]
        if missing:
            msg = f"decorated permutation frame is missing columns {missing!r}"
            raise ValueError(msg)
        rows = sorted(
            zip(df["position"], df["target"], df["decoration"], strict=True),
            key=lambda row: int(row[0]),
        )
        positions = [int(position) for position, _, _ in rows]
        if positions != list(range(1, len(rows) + 1)):
            msg = (
                f"decorated permutation frame needs each position 1..n "
                f"exactly once; got {positions!r}"
            )
            raise ValueError(msg)
        targets = tuple(int(target) for _, target, _ in rows)
        clockwise: set[int] = set()
        for position, target, decoration in rows:
            if pd.isna(decoration):
                if int(position) == int(target):
                    msg = (
                        f"fixed point {int(position)} has no decoration; every "
                        f"fixed point is colored (Postnikov Definition 13.3)"
                    )
                    raise ValueError(msg)
                continue
            if decoration not in (_CLOCKWISE, _COUNTERCLOCKWISE):
                msg = (
                    f"decoration must be {_CLOCKWISE!r}, {_COUNTERCLOCKWISE!r}, "
                    f"or null; got {decoration!r} at position {int(position)}"
                )
                raise ValueError(msg)
            if int(position) != int(target):
                msg = (
                    f"only fixed points carry a decoration (Postnikov "
                    f"Definition 13.3), but position {int(position)} maps to "
                    f"{int(target)}"
                )
                raise ValueError(msg)
            if decoration == _CLOCKWISE:
                clockwise.add(int(position))
        return DecoratedPermutation(targets, frozenset(clockwise))

    # -------------------------------------------------------- visualization
    def plot_chord_diagram(self, ax: Axes | None = None) -> Axes:
        """Draw the chord diagram onto ``ax``.

        Postnikov section 16: ``b_1, ..., b_n`` sit clockwise on a circle
        with a directed chord from ``b_i`` to ``b_pi(i)`` for each
        non-fixed ``i``; a white fixed point (clockwise here) gets a
        clockwise loop and a black one a counterclockwise loop. The loops
        are drawn outside the circle with an arrowhead showing their sense,
        clockwise ones solid and counterclockwise ones dashed.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        n = len(self.targets)
        points = unit_circle(n, phase=math.pi / 2, clockwise=True)
        boundary = unit_circle(120, phase=0.0)
        ax.plot(
            [x for x, _ in boundary] + [boundary[0][0]],
            [y for _, y in boundary] + [boundary[0][1]],
            color="0.8",
            linewidth=1,
            zorder=0,
        )
        for i, target in enumerate(self.targets, start=1):
            if target == i:
                self._draw_loop(ax, points[i - 1], clockwise=i in self.clockwise_fixed)
                continue
            ax.annotate(
                "",
                xy=points[target - 1],
                xytext=points[i - 1],
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": "0.3",
                    "shrinkA": 4,
                    "shrinkB": 4,
                    "connectionstyle": "arc3,rad=0.15",
                },
            )
        scatter_labeled(
            ax,
            points,
            [str(i) for i in range(1, n + 1)],
            [
                (-14 * x, -14 * y) if i in self.fixed_points else (14 * x, 14 * y)
                for i, (x, y) in enumerate(points, start=1)
            ],
        )
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title("Chord diagram")
        return ax

    @staticmethod
    def _draw_loop(ax: Axes, point: tuple[float, float], *, clockwise: bool) -> None:
        """Draw a small oriented loop at a boundary point, outside the circle."""
        radius = 0.16
        x, y = point
        center = (x * (1 + radius), y * (1 + radius))
        base = math.atan2(-y, -x)
        sense = -1.0 if clockwise else 1.0
        steps = 40
        angles = [base + sense * 2 * math.pi * s / steps for s in range(steps)]
        arc = [
            (center[0] + radius * math.cos(a), center[1] + radius * math.sin(a))
            for a in angles
        ]
        ax.plot(
            [p[0] for p in arc],
            [p[1] for p in arc],
            color="0.3",
            linewidth=1,
            linestyle="-" if clockwise else "--",
            zorder=1,
        )
        ax.annotate(
            "",
            xy=(x, y),
            xytext=arc[-1],
            arrowprops={
                "arrowstyle": "-|>",
                "color": "0.3",
                "shrinkA": 0,
                "shrinkB": 3,
            },
        )


def _blocks_cross(first: Collection[int], second: Collection[int]) -> bool:
    """Return whether two disjoint blocks cross: ``a < b < c < d`` interleaved."""
    labels = sorted([(x, 0) for x in first] + [(x, 1) for x in second])
    sides = [side for _, side in labels]
    collapsed = [side for side, _ in itertools.groupby(sides)]
    return len(collapsed) >= 4  # noqa: PLR2004 - the pattern a, b, c, d.


# --------------------------------------------------------------------------- #
# Canonical examples
# --------------------------------------------------------------------------- #
def postnikov_section_16_example() -> DecoratedPermutation:
    """Return Postnikov's section 16 example, in Postnikov's direction.

    ``pi = (3, 1, 5, 4, 2, 6)`` with ``col(4) = 1`` (black) and
    ``col(6) = -1`` (white). Its ``POSTNIKOV`` necklace begins
    ``I_1 = {1, 2, 6}``: the white fixed point 6 lies in every entry and
    the black fixed point 4 in none.
    """
    return DecoratedPermutation.from_colors((3, 1, 5, 4, 2, 6), {4: 1, 6: -1})


def fwz_example_7_9_2() -> DecoratedPermutation:
    """Return Fomin-Williams-Zelevinsky Example 7.9.2.

    ``(5, _2, ^3, 6, 4, 1)`` with 2 underlined and 3 overlined: three
    anti-exceedances, "namely, 1, 4, and ^3".
    """
    return DecoratedPermutation.from_underline_overline(
        (5, 2, 3, 6, 4, 1), overlined={3}, underlined={2}
    )


def fwz_definition_7_4_20_example() -> DecoratedPermutation:
    """Return the example of Fomin-Williams-Zelevinsky Definition 7.4.20.

    ``(3, 4, 5, 1, 2, ^6)``: six letters with a single overlined fixed
    point.
    """
    return DecoratedPermutation.from_underline_overline(
        (3, 4, 5, 1, 2, 6), overlined={6}
    )


def arw_section_4_3_example() -> DecoratedPermutation:
    """Return the Ardila-Rincon-Williams section 4.3 example, in their direction.

    ``(_1, 7, 9, 3, 2, ^6, 5, 10, 4, 8)``: the image of the Le-diagram of
    shape 5532 in a ``4 x 6`` box, so it has ``d = 4`` weak excedances.
    """
    return DecoratedPermutation.from_underline_overline(
        (1, 7, 9, 3, 2, 6, 5, 10, 4, 8), overlined={6}, underlined={1}
    )


def top_cell(k: int, n: int) -> DecoratedPermutation:
    """Return the top element ``i -> i + k (mod n)`` of type ``(k, n)``.

    Postnikov Lemma 17.6: the unique maximum of the circular Bruhat order,
    with all fixed points black for ``k = 0`` and all white for ``k = n``.
    It is the unique decorated permutation of its type with no alignments,
    so its cell has dimension ``k(n - k)`` (Postnikov Proposition 17.10).
    Postnikov's direction: the library-direction top cell is its
    :meth:`DecoratedPermutation.inverse`.

    Raises:
        ValueError: Unless ``0 <= k <= n``.
    """
    if not 0 <= k <= n:
        msg = f"a type (k, n) needs 0 <= k <= n; got k = {k}, n = {n}"
        raise ValueError(msg)
    targets = tuple((i + k - 1) % n + 1 for i in range(1, n + 1))
    clockwise: frozenset[int] = frozenset(range(1, n + 1) if k == n else ())
    return DecoratedPermutation(targets, clockwise)


def zero_dimensional_cell(n: int, white: Iterable[int]) -> DecoratedPermutation:
    """Return the identity of ``[n]`` with the fixed points in ``white`` counted.

    Postnikov Lemma 17.6: for ``white`` a ``k``-subset these are the
    ``binom(n, k)`` minimal elements of the circular Bruhat order, with
    ``A = k(n - k)``; the necklace is the constant necklace at ``white``.

    Raises:
        ValueError: If ``n`` is negative or ``white`` is not a subset of
            ``[n]``.
    """
    if n < 0:
        msg = f"size must be non-negative, got {n}"
        raise ValueError(msg)
    return DecoratedPermutation(tuple(range(1, n + 1)), frozenset(white))


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #
def enumerate_decorated_permutations(n: int) -> Iterator[DecoratedPermutation]:
    """Yield every decorated permutation of ``[n]``, lazily.

    Walks every permutation with every coloring of its fixed points;
    there are ``sum(n!/k!)`` of them (:func:`count_decorated_permutations`),
    so the iterator is factorial in ``n`` and single-use.

    Raises:
        ValueError: If ``n`` is negative.
    """
    if n < 0:
        msg = f"size must be non-negative, got {n}"
        raise ValueError(msg)
    for targets in itertools.permutations(range(1, n + 1)):
        fixed = [i for i, target in enumerate(targets, start=1) if target == i]
        for count in range(len(fixed) + 1):
            for clockwise in itertools.combinations(fixed, count):
                yield DecoratedPermutation(targets, frozenset(clockwise))


def count_decorated_permutations(n: int) -> int:
    """Return the number ``N_n = sum_{k=0}^{n} n!/k!`` of decorated permutations.

    Postnikov Proposition 23.2 (``N_n = n N_{n-1} + 1``, ``N_0 = 1``) in
    the closed form of Fomin-Williams-Zelevinsky Exercise 7.4.21; OEIS
    A000522.

    Raises:
        ValueError: If ``n`` is negative.
    """
    if n < 0:
        msg = f"size must be non-negative, got {n}"
        raise ValueError(msg)
    return sum(math.factorial(n) // math.factorial(k) for k in range(n + 1))


def _check_positive_type(k: int, n: int, formula: str) -> None:
    """Reject types outside ``1 <= k <= n``, where Williams's sums are empty."""
    if not 1 <= k <= n:
        msg = (
            f"{formula} is a sum over i = 0..k-1, which is empty for k = 0; "
            f"it needs 1 <= k <= n, got k = {k}, n = {n}"
        )
        raise ValueError(msg)


def count_decorated_permutations_of_type(k: int, n: int) -> int:
    r"""Return the number of decorated permutations of type ``(k, n)``.

    Williams (*Enumeration of totally positive Grassmann cells*, 2005,
    Theorem 4.1), the closed form
    ``D_{a,b} = sum_{i=0}^{a-1} (-1)^i binom(b, i) [(a-i)^i (a-i+1)^{b-i}
    - (a-i-1)^i (a-i)^{b-i}]``; the triangle is OEIS A046802.

    Raises:
        ValueError: Unless ``1 <= k <= n`` (the sum is empty at ``k = 0``,
            where the true count is 1).
    """
    _check_positive_type(k, n, "Williams's closed form D_{a,b}")
    return sum(
        (-1) ** i
        * math.comb(n, i)
        * (
            (k - i) ** i * (k - i + 1) ** (n - i)
            - (k - i - 1) ** i * (k - i) ** (n - i)
        )
        for i in range(k)
    )


type _Polynomial = list[int]


def _poly_mul(a: _Polynomial, b: _Polynomial) -> _Polynomial:
    """Multiply two coefficient lists (index = power)."""
    product = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            product[i + j] += x * y
    return product


def _poly_pow(a: _Polynomial, exponent: int) -> _Polynomial:
    """Raise a coefficient list to a non-negative power."""
    result = [1]
    for _ in range(exponent):
        result = _poly_mul(result, a)
    return result


def _q_integer(i: int) -> _Polynomial:
    """Return ``[i] = 1 + q + ... + q^{i-1}``, the zero polynomial for ``i = 0``."""
    return [1] * i if i > 0 else [0]


def _signed_sum(terms: Iterable[tuple[int, int, _Polynomial]]) -> _Polynomial:
    """Sum ``coefficient * q^shift * polynomial`` over the given terms."""
    total: dict[int, int] = {}
    for coefficient, shift, polynomial in terms:
        for power, value in enumerate(polynomial):
            total[power + shift] = total.get(power + shift, 0) + coefficient * value
    degree = max((p for p, v in total.items() if v), default=0)
    return [total.get(power, 0) for power in range(degree + 1)]


def _divide_by_q_power(polynomial: _Polynomial, power: int) -> tuple[int, ...]:
    """Divide by ``q^power``; the low coefficients must vanish."""
    low, high = polynomial[:power], polynomial[power:]
    if any(low):
        msg = f"polynomial {polynomial!r} is not divisible by q^{power}"
        raise ArithmeticError(msg)
    return tuple(high) if high else (0,)


def dimension_generating_polynomial(k: int, n: int) -> tuple[int, ...]:
    r"""Return the coefficients of Williams's ``A_{k,n}(q)``, index = power.

    Williams 2005, Main Theorem:
    ``A_{k,n}(q) = q^{-k^2} sum_{i=0}^{k-1} (-1)^i binom(n, i)
    (q^{ki} [k-i]^i [k-i+1]^{n-i} - q^{(k+1)i} [k-i-1]^i [k-i]^{n-i})``
    with ``[i] = 1 + q + ... + q^{i-1}``; the coefficient of ``q^r`` counts
    type-``(k, n)`` decorated permutations of dimension
    ``k(n - k) - A = r``.

    Raises:
        ValueError: Unless ``1 <= k <= n`` (the sum is empty at ``k = 0``,
            where the true polynomial is 1).
    """
    _check_positive_type(k, n, "Williams's A_{k,n}(q)")
    terms: list[tuple[int, int, _Polynomial]] = []
    for i in range(k):
        sign = (-1) ** i * math.comb(n, i)
        first = _poly_mul(
            _poly_pow(_q_integer(k - i), i), _poly_pow(_q_integer(k - i + 1), n - i)
        )
        second = _poly_mul(
            _poly_pow(_q_integer(k - i - 1), i), _poly_pow(_q_integer(k - i), n - i)
        )
        terms.append((sign, k * i, first))
        terms.append((-sign, (k + 1) * i, second))
    return _divide_by_q_power(_signed_sum(terms), k * k)


def q_eulerian_polynomial(k: int, n: int) -> tuple[int, ...]:
    r"""Return the coefficients of Williams's ``E_{k,n}(q)``, index = power.

    Williams 2005, section 5:
    ``E_{k,n}(q) = q^{n-k^2} sum_{i=0}^{k-1} binom(n, i) (-1)^i
    (q^{ki-i} [k-i]^n - q^{ki} [k-i-1]^n)``, the part of ``A_{k,n}(q)``
    contributed by regular decorated permutations, so that
    ``A_{k,n}(q) = sum_i binom(n, i) E_{k,n-i}(q)``. The normalized
    ``q^{k-n} E_{k,n}(q)`` interpolates between the Eulerian numbers, the
    Narayana numbers, and the binomial coefficients at ``q = 1, 0, -1``.

    Raises:
        ValueError: Unless ``1 <= k <= n``.
    """
    _check_positive_type(k, n, "Williams's E_{k,n}(q)")
    terms: list[tuple[int, int, _Polynomial]] = []
    for i in range(k):
        sign = (-1) ** i * math.comb(n, i)
        terms.append((sign, n + k * i - i, _poly_pow(_q_integer(k - i), n)))
        terms.append((-sign, n + k * i, _poly_pow(_q_integer(k - i - 1), n)))
    return _divide_by_q_power(_signed_sum(terms), k * k)
