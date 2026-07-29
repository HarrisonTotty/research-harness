"""Grassmann necklaces: cyclic sequences of Gale-minimal bases.

A Grassmann necklace of type ``(k, n)`` is a cyclic sequence
``I = (I_1, ..., I_n)`` of subsets of ``[n]`` obeying Postnikov's transition
conditions (N1)-(N2) (*Total positivity, Grassmannians, and networks*,
arXiv:math/0609764, 2006, Def. 16.1); equal entry cardinality is a
consequence of the conditions, not an axiom. This module stores the necklace
as a first-class value over a cyclically ordered ground set — sharing the
mask machinery of :mod:`research._cyclic` with :mod:`research.positroid` —
and implements the cryptomorphic dictionary of the Grassmann Necklace page:
the positroid of Oh's theorem, decorated permutations (Postnikov Lemma
16.2), bounded affine permutation windows (Lam Theorem 6.2), Oh's upper
necklace, and the juggling states of Knutson-Lam-Speyer.

Everything inherits the explicit, exponential design of the matroid module:
practical ground sets stop around 16 elements, and the enumeration helper
branches exponentially in ``n``.
"""

import functools
import itertools
from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, override

import pandas as pd

from research._bitmask import bits, fmt, indexed_ground_set, mask_from_labels
from research._cyclic import (
    check_necklace_conditions,
    decorated_necklace_masks,
    gale_geq,
    necklace_decorated,
    necklace_masks,
)
from research._plot import ensure_axes, scatter_labeled
from research.matroid import Matroid
from research.positroid import DecoratedPermutation, Positroid

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "GrassmannNecklace",
    "constant_necklace",
    "enumerate_necklaces",
    "kls_example_3_14",
    "lam_example_6_1",
    "oh_worked_example",
    "postnikov_figure_16_1",
    "uniform_necklace",
]


@dataclass(frozen=True, repr=False)
class GrassmannNecklace[T: Hashable]:
    """A Grassmann necklace over a cyclically ordered ground set.

    ``elements`` fixes the cyclic order and ``entry_masks`` stores the
    entries ``I_1, ..., I_n`` as bitmasks over element positions — the
    sequence form of Postnikov Def. 16.1 won as the representation because
    the transition conditions (N1)-(N2) validate in linear time on it and
    every derived view (Gale comparisons, transition readings, juggling
    shifts, interval ranks) is direct mask arithmetic. Build through the
    ``from_<formulation>`` constructors, which validate; calling
    ``GrassmannNecklace(...)`` directly skips validation. A necklace is a
    value: equality and hashing compare the cyclic order and the entries.

    Convention dictionary (Grassmann Necklace page, convention warning):
    :class:`research.positroid.DecoratedPermutation` views use this
    library's stored Ardila-Rincon-Williams direction — the inverse of
    Postnikov's permutation — with coloops as clockwise fixed points;
    bounded affine permutation windows, siteswaps, and juggling states use
    the Knutson-Lam-Speyer/Postnikov direction.
    """

    elements: tuple[T, ...]
    entry_masks: tuple[int, ...]

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        shown = ", ".join(fmt(mask, self.elements) for mask in self.entry_masks)
        return (
            f"GrassmannNecklace(k={self.rank}, n={len(self.elements)}, "
            f"entries=({shown}))"
        )

    # ---------------------------------------------------------- constructors
    @classmethod
    def from_entries(
        cls,
        elements: Iterable[T],
        entries: Sequence[Iterable[T]],
        *,
        validate: bool = True,
    ) -> GrassmannNecklace[T]:
        """Build a necklace from its entry sequence ``I_1, ..., I_n``.

        The primary formulation (Postnikov Def. 16.1): one entry per
        ground-set element, checked against the transition conditions
        (N1)-(N2) with indices modulo ``n``.

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            entries: The entries in cyclic order, as label collections.
            validate: Check the conditions (N1)-(N2).

        Returns:
            The validated necklace.

        Raises:
            ValueError: If the entry count or labels are wrong, or a
                transition condition fails (the message names it).
        """
        elems, index = indexed_ground_set(elements)
        masks = tuple(mask_from_labels(entry, index) for entry in entries)
        if len(masks) != len(elems):
            msg = (
                f"a Grassmann necklace has one entry per ground-set element "
                f"(Postnikov Def. 16.1); got {len(masks)} entries for "
                f"{len(elems)} elements"
            )
            raise ValueError(msg)
        if validate:
            check_necklace_conditions(elems, masks)
        return cls(elems, masks)

    @classmethod
    def from_matroid(cls, matroid: Matroid[T]) -> GrassmannNecklace[T]:
        """Return the necklace of a matroid's Gale-minimal shifted bases.

        Entry ``i`` is the lexicographically minimal (equivalently
        Gale-minimal) basis for the cyclic order starting at position ``i``;
        by Postnikov Lemma 16.3 the sequence satisfies (N1)-(N2) for *any*
        matroid — positroid or not — so no validation runs.

        Args:
            matroid: Any matroid; its stored element order is read as the
                cyclic order (a :class:`research.positroid.Positroid` works
                unchanged).

        Returns:
            The necklace ``I(M)``, certified by Lemma 16.3.
        """
        return cls(
            matroid.elements,
            necklace_masks(len(matroid.elements), matroid.independent_masks),
        )

    @classmethod
    def from_decorated_permutation(
        cls,
        elements: Iterable[T],
        decorated: DecoratedPermutation,
        *,
        validate: bool = True,
    ) -> GrassmannNecklace[T]:
        """Build the necklace of a decorated permutation.

        The inverse half of Postnikov Lemma 16.2, via the transition rule
        read backwards: ``j`` lies in ``I_k`` exactly when ``k`` is in the
        cyclic interval ``(pi(j), j]``, clockwise fixed points (coloops) in
        every entry and counterclockwise ones (loops) in none. ``decorated``
        is in this library's stored (Ardila-Rincon-Williams) direction.

        Args:
            elements: Ground-set labels; their order is the cyclic order,
                and position ``i`` of the permutation is ``elements[i-1]``.
            decorated: The decorated permutation of ``[n]``.
            validate: Re-check the derived conditions (N1)-(N2).

        Returns:
            The necklace ``I(pi)``.

        Raises:
            ValueError: If the ground-set size does not match the
                permutation.
        """
        elems, _ = indexed_ground_set(elements)
        n = len(elems)
        if len(decorated.targets) != n:
            msg = (
                f"the decorated permutation acts on [{len(decorated.targets)}] "
                f"but the ground set has {n} elements"
            )
            raise ValueError(msg)
        masks = decorated_necklace_masks(
            decorated.targets, decorated.clockwise_fixed, n
        )
        if validate:
            check_necklace_conditions(elems, masks)
        return cls(elems, masks)

    @classmethod
    def from_bounded_affine_permutation(
        cls,
        elements: Iterable[T],
        window: Sequence[int],
        *,
        validate: bool = True,
    ) -> GrassmannNecklace[T]:
        """Build the necklace of a bounded affine permutation window.

        Lam Theorem 6.2 (equivalently Knutson-Lam-Speyer Lemma 3.8):
        ``I_a = {f(b) mod n : b < a and f(b) >= a}`` under the periodic
        extension ``f(b + n) = f(b) + n``, residues taken in ``1..n``.

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            window: The window notation ``(f(1), ..., f(n))``, with
                ``a <= f(a) <= a + n`` and residues mod ``n`` a bijection
                of ``[n]`` (Knutson-Lam-Speyer section 3.2).
            validate: Re-check the derived conditions (N1)-(N2).

        Returns:
            The necklace ``I(f)``.

        Raises:
            ValueError: If the window length, bounds, or residues are wrong.
        """
        elems, _ = indexed_ground_set(elements)
        n = len(elems)
        win = tuple(window)
        if len(win) != n:
            msg = (
                f"a bounded affine permutation window has one value per "
                f"ground-set element; got {len(win)} values for {n} elements"
            )
            raise ValueError(msg)
        for a, value in enumerate(win, start=1):
            if not a <= value <= a + n:
                msg = (
                    f"a bounded affine permutation needs "
                    f"i <= f(i) <= i + n (Knutson-Lam-Speyer section 3.2); "
                    f"got f({a}) = {value} for n = {n}"
                )
                raise ValueError(msg)
        if sorted((value - 1) % n for value in win) != list(range(n)):
            msg = (
                f"a bounded affine permutation must reduce to a bijection "
                f"of [n] modulo n (Knutson-Lam-Speyer section 3.2); got "
                f"window {win!r}"
            )
            raise ValueError(msg)
        masks: list[int] = []
        for a in range(1, n + 1):
            mask = 0
            for b in range(a - n, a):
                base = (b - 1) % n
                landing = win[base] + (b - base - 1)
                if landing >= a:
                    mask |= 1 << (landing - 1) % n
            masks.append(mask)
        if validate:
            check_necklace_conditions(elems, tuple(masks))
        return cls(elems, tuple(masks))

    # ------------------------------------------------- computed properties
    @property
    def rank(self) -> int:
        """The common entry cardinality ``k``.

        Equal cardinality is a consequence of (N1)-(N2), not an axiom
        (Postnikov section 16).
        """
        return self.entry_masks[0].bit_count() if self.entry_masks else 0

    @property
    def necklace_type(self) -> tuple[int, int]:
        """The type ``(k, n)`` of the necklace (Postnikov section 17)."""
        return (self.rank, len(self.elements))

    @functools.cached_property
    def entries(self) -> tuple[frozenset[T], ...]:
        """The entries ``(I_1, ..., I_n)`` as label sets."""
        return tuple(
            frozenset(self.elements[b] for b in bits(mask)) for mask in self.entry_masks
        )

    @functools.cached_property
    def loops(self) -> frozenset[T]:
        """The elements in no entry — Postnikov's zeros.

        Under the decorated-permutation correspondence these are exactly
        the black (here counterclockwise) fixed points (Grassmann Necklace
        page, Correspondences).
        """
        union = 0
        for mask in self.entry_masks:
            union |= mask
        return frozenset(
            element
            for position, element in enumerate(self.elements)
            if not union >> position & 1
        )

    @functools.cached_property
    def coloops(self) -> frozenset[T]:
        """The elements in every entry — Postnikov's cozeros.

        Under the decorated-permutation correspondence these are exactly
        the white (here clockwise) fixed points (Grassmann Necklace page,
        Correspondences).
        """
        if not self.entry_masks:
            return frozenset()
        common = self.entry_masks[0]
        for mask in self.entry_masks:
            common &= mask
        return frozenset(self.elements[b] for b in bits(common))

    @property
    def is_constant(self) -> bool:
        """Whether the necklace is ``(I, I, ..., I)`` for a fixed entry.

        Constant necklaces are always valid (take ``j = i`` in (N1)) and
        are the necklaces of single-basis positroids (Grassmann Necklace
        page, derived vocabulary; Lam Lemma 8.3).
        """
        return len(set(self.entry_masks)) <= 1

    @functools.cached_property
    def upper_necklace(self) -> tuple[frozenset[T], ...]:
        """Oh's upper Grassmann necklace ``J_i = pi^{-1}(I_i)``.

        The companion sequence of Oh Theorem 19, with ``pi`` in Oh's
        (Postnikov) convention — whose inverse is exactly this library's
        stored direction, so the stored targets apply forward. Not itself
        claimed to satisfy (N1)-(N2); its role is that the positroid is
        equally the intersection of the cyclically shifted *dual* Schubert
        matroids ``{H : H <=_i J_i}``.
        """
        decorated = self.to_decorated_permutation()
        return tuple(
            frozenset(self.elements[decorated.targets[b] - 1] for b in bits(mask))
            for mask in self.entry_masks
        )

    @functools.cached_property
    def juggling_states(self) -> tuple[frozenset[int], ...]:
        """The Knutson-Lam-Speyer juggling states, as 1-based positions.

        ``J_r = chi^{-(r-1)}(I_r)`` for the long cycle
        ``chi = [2, 3, ..., n, 1]`` — each entry rotated back to a fixed
        window instead of being read in a shifted order (Grassmann Necklace
        page, juggling-state form; Knutson-Lam-Speyer section 3). Read as
        the scheduled landing times of ``k`` airborne balls.
        """
        n = len(self.elements)
        return tuple(
            frozenset((b - shift) % n + 1 for b in bits(mask))
            for shift, mask in enumerate(self.entry_masks)
        )

    @functools.cached_property
    def siteswap(self) -> tuple[int, ...]:
        """One period ``f(i) - i`` of the bounded affine permutation.

        The juggling siteswap of Knutson-Lam-Speyer sections 1.2 and 3.3.
        """
        return tuple(
            landing - position
            for position, landing in enumerate(
                self.to_bounded_affine_permutation(), start=1
            )
        )

    def interval_rank(self, a: int, b: int) -> int:
        """Return ``r_ab = |I_a ∩ {a, a+1, ..., b}|`` over cyclic positions.

        The comparison statistic of Postnikov Corollary 17.7; ``a`` and
        ``b`` are 1-based positions in the cyclic order, and the interval
        wraps when ``b < a``.

        Raises:
            ValueError: If a position is out of range.
        """
        n = len(self.elements)
        if not (1 <= a <= n and 1 <= b <= n):
            msg = f"positions must be in 1..{n}, got a={a}, b={b}"
            raise ValueError(msg)
        interval = 0
        for offset in range((b - a) % n + 1):
            interval |= 1 << (a - 1 + offset) % n
        return (self.entry_masks[a - 1] & interval).bit_count()

    def circular_bruhat_leq(self, other: GrassmannNecklace[T]) -> bool:
        """Return whether ``self <= other`` in the circular Bruhat order.

        Postnikov Corollary 17.7: ``pi <= sigma`` in ``CB_kn`` iff
        ``r_ab(pi) <= r_ab(sigma)`` for all ``a, b in [n]``. Quadratically
        many interval counts.

        Raises:
            ValueError: If the necklaces do not share one cyclic order and
                type — ``CB_kn`` compares necklaces of a single type.
        """
        if self.elements != other.elements or self.rank != other.rank:
            msg = (
                f"the circular Bruhat order CB_kn compares necklaces of one "
                f"type on one cyclic order (Postnikov Cor. 17.7); got types "
                f"{self.necklace_type} and {other.necklace_type}"
            )
            raise ValueError(msg)
        n = len(self.elements)
        return all(
            self.interval_rank(a, b) <= other.interval_rank(a, b)
            for a in range(1, n + 1)
            for b in range(1, n + 1)
        )

    # ------------------------------------------------- the Lam partial order
    def _lam_leq(self, other: GrassmannNecklace[T]) -> bool:
        """Compare entrywise in the shifted Gale orders (Lam section 6.3)."""
        n = len(self.elements)
        return all(
            gale_geq(other.entry_masks[a], self.entry_masks[a], a, n) for a in range(n)
        )

    def __le__(self, other: object) -> bool:
        """Lam's necklace order: ``I <= I'`` iff ``I_a <=_a I'_a`` for all ``a``.

        The partial order of Lam section 6.3, matched order-reversingly
        with bounded affine permutations by his Theorem 6.2. Necklaces on
        different ground sets or of different types are incomparable
        (``False`` both ways).
        """
        if not isinstance(other, GrassmannNecklace):
            return NotImplemented
        if self.elements != other.elements or self.rank != other.rank:
            return False
        return self._lam_leq(other)

    def __lt__(self, other: object) -> bool:
        """Strict form of Lam's necklace order (Lam section 6.3)."""
        if not isinstance(other, GrassmannNecklace):
            return NotImplemented
        return self != other and self <= other

    def __ge__(self, other: object) -> bool:
        """Reverse of Lam's necklace order (Lam section 6.3)."""
        if not isinstance(other, GrassmannNecklace):
            return NotImplemented
        return other <= self

    def __gt__(self, other: object) -> bool:
        """Strict reverse of Lam's necklace order (Lam section 6.3)."""
        if not isinstance(other, GrassmannNecklace):
            return NotImplemented
        return other < self

    # ------------------------------------------------------- transformations
    def to_positroid(self) -> Positroid[T]:
        """Return the positroid ``M(I)`` of Oh's theorem.

        The bases are ``{B : B >=_j I_j for all j}`` in the shifted Gale
        orders — equivalently the intersection of the cyclically shifted
        Schubert matroids named by the entries (Oh 2011, Theorem 6; new
        proof by Lam, Theorem 8.4). The construction certifies the result,
        so no re-validation runs.
        """
        return Positroid.from_grassmann_necklace(
            self.elements, self.entries, validate=False
        )

    def to_matroid(self) -> Matroid[T]:
        """Return ``M(I)`` as a plain matroid, forgetting the cyclic reading."""
        return self.to_positroid().to_matroid()

    def to_decorated_permutation(self) -> DecoratedPermutation:
        r"""Return the decorated permutation indexing this necklace.

        The forward half of Postnikov Lemma 16.2, read off the transitions:
        when ``I_{i+1}`` is ``(I_i \ {i}) + {j}`` the permutation sends
        ``j`` to ``i`` (this library's stored direction — the inverse of
        Postnikov's, per the page's convention warning); a coloop is a
        clockwise fixed point and a loop a counterclockwise one.

        Returns:
            The decorated permutation on positions ``1..n``.

        Raises:
            ValueError: If the transitions are inconsistent, which
                indicates the instance was built without validation.
        """
        targets, clockwise = necklace_decorated(self.entry_masks, len(self.elements))
        return DecoratedPermutation(targets, clockwise)

    def to_bounded_affine_permutation(self) -> tuple[int, ...]:
        r"""Return the bounded affine permutation window ``(f(1), ..., f(n))``.

        The inverse map of Lam Theorem 6.2: if ``a`` is not in ``I_a`` then
        ``f(a) = a``; if ``a`` is in ``I_a`` and ``I_{a+1}`` is
        ``(I_a \ {a}) + {a'}`` then ``f(a)`` is the unique ``b`` congruent
        to ``a'`` mod ``n`` with ``a < b <= a + n`` (so a coloop gets
        ``f(a) = a + n``).

        Returns:
            The window values, one per position.

        Raises:
            ValueError: If the transitions are inconsistent, which
                indicates the instance was built without validation.
        """
        n = len(self.elements)
        window: list[int] = []
        for a in range(n):
            current = self.entry_masks[a]
            following = self.entry_masks[(a + 1) % n]
            if not current >> a & 1:
                window.append(a + 1)
                continue
            base = current ^ (1 << a)
            added = following & ~base
            if added.bit_count() != 1 or base & ~following:
                msg = (
                    f"inconsistent Grassmann necklace at position {a + 1}; "
                    f"was this instance built without validation?"
                )
                raise ValueError(msg)
            landing = added.bit_length()
            window.append(landing if landing > a + 1 else landing + n)
        return tuple(window)

    def cyclic_shift(self, steps: int = 1) -> GrassmannNecklace[T]:
        """Return the necklace of the positroid with its order rotated.

        Positroids are closed under cyclic shift (Ardila-Rincon-Williams
        Lemma 3.3); the sources state no necklace-level formula, so the map
        is the composite through Oh's construction — shift the positroid,
        take its necklace (Grassmann Necklace page, closure operations).
        """
        return GrassmannNecklace.from_matroid(self.to_positroid().cyclic_shift(steps))

    # ---------------------------------------------------------- dataframes
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per entry membership.

        Columns are ``element`` and ``entry`` (nullable ``Int64``). Each
        ground-set element gets one row with a null ``entry`` (declaring
        the cyclic order), and each (entry ``i``, member) incidence gets
        one row with ``entry = i`` (1-based). A rank-0 necklace therefore
        has only ground rows. The encoding survives
        ``experiments.io.write_result`` (records-oriented JSON) when read
        back with ``pd.read_json(path, dtype=False)``.

        Returns:
            The tidy incidence frame; invert with :meth:`from_dataframe`.
        """
        element_column: list[T] = list(self.elements)
        entry_column: list[int | None] = [None] * len(self.elements)
        for i, mask in enumerate(self.entry_masks, start=1):
            for b in bits(mask):
                element_column.append(self.elements[b])
                entry_column.append(i)
        return pd.DataFrame(
            {
                "element": element_column,
                "entry": pd.array(entry_column, dtype="Int64"),
            }
        )

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> GrassmannNecklace[Hashable]:
        """Rebuild a necklace from a frame produced by :meth:`to_dataframe`.

        Tolerates the type drift of a JSON round trip (float entry ids,
        plain ``None`` markers) but not label-dtype drift — read JSON with
        ``dtype=False`` so string labels survive. A frame with no rows
        decodes to the empty necklace. Validates through
        :meth:`from_entries`.

        Args:
            df: Frame with ``element`` and ``entry`` columns, or an empty
                frame for the empty necklace.

        Returns:
            The decoded necklace.

        Raises:
            ValueError: If required columns are missing or the decoded
                entries violate (N1)-(N2).
        """
        if df.empty:
            return GrassmannNecklace.from_entries((), [])
        missing = {"element", "entry"} - set(df.columns)
        if missing:
            msg = f"dataframe is missing required columns {sorted(missing)}"
            raise ValueError(msg)
        ground_rows = df["entry"].isna()
        elements = df.loc[ground_rows, "element"].tolist()
        grouped: dict[int, set[Hashable]] = {}
        membership = df.loc[~ground_rows]
        for label, entry_id in zip(
            membership["element"], membership["entry"], strict=True
        ):
            grouped.setdefault(int(entry_id), set()).add(label)
        entries = [grouped.get(i, set()) for i in range(1, len(elements) + 1)]
        return GrassmannNecklace.from_entries(elements, entries)

    # ------------------------------------------------------- visualization
    def plot_membership_grid(self, ax: Axes | None = None) -> Axes:
        """Draw the entry-membership grid onto ``ax``.

        Row ``i`` (top down) marks the members of ``I_i`` against the
        cyclic order of the columns — the staircase of cyclically shifted
        minimal bases the necklace records (Grassmann Necklace page,
        Overview).

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        n = len(self.elements)
        if n:
            grid = [
                [1.0 if mask >> b & 1 else 0.0 for b in range(n)]
                for mask in self.entry_masks
            ]
            ax.imshow(grid, cmap="Greys", vmin=0.0, vmax=1.25)
            ax.set_xticks(range(n), [repr(e) for e in self.elements], fontsize=8)
            labels = [f"$I_{{{i}}}$" for i in range(1, n + 1)]
            ax.set_yticks(range(n), labels, fontsize=8)
        ax.set_title("Grassmann necklace membership")
        return ax

    def plot_juggling_pattern(self, ax: Axes | None = None) -> Axes:
        """Draw one period of the juggling pattern onto ``ax``.

        Knutson-Lam-Speyer's reading (Grassmann Necklace page, Overview and
        juggling-state blocks): a juggler juggling ``k`` balls, one throw
        every second, doing a pattern of period ``n``. Each throw is an arc
        from second ``t`` to its landing time ``f(t)`` in the bounded
        affine permutation window; empty seconds (loops, ``f(t) = t``) draw
        no arc and coloops arc a full period.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        # Deferred so core use never imports the plotting stack.
        from matplotlib.patches import Arc  # noqa: PLC0415

        n = len(self.elements)
        window = self.to_bounded_affine_permutation()
        for t, landing in enumerate(window, start=1):
            if landing == t:
                continue
            span = float(landing - t)
            ax.add_patch(
                Arc(
                    ((t + landing) / 2, 0.0),
                    span,
                    0.7 * span,
                    theta1=0.0,
                    theta2=180.0,
                    color="0.3",
                )
            )
        points = [(float(t), 0.0) for t in range(1, 2 * n + 1)]
        scatter_labeled(
            ax,
            points,
            [repr(self.elements[(t - 1) % n]) for t in range(1, 2 * n + 1)],
            [(0.0, -10.0)] * (2 * n),
        )
        if n:
            ax.set_ylim(-0.5, 0.35 * n + 0.5)
        ax.set_yticks([])
        siteswap = ",".join(str(digit) for digit in self.siteswap)
        ax.set_title(f"Juggling pattern (siteswap {siteswap})")
        return ax


# --------------------------------------------------------------------------- #
# Canonical examples — test fixtures per the page, exported as API
# --------------------------------------------------------------------------- #
def uniform_necklace(rank: int, n: int) -> GrassmannNecklace[int]:
    """Return the cyclic-interval necklace ``I_i = {i, ..., i+k-1}`` on ``1..n``.

    The necklace of the top cell / uniform positroid ``U_{k,n}`` — the
    unique top element of the circular Bruhat order (Grassmann Necklace
    page, Examples; Postnikov Lemma 17.6, with the page's note that his
    printed interval is a typo for the ``k``-element one). Valid by
    construction, so no validation runs.

    Args:
        rank: The entry cardinality ``k``.
        n: The ground-set size.

    Returns:
        The cyclic-interval necklace on ``(1, ..., n)``.

    Raises:
        ValueError: If ``0 <= rank <= n`` fails.
    """
    if not 0 <= rank <= n:
        msg = f"a necklace of type (k, n) needs 0 <= k <= n, got k={rank}, n={n}"
        raise ValueError(msg)
    masks: list[int] = []
    for start in range(n):
        mask = 0
        for offset in range(rank):
            mask |= 1 << (start + offset) % n
        masks.append(mask)
    return GrassmannNecklace(tuple(range(1, n + 1)), tuple(masks))


def constant_necklace(n: int, entry: Iterable[int]) -> GrassmannNecklace[int]:
    """Return the constant necklace ``(I, I, ..., I)`` on ``1..n``.

    Always a valid necklace (take ``j = i`` in (N1)); it is the necklace of
    the single-basis positroid ``{I}`` (Lam Lemma 8.3), and the constant
    necklaces are exactly the minimal elements of the circular Bruhat order
    (Postnikov Lemma 17.6). No validation runs.

    Args:
        n: The ground-set size.
        entry: The fixed entry ``I``, as labels from ``1..n``.

    Returns:
        The constant necklace on ``(1, ..., n)``.

    Raises:
        ValueError: If an entry label is outside ``1..n``.
    """
    elements = tuple(range(1, n + 1))
    mask = mask_from_labels(entry, {e: i for i, e in enumerate(elements)})
    return GrassmannNecklace(elements, (mask,) * n)


def postnikov_figure_16_1() -> GrassmannNecklace[int]:
    """Return the type-(3,6) necklace of Postnikov's Figure 16.1.

    Postnikov's example: ``pi = (3, 1, 5, 4, 2, 6)`` with ``4`` a black and
    ``6`` a white fixed point, whose necklace is ``({1,2,6}, {2,3,6},
    {1,3,6}, {1,5,6}, {1,5,6}, {1,2,6})`` (Grassmann Necklace page,
    Examples). The page's convention warning applies: this library stores
    the Ardila-Rincon-Williams direction — the inverse ``(2, 5, 1, 4, 3,
    6)`` of Postnikov's permutation — with the white (coloop) fixed point
    colored clockwise.
    """
    decorated = DecoratedPermutation((2, 5, 1, 4, 3, 6), frozenset({6}))
    return GrassmannNecklace.from_decorated_permutation(tuple(range(1, 7)), decorated)


def oh_worked_example() -> GrassmannNecklace[int]:
    """Return Oh's worked example on ``[5]``.

    The necklace ``({1,2,4}, {2,4,5}, {3,4,5}, {2,4,5}, {1,2,5})`` — note
    the repeated entry — whose positroid has the six bases ``{124, 125,
    134, 135, 245, 345}`` (Grassmann Necklace page, Examples; Oh 2011).
    """
    return GrassmannNecklace.from_entries(
        range(1, 6),
        [{1, 2, 4}, {2, 4, 5}, {3, 4, 5}, {2, 4, 5}, {1, 2, 5}],
    )


def lam_example_6_1() -> GrassmannNecklace[int]:
    """Return the necklace of Lam's Example 6.1.

    The bounded affine permutation ``f = [2, 4, 6, 5, 7, 9]`` (``k = 2``,
    ``n = 6``) with necklace ``({1,3}, {2,3}, {3,4}, {4,6}, {5,6}, {1,6})``
    (Grassmann Necklace page, Examples; Lam Theorem 6.2).
    """
    return GrassmannNecklace.from_bounded_affine_permutation(
        range(1, 7), (2, 4, 6, 5, 7, 9)
    )


def kls_example_3_14() -> GrassmannNecklace[int]:
    """Return the juggling pattern of Knutson-Lam-Speyer's Example 3.14.

    For ``n = 4``, ``k = 2``: the bounded affine permutation with window
    ``(2, 3, 5, 8)`` has siteswap ``4112`` up to rotation and
    juggling-state sequence ``({1,4}, {1,3}, {1,2}, {1,2})`` (Grassmann
    Necklace page, Examples).
    """
    return GrassmannNecklace.from_bounded_affine_permutation(range(1, 5), (2, 3, 5, 8))


def enumerate_necklaces(rank: int, n: int) -> list[GrassmannNecklace[int]]:
    """Return every Grassmann necklace of type ``(rank, n)`` on ``1..n``.

    Brute-force search over the transition conditions (N1)-(N2) — the
    direct enumeration behind the page's OEIS A046802 cross-check; the
    counts match rank-``k`` positroids on ``[n]`` by Oh's theorem
    (Grassmann Necklace page, Enumeration). Exponential in ``n``.

    Args:
        rank: The entry cardinality ``k``.
        n: The ground-set size.

    Returns:
        All necklaces of the type, pairwise distinct.

    Raises:
        ValueError: If ``0 <= rank <= n`` fails.
    """
    if not 0 <= rank <= n:
        msg = f"a necklace of type (k, n) needs 0 <= k <= n, got k={rank}, n={n}"
        raise ValueError(msg)
    elements = tuple(range(1, n + 1))
    if n == 0:
        return [GrassmannNecklace(elements, ())]
    full = (1 << n) - 1
    found: list[GrassmannNecklace[int]] = []

    def extend(prefix: list[int]) -> None:
        """Grow the entry sequence by every (N1)/(N2)-legal successor."""
        position = len(prefix) - 1
        current = prefix[position]
        if len(prefix) == n:
            first = prefix[0]
            if current >> position & 1:
                base = current ^ (1 << position)
                closes = not base & ~first and (first & ~base).bit_count() == 1
            else:
                closes = first == current
            if closes:
                found.append(GrassmannNecklace(elements, tuple(prefix)))
            return
        if not current >> position & 1:
            prefix.append(current)
            extend(prefix)
            prefix.pop()
            return
        base = current ^ (1 << position)
        for j in bits(full & ~base):
            prefix.append(base | (1 << j))
            extend(prefix)
            prefix.pop()

    for combo in itertools.combinations(range(n), rank):
        extend([sum(1 << p for p in combo)])
    return found
