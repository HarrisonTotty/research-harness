"""Positroids: the matroids of the totally nonnegative Grassmannian.

A positroid is a matroid on a cyclically ordered ground set representable by
the columns of a full-rank ``d x n`` real matrix all of whose maximal minors
are nonnegative (Postnikov, *Total positivity, Grassmannians, and networks*,
arXiv:math/0609764, 2006, Defs. 3.1-3.2; Ardila-Rincon-Williams, *Positroids
and non-crossing partitions*, arXiv:1308.2698, 2016). This module stores a
positroid as a :class:`research.matroid.Matroid` whose stored element order
carries the cyclic order, validates the positroid property through Oh's
theorem (*Positroids and Schubert matroids*, arXiv:0803.1018, 2011), and
derives the combinatorial indexings the Positroid page defines — Grassmann
necklaces and decorated permutations — as views.

Everything inherits the explicit, exponential design of the matroid module:
practical ground sets stop around 16 elements, and the enumeration helper is
factorial in ``n``.
"""

import functools
import itertools
import math
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, override

import pandas as pd

from research._bitmask import (
    down_closure,
    fmt,
    indexed_ground_set,
    mask_from_labels,
    remap,
)
from research._cyclic import (
    check_necklace_conditions,
    check_positroid,
    gale_geq,
    necklace_bases,
    necklace_masks,
    positroid_witness,
)
from research._graph import UnionFind
from research._linalg import det_q
from research._plot import ensure_axes, scatter_labeled, unit_circle
from research.matroid import Matroid

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "DecoratedPermutation",
    "Positroid",
    "enumerate_positroids",
    "is_positroid",
    "shifted_schubert_positroid",
    "uniform_positroid",
]


# --------------------------------------------------------------------------- #
# Decorated permutations
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DecoratedPermutation:
    """A bijection of ``[n]`` whose fixed points are each colored.

    ``targets[i - 1]`` is the image of position ``i`` in ``1..n``; fixed
    points in ``clockwise_fixed`` are colored clockwise and the remaining
    fixed points counterclockwise (Positroid page, decorated permutation
    block). Under the bijection with positroids, coloops correspond to
    clockwise fixed points and loops to counterclockwise ones — the page
    notes which color counts toward the rank varies by source; this library
    fixes coloop = clockwise.

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
                f"(Positroid page, decorated permutation block); got "
                f"targets {self.targets!r}"
            )
            raise ValueError(msg)
        stray = self.clockwise_fixed - self.fixed_points
        if stray:
            msg = (
                f"only fixed points carry a decoration, but positions "
                f"{sorted(stray)!r} are not fixed by {self.targets!r}"
            )
            raise ValueError(msg)

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
    def weak_excedances(self) -> frozenset[int]:
        """Positions with ``pi(i) > i``, plus the clockwise fixed points.

        Rank-``d`` positroids correspond to decorated permutations with
        exactly ``d`` weak excedances — positions ``i`` with ``pi(i) > i``
        together with fixed points of one designated color (Positroid page);
        the designated color is clockwise under this library's convention.
        """
        strict = frozenset(
            i for i, target in enumerate(self.targets, start=1) if target > i
        )
        return strict | self.clockwise_fixed

    @property
    def weak_excedance_count(self) -> int:
        """The number of weak excedances — the rank of the positroid."""
        return len(self.weak_excedances)

    @property
    def is_stabilized_interval_free(self) -> bool:
        """Whether ``pi(I) != I`` for every proper interval ``I`` of ``[n]``.

        Connected positroids correspond to stabilized-interval-free
        permutations (Ardila-Rincon-Williams Thms. 10.6-10.7; OEIS A075834).
        Quadratic in ``n``.
        """
        n = len(self.targets)
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                if b - a + 1 == n:
                    continue
                if set(self.targets[a - 1 : b]) == set(range(a, b + 1)):
                    return False
        return True


# --------------------------------------------------------------------------- #
# The positroid itself
# --------------------------------------------------------------------------- #
def _adopt[T: Hashable](base: Matroid[T]) -> Positroid[T]:
    """Reread a matroid already known to be a positroid, skipping validation.

    The transformation overrides funnel through here after a closure result
    (Ardila-Rincon-Williams Props. 3.4-3.5, Lemma 3.3) certifies the parent
    computation, so no Oh membership test is run.
    """
    return Positroid(base.elements, base.independent_masks)


@dataclass(frozen=True, eq=False, repr=False)
class Positroid[T: Hashable](Matroid[T]):
    """A positroid: a matroid whose stored element order is its cyclic order.

    The Positroid page's caveat is the design: being a positroid is not
    invariant under relabeling, so the object is a matroid *with* a cyclic
    order — exactly the ``elements`` tuple the parent class already stores.
    Every ``from_<formulation>`` constructor validates the matroid axioms
    through the parent and then the positroid property through Oh's theorem
    membership test (bases must equal ``{B : B >=_j I_j for all j}`` for the
    Grassmann necklace ``I``); calling ``Positroid(...)`` directly skips all
    validation. Equality and hashing stay label-based matroid equality, so
    the cyclic order is representation, not identity — order-dependent views
    (necklace, decorated permutation) read it from ``elements``.

    Transformations the page proves closed — duality, restriction,
    contraction, minors (Ardila-Rincon-Williams Prop. 3.5), cyclic shifts
    (Lemma 3.3), and direct sums of positroids (Prop. 3.4) — return
    ``Positroid`` instances without re-validation, citing those results.
    """

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        return (
            f"Positroid(n={len(self.elements)}, rank={self.rank()}, "
            f"bases={len(self.bases)})"
        )

    # ---------------------------------------------------------- constructors
    @classmethod
    def from_matroid(
        cls, matroid: Matroid[T], *, validate: bool = True
    ) -> Positroid[T]:
        """Adopt a matroid's stored element order as its cyclic order.

        Args:
            matroid: The matroid to reinterpret; its element order becomes
                the cyclic order.
            validate: Check the positroid property via Oh's theorem.

        Returns:
            The positroid on the same ground set and independence family.

        Raises:
            ValueError: If the matroid is not a positroid for this cyclic
                order (the message names Oh's theorem with a witness).
        """
        if validate:
            check_positroid(matroid)
        return cls(matroid.elements, matroid.independent_masks)

    @classmethod
    def from_matrix(cls, columns: Mapping[T, Sequence[Fraction | int]]) -> Positroid[T]:
        """Build a positroid from a rank-``d`` matrix with nonnegative minors.

        The page's primary definition (Ardila-Rincon-Williams phrasing;
        Postnikov Defs. 3.1-3.2): the bases are the ``d``-subsets ``I`` with
        maximal minor ``Delta_I(A) > 0``. Minors are computed exactly over
        the rationals, so no positroid re-validation is needed — the
        definition certifies the result.

        Args:
            columns: Mapping from ground-set label to its column vector; the
                mapping order fixes the cyclic order, and the shared vector
                length is the rank ``d``.

        Returns:
            The positroid of the columns.

        Raises:
            ValueError: If columns have inconsistent dimensions, some
                maximal minor is negative, or the matrix is not full rank.
        """
        labels = tuple(columns)
        dimensions = {len(v) for v in columns.values()}
        if len(dimensions) > 1:
            msg = f"columns must share one dimension, got lengths {dimensions}"
            raise ValueError(msg)
        d = dimensions.pop() if dimensions else 0
        rational = [[Fraction(a) for a in columns[label]] for label in labels]
        n = len(labels)
        basis_masks: list[int] = []
        for combo in itertools.combinations(range(n), d):
            minor = det_q([[rational[c][r] for c in combo] for r in range(d)])
            if minor < 0:
                chosen = sum(1 << c for c in combo)
                msg = (
                    f"maximal minor for {fmt(chosen, labels)} is {minor} < 0; "
                    f"a positroid representation needs every maximal minor "
                    f"nonnegative (Positroid page, Definition)"
                )
                raise ValueError(msg)
            if minor > 0:
                basis_masks.append(sum(1 << c for c in combo))
        if not basis_masks:
            msg = (
                f"the matrix must have full rank {d} (Positroid page, "
                f"Definition), but every maximal minor vanishes"
            )
            raise ValueError(msg)
        return cls(labels, down_closure(basis_masks))

    @classmethod
    def _from_necklace_masks(
        cls, elements: tuple[T, ...], masks: tuple[int, ...]
    ) -> Positroid[T]:
        """Apply Oh's construction ``B(I) = {B : B >=_j I_j for all j}``.

        Oh's theorem certifies the result is a positroid with necklace
        ``I``, so no re-validation is performed.
        """
        return cls(elements, down_closure(necklace_bases(masks, len(elements))))

    @classmethod
    def from_grassmann_necklace(
        cls,
        elements: Iterable[T],
        necklace: Sequence[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build the positroid of a Grassmann necklace via Oh's theorem.

        The bases are ``{B : B >=_j I_j for all j}`` in the shifted Gale
        orders (Positroid page, Oh's theorem block — equivalently, the
        intersection of the cyclically shifted Schubert matroids cut out by
        the necklace entries).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            necklace: The entries ``I_1, ..., I_n``, one per ground-set
                element, as label collections.
            validate: Check the necklace conditions of Postnikov section 16.

        Returns:
            The positroid ``M(I)``.

        Raises:
            ValueError: If the necklace length or labels are wrong, or a
                necklace condition fails (the message names it).
        """
        elems, index = indexed_ground_set(elements)
        masks = tuple(mask_from_labels(entry, index) for entry in necklace)
        if len(masks) != len(elems):
            msg = (
                f"a Grassmann necklace has one entry per ground-set element "
                f"(Postnikov section 16); got {len(masks)} entries for "
                f"{len(elems)} elements"
            )
            raise ValueError(msg)
        if validate:
            check_necklace_conditions(elems, masks)
        return cls._from_necklace_masks(elems, masks)

    @classmethod
    def from_decorated_permutation(
        cls,
        elements: Iterable[T],
        decorated: DecoratedPermutation,
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build the positroid of a decorated permutation.

        The Grassmann necklace is recovered from the page's transition rule
        read backwards: the element added passing from ``I_i`` to
        ``I_{i+1}`` is the ``j`` with ``pi(j) = i``, so ``j`` belongs to
        ``I_k`` exactly when ``k`` lies in the cyclic interval
        ``(pi(j), j]``; clockwise fixed points (coloops) lie in every entry
        and counterclockwise ones (loops) in none (Positroid page, Grassmann
        necklace and decorated permutation blocks). Oh's construction then
        yields the positroid.

        Args:
            elements: Ground-set labels; their order is the cyclic order,
                and position ``i`` of the permutation is ``elements[i-1]``.
            decorated: The decorated permutation of ``[n]``.
            validate: Re-check the derived necklace conditions.

        Returns:
            The positroid of the decorated permutation.

        Raises:
            ValueError: If the ground-set size does not match the
                permutation, or the derived necklace is inconsistent.
        """
        elems, _ = indexed_ground_set(elements)
        n = len(elems)
        if len(decorated.targets) != n:
            msg = (
                f"the decorated permutation acts on [{len(decorated.targets)}] "
                f"but the ground set has {n} elements"
            )
            raise ValueError(msg)
        masks = [0] * n
        for j in range(n):
            target = decorated.targets[j] - 1
            if target == j:
                if j + 1 in decorated.clockwise_fixed:
                    for k in range(n):
                        masks[k] |= 1 << j
                continue
            span = (j - target - 1) % n
            for k in range(n):
                if (k - target - 1) % n <= span:
                    masks[k] |= 1 << j
        if validate:
            check_necklace_conditions(elems, tuple(masks))
        return cls._from_necklace_masks(elems, tuple(masks))

    # ------------------------------------------- inherited-formulation gates
    # Every constructor inherited from Matroid is re-exposed so that no
    # inherited path can mint an unvalidated Positroid: each builds the
    # matroid through the parent, then applies the Oh membership test.
    @classmethod
    @override
    def from_independent_sets(
        cls,
        elements: Iterable[T],
        independent_sets: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from its independent sets.

        Matroid axioms (I1)-(I3) are checked by the parent constructor, then
        the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            independent_sets: Every independent set, as label collections.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_independent_sets(
            elements, independent_sets, validate=validate
        )
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_bases(
        cls,
        elements: Iterable[T],
        bases: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from its bases.

        Matroid axioms (B1)-(B2) are checked by the parent constructor, then
        the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            bases: The maximal independent sets.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_bases(elements, bases, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_circuits(
        cls,
        elements: Iterable[T],
        circuits: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from its circuits.

        Matroid axioms (C1)-(C3) are checked by the parent constructor, then
        the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            circuits: The minimal dependent sets.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_circuits(elements, circuits, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_rank_function(
        cls,
        elements: Iterable[T],
        rank: Callable[[frozenset[T]], int],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from a rank oracle.

        Rank axioms (R1)-(R3) are checked by the parent constructor, then
        the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            rank: Callable evaluated on every subset of the ground set.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_rank_function(elements, rank, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_closure(
        cls,
        elements: Iterable[T],
        closure: Callable[[frozenset[T]], frozenset[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from a closure oracle.

        Closure axioms (CL1)-(CL4) are checked by the parent constructor,
        then the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            closure: Callable evaluated on every subset of the ground set.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_closure(elements, closure, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_flats(
        cls,
        elements: Iterable[T],
        flats: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from its flats.

        Flat axioms (F1)-(F3) are checked by the parent constructor, then
        the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            flats: The closed sets, as label collections.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_flats(elements, flats, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_hyperplanes(
        cls,
        elements: Iterable[T],
        hyperplanes: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build a positroid from its hyperplanes.

        Hyperplane axioms (H1)-(H3) are checked by the parent constructor,
        then the positroid property via Oh's theorem (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            hyperplanes: The rank ``r(M) - 1`` flats, as label collections.
            validate: Check the matroid axioms and the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If an axiom or the positroid property fails (the
                message names the axiom or Oh's theorem).
        """
        matroid = super().from_hyperplanes(elements, hyperplanes, validate=validate)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_vectors(
        cls,
        vectors: Mapping[T, Sequence[Fraction | int]],
        *,
        field_char: int | None = None,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build the positroid of a linear matroid, if it is one.

        The parent constructor builds the linear matroid of the columns;
        the positroid property is then checked via Oh's theorem for the
        mapping order. Note the asymmetry with :meth:`from_matrix`: a
        matrix with some negative minors can still represent a positroid,
        because only the existence of *some* nonnegative representation
        matters (Positroid page, Definition).

        Args:
            vectors: Mapping from ground-set label to coordinate vector;
                the mapping order fixes the cyclic order.
            field_char: ``None`` for the rationals, or a prime ``p``.
            validate: Check the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If the input shape is wrong or the linear matroid
                is not a positroid (the message names Oh's theorem).
        """
        matroid = super().from_vectors(vectors, field_char=field_char)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_graph_edges[V: Hashable](
        cls, edges: Mapping[T, tuple[V, V]], *, validate: bool = True
    ) -> Positroid[T]:
        """Build the positroid of a graphic matroid, if it is one.

        The parent constructor builds the cycle matroid of the edges; the
        positroid property is then checked via Oh's theorem for the mapping
        order (Positroid page).

        Args:
            edges: Mapping from edge label to its endpoint pair; the
                mapping order fixes the cyclic order.
            validate: Check the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If the graphic matroid is not a positroid for this
                cyclic order (the message names Oh's theorem).
        """
        matroid = super().from_graph_edges(edges)
        return cls.from_matroid(matroid, validate=validate)

    @classmethod
    @override
    def from_transversal_system(
        cls,
        elements: Iterable[T],
        system: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Positroid[T]:
        """Build the positroid of a transversal matroid, if it is one.

        The parent constructor builds the transversal matroid of the
        system; the positroid property is then checked via Oh's theorem for
        the element order (Positroid page).

        Args:
            elements: Ground-set labels; their order is the cyclic order.
            system: The sets of the system, as label collections.
            validate: Check the positroid property.

        Returns:
            The validated positroid.

        Raises:
            ValueError: If the transversal matroid is not a positroid for
                this cyclic order (the message names Oh's theorem).
        """
        matroid = super().from_transversal_system(elements, system)
        return cls.from_matroid(matroid, validate=validate)

    # ------------------------------------------------- cryptomorphic views
    @functools.cached_property
    def _necklace_position_masks(self) -> tuple[int, ...]:
        """The Grassmann necklace as bitmasks, one entry per start position."""
        return necklace_masks(len(self.elements), self.independent_masks)

    @functools.cached_property
    def grassmann_necklace(self) -> tuple[frozenset[T], ...]:
        """The Grassmann necklace ``(I_1, ..., I_n)`` of the positroid.

        Entry ``i`` is the Gale-minimal (equivalently lexicographically
        minimal) basis with respect to the cyclically shifted order starting
        at ``elements[i-1]`` (Positroid page, Grassmann necklace block).
        Computed once and cached; ``n`` greedy passes.
        """
        return tuple(self._labels(m) for m in self._necklace_position_masks)

    def to_decorated_permutation(self) -> DecoratedPermutation:
        r"""Return the decorated permutation indexing this positroid.

        Read off the necklace transitions: when ``I_{i+1}`` is
        ``(I_i \ {i}) + {j}`` the permutation sends ``j`` to ``i``; a
        coloop is a clockwise fixed point and a loop a counterclockwise one
        (Positroid page, decorated permutation block; the color convention
        is this library's, as the page leaves it open). Its weak excedance
        count is the rank.

        Returns:
            The decorated permutation on positions ``1..n``.

        Raises:
            ValueError: If the necklace transitions are inconsistent, which
                indicates the instance was built without validation.
        """
        n = len(self.elements)
        masks = self._necklace_position_masks
        targets: list[int] = [0] * n
        clockwise: set[int] = set()
        for i in range(n):
            current = masks[i]
            following = masks[(i + 1) % n]
            if not current >> i & 1:
                targets[i] = i + 1
                continue
            base = current ^ (1 << i)
            added = following & ~base
            if added.bit_count() != 1 or base & ~following:
                msg = (
                    f"inconsistent Grassmann necklace at position {i + 1}; "
                    f"was this instance built without validation?"
                )
                raise ValueError(msg)
            j = added.bit_length() - 1
            if j == i:
                targets[i] = i + 1
                clockwise.add(i + 1)
            else:
                targets[j] = i + 1
        return DecoratedPermutation(tuple(targets), frozenset(clockwise))

    # ------------------------------------------------- computed properties
    def connected_components(self) -> frozenset[frozenset[T]]:
        """Return the connected components of the underlying matroid.

        Two elements are in one component when they share a circuit
        (Matroid page, derived vocabulary); for a positroid the components
        form a non-crossing partition of the cyclic order
        (Ardila-Rincon-Williams Thm. 7.6). Elements in no circuit (coloops)
        are singleton components.
        """
        components = UnionFind[int]()
        for position in range(len(self.elements)):
            components.add(position)
        for circuit in self.circuits:
            members = [self._index[e] for e in circuit]
            for member in members[1:]:
                components.union(members[0], member)
        return frozenset(
            frozenset(self.elements[p] for p in group) for group in components.groups()
        )

    def cyclic_rank_bounds(self) -> dict[tuple[int, int], int]:
        """Return ``{(i, j): r([i, j])}`` over all cyclic intervals of ``[n]``.

        These are the tightest right-hand sides for the cyclic-interval
        inequalities cutting out the positroid polytope
        (Ardila-Rincon-Williams Prop. 5.6): a matroid is a positroid iff its
        matroid polytope is cut out by the equality ``sum(x) = d`` together
        with inequalities over cyclic intervals. Keys are 1-based positions
        in the cyclic order; ``[i, j]`` wraps when ``j < i``.
        """
        n = len(self.elements)
        bounds: dict[tuple[int, int], int] = {}
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                length = (j - i) % n + 1
                mask = 0
                for offset in range(length):
                    mask |= 1 << (i - 1 + offset) % n
                bounds[(i, j)] = self._rank_table[mask]
        return bounds

    # ------------------------------------------------------- transformations
    def to_matroid(self) -> Matroid[T]:
        """Return the underlying matroid, forgetting the positroid reading."""
        return Matroid(self.elements, self.independent_masks)

    def with_cyclic_order(
        self, order: Iterable[T], *, validate: bool = True
    ) -> Positroid[T]:
        """Return the same labeled matroid read with a new cyclic order.

        This is the page's order-dependence caveat made executable: the
        positroid property depends on the chosen cyclic order, so an
        arbitrary reordering can fail validation while every rotation
        succeeds (Ardila-Rincon-Williams Lemma 3.3).

        Args:
            order: The ground-set labels in the new cyclic order.
            validate: Check the positroid property for the new order.

        Returns:
            The positroid with the new stored order.

        Raises:
            ValueError: If ``order`` is not a permutation of the ground set
                or the matroid is not a positroid for it.
        """
        new_elements, translate = indexed_ground_set(order)
        if frozenset(new_elements) != self.ground_set:
            msg = (
                f"the new order must be a permutation of the ground set; "
                f"got {new_elements!r}"
            )
            raise ValueError(msg)
        table = [translate[e] for e in self.elements]
        family = frozenset(remap(mask, table) for mask in self.independent_masks)
        result = Positroid(new_elements, family)
        if validate:
            check_positroid(result)
        return result

    def cyclic_shift(self, steps: int = 1) -> Positroid[T]:
        """Return the positroid with its cyclic order rotated by ``steps``.

        Positroids are closed under cyclic shift (Ardila-Rincon-Williams
        Lemma 3.3), so no re-validation is performed.
        """
        n = len(self.elements)
        k = steps % n if n else 0
        return self.with_cyclic_order(
            self.elements[k:] + self.elements[:k], validate=False
        )

    @override
    def dual(self) -> Positroid[T]:
        """Return ``M*``, a positroid for the same cyclic order.

        Closure under duality is Ardila-Rincon-Williams Prop. 3.5, so no
        re-validation is performed.
        """
        return _adopt(super().dual())

    @override
    def restrict(self, subset: Iterable[T]) -> Positroid[T]:
        """Return ``M|X``, a positroid for the inherited cyclic order.

        Closure under restriction is Ardila-Rincon-Williams Prop. 3.5, so
        no re-validation is performed.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return _adopt(super().restrict(subset))

    @override
    def delete(self, subset: Iterable[T]) -> Positroid[T]:
        """Return ``M - X``, a positroid for the inherited cyclic order.

        Deletion is restriction to the complement, so closure follows from
        Ardila-Rincon-Williams Prop. 3.5.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return _adopt(super().delete(subset))

    @override
    def contract(self, subset: Iterable[T]) -> Positroid[T]:
        """Return ``M/X``, a positroid for the inherited cyclic order.

        Closure under contraction is Ardila-Rincon-Williams Prop. 3.5, so
        no re-validation is performed.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return _adopt(super().contract(subset))

    @override
    def minor(
        self,
        *,
        deletions: Iterable[T] = (),
        contractions: Iterable[T] = (),
    ) -> Positroid[T]:
        """Return the minor ``M - D / C``, itself a positroid.

        Positroids are closed under minors (Ardila-Rincon-Williams
        Prop. 3.5), so no re-validation is performed.

        Args:
            deletions: Elements to delete (``D``).
            contractions: Elements to contract (``C``); disjoint from ``D``.

        Returns:
            The minor on the remaining elements, in inherited cyclic order.

        Raises:
            ValueError: If the two sets overlap or contain unknown labels.
        """
        return _adopt(super().minor(deletions=deletions, contractions=contractions))

    @override
    def simplification(self) -> Positroid[T]:
        """Return ``si(M)``, a positroid as a restriction of ``M``.

        Simplification restricts to loop-free parallel-class
        representatives, so closure follows from Ardila-Rincon-Williams
        Prop. 3.5.
        """
        return _adopt(super().simplification())

    @override
    def direct_sum[U: Hashable](self, other: Matroid[U]) -> Matroid[T | U]:
        """Return the direct sum, a positroid when ``other`` is one.

        Concatenating the two cyclic orders places the ground sets as
        cyclic intervals, so the sum of two positroids is a positroid
        (Ardila-Rincon-Williams Prop. 3.4) and is returned as a
        ``Positroid`` instance without re-validation; summing with a plain
        matroid falls back to the parent behavior.

        Raises:
            ValueError: If the ground sets share a label.
        """
        combined = super().direct_sum(other)
        if isinstance(other, Positroid):
            return _adopt(combined)
        return combined

    # ---------------------------------------------------------- dataframes
    @staticmethod
    @override
    def from_dataframe(df: pd.DataFrame) -> Positroid[Hashable]:
        """Rebuild a positroid from a frame produced by ``to_dataframe``.

        The inherited encoding stores the ground rows in cyclic order, so
        the order survives the round trip; decoding validates both the
        basis axioms (parent) and the positroid property (Oh's theorem).

        Args:
            df: Frame with ``element`` and ``basis`` columns, or an empty
                frame for the empty positroid.

        Returns:
            The decoded positroid.

        Raises:
            ValueError: If required columns are missing or the decoded data
                is not a positroid for the stored order.
        """
        return Positroid.from_matroid(Matroid.from_dataframe(df))

    # ------------------------------------------------------- visualization
    def _circle_positions(self) -> list[tuple[float, float]]:
        """Place the ground set clockwise on the unit circle, first on top."""
        return unit_circle(len(self.elements), phase=math.pi / 2, clockwise=True)

    def _annotate_circle(self, ax: Axes, points: list[tuple[float, float]]) -> None:
        """Draw and label the ground-set points of a circular diagram."""
        scatter_labeled(
            ax,
            points,
            [repr(element) for element in self.elements],
            [(12 * x, 12 * y) for x, y in points],
        )
        ax.set_aspect("equal")
        ax.set_axis_off()

    def plot_decorated_permutation(self, ax: Axes | None = None) -> Axes:
        """Draw the decorated permutation as a chord diagram onto ``ax``.

        Ground-set elements sit clockwise on a circle in cyclic order;
        arrows join position ``i`` to ``pi(i)``, clockwise fixed points
        (coloops) are filled markers and counterclockwise ones (loops)
        hollow — the page's clockwise/counterclockwise vocabulary drawn
        literally (Positroid page, decorated permutation block).

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        decorated = self.to_decorated_permutation()
        points = self._circle_positions()
        for i, target in enumerate(decorated.targets, start=1):
            if target == i:
                continue
            ax.annotate(
                "",
                xy=points[target - 1],
                xytext=points[i - 1],
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": "0.3",
                    "shrinkA": 8,
                    "shrinkB": 8,
                    "connectionstyle": "arc3,rad=0.2",
                },
            )
        for i in sorted(decorated.fixed_points):
            x, y = points[i - 1]
            filled = i in decorated.clockwise_fixed
            ax.scatter(
                [x],
                [y],
                s=120,
                facecolors="0.2" if filled else "none",
                edgecolors="0.2",
                zorder=3,
            )
        self._annotate_circle(ax, points)
        ax.set_title("Decorated permutation")
        return ax

    def plot_connected_components(self, ax: Axes | None = None) -> Axes:
        """Draw the non-crossing partition of components onto ``ax``.

        Ground-set elements sit clockwise on a circle in cyclic order and
        each connected component is joined into a closed chord loop; by
        Ardila-Rincon-Williams Thm. 7.6 the loops never cross.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        points = self._circle_positions()
        components = sorted(
            (
                sorted(self._index[e] for e in block)
                for block in self.connected_components()
            ),
            key=lambda block: block[0],
        )
        for block in components:
            if len(block) == 1:
                continue
            cycle = [*block, block[0]]
            ax.plot(
                [points[p][0] for p in cycle],
                [points[p][1] for p in cycle],
                linewidth=1.5,
                zorder=1,
            )
        self._annotate_circle(ax, points)
        ax.set_title("Connected components (non-crossing partition)")
        return ax


# --------------------------------------------------------------------------- #
# Canonical examples — test fixtures per the Positroid page, exported as API
# --------------------------------------------------------------------------- #
def uniform_positroid(rank: int, n: int) -> Positroid[int]:
    """Return ``U_{d,n}`` on ``range(n)`` as the positroid of the top cell.

    The uniform matroid with the standard cyclic order is the positroid of
    the totally positive Grassmannian's top cell — any totally positive
    matrix realizes it (Positroid page, Examples), so no validation is
    needed.

    Args:
        rank: The rank ``d``.
        n: The ground-set size.

    Returns:
        The uniform positroid.

    Raises:
        ValueError: If ``0 <= rank <= n`` fails.
    """
    if not 0 <= rank <= n:
        msg = f"uniform positroid needs 0 <= rank <= n, got rank={rank}, n={n}"
        raise ValueError(msg)
    family = frozenset(m for m in range(1 << n) if m.bit_count() <= rank)
    return Positroid(tuple(range(n)), family)


def shifted_schubert_positroid[T: Hashable](
    elements: Iterable[T],
    core: Iterable[T],
    position: int = 1,
    *,
    validate: bool = True,
) -> Positroid[T]:
    """Return the cyclically shifted Schubert matroid ``{B : B >=_j core}``.

    Each single Gale-order condition of Oh's theorem cuts out one shifted
    Schubert matroid, and positroids are exactly the intersections of these
    building blocks (Positroid page, Oh's theorem and Examples blocks).

    Args:
        elements: Ground-set labels; their order is the cyclic order.
        core: The ``d``-subset every basis must dominate in Gale order.
        position: The 1-based shift ``j`` starting the cyclic order.
        validate: Check the matroid axioms and the positroid property.

    Returns:
        The shifted Schubert matroid as a positroid.

    Raises:
        ValueError: If labels are unknown or duplicated, or ``position`` is
            out of range.
    """
    elems, index = indexed_ground_set(elements)
    n = len(elems)
    if not 1 <= position <= n:
        msg = f"position must be in 1..{n}, got {position}"
        raise ValueError(msg)
    core_mask = mask_from_labels(core, index)
    d = core_mask.bit_count()
    bases = [
        [elems[p] for p in combo]
        for combo in itertools.combinations(range(n), d)
        if gale_geq(sum(1 << p for p in combo), core_mask, position - 1, n)
    ]
    return Positroid.from_bases(elems, bases, validate=validate)


def enumerate_positroids(n: int) -> list[Positroid[int]]:
    """Return every positroid on the cyclically ordered ground set ``1..n``.

    Positroids on ``[n]`` biject with decorated permutations of ``[n]``, so
    the enumeration walks every permutation with every coloring of its
    fixed points; the count is ``sum(n!/k!)`` — 1, 2, 5, 16, 65, 326, ...
    (OEIS A000522; Positroid page, Enumeration). Factorial in ``n``.

    Args:
        n: Ground-set size.

    Returns:
        All positroids on ``(1, ..., n)``, pairwise distinct.

    Raises:
        ValueError: If ``n`` is negative.
    """
    if n < 0:
        msg = f"ground-set size must be non-negative, got {n}"
        raise ValueError(msg)
    elements = tuple(range(1, n + 1))
    found: list[Positroid[int]] = []
    for targets in itertools.permutations(elements):
        fixed = [i for i in elements if targets[i - 1] == i]
        for count in range(len(fixed) + 1):
            for clockwise in itertools.combinations(fixed, count):
                decorated = DecoratedPermutation(targets, frozenset(clockwise))
                found.append(
                    Positroid.from_decorated_permutation(
                        elements, decorated, validate=False
                    )
                )
    return found


def is_positroid[T: Hashable](matroid: Matroid[T]) -> bool:
    """Return whether a matroid is a positroid for its stored element order.

    Runs Oh's theorem membership test: the bases must be exactly
    ``{B : B >=_j I_j for all j}`` for the matroid's Grassmann necklace
    ``I`` (Positroid page, Oh's theorem block — "a finite, efficiently
    checkable membership test"). Remember the page's caveat: the answer
    depends on the stored order, not just the isomorphism class.
    """
    return positroid_witness(matroid) is None
