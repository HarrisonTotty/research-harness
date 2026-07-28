"""Finite matroids over explicit ground sets, following the Matroid page.

A matroid abstracts the notion of independence shared by matrix columns,
forests in a graph, and partial transversals (Whitney, *On the abstract
properties of linear dependence*, 1935; independently Nakasawa 1935-36).
This module stores exactly one primitive — the independence family, as
``int`` bitmasks over a fixed element ordering — and derives every other
standard presentation (bases, circuits, rank, closure, flats, hyperplanes)
lazily from it, per the cryptomorphism conversion table on the Matroid page.

Everything here is explicit and exponential by design: the practical
envelope is roughly ``n <= 16`` ground-set elements, and the classification
and enumeration helpers say so where they are worse. Axiom validation is on
by default in every ``from_<formulation>`` constructor and can be switched
off for large trusted inputs with ``validate=False``.
"""

import functools
import itertools
import math
from collections import Counter
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, override

import pandas as pd

from research._axioms import (
    basis_exchange_violation,
    check_basis_axioms,
    check_circuit_axioms,
    check_closure_axioms,
    check_flat_axioms,
    check_hyperplane_axioms,
    check_independence_axioms,
    check_rank_axioms,
)
from research._bitmask import (
    bits,
    down_closure,
    fmt,
    indexed_ground_set,
    mask_from_labels,
    remap,
    submasks,
)
from research._graph import has_matching, is_forest
from research._linalg import (
    RATIONALS,
    gf_scalar,
    is_prime,
    linear_independence_family,
    prime_field,
)
from research._plot import ensure_axes, scatter_labeled, unit_circle

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "Matroid",
    "empty_matroid",
    "enumerate_matroids",
    "fano_matroid",
    "free_matroid",
    "k4_matroid",
    "loopy_matroid",
    "non_fano_matroid",
    "u24",
    "uniform_matroid",
    "vamos_matroid",
]

_EXCHANGE_DISTANCE: int = 2
"""Symmetric-difference size linking two bases in the basis exchange graph."""


# --------------------------------------------------------------------------- #
# The matroid itself
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, eq=False)
class Matroid[T: Hashable]:
    """A finite matroid ``(E, I)`` stored as an explicit independence family.

    The ground set is ``elements`` (order fixes the label-to-bit bijection)
    and ``independent_masks`` holds every independent set as a bitmask. The
    ground set is stored explicitly because it is not recoverable from the
    independence family alone — a loop lies in no independent set (Matroid
    page, Definition — Independent Sets).

    Calling ``Matroid(...)`` directly skips all axiom validation; the
    supported way to build one is through the ``from_<formulation>``
    classmethods, which validate their axiom system by default. Equality and
    hashing compare the *labeled* matroid — ground-set labels plus the
    independence family as label sets — so element order is representation,
    not identity.

    Derived presentations are memoized on first use; the rank table costs
    ``O(2^n * n)`` time and ``O(2^n)`` memory, which bounds practical ground
    sets at roughly 16 elements.
    """

    elements: tuple[T, ...]
    independent_masks: frozenset[int]

    # ---------------------------------------------------------------- basics
    @functools.cached_property
    def _index(self) -> dict[T, int]:
        """Map each ground-set label to its bit position."""
        return {element: i for i, element in enumerate(self.elements)}

    @functools.cached_property
    def _full_mask(self) -> int:
        """Bitmask of the whole ground set."""
        return (1 << len(self.elements)) - 1

    def _mask(self, labels: Iterable[T]) -> int:
        """Convert labels to a bitmask, rejecting unknown labels."""
        return mask_from_labels(labels, self._index)

    def _labels(self, mask: int) -> frozenset[T]:
        """Convert a bitmask back to a frozenset of labels."""
        return frozenset(self.elements[b] for b in bits(mask))

    @property
    def ground_set(self) -> frozenset[T]:
        """The ground set ``E`` as a frozenset of labels."""
        return frozenset(self.elements)

    @property
    def size(self) -> int:
        """The number of ground-set elements ``|E|``."""
        return len(self.elements)

    @functools.cached_property
    def _canonical_family(self) -> frozenset[frozenset[T]]:
        """The independence family as label sets (the equality witness)."""
        return frozenset(self._labels(mask) for mask in self.independent_masks)

    @override
    def __eq__(self, other: object) -> bool:
        """Compare labeled matroids: same ground set, same independent sets."""
        if not isinstance(other, Matroid):
            return NotImplemented
        return (
            self.ground_set == other.ground_set
            and self._canonical_family == other._canonical_family
        )

    @override
    def __hash__(self) -> int:
        """Hash the labeled canonical form used by ``__eq__``."""
        return hash((self.ground_set, self._canonical_family))

    @override
    def __repr__(self) -> str:
        """Render compactly enough to read in a failing test."""
        return (
            f"Matroid(n={len(self.elements)}, rank={self.rank()}, "
            f"bases={len(self._basis_masks)})"
        )

    # ---------------------------------------------------------- constructors
    @classmethod
    def from_independent_sets(
        cls,
        elements: Iterable[T],
        independent_sets: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from its independent sets, checking (I1)-(I3).

        This is the primary formulation on the Matroid page (Whitney 1935,
        axioms (a)/(b), stated for general cardinalities).

        Args:
            elements: Ground-set labels in the order that fixes bit positions.
            independent_sets: Every independent set, as label collections.
            validate: Check the axioms (default). Naive (I3) checking is
                ``O(|I|^2 * n)``; disable only for large trusted input.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are unknown or duplicated, or an axiom
                fails (the message names the numbered axiom).
        """
        elems, index = indexed_ground_set(elements)
        masks = frozenset(mask_from_labels(s, index) for s in independent_sets)
        if validate:
            check_independence_axioms(elems, masks)
        return cls(elems, masks)

    @classmethod
    def from_bases(
        cls,
        elements: Iterable[T],
        bases: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from its bases, checking (B1)-(B2).

        The independence family is the downward closure of the bases
        (cryptomorphism table, B -> I). Equicardinality is a theorem and is
        not assumed of the input.

        Args:
            elements: Ground-set labels in order.
            bases: The maximal independent sets.
            validate: Check (B1) and basis exchange (B2) on the input.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are unknown or duplicated, or an axiom
                fails (the message names the numbered axiom).
        """
        elems, index = indexed_ground_set(elements)
        base_masks = frozenset(mask_from_labels(b, index) for b in bases)
        if validate:
            check_basis_axioms(elems, base_masks)
        return cls(elems, down_closure(base_masks))

    @classmethod
    def from_circuits(
        cls,
        elements: Iterable[T],
        circuits: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from its circuits, checking (C1)-(C3).

        A set is independent iff it contains no circuit (cryptomorphism
        table, C -> I); the dependent sets are propagated upward in
        ``O(2^n * n)``.

        Args:
            elements: Ground-set labels in order.
            circuits: The minimal dependent sets.
            validate: Check (C1)-(C3), including circuit elimination.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are unknown or duplicated, or an axiom
                fails (the message names the numbered axiom).
        """
        elems, index = indexed_ground_set(elements)
        circuit_masks = frozenset(mask_from_labels(c, index) for c in circuits)
        if validate:
            check_circuit_axioms(elems, circuit_masks)
        n = len(elems)
        dependent = bytearray(1 << n)
        for mask in range(1 << n):
            if mask in circuit_masks:
                dependent[mask] = 1
                continue
            for bit in bits(mask):
                if dependent[mask ^ (1 << bit)]:
                    dependent[mask] = 1
                    break
        family = frozenset(m for m in range(1 << n) if not dependent[m])
        return cls(elems, family)

    @classmethod
    def from_rank_function(
        cls,
        elements: Iterable[T],
        rank: Callable[[frozenset[T]], int],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from a rank oracle, checking (R1)-(R3).

        Validation uses Whitney's local axioms (1935, section 2): rank of the
        empty set, unit increase, and local flatness — ``O(2^n * n^2)``
        rather than ``O(4^n)`` pairwise submodularity, which Whitney derives
        as his Theorem 3. A set is independent iff ``r(X) == |X|``
        (cryptomorphism table, r -> I).

        Args:
            elements: Ground-set labels in order.
            rank: Callable evaluated on every subset of the ground set.
            validate: Check the axioms on the tabulated oracle.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are duplicated or an axiom fails (the
                message names the numbered axiom).
        """
        elems, _ = indexed_ground_set(elements)
        n = len(elems)
        table = [
            rank(frozenset(elems[b] for b in bits(mask))) for mask in range(1 << n)
        ]
        if validate:
            check_rank_axioms(elems, table)
        family = frozenset(m for m in range(1 << n) if table[m] == m.bit_count())
        return cls(elems, family)

    @classmethod
    def from_closure(
        cls,
        elements: Iterable[T],
        closure: Callable[[frozenset[T]], frozenset[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from a closure oracle, checking (CL1)-(CL4).

        (CL1)-(CL3) define a general closure operator; (CL4) is the Mac
        Lane-Steinitz exchange property that makes it a matroid closure. A
        set is independent iff no element is in the closure of the others
        (cryptomorphism table, cl -> I).

        Args:
            elements: Ground-set labels in order.
            closure: Callable evaluated on every subset; must return subsets
                of the ground set.
            validate: Check the axioms on the tabulated oracle.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are duplicated, the oracle returns unknown
                labels, or an axiom fails (the message names the axiom).
        """
        elems, index = indexed_ground_set(elements)
        n = len(elems)
        table = [
            mask_from_labels(closure(frozenset(elems[b] for b in bits(mask))), index)
            for mask in range(1 << n)
        ]
        if validate:
            check_closure_axioms(elems, table)
        family = frozenset(
            mask
            for mask in range(1 << n)
            if all(not table[mask ^ (1 << e)] >> e & 1 for e in bits(mask))
        )
        return cls(elems, family)

    @classmethod
    def from_flats(
        cls,
        elements: Iterable[T],
        flats: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from its flats, checking (F1)-(F3).

        The closure of a set is the intersection of all flats containing it
        (cryptomorphism table, F -> cl), and independence follows as in
        :meth:`from_closure`.

        Args:
            elements: Ground-set labels in order.
            flats: The closed sets, as label collections.
            validate: Check (F1), (F2), and the covering axiom (F3).

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are unknown or duplicated, or an axiom
                fails (the message names the numbered axiom).
        """
        elems, index = indexed_ground_set(elements)
        flat_masks = frozenset(mask_from_labels(f, index) for f in flats)
        if validate:
            check_flat_axioms(elems, flat_masks)
        n = len(elems)
        full = (1 << n) - 1
        table = []
        for mask in range(1 << n):
            closure_mask = full
            for flat in flat_masks:
                if mask & ~flat == 0:
                    closure_mask &= flat
            table.append(closure_mask)
        family = frozenset(
            mask
            for mask in range(1 << n)
            if all(not table[mask ^ (1 << e)] >> e & 1 for e in bits(mask))
        )
        return cls(elems, family)

    @classmethod
    def from_hyperplanes(
        cls,
        elements: Iterable[T],
        hyperplanes: Iterable[Iterable[T]],
        *,
        validate: bool = True,
    ) -> Matroid[T]:
        """Build a matroid from its hyperplanes, checking (H1)-(H3).

        Hyperplane complements are exactly the cocircuits (Matroid page,
        Definition — Hyperplanes), so the dual is built from those circuits
        and dualized back.

        Args:
            elements: Ground-set labels in order.
            hyperplanes: The rank ``r(M) - 1`` flats, as label collections.
            validate: Check (H1)-(H3) on the input.

        Returns:
            The validated matroid.

        Raises:
            ValueError: If labels are unknown or duplicated, or an axiom
                fails (the message names the numbered axiom).
        """
        elems, index = indexed_ground_set(elements)
        h_masks = frozenset(mask_from_labels(h, index) for h in hyperplanes)
        if validate:
            check_hyperplane_axioms(elems, h_masks)
        full = (1 << len(elems)) - 1
        cocircuits = [[elems[b] for b in bits(full ^ h_mask)] for h_mask in h_masks]
        # (H1)-(H3) already certify the system; the complements satisfy the
        # circuit axioms by cryptomorphism, so re-validation is skipped.
        dual = cls.from_circuits(elems, cocircuits, validate=False)
        return dual.dual()

    @classmethod
    def from_vectors(
        cls,
        vectors: Mapping[T, Sequence[Fraction | int]],
        *,
        field_char: int | None = None,
    ) -> Matroid[T]:
        """Build the linear matroid of the given vectors.

        Independence is linear independence, computed by exact Gaussian
        elimination — over the rationals when ``field_char`` is ``None``,
        else over ``GF(field_char)``. Linear matroids always satisfy the
        matroid axioms (Whitney 1935, section 1), so there is nothing to
        validate beyond the input shape. Representability is only meaningful
        relative to a field (Matroid page, terminology hazard), which is why
        the characteristic is an explicit parameter.

        Args:
            vectors: Mapping from ground-set label to coordinate vector; the
                mapping order fixes bit positions.
            field_char: ``None`` for the rationals, or a prime ``p``.

        Returns:
            The linear matroid of the columns.

        Raises:
            ValueError: If vectors have inconsistent dimensions,
                ``field_char`` is not prime, or an entry's denominator is
                divisible by ``field_char`` (no image in the field).
        """
        labels = tuple(vectors)
        dimensions = {len(v) for v in vectors.values()}
        if len(dimensions) > 1:
            msg = f"vectors must share one dimension, got lengths {dimensions}"
            raise ValueError(msg)
        if field_char is not None and not is_prime(field_char):
            msg = f"field_char must be a prime number, got {field_char}"
            raise ValueError(msg)
        if field_char is None:
            family = linear_independence_family(
                [[Fraction(a) for a in vectors[label]] for label in labels],
                RATIONALS,
            )
        else:
            family = linear_independence_family(
                [
                    [gf_scalar(a, field_char) for a in vectors[label]]
                    for label in labels
                ],
                prime_field(field_char),
            )
        return cls(labels, family)

    @classmethod
    def from_graph_edges[V: Hashable](
        cls, edges: Mapping[T, tuple[V, V]]
    ) -> Matroid[T]:
        """Build the cycle matroid ``M(G)`` of a graph given as labeled edges.

        Independent sets are the forests (acyclic edge sets) — the graph
        abstraction the Matroid page records under *Abstracts*. Self-loops
        ``(v, v)`` become matroid loops and parallel edges become parallel
        elements. Forests form a matroid for every graph, so no axiom
        validation is needed. Note that ``M(G)`` forgets graph structure:
        non-isomorphic graphs can share a cycle matroid (Whitney's
        2-isomorphism theorem, per the page).

        Args:
            edges: Mapping from edge label to its endpoint pair; the mapping
                order fixes bit positions.

        Returns:
            The graphic matroid of the edge set.
        """
        labels = tuple(edges)
        pairs = list(edges.values())
        n = len(labels)
        family = frozenset(
            mask for mask in range(1 << n) if is_forest([pairs[b] for b in bits(mask)])
        )
        return cls(labels, family)

    @classmethod
    def from_transversal_system(
        cls, elements: Iterable[T], system: Iterable[Iterable[T]]
    ) -> Matroid[T]:
        """Build the transversal matroid ``M[A]`` of a set system.

        Independent sets are the partial transversals of the system
        (Edmonds-Fulkerson 1965); the independence oracle is a bipartite
        matching between elements and the sets containing them. Partial
        transversals always form a matroid, so no axiom validation is
        needed.

        Args:
            elements: Ground-set labels in order.
            system: The sets ``A_1, ..., A_m``, as label collections.

        Returns:
            The transversal matroid of the system.

        Raises:
            ValueError: If labels are unknown or duplicated.
        """
        elems, index = indexed_ground_set(elements)
        set_masks = [mask_from_labels(a, index) for a in system]
        n = len(elems)
        independent: set[int] = {0}
        for mask in sorted(range(1 << n), key=int.bit_count):
            if mask == 0:
                continue
            hereditary = all(mask ^ (1 << b) in independent for b in bits(mask))
            if hereditary and has_matching(mask, set_masks):
                independent.add(mask)
        return cls(elems, frozenset(independent))

    # ------------------------------------------------------- rank and views
    @functools.cached_property
    def _rank_table(self) -> list[int]:
        """Rank of every subset, by DP over the independence family.

        ``r(X) = |X|`` when ``X`` is independent, else the maximum rank of a
        one-element-smaller subset (cryptomorphism table, I -> r). Costs
        ``O(2^n * n)`` time and ``O(2^n)`` memory.
        """
        n = len(self.elements)
        family = self.independent_masks
        table = [0] * (1 << n)
        for mask in range(1, 1 << n):
            if mask in family:
                table[mask] = mask.bit_count()
            else:
                table[mask] = max(table[mask ^ (1 << b)] for b in bits(mask))
        return table

    def rank(self, subset: Iterable[T] | None = None) -> int:
        """Return ``r(X)``, or the rank ``r(M)`` when no subset is given.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        mask = self._full_mask if subset is None else self._mask(subset)
        return self._rank_table[mask]

    def nullity(self, subset: Iterable[T] | None = None) -> int:
        """Return the nullity ``n(X) = |X| - r(X)`` (defaults to ``E``).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        mask = self._full_mask if subset is None else self._mask(subset)
        return mask.bit_count() - self._rank_table[mask]

    def corank(self, subset: Iterable[T] | None = None) -> int:
        """Return the dual rank ``r*(X) = |X| + r(E - X) - r(E)``.

        The formula is Whitney's (1935; Matroid page, duality theorem).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        mask = self._full_mask if subset is None else self._mask(subset)
        return (
            mask.bit_count()
            + self._rank_table[self._full_mask ^ mask]
            - self._rank_table[self._full_mask]
        )

    def _closure_mask(self, mask: int) -> int:
        """Closure of a subset mask via ``cl(X) = {e : r(X + e) = r(X)}``."""
        rank_here = self._rank_table[mask]
        closed = mask
        for e in bits(self._full_mask & ~mask):
            if self._rank_table[mask | (1 << e)] == rank_here:
                closed |= 1 << e
        return closed

    def closure(self, subset: Iterable[T]) -> frozenset[T]:
        """Return ``cl(X)`` (cryptomorphism table, r -> cl).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return self._labels(self._closure_mask(self._mask(subset)))

    def is_independent(self, subset: Iterable[T]) -> bool:
        """Return whether the subset is independent.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return self._mask(subset) in self.independent_masks

    def is_dependent(self, subset: Iterable[T]) -> bool:
        """Return whether the subset is dependent (not independent).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return not self.is_independent(subset)

    def is_spanning(self, subset: Iterable[T]) -> bool:
        """Return whether ``r(X) = r(E)`` (Derived vocabulary: spanning).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return self._rank_table[self._mask(subset)] == self._rank_table[self._full_mask]

    def is_basis(self, subset: Iterable[T]) -> bool:
        """Return whether the subset is a basis (independent and spanning).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return self._mask(subset) in self._basis_masks

    def is_circuit(self, subset: Iterable[T]) -> bool:
        """Return whether the subset is a circuit (minimal dependent set).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        return self._mask(subset) in self._circuit_masks

    @functools.cached_property
    def _basis_masks(self) -> frozenset[int]:
        """Bases as masks: the maximal (equicardinal) independent sets."""
        top = self._rank_table[self._full_mask]
        return frozenset(m for m in self.independent_masks if m.bit_count() == top)

    @functools.cached_property
    def bases(self) -> frozenset[frozenset[T]]:
        """The bases: maximal independent sets (cryptomorphism, I -> B)."""
        return frozenset(self._labels(m) for m in self._basis_masks)

    @functools.cached_property
    def _circuit_masks(self) -> frozenset[int]:
        """Circuits as masks: minimal members of the dependent family."""
        family = self.independent_masks
        return frozenset(
            mask
            for mask in range(1 << len(self.elements))
            if mask not in family and all(mask ^ (1 << b) in family for b in bits(mask))
        )

    @functools.cached_property
    def circuits(self) -> frozenset[frozenset[T]]:
        """The circuits: minimal dependent sets (cryptomorphism, I -> C)."""
        return frozenset(self._labels(m) for m in self._circuit_masks)

    @functools.cached_property
    def _flat_masks(self) -> frozenset[int]:
        """Flats as masks: the closure-fixed subsets."""
        return frozenset(
            mask
            for mask in range(1 << len(self.elements))
            if self._closure_mask(mask) == mask
        )

    @functools.cached_property
    def flats(self) -> frozenset[frozenset[T]]:
        """The flats ``{X : cl(X) = X}``; they form the lattice of flats."""
        return frozenset(self._labels(m) for m in self._flat_masks)

    @functools.cached_property
    def _hyperplane_masks(self) -> frozenset[int]:
        """Hyperplanes as masks: flats of rank ``r(M) - 1``."""
        top = self._rank_table[self._full_mask]
        return frozenset(m for m in self._flat_masks if self._rank_table[m] == top - 1)

    @functools.cached_property
    def hyperplanes(self) -> frozenset[frozenset[T]]:
        """The hyperplanes: flats of rank ``r(M) - 1``."""
        return frozenset(self._labels(m) for m in self._hyperplane_masks)

    @functools.cached_property
    def cocircuits(self) -> frozenset[frozenset[T]]:
        """The cocircuits: complements of hyperplanes (Matroid page)."""
        return frozenset(
            self._labels(self._full_mask ^ m) for m in self._hyperplane_masks
        )

    @functools.cached_property
    def loops(self) -> frozenset[T]:
        """Elements with ``r({e}) = 0``; they lie in no independent set."""
        return frozenset(
            e
            for i, e in enumerate(self.elements)
            if 1 << i not in self.independent_masks
        )

    @functools.cached_property
    def coloops(self) -> frozenset[T]:
        """Elements lying in every basis (bridges / isthmuses)."""
        common = self._full_mask
        for basis in self._basis_masks:
            common &= basis
        return self._labels(common)

    @functools.cached_property
    def parallel_classes(self) -> frozenset[frozenset[T]]:
        """The parallel classes of non-loops (singletons included).

        Non-loops ``e, f`` are parallel when ``r({e, f}) = 1``; parallelism
        is an equivalence relation on non-loops (Derived vocabulary).
        """
        loop_labels = self.loops
        classes: list[list[int]] = []
        for i, element in enumerate(self.elements):
            if element in loop_labels:
                continue
            for group in classes:
                if self._rank_table[(1 << i) | (1 << group[0])] == 1:
                    group.append(i)
                    break
            else:
                classes.append([i])
        return frozenset(
            frozenset(self.elements[i] for i in group) for group in classes
        )

    @property
    def is_simple(self) -> bool:
        """Whether the matroid has no loops and no parallel pairs."""
        return not self.loops and all(
            len(group) == 1 for group in self.parallel_classes
        )

    @functools.cached_property
    def is_connected(self) -> bool:
        """Whether every two distinct elements lie in a common circuit.

        Equivalently, the matroid is not a direct sum of two matroids on
        nonempty ground sets (Derived vocabulary). Vacuously true for ground
        sets of size at most one.
        """
        n = len(self.elements)
        if n <= 1:
            return True
        covered: set[tuple[int, int]] = set()
        for circuit in self._circuit_masks:
            members = list(bits(circuit))
            covered.update(itertools.combinations(members, 2))
        return len(covered) == n * (n - 1) // 2

    def fundamental_circuit(self, element: T, basis: Iterable[T]) -> frozenset[T]:
        """Return ``C(e, B)``: the unique circuit inside ``B + e``.

        For a basis ``B`` and ``e`` outside it, exactly one circuit lies in
        ``B + e`` and it contains ``e`` (Matroid page, unique fundamental
        circuit theorem).

        Args:
            element: An element outside the basis.
            basis: A basis of the matroid.

        Returns:
            The fundamental circuit as a frozenset of labels.

        Raises:
            ValueError: If ``basis`` is not a basis, ``element`` lies in it,
                or labels are unknown.
        """
        basis_mask = self._mask(basis)
        if basis_mask not in self._basis_masks:
            msg = f"{fmt(basis_mask, self.elements)} is not a basis"
            raise ValueError(msg)
        element_mask = self._mask((element,))
        if element_mask & basis_mask:
            msg = f"element {element!r} lies in the basis"
            raise ValueError(msg)
        target = basis_mask | element_mask
        for circuit in self._circuit_masks:
            if circuit & ~target == 0:
                return self._labels(circuit)
        msg = "no circuit found; was this instance built without validation?"
        raise ValueError(msg)

    def simplification(self) -> Matroid[T]:
        """Return ``si(M)``: drop loops, keep one element per parallel class.

        The representative kept is the earliest element (in stored order) of
        each class, making the result deterministic.
        """
        keep = [
            min(group, key=lambda e: self._index[e]) for group in self.parallel_classes
        ]
        return self.restrict(keep)

    def greedy_max_weight_independent(self, weights: Mapping[T, float]) -> frozenset[T]:
        """Return a maximum-weight independent set via the greedy algorithm.

        Elements are taken in decreasing weight order (ties broken by stored
        order), each kept when it preserves independence; non-positive
        weights are never taken, since dropping them cannot decrease the
        total in a hereditary family. On a matroid the result has maximum
        weight for every weight function — the Rado-Edmonds greedy
        characterization (Matroid page, structural theorems).

        Args:
            weights: Weight of every ground-set element.

        Returns:
            A maximum-weight independent set.

        Raises:
            KeyError: If an element is missing from ``weights``.
        """
        order = sorted(
            range(len(self.elements)),
            key=lambda i: (-weights[self.elements[i]], i),
        )
        chosen = 0
        for i in order:
            if weights[self.elements[i]] <= 0:
                break
            if chosen | (1 << i) in self.independent_masks:
                chosen |= 1 << i
        return self._labels(chosen)

    # ------------------------------------------------------------ invariants
    def independent_set_counts(self) -> tuple[int, ...]:
        """Return ``(I_0, ..., I_r)``: independent sets counted by size.

        The sequence is log-concave (Mason's conjecture, proved by
        Anari-Liu-Oveis Gharan-Vinzant and Branden-Huh, 2020; Matroid page).
        """
        counts = [0] * (self.rank() + 1)
        for mask in self.independent_masks:
            counts[mask.bit_count()] += 1
        return tuple(counts)

    def tutte_polynomial(self) -> dict[tuple[int, int], int]:
        """Return the Tutte polynomial as ``{(i, j): coefficient}``.

        Computed from the corank-nullity sum
        ``T(x, y) = sum over A of (x-1)^(r(E)-r(A)) (y-1)^(|A|-r(A))``
        (Matroid page, structural theorems; the universal
        deletion-contraction invariant). Costs ``O(2^n)`` rank lookups.

        Returns:
            Mapping from ``(x-degree, y-degree)`` to nonzero coefficient.
        """
        top = self._rank_table[self._full_mask]
        shifted: Counter[tuple[int, int]] = Counter()
        for mask in range(1 << len(self.elements)):
            rank_a = self._rank_table[mask]
            shifted[(top - rank_a, mask.bit_count() - rank_a)] += 1
        coefficients: Counter[tuple[int, int]] = Counter()
        for (i, j), count in shifted.items():
            for p in range(i + 1):
                for q in range(j + 1):
                    sign = (-1) ** (i - p + j - q)
                    term = count * math.comb(i, p) * math.comb(j, q) * sign
                    coefficients[(p, q)] += term
        return {key: c for key, c in coefficients.items() if c}

    def tutte(self, x: float, y: float) -> float:
        """Evaluate the Tutte polynomial at ``(x, y)``.

        Exact for integer arguments. ``T(1,1)`` counts bases, ``T(2,1)``
        independent sets, ``T(1,2)`` spanning sets, and ``T(2,2) = 2^|E|``
        (Matroid page, evaluations as property tests).
        """
        return sum(c * x**p * y**q for (p, q), c in self.tutte_polynomial().items())

    def characteristic_polynomial(self) -> tuple[int, ...]:
        """Return ``p_M`` coefficients, index ``k`` holding the ``lambda^k`` term.

        ``p_M(lambda) = (-1)^r(E) T_M(1 - lambda, 0)`` (Matroid page, Tutte
        polynomial block). The coefficients — Whitney numbers of the first
        kind — form a log-concave sequence in absolute value
        (Adiprasito-Huh-Katz 2018). Identically zero when the matroid has a
        loop.
        """
        top = self._rank_table[self._full_mask]
        coefficients = [0] * (top + 1)
        sign_top = (-1) ** top
        for (a, b), c in self.tutte_polynomial().items():
            if b != 0:
                continue
            for k in range(a + 1):
                coefficients[k] += sign_top * c * math.comb(a, k) * (-1) ** k
        return tuple(coefficients)

    def ingleton_holds_for(
        self,
        a: Iterable[T],
        b: Iterable[T],
        c: Iterable[T],
        d: Iterable[T],
    ) -> bool:
        """Check Ingleton's inequality on one quadruple of subsets.

        Every representable matroid satisfies it for all quadruples; its
        violation certifies non-representability, which is exactly the
        Vamos matroid's role (Matroid page, canonical examples). The check
        is one-sided: satisfying it does not imply representability.

        Args:
            a: First subset of the ground set.
            b: Second subset.
            c: Third subset.
            d: Fourth subset.

        Returns:
            Whether the inequality holds for this quadruple.

        Raises:
            ValueError: If a subset contains unknown labels.
        """
        r = self._rank_table
        ma, mb, mc, md = (self._mask(s) for s in (a, b, c, d))
        lhs = r[ma | mb] + r[ma | mc] + r[ma | md] + r[mb | mc] + r[mb | md]
        rhs = r[ma] + r[mb] + r[mc | md] + r[ma | mb | mc] + r[ma | mb | md]
        return lhs >= rhs

    # ------------------------------------------------------- transformations
    def dual(self) -> Matroid[T]:
        """Return ``M*``: bases are the complements of bases (Whitney 1935).

        Duality is an involution (``M** = M``) and swaps loops with coloops,
        circuits with cocircuits, and deletion with contraction (Matroid
        page, operations and structural theorems).
        """
        complements = (self._full_mask ^ b for b in self._basis_masks)
        return Matroid(self.elements, down_closure(complements))

    def restrict(self, subset: Iterable[T]) -> Matroid[T]:
        """Return ``M|X``: ground set ``X``, independent sets those inside it.

        Element order is inherited from the parent matroid.

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        keep_mask = self._mask(subset)
        kept = [i for i in range(len(self.elements)) if keep_mask >> i & 1]
        new_bit = {old: new for new, old in enumerate(kept)}
        family = frozenset(
            remap(mask, new_bit)
            for mask in self.independent_masks
            if mask & ~keep_mask == 0
        )
        return Matroid(tuple(self.elements[i] for i in kept), family)

    def delete(self, subset: Iterable[T]) -> Matroid[T]:
        """Return ``M - X = M|(E - X)`` (Matroid page, deletion).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        drop_mask = self._mask(subset)
        return self.restrict(self._labels(self._full_mask ^ drop_mask))

    def contract(self, subset: Iterable[T]) -> Matroid[T]:
        """Return ``M/X``: independent sets ``I`` with ``I + B_X`` independent.

        ``B_X`` is a basis of ``X``; the result does not depend on the
        choice (Matroid page, contraction, basis form).

        Raises:
            ValueError: If the subset contains unknown labels.
        """
        contract_mask = self._mask(subset)
        base_of_x = 0
        for bit in bits(contract_mask):
            if self._rank_table[base_of_x | (1 << bit)] > self._rank_table[base_of_x]:
                base_of_x |= 1 << bit
        kept = [i for i in range(len(self.elements)) if not contract_mask >> i & 1]
        new_bit = {old: new for new, old in enumerate(kept)}
        family = frozenset(
            remap(mask, new_bit)
            for mask in self.independent_masks
            if mask & contract_mask == 0
            and (mask | base_of_x) in self.independent_masks
        )
        return Matroid(tuple(self.elements[i] for i in kept), family)

    def minor(
        self,
        *,
        deletions: Iterable[T] = (),
        contractions: Iterable[T] = (),
    ) -> Matroid[T]:
        """Return the minor ``M - D / C`` in its normal form.

        Deletion and contraction commute, so every minor has this form
        (Matroid page, minors).

        Args:
            deletions: Elements to delete (``D``).
            contractions: Elements to contract (``C``); disjoint from ``D``.

        Returns:
            The minor on ``E - D - C``.

        Raises:
            ValueError: If the two sets overlap or contain unknown labels.
        """
        delete_mask = self._mask(deletions)
        contract_mask = self._mask(contractions)
        if delete_mask & contract_mask:
            msg = "deletions and contractions must be disjoint"
            raise ValueError(msg)
        return self.delete(self._labels(delete_mask)).contract(
            self._labels(contract_mask)
        )

    def direct_sum[U: Hashable](self, other: Matroid[U]) -> Matroid[T | U]:
        """Return ``M1 (+) M2`` on the disjoint union of the ground sets.

        Raises:
            ValueError: If the ground sets share a label.
        """
        overlap = self.ground_set & other.ground_set
        if overlap:
            msg = f"direct sum needs disjoint ground sets; both contain {overlap!r}"
            raise ValueError(msg)
        elements: tuple[T | U, ...] = self.elements + other.elements
        shift = len(self.elements)
        family = frozenset(
            mine | (theirs << shift)
            for mine in self.independent_masks
            for theirs in other.independent_masks
        )
        return Matroid(elements, family)

    def truncation(self, k: int) -> Matroid[T]:
        """Return ``T_k(M)``: independent sets of size at most ``k``.

        Raises:
            ValueError: If ``k`` is negative.
        """
        if k < 0:
            msg = f"truncation size must be non-negative, got {k}"
            raise ValueError(msg)
        return Matroid(
            self.elements,
            frozenset(m for m in self.independent_masks if m.bit_count() <= k),
        )

    def _aligned_family(self, other: Matroid[T]) -> frozenset[int]:
        """Remap the other matroid's family onto this element ordering.

        Raises:
            ValueError: If the ground sets differ as sets.
        """
        if self.ground_set != other.ground_set:
            msg = "both matroids must share the same ground set"
            raise ValueError(msg)
        translate = [self._index[e] for e in other.elements]
        return frozenset(remap(mask, translate) for mask in other.independent_masks)

    def union(self, other: Matroid[T]) -> Matroid[T]:
        """Return ``M1 v M2``: independent sets are unions ``I1 + I2``.

        The rank of the union obeys the Nash-Williams formula
        ``r(X) = min over Y <= X of |X - Y| + r1(Y) + r2(Y)`` (Matroid page,
        structural theorems). Solvable in polynomial time with an
        independence oracle in general; this implementation is the explicit
        ``O(|I1| * |I2|)`` product, fitting the small-n envelope.

        Raises:
            ValueError: If the ground sets differ as sets.
        """
        other_family = self._aligned_family(other)
        family = frozenset(
            mine | theirs for mine in self.independent_masks for theirs in other_family
        )
        return Matroid(self.elements, family)

    def max_common_independent(self, other: Matroid[T]) -> frozenset[T]:
        """Return a maximum common independent set of two matroids.

        Its size equals ``min over X of r1(X) + r2(E - X)`` — Edmonds'
        matroid intersection theorem (Matroid page, structural theorems).
        Intersection of two matroids is tractable; three is NP-hard (page,
        operations). Ties are broken deterministically.

        Raises:
            ValueError: If the ground sets differ as sets.
        """
        common = self.independent_masks & self._aligned_family(other)
        best = max(common, key=lambda m: (m.bit_count(), -m))
        return self._labels(best)

    # ------------------------------------------------------- classification
    def is_isomorphic_to[U: Hashable](self, other: Matroid[U]) -> bool:
        """Return whether some ground-set bijection maps ``I`` onto ``I``.

        Backtracking over label bijections constrained by per-element
        signatures (counts of independent sets of each size through the
        element). Worst case is factorial in ``n``; fine for the small
        fixtures this library targets.
        """
        if len(self.elements) != len(other.elements):
            return False
        if self.independent_set_counts() != other.independent_set_counts():
            return False

        def signatures[W: Hashable](m: Matroid[W]) -> list[tuple[int, ...]]:
            top = m.rank()
            sigs = []
            for i in range(len(m.elements)):
                counts = [0] * (top + 1)
                for mask in m.independent_masks:
                    if mask >> i & 1:
                        counts[mask.bit_count()] += 1
                sigs.append(tuple(counts))
            return sigs

        sig_self = signatures(self)
        sig_other = signatures(other)
        if Counter(sig_self) != Counter(sig_other):
            return False
        n = len(self.elements)
        other_family = set(other.independent_masks)
        used = [False] * n
        assign = [0] * n

        def backtrack(i: int) -> bool:
            if i == n:
                remapped = {remap(mask, assign) for mask in self.independent_masks}
                return remapped == other_family
            for j in range(n):
                if not used[j] and sig_other[j] == sig_self[i]:
                    used[j] = True
                    assign[i] = j
                    if backtrack(i + 1):
                        return True
                    used[j] = False
            return False

        return backtrack(0)

    def has_minor[U: Hashable](self, other: Matroid[U]) -> bool:
        """Return whether some minor of this matroid is isomorphic to ``other``.

        Searches contractions ``C`` and deletions ``D`` with ``C``
        independent and ``D`` coindependent, which lose no minors (Oxley,
        *Matroid Theory*, 2nd ed. 2011, section 3.3). Exponential in the
        ground-set size; intended for the small fixtures of the excluded
        minor characterizations (Tutte; Matroid page).
        """
        m, k = len(self.elements), len(other.elements)
        if k > m:
            return False
        target_rank = other.rank()
        full = self._full_mask
        top = self._rank_table[full]
        for contract_mask in self.independent_masks:
            if contract_mask.bit_count() > m - k:
                continue
            rest = full ^ contract_mask
            wanted = m - k - contract_mask.bit_count()
            for delete_mask in submasks(rest):
                if delete_mask.bit_count() != wanted:
                    continue
                if self._rank_table[full ^ delete_mask] != top:
                    continue
                candidate = self.minor(
                    deletions=self._labels(delete_mask),
                    contractions=self._labels(contract_mask),
                )
                if candidate.rank() == target_rank and candidate.is_isomorphic_to(
                    other
                ):
                    return True
        return False

    def is_binary(self) -> bool:
        """Return whether the matroid is binary: no ``U_{2,4}`` minor.

        Tutte's excluded-minor characterization (Matroid page, structural
        theorems). Exponential; small ground sets only.
        """
        return not self.has_minor(u24())

    def is_regular(self) -> bool:
        """Return whether the matroid is regular.

        Tutte's characterization: no ``U_{2,4}``, ``F_7``, or ``F_7*`` minor
        (Matroid page, structural theorems). Exponential; small sets only.
        """
        return (
            self.is_binary()
            and not self.has_minor(fano_matroid())
            and not self.has_minor(fano_matroid().dual())
        )

    def is_graphic(self) -> bool:
        """Return whether the matroid is graphic.

        Tutte's characterization: no ``U_{2,4}``, ``F_7``, ``F_7*``,
        ``M*(K_5)``, or ``M*(K_{3,3})`` minor (Matroid page, structural
        theorems). Exponential; small ground sets only.
        """
        return (
            self.is_regular()
            and not self.has_minor(_mstar_k5())
            and not self.has_minor(_mstar_k33())
        )

    # ---------------------------------------------------------- dataframes
    def to_dataframe(self) -> pd.DataFrame:
        """Serialize to a tidy frame, one row per incidence.

        Columns are ``element`` and ``basis`` (nullable ``Int64``). Each
        ground-set element gets one row with a null ``basis`` (declaring the
        ground set, in stored order — required because ``E`` is not
        recoverable from the bases when loops exist), and each (basis,
        element) membership gets one row, bases numbered in a deterministic
        order. A rank-0 matroid therefore has only ground rows, decoding to
        the single empty basis. The encoding survives
        ``experiments.io.write_result`` (records-oriented JSON) when read
        back with ``pd.read_json(path, dtype=False)`` — default dtype
        inference re-parses numeric-looking string labels as integers.

        Returns:
            The tidy incidence frame; invert with :meth:`from_dataframe`.
        """
        element_column: list[T] = list(self.elements)
        basis_column: list[int | None] = [None] * len(self.elements)
        for basis_id, basis_mask in enumerate(sorted(self._basis_masks)):
            for bit in bits(basis_mask):
                element_column.append(self.elements[bit])
                basis_column.append(basis_id)
        return pd.DataFrame(
            {
                "element": element_column,
                "basis": pd.array(basis_column, dtype="Int64"),
            }
        )

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> Matroid[Hashable]:
        """Rebuild a matroid from a frame produced by :meth:`to_dataframe`.

        Tolerates the type drift of a JSON round trip (float basis ids,
        plain ``None`` markers) but not label-dtype drift — read JSON with
        ``dtype=False`` so string labels survive. A frame with no rows is
        the encoding of the empty matroid and decodes to it; its JSON form
        loses the columns too, so it is accepted before the column check.
        Validates through :meth:`from_bases`.

        Args:
            df: Frame with ``element`` and ``basis`` columns, or an empty
                frame for the empty matroid.

        Returns:
            The decoded matroid.

        Raises:
            ValueError: If required columns are missing or the decoded data
                violates the basis axioms.
        """
        if df.empty:
            return Matroid.from_bases((), [set()])
        missing = {"element", "basis"} - set(df.columns)
        if missing:
            msg = f"dataframe is missing required columns {sorted(missing)}"
            raise ValueError(msg)
        ground_rows = df["basis"].isna()
        elements = df.loc[ground_rows, "element"].tolist()
        grouped: dict[int, set[Hashable]] = {}
        membership = df.loc[~ground_rows]
        for label, basis_id in zip(
            membership["element"], membership["basis"], strict=True
        ):
            grouped.setdefault(int(basis_id), set()).add(label)
        bases: list[set[Hashable]] = [grouped[i] for i in sorted(grouped)]
        if not bases:
            bases = [set()]
        return Matroid.from_bases(elements, bases)

    # ------------------------------------------------------- visualization
    def _layout_flats(
        self,
    ) -> tuple[dict[int, tuple[float, float]], list[tuple[int, int]]]:
        """Compute node positions and cover edges for the lattice of flats."""
        by_rank: dict[int, list[int]] = {}
        for flat in sorted(self._flat_masks):
            by_rank.setdefault(self._rank_table[flat], []).append(flat)
        positions: dict[int, tuple[float, float]] = {}
        for level, row in by_rank.items():
            for i, flat in enumerate(row):
                positions[flat] = (i - (len(row) - 1) / 2, float(level))
        covers = [
            (low, high)
            for low in self._flat_masks
            for high in self._flat_masks
            if low & ~high == 0 and self._rank_table[high] == self._rank_table[low] + 1
        ]
        return positions, covers

    def plot_lattice_of_flats(self, ax: Axes | None = None) -> Axes:
        """Draw the Hasse diagram of the lattice of flats onto ``ax``.

        The lattice is geometric, so covers are exactly containments with a
        rank gap of one (Matroid page, Definition — Flats). Layers are rank
        levels, bottom to top.

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        positions, covers = self._layout_flats()
        for low, high in covers:
            (x1, y1), (x2, y2) = positions[low], positions[high]
            ax.plot([x1, x2], [y1, y2], color="0.7", linewidth=1, zorder=1)
        scatter_labeled(
            ax,
            list(positions.values()),
            [fmt(flat, self.elements) for flat in positions],
        )
        ax.set_axis_off()
        ax.set_title("Lattice of flats")
        return ax

    def plot_basis_exchange_graph(self, ax: Axes | None = None) -> Axes:
        """Draw the basis exchange graph onto ``ax``.

        Nodes are bases on a circle; edges join bases whose symmetric
        difference has size two. The graph is connected (Matroid page,
        structural theorems).

        Args:
            ax: Axes to draw on; a new figure is created when omitted.

        Returns:
            The axes drawn on. Never calls ``show`` or writes files.
        """
        ax = ensure_axes(ax)
        ordered = sorted(self._basis_masks)
        positions = dict(zip(ordered, unit_circle(len(ordered)), strict=True))
        for b1, b2 in itertools.combinations(ordered, 2):
            if (b1 ^ b2).bit_count() == _EXCHANGE_DISTANCE:
                (x1, y1), (x2, y2) = positions[b1], positions[b2]
                ax.plot([x1, x2], [y1, y2], color="0.7", linewidth=1, zorder=1)
        scatter_labeled(
            ax,
            list(positions.values()),
            [fmt(basis, self.elements) for basis in ordered],
        )
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title("Basis exchange graph")
        return ax


# --------------------------------------------------------------------------- #
# Canonical examples — test fixtures per the Matroid page, exported as API
# --------------------------------------------------------------------------- #
def uniform_matroid(rank: int, n: int) -> Matroid[int]:
    """Return ``U_{r,n}`` on ``range(n)``: independent iff size at most ``r``.

    Args:
        rank: The rank ``r``.
        n: The ground-set size.

    Returns:
        The uniform matroid, satisfying ``U_{r,n}* = U_{n-r,n}``.

    Raises:
        ValueError: If ``0 <= rank <= n`` fails.
    """
    if not 0 <= rank <= n:
        msg = f"uniform matroid needs 0 <= rank <= n, got rank={rank}, n={n}"
        raise ValueError(msg)
    family = frozenset(m for m in range(1 << n) if m.bit_count() <= rank)
    return Matroid(tuple(range(n)), family)


def empty_matroid() -> Matroid[int]:
    """Return ``U_{0,0}``, the empty matroid."""
    return uniform_matroid(0, 0)


def free_matroid(n: int) -> Matroid[int]:
    """Return the free matroid ``U_{n,n}``: every subset independent.

    Raises:
        ValueError: If ``n`` is negative.
    """
    return uniform_matroid(n, n)


def loopy_matroid(n: int) -> Matroid[int]:
    """Return the loopy matroid ``U_{0,n}``: every element a loop.

    Raises:
        ValueError: If ``n`` is negative.
    """
    return uniform_matroid(0, n)


def u24() -> Matroid[int]:
    """Return ``U_{2,4}``: the excluded minor for binariness.

    The smallest non-binary matroid; rank 2 on 4 elements and self-dual
    (Matroid page, canonical examples).
    """
    return uniform_matroid(2, 4)


def k4_matroid() -> Matroid[str]:
    """Return ``M(K_4)``: the graphic matroid of the complete graph ``K_4``.

    Rank 3 on 6 elements; the standard small graphic fixture (Matroid page,
    canonical examples). Edge labels name their endpoints.
    """
    edges = {
        "12": (1, 2),
        "13": (1, 3),
        "14": (1, 4),
        "23": (2, 3),
        "24": (2, 4),
        "34": (3, 4),
    }
    return Matroid.from_graph_edges(edges)


_FANO_LINES: frozenset[frozenset[int]] = frozenset(
    frozenset(line)
    for line in [
        (1, 2, 4),
        (1, 3, 5),
        (1, 6, 7),
        (2, 3, 6),
        (2, 5, 7),
        (3, 4, 7),
        (4, 5, 6),
    ]
)
"""The seven lines of the Fano plane (Whitney 1935, section 16)."""


def _all_triples_except(lines: frozenset[frozenset[int]]) -> list[frozenset[int]]:
    """Return every 3-subset of ``{1..7}`` that is not one of ``lines``."""
    return [
        frozenset(triple)
        for triple in itertools.combinations(range(1, 8), 3)
        if frozenset(triple) not in lines
    ]


def fano_matroid() -> Matroid[int]:
    """Return the Fano matroid ``F_7``.

    Rank 3 on ``{1..7}``; bases are all 3-sets except the seven lines
    124, 135, 167, 236, 257, 347, 456 (Whitney 1935, section 16 — the lines
    of the Fano plane). Representable precisely over fields of
    characteristic 2 (Matroid page, canonical examples).
    """
    return Matroid.from_bases(range(1, 8), _all_triples_except(_FANO_LINES))


def non_fano_matroid() -> Matroid[int]:
    """Return the non-Fano matroid ``F_7^-``.

    The Fano matroid with the line 456 removed — that 3-set becomes a basis
    — representable precisely over characteristic other than 2; the pair
    with ``F_7`` is the standard characteristic-sensitivity test (Matroid
    page, canonical examples).
    """
    lines = frozenset(line for line in _FANO_LINES if line != frozenset({4, 5, 6}))
    return Matroid.from_bases(range(1, 8), _all_triples_except(lines))


def vamos_matroid() -> Matroid[int]:
    """Return the Vamos matroid ``V_8``.

    Rank 4 on ``{1..8}``: bases are all 4-sets except the five
    circuit-hyperplanes 1234, 1256, 1278, 3456, 3478 (Vamos 1968; Oxley
    2011 — the defining planes are not recorded on the Matroid page, which
    is a reported gap). Representable over no field: it violates Ingleton's
    inequality at ``({1,2}, {3,4}, {5,6}, {7,8})``, the standard certificate
    (Matroid page, canonical examples). Also non-algebraic (Ingleton-Main
    1975).
    """
    planes = frozenset(
        frozenset(p)
        for p in [(1, 2, 3, 4), (1, 2, 5, 6), (1, 2, 7, 8), (3, 4, 5, 6), (3, 4, 7, 8)]
    )
    bases = [
        frozenset(quad)
        for quad in itertools.combinations(range(1, 9), 4)
        if frozenset(quad) not in planes
    ]
    return Matroid.from_bases(range(1, 9), bases)


@functools.cache
def _mstar_k5() -> Matroid[str]:
    """Return ``M*(K_5)``, an excluded minor for graphicness (Tutte)."""
    vertices = range(1, 6)
    edges = {f"{u}{v}": (u, v) for u, v in itertools.combinations(vertices, 2)}
    return Matroid.from_graph_edges(edges).dual()


@functools.cache
def _mstar_k33() -> Matroid[str]:
    """Return ``M*(K_{3,3})``, an excluded minor for graphicness (Tutte)."""
    edges = {f"{u}{v}": (u, v) for u in (1, 2, 3) for v in (4, 5, 6)}
    return Matroid.from_graph_edges(edges).dual()


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #
def enumerate_matroids(n: int) -> list[Matroid[int]]:
    """Return all matroids on ``n`` labeled elements, up to isomorphism.

    Searches every equicardinal basis family satisfying basis exchange and
    deduplicates by isomorphism. The count for ``n = 0..8`` is
    1, 2, 4, 8, 17, 38, 98, 306, 1724 (OEIS A055545; Matroid page,
    enumeration regression) — any enumerator must reproduce this prefix.
    Doubly exponential in ``n``; practical only for ``n <= 5``.

    Args:
        n: Ground-set size.

    Returns:
        One representative per isomorphism class, on ``range(n)``.

    Raises:
        ValueError: If ``n`` is negative.
    """
    if n < 0:
        msg = f"ground-set size must be non-negative, got {n}"
        raise ValueError(msg)
    elements = tuple(range(n))
    buckets: dict[tuple[int, tuple[int, ...]], list[Matroid[int]]] = {}
    found: list[Matroid[int]] = []
    for k in range(n + 1):
        k_subsets = [
            sum(1 << i for i in combo) for combo in itertools.combinations(range(n), k)
        ]
        for pick in range(1, 1 << len(k_subsets)):
            base_masks = frozenset(k_subsets[i] for i in bits(pick))
            if basis_exchange_violation(base_masks) is not None:
                continue
            candidate = Matroid(elements, down_closure(base_masks))
            key = (k, candidate.independent_set_counts())
            bucket = buckets.setdefault(key, [])
            if any(candidate.is_isomorphic_to(seen) for seen in bucket):
                continue
            bucket.append(candidate)
            found.append(candidate)
    return found
