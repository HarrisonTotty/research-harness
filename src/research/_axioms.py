"""Matroid axiom checkers, one per axiomatization on the Matroid page.

Each checker validates one cryptomorphic presentation — independent sets,
bases, circuits, rank, closure, flats, hyperplanes — against the numbered
axioms exactly as the Matroid page numbers them, and every error message
names the first violated axiom. :mod:`research.matroid` calls these from its
``from_<formulation>`` constructors; nothing here constructs a matroid.
"""

import itertools
from collections.abc import Hashable, Sequence

from research._bitmask import bits, fmt

__all__ = [
    "basis_exchange_violation",
    "check_basis_axioms",
    "check_circuit_axioms",
    "check_closure_axioms",
    "check_flat_axioms",
    "check_hyperplane_axioms",
    "check_independence_axioms",
    "check_rank_axioms",
]


def check_independence_axioms(
    elements: tuple[Hashable, ...], masks: frozenset[int]
) -> None:
    """Check the independent-set axioms (I1)-(I3) (Whitney 1935, (a)/(b)).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    if 0 not in masks:
        msg = "(I1) violated: the empty set must be independent"
        raise ValueError(msg)
    for mask in masks:
        for bit in bits(mask):
            sub = mask ^ (1 << bit)
            if sub not in masks:
                msg = (
                    f"(I2) hereditary axiom violated: {fmt(mask, elements)} is "
                    f"independent but its subset {fmt(sub, elements)} is not"
                )
                raise ValueError(msg)
    for small in masks:
        for large in masks:
            if small.bit_count() >= large.bit_count():
                continue
            can_augment = any(
                small | (1 << bit) in masks for bit in bits(large & ~small)
            )
            if not can_augment:
                msg = (
                    f"(I3) augmentation axiom violated: {fmt(small, elements)} "
                    f"cannot be extended from {fmt(large, elements)}"
                )
                raise ValueError(msg)


def basis_exchange_violation(
    base_masks: frozenset[int],
) -> tuple[int, int, int] | None:
    """Return a witness ``(b1, b2, x)`` violating (B2), or ``None`` if none."""
    for b1 in base_masks:
        for b2 in base_masks:
            for x in bits(b1 & ~b2):
                exchanged = (
                    (b1 ^ (1 << x)) | (1 << y) in base_masks for y in bits(b2 & ~b1)
                )
                if not any(exchanged):
                    return (b1, b2, x)
    return None


def check_basis_axioms(
    elements: tuple[Hashable, ...], base_masks: frozenset[int]
) -> None:
    """Check the basis axioms (B1)-(B2).

    Equicardinality of bases is a theorem, not an axiom, and is deliberately
    not checked here (Matroid page, Definition — Bases).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    if not base_masks:
        msg = "(B1) violated: the set of bases must be nonempty"
        raise ValueError(msg)
    witness = basis_exchange_violation(base_masks)
    if witness is not None:
        b1, b2, x = witness
        msg = (
            f"(B2) basis exchange violated: no element of {fmt(b2, elements)} "
            f"can replace {elements[x]!r} in {fmt(b1, elements)}"
        )
        raise ValueError(msg)


def check_circuit_axioms(
    elements: tuple[Hashable, ...], circuit_masks: frozenset[int]
) -> None:
    """Check the circuit axioms (C1)-(C3).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    if 0 in circuit_masks:
        msg = "(C1) violated: the empty set cannot be a circuit"
        raise ValueError(msg)
    for c1 in circuit_masks:
        for c2 in circuit_masks:
            if c1 != c2 and c1 & ~c2 == 0:
                msg = (
                    f"(C2) antichain axiom violated: circuit {fmt(c1, elements)} "
                    f"is a proper subset of circuit {fmt(c2, elements)}"
                )
                raise ValueError(msg)
    for c1 in circuit_masks:
        for c2 in circuit_masks:
            if c1 == c2:
                continue
            for e in bits(c1 & c2):
                allowed = (c1 | c2) ^ (1 << e)
                if not any(c3 & ~allowed == 0 for c3 in circuit_masks):
                    msg = (
                        f"(C3) circuit elimination violated for "
                        f"{fmt(c1, elements)} and {fmt(c2, elements)} at "
                        f"{elements[e]!r}"
                    )
                    raise ValueError(msg)


def check_rank_axioms(elements: tuple[Hashable, ...], table: Sequence[int]) -> None:
    """Check rank axioms (R1)-(R3) via Whitney's local forms (1935, section 2).

    The unit-increase and local-flatness checks cost ``O(2^n * n^2)`` versus
    ``O(4^n)`` for pairwise submodularity; their failures imply failures of
    the page-numbered axioms as noted in each message (Whitney 1935, Thm. 3).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    n = len(elements)
    for mask in range(1 << n):
        if not 0 <= table[mask] <= mask.bit_count():
            msg = (
                f"(R1) violated: r({fmt(mask, elements)}) = {table[mask]} "
                f"is outside [0, {mask.bit_count()}]"
            )
            raise ValueError(msg)
    for mask in range(1 << n):
        for e in range(n):
            if mask >> e & 1:
                continue
            step = table[mask | (1 << e)] - table[mask]
            if step < 0:
                msg = (
                    f"(R2) monotonicity violated: rank drops when adding "
                    f"{elements[e]!r} to {fmt(mask, elements)}"
                )
                raise ValueError(msg)
            if step > 1:
                msg = (
                    f"(R3) submodularity (with (R1)) violated: rank grows by "
                    f"{step} when adding {elements[e]!r} to {fmt(mask, elements)}"
                )
                raise ValueError(msg)
    for mask in range(1 << n):
        outside = [e for e in range(n) if not mask >> e & 1]
        for e1, e2 in itertools.combinations(outside, 2):
            flat1 = table[mask | (1 << e1)] == table[mask]
            flat2 = table[mask | (1 << e2)] == table[mask]
            joint = table[mask | (1 << e1) | (1 << e2)]
            if flat1 and flat2 and joint != table[mask]:
                msg = (
                    f"(R3) submodularity violated: {elements[e1]!r} and "
                    f"{elements[e2]!r} each leave r({fmt(mask, elements)}) "
                    f"unchanged but together raise it"
                )
                raise ValueError(msg)


def check_closure_axioms(elements: tuple[Hashable, ...], table: Sequence[int]) -> None:
    """Check the closure axioms (CL1)-(CL4) (Mac Lane-Steinitz exchange).

    Monotonicity (CL2) is checked in single-element steps, which implies the
    general form by induction along a chain.

    Raises:
        ValueError: Naming the first violated axiom.
    """
    n = len(elements)
    for mask in range(1 << n):
        if mask & ~table[mask]:
            msg = (
                f"(CL1) extensivity violated: cl({fmt(mask, elements)}) does "
                f"not contain {fmt(mask, elements)}"
            )
            raise ValueError(msg)
    for mask in range(1 << n):
        for e in range(n):
            if mask >> e & 1:
                continue
            if table[mask] & ~table[mask | (1 << e)]:
                msg = (
                    f"(CL2) monotonicity violated between {fmt(mask, elements)} "
                    f"and {fmt(mask | (1 << e), elements)}"
                )
                raise ValueError(msg)
    for mask in range(1 << n):
        if table[table[mask]] != table[mask]:
            msg = (
                f"(CL3) idempotence violated: cl(cl({fmt(mask, elements)})) "
                f"differs from cl({fmt(mask, elements)})"
            )
            raise ValueError(msg)
    for mask in range(1 << n):
        for x in range(n):
            gained = table[mask | (1 << x)] & ~table[mask]
            for y in bits(gained):
                if not table[mask | (1 << y)] >> x & 1:
                    msg = (
                        f"(CL4) Mac Lane-Steinitz exchange violated at "
                        f"X = {fmt(mask, elements)}, x = {elements[x]!r}, "
                        f"y = {elements[y]!r}"
                    )
                    raise ValueError(msg)


def check_flat_axioms(
    elements: tuple[Hashable, ...], flat_masks: frozenset[int]
) -> None:
    """Check the flat axioms (F1)-(F3).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    full = (1 << len(elements)) - 1
    if full not in flat_masks:
        msg = "(F1) violated: the ground set must be a flat"
        raise ValueError(msg)
    for f1 in flat_masks:
        for f2 in flat_masks:
            if f1 & f2 not in flat_masks:
                msg = (
                    f"(F2) violated: the intersection of flats "
                    f"{fmt(f1, elements)} and {fmt(f2, elements)} is not a flat"
                )
                raise ValueError(msg)
    for flat in flat_masks:
        above = [
            g for g in flat_masks if g != flat and g & ~flat != 0 and flat & ~g == 0
        ]
        minimal = [
            g
            for g in above
            if not any(h != g and flat | h != flat and h & ~g == 0 for h in above)
        ]
        covered = 0
        for g in minimal:
            part = g & ~flat
            if covered & part:
                msg = (
                    f"(F3) covering axiom violated: minimal flats above "
                    f"{fmt(flat, elements)} overlap outside it"
                )
                raise ValueError(msg)
            covered |= part
        if covered != full & ~flat:
            msg = (
                f"(F3) covering axiom violated: minimal flats above "
                f"{fmt(flat, elements)} do not cover its complement"
            )
            raise ValueError(msg)


def check_hyperplane_axioms(
    elements: tuple[Hashable, ...], h_masks: frozenset[int]
) -> None:
    """Check the hyperplane axioms (H1)-(H3).

    Raises:
        ValueError: Naming the first violated axiom.
    """
    full = (1 << len(elements)) - 1
    if full in h_masks:
        msg = "(H1) violated: the ground set cannot be a hyperplane"
        raise ValueError(msg)
    for h1 in h_masks:
        for h2 in h_masks:
            if h1 != h2 and h1 & ~h2 == 0:
                msg = (
                    f"(H2) antichain axiom violated: {fmt(h1, elements)} is a "
                    f"proper subset of {fmt(h2, elements)}"
                )
                raise ValueError(msg)
    for h1 in h_masks:
        for h2 in h_masks:
            if h1 == h2:
                continue
            for e in bits(full & ~(h1 | h2)):
                required = (h1 & h2) | (1 << e)
                if not any(required & ~h3 == 0 for h3 in h_masks):
                    msg = (
                        f"(H3) violated for {fmt(h1, elements)} and "
                        f"{fmt(h2, elements)} at {elements[e]!r}"
                    )
                    raise ValueError(msg)
