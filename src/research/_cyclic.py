"""Gale orders and Grassmann-necklace machinery for positroids.

Everything here works on bitmasks over a cyclically ordered ground set: the
shifted Gale order comparisons, the greedy Grassmann necklace of an
independence family, Oh's membership test for the positroid property, and
the necklace consistency conditions of Postnikov section 16.
:mod:`research.positroid` builds its constructors and views on these.
"""

import itertools
from collections.abc import Hashable, Sequence

from research._bitmask import bits, fmt
from research.matroid import Matroid

__all__ = [
    "check_necklace_conditions",
    "check_positroid",
    "gale_geq",
    "necklace_bases",
    "necklace_masks",
    "positroid_witness",
]


def gale_geq(candidate: int, minimum: int, start: int, n: int) -> bool:
    """Return whether ``candidate >= minimum`` in the Gale order ``<=_start``.

    Both masks must have equal cardinality. Positions are compared in the
    cyclically shifted order ``start <_start start+1 <_start ...`` after
    sorting both sets (Positroid page, Grassmann necklace block: the Gale
    order on ``d``-subsets).
    """
    shifted_candidate = sorted((p - start) % n for p in bits(candidate))
    shifted_minimum = sorted((p - start) % n for p in bits(minimum))
    return all(c >= m for c, m in zip(shifted_candidate, shifted_minimum, strict=True))


def necklace_masks(n: int, family: frozenset[int]) -> tuple[int, ...]:
    """Return the Grassmann necklace of an independence family, as masks.

    Entry ``i`` is the ``<=_i``-minimal basis, computed greedily in the
    cyclically shifted order — the lexicographically minimal basis w.r.t.
    ``<_i``, which the Positroid page identifies with the Gale-minimal one
    (Grassmann necklace block, "From a matroid").
    """
    masks: list[int] = []
    for start in range(n):
        current = 0
        for offset in range(n):
            e = (start + offset) % n
            if current | (1 << e) in family:
                current |= 1 << e
        masks.append(current)
    return tuple(masks)


def necklace_bases(masks: Sequence[int], n: int) -> list[int]:
    """Return Oh's construction ``B(I) = {B : B >=_j I_j for all j}``.

    The ``d``-subsets dominating every necklace entry in its shifted Gale
    order, as masks in combination order (Positroid page, Oh's theorem
    block).
    """
    d = masks[0].bit_count() if masks else 0
    found: list[int] = []
    for combo in itertools.combinations(range(n), d):
        mask = sum(1 << p for p in combo)
        if all(gale_geq(mask, masks[s], s, n) for s in range(n)):
            found.append(mask)
    return found


def positroid_witness[T: Hashable](matroid: Matroid[T]) -> int | None:
    """Return a mask witnessing failure of Oh's membership test, or ``None``.

    Oh's theorem: a matroid is a positroid for its stored cyclic order iff
    its bases are exactly ``{B : B >=_j I_j for all j}`` where ``I`` is its
    Grassmann necklace (Positroid page, Oh's theorem block). The witness is
    a ``d``-subset on the wrong side of that equality.
    """
    n = len(matroid.elements)
    family = matroid.independent_masks
    necklace = necklace_masks(n, family)
    d = matroid.rank()
    basis_masks = {m for m in family if m.bit_count() == d}
    members = set(necklace_bases(necklace, n))
    for combo in itertools.combinations(range(n), d):
        mask = sum(1 << p for p in combo)
        if (mask in members) != (mask in basis_masks):
            return mask
    return None


def check_positroid[T: Hashable](matroid: Matroid[T]) -> None:
    """Check the positroid property via Oh's theorem.

    Raises:
        ValueError: Naming Oh's theorem, with the witness subset.
    """
    witness = positroid_witness(matroid)
    if witness is None:
        return
    formatted = fmt(witness, matroid.elements)
    if witness in matroid.independent_masks:
        msg = (
            f"not a positroid for this cyclic order: basis {formatted} "
            f"fails B >=_j I_j against the Grassmann necklace, violating "
            f"Oh's theorem"
        )
    else:
        msg = (
            f"not a positroid for this cyclic order: by Oh's theorem the "
            f"bases must be exactly {{B : B >=_j I_j for all j}}, but "
            f"{formatted} satisfies every cyclically shifted Schubert "
            f"condition without being a basis"
        )
    raise ValueError(msg)


def check_necklace_conditions(
    elements: tuple[Hashable, ...], masks: Sequence[int]
) -> None:
    r"""Check the Grassmann necklace conditions (Postnikov section 16).

    Indices modulo ``n``: if ``i`` is in ``I_i`` then ``I_{i+1}`` must be
    ``(I_i \ {i}) + {j}`` for some ``j``; otherwise ``I_{i+1} = I_i``.

    Raises:
        ValueError: Naming the violated condition.
    """
    n = len(elements)
    for i in range(n):
        current = masks[i]
        following = masks[(i + 1) % n]
        if current >> i & 1:
            base = current ^ (1 << i)
            added = following & ~base
            if base & ~following or added.bit_count() != 1:
                msg = (
                    f"Grassmann necklace condition violated (Postnikov "
                    f"section 16): {elements[i]!r} is in "
                    f"I_{i + 1} = {fmt(current, elements)}, so the next "
                    f"entry must be (I_i \\ {{i}}) + {{j}}, got "
                    f"{fmt(following, elements)}"
                )
                raise ValueError(msg)
        elif following != current:
            msg = (
                f"Grassmann necklace condition violated (Postnikov "
                f"section 16): {elements[i]!r} is not in "
                f"I_{i + 1} = {fmt(current, elements)}, so the next entry "
                f"must equal it, got {fmt(following, elements)}"
            )
            raise ValueError(msg)
