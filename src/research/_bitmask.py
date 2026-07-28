"""Shared bitmask-subset helpers for the combinatorial structure modules.

Subsets of an ordered ground set are stored as ``int`` bitmasks throughout
``research`` (Matroid page, implementation notes); this internal module
holds the label/mask plumbing shared by :mod:`research.matroid` and
:mod:`research.positroid` so each structure module stays about one concern.
"""

from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence

__all__ = [
    "bits",
    "down_closure",
    "fmt",
    "indexed_ground_set",
    "mask_from_labels",
    "remap",
    "require_distinct",
    "submasks",
]


def bits(mask: int) -> Iterator[int]:
    """Yield the set bit positions of ``mask`` in ascending order."""
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def submasks(mask: int) -> Iterator[int]:
    """Yield every submask of ``mask``, including ``0`` and ``mask`` itself."""
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & mask


def down_closure(masks: Iterable[int]) -> frozenset[int]:
    """Return the downward closure (all submasks) of the given masks."""
    closed: set[int] = set()
    for mask in masks:
        closed.update(submasks(mask))
    return frozenset(closed)


def fmt(mask: int, elements: tuple[Hashable, ...]) -> str:
    """Render a bitmask as a readable set of ground-set labels."""
    if mask == 0:
        return "{}"
    return "{" + ", ".join(repr(elements[b]) for b in bits(mask)) + "}"


def mask_from_labels[T: Hashable](labels: Iterable[T], index: Mapping[T, int]) -> int:
    """Convert a collection of labels to a bitmask.

    Args:
        labels: Labels to convert; duplicates are collapsed.
        index: Label-to-bit-position mapping for the ground set.

    Returns:
        The bitmask with one bit per distinct label.

    Raises:
        ValueError: If a label is not in the ground set.
    """
    mask = 0
    for label in labels:
        position = index.get(label)
        if position is None:
            msg = f"element {label!r} is not in the ground set"
            raise ValueError(msg)
        mask |= 1 << position
    return mask


def remap(mask: int, table: Mapping[int, int] | Sequence[int]) -> int:
    """Return the mask with each set bit ``b`` moved to position ``table[b]``.

    Every set bit of ``mask`` must have an entry in ``table``; targets must
    be distinct for the result to preserve cardinality.
    """
    return sum(1 << table[b] for b in bits(mask))


def indexed_ground_set[T: Hashable](
    elements: Iterable[T],
) -> tuple[tuple[T, ...], dict[T, int]]:
    """Return the ground tuple and its label-to-bit index.

    The shared preamble of every ``from_<formulation>`` constructor: fix the
    element order, reject duplicates, and index labels by bit position.

    Raises:
        ValueError: If two elements compare equal.
    """
    elems = tuple(elements)
    require_distinct(elems)
    return elems, {e: i for i, e in enumerate(elems)}


def require_distinct(elements: tuple[Hashable, ...]) -> None:
    """Reject ground sets with repeated labels.

    Raises:
        ValueError: If two elements compare equal.
    """
    if len(set(elements)) != len(elements):
        msg = f"ground set labels must be distinct, got {elements!r}"
        raise ValueError(msg)
