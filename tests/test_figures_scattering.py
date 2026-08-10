"""Tests for the colour-ordered tree enumeration behind the Feynman gallery.

The gallery draws diagrams produced by
:func:`figures.scattering_amplitudes.dissections`, and the growth curve plots
counts produced by :func:`figures.scattering_amplitudes.dissection_count`.
Both come from one recursion, so the published sequence is pinned here rather
than in the module: if the enumerator drifts, the curve and the gallery drift
together and only this file notices.
"""

import itertools

import pytest

from figures.scattering_amplitudes import dissection_count, dissections

PUBLISHED_COUNTS: dict[int, int] = {
    3: 1,
    4: 3,
    5: 10,
    6: 38,
    7: 154,
    8: 654,
    9: 2871,
    10: 12925,
}
"""Colour-ordered tree diagrams on ``n`` legs, ``n`` mapping to the count.

Elvang and Huang, *Scattering Amplitudes in Gauge Theory and Gravity* §2.6 for
``n <= 7``; the continuation is recorded as verified on the knowledge graph's
**BCFW Recursion** page.
"""


@pytest.mark.parametrize(("legs", "expected"), sorted(PUBLISHED_COUNTS.items()))
def test_dissection_count_matches_published_sequence(legs: int, expected: int) -> None:
    """The counted dissections reproduce the published diagram counts."""
    assert dissection_count(legs) == expected


@pytest.mark.parametrize("legs", [3, 4, 5, 6, 7])
def test_enumeration_agrees_with_the_count(legs: int) -> None:
    """Materializing the dissections yields exactly as many as counting does."""
    assert len(dissections(tuple(range(legs)))) == dissection_count(legs)


@pytest.mark.parametrize("legs", [4, 5, 6, 7])
def test_every_dissection_is_distinct(legs: int) -> None:
    """No dissection is enumerated twice, as a set of unordered cells."""
    found = dissections(tuple(range(legs)))
    assert len({frozenset(found_cells) for found_cells in found}) == len(found)


@pytest.mark.parametrize("legs", [4, 5, 6, 7])
def test_cells_are_triangles_or_quadrilaterals(legs: int) -> None:
    """Every cell comes from a cubic or quartic vertex, and nothing else."""
    sizes = {len(cell) for found in dissections(tuple(range(legs))) for cell in found}
    assert sizes == {3, 4}


@pytest.mark.parametrize("legs", [4, 5, 6, 7])
def test_every_polygon_edge_lies_in_exactly_one_cell(legs: int) -> None:
    """Each boundary edge has one home, so each external leg has one vertex.

    :func:`figures.scattering_amplitudes._diagram` attaches a leg to
    ``sides[edge][0]``; that is only well defined because the choice is
    unique.
    """
    boundary = [frozenset((i, (i + 1) % legs)) for i in range(legs)]
    for found in dissections(tuple(range(legs))):
        counts = dict.fromkeys(boundary, 0)
        for cell in found:
            for a, b in itertools.pairwise((*cell, cell[0])):
                edge = frozenset((a, b))
                if edge in counts:
                    counts[edge] += 1
        assert set(counts.values()) == {1}


@pytest.mark.parametrize("legs", [4, 5, 6, 7])
def test_internal_lines_join_exactly_two_cells(legs: int) -> None:
    """Every diagonal is shared by two cells, so every internal line has ends.

    The drawing routine reads an internal line off any side belonging to two
    cells and ignores the rest; a side in three cells would silently lose one.
    """
    boundary = {frozenset((i, (i + 1) % legs)) for i in range(legs)}
    for found in dissections(tuple(range(legs))):
        shared: dict[frozenset[int], int] = {}
        for cell in found:
            for a, b in itertools.pairwise((*cell, cell[0])):
                edge = frozenset((a, b))
                shared[edge] = shared.get(edge, 0) + 1
        diagonals = {e: n for e, n in shared.items() if e not in boundary}
        assert set(diagonals.values()) <= {2}


def test_four_particle_diagrams_are_the_three_textbook_topologies() -> None:
    """Two cubic exchanges and one contact term, and nothing else."""
    found = dissections((0, 1, 2, 3))
    shapes = sorted(sorted(len(cell) for cell in cells) for cells in found)
    assert shapes == [[3, 3], [3, 3], [4]]
