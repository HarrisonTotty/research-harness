"""Tests for research.strand_dynamics.

Strand dynamics has no Logseq page — the scheme is original to this project
— so the oracles are, in order of independence:

* hand-computed fixtures small enough to trace on paper (the single bridge
  of Gr(1,2), lollipops);
* the library's own trip tracer: loads and trip permutations of arbitrary
  colorings are recomputed from ``PlabicGraph.trips`` of the recolored graph,
  which shares no code with the automaton's dart tables;
* the published tables of "Finding 020" of the predecessor repository
  (``docs/findings/020-strand-gradient-landscape.md`` there), produced by an
  independent implementation of the same rule on the same bridge graphs.

The orbit laws (a repeat is a true period; the trip-permutation period
divides it) are properties of any deterministic map and are tested as such.
"""

import itertools
from collections.abc import Hashable, Iterator

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from research import decorated_permutation as dp
from research import plabic_graph as pg
from research.le_diagram import LeDiagram
from research.strand_dynamics import Orbit, StrandAutomaton

BLACK, WHITE = pg.BLACK, pg.WHITE
Graph = pg.PlabicGraph
type AnyGraph = pg.PlabicGraph[Hashable]

MAX_STEPS = 10_000


# --------------------------------------------------------------------------- #
# Helpers and strategies
# --------------------------------------------------------------------------- #
def _derangements(k: int, n: int) -> Iterator[dp.DecoratedPermutation]:
    """Yield the fixed-point-free decorated permutations of type ``(k, n)``."""
    for targets in itertools.permutations(range(1, n + 1)):
        if all(target != i for i, target in enumerate(targets, start=1)):
            decorated = dp.DecoratedPermutation(targets)
            if decorated.anti_exceedance_count == k:
                yield decorated


def _orbits(k: int, n: int) -> list[Orbit]:
    orbits: list[Orbit] = []
    for decorated in _derangements(k, n):
        graph = Graph.from_bridge_decomposition(decorated)
        orbit = StrandAutomaton(graph).orbit(max_steps=MAX_STEPS)
        assert orbit is not None
        orbits.append(orbit)
    return orbits


def _is_noncrossing_matching(decorated: dp.DecoratedPermutation) -> bool:
    arcs = {(i, t) for i, t in enumerate(decorated.targets, start=1) if i < t}
    involution = all(
        decorated.targets[t - 1] == i for i, t in enumerate(decorated.targets, start=1)
    )
    return involution and not any(
        a < c < b < d for (a, b), (c, d) in itertools.permutations(arcs, 2)
    )


def _reference_loads(graph: AnyGraph) -> tuple[int, ...]:
    """Recompute the loads from the library's trips (edge e: darts 2e, 2e+1)."""
    load = dict.fromkeys(graph.internal_vertices, 0)
    for trip in graph.trips()[: graph.size]:
        visited = {graph.edges[dart >> 1][1 - (dart & 1)] for dart in trip}
        for v in visited & load.keys():
            load[v] += len(trip)
    return tuple(load[v] for v in graph.internal_vertices)


@st.composite
def colored_automata(draw: st.DrawFn) -> tuple[StrandAutomaton[Hashable], list[int]]:
    """Draw a bridge or Le-graph of a random cell with a random coloring."""
    n = draw(st.integers(min_value=1, max_value=6))
    targets = tuple(draw(st.permutations(range(1, n + 1))))
    fixed = [i for i, target in enumerate(targets, start=1) if target == i]
    clockwise = (
        frozenset(draw(st.sets(st.sampled_from(fixed)))) if fixed else frozenset()
    )
    decorated = dp.DecoratedPermutation(targets, clockwise)
    if draw(st.booleans()):
        graph = Graph.from_bridge_decomposition(decorated)
    else:
        graph = LeDiagram.from_decorated_permutation(decorated).to_plabic_graph()
    state = draw(
        st.lists(
            st.sampled_from([BLACK, WHITE]),
            min_size=len(graph.internal_vertices),
            max_size=len(graph.internal_vertices),
        )
    )
    return StrandAutomaton(graph), state


# --------------------------------------------------------------------------- #
# Hand-computed fixtures
# --------------------------------------------------------------------------- #
class TestSingleBridge:
    """Gr(1,2): one bridge; both trips have 3 darts and visit both vertices."""

    def test_both_vertices_carry_the_load_of_both_trips(self):
        automaton = StrandAutomaton(Graph.from_bridge_decomposition(dp.top_cell(1, 2)))
        assert automaton.loads() == (6, 6)

    def test_equal_loads_do_not_flip(self):
        automaton = StrandAutomaton(Graph.from_bridge_decomposition(dp.top_cell(1, 2)))
        assert automaton.step() == automaton.initial_state == (WHITE, BLACK)

    def test_the_orbit_is_a_fixed_point(self):
        automaton = StrandAutomaton(Graph.from_bridge_decomposition(dp.top_cell(1, 2)))
        assert automaton.orbit(max_steps=1) == Orbit(0, 1, 1, 1)

    def test_color_is_immaterial_at_a_degree_two_vertex(self):
        """Left and right of the arrival edge are the same, other, edge."""
        automaton = StrandAutomaton(Graph.from_bridge_decomposition(dp.top_cell(1, 2)))
        assert automaton.trip_permutation([WHITE, WHITE]) == (2, 1)
        assert automaton.loads([BLACK, WHITE]) == (6, 6)


class TestLollipops:
    def test_a_lollipop_carries_its_own_u_turn_and_never_flips(self):
        automaton = StrandAutomaton(pg.lollipop_graph([BLACK, WHITE, BLACK]))
        assert automaton.loads() == (2, 2, 2)
        assert automaton.orbit(max_steps=1) == Orbit(0, 1, 1, 1)


# --------------------------------------------------------------------------- #
# The library's trips as an independent oracle
# --------------------------------------------------------------------------- #
class TestAgainstLibraryTrips:
    @given(colored_automata())
    @settings(max_examples=150, deadline=None)
    def test_trip_permutation_is_that_of_the_recolored_graph(self, drawn):
        automaton, state = drawn
        recolored = automaton.recolored(state)
        assert automaton.trip_permutation(state) == recolored.trip_permutation

    @given(colored_automata())
    @settings(max_examples=150, deadline=None)
    def test_loads_sum_the_lengths_of_the_visiting_trips(self, drawn):
        automaton, state = drawn
        assert automaton.loads(state) == _reference_loads(automaton.recolored(state))

    @given(colored_automata())
    @settings(max_examples=150, deadline=None)
    def test_a_vertex_flips_iff_its_load_exceeds_its_neighbors_mean(self, drawn):
        automaton, state = drawn
        graph = automaton.recolored(state)
        load = dict(zip(graph.internal_vertices, automaton.loads(state), strict=True))
        expected = []
        for v, color in graph.colors:
            ends = [graph.edges[d >> 1][1 - (d & 1)] for d in dict(graph.rotations)[v]]
            around = [load[u] for u in ends if u in load]
            flips = bool(around) and load[v] > sum(around) / len(around)
            expected.append(-color if flips else color)
        assert automaton.step(state) == tuple(expected)


# --------------------------------------------------------------------------- #
# Orbit laws
# --------------------------------------------------------------------------- #
class TestOrbit:
    @given(colored_automata())
    @settings(max_examples=100, deadline=None)
    def test_the_orbit_describes_the_iterated_step(self, drawn):
        automaton, state = drawn
        orbit = automaton.orbit(state, max_steps=MAX_STEPS)
        assert orbit is not None
        trajectory = [tuple(state)]
        for _ in range(orbit.transient + orbit.period):
            trajectory.append(automaton.step(trajectory[-1]))
        assert trajectory[-1] == trajectory[orbit.transient]
        assert len(set(trajectory[:-1])) == len(trajectory) - 1
        cycle = [
            automaton.trip_permutation(s) for s in trajectory[orbit.transient : -1]
        ]
        assert orbit.cycle_permutations == len(set(cycle))
        assert orbit.period % orbit.permutation_period == 0
        shifted = cycle[orbit.permutation_period :] + cycle[: orbit.permutation_period]
        assert shifted == cycle

    def test_the_step_budget_counts_rule_applications(self):
        """Gr(1,3) has period 2 from the start: two steps to see the repeat."""
        automaton = StrandAutomaton(Graph.from_bridge_decomposition(dp.top_cell(1, 3)))
        assert automaton.orbit(max_steps=1) is None
        orbit = automaton.orbit(max_steps=2)
        assert orbit is not None
        assert (orbit.transient, orbit.period) == (0, 2)

    def test_a_negative_budget_is_rejected(self):
        automaton = StrandAutomaton(pg.lollipop_graph([BLACK]))
        with pytest.raises(ValueError, match="max_steps must be non-negative"):
            automaton.orbit(max_steps=-1)

    @pytest.mark.parametrize("state", [[BLACK], [BLACK, 0], [BLACK, WHITE, WHITE]])
    def test_a_state_must_color_every_internal_vertex(self, state):
        automaton = StrandAutomaton(pg.lollipop_graph([BLACK, WHITE]))
        with pytest.raises(ValueError, match="each of the 2 internal vertices"):
            automaton.step(state)


# --------------------------------------------------------------------------- #
# Finding 020 of the predecessor repository (bridge graphs, n <= 6)
# --------------------------------------------------------------------------- #
class TestFinding020:
    @pytest.mark.parametrize(
        (
            "k",
            "n",
            "cells",
            "period_two",
            "longer",
            "fixed",
            "max_period",
            "max_transient",
        ),
        [
            (1, 2, 1, 0, 0, 1, 1, 0),
            (2, 4, 7, 5, 0, 2, 2, 2),
            (2, 5, 21, 20, 1, 0, 4, 2),
            (3, 5, 21, 17, 4, 0, 6, 2),
            (2, 6, 51, 46, 5, 0, 10, 10),
            (3, 6, 161, 130, 26, 5, 26, 18),
            (4, 6, 51, 44, 7, 0, 8, 14),
        ],
    )
    def test_per_grassmannian_summary_row(
        self, k, n, cells, period_two, longer, fixed, max_period, max_transient
    ):
        orbits = _orbits(k, n)
        assert len(orbits) == cells
        assert sum(o.period == 2 for o in orbits) == period_two
        assert sum(o.period > 2 for o in orbits) == longer
        assert sum(o.period == 1 for o in orbits) == fixed
        assert max(o.period for o in orbits) == max_period
        assert max(o.transient for o in orbits) == max_transient

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_extreme_ranks_have_period_two_and_no_transient(self, n):
        for k in (1, n - 1):
            assert {(o.transient, o.period) for o in _orbits(k, n)} == {(0, 2)}

    @pytest.mark.parametrize("m", [1, 2, 3])
    def test_fixed_points_of_gr_m_2m_are_the_noncrossing_matchings(self, m):
        fixed = {
            decorated
            for decorated in _derangements(m, 2 * m)
            if StrandAutomaton(Graph.from_bridge_decomposition(decorated)).orbit(
                max_steps=1
            )
            is not None
        }
        matchings = {d for d in _derangements(m, 2 * m) if _is_noncrossing_matching(d)}
        assert fixed == matchings
