"""Tests for research.matroid, transcribed from the Logseq Matroid page.

Fixtures assert exactly what the page's canonical-example blocks certify;
property tests transcribe the structural-theorem blocks (attribution in each
docstring-free test name and comment); round-trip laws come from the API
contract; and every numbered axiom has a rejection test naming it.
"""

import itertools
from collections import Counter
from fractions import Fraction

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from experiments import io
from research import matroid as mt

matplotlib.use("Agg")


# --------------------------------------------------------------------------- #
# Strategies: small random matroids through the public constructors
# --------------------------------------------------------------------------- #
def _draw_matroid_on(draw: st.DrawFn, labels: tuple[int, ...]) -> mt.Matroid[int]:
    n = len(labels)
    kind = draw(st.sampled_from(["uniform", "graphic", "linear", "transversal"]))
    if kind == "uniform":
        r = draw(st.integers(0, n))
        return mt.Matroid.from_transversal_system(labels, [labels] * r)
    if kind == "graphic":
        vertex_count = draw(st.integers(1, 4))
        pairs = [
            (
                draw(st.integers(1, vertex_count)),
                draw(st.integers(1, vertex_count)),
            )
            for _ in labels
        ]
        return mt.Matroid.from_graph_edges(dict(zip(labels, pairs, strict=True)))
    if kind == "linear":
        p = draw(st.sampled_from([2, 3, None]))
        dim = draw(st.integers(1, 3))
        vectors = {
            label: tuple(draw(st.integers(-2, 2)) for _ in range(dim))
            for label in labels
        }
        return mt.Matroid.from_vectors(vectors, field_char=p)
    if n == 0:
        return mt.Matroid.from_transversal_system(labels, [])
    system = draw(
        st.lists(st.sets(st.sampled_from(list(labels))), min_size=0, max_size=3)
    )
    return mt.Matroid.from_transversal_system(labels, system)


@st.composite
def matroids(draw: st.DrawFn) -> mt.Matroid[int]:
    n = draw(st.integers(0, 5))
    m = _draw_matroid_on(draw, tuple(range(n)))
    op = draw(st.sampled_from(["none", "dual", "truncate"]))
    if op == "dual":
        return m.dual()
    if op == "truncate":
        return m.truncation(draw(st.integers(0, m.rank())))
    return m


@st.composite
def matroid_pairs(draw: st.DrawFn) -> tuple[mt.Matroid[int], mt.Matroid[int]]:
    n = draw(st.integers(0, 4))
    labels = tuple(range(n))
    return _draw_matroid_on(draw, labels), _draw_matroid_on(draw, labels)


def _all_subsets(m: mt.Matroid[int]) -> list[frozenset[int]]:
    elems = list(m.elements)
    return [
        frozenset(combo)
        for size in range(len(elems) + 1)
        for combo in itertools.combinations(elems, size)
    ]


def _label_family[T](m: mt.Matroid[T]) -> set[frozenset[T]]:
    elems = m.elements
    return {
        frozenset(elems[i] for i in range(len(elems)) if mask >> i & 1)
        for mask in m.independent_masks
    }


_FANO_VECTORS = {
    1: (0, 0, 1),
    2: (0, 1, 0),
    3: (1, 0, 0),
    4: (0, 1, 1),
    5: (1, 0, 1),
    6: (1, 1, 0),
    7: (1, 1, 1),
}


# --------------------------------------------------------------------------- #
# Canonical examples: each asserts what the page says the example certifies
# --------------------------------------------------------------------------- #
class TestCanonicalExamples:
    def test_u24_is_not_binary(self):
        assert mt.u24().is_binary() is False

    def test_u24_is_excluded_minor_every_proper_minor_is_binary(self):
        m = mt.u24()
        for e in m.elements:
            assert m.delete([e]).is_binary() is True
            assert m.contract([e]).is_binary() is True

    def test_u24_is_self_dual(self):
        assert mt.u24().dual() == mt.u24()

    def test_empty_matroid_is_trivial(self):
        m = mt.empty_matroid()
        assert m.size == 0
        assert m.rank() == 0
        assert m.bases == {frozenset()}
        assert m.tutte_polynomial() == {(0, 0): 1}

    def test_free_matroid_has_every_subset_independent(self):
        m = mt.free_matroid(3)
        assert _label_family(m) == set(_all_subsets(m))

    def test_loopy_matroid_makes_every_element_a_loop(self):
        m = mt.loopy_matroid(3)
        assert m.loops == {0, 1, 2}
        assert m.rank() == 0

    def test_k4_has_rank_3_on_6_elements(self):
        m = mt.k4_matroid()
        assert m.rank() == 3
        assert m.size == 6

    def test_k4_bases_are_the_16_spanning_trees(self):
        assert len(mt.k4_matroid().bases) == 16

    def test_k4_circuits_are_triangles_and_quadrilaterals(self):
        sizes = Counter(len(c) for c in mt.k4_matroid().circuits)
        assert sizes == {3: 4, 4: 3}

    def test_k4_is_graphic_hence_regular_and_binary(self):
        m = mt.k4_matroid()
        assert m.is_binary() is True
        assert m.is_regular() is True
        assert m.is_graphic() is True

    def test_fano_bases_are_all_triples_except_whitneys_seven_lines(self):
        lines = {
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
        }
        expected = {
            frozenset(t)
            for t in itertools.combinations(range(1, 8), 3)
            if frozenset(t) not in lines
        }
        assert mt.fano_matroid().bases == expected

    def test_fano_is_representable_over_gf2(self):
        linear = mt.Matroid.from_vectors(_FANO_VECTORS, field_char=2)
        assert linear == mt.fano_matroid()

    def test_fano_vectors_over_rationals_give_the_non_fano(self):
        linear = mt.Matroid.from_vectors(_FANO_VECTORS)
        assert linear == mt.non_fano_matroid()

    def test_fano_is_binary_but_not_regular(self):
        assert mt.fano_matroid().is_binary() is True
        assert mt.fano_matroid().is_regular() is False

    def test_vamos_has_rank_4_on_8_elements(self):
        m = mt.vamos_matroid()
        assert m.rank() == 4
        assert m.size == 8
        assert len(m.bases) == 65

    def test_vamos_violates_ingleton_at_the_standard_quadruple(self):
        certificate = mt.vamos_matroid().ingleton_holds_for(
            {1, 2}, {3, 4}, {5, 6}, {7, 8}
        )
        assert certificate is False

    def test_uniform_matroids_are_transversal(self):
        built = mt.Matroid.from_transversal_system(range(4), [range(4)] * 2)
        assert built == mt.uniform_matroid(2, 4)

    @pytest.mark.parametrize(
        ("n", "expected"), [(0, 1), (1, 2), (2, 4), (3, 8), (4, 17), (5, 38)]
    )
    def test_enumeration_matches_oeis_a055545_prefix(self, n, expected):
        assert len(mt.enumerate_matroids(n)) == expected


# --------------------------------------------------------------------------- #
# Structural theorems as property tests (the page's property-test oracle)
# --------------------------------------------------------------------------- #
class TestStructuralTheorems:
    @settings(max_examples=30, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_rado_edmonds_greedy_is_optimal(self, m, data):
        # Rado-Edmonds greedy characterization (Matroid page).
        weights = {
            e: data.draw(st.integers(-3, 6), label=f"w[{e}]") for e in m.elements
        }
        greedy = m.greedy_max_weight_independent(weights)
        best = max(sum(weights[e] for e in members) for members in _label_family(m))
        assert sum(weights[e] for e in greedy) == best

    def test_greedy_fails_on_a_hereditary_non_matroid(self):
        # The iff direction of Rado-Edmonds, on the fixed witness family
        # {<=ab, c} where greedy takes c (weight 3) over ab (total 4).
        family = [set(), {"a"}, {"b"}, {"c"}, {"a", "b"}]
        with pytest.raises(ValueError, match=r"\(I3\)"):
            mt.Matroid.from_independent_sets("abc", family)
        weights = {"a": 2, "b": 2, "c": 3}
        picked: set[str] = set()
        for e in sorted(weights, key=lambda e: -weights[e]):
            if any(picked | {e} <= s for s in family):
                picked |= {e}
        best = max(sum(weights[e] for e in s) for s in family)
        assert sum(weights[e] for e in picked) < best

    @settings(max_examples=40, deadline=None)
    @given(m=matroids())
    def test_duality_is_an_involution(self, m):
        # Whitney 1935: M** = M.
        assert m.dual().dual() == m

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_dual_rank_formula(self, m):
        # r*(X) = |X| + r(E - X) - r(E) (Whitney 1935).
        dual = m.dual()
        for subset in _all_subsets(m):
            assert dual.rank(subset) == m.corank(subset)

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_duality_swaps_loops_and_coloops(self, m):
        assert m.dual().loops == m.coloops
        assert m.dual().coloops == m.loops

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_duality_swaps_circuits_and_cocircuits(self, m):
        assert m.dual().circuits == m.cocircuits
        assert m.dual().cocircuits == m.circuits

    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_duality_swaps_deletion_and_contraction(self, m, data):
        if m.size == 0:
            return
        e = data.draw(st.sampled_from(list(m.elements)), label="e")
        assert m.delete([e]).dual() == m.dual().contract([e])
        assert m.contract([e]).dual() == m.dual().delete([e])

    @given(r=st.integers(0, 4), extra=st.integers(0, 3))
    def test_uniform_dual_law(self, r, extra):
        # U_{r,n}* = U_{n-r,n} (Matroid page, canonical examples).
        n = r + extra
        assert mt.uniform_matroid(r, n).dual() == mt.uniform_matroid(n - r, n)

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_circuit_cocircuit_orthogonality(self, m):
        # |C intersect C*| != 1 (Matroid page).
        for circuit in m.circuits:
            for cocircuit in m.cocircuits:
                assert len(circuit & cocircuit) != 1

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_unique_fundamental_circuit(self, m):
        if not m.bases:
            return
        basis = min(m.bases, key=sorted)
        for e in set(m.elements) - basis:
            inside = [c for c in m.circuits if c <= basis | {e}]
            assert len(inside) == 1
            assert e in inside[0]
            assert m.fundamental_circuit(e, basis) == inside[0]

    @settings(max_examples=40, deadline=None)
    @given(m=matroids())
    def test_all_bases_are_equicardinal(self, m):
        # A theorem, not an axiom (Matroid page, Definition - Bases).
        assert len({len(b) for b in m.bases}) == 1

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_basis_exchange_graph_is_connected(self, m):
        bases = sorted(m.bases, key=sorted)
        reached = {0}
        frontier = [0]
        while frontier:
            i = frontier.pop()
            for j in range(len(bases)):
                if j not in reached and len(bases[i] ^ bases[j]) == 2:
                    reached.add(j)
                    frontier.append(j)
        assert reached == set(range(len(bases)))

    @settings(max_examples=20, deadline=None)
    @given(m=matroids())
    def test_symmetric_exchange_brylawski(self, m):
        for b1, b2 in itertools.product(sorted(m.bases, key=sorted), repeat=2):
            for x in b1 - b2:
                assert any(
                    m.is_basis((b1 - {x}) | {y}) and m.is_basis((b2 - {y}) | {x})
                    for y in b2 - b1
                )

    @settings(max_examples=20, deadline=None)
    @given(m=matroids())
    def test_strong_circuit_elimination(self, m):
        circuits = sorted(m.circuits, key=sorted)
        for c1, c2 in itertools.product(circuits, repeat=2):
            if c1 == c2:
                continue
            for e in c1 & c2:
                for f in c1 - c2:
                    assert any(f in c3 and c3 <= (c1 | c2) - {e} for c3 in circuits)

    @settings(max_examples=25, deadline=None)
    @given(pair=matroid_pairs())
    def test_edmonds_intersection_min_max(self, pair):
        # Edmonds' matroid intersection theorem (Matroid page).
        m1, m2 = pair
        common = m1.max_common_independent(m2)
        ground = set(m1.elements)
        certificate = min(m1.rank(x) + m2.rank(ground - x) for x in _all_subsets(m1))
        assert len(common) == certificate

    @settings(max_examples=15, deadline=None)
    @given(pair=matroid_pairs())
    def test_nash_williams_union_rank_formula(self, pair):
        m1, m2 = pair
        union = m1.union(m2)
        for x in _all_subsets(m1):
            expected = min(
                len(x - y) + m1.rank(y) + m2.rank(y)
                for y in (
                    frozenset(c)
                    for size in range(len(x) + 1)
                    for c in itertools.combinations(sorted(x), size)
                )
            )
            assert union.rank(x) == expected

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_base_packing_certificate_tutte_nash_williams(self, m):
        # M has 2 disjoint bases iff 2 r(X) + |E - X| >= 2 r(E) for all X.
        ground = set(m.elements)
        certificate = all(
            2 * m.rank(x) + len(ground - x) >= 2 * m.rank() for x in _all_subsets(m)
        )
        packs = m.union(m).rank() == 2 * m.rank()
        assert packs == certificate

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_covering_certificate_edmonds(self, m):
        # E is a union of 2 independent sets iff |X| <= 2 r(X) for all X.
        certificate = all(len(x) <= 2 * m.rank(x) for x in _all_subsets(m))
        covers = m.union(m).is_independent(m.elements)
        assert covers == certificate

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_tutte_evaluations_count_bases_independents_spanning(self, m):
        family = _label_family(m)
        spanning = sum(1 for x in _all_subsets(m) if m.is_spanning(x))
        assert m.tutte(1, 1) == len(m.bases)
        assert m.tutte(2, 1) == len(family)
        assert m.tutte(1, 2) == spanning
        assert m.tutte(2, 2) == 2**m.size

    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_tutte_deletion_contraction_recurrence(self, m, data):
        if m.size == 0:
            return
        e = data.draw(st.sampled_from(list(m.elements)), label="e")
        t = Counter(m.tutte_polynomial())
        if e in m.loops:
            shifted = Counter(
                {
                    (p, q + 1): c
                    for (p, q), c in m.delete([e]).tutte_polynomial().items()
                }
            )
            assert t == shifted
        elif e in m.coloops:
            shifted = Counter(
                {
                    (p + 1, q): c
                    for (p, q), c in m.contract([e]).tutte_polynomial().items()
                }
            )
            assert t == shifted
        else:
            total = Counter(m.delete([e]).tutte_polynomial())
            total.update(m.contract([e]).tutte_polynomial())
            assert t == +total

    @settings(max_examples=25, deadline=None)
    @given(m=matroids())
    def test_tutte_duality_swaps_variables(self, m):
        transposed = {(q, p): c for (p, q), c in m.tutte_polynomial().items()}
        assert m.dual().tutte_polynomial() == transposed

    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), lam=st.integers(-2, 3))
    def test_characteristic_polynomial_is_signed_tutte_specialization(self, m, lam):
        coefficients = m.characteristic_polynomial()
        value = sum(c * lam**k for k, c in enumerate(coefficients))
        assert value == (-1) ** m.rank() * m.tutte(1 - lam, 0)

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_characteristic_coefficients_are_log_concave(self, m):
        # Adiprasito-Huh-Katz 2018 (Heron-Rota-Welsh conjecture).
        w = [abs(c) for c in m.characteristic_polynomial()]
        for k in range(1, len(w) - 1):
            assert w[k] * w[k] >= w[k - 1] * w[k + 1]

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_independent_set_counts_are_log_concave(self, m):
        # Mason's conjecture, strongest form (ALOGV / Branden-Huh 2020).
        counts = m.independent_set_counts()
        for k in range(1, len(counts) - 1):
            assert counts[k] * counts[k] >= counts[k - 1] * counts[k + 1]

    @settings(max_examples=25, deadline=None)
    @given(data=st.data())
    def test_ingleton_holds_for_representable_matroids(self, data):
        # One-sided, as the page states: representable => Ingleton.
        n = data.draw(st.integers(0, 5), label="n")
        vectors = {
            i: tuple(data.draw(st.integers(0, 2)) for _ in range(3)) for i in range(n)
        }
        m = mt.Matroid.from_vectors(vectors, field_char=3)
        quadruple = [
            data.draw(st.sets(st.sampled_from(range(n))), label=f"S{j}") if n else set()
            for j in range(4)
        ]
        assert m.ingleton_holds_for(*quadruple) is True

    @settings(max_examples=20, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_deletion_and_contraction_commute(self, m, data):
        if m.size < 2:
            return
        e, f = data.draw(st.permutations(list(m.elements)), label="order")[:2]
        assert m.delete([e]).contract([f]) == m.contract([f]).delete([e])
        assert m.minor(deletions=[e], contractions=[f]) == m.delete([e]).contract([f])

    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_contraction_rank_formula(self, m, data):
        # r_{M/X}(Y) = r(Y + X) - r(X) (Matroid page, contraction).
        elems = list(m.elements)
        x = data.draw(st.sets(st.sampled_from(elems)) if elems else st.just(set()))
        contracted = m.contract(x)
        for y in _all_subsets(contracted):
            assert contracted.rank(y) == m.rank(y | set(x)) - m.rank(x)


# --------------------------------------------------------------------------- #
# Cryptomorphism and serialization round trips
# --------------------------------------------------------------------------- #
class TestRoundTrips:
    @settings(max_examples=20, deadline=None)
    @given(m=matroids())
    def test_all_seven_presentations_rebuild_the_same_matroid(self, m):
        elems = m.elements
        assert mt.Matroid.from_independent_sets(elems, _label_family(m)) == m
        assert mt.Matroid.from_bases(elems, m.bases) == m
        assert mt.Matroid.from_circuits(elems, m.circuits) == m
        assert mt.Matroid.from_rank_function(elems, m.rank) == m
        assert mt.Matroid.from_closure(elems, m.closure) == m
        assert mt.Matroid.from_flats(elems, m.flats) == m
        assert mt.Matroid.from_hyperplanes(elems, m.hyperplanes) == m

    @settings(max_examples=30, deadline=None)
    @given(m=matroids())
    def test_dataframe_round_trip(self, m):
        assert mt.Matroid.from_dataframe(m.to_dataframe()) == m

    def test_dataframe_shape_for_fano(self):
        frame = mt.fano_matroid().to_dataframe()
        assert list(frame.columns) == ["element", "basis"]
        assert len(frame) == 7 + 28 * 3
        assert frame["basis"].isna().sum() == 7

    def test_dataframe_survives_experiment_io(self, tmp_path):
        fano = mt.fano_matroid()
        path = io.write_result(fano.to_dataframe(), tmp_path / "fano.json")
        decoded = mt.Matroid.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == fano

    def test_experiment_io_keeps_numeric_looking_string_labels(self, tmp_path):
        k4 = mt.k4_matroid()
        path = io.write_result(k4.to_dataframe(), tmp_path / "k4.json")
        decoded = mt.Matroid.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == k4
        assert decoded.ground_set == k4.ground_set

    def test_empty_matroid_survives_experiment_io(self, tmp_path):
        empty = mt.empty_matroid()
        path = io.write_result(empty.to_dataframe(), tmp_path / "empty.json")
        decoded = mt.Matroid.from_dataframe(pd.read_json(path, dtype=False))
        assert decoded == empty

    def test_equality_ignores_element_order(self):
        subsets = [set(), {"a"}, {"b"}, {"a", "b"}]
        forward = mt.Matroid.from_independent_sets("ab", subsets)
        backward = mt.Matroid.from_independent_sets("ba", subsets)
        assert forward == backward
        assert hash(forward) == hash(backward)


# --------------------------------------------------------------------------- #
# Operations and derived vocabulary
# --------------------------------------------------------------------------- #
class TestOperations:
    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), k=st.integers(0, 5))
    def test_truncation_keeps_small_independent_sets(self, m, k):
        truncated = m.truncation(k)
        assert _label_family(truncated) == {s for s in _label_family(m) if len(s) <= k}

    @settings(max_examples=25, deadline=None)
    @given(m=matroids(), data=st.data())
    def test_restriction_keeps_inside_independent_sets(self, m, data):
        elems = list(m.elements)
        x = data.draw(st.sets(st.sampled_from(elems)) if elems else st.just(set()))
        restricted = m.restrict(x)
        assert restricted.ground_set == frozenset(x)
        assert _label_family(restricted) == {
            s for s in _label_family(m) if s <= frozenset(x)
        }

    @settings(max_examples=20, deadline=None)
    @given(pair=matroid_pairs())
    def test_direct_sum_of_nonempty_parts_is_disconnected(self, pair):
        m1, m2 = pair
        if m1.size == 0:
            return
        shifted = mt.Matroid.from_independent_sets(
            [f"b{e}" for e in m2.elements],
            [[f"b{e}" for e in s] for s in _label_family(m2)],
            validate=False,
        )
        total = m1.direct_sum(shifted)
        assert total.size == m1.size + m2.size
        assert total.is_connected is False

    def test_direct_sum_rejects_overlapping_ground_sets(self):
        with pytest.raises(ValueError, match="disjoint ground sets"):
            mt.u24().direct_sum(mt.uniform_matroid(1, 2))

    def test_union_rejects_different_ground_sets(self):
        with pytest.raises(ValueError, match="same ground set"):
            mt.uniform_matroid(1, 2).union(mt.uniform_matroid(1, 3))

    def test_minor_rejects_overlapping_sets(self):
        with pytest.raises(ValueError, match="disjoint"):
            mt.u24().minor(deletions=[0], contractions=[0])

    def test_graph_loops_and_parallel_edges(self):
        m = mt.Matroid.from_graph_edges(
            {"e1": (1, 2), "e2": (1, 2), "e3": (2, 3), "ell": (1, 1)}
        )
        assert m.loops == {"ell"}
        assert {"e1", "e2"} in m.parallel_classes
        assert m.is_simple is False
        simplified = m.simplification()
        assert simplified.is_simple is True
        assert simplified.size == 2

    def test_simplification_keeps_one_element_per_parallel_class(self):
        m = mt.Matroid.from_graph_edges(
            {"e1": (1, 2), "e2": (1, 2), "e3": (2, 3), "ell": (1, 1)}
        )
        assert m.simplification().ground_set == {"e1", "e3"}

    def test_fundamental_circuit_in_k4(self):
        m = mt.k4_matroid()
        assert m.fundamental_circuit("23", {"12", "13", "14"}) == {
            "12",
            "13",
            "23",
        }

    def test_fundamental_circuit_rejects_non_basis(self):
        with pytest.raises(ValueError, match="not a basis"):
            mt.k4_matroid().fundamental_circuit("23", {"12", "13"})

    def test_fundamental_circuit_rejects_element_inside_basis(self):
        with pytest.raises(ValueError, match="lies in the basis"):
            mt.k4_matroid().fundamental_circuit("12", {"12", "13", "14"})

    def test_u24_flats_and_hyperplanes(self):
        m = mt.u24()
        singletons = {frozenset({e}) for e in range(4)}
        assert m.flats == {frozenset()} | singletons | {frozenset(range(4))}
        assert m.hyperplanes == singletons
        assert m.cocircuits == {frozenset(range(4)) - s for s in singletons}

    def test_u24_tutte_polynomial_exact(self):
        assert mt.u24().tutte_polynomial() == {
            (2, 0): 1,
            (1, 0): 2,
            (0, 1): 2,
            (0, 2): 1,
        }

    def test_u24_characteristic_polynomial_exact(self):
        assert mt.u24().characteristic_polynomial() == (3, -4, 1)

    def test_greedy_requires_a_weight_for_every_element(self):
        with pytest.raises(KeyError):
            mt.u24().greedy_max_weight_independent({0: 1})

    def test_repr_is_compact(self):
        assert repr(mt.u24()) == "Matroid(n=4, rank=2, bases=6)"


# --------------------------------------------------------------------------- #
# Axiom violations: every numbered axiom rejects with its name
# --------------------------------------------------------------------------- #
class TestAxiomViolations:
    def test_i1_empty_set_must_be_independent(self):
        with pytest.raises(ValueError, match=r"\(I1\)"):
            mt.Matroid.from_independent_sets("a", [{"a"}])

    def test_i2_hereditary(self):
        with pytest.raises(ValueError, match=r"\(I2\)"):
            mt.Matroid.from_independent_sets("ab", [set(), {"a", "b"}])

    def test_i3_augmentation(self):
        family = [set(), {"a"}, {"b"}, {"c"}, {"b", "c"}]
        with pytest.raises(ValueError, match=r"\(I3\)"):
            mt.Matroid.from_independent_sets("abc", family)

    def test_b1_nonempty(self):
        with pytest.raises(ValueError, match=r"\(B1\)"):
            mt.Matroid.from_bases("ab", [])

    def test_b2_basis_exchange(self):
        with pytest.raises(ValueError, match=r"\(B2\)"):
            mt.Matroid.from_bases("abcd", [{"a", "b"}, {"c", "d"}])

    def test_c1_no_empty_circuit(self):
        with pytest.raises(ValueError, match=r"\(C1\)"):
            mt.Matroid.from_circuits("a", [set()])

    def test_c2_antichain(self):
        with pytest.raises(ValueError, match=r"\(C2\)"):
            mt.Matroid.from_circuits("ab", [{"a"}, {"a", "b"}])

    def test_c3_circuit_elimination(self):
        with pytest.raises(ValueError, match=r"\(C3\)"):
            mt.Matroid.from_circuits("abc", [{"a", "b"}, {"b", "c"}])

    def test_r1_bounded_by_cardinality(self):
        with pytest.raises(ValueError, match=r"\(R1\)"):
            mt.Matroid.from_rank_function("a", lambda s: 2 * len(s))

    def test_r2_monotone(self):
        table = {
            frozenset(): 0,
            frozenset("a"): 1,
            frozenset("b"): 1,
            frozenset("ab"): 0,
        }
        with pytest.raises(ValueError, match=r"\(R2\)"):
            mt.Matroid.from_rank_function("ab", lambda s: table[frozenset(s)])

    def test_r3_submodular(self):
        table = {
            frozenset(): 0,
            frozenset("a"): 1,
            frozenset("b"): 1,
            frozenset("c"): 1,
            frozenset("ab"): 1,
            frozenset("ac"): 1,
            frozenset("bc"): 2,
            frozenset("abc"): 2,
        }
        with pytest.raises(ValueError, match=r"\(R3\)"):
            mt.Matroid.from_rank_function("abc", lambda s: table[frozenset(s)])

    def test_cl1_extensive(self):
        with pytest.raises(ValueError, match=r"\(CL1\)"):
            mt.Matroid.from_closure("a", lambda s: frozenset())

    def test_cl2_monotone(self):
        table = {
            frozenset(): frozenset("a"),
            frozenset("a"): frozenset("a"),
            frozenset("b"): frozenset("b"),
            frozenset("ab"): frozenset("ab"),
        }
        with pytest.raises(ValueError, match=r"\(CL2\)"):
            mt.Matroid.from_closure("ab", lambda s: table[frozenset(s)])

    def test_cl3_idempotent(self):
        table = {
            frozenset(): frozenset("a"),
            frozenset("a"): frozenset("ab"),
            frozenset("b"): frozenset("ab"),
            frozenset("ab"): frozenset("ab"),
        }
        with pytest.raises(ValueError, match=r"\(CL3\)"):
            mt.Matroid.from_closure("ab", lambda s: table[frozenset(s)])

    def test_cl4_mac_lane_steinitz_exchange(self):
        table = {
            frozenset(): frozenset(),
            frozenset("a"): frozenset("ab"),
            frozenset("b"): frozenset("b"),
            frozenset("ab"): frozenset("ab"),
        }
        with pytest.raises(ValueError, match=r"\(CL4\)"):
            mt.Matroid.from_closure("ab", lambda s: table[frozenset(s)])

    def test_f1_ground_set_is_a_flat(self):
        with pytest.raises(ValueError, match=r"\(F1\)"):
            mt.Matroid.from_flats("a", [set()])

    def test_f2_intersection_closed(self):
        with pytest.raises(ValueError, match=r"\(F2\)"):
            mt.Matroid.from_flats("abc", [{"a", "b", "c"}, {"a"}, {"b"}])

    def test_f3_covering(self):
        with pytest.raises(ValueError, match=r"\(F3\)"):
            mt.Matroid.from_flats("ab", [{"a", "b"}, set(), {"a"}])

    def test_h1_ground_set_is_not_a_hyperplane(self):
        with pytest.raises(ValueError, match=r"\(H1\)"):
            mt.Matroid.from_hyperplanes("a", [{"a"}])

    def test_h2_antichain(self):
        with pytest.raises(ValueError, match=r"\(H2\)"):
            mt.Matroid.from_hyperplanes("abc", [{"a"}, {"a", "b"}])

    def test_h3_hyperplane_axiom(self):
        with pytest.raises(ValueError, match=r"\(H3\)"):
            mt.Matroid.from_hyperplanes("abc", [{"a"}, {"b"}])


# --------------------------------------------------------------------------- #
# Boundary validation beyond the axioms
# --------------------------------------------------------------------------- #
class TestBoundaryValidation:
    def test_duplicate_ground_set_labels_rejected(self):
        with pytest.raises(ValueError, match="distinct"):
            mt.Matroid.from_independent_sets("aa", [set()])

    def test_unknown_label_rejected(self):
        with pytest.raises(ValueError, match="not in the ground set"):
            mt.Matroid.from_independent_sets("a", [set(), {"z"}])

    def test_from_vectors_rejects_ragged_dimensions(self):
        with pytest.raises(ValueError, match="dimension"):
            mt.Matroid.from_vectors({1: (1, 0), 2: (1,)})

    def test_from_vectors_rejects_composite_characteristic(self):
        with pytest.raises(ValueError, match="prime"):
            mt.Matroid.from_vectors({1: (1,)}, field_char=4)

    def test_from_vectors_reduces_fractions_mod_p(self):
        vectors: dict[str, tuple[Fraction | int, ...]] = {
            "a": (Fraction(1, 2),),
            "b": (3,),
        }
        m = mt.Matroid.from_vectors(vectors, field_char=5)
        assert m.rank(["a"]) == 1
        assert m.rank(["a", "b"]) == 1

    def test_from_vectors_rejects_fraction_without_gf_image(self):
        with pytest.raises(ValueError, match=r"no image in GF\(3\)"):
            mt.Matroid.from_vectors({"a": (Fraction(1, 3),)}, field_char=3)

    def test_uniform_matroid_rejects_bad_rank(self):
        with pytest.raises(ValueError, match="0 <= rank <= n"):
            mt.uniform_matroid(3, 2)

    def test_truncation_rejects_negative_size(self):
        with pytest.raises(ValueError, match="non-negative"):
            mt.u24().truncation(-1)

    def test_enumerate_rejects_negative_size(self):
        with pytest.raises(ValueError, match="non-negative"):
            mt.enumerate_matroids(-1)

    def test_from_dataframe_rejects_missing_columns(self):
        with pytest.raises(ValueError, match="missing required columns"):
            mt.Matroid.from_dataframe(pd.DataFrame({"element": [1]}))


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
class TestPlots:
    def test_plot_lattice_of_flats_draws_on_provided_axes(self):
        _, ax = plt.subplots()
        try:
            returned = mt.u24().plot_lattice_of_flats(ax)
            assert returned is ax
            assert len(ax.collections) == 1
            assert len(ax.texts) == 6
        finally:
            plt.close("all")

    def test_plot_lattice_of_flats_creates_axes_when_omitted(self):
        try:
            ax = mt.u24().plot_lattice_of_flats()
            assert ax.get_title() == "Lattice of flats"
        finally:
            plt.close("all")

    def test_plot_basis_exchange_graph_draws_all_bases(self):
        _, ax = plt.subplots()
        try:
            returned = mt.u24().plot_basis_exchange_graph(ax)
            assert returned is ax
            assert len(ax.texts) == 6
        finally:
            plt.close("all")
