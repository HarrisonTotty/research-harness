/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Algebra.MvPolynomial.Eval
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Combinatorics.Matroid.Constructions
import Mathlib.Combinatorics.Quiver.Path.Weight
import Mathlib.LinearAlgebra.Matrix.Transvection
import Mathlib.Logic.Equiv.Fin.Rotate
import Mathlib.RingTheory.PowerSeries.Basic
import Mathlib.Tactic
import Theorems.Positroid

/-!
# The boundary measurement map

Postnikov's boundary measurement map sends a planar directed network in a disk, with positive
weights on its edges, to a point of the Grassmannian: for each boundary source `bᵢ` and boundary
sink `bⱼ` it records a signed weighted sum over the directed paths from `bᵢ` to `bⱼ`, and reads
those numbers as ratios of Plücker coordinates (Postnikov, *Total positivity, Grassmannians, and
networks*, 2006, §4). Mathlib has none of the planar-network side of this theory — no planar
graphs, no winding index, no flows, no plabic graphs — so this module formalizes the pieces of
the theory that live on objects Mathlib (or this repository) already has: matrices, their
maximal minors, quiver path weights, polynomials, and matroids.

## Main statements

* `MvPolynomial.eval_nonneg` / `MvPolynomial.eval_pos` — a polynomial with nonnegative
  coefficients evaluates nonnegatively at nonnegative points, and (when nonzero) positively at
  positive points. This is the evaluation half of Postnikov's Definition 4.4: a subtraction-free
  rational expression (a quotient of polynomials with positive coefficients) specializes at
  positive weights to a well-defined nonnegative real, because its denominator cannot vanish.
* `PowerSeries.mk_neg_one_pow_mul_one_add_eq_one` — the alternating geometric series
  `1 - X + X² - ⋯` sums against `1 + X` to `1`. This is the formal content of Postnikov's
  Example 4.5, the smallest witness that a cyclic network needs the subtraction-free rational
  form: the raw series `M₁₂ = xyt - xyzyt + ⋯` is infinite and alternating, the answer
  `xyt/(1 + yz)` finite.
* `Quiver.Path.weight_conj` / `Quiver.Path.weight_conj_of_eq_one` — reweighting the edges of a
  quiver by a vertex potential, `w' e = t a * w e * (t b)⁻¹` for `e : a ⟶ b`, rescales the
  weight of a path `p : Path i j` by `t i * ⬝ * (t j)⁻¹`; when the potential is `1` at both
  endpoints the path weight is unchanged. This is Postnikov's gauge transformation (eq. (4.2))
  and its invariance claim, at the level of a single path term of the boundary measurement.
* `Matrix.maximalMinor_pair` — a maximal minor of a two-row matrix on a sorted column pair is
  the familiar `2 × 2` determinant.
* `Matrix.plucker_fin_two` — the three-term Grassmann–Plücker relation
  `Δ_{ij}Δ_{kl} - Δ_{ik}Δ_{jl} + Δ_{il}Δ_{jk} = 0` for the maximal minors of any two-row
  matrix. Postnikov notes (after his Lemma 5.1) that the `r = 2` case of his Plücker identity
  is equivalent to this relation.
* `Matrix.maximalMinor_mul_transvection_of_notMem` /
  `Matrix.maximalMinor_mul_transvection_of_mem_of_mem` /
  `Matrix.maximalMinor_mul_transvection_of_notMem_of_mem` — how right multiplication by the
  transvection `xᵢ(a) = 1 + a·Eᵢⱼ` (for adjacent columns `j = i + 1`) acts on maximal minors:
  minors not selecting column `j` (or selecting both `i` and `j`) are unchanged, and
  `Δ_I ↦ Δ_I + a·Δ_{(I∖{j})∪{i}}` when `j ∈ I ∌ i`. This is the matrix side of Lam's
  Lemma 7.6 (*Totally nonnegative Grassmannian and Grassmann polytopes*, 2015): adding a bridge
  to a network multiplies its boundary measurement point by a Chevalley generator.
* `Matroid.isPositroid_uniqueBaseOn` — for any `I ⊆ [n]`, the matroid on `[n]` whose unique
  base is `I` is a positroid, realized by the 0/1 matrix with identity in columns `I`. These
  are the torus-fixed points of the totally nonnegative Grassmannian — the zero-dimensional
  cells that Lam's lollipop example (Example 4.2) certifies.

## Canonical examples

* `BoundaryMeasurement.squareMatrix` — Lam's `Gr(2,4)` square-graph fixture (Example 4.3): the
  matrix `!![b, 1, c, 0; -a, 0, d, 1]` whose six maximal minors are `a, ac+bd, b, d, 1, c`,
  with the Plücker relation, positivity of all minors at positive weights (the top cell), full
  rank, and the resulting base family of the uniform matroid `U_{2,4}` (reusing
  `Matroid.unifOn`). The entry `-a` is Postnikov's `(-1)^s` sign from Definition 4.6.
* `BoundaryMeasurement.squareTrip` — the trip permutation `(3,4,1,2)` of the square graph
  (0-indexed: `i ↦ i + 2`) as a `DecoratedPermutation 4` with two weak excedances, matching
  the rank of the cell.
* `BoundaryMeasurement.lollipopMatrix` — Lam's lollipop fixture (Example 4.2): the matrix
  realizing `span(e₃, e₄) ∈ Gr(2,4)`, whose only nonvanishing maximal minor is at `{3,4}`
  (0-indexed `{2,3}`); with `Matroid.isPositroid_uniqueBaseOn` this certifies the
  zero-dimensional cells.
* `BoundaryMeasurement.fourVertexMatrix` — Postnikov's Example 4.7 boundary measurement matrix
  for `I = {1,3}` in `[4]`, certifying the placement of the `(-1)^s` signs: exactly one entry
  is negated and every exchange minor `Δ_{(I∖{i})∪{j}}` comes out to the *unsigned*
  measurement `M_{ij}`, together with the companion minor
  `Δ_{24} = M₁₂M₃₄ + M₁₄M₃₂` (the smallest instance of Postnikov's Proposition 5.2).
  Note: Postnikov's printed Example 4.7 contains an uncorrected typo — it prints
  `M₁₄ = Δ₂₄/Δ₁₃` where his own Definition 4.6 and §5 force `M₁₄ = Δ₃₄/Δ₁₃`; the statements
  here prove the corrected relations (see the coverage notes for the evidence trail).
* `BoundaryMeasurement.talaskaMatrix` — Talaska's `n = 5` analogue (Example 2.7 of *A formula
  for Plücker coordinates associated with a planar network*, 2008), certifying the `(-1)^s`
  rule at `n = 5`: one negated entry, at `(1,5)`, and `Δ_{45} = M₁₅`.
* `BoundaryMeasurement.fiveBasisMatroid` — Postnikov's Example 11.9: the rank-2 positroid on
  `[4]` whose bases are all `2`-subsets except `{2,3}` (0-indexed `{1,2}`), constructed from
  its basis list and proved to be a positroid via the realizing matrix
  `!![1, 1, 1, 0; 0, 1, 1, 1]`.

## Implementation notes

Mathlib (pinned revision) covers the substrates used here — `Matrix.transvection` with its
column action, `Quiver.Path.weight`, `MvPolynomial.eval`, `PowerSeries`, geometric sums
(`geom_sum_mul_neg`), `Matroid.uniqueBaseOn`, `Finset.orderEmbOfFin` — but has none of:
planar graphs or disk embeddings, a winding index, network flows, perfect orientations,
plabic graphs, the Lindström–Gessel–Viennot lemma, total positivity of matrices, Plücker
relations, or a subspace Grassmannian over `ℝ` (its `Module.Grassmannian` uses the
quotient-rank convention and has no topology). Consequently the map itself — and with it
Postnikov's Lemma 4.3, Theorems 4.8/4.10/4.11/6.5/10.1/12.1/12.7, Talaska's flow formula,
Lam's matching formula, the twist map, and the Marsh–Rietsch comparison — is out of reach of a
faithful formalization today and is recorded as backlog in the research notes, not stated here.

`Matrix.maximalMinor`, `Matroid.IsPositroid`, `Matroid.unifOn`, and `DecoratedPermutation`
come from `Theorems.Positroid`. The transvection lemmas take the exchanged-set cardinality
`(insert i (I.erase j)).card = d` as a hypothesis rather than deriving it, matching
`Matrix.maximalMinor`'s dependent signature at use sites. The adjacency hypothesis
`(i : ℕ) + 1 = j` is stated on natural values to avoid the wrap-around of `Fin` addition; it
is what makes the exchange minor appear with coefficient `+a` (a non-adjacent exchange picks
up a sign and is not stated here).

Concrete `Finset` minors are evaluated through `Matrix.maximalMinor_pair` and
`Finset.orderEmbOfFin_unique`, never through `decide` on `Finset.sort` (which does not reduce
in the kernel; see `Theorems.Positroid`).

## References

* [A. Postnikov, *Total positivity, Grassmannians, and networks*][postnikov2006],
  arXiv:math/0609764. §4 (Definitions 4.1/4.4/4.6, Examples 4.5/4.7), §5, §11 (Example 11.9),
  eq. (4.2) (gauge transformations).
* [T. Lam, *Totally nonnegative Grassmannian and Grassmann polytopes*][lam2015],
  arXiv:1506.00603. Examples 4.2/4.3, §7 (trips, Lemma 7.6).
* [K. Talaska, *A formula for Plücker coordinates associated with a planar network*][talaska2008],
  arXiv:0801.4822. Example 2.7.

## Tags

boundary measurement, total positivity, positroid, transvection, gauge transformation,
Plücker relation, planar network
-/

open Finset

/-! ### Positive-coefficient polynomials evaluate positively

The evaluation half of Postnikov's Definition 4.4: subtraction-freeness plus `xₑ > 0` means no
denominator vanishes, so boundary measurements are well-defined nonnegative reals. -/

namespace MvPolynomial

/-- A multivariate polynomial with nonnegative coefficients evaluates to a nonnegative value at
a nonnegative point. Together with `MvPolynomial.eval_pos` this is why a subtraction-free
rational expression — "a quotient of two polynomial expressions with positive coefficients"
(Postnikov, 2006, §4) — specializes at positive weights to a well-defined nonnegative real
(Postnikov, 2006, Definition 4.4). -/
theorem eval_nonneg {σ R : Type*} [CommSemiring R] [PartialOrder R] [IsOrderedRing R]
    {p : MvPolynomial σ R} (hp : ∀ d, 0 ≤ p.coeff d) {x : σ → R} (hx : ∀ i, 0 ≤ x i) :
    0 ≤ eval x p := by
  rw [eval_eq]
  exact Finset.sum_nonneg fun d _ ↦
    mul_nonneg (hp d) (Finset.prod_nonneg fun i _ ↦ pow_nonneg (hx i) _)

/-- A nonzero multivariate polynomial with nonnegative coefficients evaluates to a strictly
positive value at a strictly positive point. This is the denominator half of Postnikov's
Definition 4.4 (2006): the denominator of a subtraction-free expression cannot vanish at
positive edge weights. -/
theorem eval_pos {σ R : Type*} [CommSemiring R] [LinearOrder R] [IsStrictOrderedRing R]
    {p : MvPolynomial σ R} (hp : ∀ d, 0 ≤ p.coeff d) (hp0 : p ≠ 0) {x : σ → R}
    (hx : ∀ i, 0 < x i) : 0 < eval x p := by
  rw [eval_eq]
  obtain ⟨d₀, hd₀⟩ : p.support.Nonempty :=
    Finset.nonempty_iff_ne_empty.2 fun h ↦ hp0 (support_eq_empty.1 h)
  refine Finset.sum_pos'
    (fun d _ ↦ mul_nonneg (hp d) (Finset.prod_nonneg fun i _ ↦ pow_nonneg (hx i).le _))
    ⟨d₀, hd₀, ?_⟩
  have hc : 0 < p.coeff d₀ :=
    lt_of_le_of_ne (hp d₀) (Ne.symm (mem_support_iff.1 hd₀))
  exact mul_pos hc (Finset.prod_pos fun i _ ↦ pow_pos (hx i) _)

end MvPolynomial

/-! ### The alternating geometric series

Postnikov's Example 4.5: in a network with one cycle of weight `u = yz`, the formal boundary
measurement is the alternating series `xyt·(1 - u + u² - ⋯)`, which sums against `1 + u` to
`xyt`. The finite telescoping identity the research page verifies numerically is Mathlib's
`geom_sum_mul_neg`, pinned by the `example` below. -/

namespace PowerSeries

/-- The alternating geometric series `1 - X + X² - ⋯` is the inverse of `1 + X`: formally,
`(∑ₙ (-1)ⁿ Xⁿ) * (1 + X) = 1`. This is the summed form of the boundary measurement series of
Postnikov's Example 4.5 (2006), the smallest witness that a cyclic network's boundary
measurement is an infinite alternating series with a finite subtraction-free value. -/
theorem mk_neg_one_pow_mul_one_add_eq_one {R : Type*} [Ring R] :
    (mk fun n ↦ (-1 : R) ^ n) * (1 + X) = 1 := by
  ext n
  cases n with
  | zero => simp [mul_add]
  | succ n =>
      rw [mul_add, mul_one, map_add, coeff_succ_mul_X, coeff_mk, coeff_mk, coeff_one]
      simp [pow_succ]

end PowerSeries

/-- Postnikov's telescoping certificate for Example 4.5 (2006), as the research page records
it: `(∑_{m=0}^{M} (-u)^m)(1 + u) = 1 - (-u)^{M+1}`. This is Mathlib's `geom_sum_mul_neg`
specialized at `x := -u`; pinned here so a change in the library form is caught. -/
example (u : ℝ) (M : ℕ) :
    (∑ m ∈ Finset.range (M + 1), (-u) ^ m) * (1 + u) = 1 - (-u) ^ (M + 1) := by
  simpa [sub_neg_eq_add] using geom_sum_mul_neg (-u) (M + 1)

/-- At unit weights `x = y = z = t = 1`, Postnikov's Example 4.5 (2006) boundary measurement
`M₁₂ = xyt/(1 + yz)` specializes to `1/2` — a positive rational, not a formal series. -/
example : (1 : ℚ) * 1 * 1 / (1 + 1 * 1) = 1 / 2 := by norm_num

/-! ### Gauge transformations act on path weights by a vertex potential

Postnikov's eq. (4.2): pick `t_v` at each vertex and replace the weight of each edge
`e : u ⟶ v` by `x'_e = t_u · x_e · t_v⁻¹`. The weight of a path from `i` to `j` then changes
by `t_i · ⬝ · t_j⁻¹` — so if the potential is trivial at the boundary, every path weight, and
hence every term of a boundary measurement, is unchanged. -/

namespace Quiver.Path

universe v u

variable {V : Type u} [Quiver.{v} V] {G : Type*} [Group G]

/-- Reweighting the edges of a quiver by a vertex potential `t` — replacing the weight of
`e : a ⟶ b` by `t a * w e * (t b)⁻¹` — rescales the weight of a path `p : Path i j` to
`t i * p.weight w * (t j)⁻¹`: the interior potentials telescope away. This is Postnikov's
gauge transformation of a network (2006, eq. (4.2)) acting on a single path. Stated for an
arbitrary group so that it applies to positive real weights (`{x : ℝ // 0 < x}`) as well as to
units of a ring. -/
theorem weight_conj (t : V → G) (w : ∀ {a b : V}, (a ⟶ b) → G) {i j : V} (p : Path i j) :
    p.weight (fun {a b} e ↦ t a * w e * (t b)⁻¹) = t i * p.weight w * (t j)⁻¹ := by
  induction p with
  | nil => simp
  | cons q e ih =>
      rw [weight_cons, weight_cons, ih]
      group

/-- A gauge transformation whose potential is trivial at both endpoints of a path preserves the
path's weight. For boundary-to-boundary paths with `t = 1` on boundary vertices this is
Postnikov's claim (2006, after eq. (4.2)) that gauge transformations preserve all boundary
measurements, term by term. -/
theorem weight_conj_of_eq_one (t : V → G) (w : ∀ {a b : V}, (a ⟶ b) → G) {i j : V}
    (p : Path i j) (hi : t i = 1) (hj : t j = 1) :
    p.weight (fun {a b} e ↦ t a * w e * (t b)⁻¹) = p.weight w := by
  rw [weight_conj, hi, hj]
  simp

end Quiver.Path

/-! ### Maximal minors of two-row matrices and the three-term Plücker relation -/

namespace Matrix

variable {R : Type*} [CommRing R] {n : ℕ}

/-- The maximal minor of a two-row matrix on the column pair `{i, j}` with `i < j` is the
`2 × 2` determinant `A 0 i * A 1 j - A 0 j * A 1 i`. This evaluation is how the concrete
boundary measurement fixtures below are computed. -/
theorem maximalMinor_pair (A : Matrix (Fin 2) (Fin n) R) {i j : Fin n} (hij : i < j) :
    A.maximalMinor {i, j} (Finset.card_pair_eq_two_iff.2 hij.ne) =
      A 0 i * A 1 j - A 0 j * A 1 i := by
  have hemb : ⇑(({i, j} : Finset (Fin n)).orderEmbOfFin
      (Finset.card_pair_eq_two_iff.2 hij.ne)) = ![i, j] := by
    symm
    refine Finset.orderEmbOfFin_unique _ (fun x ↦ ?_) ?_
    · fin_cases x <;> simp
    · intro x y hxy
      fin_cases x <;> fin_cases y <;>
        first
          | exact absurd hxy (by decide)
          | simpa using hij
  unfold maximalMinor
  rw [det_fin_two]
  simp [hemb]

/-- The three-term Grassmann–Plücker relation for a two-row matrix: for columns
`i < j < k < l`, `Δ_{ij}Δ_{kl} - Δ_{ik}Δ_{jl} + Δ_{il}Δ_{jk} = 0`. Postnikov (2006, after
Lemma 5.1) notes that the `r = 2` case of his Plücker identity is equivalent to a three-term
Grassmann–Plücker relation; this is that relation for the rows of a `2 × n` matrix. -/
theorem plucker_fin_two (A : Matrix (Fin 2) (Fin n) R) {i j k l : Fin n} (hij : i < j)
    (hjk : j < k) (hkl : k < l) :
    A.maximalMinor {i, j} (Finset.card_pair_eq_two_iff.2 hij.ne) *
        A.maximalMinor {k, l} (Finset.card_pair_eq_two_iff.2 hkl.ne) -
      A.maximalMinor {i, k} (Finset.card_pair_eq_two_iff.2 (hij.trans hjk).ne) *
        A.maximalMinor {j, l} (Finset.card_pair_eq_two_iff.2 (hjk.trans hkl).ne) +
      A.maximalMinor {i, l} (Finset.card_pair_eq_two_iff.2 ((hij.trans hjk).trans hkl).ne) *
        A.maximalMinor {j, k} (Finset.card_pair_eq_two_iff.2 hjk.ne) = 0 := by
  rw [A.maximalMinor_pair hij, A.maximalMinor_pair hkl, A.maximalMinor_pair (hij.trans hjk),
    A.maximalMinor_pair (hjk.trans hkl), A.maximalMinor_pair ((hij.trans hjk).trans hkl),
    A.maximalMinor_pair hjk]
  ring

/-! ### Bridges act on maximal minors as Chevalley generators

The matrix side of Lam's Lemma 7.6 (2015): appending a bridge of weight `a` at `i` multiplies
the boundary measurement point on the right by the transvection `xᵢ(a) = 1 + a·E_{i,i+1}`,
which fixes every maximal minor not selecting column `i + 1` (or selecting both columns) and
sends `Δ_I ↦ Δ_I + a·Δ_{(I∖{i+1})∪{i}}` when `i + 1 ∈ I ∌ i`. -/

variable {d : ℕ}

/-- Right multiplication by a transvection `1 + a·Eᵢⱼ` does not change maximal minors that do
not select column `j`: only column `j` of the product differs from `A`. Together with
`Matrix.maximalMinor_mul_transvection_of_mem_of_mem` and
`Matrix.maximalMinor_mul_transvection_of_notMem_of_mem`, this is the matrix half of Lam's
Lemma 7.6 (2015) on adding a bridge to a network. -/
theorem maximalMinor_mul_transvection_of_notMem (A : Matrix (Fin d) (Fin n) R) (i : Fin n)
    (a : R) {I : Finset (Fin n)} (hI : I.card = d) {j : Fin n} (hj : j ∉ I) :
    (A * transvection i j a).maximalMinor I hI = A.maximalMinor I hI := by
  unfold maximalMinor
  congr 1
  ext r p
  simp only [submatrix_apply, id_eq]
  exact mul_transvection_apply_of_ne i j r (I.orderEmbOfFin hI p)
    (fun h ↦ hj (h ▸ I.orderEmbOfFin_mem hI p)) a A

/-- Right multiplication by a transvection `1 + a·Eᵢⱼ` does not change maximal minors that
select both column `i` and column `j`: inside the minor, the transvection adds a multiple of
one selected column to another, which fixes the determinant. Part of the matrix half of Lam's
Lemma 7.6 (2015). -/
theorem maximalMinor_mul_transvection_of_mem_of_mem (A : Matrix (Fin d) (Fin n) R)
    {i j : Fin n} (hij : i ≠ j) (a : R) {I : Finset (Fin n)} (hI : I.card = d) (hi : i ∈ I)
    (hj : j ∈ I) :
    (A * transvection i j a).maximalMinor I hI = A.maximalMinor I hI := by
  unfold maximalMinor
  obtain ⟨q, hq⟩ : ∃ q, I.orderEmbOfFin hI q = i := by
    have hmem : i ∈ Set.range (I.orderEmbOfFin hI) := by
      rw [Finset.range_orderEmbOfFin]
      exact Finset.mem_coe.2 hi
    exact hmem
  obtain ⟨p, hp⟩ : ∃ p, I.orderEmbOfFin hI p = j := by
    have hmem : j ∈ Set.range (I.orderEmbOfFin hI) := by
      rw [Finset.range_orderEmbOfFin]
      exact Finset.mem_coe.2 hj
    exact hmem
  have hpq : p ≠ q := fun h ↦ hij (by rw [← hq, ← hp, h])
  have hsub : (A * transvection i j a).submatrix id (I.orderEmbOfFin hI) =
      (A.submatrix id (I.orderEmbOfFin hI)).updateCol p
        (fun r ↦ A.submatrix id (I.orderEmbOfFin hI) r p +
          a • A.submatrix id (I.orderEmbOfFin hI) r q) := by
    ext r p'
    by_cases hp' : p' = p
    · subst hp'
      simp only [updateCol_self, submatrix_apply, id_eq, hp, hq, smul_eq_mul]
      exact mul_transvection_apply_same i j r a A
    · simp only [updateCol_ne hp', submatrix_apply, id_eq]
      exact mul_transvection_apply_of_ne i j r (I.orderEmbOfFin hI p')
        (fun h ↦ hp' ((I.orderEmbOfFin hI).injective (h.trans hp.symm))) a A
  rw [hsub, det_updateCol_add_smul_self _ hpq a]

/-- The exchange case of Lam's Lemma 7.6 (2015): for adjacent columns `j = i + 1`, right
multiplication by the transvection `1 + a·Eᵢⱼ` sends `Δ_I` to
`Δ_I + a·Δ_{(I∖{j})∪{i}}` whenever `j ∈ I` and `i ∉ I`. Adjacency is what makes the exchanged
set sort with `i` in the slot `j` vacated, so the second minor appears with coefficient `+a`
and no sign; a non-adjacent exchange picks up a sign and is not stated here. The hypothesis
`(i : ℕ) + 1 = j` is on natural values so that no `Fin` wrap-around can occur, and the
cardinality of the exchanged set is taken as a hypothesis to match `Matrix.maximalMinor`'s
dependent signature. -/
theorem maximalMinor_mul_transvection_of_notMem_of_mem (A : Matrix (Fin d) (Fin n) R)
    {i j : Fin n} (hadj : (i : ℕ) + 1 = j) (a : R) {I : Finset (Fin n)} (hI : I.card = d)
    (hi : i ∉ I) (hj : j ∈ I) (hI' : (insert i (I.erase j)).card = d) :
    (A * transvection i j a).maximalMinor I hI =
      A.maximalMinor I hI + a * A.maximalMinor (insert i (I.erase j)) hI' := by
  unfold maximalMinor
  obtain ⟨p, hp⟩ : ∃ p, I.orderEmbOfFin hI p = j := by
    have hmem : j ∈ Set.range (I.orderEmbOfFin hI) := by
      rw [Finset.range_orderEmbOfFin]
      exact Finset.mem_coe.2 hj
    exact hmem
  have hmemI : ∀ x : Fin d, I.orderEmbOfFin hI x ∈ I := I.orderEmbOfFin_mem hI
  have hsub : (A * transvection i j a).submatrix id (I.orderEmbOfFin hI) =
      (A.submatrix id (I.orderEmbOfFin hI)).updateCol p
        (fun r ↦ A.submatrix id (I.orderEmbOfFin hI) r p + a * A r i) := by
    ext r p'
    by_cases hp' : p' = p
    · subst hp'
      simp only [updateCol_self, submatrix_apply, id_eq, hp]
      exact mul_transvection_apply_same i j r a A
    · simp only [updateCol_ne hp', submatrix_apply, id_eq]
      exact mul_transvection_apply_of_ne i j r (I.orderEmbOfFin hI p')
        (fun h ↦ hp' ((I.orderEmbOfFin hI).injective (h.trans hp.symm))) a A
  have hmem' : ∀ x : Fin d,
      Function.update (⇑(I.orderEmbOfFin hI)) p i x ∈ insert i (I.erase j) := by
    intro x
    by_cases hx : x = p
    · subst hx
      rw [Function.update_self]
      exact Finset.mem_insert_self i _
    · rw [Function.update_of_ne hx]
      exact Finset.mem_insert_of_mem (Finset.mem_erase.2
        ⟨fun h ↦ hx ((I.orderEmbOfFin hI).injective (h.trans hp.symm)), hmemI x⟩)
  have hmono : StrictMono (Function.update (⇑(I.orderEmbOfFin hI)) p i) := by
    have hlt : ∀ x : Fin d, x < p → I.orderEmbOfFin hI x < i := by
      intro x hx
      have h1 : ((I.orderEmbOfFin hI) x : ℕ) < (j : ℕ) :=
        Fin.lt_def.1 (hp ▸ (I.orderEmbOfFin hI).strictMono hx)
      have h2 : ((I.orderEmbOfFin hI) x : ℕ) ≠ (i : ℕ) :=
        fun h ↦ hi (Fin.val_injective h ▸ hmemI x)
      exact Fin.lt_def.2 (by omega)
    have hgt : ∀ x : Fin d, p < x → i < I.orderEmbOfFin hI x := by
      intro x hx
      have h1 : (j : ℕ) < ((I.orderEmbOfFin hI) x : ℕ) :=
        Fin.lt_def.1 (hp ▸ (I.orderEmbOfFin hI).strictMono hx)
      exact Fin.lt_def.2 (by omega)
    intro x y hxy
    rcases eq_or_ne x p with hx | hx <;> rcases eq_or_ne y p with hy | hy
    · rw [hx, hy] at hxy
      exact absurd hxy (lt_irrefl p)
    · subst hx
      rw [Function.update_self, Function.update_of_ne hy]
      exact hgt y hxy
    · subst hy
      rw [Function.update_self, Function.update_of_ne hx]
      exact hlt x hxy
    · rw [Function.update_of_ne hx, Function.update_of_ne hy]
      exact (I.orderEmbOfFin hI).strictMono hxy
  have hemb' : ⇑((insert i (I.erase j)).orderEmbOfFin hI') =
      Function.update (⇑(I.orderEmbOfFin hI)) p i :=
    (Finset.orderEmbOfFin_unique hI' hmem' hmono).symm
  have hsplitfun : (fun r ↦ A.submatrix id (I.orderEmbOfFin hI) r p + a * A r i) =
      (fun r ↦ A.submatrix id (I.orderEmbOfFin hI) r p) + a • (fun r ↦ A r i) := rfl
  have hcolself : (A.submatrix id (I.orderEmbOfFin hI)).updateCol p
      (fun r ↦ A.submatrix id (I.orderEmbOfFin hI) r p) =
        A.submatrix id (I.orderEmbOfFin hI) := by
    ext r p'
    by_cases hp' : p' = p
    · subst hp'
      rw [updateCol_self]
    · rw [updateCol_ne hp']
  have hcolexch : (A.submatrix id (I.orderEmbOfFin hI)).updateCol p (fun r ↦ A r i) =
      A.submatrix id ((insert i (I.erase j)).orderEmbOfFin hI') := by
    ext r p'
    rw [submatrix_apply, id_eq, hemb']
    by_cases hp' : p' = p
    · subst hp'
      rw [updateCol_self, Function.update_self]
    · rw [updateCol_ne hp', Function.update_of_ne hp', submatrix_apply, id_eq]
  rw [hsub, hsplitfun, det_updateCol_add, det_updateCol_smul, hcolself, hcolexch]

end Matrix

/-! ### Torus-fixed points are positroids -/

namespace Matroid

/-- For any `I ⊆ [n]`, the matroid on `[n]` whose unique base is `I` is a positroid: it is
realized by the 0/1 matrix with the identity in columns `I` and zeros elsewhere, whose only
nonvanishing maximal minor is at `I`. These matroids index the zero-dimensional cells of the
totally nonnegative Grassmannian — its torus-fixed points, `span(eᵢ : i ∈ I)` — which Lam's
lollipop networks parameterize (Lam, 2015, Example 4.2). -/
theorem isPositroid_uniqueBaseOn {n : ℕ} (I : Finset (Fin n)) :
    (uniqueBaseOn (↑I) (Set.univ : Set (Fin n))).IsPositroid := by
  classical
  set A : Matrix (Fin I.card) (Fin n) ℝ :=
    Matrix.of fun r c ↦ if c = I.orderEmbOfFin rfl r then 1 else 0 with hA
  have hself : ∀ hJ : I.card = I.card, A.maximalMinor I hJ = 1 := by
    intro hJ
    have hone : A.submatrix id (I.orderEmbOfFin hJ) = 1 := by
      ext r q
      rw [Matrix.submatrix_apply, id_eq, hA, Matrix.of_apply, Matrix.one_apply]
      by_cases h : r = q
      · subst h
        rw [if_pos rfl, if_pos rfl]
      · rw [if_neg h, if_neg fun hc ↦ h ((I.orderEmbOfFin rfl).injective hc).symm]
    unfold Matrix.maximalMinor
    rw [hone, Matrix.det_one]
  have hother : ∀ (J : Finset (Fin n)) (hJ : J.card = I.card), J ≠ I →
      A.maximalMinor J hJ = 0 := by
    intro J hJ hne
    obtain ⟨c₀, hc₀J, hc₀I⟩ : ∃ c₀ ∈ J, c₀ ∉ I := by
      by_contra h
      push Not at h
      exact hne (Finset.eq_of_subset_of_card_le h (le_of_eq hJ.symm))
    obtain ⟨q₀, hq₀⟩ : ∃ q₀, J.orderEmbOfFin hJ q₀ = c₀ := by
      have hmem : c₀ ∈ Set.range (J.orderEmbOfFin hJ) := by
        rw [Finset.range_orderEmbOfFin]
        exact Finset.mem_coe.2 hc₀J
      exact hmem
    unfold Matrix.maximalMinor
    refine Matrix.det_eq_zero_of_column_eq_zero q₀ fun r ↦ ?_
    rw [Matrix.submatrix_apply, id_eq, hq₀, hA, Matrix.of_apply,
      if_neg fun h ↦ hc₀I (by rw [h]; exact I.orderEmbOfFin_mem rfl r)]
  have hrank : A.rank = I.card := by
    have hunit : IsUnit (A.submatrix id (I.orderEmbOfFin rfl)) := by
      rw [Matrix.isUnit_iff_isUnit_det]
      rw [show (A.submatrix id (I.orderEmbOfFin rfl)).det = A.maximalMinor I rfl from rfl,
        hself rfl]
      exact isUnit_one
    have hle : A.rank ≤ I.card := by simpa using A.rank_le_card_height
    have hge : I.card ≤ A.rank := by
      have h := Matrix.rank_submatrix_le A id (I.orderEmbOfFin rfl)
      rwa [Matrix.rank_of_isUnit _ hunit, Fintype.card_fin] at h
    exact le_antisymm hle hge
  refine ⟨uniqueBaseOn_ground, I.card, A, hrank, fun J hJ ↦ ?_, fun B ↦ ?_⟩
  · rcases eq_or_ne J I with h | h
    · subst h
      rw [hself hJ]
      norm_num
    · rw [hother J hJ h]
  · rw [uniqueBaseOn_isBase_iff (Set.subset_univ _)]
    constructor
    · intro hB
      have hBI : B = I := Finset.coe_injective hB
      subst hBI
      exact ⟨rfl, by rw [hself rfl]; norm_num⟩
    · rintro ⟨hB, hpos⟩
      rcases eq_or_ne B I with h | h
      · rw [h]
      · rw [hother B hB h] at hpos
        exact absurd hpos (lt_irrefl 0)

end Matroid

namespace BoundaryMeasurement

/-! ### The `Gr(2,4)` square-graph fixture (Lam, Example 4.3) -/

variable {R : Type*} [CommRing R]

/-- Lam's square-graph boundary measurement matrix (2015, Example 4.3): the point of `Gr(2,4)`
measured from the square network with edge weights `a, b, c, d`, written with the identity in
the source columns `{2, 4}` (0-indexed `{1, 3}`). The entry `-a` is Postnikov's `(-1)^s` sign
(2006, Definition 4.6): `s = 1` because one source, `2`, lies strictly between source `4` and
sink `1` — so this single matrix exhibits Lam's sign-free convention and Postnikov's signed
convention agreeing. -/
def squareMatrix (a b c d : R) : Matrix (Fin 2) (Fin 4) R :=
  !![b, 1, c, 0; -a, 0, d, 1]

/-- `Δ₁₂ = a` for the square-graph fixture (Lam, 2015, Example 4.3; sets 0-indexed here, so
this is the minor on columns `{0, 1}`). -/
theorem maximalMinor_squareMatrix_01 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {0, 1} (by decide) = a := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 1)]
  simp [squareMatrix]

/-- `Δ₁₃ = ac + bd` for the square-graph fixture (Lam, 2015, Example 4.3): the two-path minor,
sum over the two vertex-disjoint path families. -/
theorem maximalMinor_squareMatrix_02 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {0, 2} (by decide) = a * c + b * d := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 2)]
  simp [squareMatrix]
  ring

/-- `Δ₁₄ = b` for the square-graph fixture (Lam, 2015, Example 4.3). -/
theorem maximalMinor_squareMatrix_03 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {0, 3} (by decide) = b := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 3)]
  simp [squareMatrix]

/-- `Δ₂₃ = d` for the square-graph fixture (Lam, 2015, Example 4.3). -/
theorem maximalMinor_squareMatrix_12 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {1, 2} (by decide) = d := by
  rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 2)]
  simp [squareMatrix]

/-- `Δ₂₄ = 1` for the square-graph fixture (Lam, 2015, Example 4.3): the source columns carry
the identity, so this minor is `1` for every choice of weights. -/
theorem maximalMinor_squareMatrix_13 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {1, 3} (by decide) = 1 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 3)]
  simp [squareMatrix]

/-- `Δ₃₄ = c` for the square-graph fixture (Lam, 2015, Example 4.3). -/
theorem maximalMinor_squareMatrix_23 (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {2, 3} (by decide) = c := by
  rw [Matrix.maximalMinor_pair _ (by decide : (2 : Fin 4) < 3)]
  simp [squareMatrix]

/-- The square-graph fixture satisfies the Plücker relation
`Δ₁₂Δ₃₄ - Δ₁₃Δ₂₄ + Δ₁₄Δ₂₃ = 0` — the page's `ac - (ac + bd) + bd = 0` certificate that the
boundary measurements really are the Plücker coordinates of a point of the Grassmannian
(Lam, 2015, Example 4.3). An instance of `Matrix.plucker_fin_two`. -/
theorem plucker_squareMatrix (a b c d : R) :
    (squareMatrix a b c d).maximalMinor {0, 1} (by decide) *
        (squareMatrix a b c d).maximalMinor {2, 3} (by decide) -
      (squareMatrix a b c d).maximalMinor {0, 2} (by decide) *
        (squareMatrix a b c d).maximalMinor {1, 3} (by decide) +
      (squareMatrix a b c d).maximalMinor {0, 3} (by decide) *
        (squareMatrix a b c d).maximalMinor {1, 2} (by decide) = 0 := by
  exact Matrix.plucker_fin_two _ (by decide) (by decide) (by decide)

section Ordered

variable {R : Type*} [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]

/-- At strictly positive weights, every maximal minor of the square-graph fixture is strictly
positive: the network parameterizes the top cell of the totally nonnegative Grassmannian —
all `binom(4)(2) = 6` bases, the uniform matroid `U_{2,4}` (Lam, 2015, Example 4.3; verified
against the repository's Python implementation on the research page). -/
theorem maximalMinor_squareMatrix_pos {a b c d : R} (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hd : 0 < d) {I : Finset (Fin 4)} (hI : I.card = 2) :
    0 < (squareMatrix a b c d).maximalMinor I hI := by
  have hcases : ∀ J : Finset (Fin 4), J.card = 2 →
      J = {0, 1} ∨ J = {0, 2} ∨ J = {0, 3} ∨ J = {1, 2} ∨ J = {1, 3} ∨ J = {2, 3} := by
    decide
  rcases hcases I hI with h | h | h | h | h | h <;> subst h
  · rw [maximalMinor_squareMatrix_01]
    exact ha
  · rw [maximalMinor_squareMatrix_02]
    exact add_pos (mul_pos ha hc) (mul_pos hb hd)
  · rw [maximalMinor_squareMatrix_03]
    exact hb
  · rw [maximalMinor_squareMatrix_12]
    exact hd
  · rw [maximalMinor_squareMatrix_13]
    exact one_pos
  · rw [maximalMinor_squareMatrix_23]
    exact hc

end Ordered

/-- The square-graph fixture has full rank `2` for every choice of weights: the source-column
minor `Δ₂₄ = 1` never vanishes (Lam, 2015, Example 4.3). -/
theorem rank_squareMatrix (a b c d : ℝ) : (squareMatrix a b c d).rank = 2 := by
  have hI : ({1, 3} : Finset (Fin 4)).card = 2 := by decide
  have hunit : IsUnit ((squareMatrix a b c d).submatrix id
      (({1, 3} : Finset (Fin 4)).orderEmbOfFin hI)) := by
    rw [Matrix.isUnit_iff_isUnit_det]
    have hdet : ((squareMatrix a b c d).submatrix id
        (({1, 3} : Finset (Fin 4)).orderEmbOfFin hI)).det = 1 :=
      maximalMinor_squareMatrix_13 a b c d
    rw [hdet]
    exact isUnit_one
  have hle : (squareMatrix a b c d).rank ≤ 2 := by
    simpa using (squareMatrix a b c d).rank_le_card_height
  have hge : 2 ≤ (squareMatrix a b c d).rank := by
    have h := Matrix.rank_submatrix_le (squareMatrix a b c d) id
      (({1, 3} : Finset (Fin 4)).orderEmbOfFin hI)
    rwa [Matrix.rank_of_isUnit _ hunit, Fintype.card_fin] at h
  exact le_antisymm hle hge

/-- At strictly positive weights the square-graph fixture realizes the uniform matroid
`U_{2,4}`: the bases of `Matroid.unifOn` on `[4]` of rank `2` are exactly the column pairs
with positive minor. This is the fixture's "it lands in the top cell" certificate (Lam, 2015,
Example 4.3), stated against the repository's positroid vocabulary. -/
theorem unifOn_isBase_iff_maximalMinor_squareMatrix_pos {a b c d : ℝ} (ha : 0 < a)
    (hb : 0 < b) (hc : 0 < c) (hd : 0 < d) {B : Finset (Fin 4)} :
    (Matroid.unifOn (Set.univ : Set (Fin 4)) 2).IsBase ↑B ↔
      ∃ hB : B.card = 2, 0 < (squareMatrix a b c d).maximalMinor B hB := by
  rw [Matroid.unifOn_isBase_iff (by norm_num)]
  exact ⟨fun hB ↦ ⟨hB, maximalMinor_squareMatrix_pos ha hb hc hd hB⟩, fun ⟨hB, _⟩ ↦ hB⟩

/-- The trip permutation of the square graph (Lam, 2015, §7.1): `π_G = (3,4,1,2)` 1-indexed,
that is `i ↦ i + 2` on `Fin 4`, with no fixed points and hence an empty decoration. Under
Postnikov's bijection this decorated permutation names the top cell of `Gr(2,4)^{tnn}` — the
research page verifies it equals the decorated permutation the repository's Python
implementation computes for the same fixture. -/
def squareTrip : DecoratedPermutation 4 where
  toPerm := finRotate 4 ^ 2
  clockwise := ∅
  apply_eq_self_of_mem_clockwise := by simp

/-- The square-graph trip permutation sends `i` to `i + 2` (0-indexed form of Lam's
`(3,4,1,2)`, 2015, §7.1). -/
theorem squareTrip_toPerm_apply (i : Fin 4) : squareTrip.toPerm i = i + 2 := by
  revert i
  decide

/-- The square-graph trip permutation has exactly two weak excedances, matching the rank `2`
of the cell it names (Postnikov's bijection; Lam, 2015, §7.1). -/
theorem card_weakExcedances_squareTrip : squareTrip.weakExcedances.card = 2 := by
  decide

/-! ### The lollipop fixture (Lam, Example 4.2) -/

/-- Lam's lollipop network point (2015, Example 4.2): the matrix realizing
`span(e₃, e₄) ∈ Gr(2,4)` (0-indexed: columns `2` and `3`), a torus-fixed point. The lollipop
graph has a single almost perfect matching, whose boundary subset is `{3, 4}`. -/
def lollipopMatrix : Matrix (Fin 2) (Fin 4) ℝ :=
  !![0, 0, 1, 0; 0, 0, 0, 1]

/-- The lollipop point's minor on its base `{3, 4}` (0-indexed `{2, 3}`) is `1` (Lam, 2015,
Example 4.2). -/
theorem maximalMinor_lollipopMatrix_self :
    lollipopMatrix.maximalMinor {2, 3} (by decide) = 1 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (2 : Fin 4) < 3)]
  simp [lollipopMatrix]

/-- Every other maximal minor of the lollipop point vanishes: the point is the torus-fixed
point `span(e₃, e₄)`, and its cell is zero-dimensional (Lam, 2015, Example 4.2). -/
theorem maximalMinor_lollipopMatrix_of_ne {I : Finset (Fin 4)} (hI : I.card = 2)
    (hne : I ≠ {2, 3}) : lollipopMatrix.maximalMinor I hI = 0 := by
  have hcases : ∀ J : Finset (Fin 4), J.card = 2 →
      J = {0, 1} ∨ J = {0, 2} ∨ J = {0, 3} ∨ J = {1, 2} ∨ J = {1, 3} ∨ J = {2, 3} := by
    decide
  rcases hcases I hI with h | h | h | h | h | h <;> subst h
  · rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 1)]
    simp [lollipopMatrix]
  · rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 2)]
    simp [lollipopMatrix]
  · rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 3)]
    simp [lollipopMatrix]
  · rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 2)]
    simp [lollipopMatrix]
  · rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 3)]
    simp [lollipopMatrix]
  · exact absurd rfl hne

/-- The lollipop's matroid — unique base `{3, 4}` (0-indexed `{2, 3}`) — is a positroid: the
zero-dimensional-cell certificate of Lam's Example 4.2 (2015), as an instance of
`Matroid.isPositroid_uniqueBaseOn`. -/
theorem isPositroid_lollipop :
    (Matroid.uniqueBaseOn (↑({2, 3} : Finset (Fin 4))) (Set.univ : Set (Fin 4))).IsPositroid :=
  Matroid.isPositroid_uniqueBaseOn _

/-! ### Postnikov's four-vertex matrix (Example 4.7) -/

/-- Postnikov's Example 4.7 (2006): the boundary measurement matrix of a network on four
boundary vertices with source set `I = {1, 3}` and sink set `{2, 4}`. Exactly one entry is
negated — `-M₁₄`, at row `1`, column `4` — because `s = 1`: the source `3` lies strictly
between `1` and `4` (Definition 4.6's `(-1)^s` rule). -/
def fourVertexMatrix (m12 m14 m32 m34 : R) : Matrix (Fin 2) (Fin 4) R :=
  !![1, m12, 0, -m14; 0, m32, 1, m34]

/-- `Δ₁₃ = 1` for Postnikov's Example 4.7 (2006): the source columns carry the identity. -/
theorem maximalMinor_fourVertexMatrix_02 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {0, 2} (by decide) = 1 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 2)]
  simp [fourVertexMatrix]

/-- `Δ₂₃ = M₁₂` for Postnikov's Example 4.7 (2006): the exchange minor `(I∖{1})∪{2}` recovers
the unsigned boundary measurement, as Definition 4.6's sign rule guarantees. -/
theorem maximalMinor_fourVertexMatrix_12 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {1, 2} (by decide) = m12 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 2)]
  simp [fourVertexMatrix]

/-- `Δ₁₂ = M₃₂` for Postnikov's Example 4.7 (2006): the exchange minor `(I∖{3})∪{2}`. -/
theorem maximalMinor_fourVertexMatrix_01 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {0, 1} (by decide) = m32 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 1)]
  simp [fourVertexMatrix]

/-- `Δ₃₄ = M₁₄` for Postnikov's Example 4.7 (2006): the sign-critical exchange minor
`(I∖{1})∪{4}` — the negated entry `-M₁₄` is exactly what makes this minor come out to the
*unsigned* `M₁₄`. Note: Postnikov's printed example carries an uncorrected typo here
(`M₁₄ = Δ₂₄/Δ₁₃` in the source); his own Definition 4.6 and §5 force `Δ₃₄`, which is what
this theorem proves. -/
theorem maximalMinor_fourVertexMatrix_23 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {2, 3} (by decide) = m14 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (2 : Fin 4) < 3)]
  simp [fourVertexMatrix]

/-- `Δ₁₄ = M₃₄` for Postnikov's Example 4.7 (2006): the exchange minor `(I∖{3})∪{4}`. -/
theorem maximalMinor_fourVertexMatrix_03 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {0, 3} (by decide) = m34 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 3)]
  simp [fourVertexMatrix]

/-- The companion minor of Postnikov's Example 4.7: `Δ₂₄ = M₁₂M₃₄ + M₁₄M₃₂`, "because both
bijections `π : {1,3} → {2,4}` have just one misalignment and no crossings" (Postnikov, 2006,
§5) — the smallest instance of his Proposition 5.2 expressing minors as immanants in the
boundary measurements. -/
theorem maximalMinor_fourVertexMatrix_13 (m12 m14 m32 m34 : R) :
    (fourVertexMatrix m12 m14 m32 m34).maximalMinor {1, 3} (by decide) =
      m12 * m34 + m14 * m32 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 3)]
  simp [fourVertexMatrix]

/-! ### Talaska's `n = 5` sign fixture (Example 2.7) -/

/-- Talaska's Example 2.7 (2008): the boundary measurement matrix for her Figure 1 network,
with source set `I = {1, 4}` in `[5]`. Again exactly one entry is negated — `-M₁₅`, at row
`1`, column `5` — since the source `4` lies strictly between `1` and `5`. -/
def talaskaMatrix (m12 m13 m15 m42 m43 m45 : R) : Matrix (Fin 2) (Fin 5) R :=
  !![1, m12, m13, 0, -m15; 0, m42, m43, 1, m45]

/-- `Δ₁₄ = 1` for Talaska's Example 2.7 (2008): the source columns carry the identity. -/
theorem maximalMinor_talaskaMatrix_03 (m12 m13 m15 m42 m43 m45 : R) :
    (talaskaMatrix m12 m13 m15 m42 m43 m45).maximalMinor {0, 3} (by decide) = 1 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 5) < 3)]
  simp [talaskaMatrix]

/-- `Δ₄₅ = M₁₅` for Talaska's Example 2.7 (2008): the sign-critical exchange minor at `n = 5` —
the negated entry `-M₁₅` makes the minor on `(I∖{1})∪{5}` come out to the unsigned `M₁₅`,
certifying the `(-1)^s` rule beyond `n = 4`. -/
theorem maximalMinor_talaskaMatrix_34 (m12 m13 m15 m42 m43 m45 : R) :
    (talaskaMatrix m12 m13 m15 m42 m43 m45).maximalMinor {3, 4} (by decide) = m15 := by
  rw [Matrix.maximalMinor_pair _ (by decide : (3 : Fin 5) < 4)]
  simp [talaskaMatrix]

/-! ### Postnikov's five-basis positroid (Example 11.9) -/

/-- The basis family of Postnikov's Example 11.9 (2006): the five source sets of the perfect
orientations of his plabic graph — every `2`-subset of `[4]` except `{2, 3}` (0-indexed:
except `{1, 2}`). -/
def fiveBasisBases : Finset (Finset (Fin 4)) :=
  {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {2, 3}}

/-- The matroid of Postnikov's Example 11.9 (2006): the rank-2 matroid on `[4]` whose bases
are the five `2`-subsets other than `{2, 3}` (0-indexed `{1, 2}`) — the matroid `𝓜_G` of his
five-perfect-orientation plabic graph, and the smallest non-uniform positroid fixture. -/
def fiveBasisMatroid : Matroid (Fin 4) :=
  (IndepMatroid.ofFinset (Set.univ : Set (Fin 4)) (fun J ↦ ∃ B ∈ fiveBasisBases, J ⊆ B)
    ⟨{0, 1}, by decide, Finset.empty_subset _⟩
    (fun I J hJ hIJ ↦ by
      obtain ⟨B, hB, hJB⟩ := hJ
      exact ⟨B, hB, hIJ.trans hJB⟩)
    (by decide)
    (fun I _ ↦ Set.subset_univ _)).matroid

/-- The bases of `fiveBasisMatroid` are exactly the five listed `2`-subsets (Postnikov, 2006,
Example 11.9). -/
theorem fiveBasisMatroid_isBase_iff {B : Finset (Fin 4)} :
    fiveBasisMatroid.IsBase ↑B ↔ B ∈ fiveBasisBases := by
  have hindep : ∀ J : Finset (Fin 4),
      fiveBasisMatroid.Indep ↑J ↔ ∃ B' ∈ fiveBasisBases, J ⊆ B' := by
    intro J
    simp [fiveBasisMatroid]
  have hcard2 : ∀ B' ∈ fiveBasisBases, Finset.card B' = 2 := by decide
  rw [Matroid.isBase_iff_maximal_indep]
  constructor
  · rintro ⟨hind, hmax⟩
    obtain ⟨B', hB', hBB'⟩ := (hindep B).1 hind
    have hsub : ↑B' ⊆ (↑B : Set (Fin 4)) :=
      hmax ((hindep B').2 ⟨B', hB', Finset.Subset.refl B'⟩) (Finset.coe_subset.2 hBB')
    have hBeq : B = B' := Finset.Subset.antisymm hBB' (Finset.coe_subset.1 hsub)
    rw [hBeq]
    exact hB'
  · intro hB
    refine ⟨(hindep B).2 ⟨B, hB, Finset.Subset.refl B⟩, ?_⟩
    intro Y hY hBY
    simp only [fiveBasisMatroid] at hY
    rw [IndepMatroid.matroid_indep_iff, IndepMatroid.ofFinset_indep'] at hY
    intro x hxY
    by_contra hxB
    have hsubY : ↑(insert x B) ⊆ Y := by
      rw [Finset.coe_insert]
      exact Set.insert_subset_iff.2 ⟨hxY, hBY⟩
    obtain ⟨B'', hB'', hins⟩ := hY _ hsubY
    have h3 : (insert x B).card = 3 := by
      rw [Finset.card_insert_of_notMem (by simpa using hxB), hcard2 B hB]
    have hle := Finset.card_le_card hins
    rw [h3, hcard2 B'' hB''] at hle
    omega

/-- A realizing matrix for Postnikov's Example 11.9: columns `(1,0), (1,1), (1,1), (0,1)`.
Its maximal minor vanishes exactly on the column pair `{2, 3}` (0-indexed `{1, 2}`), where the
two parallel columns meet, and equals `1` on the other five pairs. -/
def fiveBasisMatrix : Matrix (Fin 2) (Fin 4) ℝ :=
  !![1, 1, 1, 0; 0, 1, 1, 1]

/-- Postnikov's Example 11.9 matroid is a positroid (2006, Proposition 11.7 predicts this for
the matroid of any perfectly orientable plabic graph; here it is certified directly by the
realizing matrix `fiveBasisMatrix`). The smallest non-uniform companion to
`Matroid.isPositroid_unifOn`. -/
theorem isPositroid_fiveBasisMatroid : fiveBasisMatroid.IsPositroid := by
  have h01 : fiveBasisMatrix.maximalMinor {0, 1} (by decide) = 1 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 1)]
    simp [fiveBasisMatrix]
  have h02 : fiveBasisMatrix.maximalMinor {0, 2} (by decide) = 1 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 2)]
    simp [fiveBasisMatrix]
  have h03 : fiveBasisMatrix.maximalMinor {0, 3} (by decide) = 1 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (0 : Fin 4) < 3)]
    simp [fiveBasisMatrix]
  have h12 : fiveBasisMatrix.maximalMinor {1, 2} (by decide) = 0 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 2)]
    simp [fiveBasisMatrix]
  have h13 : fiveBasisMatrix.maximalMinor {1, 3} (by decide) = 1 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (1 : Fin 4) < 3)]
    simp [fiveBasisMatrix]
  have h23 : fiveBasisMatrix.maximalMinor {2, 3} (by decide) = 1 := by
    rw [Matrix.maximalMinor_pair _ (by decide : (2 : Fin 4) < 3)]
    simp [fiveBasisMatrix]
  have hcases : ∀ J : Finset (Fin 4), J.card = 2 →
      J = {0, 1} ∨ J = {0, 2} ∨ J = {0, 3} ∨ J = {1, 2} ∨ J = {1, 3} ∨ J = {2, 3} := by
    decide
  have hrank : fiveBasisMatrix.rank = 2 := by
    have hI : ({0, 1} : Finset (Fin 4)).card = 2 := by decide
    have hunit : IsUnit (fiveBasisMatrix.submatrix id
        (({0, 1} : Finset (Fin 4)).orderEmbOfFin hI)) := by
      rw [Matrix.isUnit_iff_isUnit_det]
      have hdet : (fiveBasisMatrix.submatrix id
          (({0, 1} : Finset (Fin 4)).orderEmbOfFin hI)).det = 1 := h01
      rw [hdet]
      exact isUnit_one
    have hle : fiveBasisMatrix.rank ≤ 2 := by simpa using fiveBasisMatrix.rank_le_card_height
    have hge : 2 ≤ fiveBasisMatrix.rank := by
      have h := Matrix.rank_submatrix_le fiveBasisMatrix id
        (({0, 1} : Finset (Fin 4)).orderEmbOfFin hI)
      rwa [Matrix.rank_of_isUnit _ hunit, Fintype.card_fin] at h
    exact le_antisymm hle hge
  refine ⟨by simp [fiveBasisMatroid], 2, fiveBasisMatrix, hrank, fun J hJ ↦ ?_, fun B ↦ ?_⟩
  · rcases hcases J hJ with h | h | h | h | h | h <;> subst h
    · rw [h01]
      norm_num
    · rw [h02]
      norm_num
    · rw [h03]
      norm_num
    · rw [h12]
    · rw [h13]
      norm_num
    · rw [h23]
      norm_num
  · rw [fiveBasisMatroid_isBase_iff]
    constructor
    · intro hB
      have hmem : B = {0, 1} ∨ B = {0, 2} ∨ B = {0, 3} ∨ B = {1, 3} ∨ B = {2, 3} :=
        (by decide : ∀ B' : Finset (Fin 4), B' ∈ fiveBasisBases →
          B' = {0, 1} ∨ B' = {0, 2} ∨ B' = {0, 3} ∨ B' = {1, 3} ∨ B' = {2, 3}) B hB
      rcases hmem with h | h | h | h | h <;> subst h
      · exact ⟨by decide, by rw [h01]; norm_num⟩
      · exact ⟨by decide, by rw [h02]; norm_num⟩
      · exact ⟨by decide, by rw [h03]; norm_num⟩
      · exact ⟨by decide, by rw [h13]; norm_num⟩
      · exact ⟨by decide, by rw [h23]; norm_num⟩
    · rintro ⟨hB, hpos⟩
      rcases hcases B hB with h | h | h | h | h | h <;> subst h
      · decide
      · decide
      · decide
      · rw [h12] at hpos
        exact absurd hpos (lt_irrefl 0)
      · decide
      · decide

end BoundaryMeasurement
