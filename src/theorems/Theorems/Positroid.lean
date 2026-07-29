/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Combinatorics.Matroid.Basic
import Mathlib.Combinatorics.Matroid.IndepAxioms
import Mathlib.Combinatorics.Young.YoungDiagram
import Mathlib.Data.Finset.Powerset
import Mathlib.Data.Finset.Sort
import Mathlib.Data.Real.Basic
import Mathlib.GroupTheory.Perm.Basic
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.Tactic

/-!
# Positroids

A *positroid* of rank `d` on the ordered ground set `[n]` is a matroid that can be represented
by the columns of a full-rank `d × n` real matrix all of whose maximal minors are nonnegative;
its bases are the `d`-subsets whose minor is strictly positive (Ardila–Rincón–Williams,
*Positroids and non-crossing partitions*, 2016, following Postnikov, *Total positivity,
Grassmannians, and networks*, 2006). This module formalizes the positroid predicate on
`Matroid (Fin n)` together with the combinatorial indexing objects the theory is built from:
the Gale order and its cyclic shifts, Grassmann necklaces, decorated permutations, and
Le-diagrams.

Being a positroid genuinely depends on the linear order of the ground set (it is invariant
only under cyclic shifts), so the definition is stated for matroids on `Fin n` and inherits
that type's order.

## Main definitions

* `Matrix.maximalMinor` — the `d × d` minor of a `d × n` matrix on a column set of size `d`,
  taken with the columns in increasing order.
* `Matroid.unifOn` — the uniform matroid of rank `k` on a ground set `E`, whose independent
  sets are the subsets of `E` of size at most `k`.
* `Matroid.IsPositroid` — the Ardila–Rincón–Williams positroid predicate on `Matroid (Fin n)`.
* `Finset.GaleLE` — the Gale order: `S` is below `T` when they have equal size and the `j`-th
  smallest element of `S` is at most the `j`-th smallest element of `T` for every `j`.
* `Finset.GaleShiftLE` — the Gale order induced by the cyclically shifted order `<ᵢ` on `Fin n`.
* `GrassmannNecklace` — Postnikov's Grassmann necklaces of type `(d, n)`.
* `GrassmannNecklace.ohFamily` — the family `𝓑(𝓘) = {B : Iⱼ ≤ⱼ B for all j}` from Oh's
  theorem, the candidate basis family of the positroid indexed by a necklace.
* `DecoratedPermutation` — permutations of `Fin n` with fixed points colored clockwise or
  counterclockwise, and their weak excedances.
* `LeDiagram` — 0/1 fillings of a Young diagram satisfying the Γ-property.

## Main statements

* `Matrix.det_vandermonde_pos` — a Vandermonde determinant with strictly increasing nodes is
  positive.
* `Matroid.unifOn_isBase_iff` — the bases of the rank-`k` uniform matroid on `Fin n` (with
  `k ≤ n`) are exactly the `k`-subsets.
* `Matroid.isPositroid_unifOn` — the uniform matroid `U_{d,n}` with the standard order is a
  positroid, realized by a totally positive Vandermonde-power matrix. This is the "top cell"
  example of the theory.
* `Finset.GaleLE.trans` / `Finset.GaleLE.antisymm` — the Gale order is a partial order.

## Implementation notes

Mathlib's matroid library supplies the matroid itself (`Matroid`, its bases, duals, minors,
loops, and sums), `Finset.orderEmbOfFin` supplies sorted access to a finite set, and
`Matrix.vandermonde` supplies the determinant formula used for the top-cell example.
Everything positroid-specific is new here: Mathlib currently has no matroid representability,
no uniform matroids, no total positivity, no Gale order, and none of the cryptomorphic
indexing objects. `Matroid.unifOn` fills the uniform-matroid gap following the naming used by
more recent Mathlib revisions so it can be deleted on upgrade; unlike upstream it is built
from `IndepMatroid.ofFinset` and therefore assumes `DecidableEq`.

`Finset.GaleLE S T` is stated as an existential over the proof that the cardinalities agree,
so that it is a total (and decidable) relation: sets of different sizes are incomparable.
For the shifted order `<ᵢ` on `Fin n` we use cyclic subtraction — `a ≤ᵢ b ↔ a - i ≤ b - i` —
so `Finset.GaleShiftLE` compares images under `· - i` rather than introducing a rotated
`LinearOrder` instance.

For `DecoratedPermutation.weakExcedances` the page-level convention choice (which fixed-point
color counts toward the rank) is resolved the same way as this repository's Python
implementation: clockwise fixed points count.

Although `Finset.GaleLE` is decidable, `decide` cannot evaluate it on concrete sets:
`Finset.sort` is defined by well-founded recursion and does not reduce in the kernel.
Concrete instances are proved through the `Finset.sort_insert` and `Finset.sort_singleton`
simp lemmas instead, as in the example following `Finset.GaleShiftLE`.

The theorems of the source literature that need machinery Mathlib lacks — Oh's theorem that
`ohFamily` is the basis family of a positroid, the necklace/decorated-permutation bijections,
closure of positroids under duality, minors, cyclic shifts and direct sums, and the
non-crossing partition decomposition — are recorded as backlog in the research notes rather
than stated here.

## References

* [A. Postnikov, *Total positivity, Grassmannians, and networks*][postnikov2006],
  arXiv:math/0609764. Grassmann necklaces (§16), decorated permutations, Le-diagrams (§6).
* [S. Oh, *Positroids and Schubert matroids*][oh2011], arXiv:0803.1018. The membership test
  `𝓑(𝓘)` and the Gale-order characterization of positroids.
* [F. Ardila, F. Rincón, L. Williams, *Positroids and non-crossing partitions*][arw2016],
  arXiv:1308.2698. The matrix definition of positroids used here (their §2–3).

## Tags

positroid, matroid, total positivity, Grassmann necklace, decorated permutation, Gale order
-/

open Finset

namespace Matrix

/-- The maximal minor of a `d × n` matrix `A` at a set `I` of `d` columns: the determinant of
the square submatrix built from the columns in `I`, taken in increasing order. For a positroid
this is Postnikov's Plücker coordinate `Δ_I(A)`. -/
def maximalMinor {R : Type*} [CommRing R] {d n : ℕ} (A : Matrix (Fin d) (Fin n) R)
    (I : Finset (Fin n)) (hI : I.card = d) : R :=
  (A.submatrix id (I.orderEmbOfFin hI)).det

/-- A Vandermonde determinant with strictly increasing nodes is positive. This is the reason
Vandermonde matrices realize the top cell of the totally nonnegative Grassmannian (Postnikov,
2006). -/
theorem det_vandermonde_pos {R : Type*} [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]
    {m : ℕ} {v : Fin m → R} (hv : StrictMono v) : 0 < (vandermonde v).det := by
  rw [det_vandermonde]
  refine Finset.prod_pos fun i _ ↦ Finset.prod_pos fun j hj ↦ ?_
  exact sub_pos.2 (hv (Finset.mem_Ioi.1 hj))

end Matrix

namespace Finset

variable {α : Type*} [LinearOrder α]

/-- The Gale order on finite sets in a linear order (Gale, 1968; used by Postnikov and Oh to
index positroids): `S.GaleLE T` holds when `S` and `T` have the same size and, listing both in
increasing order, the `j`-th element of `S` is at most the `j`-th element of `T` for every
`j`. Sets of different sizes are incomparable. -/
def GaleLE (S T : Finset α) : Prop :=
  ∃ h : T.card = S.card, ∀ j : Fin S.card, S.orderEmbOfFin rfl j ≤ T.orderEmbOfFin h j

instance (S T : Finset α) : Decidable (S.GaleLE T) :=
  inferInstanceAs (Decidable (∃ _ : T.card = S.card, ∀ _ : Fin S.card, _ ≤ _))

/-- The Gale order relates only sets of equal size. -/
theorem GaleLE.card_eq {S T : Finset α} (h : S.GaleLE T) : T.card = S.card := h.1

/-- `Finset.orderEmbOfFin` does not depend on which cardinality proof indexes it: the two
enumerations agree up to an index cast. -/
theorem orderEmbOfFin_cast {S : Finset α} {k k' : ℕ} (h : S.card = k) (h' : S.card = k')
    (j : Fin k) : S.orderEmbOfFin h j = S.orderEmbOfFin h' (Fin.cast (h.symm.trans h') j) := by
  subst h h'
  rfl

/-- The Gale order is reflexive. -/
theorem GaleLE.refl (S : Finset α) : S.GaleLE S := ⟨rfl, fun _ ↦ le_rfl⟩

/-- The Gale order is transitive. -/
theorem GaleLE.trans {S T U : Finset α} (hST : S.GaleLE T) (hTU : T.GaleLE U) : S.GaleLE U := by
  obtain ⟨hTS, pST⟩ := hST
  obtain ⟨hUT, pTU⟩ := hTU
  refine ⟨hUT.trans hTS, fun j ↦ (pST j).trans ?_⟩
  calc T.orderEmbOfFin hTS j
      = T.orderEmbOfFin rfl (Fin.cast hTS.symm j) := orderEmbOfFin_cast hTS rfl j
    _ ≤ U.orderEmbOfFin hUT (Fin.cast hTS.symm j) := pTU _
    _ = U.orderEmbOfFin (hUT.trans hTS) j := orderEmbOfFin_cast hUT (hUT.trans hTS) _

/-- The Gale order is antisymmetric. -/
theorem GaleLE.antisymm {S T : Finset α} (hST : S.GaleLE T) (hTS : T.GaleLE S) : S = T := by
  obtain ⟨hTS', pST⟩ := hST
  obtain ⟨hST', pTS⟩ := hTS
  have key : ∀ j : Fin S.card, S.orderEmbOfFin rfl j = T.orderEmbOfFin hTS' j := by
    intro j
    refine le_antisymm (pST j) ?_
    calc T.orderEmbOfFin hTS' j
        = T.orderEmbOfFin rfl (Fin.cast hTS'.symm j) := orderEmbOfFin_cast hTS' rfl j
      _ ≤ S.orderEmbOfFin hST' (Fin.cast hTS'.symm j) := pTS _
      _ = S.orderEmbOfFin rfl j := orderEmbOfFin_cast hST' rfl _
  apply Finset.coe_injective
  calc (S : Set α)
      = Set.range (S.orderEmbOfFin rfl) := (range_orderEmbOfFin S rfl).symm
    _ = Set.range (T.orderEmbOfFin hTS') := congrArg Set.range (funext key)
    _ = (T : Set α) := range_orderEmbOfFin T hTS'

/-- The Gale order on subsets of `Fin n` induced by the cyclically shifted order
`i <ᵢ i+1 <ᵢ ⋯ <ᵢ i-1` (Postnikov, 2006, §16). Comparing two elements in the shifted order is
comparing their cyclic differences from `i`, so the sorted lists relative to `<ᵢ` are the
sorted lists of the images under `· - i`. -/
def GaleShiftLE {n : ℕ} (i : Fin n) (S T : Finset (Fin n)) : Prop :=
  GaleLE (S.image fun a ↦ a - i) (T.image fun a ↦ a - i)

instance {n : ℕ} (i : Fin n) (S T : Finset (Fin n)) : Decidable (GaleShiftLE i S T) :=
  inferInstanceAs (Decidable (GaleLE _ _))

/-- `{1, 2} ≤ {1, 3}` in the Gale order: sorted componentwise, `1 ≤ 1` and `2 ≤ 3`. -/
example : ({1, 2} : Finset ℕ).GaleLE {1, 3} := by
  refine ⟨by decide, fun j ↦ ?_⟩
  fin_cases j <;>
    simp [Finset.orderEmbOfFin_apply, Finset.sort_insert, Finset.sort_singleton]

end Finset

namespace Matroid

/-- The uniform matroid of rank `k` on the ground set `E`: every subset of `E` with at most
`k` elements is independent. Follows the naming of more recent Mathlib revisions (which this
declaration should be replaced by on upgrade), but is constructed from
`IndepMatroid.ofFinset`, whence the `DecidableEq` assumption. -/
def unifOn {α : Type*} [DecidableEq α] (E : Set α) (k : ℕ) : Matroid α :=
  (IndepMatroid.ofFinset E (fun I ↦ ↑I ⊆ E ∧ I.card ≤ k)
    (by simp)
    (fun I J hJ hIJ ↦ ⟨(Finset.coe_subset.2 hIJ).trans hJ.1,
      (Finset.card_le_card hIJ).trans hJ.2⟩)
    (fun I J hI hJ hIJ ↦ by
      obtain ⟨e, heJ, heI⟩ := Finset.exists_mem_notMem_of_card_lt_card hIJ
      refine ⟨e, heJ, heI, ?_, ?_⟩
      · rw [Finset.coe_insert]
        exact Set.insert_subset_iff.2 ⟨hJ.1 (Finset.mem_coe.2 heJ), hI.1⟩
      · exact (Finset.card_insert_le e I).trans (Nat.succ_le_of_lt (hIJ.trans_le hJ.2)))
    (fun I hI ↦ hI.1)).matroid

/-- A finite set is independent in `Matroid.unifOn E k` iff it lies in `E` and has at most `k`
elements. -/
@[simp]
theorem unifOn_indep_finset_iff {α : Type*} [DecidableEq α] {E : Set α} {k : ℕ}
    {I : Finset α} : (unifOn E k).Indep ↑I ↔ ↑I ⊆ E ∧ I.card ≤ k := by
  simp [unifOn]

/-- The bases of the rank-`k` uniform matroid on all of `Fin n`, for `k ≤ n`, are exactly the
`k`-element subsets. -/
theorem unifOn_isBase_iff {n k : ℕ} (hkn : k ≤ n) {B : Finset (Fin n)} :
    (unifOn (Set.univ : Set (Fin n)) k).IsBase ↑B ↔ B.card = k := by
  rw [Matroid.isBase_iff_maximal_indep]
  constructor
  · rintro ⟨hind, hmax⟩
    rw [unifOn_indep_finset_iff] at hind
    by_contra hne
    obtain ⟨C, hBC, hC⟩ :=
      Finset.exists_superset_card_eq hind.2 (by simpa using hkn)
    have hCB : C ⊆ B := by
      exact_mod_cast hmax (unifOn_indep_finset_iff.2 ⟨by simp, hC.le⟩)
        (Finset.coe_subset.2 hBC)
    have hCeq : C = B := Finset.Subset.antisymm hCB hBC
    exact hne (hCeq ▸ hC)
  · intro hB
    refine ⟨unifOn_indep_finset_iff.2 ⟨by simp, hB.le⟩, ?_⟩
    intro Y hY hBY
    simp only [unifOn] at hY
    rw [IndepMatroid.matroid_indep_iff, IndepMatroid.ofFinset_indep'] at hY
    intro x hxY
    by_contra hxB
    have hsub : ↑(insert x B) ⊆ Y := by
      rw [Finset.coe_insert]
      exact Set.insert_subset_iff.2 ⟨hxY, hBY⟩
    have hcard := (hY _ hsub).2
    rw [Finset.card_insert_of_notMem (by simpa using hxB), hB] at hcard
    exact Nat.not_succ_le_self k hcard

/-- A matroid on the ordered ground set `Fin n` is a *positroid* if it is represented by the
columns of a full-rank real `d × n` matrix all of whose maximal minors are nonnegative, its
bases being the `d`-subsets with strictly positive minor (Ardila–Rincón–Williams, 2016,
following Postnikov, 2006). The property depends on the order of `Fin n`: it is preserved by
cyclic shifts of the ground set but not by arbitrary relabelings. -/
def IsPositroid {n : ℕ} (M : Matroid (Fin n)) : Prop :=
  M.E = Set.univ ∧
  ∃ (d : ℕ) (A : Matrix (Fin d) (Fin n) ℝ), A.rank = d ∧
    (∀ (I : Finset (Fin n)) (hI : I.card = d), 0 ≤ A.maximalMinor I hI) ∧
    ∀ B : Finset (Fin n), (M.IsBase ↑B ↔ ∃ hB : B.card = d, 0 < A.maximalMinor B hB)

/-- The uniform matroid `U_{d,n}` with the standard order on `[n]` is a positroid: the matrix
whose `(j, i)` entry is `(i+1)^j` is totally positive (its maximal minors are Vandermonde
determinants with strictly increasing nodes), so it realizes the "top cell" of the totally
nonnegative Grassmannian (Postnikov, 2006; Ardila–Rincón–Williams, 2016). -/
theorem isPositroid_unifOn {n d : ℕ} (hdn : d ≤ n) :
    (unifOn (Set.univ : Set (Fin n)) d).IsPositroid := by
  set A : Matrix (Fin d) (Fin n) ℝ := Matrix.of fun j i ↦ ((i : ℕ) + 1 : ℝ) ^ (j : ℕ) with hA
  have hminor : ∀ (I : Finset (Fin n)) (hI : I.card = d), 0 < A.maximalMinor I hI := by
    intro I hI
    have hmono : StrictMono fun p : Fin d ↦ (((I.orderEmbOfFin hI p : Fin n) : ℕ) + 1 : ℝ) := by
      intro p q hpq
      have hlt : ((I.orderEmbOfFin hI p : Fin n) : ℕ) < ((I.orderEmbOfFin hI q : Fin n) : ℕ) :=
        (I.orderEmbOfFin hI).strictMono hpq
      dsimp only
      exact_mod_cast Nat.add_lt_add_right hlt 1
    have hcalc : A.maximalMinor I hI =
        (Matrix.vandermonde fun p : Fin d ↦
          (((I.orderEmbOfFin hI p : Fin n) : ℕ) + 1 : ℝ)).det := by
      rw [Matrix.maximalMinor, ← Matrix.det_transpose]
      congr 1
    rw [hcalc]
    exact Matrix.det_vandermonde_pos hmono
  have hrank : A.rank = d := by
    obtain ⟨I₀, -, hI₀⟩ :=
      Finset.exists_superset_card_eq (s := (∅ : Finset (Fin n))) (by simp) (by simpa using hdn)
    have hle : A.rank ≤ d := by simpa using A.rank_le_card_height
    have hunit : IsUnit (A.submatrix id (I₀.orderEmbOfFin hI₀)) := by
      rw [Matrix.isUnit_iff_isUnit_det]
      exact (hminor I₀ hI₀).ne'.isUnit
    have hsub : (A.submatrix id (I₀.orderEmbOfFin hI₀)).rank = d := by
      rw [Matrix.rank_of_isUnit _ hunit, Fintype.card_fin]
    have hge : d ≤ A.rank := by
      have h := Matrix.rank_submatrix_le A id (I₀.orderEmbOfFin hI₀)
      rwa [hsub] at h
    exact le_antisymm hle hge
  refine ⟨by simp [unifOn], d, A, hrank, fun I hI ↦ (hminor I hI).le, fun B ↦ ?_⟩
  rw [unifOn_isBase_iff hdn]
  exact ⟨fun hB ↦ ⟨hB, hminor B hB⟩, fun ⟨hB, _⟩ ↦ hB⟩

end Matroid

/-- A *Grassmann necklace* of type `(d, n)` (Postnikov, 2006, §16): a cyclic sequence
`I₁, …, Iₙ` of `d`-subsets of `[n]` such that, indices modulo `n`, the set `I_{i+1}` is
obtained from `Iᵢ` by exchanging `i` for some element `j` when `i ∈ Iᵢ`, and equals `Iᵢ`
otherwise. The exchanged element `j` may be `i` itself, which allows `I_{i+1} = Iᵢ` in the
membership case as well. -/
structure GrassmannNecklace (n d : ℕ) [NeZero n] where
  /-- The `i`-th subset of the necklace. -/
  toFun : Fin n → Finset (Fin n)
  /-- Every entry of the necklace has exactly `d` elements. -/
  card_toFun : ∀ i : Fin n, (toFun i).card = d
  /-- When `i` is not in the `i`-th entry, the next entry is unchanged. -/
  toFun_succ_of_notMem : ∀ i : Fin n, i ∉ toFun i → toFun (i + 1) = toFun i
  /-- When `i` is in the `i`-th entry, the next entry exchanges `i` for a single element. -/
  exists_toFun_succ_of_mem : ∀ i : Fin n, i ∈ toFun i →
    ∃ j : Fin n, toFun (i + 1) = insert j ((toFun i).erase i)

namespace GrassmannNecklace

/-- The candidate basis family `𝓑(𝓘) = {B : Iⱼ ≤ⱼ B for all j}` attached to a Grassmann
necklace, where `≤ⱼ` is the Gale order for the cyclic shift by `j`. By Oh's theorem (Oh,
2011, proving Postnikov's conjecture) this is exactly the set of bases of the positroid
indexed by the necklace; that theorem is not yet formalized here. -/
def ohFamily {n d : ℕ} [NeZero n] (N : GrassmannNecklace n d) : Finset (Finset (Fin n)) :=
  (Finset.powersetCard d Finset.univ).filter fun B ↦
    ∀ j : Fin n, Finset.GaleShiftLE j (N.toFun j) B

/-- The Grassmann necklace of the uniform matroid `U_{2,4}`: the `i`-th entry is the cyclic
interval `{i, i+1}` (the Gale-minimal basis for the shifted order at `i`). -/
def uniformTwoFour : GrassmannNecklace 4 2 where
  toFun i := {i, i + 1}
  card_toFun := by decide
  toFun_succ_of_notMem := by decide
  exists_toFun_succ_of_mem := by decide

end GrassmannNecklace

/-- A *decorated permutation* of `[n]` (Postnikov, 2006): a permutation of `Fin n` together
with a coloring of each fixed point as either clockwise or counterclockwise, recorded here as
the set of clockwise fixed points. Under Postnikov's bijection with positroids, coloops
correspond to fixed points of one color and loops to fixed points of the other — which is why
the decoration is needed. -/
structure DecoratedPermutation (n : ℕ) where
  /-- The underlying permutation. -/
  toPerm : Equiv.Perm (Fin n)
  /-- The set of fixed points colored clockwise. -/
  clockwise : Finset (Fin n)
  /-- Only fixed points may be colored clockwise. -/
  apply_eq_self_of_mem_clockwise : ∀ i ∈ clockwise, toPerm i = i

namespace DecoratedPermutation

/-- The weak excedances of a decorated permutation: the positions `i` with `i < π i`,
together with the clockwise fixed points. Which fixed-point color counts varies by source;
this development (like the repository's Python implementation) counts the clockwise ones.
Under Postnikov's bijection the number of weak excedances is the rank of the positroid. -/
def weakExcedances {n : ℕ} (σ : DecoratedPermutation n) : Finset (Fin n) :=
  (Finset.univ.filter fun i ↦ i < σ.toPerm i) ∪ σ.clockwise

/-- The decorated permutation exchanging `0` and `1` in `[3]`, with `2` a clockwise fixed
point. -/
def swapZeroOne : DecoratedPermutation 3 where
  toPerm := Equiv.swap 0 1
  clockwise := {2}
  apply_eq_self_of_mem_clockwise := by decide

/-- `swapZeroOne` has two weak excedances (the excedance at `0` and the clockwise fixed point
`2`), so it indexes a rank-`2` positroid on `[3]`. -/
theorem card_weakExcedances_swapZeroOne : swapZeroOne.weakExcedances.card = 2 := by
  decide

end DecoratedPermutation

/-- A *Le-diagram* (Γ-diagram, Postnikov, 2006, §6): a filling of the cells of a Young
diagram with `0`s and `1`s such that no `0` has simultaneously a `1` above it in its column
and a `1` to its left in its row. The condition is stated positively: whenever a cell of the
diagram has a `1` above it and a `1` to its left, the cell itself is a `1`. Coordinates are
matrix-style `(row, column)` pairs, matching `YoungDiagram`. -/
structure LeDiagram (μ : YoungDiagram) where
  /-- The filling; `true` encodes a `1`. -/
  toFun : ℕ → ℕ → Bool
  /-- Cells outside the diagram are `0`. -/
  toFun_eq_false_of_notMem : ∀ i j : ℕ, (i, j) ∉ μ → toFun i j = false
  /-- The Γ-property: a cell with a `1` above it in its column and a `1` to its left in its
  row is itself a `1`. -/
  le_property : ∀ ⦃i₁ i₂ j₁ j₂ : ℕ⦄, i₁ < i₂ → j₁ < j₂ → (i₂, j₂) ∈ μ →
    toFun i₁ j₂ = true → toFun i₂ j₁ = true → toFun i₂ j₂ = true
