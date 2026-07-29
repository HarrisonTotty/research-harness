/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Combinatorics.Matroid.Circuit
import Mathlib.Combinatorics.Matroid.Closure
import Mathlib.Combinatorics.Matroid.Constructions
import Mathlib.Combinatorics.Matroid.Dual
import Mathlib.Combinatorics.Matroid.Loop
import Mathlib.Combinatorics.Matroid.Minor.Contract
import Mathlib.Combinatorics.Matroid.Rank.ENat
import Mathlib.Combinatorics.Matroid.Rank.Finite
import Mathlib.Tactic
import Theorems.Positroid

/-!
# Matroids

A *matroid* is a combinatorial structure abstracting the notion of independence, introduced by
Whitney (*On the abstract properties of linear dependence*, 1935) and independently by Nakasawa
(1935–1936). Mathlib's `Matroid` (under `Mathlib/Combinatorics/Matroid/`) is the library of
record: it stores a possibly infinite ground set `M.E` with a base predicate `M.IsBase` as the
primitive, and already provides independence (`Matroid.Indep`), circuits (`Matroid.IsCircuit`),
closure (`Matroid.closure`), flats (`Matroid.IsFlat`), loops and coloops, duality (`M✶`),
restriction, deletion, contraction, minors, direct sums, and the `ℕ∞`-valued rank `Matroid.eRk`.

This module fills the gaps between Mathlib's coverage and the standard finite-matroid
vocabulary: hyperplanes and their identification with cocircuit complements, parallel elements
and simple matroids, nullity, truncation, the contraction rank identity, and rank/duality
lemmas for the uniform matroids `Matroid.unifOn` defined in `Theorems.Positroid`.

## Main definitions

* `Matroid.IsHyperplane` — a hyperplane, as a maximal nonspanning subset of the ground set.
* `Matroid.Parallel` — parallel elements: nonloops with equal closures.
* `Matroid.Simple` — simple matroids (no loops and no distinct parallel pairs).
* `Matroid.nullity` — the `ℕ∞`-valued nullity of a set, via the dual of the restriction.
* `Matroid.truncateTo` — the truncation of a matroid to rank at most `k`.

## Main statements

* `Matroid.isHyperplane_iff_compl_isCocircuit` — hyperplanes are exactly the complements of
  cocircuits within the ground set.
* `Matroid.IsHyperplane.exists_isHyperplane_insert_inter_subset` — the hyperplane "elimination"
  axiom (H3), derived from circuit elimination in the dual.
* `Matroid.parallel_iff_isCircuit_pair` — distinct elements are parallel iff they form a
  two-element circuit.
* `Matroid.eRk_add_nullity_eq_encard` — the subtraction-free form of `n(X) = |X| - r(X)`.
* `Matroid.contract_eRk_add_eRk_eq_eRk_union` — the subtraction-free contraction rank formula
  `r_{M/C}(Y) + r(C) = r(Y ∪ C)`.
* `Matroid.mem_closure_iff_eRk_insert_eq` — the rank characterization of closure membership.
* `Matroid.IsColoop.eRk_ground_sdiff_add_one_eq_eRank` — the coloop rank identity
  `r(E \ {e}) + 1 = r(E)`, subtraction-free.
* `Matroid.unifOn_dual_eq` — uniform matroid duality `U_{k,n}✶ = U_{l,n}` for `k + l = n`,
  with the self-duality of `U_{2,4}` as the corollary `Matroid.unifOn_two_four_dual_eq`.

## Implementation notes

The specification for this module is the research graph's *Matroid* page; each declaration's
docstring records the claim it formalizes. Claims already covered by Mathlib at the pinned
revision are used, not restated — in particular the independence axioms (I1)–(I3)
(`IndepMatroid.ofFinite`), base exchange (`Matroid.isBase_exchange`), equicardinality of bases
(`Matroid.IsBase.encard_eq_encard_of_isBase`), the circuit axioms and (strong) circuit
elimination (`Matroid.IsCircuit.elimination`, `Matroid.IsCircuit.strong_elimination`), the rank
axioms (R1)–(R3) (`Matroid.eRk_le_encard`, `Matroid.eRk_mono`, `Matroid.eRk_submod`), the
closure axioms (CL1)–(CL4) (`Matroid.subset_closure`, `Matroid.closure_mono`,
`Matroid.closure_closure`, `Matroid.closure_exchange`), circuit–cocircuit orthogonality
(`Matroid.IsCircuit.inter_isCocircuit_ne_singleton`), the unique fundamental circuit
(`Matroid.fundCircuit`), and the duality, minor, and direct-sum calculus.

`IsHyperplane` is defined as a maximal nonspanning subset rather than as a corank-one flat:
the pinned Mathlib has no flat-order machinery, and this form both composes directly with
`Matroid.isCocircuit_iff_minimal_compl_nonspanning` and is provably a flat
(`IsHyperplane.isFlat`). `Parallel` is defined by equality of closures rather than by the
page's rank equation so that it composes with Mathlib's closure API; the page's definition is
recovered by `parallel_iff_eRk_pair_eq_one`. `nullity` follows the dual-of-restriction design
so that it needs no finiteness hypothesis. Statements involving the page's subtractions
(`r(E) - 1`, `|X| - r(X)`, `r(Y ∪ X) - r(X)`) are stated subtraction-free by adding on the
other side, per the repository's statement-hygiene rules.

## References

* [H. Whitney, *On the abstract properties of linear dependence*][whitney1935],
  American Journal of Mathematics **57** (1935), 509–533.
* [J. Oxley, *Matroid Theory*, 2nd ed.][oxley2011], Oxford, 2011.

## Tags

matroid, hyperplane, cocircuit, parallel, simple matroid, nullity, truncation, uniform matroid
-/

open Set

namespace Matroid

variable {α : Type*} {M : Matroid α} {C F₁ F₂ H H₁ H₂ K X Y I : Set α} {e f g : α}

/-! ### Flats -/

/-- The intersection of two flats is a flat. This is the binary case of the page's flat axiom
(F2) (Whitney, 1935; Oxley, 2011, §1.7); Mathlib has the indexed version
`Matroid.IsFlat.iInter`, which requires a nonempty index type. -/
theorem IsFlat.inter (h₁ : M.IsFlat F₁) (h₂ : M.IsFlat F₂) : M.IsFlat (F₁ ∩ F₂) := by
  rw [inter_eq_iInter]
  exact IsFlat.iInter fun b ↦ b.rec h₂ h₁

/-! ### Hyperplanes -/

/-- A *hyperplane* of a matroid is a maximal nonspanning subset of the ground set — a flat of
rank `r(M) - 1` (Oxley, 2011, §1.4). Since the pinned Mathlib has no flat-order machinery,
maximal-nonspanning is taken as the definition; `Matroid.IsHyperplane.isFlat` recovers
flatness, and `Matroid.isHyperplane_iff_compl_isCocircuit` the cocircuit-complement
characterization. -/
def IsHyperplane (M : Matroid α) (H : Set α) : Prop :=
  Maximal (fun X ↦ ¬ M.Spanning X ∧ X ⊆ M.E) H

/-- A hyperplane is not spanning. -/
theorem IsHyperplane.not_spanning (hH : M.IsHyperplane H) : ¬ M.Spanning H :=
  hH.prop.1

/-- A hyperplane is contained in the ground set. -/
theorem IsHyperplane.subset_ground (hH : M.IsHyperplane H) : H ⊆ M.E :=
  hH.prop.2

/-- The ground set is not a hyperplane. This is the page's hyperplane axiom (H1). -/
theorem ground_not_isHyperplane (M : Matroid α) : ¬ M.IsHyperplane M.E :=
  fun h ↦ h.not_spanning M.ground_spanning

/-- Distinct hyperplanes are incomparable: the set of hyperplanes is an antichain. This is the
page's hyperplane axiom (H2). -/
theorem IsHyperplane.eq_of_subset_isHyperplane (h₁ : M.IsHyperplane H₁) (h₂ : M.IsHyperplane H₂)
    (h : H₁ ⊆ H₂) : H₁ = H₂ :=
  h.antisymm (h₁.2 h₂.prop h)

/-- A hyperplane is a flat: the page's description of hyperplanes as "flats of rank
`r(M) - 1`" (Oxley, 2011, §1.4). -/
theorem IsHyperplane.isFlat (hH : M.IsHyperplane H) : M.IsFlat H := by
  have hHE : H ⊆ M.E := hH.subset_ground
  rw [isFlat_iff_closure_eq]
  refine (hH.2 ⟨fun hsp ↦ hH.not_spanning ?_, M.closure_subset_ground H⟩
      (M.subset_closure H hHE)).antisymm (M.subset_closure H hHE)
  rw [spanning_iff_closure_eq (M.closure_subset_ground H), closure_closure] at hsp
  rw [spanning_iff_closure_eq hHE]
  exact hsp

/-- A subset of the ground set is a hyperplane iff its complement is a cocircuit. This is the
page's claim that "hyperplanes are exactly the complements of the cocircuits" (Oxley, 2011,
Proposition 2.1.6). -/
theorem isHyperplane_iff_compl_isCocircuit (hH : H ⊆ M.E) :
    M.IsHyperplane H ↔ M.IsCocircuit (M.E \ H) := by
  rw [isCocircuit_iff_minimal_compl_nonspanning]
  constructor
  · rintro ⟨⟨hsp, -⟩, hmax⟩
    refine ⟨show ¬ M.Spanning (M.E \ (M.E \ H)) by rwa [sdiff_sdiff_cancel_left hH],
      fun X hX hXH ↦ ?_⟩
    have hHX : H ⊆ M.E \ X := subset_sdiff.2 ⟨hH, (subset_sdiff.1 hXH).2.symm⟩
    exact sdiff_subset_comm.1 (hmax ⟨hX, sdiff_subset⟩ hHX)
  · rintro ⟨hsp, hmin⟩
    rw [sdiff_sdiff_cancel_left hH] at hsp
    refine ⟨⟨hsp, hH⟩, fun X hX hHX ↦ ?_⟩
    have h' : M.E \ H ⊆ M.E \ X :=
      hmin (show ¬ M.Spanning (M.E \ (M.E \ X)) by
        rw [sdiff_sdiff_cancel_left hX.2]; exact hX.1) (sdiff_subset_sdiff_right hHX)
    intro x hxX
    by_contra hxH
    exact (h' ⟨hX.2 hxX, hxH⟩).2 hxX

/-- The complement of a hyperplane is a cocircuit. -/
theorem IsHyperplane.compl_isCocircuit (hH : M.IsHyperplane H) : M.IsCocircuit (M.E \ H) :=
  (isHyperplane_iff_compl_isCocircuit hH.subset_ground).1 hH

/-- The complement of a cocircuit is a hyperplane. -/
theorem IsCocircuit.compl_isHyperplane (hK : M.IsCocircuit K) : M.IsHyperplane (M.E \ K) := by
  have hKE : K ⊆ M.E := by simpa using hK.subset_ground
  rw [isHyperplane_iff_compl_isCocircuit sdiff_subset, sdiff_sdiff_cancel_left hKE]
  exact hK

/-- The hyperplane "elimination" axiom: for distinct hyperplanes `H₁, H₂` and an element `e`
of the ground set, some hyperplane contains `e` together with `H₁ ∩ H₂`. This is the page's
hyperplane axiom (H3), proved by dualizing circuit elimination. The page states the axiom for
`e ∈ E \ (H₁ ∪ H₂)`; the proof does not use `e ∉ H₁ ∪ H₂`, so the statement omits it — the
axiom's usual form is the special case. -/
theorem IsHyperplane.exists_isHyperplane_insert_inter_subset (h₁ : M.IsHyperplane H₁)
    (h₂ : M.IsHyperplane H₂) (hne : H₁ ≠ H₂) (heE : e ∈ M.E) :
    ∃ H₃, M.IsHyperplane H₃ ∧ insert e (H₁ ∩ H₂) ⊆ H₃ := by
  have hK₁ : M.IsCocircuit (M.E \ H₁) := h₁.compl_isCocircuit
  have hK₂ : M.IsCocircuit (M.E \ H₂) := h₂.compl_isCocircuit
  have hne' : M.E \ H₁ ≠ M.E \ H₂ := fun h ↦ hne <| by
    rw [← sdiff_sdiff_cancel_left h₁.subset_ground, h, sdiff_sdiff_cancel_left h₂.subset_ground]
  obtain ⟨K, hKsub, hK⟩ := hK₁.elimination hK₂ hne' e
  refine ⟨M.E \ K, IsCocircuit.compl_isHyperplane hK, insert_subset ?_ ?_⟩
  · exact ⟨heE, fun heK ↦ (hKsub heK).2 (mem_singleton e)⟩
  · intro x hx
    refine ⟨h₁.subset_ground hx.1, fun hxK ↦ ?_⟩
    obtain ⟨h | h, -⟩ := hKsub hxK
    · exact h.2 hx.1
    · exact h.2 hx.2

/-! ### Closure and rank -/

/-- An element of the ground set lies in the closure of a finite-rank set `X` iff adding it
does not increase the rank. This is the page's cryptomorphism `r → cl`:
`cl(X) = {e ∈ E : r(X ∪ {e}) = r(X)}` (Whitney, 1935, §4). The finite-rank hypothesis is
necessary: adding any element to an infinite-rank set preserves `eRk`. -/
theorem mem_closure_iff_eRk_insert_eq (hfin : M.IsRkFinite X) (hX : X ⊆ M.E) (he : e ∈ M.E) :
    e ∈ M.closure X ↔ M.eRk (insert e X) = M.eRk X := by
  constructor
  · intro h
    have hsub : insert e X ⊆ M.closure X := insert_subset h (M.subset_closure X hX)
    have hle : M.eRk (insert e X) ≤ M.eRk (M.closure X) := M.eRk_mono hsub
    rw [eRk_closure_eq] at hle
    exact hle.antisymm (M.eRk_mono (subset_insert e X))
  · intro h
    by_contra hcl
    have hplus := M.eRk_insert_eq_add_one (X := X) (e := e) ⟨he, hcl⟩
    rw [h] at hplus
    have h1 : M.eRk X + 1 ≤ M.eRk X := hplus.symm.le
    rw [ENat.add_one_le_iff hfin.eRk_lt_top.ne] at h1
    exact lt_irrefl _ h1

/-! ### Nullity -/

/-- The *nullity* of a set `X` is the corank of the restriction to `X` — for finite sets,
`n(X) = |X| - r(X)` (Whitney, 1935, §2, where it is written `n(N)`). Defined via the dual so
that no finiteness hypothesis is needed; the subtraction-free identity is
`Matroid.eRk_add_nullity_eq_encard`. -/
noncomputable def nullity (M : Matroid α) (X : Set α) : ℕ∞ :=
  (M ↾ X)✶.eRank

/-- Rank plus nullity is cardinality: the subtraction-free form of the page's definition
`n(X) = |X| - r(X)` (Whitney, 1935, §2). -/
theorem eRk_add_nullity_eq_encard (M : Matroid α) (X : Set α) :
    M.eRk X + M.nullity X = X.encard := by
  show (M ↾ X).eRank + (M ↾ X)✶.eRank = X.encard
  simpa using (M ↾ X).eRank_add_eRank_dual

/-- Independent sets are exactly the sets of nullity zero; this direction validates the
definition of `Matroid.nullity` against Whitney's. -/
theorem Indep.nullity_eq_zero (hI : M.Indep I) : M.nullity I = 0 := by
  rw [nullity, hI.restrict_eq_freeOn, freeOn_dual_eq, eRank_loopyOn]

/-! ### Coloops -/

/-- Removing a coloop from the ground set drops the rank by exactly one: the subtraction-free
form of the page's characterization `r(E \ {e}) = r(E) - 1` (Oxley, 2011, §1.5). Stated in
`ℕ∞`, so no finiteness hypothesis is needed. -/
theorem IsColoop.eRk_ground_sdiff_add_one_eq_eRank (he : M.IsColoop e) :
    M.eRk (M.E \ {e}) + 1 = M.eRank := by
  have hloop : M✶.IsLoop e := dual_isLoop_iff_isColoop.2 he
  have heE : e ∈ M.E := by simpa using hloop.mem_ground
  have h := M.eRk_dual_add_eRank (X := {e}) (singleton_subset_iff.2 heE)
  rw [hloop.eRk_eq, zero_add, encard_singleton] at h
  exact h.symm

/-- In a matroid of finite rank, the coloops are exactly the elements whose removal drops the
rank: the page's "equivalently" for coloops, as a biconditional (Oxley, 2011, §1.5). Finite
rank is necessary — in an infinite-rank matroid the rank equation holds vacuously at every
element of the ground set. -/
theorem isColoop_iff_eRk_ground_sdiff_add_one_eq_eRank [M.RankFinite] (heE : e ∈ M.E) :
    M.IsColoop e ↔ M.eRk (M.E \ {e}) + 1 = M.eRank := by
  refine ⟨IsColoop.eRk_ground_sdiff_add_one_eq_eRank, fun h ↦ ?_⟩
  have hd := M.eRk_dual_add_eRank (X := {e}) (singleton_subset_iff.2 heE)
  rw [encard_singleton, h] at hd
  have h0 : M✶.eRk {e} + M.eRank = 0 + M.eRank := by rwa [zero_add]
  have hzero : M✶.eRk {e} = 0 :=
    WithTop.add_right_cancel (M.eRank_ne_top_iff.2 inferInstance) h0
  have hloops : ({e} : Set α) ⊆ M✶.loops :=
    (eRk_eq_zero_iff (by simpa using singleton_subset_iff.2 heE)).1 hzero
  exact dual_isLoop_iff_isColoop.1 (hloops (mem_singleton e))

/-! ### Parallel elements and simple matroids -/

/-- Two elements are *parallel* if both are nonloops and they have the same closure
(Oxley, 2011, §1.1). The page defines parallel elements as nonloops with `r({e, f}) = 1`;
that characterization is `Matroid.parallel_iff_eRk_pair_eq_one`, and the closure form is
taken as the definition so that the relation composes with Mathlib's closure API. -/
def Parallel (M : Matroid α) (e f : α) : Prop :=
  M.IsNonloop e ∧ M.IsNonloop f ∧ M.closure {e} = M.closure {f}

/-- The left element of a parallel pair is a nonloop. -/
theorem Parallel.isNonloop_left (h : M.Parallel e f) : M.IsNonloop e :=
  h.1

/-- The right element of a parallel pair is a nonloop. -/
theorem Parallel.isNonloop_right (h : M.Parallel e f) : M.IsNonloop f :=
  h.2.1

/-- Parallel elements have equal closures. -/
theorem Parallel.closure_eq_closure (h : M.Parallel e f) : M.closure {e} = M.closure {f} :=
  h.2.2

/-- Every nonloop is parallel to itself: parallelism is reflexive on nonloops, part of the
page's claim that "parallelism is an equivalence relation on non-loops". -/
theorem IsNonloop.parallel_self (he : M.IsNonloop e) : M.Parallel e e :=
  ⟨he, he, rfl⟩

/-- Parallelism is symmetric, part of the page's claim that "parallelism is an equivalence
relation on non-loops". -/
theorem Parallel.symm (h : M.Parallel e f) : M.Parallel f e :=
  ⟨h.2.1, h.1, h.2.2.symm⟩

/-- Parallelism is transitive, part of the page's claim that "parallelism is an equivalence
relation on non-loops". -/
theorem Parallel.trans (hef : M.Parallel e f) (hfg : M.Parallel f g) : M.Parallel e g :=
  ⟨hef.1, hfg.2.1, hef.2.2.trans hfg.2.2⟩

/-- Distinct elements are parallel iff they form a two-element circuit: the page's
"equivalently `{e, f}` a circuit" (Oxley, 2011, §1.1). -/
theorem parallel_iff_isCircuit_pair (hef : e ≠ f) : M.Parallel e f ↔ M.IsCircuit {e, f} := by
  constructor
  · rintro ⟨he, hf, hcl⟩
    exact (he.closure_eq_closure_iff_isCircuit_of_ne hef).1 hcl
  · intro h
    have he : M.IsNonloop e := h.isNonloop_of_mem (nontrivial_pair hef) (mem_insert e {f})
    have hf : M.IsNonloop f := h.isNonloop_of_mem (nontrivial_pair hef) (by simp)
    exact ⟨he, hf, (he.closure_eq_closure_iff_isCircuit_of_ne hef).2 h⟩

/-- Elements are parallel iff both are nonloops and their pair has rank one: the page's
definition "non-loops `e, f` with `r({e, f}) = 1`". Also covers `e = f`, where it says that a
nonloop is a rank-one singleton. -/
theorem parallel_iff_eRk_pair_eq_one :
    M.Parallel e f ↔ M.IsNonloop e ∧ M.IsNonloop f ∧ M.eRk {e, f} = 1 := by
  obtain rfl | hne := eq_or_ne e f
  · rw [pair_eq_singleton, eRk_singleton_eq_one_iff]
    exact ⟨fun h ↦ ⟨h.1, h.1, h.1⟩, fun h ↦ h.1.parallel_self⟩
  constructor
  · rintro ⟨he, hf, hcl⟩
    refine ⟨he, hf, ?_⟩
    have hsub : ({e, f} : Set α) ⊆ M.closure {e} := insert_subset
      (M.mem_closure_self e he.mem_ground)
      (singleton_subset_iff.2 (hcl ▸ M.mem_closure_self f hf.mem_ground))
    have hle : M.eRk {e, f} ≤ M.eRk (M.closure {e}) := M.eRk_mono hsub
    rw [eRk_closure_eq, eRk_singleton_eq_one_iff.2 he] at hle
    have hge : M.eRk {e} ≤ M.eRk {e, f} := M.eRk_mono (singleton_subset_iff.2 (mem_insert e {f}))
    rw [eRk_singleton_eq_one_iff.2 he] at hge
    exact hle.antisymm hge
  · rintro ⟨he, hf, hr⟩
    refine ⟨he, hf, ?_⟩
    have hfe : f ∈ M.closure {e} := by
      by_contra hcl
      have hplus := M.eRk_insert_eq_add_one (X := {e}) (e := f) ⟨hf.mem_ground, hcl⟩
      rw [eRk_singleton_eq_one_iff.2 he, pair_comm f e, hr] at hplus
      norm_num at hplus
    exact (hf.closure_eq_of_mem_closure hfe).symm

/-- A matroid is *simple* (a "combinatorial geometry" in the Crapo–Rota sense) if it has no
loops and no parallel pairs of distinct elements — equivalently, if parallelism on the ground
set is the identity (Oxley, 2011, §1.1). Loops are excluded because a ground-set element `e`
with `M.Parallel e e` is a nonloop; see `Matroid.Simple.loopless`. -/
class Simple (M : Matroid α) : Prop where
  /-- On the ground set, parallelism coincides with equality. -/
  parallel_iff_eq : ∀ ⦃e f : α⦄, e ∈ M.E → (M.Parallel e f ↔ e = f)

/-- In a simple matroid, parallel elements are equal: "no parallel pairs". -/
theorem Simple.eq_of_parallel [M.Simple] (h : M.Parallel e f) : e = f :=
  (Simple.parallel_iff_eq h.isNonloop_left.mem_ground).1 h

/-- In a simple matroid, every ground-set element is a nonloop. -/
theorem Simple.isNonloop [M.Simple] (he : e ∈ M.E) : M.IsNonloop e :=
  ((Simple.parallel_iff_eq he).2 rfl).isNonloop_left

/-- A simple matroid has no loops: the "no loops" half of the page's definition of a simple
matroid. -/
instance Simple.loopless [M.Simple] : M.Loopless := by
  rw [loopless_iff_forall_isNonloop]
  exact fun e he ↦ Simple.isNonloop he

/-! ### Truncation -/

/-- The *truncation* of a matroid to rank `k`: the matroid on the same ground set whose
independent sets are the independent sets of `M` with at most `k` elements — the page's
operation "`T_k(M)` keeping independent sets of size `≤ k`". Built from
`IndepMatroid.ofBddAugment`, so `M` itself may be infinite. -/
def truncateTo (M : Matroid α) (k : ℕ) : Matroid α :=
  (IndepMatroid.ofBddAugment (E := M.E)
    (Indep := fun I ↦ M.Indep I ∧ I.encard ≤ (k : ℕ∞))
    (indep_empty := ⟨M.empty_indep, by simp⟩)
    (indep_subset := fun I J hJ hIJ ↦ ⟨hJ.1.subset hIJ, (encard_mono hIJ).trans hJ.2⟩)
    (indep_aug := by
      rintro I J ⟨hI, hIk⟩ ⟨hJ, hJk⟩ hIJ
      obtain ⟨e, he, hi⟩ := hI.exists_insert_of_encard_lt hJ hIJ
      exact ⟨e, he.1, he.2, hi,
        (encard_insert_le I e).trans ((Order.add_one_le_of_lt hIJ).trans hJk)⟩)
    (indep_bdd := ⟨k, fun I hI ↦ hI.2⟩)
    (subset_ground := fun I hI ↦ hI.1.subset_ground)).matroid

/-- The truncation has the same ground set as the original matroid. -/
@[simp]
theorem truncateTo_ground_eq (M : Matroid α) (k : ℕ) : (M.truncateTo k).E = M.E := by
  simp [truncateTo]

/-- A set is independent in the truncation `T_k(M)` iff it is independent in `M` and has at
most `k` elements: the page's "keeping independent sets of size `≤ k`". -/
@[simp]
theorem truncateTo_indep_iff {k : ℕ} :
    (M.truncateTo k).Indep I ↔ M.Indep I ∧ I.encard ≤ (k : ℕ∞) := by
  simp [truncateTo]

/-! ### Contraction and rank -/

/-- The contraction rank formula: the subtraction-free form of the page's
`r_{M/X}(Y) = r(Y ∪ X) - r(X)` (Whitney, 1935, §9; Oxley, 2011, Proposition 3.1.6), stated
for sets `Y` in the ground set of the contraction. -/
theorem contract_eRk_add_eRk_eq_eRk_union (hC : C ⊆ M.E) (hY : Y ⊆ M.E \ C) :
    (M ／ C).eRk Y + M.eRk C = M.eRk (Y ∪ C) := by
  obtain ⟨I, hI⟩ := M.exists_isBasis C hC
  obtain ⟨J, hJ⟩ := (M ／ C).exists_isBasis Y (by rwa [contract_ground])
  have h1 : (M ／ C).eRk Y = J.encard := hJ.eRk_eq_encard
  have h2 : M.eRk C = I.encard := hI.eRk_eq_encard
  have hJ' : (M ／ I).IsBasis J Y := by
    rw [hI.contract_eq_contract_delete, delete_isBasis_iff] at hJ
    exact hJ.1
  have hbasis : M.IsBasis (J ∪ I) (Y ∪ I) :=
    hI.indep.union_isBasis_union_of_contract_isBasis hJ'
  have hcl : M.closure (Y ∪ C) = M.closure (Y ∪ I) := by
    rw [← closure_union_closure_right_eq, ← hI.closure_eq_closure,
      closure_union_closure_right_eq]
  have hrk : M.eRk (Y ∪ C) = M.eRk (Y ∪ I) := by
    rw [← eRk_closure_eq, hcl, eRk_closure_eq]
  have hdisj : Disjoint J I :=
    ((subset_sdiff.1 hY).2.mono hJ.subset hI.subset)
  rw [h1, h2, hrk, hbasis.eRk_eq_encard, encard_union_eq hdisj]

/-! ### Uniform matroids -/

section Uniform

variable [DecidableEq α] {k l n : ℕ}

/-- The ground set of the uniform matroid `Matroid.unifOn E k` is `E`. -/
@[simp]
theorem unifOn_ground_eq (E : Set α) (k : ℕ) : (unifOn E k).E = E := by
  simp [unifOn]

/-- A set of any cardinality is independent in `Matroid.unifOn E k` iff it lies in `E` and has
at most `k` elements: the page's "independent sets are all subsets of size `≤ r`", extended
from the finite characterization `Matroid.unifOn_indep_finset_iff` of `Theorems.Positroid`. -/
theorem unifOn_indep_iff {E I : Set α} {k : ℕ} :
    (unifOn E k).Indep I ↔ I ⊆ E ∧ I.encard ≤ (k : ℕ∞) := by
  rw [unifOn, IndepMatroid.matroid_indep_iff, IndepMatroid.ofFinset_indep']
  constructor
  · intro h
    refine ⟨fun x hx ↦ (h {x} (by simpa using hx)).1 (by simp), ?_⟩
    by_contra hlt
    rw [not_le] at hlt
    obtain ⟨t, hts, htcard⟩ :=
      exists_subset_encard_eq (Order.add_one_le_of_lt hlt)
    have htfin : t.Finite := finite_of_encard_eq_coe (k := k + 1) (by exact_mod_cast htcard)
    have hcard := (h htfin.toFinset (by simpa using hts)).2
    have henc : ((k : ℕ∞) + 1) ≤ (k : ℕ∞) := by
      rw [← htcard, htfin.encard_eq_coe_toFinset_card]
      exact_mod_cast hcard
    rw [ENat.add_one_le_iff (by simp : (k : ℕ∞) ≠ ⊤)] at henc
    exact lt_irrefl _ henc
  · rintro ⟨hIE, hIk⟩ J hJ
    refine ⟨hJ.trans hIE, ?_⟩
    have h := (encard_mono hJ).trans hIk
    rwa [encard_coe_eq_coe_finsetCard, Nat.cast_le] at h

/-- The rank-zero uniform matroid is the loopy matroid: the page's "`U_{0,n}` the loopy
matroid (every element a loop)". -/
theorem unifOn_zero_eq_loopyOn (E : Set α) : unifOn E 0 = loopyOn E := by
  refine ext_indep (by simp) fun I hI ↦ ?_
  rw [unifOn_indep_iff, loopyOn_indep_iff]
  constructor
  · rintro ⟨-, h⟩
    rwa [Nat.cast_zero, le_zero_iff, encard_eq_zero] at h
  · rintro rfl
    simp

/-- The full-rank uniform matroid is the free matroid: the page's "`U_{n,n}` is the free
matroid". -/
theorem unifOn_eq_freeOn {E : Set α} {k : ℕ} (hE : E.encard ≤ (k : ℕ∞)) :
    unifOn E k = freeOn E := by
  refine ext_indep (by simp) fun I hI ↦ ?_
  rw [unifOn_indep_iff, freeOn_indep_iff]
  exact and_iff_left_of_imp fun h ↦ (encard_mono h).trans hE

end Uniform

section UniformFin

variable {k l n : ℕ}

/-- The uniform matroid `U_{k,n}` has rank `k` for `k ≤ n`: the page's "`U_{r,n}` … rank `r`
on `n` elements". -/
theorem unifOn_eRank_eq (hkn : k ≤ n) : (unifOn (univ : Set (Fin n)) k).eRank = (k : ℕ∞) := by
  obtain ⟨B, -, hB⟩ := Finset.exists_subset_card_eq
    (show k ≤ (Finset.univ : Finset (Fin n)).card by simpa using hkn)
  have hbase : (unifOn (univ : Set (Fin n)) k).IsBase ↑B := (unifOn_isBase_iff hkn).2 hB
  rw [← hbase.eRk_eq_eRank, hbase.indep.eRk_eq_encard, encard_coe_eq_coe_finsetCard, hB]

/-- Uniform matroid duality: `U_{k,n}✶ = U_{l,n}` when `k + l = n`. This is the page's
"`U_{r,n}✶ = U_{n-r,n}` — a clean duality test", stated additively to avoid natural
subtraction. -/
theorem unifOn_dual_eq (hkl : k + l = n) :
    (unifOn (univ : Set (Fin n)) k)✶ = unifOn univ l := by
  have hkn : k ≤ n := hkl ▸ Nat.le_add_right k l
  have hln : l ≤ n := hkl ▸ Nat.le_add_left l k
  refine ext_isBase (by simp) fun B hB ↦ ?_
  obtain ⟨B', rfl⟩ : ∃ B' : Finset (Fin n), ↑B' = B :=
    ⟨B.toFinite.toFinset, B.toFinite.coe_toFinset⟩
  rw [dual_isBase_iff (show (↑B' : Set (Fin n)) ⊆ (unifOn (univ : Set (Fin n)) k).E by simp)]
  have hcompl : (unifOn (univ : Set (Fin n)) k).E \ ↑B' = ↑(B'ᶜ) := by
    rw [unifOn_ground_eq, Finset.coe_compl, ← compl_eq_univ_sdiff]
  rw [hcompl, unifOn_isBase_iff hkn, unifOn_isBase_iff hln, Finset.card_compl,
    Fintype.card_fin]
  have hle : B'.card ≤ n := by simpa using B'.card_le_univ
  omega

/-- The uniform matroid `U_{2,4}` is self-dual: part of the page's canonical example
"`U_{2,4}` — … rank 2 on 4 elements, self-dual". -/
theorem unifOn_two_four_dual_eq : (unifOn (univ : Set (Fin 4)) 2)✶ = unifOn univ 2 :=
  unifOn_dual_eq rfl

end UniformFin

end Matroid
