/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic
import Theorems.Positroid

/-!
# Grassmann necklaces

A *Grassmann necklace* of type `(d, n)` (Postnikov, *Total positivity, Grassmannians, and
networks*, 2006, Definition 16.1) is a cyclic sequence `I₁, …, Iₙ` of `d`-element subsets of
`[n]` in which each entry is obtained from the previous one by the smallest possible cyclic
update: passing from `Iᵢ` to `Iᵢ₊₁`, the element `i` is dropped (if present) and at most one
new element enters. The structure `GrassmannNecklace` itself is defined in
`Theorems.Positroid`, together with the Gale orders `Finset.GaleLE` / `Finset.GaleShiftLE` and
the candidate basis family `GrassmannNecklace.ohFamily`. This module develops the necklace
theory on top of those definitions.

## Main definitions

* `GrassmannNecklace.ofTransitions` — build a necklace from the transition conditions alone,
  with the cardinality of every entry *derived* rather than assumed.
* `GrassmannNecklace.ofEraseSubsetSucc` — build a necklace from the containment form of the
  axioms (Oh–Postnikov–Speyer, 2015, Definition 4.1).
* `GrassmannNecklace.const` — the constant necklace `(I, I, …, I)`.
* `GrassmannNecklace.instPartialOrder` — Lam's partial order on necklaces (Lam, 2014, §6.3):
  `𝓘 ≤ 𝓘'` iff `Iₐ ≤ₐ I'ₐ` in the shifted Gale order for every `a`.
* `GrassmannNecklace.cyclicShift` — the cyclic shift of a necklace, rotating both the entries
  and their indexing.
* `GrassmannNecklace.jugglingState` — the juggling states of Knutson–Lam–Speyer (2013, §3):
  the entries rotated back to a fixed window.
* `GrassmannNecklace.uniform` — the cyclic-interval necklace `Iᵢ = {i, i+1, …, i+d-1}` of the
  uniform matroid (the "top cell" necklace).
* `GrassmannNecklace.postnikovFigure16_1`, `GrassmannNecklace.ohWorkedExample`,
  `GrassmannNecklace.nonNecklace` — the canonical fixtures of the research page, 0-indexed.

## Main statements

* `GrassmannNecklace.card_eq_card_of_transitions` — the transition conditions alone force all
  entries to have the same cardinality (Postnikov's remark after Definition 16.1: "In
  particular, we have `|I₁| = ⋯ = |Iₙ|`").
* `GrassmannNecklace.instUniqueOfRankZero` / `instUniqueOfFullRank` — for `d = 0` and `d = n`
  there is exactly one necklace.
* `GrassmannNecklace.erase_subset_toFun_succ` — the containment form of the transition
  conditions holds in every necklace.
* `GrassmannNecklace.image_sub_one_erase_zero_jugglingState_subset` — the juggling-state
  condition `J_{r+1} ⊇ (J_r \ {0}) - 1` of Knutson–Lam–Speyer (2013, §3.3).
* `Finset.galeLE_filter_val_lt` — the initial segment `{0, …, d-1}` is Gale-below every
  `d`-element subset of `Fin n`.
* `GrassmannNecklace.isBot_uniform` — the cyclic-interval necklace is the least element of
  Lam's order; equivalently, the unique top element of Postnikov's circular Bruhat order
  (Postnikov, 2006, Lemma 17.6 — see the implementation notes on order direction).
* `GrassmannNecklace.ohFamily_uniform` — the basis family of the cyclic-interval necklace is
  the full family of `d`-subsets, the bases of the uniform matroid `U_{d,n}` (OPS, 2015, §4).
* `GrassmannNecklace.toFun_ne_nonNecklace` — the page's rejection fixture is not the entry
  function of any Grassmann necklace.

## Implementation notes

Everything here is 0-indexed: `Fin n` plays the role of Postnikov's `[n] = {1, …, n}`, so all
fixture data from the literature appears with every element shifted down by one, and the
transition at index `i` relates `toFun i` to `toFun (i + 1)` with wraparound provided by the
group structure on `Fin n`.

Postnikov's Definition 16.1 does not assume the entries have equal cardinality — he remarks
that it follows. The structure `GrassmannNecklace` of `Theorems.Positroid` carries the
cardinality as a field for convenience; `card_eq_card_of_transitions` proves the remark, and
the constructor `ofTransitions` discharges the field from the transition conditions alone, so
nothing is lost by storing it.

On order direction: Lam (2014, §6.3) orders necklaces by `𝓘 ≤ 𝓘' ↔ ∀ a, Iₐ ≤ₐ I'ₐ`, and
orders bounded affine permutations by the *dual* of the Bruhat order, making his necklace
order the dual of Postnikov's circular Bruhat order on the same data. Postnikov's Lemma 17.6
("`CB_{kn}` has a unique top element", the cyclic-interval necklace) therefore appears here as
`isBot_uniform`: in the entrywise Gale order the cyclic-interval entries are minimal, since
`{a, a+1, …, a+d-1}` is the Gale-least `d`-subset for the order shifted to start at `a`.

`GrassmannNecklace.cyclicShift` proves only that the rotated sequence is again a necklace.
That it agrees with the necklace of the cyclically shifted positroid (Ardila–Rincón–Williams,
2016, Lemma 3.3) is not formalized: it needs the necklace-of-a-matroid construction
(Postnikov's Lemma 16.3), which is recorded as backlog together with Lemma 16.2, Oh's
theorem, and the geometric statements — see the research notes.

`Finset.GaleLE` is decidable but not kernel-evaluable (`Finset.sort` is defined by
well-founded recursion), so `ohFamily` computations on concrete necklaces cannot close with
`decide`; `ohFamily_uniform` is proved structurally instead, and the corresponding equality
for `ohWorkedExample` is backlog.

## References

* [A. Postnikov, *Total positivity, Grassmannians, and networks*][postnikov2006],
  arXiv:math/0609764. Definition 16.1, the cardinality remark, §17, Lemma 17.6, Figure 16.1.
* [S. Oh, *Positroids and Schubert matroids*][oh2011], arXiv:0803.1018. The worked example
  and the basis family `𝓑(𝓘)`.
* [A. Knutson, T. Lam, D. Speyer, *Positroid varieties: juggling and geometry*][kls2013],
  arXiv:0903.3694. Juggling states (§3).
* [S. Oh, A. Postnikov, D. Speyer, *Weak separation and plabic graphs*][ops2015],
  arXiv:1109.4434. The containment form (Definition 4.1) and the top-cell necklace (§4).
* [T. Lam, *Totally nonnegative Grassmannian and Grassmann polytopes*][lam2014],
  arXiv:1506.00603. The partial order on necklaces (§6.3).

## Tags

Grassmann necklace, positroid, Gale order, cyclic shift, juggling state, total positivity
-/

open Finset

namespace Finset

/-- The shifted Gale order is reflexive. -/
theorem GaleShiftLE.refl {n : ℕ} (i : Fin n) (S : Finset (Fin n)) : GaleShiftLE i S S :=
  GaleLE.refl _

/-- The shifted Gale order is transitive. -/
theorem GaleShiftLE.trans {n : ℕ} {i : Fin n} {S T U : Finset (Fin n)}
    (hST : GaleShiftLE i S T) (hTU : GaleShiftLE i T U) : GaleShiftLE i S U :=
  GaleLE.trans hST hTU

/-- The shifted Gale order is antisymmetric: comparability both ways forces equality of the
underlying sets, since shifting by `· - i` is injective. -/
theorem GaleShiftLE.antisymm {n : ℕ} {i : Fin n} {S T : Finset (Fin n)}
    (hST : GaleShiftLE i S T) (hTS : GaleShiftLE i T S) : S = T := by
  haveI : NeZero n := ⟨i.pos.ne'⟩
  exact Finset.image_injective (Equiv.subRight i).injective (GaleLE.antisymm hST hTS)

/-- A strictly monotone function between `Fin` types grows at least linearly: the `j`-th
value is at least `j`. -/
private theorem le_val_of_strictMono {n k : ℕ} {f : Fin k → Fin n} (hf : StrictMono f)
    (j : Fin k) : (j : ℕ) ≤ (f j : ℕ) := by
  suffices h : ∀ (m : ℕ) (hm : m < k), m ≤ (f ⟨m, hm⟩ : ℕ) by
    simpa using h j.val j.isLt
  intro m
  induction m with
  | zero => exact fun hm ↦ Nat.zero_le _
  | succ p ih =>
      intro hm
      have hp : p < k := Nat.lt_of_succ_lt hm
      have h1 : p ≤ (f ⟨p, hp⟩ : ℕ) := ih hp
      have h2 : (f ⟨p, hp⟩ : ℕ) < (f ⟨p + 1, hm⟩ : ℕ) := hf (by simp [Fin.mk_lt_mk])
      omega

/-- The initial segment `{0, …, d-1}` of `Fin n`, as a filter, has exactly `d` elements. -/
theorem card_filter_val_lt {n d : ℕ} (hdn : d ≤ n) :
    (univ.filter fun x : Fin n ↦ (x : ℕ) < d).card = d := by
  have h : (univ.filter fun x : Fin n ↦ (x : ℕ) < d)
      = (univ : Finset (Fin d)).map (Fin.castLEEmb hdn) := by
    ext x
    simp only [mem_filter, mem_univ, true_and, Finset.mem_map]
    constructor
    · intro hx
      refine ⟨⟨x.val, hx⟩, ?_⟩
      ext
      simp [Fin.castLEEmb]
    · rintro ⟨y, -, rfl⟩
      simp [Fin.castLEEmb]
  rw [h, Finset.card_map, Finset.card_univ, Fintype.card_fin]

/-- The initial segment `{0, …, d-1}` of `Fin n` is below every `d`-element subset in the
Gale order: listing a `d`-subset in increasing order, its `j`-th element is at least `j`.
This is why the cyclic-interval entries of the top-cell necklace are Gale-minimal
(Postnikov, 2006, §17; Oh–Postnikov–Speyer, 2015, §4). -/
theorem galeLE_filter_val_lt {n d : ℕ} (hdn : d ≤ n) {S : Finset (Fin n)} (hS : S.card = d) :
    GaleLE (univ.filter fun x : Fin n ↦ (x : ℕ) < d) S := by
  have hDcard : (univ.filter fun x : Fin n ↦ (x : ℕ) < d).card = d := card_filter_val_lt hdn
  refine ⟨hS.trans hDcard.symm, fun j ↦ ?_⟩
  have hf : StrictMono fun x : Fin (univ.filter fun x : Fin n ↦ (x : ℕ) < d).card ↦
      (⟨(x : ℕ), lt_of_lt_of_le (lt_of_lt_of_le x.isLt hDcard.le) hdn⟩ : Fin n) := by
    intro a b hab
    exact Fin.mk_lt_mk.mpr (Fin.lt_def.mp hab)
  have hmem : ∀ x : Fin (univ.filter fun x : Fin n ↦ (x : ℕ) < d).card,
      (OrderEmbedding.ofStrictMono _ hf) x ∈ (univ.filter fun x : Fin n ↦ (x : ℕ) < d) := by
    intro x
    simp only [mem_filter, mem_univ, true_and]
    show ((⟨(x : ℕ), _⟩ : Fin n) : ℕ) < d
    exact lt_of_lt_of_le x.isLt hDcard.le
  have huniq := Finset.orderEmbOfFin_unique' rfl hmem
  have hstep : (univ.filter fun x : Fin n ↦ (x : ℕ) < d).orderEmbOfFin rfl j
      = (OrderEmbedding.ofStrictMono _ hf) j := (DFunLike.congr_fun huniq j).symm
  rw [Fin.le_def, hstep]
  show (j : ℕ) ≤ _
  exact le_val_of_strictMono (S.orderEmbOfFin _).strictMono j

end Finset

namespace GrassmannNecklace

variable {n d : ℕ} [NeZero n]

/-- Two Grassmann necklaces with the same entries are equal. -/
@[ext]
theorem ext {N M : GrassmannNecklace n d} (h : ∀ i : Fin n, N.toFun i = M.toFun i) : N = M := by
  cases N
  cases M
  simp only [mk.injEq]
  exact funext h

/-! ### The cardinality of the entries is determined by the transitions -/

/-- The transition conditions (N1)–(N2) alone force all entries of a Grassmann necklace to
have the same cardinality: Postnikov's remark after Definition 16.1, "In particular, we have
`|I₁| = ⋯ = |Iₙ|`" (Postnikov, 2006). Each transition weakly decreases the cardinality, and
following the transitions all the way around the cycle returns to the starting entry, so no
transition can decrease it. -/
theorem card_eq_card_of_transitions {I : Fin n → Finset (Fin n)}
    (h_mem : ∀ i : Fin n, i ∈ I i → ∃ j : Fin n, I (i + 1) = insert j ((I i).erase i))
    (h_notMem : ∀ i : Fin n, i ∉ I i → I (i + 1) = I i) (i j : Fin n) :
    (I i).card = (I j).card := by
  have step : ∀ a : Fin n, (I (a + 1)).card ≤ (I a).card := by
    intro a
    by_cases ha : a ∈ I a
    · obtain ⟨b, hb⟩ := h_mem a ha
      calc (I (a + 1)).card = (insert b ((I a).erase a)).card := by rw [hb]
        _ ≤ ((I a).erase a).card + 1 := Finset.card_insert_le _ _
        _ = (I a).card := Finset.card_erase_add_one ha
    · rw [h_notMem a ha]
  have key : ∀ (m : ℕ) (a x : Fin n), (x : ℕ) = m → (I (a + x)).card ≤ (I a).card := by
    intro m
    induction m with
    | zero =>
        intro a x hx
        obtain rfl : x = 0 := Fin.ext (by simp [hx])
        simp
    | succ p ih =>
        intro a x hx
        have hx0 : x ≠ 0 := by
          intro h
          rw [h] at hx
          simp at hx
        have hy : ((x - 1 : Fin n) : ℕ) = p := by
          rw [Fin.val_sub_one_of_ne_zero hx0, hx]
          omega
        have hsplit : a + x = (a + (x - 1)) + 1 := by
          rw [add_assoc, sub_add_cancel]
        calc (I (a + x)).card = (I ((a + (x - 1)) + 1)).card := by rw [hsplit]
          _ ≤ (I (a + (x - 1))).card := step _
          _ ≤ (I a).card := ih a (x - 1) hy
  have le_all : ∀ a b : Fin n, (I b).card ≤ (I a).card := by
    intro a b
    have hab : a + (b - a) = b := by abel
    calc (I b).card = (I (a + (b - a))).card := by rw [hab]
      _ ≤ (I a).card := key _ a (b - a) rfl
  exact le_antisymm (le_all j i) (le_all i j)

/-- Build a Grassmann necklace from the transition conditions (N1)–(N2) alone, deriving the
cardinality of every entry from the cardinality of the entry at `0` via
`card_eq_card_of_transitions`. This is Postnikov's Definition 16.1 as printed — with equal
cardinality a consequence, not an axiom (Postnikov, 2006). -/
def ofTransitions (I : Fin n → Finset (Fin n))
    (h_mem : ∀ i : Fin n, i ∈ I i → ∃ j : Fin n, I (i + 1) = insert j ((I i).erase i))
    (h_notMem : ∀ i : Fin n, i ∉ I i → I (i + 1) = I i) (hd : (I 0).card = d) :
    GrassmannNecklace n d where
  toFun := I
  card_toFun i := (card_eq_card_of_transitions h_mem h_notMem i 0).trans hd
  toFun_succ_of_notMem := h_notMem
  exists_toFun_succ_of_mem := h_mem

/-- The entries of `ofTransitions` are the given sequence. -/
@[simp]
theorem ofTransitions_toFun {I : Fin n → Finset (Fin n)}
    (h_mem : ∀ i : Fin n, i ∈ I i → ∃ j : Fin n, I (i + 1) = insert j ((I i).erase i))
    (h_notMem : ∀ i : Fin n, i ∉ I i → I (i + 1) = I i) (hd : (I 0).card = d) :
    (ofTransitions I h_mem h_notMem hd).toFun = I :=
  rfl

/-! ### The containment form of the axioms -/

/-- In every Grassmann necklace, each entry with its index removed is contained in the next
entry. This is the forward half of the containment form of the axioms (Oh–Postnikov–Speyer,
2015, Definition 4.1). -/
theorem erase_subset_toFun_succ (N : GrassmannNecklace n d) (i : Fin n) :
    (N.toFun i).erase i ⊆ N.toFun (i + 1) := by
  by_cases hi : i ∈ N.toFun i
  · obtain ⟨j, hj⟩ := N.exists_toFun_succ_of_mem i hi
    rw [hj]
    exact Finset.subset_insert _ _
  · rw [N.toFun_succ_of_notMem i hi]
    exact Finset.erase_subset _ _

/-- Build a Grassmann necklace from the containment form of the axioms: a sequence of
`d`-element subsets with `Iᵢ₊₁ ⊇ Iᵢ \ {i}` for all `i`, and `Iᵢ₊₁ = Iᵢ` whenever `i ∉ Iᵢ`
(Oh–Postnikov–Speyer, 2015, Definition 4.1). The exchange condition (N1) is recovered by a
cardinality count: `Iᵢ \ {i}` has one element fewer than `Iᵢ₊₁`, so exactly one element
enters. -/
def ofEraseSubsetSucc (I : Fin n → Finset (Fin n)) (hcard : ∀ i : Fin n, (I i).card = d)
    (hsub : ∀ i : Fin n, (I i).erase i ⊆ I (i + 1))
    (h_notMem : ∀ i : Fin n, i ∉ I i → I (i + 1) = I i) : GrassmannNecklace n d where
  toFun := I
  card_toFun := hcard
  toFun_succ_of_notMem := h_notMem
  exists_toFun_succ_of_mem i hi := by
    have hlt : ((I i).erase i).card < (I (i + 1)).card := by
      rw [hcard (i + 1)]
      calc ((I i).erase i).card < (I i).card := Finset.card_erase_lt_of_mem hi
        _ = d := hcard i
    obtain ⟨j, hj1, hj2⟩ := Finset.exists_mem_notMem_of_card_lt_card hlt
    refine ⟨j, (Finset.eq_of_subset_of_card_le (Finset.insert_subset hj1 (hsub i)) ?_).symm⟩
    rw [Finset.card_insert_of_notMem hj2, Finset.card_erase_add_one hi, hcard i, hcard (i + 1)]

/-! ### Constant necklaces and the extreme ranks -/

/-- The constant necklace `(I, I, …, I)` at a `d`-element subset `I`: always a valid
Grassmann necklace, taking `j = i` in the exchange condition (the page's "constant necklace";
by Lam, 2014, Lemma 8.3 it is the necklace of the single-basis positroid `{I}`, a statement
not yet formalized here). -/
def const (I : Finset (Fin n)) (hI : I.card = d) : GrassmannNecklace n d where
  toFun _ := I
  card_toFun _ := hI
  toFun_succ_of_notMem _ _ := rfl
  exists_toFun_succ_of_mem i hi := ⟨i, (Finset.insert_erase hi).symm⟩

/-- The entries of the constant necklace. -/
@[simp]
theorem const_toFun (I : Finset (Fin n)) (hI : I.card = d) (i : Fin n) :
    (const I hI).toFun i = I :=
  rfl

/-- There is exactly one Grassmann necklace of type `(0, n)`: the constant necklace at `∅`
(the page's "`k = 0` gives the single necklace `(∅, …, ∅)`"). -/
instance instUniqueOfRankZero : Unique (GrassmannNecklace n 0) where
  default := const ∅ Finset.card_empty
  uniq N := ext fun i ↦ Finset.card_eq_zero.mp (N.card_toFun i)

/-- There is exactly one Grassmann necklace of type `(n, n)`: the constant necklace at the
full ground set (the page's "`k = n` the single necklace `([n], …, [n])`"). -/
instance instUniqueOfFullRank : Unique (GrassmannNecklace n n) where
  default := const Finset.univ (by simp)
  uniq N := ext fun i ↦
    Finset.eq_univ_of_card _ ((N.card_toFun i).trans (Fintype.card_fin n).symm)

/-! ### Lam's partial order -/

/-- Lam's partial order on Grassmann necklaces of type `(d, n)` (Lam, 2014, §6.3):
`𝓘 ≤ 𝓘'` iff every entry satisfies `Iₐ ≤ₐ I'ₐ` in the Gale order shifted to start at `a`. -/
instance instPartialOrder : PartialOrder (GrassmannNecklace n d) where
  le N M := ∀ a : Fin n, Finset.GaleShiftLE a (N.toFun a) (M.toFun a)
  le_refl N a := Finset.GaleShiftLE.refl a (N.toFun a)
  le_trans N M P hNM hMP a := (hNM a).trans (hMP a)
  le_antisymm N M hNM hMN := ext fun a ↦ Finset.GaleShiftLE.antisymm (hNM a) (hMN a)

/-- Unfolding lemma for Lam's order on necklaces. -/
theorem le_def {N M : GrassmannNecklace n d} :
    N ≤ M ↔ ∀ a : Fin n, Finset.GaleShiftLE a (N.toFun a) (M.toFun a) :=
  Iff.rfl

/-! ### Cyclic shift -/

/-- The cyclic shift of a Grassmann necklace: entry `i` of the shifted necklace is entry
`i - 1` of the original with every element increased by one. This realizes, at the necklace
level, the cyclic-shift closure of the theory (Ardila–Rincón–Williams, 2016, Lemma 3.3 for
positroids); only necklace-hood of the shifted sequence is proved here — see the
implementation notes. -/
def cyclicShift (N : GrassmannNecklace n d) : GrassmannNecklace n d where
  toFun i := (N.toFun (i - 1)).image (· + 1)
  card_toFun i :=
    (Finset.card_image_of_injective _ (add_left_injective 1)).trans (N.card_toFun _)
  toFun_succ_of_notMem i hi := by
    have hmem : i - 1 ∉ N.toFun (i - 1) := by
      intro h
      have h' := Finset.mem_image_of_mem (· + 1) h
      rw [sub_add_cancel] at h'
      exact hi h'
    have h2 := N.toFun_succ_of_notMem (i - 1) hmem
    rw [sub_add_cancel] at h2
    rw [add_sub_cancel_right, h2]
  exists_toFun_succ_of_mem i hi := by
    have hmem : i - 1 ∈ N.toFun (i - 1) := by
      obtain ⟨x, hx, hxi⟩ := Finset.mem_image.mp hi
      obtain rfl : x = i - 1 := eq_sub_of_add_eq hxi
      exact hx
    obtain ⟨j, hj⟩ := N.exists_toFun_succ_of_mem (i - 1) hmem
    rw [sub_add_cancel] at hj
    refine ⟨j + 1, ?_⟩
    rw [add_sub_cancel_right, hj, Finset.image_insert,
      Finset.image_erase (add_left_injective 1), sub_add_cancel]

/-- The entries of the cyclic shift. -/
@[simp]
theorem cyclicShift_toFun (N : GrassmannNecklace n d) (i : Fin n) :
    N.cyclicShift.toFun i = (N.toFun (i - 1)).image (· + 1) :=
  rfl

/-! ### Juggling states -/

/-- The juggling states of a Grassmann necklace (Knutson–Lam–Speyer, 2013, §3): the `r`-th
state is the `r`-th entry rotated back to a fixed window, `J_r = χ^{-r}(I_r)` for the long
cycle `χ`. Each state is read as the set of scheduled landing times of `d` airborne balls. -/
def jugglingState (N : GrassmannNecklace n d) (r : Fin n) : Finset (Fin n) :=
  (N.toFun r).image (· - r)

/-- The `0`-th juggling state is the `0`-th entry. -/
@[simp]
theorem jugglingState_zero (N : GrassmannNecklace n d) : N.jugglingState 0 = N.toFun 0 := by
  simp [jugglingState]

/-- The juggling-state condition of Knutson–Lam–Speyer (2013, §3.3): each state, with the
landing time `0` removed and every remaining time decreased by one, is contained in the next
state — throwing at most one new ball each second. (KLS state this 1-indexed as
`J_{i+1} ⊇ (J_i \ {1}) - 1`.) -/
theorem image_sub_one_erase_zero_jugglingState_subset (N : GrassmannNecklace n d) (r : Fin n) :
    ((N.jugglingState r).erase 0).image (· - 1) ⊆ N.jugglingState (r + 1) := by
  have hinj : Function.Injective (· - r : Fin n → Fin n) := by
    intro a b h
    simpa [sub_add_cancel] using congrArg (· + r) h
  have h0 : ((N.toFun r).image (· - r)).erase 0 = ((N.toFun r).erase r).image (· - r) := by
    rw [Finset.image_erase hinj, sub_self]
  rw [jugglingState, jugglingState, h0, Finset.image_image]
  have hfun : ((· - (1 : Fin n)) ∘ (· - r)) = (· - (r + 1)) := by
    funext x
    show x - r - 1 = x - (r + 1)
    rw [sub_sub]
  rw [hfun]
  exact Finset.image_subset_image (N.erase_subset_toFun_succ r)

/-! ### The cyclic-interval necklace of the top cell -/

/-- A cyclic-interval entry has exactly `d` elements: it is the image of the initial segment
under the injective translation `· + i`. -/
private theorem card_uniformEntry (hdn : d ≤ n) (i : Fin n) :
    (univ.filter fun j : Fin n ↦ ((j - i : Fin n) : ℕ) < d).card = d := by
  have himg : (univ.filter fun j : Fin n ↦ ((j - i : Fin n) : ℕ) < d)
      = (univ.filter fun x : Fin n ↦ (x : ℕ) < d).image (· + i) := by
    ext j
    simp only [mem_filter, mem_univ, true_and, Finset.mem_image]
    constructor
    · intro hj
      exact ⟨j - i, hj, sub_add_cancel j i⟩
    · rintro ⟨x, hx, rfl⟩
      simpa [add_sub_cancel_right] using hx
  rw [himg, Finset.card_image_of_injective _ (add_left_injective i),
    Finset.card_filter_val_lt hdn]

/-- The cyclic-interval necklace `Iᵢ = {i, i+1, …, i+d-1}` (Postnikov, 2006, §17;
Oh–Postnikov–Speyer, 2015, §4): the Grassmann necklace of the uniform matroid `U_{d,n}`,
indexing the top cell of the totally nonnegative Grassmannian. Entry `i` collects the `j`
whose cyclic distance from `i` is less than `d`. -/
def uniform (hdn : d ≤ n) : GrassmannNecklace n d where
  toFun i := univ.filter fun j : Fin n ↦ ((j - i : Fin n) : ℕ) < d
  card_toFun i := card_uniformEntry hdn i
  toFun_succ_of_notMem i hi := by
    have hd0 : d = 0 := by
      have h : ¬ 0 < d := by
        intro h
        exact hi (Finset.mem_filter.mpr ⟨Finset.mem_univ i, by simpa [sub_self] using h⟩)
      omega
    subst hd0
    simp
  exists_toFun_succ_of_mem i hi := by
    have hd0 : 0 < d := by simpa [sub_self] using (Finset.mem_filter.mp hi).2
    rcases eq_or_lt_of_le hdn with hdn' | hdn'
    · refine ⟨i, ?_⟩
      have huniv : ∀ a : Fin n,
          (univ.filter fun j' : Fin n ↦ ((j' - a : Fin n) : ℕ) < d) = univ := by
        intro a
        apply Finset.filter_true_of_mem
        intro j' _
        rw [hdn']
        exact (j' - a).isLt
      rw [huniv, huniv, Finset.insert_erase (mem_univ i)]
    · refine ⟨i + ⟨d, hdn'⟩, ?_⟩
      ext j
      simp only [mem_filter, mem_univ, true_and, Finset.mem_insert, Finset.mem_erase, ne_eq]
      rw [show j - (i + 1) = j - i - 1 from by rw [sub_sub],
        show (j = i + ⟨d, hdn'⟩) ↔ j - i = ⟨d, hdn'⟩ from by rw [sub_eq_iff_eq_add, add_comm],
        show (j = i) ↔ j - i = 0 from by rw [sub_eq_zero]]
      by_cases hx0 : j - i = 0
      · rw [hx0]
        obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (NeZero.ne n)
        constructor
        · intro h
          simp at h
          omega
        · rintro (h | ⟨h, -⟩)
          · have hd' := congrArg Fin.val h
            simp at hd'
            omega
          · exact absurd rfl h
      · have hval : ((j - i - 1 : Fin n) : ℕ) = ((j - i : Fin n) : ℕ) - 1 :=
          Fin.val_sub_one_of_ne_zero hx0
        have hne : ((j - i : Fin n) : ℕ) ≠ 0 := by
          intro h
          exact hx0 (Fin.ext (by simp [h]))
        rw [hval]
        simp only [Fin.ext_iff, Fin.val_zero]
        omega

/-- Membership in the cyclic-interval necklace: `j` lies in entry `i` iff its cyclic distance
from `i` is less than `d`. -/
@[simp]
theorem mem_uniform_toFun {hdn : d ≤ n} {i j : Fin n} :
    j ∈ (uniform hdn).toFun i ↔ ((j - i : Fin n) : ℕ) < d := by
  simp [uniform]

/-- The general cyclic-interval necklace specializes to the `U_{2,4}` fixture of
`Theorems.Positroid`. -/
theorem uniform_eq_uniformTwoFour : uniform (by norm_num : 2 ≤ 4) = uniformTwoFour := by
  apply ext
  decide

/-- The Gale-shift comparison at `a` between a cyclic-interval entry and any `d`-element
subset: shifting back by `a` turns the entry into the initial segment, which is Gale-least. -/
private theorem galeShiftLE_uniformEntry (hdn : d ≤ n) (a : Fin n) {S : Finset (Fin n)}
    (hS : S.card = d) :
    Finset.GaleShiftLE a (univ.filter fun j : Fin n ↦ ((j - a : Fin n) : ℕ) < d) S := by
  have himg : ((univ.filter fun j : Fin n ↦ ((j - a : Fin n) : ℕ) < d).image (· - a))
      = univ.filter fun x : Fin n ↦ (x : ℕ) < d := by
    ext x
    simp only [Finset.mem_image, mem_filter, mem_univ, true_and]
    constructor
    · rintro ⟨j, hj, rfl⟩
      exact hj
    · intro hx
      exact ⟨x + a, by simpa [add_sub_cancel_right] using hx, add_sub_cancel_right x a⟩
  have hinj : Function.Injective (· - a : Fin n → Fin n) := by
    intro x y h
    simpa [sub_add_cancel] using congrArg (· + a) h
  show Finset.GaleLE _ _
  rw [himg]
  exact Finset.galeLE_filter_val_lt hdn ((Finset.card_image_of_injective _ hinj).trans hS)

/-- The cyclic-interval necklace is the least element of Lam's order: every necklace lies
above it entrywise, because the cyclic interval starting at `a` is the Gale-least `d`-subset
for the order shifted to `a`. Under the order-direction dictionary of the implementation
notes this is Postnikov's Lemma 17.6: the circular Bruhat order has a unique *top* element,
the cyclic-interval necklace of the top cell (Postnikov, 2006). -/
theorem isBot_uniform (hdn : d ≤ n) : IsBot (uniform (n := n) hdn) := by
  intro N a
  exact galeShiftLE_uniformEntry hdn a (N.card_toFun a)

/-- The basis family of the cyclic-interval necklace is the full family of `d`-subsets — the
bases of the uniform matroid `U_{d,n}`, whose cell is the whole totally nonnegative
Grassmannian ("Then `𝓜(𝓘) = ([n] choose k)`", Oh–Postnikov–Speyer, 2015, §4). -/
theorem ohFamily_uniform (hdn : d ≤ n) :
    (uniform (n := n) hdn).ohFamily = Finset.powersetCard d Finset.univ := by
  unfold ohFamily
  refine Finset.filter_eq_self.mpr ?_
  intro B hB j
  exact galeShiftLE_uniformEntry hdn j (Finset.mem_powersetCard_univ.mp hB)

/-! ### Fixtures from the literature (0-indexed) -/

/-- Postnikov's Figure 16.1 (Postnikov, 2006): the necklace of the decorated permutation
`π = (3,1,5,4,2,6)` with `4` a black and `6` a white fixed point — 1-indexed entries
`({1,2,6}, {2,3,6}, {1,3,6}, {1,5,6}, {1,5,6}, {1,2,6})`, stored here 0-indexed. -/
def postnikovFigure16_1 : GrassmannNecklace 6 3 where
  toFun := ![{0, 1, 5}, {1, 2, 5}, {0, 2, 5}, {0, 4, 5}, {0, 4, 5}, {0, 1, 5}]
  card_toFun := by decide
  toFun_succ_of_notMem := by decide
  exists_toFun_succ_of_mem := by decide

/-- In Postnikov's Figure 16.1, the white fixed point `6` (here `5`) lies in **every** entry
of the necklace: white fixed points are the coloops of the positroid (Postnikov, 2006,
§16–17). -/
theorem mem_toFun_postnikovFigure16_1 : ∀ i : Fin 6, (5 : Fin 6) ∈ postnikovFigure16_1.toFun i :=
  by decide

/-- In Postnikov's Figure 16.1, the black fixed point `4` (here `3`) lies in **no** entry of
the necklace: black fixed points are the loops of the positroid (Postnikov, 2006, §16–17). -/
theorem notMem_toFun_postnikovFigure16_1 :
    ∀ i : Fin 6, (3 : Fin 6) ∉ postnikovFigure16_1.toFun i := by
  decide

/-- Oh's worked example (Oh, 2011): the 1-indexed necklace
`({1,2,4}, {2,4,5}, {3,4,5}, {2,4,5}, {1,2,5})` on `[5]`, stored here 0-indexed. Its basis
family is the rank-3 positroid with the six bases `{124, 125, 134, 135, 245, 345}`; that
`ohFamily` computation is recorded as backlog, since the Gale order does not kernel-evaluate. -/
def ohWorkedExample : GrassmannNecklace 5 3 where
  toFun := ![{0, 1, 3}, {1, 3, 4}, {2, 3, 4}, {1, 3, 4}, {0, 1, 4}]
  card_toFun := by decide
  toFun_succ_of_notMem := by decide
  exists_toFun_succ_of_mem := by decide

/-- Oh's worked example has a repeated entry — the feature the page says it certifies. -/
example : ohWorkedExample.toFun 1 = ohWorkedExample.toFun 3 := by decide

/-- The page's non-example, 0-indexed: the 1-indexed sequence `({1,2}, {1,3}, {1,3})` on
`[3]`, which is *not* a Grassmann necklace. -/
def nonNecklace : Fin 3 → Finset (Fin 3) :=
  ![{0, 1}, {0, 2}, {0, 2}]

/-- The non-example violates the exchange condition (N1) at index `0`: `0 ∈ I₀`, but `I₁` is
not of the form `(I₀ \ {0}) ∪ {j}` for any `j`. -/
theorem nonNecklace_not_exists_insert_erase :
    ¬∃ j : Fin 3, nonNecklace (0 + 1) = insert j ((nonNecklace 0).erase 0) := by
  decide

/-- The non-example satisfies the fixed-entry condition (N2) at index `1`: `1 ∉ I₁` and the
next entry is unchanged. The violation of (N1) at index `0` is one-sided. -/
theorem nonNecklace_notMem_and_succ_eq :
    (1 : Fin 3) ∉ nonNecklace 1 ∧ nonNecklace (1 + 1) = nonNecklace 1 := by
  decide

/-- The non-example satisfies the exchange condition (N1) at index `2` (wrapping around to
entry `0`). Together with `nonNecklace_notMem_and_succ_eq`, the failure is exactly at
index `0`. -/
theorem nonNecklace_exists_insert_erase :
    (2 : Fin 3) ∈ nonNecklace 2 ∧
      ∃ j : Fin 3, nonNecklace (2 + 1) = insert j ((nonNecklace 2).erase 2) := by
  decide

/-- No Grassmann necklace of type `(2, 3)` has the non-example as its entry function: the
transition condition (N1) genuinely constrains beyond equal cardinality. -/
theorem toFun_ne_nonNecklace (N : GrassmannNecklace 3 2) : N.toFun ≠ nonNecklace := by
  intro h
  have h0 : (0 : Fin 3) ∈ N.toFun 0 := by
    rw [h]
    decide
  obtain ⟨j, hj⟩ := N.exists_toFun_succ_of_mem 0 h0
  rw [h] at hj
  exact nonNecklace_not_exists_insert_erase ⟨j, hj⟩

end GrassmannNecklace
