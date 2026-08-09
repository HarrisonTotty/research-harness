/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Analysis.Calculus.LocalExtr.Rolle
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Gaussian.GaussianIntegral
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.Convex.Function
import Mathlib.LinearAlgebra.Matrix.Circulant
import Mathlib.LinearAlgebra.Matrix.Symmetric
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.Topology.Instances.Matrix
import Mathlib.Tactic
import Theorems.Positroid

/-!
# Totally positive kernels

A *totally positive kernel* is a function `K : X → Y → ℝ` on linearly ordered sets whose every
"minor" — the determinant `det [K (x i) (y j)]` built from strictly increasing arguments
`x 1 < ⋯ < x m`, `y 1 < ⋯ < y m` — is nonnegative (Karlin, *Total positivity, absorption
probabilities and applications*, Trans. AMS 111 (1964); *Total Positivity* Vol. I, 1968). Since
Mathlib's `Matrix X Y R` is literally the function type `X → Y → R` with arbitrary index types,
kernels and matrices are captured by a single predicate on `Matrix`, quantifying over strictly
monotone samplings — mirroring the statement shape of `Matrix.IsTotallyUnimodular`, with
`StrictMono` in place of `Function.Injective` because total positivity is not invariant under
row or column permutation.

## Naming convention

This module uses the modern matrix-school convention (Fomin–Zelevinsky;
Belton–Guillot–Khare–Putinar 2023, henceforth BGKP): `IsTotallyNonneg` means all minors are
nonnegative and `IsTotallyPos` means all minors are strictly positive. In the
Schoenberg/Karlin/Ando/Pinkus school these are called "totally positive (TP)" and "strictly
totally positive (STP)" respectively; Karlin's `TP_r` is `IsTotallyNonnegOfOrder · r` here.

## Main definitions

* `Matrix.IsTotallyNonneg` — every minor sampled at strictly increasing row and column
  arguments is nonnegative (Karlin's `TP_∞`).
* `Matrix.IsTotallyPos` — every such minor is strictly positive (Karlin's `STP_∞`).
* `Matrix.IsTotallyNonnegOfOrder`, `Matrix.IsTotallyPosOfOrder` — the order-`r` variants
  (Karlin's `TP_r`, `STP_r`), quantifying over minors of size at most `r`.
* `Matrix.IsSignRegularOfOrder` — minors of each size `k ≤ r` have the sign prescribed by
  `ε k` (Karlin's `SR_r`; Schoenberg's 1930 "minorendefinit").
* `Matrix.hankel` — the Hankel kernel `(i, j) ↦ v (i + j)`, the additive companion of
  `Matrix.circulant` (which, over a type with subtraction such as `ℝ` or `ℤ`, is literally the
  Toeplitz kernel `(i, j) ↦ v (i - j)`).
* `Matrix.IsOscillatory` — a totally nonnegative square matrix some positive power of which is
  totally positive (Gantmacher–Krein 1937).
* `Function.IsPolyaFrequency` — Pólya frequency functions in the modern formulation (BGKP
  Def. 8.1): nonnegative, Lebesgue integrable, nonzero at two or more points, with totally
  nonnegative Toeplitz kernel.

## Main statements

* `TotallyPositiveKernel.isTotallyPos_expMulKernel` — the exponential kernel
  `(x, y) ↦ exp (x * y)` on `ℝ × ℝ` is totally positive: the fundamental example of the
  theory (Johnson–Richards 2024; BGKP p. 8).
* `TotallyPositiveKernel.isTotallyPos_rpowKernel` — the generalized Vandermonde kernel
  `(x, y) ↦ x ^ y` on `(0, ∞) × ℝ` is totally positive (Gantmacher, *Theory of Matrices*,
  Ch. XIII §8).
* `Matrix.isTotallyPos_vandermonde` — a Vandermonde matrix with positive strictly increasing
  nodes is totally positive (all minors, refining `Matrix.det_vandermonde_pos`).
* `TotallyPositiveKernel.isTotallyPos_circulant_gaussian` — the Gaussian Toeplitz kernel
  `(x, y) ↦ exp (-(x - y)²)` is totally positive, by Schoenberg's 1947 reduction to the
  exponential kernel.
* `TotallyPositiveKernel.isPolyaFrequency_gaussian` and
  `TotallyPositiveKernel.isPolyaFrequency_oneSidedExp` — the Gaussian is a Pólya frequency
  function (Schoenberg 1947, eq. (4)), and so is the one-sided exponential
  `x ↦ if 0 ≤ x then exp (-x) else 0` (Schoenberg 1947, eq. (3)), the canonical
  discontinuous example.
* `TotallyPositiveKernel.isTotallyNonneg_indicatorKernel` and
  `TotallyPositiveKernel.det_submatrix_indicatorKernel_mem` — the indicator kernel
  `(x, y) ↦ if x ≤ y then 1 else 0` is totally nonnegative with every minor equal to `0` or
  `1`; `TotallyPositiveKernel.not_isTotallyPos_indicatorKernel` records that it is not
  totally positive.
* `Matrix.isTotallyNonnegOfOrder_two_circulant_iff` — for a positive function `Λ`, the
  Toeplitz kernel of `Λ` is totally nonnegative of order 2 iff `Λ` satisfies the
  multiplicative concavity inequality `Λ a' * Λ b' ≤ Λ a * Λ b` for nested same-sum pairs —
  the quantifier-explicit form of "`PF₂ ⇔ log-concave`" (Karlin 1968, Thms 4.1.8–4.1.9).
  `Matrix.isTotallyNonnegOfOrder_two_circulant_of_concaveOn_log` derives order-2 total
  nonnegativity from log-concavity proper.

## Implementation notes

* The predicates quantify over `k : ℕ` including `k = 0`, where the empty minor is `1`;
  Karlin quantifies over `1 ≤ m ≤ r` only. Over `ℝ` (indeed over any ordered ring) the two
  forms agree, and including `k = 0` matches `Matrix.IsTotallyUnimodular` and removes side
  conditions from the API.
* The definitions assume only `[Preorder X] [Preorder Y]` and `[CommRing R] [PartialOrder R]`
  — the weakest classes under which they are meaningful. Karlin's linearly ordered domains
  and real values are the intended specialization, used by all example theorems.
* Toeplitz kernels are deliberately *not* redefined: `Matrix.circulant` requires only
  `[Sub n]` on its index type, so `Matrix.circulant Λ : Matrix ℝ ℝ ℝ` is the Toeplitz kernel
  `(x, y) ↦ Λ (x - y)` and `Matrix.circulant a : Matrix ℤ ℤ ℝ` is the Toeplitz matrix of a
  bi-infinite sequence. `Matrix.hankel` fills the genuinely missing additive companion.
* Of the structural theory, only the order-2/log-concavity circle is formalized here. The
  variation-diminishing theorem, Fekete's criterion, Whitney density, the Gantmacher–Krein
  spectral theorem, the Schoenberg representation, the basic composition formula (blocked on
  the Cauchy–Binet formula, absent from Mathlib), and the Hankel/moment correspondence are
  recorded as backlog in the research notes.
* In `Matrix.isTotallyNonnegOfOrder_two_circulant_iff`, the interesting direction of
  "`PF₂ ⇔ log-concave`" that upgrades the pair inequality to `ConcaveOn` requires
  measurability regularization (Sierpiński's theorem on measurable midpoint-convex
  functions), which Mathlib does not yet have; the equivalence is therefore stated in the
  regularity-free pair form, with the `ConcaveOn` implication proved separately.

## References

* [S. Karlin, *Total positivity, absorption probabilities and applications*][karlin1964],
  Trans. AMS 111 (1964), 33–107. The `TP_r`/`SR_r` definitions (§0–§1).
* [I. J. Schoenberg, *On totally positive functions, Laplace transforms and total
  positivity*][schoenberg1947], Proc. Nat. Acad. Sci. 33 (1947), 11–17. The Pólya frequency
  axioms and the Gaussian/one-sided-exponential examples.
* [A. Belton, D. Guillot, A. Khare, M. Putinar, *Totally positive kernels, Pólya
  frequency functions, and their transforms*][bgkp2023], J. Analyse Math. 150 (2023).
  Modern kernel formulation (§1.1) and PF definition (Def. 8.1).
* [F. R. Gantmacher, M. G. Krein, *Sur les matrices complètement non négatives et
  oscillatoires*][gantmacherKrein1937], Compositio Math. 4 (1937), 445–476. Oscillation
  matrices.

## Tags

totally positive kernel, total positivity, totally nonnegative, Polya frequency function,
Toeplitz kernel, Hankel kernel, oscillation matrix, sign-regular, exponential kernel
-/

open Finset

namespace Matrix

section Defs

variable {X Y : Type*} [Preorder X] [Preorder Y] {R : Type*} [CommRing R] [PartialOrder R]

/-- A kernel (equivalently, a matrix with arbitrary ordered index types) is **totally
nonnegative** when every minor sampled at strictly increasing row and column arguments is
nonnegative. This is Karlin's "totally positive of order `∞`" (`TP_∞`, Karlin 1964, §0) in
the modern matrix-school naming (BGKP §1.1, "totally non-negative"); the statement shape
mirrors `Matrix.IsTotallyUnimodular` with `StrictMono` samplings. -/
def IsTotallyNonneg (A : Matrix X Y R) : Prop :=
  ∀ (k : ℕ) (f : Fin k → X) (g : Fin k → Y), StrictMono f → StrictMono g →
    0 ≤ (A.submatrix f g).det

/-- A kernel is **totally positive** when every minor sampled at strictly increasing row and
column arguments is strictly positive. This is Karlin's "strictly totally positive"
(`STP_∞`, Karlin 1964, §1) in the modern matrix-school naming (BGKP §1.1, "totally
positive"). -/
def IsTotallyPos (A : Matrix X Y R) : Prop :=
  ∀ (k : ℕ) (f : Fin k → X) (g : Fin k → Y), StrictMono f → StrictMono g →
    0 < (A.submatrix f g).det

/-- A kernel is **totally nonnegative of order `r`** when every minor of size at most `r`
sampled at strictly increasing arguments is nonnegative. This is Karlin's `TP_r` (Karlin
1964, §0, display (0.1)) in the modern naming (BGKP Def. 2.2, `TN_p`). -/
def IsTotallyNonnegOfOrder (A : Matrix X Y R) (r : ℕ) : Prop :=
  ∀ k ≤ r, ∀ (f : Fin k → X) (g : Fin k → Y), StrictMono f → StrictMono g →
    0 ≤ (A.submatrix f g).det

/-- A kernel is **totally positive of order `r`** when every minor of size at most `r`
sampled at strictly increasing arguments is strictly positive. This is Karlin's `STP_r`
(Karlin 1964, §1) in the modern naming (BGKP Def. 2.2, `TP_p`). -/
def IsTotallyPosOfOrder (A : Matrix X Y R) (r : ℕ) : Prop :=
  ∀ k ≤ r, ∀ (f : Fin k → X) (g : Fin k → Y), StrictMono f → StrictMono g →
    0 < (A.submatrix f g).det

/-- A kernel is **sign-regular of order `r`** with sign sequence `ε` when every minor of size
`k ≤ r` sampled at strictly increasing arguments has the sign prescribed by `ε k`, i.e.
`0 ≤ ε k * det`. With `ε` valued in `{±1}` this is Karlin's `SR_r` (Karlin 1964, §0;
Schoenberg 1930, "minorendefinit"); the single-size condition (Karlin's sign consistency
`SC_m`) is the case `r = m` restricted to `k = m`. Total nonnegativity of order `r` is the
case `ε ≡ 1`, see `Matrix.isSignRegularOfOrder_one_iff`. -/
def IsSignRegularOfOrder (A : Matrix X Y R) (ε : ℕ → R) (r : ℕ) : Prop :=
  ∀ k ≤ r, ∀ (f : Fin k → X) (g : Fin k → Y), StrictMono f → StrictMono g →
    0 ≤ ε k * (A.submatrix f g).det

/-- Sign regularity with the constant sign sequence `1` is exactly total nonnegativity of the
same order (Karlin 1964, §0: total positivity is sign regularity with every `ε_m = +1`). -/
theorem isSignRegularOfOrder_one_iff {A : Matrix X Y R} {r : ℕ} :
    A.IsSignRegularOfOrder 1 r ↔ A.IsTotallyNonnegOfOrder r := by
  simp only [IsSignRegularOfOrder, IsTotallyNonnegOfOrder, Pi.one_apply, one_mul]

/-- A totally positive kernel is totally nonnegative. -/
theorem IsTotallyPos.isTotallyNonneg {A : Matrix X Y R} (h : A.IsTotallyPos) :
    A.IsTotallyNonneg :=
  fun k f g hf hg => (h k f g hf hg).le

/-- A kernel totally positive of order `r` is totally nonnegative of order `r`. -/
theorem IsTotallyPosOfOrder.isTotallyNonnegOfOrder {A : Matrix X Y R} {r : ℕ}
    (h : A.IsTotallyPosOfOrder r) : A.IsTotallyNonnegOfOrder r :=
  fun k hk f g hf hg => (h k hk f g hf hg).le

/-- A totally nonnegative kernel is totally nonnegative of every order. -/
theorem IsTotallyNonneg.isTotallyNonnegOfOrder {A : Matrix X Y R} (h : A.IsTotallyNonneg)
    (r : ℕ) : A.IsTotallyNonnegOfOrder r :=
  fun k _ f g hf hg => h k f g hf hg

/-- A totally positive kernel is totally positive of every order. -/
theorem IsTotallyPos.isTotallyPosOfOrder {A : Matrix X Y R} (h : A.IsTotallyPos) (r : ℕ) :
    A.IsTotallyPosOfOrder r :=
  fun k _ f g hf hg => h k f g hf hg

/-- Total nonnegativity of order `r` restricts to any smaller order (Karlin 1964: the
defining inequalities for `TP_r` include those for `TP_{r'}`, `r' ≤ r`). -/
theorem IsTotallyNonnegOfOrder.of_le {A : Matrix X Y R} {r r' : ℕ}
    (h : A.IsTotallyNonnegOfOrder r) (hr : r' ≤ r) : A.IsTotallyNonnegOfOrder r' :=
  fun k hk f g hf hg => h k (hk.trans hr) f g hf hg

/-- Total nonnegativity is the conjunction of total nonnegativity of all finite orders
(Karlin 1964, §0: "if the subscript `∞` is written … the property in question will be
understood to hold for all values of `r`"). -/
theorem isTotallyNonneg_iff_forall_isTotallyNonnegOfOrder {A : Matrix X Y R} :
    A.IsTotallyNonneg ↔ ∀ r : ℕ, A.IsTotallyNonnegOfOrder r :=
  ⟨fun h r => h.isTotallyNonnegOfOrder r, fun h k f g hf hg => h k k le_rfl f g hf hg⟩

/-- Entries of a totally nonnegative kernel are nonnegative (the `1 × 1` minors). -/
theorem IsTotallyNonneg.entry_nonneg {A : Matrix X Y R} (h : A.IsTotallyNonneg) (x : X)
    (y : Y) : 0 ≤ A x y := by
  simpa [det_fin_one] using
    h 1 (fun _ => x) (fun _ => y) (Subsingleton.strictMono _) (Subsingleton.strictMono _)

/-- Entries of a totally positive kernel are positive (the `1 × 1` minors). -/
theorem IsTotallyPos.entry_pos {A : Matrix X Y R} (h : A.IsTotallyPos) (x : X) (y : Y) :
    0 < A x y := by
  simpa [det_fin_one] using
    h 1 (fun _ => x) (fun _ => y) (Subsingleton.strictMono _) (Subsingleton.strictMono _)

/-- Restricting a totally nonnegative kernel along strictly monotone maps of both variables
preserves total nonnegativity — Karlin's restriction principle: any restriction of a `TP_r`
kernel to `X' × Y'` with `X' ⊆ X`, `Y' ⊆ Y` is `TP_r`, "immediate from the definition's
quantifier over finite increasing tuples". -/
theorem IsTotallyNonneg.submatrix {X' Y' : Type*} [Preorder X'] [Preorder Y']
    {A : Matrix X Y R} (h : A.IsTotallyNonneg) {f : X' → X} {g : Y' → Y} (hf : StrictMono f)
    (hg : StrictMono g) : (A.submatrix f g).IsTotallyNonneg := by
  intro k r c hr hc
  rw [submatrix_submatrix]
  exact h k _ _ (hf.comp hr) (hg.comp hc)

/-- Restricting a totally positive kernel along strictly monotone maps of both variables
preserves total positivity. -/
theorem IsTotallyPos.submatrix {X' Y' : Type*} [Preorder X'] [Preorder Y']
    {A : Matrix X Y R} (h : A.IsTotallyPos) {f : X' → X} {g : Y' → Y} (hf : StrictMono f)
    (hg : StrictMono g) : (A.submatrix f g).IsTotallyPos := by
  intro k r c hr hc
  rw [submatrix_submatrix]
  exact h k _ _ (hf.comp hr) (hg.comp hc)

/-- The transpose of a totally nonnegative kernel is totally nonnegative: the definition is
symmetric in the two variables since `det` is invariant under transposition. -/
theorem IsTotallyNonneg.transpose {A : Matrix X Y R} (h : A.IsTotallyNonneg) :
    Aᵀ.IsTotallyNonneg := by
  intro k f g hf hg
  rw [← transpose_submatrix, det_transpose]
  exact h k g f hg hf

/-- The transpose of a totally positive kernel is totally positive. -/
theorem IsTotallyPos.transpose {A : Matrix X Y R} (h : A.IsTotallyPos) :
    Aᵀ.IsTotallyPos := by
  intro k f g hf hg
  rw [← transpose_submatrix, det_transpose]
  exact h k g f hg hf

/-- Maximal minors of a totally nonnegative matrix are nonnegative — the bridge to the
Plücker coordinates `Δ_I` of the positroid theory (`Matrix.maximalMinor`). -/
theorem IsTotallyNonneg.maximalMinor_nonneg {d n : ℕ} {A : Matrix (Fin d) (Fin n) R}
    (hA : A.IsTotallyNonneg) {I : Finset (Fin n)} (hI : I.card = d) :
    0 ≤ A.maximalMinor I hI :=
  hA d id (I.orderEmbOfFin hI) strictMono_id (I.orderEmbOfFin hI).strictMono

end Defs

section Hankel

variable {α β : Type*}

/-- The **Hankel kernel** of a function `v`: the kernel `(i, j) ↦ v (i + j)` (BGKP §1.2,
`H_f(x, y) = f(x + y)`). Over an index type with subtraction the multiplicative companion —
the Toeplitz kernel `(i, j) ↦ v (i - j)` — is `Matrix.circulant`. -/
def hankel [Add α] (v : α → β) : Matrix α α β :=
  of fun i j => v (i + j)

@[simp]
theorem hankel_apply [Add α] (v : α → β) (i j : α) : hankel v i j = v (i + j) :=
  rfl

/-- A Hankel kernel over a commutative additive index type is a symmetric matrix. -/
theorem isSymm_hankel [AddCommMagma α] (v : α → β) : (hankel v).IsSymm := by
  ext i j
  show hankel v j i = hankel v i j
  simp only [hankel_apply]
  exact congrArg v (add_comm j i)

end Hankel

section Oscillatory

variable {n : Type*} [LinearOrder n] [Fintype n] [DecidableEq n] {R : Type*} [CommRing R]
    [PartialOrder R]

/-- An **oscillation matrix** (matrice oscillatoire) is a totally nonnegative square matrix
some positive power of which is totally positive (Gantmacher–Krein 1937: "Une matrice A
complètement non négative est dite oscillatoire, si une certaine puissance `A^r` … est
complètement positive"). -/
def IsOscillatory (A : Matrix n n R) : Prop :=
  A.IsTotallyNonneg ∧ ∃ k : ℕ, 0 < k ∧ (A ^ k).IsTotallyPos

/-- A totally positive square matrix is an oscillation matrix (take the first power). -/
theorem IsTotallyPos.isOscillatory {A : Matrix n n R} (h : A.IsTotallyPos) :
    A.IsOscillatory :=
  ⟨h.isTotallyNonneg, 1, one_pos, by simpa [pow_one] using h⟩

end Oscillatory

end Matrix

namespace Function

/-- A **Pólya frequency function** in the modern formulation (BGKP Def. 8.1, after
Schoenberg 1947): a nonnegative Lebesgue-integrable function on `ℝ`, nonzero at two or more
points, whose Toeplitz kernel `(x, y) ↦ Λ (x - y)` (that is, `Matrix.circulant Λ`) is
totally nonnegative. Schoenberg's original axioms (α)(β)(γ) require measurability in place
of integrability and add the normalization split; the two formulations agree on the
integrable class Schoenberg names "Pólya frequency functions". -/
structure IsPolyaFrequency (Λ : ℝ → ℝ) : Prop where
  /-- A Pólya frequency function is nonnegative. -/
  nonneg : ∀ x : ℝ, 0 ≤ Λ x
  /-- A Pólya frequency function is Lebesgue integrable. -/
  integrable : MeasureTheory.Integrable Λ
  /-- A Pólya frequency function is nonzero at two or more points. -/
  exists_pair_ne_zero : ∃ x y : ℝ, x ≠ y ∧ Λ x ≠ 0 ∧ Λ y ≠ 0
  /-- The Toeplitz kernel of a Pólya frequency function is totally nonnegative. -/
  isTotallyNonneg : (Matrix.circulant Λ).IsTotallyNonneg

end Function

namespace Matrix

section LogConcavity

/-- A map out of `Fin 2` is strictly monotone as soon as it increases from `0` to `1`. Not private,
because `Theorems.LoewnerWhitney` builds `2 × 2` minors with it. -/
theorem strictMono_fin_two {α : Type*} [Preorder α] {f : Fin 2 → α} (h : f 0 < f 1) :
    StrictMono f := by
  intro i j hij
  fin_cases i <;> fin_cases j <;> simp_all

/-- The `2 × 2` minors of a Toeplitz kernel, expanded. -/
private theorem det_fin_two_circulant (Λ : ℝ → ℝ) (x y : Fin 2 → ℝ) :
    ((circulant Λ).submatrix x y).det
      = Λ (x 0 - y 0) * Λ (x 1 - y 1) - Λ (x 0 - y 1) * Λ (x 1 - y 0) := by
  rw [det_fin_two]
  rfl

/-- **`PF₂` as a multiplicative concavity inequality** (quantifier-explicit form of Karlin
1968, Thms 4.1.8–4.1.9: a translation kernel is `PF₂` iff the function is log-concave). For
a positive function `Λ`, the Toeplitz kernel `Matrix.circulant Λ` is totally nonnegative of
order 2 iff `Λ a' * Λ b' ≤ Λ a * Λ b` whenever the pair `(a, b)` is nested in the pair
`(a', b')` with the same sum. For measurable `Λ` this inequality is equivalent to
log-concavity of `Λ` (Sierpiński regularization — not yet formalizable in Mathlib); the
implication from log-concavity proper is
`Matrix.isTotallyNonnegOfOrder_two_circulant_of_concaveOn_log`. -/
theorem isTotallyNonnegOfOrder_two_circulant_iff {Λ : ℝ → ℝ} (hΛ : ∀ x : ℝ, 0 < Λ x) :
    (circulant Λ).IsTotallyNonnegOfOrder 2 ↔
      ∀ a b a' b' : ℝ, a' ≤ a → a ≤ b → b ≤ b' → a + b = a' + b' →
        Λ a' * Λ b' ≤ Λ a * Λ b := by
  constructor
  · intro h a b a' b' ha' hab hb' hsum
    rcases eq_or_lt_of_le ha' with rfl | ha'
    · obtain rfl : b = b' := by linarith
      exact le_rfl
    · have hdet := h 2 le_rfl ![0, b' - a] ![-a, -a']
        (strictMono_fin_two (by norm_num; linarith))
        (strictMono_fin_two (by norm_num; linarith))
      rw [det_fin_two_circulant] at hdet
      simp only [Matrix.cons_val_zero, Matrix.cons_val_one] at hdet
      have e00 : (0 : ℝ) - -a = a := by ring
      have e01 : (0 : ℝ) - -a' = a' := by ring
      have e10 : b' - a - -a = b' := by ring
      have e11 : b' - a - -a' = b := by linarith
      rw [e00, e01, e10, e11] at hdet
      linarith
  · rintro h (_ | _ | (_ | k)) hk f g hf hg
    · simp [det_fin_zero]
    · simpa [det_fin_one, submatrix_apply, circulant_apply] using (hΛ (f 0 - g 0)).le
    · rw [det_fin_two_circulant, sub_nonneg]
      have hf01 : f 0 < f 1 := hf (by decide)
      have hg01 : g 0 < g 1 := hg (by decide)
      rcases le_total (f 0 - g 0) (f 1 - g 1) with hpq | hpq
      · exact h (f 0 - g 0) (f 1 - g 1) (f 0 - g 1) (f 1 - g 0)
          (by linarith) hpq (by linarith) (by ring)
      · calc Λ (f 0 - g 1) * Λ (f 1 - g 0)
            ≤ Λ (f 1 - g 1) * Λ (f 0 - g 0) :=
              h (f 1 - g 1) (f 0 - g 0) (f 0 - g 1) (f 1 - g 0)
                (by linarith) hpq (by linarith) (by ring)
          _ = Λ (f 0 - g 0) * Λ (f 1 - g 1) := mul_comm _ _
    · exact absurd hk (by omega)

/-- **Log-concavity implies `PF₂`** (one direction of Karlin 1968, Thms 4.1.8–4.1.9): if a
positive function `Λ` is log-concave on `ℝ`, its Toeplitz kernel is totally nonnegative of
order 2. -/
theorem isTotallyNonnegOfOrder_two_circulant_of_concaveOn_log {Λ : ℝ → ℝ}
    (hΛ : ∀ x : ℝ, 0 < Λ x)
    (hlog : ConcaveOn ℝ Set.univ fun x : ℝ => Real.log (Λ x)) :
    (circulant Λ).IsTotallyNonnegOfOrder 2 := by
  rw [isTotallyNonnegOfOrder_two_circulant_iff hΛ]
  intro a b a' b' ha' hab hb' hsum
  rcases eq_or_lt_of_le ((ha'.trans hab).trans hb') with heq | hlt
  · obtain rfl : a' = a := by linarith
    obtain rfl : b = b' := by linarith
    exact le_rfl
  · have hd : (0 : ℝ) < b' - a' := by linarith
    set θ : ℝ := (b' - a) / (b' - a') with hθ
    have hθmul : θ * (b' - a') = b' - a := by
      rw [hθ]
      field_simp
    have hθ0 : 0 ≤ θ := div_nonneg (by linarith) hd.le
    have hθ1 : θ ≤ 1 := (div_le_one hd).2 (by linarith)
    have hθ1' : 0 ≤ 1 - θ := by linarith
    have hc1 : θ * a' + (1 - θ) * b' = a := by linear_combination -hθmul
    have hc2 : (1 - θ) * a' + θ * b' = b := by linear_combination hθmul - hsum
    have h1 : θ * Real.log (Λ a') + (1 - θ) * Real.log (Λ b') ≤ Real.log (Λ a) := by
      have := hlog.2 (Set.mem_univ a') (Set.mem_univ b') hθ0 hθ1' (by ring)
      simpa [smul_eq_mul, hc1] using this
    have h2 : (1 - θ) * Real.log (Λ a') + θ * Real.log (Λ b') ≤ Real.log (Λ b) := by
      have := hlog.2 (Set.mem_univ a') (Set.mem_univ b') hθ1' hθ0 (by ring)
      simpa [smul_eq_mul, hc2] using this
    have hlog_le :
        Real.log (Λ a') + Real.log (Λ b') ≤ Real.log (Λ a) + Real.log (Λ b) := by
      linarith
    have hexp := Real.exp_le_exp.2 hlog_le
    rwa [Real.exp_add, Real.exp_add, Real.exp_log (hΛ a'), Real.exp_log (hΛ b'),
      Real.exp_log (hΛ a), Real.exp_log (hΛ b)] at hexp

end LogConcavity

end Matrix

namespace TotallyPositiveKernel

open Matrix

/-! ### The exponential kernel

The fundamental example (Johnson–Richards 2024: "the fundamental example of a totally
positive kernel"; BGKP p. 8): `(x, y) ↦ exp (x * y)` is totally positive on `ℝ × ℝ`. The
proof of `isTotallyPos_expMulKernel` follows the classical argument: an exponential
polynomial `t ↦ ∑ c j * exp (a j * t)` with `k` distinct frequencies and a nonzero
coefficient vector has fewer than `k` zeros (induction on `k` via Rolle's theorem), so the
sampled determinants never vanish; a segment path to the Vandermonde evaluation point
`x i = i`, `y j = j` plus the intermediate value theorem pins the sign. -/

/-- The exponential kernel `(x, y) ↦ exp (x * y)` on `ℝ × ℝ` (Karlin 1968; Johnson–Richards
2024, "the fundamental example of a totally positive kernel"). -/
noncomputable def expMulKernel : Matrix ℝ ℝ ℝ :=
  of fun x y : ℝ => Real.exp (x * y)

@[simp]
theorem expMulKernel_apply (x y : ℝ) : expMulKernel x y = Real.exp (x * y) :=
  rfl

/-- An exponential polynomial `t ↦ ∑ j, c j * exp (a j * t)` with `n` distinct frequencies
`a j` that vanishes at `n` strictly increasing points has all coefficients zero. This is the
zero-counting core of the total positivity of the exponential kernel — the contrapositive of
"fewer than `n` zeros" — proved by induction on the number of terms via Rolle's theorem. -/
private theorem eq_zero_of_sum_mul_exp_eq_zero :
    ∀ (n : ℕ) (a c t : Fin n → ℝ), Function.Injective a → StrictMono t →
      (∀ i : Fin n, ∑ j : Fin n, c j * Real.exp (a j * t i) = 0) → c = 0 := by
  intro n
  induction n with
  | zero =>
    intro a c t _ _ _
    funext j
    exact j.elim0
  | succ n ih =>
    intro a c t ha ht hzero
    have hgderiv : ∀ u : ℝ,
        HasDerivAt (fun s : ℝ => ∑ j : Fin (n + 1), c j * Real.exp ((a j - a 0) * s))
          (∑ j : Fin (n + 1), c j * (a j - a 0) * Real.exp ((a j - a 0) * u)) u := by
      intro u
      have hterm : ∀ j : Fin (n + 1),
          HasDerivAt (fun s : ℝ => c j * Real.exp ((a j - a 0) * s))
            (c j * (a j - a 0) * Real.exp ((a j - a 0) * u)) u := by
        intro j
        have hid : HasDerivAt (fun s : ℝ => s) (1 : ℝ) u := hasDerivAt_id u
        have h1 : HasDerivAt (fun s : ℝ => (a j - a 0) * s) (a j - a 0) u := by
          simpa using hid.const_mul (a j - a 0)
        have h2 := (h1.exp).const_mul (c j)
        have heq : c j * (Real.exp ((a j - a 0) * u) * (a j - a 0))
            = c j * (a j - a 0) * Real.exp ((a j - a 0) * u) := by ring
        rwa [heq] at h2
      exact HasDerivAt.fun_sum fun j _ => hterm j
    have hgzero : ∀ i : Fin (n + 1),
        (∑ j : Fin (n + 1), c j * Real.exp ((a j - a 0) * t i)) = 0 := by
      intro i
      have hfac : ∀ j : Fin (n + 1),
          c j * Real.exp ((a j - a 0) * t i)
            = Real.exp (-(a 0 * t i)) * (c j * Real.exp (a j * t i)) := by
        intro j
        rw [show (a j - a 0) * t i = a j * t i + -(a 0 * t i) from by ring, Real.exp_add]
        ring
      calc (∑ j : Fin (n + 1), c j * Real.exp ((a j - a 0) * t i))
          = ∑ j : Fin (n + 1), Real.exp (-(a 0 * t i)) * (c j * Real.exp (a j * t i)) :=
            Finset.sum_congr rfl fun j _ => hfac j
        _ = Real.exp (-(a 0 * t i)) * ∑ j : Fin (n + 1), c j * Real.exp (a j * t i) :=
            (Finset.mul_sum _ _ _).symm
        _ = 0 := by rw [hzero i, mul_zero]
    have hcont : Continuous fun s : ℝ => ∑ j : Fin (n + 1), c j * Real.exp ((a j - a 0) * s) :=
      continuous_finsetSum _ fun j _ => by fun_prop
    have hrolle : ∀ i : Fin n, ∃ u ∈ Set.Ioo (t i.castSucc) (t i.succ),
        deriv (fun s : ℝ => ∑ j : Fin (n + 1), c j * Real.exp ((a j - a 0) * s)) u = 0 := by
      intro i
      refine exists_deriv_eq_zero (ht (Fin.castSucc_lt_succ (i := i))) hcont.continuousOn ?_
      rw [hgzero i.castSucc, hgzero i.succ]
    choose s hs hs' using hrolle
    have hsmono : StrictMono s := by
      intro i i' hii
      have h1 : s i < t i.succ := (hs i).2
      have h2 : t i'.castSucc < s i' := (hs i').1
      have h3 : t i.succ ≤ t i'.castSucc := by
        refine ht.monotone ?_
        simp only [Fin.le_def, Fin.val_succ, Fin.val_castSucc]
        omega
      linarith
    have hG : ∀ i : Fin n,
        ∑ j : Fin n, c j.succ * (a j.succ - a 0) * Real.exp ((a j.succ - a 0) * s i) = 0 := by
      intro i
      have hderiv := (hgderiv (s i)).deriv
      have hzero' :
          ∑ j : Fin (n + 1), c j * (a j - a 0) * Real.exp ((a j - a 0) * s i) = 0 := by
        rw [← hderiv]
        exact hs' i
      rw [Fin.sum_univ_succ] at hzero'
      simpa using hzero'
    have hc' : (fun j : Fin n => c j.succ * (a j.succ - a 0)) = 0 := by
      refine ih (fun j => a j.succ - a 0) (fun j => c j.succ * (a j.succ - a 0)) s ?_ hsmono hG
      intro j₁ j₂ hj
      exact Fin.succ_injective n (ha (sub_left_inj.mp hj))
    have hcsucc : ∀ j : Fin n, c j.succ = 0 := by
      intro j
      have h0 : c j.succ * (a j.succ - a 0) = 0 := congrFun hc' j
      have hne : a j.succ - a 0 ≠ 0 := by
        intro hcontra
        exact Fin.succ_ne_zero j (ha (sub_eq_zero.mp hcontra))
      exact (mul_eq_zero.mp h0).resolve_right hne
    have hc0 : c 0 = 0 := by
      have h0 := hzero 0
      rw [Fin.sum_univ_succ] at h0
      have hrest : ∑ j : Fin n, c j.succ * Real.exp (a j.succ * t 0) = 0 :=
        Finset.sum_eq_zero fun j _ => by rw [hcsucc j, zero_mul]
      rw [hrest, add_zero] at h0
      exact (mul_eq_zero.mp h0).resolve_right (Real.exp_ne_zero _)
    funext j
    show c j = 0
    exact Fin.cases hc0 hcsucc j

/-- Minors of the exponential kernel never vanish: the columns `t ↦ exp (y j * t)` are
linearly independent over distinct frequencies. -/
private theorem det_submatrix_expMulKernel_ne_zero {k : ℕ} {x y : Fin k → ℝ}
    (hx : StrictMono x) (hy : Function.Injective y) :
    (expMulKernel.submatrix x y).det ≠ 0 := by
  intro hdet
  obtain ⟨v, hv, hmul⟩ := Matrix.exists_mulVec_eq_zero_iff.2 hdet
  refine hv (eq_zero_of_sum_mul_exp_eq_zero k y v x hy hx fun i => ?_)
  have hi : ∑ j : Fin k, (expMulKernel.submatrix x y) i j * v j = 0 := congrFun hmul i
  rw [← hi]
  exact Finset.sum_congr rfl fun j _ => by
    simp only [Matrix.submatrix_apply, expMulKernel_apply]
    rw [mul_comm (v j), mul_comm (x i)]

/-- A convex combination of a strictly monotone tuple with the standard increasing tuple is
strictly monotone. This is the deformation path along which the sign of an exponential-kernel
minor is transported to the Vandermonde evaluation point. -/
private theorem strictMono_interp {k : ℕ} {x : Fin k → ℝ} (hx : StrictMono x) {τ : ℝ}
    (hτ0 : 0 ≤ τ) (hτ1 : τ ≤ 1) :
    StrictMono fun i : Fin k => (1 - τ) * x i + τ * ((i : ℕ) : ℝ) := by
  intro i j hij
  have h1 : x i < x j := hx hij
  have h2 : ((i : ℕ) : ℝ) < ((j : ℕ) : ℝ) := by exact_mod_cast Fin.lt_def.mp hij
  rcases eq_or_lt_of_le hτ0 with rfl | hτpos
  · simpa using h1
  · have h3 : 0 < τ * (((j : ℕ) : ℝ) - ((i : ℕ) : ℝ)) := mul_pos hτpos (by linarith)
    have h4 : 0 ≤ (1 - τ) * (x j - x i) := mul_nonneg (by linarith) (by linarith)
    nlinarith

/-- **The exponential kernel is totally positive** (Johnson–Richards 2024; BGKP p. 8:
`e^{xy}` is TP on any `I' ⊆ ℝ` in their strict sense — `STP_∞` in Karlin's terms). The page
certifies it as a totally positive kernel on the full plane and the target of the Gaussian's
reduction identity. -/
theorem isTotallyPos_expMulKernel : expMulKernel.IsTotallyPos := by
  intro k x y hx hy
  set X : ℝ → Fin k → ℝ := fun τ i => (1 - τ) * x i + τ * ((i : ℕ) : ℝ) with hX
  set Y : ℝ → Fin k → ℝ := fun τ j => (1 - τ) * y j + τ * ((j : ℕ) : ℝ) with hY
  set h : ℝ → ℝ := fun τ => (expMulKernel.submatrix (X τ) (Y τ)).det with hh
  have hne : ∀ τ ∈ Set.Icc (0 : ℝ) 1, h τ ≠ 0 := by
    intro τ hτ
    exact det_submatrix_expMulKernel_ne_zero (strictMono_interp hx hτ.1 hτ.2)
      (strictMono_interp hy hτ.1 hτ.2).injective
  have hcont : Continuous h := by
    refine Continuous.matrix_det (continuous_matrix fun i j => ?_)
    simp only [hX, hY, Matrix.submatrix_apply, expMulKernel_apply]
    fun_prop
  have h1pos : 0 < h 1 := by
    have hmat : expMulKernel.submatrix (X 1) (Y 1)
        = Matrix.vandermonde fun i : Fin k => Real.exp ((i : ℕ) : ℝ) := by
      ext i j
      simp only [hX, hY, Matrix.submatrix_apply, expMulKernel_apply, Matrix.vandermonde_apply]
      rw [← Real.rpow_natCast (Real.exp ((i : ℕ) : ℝ)) (j : ℕ),
        Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
      norm_num
    rw [hh]
    simp only
    rw [hmat]
    exact Matrix.det_vandermonde_pos fun i j hij => by
      exact Real.exp_lt_exp.2 (by exact_mod_cast Fin.lt_def.mp hij)
  have h0 : h 0 = (expMulKernel.submatrix x y).det := by
    simp [hh, hX, hY]
  rw [← h0]
  rcases lt_trichotomy (h 0) 0 with hlt | heq | hgt
  · obtain ⟨τ, hτmem, hτ0⟩ :=
      intermediate_value_Icc (by norm_num : (0 : ℝ) ≤ 1) hcont.continuousOn
        (Set.mem_Icc.2 ⟨hlt.le, h1pos.le⟩)
    exact absurd hτ0 (hne τ hτmem)
  · exact absurd heq (hne 0 (by norm_num))
  · exact hgt

/-! ### The generalized Vandermonde kernel -/

/-- The generalized Vandermonde kernel `(x, y) ↦ x ^ y` on `(0, ∞) × ℝ` (Gantmacher,
*Theory of Matrices*, Ch. XIII §8, via BGKP Ex. 2.3). -/
noncomputable def rpowKernel : Matrix {x : ℝ // 0 < x} ℝ ℝ :=
  of fun x y => (x : ℝ) ^ y

/-- **The generalized Vandermonde kernel is totally positive**: the matrices `(u i ^ α j)`
with `0 < u 1 < ⋯ < u n` and `α 1 < ⋯ < α n` have all minors strictly positive (Gantmacher,
*Theory of Matrices*, Ch. XIII §8, via BGKP Ex. 2.3). Via `x ^ y = exp (log x * y)` this is
the exponential kernel sampled at strictly increasing arguments. The page certifies it as
total positivity without symmetry. -/
theorem isTotallyPos_rpowKernel : rpowKernel.IsTotallyPos := by
  intro k f g hf hg
  have heq : rpowKernel.submatrix f g
      = expMulKernel.submatrix (fun i => Real.log (f i : ℝ)) g := by
    ext i j
    simp only [Matrix.submatrix_apply, rpowKernel, expMulKernel_apply, Matrix.of_apply]
    exact Real.rpow_def_of_pos (f i).2 (g j)
  rw [heq]
  refine isTotallyPos_expMulKernel k _ g (fun i j hij => ?_) hg
  exact Real.log_lt_log (f i).2 (hf hij)

/-! ### The Gaussian -/

/-- The Gaussian `x ↦ exp (-x²)`, the flagship Pólya frequency function (Schoenberg 1947,
eq. (4)). -/
noncomputable def gaussian : ℝ → ℝ :=
  fun x : ℝ => Real.exp (-x ^ 2)

@[simp]
theorem gaussian_apply (x : ℝ) : gaussian x = Real.exp (-x ^ 2) :=
  rfl

/-- **The Gaussian Toeplitz kernel `(x, y) ↦ exp (-(x - y)²)` is totally positive**, via
Schoenberg's reduction (Schoenberg 1947, eq. (4)):
`det [exp (-(x i - y j)²)] = exp (-∑ x i²) * exp (-∑ y j²) * det [exp (2 * x i * y j)]`,
inheriting positivity from the exponential kernel. -/
theorem isTotallyPos_circulant_gaussian : (circulant gaussian).IsTotallyPos := by
  intro k x y hx hy
  have hsplit : (circulant gaussian).submatrix x y
      = of fun i j : Fin k => Real.exp (-(x i) ^ 2) *
          (of fun i' j' : Fin k => Real.exp (-(y j') ^ 2) *
            (expMulKernel.submatrix (fun i'' => 2 * x i'') y) i' j') i j := by
    ext i j
    simp only [Matrix.submatrix_apply, circulant_apply, Matrix.of_apply, gaussian_apply,
      expMulKernel_apply]
    rw [← Real.exp_add, ← Real.exp_add]
    congr 1
    ring
  rw [hsplit, det_mul_column, det_mul_row]
  have hB : 0 < (expMulKernel.submatrix (fun i => 2 * x i) y).det :=
    isTotallyPos_expMulKernel k _ y (fun i j hij => by linarith [hx hij]) hy
  have hu : 0 < ∏ i : Fin k, Real.exp (-(x i) ^ 2) :=
    Finset.prod_pos fun i _ => Real.exp_pos _
  have hw : 0 < ∏ j : Fin k, Real.exp (-(y j) ^ 2) :=
    Finset.prod_pos fun j _ => Real.exp_pos _
  exact mul_pos hu (mul_pos hw hB)

/-- **The Gaussian is a Pólya frequency function** (Schoenberg 1947, eq. (4); BGKP p. 38). -/
theorem isPolyaFrequency_gaussian : Function.IsPolyaFrequency gaussian where
  nonneg x := (Real.exp_pos _).le
  integrable := by
    show MeasureTheory.Integrable fun x : ℝ => Real.exp (-x ^ 2)
    simpa using integrable_exp_neg_mul_sq (b := 1) one_pos
  exists_pair_ne_zero := ⟨0, 1, by norm_num, Real.exp_ne_zero _, Real.exp_ne_zero _⟩
  isTotallyNonneg := isTotallyPos_circulant_gaussian.isTotallyNonneg

/-! ### The indicator kernel -/

/-- The indicator kernel `(x, y) ↦ if x ≤ y then 1 else 0` on a linear order — the
degenerate-rank totally nonnegative example, sharing the two-sided product shape of
Sturm–Liouville Green's functions (Karlin 1964, §8). -/
def indicatorKernel (α : Type*) [LinearOrder α] : Matrix α α ℝ :=
  of fun x y : α => if x ≤ y then 1 else 0

@[simp]
theorem indicatorKernel_apply {α : Type*} [LinearOrder α] (x y : α) :
    indicatorKernel α x y = if x ≤ y then 1 else 0 :=
  rfl

/-- **Every minor of the indicator kernel is `0` or `1`** (page: "TP with minors in
`{0, 1}`", locally verified in exact integer arithmetic). The proof subtracts consecutive
rows, turning the sampled matrix into one whose rows are indicators of pairwise disjoint,
ordered intervals; only the identity permutation survives in the Leibniz expansion. -/
theorem det_submatrix_indicatorKernel_mem {α : Type*} [LinearOrder α] {k : ℕ}
    {x y : Fin k → α} (hx : StrictMono x) (hy : StrictMono y) :
    ((indicatorKernel α).submatrix x y).det ∈ ({0, 1} : Set ℝ) := by
  classical
  have key : ((indicatorKernel α).submatrix x y).det = 0
      ∨ ((indicatorKernel α).submatrix x y).det = 1 := by
    by_cases hB : ∃ i j : Fin k, j < i ∧ x i ≤ y j
    · left
      set S : Finset (Fin k) :=
        Finset.univ.filter (fun j : Fin k => ∃ i : Fin k, j < i ∧ x i ≤ y j) with hS
      have hSne : S.Nonempty := by
        obtain ⟨i, j, hji, hxy⟩ := hB
        refine ⟨j, ?_⟩
        simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and]
        exact ⟨i, hji, hxy⟩
      have hjS : S.min' hSne ∈ S := S.min'_mem hSne
      set j : Fin k := S.min' hSne with hjdef
      obtain ⟨i, hji, hxyj⟩ : ∃ i : Fin k, j < i ∧ x i ≤ y j := by
        simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and] at hjS
        exact hjS
      have hmin : ∀ m : Fin k, m < j → ∀ i' : Fin k, m < i' → ¬ x i' ≤ y m := by
        intro m hm i' hi' hle
        have hmS : m ∈ S := by
          simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and]
          exact ⟨i', hi', hle⟩
        exact absurd (S.min'_le m hmS) (not_le.mpr hm)
      have hjval : (j : ℕ) < (i : ℕ) := hji
      have hik : (i : ℕ) < k := i.isLt
      have hjk : (j : ℕ) + 1 < k := by omega
      set j' : Fin k := ⟨(j : ℕ) + 1, hjk⟩ with hj'def
      have hjj' : j < j' := by
        simp only [hj'def, Fin.lt_def]
        omega
      have hj'i : j' ≤ i := by
        simp only [hj'def, Fin.le_def]
        omega
      have hne : j ≠ j' := ne_of_lt hjj'
      refine det_zero_of_row_eq hne ?_
      funext m
      simp only [Matrix.submatrix_apply, indicatorKernel_apply]
      by_cases hm : j ≤ m
      · have hxj' : x j' ≤ y m :=
          le_trans (hx.monotone hj'i) (le_trans hxyj (hy.monotone hm))
        rw [if_pos (le_trans (hx.monotone hjj'.le) hxj'), if_pos hxj']
      · replace hm : m < j := not_le.mp hm
        rw [if_neg (hmin m hm j hm), if_neg (hmin m hm j' (hm.trans hjj'))]
    · have hB' : ∀ i j : Fin k, j < i → ¬ x i ≤ y j := fun i j hij hle =>
        hB ⟨i, j, hij, hle⟩
      have htri : ((indicatorKernel α).submatrix x y).BlockTriangular id := by
        intro i j hij
        simp only [Matrix.submatrix_apply, indicatorKernel_apply]
        exact if_neg (hB' i j hij)
      rw [det_of_upperTriangular htri]
      refine Finset.prod_induction _ (fun r : ℝ => r = 0 ∨ r = 1) ?_ ?_ ?_
      · rintro p q (rfl | rfl) (rfl | rfl) <;> simp
      · exact Or.inr rfl
      · intro i _
        simp only [Matrix.submatrix_apply, indicatorKernel_apply]
        by_cases hc : x i ≤ y i
        · exact Or.inr (if_pos hc)
        · exact Or.inl (if_neg hc)
  rcases key with h | h <;> simp [h]

/-- **The indicator kernel is totally nonnegative** (Karlin 1964, §8: the Green's-function
shape; minors are `0` or `1`). -/
theorem isTotallyNonneg_indicatorKernel (α : Type*) [LinearOrder α] :
    (indicatorKernel α).IsTotallyNonneg := by
  intro k x y hx hy
  have h := det_submatrix_indicatorKernel_mem hx hy
  simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at h
  rcases h with h | h
  · exact h.ge
  · rw [h]
    exact zero_le_one

/-- **The indicator kernel is not totally positive**: on `ℝ` the `2 × 2` minor at
`x = (0, 1)`, `y = (2, 3)` is `0`. Together with `isTotallyNonneg_indicatorKernel` this
certifies "totally positive ≠ strictly totally positive" (in the page's Karlin terms). -/
theorem not_isTotallyPos_indicatorKernel : ¬ (indicatorKernel ℝ).IsTotallyPos := by
  intro h
  have hpos := h 2 ![0, 1] ![2, 3] (strictMono_fin_two (by norm_num))
    (strictMono_fin_two (by norm_num))
  rw [det_fin_two] at hpos
  norm_num [Matrix.submatrix_apply, indicatorKernel_apply] at hpos

/-! ### The one-sided exponential -/

/-- The one-sided exponential `x ↦ if 0 ≤ x then exp (-x) else 0`, Schoenberg's canonical
discontinuous Pólya frequency function (Schoenberg 1947, eq. (3)). -/
noncomputable def oneSidedExp : ℝ → ℝ :=
  fun x : ℝ => if 0 ≤ x then Real.exp (-x) else 0

/-- **The one-sided exponential is a Pólya frequency function** (Schoenberg 1947, eq. (3)).
Total nonnegativity of its Toeplitz kernel reduces to the indicator kernel by factoring
`exp (-(x - y)) = exp (-x) * exp y` out of rows and columns. Schoenberg proves its affine
images are the *only* discontinuous totally positive functions; the discontinuity at `0`
makes it the boundary case of that classification. -/
theorem isPolyaFrequency_oneSidedExp : Function.IsPolyaFrequency oneSidedExp where
  nonneg x := by
    simp only [oneSidedExp]
    split
    · exact (Real.exp_pos _).le
    · exact le_rfl
  integrable := by
    have hset : oneSidedExp = Set.indicator (Set.Ici 0) fun x : ℝ => Real.exp (-x) := by
      funext x
      simp only [oneSidedExp, Set.indicator_apply, Set.mem_Ici]
    rw [hset]
    exact (MeasureTheory.integrable_indicator_iff measurableSet_Ici).2
      (Iff.mpr integrableOn_Ici_iff_integrableOn_Ioi (integrableOn_exp_neg_Ioi 0))
  exists_pair_ne_zero := by
    refine ⟨0, 1, by norm_num, ?_, ?_⟩
    · simp [oneSidedExp]
    · simp [oneSidedExp]
  isTotallyNonneg := by
    intro k x y hx hy
    have hsplit : (circulant oneSidedExp).submatrix x y
        = of fun i j : Fin k => Real.exp (-(x i)) *
            ((of fun i' j' : Fin k =>
              Real.exp (y j') * (((indicatorKernel ℝ).submatrix y x)ᵀ i' j')) i j) := by
      ext i j
      simp only [Matrix.submatrix_apply, circulant_apply, Matrix.of_apply,
        Matrix.transpose_apply, oneSidedExp, indicatorKernel_apply]
      by_cases hc : y j ≤ x i
      · rw [if_pos (by linarith : (0 : ℝ) ≤ x i - y j), if_pos hc, mul_one, ← Real.exp_add]
        congr 1
        ring
      · rw [if_neg (by linarith : ¬ (0 : ℝ) ≤ x i - y j), if_neg hc]
        ring
    rw [hsplit, det_mul_column, det_mul_row, det_transpose]
    have hN : 0 ≤ ((indicatorKernel ℝ).submatrix y x).det :=
      isTotallyNonneg_indicatorKernel ℝ k y x hy hx
    have h1 : 0 < ∏ i : Fin k, Real.exp (-(x i)) :=
      Finset.prod_pos fun i _ => Real.exp_pos _
    have h2 : 0 < ∏ j : Fin k, Real.exp (y j) :=
      Finset.prod_pos fun j _ => Real.exp_pos _
    exact mul_nonneg h1.le (mul_nonneg h2.le hN)

/-! ### The trivial totally nonnegative kernel -/

/-- Minors of size at least 2 of the Toeplitz kernel of `t ↦ exp (a * t + b)` vanish: the
kernel has rank one (Schoenberg 1947, eq. (2): `e^{ax+b}` is trivially totally positive).
Stated for arbitrary samplings — monotonicity is not needed for the vanishing. -/
theorem det_submatrix_circulant_exp_affine {a b : ℝ} {k : ℕ} (hk : 2 ≤ k)
    (x y : Fin k → ℝ) :
    ((circulant fun t : ℝ => Real.exp (a * t + b)).submatrix x y).det = 0 := by
  have hsplit : (circulant fun t : ℝ => Real.exp (a * t + b)).submatrix x y
      = of fun i j : Fin k =>
          Real.exp (a * x i + b) * (of fun _ j' : Fin k => Real.exp (-(a * y j'))) i j := by
    ext i j
    simp only [submatrix_apply, circulant_apply, of_apply]
    rw [← Real.exp_add]
    exact congrArg Real.exp (by ring)
  rw [hsplit, det_mul_column]
  have h01 : (⟨0, by omega⟩ : Fin k) ≠ ⟨1, by omega⟩ := by
    simp [Fin.ext_iff]
  rw [det_zero_of_row_eq h01 rfl, mul_zero]

/-- **The Toeplitz kernel of `t ↦ exp (a * t + b)` is totally nonnegative**, with all minors
of size at least 2 equal to zero (Schoenberg 1947, eq. (2)): the degenerate stratum of the
definition. It is a totally positive function in Schoenberg's axioms yet, being monotone
(hence not integrable), it is not a Pólya frequency function — `Function.IsPolyaFrequency`
rejects it through the integrability field. -/
theorem isTotallyNonneg_circulant_exp_affine (a b : ℝ) :
    (circulant fun t : ℝ => Real.exp (a * t + b)).IsTotallyNonneg := by
  rintro (_ | _ | k) f g hf hg
  · simp [det_fin_zero]
  · simpa [det_fin_one, circulant_apply] using (Real.exp_pos (a * (f 0 - g 0) + b)).le
  · rw [det_submatrix_circulant_exp_affine (by omega)]

end TotallyPositiveKernel

namespace Matrix

/-- **Vandermonde matrices with positive strictly increasing nodes are totally positive**:
every minor (not just the full determinant of `Matrix.det_vandermonde_pos`) is strictly
positive, since `v i ^ (j : ℕ)` is the generalized Vandermonde kernel sampled at natural
exponents (Gantmacher, *Theory of Matrices*, Ch. XIII §8: "the classical Vandermonde as the
integer-exponent specialization"). -/
theorem isTotallyPos_vandermonde {n : ℕ} {v : Fin n → ℝ} (hv : StrictMono v)
    (h0 : ∀ i : Fin n, 0 < v i) : (vandermonde v).IsTotallyPos := by
  have heq : vandermonde v
      = TotallyPositiveKernel.rpowKernel.submatrix
          (fun i => (⟨v i, h0 i⟩ : {x : ℝ // 0 < x})) fun j : Fin n => ((j : ℕ) : ℝ) := by
    ext i j
    simp only [vandermonde_apply, Matrix.submatrix_apply, TotallyPositiveKernel.rpowKernel,
      Matrix.of_apply]
    exact (Real.rpow_natCast (v i) (j : ℕ)).symm
  rw [heq]
  refine TotallyPositiveKernel.isTotallyPos_rpowKernel.submatrix (fun i j hij => ?_)
    fun i j hij => ?_
  · exact hv hij
  · exact_mod_cast Fin.lt_def.mp hij

end Matrix
