/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Mathlib.Tactic

/-!
# Scaffolding placeholder

A minimal, self-contained example that exercises the toolchain and demonstrates
the file conventions (copyright header, module docstring, documented
declarations, a warning-free proof). Replace it with real material as the
library grows.

## Main definitions

* `ResearchHarness.double` — the sum of a natural number with itself.

## Main statements

* `ResearchHarness.double_eq_two_mul` — `double n` equals `2 * n`.
-/

namespace ResearchHarness

/-- `double n` is `n` added to itself. -/
def double (n : ℕ) : ℕ := n + n

/-- Doubling a natural number is the same as multiplying it by two. -/
theorem double_eq_two_mul (n : ℕ) : double n = 2 * n := by
  simp only [double, two_mul]

end ResearchHarness
