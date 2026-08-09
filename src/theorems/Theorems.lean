/-
Copyright (c) 2026 Harrison Totty. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Harrison Totty
-/
import Theorems.Basic
import Theorems.BoundaryMeasurement
import Theorems.GrassmannNecklace
import Theorems.LoewnerWhitney
import Theorems.Matroid
import Theorems.Positroid
import Theorems.TotallyPositiveKernel

/-!
# Research Harness — Theorems

Root module of the repository's Lean library. It imports every submodule so that
`lake build` (via `just lean-build`) elaborates the whole library from this
single target.

Add each new file under `Theorems/` and import it here to keep it in the build.

## Main definitions

None; this module only aggregates imports.

## Main statements

None; this module only aggregates imports.
-/
