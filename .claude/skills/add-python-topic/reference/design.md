# Design contract for topic implementations

## Page-section → code-artifact mapping

The Logseq page template and this contract are two halves of one pipeline.
`spec.md` must contain this table, filled in for the topic:

| Page section                            | Code artifact                                        |
| --------------------------------------- | ---------------------------------------------------- |
| Definition — primary formulation        | Internal representation + boundary validation        |
| Definition — equivalent axiomatizations | `from_<formulation>` classmethod constructors        |
| Abstracts / generalizes (intro links)   | Transformation targets — or backlog entries          |
| Derived vocabulary                      | Computed properties / methods                        |
| Operations and constructions            | Transformation methods (`to_*`, `from_*`, `dual`, …) |
| Structural theorems                     | Hypothesis property tests (see testing.md)           |
| Canonical examples                      | Named example constructors + test fixtures           |
| Open questions                          | Candidate experiments (report only, do not build)    |

Section names follow the page template of `add-logseq-topic`; the equivalent
axiomatizations are sibling blocks under **Definition**, not a section of
their own.

## Choosing the internal representation

When the page gives several equivalent axiomatizations, pick **one** as the
stored representation and derive the rest as views. Prefer the formulation
that makes boundary validation cheapest and the most derived vocabulary
directly computable; say in `spec.md` which one won and why. Every other
standard formulation still gets a `from_<formulation>` constructor (which
converts and validates) and, where useful, a corresponding accessor — that
is the code-side cryptomorphism table.

## API shape

- One module per topic: `src/research/<topic>.py`; split into a package only
  when a single module genuinely becomes a grab bag.
- The object is a frozen dataclass over immutable collections (`frozenset`,
  `tuple`). Mathematical objects are values, not identities.
- `__init__` stays cheap; axiom checking lives in the classmethod
  constructors, the only supported way to build the object. A dataclass
  `__init__` cannot actually be hidden, so say in the class docstring that
  calling it directly skips validation. Validation
  errors state which numbered axiom (as the page numbers them) failed and on
  what input.
- Expensive derived invariants may use `functools.cached_property` (works on
  frozen, non-slotted dataclasses); document the cost in the docstring.
- Implement `__repr__` compactly enough to read in a failing test.
- Every formula or algorithm docstring cites the source the page attributes
  it to (author, year, theorem/section number) — the page already carries
  the attribution; carry it into the code.

## The four required capabilities

**1. Computed properties.** Each derived-vocabulary term and each numeric
invariant on the page becomes a method or property with the same name the
page uses (snake-cased). Preserve claim strength: a quantity the page defines
only for a special case raises `ValueError` outside it rather than guessing.

**2. Transformations.** `to_<other>()` / `from_<other>()` methods for each
operation and each related structure on the page. Implement only
transformations whose target exists — a stdlib/collections type, an existing
`src/research` module, or the object's own type (duals, minors,
restrictions). Targets whose structure has no implementation yet go in the
spec's **transformation backlog**, not into `NotImplementedError` stubs; the
backlog mirrors the graph's red links and is reported in Step 8.

**3. DataFrame serialization.** `to_dataframe()` returns a canonical tidy
encoding (one row per element/relation; document the columns and orientation
in the docstring) and `from_dataframe()` inverts it exactly. Keep the
encoding compatible with `experiments.io.write_result` (`records`-oriented
JSON) so experiment results over these objects serialize for free.

**4. Visualization.** At least one graphical form: `plot_<form>()` methods
that draw onto a provided `matplotlib.axes.Axes` (creating one only when the
caller passes none) and return it — never call `show()` or write files. Add
`matplotlib` via `uv add matplotlib` if the project lacks it. Choose forms
the page itself suggests (a lattice of flats, a bipartite incidence graph, a
geometric representation); a text/ASCII rendering is a welcome extra, not a
substitute. Keep plotting at the edge: core methods must not import or
depend on matplotlib.
