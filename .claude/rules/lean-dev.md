---
paths:
  - "**/*.lean"
---

# Lean Development

Target is **Lean 4**, current stable series, with **Mathlib** conventions. The
repository's toolchain file is the authoritative version — never work around a
mismatch by rewriting code to suit a different toolchain.

Follow Mathlib's conventions even in modules that do not depend on Mathlib, so
material can migrate in either direction without a rewrite.

## Naming

- Theorems, lemmas, and any term whose type is a `Prop` use `snake_case`.
- Types, structures, classes, and `Prop`-valued definitions use
  `UpperCamelCase`. All other terms use `lowerCamelCase`.
- Treat acronyms as a unit, cased by what the first character would be.
- Name a theorem after its statement, not its proof. Describe the conclusion
  first, then hypotheses in order of appearance, joined by `_of_`.
- Use the standard vocabulary of name fragments for logical connectives, set
  operations, algebraic operations, and relation properties. Prefer the short
  established forms — `pos`, `neg`, `nonneg`, `nonpos` — over spelling out a
  comparison against zero.
- Use `le` and `lt` for the natural argument order; reserve `ge` and `gt` for
  statements whose arguments are genuinely swapped.
- Place declarations in the namespace of their principal subject so dot notation
  works at the call site, and name extensionality and injectivity results with
  the conventional suffixes.
- A primed name marks a variant of an existing result, not a second attempt.
  Explain in the docstring how it differs from the unprimed version.

## Documentation

- Every file opens with the copyright header, then imports, then a module
  docstring delimited by `/-!` and `-/` on their own lines.
- The module docstring starts with a top-level Markdown header and a summary,
  then covers main definitions, main statements, notation introduced,
  implementation notes, references, and tags. Omit a section only when it is
  genuinely empty; notation must be documented whenever any is introduced.
- Every definition and every major theorem carries a doc comment. Docstrings on
  minor lemmas are encouraged whenever they carry mathematical content or are
  used outside their own file.
- Write docstrings in complete sentences ending in periods. Reference other
  declarations in backticks so they render as links.
- A docstring conveys mathematical meaning; it may simplify away implementation
  detail, but it must never contradict the statement.
- Record why a definition takes the form it does — especially a choice made for
  definitional convenience — in the implementation notes rather than leaving the
  reader to reverse-engineer it.
- For research code, cite the paper, section, or equation a definition or result
  formalizes, and state explicitly where the formalization diverges from the
  source.

## Formatting

- Keep lines within 100 characters.
- Indent proof bodies two spaces. Indent continuation lines of a multi-line
  statement four spaces, so the proof body remains visually distinct.
- Put spaces on both sides of `:`, `:=`, and infix operators, and break the line
  after the operator rather than before it.
- Never leave `by` alone on a line — end the statement line with `:= by` and
  indent the tactic block beneath it.
- Introduce each new goal with a focusing dot and indent its block. Structure
  every branching proof this way rather than letting goals run together.
- Write anonymous functions with `fun` and `↦`. The `λ` spelling is disallowed.
- Use `<|` rather than `$`.
- Put a space after `←` in rewrite and simp arguments, and after binders.
- Write binder types explicitly even where they could be inferred.
- Construct structures and instances with `where` rather than enclosing braces.

## Definitions

- Disable automatic implicits at the file or project level and declare variables
  explicitly. Auto-bound implicits silently turn a typo into a new parameter.
- Use `Type*` and `Sort*` for universe polymorphism rather than an underscore,
  which can be unexpectedly specialized by later code.
- Prefer structural recursion. Annotate the termination argument explicitly even
  when it would be found automatically: it speeds elaboration, documents the
  argument, and prevents silently falling back to well-founded recursion.
- Reach for well-founded recursion only when the recursion is genuinely not
  structural. Its results are not definitionally equal to their return values
  and reduce slowly, so downstream proofs must rely on the equation lemmas.
- Avoid `partial` and `unsafe` in anything a proof depends on. `partial` gives
  up the equational theory of the definition, and `unsafe` propagates to every
  declaration that refers to it.
- Keep `Prop` for mathematical statements and `Bool` for computation. Decide the
  boundary deliberately and convert explicitly rather than letting coercions
  accumulate.
- Use `let` rather than `have` when introducing data — `have` discards the value
  and leaves the term unusable in later proofs.
- Use a subtype when a witness must be extracted and used; an existential
  statement is proof-irrelevant and its witness cannot be recovered
  constructively.
- State isomorphism of types with an equivalence, never with equality of types,
  which is badly behaved and usually neither provable nor disprovable.
- Mark a definition `noncomputable` deliberately rather than by accident; if the
  computational content matters to an experiment, choose a representation that
  keeps it.
- Avoid `Float` anywhere a proof is intended. Use rationals or reals for
  mathematics, and confine floating point to the boundary where results leave
  Lean.

## Typeclasses and Generality

- Assume the weakest typeclass that supports the proof. Strengthening an
  assumption to make one step easier narrows every downstream application.
- Never add an instance argument for a structure that already has a global
  instance. Doing so silently generalizes the statement to an arbitrary
  structure and can make it vacuous or simply false.
- Group shared hypotheses in `variable` blocks rather than repeating them, and
  keep those blocks small enough that a reader can see what is in scope.
- Prefer instance-implicit arguments for class assumptions and strict-implicit
  binders where an ordinary implicit would leave unsolved metavariables when the
  function is mentioned without arguments.
- Generalize a result when the generalization is natural and the proof is
  unchanged. Do not contort a proof to reach a generality nothing needs.

## Statement Hygiene

- Total functions in Lean return junk values outside their mathematical domain:
  division by zero is zero, natural subtraction truncates at zero, and functions
  like square root, logarithm, derivative, and infinite sum are defined
  arbitrarily where they are undefined.
- Because of this, a statement can be trivially true or quietly false rather
  than ill-formed. Add the hypotheses that constrain the argument to the real
  domain, and never read a proved theorem as evidence that the edge case was
  handled.
- Avoid natural subtraction in statements. Reformulate to add on the other side
  of the relation instead, and where subtraction is unavoidable, carry the
  hypothesis that keeps the expression out of the truncation regime.
- Prefer `<` and `≤` over `>` and `≥` in statements; the library is oriented
  that way and the flipped forms fail to match.
- Be explicit about numeric types in statements involving literals, division, or
  coercion — the type determines which operation is meant.
- Keep coercions at the boundary of a statement rather than scattered through
  it, and normalize them with the dedicated cast tactics instead of by hand.
- State membership and subset hypotheses with explicit element and proof
  arguments rather than coercing a set to a type, which complicates rewriting
  and induction.
- Read every statement back as mathematics before proving it. A proof of the
  wrong statement is the most expensive failure mode in this repository.

## Proofs

- Never commit `sorry`. A file containing one is unfinished work, not a partial
  result, and the warning it emits must never be suppressed.
- Every `simp` must close its goal. Mid-proof, use `simp only` with an explicit
  lemma list, fold the simplification into a `have` or `suffices`, or finish
  with `simpa` — a bare non-terminal `simp` breaks as soon as the simp set
  changes.
- Replace an exploratory `simp` with the explicit lemma list it reports before
  committing.
- Rewrite under binders with the dedicated tactics or conversion mode; ordinary
  rewriting cannot reach inside a bound variable.
- Do not depend on automatically generated hypothesis names. Name what you
  introduce, and select goals by their case name rather than positionally where
  the proof branches.
- Prefer a structured proof — explicit intermediate statements, named
  hypotheses, focusing dots — over a long opaque tactic chain. Intermediate
  statements are where a reader checks the argument.
- Use term mode where it is shorter and clearer than tactic mode, particularly
  for applying a lemma directly or constructing a structure.
- Extract a reusable step into its own lemma with a proper name and docstring
  instead of repeating it. Private helper lemmas are preferable to duplication.
- Prefer decision procedures and automation for the goals they are designed for,
  but do not paper over a failing step with heavier automation — find out why it
  fails first.
- When a proof breaks after a library update, fix the argument rather than
  pinning around it.

## Simp Set Discipline

- Tag a lemma `@[simp]` only when it rewrites toward a genuinely simpler form,
  and when the result is the library's normal form for that concept.
- Ensure the left-hand side's arguments are themselves in normal form, since
  simplification works from the inside out.
- Never tag a lemma whose rewrite can loop or make the term grow, and never tag
  one that merely happens to be useful in a particular proof — invoke it by name
  instead.
- Aim for confluence: the outcome should not depend on the order rewrites are
  attempted.
- Adding a simp lemma changes every downstream proof. Treat it as a change to a
  shared interface, and check the library still builds.

## Trust and Verification

- Confirm that finished results depend only on the standard axioms by printing
  the axiom dependencies of the top-level theorems. This is the check that
  catches an unnoticed `sorry`.
- Treat `decide` as a proof and native evaluation as an extension of the trusted
  base: native computation is recorded as an extra axiom per computation and is
  trusted rather than kernel-checked. Never use it to close a result the
  repository presents as proved.
- Run the library linters before considering a file finished, including the
  checks for missing docstrings and simp-normal-form violations.
- Treat every warning as a defect. Unused-variable warnings in particular
  usually mean the statement does not say what was intended.
- Use `example` blocks to pin down intended behavior of a definition and to
  guard against a redefinition silently changing its meaning.
- Keep definitions, the results proved about them, and the code that executes
  them in agreement — a definition that has drifted from the experiment it
  models invalidates both.
