# Testing contract for topic implementations

The Logseq page was written to be executable: its theorems section is titled
"use these as the property-test oracle" and its examples "test fixtures".
This file says how to transcribe them. Tests live in `tests/test_<topic>.py`.

## Fixtures from canonical examples

- One named constructor per canonical example (e.g. `fano_matroid()`),
  exported by the topic module — canonical examples are library API, not
  test-only helpers — matching the page's data exactly.
- Each example block on the page says what the example certifies ("smallest
  X that is not Y"). Write exactly that assertion as a test — the
  certification is the test's reason to exist, and its name should say so.
- Counterexamples test the failing direction explicitly: assert the property
  they violate *and* the property they retain, so the one-sidedness the page
  records is enforced, not just remembered.

## Property tests from structural theorems

- One Hypothesis property test per structural theorem, named after the
  theorem, with the attribution (name, year) in the docstring.
- Write one shared Hypothesis strategy that generates small random valid
  instances through the public constructors, and derive constrained variants
  from it. Keep generated sizes small — these laws are combinatorial and
  explode quickly; a bounded ground set that runs in milliseconds beats a
  general one that times out.
- Transcribe claim strength faithfully. A one-sided implication on the page
  is tested one-sided; testing the converse of a theorem the page does not
  state is asserting something the graph never claimed.
- Quantitative theorems (bounds, counts) assert the exact inequality or
  value; for numeric properties compare exactly when the mathematics is
  exact, and with explicit chosen tolerances only when it is not.

## Round-trip laws

Beyond the page's theorems, the API contract itself yields properties —
test each that applies:

- `from_dataframe(x.to_dataframe()) == x` for generated instances.
- `from_<formulation>` constructors agree: building the same object through
  two axiomatizations yields equal values.
- Involutions and inverses the page records (`x.dual().dual() == x`,
  restriction/extension pairs) hold on generated instances.
- Invalid inputs: each numbered axiom has a test feeding a violating input
  to the constructor and asserting `ValueError` with the axiom named in the
  message.

## Oracle discipline

A failing property test means one of three things: the implementation is
wrong, the test mistranscribed the theorem, or the page's claim is false.
Diagnose which — shrink the failing case by hand, check it against the
page's verbatim statement — before editing anything. Never weaken a test to
make it pass: if the page is wrong, that is a research finding; stop and
report it (SKILL.md, Step 7). Standard suite hygiene (determinism, seeded
randomness, tmp dirs) is covered by the house Python rules and applies here
unchanged.
