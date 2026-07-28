# Theorems

Proven results, with their proofs. This section is the prose counterpart to the
Lean library: it carries the statement, the argument, and the context a reader
needs, while the machine-checked version lives in `src/theorems/` and is
rendered in the [Lean API](../reference/lean/index.md).

A page in this section should carry:

- the statement, in the same form as its Lean counterpart;
- a proof written for a human reader, not a transcription of the tactic script;
- the name of the Lean declaration that discharges it, so the two can be checked
  against each other;
- what it depends on, and what it is used to prove.

A result that is stated but not yet proven belongs in
[Conjectures](../conj/index.md) until the proof closes.
