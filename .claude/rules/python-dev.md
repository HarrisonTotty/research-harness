---
paths:
  - "**/*.py"
  - "**/*.pyi"
---

# Python Development

Target runtime is **Python 3.14**. Write for that version — do not add
compatibility shims, version guards, or backports for older interpreters.

## Typing

- Annotate every public function signature, method signature, and module-level
  constant. Local variables only need annotations where inference is ambiguous.
- Annotations are deferred and evaluated lazily as of 3.14, so forward
  references resolve without quoting. `from __future__ import annotations` is
  still supported but opts back into the older stringizing behavior — do not add
  it to new modules.
- Prefer built-in generics over the deprecated `typing` aliases, and the union
  operator over `Union`. Write `X | None`, not `Optional[X]`.
- Declare type aliases with the `type` statement rather than the deprecated
  explicit alias annotation; it supports forward references natively.
- Use inline type-parameter syntax for generic functions, classes, and aliases
  rather than explicitly constructed type variables.
- Accept the widest reasonable abstract type in parameters and return the most
  specific concrete type. Prefer abstract collection protocols for inputs.
- Mark parameters that must not be mutated with immutable or read-only types
  rather than relying on convention.
- Use structural typing (protocols) for duck-typed interfaces instead of forcing
  callers into a nominal hierarchy. Runtime-checkable protocols verify only that
  members exist, not their signatures — a passing check is not conformance.
- Use `object` for a value that could be any type but must be narrowed before
  use; reserve `Any` for genuinely dynamic boundaries with untyped data.
- Use literal, enum, and typed-mapping types to make invalid states
  unrepresentable instead of validating stringly-typed arguments at runtime.
- For custom type predicates prefer `TypeIs` over `TypeGuard`: it narrows on the
  false branch too, and requires the narrowed type to be a subtype of the input.
  Use `TypeGuard` only when that subtype relationship genuinely does not hold.
- Annotate functions that never return normally with `Never` rather than the
  older equivalent spelling.
- Type checking must pass with strict settings. Silence a diagnostic only with a
  narrow, specific ignore comment plus a one-line justification.

## Docstrings

- Every module, public class, function, and method carries a docstring. Private
  helpers need one whenever the intent is not obvious from the signature.
- Follow PEP 257. Always use triple double quotes, and raw triple double quotes
  if the text contains backslashes.
- Open with a one-line summary on the same line as the opening quotes, phrased
  as a command — "Return the ..." not "Returns the ..." — ending in a period.
- One-line docstrings keep the closing quotes on that line, with no blank line
  before or after. Multi-line docstrings put a blank line after the summary and
  the closing quotes on a line of their own.
- Insert a blank line after every class docstring.
- Document arguments, return values, side effects, exceptions raised, and any
  restriction on when the callable may be invoked. Use one structured style
  consistently across the repository; do not mix styles.
- Never repeat the signature or restate types already present in it. Document
  meaning, units, valid ranges, and ownership instead.
- For research code, cite the source of a formula or algorithm — paper, section,
  or equation number — in the docstring of the implementing function.
- Module docstrings state what the module is for and how its pieces fit
  together, not a list of its contents.

## Errors and Exceptions

- Raise the most specific built-in exception that fits. Define a module-level
  base exception and derive from it only when callers need to catch your errors
  as a group.
- Never write a bare `except:` — it swallows interpreter exit and keyboard
  interrupts. Catch specific exception types; fall back to `except Exception:`
  only at a top-level boundary that logs and re-raises.
- Keep the `try` clause to the minimum code that can raise, so it cannot mask
  failures from unrelated statements.
- Multiple exception types no longer need parentheses when there is no `as`
  clause.
- Never swallow an exception silently. Handle it, re-raise it, or wrap it with
  explicit chaining that preserves the cause. Suppress chaining explicitly when
  the cause is genuinely irrelevant, so the omission reads as a decision.
- Use exception groups and `except*` when multiple independent failures can
  occur, particularly in concurrent code.
- Add context to error messages: what was attempted, with which inputs. Messages
  should be actionable without a debugger.
- Never `return`, `break`, or `continue` out of a `finally` block — it discards
  in-flight exceptions and raises a syntax warning as of 3.14.
- Validate inputs at API boundaries and trust them internally. Do not
  defensively re-validate on every internal call.
- Use assertions only for invariants that indicate programmer error. Never use
  them to validate external input or enforce runtime contracts.

## Data Modeling

- Prefer immutable structures. Use frozen dataclasses or named tuples for value
  objects; reserve mutable classes for things with genuine identity.
- Never use a mutable default argument. Use a sentinel and construct the default
  inside the function body.
- Use enums for closed sets of values rather than string or integer constants.
- Compare against `None` and other singletons with `is` / `is not`, never with
  equality operators. Test types with instance checks, not by comparing types.
- Keep `__init__` cheap and side-effect free; put expensive construction behind
  explicit classmethod constructors.
- Implement `__repr__` for any class appearing in logs, test failures, or
  interactive sessions.
- Use keyword-only parameters for optional and boolean flags so call sites stay
  readable, and positional-only parameters where the name is an implementation
  detail.
- Be consistent about returns: if any path returns a value, every path does, and
  bare exits are written as an explicit `return None`.
- Use `match` for structural dispatch over heterogeneous data; prefer plain
  conditionals for simple value comparisons.

## Iteration and Comprehensions

- Use a comprehension when it builds a collection from a single expression with
  at most one condition and fits readably on a few lines. Otherwise write a loop.
- Never write a comprehension for its side effects; use a `for` loop.
- Avoid nesting beyond two clauses — it reads worse than the loop it replaces.
- Prefer generator expressions over list comprehensions when the result is
  consumed once, especially for large or streaming data.
- Iterate directly over objects; use enumeration and zipping helpers rather than
  indexing by range. Require strict length agreement when zipping sequences that
  must be parallel.
- Do not mutate a collection while iterating it. Build a new one.
- Return iterators from library functions producing large sequences, and make
  laziness explicit in the docstring so callers know the result is single-use.

## Concurrency

- Choose by workload: async for I/O-bound concurrency, threads for I/O and for
  extension calls that release the interpreter lock, separate processes or
  interpreters for CPU-bound work. Threads do not speed up pure-Python
  computation on a default build.
- Do not hand-roll thread or process management. Use the standard executor and
  task-group abstractions, which propagate exceptions and guarantee cleanup.
- Prefer structured concurrency: scope every task to a task group so a failure
  cancels its siblings and no task outlives its parent.
- Never block the event loop with synchronous I/O or CPU-heavy work. Offload it.
- Always await or explicitly cancel tasks you create. The loop holds only weak
  references, so retain a strong reference to any fire-and-forget task and drop
  it on completion, or it may be collected mid-execution.
- Treat cancellation as expected: let it propagate, never swallow it, and use
  shielding only for cleanup that must complete. Suppressing it breaks the
  timeout and task-group machinery built on top of it.
- Prefer message passing over shared mutable state. Where state must be shared,
  guard it with an explicit lock and document the invariant the lock protects.
- Acquire locks in a consistent global order and hold them for the shortest
  possible span. Use context managers, never manual acquire/release.
- Do not rely on the interpreter lock, or on the internal locking of built-in
  containers, to make a sequence of operations atomic. Neither is a guarantee,
  and free-threaded builds are officially supported as of 3.14.
- Set explicit timeouts on every network and inter-process wait.
- Keep concurrency at the edges of the codebase. Numerical and model code should
  be pure and synchronous so it stays testable.

## Resource Management

- Manage every file, socket, lock, and connection with a context manager.
- Write custom context managers for paired setup/teardown rather than relying on
  callers to remember cleanup. Use the exit-stack pattern for a dynamic number.
- Always specify encoding explicitly when opening text files.
- Use the path abstraction rather than string manipulation for filesystem paths.

## Testing

- Every public function gets tests. Every bug fix gets a regression test that
  fails before the fix.
- One behavior per test, named for the expected behavior and the condition under
  which it holds.
- Structure tests as arrange/act/assert with no branching. A test containing
  conditionals should be several tests or a parameterized one.
- Parameterize over cases rather than looping inside a test, so each case
  reports independently.
- Assert on specific values and specific exception types with matched messages.
  Never assert merely that something is truthy or that nothing was raised.
- Tests must be deterministic and order-independent. Seed every random generator
  explicitly and inject clocks rather than reading wall time.
- Use temporary directories for filesystem tests; never write into the source
  tree or depend on files left behind by another test.
- For numerical code, assert with explicit absolute and relative tolerances
  chosen for the problem rather than copied. Test boundary conditions, empty
  inputs, and known analytic solutions.
- Prefer property-based tests wherever an invariant can be stated — round-trips,
  symmetries, and conservation laws.
- Fake only what you own or what crosses a process boundary. Heavy mocking of
  internal calls tests the implementation rather than the behavior.
- Keep fixtures small and composable, with the narrowest scope that works.

## Modules and Safety

- One concern per module. Split a module before it becomes a grab bag.
- Use absolute imports, never wildcard imports. Import modules rather than
  individual names when the source matters for readability.
- Keep imports at module top level. Defer an import into a function only to
  break a genuine cycle or to avoid an expensive optional dependency, with a
  comment saying which.
- Declare the public surface explicitly and prefix internal names with a single
  underscore.
- Guard script entry points so importing a module never executes work, and do no
  I/O, network calls, or expensive computation at import time.
- Avoid mutable module-level state. Pass configuration explicitly.
- Pass logging arguments for lazy interpolation instead of pre-formatting the
  message, and log tracebacks from within an exception handler using the
  dedicated exception-logging call. Never log secrets, credentials, or full
  input payloads.
- Build SQL, shell commands, and markup with template strings or a parameterized
  interface, never by interpolating values into an f-string.
- Never use `eval` or `exec` on data from outside the repository, and never
  unpickle data from an untrusted or tamperable source — it executes arbitrary
  code. Use a structured, schema-checked format for untrusted input, and sign
  payloads that must survive an untrusted channel.
