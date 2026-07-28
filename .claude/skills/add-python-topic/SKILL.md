---
name: add-python-topic
description: Implements a mathematical object from the Logseq knowledge graph as a Python structure in src/research — visualizable, DataFrame-serializable, with transformations to/from other objects and computed properties — plus tests derived from the page's theorems and examples. Use when the user asks to implement, code up, or build a Python representation of a topic that has a Logseq page.
argument-hint: [topic]
allowed-tools: Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Implement a Logseq topic in Python

Implement **$ARGUMENTS** as a Python structure in `src/research`. The Logseq
page is the specification: its definition fixes the representation, its
derived vocabulary becomes computed properties, its operations become
transformation methods, its theorems become the property-test oracle, and its
canonical examples become test fixtures. The finished object must be
visualizable, serializable to a pandas DataFrame, convertible to and from
related objects, and able to compute its own useful properties.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Read the topic's Logseq page
- [ ] Step 2: Survey src/research and the toolchain
- [ ] Step 3: Write the design spec
- [ ] Step 4: Implement the module
- [ ] Step 5: Write the tests
- [ ] Step 6: Run the full gate
- [ ] Step 7: Close the loop in Logseq and report
```

## Working file

Before designing, create `spec.md` in the session scratchpad directory. Step 3
writes it; Steps 4–7 read from it. It is the ground truth that survives
context compaction: wherever conversation memory and the file disagree, trust
the file.

## Step 1: Read the topic's Logseq page

Fetch the page with `get_page` (titles are usually singular; fall back to
`search_logseq` with singular/plural variants). The page must actually
specify the object: a populated **Definition** section, **structural
theorems**, and **canonical examples**. If the page is missing or too thin to
implement from, stop and tell the user to run `/add-logseq-topic` first —
do not substitute your own research; the graph is the specification.

Also read the pages linked from the topic's **Operations and constructions**
and **Abstracts / generalizes** blocks — those are the candidate targets for
transformation methods. Note which linked structures are red links (no page
yet): they cannot be transformation targets today.

Read the page's existing **Implementation notes — Python** section as well:
the research skill writes it as implementation hints (natural
representations, which theorems become property tests). It is design input
for Step 3; Step 7 replaces its speculation with what was actually built.

## Step 2: Survey src/research and the toolchain

Read what exists in `src/research` — module layout, naming, and any objects
this one should transform to or from — and check `pyproject.toml` for
available dependencies. Existing modules override anything a generic design
would suggest: extend the established shape rather than introducing a second
style. Visualization needs a plotting library; if none is present, plan to
add `matplotlib` with `uv add matplotlib` in Step 4.

## Step 3: Write the design spec

Write `spec.md` following [reference/design.md](reference/design.md). It
must contain the page-section → code-artifact mapping table, the chosen
internal representation (and why that axiomatization won), the constructor
list, and the full method inventory across the four required capabilities —
properties, transformations, DataFrame serialization, visualization — each
entry citing the page block (definition, theorem, or example) it implements.
A method with no supporting block on the page is not ready to implement;
either take it back to the graph via `/add-logseq-topic` or drop it.

## Step 4: Implement the module

Implement `src/research/<topic>.py` from `spec.md`, following the design
contract in [reference/design.md](reference/design.md). The house Python
rules apply as always; the skill-specific obligations are the docstring
citations (every formula and algorithm cites the source the Logseq page
attributes it to) and the four capabilities. Re-read the relevant section of
`spec.md` before implementing it rather than reciting it from memory.

## Step 5: Write the tests

Write `tests/test_<topic>.py` following
[reference/testing.md](reference/testing.md): fixtures transcribed from the
page's canonical examples (each asserting exactly what the page says the
example certifies), Hypothesis property tests transcribed from the
structural theorems, the round-trip laws (DataFrame, constructors,
involutions), and the axiom-violation tests (each numbered axiom rejected
with a `ValueError` that names it).

## Step 6: Run the full gate

Run `just check` and fix findings until it passes clean. Treat a failing
property test as a three-way question — implementation bug, mistranscribed
test, or a wrong claim on the page — and diagnose which before editing
anything. If the page's claim is wrong, do not silently fix or weaken the
test: report the discrepancy to the user, since the page (and possibly its
sources) needs correcting. The implementation is an executable audit of the
graph.

## Step 7: Close the loop in Logseq and report

Update the page's **Implementation notes — Python** section with the real
module path, class name, constructor names, and which theorems became which
property tests. Respect the Logseq MCP sharp edges documented in
[../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md):
`get_block` before `update_block`, preserve UUID link refs verbatim, write
serially, and on a write timeout re-read before retrying.

Then report to the user: what was built and where, the test inventory
(which theorem backs which test), any page discrepancies found in Step 6,
the transformation backlog (targets skipped because their structure has no
page or no implementation yet — each is a candidate for a future
`/add-logseq-topic` or `/add-python-topic` run), and any **Open questions**
from the page that the new object makes cheap to explore as experiments in
`src/experiments`.
