---
name: add-lean-topic
description: Formalizes a mathematical object from the Logseq knowledge graph in Lean 4 — reusing Mathlib wherever it already covers the topic and extending it in src/theorems where definitions or theorems are missing. Use when the user asks to formalize, prove, or build a Lean or Mathlib representation of a topic that has a Logseq page.
argument-hint: [topic]
allowed-tools: Task, Read, Grep, Glob, Write, Edit, Bash, mcp__logseq
---

# Formalize a Logseq topic in Lean

Formalize **$ARGUMENTS** in `src/theorems`. The Logseq page is the
specification, but Mathlib is the library of record: every claim on the page
is first mapped to an existing Mathlib declaration, and only genuine gaps —
missing definitions, missing theorems — are formalized here, in Mathlib's
conventions, so they could migrate upstream without a rewrite. Never redefine
what Mathlib already has: a second definition of the same object connects to
nothing downstream and forks the library.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Read the topic's Logseq page
- [ ] Step 2: Survey Mathlib and src/theorems
- [ ] Step 3: Write the coverage map
- [ ] Step 4: State before proving
- [ ] Step 5: Audit the statements
- [ ] Step 6: Prove
- [ ] Step 7: Run the full gate
- [ ] Step 8: Close the loop in Logseq and report
```

## Working file

Before surveying, create `coverage.md` in the session scratchpad directory.
Steps 2–3 write it; Steps 4–8 read from it. It is the ground truth that
survives context compaction: wherever conversation memory and the file
disagree, trust the file.

## Step 1: Read the topic's Logseq page

Fetch the page with `get_page` (titles are usually singular; fall back to
`search_logseq` with singular/plural variants). The page must actually
specify the object: a populated **Definition** section, **structural
theorems**, and **canonical examples**. If the page is missing or too thin to
formalize from, stop and tell the user to run `/add-logseq-topic` first — do
not substitute your own research; the graph is the specification.

Read the page's existing **Implementation notes — Lean** section as well: the
research skill writes it as speculation about what Mathlib provides. Treat
its module paths and declaration names as scout leads, not facts — Step 2
verifies them. Step 8 replaces the speculation with what was actually found
and built.

## Step 2: Survey Mathlib and src/theorems

Confirm the pinned Mathlib checkout exists at `.lake/packages/mathlib/`
(run `just lean-update` if it does not — it also fetches the build cache,
without which everything downstream is unusably slow).

**Keep Mathlib source out of this conversation.** Dispatch `mathlib-scout`
subagents — one per claim cluster (definition and axiomatizations; derived
vocabulary; operations and constructions; structural theorems; canonical
examples, which Mathlib often provides as named constructions), launched in
parallel. Give each the topic, its claims quoted verbatim from the page, and
any leads from the page's Lean notes; each returns a disposition per claim —
covered, partial, or missing — with locally verified declaration names and
module paths. Append scout output to `coverage.md` as it returns.

Also read what exists in `src/theorems/Theorems/` — prior topic modules are
targets for cross-topic theorems and set the local style.

## Step 3: Write the coverage map

Complete `coverage.md` following [reference/design.md](reference/design.md):
for every page claim, the verbatim claim text, the scout's disposition, and
the decision — **reuse** (cite the Mathlib declaration; nothing to build),
**extend** (a planned repo declaration about Mathlib's object), **define**
(Mathlib lacks the object itself), or **backlog** (no tractable plan; report,
do not build). Where the page's primary axiomatization differs from the one
Mathlib chose, work with Mathlib's — state missing cryptomorphisms as
theorems or constructors; do not re-derive the object. A claim with no row in
the map is not ready to formalize.

## Step 4: State before proving

Write `src/theorems/Theorems/<Topic>.lean` from `coverage.md`, following
[reference/proving.md](reference/proving.md): copyright header, module
docstring, minimal imports, then every planned declaration with its full
statement, docstring, and attribution — with `sorry` for each proof body.
`sorry` is scaffolding that exists only between here and the end of Step 6;
it is never committed. Add the module to the imports in
`src/theorems/Theorems.lean` and confirm the statements elaborate with
`lake build Theorems.<Topic>`. Re-read the relevant rows of `coverage.md`
before stating them rather than reciting them from memory.

## Step 5: Audit the statements

Do not spend proof effort on unaudited statements. Dispatch a
`statement-auditor` subagent with the paths to `coverage.md` and the new
module; it returns a severity-ordered findings list, or `CLEAN`. Fix every
finding — a meaning finding usually means editing the statement, not
defending it — then re-dispatch. Repeat until it returns `CLEAN`: only a
clean audit unlocks Step 6, because a proved wrong statement is worth less
than no statement at all.

## Step 6: Prove

Replace every `sorry`, iterating with `lake build Theorems.<Topic>` for fast
feedback. Before investing in a hard proof, hunt for counterexamples with
`plausible` and check the statement against the page's canonical examples.
Treat a stuck proof as a three-way question — the statement is subtly wrong,
the proof is genuinely hard, or the page's claim is false — and diagnose
which before pushing harder. A false page claim (especially a `plausible`
counterexample) is a research finding: report it, do not weaken the statement
to something provable but unclaimed. A genuinely hard proof gets its
statement removed and recorded as backlog in `coverage.md` — never a
committed `sorry`. The formalization is an executable audit of the graph.

## Step 7: Run the full gate

Run `just lean-check` and fix findings until it passes clean, treating every
warning as a defect. Then audit axioms per
[reference/proving.md](reference/proving.md): each top-level theorem must
depend on nothing beyond `propext`, `Classical.choice`, and `Quot.sound` —
`sorryAx` means unfinished work escaped Step 6, and native-evaluation axioms
mean a result the repository would present as proved is merely trusted.

## Step 8: Close the loop in Logseq and report

Update the page's **Implementation notes — Lean** section with what is now
verified fact: the Mathlib declarations and module paths that cover the
topic, the repo module's path and its declaration names, and what remains
missing. Respect the Logseq MCP sharp edges documented in
[../add-logseq-topic/reference/logseq-api.md](../add-logseq-topic/reference/logseq-api.md):
`get_block` before `update_block`, preserve UUID link refs verbatim, write
serially, and on a write timeout re-read before retrying.

Then report to the user: the coverage map summary (reused / extended /
defined / backlog counts and highlights), each new declaration and the page
claim it formalizes, upstream candidates (extensions in a Mathlib namespace
that Mathlib itself might accept), the formalization backlog, any page
discrepancies found in Steps 5–6, and any page **Open questions** that the
formalization makes precise enough to attack.
