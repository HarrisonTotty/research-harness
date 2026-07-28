---
name: add-logseq-topic
description: Researches a topic from official and original sources, then creates a comprehensive Logseq page containing its full formal definition, properties, examples, and references. Use when the user asks to add, create, or research a topic page in Logseq, or to expand the knowledge graph with a new concept.
argument-hint: [topic]
allowed-tools: WebSearch, WebFetch, mcp__logseq
---

# Add a researched topic page to Logseq

Research **$ARGUMENTS** and publish it as a comprehensive Logseq page. The
finished page must contain the topic's full formal definition (all standard
equivalent formulations), its major properties and theorems with attribution,
canonical examples, and references that prioritize original sources.

Copy this checklist and track progress through it:

```
Progress:
- [ ] Step 1: Check the graph for existing coverage
- [ ] Step 2: Research from primary sources
- [ ] Step 3: Draft the page structure
- [ ] Step 4: Audit the draft against the fact map
- [ ] Step 5: Write the page to Logseq
- [ ] Step 6: Verify the published page
- [ ] Step 7: Report with sources
```

## Working files

Before researching, create two files in the session scratchpad directory:

- `sources.md` — the running fact → source map. Step 2 appends to it;
  Steps 3, 4, and 7 read from it.
- `draft.md` — the complete page draft. Step 3 writes it; Steps 4–6
  read from it.

These files are the ground truth that survives context compaction. Keep
them current as you work; wherever conversation memory and file content
disagree, trust the files.

## Step 1: Check the graph for existing coverage

Search before creating: `search_logseq` for the topic name **and** its
singular/plural variants and synonyms (pages are usually titled in the
singular). If a page already exists, switch to reviewing and extending it
instead of creating a duplicate. Also collect the exact titles of related
existing pages so the new page can link to them.

## Step 2: Research from primary sources

Follow [reference/research.md](reference/research.md) for the
source-priority ladder and verification rules. The core rules:

- Prioritize the **original source** (the paper or spec that introduced the
  concept) and **official references** (standards documents, authoritative
  databases such as OEIS/DLMF, maintained official docs). Use encyclopedias
  only for orientation and leads — verify anything you keep against a
  higher-tier source.
- Every theorem or named result gets an attribution (who, year). Every
  numeric claim gets checked against an authoritative database or computed
  directly. Do not publish a claim you could check but did not.

**Keep raw sources out of this conversation.** Dispatch one
`source-reader` subagent per source (or per subtopic when one source
covers several), launched in parallel; each returns distilled facts with
provenance, pre-formatted for `sources.md`. Append each agent's output to
`sources.md` as it returns. See "Delegating the reading" in research.md
for the acceptance rules. Fetch directly yourself only for quick lookups
(a single OEIS entry, resolving a DOI) — never for papers or long
documents.

## Step 3: Draft the page structure

Write the complete draft to `draft.md`, following the section template in
[reference/page-template.md](reference/page-template.md). It encodes the house
style: formal definition with all standard equivalent axiomatizations, a
theorems section written to double as a property-test oracle, examples written
as test fixtures, and implementation notes tied to this repository (Python
research library, Lean/Mathlib). Build the draft only from facts recorded in
`sources.md` — if a fact is not in the map, it is not ready to publish. Any
claim about this repo — module paths, API names, toolchain versions — must be
verified against the working tree before it goes in the draft.

## Step 4: Audit the draft against the fact map

Do not publish an unaudited draft. Dispatch a `draft-auditor` subagent
with the paths to `sources.md` and `draft.md` and the topic focus; it
returns a severity-ordered findings list, or `CLEAN`. Fix every finding:
edit `draft.md` for drift, attribution, and consistency problems; for an
unsupported claim, either re-research it (a new `source-reader` dispatch,
appending to `sources.md`) or cut it from the draft. Then re-dispatch the
auditor. Repeat until it returns `CLEAN` — only a clean audit unlocks
Step 5.

## Step 5: Write the page to Logseq

Follow [reference/logseq-api.md](reference/logseq-api.md) exactly — the
Logseq MCP tools have sharp edges (UUID link refs, append-only ordering,
transient timeouts). The short version: create the page, then transcribe
`draft.md` one section at a time (`create_block` for the section parent,
`create_blocks` for its children), strictly serially, in final display
order. Re-read the relevant section of `draft.md` before writing it rather
than reciting it from memory.

## Step 6: Verify the published page

Re-read the whole page with `get_page` and diff it against the `draft.md`
file — read the file; do not verify against your memory of what you
drafted: every section present and correctly nested, math delimiters
intact, `[[...]]` links spelled exactly as drafted. Run
`find_missing_pages` and confirm each unresolved link is an intentional red
link (future work), not a typo of an existing page's title. Fix
discrepancies with `update_block` / `delete_block`, then re-read once more.
Only finish when the published page matches `draft.md`.

## Step 7: Report with sources

Summarize for the user: what was created, which claims were verified and how,
and a Sources list built from `sources.md`, with markdown links to the URLs
actually used. Also list
the red links the new page introduced — each is a candidate for a future
`/add-logseq-topic` run — and note any red-link targets that other pages
already point to, since shared targets are the highest-value next topics.
