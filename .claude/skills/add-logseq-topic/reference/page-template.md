# Page structure template

The house style, modeled on the existing **Matroid** page. Adapt section
titles to the topic's field (a protocol page has no "cryptomorphism table"),
but keep the ordering and the intent of each section: definition first,
theorems as a test oracle, examples as fixtures, then repo-specific
implementation notes.

## Page skeleton

```
- <intro block: one-paragraph plain-language statement of what the topic is>
  - **Also known as:** <synonyms and historical names>
  - **Type:** <e.g. mathematical structure, algorithm, protocol>
  - **Introduced:** <person(s), year, original publication>
  - **Abstracts / generalizes:** <the concrete notions it unifies, as [[links]]>
- ## Overview
  - <history and motivation: who introduced it, what problem it solved,
     independent discoveries, how the modern formulation emerged>
- ## Definition
  - <primary formal definition: ground set/domain stated, axioms numbered>
  - <each standard equivalent definition as a sibling block, one per
     axiomatization: independence, bases, circuits, rank, closure, ...>
  - <equivalence/cryptomorphism table mapping between formulations, if the
     topic has several>
- ## Derived vocabulary
  - <one block per term: formal definition in a sentence, [[link]] if the
     term has (or deserves) its own page>
- ## Operations and constructions
  - <duality, minors, sums, products, restrictions — whatever applies>
- ## Structural theorems — use these as the property-test oracle
  - <one block per theorem: precise statement in KaTeX, attribution
     (name, year), children for corollaries or the quantitative form.
     Written so each can be transcribed into a property-based test.>
- ## Canonical examples — test fixtures
  - <one block per example: the standard small objects, extremal cases, and
     counterexamples. Say what each one certifies (e.g. "smallest X that is
     not Y") — that is what makes it a useful fixture.>
- ## Implementation notes — Python
  - <how the topic maps onto `src/research`: existing modules to build on,
     natural representations, which theorems become property tests>
- ## Implementation notes — Lean
  - <what Mathlib provides: exact module paths and declaration names,
     verified against `.lake/packages/mathlib`; what is missing; how the
     repo's `src/theorems` could use or extend it>
- ## Open questions
  - <known open problems, plus conjectures worth exploring in this sandbox>
- ## References
  - <original sources first, then official references, then textbooks and
     surveys; each as a markdown link with author, title, year>
```

## Formatting conventions

- Math in KaTeX: `$...$` inline, `$$...$$` for display blocks.
- One fact per block; supporting detail goes in child blocks, not run-on
  paragraphs. Blocks should survive being read (and transcluded) alone.
- Cross-link concepts with `[[Page Name]]` — link on first meaningful use.
  **Deliberately create red links** to pages that do not exist yet: every
  term that deserves its own page (derived vocabulary, named theorems,
  related structures) should be linked even though the target is missing.
  Red links are the graph's backlog — `find_missing_pages` turns them into
  a queue of future topics, and the links resolve retroactively once the
  page is written. Do not link generic words that will never warrant a
  page; that turns the backlog into noise.
- Bold the lead term of a block (`**Rank function.**`) so the outline is
  scannable when collapsed.
- Counterexamples state their direction explicitly ("violates X, which every
  Y satisfies") — never leave a one-sided test reading as two-sided.
