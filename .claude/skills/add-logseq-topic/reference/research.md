# Research and verification rules

## Delegating the reading

Raw source text (papers, book chapters, long reference pages) must not enter
the main conversation — it crowds out working context and rots under
compaction. The **`source-reader`** subagent (defined in `.claude/agents/`)
owns the reading: given one source and a topic focus, it returns distilled
facts with provenance, formatted for direct appending to `sources.md`.

- Dispatch one `source-reader` per source (or per subtopic when one source
  covers several). Give each the exact URL/DOI and the topic focus, and
  run independent dispatches in parallel.
- Append each agent's output to `sources.md` (see SKILL.md, "Working
  files") as soon as it returns.
- Accept nothing without provenance. If a definition or formal theorem
  statement comes back paraphrased without its verbatim quote, re-dispatch
  for the exact text.
- Delegation does not waive verification: numeric claims reported by a
  subagent still get checked against an authoritative database or
  recomputed directly by you.
- Fetch directly only for quick, small lookups — a single OEIS entry,
  resolving a DOI, confirming a title.

## Source priority

Work down this ladder; cite the highest tier that actually supports the claim.

1. **Original sources.** The paper, book, spec, or RFC that introduced the
   concept. Locate via DOI, arXiv, journal archive, or publisher. For
   historical works, a scanned original or a faithful republication counts;
   note the year and venue.
2. **Official references.** Standards documents, maintained official
   documentation, and authoritative curated databases: OEIS (integer
   sequences), DLMF/NIST (special functions and constants), IETF RFCs,
   ISO/ECMA standards, language references, upstream project docs.
3. **Canonical textbooks and monographs.** The standard graduate reference
   for the field (e.g. Oxley for matroid theory, Diestel for graph theory).
   Use for consolidated statements of definitions and theorems.
4. **Peer-reviewed surveys.** Good for a map of the area, attribution
   chains, and open problems.
5. **Encyclopedias** (Wikipedia, nLab, Encyclopedia of Mathematics,
   MathWorld). Orientation and leads only. Never cite as the sole support
   for a claim that a higher tier could confirm — climb the ladder via
   their references instead.

## Verification requirements

- **Definitions** are stated in full formal precision: quantifiers explicit,
  axioms individually numbered, ground set / domain / types stated. When a
  concept has several standard equivalent definitions (cryptomorphisms,
  alternative axiomatizations), include each one and say what connects them.
- **Theorems** carry attribution: name (if named), discoverer(s), year.
  If two sources disagree on attribution, note the discrepancy rather than
  picking silently.
- **Numeric data** (enumeration counts, sequence values, constants) is
  checked against an authoritative database or recomputed directly. Small
  finite claims (e.g. "these 7 triples cover every pair exactly once") are
  cheap to verify by hand or with a few lines of Python — do it.
- **One-sided implications** are flagged as such. Do not let a
  characterization claim drift into a biconditional the source does not
  support.
- **Claims about this repository** (Mathlib module paths, API names,
  toolchain pins, library capabilities) are verified against the working
  tree — read the actual files under `.lake/packages/` or `src/` — never
  from memory of what an upstream library "should" contain.
- The running fact → source map lives in the `sources.md` working file, not
  in conversation memory. Every entry in the final References section must
  have supported at least one claim on the page, and every non-obvious
  claim must trace to an entry.
- The claim → fact tracing is not self-certified: before publication, the
  **`draft-auditor`** subagent (defined in `.claude/agents/`) reads
  `sources.md` and `draft.md` with fresh eyes and reports unsupported
  claims, claim-strength drift, attribution mismatches, unused facts, and
  internal inconsistencies. Publication is gated on a `CLEAN` audit
  (SKILL.md, Step 4).
