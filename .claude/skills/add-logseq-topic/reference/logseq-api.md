# Logseq MCP tooling notes

Hard-won operational rules for the `logseq` MCP server. Violating these
produces broken links, duplicated blocks, or timeout cascades.

## Contents

- Finding pages
- Writing a new page
- UUID link refs when editing existing blocks
- Ordering and nesting
- Timeouts and retries
- Final verification

## Finding pages

- `get_page` requires the exact title. Pages are usually titled in the
  **singular** ("Matroid", not "Matroids"). When a lookup misses, fall back
  to `search_logseq` before concluding the page does not exist.

## Writing a new page

- Create the page with `create_page`, then write content **one section at a
  time**: `create_block` for the section heading, then `create_blocks` for
  that section's children in a single batch.
- Do not send one giant deeply-nested `create_blocks` call for the whole
  page — large hierarchical writes are the most common trigger for the
  30-second timeout. Section-sized batches land reliably.
- Issue MCP **write calls strictly serially**. Parallel writes make the
  Logseq API time out.

## UUID link refs when editing existing blocks

- Logseq stores page links inside existing blocks as UUID refs
  (`[[6a67da97-...]]`), not as `[[Page Name]]`.
- Before `update_block`, always `get_block` and reuse the UUID refs
  verbatim in the replacement content — rewriting them as page names breaks
  the links.
- In **newly created** blocks, plain `[[Page Name]]` is correct; Logseq
  resolves it (and unresolved names become useful red links).

## Ordering and nesting

- New blocks **append at the end** of a parent's children; there is no
  positional insert. Write sections and children in final display order.
- To place content under an existing section, pass that section block's
  UUID as the parent — capture UUIDs from each create call's response
  rather than re-fetching the page.

## Timeouts and retries

- `Request to Logseq timed out after 30000ms` is transient — it happens
  when the Logseq window is backgrounded or the app is busy digesting a
  prior write. Timed-out **reads** are safe to retry as-is.
- A timed-out **write may still have succeeded**. Before retrying any
  write, re-read (`get_block` / `get_page`) to check whether the content
  already landed — blind retries create duplicate blocks. Retry only what
  is confirmed missing.

## Final verification

- After all writes, `get_page` the full page and diff it against the draft:
  section order, nesting depth, math delimiters, and link refs. Fix with
  `update_block` / `delete_block`, then re-read once more.
- Run `find_missing_pages` to list unresolved `[[...]]` targets: intentional
  red links are fine, but a near-miss of an existing page's title (plural,
  typo, casing) is a broken cross-link — fix the spelling.
