# Skills, agents & library for the Publishing phase (figures first)

Brainstorming for tooling that covers the back end of the research process:
turning finished artifacts — `docs/theorems/`, `docs/results/`, `docs/conj/`
pages — into blog posts and papers. Publishing is a wide phase (figure design,
post drafting, paper scaffolding); this doc scopes the **figure slice**, which
both downstream deliverables consume and which is where the repo currently has
the least support. Post/paper drafting is sketched only far enough to make
sure the figure tooling feeds it cleanly.

Two gaps motivate this:

1. **No identification scaffolding.** Nothing helps decide *which* figures a
   results doc or theorem doc warrants. The exploration design deliberately
   kept figure generation out of `explore-results` (see the decision at the
   end of `experimentation-exploration-skills.md`) with the note that keepers
   "graduate to `src/figures`" — but the graduation path was never built.
2. **No standard visualization library.** `figures.style` (validated blog
   palette + rcParams) and `figures.cli` (the `@figure` harness) exist and are
   solid, but everything between them and a finished figure is bespoke: all
   six existing modules hand-roll the same `_box` / `_arrow` / `_edge_label` /
   panel primitives, and nothing at all exists for data-driven charts from
   experiment DataFrames.

## House patterns to reuse

1. **Checklist skills with scratchpad working files** surviving compaction;
   interactive phases structured around explicit user checkpoints.
2. **Fresh-eyes auditor agents** gating publication-bound artifacts
   (findings-list output, severity-ordered, `CLEAN` sentinel, loop until
   clean).
3. **Scout agents keeping bulk content out of the main context** — for
   figures, `data-profiler` already covers the data side.
4. **The docstring as the figure's argument.** Every existing figure module
   opens with prose stating what claim the figure carries and why each element
   sits where it does (`matroid_page.py` is the exemplar). This is the
   publishing analog of pre-registration: the claim is stated before the
   reader sees the picture, and an auditor can check the picture against it.

## Current state (what already exists)

- `figures.style` — blog-exact palette with documented obligations: the
  categorical set fails colorblind-separation gates, so color must never
  carry series identity alone; indianred/slate carry the blog's
  falsified/confirmed callout semantics; sequential/diverging/ordinal ramps
  are validated.
- `figures.cli` — `@figure` decorator: standard options, `FigureContext.save`
  into `docs/fig/<name>/`, blog style applied automatically, stable stems
  because posts reference outputs by path.
- Six bespoke diagram modules (`process_*`, `matroid_page`, `palette_demo`),
  each duplicating box/arrow/panel plumbing (~100–300 lines of it per module).
- `src/research` structure plot methods (`plot_lattice_of_flats`,
  `plot_decorated_permutation`, ...) with the shared `_plot.py` contract:
  draw onto given axes, never `show`.
- The `dataviz` skill (form heuristic, mark specs, palette validator) — used
  to build `figures.style`; remains the reference for chart-form choices.

## Design principles specific to this phase

- **A figure is an argument, not an illustration.** Identification means
  extracting from a doc the small set of claims that are better carried
  visually than in prose, then choosing the form that carries each claim —
  not decorating every section. A candidate figure is specified as:
  *claim carried* → *form* → *data source* → *encoding plan*.
- **Style obligations become API, not convention.** The palette's documented
  hazards (no identity by color alone, ≤3 series in scatter, assign slots in
  order, red = falsified) are currently docstring text that every new module
  must re-read and re-honor. Chart builders should enforce them by
  construction: secondary encoding applied automatically, series caps raised
  as errors, prediction outcomes colored through named semantic constants.
- **Figures regenerate from checked-in sources.** A data figure reads
  `data/results/<name>.<ts>.json` (path recorded in the module); a structure
  figure derives from `src/research` objects; a diagram derives from constants
  in the module. Nothing is hand-edited post-render; re-running overwrites in
  place. The regeneration command appears in the module docstring.
- **Raw DataFrames stay out of the main context.** Figure design consults the
  run through `data-profiler` or scripted pandas; only the implementation
  reads it, in code.
- **Dual targets, one source.** Blog (SVG, ~700 px column) and paper
  (PDF/PNG, LaTeX column widths) want different sizes/fonts of the *same*
  figure. Better one module rendering both than parallel modules drifting
  apart — the harness already supports multi-format saves; sizing presets are
  the missing piece.

## Proposed artifact chain

```
docs/results/<name>.<ts>.md | docs/theorems/<name>.md | docs/conj/<name>.md
  → figure plan                          (design-figures, interactive)
  → src/figures/<name>.py               (add-figure, audited against the plan)
  → docs/fig/<name>/*.{svg,png,pdf}     (regenerable; referenced by path)
  → blog post / paper draft             (later slice, consumes docs/fig)
```

## Library

### 1. `figures.diagram` — extract the duplicated primitives

Mechanical consolidation of what the six modules already share, keeping the
canvas-inches idiom they all use:

- fixed-frame canvas bootstrap (figure sized to the drawing, full-figure axes,
  equal aspect, cropped by construction — the `_render` preamble every module
  repeats);
- rounded `_box` panels; titled list panels; header bands;
- edge-anchored arrows: the `_box_exit` / `_edge_endpoints` geometry (arrows
  meet box *edges* at any approach angle, bowed arcs handled) plus edge
  labels;
- reveal sequencing: the build-up pattern from `matroid_page` (render step
  1..N into an identical frame so slides stack without jitter) as a helper.

Migration of the existing modules is verifiable by image diff — outputs must
be pixel-identical (or reviewed-equivalent) after the refactor. This layer
has immediate value independent of everything else in this doc.

### 2. `figures.charts` — data-driven builders over experiment results

The genuinely new layer: builders consuming tidy DataFrames from
`data/results` runs. Initial set, driven by what results docs will contain:

- **sweep lines** — metric vs. swept axis, one series per arm, replication
  spread as bands; secondary encoding (markers + direct labels) applied
  automatically.
- **prediction outcomes** — the pre-registration verdict figure: per-region
  predicted vs. observed, colored by confirmed (slate) / contradicted
  (indianred) / no-prediction, via semantic constants, matching the blog's
  callout colors.
- **parameter-grid heatmap** — metric over a 2-D sweep, `blog-sequential`
  (or diverging around a meaningful zero).
- **facet grid** — small multiples over a third axis, the escape hatch the
  palette's series cap forces.
- **run comparison** — outcome deltas between two timestamps, mirroring
  `explore-results`' compare mode.

Builders take a DataFrame plus explicit column semantics (value/axis/series
column names), draw onto provided axes (the `_plot.py` contract), and raise
on obligation violations rather than documenting them. They do **not** hide
matplotlib — a figure module composes builders and then adjusts, because
publication figures always need final-mile control.

Growth is demand-driven: build a chart form when the first results doc needs
it, generalizing from that concrete use. No `data/results` runs exist yet, so
committing to a large chart API now would be speculation.

### 3. Structure plots stay where they are

`src/research` plot methods remain structure-owned (they are analysis tools
too). Figure modules compose them under the blog style — which they already
inherit via rcParams — adding only publication chrome. No adapter layer
needed until one is actually missed.

## Skills

### 1. `design-figures` (interactive)

The identification scaffolding. Input: a results doc, theorem doc, conjecture
doc, or a post outline; the skill reads it plus its neighbors (design doc,
Lean sources, the run's profile via `data-profiler` — never the raw frame).

Steps sketch:

1. Gather: the doc's claims, its pre-registered-prediction outcomes or lemma
   DAG, what structures involved already have `plot_*` methods, what figures
   already exist in `docs/fig`.
2. Propose candidates — each stated as *claim carried / form / data source /
   encoding plan / where it sits in the post or paper* — deliberately few,
   ranked, with "not worth a figure" rationale for prose-sufficient claims.
   A doc-type playbook seeds the proposals:
   - *results doc*: the prediction-verdict figure; one sweep figure per
     headline trend; a surprise figure only if it survived bug-first triage.
   - *theorem doc*: the object itself (via structure plots); the lemma DAG /
     proof structure; a scope figure (proven range vs. experimentally tested
     range) — the honest picture of what the theorem does and does not cover.
   - *conjecture doc*: evidence plot with tested region marked; counterexample
     search coverage map.
3. Checkpoint with the user: select, refine forms, fix sizing targets.
4. Record the plan and stub follow-ups for unselected candidates.

The skill's value is the checklist of figure-design failure modes: form
mismatched to claim (the `dataviz` heuristic), series overload vs. this
palette, missing uncertainty representation, unlabeled units/axes, a figure
implying a claim the data can't support (drift, visually), blog/paper sizing
picked too late.

### 2. `add-figure`

Implements `src/figures/<name>.py` from the plan, following the harness
conventions (`@figure`, argument-docstring, stable stems, constants over
magic numbers, `figures.diagram`/`figures.charts` before bespoke drawing).

Steps sketch: read plan → survey conventions → implement → pilot render to
the scratchpad and *look at it* (Read the PNG) → iterate → audit gate
(`figure-auditor`) → final render into `docs/fig/<name>/` → report with the
`just figure <name>` invocation and the paths written.

## Agents

### `figure-auditor` (gates `add-figure`)

Fresh-eyes audit of the rendered figure + module against the plan and the
source doc. Given the plan, the module, the rendered outputs, and (for data
figures) the result path; checks:

- **[claim]** the figure supports exactly the claim its docstring states —
  not a stronger one (the visual analog of claim-strength drift).
- **[num]** data-derived marks are recomputable from the referenced result
  file (recompute via pandas; compare against extracted mark values or the
  module's computation, not by eye alone).
- **[style]** palette obligations honored: secondary encoding present
  wherever color distinguishes series, series caps respected, slots assigned
  in order, semantic red/slate not repurposed.
- **[decode]** every encoded channel is decodable — legend or direct labels,
  axis units, colorbar for ramps.
- **[repro]** regenerates deterministically from checked-in sources; stems
  stable; regeneration command in the docstring; no hand-placed value that
  should be computed.

Findings-list output, `CLEAN` sentinel, same contract as the other auditors.

## Later slices (out of scope here, noted for fit)

- **`draft-post`** — from a theorem/results doc plus its `docs/fig` figures
  to a blog post draft; needs a decision on where drafts live (this repo vs.
  the blog repo) before it can be designed.
- **Paper scaffolding** — LaTeX (or Typst) skeleton consuming the same
  figures as PDF; brings citation management (`docs/ref` → BibTeX) with it.
- Both consume the figure tooling above unchanged; nothing here should need
  rework when they arrive.

## Build order

1. **`figures.diagram`** — pure consolidation of existing duplication,
   verifiable against current outputs, and every later figure builds on it.
2. **`design-figures` + the plan convention** — useful immediately for the
   existing process/page figures' successors and for the first theorem doc;
   establishes the claim-first discipline everything else audits against.
3. **`add-figure` + `figure-auditor`** — the graduation path the exploration
   phase left open. Degrades gracefully: with no plan doc, the auditor audits
   against the docstring argument alone.
4. **`figures.charts`** — grown builder-by-builder as real results docs
   arrive; the first prediction-verdict and sweep figures extracted into the
   library from their first concrete use.

## Open questions

- **Where does the figure plan live?** Options: (a) the module docstring only
  — the house "docstring as argument" pattern, no new artifact, but then the
  plan exists only after implementation; (b) `docs/fig/<name>/plan.md`
  alongside the outputs; (c) a design section appended to the *source* doc
  (results/theorem doc gets a "Figures" section). Leaning (c): it keeps the
  claim-to-figure mapping in the doc being published from, and `add-figure`
  copies the argument into the docstring where the auditor checks it. -
  DECISION:

- **Should identification be a closing step of `explore-results` /
  `attack-conjecture` instead of a separate skill?** The earlier decision
  kept figure *generation* out of exploration, but a cheap "candidate
  figures" list in the results/theorem doc (claims only, no implementation)
  may beat re-reading the doc later in `design-figures`. Separate skill still
  owns the interactive refinement either way. - DECISION:

- **Does `figures.charts` bind to the experiment output schema or stay
  schema-agnostic?** Binding (reading `.meta.json`, knowing arm/seed/axis
  conventions) makes the prediction-verdict figure nearly automatic;
  agnostic (tidy frame + column names) keeps the library usable for
  non-experiment data. Leaning agnostic core + a thin experiment-aware
  wrapper. - DECISION:

- **Blog sizing contract.** Is there a fixed content-column width (px) and
  font-size floor the blog imposes, worth encoding as a sizing preset next to
  `figures.style`? Requires numbers from the blog repo. Same question for the
  eventual paper class. - DECISION:

- **SVG font handling.** SVGs embed font *names*; a reader without Roboto
  Mono falls back. Accept fallback, embed subset fonts, or render SVG text as
  paths (`svg.fonttype: 'path'`, larger files, text no longer selectable)?
  Affects `figures.style`, not per-figure code. - DECISION:
