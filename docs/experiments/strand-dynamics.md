# Strand dynamics on canonical plabic graphs

Strand dynamics is a cellular automaton on the vertex colors of a plabic
graph whose update rule is driven by the graph's own strands (trips). This
experiment runs it from the canonical coloring of **every fixed-point-free
positroid cell** of every $\mathrm{Gr}(k,n)$, $2 \le n \le 8$, on **two
canonical graphs per cell** — the bridge graph and the Le-graph — and
records the transient and period of each trajectory. It exists to port the
predecessor repository's `strand_gradient` sweep ("Finding 020") into this
harness, to check that the port reproduces it, and to ask the first new
question the port makes cheap: which of Finding 020's observations belong
to the cell, and which to the graph chosen to represent it.

## Hypothesis

On **bridge graphs**, (H1) the dynamics reproduce the per-Grassmannian
summary table of Finding 020 exactly, and (H2) the cells of
$\mathrm{Gr}(m,2m)$ whose canonical coloring is a fixed point of the rule
are exactly the non-crossing perfect matchings of $[2m]$ — so there are
$C_m$ of them — while no fixed-point-free cell of any other
$\mathrm{Gr}(k,n)$ is a fixed point.

## Provenance

- The rule and the sweep: `~/gh/research` at `f7f0f01` —
  `src/research/plabic_automata/automaton.py` (`StrandGradientRule`),
  `src/research/experiments/plabic_automata_cli.py` (`run_automaton`,
  `sweep`), and the write-up
  `docs/findings/020-strand-gradient-landscape.md` (2026-03-27), which
  states H2 as a conjecture verified through $m = 5$.
- The rule is original to this project; there is no Logseq page and no
  literature prediction for it. The bridge construction is the BCFW-bridge
  decomposition (Arkani-Hamed et al., arXiv:1212.5605;
  Fomin-Williams-Zelevinsky section 7.10), the Le-graph is Postnikov
  section 20.
- Port: `research.strand_dynamics.StrandAutomaton` and
  `research.plabic_graph.PlabicGraph.from_bridge_decomposition`. Before
  this doc was written the port was checked cell-by-cell against the old
  implementation: all 2,175 fixed-point-free cells with $n \le 7$ give
  graphs identical to the old ones and identical `(transient, period)`.
  H1 therefore carries real risk only at $n = 8$ (and $n = 9$ if run),
  where the old implementation was not re-run and the comparison is
  against Finding 020's published table.

## Falsification criterion

- **H1** is falsified by any bridge-arm row of the summary table below
  whose cell count, period-2 count, period-$>2$ count, period-1 count,
  maximum period, or maximum transient differs from Finding 020's value,
  or whose mean period differs after rounding to two decimals.
- **H2** is falsified by any bridge-arm cell with `period == 1` and
  `transient == 0` that is not a non-crossing perfect matching, or any
  non-crossing perfect matching of $[2m]$, $2m \le n_{\max}$, whose row
  does not have `period == 1` and `transient == 0`.

A bridge-arm row with `period == 1` and `transient > 0` falsifies neither
clause but is reported: Finding 020 counts "fixed points" by period and
its oracle has none with a transient for $n \le 7$.

## Parameter space

| Axis | Range | Granularity | Rationale |
| ---- | ----- | ----------- | --------- |
| `n` | 2..8 (default); 9 opt-in via `--n-max 9` | every integer | Finding 020 covers 2..9. $n \le 8$ is 17,008 cells and about two minutes; $n = 9$ adds 133,496 cells and is dominated by Le-graph construction (about 6 ms per cell), so it is opt-in. H2 needs even $n$: 2, 4, 6, 8 give $m = 1..4$. |
| cell | every fixed-point-free decorated permutation of $[n]$; $k$ is its anti-exceedance count (Postnikov's type), so $1 \le k \le n-1$ | exhaustive | Matches Finding 020, whose "all positroid cells" are the derangements (7 cells in $\mathrm{Gr}(2,4)$, not 33). A fixed point of the permutation is a lollipop, which has no internal neighbor, never flips, and shares no strand with the rest of the graph; cells with fixed points are left to a follow-up rather than mixed into the replication. |
| `graph` | `bridge`, `le` | both per cell | `bridge` is the replication arm. `le` is the library's canonical graph; the rule depends on the topology, and the two graphs differ for every cell with $n \le 7$ (different vertex counts; bridge graphs keep degree-2 vertices). |
| `max_steps` | 100,000 | fixed | Finding 020's budget; its longest transient + period for $n \le 9$ is under 30,000. |

Fixed, not swept: the rule's `sensitivity` is 1 (the only value Finding 020
used; the port compares in exact integers and does not expose it). The
initial coloring is the graph's canonical one — no random colorings.

## Controls and baselines

The `bridge` arm is the control: it is the old experiment, and its
expected output is known in advance (the table below). Any surprise in the
`le` arm is only interpretable once the bridge arm has matched, which is
why both arms are computed in the same run from the same cell enumeration.

## Replication and seeds

The experiment is deterministic: exhaustive enumeration, canonical initial
colorings, a deterministic rule. There is no randomness, hence no seed and
a replication count of 1. The metadata records this explicitly
(`"seeds": []`) so the absence reads as a decision.

## Output schema

One row per (cell, graph).

| Column | Dtype | Meaning / units |
| ------ | ----- | --------------- |
| `n` | int | Number of boundary vertices. |
| `k` | int | Anti-exceedance count of the permutation (Postnikov type). |
| `permutation` | str | One-line images $\pi(1),\dots,\pi(n)$, comma-separated, Postnikov's direction (the trip from $b_i$ ends at $b_{\pi(i)}$); the join key between arms. |
| `graph` | str | `bridge` or `le`. |
| `dimension` | int | Dimension of the positroid cell. |
| `internal_vertices` | int | Number of internal vertices of the graph (the automaton's state has this many bits). |
| `edges` | int | Number of edges of the graph. |
| `noncrossing_matching` | bool | Whether the permutation is an involution whose arcs pairwise do not cross (it is fixed-point-free by construction). |
| `converged` | bool | Whether a coloring repeated within `max_steps` steps. |
| `transient` | int, null if not converged | Steps before the trajectory enters its cycle. |
| `period` | int, null if not converged | Length of the cycle of colorings. |
| `permutation_period` | int, null if not converged | Minimal period of the trip-permutation sequence around the cycle (a true period, unlike Finding 020's first-recurrence gap). |
| `cycle_permutations` | int, null if not converged | Number of distinct trip permutations on the cycle. |

## Pre-registered predictions

All predictions are for the `bridge` arm and are transcribed from Finding
020's "Per-Grassmannian Summary" (source for every row:
`~/gh/research/docs/findings/020-strand-gradient-landscape.md`). "Fixed"
is the number of period-1 cells.

| Region | Prediction | Source |
| ------ | ---------- | ------ |
| every bridge row | `converged` is true | Finding 020, "Overall Convergence" |
| $\mathrm{Gr}(1,2)$ | 1 cell, period 1, transient 0 | Finding 020 table |
| $\mathrm{Gr}(1,n)$ and $\mathrm{Gr}(n-1,n)$, $n \ge 3$ | 1 cell each, period 2, transient 0 | Finding 020 table and the note under it |
| $\mathrm{Gr}(2,4)$ | 7 cells: 5 period-2, 0 longer, 2 fixed; mean 1.71, max period 2, max transient 2 | Finding 020 table |
| $\mathrm{Gr}(2,5)$ | 21: 20 / 1 / 0; mean 2.10, max 4, transient 2 | Finding 020 table |
| $\mathrm{Gr}(3,5)$ | 21: 17 / 4 / 0; mean 2.48, max 6, transient 2 | Finding 020 table |
| $\mathrm{Gr}(2,6)$ | 51: 46 / 5 / 0; mean 2.41, max 10, transient 10 | Finding 020 table |
| $\mathrm{Gr}(3,6)$ | 161: 130 / 26 / 5; mean 2.63, max 26, transient 18 | Finding 020 table |
| $\mathrm{Gr}(4,6)$ | 51: 44 / 7 / 0; mean 2.43, max 8, transient 14 | Finding 020 table |
| $\mathrm{Gr}(2,7)$ | 113: 90 / 23 / 0; mean 2.59, max 10, transient 80 | Finding 020 table |
| $\mathrm{Gr}(3,7)$ | 813: 650 / 163 / 0; mean 3.78, max 169, transient 112 | Finding 020 table |
| $\mathrm{Gr}(4,7)$ | 813: 629 / 184 / 0; mean 3.42, max 126, transient 101 | Finding 020 table |
| $\mathrm{Gr}(5,7)$ | 113: 101 / 12 / 0; mean 2.36, max 10, transient 29 | Finding 020 table |
| $\mathrm{Gr}(2,8)$ | 239: 186 / 53 / 0; mean 3.53, max 48, transient 54 | Finding 020 table |
| $\mathrm{Gr}(3,8)$ | 3,361: 2,395 / 966 / 0; mean 6.81, max 642, transient 680 | Finding 020 table |
| $\mathrm{Gr}(4,8)$ | 7,631: 5,240 / 2,377 / 14; mean 7.80, max 2,654, transient 1,509 | Finding 020 table |
| $\mathrm{Gr}(5,8)$ | 3,361: 2,346 / 1,015 / 0; mean 5.77, max 896, transient 786 | Finding 020 table |
| $\mathrm{Gr}(6,8)$ | 239: 194 / 45 / 0; mean 2.83, max 18, transient 85 | Finding 020 table |
| $\mathrm{Gr}(2,9)$ (opt-in) | 493: 359 / 134 / 0; mean 5.90, max 258, transient 90 | Finding 020 table |
| $\mathrm{Gr}(3,9)$ (opt-in) | 12,421: 7,786 / 4,635 / 0; mean 14.62, max 2,614, transient 2,470 | Finding 020 table |
| $\mathrm{Gr}(4,9)$ (opt-in) | 53,833: 32,633 / 21,200 / 0; mean 26.85, max 17,790, transient 10,784 | Finding 020 table |
| $\mathrm{Gr}(5,9)$ (opt-in) | 53,833: 32,495 / 21,338 / 0; mean 23.70, max 11,558, transient 14,051 | Finding 020 table |
| $\mathrm{Gr}(6,9)$ (opt-in) | 12,421: 7,973 / 4,448 / 0; mean 12.13, max 3,420, transient 2,780 | Finding 020 table |
| $\mathrm{Gr}(7,9)$ (opt-in) | 493: 380 / 113 / 0; mean 3.69, max 100, transient 105 | Finding 020 table |
| $\mathrm{Gr}(m,2m)$, bridge | fixed points = non-crossing perfect matchings: 1, 2, 5, 14 for $m = 1..4$ | Finding 020, "Fixed Points Are Counted by Catalan Numbers" |
| bridge, $k$ vs $n-k$ | period distributions of $\mathrm{Gr}(k,n)$ and $\mathrm{Gr}(n-k,n)$ differ for every pair with $2 \le k < n/2$ (the table's rows already differ from $n = 5$ on) | Finding 020, "The k vs n-k Asymmetry" |
| **every `le` row** | **no prediction** — nothing in the literature or in Finding 020 speaks to Le-graphs: not convergence within the budget, not the period distribution, not whether any fixed points exist or are Catalan-counted, not whether period agrees with the bridge arm cell-by-cell | — |
| `permutation_period`, `cycle_permutations`, both arms | **no prediction** — Finding 020 measured a first-recurrence gap it itself flags as an artifact; the true period is new | — |

Disclosure: sizing this sweep ran both arms for $n \le 8$ before this doc
was written. The only `le`-arm outcomes seen were: no unconverged cells
for $n \le 7$; maximum period 26 ($n=6$) and 84 ($n=7$); maximum transient
21 and 120. No `le`-arm fixed-point count, period distribution, or
cell-by-cell comparison was looked at, and no $n = 8$ outcome of either
arm.

## Analysis plan

1. **Integrity.** Row count is $2 \sum_n D_n$ with $D_n$ the derangement
   numbers (34,016 for $n = 2..8$); each `(n, permutation)` appears once
   per arm; every row converged; `permutation_period` divides `period`;
   `k` ranges over $1..n-1$.
2. **H1.** Group the bridge arm by `(n, k)`: cell count, counts of
   `period == 2`, `period > 2`, `period == 1`, mean and max `period`, max
   `transient`; compare against the table above.
3. **H2.** On the bridge arm, compare the set `period == 1 & transient ==
   0` with the set `noncrossing_matching`, per $n$; list any
   `period == 1` row with a transient, and any outside $k = n/2$.
4. **Representative dependence (exploratory).** Join the arms on
   `(n, permutation)`: fraction of cells with equal `period`, with equal
   `(transient, period)`; the `le`-arm fixed-point set against
   `noncrossing_matching`; the `le`-arm version of the H1 table; how
   `period` scales with `internal_vertices` in each arm.
5. **Duality (exploratory).** Per arm, compare the period distributions
   of `(k, n)` and `(n-k, n)`.
6. **Trip permutations (exploratory).** Distribution of `period /
   permutation_period` and of `cycle_permutations`; the fraction of cells
   whose cycle holds a single trip permutation.

## Implementation notes

- **Module:** `src/experiments/strand_dynamics.py` (CLI name
  `strand-dynamics`), composing `research.strand_dynamics.StrandAutomaton`,
  `PlabicGraph.from_bridge_decomposition`, and
  `LeDiagram.from_decorated_permutation(pi.inverse()).to_plabic_graph()` —
  the inverse because the Le-graph's trip permutation is the inverse of its
  diagram's decorated permutation. Every row is guarded by a check that the
  built graph's trip permutation equals the row's `permutation`.
- **Pilot** ($n = 2..6$, both arms, into the session scratchpad): 642 rows
  = $2 \times (1 + 2 + 9 + 44 + 265)$ in 0.9 s; schema and dtypes as
  designed; all rows converged. The bridge arm's `(n, k)` summary equals
  the prediction table above for every row with $n \le 6$ — expected, since
  the port had already been checked cell-by-cell through $n = 7$.
- **Runtime** (single process, measured while sizing): $n = 7$ about 8 s,
  $n = 8$ about 105 s, of which about 90 s is Le-graph construction (6 ms
  per cell) and under 10 s the bridge arm. $n = 9$ was not piloted; expect
  a quarter of an hour of Le-graph construction alone plus the dynamics,
  which Finding 020 reports as heavy-tailed at $n = 9$ (periods to 17,790).
  There is no replication variance: the run is deterministic.
- **Memory:** `StrandAutomaton.orbit` keeps every visited coloring until
  the trajectory closes; at $n \le 9$ that is at most tens of thousands of
  short tuples per cell, released between cells.
- **Full sweep:** `just experiment strand-dynamics` (the designed
  $n = 2..8$, both arms, 34,016 rows, about two minutes). With Finding
  020's full range: `just experiment strand-dynamics --n-max 9`.
