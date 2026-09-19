# Strand dynamics on canonical plabic graphs: 20260919T140150Z run

On bridge graphs the port reproduces every row of Finding 020's
per-Grassmannian table for $2 \le n \le 8$, including the $n = 8$ rows that
had not been checked before (H1 confirmed), and the fixed points of the rule
are exactly the non-crossing perfect matchings — 1, 2, 5, 14 cells for
$m = 1..4$ (H2 confirmed). In the exploratory Le-graph arm — no prediction
was registered for it — the same Catalan fixed points recur unchanged, but
they are a property of the two canonical constructions, not of the cell: a
fixed point is a coloring whose load is constant on each component (an
elementary argument, not machine-checked, that every graph of the run
obeys), and the symmetric square graph of the top cell of $\mathrm{Gr}(2,4)$, not a
non-crossing matching, is one (see
[Critique addendum](#critique-addendum-fixed-points-are-constant-load-colorings)).
Finding 020's $k$ versus $n-k$ asymmetry appears to belong to the bridge
graph: on Le-graphs every measured quantity is invariant under
$\pi \mapsto w_0 \pi w_0$, which exchanges $k$ and $n-k$ — the
Le-construction is $w_0$-equivariant, the bridge construction is anchored
at label 1, and by the rule's mirror symmetry (checked on every graph of the
run, not proved) anchoring it at label $n$ would reverse the asymmetry.
Periods themselves are representative-dependent: the two graphs of a cell
agree on the period in 46% of cells.

## Run

- Invocation: `just experiment strand-dynamics`
- Design doc: [Strand dynamics on canonical plabic graphs](../experiments/strand-dynamics.md)
- Raw artifacts: `data/results/strand-dynamics.20260919T140150Z.json` and
  `.meta.json`.

| Field | Value |
| ----- | ----- |
| `git_commit` | `d4312c9d82e31592b260b977152885c6d751e68c-dirty` (see [Data integrity](#data-integrity): observed clean in-session; not verifiable from the artifacts) |
| `n_min`, `n_max` | 2, 8 |
| `graphs` | `bridge`, `le` |
| `max_steps` | 100,000 |
| `seeds` | `[]` (`deterministic: true`) |
| rows | 34,016 (0 unconverged) |
| elapsed | 112.262 s |
| versions | Python 3.14.7, pandas 3.0.5, research-harness 0.1.0 |

Every conclusion below is claimed only for fixed-point-free cells with
$2 \le n \le 8$, the canonical initial coloring, and these two graphs. The
run is an exhaustive, deterministic enumeration, so within that range there
is no sampling or seed uncertainty; all uncertainty is about what happens
outside it. The critique addenda step outside that scope in two labeled
places: two hand-built graphs (the contracted square and the star) and the
recolorings E1, E2 and N of the run's graphs.

## Data integrity

The profiler passed every check; two caveats, neither touching a data
value, come first.

- **`git_commit` carries `-dirty`; the tree was observed clean in-session,
  which the artifacts cannot verify.** The
  experiment computes the commit after it writes the result file, and its
  dirtiness test counts untracked files, so the run marks its own tree
  dirty. `git status` was clean immediately before the run (the session's
  starting snapshot) and immediately after it listed only this run's two
  artifacts, which would put the run at `d4312c9`; but that is a session
  observation, and neither the sidecar nor the repository can now verify
  it — see [Suspected bugs](#suspected-bugs).
  No data value is affected.
- **The analysis environment differs from the run's.** After the run,
  `scipy>=1.18.1` and `scipy-stubs` were added to `pyproject.toml` and
  `uv.lock` (uncommitted at the time of writing; both files were modified
  after the result was written). The experiment does not import scipy. The
  Spearman correlations below are the only numbers that depend on it, since
  pandas delegates `method="spearman"` to scipy; they were first computed
  without scipy as Pearson correlations of average ranks, with identical
  output.
- Passed: 34,016 rows, equal to $2 \sum_n D_n$ with per-arm counts
  1, 2, 9, 44, 265, 1,854, 14,833; the 13 designed columns with the designed
  dtypes and no nulls; each `(n, permutation)` once per arm, and the
  permutation sets equal an independently enumerated set of derangements;
  all rows converged, with the largest `transient + period` at 2,987 against
  a budget of 100,000; `permutation_period` divides `period`; `k` equals a
  recomputed anti-exceedance count and lies in $1..n-1$; `k`, `dimension`
  and `noncrossing_matching` agree between arms; the non-crossing matching
  counts are 1, 0, 2, 0, 5, 0, 14 for $n = 2..8$.

## Prediction outcomes

All predictions are for the bridge arm. "Cells: a / b / c" is period-2 /
period-$>2$ / period-1 counts, followed by mean period, maximum period and
maximum transient.

| Region | Predicted (source: Finding 020) | Observed | Outcome |
| ------ | ------------------------------- | -------- | ------- |
| every bridge row | converged | 17,008 of 17,008 converged | confirmed |
| $\mathrm{Gr}(1,2)$ | 1 cell, period 1, transient 0 | 1 cell, period 1, transient 0 | confirmed |
| $\mathrm{Gr}(1,n)$, $\mathrm{Gr}(n-1,n)$, $3 \le n \le 8$ | 1 cell each, period 2, transient 0 | 12 rows, each 1 cell, period 2, transient 0 | confirmed |
| $\mathrm{Gr}(2,4)$ | 7: 5 / 0 / 2; 1.71, 2, 2 | 7: 5 / 0 / 2; 1.71, 2, 2 | confirmed |
| $\mathrm{Gr}(2,5)$ | 21: 20 / 1 / 0; 2.10, 4, 2 | 21: 20 / 1 / 0; 2.10, 4, 2 | confirmed |
| $\mathrm{Gr}(3,5)$ | 21: 17 / 4 / 0; 2.48, 6, 2 | 21: 17 / 4 / 0; 2.48, 6, 2 | confirmed |
| $\mathrm{Gr}(2,6)$ | 51: 46 / 5 / 0; 2.41, 10, 10 | 51: 46 / 5 / 0; 2.41, 10, 10 | confirmed |
| $\mathrm{Gr}(3,6)$ | 161: 130 / 26 / 5; 2.63, 26, 18 | 161: 130 / 26 / 5; 2.63, 26, 18 | confirmed |
| $\mathrm{Gr}(4,6)$ | 51: 44 / 7 / 0; 2.43, 8, 14 | 51: 44 / 7 / 0; 2.43, 8, 14 | confirmed |
| $\mathrm{Gr}(2,7)$ | 113: 90 / 23 / 0; 2.59, 10, 80 | 113: 90 / 23 / 0; 2.59, 10, 80 | confirmed |
| $\mathrm{Gr}(3,7)$ | 813: 650 / 163 / 0; 3.78, 169, 112 | 813: 650 / 163 / 0; 3.78, 169, 112 | confirmed |
| $\mathrm{Gr}(4,7)$ | 813: 629 / 184 / 0; 3.42, 126, 101 | 813: 629 / 184 / 0; 3.42, 126, 101 | confirmed |
| $\mathrm{Gr}(5,7)$ | 113: 101 / 12 / 0; 2.36, 10, 29 | 113: 101 / 12 / 0; 2.36, 10, 29 | confirmed |
| $\mathrm{Gr}(2,8)$ | 239: 186 / 53 / 0; 3.53, 48, 54 | 239: 186 / 53 / 0; 3.53, 48, 54 | confirmed |
| $\mathrm{Gr}(3,8)$ | 3,361: 2,395 / 966 / 0; 6.81, 642, 680 | 3,361: 2,395 / 966 / 0; 6.81, 642, 680 | confirmed |
| $\mathrm{Gr}(4,8)$ | 7,631: 5,240 / 2,377 / 14; 7.80, 2,654, 1,509 | 7,631: 5,240 / 2,377 / 14; 7.80, 2,654, 1,509 | confirmed |
| $\mathrm{Gr}(5,8)$ | 3,361: 2,346 / 1,015 / 0; 5.77, 896, 786 | 3,361: 2,346 / 1,015 / 0; 5.77, 896, 786 | confirmed |
| $\mathrm{Gr}(6,8)$ | 239: 194 / 45 / 0; 2.83, 18, 85 | 239: 194 / 45 / 0; 2.83, 18, 85 | confirmed |
| $\mathrm{Gr}(k,9)$ (opt-in rows) | six rows of Finding 020 | not run (`n_max` = 8) | not tested |
| $\mathrm{Gr}(m,2m)$, bridge | fixed points = non-crossing perfect matchings: 1, 2, 5, 14 | set equality for every $n$; counts 1, 2, 5, 14 | confirmed |
| bridge, $k$ vs $n-k$ | period distributions differ for every pair $2 \le k < n/2$ | differ for all 6 pairs ($n = 5..8$) | confirmed |
| every `le` row | no prediction | see Findings and Unexpected observations | no prediction existed |
| `permutation_period`, `cycle_permutations` | no prediction | see Unexpected observations | no prediction existed |

## Findings

### H1: the port reproduces Finding 020 through $n = 8$

All 28 bridge-arm `(n, k)` rows equal Finding 020's table in all seven
compared statistics, with 0 mismatches. The mean is compared after rounding
to two decimals, as pre-registered; for example $\mathrm{Gr}(4,8)$ has an
unrounded mean of 7.8004 against a published 7.80.

This shows agreement of the port with the old implementation's published
$n = 8$ summary, which the design doc identified as the only part of H1
carrying real risk. Alternatives considered:

- *Implementation.* The comparison is of seven aggregates per row, not of
  cells, so compensating cell-level differences at $n = 8$ are not excluded
  in principle. Matching a maximum period of 2,654 and a maximum transient
  of 1,509 by coincidence is implausible, and the port was identical
  cell-by-cell for $n \le 7$ before the run.
- *Numerics.* Counts and maxima are exact integers; only the mean involves
  rounding, and it is the pre-registered comparison.
- *Selection, seed, grid.* Not applicable: every row was pre-registered and
  the enumeration is exhaustive and deterministic.

It would be falsified by a cell-by-cell comparison at $n = 8$ against the
old implementation finding a differing `(transient, period)`. The six
$n = 9$ rows remain untested.

*Critique addendum: that comparison was run and found no difference.* The
old implementation (`~/gh/research` at `f7f0f01`, clean tree; `cmd_sweep`'s
path — `automaton_from_permutation(..., "strand_gradient")` at sensitivity
1, canonical coloring, `run_automaton` with 100,000 steps) was run on all
14,833 derangements of $[8]$ and joined to this run's bridge rows on the
permutation. `k`, `transient`, `period`, `internal_vertices` and `edges`
each agree in 14,833 of 14,833 cells, with no cell missing on either side.
Controls on the join: the same mapping (old 0-indexed images plus one, no
inversion) agrees on all five columns in 2,175 of 2,175 cells for
$n \le 7$, and the inverse mapping agrees on `(k, transient, period)` in 0
of 44 cells at $n = 5$ and 2,333 of 14,833 at $n = 8$, so a direction error
would have been visible. (The $n = 5$ zero is carried by `k` alone —
inversion sends $k$ to $n - k$, which never matches at odd $n$ — and
`(transient, period)` by itself still agrees in 36 of 44; the $n = 8$ figure
is the discriminating one, with `(transient, period)` alone agreeing in
4,731 of 14,833.) The "compensating cell-level differences"
alternative above is thereby excluded for the bridge arm at $n \le 8$. Not
checked: graph identity at $n = 8$ beyond the vertex and edge counts.
(That a rerun of the old code reproduces Finding 020's published $n = 8$
rows was not computed directly, but follows: this run's aggregates equal
the table and its cells equal the old code's.) Scripts
and outputs are in the session scratchpad under `recompute/h1_n8/`.

### H2: bridge fixed points are the non-crossing perfect matchings

On the bridge arm the set of cells with `period == 1` and `transient == 0`
equals the set of non-crossing perfect matchings for every $n$ from 2 to 8:
1, 0, 2, 0, 5, 0, 14 cells, with both set differences empty. No bridge row
has period 1 with a positive transient, and no period-1 row lies outside
$k = n/2$.

This shows H2 for $m \le 4$ in this run; Finding 020 reports $m = 5$, which
this run does not reach. Four even values of $n$ cannot distinguish the
Catalan numbers from other sequences beginning 1, 2, 5, 14 — the claim with
content is the set equality, which is checked exactly. Alternative
considered — *implementation*: `noncrossing_matching` is computed from the
permutation alone and was recomputed independently by the profiler with 0
mismatches, so the equality is not an artifact of the flag sharing code
with the dynamics. It would be falsified by a fixed point at $n \ge 9$ that
is not a non-crossing matching; the other direction cannot fail as long as
the canonical graph of a non-crossing matching is a union of paths, which
is observed for $n \le 8$ and expected from the literature inference below,
not proved (next section).

#### Critique addendum: fixed points are constant-load colorings

Added in the critique session, after the audit; scripts and outputs are in
the session scratchpad (`recompute/constant_load.py`, `square_probe.py`,
`star_probe.py`), and everything here is post hoc.

*Reformulation.* A vertex of maximum load within a connected component of
internal vertices flips unless every neighbor ties it, so a coloring is a
fixed point iff its load is constant on each component (loads are
non-negative, and a component whose maximum is 0 is already constant, so
positivity is not needed). That maximum-principle argument is not machine-checked; the data
agree with it without exception. Rebuilding both graphs of all 17,008
cells: the initial load is constant on every component in exactly the 22
fixed cells of each arm and in none of the other 16,986; a recomputed
one-step fixedness equals the run's `period == 1 and transient == 0` in
17,008 of 17,008 rows per arm; no graph has a non-positive load. In both
arms those 22 graphs are exactly the ones with no vertex of degree at least
3 (unions of paths), and exactly the non-crossing matchings. Component by
component (`recheck/percomp/`): each arm has 2,602 components without a
vertex of degree at least 3, all of constant load, and 17,948 with one, none
of constant load. On the Le arm every constant-load component is a single
vertex with no internal neighbor, so the 22 Le fixed points are fixed
through the rule's never-flips clause; on the bridge arm none is.

*What H2 then says.* "Non-crossing matching implies fixed" is immediate:
the graph is a union of boundary-to-boundary paths, and both trips of a
path visit all of its vertices. The Catalan count is a corollary of "union
of paths". The content is the converse, in graph-level form: on bridge and
Le graphs with $n \le 8$, a component containing a vertex of degree at
least 3 never has constant canonical load.

*It is not a statement about cells.* The contracted square graph of the top
cell of $\mathrm{Gr}(2,4)$ — trip permutation `3,4,1,2`, a crossing
matching, `is_reduced()` true — has loads (12, 12, 12, 12) and is a fixed
point (transient 0, period 1) in both of its colorings, although both
canonical graphs of that cell have period 2 in the run's data. Every vertex of the square has
internal neighbors, so the rule's "no internal neighbor never flips" clause
plays no part. The clause gives a second family: the star with one
$n$-valent internal vertex, a reduced graph of the cell of
$\mathrm{Gr}(1,n)$ or $\mathrm{Gr}(n-1,n)$, is fixed for $n = 3, 4, 5$
(loads $2n$), against period 2 on both canonical graphs. The canonical
constructions break a symmetry the square keeps; which representatives of a
cell have constant load is open.

### Bridge-arm $k$ versus $n-k$ asymmetry

The period distributions of $\mathrm{Gr}(k,n)$ and $\mathrm{Gr}(n-k,n)$
differ on the bridge arm for all 6 pairs with $2 \le k < n/2$,
$5 \le n \le 8$, as predicted. The Le arm changes how this should be read:
see the first unexpected observation.

## Unexpected observations

Everything in this section is exploratory: the design doc registered no
prediction for the Le arm or for the trip-permutation columns, and the
regions below were singled out after seeing the data.

### The Le arm is exactly symmetric under $\pi \mapsto w_0 \pi w_0$

Let $w_0 : i \mapsto n + 1 - i$. Conjugation by $w_0$ moves 16,740 of the
17,008 cells and sends $k$ to $n - k$ in all 17,008. On the Le arm it
preserves `internal_vertices`, `edges`, `transient`, `period`,
`permutation_period` and `cycle_permutations` in 17,008 of 17,008 cells.
Consequently the Le-arm summary rows of $(n, k)$ and $(n, n-k)$ are
identical and the joint `(transient, period)` distributions are equal for
all 6 pairs. On the bridge arm the same map preserves the vertex and edge
counts in all cells but the period in only 10,036 and the transient in
7,118.

Inversion and rotation are not symmetries of either arm: on the Le arm the
period is preserved in 7,200 (inverse) and 7,316 (rotation) of 17,008
cells, on the bridge arm in 9,848 and 10,364.

This suggests that Finding 020's "$k$ vs $n-k$ asymmetry" is a property of
the bridge graph rather than of the positroid cell, for $n \le 8$.

Mechanism, tested here for $n \le 7$ (graphs rebuilt from the experiment
module; $n = 8$ was tested later, in the critique addendum below). A boundary-anchored isomorphism test of rotation
systems — $b_i \mapsto b_{n+1-i}$, counterclockwise successor matched to
clockwise successor — finds the Le-graph of $w_0 \pi w_0$ to be the mirror
image of the Le-graph of $\pi$ with the colors of its vertices of degree at
least 3 swapped and its degree-2 vertices left as they are, in 2,175 of
2,175 cells. The bridge arm is the discriminating control: the same test
succeeds for 0 of 2,175 bridge graphs, and the variant swapping every color
for 59 of 2,175 (and for 0 Le-graphs). Those 59 bridge cells are a positive
control for the argument that follows: they are exactly the bridge cells
whose graphs are mirror-isomorphic with colors ignored altogether (every
other bridge cell fails on topology), and the run's `(transient, period)` equals that of the
$w_0$-conjugate in all of them. Of the 2,175 cells with $n \le 7$, 35 are
self-conjugate ($w_0 \pi w_0 = \pi$) and compare a row with itself; 11 of
the 59 are among them, so the control carries 48 of 48 non-trivial
comparisons — as the Le-arm isomorphism carries 2,140 of 2,140. The
isomorphism is sufficient for the invariance, not necessary: the bridge
`(transient, period)` is $w_0$-invariant in 1,113 of the 2,175 cells,
1,078 of the 2,140 that are not self-conjugate. Under the rule, a mirror image with
swapped colors has the same trips as vertex sets and a degree-2 vertex
routes a trip identically in either color, so loads and hence flips would
correspond; that last step is an argument, not machine-checked, but the
17,008 of 17,008 data-level agreement is what it predicts.

*Critique addendum: the mechanism at $n = 8$, and the rule's symmetries
tested directly* (post hoc; scratchpad `recompute/mirror_n8/` and
`recompute/equivariance/`). The isomorphism test was first re-read and
found to be a genuine boundary-anchored rotation-system isomorphism (the
dart map is forced from the boundary, rejected on any conflict, degree or
color mismatch, and required total and injective); it reproduces every
$n \le 7$ count above. At $n = 8$ (14,833 cells, 233 self-conjugate): the
Le-graph of $w_0 \pi w_0$ is the mirror image with the colors of degree
$\ge 3$ swapped in 14,833 of 14,833 cells (14,600 of 14,600 that are not
self-conjugate) and with every color swapped in 0; bridge graphs pass the
first test in 0 cells and the full swap in 141 (112 not self-conjugate),
which are again exactly the cells passing with colors ignored, and the
run's `(transient, period)` equals the conjugate's in 141 of 141 (112 of
112). The second falsifier below — as first written, a failed isomorphism
at $n \ge 8$ — is thereby tested at $n = 8$ and survived; it now reads
$n \ge 9$. Sufficiency without necessity
persists: 4,889 of the 14,833 bridge cells are $w_0$-invariant.

The step flagged above as "an argument, not machine-checked" was then
tested on every graph of the run, 17,008 per arm, with the orbit of each
graph first recomputed and found equal to the recorded row in 34,016 of
34,016. Three recolorings preserve the orbit (all four recorded fields,
plus per-vertex loads and states over the first 40 steps) in every graph of
both arms: (E1) mirror image with every color swapped, which sends the trip
permutation to $w_0 \pi w_0$; (E2) negating the colors of the degree-2
vertices, which keeps $\pi$; and (N) negating every color, which sends
$\pi$ to $\pi^{-1}$ and negates every state with loads unchanged. Controls
that do break the orbit: negating only the first vertex of degree $\ge 3$
preserves it in 6,841 of 17,008 bridge graphs and 5,988 Le-graphs, negating
every other one in 4,033 and 3,601 (each count includes the 22 graphs per
arm with no such vertex, which these controls leave unchanged). Off the
canonical families, a recheck script (`sdcheck/check.py` in the scratchpad)
found E1, E2, N and the constant-load characterization to hold on up to 64
colorings each of a handful of small graphs — the square, the stars, and a
non-reduced graph with a loop and a double edge, including colorings with a
zero load; that is a smoke test, not the sweep the follow-up below asks
for. The Le-arm
mechanism is the
isomorphism composed with E1 and E2, so for $n \le 8$ every link of it is
now checked exhaustively — by computation, not proof. Two controls tried
first in the critique session, "mirror without swap" and "swap degree
$\ge 3$ without mirroring", are non-discriminating for the orbit (they are
E1 and E2 composed with N) and change only the trip permutation.

*What "belongs to the bridge graph" means.* By E1, mirroring and
color-swapping the bridge graph of $w_0 \pi w_0$ gives a graph of the cell
$\pi$ with the conjugate's dynamics: a bridge decomposition anchored at
label $n$ instead of label 1 (an identification made by definition — the
library has no such constructor to check against). Its summary row for
$\mathrm{Gr}(k,n)$ is then the recorded row of $\mathrm{Gr}(n-k,n)$; for
example $\mathrm{Gr}(2,8)$ becomes 239: 194 / 45 / 0; 2.83, 18, 85. That
row identity follows from E1 and from conjugation mapping the cells of
$\mathrm{Gr}(k,n)$ one-to-one onto those of $\mathrm{Gr}(n-k,n)$; the
numerical comparison of the 6 pairs run in the critique session looks up
the conjugate's recorded row and cannot fail, so it is not independent
evidence — E1 on all 34,016 graphs is. Granting E1 and the identification,
the direction of Finding 020's asymmetry is set by the anchoring convention
of the decomposition, and reverses with it. The interpretation the data support is that the
Le-construction is $w_0$-equivariant and the bridge construction is not —
not that Le-graphs see a symmetry of the cell: periods are
representative-dependent (below), so no dynamical quantity measured here is
a property of the cell.

Disclosure: this test replaced a weaker probe chosen after seeing the
data. Comparing multisets of (color, degree) pairs, a full color swap
matched 0 of 2,175 Le cells, no swap 169, and a swap restricted to degree
at least 3 matched 2,175 — but that last probe does not discriminate: it
also passes on all 2,175 bridge cells, where the invariance fails. It is
reported only as the path to the isomorphism test, which does discriminate.

Literature disposition: **[consistent]** in part, otherwise
**[no-coverage]**; nothing recorded contradicts it. The Logseq page
*Le-Diagram* records (Williams 2005, section 4) that reflecting a
Le-diagram of type $(k,n)$ over the main diagonal gives a Le-diagram of type
$(n-k,n)$ of the same rank, and (Postnikov section 20) that the Le-graph
trip rule is symmetric about the axis $x + y = 0$ — a $k \leftrightarrow
n-k$ involution of the right shape, though the graph does not record which
decorated permutation the transposed diagram carries. For the bridge arm,
the page *BCFW Recursion* (Arkani-Hamed et al., arXiv:1212.5605) defines the
decomposition by the lexicographically first admissible pair $a < c$, a
choice anchored at label 1 that nothing recorded says commutes with $w_0$.
No page or paper in `docs/ref` records the effect of mirroring or
color-swapping a plabic graph on its trip permutation, nor any statement
about the bridge graph and duality. The checker could not re-verify these
Logseq claims against their source PDFs, which are not in `docs/ref`. On
direction: the `permutation` column is the graph's trip permutation
(guarded per row), and the data separate the candidates — conjugation by
$w_0$ is preserved in every Le cell while inversion is not, so
$\pi \mapsto w_0 \pi^{-1} w_0$ is not a Le-arm symmetry.

It would be falsified by a Le-arm cell at $n = 9$ whose
`(transient, period)` differs from that of its $w_0$-conjugate, or by a
pair of Le-graphs at $n \ge 9$ that fail the mirror-and-swap isomorphism
test ($n = 8$ was tested in the critique session and passed; see the
addendum above).

### Le-arm fixed points are also the non-crossing perfect matchings

On the Le arm the cells with `period == 1` and `transient == 0` are again
exactly the non-crossing perfect matchings: 1, 0, 2, 0, 5, 0, 14 for
$n = 2..8$, with both set differences empty for every $n$, and no period-1
row with a positive transient. The two arms therefore have the same
fixed-point set although their graphs never share a vertex count. The
audited draft read this as suggesting that H2 is about the cell, not the
bridge graph, and named a third representative as the falsifier. **That
reading is withdrawn**: the critique session found the falsifier — the
square graph of the top cell of $\mathrm{Gr}(2,4)$ is a fixed point of a
cell that is not a non-crossing matching (see
[Critique addendum](#critique-addendum-fixed-points-are-constant-load-colorings)).
What the two arms share is that both constructions yield a union of paths
exactly on the non-crossing matchings and a non-constant load everywhere
else.

Literature disposition: **[consistent]** for one direction,
**[no-coverage]** for the converse. The Logseq pages *Decorated
Permutation* and *Positroid* record (Ardila-Rincón-Williams Theorem 7.6,
Proposition 7.8, Corollary 7.9) that a positroid's connected components are
the blocks of the finest non-crossing partition containing every pair
$\{i, \pi(i)\}$; for a non-crossing matching these are the $m$ pairs, and
for a crossing matching the crossing pairs merge. The page *Hypersimplex*
records (Łukowski-Parisi-Williams Propositions 3.15-3.16) that a cell has
dimension $n - c$, with $c$ its number of components, iff every reduced
plabic graph of it is a forest. A non-crossing matching cell has $n = 2m$,
$c = m$ and dimension $m$ — in this run the 14 such cells at $n = 8$ have
dimension 4 — so each of its reduced graphs is a forest whose components
are paths between two boundary vertices. On a path every internal vertex
lies on the same two trips, all loads are equal, and no load strictly
exceeds its neighbors' mean. That inference (the checker's, from recorded
statements; not machine-checked) explains "non-crossing matching implies
fixed point" for any reduced graph. Nothing recorded bears on the converse,
which is purely empirical on the two canonical graphs and false for reduced
graphs in general (the square above).

### Periods depend on the representative

Joining the arms on `(n, permutation)`, the two graphs give the same period
in 7,884 of 17,008 cells (46.35%) and the same `(transient, period)` in
3,271 (19.23%). Agreement on the period falls with $n$: 100% for
$n \le 4$, then 79.55%, 69.81%, 58.47% and 44.28% for $n = 5..8$. Most of
the agreement is the shared period-2 bulk: among the 9,533 cells with a
period above 2 in either arm, the periods agree in 4.29%. A period above 2
occurs in both arms for 3,055 cells, on the bridge arm only for 1,826 and
on the Le arm only for 4,652. The Spearman correlation of the two periods
is 0.421, 0.328, 0.282 and 0.243 for $n = 5..8$.

The Le-arm version of the H1 table (analysis plan step 4), in the same
format; rows $(n, k)$ and $(n, n-k)$ coincide, so each pair is listed once.

| Region (Le arm) | Cells: period-2 / period-$>2$ / period-1; mean, max period, max transient |
| --------------- | -------------------------------------------------------------------------- |
| $\mathrm{Gr}(1,2)$ | 1: 0 / 0 / 1; 1.00, 1, 0 |
| $\mathrm{Gr}(1,n)$, $\mathrm{Gr}(n-1,n)$, $3 \le n \le 8$ | 1: 1 / 0 / 0; 2.00, 2, max transient 0 for $n = 3, 4, 5, 7$ and 1 for $n = 6, 8$ |
| $\mathrm{Gr}(2,4)$ | 7: 5 / 0 / 2; 1.71, 2, 0 |
| $\mathrm{Gr}(2,5)$, $\mathrm{Gr}(3,5)$ | 21: 17 / 4 / 0; 2.62, 8, 4 |
| $\mathrm{Gr}(2,6)$, $\mathrm{Gr}(4,6)$ | 51: 36 / 15 / 0; 3.84, 26, 9 |
| $\mathrm{Gr}(3,6)$ | 161: 121 / 35 / 5; 3.30, 26, 21 |
| $\mathrm{Gr}(2,7)$, $\mathrm{Gr}(5,7)$ | 113: 73 / 40 / 0; 5.51, 42, 35 |
| $\mathrm{Gr}(3,7)$, $\mathrm{Gr}(4,7)$ | 813: 528 / 285 / 0; 5.51, 84, 120 |
| $\mathrm{Gr}(2,8)$, $\mathrm{Gr}(6,8)$ | 239: 127 / 112 / 0; 8.74, 188, 58 |
| $\mathrm{Gr}(3,8)$, $\mathrm{Gr}(5,8)$ | 3,361: 1,748 / 1,613 / 0; 11.90, 618, 503 |
| $\mathrm{Gr}(4,8)$ | 7,631: 4,083 / 3,534 / 14; 13.61, 1,003, 1,195 |

Even where the periods agree the transients need not: the single cell of
$\mathrm{Gr}(1,n)$ and of $\mathrm{Gr}(n-1,n)$ has transient 1 on the Le
arm at $n = 6$ and $n = 8$ against 0 on the bridge arm, and
$\mathrm{Gr}(2,4)$ has maximum transient 0 on the Le arm against 2.

The period is therefore not a cell invariant under this rule: the summary
statistics of either arm describe the (cell, graph) pair. The falling
agreement over four values of $n$ is a trend on few points, not a fitted
law.

Literature disposition: **[consistent]**. The Logseq page *Plabic Graph*
records (Postnikov Lemma 13.1, Theorem 13.4) that reduced graphs are
move-equivalent iff they share the decorated trip permutation, and every
recorded invariant of a move class is static (trip permutation, positroid,
dimension, face count). The moves themselves change the automaton's state
space — the square move switches four colors, (M2) and (M3) add or remove
vertices — so no recorded statement would have led one to expect a
representative-independent period.

### The Le arm is slower on average; the bridge arm has the longer extremes

Le-graphs have strictly fewer internal vertices than bridge graphs in all
17,008 cells (at $n = 8$: 4 to 25 against 8 to 32), yet the mean period is
higher on the Le arm for every $n \ge 5$ — 12.674 against 6.965 at
$n = 8$ — and the $n = 8$ period quantiles at 0.5 / 0.9 / 0.99 / 0.999 are
2 / 26 / 176 / 450 (Le) against 2 / 8 / 101.36 / 379.06 (bridge). The
extremes go the other way: maximum period 2,654 (bridge) against 1,003
(Le), maximum transient 1,509 against 1,195.

At $n = 8$ — the only $n$ examined for this — the period is
rank-correlated with the size of the graph within each arm: the Spearman
correlation of `period` with `internal_vertices` is 0.425 (bridge) and
0.507 (Le). The relation is not monotone. The median period is 1 at dimension 4
(the 14 non-crossing matching cells) and flat at 2 over dimensions 5 to 9
on the bridge arm and 5 to 8 on the Le arm, then
rises to 122 at dimension 15 (bridge, 10 cells) and 173 at dimension 14
(Le, 52 cells), and falls back at the top: the Le median is 38 at dimension
15, and by vertex count it is 182 at 21 vertices, then 126, 36, 40 and 18
at 22 to 25. The unique top-dimensional cell of $\mathrm{Gr}(4,8)$
(dimension 16) has period 2 on its bridge graph and 18 on its Le-graph. *Confound:* on the bridge arm
`internal_vertices` equals twice `dimension` in all 17,008 cells, so the
two axes cannot be separated there; on the Le arm up to 4 vertex counts
share a dimension, and at $n = 8$ the period's Spearman correlation with
`dimension` is 0.516.

*Frozen vertices (critique addendum).* 2,310 of the 17,008 Le-graphs have
at least one internal vertex with no internal neighbor, which the rule
never flips; no bridge graph has one. On those Le-graphs the automaton's
effective state space is smaller than `internal_vertices` indicates, which
the vertex-count comparisons above do not correct for. How the Le-arm
statistics split between graphs with and without a frozen vertex was not
computed.

Bug-first triage of "fewer vertices, longer periods". Literature
disposition: **[no-coverage]** — the rule is original, and the recorded
invariants of a move class (see the previous observation) include neither
vertex counts nor anything dynamical, so nothing predicts which graph
should be slower. Because a defect confined to the Le arm would produce
exactly this pattern and the Le arm has no oracle, the bug classes were
checked: (i) *wrong cell* — every row passed the experiment's guard that
the built graph's trip permutation equals `permutation`; (ii) *wrong
dynamics on Le-graphs* — a naive re-implementation of the rule, which
traces trips with the graph library's own `PlabicGraph.trips()` on the
recolored graph instead of the automaton's dart tables, reproduces the
run's `(transient, period)` in 321 of 321 cells on each arm for
$n = 2..6$ (the arms first differ at $n = 4$ in transient and $n = 5$ in
period); (iii) *drift since sizing* — the Le-arm maxima for $n \le 7$
equal those disclosed in the design doc. No bug found; $n = 7, 8$ were not
re-derived independently.

Outliers: the bridge maximum period 2,654 is the cell `5,6,7,8,2,1,3,4`
($k = 4$, dimension 15, 30 vertices, transient 333, 1,962 distinct trip
permutations on its cycle); the bridge maximum transient 1,509 is
`5,6,7,8,1,3,2,4`, which then settles to period 2. The Le maximum period
1,003 is attained by `5,6,7,8,1,2,4,3` and its $w_0$-conjugate
`6,5,7,8,1,2,3,4` (transient 865), the Le maximum transient 1,195 by
`5,6,7,8,2,1,4,3` and its conjugate.

### Trip permutations almost always move, and with the full period

`period / permutation_period` is 1 in 16,991 bridge cells and 16,966 Le
cells, and 2 in the remaining 17 and 42; no other ratio occurs. So for
$n \le 8$ the trip-permutation sequence has the coloring's period except in
59 rows where it has half of it. A cycle holds a single trip permutation in
36 bridge cells and 63 Le cells; excluding fixed points, in 14 of 16,986
and 41 of 16,986. The number of distinct trip permutations on the cycle
equals `permutation_period` in 94.71% (bridge) and 87.84% (Le) of cells,
and reaches 1,962 and 856. Odd periods above 1 occur in 159 bridge cells
and 859 Le cells, so the period-2 bulk does not reflect a global parity
constraint.

The trip permutation therefore changes along the cycle of almost every
non-fixed trajectory, in both arms. Whether the trajectory leaves the
positroid cell is a stronger statement the data cannot carry: nothing
checks that the recolored graphs are reduced, and the trip permutation of a
non-reduced graph does not determine a positroid. Literature disposition: **[no-coverage]** — the
checker confirmed that no Logseq page or block describes this rule
("strand dynamics", "strand gradient", "plabic automaton"), as the design
doc states.

## Suspected bugs

One, in the metadata only. `src/experiments/strand_dynamics.py`, as of the
run's commit `d4312c9`, evaluates
`_git_commit()` while building the metadata mapping, after
`ctx.write_result(frame)` has created an untracked file under
`data/results/`; `git status --porcelain` then reports that file and the
commit is suffixed `-dirty` on every run from a clean tree. Evidence,
checkable in the source at `d4312c9` (the line numbers are that commit's;
the fix below moved them): `ctx.write_result(frame)` (line 198) precedes
`_git_commit()` (line 212), the dirtiness test is a plain
`git status --porcelain`, and `data/results/` is not ignored. In-session,
`git status` was seen clean before the run and listing only the run's two
artifacts directly after it; the artifacts cannot verify that observation.
The defect makes every run look irreproducible; it does
not touch the data. Fix: capture the commit before writing, or ignore
untracked files (`--untracked-files=no`) in the dirtiness test.

*Fixed after this run* (critique session, uncommitted at the time of
writing): the experiment now reads the commit before the sweep starts, so
untracked source files still count as dirty but the run's own artifacts do
not; a regression test pins the ordering. This run's sidecar keeps its
`-dirty` suffix.

No suspected defect in the dynamics, the graph constructions, or the
schema: every row passed the experiment's own trip-permutation guard, and
the bridge arm matched its oracle.

## Follow-up experiments

- **$n = 9$, both arms** (`--n-max 9`) — a replication chore, ranked last
  after the critique session. It completes H1's six remaining rows, but H1
  is now settled cell-by-cell through $n = 8$ and nothing suggests the port
  diverges at $n = 9$; the Le-arm $w_0$ invariance has a mechanism checked
  exhaustively for $n \le 8$, so one more $n$ is weak evidence next to
  proving the lemmas or testing the rule's symmetries off the canonical
  families; and the absence of fixed points at odd $n$ cannot fail if the
  constant-load reading holds and odd-$n$ canonical graphs always contain a
  vertex of degree at least 3. A Le cell whose `(transient, period)`
  differs from its conjugate's would still kill the invariance.
- **The Le-arm invariance as lemmas** (the $n = 8$ half of this follow-up
  was run in the critique session: the isomorphism holds for all 17,008
  cells with $n \le 8$, and the rule's symmetries E1, E2 and N hold on every
  graph of the run). What remains is proof, not more sweeping: (a) the
  Le-graph of $w_0 \pi w_0$ is the mirror image of the Le-graph of $\pi$
  with the colors of degree $\ge 3$ swapped — a statement about
  Le-diagrams, where Williams's transposition is the natural starting
  point; (b) the rule's orbit is invariant under mirror image with full
  color swap, under global negation, and under recoloring degree-2
  vertices — statements about the rule on any plabic graph, each killed by
  one graph, of any kind, on which the orbit changes. Testing (b) on graphs
  outside the two canonical families (random planar graphs, non-reduced
  ones) has more discriminating power than $n = 9$.
- **Constant load within a move class** (replaces "a third
  representative", whose fixed-point half the critique session resolved:
  the square graph of `3,4,1,2` is a fixed point, so fixed points are not a
  cell property, and the non-crossing matching half cannot fail on any
  representative that is a union of paths). Enumerate
  the graphs reachable from the canonical ones by square moves and
  (un)contractions and record which have constant load on every component.
  The forming guess — they are the representatives with a symmetry the
  canonical constructions lack — is killed by one asymmetric constant-load
  graph with a vertex of degree at least 3. The same sweep measures how far
  the period moves within a move-equivalence class.
- **$\mathrm{Gr}(5,10)$ restricted to involutions.** H2 predicts exactly 42
  fixed points among the fixed-point-free involutions of $[10]$, on both
  arms; the restriction keeps the sweep to 945 cells per arm. Rescoped in
  the critique session: the 42 non-crossing matchings cannot fail to be
  fixed if their canonical graphs are unions of paths (true of all 22 such
  graphs with $n \le 8$; extrapolated to $n = 10$), so the main live question is
  whether any of the 903 crossing matchings has constant canonical load.
  The restriction also leaves H2's other half untested at $n = 10$ — a
  fixed point that is not an involution — which needs the full sweep.
- **Cells with fixed points.** The design doc argues lollipops are inert;
  a sweep over all decorated permutations would test that the statistics
  of a cell with fixed points equal those of the cell with them removed.

## Conjecture links

None yet — `docs/conj/` has no pages. Two statements from this run are
candidate seeds. Apart from the bridge half of the first, which is the
pre-registered H2, both come from regions with no pre-registered prediction
or from post hoc analysis. First, the fixed-point characterization, rescoped by the critique
session: a lemma (a coloring is fixed iff its load is constant on each
component of internal vertices) plus a construction-level statement (on
bridge and Le graphs the canonical load is constant on every component iff
the graph is a union of paths, iff the cell is a non-crossing perfect
matching). It is not a statement about cells — see the
[critique addendum](#critique-addendum-fixed-points-are-constant-load-colorings)
— and whether a statement about two particular constructions deserves a
conjecture page, as opposed to the lemma alone, is open. Second, the Le-arm
invariance of the orbit shape under $\pi \mapsto w_0 \pi w_0$, which the
critique session split into two provable-looking pieces: the Le-graph
mirror statement, and the rule's three orbit symmetries (mirror image with
full color swap, global negation, degree-2 recoloring), each checked on
every graph of this run.
