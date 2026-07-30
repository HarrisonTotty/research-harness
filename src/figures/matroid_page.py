"""The knowledge graph's Matroid page, as slide graphics.

Renders `docs/conj`-adjacent knowledge-graph content rather than a process:
the Logseq page **Matroid** is the graph's largest mathematical entry, and
this module turns it into two figures sized for a presentation slide.

The first figure, ``matroid-page``, is the page's own argument. Its centre is
the fact the page is organised around — a matroid has seven equivalent
definitions, and choosing between them is the first design decision any
library makes. That is drawn as a hub and ring: independence sits in the
middle because it is the primitive both stacks actually store, and the six
other axiom systems ring it, each carrying its own numbered axioms. Arrows
show only *which* conversions the page records; the formulas themselves live
in the conversion-table band below the wheel, where there is room to read
them. Three rim arrows close the ring where the table's conversions do not pass
through independence at all — rank and the flats both fix the closure, and the
hyperplanes are the cocircuit complements.

Around the wheel sit the page's other sections, in the positions their content
warrants: the vocabulary and operations the definitions generate on the left,
the structural theorems and canonical fixtures on the right — those two are
the page's property-test oracle and its test data, which is why they are
adjacent. The band along the bottom is what the page was written for: one
Python module and one Lean module, each reporting what it reused, what it had
to build, and what it left behind.

The second figure, ``matroid-page-links``, is the page's neighbourhood in the
graph — 64 forward links clustered by the role they play on the page, and the
6 pages that cite Matroid back. The clusters are what makes the count legible:
a flat list of 64 titles says nothing, whereas the split shows that a third of
them are classes and fixtures and only seven are people or tooling. Backlink
arrows run inward, so the page reads as a hub for the positroid/Grassmannian
line of work rather than a leaf.

Both graphics carry no title chrome beyond the page's own definition — the
slide supplies the rest — and are saved cropped to the drawing itself, so they
can be placed and scaled freely.

Besides the complete breakdown, the module writes a build-up sequence — one
image per reveal step, each adding the next layer — for walking an audience
through it. Every step is rendered into the same frame as the whole, so the
images can be stacked on one slide without anything shifting between them.

Regenerate with ``just figure matroid-page``.
"""

import math
from dataclasses import dataclass
from enum import Enum, auto

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from figures import style
from figures.cli import FigureContext, figure

_PAD: float = 0.2
"""Border kept around the drawing on every side, in inches."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_MARGIN_L: float = 0.5
"""Left edge of every full-width element."""

_MARGIN_R: float = 15.5
"""Right edge of every full-width element."""

# --------------------------------------------------------------------------
# Breakdown figure: vertical bands, bottom-up.
# --------------------------------------------------------------------------

_FOOTER_BOTTOM: float = 0.35
"""Bottom edge of the open-questions strip, and of the whole drawing."""

_FOOTER_TOP: float = 0.88
"""Top edge of the open-questions strip."""

_IMPL_BOTTOM: float = 1.06
"""Bottom edge of the two implementation panels."""

_IMPL_TOP: float = 2.78
"""Top edge of the two implementation panels."""

_TABLE_BOTTOM: float = 2.90
"""Bottom edge of the cryptomorphism conversion table."""

_TABLE_TOP: float = 4.05
"""Top edge of the cryptomorphism conversion table."""

_MAIN_BOTTOM: float = 4.18
"""Bottom edge of the wheel and its flanking section panels."""

_MAIN_TOP: float = 8.73
"""Top edge of the wheel and its flanking section panels."""

_HEADER_BOTTOM: float = 8.91
"""Bottom edge of the header band."""

_HEADER_TOP: float = 9.72
"""Top edge of the header band, and of the whole drawing."""

_COL_W: float = 3.85
"""Width of the section panel columns flanking the wheel."""

_COL_L: float = _MARGIN_L
"""Left edge of the left-hand section column."""

_COL_R: float = _MARGIN_R - _COL_W
"""Left edge of the right-hand section column."""

_PANEL_GAP: float = 0.25
"""Vertical gap between the two panels stacked in a section column."""

_UPPER_H: float = 2.15
"""Height of the upper panel in each section column."""

_LOWER_H: float = _MAIN_TOP - _MAIN_BOTTOM - _PANEL_GAP - _UPPER_H
"""Height of the lower panel in each section column."""

_UPPER_BOTTOM: float = _MAIN_TOP - _UPPER_H
"""Bottom edge of the upper panel in each section column."""

_PANEL_PAD: float = 0.12
"""Inset from a panel's edge to its text, in inches."""

_LINE_STEP: float = 0.135
"""Vertical distance between consecutive body lines inside a panel."""

_PANEL_FIRST_LINE: float = 0.36
"""Drop from a panel's top edge to its first body line, clearing the title."""

_WHEEL_CX: float = 8.0
"""Horizontal centre of the cryptomorphism wheel."""

_WHEEL_CY: float = (_MAIN_BOTTOM + _MAIN_TOP) / 2.0 - 0.12
"""Vertical centre of the cryptomorphism wheel, nudged below the band's own
centre to open a strip along the top for the wheel's caption."""

_RING_RX: float = 2.95
"""Horizontal radius of the ring the six derived systems sit on. The ring is
an ellipse rather than a circle because the drawing is wider than it is tall,
and a circle would waste the horizontal room the chips need."""

_RING_RY: float = 1.62
"""Vertical radius of the ring the six derived systems sit on."""

_CHIP_W: float = 1.80
"""Width of a ring chip."""

_CHIP_H: float = 0.80
"""Height of a ring chip; it carries a badge line and two axiom lines."""

_HUB_W: float = 2.40
"""Width of the independence hub at the wheel's centre."""

_HUB_H: float = 1.30
"""Height of the independence hub."""

_EDGE_GAP: float = 0.05
"""Clearance left between an arrowhead and the box it points at, so the head
never lands on a rounded corner and read as part of the outline."""

_SPOKE_COLOR: str = style.CATEGORICAL[0]
"""Slate — carries the conversions that pass through independence, which is
every conversion the page records except the three on the rim."""

_RIM_COLOR: str = style.CATEGORICAL[4]
"""Purple — carries the three conversions that bypass independence entirely
(rank to closure, flats to closure, hyperplanes to cocircuits), so the rim
reads as a different kind of edge from the spokes."""

_RIM_BOW: float = 0.3
"""How far a rim conversion bows off its chord, in chord lengths."""

_BUILT_COLOR: str = style.CATEGORICAL[1]
"""Teal — marks what the two implementation stacks actually built."""


def _radians(degrees: float) -> float:
    """Return ``degrees`` in radians."""
    return degrees * math.pi / 180.0


def _box_exit(half_w: float, half_h: float, degrees: float) -> tuple[float, float]:
    """Return where a ray leaves an axis-aligned box centred on the origin.

    Args:
        half_w: Half the box's width.
        half_h: Half the box's height.
        degrees: Direction of the ray, measured counter-clockwise from east.

    Returns:
        The exit point as an offset from the box's centre, so a caller adds it
        to the centre to get canvas coordinates. Used to start and end every
        spoke on a box edge rather than at its centre, which keeps arrowheads
        off the boxes at any ring angle.
    """
    angle = _radians(degrees)
    dx, dy = math.cos(angle), math.sin(angle)
    scale = min(
        half_w / abs(dx) if dx else math.inf,
        half_h / abs(dy) if dy else math.inf,
    )
    return scale * dx, scale * dy


def _edge_endpoints(
    source: tuple[float, float],
    target: tuple[float, float],
    source_half: tuple[float, float],
    target_half: tuple[float, float],
    rad: float = 0.0,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return where an arrow between two boxes meets each box's edge.

    Args:
        source: Centre of the box the arrow leaves.
        target: Centre of the box the arrow points at.
        source_half: Half-width and half-height of the source box.
        target_half: Half-width and half-height of the target box.
        rad: Bow of the arc, matching the ``rad`` of :func:`_arrow`.

    Returns:
        The start and end points, as canvas coordinates.

    Both ends are taken along the direction the drawn curve actually travels,
    not along a ring angle: chips sit on an ellipse, so the bearing from the
    hub to a chip is not the angle that placed it, and reusing the placement
    angle lands arrowheads on box corners instead of edges.
    """
    chord = math.degrees(math.atan2(target[1] - source[1], target[0] - source[0]))
    # `arc3` puts its control point one `rad`-scaled chord-length off the
    # midpoint, along the chord's *clockwise* normal, so the curve leaves each
    # end `atan(2 * rad)` off the chord in that direction.
    skew = math.degrees(math.atan(2.0 * rad))
    source_dx, source_dy = _box_exit(source_half[0], source_half[1], chord - skew)
    target_dx, target_dy = _box_exit(
        target_half[0], target_half[1], chord + 180.0 + skew
    )
    return (
        (source[0] + source_dx, source[1] + source_dy),
        (target[0] + target_dx, target[1] + target_dy),
    )


class _Spoke(Enum):
    """How a system's chip is joined to the independence hub.

    The page's conversion table is not symmetric. Most systems convert both to
    and from independence, but the closure is only ever converted *to* it — the
    table reaches the closure itself from the rank and from the flats, never
    from independence — so its spoke carries a single arrowhead.
    """

    NONE = auto()
    INWARD = auto()
    BOTH = auto()


@dataclass(frozen=True, slots=True)
class _System:
    """One of the six axiom systems ringing the independence hub.

    ``angle`` places the chip on the ring, in degrees counter-clockwise from
    east. The order around the ring is chosen so that every rim conversion
    connects neighbours; see :data:`_SYSTEMS`.
    """

    symbol: str
    name: str
    axioms: tuple[str, str]
    angle: float
    spoke: _Spoke


_SYSTEMS: tuple[_System, ...] = (
    _System(
        symbol=r"$\mathcal{B}$",
        name="BASES",
        axioms=("(B1) at least one basis", "(B2) basis exchange"),
        angle=90.0,
        spoke=_Spoke.BOTH,
    ),
    _System(
        symbol=r"$\mathcal{C}$",
        name="CIRCUITS",
        axioms=("(C1) no empty circuit", "(C2) antichain (C3) elim."),
        angle=30.0,
        spoke=_Spoke.BOTH,
    ),
    _System(
        symbol=r"$\mathcal{H}$",
        name="HYPERPLANES",
        axioms=("(H1) E is not one", "(H2) antichain (H3) elim."),
        angle=330.0,
        spoke=_Spoke.NONE,
    ),
    _System(
        symbol=r"$\mathcal{F}$",
        name="FLATS",
        axioms=("(F1) E is a flat", "(F2) meet-closed (F3) cover"),
        angle=270.0,
        spoke=_Spoke.NONE,
    ),
    _System(
        symbol=r"$\mathrm{cl}$",
        name="CLOSURE",
        axioms=("(CL1-3) closure operator", "(CL4) Mac Lane-Steinitz"),
        angle=210.0,
        spoke=_Spoke.INWARD,
    ),
    _System(
        symbol=r"$r$",
        name="RANK",
        axioms=("(R1) r(X) <= |X| (R2) mono", "(R3) submodular"),
        angle=150.0,
        spoke=_Spoke.BOTH,
    ),
)
"""The six derived systems, clockwise from the top.

The ring order is load-bearing: bases, circuits, hyperplanes, flats, closure,
rank puts rank next to closure and hyperplanes next to circuits, which are
exactly the two conversions that do not pass through independence. Both rim
arrows therefore join adjacent chips instead of cutting across the wheel.
"""


@dataclass(frozen=True, slots=True)
class _Rim:
    """A conversion drawn on the rim, between two neighbouring chips.

    The caption is anchored explicitly rather than derived from the arrow's
    midpoint: the three rim edges sit in three differently shaped pockets of
    free space, and only a hand-placed label clears both the chips and the
    section columns on every one of them.
    """

    source: str
    target: str
    tag: tuple[str, str]
    label: tuple[float, float]
    ha: str


_RIM: tuple[_Rim, ...] = (
    _Rim("RANK", "CLOSURE", ("closure", "from rank"), (5.78, 6.42), "left"),
    _Rim("FLATS", "CLOSURE", ("closure =", "meet of flats"), (6.62, 4.95), "center"),
    _Rim(
        "HYPERPLANES",
        "CIRCUITS",
        ("cocircuits are", "complements"),
        (10.22, 6.42),
        "right",
    ),
)
"""Conversions drawn on the rim, inside the ring.

The tags are deliberately short — the formulas are in the conversion table,
which is the only place on the drawing wide enough to hold them.
"""

_HUB_LINES: tuple[str, ...] = (
    "(I1) the empty set is independent",
    "(I2) hereditary      - Whitney (a)",
    "(I3) augmentation    - Whitney (b)",
    "",
    "the stored primitive - but note E is",
    "not recoverable from the family alone",
)
"""Body of the independence hub, below its badge line."""

_HUB_AXIOM_LINES: int = 3
"""How many of :data:`_HUB_LINES` are the numbered axioms; those are set in
ink, and the commentary below them is set back in the muted grey."""

_CONVERSIONS: tuple[tuple[str, str], ...] = (
    (
        r"$\mathcal{I} \rightarrow \mathcal{B}$",
        "bases are the maximal independent sets",
    ),
    (
        r"$\mathcal{I} \rightarrow \mathcal{C}$",
        "circuits are the minimal dependent sets",
    ),
    (
        r"$\mathcal{I} \rightarrow r$",
        "r(X) = size of a largest independent subset",
    ),
    (
        r"$\mathcal{B} \rightarrow \mathcal{I}$",
        "independent = contained in some basis",
    ),
    (
        r"$\mathcal{C} \rightarrow \mathcal{I}$",
        "independent = contains no circuit",
    ),
    (
        r"$r \rightarrow \mathcal{I}$",
        "independent iff r(X) equals |X|",
    ),
    (
        r"$r \rightarrow \mathrm{cl}$",
        "cl(X) = every e with r(X + e) = r(X)",
    ),
    (
        r"$\mathrm{cl} \rightarrow \mathcal{I}$",
        "independent iff no e lies in cl(X - e)",
    ),
    (
        r"$\mathcal{F} \rightarrow \mathrm{cl}$",
        "cl(X) = the meet of all flats above X",
    ),
)
"""The page's nine conversion rows, read down each column of three in turn."""

_CONVERSION_NOTE: str = (
    "and two more by duality:   H -> C*  cocircuits are the hyperplane complements   "
    "*   B -> B*  the dual's bases are the complements of the bases"
)
"""Trailing row of the conversion table: the two duality-flavoured rows."""

_CONVERSION_RULE: str = (
    "store exactly one primitive and derive the other six lazily - keeping several "
    "in sync is the correctness hazard; round-trip every pair as a regression test"
)
"""The design consequence the page draws from the table, and the reason the
table is on the slide at all."""


@dataclass(frozen=True, slots=True)
class _Panel:
    """A titled list of body lines occupying one rectangle of the layout."""

    title: str
    lines: tuple[str, ...]


_VOCABULARY: _Panel = _Panel(
    title="DERIVED VOCABULARY",
    lines=(
        "dependent    not independent",
        "spanning     r(X) = r(E)",
        "basis        independent and spanning",
        "circuit      minimal dependent",
        "cocircuit    circuit of M* = E - a hyperplane",
        "loop         r({e}) = 0; lies in no basis",
        "coloop       lies in every basis; loop of M*",
        "parallel     non-loops with r({e,f}) = 1",
        "simple       no loops, no parallel pair",
        "nullity      n(X) = |X| - r(X)",
        "corank       r*(X) = |X| + r(E-X) - r(E)",
        "connected    every pair shares a circuit",
        "C(e,B)       the unique circuit in B + e",
    ),
)
"""Left column, upper panel: what the seven definitions name."""

_OPERATIONS: _Panel = _Panel(
    title="OPERATIONS AND CONSTRUCTIONS",
    lines=(
        "restriction   M|X - the independent sets inside X",
        "deletion      M\\X = M|(E-X)",
        "contraction   M/X - rank r(Y + X) - r(X)",
        "minor         M\\D/C - deletion and contraction commute",
        "duality       the dual's bases are the complements;",
        "              swaps \\ with /, loop/coloop, C with C*",
        "direct sum    M1 (+) M2 on disjoint ground sets",
        "truncation    T_k keeps independent sets of size <= k",
        "extensions    free, principal; series and parallel",
        "              connection; the 2-sum",
        "union, inter. poly-time given an independence oracle -",
        "              but 3-matroid intersection is NP-hard",
    ),
)
"""Left column, lower panel: what the definitions let you build."""

_THEOREMS: _Panel = _Panel(
    title="STRUCTURAL THEOREMS - THE PROPERTY-TEST ORACLE",
    lines=(
        "Rado-Edmonds  a hereditary system is a matroid iff",
        "              greedy is optimal for every weight fn.",
        "duality       M** = M, an involution",
        "orthogonality |C and C* meet| is never exactly 1",
        "fundamental   I + e dependent => exactly one circuit",
        "bases         equicardinal; exchange graph connected",
        "Edmonds       max |I| = min r1(X) + r2(E - X)",
        "Nash-Williams union rank; base packing and covering",
        "              as min-max certificates",
        "excluded min. binary iff no U(2,4) minor; regular",
        "              also excludes F7 and its dual",
        "Tutte T(x,y)  universal deletion-contraction invariant",
        "AHK 2018      characteristic polynomial log-concave",
    ),
)
"""Right column, upper panel: the oracle every implementation tests against."""

_FIXTURES: _Panel = _Panel(
    title="CANONICAL EXAMPLES - THE TEST FIXTURES",
    lines=(
        "U(r,n)      independent iff |I| <= r;",
        "            duality law U(r,n)* = U(n-r,n)",
        "U(n,n) free           U(0,n) loopy",
        "U(2,4)      smallest non-binary; self-dual",
        "M(K4)       graphic, rank 3 on 6 elements",
        "F7   Fano   representable iff characteristic 2",
        "F7-  non-Fano         iff characteristic not 2",
        "V8   Vamos  no field at all - it violates the",
        "            Ingleton inequality",
        "R10, Pappus, non-Pappus  round out the set",
        "M[A]        transversal; a matching oracle",
        "count       1,2,4,8,17,38,98,306,1724 (A055545)",
    ),
)
"""Right column, lower panel: the small matroids the tests run on."""

_PYTHON: _Panel = _Panel(
    title="PYTHON  -  src/research/matroid.py",
    lines=(
        "frozen dataclass Matroid[T] - an ordered ground set plus the "
        "independence family as int bitmasks;",
        "every other presentation is derived lazily over an O(2^n * n) rank "
        "table, so n is about 16 in practice.",
        "axiom ctors  from_independent_sets / bases / circuits / rank_function "
        "/ closure / flats / hyperplanes;",
        "             each failure raises ValueError naming the numbered axiom "
        "from this page",
        "structural   from_vectors (over Q or GF(p) - the field is always an "
        "explicit parameter),",
        "             from_graph_edges, from_transversal_system, from_dataframe",
        "fixtures     uniform, free, loopy, empty, u24, k4, fano, non_fano, "
        "vamos; enumerate_matroids reproduces",
        "             the A055545 prefix to n = 5.  Also plot_lattice_of_flats "
        "and plot_basis_exchange_graph",
        "tests        107 green under just check - every theorem on the left "
        "mapped to a property test",
    ),
)
"""Implementation band, left: what the Python stack stores and checks."""

_LEAN: _Panel = _Panel(
    title="LEAN  -  src/theorems/Theorems/Matroid.lean",
    lines=(
        "extends Mathlib's base-primitive Matroid (E, IsBase, Indep) rather "
        "than redefining it; rank is",
        "N-infinity valued, so infinite ground sets stay in scope. "
        "Kernel-checked, zero sorry.",
        "reused       about 30 page claims - every axiom system as a theorem, "
        "base equicardinality, strong",
        "             circuit elimination, orthogonality, the duality / minor "
        "/ direct-sum calculus,",
        "             and 7 of the 9 conversion rows",
        "added        IsHyperplane (= the cocircuit complements), Parallel, "
        "Simple, nullity, truncateTo,",
        "             the contraction and coloop rank identities, "
        "IsFlat.inter, and the uniform-matroid API",
        "backlog      Tutte polynomial (wants an N-valued rank), union and "
        "intersection, connectivity,",
        "             representability, graphic and transversal matroids, the "
        "F7 / V8 / R10 fixtures",
    ),
)
"""Implementation band, right: what the Lean stack reused, added, and deferred."""

_OPEN_QUESTIONS: tuple[tuple[str, str], ...] = (
    (
        "finite-only, or finite plus finitary-infinite?  Python",
        "wants finite; Lean inherits Mathlib's infinite definition",
    ),
    (
        "which primitive is canonical in Python - the independence",
        "oracle or the rank oracle?  It is baked into every signature",
    ),
    (
        "is a decidable Lean mirror in scope, or does Lean stay",
        "purely proof-side with all computation living in Python?",
    ),
)
"""The three questions the page says to settle before writing more code."""

_STEP_HUB: int = 1
"""Reveal step that states the definition."""

_STEP_WHEEL: int = 2
"""Reveal step that adds the six other axiom systems."""

_STEP_TABLE: int = 3
"""Reveal step that adds the conversion table."""

_STEP_LEFT: int = 4
"""Reveal step that adds the vocabulary and operations the definitions generate."""

_STEP_RIGHT: int = 5
"""Reveal step that adds the theorem oracle and the fixtures."""

_STEP_IMPL: int = 6
"""Reveal step that adds the two implementation stacks."""

_STEP_OPEN: int = 7
"""Reveal step that adds the open questions; the last one is the whole figure."""

_TOTAL_STEPS: int = _STEP_OPEN
"""Number of images in the build-up sequence."""


def _box(
    ax: Axes,
    rect: tuple[float, float, float, float],
    *,
    facecolor: str,
    edgecolor: str = "none",
    linewidth: float = 1.0,
    linestyle: str = "solid",
) -> None:
    """Draw one rounded panel from ``rect`` as ``(x, y, width, height)``.

    The position is the panel's lower-left corner, in canvas inches.
    """
    x, y, width, height = rect
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle=f"round,pad=0,rounding_size={_CORNER}",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            linestyle=linestyle,
        )
    )


def _arrow(
    ax: Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    rad: float = 0.0,
    linewidth: float = 1.2,
    two_way: bool = False,
) -> None:
    """Draw an arrow from ``start`` to ``end``, bowed by ``rad``.

    Set ``two_way`` for a conversion that runs in both directions, which the
    page's table records for every spoke of the wheel.
    """
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="<|-|>" if two_way else "-|>",
            mutation_scale=11,
            linewidth=linewidth,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=0,
            shrinkB=0,
        )
    )


def _panel(ax: Axes, rect: tuple[float, float, float, float], panel: _Panel) -> None:
    """Draw a titled list panel into ``rect`` as ``(x, y, width, height)``.

    Body lines are pre-wrapped in the panel's own definition rather than
    flowed: at this size the wrap points are a layout decision, and the
    left-hand keyword column only lines up because they are fixed.
    """
    x, y, _, height = rect
    _box(
        ax,
        rect,
        facecolor=style.PAPER,
        edgecolor=style.MIST,
        linewidth=1.0,
    )
    ax.text(
        x + _PANEL_PAD,
        y + height - 0.17,
        panel.title,
        fontsize=8.0,
        fontweight="bold",
        color=style.SLATE,
        va="center",
    )
    for index, line in enumerate(panel.lines):
        ax.text(
            x + _PANEL_PAD,
            y + height - _PANEL_FIRST_LINE - index * _LINE_STEP,
            line,
            fontsize=7.0,
            color=style.INK,
            va="center",
        )


def _header(ax: Axes) -> None:
    """Draw the page's identity and its one-sentence definition."""
    ax.text(
        _MARGIN_L,
        _HEADER_BOTTOM + 0.50,
        "Matroid",
        fontsize=23,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    ax.text(
        _MARGIN_L,
        _HEADER_BOTTOM + 0.14,
        "a Logseq page, broken down",
        fontsize=8.0,
        color=style.MIST,
        va="center",
    )
    ax.plot(
        [3.15, 3.15],
        [_HEADER_BOTTOM + 0.04, _HEADER_TOP - 0.04],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )
    definition = (
        "a combinatorial structure that abstracts independence - exactly the "
        "properties of dependence shared by the columns of a matrix,",
        "the forests of a graph, algebraic independence in a field extension, "
        "and the partial transversals of a set system",
    )
    for index, line in enumerate(definition):
        ax.text(
            3.40,
            _HEADER_TOP - 0.20 - index * 0.24,
            line,
            fontsize=9.0,
            color=style.INK,
            va="center",
        )
    ax.text(
        3.40,
        _HEADER_BOTTOM + 0.14,
        "Whitney 1935  *  Nakasawa 1935-36, independently  *  Oxley, Matroid "
        "Theory 2nd ed. 2011  *  7 cryptomorphic axiom systems  "
        f"*  {_FORWARD_COUNT} forward links",
        fontsize=7.5,
        color=style.MIST,
        va="center",
    )


def _hub(ax: Axes) -> None:
    """Draw the independence hub at the centre of the wheel."""
    left = _WHEEL_CX - _HUB_W / 2.0
    bottom = _WHEEL_CY - _HUB_H / 2.0
    _box(
        ax,
        (left, bottom, _HUB_W, _HUB_H),
        facecolor=style.PAPER,
        edgecolor=style.SLATE,
        linewidth=1.6,
    )
    ax.text(
        left + 0.13,
        bottom + _HUB_H - 0.20,
        r"$\mathcal{I}$",
        fontsize=13,
        color=style.SLATE,
        va="center",
    )
    ax.text(
        left + 0.36,
        bottom + _HUB_H - 0.21,
        "INDEPENDENT SETS",
        fontsize=9.0,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    for index, line in enumerate(_HUB_LINES):
        ax.text(
            left + 0.13,
            bottom + _HUB_H - 0.44 - index * 0.155,
            line,
            fontsize=6.8,
            color=style.INK if index < _HUB_AXIOM_LINES else style.MIST,
            va="center",
        )


def _system_center(system: _System) -> tuple[float, float]:
    """Return the canvas centre of ``system``'s chip on the ring."""
    angle = _radians(system.angle)
    return (
        _WHEEL_CX + _RING_RX * math.cos(angle),
        _WHEEL_CY + _RING_RY * math.sin(angle),
    )


def _chip(ax: Axes, system: _System) -> None:
    """Draw one axis-system chip on the ring, with its numbered axioms."""
    cx, cy = _system_center(system)
    left = cx - _CHIP_W / 2.0
    top = cy + _CHIP_H / 2.0
    _box(
        ax,
        (left, cy - _CHIP_H / 2.0, _CHIP_W, _CHIP_H),
        facecolor=style.PAPER,
        edgecolor=style.SLATE,
        linewidth=1.2,
    )
    ax.text(
        left + 0.11,
        top - 0.19,
        system.symbol,
        fontsize=11,
        color=style.SLATE,
        va="center",
    )
    ax.text(
        left + 0.36,
        top - 0.20,
        system.name,
        fontsize=8.5,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    for index, line in enumerate(system.axioms):
        ax.text(
            left + 0.11,
            top - 0.44 - index * 0.185,
            line,
            fontsize=6.6,
            color=style.INK,
            va="center",
        )


def _spokes(ax: Axes) -> None:
    """Draw the conversions between independence and the ring.

    Only four of the six systems take a spoke: the page's table converts flats
    and hyperplanes through the closure and the circuits respectively, never
    directly to independence, and drawing a spoke there would claim a
    conversion the page does not record. Each spoke's heads are governed by the
    same rule — see :class:`_Spoke` for why the closure takes only one.
    """
    hub = (_WHEEL_CX, _WHEEL_CY)
    hub_half = (_HUB_W / 2.0 + _EDGE_GAP, _HUB_H / 2.0 + _EDGE_GAP)
    chip_half = (_CHIP_W / 2.0 + _EDGE_GAP, _CHIP_H / 2.0 + _EDGE_GAP)
    for system in _SYSTEMS:
        if system.spoke is _Spoke.NONE:
            continue
        chip = _system_center(system)
        if system.spoke is _Spoke.INWARD:
            start, end = _edge_endpoints(chip, hub, chip_half, hub_half)
        else:
            start, end = _edge_endpoints(hub, chip, hub_half, chip_half)
        _arrow(
            ax,
            start,
            end,
            color=_SPOKE_COLOR,
            two_way=system.spoke is _Spoke.BOTH,
        )


def _rim_bow(source: tuple[float, float], target: tuple[float, float]) -> float:
    """Return the bow that bends a rim edge away from the wheel's centre.

    Which sign of ``rad`` points outward depends on the direction of travel,
    not on which side of the wheel the edge sits: ``arc3`` always bulges along
    the chord's clockwise normal, so an edge running up the right of the ring
    and one running down the left bow outward at the *same* sign.
    """
    dx, dy = target[0] - source[0], target[1] - source[1]
    midpoint_x = (source[0] + target[0]) / 2.0
    midpoint_y = (source[1] + target[1]) / 2.0
    away = dy * (midpoint_x - _WHEEL_CX) - dx * (midpoint_y - _WHEEL_CY)
    return _RIM_BOW if away > 0 else -_RIM_BOW


def _rim(ax: Axes) -> None:
    """Draw the conversions that never pass through independence."""
    by_name = {system.name: system for system in _SYSTEMS}
    for edge in _RIM:
        source = _system_center(by_name[edge.source])
        target = _system_center(by_name[edge.target])
        half = (_CHIP_W / 2.0 + _EDGE_GAP, _CHIP_H / 2.0 + _EDGE_GAP)
        # Bow away from the wheel's centre, so the rim never crosses a spoke.
        outward = _rim_bow(source, target)
        start, end = _edge_endpoints(source, target, half, half, rad=outward)
        _arrow(ax, start, end, color=_RIM_COLOR, rad=outward)
        label_x, label_y = edge.label
        for index, line in enumerate(edge.tag):
            ax.text(
                label_x,
                label_y - index * 0.155,
                line,
                fontsize=6.5,
                color=_RIM_COLOR,
                ha=edge.ha,
                va="center",
            )


def _wheel_caption(ax: Axes) -> None:
    """Name the wheel, above the bases chip."""
    ax.text(
        _WHEEL_CX,
        _MAIN_TOP - 0.06,
        "SEVEN EQUIVALENT DEFINITIONS - pick one, derive the rest",
        fontsize=8.0,
        fontweight="bold",
        color=style.SLATE,
        ha="center",
        va="top",
    )


def _conversion_table(ax: Axes) -> None:
    """Draw the cryptomorphism conversion table beneath the wheel."""
    width = _MARGIN_R - _MARGIN_L
    height = _TABLE_TOP - _TABLE_BOTTOM
    _box(
        ax,
        (_MARGIN_L, _TABLE_BOTTOM, width, height),
        facecolor=style.PAPER,
        edgecolor=style.MIST,
        linewidth=1.0,
    )
    ax.text(
        _MARGIN_L + _PANEL_PAD,
        _TABLE_TOP - 0.17,
        "CRYPTOMORPHISM - THE CONVERSION TABLE",
        fontsize=8.0,
        fontweight="bold",
        color=_SPOKE_COLOR,
        va="center",
    )
    ax.text(
        _MARGIN_L + 3.45,
        _TABLE_TOP - 0.17,
        "each row is the definition of one primitive in terms of another",
        fontsize=7.5,
        color=style.MIST,
        va="center",
    )
    for index, (key, gloss) in enumerate(_CONVERSIONS):
        column, row = divmod(index, 3)
        x = _MARGIN_L + _PANEL_PAD + column * 5.0
        y = _TABLE_TOP - 0.40 - row * 0.175
        ax.text(x, y, key, fontsize=8.0, color=_SPOKE_COLOR, va="center")
        ax.text(x + 0.72, y, gloss, fontsize=7.0, color=style.INK, va="center")
    ax.text(
        _MARGIN_L + _PANEL_PAD,
        _TABLE_BOTTOM + 0.25,
        _CONVERSION_NOTE,
        fontsize=7.0,
        color=style.MIST,
        va="center",
    )
    ax.text(
        _MARGIN_L + _PANEL_PAD,
        _TABLE_BOTTOM + 0.11,
        _CONVERSION_RULE,
        fontsize=7.0,
        color=_RIM_COLOR,
        va="center",
    )


def _implementations(ax: Axes) -> None:
    """Draw the Python and Lean panels, and the tie from each to the wheel."""
    height = _IMPL_TOP - _IMPL_BOTTOM
    width = 7.35
    for index, panel in enumerate((_PYTHON, _LEAN)):
        x = _MARGIN_L + index * (width + 0.30)
        _box(
            ax,
            (x, _IMPL_BOTTOM, width, height),
            facecolor=style.PAPER,
            edgecolor=_BUILT_COLOR,
            linewidth=1.3,
        )
        ax.text(
            x + _PANEL_PAD,
            _IMPL_TOP - 0.18,
            panel.title,
            fontsize=8.5,
            fontweight="bold",
            color=_BUILT_COLOR,
            va="center",
        )
        ax.text(
            x + width - _PANEL_PAD,
            _IMPL_TOP - 0.18,
            "as built, 2026-07-28",
            fontsize=7.0,
            color=style.MIST,
            ha="right",
            va="center",
        )
        for line_index, line in enumerate(panel.lines):
            ax.text(
                x + _PANEL_PAD,
                _IMPL_TOP - 0.40 - line_index * 0.155,
                line,
                fontsize=7.0,
                color=style.INK,
                va="center",
            )


def _open_questions(ax: Axes) -> None:
    """Draw the strip of questions the page leaves open."""
    ax.text(
        _MARGIN_L,
        (_FOOTER_BOTTOM + _FOOTER_TOP) / 2.0,
        "OPEN\nQUESTIONS",
        fontsize=7.5,
        fontweight="bold",
        color=style.MIST,
        va="center",
        linespacing=1.4,
    )
    for index, lines in enumerate(_OPEN_QUESTIONS):
        x = 1.75 + index * 4.62
        for line_index, line in enumerate(lines):
            ax.text(
                x,
                _FOOTER_TOP - 0.14 - line_index * 0.155,
                line,
                fontsize=6.8,
                color=style.INK,
                va="center",
            )


def _render(step: int) -> Figure:
    """Return the breakdown revealed up to ``step``, in the full-figure frame.

    The frame is fixed by the layout constants rather than by what is drawn,
    so every step lands in an identically sized image and nothing already on
    screen moves as the next layer arrives.
    """
    # The figure is the drawing plus its border, and the axes spans the figure
    # with one drawing unit per inch, so the file is cropped to the graphic by
    # construction — a tight bbox would not crop it, since a full-figure axes
    # reports its own extent rather than the artists inside it.
    left, right = _MARGIN_L - _PAD, _MARGIN_R + _PAD
    bottom, top = _FOOTER_BOTTOM - _PAD, _HEADER_TOP + _PAD

    fig = plt.figure(figsize=(right - left, top - bottom))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_axis_off()

    _header(ax)
    if step >= _STEP_HUB:
        _hub(ax)
    if step >= _STEP_WHEEL:
        _wheel_caption(ax)
        for system in _SYSTEMS:
            _chip(ax, system)
        _spokes(ax)
        _rim(ax)
    if step >= _STEP_TABLE:
        _conversion_table(ax)
    if step >= _STEP_LEFT:
        _panel(ax, (_COL_L, _UPPER_BOTTOM, _COL_W, _UPPER_H), _VOCABULARY)
        _panel(ax, (_COL_L, _MAIN_BOTTOM, _COL_W, _LOWER_H), _OPERATIONS)
    if step >= _STEP_RIGHT:
        _panel(ax, (_COL_R, _UPPER_BOTTOM, _COL_W, _UPPER_H), _THEOREMS)
        _panel(ax, (_COL_R, _MAIN_BOTTOM, _COL_W, _LOWER_H), _FIXTURES)
    if step >= _STEP_IMPL:
        _implementations(ax)
    if step >= _STEP_OPEN:
        _open_questions(ax)
    return fig


# --------------------------------------------------------------------------
# Link figure: the page's neighbourhood in the graph.
# --------------------------------------------------------------------------

_LINK_ROW_H: float = 1.83
"""Height of one cluster cell in the link graph; sized to the largest cell's
eight rows, so every cell is the same height and the grid reads as a grid."""

_LINK_ROW_GAP: float = 0.72
"""Vertical gap between cluster rows, holding the arrows out of the hub."""

_LINK_COL_GAP: float = 0.50
"""Horizontal gap between cluster columns."""

_LINK_COL_W: float = (_MARGIN_R - _MARGIN_L - 2 * _LINK_COL_GAP) / 3.0
"""Width of one cluster cell."""

_LINK_BOTTOM: float = 1.00
"""Bottom edge of the lowest cluster row."""

_LINK_MID: float = _LINK_BOTTOM + _LINK_ROW_H + _LINK_ROW_GAP
"""Bottom edge of the middle cluster row, which also holds the hub."""

_LINK_TOP_ROW: float = _LINK_MID + _LINK_ROW_H + _LINK_ROW_GAP
"""Bottom edge of the highest cluster row."""

_LINK_HEADER_BOTTOM: float = _LINK_TOP_ROW + _LINK_ROW_H + 0.18
"""Bottom edge of the link figure's header band."""

_LINK_HEADER_TOP: float = _LINK_HEADER_BOTTOM + 0.81
"""Top edge of the link figure's header band, and of that drawing."""

_LINK_LEGEND_Y: float = 0.55
"""Baseline of the link figure's legend, below the lowest row."""

_LINK_HUB_W: float = 2.50
"""Width of the Matroid hub in the link graph."""

_LINK_HUB_H: float = 1.15
"""Height of the Matroid hub in the link graph."""

_BACKLINK_COLOR: str = style.CATEGORICAL[4]
"""Purple — carries the six pages that cite Matroid, whose arrows run inward.
Direction and the cluster's own caption already distinguish them; the colour
only reinforces it, since this palette cannot carry identity alone."""


@dataclass(frozen=True, slots=True)
class _Cluster:
    """One group of linked pages, placed in a cell of the three-by-three grid.

    ``column`` and ``row`` index the grid from the top left; the centre cell
    is the hub and never carries a cluster. ``columns`` is how many columns
    the titles are laid out in inside the cell.
    """

    title: str
    column: int
    row: int
    pages: tuple[str, ...]
    columns: int = 1
    inbound: bool = False


_CLUSTERS: tuple[_Cluster, ...] = (
    _Cluster(
        title="AXIOM SYSTEMS AND PRIMITIVES",
        column=0,
        row=0,
        columns=2,
        pages=(
            "Independent Set",
            "Matroid Basis",
            "Matroid Circuit",
            "Rank Function",
            "Matroid Closure",
            "Matroid Flat",
            "Matroid Hyperplane",
            "Ground Set",
            "Cryptomorphism",
            "Basis Exchange",
            "Circuit Elimination",
            "Submodularity",
            "Mac Lane–Steinitz Exchange",
        ),
    ),
    _Cluster(
        title="DERIVED VOCABULARY",
        column=1,
        row=0,
        pages=(
            "Cocircuit",
            "Coloop",
            "Matroid Loop",
            "Parallel Elements",
            "Simple Matroid",
            "Nullity",
            "Fundamental Circuit",
            "Matroid Connectivity",
        ),
    ),
    _Cluster(
        title="OPERATIONS AND CONSTRUCTIONS",
        column=2,
        row=0,
        columns=2,
        pages=(
            "Matroid Restriction",
            "Matroid Deletion",
            "Matroid Contraction",
            "Matroid Minor",
            "Matroid Duality",
            "Direct Sum",
            "Truncation",
            "Matroid Union",
            "Matroid Intersection",
            "Lattice of Flats",
        ),
    ),
    _Cluster(
        title="WHAT IT ABSTRACTS, AND WHERE IT SITS",
        column=0,
        row=1,
        pages=(
            "Linear Independence",
            "Algebraic Independence",
            "Graph",
            "Transversal",
            "Mathematical Structure",
            "Combinatorics",
            "Matroid Theory",
        ),
    ),
    _Cluster(
        title="STRUCTURAL THEOREMS",
        column=2,
        row=1,
        pages=(
            "Rado–Edmonds Theorem",
            "Greedy Algorithm",
            "Circuit–Cocircuit Orthogonality",
            "Excluded Minor",
            "Tutte Polynomial",
        ),
    ),
    _Cluster(
        title="PEOPLE, SOURCES, AND TOOLING",
        column=0,
        row=2,
        pages=(
            "Hassler Whitney",
            "Takeo Nakasawa",
            "Saunders MacLane",
            "James Oxley",
            "W. T. Tutte",
            "Mathlib",
            "SageMath",
        ),
    ),
    _Cluster(
        title="CLASSES AND CANONICAL EXAMPLES",
        column=1,
        row=2,
        columns=2,
        pages=(
            "Uniform Matroid",
            "Free Matroid",
            "Binary Matroid",
            "Regular Matroid",
            "Graphic Matroid",
            "Transversal Matroid",
            "Representable Matroid",
            "Matroid Representability",
            "Fano Matroid",
            "Non-Fano Matroid",
            "Vámos Matroid",
            "Pappus Matroid",
            "Rota's Conjecture",
            "Whitney's 2-Isomorphism Theorem",
        ),
    ),
    _Cluster(
        title="CITED BY - PAGES THAT LINK HERE",
        column=2,
        row=2,
        inbound=True,
        pages=(
            "Positroid",
            "Grassmannian",
            "Grassmann Necklace",
            "Plabic Graph",
            "Amplituhedron",
            "Totally Nonnegative Grassmannian",
        ),
    ),
)
"""The page's 64 forward links grouped by the role they play on it, plus the
6 backlinks as a cluster of their own.

The grouping is the whole point of the figure: a flat list of 64 titles is
unreadable on a slide, whereas the split shows at a glance that classes and
fixtures are the largest group and that only seven links are people or
tooling. Every title here is the page's actual link target, so the counts add
back to 64 and 6.
"""

_FORWARD_COUNT: int = sum(len(c.pages) for c in _CLUSTERS if not c.inbound)
"""Forward links drawn, checked against the graph's own count in the header."""

_BACKLINK_COUNT: int = sum(len(c.pages) for c in _CLUSTERS if c.inbound)
"""Backlinks drawn."""


def _cluster_rect(cluster: _Cluster) -> tuple[float, float, float, float]:
    """Return ``(x, y, width, height)`` for ``cluster``'s cell."""
    x = _MARGIN_L + cluster.column * (_LINK_COL_W + _LINK_COL_GAP)
    y = (_LINK_TOP_ROW, _LINK_MID, _LINK_BOTTOM)[cluster.row]
    return x, y, _LINK_COL_W, _LINK_ROW_H


def _cluster(ax: Axes, cluster: _Cluster) -> None:
    """Draw one cluster cell: its caption, count, and page titles."""
    x, y, width, height = _cluster_rect(cluster)
    color = _BACKLINK_COLOR if cluster.inbound else style.SLATE
    _box(
        ax,
        (x, y, width, height),
        facecolor=style.PAPER,
        edgecolor=color if cluster.inbound else style.MIST,
        linewidth=1.2 if cluster.inbound else 1.0,
    )
    ax.text(
        x + _PANEL_PAD,
        y + height - 0.17,
        cluster.title,
        fontsize=7.8,
        fontweight="bold",
        color=color,
        va="center",
    )
    ax.text(
        x + width - _PANEL_PAD,
        y + height - 0.17,
        str(len(cluster.pages)),
        fontsize=8.5,
        fontweight="bold",
        color=style.MIST,
        ha="right",
        va="center",
    )
    per_column = -(-len(cluster.pages) // cluster.columns)  # ceiling division
    sub_width = (width - 2 * _PANEL_PAD) / cluster.columns
    for index, page in enumerate(cluster.pages):
        column, row = divmod(index, per_column)
        ax.text(
            x + _PANEL_PAD + column * sub_width,
            y + height - 0.40 - row * 0.185,
            page,
            fontsize=7.0,
            color=style.INK,
            va="center",
        )


def _link_hub(ax: Axes) -> tuple[float, float]:
    """Draw the Matroid hub of the link graph, and return its centre."""
    cx = _MARGIN_L + (_MARGIN_R - _MARGIN_L) / 2.0
    cy = _LINK_MID + _LINK_ROW_H / 2.0
    _box(
        ax,
        (cx - _LINK_HUB_W / 2.0, cy - _LINK_HUB_H / 2.0, _LINK_HUB_W, _LINK_HUB_H),
        facecolor=style.PAPER,
        edgecolor=style.SLATE,
        linewidth=1.8,
    )
    ax.text(
        cx,
        cy + 0.26,
        "Matroid",
        fontsize=16,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        cx,
        cy - 0.06,
        f"{_FORWARD_COUNT} forward links",
        fontsize=8.0,
        color=style.SLATE,
        ha="center",
        va="center",
    )
    ax.text(
        cx,
        cy - 0.28,
        f"{_BACKLINK_COUNT} backlinks",
        fontsize=8.0,
        color=_BACKLINK_COLOR,
        ha="center",
        va="center",
    )
    return cx, cy


def _link_edges(ax: Axes, hub: tuple[float, float]) -> None:
    """Draw one arrow between the hub and each cluster cell.

    Corner cells are met at their nearest corner and side cells at the middle
    of their facing edge, so no arrow crosses a cell it does not belong to.
    """
    cx, cy = hub
    half_w, half_h = _LINK_HUB_W / 2.0, _LINK_HUB_H / 2.0
    for cluster in _CLUSTERS:
        x, y, width, height = _cluster_rect(cluster)
        if cluster.column == 1:  # directly above or below the hub
            start = (cx, cy + half_h if cluster.row == 0 else cy - half_h)
            end = (cx, y if cluster.row == 0 else y + height)
        elif cluster.row == 1:  # directly left or right of the hub
            on_left = cluster.column == 0
            start = (cx - half_w if on_left else cx + half_w, cy)
            end = (x + width if on_left else x, cy)
        else:  # a corner cell, met diagonally
            on_left = cluster.column == 0
            above = cluster.row == 0
            start = (
                cx - half_w + 0.10 if on_left else cx + half_w - 0.10,
                cy + half_h if above else cy - half_h,
            )
            end = (
                x + width - 0.15 if on_left else x + 0.15,
                y - 0.04 if above else y + height + 0.04,
            )
        if cluster.inbound:
            start, end = end, start
        _arrow(
            ax,
            start,
            end,
            color=_BACKLINK_COLOR if cluster.inbound else style.SLATE,
            linewidth=1.3,
        )


def _link_header(ax: Axes) -> None:
    """Draw the link figure's title and the count it is accounting for."""
    ax.text(
        _MARGIN_L,
        _LINK_HEADER_BOTTOM + 0.52,
        "Matroid - the page's neighbourhood",
        fontsize=20,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    ax.text(
        _MARGIN_L,
        _LINK_HEADER_BOTTOM + 0.15,
        f"every page the Logseq entry links to, grouped by the role it plays there - "
        f"{_FORWARD_COUNT} outward, {_BACKLINK_COUNT} inward",
        fontsize=8.5,
        color=style.MIST,
        va="center",
    )


def _link_legend(ax: Axes) -> None:
    """Caption the two edge directions, below the lowest cluster row."""
    ax.text(
        _MARGIN_L,
        _LINK_LEGEND_Y,
        "outward - concepts this page defines, uses, or cites",
        fontsize=7.5,
        color=style.SLATE,
        va="center",
    )
    ax.text(
        _MARGIN_L + 5.6,
        _LINK_LEGEND_Y,
        "inward - the positroid and Grassmannian line of work builds on this page",
        fontsize=7.5,
        color=_BACKLINK_COLOR,
        va="center",
    )


def _render_links() -> Figure:
    """Return the link-graph figure, cropped to the drawing."""
    left, right = _MARGIN_L - _PAD, _MARGIN_R + _PAD
    bottom, top = _LINK_LEGEND_Y - 0.20 - _PAD, _LINK_HEADER_TOP + _PAD

    fig = plt.figure(figsize=(right - left, top - bottom))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_axis_off()

    _link_header(ax)
    hub = _link_hub(ax)
    _link_edges(ax, hub)
    for cluster in _CLUSTERS:
        _cluster(ax, cluster)
    _link_legend(ax)
    return fig


@figure(name="matroid-page")
def matroid_page(ctx: FigureContext) -> None:
    """Render the Matroid page breakdown, its build-up, and its link graph.

    Writes ``matroid-page`` (the complete breakdown),
    ``matroid-page-NN`` for each reveal step in presentation order, and
    ``matroid-page-links`` (the knowledge-graph neighbourhood).
    """
    # The house style saves with a tight bbox, which would re-crop the frames
    # that were just sized deliberately — and would crop each step to its own
    # contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"matroid-page-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "matroid-page")
        plt.close(fig)

        fig = _render_links()
        ctx.save(fig, "matroid-page-links")
        plt.close(fig)
