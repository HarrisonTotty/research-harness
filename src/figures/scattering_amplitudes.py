r"""Two matched slide graphics: how an amplitude is computed, before and after.

Renders the same question — *where does a scattering amplitude come from?* —
answered twice, as a pair of 16:9 diagrams meant for consecutive slides. The
first, ``scattering-feynman``, is the textbook route: fix a Lagrangian, draw
every topology, evaluate each one off-shell, and watch the gauge dependence
cancel between terms. The second, ``scattering-amplituhedron``, is the route of
Arkani-Hamed and Trnka: feed one point of the positive Grassmannian and one
positive kinematic matrix into a single linear map, and read the amplitude off
the boundary of the region that comes out.

Both figures share one frame, one four-stage band, one panel band, and one
footer strip, in identical positions — so advancing from the first slide to the
second replaces the mechanism in place while the scaffolding holds still. That
substitution *is* the argument, and it is why these are two figures rather than
a split one.

Each carries a computed panel rather than a cartoon:

* the Feynman figure plots the enumeration cost — colour-ordered tree diagrams
  for :math:`n` gluons, counted as dissections of a convex :math:`n`-gon into
  triangles and quadrilaterals, which run ``1, 3, 10, 38, 154, 654, 2871,
  12925`` for :math:`n = 3, \dots, 10` (Elvang and Huang §2.6, verified in the
  knowledge graph's **BCFW Recursion** page), annotated at the one head-to-head
  comparison that page certifies: the split-helicity 6-gluon NMHV amplitude,
  under a :math:`[1,2\rangle` shift, costs 38 diagrams and 2 BCFW terms;
* the amplituhedron figure draws the map itself at :math:`k = 1`, :math:`m = 2`,
  :math:`n = 6`, where the amplituhedron genuinely *is* the interior of a convex
  hexagon. Six vertices in counterclockwise convex position are a positive
  :math:`Z` (every ordered :math:`3 \times 3` minor is twice a positive signed
  area), the twenty 2-dimensional cells of :math:`Gr^{\geq 0}_{1,6}` are the
  three-element supports, and the four of the form :math:`\{1, i, i{+}1\}` map
  to the four triangles of the fan from :math:`Z_1` — matching the proved
  :math:`m = 2` tile count :math:`\binom{n-2}{k} = 4`. Cells and tiles are
  drawn in the same four ordinal colors, so the correspondence is visible
  rather than asserted.

The footer strip is the punchline both panels build to, taken from §11 of
Arkani-Hamed and Trnka: locality and unitarity are *inputs* on the left and
*consequences of positivity* on the right, and gauge redundancy — which infests
every intermediate step of the first route — never appears in the second.

Neither graphic carries a title or summary chrome; the slide supplies those.
Both are saved cropped to the drawing itself so they can be placed and scaled
freely, and both write a build-up sequence — one image per reveal step, each
adding the next piece — into the same frame as the whole, so the steps can be
stacked on one slide without anything shifting between them::

    1  stage 1        5  panel, left
    2  stage 2        6  panel, right
    3  stage 3        7  footer strip
    4  stage 4

Regenerate with ``just figure scattering-amplitudes``.
"""

import math
from dataclasses import dataclass

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon

from figures import style
from figures.cli import FigureContext, figure

# --------------------------------------------------------------------------
# Shared frame — both figures are laid out on exactly these coordinates.
# --------------------------------------------------------------------------

_PAD: float = 0.2
"""Border kept around the drawing on every side, in inches."""

_MARGIN_L: float = 0.5
"""Left edge of every full-width element."""

_MARGIN_R: float = 16.5
"""Right edge of every full-width element."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_CONTENT_TOP: float = 9.20
"""Top of the drawing: the cap height of the eyebrow line."""

_EYEBROW_Y: float = 8.95
"""Vertical center of the one-line mechanism summary above the stage row."""

_STAGE_TOP: float = 8.45
"""Top edge of the stage row."""

_STAGE_H: float = 2.70
"""Height of a pipeline stage box."""

_STAGE_BOTTOM: float = _STAGE_TOP - _STAGE_H
"""Bottom edge of the stage row."""

_STAGE_GAP: float = 0.5
"""Horizontal gap between consecutive stage boxes, holding the flow arrow."""

_STAGE_W: float = (_MARGIN_R - _MARGIN_L - 3.0 * _STAGE_GAP) / 4.0
"""Width of a pipeline stage box: four boxes and three gaps span the margins."""

_PANEL_TOP: float = 5.50
"""Top edge of the band holding each figure's computed panel."""

_PANEL_BOTTOM: float = 1.85
"""Bottom edge of the panel band."""

_STRIP_TOP: float = 1.55
"""Top edge of the input-versus-consequence footer strip."""

_STRIP_BOTTOM: float = 0.55
"""Bottom edge of the footer strip, and of the whole drawing."""

_CONTENT_BOTTOM: float = _STRIP_BOTTOM
"""Bottom of the drawing."""

_FRAME_LEFT: float = _MARGIN_L - _PAD
"""Left limit of the full-figure axes, in inches."""

_FRAME_RIGHT: float = _MARGIN_R + _PAD
"""Right limit of the full-figure axes, in inches."""

_FRAME_BOTTOM: float = _CONTENT_BOTTOM - _PAD
"""Bottom limit of the full-figure axes, in inches."""

_FRAME_TOP: float = _CONTENT_TOP + _PAD
"""Top limit of the full-figure axes, in inches."""

_STRIP_GAP: float = 0.4
"""Horizontal gap between the three footer chips."""

_STRIP_W: float = (_MARGIN_R - _MARGIN_L - 2.0 * _STRIP_GAP) / 3.0
"""Width of one footer chip."""

_TINT: str = style.SEQUENTIAL_ANCHORS[0]
"""Lightest step of the blog's slate ramp; fill for the footer chips."""

_COST_COLOR: str = style.CATEGORICAL[5]
"""Indianred — marks the redundancy the Feynman route pays for."""

_MAP_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the amplituhedron map, distinct from the stage flow."""

_ORDINAL_TEXT: tuple[str, ...] = (
    style.INK,
    style.PAPER,
    style.PAPER,
    style.PAPER,
)
"""Label color for a mark filled with the matching :data:`style.ORDINAL_STEPS`
step. The ramp's light end is too pale to carry white type, so only the lower
three steps get it."""

_TOTAL_STEPS: int = 7
"""Number of images in each build-up sequence; the last is the whole diagram."""

_STEP_PANEL_LEFT: int = 5
"""Reveal step adding the left half of the panel band."""

_STEP_PANEL_RIGHT: int = 6
"""Reveal step adding the right half of the panel band."""

_STEP_STRIP: int = 7
"""Reveal step adding the footer strip that contrasts the two routes."""


@dataclass(frozen=True, slots=True)
class _Stage:
    """One box in a four-stage mechanism band.

    Body lines are pre-wrapped rather than flowed: at this size the wrap
    points are a layout decision, not something to leave to a text engine.
    """

    number: int
    title: str
    body: tuple[str, ...]
    tag: str


@dataclass(frozen=True, slots=True)
class _Note:
    """One lead/detail pair in a panel's right-hand reading column."""

    lead: str
    detail: str
    accent: bool = False


@dataclass(frozen=True, slots=True)
class _Chip:
    """One footer chip: a property, and how the route obtains it."""

    label: str
    body: tuple[str, str]


# --------------------------------------------------------------------------
# Figure one: the Feynman expansion.
# --------------------------------------------------------------------------

_FEYNMAN_EYEBROW: str = (
    "THE FEYNMAN EXPANSION   ·   build the answer out of parts, "
    "then watch the parts cancel"
)
"""One-line mechanism summary above the Feynman stage row."""

_FEYNMAN_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="LAGRANGIAN",
        body=(
            "fix the Yang-Mills action and",
            "a gauge, then read off the",
            "three- and four-gluon vertices",
            "and the propagator",
        ),
        tag="the gauge choice enters here",
    ),
    _Stage(
        number=2,
        title="ENUMERATE",
        body=(
            "draw every tree topology on n",
            "colour-ordered legs: the",
            "dissections of an n-gon into",
            "triangles and quadrilaterals",
        ),
        tag="1, 3, 10, 38, 154, 654, ...",
    ),
    _Stage(
        number=3,
        title="EVALUATE",
        body=(
            "contract each diagram into a",
            "rational function; internal",
            "lines go off-shell and every",
            "term depends on the gauge",
        ),
        tag="no single term is physical",
    ),
    _Stage(
        number=4,
        title="CANCEL",
        body=(
            "add them up: the gauge and",
            "off-shell dependence cancels",
            "between terms, collapsing the",
            "sum to a handful of pieces",
        ),
        tag="38 diagrams in, 2 terms out",
    ),
)
"""The textbook pipeline, left to right."""

_GLUONS: tuple[int, ...] = (3, 4, 5, 6, 7, 8, 9, 10)
"""Number of external gluons for each plotted diagram count."""

_TREE_DIAGRAMS: tuple[int, ...] = (1, 3, 10, 38, 154, 654, 2871, 12925)
"""Colour-ordered tree Feynman diagrams for :data:`_GLUONS` gluons.

Dissections of a convex :math:`n`-gon into triangles and quadrilaterals, one
per topology built from the cubic and quartic Yang-Mills vertices. The first
five entries are Elvang-Huang §2.6; the sequence and its continuation are
recorded as verified on the knowledge graph's **BCFW Recursion** page.
"""

_HIGHLIGHT_GLUONS: int = 6
r"""The leg count called out on the curve.

The one head-to-head comparison the knowledge graph certifies: 38 diagrams
against 2 BCFW terms, for the split-helicity amplitude
:math:`A_6[1^-2^-3^-4^+5^+6^+]` under a :math:`[1,2\rangle` shift. The shift is
part of the claim, not decoration - the same amplitude under :math:`[2,1\rangle`
has three terms, so a term count is never a property of the amplitude alone.
"""

_FEYNMAN_NOTES: tuple[_Note, ...] = (
    _Note(
        lead="every term needs a gauge",
        detail="the amplitude it adds up to does not",
    ),
    _Note(
        lead="every internal line is off-shell",
        detail="the amplitude is on-shell throughout",
    ),
    _Note(
        lead=r"terms grow as $z,\ z^2,\ z^3$ at large complex momenta",
        detail=r"a valid shift makes their sum fall as $1/z$: pure cancellation",
    ),
    _Note(
        lead=r"n = 6, split helicity, $[1,2\rangle$ shift",
        detail="38 diagrams carry exactly what 2 BCFW terms carry",
        accent=True,
    ),
)
"""The reading column beside the growth curve."""

_FEYNMAN_CHIPS: tuple[_Chip, ...] = (
    _Chip(
        label="LOCALITY",
        body=(
            "put in by hand: one propagator per",
            "internal line, one vertex per interaction",
        ),
    ),
    _Chip(
        label="UNITARITY",
        body=(
            "put in by hand: a sum over the",
            "intermediate states the theory allows",
        ),
    ),
    _Chip(
        label="GAUGE REDUNDANCY",
        body=(
            "in every single term, and gone from",
            "the answer they add up to",
        ),
    ),
)
"""What the Feynman route must assume, and what it pays for assuming it."""


# --------------------------------------------------------------------------
# Figure two: the amplituhedron.
# --------------------------------------------------------------------------

_AMPLITUHEDRON_EYEBROW: str = (
    "THE AMPLITUHEDRON   ·   build one geometry, then read the answer off its boundary"
)
"""One-line mechanism summary above the amplituhedron stage row."""

_AMPLITUHEDRON_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="POSITIVE DATA",
        body=(
            r"one point $C$ of the positive",
            r"Grassmannian $Gr^{\geq 0}_{k,n}$, and",
            r"kinematic data $Z$ whose ordered",
            r"maximal minors are all $> 0$",
        ),
        tag=r"$Z$ = the n bosonized twistors",
    ),
    _Stage(
        number=2,
        title="ONE MAP",
        body=(
            r"push $C$ through $Z$: $Y = C \cdot Z$.",
            r"as $C$ sweeps the whole positive",
            r"Grassmannian, $Y$ sweeps out one",
            r"region $\mathcal{A}_{n,k,4}(Z)$",
        ),
        tag=r"$\dim = 4k$ — nothing is summed",
    ),
    _Stage(
        number=3,
        title="TILE",
        body=(
            "the BCFW cells of the positive",
            "Grassmannian land on disjoint",
            "tiles whose closures cover that",
            "region, each covered once",
        ),
        tag=r"$N(n{-}3,\, k{+}1)$ of them",
    ),
    _Stage(
        number=4,
        title="CANONICAL FORM",
        body=(
            r"the region carries one form $\Omega$",
            "with a logarithmic pole on every",
            "boundary and nowhere else;",
            "the amplitude is that form",
        ),
        tag="one geometry, one answer",
    ),
)
"""The positive-geometry pipeline, left to right."""

_N_PARTICLES: int = 6
"""Legs in the worked example drawn in the panel."""

_CELLS_2D: tuple[tuple[int, ...], ...] = tuple(
    (a, b, c)
    for a in range(1, _N_PARTICLES + 1)
    for b in range(a + 1, _N_PARTICLES + 1)
    for c in range(b + 1, _N_PARTICLES + 1)
)
"""The 2-dimensional cells of :math:`Gr^{\\geq 0}_{1,6}`, in lexicographic order.

A cell of :math:`Gr^{\\geq 0}_{1,n}` is fixed by the support of its single row,
and has dimension one less than that support's size — so the 2-dimensional
cells are exactly the three-element subsets of :math:`[n]`.
"""

_TILING_CELLS: tuple[tuple[int, ...], ...] = tuple(
    (1, i, i + 1) for i in range(2, _N_PARTICLES)
)
"""The four cells whose images tile the hexagon: the fan from :math:`Z_1`."""

_CELL_COLS: int = 5
"""Columns in the grid of 2-dimensional cells."""

_AMPLITUHEDRON_NOTES: tuple[_Note, ...] = (
    _Note(
        lead=r"$\dim \mathcal{A}_{6,1,2} = km = 2$",
        detail="one region, not a sum over topologies",
    ),
    _Note(
        lead=r"$\binom{n-2}{k} = 4$ tiles",
        detail="and every tiling of it has the same four",
    ),
    _Note(
        lead="the six edges are the poles",
        detail=r"$\langle Y\, Z_i\, Z_{i+1} \rangle = 0$, and nothing else is",
    ),
    _Note(
        lead=r"$\Omega$ is one dlog form per tile",
        detail=r"at $m = 4$: $\dim = 4k$, $N(n{-}3,\,k{+}1)$ BCFW tiles",
        accent=True,
    ),
)
"""The reading column beside the worked hexagon."""

_AMPLITUHEDRON_CHIPS: tuple[_Chip, ...] = (
    _Chip(
        label="LOCALITY",
        body=(
            "falls out: the only boundaries sit at",
            r"$\langle Y\, Z_i\, Z_{i+1}\, Z_j\, Z_{j+1} \rangle = 0$",
        ),
    ),
    _Chip(
        label="UNITARITY",
        body=(
            "falls out: on a boundary, positivity",
            r"forces $C$ to split, $k_L + k_R = k - 1$",
        ),
    ),
    _Chip(
        label="GAUGE REDUNDANCY",
        body=(
            "never appears — the construction has",
            "no off-shell step for it to live in",
        ),
    ),
)
"""What positivity hands back for free, in the same slots as the left figure."""


# --------------------------------------------------------------------------
# Shared drawing primitives.
# --------------------------------------------------------------------------


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
    linewidth: float = 1.3,
) -> None:
    """Draw a single straight arrow from ``start`` to ``end``."""
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=linewidth,
            color=color,
            shrinkA=0,
            shrinkB=0,
        )
    )


def _stage_x(index: int) -> float:
    """Return the left edge of the stage box at ``index``."""
    return _MARGIN_L + index * (_STAGE_W + _STAGE_GAP)


def _stage(ax: Axes, index: int, stage: _Stage) -> None:
    """Draw the stage box at ``index``: badge, title, body, and tag line."""
    x = _stage_x(index)
    _box(
        ax,
        (x, _STAGE_BOTTOM, _STAGE_W, _STAGE_H),
        facecolor=style.PAPER,
        edgecolor=style.SLATE,
        linewidth=1.4,
    )

    badge_y = _STAGE_TOP - 0.40
    ax.add_patch(Circle((x + 0.38, badge_y), 0.22, facecolor=style.SLATE, lw=0))
    ax.text(
        x + 0.38,
        badge_y,
        str(stage.number),
        fontsize=12,
        fontweight="bold",
        color=style.PAPER,
        ha="center",
        va="center",
    )
    ax.text(
        x + 0.72,
        badge_y,
        stage.title,
        fontsize=17,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    ax.plot(
        [x + 0.2, x + _STAGE_W - 0.2],
        [_STAGE_TOP - 0.78, _STAGE_TOP - 0.78],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )

    for line_index, line in enumerate(stage.body):
        ax.text(
            x + 0.2,
            _STAGE_TOP - 1.12 - line_index * 0.34,
            line,
            fontsize=12,
            color=style.INK,
            va="center",
        )
    ax.text(
        x + 0.2,
        _STAGE_BOTTOM + 0.26,
        stage.tag,
        fontsize=12,
        color=style.MIST,
        va="center",
    )


def _stages(ax: Axes, stages: tuple[_Stage, ...], step: int) -> None:
    """Draw every stage revealed by ``step``, with the arrows between them."""
    mid_y = _STAGE_BOTTOM + _STAGE_H / 2.0
    for index, stage in enumerate(stages):
        if step < index + 1:
            continue
        _stage(ax, index, stage)
        if index > 0:
            _arrow(
                ax,
                (_stage_x(index - 1) + _STAGE_W, mid_y),
                (_stage_x(index), mid_y),
                color=style.SLATE,
            )


def _eyebrow(ax: Axes, text: str) -> None:
    """Draw the one-line mechanism summary above the stage row."""
    ax.text(
        _MARGIN_L,
        _EYEBROW_Y,
        text,
        fontsize=13,
        fontweight="bold",
        color=style.MIST,
        va="center",
    )


def _notes(
    ax: Axes,
    left: float,
    header: str,
    notes: tuple[_Note, ...],
    *,
    accent: str,
) -> None:
    """Draw a panel's reading column, headed and ruled, from ``left``.

    Args:
        ax: The full-figure axes.
        left: Left edge of the column, in drawing inches.
        header: Column heading, drawn above a rule.
        notes: Lead/detail pairs, top to bottom.
        accent: Color for the notes flagged ``accent``. The blog palette reads
            indianred as "falsified", so only the cost figure may use it; the
            amplituhedron column accents in the teal that already carries its
            map, and the difference is itself part of the contrast.
    """
    ax.text(
        left,
        _PANEL_TOP - 0.25,
        header,
        fontsize=14,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    ax.plot(
        [left, _MARGIN_R],
        [_PANEL_TOP - 0.50, _PANEL_TOP - 0.50],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )
    for index, note in enumerate(notes):
        top = _PANEL_TOP - 0.90 - index * 0.72
        ax.text(
            left,
            top,
            note.lead,
            fontsize=13,
            fontweight="bold",
            color=accent if note.accent else style.INK,
            va="center",
        )
        ax.text(
            left,
            top - 0.32,
            note.detail,
            fontsize=12,
            color=accent if note.accent else style.MIST,
            va="center",
        )


def _strip(ax: Axes, chips: tuple[_Chip, ...], step: int) -> None:
    """Draw the three-chip footer contrasting inputs against consequences."""
    if step < _STEP_STRIP:
        return
    for index, chip in enumerate(chips):
        x = _MARGIN_L + index * (_STRIP_W + _STRIP_GAP)
        _box(
            ax,
            (x, _STRIP_BOTTOM, _STRIP_W, _STRIP_TOP - _STRIP_BOTTOM),
            facecolor=_TINT,
        )
        ax.text(
            x + 0.22,
            _STRIP_TOP - 0.26,
            chip.label,
            fontsize=12,
            fontweight="bold",
            color=style.INK,
            va="center",
        )
        for line_index, line in enumerate(chip.body):
            ax.text(
                x + 0.22,
                _STRIP_TOP - 0.56 - line_index * 0.30,
                line,
                fontsize=11,
                color=style.SLATE,
                va="center",
            )


def _new_figure() -> tuple[Figure, Axes]:
    """Return the shared frame: a figure spanning one drawing inch per unit.

    The frame is fixed by the layout constants rather than by what is drawn,
    so every reveal step of both figures lands in an identically sized image,
    nothing moves as the next piece arrives, and the two finished graphics
    swap in place on consecutive slides.
    """
    # The axes spans the whole figure at one drawing unit per inch, so the file
    # is cropped to the graphic by construction — a tight bbox would not crop
    # it, since a full-figure axes reports its own extent rather than the
    # artists inside it.
    fig = plt.figure(figsize=(_FRAME_RIGHT - _FRAME_LEFT, _FRAME_TOP - _FRAME_BOTTOM))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(_FRAME_LEFT, _FRAME_RIGHT)
    ax.set_ylim(_FRAME_BOTTOM, _FRAME_TOP)
    ax.set_aspect("equal")
    ax.set_axis_off()
    return fig, ax


def _inset(fig: Figure, rect: tuple[float, float, float, float]) -> Axes:
    """Add a plotting axes at ``rect``, given in drawing inches.

    Args:
        fig: The figure built by :func:`_new_figure`.
        rect: ``(x, y, width, height)`` of the axes box, lower-left origin, in
            the same inch coordinates as every other element.

    Returns:
        The new axes, positioned by converting ``rect`` into the figure
        fractions ``Figure.add_axes`` expects.
    """
    x, y, width, height = rect
    span_x = _FRAME_RIGHT - _FRAME_LEFT
    span_y = _FRAME_TOP - _FRAME_BOTTOM
    return fig.add_axes(
        (
            (x - _FRAME_LEFT) / span_x,
            (y - _FRAME_BOTTOM) / span_y,
            width / span_x,
            height / span_y,
        )
    )


# --------------------------------------------------------------------------
# Figure one: panel.
# --------------------------------------------------------------------------


def _growth_panel(fig: Figure) -> None:
    """Plot the enumeration cost: tree diagrams against the number of gluons.

    The vertical scale is logarithmic because the sequence is — the point of
    the panel is that the number of parts to write down outruns the answer's
    length, not the particular values along the way.
    """
    ax = _inset(fig, (1.45, 2.55, 6.75, 2.55))
    ax.semilogy(
        _GLUONS,
        _TREE_DIAGRAMS,
        color=style.SLATE,
        marker="o",
        markerfacecolor=style.PAPER,
        markeredgewidth=2.0,
    )

    index = _GLUONS.index(_HIGHLIGHT_GLUONS)
    ax.plot(
        [_HIGHLIGHT_GLUONS],
        [_TREE_DIAGRAMS[index]],
        marker="o",
        color=_COST_COLOR,
        markersize=9,
    )
    ax.annotate(
        "38 diagrams",
        xy=(_HIGHLIGHT_GLUONS, _TREE_DIAGRAMS[index]),
        xytext=(6.35, 2.5),
        fontsize=12,
        color=_COST_COLOR,
        arrowprops={"arrowstyle": "-", "color": _COST_COLOR, "linewidth": 1.0},
    )
    ax.annotate(
        f"{_TREE_DIAGRAMS[-1]:,}",
        xy=(_GLUONS[-1], _TREE_DIAGRAMS[-1]),
        xytext=(-6, 4),
        textcoords="offset points",
        fontsize=12,
        color=style.SLATE,
        ha="right",
    )

    ax.set_title("the enumeration cost", fontsize=14)
    ax.set_xlabel("external gluons n", fontsize=12)
    ax.set_ylabel("colour-ordered tree diagrams", fontsize=12)
    ax.set_xticks(list(_GLUONS))
    ax.set_xlim(2.6, 10.4)
    ax.tick_params(labelsize=11)


# --------------------------------------------------------------------------
# Figure two: panel.
# --------------------------------------------------------------------------

_GRID_L: float = 0.70
"""Left edge of the grid of 2-dimensional cells."""

_GRID_R: float = 4.80
"""Right edge of the grid of 2-dimensional cells."""

_GRID_TOP: float = 4.55
"""Top edge of the grid of 2-dimensional cells."""

_GRID_BOTTOM: float = 2.85
"""Bottom edge of the grid of 2-dimensional cells."""

_HEX_TOP: float = _PANEL_TOP - 0.80
"""Top of the space the hexagon and its vertex labels may occupy, below the
panel's two caption lines."""

_HEX_BOTTOM: float = _PANEL_BOTTOM - 0.15
"""Bottom of that space: the panel floor, plus the gutter above the footer."""

_HEX_LABEL_RADIUS: float = (_HEX_TOP - _HEX_BOTTOM) / 2.0 - 0.10
"""Radius at which the vertex labels sit, clear of the polygon. Sized so the
labelled polygon exactly fills the space left between the captions and the
footer strip — the band is the binding constraint, not the drawing."""

_HEX_RADIUS: float = _HEX_LABEL_RADIUS * 0.84
"""Circumradius of the drawn hexagon, inset within its ring of labels."""

_HEX_CENTER: tuple[float, float] = (9.0, (_HEX_TOP + _HEX_BOTTOM) / 2.0)
"""Center of the drawn hexagon, in drawing inches."""

_MAP_Y: float = 3.55
"""Height of the map arrow, between the two panel halves it joins."""


def _hex_vertex(a: int) -> tuple[float, float]:
    r"""Return the position of :math:`Z_a`, for ``a`` in ``1..6``.

    The six vertices run counterclockwise from the upper left, which is what
    makes the configuration a *positive* :math:`Z`: each ordered
    :math:`3 \times 3` minor of :math:`(1, x_a, y_a)` is twice the signed area
    of the triangle on those three vertices, and every ordered triple of points
    in counterclockwise convex position is positively oriented.
    """
    angle = math.radians(150.0 + (a - 1) * 60.0)
    return (
        _HEX_CENTER[0] + _HEX_RADIUS * math.cos(angle),
        _HEX_CENTER[1] + _HEX_RADIUS * math.sin(angle),
    )


def _cell_grid(ax: Axes) -> None:
    """Draw the 2-dimensional cells of the positive Grassmannian as chips.

    Twenty chips, one per three-element support, with the four that tile filled
    in the ordinal ramp — the same four colors the hexagon's tiles carry, so
    the correspondence reads without a legend line joining them.
    """
    ax.text(
        (_GRID_L + _GRID_R) / 2.0,
        _PANEL_TOP - 0.22,
        r"$Gr^{\geq 0}_{1,6}$ — 63 positroid cells",
        fontsize=13,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        (_GRID_L + _GRID_R) / 2.0,
        _PANEL_TOP - 0.55,
        "twenty of them are 2-dimensional",
        fontsize=12,
        color=style.MIST,
        ha="center",
        va="center",
    )

    col_pitch = (_GRID_R - _GRID_L) / _CELL_COLS
    rows = math.ceil(len(_CELLS_2D) / _CELL_COLS)
    row_pitch = (_GRID_TOP - _GRID_BOTTOM) / rows
    chip_w = col_pitch - 0.14
    chip_h = row_pitch - 0.13

    for index, cell in enumerate(_CELLS_2D):
        row, col = divmod(index, _CELL_COLS)
        x = _GRID_L + col * col_pitch
        y = _GRID_TOP - (row + 1) * row_pitch + (row_pitch - chip_h)
        if cell in _TILING_CELLS:
            tile = _TILING_CELLS.index(cell)
            _box(ax, (x, y, chip_w, chip_h), facecolor=style.ORDINAL_STEPS[tile])
            ax.text(
                x + chip_w / 2.0,
                y + chip_h / 2.0,
                "".join(str(part) for part in cell),
                fontsize=11,
                fontweight="bold",
                color=_ORDINAL_TEXT[tile],
                ha="center",
                va="center",
            )
        else:
            _box(
                ax,
                (x, y, chip_w, chip_h),
                facecolor=style.PAPER,
                edgecolor=style.PARCHMENT,
                linewidth=1.0,
            )
            ax.text(
                x + chip_w / 2.0,
                y + chip_h / 2.0,
                "".join(str(part) for part in cell),
                fontsize=11,
                color=style.MIST,
                ha="center",
                va="center",
            )

    ax.text(
        (_GRID_L + _GRID_R) / 2.0,
        _GRID_BOTTOM - 0.35,
        r"the four cells $\{1,\, i,\, i{+}1\}$ tile the image",
        fontsize=12,
        color=style.INK,
        ha="center",
        va="center",
    )


def _map_arrow(ax: Axes) -> None:
    """Draw the amplituhedron map between the cell grid and the hexagon."""
    _arrow(
        ax,
        (5.25, _MAP_Y),
        (6.95, _MAP_Y),
        color=_MAP_COLOR,
        linewidth=1.6,
    )
    ax.text(
        6.10,
        _MAP_Y + 0.30,
        r"$Y = C \cdot Z$",
        fontsize=15,
        fontweight="bold",
        color=_MAP_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        6.10,
        _MAP_Y - 0.30,
        r"$Z$ positive",
        fontsize=12,
        color=style.MIST,
        ha="center",
        va="center",
    )


def _hexagon(ax: Axes) -> None:
    """Draw the :math:`k=1`, :math:`m=2`, :math:`n=6` amplituhedron and its fan.

    At :math:`k = 1` the amplituhedron is a genuine projective polytope, so the
    picture is the object rather than a schematic of it: six positive columns
    of :math:`Z` are six vertices in convex position, and the four triangles of
    the fan from :math:`Z_1` are the images of the four highlighted cells.
    """
    ax.text(
        _HEX_CENTER[0],
        _PANEL_TOP - 0.22,
        r"$\mathcal{A}_{6,1,2}(Z)$ — the hexagon",
        fontsize=13,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        _HEX_CENTER[0],
        _PANEL_TOP - 0.55,
        "four tiles, six boundary poles",
        fontsize=12,
        color=style.MIST,
        ha="center",
        va="center",
    )

    for index, cell in enumerate(_TILING_CELLS):
        ax.add_patch(
            Polygon(
                [_hex_vertex(a) for a in cell],
                closed=True,
                facecolor=style.ORDINAL_STEPS[index],
                edgecolor=style.PAPER,
                linewidth=1.4,
            )
        )

    vertices = [_hex_vertex(a) for a in range(1, _N_PARTICLES + 1)]
    ax.add_patch(
        Polygon(
            vertices,
            closed=True,
            facecolor="none",
            edgecolor=style.INK,
            linewidth=1.8,
        )
    )
    for a, (vx, vy) in enumerate(vertices, start=1):
        ax.plot([vx], [vy], marker="o", color=style.INK, markersize=6)
        angle = math.radians(150.0 + (a - 1) * 60.0)
        ax.text(
            _HEX_CENTER[0] + _HEX_LABEL_RADIUS * math.cos(angle),
            _HEX_CENTER[1] + _HEX_LABEL_RADIUS * math.sin(angle),
            f"$Z_{a}$",
            fontsize=13,
            color=style.INK,
            ha="center",
            va="center",
        )


# --------------------------------------------------------------------------
# Renderers.
# --------------------------------------------------------------------------


def _render_feynman(step: int) -> Figure:
    """Return the Feynman-expansion diagram revealed up to ``step``."""
    fig, ax = _new_figure()
    _eyebrow(ax, _FEYNMAN_EYEBROW)
    _stages(ax, _FEYNMAN_STAGES, step)
    if step >= _STEP_PANEL_RIGHT:
        _notes(ax, 9.1, "WHAT THE SUM COSTS", _FEYNMAN_NOTES, accent=_COST_COLOR)
    _strip(ax, _FEYNMAN_CHIPS, step)
    # Added last: the plotting axes must sit above the full-figure axes, and
    # `add_axes` stacks in call order.
    if step >= _STEP_PANEL_LEFT:
        _growth_panel(fig)
    return fig


def _render_amplituhedron(step: int) -> Figure:
    """Return the amplituhedron diagram revealed up to ``step``."""
    fig, ax = _new_figure()
    _eyebrow(ax, _AMPLITUHEDRON_EYEBROW)
    _stages(ax, _AMPLITUHEDRON_STAGES, step)
    if step >= _STEP_PANEL_LEFT:
        _cell_grid(ax)
    if step >= _STEP_PANEL_RIGHT:
        _map_arrow(ax)
        _hexagon(ax)
        _notes(
            ax,
            11.3,
            "WHAT THE PICTURE SHOWS",
            _AMPLITUHEDRON_NOTES,
            accent=_MAP_COLOR,
        )
    _strip(ax, _AMPLITUHEDRON_CHIPS, step)
    return fig


@figure(name="scattering-amplitudes")
def scattering_amplitudes(ctx: FigureContext) -> None:
    """Render both amplitude-mechanism diagrams, whole and as build-up sequences.

    Writes ``scattering-feynman`` and ``scattering-amplituhedron`` (the complete
    graphics) plus ``<stem>-NN`` for each reveal step, numbered in presentation
    order. The two finished images share a frame, so they can be placed at the
    same position on consecutive slides.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter and the two figures would
    # stop matching. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for stem, render in (
            ("scattering-feynman", _render_feynman),
            ("scattering-amplituhedron", _render_amplituhedron),
        ):
            for step in range(1, _TOTAL_STEPS + 1):
                fig = render(step)
                ctx.save(fig, f"{stem}-{step:02d}")
                plt.close(fig)

            fig = render(_TOTAL_STEPS)
            ctx.save(fig, stem)
            plt.close(fig)
