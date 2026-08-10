r"""Two matched slide graphics: where a collision calculation's answer comes from.

Renders the same question — *where does the answer to a particle-collision
calculation come from?* — answered twice, as a pair of diagrams meant for
consecutive slides. The first, ``scattering-feynman``, is the textbook route:
fix the rules, draw every distinct way the particles can meet, evaluate each
drawing, and watch nearly all of it cancel. The second,
``scattering-amplituhedron``, is the route of Arkani-Hamed and Trnka: feed one
nonnegative matrix and one kinematic matrix into a single product, and read the
answer off the boundary of the region that comes out.

Both figures are written for readers who work with matrices, determinants, and
exponential blow-up rather than with quantum field theory, so the physics
vocabulary is glossed in place: a *gauge* is an arbitrary bookkeeping choice and
an *internal line* is a quantity that is never observed.

Both share one frame, one four-stage band, and one panel band, in identical
positions — so advancing from the first slide to the second replaces the
mechanism in place while the scaffolding holds still. That substitution *is* the
argument, and it is why these are two figures rather than a split one.

Each carries computed content rather than a cartoon:

* the Feynman figure draws real diagrams. Colour-ordered tree topologies on
  :math:`n` legs are in bijection with the dissections of a convex
  :math:`n`-gon into triangles and quadrilaterals (the cubic and quartic
  Yang-Mills vertices); :func:`dissections` enumerates them and
  :func:`_diagram` draws each one's planar dual, so both galleries are provably
  complete: all three 4-particle diagrams, and all 38 of the 6-particle ones.
  The counts run ``1, 3, 10, 38, 154, 654, 2871, 12925`` for
  :math:`n = 3, \dots, 10` (Elvang and Huang §2.6, verified in the knowledge
  graph's **BCFW Recursion** page and reproduced by
  :func:`dissection_count`), plotted beside the gallery and annotated at the
  one head-to-head comparison that page certifies: the split-helicity 6-gluon
  NMHV amplitude, under a :math:`[1,2\rangle` shift, costs 38 diagrams and 2
  BCFW terms;
* the amplituhedron figure draws the map itself at :math:`k = 1`, :math:`m = 2`,
  :math:`n = 6`, where the amplituhedron genuinely *is* the interior of a convex
  hexagon. Six vertices in counterclockwise convex position are a positive
  :math:`Z` (every ordered :math:`3 \times 3` minor is twice a positive signed
  area), the twenty 2-dimensional cells of :math:`Gr^{\geq 0}_{1,6}` are the
  three-element supports — at :math:`k = 1`, simply which three of the six
  entries of the single row are nonzero — and the four of the form
  :math:`\{1, i, i{+}1\}` map to the four triangles of the fan from
  :math:`Z_1`, matching the proved :math:`m = 2` tile count
  :math:`\binom{n-2}{k} = 4`. Cells and tiles are drawn in the same four
  ordinal colors, so the correspondence is visible rather than asserted.

Neither graphic carries a title, subtitle, or summary chrome of any kind: these
are placed in a slide deck, and the slide supplies all of that. Both are saved
cropped to the drawing itself so they can be placed and scaled freely, and both
write a build-up sequence — one image per reveal step, each adding the next
piece — into the same frame as the whole, so the steps can be stacked on one
slide without anything shifting between them::

    1  stage 1        5  panel, left
    2  stage 2        6  panel, right
    3  stage 3
    4  stage 4

Regenerate with ``just figure scattering-amplitudes``.
"""

import functools
import itertools
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

_PANEL_BOTTOM: float = 2.20
"""Bottom edge of the band holding each figure's computed panel."""

_PANEL_TOP: float = 6.10
"""Top edge of the panel band."""

_STAGE_BOTTOM: float = 6.45
"""Bottom edge of the stage row."""

_STAGE_H: float = 2.70
"""Height of a pipeline stage box."""

_STAGE_TOP: float = _STAGE_BOTTOM + _STAGE_H
"""Top edge of the stage row, and the top of the drawing.

Nothing sits above it: the deck supplies the title, so the figure starts at
its first piece of content.
"""

_STAGE_GAP: float = 0.5
"""Horizontal gap between consecutive stage boxes, holding the flow arrow."""

_STAGE_W: float = (_MARGIN_R - _MARGIN_L - 3.0 * _STAGE_GAP) / 4.0
"""Width of a pipeline stage box: four boxes and three gaps span the margins."""

_CONTENT_TOP: float = _STAGE_TOP
"""Top of the drawing."""

_CONTENT_BOTTOM: float = 1.85
"""Bottom of the drawing.

Set by the thumbnail block, which hangs below :data:`_PANEL_BOTTOM` and is the
lowest thing either figure draws. The amplituhedron figure stops higher and
carries the difference as white space — the frame is shared, so one of the two
always has slack somewhere.
"""

_FRAME_LEFT: float = _MARGIN_L - _PAD
"""Left limit of the full-figure axes, in inches."""

_FRAME_RIGHT: float = _MARGIN_R + _PAD
"""Right limit of the full-figure axes, in inches."""

_FRAME_BOTTOM: float = _CONTENT_BOTTOM - _PAD
"""Bottom limit of the full-figure axes, in inches."""

_FRAME_TOP: float = _CONTENT_TOP + _PAD
"""Top limit of the full-figure axes, in inches."""

_NOTES_L: float = 11.95
"""Left edge of the reading column, the same in both figures."""

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

_TOTAL_STEPS: int = 6
"""Number of images in each build-up sequence; the last is the whole diagram."""

_STEP_PANEL_LEFT: int = 5
"""Reveal step adding the left half of the panel band."""

_STEP_PANEL_RIGHT: int = 6
"""Reveal step adding the right half of the panel band."""


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


# --------------------------------------------------------------------------
# Figure one: the Feynman expansion.
# --------------------------------------------------------------------------

_FEYNMAN_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="THE RULES",
        body=(
            "the theory says which particles",
            "can meet, how strongly, and how",
            "they travel in between. it also",
            "makes one arbitrary choice.",
        ),
        tag="that choice is called a gauge",
    ),
    _Stage(
        number=2,
        title="ENUMERATE",
        body=(
            "draw every distinct way the",
            "particles can meet and split on",
            "the way from the inputs to the",
            "outputs. one drawing per way.",
        ),
        tag="1, 3, 10, 38, 154, 654, ...",
    ),
    _Stage(
        number=3,
        title="EVALUATE",
        body=(
            "turn each drawing into one",
            "algebraic term. its internal",
            "lines are not real particles,",
            "and carry the arbitrary choice.",
        ),
        tag="no single term is measurable",
    ),
    _Stage(
        number=4,
        title="CANCEL",
        body=(
            "add the terms up. nearly all of",
            "it cancels: the arbitrary choice",
            "vanishes, and a huge sum",
            "collapses to a few pieces.",
        ),
        # Deliberately not "38 drawings in, 2 pieces out": 2 is the BCFW term
        # count, which this route does not produce. Claiming it here would put
        # the second method's payoff on the first method's last stage.
        tag="almost none of the work survives",
    ),
)
"""The textbook pipeline, left to right."""

_PARTICLE_COUNTS: tuple[int, ...] = (3, 4, 5, 6, 7, 8, 9, 10)
"""Number of external particles for each plotted diagram count."""

_HIGHLIGHT_PARTICLES: int = 6
r"""The leg count called out on the curve.

The one head-to-head comparison the knowledge graph certifies: 38 diagrams
against 2 BCFW terms, for the split-helicity amplitude
:math:`A_6[1^-2^-3^-4^+5^+6^+]` under a :math:`[1,2\rangle` shift. The shift is
part of the claim, not decoration - the same amplitude under :math:`[2,1\rangle`
has three terms, so a term count is never a property of the amplitude alone.
"""

_GALLERY_FULL_N: int = 4
"""Particle count whose diagrams are drawn in full, large and labelled."""

_GALLERY_SAMPLE_N: int = 6
"""Particle count drawn beneath it in full, at thumbnail size."""

_FEYNMAN_NOTES: tuple[_Note, ...] = (
    _Note(
        lead="the parts are messier than the whole",
        detail="no one drawing is measurable; their sum is",
    ),
    _Note(
        lead="every term needs an arbitrary choice",
        detail="the answer they add up to does not",
    ),
    _Note(
        lead="catastrophic cancellation, by design",
        # Naive power counting gives individual diagrams z, z^2, or z^3
        # depending on the helicity case, while an adjacent shift of the sum
        # falls as 1/z (Arkani-Hamed and Kaplan 2008, recorded on the
        # knowledge graph's **BCFW Recursion** page). Written as a range: the
        # three powers are three separate cases, not three terms of one sum.
        detail=r"terms diverge up to $z^3$; the sum falls as $1/z$",
    ),
    _Note(
        lead="6 particles: 38 drawings, 2 terms",
        detail="same answer, from a method with no sum",
        accent=True,
    ),
)
"""The reading column beside the gallery and the growth curve."""

# --------------------------------------------------------------------------
# Figure two: the amplituhedron.
# --------------------------------------------------------------------------

_AMPLITUHEDRON_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="POSITIVE DATA",
        body=(
            r"take a $k \times n$ matrix $C$ whose",
            r"ordered $k \times k$ determinants are",
            r"all $\geq 0$, plus a matrix $Z$ built",
            "from the particles' momenta.",
        ),
        tag="positivity is the whole input",
    ),
    _Stage(
        number=2,
        title="ONE MAP",
        body=(
            r"multiply: $Y = C \cdot Z$. let $C$ range",
            "over every such matrix and the",
            r"outputs $Y$ fill in one connected",
            "region: a single shape.",
        ),
        tag="a matrix product, nothing more",
    ),
    _Stage(
        number=3,
        title="TILE",
        body=(
            "that region is covered exactly",
            "once by finitely many tiles, the",
            "images of finitely many pieces",
            r"of the space of matrices $C$.",
        ),
        tag="the tiles are the answer's terms",
    ),
    _Stage(
        number=4,
        title="READ IT OFF",
        body=(
            "a region like this carries one",
            "distinguished function, fixed by",
            "where its boundaries are.",
            "that function is the answer.",
        ),
        tag="one shape, one answer",
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
cells are exactly the three-element subsets of :math:`[n]`. At :math:`k = 1`
that reads without any Grassmannian vocabulary at all: :math:`C` is one
nonnegative row of six numbers, and a piece of that space is fixed by which
three of the six entries are nonzero.
"""

_TILING_CELLS: tuple[tuple[int, ...], ...] = tuple(
    (1, i, i + 1) for i in range(2, _N_PARTICLES)
)
"""The four cells whose images tile the hexagon: the fan from :math:`Z_1`."""

_CELL_COLS: int = 5
"""Columns in the grid of 2-dimensional cells."""

_AMPLITUHEDRON_NOTES: tuple[_Note, ...] = (
    _Note(
        lead="one region, not a sum of terms",
        detail="the answer is a property of a single shape",
    ),
    _Note(
        lead="four tiles, however you tile it",
        detail=r"$\binom{n-2}{k} = 4$ here; the count is intrinsic",
    ),
    _Note(
        lead="the six edges are the only poles",
        detail="exactly one function has just these poles",
    ),
    _Note(
        lead="same answer, no cancellation",
        detail="nothing large is built only to cancel away",
        accent=True,
    ),
)
"""The reading column beside the worked hexagon."""


# --------------------------------------------------------------------------
# Colour-ordered tree diagrams, enumerated.
# --------------------------------------------------------------------------

_MIN_SIDES: int = 3
"""Fewest sides a polygon can have; below it there is nothing left to cut."""

_CELL_SIDES: tuple[int, ...] = (3, 4)
"""Cell sizes a dissection may use: the cubic and quartic Yang-Mills vertices."""

_CELLS_PER_DIAGONAL: int = 2
"""Cells meeting along an internal line. A polygon side meets exactly one."""

type _Cell = tuple[int, ...]
"""One cell of a polygon dissection: its vertices, in boundary order."""

type _Dissection = tuple[_Cell, ...]
"""One dissection of a convex polygon into triangles and quadrilaterals.

Dual to one colour-ordered tree diagram: a cell is a vertex of the diagram, a
polygon edge is an external leg, and a diagonal is an internal line.
"""


def _cell_choices(last: int) -> list[tuple[int, ...]]:
    """Return the interior index picks for the cell on a polygon's root edge.

    The cell containing the root edge runs from index ``0`` to index ``last``,
    so a cell of ``s`` vertices contributes ``s - 2`` intermediate indices
    drawn in order from ``1 .. last - 1``, once per allowed ``s`` in
    :data:`_CELL_SIDES`.
    """
    return [
        interior
        for cell_sides in _CELL_SIDES
        for interior in itertools.combinations(range(1, last), cell_sides - 2)
    ]


def dissections(vertices: tuple[int, ...]) -> list[_Dissection]:
    """Return every dissection of a convex polygon into triangles and quads.

    Args:
        vertices: The polygon's vertices in boundary order. The pair
            ``(vertices[0], vertices[-1])`` is the root edge: the recursion
            splits on which cell contains it, then recurses into the
            sub-polygons cut off by that cell's other sides.

    Returns:
        One entry per dissection, each a tuple of cells. Cells appear in
        recursion order, which is unspecified but deterministic.
    """
    if len(vertices) < _MIN_SIDES:
        return [()]
    last = len(vertices) - 1
    found: list[_Dissection] = []
    for interior in _cell_choices(last):
        picks = (0, *interior, last)
        cell = tuple(vertices[index] for index in picks)
        sub_polygons = [
            vertices[a : b + 1] for a, b in itertools.pairwise(picks) if b > a + 1
        ]
        found.extend(
            (cell, *itertools.chain.from_iterable(combination))
            for combination in itertools.product(
                *(dissections(sub) for sub in sub_polygons)
            )
        )
    return found


@functools.cache
def dissection_count(sides: int) -> int:
    """Return how many such dissections a convex ``sides``-gon admits.

    The same recursion as :func:`dissections`, counted rather than
    materialized: the number depends only on the polygon's size, so it
    memoizes on ``sides`` and stays cheap out to the largest leg count
    plotted. A polygon of fewer than three sides is a bare edge, which admits
    exactly one (empty) dissection.
    """
    if sides < _MIN_SIDES:
        return 1
    last = sides - 1
    total = 0
    for interior in _cell_choices(last):
        picks = (0, *interior, last)
        ways = 1
        for a, b in itertools.pairwise(picks):
            if b > a + 1:
                ways *= dissection_count(b - a + 1)
        total += ways
    return total


_TREE_DIAGRAMS: tuple[int, ...] = tuple(
    dissection_count(count) for count in _PARTICLE_COUNTS
)
"""Colour-ordered tree Feynman diagrams for :data:`_PARTICLE_COUNTS` particles.

One per topology built from the cubic and quartic Yang-Mills vertices, so one
per dissection of a convex :math:`n`-gon into triangles and quadrilaterals.
Computed rather than transcribed, so the plotted curve and the drawn gallery
cannot disagree; ``tests/test_figures_scattering.py`` pins the result against
the published sequence ``1, 3, 10, 38, 154, 654, 2871, 12925`` (the first five
are Elvang-Huang §2.6, and the whole sequence is recorded as verified on the
knowledge graph's **BCFW Recursion** page). Memoized and linear in the largest
leg count, so evaluating it at import costs nothing.
"""


def _leg_angle(index: int, legs: int) -> float:
    """Return the outward direction of external leg ``index``, in radians.

    Legs run counterclockwise from the upper left, which puts the incoming
    half of the process on the left and the outgoing half on the right — the
    orientation colour-ordered diagrams are conventionally drawn in. Leg
    ``index`` is the polygon edge from vertex ``index`` to vertex
    ``index + 1``, so its direction bisects those two vertices'.
    """
    return math.radians(90.0 + 360.0 * (index + 0.5) / legs)


def _polygon_vertex(index: int, legs: int) -> tuple[float, float]:
    """Return polygon vertex ``index`` on the unit circle, in diagram units."""
    angle = math.radians(90.0 + 360.0 * index / legs)
    return (math.cos(angle), math.sin(angle))


def _centroid(points: list[tuple[float, float]]) -> tuple[float, float]:
    """Return the arithmetic mean of ``points``."""
    count = len(points)
    return (
        sum(x for x, _ in points) / count,
        sum(y for _, y in points) / count,
    )


_WAVE_SAMPLES: int = 96
"""Points used to draw one wavy line segment."""

_WAVE_LENGTH: float = 0.30
"""Target wavelength of a drawn wave, in diagram units."""

_WAVE_AMPLITUDE: float = 0.05
"""Peak transverse offset of a drawn wave, in diagram units."""


def _wave(
    ax: Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    origin: tuple[float, float],
    scale: float,
    color: str,
    linewidth: float,
) -> None:
    """Draw one wavy line between two points of a diagram's local frame.

    Gluon lines are conventionally drawn as coils; a transverse sine is the
    cheap stand-in, and at this size reads the same. A whole number of periods
    is used so both endpoints sit exactly on the straight line between
    ``start`` and ``end``, and the vertices they meet stay sharp.

    Args:
        ax: The full-figure axes.
        start: Segment start, in diagram units.
        end: Segment end, in diagram units.
        origin: Where the diagram's center falls, in drawing inches.
        scale: Diagram units per drawing inch.
        color: Stroke color.
        linewidth: Stroke width in points.
    """
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return
    periods = max(1, round(length / _WAVE_LENGTH))
    ux, uy = dx / length, dy / length
    xs: list[float] = []
    ys: list[float] = []
    for sample in range(_WAVE_SAMPLES + 1):
        t = sample / _WAVE_SAMPLES
        offset = _WAVE_AMPLITUDE * math.sin(2.0 * math.pi * periods * t)
        xs.append(origin[0] + scale * (start[0] + dx * t - uy * offset))
        ys.append(origin[1] + scale * (start[1] + dy * t + ux * offset))
    ax.plot(xs, ys, color=color, linewidth=linewidth, solid_capstyle="round", zorder=2)


_LEG_RADIUS: float = 1.34
"""Radius at which an external leg ends, in diagram units."""

_LEG_LABEL_RADIUS: float = 1.60
"""Radius at which an external leg's number sits, in diagram units."""

_DIAGRAM_EXTENT_LABELLED: float = _LEG_LABEL_RADIUS + 0.10
"""Half-width of the box a labelled diagram occupies, in diagram units."""

_DIAGRAM_EXTENT_BARE: float = _LEG_RADIUS + 0.06
"""The same for an unlabelled one, which needs no room for leg numbers.

Kept separate so thumbnails scale to the drawing rather than to the empty
ring where the numbers would have gone.
"""


@dataclass(frozen=True, slots=True)
class _DiagramStyle:
    """How one tree diagram is drawn: everything that varies with its size.

    The gallery draws diagrams at exactly two sizes, so these travel together
    as a pair of presets rather than as loose arguments at each call.
    """

    height: float
    """Full height of the box the diagram occupies, in drawing inches."""

    wavy: bool
    """Draw lines as waves rather than straight segments.

    Reads as a Feynman diagram at gallery size and as noise at thumbnail size.
    """

    labels: bool
    """Number the external legs."""

    linewidth: float
    """Stroke width in points."""

    dot_size: float
    """Marker size of an interaction vertex, in points."""

    label_size: float
    """Font size of a leg number."""


def _diagram(
    ax: Axes,
    dissection: _Dissection,
    legs: int,
    *,
    center: tuple[float, float],
    spec: _DiagramStyle,
) -> None:
    """Draw one colour-ordered tree diagram: the planar dual of ``dissection``.

    Each cell becomes an interaction vertex at its centroid, each polygon edge
    an external leg leaving through that edge's midpoint and continuing
    radially outward, and each diagonal an internal line joining the two cells
    that share it. Routing legs through their own edge's midpoint keeps every
    stroke inside the cell it belongs to until it leaves the polygon, so the
    drawing is planar by construction — no topology-specific layout needed.

    External legs are drawn in slate and internal lines in indianred, because
    the distinction is the whole point of the first figure: the internal lines
    are the parts that are never observed and that later cancel.

    Args:
        ax: The full-figure axes.
        dissection: The dissection to dualize.
        legs: Number of external legs, i.e. sides of the dissected polygon.
        center: Where the diagram's center falls, in drawing inches.
        spec: Size-dependent drawing settings.
    """
    extent = _DIAGRAM_EXTENT_LABELLED if spec.labels else _DIAGRAM_EXTENT_BARE
    scale = spec.height / (2.0 * extent)
    vertices = [_polygon_vertex(index, legs) for index in range(legs)]
    centroids = {cell: _centroid([vertices[v] for v in cell]) for cell in dissection}

    sides: dict[frozenset[int], list[_Cell]] = {}
    for cell in dissection:
        for a, b in itertools.pairwise((*cell, cell[0])):
            sides.setdefault(frozenset((a, b)), []).append(cell)

    def stroke(
        start: tuple[float, float], end: tuple[float, float], color: str
    ) -> None:
        """Draw one line of the diagram, wavy or straight as configured."""
        if spec.wavy:
            _wave(
                ax,
                start,
                end,
                origin=center,
                scale=scale,
                color=color,
                linewidth=spec.linewidth,
            )
        else:
            ax.plot(
                [center[0] + scale * start[0], center[0] + scale * end[0]],
                [center[1] + scale * start[1], center[1] + scale * end[1]],
                color=color,
                linewidth=spec.linewidth,
                solid_capstyle="round",
                zorder=2,
            )

    for index in range(legs):
        following = (index + 1) % legs
        cell = sides[frozenset((index, following))][0]
        midpoint = _centroid([vertices[index], vertices[following]])
        angle = _leg_angle(index, legs)
        outward = (_LEG_RADIUS * math.cos(angle), _LEG_RADIUS * math.sin(angle))
        stroke(centroids[cell], midpoint, style.SLATE)
        stroke(midpoint, outward, style.SLATE)
        if spec.labels:
            ax.text(
                center[0] + scale * _LEG_LABEL_RADIUS * math.cos(angle),
                center[1] + scale * _LEG_LABEL_RADIUS * math.sin(angle),
                str(index + 1),
                fontsize=spec.label_size,
                color=style.INK,
                ha="center",
                va="center",
            )

    for shared in sides.values():
        if len(shared) == _CELLS_PER_DIAGONAL:
            stroke(centroids[shared[0]], centroids[shared[1]], _COST_COLOR)

    for point in centroids.values():
        ax.plot(
            [center[0] + scale * point[0]],
            [center[1] + scale * point[1]],
            marker="o",
            color=style.INK,
            markersize=spec.dot_size,
            zorder=3,
        )


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


def _caption(
    ax: Axes,
    x: float,
    y: float,
    lead: str,
    detail: str | None = None,
    *,
    ha: str = "left",
) -> None:
    """Label a piece of the panel band: a bold line, optionally a muted one.

    These are content labels, not figure titles — the deck supplies the title
    and the subtitle, so nothing here restates what the slide already says.
    """
    ax.text(
        x, y, lead, fontsize=13, fontweight="bold", color=style.INK, ha=ha, va="center"
    )
    if detail is not None:
        ax.text(x, y - 0.30, detail, fontsize=11, color=style.MIST, ha=ha, va="center")


def _notes(
    ax: Axes,
    header: str,
    notes: tuple[_Note, ...],
    *,
    accent: str,
) -> None:
    """Draw a panel's reading column, headed and ruled, at :data:`_NOTES_L`.

    Args:
        ax: The full-figure axes.
        header: Column heading, drawn above a rule.
        notes: Lead/detail pairs, top to bottom.
        accent: Color for the notes flagged ``accent``. The blog palette reads
            indianred as "falsified", so only the cost figure may use it; the
            amplituhedron column accents in the teal that already carries its
            map, and the difference is itself part of the contrast.
    """
    ax.text(
        _NOTES_L,
        _PANEL_TOP - 0.10,
        header,
        fontsize=14,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    ax.plot(
        [_NOTES_L, _MARGIN_R],
        [_PANEL_TOP - 0.35, _PANEL_TOP - 0.35],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )
    for index, note in enumerate(notes):
        top = _PANEL_TOP - 0.80 - index * 0.80
        ax.text(
            _NOTES_L,
            top,
            note.lead,
            fontsize=13,
            fontweight="bold",
            color=accent if note.accent else style.INK,
            va="center",
        )
        ax.text(
            _NOTES_L,
            top - 0.32,
            note.detail,
            fontsize=12,
            color=accent if note.accent else style.MIST,
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

_GALLERY_L: float = _MARGIN_L
"""Left edge of the diagram gallery."""

_GALLERY_R: float = 7.30
"""Right edge of the gallery.

Set by where the growth panel's label begins, not by the page margin: the two
share the band, and the thumbnail block is the piece that has to give way.
"""

_GALLERY_HEADING_GAP: float = 0.28
"""Clearance between a row's heading and the top of its diagram boxes."""

_FULL_STYLE: _DiagramStyle = _DiagramStyle(
    height=1.70,
    wavy=True,
    labels=True,
    linewidth=1.5,
    dot_size=4.5,
    label_size=10.0,
)
"""The gallery's large, labelled, wavy diagrams: one particle process, legible."""

_GALLERY_FULL_Y: float = (
    _PANEL_TOP - 0.10 - _FULL_STYLE.height / 2.0 - _GALLERY_HEADING_GAP
)
"""Vertical center of the row of fully drawn diagrams.

Hung from the panel heading rather than set outright, so the row sits the same
:data:`_GALLERY_HEADING_GAP` below its label as the thumbnail block does below
its own. The heading is pinned to :data:`_PANEL_TOP` because it aligns with the
notes column opposite, so the row has to move to close the gap, not the label.
"""

_GALLERY_FULL_PITCH: float = 2.0
"""Horizontal spacing between fully drawn diagrams."""

_GALLERY_SAMPLE_TOP: float = 3.45
"""Top edge of the thumbnail block, clear of the labelled row's leg numbers."""

_GALLERY_SAMPLE_BOTTOM: float = _CONTENT_BOTTOM
"""Bottom edge of the thumbnail block: the floor of the drawing itself.

The block is the lowest element in either figure, so nothing below it competes
for the space — it simply runs out at the frame.
"""

_GALLERY_SAMPLE_COLS: int = 13
"""Thumbnails per row.

Chosen to leave the last row nearly full: 13 splits the 38 as 13/13/12, where
14 would strand 10 alone.
"""

_GALLERY_SAMPLE_ROWS: int = math.ceil(
    dissection_count(_GALLERY_SAMPLE_N) / _GALLERY_SAMPLE_COLS
)
"""Rows needed to hold every diagram at :data:`_GALLERY_SAMPLE_COLS` per row."""

_GALLERY_SAMPLE_PITCH: float = min(
    (_GALLERY_R - _GALLERY_L) / _GALLERY_SAMPLE_COLS,
    (_GALLERY_SAMPLE_TOP - _GALLERY_SAMPLE_BOTTOM) / _GALLERY_SAMPLE_ROWS,
)
"""Center-to-center spacing of thumbnails, equal on both axes.

The block shows every diagram, so the cell is whatever divides the space
available — never a tuned number. Taking the smaller of the two axes keeps the
cell square, which matters because the diagrams are drawn on a circle and would
otherwise collide along whichever axis was tighter.
"""

_GALLERY_SAMPLE_PACKING: float = 0.90
"""Fraction of a cell a thumbnail fills, leaving a gutter between neighbours.

Colour-ordered legs run almost to the edge of a diagram's box, so thumbnails
drawn at the full cell size would touch and read as one continuous mesh.
"""


_THUMB_STYLE: _DiagramStyle = _DiagramStyle(
    height=_GALLERY_SAMPLE_PITCH * _GALLERY_SAMPLE_PACKING,
    wavy=False,
    labels=False,
    linewidth=0.9,
    dot_size=2.0,
    label_size=8.0,
)
"""The block beneath it: the same object, small enough that volume is the point."""


def _gallery(ax: Axes) -> None:
    """Draw every 4-particle diagram, then every 6-particle one.

    The first row is complete and labelled, so a reader who has never seen a
    Feynman diagram can read one: two of the three drawings route the process
    through an internal line, splitting the legs as ``(1,4)(2,3)`` and
    ``(1,2)(3,4)`` — the only two splits a fixed cyclic order admits — and the
    third has all four meet at a single vertex. The block below is the same
    object at six particles, drawn in full at thumbnail size: no thumbnail is
    meant to be read, and the eye is meant to give up counting, which is the
    argument the panel is making.

    Both are exhaustive, so nothing here is a curated sample — the block is
    what the enumerator returns, in the order it returns it.

    No colour key is drawn. Slate still means an external leg and indianred an
    internal line (see :func:`_diagram`), but the presenter says so aloud, so
    printing it here would only duplicate the narration.
    """
    _caption(
        ax,
        _GALLERY_L,
        _PANEL_TOP - 0.10,
        f"{_GALLERY_FULL_N} particles: all "
        f"{dissection_count(_GALLERY_FULL_N)} drawings",
    )
    diagrams = dissections(tuple(range(_GALLERY_FULL_N)))
    for index, dissection in enumerate(diagrams):
        _diagram(
            ax,
            dissection,
            _GALLERY_FULL_N,
            center=(
                _GALLERY_L + _GALLERY_FULL_PITCH * (index + 0.5),
                _GALLERY_FULL_Y,
            ),
            spec=_FULL_STYLE,
        )

    found = dissections(tuple(range(_GALLERY_SAMPLE_N)))
    _caption(
        ax,
        _GALLERY_L,
        _GALLERY_SAMPLE_TOP + _GALLERY_HEADING_GAP,
        f"{_GALLERY_SAMPLE_N} particles: all {len(found)} drawings",
    )
    for index, dissection in enumerate(found):
        row, column = divmod(index, _GALLERY_SAMPLE_COLS)
        _diagram(
            ax,
            dissection,
            _GALLERY_SAMPLE_N,
            center=(
                _GALLERY_L + _GALLERY_SAMPLE_PITCH * (column + 0.5),
                _GALLERY_SAMPLE_TOP - _GALLERY_SAMPLE_PITCH * (row + 0.5),
            ),
            spec=_THUMB_STYLE,
        )


def _growth_panel(fig: Figure, ax: Axes) -> None:
    """Plot the enumeration cost: drawings against the number of particles.

    The vertical scale is logarithmic because the sequence is — the point of
    the panel is that the number of parts to write down outruns the answer's
    length, not the particular values along the way.
    """
    _caption(
        ax,
        7.55,
        _PANEL_TOP - 0.10,
        "and it keeps going",
        "drawings per particle count, log scale",
    )
    plot = _inset(fig, (7.90, 2.75, 3.45, 2.35))
    plot.semilogy(
        _PARTICLE_COUNTS,
        _TREE_DIAGRAMS,
        color=style.SLATE,
        marker="o",
        markerfacecolor=style.PAPER,
        markeredgewidth=2.0,
    )

    index = _PARTICLE_COUNTS.index(_HIGHLIGHT_PARTICLES)
    plot.plot(
        [_HIGHLIGHT_PARTICLES],
        [_TREE_DIAGRAMS[index]],
        marker="o",
        color=_COST_COLOR,
        markersize=9,
    )
    plot.annotate(
        "38",
        xy=(_HIGHLIGHT_PARTICLES, _TREE_DIAGRAMS[index]),
        xytext=(6.5, 2.2),
        fontsize=12,
        color=_COST_COLOR,
        arrowprops={"arrowstyle": "-", "color": _COST_COLOR, "linewidth": 1.0},
    )
    plot.annotate(
        f"{_TREE_DIAGRAMS[-1]:,}",
        xy=(_PARTICLE_COUNTS[-1], _TREE_DIAGRAMS[-1]),
        xytext=(-6, 4),
        textcoords="offset points",
        fontsize=12,
        color=style.SLATE,
        ha="right",
    )

    plot.set_xlabel("particles", fontsize=12)
    plot.set_xticks(list(_PARTICLE_COUNTS))
    plot.set_xlim(2.6, 10.4)
    plot.tick_params(labelsize=11)


# --------------------------------------------------------------------------
# Figure two: panel.
# --------------------------------------------------------------------------

_GRID_L: float = 0.60
"""Left edge of the grid of 2-dimensional cells."""

_GRID_R: float = 4.90
"""Right edge of the grid of 2-dimensional cells."""

_GRID_TOP: float = 5.10
"""Top edge of the grid of 2-dimensional cells."""

_GRID_BOTTOM: float = 2.95
"""Bottom edge of the grid of 2-dimensional cells."""

_HEX_TOP: float = _PANEL_TOP - 0.85
"""Top of the space the hexagon and its vertex labels may occupy, below the
panel's two caption lines."""

_HEX_BOTTOM: float = _PANEL_BOTTOM - 0.05
"""Bottom of that space: the panel floor, plus the gutter above the footer."""

_HEX_LABEL_RADIUS: float = (_HEX_TOP - _HEX_BOTTOM) / 2.0 - 0.10
"""Radius at which the vertex labels sit, clear of the polygon. Sized so the
labelled polygon exactly fills the space left between the captions and the
footer strip — the band is the binding constraint, not the drawing."""

_HEX_RADIUS: float = _HEX_LABEL_RADIUS * 0.84
"""Circumradius of the drawn hexagon, inset within its ring of labels."""

_HEX_CENTER: tuple[float, float] = (8.95, (_HEX_TOP + _HEX_BOTTOM) / 2.0)
"""Center of the drawn hexagon, in drawing inches."""

_MAP_Y: float = 3.70
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

    At :math:`k = 1` the caption can say what a "cell" is without leaving
    linear algebra: :math:`C` is a single nonnegative row, and a cell records
    which of its entries are nonzero.
    """
    _caption(
        ax,
        (_GRID_L + _GRID_R) / 2.0,
        _PANEL_TOP - 0.10,
        r"$C$ is one nonnegative row of 6",
        "a piece of that space: which 3 entries are nonzero",
        ha="center",
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
        "four of the twenty pieces tile the shape",
        fontsize=12,
        color=style.INK,
        ha="center",
        va="center",
    )


def _map_arrow(ax: Axes) -> None:
    """Draw the amplituhedron map between the cell grid and the hexagon."""
    _arrow(
        ax,
        (5.30, _MAP_Y),
        (6.85, _MAP_Y),
        color=_MAP_COLOR,
        linewidth=1.6,
    )
    ax.text(
        6.075,
        _MAP_Y + 0.30,
        r"$Y = C \cdot Z$",
        fontsize=15,
        fontweight="bold",
        color=_MAP_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        6.075,
        _MAP_Y - 0.30,
        "one matrix product",
        fontsize=11,
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
    _caption(
        ax,
        _HEX_CENTER[0],
        _PANEL_TOP - 0.10,
        "the region the outputs fill",
        "four tiles, six edges, no sum",
        ha="center",
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
    _stages(ax, _FEYNMAN_STAGES, step)
    if step >= _STEP_PANEL_LEFT:
        _gallery(ax)
    if step >= _STEP_PANEL_RIGHT:
        _notes(ax, "WHAT THIS COSTS", _FEYNMAN_NOTES, accent=_COST_COLOR)
    # Added last: the plotting axes must sit above the full-figure axes, and
    # `add_axes` stacks in call order.
    if step >= _STEP_PANEL_RIGHT:
        _growth_panel(fig, ax)
    return fig


def _render_amplituhedron(step: int) -> Figure:
    """Return the amplituhedron diagram revealed up to ``step``."""
    fig, ax = _new_figure()
    _stages(ax, _AMPLITUHEDRON_STAGES, step)
    if step >= _STEP_PANEL_LEFT:
        _cell_grid(ax)
    if step >= _STEP_PANEL_RIGHT:
        _map_arrow(ax)
        _hexagon(ax)
        _notes(ax, "WHAT THIS BUYS", _AMPLITUHEDRON_NOTES, accent=_MAP_COLOR)
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
