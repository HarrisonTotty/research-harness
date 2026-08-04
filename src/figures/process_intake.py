"""Meta-figure: the harness's Intake & Refinement phase, as a slide graphic.

Renders the first phase of the research process described in the repository
README — how papers and primary sources become a de-duplicated, cross-linked
concept graph in Logseq — as a single 16:9 diagram sized for a presentation
slide. The pipeline stages, their subagents, and their working files are the
ones the ``add-logseq-topic`` skill actually runs; the two feedback edges (the
audit gate on publication, and red links seeding the next intake) are drawn
because they are what make the phase a loop rather than a line.

The graphic carries no title, standfirst, or summary chrome — the slide it goes
on supplies those — and is saved cropped to the drawing itself rather than to
the 16:9 frame, so it can be placed and scaled freely on the slide.

Besides the complete diagram, the module writes a build-up sequence — one image
per reveal step, each adding the next piece of the process — for walking an
audience through it a stage at a time. A step reveals a stage together with
whatever feeds it, so step *n* is stage *n* until the last step closes the
loop. Every step is rendered into the same frame as the whole, so the images
can be stacked on one slide (or advanced through) without anything shifting
between them.

First of a series of meta-figures, one per README process phase. Regenerate
with ``just figure process-intake``.
"""

from dataclasses import dataclass

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path

from figures import style
from figures.cli import FigureContext, figure

_PAD: float = 0.2
"""Border kept around the drawing on every side, in inches."""

_MARGIN_L: float = 0.6
"""Left edge of every full-width element."""

_MARGIN_R: float = 17.4
"""Right edge of every full-width element."""

_STAGE_W: float = 3.0
"""Width of a pipeline stage box."""

_STAGE_H: float = 2.5
"""Height of a pipeline stage box."""

_STAGE_GAP: float = 0.45
"""Horizontal gap between consecutive stage boxes, holding the flow arrow."""

_STAGE_TOP: float = 6.0
"""Top edge of the stage row."""

_STAGE_BOTTOM: float = _STAGE_TOP - _STAGE_H
"""Bottom edge of the stage row."""

_ARC_CLEARANCE: float = 0.55
"""Headroom above the stage row for the bow of the audit arc."""

_FILE_TOP: float = 3.0
"""Top edge of the working-file chips, below the stage row."""

_FILE_H: float = 0.9
"""Height of a working-file chip."""

_RETURN_Y: float = 1.7
"""Height of the horizontal run of the red-link return path."""

_CONTENT_TOP: float = _STAGE_TOP + _ARC_CLEARANCE
"""Top of the drawing: the apex of the audit arc above the stage row."""

_CONTENT_BOTTOM: float = _RETURN_Y
"""Bottom of the drawing: the horizontal run of the return path."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_AUDIT_COLOR: str = style.CATEGORICAL[5]
"""Indianred — the blog's falsification hue, carrying the audit gate."""

_RETURN_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the red-link loop, distinct from the audit loop."""


@dataclass(frozen=True, slots=True)
class _Stage:
    """One box in the intake pipeline.

    Body lines are pre-wrapped rather than flowed: at this size the wrap
    points are a layout decision, not something to leave to a text engine.
    """

    number: int
    title: str
    body: tuple[str, ...]
    tool: str


_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="SURVEY",
        body=(
            "search the graph for the",
            "topic and its synonyms —",
            "extend, never fork",
        ),
        tool="search_logseq",
    ),
    _Stage(
        number=2,
        title="RESEARCH",
        body=(
            "one subagent per source;",
            "original sources first,",
            "every claim attributed",
        ),
        tool="source-reader agents",
    ),
    _Stage(
        number=3,
        title="DRAFT",
        body=(
            "build the page from the",
            "fact map alone: theorems,",
            "examples, attributions",
        ),
        tool="house page template",
    ),
    _Stage(
        number=4,
        title="AUDIT",
        body=(
            "a second agent diffs draft",
            "against facts; unsupported",
            "claims never ship",
        ),
        tool="draft-auditor · until CLEAN",
    ),
    _Stage(
        number=5,
        title="PUBLISH",
        body=(
            "transcribe, then re-read",
            "the live page against the",
            "draft; resolve every link",
        ),
        tool="logseq mcp · link check",
    ),
)
"""The pipeline, left to right."""

_STEP_FIRST_STAGE: int = 1
"""Reveal step that adds SURVEY; each later stage follows one step behind, so
a step number below :data:`_STEP_RETURN` is also the stage's badge number."""

_STEP_RETURN: int = _STEP_FIRST_STAGE + len(_STAGES)  # one past the last stage
"""Reveal step that closes the loop with the red-link return path."""

_TOTAL_STEPS: int = _STEP_RETURN
"""Number of images in the build-up sequence; the last one is the whole
diagram."""


def _stage_step(index: int) -> int:
    """Return the reveal step at which the stage at ``index`` appears."""
    return _STEP_FIRST_STAGE + index


def _stage_x(index: int) -> float:
    """Return the left edge of the stage box at ``index``."""
    return _MARGIN_L + index * (_STAGE_W + _STAGE_GAP)


def _stage_cx(index: int) -> float:
    """Return the horizontal center of the stage box at ``index``."""
    return _stage_x(index) + _STAGE_W / 2.0


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
    linewidth: float = 1.3,
) -> None:
    """Draw a single arrow from ``start`` to ``end``, bowed by ``rad``."""
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=linewidth,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=0,
            shrinkB=0,
        )
    )


def _stage(ax: Axes, index: int, stage: _Stage) -> None:
    """Draw the stage box at ``index``: badge, title, body, and tooling tag."""
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
        fontsize=18,
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
            _STAGE_TOP - 1.10 - line_index * 0.34,
            line,
            fontsize=12,
            color=style.INK,
            va="center",
        )
    ax.text(
        x + 0.2,
        _STAGE_BOTTOM + 0.30,
        stage.tool,
        fontsize=12,
        color=style.MIST,
        va="center",
    )


def _flow(ax: Axes, step: int) -> None:
    """Draw the left-to-right arrows between consecutive revealed stages."""
    mid_y = _STAGE_BOTTOM + _STAGE_H / 2.0
    for index in range(len(_STAGES) - 1):
        if step < _stage_step(index + 1):
            continue
        _arrow(
            ax,
            (_stage_x(index) + _STAGE_W, mid_y),
            (_stage_x(index + 1), mid_y),
            color=style.SLATE,
        )


def _audit_loop(ax: Axes, step: int) -> None:
    """Draw the audit's return edge — the gate that unlocks publication."""
    if step < _stage_step(3):
        return
    _arrow(
        ax,
        (_stage_cx(3), _STAGE_TOP),
        (_stage_cx(2), _STAGE_TOP),
        color=_AUDIT_COLOR,
        rad=0.28,
        linewidth=1.2,
    )


def _working_files(ax: Axes, step: int) -> None:
    """Draw the scratchpad files the middle stages write and read.

    A file appears with the stage that first writes it.
    """
    files = (
        (1, "sources.md", "every fact, with source"),
        (2, "draft.md", "the page, ready to ship"),
    )
    for index, name, gloss in files:
        if step < _stage_step(index):
            continue
        _box(
            ax,
            (_stage_x(index), _FILE_TOP - _FILE_H, _STAGE_W, _FILE_H),
            facecolor=style.PAPER,
            edgecolor=style.MIST,
            linewidth=1.0,
            linestyle="dashed",
        )
        ax.plot(
            [_stage_cx(index), _stage_cx(index)],
            [_STAGE_BOTTOM, _FILE_TOP],
            color=style.MIST,
            linewidth=0.9,
            linestyle=(0, (2, 2)),
        )
        ax.text(
            _stage_cx(index),
            _FILE_TOP - 0.28,
            name,
            fontsize=12,
            fontweight="bold",
            color=style.INK,
            ha="center",
            va="center",
        )
        ax.text(
            _stage_cx(index),
            _FILE_TOP - 0.62,
            gloss,
            fontsize=12,
            color=style.MIST,
            ha="center",
            va="center",
        )


def _return_loop(ax: Axes, step: int) -> None:
    """Draw the red-link path from the published page back to a new survey."""
    if step < _STEP_RETURN:
        return
    vertices = [
        (_stage_cx(4), _STAGE_BOTTOM),
        (_stage_cx(4), _RETURN_Y),
        (_stage_cx(0), _RETURN_Y),
        (_stage_cx(0), _STAGE_BOTTOM),
    ]
    ax.add_patch(
        FancyArrowPatch(
            path=Path(vertices, [Path.MOVETO, *[Path.LINETO] * 3]),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.2,
            color=_RETURN_COLOR,
            fill=False,
        )
    )


def _render(step: int) -> Figure:
    """Return the diagram revealed up to ``step``, in the full-diagram frame.

    The frame is fixed by the layout constants rather than by what is drawn, so
    every step in the sequence lands in an identically sized image and the
    pieces already on screen never move as the next one arrives.
    """
    # The figure is the drawing plus its border, and the axes spans the figure
    # with one drawing unit per inch, so the file is cropped to the graphic by
    # construction — a tight bbox would not crop it, since a full-figure axes
    # reports its own extent rather than the artists inside it.
    left = _MARGIN_L - _PAD
    right = _MARGIN_R + _PAD
    bottom = _CONTENT_BOTTOM - _PAD
    top = _CONTENT_TOP + _PAD

    fig = plt.figure(figsize=(right - left, top - bottom))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_axis_off()

    for index, stage in enumerate(_STAGES):
        if step >= _stage_step(index):
            _stage(ax, index, stage)
    _flow(ax, step)
    _audit_loop(ax, step)
    _working_files(ax, step)
    _return_loop(ax, step)
    return fig


@figure(name="process-intake")
def process_intake(ctx: FigureContext) -> None:
    """Render the Intake & Refinement diagram, whole and as a build-up sequence.

    Writes ``process-intake`` (the complete graphic) plus ``process-intake-NN``
    for each reveal step, numbered in presentation order.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"process-intake-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "process-intake")
        plt.close(fig)
