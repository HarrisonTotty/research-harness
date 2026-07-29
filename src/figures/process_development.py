"""Meta-figure: the harness's Development phase, as a slide graphic.

Renders the second phase of the research process described in the repository
README — how a page in the concept graph becomes a tested Python structure and
a proved Lean module — as a single diagram sized for a presentation slide. The
stages are the ones the ``add-python-topic`` and ``add-lean-topic`` skills
actually run, and the two skills are drawn as halves of one box per stage
rather than as two pipelines, because they are step-for-step the same process:
survey what exists, write the plan file, build, audit against the plan, gate.
The two feedback edges (the audit gate on proving and on the test run, and the
implementation notes written back to the page) are drawn because they are what
make the phase a loop rather than a line.

The graphic carries no title, standfirst, or summary chrome — the slide it goes
on supplies those — and is saved cropped to the drawing itself, so it can be
placed and scaled freely on the slide. Its geometry is the sibling of
:mod:`figures.process_intake`: same width, same stage width, same reveal
mechanics, so the two figures sit together in one deck.

Besides the complete diagram, the module writes a build-up sequence — one image
per reveal step, each adding the next piece of the process — for walking an
audience through it a stage at a time. A step reveals a stage together with
whatever feeds it, so step *n* is stage *n* until the last step closes the
loop. Every step is rendered into the same frame as the whole, so the images
can be stacked on one slide (or advanced through) without anything shifting
between them.

Second of a series of meta-figures, one per README process phase. Regenerate
with ``just figure process-development``.
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

_MARGIN_R: float = 15.4
"""Right edge of every full-width element."""

_STAGE_W: float = 2.6
"""Width of a pipeline stage box."""

_STAGE_H: float = 3.8
"""Height of a pipeline stage box; it carries two language tracks."""

_STAGE_GAP: float = 0.45
"""Horizontal gap between consecutive stage boxes, holding the flow arrow."""

_STAGE_TOP: float = 6.0
"""Top edge of the stage row."""

_STAGE_BOTTOM: float = _STAGE_TOP - _STAGE_H
"""Bottom edge of the stage row."""

_INPUT_BOTTOM: float = 6.55
"""Bottom edge of the input chip, above the stage row."""

_INPUT_H: float = 0.8
"""Height of the input chip."""

_PYTHON_LABEL_Y: float = _STAGE_TOP - 0.98
"""Baseline of the Python track's label inside a stage box."""

_LEAN_LABEL_Y: float = _STAGE_TOP - 2.44
"""Baseline of the Lean track's label inside a stage box."""

_TRACK_DIVIDER_Y: float = _STAGE_TOP - 2.12
"""Height of the hairline separating the two language tracks."""

_LINE_STEP: float = 0.28
"""Vertical distance between consecutive body lines."""

_FILE_TOP: float = 1.9
"""Top edge of the working-file chips, below the stage row."""

_FILE_H: float = 0.7
"""Height of a working-file chip."""

_RETURN_Y: float = 0.75
"""Height of the horizontal run of the write-back path."""

_SPINE_X: float = _MARGIN_L - 0.32
"""Left-hand column the write-back path climbs, clear of every box."""

_CONTENT_TOP: float = _INPUT_BOTTOM + _INPUT_H
"""Top of the drawing: the upper edge of the input chip."""

_CONTENT_BOTTOM: float = _RETURN_Y
"""Bottom of the drawing: the horizontal run of the write-back path."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_AUDIT_COLOR: str = style.CATEGORICAL[5]
"""Indianred — the blog's falsification hue, carrying the audit gate."""

_RETURN_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the write-back loop, distinct from the audit loop."""


@dataclass(frozen=True, slots=True)
class _Stage:
    """One box in the development pipeline, split across the two stacks.

    Body lines are pre-wrapped rather than flowed: at this size the wrap
    points are a layout decision, not something to leave to a text engine.
    """

    number: int
    title: str
    python: tuple[str, ...]
    lean: tuple[str, ...]
    tool: str


_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="SURVEY",
        python=(
            "read what src/research already",
            "has and extend its shape —",
            "never open a second style",
        ),
        lean=(
            "map every claim onto Mathlib",
            "first — parallel scout agents;",
            "never redefine what it has",
        ),
        tool="page.md · mathlib-scout",
    ),
    _Stage(
        number=2,
        title="MAP",
        python=(
            "spec.md maps page section to",
            "code artifact: representation,",
            "constructors, method by method",
        ),
        lean=(
            "coverage.md rules on each",
            "claim: reuse · extend ·",
            "define · backlog",
        ),
        tool="every entry cites its block",
    ),
    _Stage(
        number=3,
        title="BUILD",
        python=(
            "implement the module, then the",
            "suite: examples become",
            "fixtures, theorems properties",
        ),
        lean=(
            "state every declaration with",
            "docstring and attribution;",
            "proof bodies stay sorry",
        ),
        tool="docstrings cite the page",
    ),
    _Stage(
        number=4,
        title="AUDIT",
        python=(
            "test-auditor diffs the suite",
            "against page.md — one reading",
            "wrote both; this is a second",
        ),
        lean=(
            "statement-auditor diffs every",
            "statement against the map —",
            "prove nothing unaudited",
        ),
        tool="re-dispatch until CLEAN",
    ),
    _Stage(
        number=5,
        title="GATE",
        python=(
            "just check — a red test is a",
            "bug, a bad transcription, or",
            "a false claim on the page",
        ),
        lean=(
            "discharge every sorry, then",
            "just lean-check and audit the",
            "axioms: no sorryAx escapes",
        ),
        tool="the build audits the graph",
    ),
)
"""The pipeline, left to right."""

_STEP_FIRST_STAGE: int = 1
"""Reveal step that adds SURVEY; each later stage follows one step behind, so
a step number below :data:`_STEP_RETURN` is also the stage's badge number."""

_STEP_RETURN: int = _STEP_FIRST_STAGE + len(_STAGES)  # one past the last stage
"""Reveal step that closes the loop with the write-back path."""

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


def _rule(ax: Axes, index: int, y: float) -> None:
    """Draw a hairline across the stage box at ``index``, at height ``y``."""
    x = _stage_x(index)
    ax.plot(
        [x + 0.2, x + _STAGE_W - 0.2],
        [y, y],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )


def _page_input(ax: Axes) -> None:
    """Draw the page feed above SURVEY, and its arrow into the stage."""
    _box(
        ax,
        (_stage_x(0), _INPUT_BOTTOM, _STAGE_W, _INPUT_H),
        facecolor=style.PAPER,
        edgecolor=style.MIST,
        linewidth=1.0,
    )
    ax.text(
        _stage_cx(0),
        7.13,
        "the Logseq page",
        fontsize=9.5,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        _stage_cx(0),
        6.83,
        "definition, theorems, examples",
        fontsize=8.5,
        color=style.MIST,
        ha="center",
        va="center",
    )
    _arrow(
        ax,
        (_stage_cx(0), _INPUT_BOTTOM),
        (_stage_cx(0), _STAGE_TOP),
        color=style.SLATE,
    )


def _track(
    ax: Axes, index: int, label: str, lines: tuple[str, ...], top: float
) -> None:
    """Draw one language track inside a stage box, downward from ``top``."""
    x = _stage_x(index)
    ax.text(
        x + 0.2,
        top,
        label,
        fontsize=8.5,
        fontweight="bold",
        color=style.SLATE,
        va="center",
    )
    for line_index, line in enumerate(lines):
        ax.text(
            x + 0.2,
            top - 0.30 - line_index * _LINE_STEP,
            line,
            fontsize=8.5,
            color=style.INK,
            va="center",
        )


def _stage(ax: Axes, index: int, stage: _Stage) -> None:
    """Draw the stage box at ``index``: badge, title, both tracks, tooling tag."""
    x = _stage_x(index)
    _box(
        ax,
        (x, _STAGE_BOTTOM, _STAGE_W, _STAGE_H),
        facecolor=style.PAPER,
        edgecolor=style.SLATE,
        linewidth=1.4,
    )

    badge_y = _STAGE_TOP - 0.36
    ax.add_patch(Circle((x + 0.34, badge_y), 0.17, facecolor=style.SLATE, lw=0))
    ax.text(
        x + 0.34,
        badge_y,
        str(stage.number),
        fontsize=9.5,
        fontweight="bold",
        color=style.PAPER,
        ha="center",
        va="center",
    )
    ax.text(
        x + 0.64,
        badge_y,
        stage.title,
        fontsize=12,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    _rule(ax, index, _STAGE_TOP - 0.66)

    _track(ax, index, "PYTHON", stage.python, _PYTHON_LABEL_Y)
    _rule(ax, index, _TRACK_DIVIDER_Y)
    _track(ax, index, "LEAN", stage.lean, _LEAN_LABEL_Y)

    ax.text(
        x + 0.2,
        _STAGE_BOTTOM + 0.22,
        stage.tool,
        fontsize=8,
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
    """Draw the audit's return edge — the gate that unlocks the run."""
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
    """Draw the scratchpad files the early stages write and later stages read.

    A file appears with the stage that first writes it.
    """
    files = (
        (0, "page.md", "the page, verbatim"),
        (1, "spec.md · coverage.md", "the plan, one row per claim"),
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
            _FILE_TOP - 0.24,
            name,
            fontsize=9,
            fontweight="bold",
            color=style.INK,
            ha="center",
            va="center",
        )
        ax.text(
            _stage_cx(index),
            _FILE_TOP - 0.48,
            gloss,
            fontsize=8,
            color=style.MIST,
            ha="center",
            va="center",
        )


def _return_loop(ax: Axes, step: int) -> None:
    """Draw the path from a passing gate back to the page it was built from.

    It climbs the left spine rather than the stage centers, so it clears the
    working-file chips, and it lands on the page itself: what the gate proves
    is what the implementation notes on the page then claim.
    """
    if step < _STEP_RETURN:
        return
    vertices = [
        (_stage_cx(4), _STAGE_BOTTOM),
        (_stage_cx(4), _RETURN_Y),
        (_SPINE_X, _RETURN_Y),
        (_SPINE_X, _INPUT_BOTTOM + _INPUT_H / 2.0),
        (_stage_x(0), _INPUT_BOTTOM + _INPUT_H / 2.0),
    ]
    ax.add_patch(
        FancyArrowPatch(
            path=Path(vertices, [Path.MOVETO, *[Path.LINETO] * 4]),
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
    left = _SPINE_X - _PAD
    right = _MARGIN_R + _PAD
    bottom = _CONTENT_BOTTOM - _PAD
    top = _CONTENT_TOP + _PAD

    fig = plt.figure(figsize=(right - left, top - bottom))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # The page arrives with the stage that reads it.
    if step >= _stage_step(0):
        _page_input(ax)
    for index, stage in enumerate(_STAGES):
        if step >= _stage_step(index):
            _stage(ax, index, stage)
    _flow(ax, step)
    _audit_loop(ax, step)
    _working_files(ax, step)
    _return_loop(ax, step)
    return fig


@figure(name="process-development")
def process_development(ctx: FigureContext) -> None:
    """Render the Development diagram, whole and as a build-up sequence.

    Writes ``process-development`` (the complete graphic) plus
    ``process-development-NN`` for each reveal step, numbered in presentation
    order.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"process-development-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "process-development")
        plt.close(fig)
