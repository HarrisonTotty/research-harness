"""Meta-figure: the harness's Experimentation & Exploration phases, as a slide graphic.

Renders the middle of the research process described in the repository README —
how a hypothesis becomes a pre-registered design, a runnable sweep, and finally
an audited results document — as a single diagram sized for a presentation
slide. The two README sections are drawn as one figure because they are one
chain: the ``design-experiment`` / ``add-experiment`` stages exist to make the
``explore-results`` stages honest, and the sweep the user invokes by hand is
the hinge between them.

Three edges carry the argument, and are the reason this figure is not just a
row of boxes:

* the **pre-registration span** (purple) runs from the frozen design doc to
  ANALYZE, because predictions written before the code exists are what make an
  unexpected result well-defined rather than a story told afterwards;
* the **bug-first return** (indianred) runs from INTERPRET all the way back to
  IMPLEMENT, because a contradiction with the literature is a bug in the
  experiment until proven otherwise;
* the **write-back loop** (teal) runs from the published results doc to the
  next hypothesis, which is what makes the phase a cycle rather than a line.

The two audit gates are drawn as loops on the stages that own them —
``experiment-auditor`` on IMPLEMENT, ``results-auditor`` on INTERPRET — rather
than as stages of their own, because here they sit inside a skill instead of
between two of them. Artifacts hanging below the row are drawn solid when they
are durable (the design doc, the run, the results doc) and dashed when they
live and die in the session scratchpad (the pilot, the analysis scripts): the
alternation is itself part of the process being described.

The graphic carries no title, standfirst, or summary chrome — the slide it goes
on supplies those — and is saved cropped to the drawing itself, so it can be
placed and scaled freely on the slide. Its geometry is the sibling of
:mod:`figures.process_intake` and :mod:`figures.process_development`: same
width, same stage width, same box heights, same reveal mechanics, so the three
figures sit together in one deck.

Besides the complete diagram, the module writes a build-up sequence — one image
per reveal step, each adding the next piece of the process — for walking an
audience through it a piece at a time. A step reveals a single element rather
than a whole stage: the box, then the gate or return edge it owns, then the
artifact it leaves behind::

    1  the hypothesis          9  ANALYZE
    2  DESIGN                 10  pre-registration span
    3  the design doc         11  the analysis scratchpad
    4  IMPLEMENT              12  INTERPRET
    5  experiment-auditor     13  results-auditor
    6  the pilot              14  bug-first return
    7  RUN                    15  the results doc
    8  the run                16  write-back

Every step is rendered into the same frame as the whole, so the images can be
stacked on one slide (or advanced through) without anything shifting between
them.

Third of a series of meta-figures, one per README process phase. Regenerate
with ``just figure process-experimentation``.
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

_STAGE_H: float = 3.8
"""Height of a pipeline stage box; it carries two beats of its stage."""

_STAGE_GAP: float = 0.45
"""Horizontal gap between consecutive stage boxes, holding the flow arrow."""

_STAGE_TOP: float = 6.0
"""Top edge of the stage row."""

_STAGE_BOTTOM: float = _STAGE_TOP - _STAGE_H
"""Bottom edge of the stage row."""

_INPUT_BOTTOM: float = 6.55
"""Bottom edge of the input chip, above the stage row."""

_INPUT_H: float = 0.9
"""Height of the input chip."""

_FIRST_LABEL_Y: float = _STAGE_TOP - 1.10
"""Baseline of the first beat's label inside a stage box."""

_SECOND_LABEL_Y: float = _STAGE_TOP - 2.40
"""Baseline of the second beat's label inside a stage box."""

_BEAT_DIVIDER_Y: float = _STAGE_TOP - 2.10
"""Height of the hairline separating a stage's two beats."""

_LINE_STEP: float = 0.34
"""Vertical distance between consecutive body lines."""

_GATE_SPAN: float = 0.45
"""Half-width of an audit gate's loop, measured from the stage's center."""

_BUG_Y: float = 6.95
"""Height of the horizontal run of the bug-first return path."""

_BUG_X_OFFSET: float = 0.3
"""Inset of the bug-first path's verticals from a stage box's left edge; it
keeps them clear of the audit loops, which are centered on the box."""

_FILE_TOP: float = 1.9
"""Top edge of the artifact chips, below the stage row."""

_FILE_H: float = 0.9
"""Height of an artifact chip."""

_PREREG_Y: float = 0.65
"""Height of the horizontal run of the pre-registration span, below the chips."""

_PREREG_ENTRY_Y: float = _STAGE_BOTTOM + 0.5
"""Height at which the pre-registration span enters ANALYZE, below the flow
arrow that shares the same gap."""

_RETURN_Y: float = 0.25
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
"""Indianred — the blog's falsification hue, carrying both audit gates and the
bug-first return; every edge in it says "something failed, go back"."""

_RETURN_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the write-back loop, distinct from the audit edges."""

_PREREG_COLOR: str = style.CATEGORICAL[4]
"""Purple — carries the pre-registration span, the only edge that moves a
commitment forward in time rather than a failure backward."""


@dataclass(frozen=True, slots=True)
class _Beat:
    """One half of a stage box: a labelled move within that stage.

    Body lines are pre-wrapped rather than flowed: at this size the wrap
    points are a layout decision, not something to leave to a text engine.
    """

    label: str
    lines: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _Stage:
    """One box in the pipeline, split into the two beats it runs in order."""

    number: int
    title: str
    first: _Beat
    second: _Beat
    tool: str


@dataclass(frozen=True, slots=True)
class _Chip:
    """An artifact hanging below the stage that writes it.

    ``durable`` distinguishes what survives the session — the design doc, the
    run, the results doc — from what stays in the scratchpad; the two are drawn
    with solid and dashed edges respectively.
    """

    index: int
    step: int
    name: str
    gloss: str
    durable: bool


_STAGES: tuple[_Stage, ...] = (
    _Stage(
        number=1,
        title="DESIGN",
        first=_Beat(
            label="CHECKPOINTS",
            lines=(
                "hypothesis · controls ·",
                "seeds · schema, one by one",
            ),
        ),
        second=_Beat(
            label="PRE-REGISTER",
            lines=(
                "predict per region, cited;",
                "mark the unpredicted",
            ),
        ),
        tool="/design-experiment",
    ),
    _Stage(
        number=2,
        title="IMPLEMENT",
        first=_Beat(
            label="BUILD",
            lines=(
                "one module whose defaults",
                "are the designed sweep",
            ),
        ),
        second=_Beat(
            label="PILOT",
            lines=(
                "a few cells only — check",
                "schema, runtime, variance",
            ),
        ),
        tool="/add-experiment · pilot",
    ),
    _Stage(
        number=3,
        title="RUN",
        first=_Beat(
            label="SWEEP",
            lines=(
                "you launch the full grid —",
                "a clean audit unlocks it",
            ),
        ),
        second=_Beat(
            label="ARTIFACTS",
            lines=(
                "a timestamped frame beside",
                "meta.json: seeds, commit",
            ),
        ),
        tool="just experiment <name>",
    ),
    _Stage(
        number=4,
        title="ANALYZE",
        first=_Beat(
            label="PROFILE",
            lines=(
                "profile the frame against",
                "the expected grid first",
            ),
        ),
        second=_Beat(
            label="SCRIPT",
            lines=(
                "grade every prediction,",
                "never from raw frames",
            ),
        ),
        tool="/explore-results · profiler",
    ),
    _Stage(
        number=5,
        title="INTERPRET",
        first=_Beat(
            label="TRIAGE",
            lines=(
                "a contradiction with the",
                "literature is a bug first",
            ),
        ),
        second=_Beat(
            label="WRITE",
            lines=(
                "draft from the numbers,",
                "audit, publish, write back",
            ),
        ),
        tool="literature-checker · Logseq",
    ),
)
"""The pipeline, left to right."""

_CHIPS: tuple[_Chip, ...] = (
    _Chip(
        index=0,
        step=3,
        name="docs/experiments/<name>.md",
        gloss="frozen at the first run",
        durable=True,
    ),
    _Chip(
        index=1,
        step=6,
        name="pilot/",
        gloss="scratchpad only, never a run",
        durable=False,
    ),
    _Chip(
        index=2,
        step=8,
        name="data/results/<name>.<ts>",
        gloss="the frame and its metadata",
        durable=True,
    ),
    _Chip(
        index=3,
        step=11,
        name="analysis/ · findings.md",
        gloss="scripts, outputs, verdicts",
        durable=False,
    ),
    _Chip(
        index=4,
        step=15,
        name="docs/results/<name>.<ts>.md",
        gloss="every number recomputable",
        durable=True,
    ),
)
"""One artifact per stage, in stage order. Each follows the stage that writes it
by a step, and where the stage owns a gate or a return edge, follows that too."""

_GATES: tuple[tuple[int, int], ...] = ((1, 5), (4, 13))
"""Stages that loop on their own auditor, each paired with the reveal step that
adds the loop — the step after the stage it holds."""

_STEP_INPUT: int = 1
"""Reveal step that puts the hypothesis on screen, before anything acts on it."""

_STAGE_STEPS: tuple[int, ...] = (2, 4, 7, 9, 12)
"""Reveal step at which each stage box appears, in pipeline order. The gaps are
uneven because a stage is followed by whatever it owns — an audit gate, a return
edge, an artifact — one element per step."""

_STEP_PREREG: int = 10
"""Reveal step that adds the pre-registration span, once ANALYZE exists to
receive the predictions the design doc froze."""

_STEP_BUG_RETURN: int = 14
"""Reveal step that adds the bug-first return, after the triage gate that sends
work back down it."""

_STEP_RETURN: int = 16
"""Reveal step that closes the loop with the write-back path."""

_TOTAL_STEPS: int = _STEP_RETURN
"""Number of images in the build-up sequence; the last one is the whole
diagram."""


def _stage_step(index: int) -> int:
    """Return the reveal step at which the stage at ``index`` appears."""
    return _STAGE_STEPS[index]


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


def _routed_arrow(
    ax: Axes,
    vertices: list[tuple[float, float]],
    *,
    color: str,
    linewidth: float = 1.2,
) -> None:
    """Draw an arrow along the straight run through ``vertices``.

    The head lands on the final vertex, so a route is written in the order it
    is travelled.
    """
    ax.add_patch(
        FancyArrowPatch(
            path=Path(vertices, [Path.MOVETO, *[Path.LINETO] * (len(vertices) - 1)]),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=linewidth,
            color=color,
            fill=False,
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


def _hypothesis_input(ax: Axes) -> None:
    """Draw the hypothesis feed above DESIGN."""
    _box(
        ax,
        (_stage_x(0), _INPUT_BOTTOM, _STAGE_W, _INPUT_H),
        facecolor=style.PAPER,
        edgecolor=style.MIST,
        linewidth=1.0,
    )
    ax.text(
        _stage_cx(0),
        _INPUT_BOTTOM + 0.62,
        "the hypothesis",
        fontsize=12,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        _stage_cx(0),
        _INPUT_BOTTOM + 0.28,
        "question, claim, follow-up",
        fontsize=12,
        color=style.MIST,
        ha="center",
        va="center",
    )


def _hypothesis_arrow(ax: Axes) -> None:
    """Draw the arrow from the hypothesis down into DESIGN.

    It arrives with the box it lands on rather than with the chip it leaves, so
    no step ever shows an arrow into empty space.
    """
    _arrow(
        ax,
        (_stage_cx(0), _INPUT_BOTTOM),
        (_stage_cx(0), _STAGE_TOP),
        color=style.SLATE,
    )


def _beat(ax: Axes, index: int, beat: _Beat, top: float) -> None:
    """Draw one beat inside a stage box, downward from ``top``."""
    x = _stage_x(index)
    ax.text(
        x + 0.2,
        top,
        beat.label,
        fontsize=12,
        fontweight="bold",
        color=style.SLATE,
        va="center",
    )
    for line_index, line in enumerate(beat.lines):
        ax.text(
            x + 0.2,
            top - 0.38 - line_index * _LINE_STEP,
            line,
            fontsize=12,
            color=style.INK,
            va="center",
        )


def _stage(ax: Axes, index: int, stage: _Stage) -> None:
    """Draw the stage box at ``index``: badge, title, both beats, tooling tag."""
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
    _rule(ax, index, _STAGE_TOP - 0.78)

    _beat(ax, index, stage.first, _FIRST_LABEL_Y)
    _rule(ax, index, _BEAT_DIVIDER_Y)
    _beat(ax, index, stage.second, _SECOND_LABEL_Y)

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


def _audit_gates(ax: Axes, step: int) -> None:
    """Draw the loop each auditor holds its own stage in until it says CLEAN.

    The gates sit on their stages rather than between stages because in this
    phase both auditors are dispatched from inside a skill: nothing downstream
    of them starts, and nothing upstream is revisited, until the audit passes.
    """
    for index, gate_step in _GATES:
        if step < gate_step:
            continue
        center = _stage_cx(index)
        _arrow(
            ax,
            (center + _GATE_SPAN, _STAGE_TOP),
            (center - _GATE_SPAN, _STAGE_TOP),
            color=_AUDIT_COLOR,
            rad=0.85,
            linewidth=1.2,
        )


def _artifacts(ax: Axes, step: int) -> None:
    """Draw the artifact chips, each a step behind the stage that writes it."""
    for chip in _CHIPS:
        if step < chip.step:
            continue
        _box(
            ax,
            (_stage_x(chip.index), _FILE_TOP - _FILE_H, _STAGE_W, _FILE_H),
            facecolor=style.PAPER,
            edgecolor=style.MIST,
            linewidth=1.0,
            linestyle="solid" if chip.durable else "dashed",
        )
        ax.plot(
            [_stage_cx(chip.index), _stage_cx(chip.index)],
            [_STAGE_BOTTOM, _FILE_TOP],
            color=style.MIST,
            linewidth=0.9,
            linestyle=(0, (2, 2)),
        )
        ax.text(
            _stage_cx(chip.index),
            _FILE_TOP - 0.28,
            chip.name,
            fontsize=12,
            fontweight="bold",
            color=style.INK,
            ha="center",
            va="center",
        )
        ax.text(
            _stage_cx(chip.index),
            _FILE_TOP - 0.62,
            chip.gloss,
            fontsize=12,
            color=style.MIST,
            ha="center",
            va="center",
        )


def _prereg_span(ax: Axes, step: int) -> None:
    """Draw the path from the frozen design doc to the stage that grades it.

    It runs under the chip row and enters ANALYZE from the gap on its left, so
    it crosses nothing: the predictions travel the length of the phase
    untouched, which is the point of writing them down first.
    """
    if step < _STEP_PREREG:
        return
    gap_x = _stage_x(3) - _STAGE_GAP / 2.0
    _routed_arrow(
        ax,
        [
            (_stage_cx(0), _FILE_TOP - _FILE_H),
            (_stage_cx(0), _PREREG_Y),
            (gap_x, _PREREG_Y),
            (gap_x, _PREREG_ENTRY_Y),
            (_stage_x(3), _PREREG_ENTRY_Y),
        ],
        color=_PREREG_COLOR,
    )


def _bug_first_return(ax: Axes, step: int) -> None:
    """Draw the return from a literature contradiction to the experiment."""
    if step < _STEP_BUG_RETURN:
        return
    start_x = _stage_x(4) + _BUG_X_OFFSET
    end_x = _stage_x(1) + _BUG_X_OFFSET
    _routed_arrow(
        ax,
        [
            (start_x, _STAGE_TOP),
            (start_x, _BUG_Y),
            (end_x, _BUG_Y),
            (end_x, _STAGE_TOP),
        ],
        color=_AUDIT_COLOR,
    )


def _return_loop(ax: Axes, step: int) -> None:
    """Draw the path from the published results doc back to a new hypothesis.

    It climbs the left spine rather than the stage centers, so it clears the
    artifact chips, and it lands on the hypothesis itself: the follow-ups a
    results doc proposes are what the next design starts from.
    """
    if step < _STEP_RETURN:
        return
    _routed_arrow(
        ax,
        [
            (_stage_cx(4), _FILE_TOP - _FILE_H),
            (_stage_cx(4), _RETURN_Y),
            (_SPINE_X, _RETURN_Y),
            (_SPINE_X, _INPUT_BOTTOM + _INPUT_H / 2.0),
            (_stage_x(0), _INPUT_BOTTOM + _INPUT_H / 2.0),
        ],
        color=_RETURN_COLOR,
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

    # The hypothesis opens the sequence; the arrow into DESIGN waits for DESIGN.
    if step >= _STEP_INPUT:
        _hypothesis_input(ax)
    if step >= _stage_step(0):
        _hypothesis_arrow(ax)
    for index, stage in enumerate(_STAGES):
        if step >= _stage_step(index):
            _stage(ax, index, stage)
    _flow(ax, step)
    _audit_gates(ax, step)
    _artifacts(ax, step)
    _prereg_span(ax, step)
    _bug_first_return(ax, step)
    _return_loop(ax, step)
    return fig


@figure(name="process-experimentation")
def process_experimentation(ctx: FigureContext) -> None:
    """Render the Experimentation & Exploration diagram, whole and as a sequence.

    Writes ``process-experimentation`` (the complete graphic) plus
    ``process-experimentation-NN`` for each reveal step, numbered in
    presentation order.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"process-experimentation-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "process-experimentation")
        plt.close(fig)
