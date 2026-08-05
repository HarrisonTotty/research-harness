"""Meta-figure: the harness's Conjecture & Theorem Proving phase, as a slide graphic.

Renders the last research stage described in the repository README — how an
informal conjecture forming out of the results becomes a stress-tested formal
statement, a committed Lean ``Prop``, and finally a proved theorem — as a
single diagram sized for a presentation slide. The stages are the ones the
``add-conjecture`` and ``attack-conjecture`` skills actually run, drawn as one
row rather than two pipelines because they are one chain: the proposition
registered by the first is the literal term the second proves.

Four edges carry the argument, and are the reason this figure is not just a
row of boxes:

* the **refute → refine loop** (indianred) holds STATE until a hunt over
  ranges *beyond* what the experiments swept comes back without a
  counterexample, because the evidence that motivated a conjecture cannot also
  be its stress test;
* the **audit gates** (indianred) hold REGISTER, ATTACK, and PROMOTE until
  their auditor says ``CLEAN`` — drawn in the same hue and the same shape as
  the refute loop, since both say the same thing: nothing leaves this stage
  until it survives a check it did not write;
* the **statement span** (purple) runs from the registered ``Prop`` def to
  PROMOTE, because the def is stated once and referenced literally by the
  proving theorem — the proposition the evidence supported and the one the
  proof establishes are the same term, not two paraphrases of it;
* the **return loop** (teal) runs from what closed back to the graph, which is
  what makes the phase a cycle rather than a terminus: a refuted conjecture
  leaves a counterexample that often matters more than the conjecture did.

Artifacts hanging below the row are drawn solid when they are durable (the
conjecture page, the theorem page) and dashed when they live and die in the
session scratchpad (the working file, the hunt harness, the attack plan): the
alternation is itself part of the process being described. ``hunt/`` hangs
beneath the stage that keeps re-entering it, because it is the one scratchpad
artifact carried across the refute loop's dispatches — each hunter is told what
is already there, so a re-hunt widens the coverage instead of repeating it.

The graphic carries no title, standfirst, or summary chrome — the slide it goes
on supplies those — and is saved cropped to the drawing itself, so it can be
placed and scaled freely on the slide. Its geometry is the sibling of
:mod:`figures.process_intake`, :mod:`figures.process_development`, and
:mod:`figures.process_experimentation`: the same stage grid and the same reveal
mechanics throughout, and the same width and box heights as the other two-beat
figures (intake's row is shorter, one beat per stage), so the four figures sit
together in one deck.

Besides the complete diagram, the module writes a build-up sequence — one image
per reveal step, each adding the next piece of the process — for walking an
audience through it a piece at a time. A step reveals a single element rather
than a whole stage: the box, then the gate or return edge it owns, then the
artifact it leaves behind::

    1  the informal conjecture   10  ATTACK
    2  GROUND                    11  attack audit gate
    3  the working file          12  refutation return
    4  STATE                     13  the attack plan
    5  the falsification hunt    14  PROMOTE
    6  the hunt harness          15  promotion audit gate
    7  REGISTER                  16  statement span
    8  register audit gate       17  the theorem page
    9  the conjecture page       18  return to the graph

Every step is rendered into the same frame as the whole, so the images can be
stacked on one slide (or advanced through) without anything shifting between
them.

Fourth of a series of meta-figures, one per README process phase. Regenerate
with ``just figure process-theorem-proving``.
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
"""Half-width of a gate's loop, measured from the stage's center."""

_REFUTE_Y: float = 6.95
"""Height of the horizontal run of the mid-attack refutation path."""

_REFUTE_X_OFFSET: float = 0.3
"""Inset of the refutation path's verticals from a stage box's left edge; it
keeps them clear of the gate loops, which are centered on the box."""

_FILE_TOP: float = 1.9
"""Top edge of the artifact chips, below the stage row."""

_FILE_H: float = 0.9
"""Height of an artifact chip."""

_STATEMENT_Y: float = 0.65
"""Height of the horizontal run of the statement span, below the chips."""

_STATEMENT_ENTRY_Y: float = _STAGE_BOTTOM + 0.5
"""Height at which the statement span enters PROMOTE, below the flow arrow
that shares the same gap."""

_RETURN_Y: float = 0.25
"""Height of the horizontal run of the return path."""

_SPINE_X: float = _MARGIN_L - 0.32
"""Left-hand column the return path climbs, clear of every box."""

_CONTENT_TOP: float = _INPUT_BOTTOM + _INPUT_H
"""Top of the drawing: the upper edge of the input chip."""

_CONTENT_BOTTOM: float = _RETURN_Y
"""Bottom of the drawing: the horizontal run of the return path."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_AUDIT_COLOR: str = style.CATEGORICAL[5]
"""Indianred — the blog's falsification hue, carrying every gate loop and the
mid-attack refutation; each edge in it says "it did not survive, go back"."""

_RETURN_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the return loop, distinct from the falsification edges."""

_STATEMENT_COLOR: str = style.CATEGORICAL[4]
"""Purple — carries the statement span, the only edge that moves a commitment
forward in time rather than a failure backward."""


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

    ``durable`` distinguishes what survives the session — the conjecture page,
    the theorem page — from what stays in the scratchpad; the two are drawn
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
        title="GROUND",
        first=_Beat(
            label="EVIDENCE",
            lines=(
                "quote note, results, and",
                "the ranges actually swept",
            ),
        ),
        second=_Beat(
            label="DISPOSITION",
            lines=(
                "known-true is intake,",
                "known-false is a bug",
            ),
        ),
        tool="/add-conjecture · the note",
    ),
    _Stage(
        number=2,
        title="STATE",
        first=_Beat(
            label="CHECKPOINTS",
            lines=(
                "every quantifier and side",
                "condition made explicit",
            ),
        ),
        second=_Beat(
            label="FALSIFY",
            lines=(
                "hunt counterexamples",
                "beyond the swept ranges",
            ),
        ),
        tool="counterexample-hunter",
    ),
    _Stage(
        number=3,
        title="REGISTER",
        first=_Beat(
            label="LEAN",
            lines=(
                "commit the statement as",
                "a Prop def in Lean",
            ),
        ),
        second=_Beat(
            label="DOC",
            lines=(
                "the page carries evidence,",
                "coverage, attack plan",
            ),
        ),
        tool="Prop def · statement audit",
    ),
    _Stage(
        number=4,
        title="ATTACK",
        first=_Beat(
            label="PLAN",
            lines=(
                "one lemma DAG — take the",
                "riskiest lemma first",
            ),
        ),
        second=_Beat(
            label="PROVE",
            lines=(
                "state, audit, prove per",
                "lemma, under effort caps",
            ),
        ),
        tool="/attack-conjecture",
    ),
    _Stage(
        number=5,
        title="PROMOTE",
        first=_Beat(
            label="RECORD",
            lines=(
                "every status change is",
                "written through live",
            ),
        ),
        second=_Beat(
            label="CLOSE",
            lines=(
                "audit axioms; prose comes",
                "from the Lean proof",
            ),
        ),
        tool="docs/conj → docs/theorems",
    ),
)
"""The pipeline, left to right."""

_CHIPS: tuple[_Chip, ...] = (
    _Chip(
        index=0,
        step=3,
        name="conjecture.md",
        gloss="evidence map, dispositions",
        durable=False,
    ),
    _Chip(
        index=1,
        step=6,
        name="hunt/",
        gloss="kept across every re-hunt",
        durable=False,
    ),
    _Chip(
        index=2,
        step=9,
        name="docs/conj/<name>.md",
        gloss="beside its Prop def in Lean",
        durable=True,
    ),
    _Chip(
        index=3,
        step=13,
        name="attack.md",
        gloss="the DAG, session-scoped",
        durable=False,
    ),
    _Chip(
        index=4,
        step=17,
        name="docs/theorems/<name>.md",
        gloss="written from the Lean proof",
        durable=True,
    ),
)
"""One artifact per stage, in stage order. Each follows the stage that writes it
by a step, and where the stage owns a gate or a return edge, follows that too."""

_GATES: tuple[tuple[int, int], ...] = ((1, 5), (2, 8), (3, 11), (4, 15))
"""Stages held in a loop until a check they did not write comes back, each
paired with the reveal step that adds the loop — the step after the stage it
holds. GROUND has none: it is the only stage whose output is a reading of what
already exists rather than a claim of its own."""

_STEP_INPUT: int = 1
"""Reveal step that puts the informal conjecture on screen, before anything
acts on it."""

_STAGE_STEPS: tuple[int, ...] = (2, 4, 7, 10, 14)
"""Reveal step at which each stage box appears, in pipeline order. The gaps are
uneven because a stage is followed by whatever it owns — a gate loop, a return
edge, an artifact — one element per step."""

_STEP_REFUTATION: int = 12
"""Reveal step that adds the mid-attack refutation return, after the gate on
ATTACK that sends work back down it."""

_STEP_STATEMENT: int = 16
"""Reveal step that adds the statement span, once PROMOTE exists to receive the
``Prop`` def REGISTER committed."""

_STEP_RETURN: int = 18
"""Reveal step that closes the loop with the return path."""

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


def _conjecture_input(ax: Axes) -> None:
    """Draw the informal-conjecture feed above GROUND."""
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
        "the informal conjecture",
        fontsize=12,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        _stage_cx(0),
        _INPUT_BOTTOM + 0.28,
        "left by /critique-results",
        fontsize=12,
        color=style.MIST,
        ha="center",
        va="center",
    )


def _conjecture_arrow(ax: Axes) -> None:
    """Draw the arrow from the informal conjecture down into GROUND.

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


def _gates(ax: Axes, step: int) -> None:
    """Draw the loop each stage is held in until its check comes back.

    The hunt loop on STATE and the three auditor loops share one shape and one
    hue because they are one move: the work does not leave the stage until
    something that did not write it — a counterexample search, a fresh-eyes
    auditor — fails to break it.
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


def _statement_span(ax: Axes, step: int) -> None:
    """Draw the path from the registered Prop def to the stage that proves it.

    It runs under the chip row and enters PROMOTE from the gap on its left, so
    it crosses nothing: the statement travels the length of the attack
    untouched, which is the point of committing it as a def instead of
    restating it at proof time.
    """
    if step < _STEP_STATEMENT:
        return
    gap_x = _stage_x(4) - _STAGE_GAP / 2.0
    _routed_arrow(
        ax,
        [
            (_stage_cx(2), _FILE_TOP - _FILE_H),
            (_stage_cx(2), _STATEMENT_Y),
            (gap_x, _STATEMENT_Y),
            (gap_x, _STATEMENT_ENTRY_Y),
            (_stage_x(4), _STATEMENT_ENTRY_Y),
        ],
        color=_STATEMENT_COLOR,
    )


def _refutation_return(ax: Axes, step: int) -> None:
    """Draw the return from a mid-attack refutation to the statement."""
    if step < _STEP_REFUTATION:
        return
    start_x = _stage_x(3) + _REFUTE_X_OFFSET
    end_x = _stage_x(1) + _REFUTE_X_OFFSET
    _routed_arrow(
        ax,
        [
            (start_x, _STAGE_TOP),
            (start_x, _REFUTE_Y),
            (end_x, _REFUTE_Y),
            (end_x, _STAGE_TOP),
        ],
        color=_AUDIT_COLOR,
    )


def _return_loop(ax: Axes, step: int) -> None:
    """Draw the path from what closed back to the next informal conjecture.

    It climbs the left spine rather than the stage centers, so it clears the
    artifact chips, and it lands on the informal conjecture itself: a proved
    theorem and a refuted statement both go back to the graph, and the
    counterexample is often worth more than the conjecture was.
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

    # The conjecture opens the sequence; the arrow into GROUND waits for GROUND.
    if step >= _STEP_INPUT:
        _conjecture_input(ax)
    if step >= _stage_step(0):
        _conjecture_arrow(ax)
    for index, stage in enumerate(_STAGES):
        if step >= _stage_step(index):
            _stage(ax, index, stage)
    _flow(ax, step)
    _gates(ax, step)
    _artifacts(ax, step)
    _statement_span(ax, step)
    _refutation_return(ax, step)
    _return_loop(ax, step)
    return fig


@figure(name="process-theorem-proving")
def process_theorem_proving(ctx: FigureContext) -> None:
    """Render the Conjecture & Theorem Proving diagram, whole and as a sequence.

    Writes ``process-theorem-proving`` (the complete graphic) plus
    ``process-theorem-proving-NN`` for each reveal step, numbered in
    presentation order.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"process-theorem-proving-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "process-theorem-proving")
        plt.close(fig)
