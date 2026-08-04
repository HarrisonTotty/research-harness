"""Meta-figure: the whole harness process in five phases, as a slide graphic.

Renders the research process described in the repository README end to end —
how a paper becomes a concept page, a tested implementation, a graded
experiment, a proved theorem, and finally something published — as a single
diagram sized for a presentation slide. It is the overview the four
phase-level meta-figures hang beneath: where those draw the stages *inside* a
phase, this one draws a phase per box and keeps only what a phase is, what you
type to run it, and the audit that has to come back ``CLEAN`` before it ends.

The README's own sectioning is collapsed slightly to get five: Experimentation,
Exploration, and Feedback are one box (EXPLORATION), because designing a sweep,
grading it, and arguing about the reading are one loop around a single run,
and Figure Generation joins the docs build and the post as PUBLISHING.

Four edges carry the argument, and are the reason this is not just a row of
boxes:

* the **bug-first return** (indianred) runs from EXPLORATION back to
  DEVELOPMENT, because a contradiction with the recorded literature is a bug
  in the code until proven otherwise;
* the **refutation return** (indianred) runs from THEOREM PROVING back to
  EXPLORATION, because a counterexample found while stating a conjecture
  re-opens the results that suggested it;
* the **write-back rail** (purple) runs from every phase that builds something
  back into the graph, because each skill ends by recording what was actually
  built or found on the relevant Logseq page — the graph is the spine, not the
  first phase's output;
* the **publication loop** (teal) runs from what ships back to the next intake,
  which is what makes the process a cycle rather than a line.

The graphic carries no title, standfirst, or summary chrome — the slide it goes
on supplies those — and is saved cropped to the drawing itself, so it can be
placed and scaled freely on the slide. Its geometry is the sibling of
:mod:`figures.process_intake`, :mod:`figures.process_development`, and
:mod:`figures.process_experimentation`: same width, same stage width, same gap,
same reveal mechanics, so all of them sit together in one deck.

Besides the complete diagram, the module writes a build-up sequence — one image
per reveal step, each adding the next piece of the process — for walking an
audience through it a phase at a time. A step reveals a phase together with
whatever feeds it, so step *n* is phase *n* until the last step closes the
loop. Every step is rendered into the same frame as the whole, so the images
can be stacked on one slide (or advanced through) without anything shifting
between them.

Regenerate with ``just figure process-overview``.
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
"""Width of a phase box."""

_STAGE_H: float = 4.09
"""Height of a phase box; it carries prose, commands, and the audit gate."""

_STAGE_GAP: float = 0.45
"""Horizontal gap between consecutive phase boxes, holding the flow arrow."""

_STAGE_TOP: float = 8.0
"""Top edge of the phase row."""

_STAGE_BOTTOM: float = _STAGE_TOP - _STAGE_H
"""Bottom edge of the phase row."""

_HEAD_RULE_Y: float = _STAGE_TOP - 0.66
"""Height of the hairline under a phase's title."""

_PROSE_TOP_Y: float = _STAGE_TOP - 0.98
"""Baseline of the first prose line inside a phase box."""

_LINE_STEP: float = 0.28
"""Vertical distance between consecutive body lines."""

_GATE_RULE_Y: float = _STAGE_BOTTOM + 0.81
"""Height of the hairline separating the audit gate from the commands."""

_COMMAND_CENTER_Y: float = _GATE_RULE_Y + 0.68
"""Middle of the band holding a phase's commands. The list is centered on it
rather than anchored to either end: the phases run from one command to four,
and a centered block leaves that difference as even margin instead of a hole
under the prose."""

_GATE_TOP_Y: float = _GATE_RULE_Y - 0.27
"""Baseline of the first audit-gate line."""

_INPUT_BOTTOM: float = 8.55
"""Bottom edge of the input and output chips, above the phase row."""

_INPUT_H: float = 0.8
"""Height of the input and output chips."""

_BUG_Y: float = 8.62
"""Height of the horizontal run of the bug-first return path."""

_BUG_X_OFFSET: float = 0.3
"""Inset of the bug-first path's verticals from a phase box's left edge."""

_REFUTE_Y: float = 9.16
"""Height of the horizontal run of the refutation return path, above the
bug-first path it nests around."""

_REFUTE_X_OFFSET: float = 0.65
"""Inset of the refutation path's verticals; wider than
:data:`_BUG_X_OFFSET` so the two paths never share a vertical where they meet
on the same box."""

_FILE_TOP: float = _STAGE_BOTTOM - 0.3
"""Top edge of the artifact chips, below the phase row."""

_FILE_H: float = 0.7
"""Height of an artifact chip."""

_WRITEBACK_Y: float = _FILE_TOP - _FILE_H - 0.7
"""Height of the horizontal run of the write-back rail, below the chips."""

_RETURN_Y: float = _WRITEBACK_Y - 0.7
"""Height of the horizontal run of the publication loop."""

_SPINE_X: float = _MARGIN_L - 0.32
"""Left-hand column the publication loop climbs, clear of every box."""

_CONTENT_TOP: float = _INPUT_BOTTOM + _INPUT_H
"""Top of the drawing: the upper edge of the input and output chips."""

_CONTENT_BOTTOM: float = _RETURN_Y
"""Bottom of the drawing: the horizontal run of the publication loop."""

_CORNER: float = 0.12
"""Corner rounding applied to every box, in inches."""

_AUDIT_COLOR: str = style.CATEGORICAL[5]
"""Indianred — the blog's falsification hue, carrying both return paths and
the audit gate printed in every phase box."""

_RETURN_COLOR: str = style.CATEGORICAL[1]
"""Teal — carries the publication loop, distinct from the audit edges."""

_WRITEBACK_COLOR: str = style.CATEGORICAL[4]
"""Purple — carries the write-back rail into the graph, the only edge that
runs in every phase rather than between two of them."""


@dataclass(frozen=True, slots=True)
class _Phase:
    """One box in the process: a phase, what runs it, and what gates it.

    Body lines are pre-wrapped rather than flowed: at this size the wrap
    points are a layout decision, not something to leave to a text engine.
    """

    number: int
    title: str
    prose: tuple[str, ...]
    commands: tuple[str, ...]
    gate: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _Chip:
    """The durable artifact a phase leaves behind, drawn under its box."""

    index: int
    name: str
    gloss: str


_PHASES: tuple[_Phase, ...] = (
    _Phase(
        number=1,
        title="INTAKE",
        prose=(
            "papers and primary sources",
            "become one page per concept —",
            "read by a subagent each, so",
            "every fact keeps its source",
        ),
        commands=("/add-logseq-topic",),
        gate=("draft-auditor · until CLEAN", "link check on the live page"),
    ),
    _Phase(
        number=2,
        title="DEVELOPMENT",
        prose=(
            "each page becomes a tested",
            "Python structure and a proved",
            "Lean module — the page is the",
            "spec, Mathlib the library",
        ),
        commands=("/add-python-topic", "/add-lean-topic"),
        gate=("test-auditor · just check", "statement-auditor · no sorry"),
    ),
    _Phase(
        number=3,
        title="EXPLORATION",
        prose=(
            "predictions are registered",
            "before the sweep exists; the",
            "data is graded against them,",
            "then argued over line by line",
        ),
        commands=(
            "/design-experiment",
            "/add-experiment",
            "/explore-results",
            "/critique-results",
        ),
        gate=("experiment-auditor · pilot", "results-auditor · until CLEAN"),
    ),
    _Phase(
        number=4,
        title="THEOREM PROVING",
        prose=(
            "the informal conjecture is",
            "stated formally, hunted for",
            "counterexamples beyond the",
            "swept ranges, then proved",
        ),
        commands=("/add-conjecture", "/attack-conjecture"),
        gate=("conjecture-auditor · CLEAN", "proof-auditor · axiom audit"),
    ),
    _Phase(
        number=5,
        title="PUBLISHING",
        prose=(
            "what survived is rebuilt as",
            "parameterized figures beside",
            "a docs site whose API refs",
            "are generated, never written",
        ),
        commands=("just figure <name>", "just docs"),
        gate=("docs-check · zero warnings", "no scratch plot ever ships"),
    ),
)
"""The process, left to right."""

_CHIPS: tuple[_Chip, ...] = (
    _Chip(
        index=0,
        name="logseq graph · docs/ref",
        gloss="concepts, links, and the PDFs",
    ),
    _Chip(
        index=1,
        name="src/research · src/theorems",
        gloss="tested Python, proved Lean",
    ),
    _Chip(
        index=2,
        name="data/results · docs/results",
        gloss="frames, and what they meant",
    ),
    _Chip(
        index=3,
        name="docs/conj → docs/theorems",
        gloss="what is open, and what closed",
    ),
    _Chip(
        index=4,
        name="docs/fig · site/ · the post",
        gloss="figures, docs, and the writeup",
    ),
)
"""One durable artifact per phase, in phase order."""

_WRITEBACK_SOURCE: int = 3
"""Rightmost phase feeding the write-back rail; the rail starts here and runs
left into the graph."""

_WRITEBACK_FEEDERS: tuple[int, ...] = (1, 2)
"""Phases whose artifacts join the rail on its way past them."""

_STEP_FIRST_PHASE: int = 1
"""Reveal step that adds INTAKE; each later phase follows one step behind, so
a step number below :data:`_STEP_RETURN` is also the phase's badge number."""

_STEP_RETURN: int = _STEP_FIRST_PHASE + len(_PHASES)  # one past the last phase
"""Reveal step that closes the loop with the publication path."""

_TOTAL_STEPS: int = _STEP_RETURN
"""Number of images in the build-up sequence; the last one is the whole
diagram."""


def _phase_step(index: int) -> int:
    """Return the reveal step at which the phase at ``index`` appears."""
    return _STEP_FIRST_PHASE + index


def _phase_x(index: int) -> float:
    """Return the left edge of the phase box at ``index``."""
    return _MARGIN_L + index * (_STAGE_W + _STAGE_GAP)


def _phase_cx(index: int) -> float:
    """Return the horizontal center of the phase box at ``index``."""
    return _phase_x(index) + _STAGE_W / 2.0


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
    """Draw a hairline across the phase box at ``index``, at height ``y``."""
    x = _phase_x(index)
    ax.plot(
        [x + 0.2, x + _STAGE_W - 0.2],
        [y, y],
        color=style.PARCHMENT,
        linewidth=1.0,
        solid_capstyle="butt",
    )


def _chip(ax: Axes, index: int, title: str, gloss: str, bottom: float) -> None:
    """Draw a titled chip of :data:`_INPUT_H` height across a phase's column."""
    _box(
        ax,
        (_phase_x(index), bottom, _STAGE_W, _INPUT_H),
        facecolor=style.PAPER,
        edgecolor=style.MIST,
        linewidth=1.0,
    )
    ax.text(
        _phase_cx(index),
        bottom + 0.58,
        title,
        fontsize=9.5,
        fontweight="bold",
        color=style.INK,
        ha="center",
        va="center",
    )
    ax.text(
        _phase_cx(index),
        bottom + 0.28,
        gloss,
        fontsize=8.5,
        color=style.MIST,
        ha="center",
        va="center",
    )


def _source_input(ax: Axes) -> None:
    """Draw what starts the process, above INTAKE, and its arrow into it."""
    _chip(ax, 0, "the source", "a paper, a spec, a curiosity", _INPUT_BOTTOM)
    _arrow(
        ax,
        (_phase_cx(0), _INPUT_BOTTOM),
        (_phase_cx(0), _STAGE_TOP),
        color=style.SLATE,
    )


def _published_output(ax: Axes) -> None:
    """Draw what ships, above PUBLISHING, and the arrow up into it."""
    _chip(ax, 4, "the writeup", "a post, a paper, a talk", _INPUT_BOTTOM)
    _arrow(
        ax,
        (_phase_cx(4), _STAGE_TOP),
        (_phase_cx(4), _INPUT_BOTTOM),
        color=style.SLATE,
    )


def _phase(ax: Axes, index: int, phase: _Phase) -> None:
    """Draw the phase box at ``index``: badge, title, prose, commands, gate."""
    x = _phase_x(index)
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
        str(phase.number),
        fontsize=9.5,
        fontweight="bold",
        color=style.PAPER,
        ha="center",
        va="center",
    )
    ax.text(
        x + 0.64,
        badge_y,
        phase.title,
        fontsize=12,
        fontweight="bold",
        color=style.INK,
        va="center",
    )
    _rule(ax, index, _HEAD_RULE_Y)

    for line_index, line in enumerate(phase.prose):
        ax.text(
            x + 0.2,
            _PROSE_TOP_Y - line_index * _LINE_STEP,
            line,
            fontsize=8.5,
            color=style.INK,
            va="center",
        )

    top = _COMMAND_CENTER_Y + (len(phase.commands) - 1) * _LINE_STEP / 2.0
    for line_index, command in enumerate(phase.commands):
        ax.text(
            x + 0.2,
            top - line_index * _LINE_STEP,
            command,
            fontsize=8.5,
            fontweight="bold",
            color=style.SLATE,
            va="center",
        )

    _rule(ax, index, _GATE_RULE_Y)
    for line_index, line in enumerate(phase.gate):
        ax.text(
            x + 0.2,
            _GATE_TOP_Y - line_index * _LINE_STEP,
            line,
            fontsize=8,
            color=_AUDIT_COLOR,
            va="center",
        )


def _flow(ax: Axes, step: int) -> None:
    """Draw the left-to-right arrows between consecutive revealed phases."""
    mid_y = _STAGE_BOTTOM + _STAGE_H / 2.0
    for index in range(len(_PHASES) - 1):
        if step < _phase_step(index + 1):
            continue
        _arrow(
            ax,
            (_phase_x(index) + _STAGE_W, mid_y),
            (_phase_x(index + 1), mid_y),
            color=style.SLATE,
        )


def _artifacts(ax: Axes, step: int) -> None:
    """Draw the artifact chips, each with the phase that leaves it behind."""
    for chip in _CHIPS:
        if step < _phase_step(chip.index):
            continue
        _box(
            ax,
            (_phase_x(chip.index), _FILE_TOP - _FILE_H, _STAGE_W, _FILE_H),
            facecolor=style.PAPER,
            edgecolor=style.MIST,
            linewidth=1.0,
        )
        ax.plot(
            [_phase_cx(chip.index), _phase_cx(chip.index)],
            [_STAGE_BOTTOM, _FILE_TOP],
            color=style.MIST,
            linewidth=0.9,
            linestyle=(0, (2, 2)),
        )
        ax.text(
            _phase_cx(chip.index),
            _FILE_TOP - 0.24,
            chip.name,
            fontsize=9,
            fontweight="bold",
            color=style.INK,
            ha="center",
            va="center",
        )
        ax.text(
            _phase_cx(chip.index),
            _FILE_TOP - 0.48,
            chip.gloss,
            fontsize=8,
            color=style.MIST,
            ha="center",
            va="center",
        )


def _bug_first_return(ax: Axes, step: int) -> None:
    """Draw the return from a literature contradiction to the code."""
    if step < _phase_step(2):
        return
    start_x = _phase_x(2) + _BUG_X_OFFSET
    end_x = _phase_x(1) + _BUG_X_OFFSET
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


def _refutation_return(ax: Axes, step: int) -> None:
    """Draw the return from a counterexample to the results that suggested it.

    It nests above the bug-first path and lands on a different column of the
    box they share, so the two never run along the same vertical.
    """
    if step < _phase_step(3):
        return
    start_x = _phase_x(3) + _REFUTE_X_OFFSET
    end_x = _phase_x(2) + _REFUTE_X_OFFSET
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


def _writeback_rail(ax: Axes, step: int) -> None:
    """Draw the rail every building phase writes its outcome back along.

    It runs under the artifact chips and lands on the graph rather than on a
    phase box: what a skill records at the end of a run is a page, and the
    next phase to read that page may be any of them.
    """
    if step < _phase_step(_WRITEBACK_SOURCE):
        return
    chip_bottom = _FILE_TOP - _FILE_H
    for index in _WRITEBACK_FEEDERS:
        ax.plot(
            [_phase_cx(index), _phase_cx(index)],
            [chip_bottom, _WRITEBACK_Y],
            color=_WRITEBACK_COLOR,
            linewidth=1.2,
        )
    _routed_arrow(
        ax,
        [
            (_phase_cx(_WRITEBACK_SOURCE), chip_bottom),
            (_phase_cx(_WRITEBACK_SOURCE), _WRITEBACK_Y),
            (_phase_cx(0), _WRITEBACK_Y),
            (_phase_cx(0), chip_bottom),
        ],
        color=_WRITEBACK_COLOR,
    )


def _return_loop(ax: Axes, step: int) -> None:
    """Draw the path from what was published back to the next intake.

    It climbs the left spine rather than the phase centers, so it clears the
    artifact chips, and it lands on the source itself: the terms a writeup
    needs and the questions it leaves open are the next thing read.
    """
    if step < _STEP_RETURN:
        return
    _routed_arrow(
        ax,
        [
            (_phase_cx(4), _FILE_TOP - _FILE_H),
            (_phase_cx(4), _RETURN_Y),
            (_SPINE_X, _RETURN_Y),
            (_SPINE_X, _INPUT_BOTTOM + _INPUT_H / 2.0),
            (_phase_x(0), _INPUT_BOTTOM + _INPUT_H / 2.0),
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

    # The source arrives with the phase that turns it into a page, and the
    # writeup with the phase that produces it.
    if step >= _phase_step(0):
        _source_input(ax)
    if step >= _phase_step(4):
        _published_output(ax)
    for index, phase in enumerate(_PHASES):
        if step >= _phase_step(index):
            _phase(ax, index, phase)
    _flow(ax, step)
    _artifacts(ax, step)
    _bug_first_return(ax, step)
    _refutation_return(ax, step)
    _writeback_rail(ax, step)
    _return_loop(ax, step)
    return fig


@figure(name="process-overview")
def process_overview(ctx: FigureContext) -> None:
    """Render the whole-process diagram, whole and as a build-up sequence.

    Writes ``process-overview`` (the complete graphic) plus
    ``process-overview-NN`` for each reveal step, numbered in presentation
    order.
    """
    # The house style saves with a tight bbox, which would re-crop and re-pad
    # the frame that was just sized deliberately — and would crop each step to
    # its own contents, so the sequence would jitter. Write them as laid out.
    with mpl.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        for step in range(1, _TOTAL_STEPS + 1):
            fig = _render(step)
            ctx.save(fig, f"process-overview-{step:02d}")
            plt.close(fig)

        fig = _render(_TOTAL_STEPS)
        ctx.save(fig, "process-overview")
        plt.close(fig)
