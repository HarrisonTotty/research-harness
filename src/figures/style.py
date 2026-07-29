"""Blog-standard color scheme and matplotlib styling for figures.

Every color here is taken verbatim from the palette of
``harrisontotty.github.io`` (``assets/style.css``) — the chrome roles
(:data:`INK`, :data:`PAPER`, :data:`PARCHMENT`, :data:`MIST`, :data:`SLATE`),
the :data:`CATEGORICAL` series slots, and the poles of both colormaps are the
blog's exact hexes; ramp tints between those stops are interpolated in OKLCH
on the blog hues. Exact brand match is the deliberate priority: the blog's
accents are muted enough that the series palette does **not** clear the usual
colorblind-separation gates (see :data:`CATEGORICAL` for the numbers and the
obligations that follow), whereas contrast against :data:`PAPER` (WCAG >= 3:1)
does hold for every slot.

Use :func:`context` (or :func:`apply`) to activate the style;
:func:`figures.cli.figure` does so automatically around every figure body.
"""

from collections.abc import Generator
from contextlib import contextmanager

import matplotlib as mpl
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.typing import RcKeyType

INK: str = "#1d1b1b"
"""Near-black text color of the blog; primary ink for titles and labels."""

PAPER: str = "#ffffff"
"""The blog's page background; figure and axes surface."""

PARCHMENT: str = "#ddd8cf"
"""The blog's warm parchment (``pre`` background); gridline color."""

MIST: str = "#abafb6"
"""The blog's muted gray-blue (rules, footers); axis spines and ticks."""

SLATE: str = "#60728d"
"""The blog's slate-blue link color; hue anchor of the sequential ramp."""

CATEGORICAL: tuple[str, ...] = (
    "#60728d",  # slate (links)
    "#608d7b",  # teal (numeric/constant literals)
    "#8d7260",  # brown (comments)
    "#698d60",  # green (decorators)
    "#84608d",  # purple (entities/symbols)
    "#cd5c5c",  # indianred (errors, falsified callouts)
)
"""Series colors in fixed assignment order — assign in sequence, never cycle.

These are the blog's accent hexes verbatim, ordered slate-first (the blog's
signature link color) with the same-family neighbors separated as far as the
set allows. Exact brand match wins over palette validation here, and that has
a real cost: the set fails the categorical colorblind gates (worst adjacent
pair under simulated deuteranopia is indianred/green at OKLab ΔE 2.8, and the
worst adjacent normal-vision pair is teal/brown at ΔE 8.9, both far below the
8/15 floors). Color therefore must never carry series identity alone — pair
it with secondary encoding: a legend plus direct labels, distinct markers or
linestyles, or facets. Keep scatter/bubble figures to two or three series at
most, since there *every* pair can collide, not just adjacent ones.
"""

SEQUENTIAL_ANCHORS: tuple[str, ...] = (
    "#ebeff4",
    "#b5becd",
    "#8190a6",
    "#60728d",
    "#4c6181",
    "#385074",
    "#314769",
    "#263a59",
    "#1c2e49",
)
"""One-hue magnitude ramp, light to dark, through the blog's exact blues.

Interpolation anchors for :data:`SEQUENTIAL_CMAP`: a near-white tint falling
monotonically in lightness through :data:`SLATE` (``#60728d``) and the blog's
deep blue (``#385074``, its type/class syntax color) into a darker tail, all
on the slate hue angle at the blog's own muted chroma. The light end recedes
toward the surface by design — it means "near zero".
"""

ORDINAL_STEPS: tuple[str, ...] = (
    "#8190a6",
    "#60728d",
    "#385074",
    "#1c2e49",
)
"""Discrete ordered-category steps (funnel stages, tiers), light to dark.

Four steps from the slate ramp — the middle two are the blog's exact
``#60728d`` and ``#385074`` — validated as an ordinal ramp: monotone
lightness, adjacent ΔL >= 0.06, and a light end at 3.24:1 on :data:`PAPER`
so every step reads as a mark. Use slices from the dark end when fewer steps
are needed; more than four would squeeze the gaps below the visible floor.
"""

DIVERGING_ANCHORS: tuple[str, ...] = (
    "#cd5c5c",
    "#d9817e",
    "#e3a4a1",
    "#ebc6c4",
    "#eeebe5",
    "#c4ccd7",
    "#a2adbe",
    "#808fa5",
    "#60728d",
)
"""Polarity ramp: exact blog indianred pole to exact slate-blue pole.

Interpolation anchors for :data:`DIVERGING_CMAP`. The poles are the blog's
``#cd5c5c`` and ``#60728d`` verbatim (matching the falsified/confirmed
callout semantics of the blog, red = negative); each arm fades toward a
near-neutral parchment-hue midpoint so "zero" reads as nothing.
"""

SEQUENTIAL_CMAP: LinearSegmentedColormap = LinearSegmentedColormap.from_list(
    "blog-sequential", SEQUENTIAL_ANCHORS
)
"""Continuous magnitude colormap built from :data:`SEQUENTIAL_ANCHORS`."""

DIVERGING_CMAP: LinearSegmentedColormap = LinearSegmentedColormap.from_list(
    "blog-diverging", DIVERGING_ANCHORS
)
"""Continuous polarity colormap built from :data:`DIVERGING_ANCHORS`."""

CATEGORICAL_CMAP: ListedColormap = ListedColormap(CATEGORICAL, name="blog-categorical")
"""Discrete colormap over the :data:`CATEGORICAL` slots, in slot order."""

_MONOSPACE_STACK: list[str] = [
    "Roboto Mono",
    "DejaVu Sans Mono",
    "Consolas",
    "Courier New",
]
"""Monospace preference order matching the blog's font stack."""


def register_colormaps() -> None:
    """Register the blog colormaps with matplotlib, once.

    Makes ``"blog-sequential"``, ``"blog-diverging"``, and
    ``"blog-categorical"`` resolvable by name (e.g. via ``cmap=`` strings or
    ``image.cmap``). Safe to call repeatedly; already-registered names are
    left untouched.
    """
    for cmap in (SEQUENTIAL_CMAP, DIVERGING_CMAP, CATEGORICAL_CMAP):
        if cmap.name not in mpl.colormaps:
            mpl.colormaps.register(cmap)


def rc_params() -> dict[RcKeyType, object]:
    """Return the rcParams overrides that realize the blog style.

    The values pair the blog's chrome colors with the standard mark specs:
    monospace type, recessive spines and horizontal-only gridlines, 2pt lines,
    frameless legends, and the :data:`CATEGORICAL` property cycle. Feed the
    result to :func:`matplotlib.rc_context` or ``rcParams.update`` — or use
    :func:`context` / :func:`apply`, which also register the colormaps.
    """
    return {
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "savefig.bbox": "tight",
        "font.family": "monospace",
        "font.monospace": _MONOSPACE_STACK,
        "font.size": 11.0,
        "text.color": INK,
        "axes.facecolor": PAPER,
        "axes.edgecolor": MIST,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.titlelocation": "left",
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "axes.prop_cycle": cycler(color=list(CATEGORICAL)),
        "grid.color": PARCHMENT,
        "grid.linewidth": 0.8,
        "xtick.color": MIST,
        "ytick.color": MIST,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "lines.linewidth": 2.0,
        "lines.markersize": 7.0,
        "legend.frameon": False,
        "figure.titleweight": "bold",
        "image.cmap": "blog-sequential",
    }


@contextmanager
def context() -> Generator[None]:
    """Activate the blog style within a ``with`` block, then restore.

    Registers the colormaps and applies :func:`rc_params` through
    :func:`matplotlib.rc_context`, so any rcParams the block changes are
    rolled back on exit. Colormap registration persists (it is idempotent and
    additive, not stateful styling).
    """
    register_colormaps()
    with mpl.rc_context(rc_params()):
        yield


def apply() -> None:
    """Activate the blog style globally for the running process.

    Prefer :func:`context` where scoping matters (tests, notebooks that mix
    styles); this variant suits scripts that render a single figure.
    """
    register_colormaps()
    mpl.rcParams.update(rc_params())
