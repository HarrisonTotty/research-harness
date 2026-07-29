"""Demonstration sheet of the blog figure style.

Renders one figure with four panels — multi-series lines, a single-series bar
chart, a sequential heatmap, and a diverging heatmap — exercising the
categorical slots, both colormaps, and the chrome defaults from
:mod:`figures.style`. Regenerate after any style change to eyeball the result:
``just figure palette-demo``.
"""

import math

from matplotlib import pyplot as plt
from matplotlib.axes import Axes

from figures import style
from figures.cli import FigureContext, figure

_SERIES: tuple[tuple[str, float], ...] = (
    ("alpha", 1.0),
    ("beta", 1.9),
    ("gamma", 2.8),
    ("delta", 3.7),
)
"""Demo line series: label and angular frequency of a damped cosine."""


def _damped(omega: float, x: float) -> float:
    """Return the damped cosine ``exp(-x/3) * cos(omega * x)``."""
    return math.exp(-x / 3.0) * math.cos(omega * x)


def _lines_panel(ax: Axes) -> None:
    """Draw the multi-series line panel.

    The series all decay toward zero, so end-of-line direct labels would
    collide; the legend carries identity instead.
    """
    xs = [i * 6.0 / 120 for i in range(121)]
    for label, omega in _SERIES:
        ax.plot(xs, [_damped(omega, x) for x in xs], label=label)
    ax.set_title("categorical · lines")
    ax.set_xlim(0, 6.0)
    ax.legend(loc="upper right", fontsize=9)


def _bars_panel(ax: Axes) -> None:
    """Draw the single-series bar panel in the slot-1 hue."""
    labels = ["A", "B", "C", "D", "E"]
    values = [3.0, 7.0, 4.0, 6.0, 5.0]
    ax.bar(labels, values, color=style.CATEGORICAL[0], width=0.6)
    ax.set_title("nominal · bars")


def _sequential_panel(ax: Axes) -> None:
    """Draw the sequential heatmap panel (a smooth non-negative field)."""
    n = 40
    field = [
        [
            math.exp(-(((col / n) - 0.6) ** 2 + ((row / n) - 0.4) ** 2) * 6.0)
            for col in range(n)
        ]
        for row in range(n)
    ]
    image = ax.imshow(field, cmap=style.SEQUENTIAL_CMAP, origin="lower")
    ax.grid(visible=False)
    ax.set_title("sequential · heatmap")
    colorbar = ax.figure.colorbar(image, ax=ax)
    colorbar.outline.set_edgecolor(style.MIST)


def _diverging_panel(ax: Axes) -> None:
    """Draw the diverging heatmap panel (a signed field about zero)."""
    n = 40
    field = [
        [
            math.sin(2.0 * math.pi * col / n) * math.cos(2.0 * math.pi * row / n)
            for col in range(n)
        ]
        for row in range(n)
    ]
    image = ax.imshow(
        field, cmap=style.DIVERGING_CMAP, origin="lower", vmin=-1.0, vmax=1.0
    )
    ax.grid(visible=False)
    ax.set_title("diverging · heatmap")
    colorbar = ax.figure.colorbar(image, ax=ax)
    colorbar.outline.set_edgecolor(style.MIST)


@figure(name="palette-demo")
def palette_demo(ctx: FigureContext) -> None:
    """Render the four-panel style sheet and save it as ``palette-demo``."""
    fig = plt.figure(figsize=(9.5, 7.0), layout="constrained")
    _lines_panel(fig.add_subplot(2, 2, 1))
    _bars_panel(fig.add_subplot(2, 2, 2))
    _sequential_panel(fig.add_subplot(2, 2, 3))
    _diverging_panel(fig.add_subplot(2, 2, 4))
    fig.suptitle("blog figure style")
    ctx.save(fig, "palette-demo")
    plt.close(fig)
