"""Shared matplotlib plumbing for the structure modules' plot methods.

Every plot method in :mod:`research.matroid` and :mod:`research.positroid`
follows the same contract — draw onto given axes or a fresh figure, never
call ``show`` — so the axes bootstrap, the labeled scatter, and the circular
layouts live here once.
"""

import math
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "ensure_axes",
    "scatter_labeled",
    "unit_circle",
]


def ensure_axes(ax: Axes | None) -> Axes:
    """Return ``ax``, or the axes of a fresh figure when none is given."""
    if ax is not None:
        return ax
    # Deferred so core use never imports the plotting stack.
    import matplotlib.pyplot as plt  # noqa: PLC0415

    _, fresh = plt.subplots()
    return fresh


def scatter_labeled(
    ax: Axes,
    points: Sequence[tuple[float, float]],
    texts: Iterable[str],
    offsets: Iterable[tuple[float, float]] | None = None,
) -> None:
    """Scatter the points and annotate each with its text.

    Offsets are in display points, one per point, letting circular layouts
    push labels radially outward; the default nudges every label upward.
    """
    if offsets is None:
        offsets = [(0.0, 6.0)] * len(points)
    if points:
        ax.scatter([x for x, _ in points], [y for _, y in points], zorder=2)
    for text, (x, y), offset in zip(texts, points, offsets, strict=True):
        ax.annotate(
            text,
            (x, y),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            va="center",
            fontsize=8,
        )


def unit_circle(
    n: int, *, phase: float = 0.0, clockwise: bool = False
) -> list[tuple[float, float]]:
    """Return ``n`` evenly spaced points on the unit circle.

    Args:
        n: Number of points; zero yields an empty layout.
        phase: Angle of the first point, in radians.
        clockwise: Walk clockwise instead of the mathematical direction.
    """
    orientation = -1.0 if clockwise else 1.0
    return [
        (
            math.cos(phase + orientation * 2 * math.pi * i / n),
            math.sin(phase + orientation * 2 * math.pi * i / n),
        )
        for i in range(n)
    ]
