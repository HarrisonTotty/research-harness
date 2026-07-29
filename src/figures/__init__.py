"""Figure generation for papers and blog posts.

Besides the individual figure modules, this package provides the shared
harness they build on: the :func:`figure` command decorator and its
:class:`FigureContext`, plus :mod:`figures.style` — the blog-derived color
scheme and matplotlib styling every figure is rendered under.
"""

from figures.cli import (
    DEFAULT_DPI,
    DEFAULT_FIGURES_DIR,
    DEFAULT_FORMATS,
    FigureContext,
    FigureFormat,
    figure,
)

__all__ = [
    "DEFAULT_DPI",
    "DEFAULT_FIGURES_DIR",
    "DEFAULT_FORMATS",
    "FigureContext",
    "FigureFormat",
    "figure",
]
