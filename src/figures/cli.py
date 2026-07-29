"""Standardized command-line surface for figure generators.

The :func:`figure` decorator turns a plain function into a :mod:`click`
command carrying the five options every figure shares — ``--log-level``,
``--log-file``, ``--out-dir``, ``--format``, and ``--dpi`` — and hands the
wrapped function a :class:`FigureContext` as its first argument. The context
bundles the configured logger and run timestamp with the resolved output
directory and rendering settings, so a figure body only declares its own
parameters, builds its matplotlib figures under the blog style (applied
automatically via :func:`figures.style.context`), and calls
:meth:`FigureContext.save`.
"""

import functools
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from logging import Logger
from pathlib import Path

import click
from matplotlib.figure import Figure

from experiments.logging import LogLevel, configure_logging
from figures import style

DEFAULT_FIGURES_DIR: Path = Path("docs/fig")
"""Directory holding generated figures, relative to the repository root."""


class FigureFormat(StrEnum):
    """Selectable output file formats for a saved figure."""

    PNG = "png"
    SVG = "svg"
    PDF = "pdf"


DEFAULT_FORMATS: tuple[FigureFormat, ...] = (FigureFormat.PNG, FigureFormat.SVG)
"""Formats written when ``--format`` is not given: raster for previews and
papers, vector for blog embedding."""

DEFAULT_DPI: int = 300
"""Raster resolution used when ``--dpi`` is not given."""

_LOG_LEVEL_VALUES: list[str] = [level.value for level in LogLevel]
"""Accepted ``--log-level`` strings, in declaration order."""

_FORMAT_VALUES: list[str] = [fmt.value for fmt in FigureFormat]
"""Accepted ``--format`` strings, in declaration order."""

_TIMESTAMP_FORMAT: str = "%Y%m%dT%H%M%SZ"
"""Compact, filesystem-safe UTC stamp recorded for the run."""

# `...` is an explicit `Any` parameter list: a figure body declares its own
# click options, so its signature is a genuinely dynamic boundary.
type FigureFunc = Callable[..., object]  # type: ignore[explicit-any]
"""A figure body: a context plus whatever options it declares itself."""


def _default_clock() -> datetime:
    """Return the current UTC time, the default source of a run's timestamp."""
    return datetime.now(UTC)


def _format_timestamp(moment: datetime) -> str:
    """Render ``moment`` as a compact UTC stamp."""
    return moment.astimezone(UTC).strftime(_TIMESTAMP_FORMAT)


@dataclass(frozen=True, slots=True)
class FigureContext:
    """Resolved run environment handed to a figure body.

    Instances are constructed by :func:`figure` from the standard
    command-line options; figure bodies never build one directly.
    """

    name: str
    logger: Logger
    timestamp: str
    output_dir: Path
    formats: tuple[FigureFormat, ...]
    dpi: int

    def save(self, fig: Figure, stem: str) -> list[Path]:
        """Write ``fig`` into the output directory, once per format.

        Args:
            fig: The rendered matplotlib figure.
            stem: File name without extension; each format appends its own
                suffix (``<stem>.png``, ``<stem>.svg``, ...). Use stable stems
                — figures are referenced by path from posts and papers, so
                re-running a generator overwrites in place.

        Returns:
            The paths written, in :attr:`formats` order.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        paths: list[Path] = []
        for fmt in self.formats:
            path = self.output_dir / f"{stem}.{fmt.value}"
            fig.savefig(path, format=fmt.value, dpi=self.dpi)
            self.logger.info("wrote %s", path)
            paths.append(path)
        return paths


def _standard_options[**P, R](func: Callable[P, R]) -> Callable[P, R]:
    """Attach the five shared figure options to ``func``."""
    options = [
        click.option(
            "-l",
            "--log-level",
            type=click.Choice(_LOG_LEVEL_VALUES, case_sensitive=False),
            default=LogLevel.INFO.value,
            show_default=True,
            help="Minimum severity of log records emitted.",
        ),
        click.option(
            "-L",
            "--log-file",
            type=click.Path(dir_okay=False, path_type=Path),
            default=None,
            help="Also write log records to this file, in addition to stderr.",
        ),
        click.option(
            "-d",
            "--out-dir",
            type=click.Path(file_okay=False, path_type=Path),
            default=None,
            help="Override the output directory (default: docs/fig/<name>).",
        ),
        click.option(
            "-f",
            "--format",
            "formats",
            type=click.Choice(_FORMAT_VALUES, case_sensitive=False),
            multiple=True,
            help="Output format; repeat for several. [default: png, svg]",
        ),
        click.option(
            "--dpi",
            type=click.IntRange(min=1),
            default=DEFAULT_DPI,
            show_default=True,
            help="Raster resolution in dots per inch.",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def figure(
    name: str | None = None,
    *,
    figures_dir: Path = DEFAULT_FIGURES_DIR,
    clock: Callable[[], datetime] = _default_clock,
) -> Callable[[FigureFunc], click.Command]:
    """Turn a figure-generating function into a :mod:`click` command.

    The decorated function is invoked with a freshly built
    :class:`FigureContext` as its first positional argument, followed by
    keyword arguments for any :mod:`click` options the function declares
    itself. The five standard options are consumed here and do not reach the
    function, and the body runs inside :func:`figures.style.context` so every
    figure it builds wears the blog style without further setup.

    Args:
        name: Command name and stem of the default output directory; defaults
            to the function's own name with underscores turned into hyphens.
        figures_dir: Directory whose ``<name>`` subdirectory is the default
            output destination.
        clock: Source of the run timestamp; called once per invocation.
            Injected for testing, and defaults to the current UTC time.

    Returns:
        A decorator producing the configured :class:`click.Command`.
    """

    def decorator(func: FigureFunc) -> click.Command:
        fig_name = name if name is not None else func.__name__.replace("_", "-")

        @functools.wraps(func)
        def callback(
            *,
            log_level: str,
            log_file: Path | None,
            out_dir: Path | None,
            formats: tuple[str, ...],
            dpi: int,
            **kwargs: object,
        ) -> object:
            logger = configure_logging(level=LogLevel(log_level), log_file=log_file)
            context = FigureContext(
                name=fig_name,
                logger=logger,
                timestamp=_format_timestamp(clock()),
                output_dir=out_dir if out_dir is not None else figures_dir / fig_name,
                formats=tuple(FigureFormat(fmt) for fmt in formats) or DEFAULT_FORMATS,
                dpi=dpi,
            )
            with style.context():
                return func(context, **kwargs)

        # ``functools.wraps`` copied any click options the figure declared onto
        # ``callback``; copy the list so appending the standard options does
        # not mutate the original function's parameters.
        callback.__click_params__ = list(  # type: ignore[attr-defined]  # click's private param store
            getattr(func, "__click_params__", [])
        )
        return click.command(name=fig_name)(_standard_options(callback))

    return decorator
