"""Standardized command-line surface for experiments.

The :func:`experiment` decorator turns a plain function into a :mod:`click`
command carrying the four options every experiment shares — ``--log-level``,
``--log-file``, ``--out``, and ``--meta-out`` — and hands the wrapped function
an :class:`ExperimentContext` as its first argument. The context bundles the
configured logger and the run timestamp with pre-resolved JSON result and
metadata destinations, so an experiment body only declares its own parameters
and calls :meth:`ExperimentContext.write_result` /
:meth:`ExperimentContext.write_metadata`.
"""

import functools
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from logging import Logger
from pathlib import Path

import click
import pandas as pd

from experiments import io
from experiments.logging import LogLevel, configure_logging

_LOG_LEVEL_VALUES: list[str] = [level.value for level in LogLevel]
"""Accepted ``--log-level`` strings, in declaration order."""

_TIMESTAMP_FORMAT: str = "%Y%m%dT%H%M%SZ"
"""Compact, filesystem-safe UTC stamp used in default output file names."""


def _default_clock() -> datetime:
    """Return the current UTC time, the default source of a run's timestamp."""
    return datetime.now(UTC)


def _format_timestamp(moment: datetime) -> str:
    """Render ``moment`` as a compact UTC stamp for use in a file name."""
    return moment.astimezone(UTC).strftime(_TIMESTAMP_FORMAT)


@dataclass(frozen=True, slots=True)
class ExperimentContext:
    """Resolved run environment handed to an experiment body.

    Instances are constructed by :func:`experiment` from the standard
    command-line options; experiments never build one directly.
    """

    name: str
    logger: Logger
    timestamp: str
    result_path: Path
    metadata_path: Path

    def write_result(self, data: pd.DataFrame) -> Path:
        """Write ``data`` to the resolved result path and return that path."""
        return io.write_result(data, self.result_path)

    def write_metadata(self, metadata: Mapping[str, object]) -> Path:
        """Write ``metadata`` to the resolved metadata path and return it."""
        return io.write_metadata(metadata, self.metadata_path)


def _standard_options[**P, R](func: Callable[P, R]) -> Callable[P, R]:
    """Attach the four shared experiment options to ``func``."""
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
            "-o",
            "--out",
            type=click.Path(dir_okay=False, path_type=Path),
            default=None,
            help="Override the result destination (must be a .json file).",
        ),
        click.option(
            "-O",
            "--meta-out",
            type=click.Path(dir_okay=False, path_type=Path),
            default=None,
            help="Override the run-metadata destination (must be a .json file).",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def _require_json(path: Path | None, param_hint: str) -> None:
    """Reject a non-JSON override path with a click usage error.

    Args:
        path: The overridden destination, or ``None`` when not overridden.
        param_hint: Option spelling to name in the error message.

    Raises:
        click.BadParameter: If ``path`` is given but is not a ``.json`` file.
    """
    if path is not None and path.suffix.lower() != io.RESULT_SUFFIX:
        msg = "experiment output is written only as JSON; use a .json path"
        raise click.BadParameter(msg, param_hint=param_hint)


def experiment(
    name: str | None = None,
    *,
    results_dir: Path = io.DEFAULT_RESULTS_DIR,
    clock: Callable[[], datetime] = _default_clock,
) -> Callable[[Callable[..., object]], click.Command]:
    """Turn an experiment function into a :mod:`click` command.

    The decorated function is invoked with a freshly built
    :class:`ExperimentContext` as its first positional argument, followed by
    keyword arguments for any :mod:`click` options the function declares itself.
    The four standard options are consumed here and do not reach the function.

    Args:
        name: Command name and stem for default output paths; defaults to the
            function's own name with underscores turned into hyphens.
        results_dir: Directory holding the default, timestamped result and
            metadata paths.
        clock: Source of the run timestamp; called once per invocation. Injected
            for testing, and defaults to the current UTC time.

    Returns:
        A decorator producing the configured :class:`click.Command`.
    """

    def decorator(func: Callable[..., object]) -> click.Command:
        exp_name = name if name is not None else func.__name__.replace("_", "-")

        @functools.wraps(func)
        def callback(
            *,
            log_level: str,
            log_file: Path | None,
            out: Path | None,
            meta_out: Path | None,
            **kwargs: object,
        ) -> object:
            _require_json(out, "'-o' / '--out'")
            _require_json(meta_out, "'-O' / '--meta-out'")
            logger = configure_logging(level=LogLevel(log_level), log_file=log_file)
            timestamp = _format_timestamp(clock())
            result_path = (
                out
                if out is not None
                else io.default_result_path(exp_name, timestamp, results_dir)
            )
            metadata_path = (
                meta_out
                if meta_out is not None
                else io.default_metadata_path(exp_name, timestamp, results_dir)
            )
            context = ExperimentContext(
                name=exp_name,
                logger=logger,
                timestamp=timestamp,
                result_path=result_path,
                metadata_path=metadata_path,
            )
            return func(context, **kwargs)

        # ``functools.wraps`` copied any click options the experiment declared
        # onto ``callback``; copy the list so appending the standard options does
        # not mutate the original function's parameters.
        callback.__click_params__ = list(  # type: ignore[attr-defined]  # click's private param store
            getattr(func, "__click_params__", [])
        )
        return click.command(name=exp_name)(_standard_options(callback))

    return decorator
