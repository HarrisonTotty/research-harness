"""Logging configuration shared by every experiment entry point.

Experiments are applications, so configuring the root logger here is
appropriate: it captures both experiment code and anything the :mod:`research`
library logs. :func:`configure_logging` is idempotent — it only removes and
replaces the handlers it previously installed, tagged with
:data:`_MANAGED_ATTR`, so repeated calls (and a test harness's own handlers)
are left undisturbed.
"""

import logging
import sys
from enum import StrEnum
from pathlib import Path
from typing import TextIO

_MANAGED_ATTR: str = "_research_managed"
"""Attribute stamped on handlers this module installs, so it can replace them."""

_LOG_FORMAT: str = "%(asctime)s %(levelname)-8s %(name)s: %(message)s"
"""Format applied to both the stream and file handlers."""


class LogLevel(StrEnum):
    """Selectable logging verbosity, named for readable command-line values."""

    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"

    @property
    def numeric(self) -> int:
        """Return the :mod:`logging` integer threshold for this level."""
        return logging.getLevelNamesMapping()[self.name]


def _managed_handlers(logger: logging.Logger) -> list[logging.Handler]:
    """Return the handlers on ``logger`` that this module installed."""
    return [h for h in logger.handlers if getattr(h, _MANAGED_ATTR, False)]


def configure_logging(
    *,
    level: LogLevel = LogLevel.INFO,
    log_file: Path | None = None,
    stream: TextIO | None = None,
) -> logging.Logger:
    """Configure and return the root logger for an experiment run.

    Installs a stream handler (writing to ``stream``, defaulting to standard
    error) and, when ``log_file`` is given, an additional file handler that
    mirrors the same records to that path — its parent directories are created
    if missing. Any handlers installed by a previous call are removed first, so
    invoking this repeatedly does not duplicate output.

    Args:
        level: Threshold below which records are dropped.
        log_file: Optional path to also write log records to, in addition to the
            stream. Encoded as UTF-8.
        stream: Destination for the stream handler; standard error when omitted.

    Returns:
        The configured root logger.
    """
    root = logging.getLogger()
    root.setLevel(level.numeric)
    for handler in _managed_handlers(root):
        root.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter(_LOG_FORMAT)

    stream_handler = logging.StreamHandler(stream if stream is not None else sys.stderr)
    stream_handler.setFormatter(formatter)
    setattr(stream_handler, _MANAGED_ATTR, True)
    root.addHandler(stream_handler)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)
        setattr(file_handler, _MANAGED_ATTR, True)
        root.addHandler(file_handler)

    return root
