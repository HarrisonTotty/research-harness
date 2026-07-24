"""Research experiments composed from the :mod:`research` library.

Besides the individual experiment modules, this package provides the shared
harness they build on: the :func:`experiment` command decorator and its
:class:`ExperimentContext`, logging configuration, and result/metadata I/O.
"""

from experiments.cli import ExperimentContext, experiment
from experiments.io import (
    default_metadata_path,
    default_result_path,
    require_json_destination,
    write_metadata,
    write_result,
)
from experiments.logging import LogLevel, configure_logging

__all__ = [
    "ExperimentContext",
    "LogLevel",
    "configure_logging",
    "default_metadata_path",
    "default_result_path",
    "experiment",
    "require_json_destination",
    "write_metadata",
    "write_result",
]
