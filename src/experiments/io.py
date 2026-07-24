"""Result and metadata serialization for experiments.

Every experiment run writes two UTF-8 JSON files under
:data:`DEFAULT_RESULTS_DIR`: the result — a serialized :class:`pandas.DataFrame`
— and a metadata sidecar recording the run's parameters and provenance. JSON is
the only supported format, so both destinations must carry a ``.json`` suffix.
Default file names are stamped with the run timestamp
(``{name}.{timestamp}.json`` and ``{name}.{timestamp}.meta.json``) so repeated
runs do not overwrite one another, matching the repository README.
"""

import json
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

DEFAULT_RESULTS_DIR: Path = Path("data/results")
"""Directory holding raw experiment results, relative to the repository root."""

RESULT_SUFFIX: str = ".json"
"""Required extension for every result and metadata destination."""


def default_result_path(
    name: str, timestamp: str, results_dir: Path = DEFAULT_RESULTS_DIR
) -> Path:
    """Return the default, timestamped result path for experiment ``name``."""
    return results_dir / f"{name}.{timestamp}{RESULT_SUFFIX}"


def default_metadata_path(
    name: str, timestamp: str, results_dir: Path = DEFAULT_RESULTS_DIR
) -> Path:
    """Return the default, timestamped metadata path for experiment ``name``."""
    return results_dir / f"{name}.{timestamp}.meta.json"


def require_json_destination(path: Path) -> None:
    """Validate that ``path`` names a JSON file.

    Args:
        path: Candidate result or metadata destination.

    Raises:
        ValueError: If ``path`` does not end in :data:`RESULT_SUFFIX`.
    """
    if path.suffix.lower() != RESULT_SUFFIX:
        msg = (
            f"output destination must be a {RESULT_SUFFIX!r} file, got {path.name!r}; "
            "experiment results and metadata are written only as JSON"
        )
        raise ValueError(msg)


def write_result(data: pd.DataFrame, path: Path) -> Path:
    """Write a result frame to ``path`` as JSON, creating parent directories.

    The frame is serialized with the ``records`` orientation.

    Args:
        data: The result frame to serialize.
        path: Destination path; must be a ``.json`` file. Parent directories are
            created if missing.

    Returns:
        The path written to.

    Raises:
        ValueError: If ``path`` is not a ``.json`` file.
    """
    require_json_destination(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_json(path, orient="records", indent=2)
    return path


def write_metadata(metadata: Mapping[str, object], path: Path) -> Path:
    """Write run metadata to ``path`` as pretty-printed UTF-8 JSON.

    Keys are sorted for reproducible output, and values without a native JSON
    encoding (paths, timestamps) are rendered via :func:`str`. Parent
    directories are created if missing.

    Args:
        metadata: Mapping of metadata fields to serialize.
        path: Destination path; must be a ``.json`` file. Parent directories are
            created if missing.

    Returns:
        The path written to.

    Raises:
        ValueError: If ``path`` is not a ``.json`` file.
    """
    require_json_destination(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    return path
