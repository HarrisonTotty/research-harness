"""Tests for experiment result and metadata I/O."""

import json
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from experiments.io import (
    default_metadata_path,
    default_result_path,
    require_json_destination,
    write_metadata,
    write_result,
)


@pytest.fixture
def frame():
    return pd.DataFrame({"n": [1, 2, 3], "label": ["a", "b", "c"]})


def test_default_result_path_stamps_the_name_with_the_timestamp():
    path = default_result_path("sweep", "20260724T120000Z")
    assert path == Path("data/results/sweep.20260724T120000Z.json")


def test_default_metadata_path_is_a_timestamped_json_sidecar():
    path = default_metadata_path("sweep", "20260724T120000Z")
    assert path == Path("data/results/sweep.20260724T120000Z.meta.json")


@pytest.mark.parametrize("name", ["r.csv", "r.parquet", "r.txt", "r"])
def test_require_json_destination_rejects_non_json(name):
    with pytest.raises(ValueError, match="written only as JSON"):
        require_json_destination(Path(name))


def test_write_result_roundtrips_json(frame, tmp_path):
    path = write_result(frame, tmp_path / "out" / "r.json")
    assert_frame_equal(pd.read_json(path), frame)


def test_write_result_rejects_a_non_json_destination(frame, tmp_path):
    with pytest.raises(ValueError, match="written only as JSON"):
        write_result(frame, tmp_path / "r.csv")


def test_write_metadata_writes_sorted_pretty_json(tmp_path):
    path = write_metadata({"b": 2, "a": 1}, tmp_path / "meta.json")
    text = path.read_text(encoding="utf-8")
    assert json.loads(text) == {"a": 1, "b": 2}
    assert text.index('"a"') < text.index('"b"')
    assert text.endswith("\n")


def test_write_metadata_stringifies_non_json_values(tmp_path):
    path = write_metadata({"where": Path("data/results")}, tmp_path / "meta.json")
    assert json.loads(path.read_text(encoding="utf-8")) == {"where": "data/results"}


def test_write_metadata_rejects_a_non_json_destination(tmp_path):
    with pytest.raises(ValueError, match="written only as JSON"):
        write_metadata({"a": 1}, tmp_path / "meta.yaml")
