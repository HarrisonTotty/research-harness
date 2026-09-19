"""Tests for the strand-dynamics experiment (docs/experiments/strand-dynamics.md)."""

import json
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner, Result

from experiments import strand_dynamics as exp
from research.decorated_permutation import DecoratedPermutation

COLUMNS = [
    "n",
    "k",
    "permutation",
    "graph",
    "dimension",
    "internal_vertices",
    "edges",
    "noncrossing_matching",
    "converged",
    "transient",
    "period",
    "permutation_period",
    "cycle_permutations",
]


def _run(tmp_path: Path, *args: str) -> tuple[Result, Path, Path]:
    out, meta = tmp_path / "result.json", tmp_path / "result.meta.json"
    outcome = CliRunner().invoke(
        exp.strand_dynamics, [*args, "--out", str(out), "--meta-out", str(meta)]
    )
    return outcome, out, meta


def test_emits_one_row_per_cell_and_graph_with_the_designed_columns(tmp_path):
    outcome, out, _ = _run(tmp_path, "--n-max", "4")

    frame = pd.read_json(out, dtype={"permutation": str})
    assert outcome.exit_code == 0
    assert list(frame.columns) == COLUMNS
    assert len(frame) == 2 * (1 + 2 + 9)
    assert not frame.duplicated(["n", "permutation", "graph"]).any()
    assert frame.groupby("graph").size().to_dict() == {"bridge": 12, "le": 12}


def test_bridge_fixed_points_of_gr_2_4_are_the_two_noncrossing_matchings(tmp_path):
    _, out, _ = _run(tmp_path, "--n-min", "4", "--n-max", "4", "--graph", "bridge")

    frame = pd.read_json(out, dtype={"permutation": str})
    fixed = frame[(frame["period"] == 1) & (frame["transient"] == 0)]
    assert sorted(fixed["permutation"]) == ["2,1,4,3", "4,3,2,1"]
    assert sorted(frame[frame["noncrossing_matching"]]["permutation"]) == [
        "2,1,4,3",
        "4,3,2,1",
    ]
    assert set(fixed["k"]) == {2}


def test_an_exhausted_budget_is_recorded_as_an_unconverged_row_with_nulls(tmp_path):
    outcome, out, meta = _run(
        tmp_path,
        "--n-min",
        "3",
        "--n-max",
        "3",
        "--graph",
        "bridge",
        "--max-steps",
        "1",
    )

    frame = pd.read_json(out, dtype={"permutation": str})
    assert outcome.exit_code == 0
    assert frame["converged"].tolist() == [False, False]
    assert frame["period"].isna().all()
    assert json.loads(meta.read_text(encoding="utf-8"))["unconverged_rows"] == 2


def test_metadata_records_every_parameter_and_the_absence_of_seeds(tmp_path):
    _, _, meta = _run(tmp_path, "--n-max", "3")

    metadata = json.loads(meta.read_text(encoding="utf-8"))
    assert metadata["parameters"] == {
        "n_min": 2,
        "n_max": 3,
        "graphs": ["bridge", "le"],
        "max_steps": 100_000,
    }
    assert metadata["seeds"] == []
    assert metadata["rows"] == 6
    assert set(metadata["versions"]) == {"python", "pandas", "research-harness"}
    assert "git_commit" in metadata


def test_an_empty_range_is_a_usage_error(tmp_path):
    outcome, out, _ = _run(tmp_path, "--n-min", "5", "--n-max", "4")

    assert outcome.exit_code == 2
    assert "--n-min (5) must not exceed --n-max (4)" in outcome.output
    assert not out.exists()


@pytest.mark.parametrize(
    ("targets", "expected"),
    [
        ((2, 1), True),
        ((2, 1, 4, 3), True),
        ((4, 3, 2, 1), True),
        ((3, 4, 1, 2), False),
        ((2, 3, 1), False),
        ((2, 3, 4, 1), False),
    ],
)
def test_noncrossing_matchings_are_the_involutions_without_crossing_arcs(
    targets, expected
):
    # exercised through the module's private helper: it defines a schema column
    decorated = DecoratedPermutation(targets)
    assert exp._is_noncrossing_matching(decorated) is expected
