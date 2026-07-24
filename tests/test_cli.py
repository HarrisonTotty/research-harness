"""Tests for the experiment command decorator and its context."""

import json
import logging
from datetime import UTC, datetime

import click
import pandas as pd
import pytest
from click.testing import CliRunner

from experiments.cli import ExperimentContext, experiment


@pytest.fixture
def runner():
    return CliRunner()


def _frame() -> pd.DataFrame:
    return pd.DataFrame({"n": [1, 2]})


def _fixed_clock() -> datetime:
    return datetime(2026, 7, 24, 12, 0, 0, tzinfo=UTC)


_FIXED_STAMP = "20260724T120000Z"


def test_experiment_passes_a_context_and_declared_options(runner):
    seen: dict[str, object] = {}

    @experiment(name="demo")
    @click.option("--n", type=int, default=7)
    def run(ctx: ExperimentContext, n: int) -> None:
        seen["ctx"] = ctx
        seen["n"] = n

    result = runner.invoke(run, ["--n", "3"], catch_exceptions=False)

    assert result.exit_code == 0
    assert seen["n"] == 3
    assert isinstance(seen["ctx"], ExperimentContext)
    assert seen["ctx"].name == "demo"


def test_experiment_name_defaults_to_hyphenated_function_name(runner):
    seen: dict[str, object] = {}

    @experiment()
    def my_sweep(ctx: ExperimentContext) -> None:
        seen["name"] = ctx.name

    runner.invoke(my_sweep, [], catch_exceptions=False)

    assert seen["name"] == "my-sweep"


def test_experiment_produces_a_click_command():
    @experiment(name="demo")
    def run(ctx: ExperimentContext) -> None:
        pass

    assert isinstance(run, click.Command)
    assert run.name == "demo"


def test_context_exposes_the_run_timestamp(runner):
    seen: dict[str, object] = {}

    @experiment(name="demo", clock=_fixed_clock)
    def run(ctx: ExperimentContext) -> None:
        seen["timestamp"] = ctx.timestamp

    runner.invoke(run, [], catch_exceptions=False)

    assert seen["timestamp"] == _FIXED_STAMP


def test_write_result_defaults_to_a_timestamped_json(runner, tmp_path):
    @experiment(name="demo", results_dir=tmp_path, clock=_fixed_clock)
    def run(ctx: ExperimentContext) -> None:
        ctx.write_result(_frame())

    runner.invoke(run, [], catch_exceptions=False)

    written = tmp_path / f"demo.{_FIXED_STAMP}.json"
    assert written.exists()
    pd.testing.assert_frame_equal(pd.read_json(written), _frame())


def test_metadata_defaults_to_a_timestamped_sidecar(runner, tmp_path):
    @experiment(name="demo", results_dir=tmp_path, clock=_fixed_clock)
    def run(ctx: ExperimentContext) -> None:
        ctx.write_metadata({"params": {"n": 2}})

    runner.invoke(run, [], catch_exceptions=False)

    sidecar = tmp_path / f"demo.{_FIXED_STAMP}.meta.json"
    assert json.loads(sidecar.read_text(encoding="utf-8")) == {"params": {"n": 2}}


def test_out_option_overrides_the_result_path(runner, tmp_path):
    @experiment(name="demo", results_dir=tmp_path)
    def run(ctx: ExperimentContext) -> None:
        ctx.write_result(_frame())

    out = tmp_path / "custom.json"
    runner.invoke(run, ["--out", str(out)], catch_exceptions=False)

    assert out.exists()
    pd.testing.assert_frame_equal(pd.read_json(out), _frame())


def test_out_option_rejects_a_non_json_path(runner, tmp_path):
    @experiment(name="demo", results_dir=tmp_path)
    def run(ctx: ExperimentContext) -> None:
        ctx.write_result(_frame())

    result = runner.invoke(run, ["--out", str(tmp_path / "custom.csv")])

    assert result.exit_code == 2
    assert "JSON" in result.output


def test_meta_out_option_overrides_the_metadata_path(runner, tmp_path):
    @experiment(name="demo", results_dir=tmp_path)
    def run(ctx: ExperimentContext) -> None:
        ctx.write_metadata({"params": {"n": 2}})

    meta = tmp_path / "custom.meta.json"
    runner.invoke(run, ["--meta-out", str(meta)], catch_exceptions=False)

    assert json.loads(meta.read_text(encoding="utf-8")) == {"params": {"n": 2}}


def test_meta_out_option_rejects_a_non_json_path(runner, tmp_path):
    @experiment(name="demo")
    def run(ctx: ExperimentContext) -> None:
        pass

    result = runner.invoke(run, ["--meta-out", str(tmp_path / "meta.yaml")])

    assert result.exit_code == 2
    assert "JSON" in result.output


def test_log_level_option_sets_the_root_threshold(runner):
    @experiment(name="demo")
    def run(ctx: ExperimentContext) -> None:
        pass

    runner.invoke(run, ["--log-level", "debug"], catch_exceptions=False)

    assert logging.getLogger().level == logging.DEBUG


def test_invalid_log_level_is_rejected(runner):
    @experiment(name="demo")
    def run(ctx: ExperimentContext) -> None:
        pass

    result = runner.invoke(run, ["--log-level", "verbose"])

    assert result.exit_code == 2
    assert "verbose" in result.output
