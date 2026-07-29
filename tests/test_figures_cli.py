"""Tests for the figure command decorator and its context."""

import logging
from datetime import UTC, datetime

import click
import matplotlib
import pytest
from click.testing import CliRunner
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from figures.cli import FigureContext, FigureFormat, figure

matplotlib.use("Agg")


@pytest.fixture
def runner():
    return CliRunner()


def _fixed_clock() -> datetime:
    return datetime(2026, 7, 24, 12, 0, 0, tzinfo=UTC)


_FIXED_STAMP = "20260724T120000Z"


def _blank_figure() -> Figure:
    fig = plt.figure(figsize=(1, 1))
    fig.add_subplot().plot([0, 1], [0, 1])
    return fig


def test_figure_passes_a_context_and_declared_options(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path)
    @click.option("--n", type=int, default=7)
    def run(ctx: FigureContext, n: int) -> None:
        seen["ctx"] = ctx
        seen["n"] = n

    result = runner.invoke(run, ["--n", "3"], catch_exceptions=False)

    assert result.exit_code == 0
    assert seen["n"] == 3
    assert isinstance(seen["ctx"], FigureContext)
    assert seen["ctx"].name == "demo"


def test_figure_name_defaults_to_hyphenated_function_name(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(figures_dir=tmp_path)
    def my_sheet(ctx: FigureContext) -> None:
        seen["name"] = ctx.name

    runner.invoke(my_sheet, [], catch_exceptions=False)

    assert seen["name"] == "my-sheet"


def test_figure_produces_a_click_command(tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        pass

    assert isinstance(run, click.Command)
    assert run.name == "demo"


def test_context_exposes_the_run_timestamp(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path, clock=_fixed_clock)
    def run(ctx: FigureContext) -> None:
        seen["timestamp"] = ctx.timestamp

    runner.invoke(run, [], catch_exceptions=False)

    assert seen["timestamp"] == _FIXED_STAMP


def test_output_dir_defaults_to_the_figure_name(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        seen["output_dir"] = ctx.output_dir

    runner.invoke(run, [], catch_exceptions=False)

    assert seen["output_dir"] == tmp_path / "demo"


def test_save_writes_one_file_per_default_format(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        fig = _blank_figure()
        ctx.save(fig, "sheet")
        plt.close(fig)

    runner.invoke(run, [], catch_exceptions=False)

    assert (tmp_path / "demo" / "sheet.png").exists()
    assert (tmp_path / "demo" / "sheet.svg").exists()


def test_format_option_selects_the_output_formats(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        fig = _blank_figure()
        ctx.save(fig, "sheet")
        plt.close(fig)

    runner.invoke(run, ["--format", "pdf"], catch_exceptions=False)

    assert (tmp_path / "demo" / "sheet.pdf").exists()
    assert not (tmp_path / "demo" / "sheet.png").exists()


def test_invalid_format_is_rejected(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        pass

    result = runner.invoke(run, ["--format", "webp"])

    assert result.exit_code == 2
    assert "webp" in result.output


def test_out_dir_option_overrides_the_output_directory(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path / "unused")
    def run(ctx: FigureContext) -> None:
        fig = _blank_figure()
        ctx.save(fig, "sheet")
        plt.close(fig)

    override = tmp_path / "elsewhere"
    runner.invoke(run, ["--out-dir", str(override)], catch_exceptions=False)

    assert (override / "sheet.png").exists()
    assert not (tmp_path / "unused").exists()


def test_save_returns_the_written_paths_in_format_order(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        fig = _blank_figure()
        seen["paths"] = ctx.save(fig, "sheet")
        plt.close(fig)

    runner.invoke(run, ["-f", "svg", "-f", "png"], catch_exceptions=False)

    assert seen["paths"] == [
        tmp_path / "demo" / "sheet.svg",
        tmp_path / "demo" / "sheet.png",
    ]


def test_body_runs_under_the_blog_style(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        seen["grid.color"] = matplotlib.rcParams["grid.color"]

    runner.invoke(run, [], catch_exceptions=False)

    assert seen["grid.color"] == "#ddd8cf"


def test_context_formats_are_typed_figure_formats(runner, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        seen["formats"] = ctx.formats

    runner.invoke(run, ["-f", "pdf", "-f", "svg"], catch_exceptions=False)

    assert seen["formats"] == (FigureFormat.PDF, FigureFormat.SVG)


def test_log_level_option_sets_the_root_threshold(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        pass

    runner.invoke(run, ["--log-level", "debug"], catch_exceptions=False)

    assert logging.getLogger().level == logging.DEBUG


def test_invalid_log_level_is_rejected(runner, tmp_path):
    @figure(name="demo", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        pass

    result = runner.invoke(run, ["--log-level", "verbose"])

    assert result.exit_code == 2
    assert "verbose" in result.output
