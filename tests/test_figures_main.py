"""Tests for the figure launcher."""

import sys
import types

import click
import pytest

from figures import FigureContext, figure
from figures.__main__ import main


def _register(
    monkeypatch: pytest.MonkeyPatch, name: str, command: click.Command
) -> None:
    """Install ``command`` as the sole command of a fake ``figures.<name>``."""
    module = types.ModuleType(f"figures.{name}")
    setattr(module, "command", command)  # noqa: B010 - module has no static attr
    monkeypatch.setitem(sys.modules, f"figures.{name}", module)


def test_main_returns_usage_code_when_no_figure_named():
    assert main([]) == 2


def test_main_reports_unknown_figure_with_its_name(capsys):
    code = main(["ghost", "--dpi", "72"])

    assert code == 2
    assert "ghost" in capsys.readouterr().err


def test_main_dispatches_to_the_named_figure(monkeypatch, tmp_path):
    seen: dict[str, object] = {}

    @figure(name="probe", figures_dir=tmp_path)
    @click.option("--n", type=int, default=1)
    def run(ctx: FigureContext, n: int) -> None:
        seen["n"] = n

    _register(monkeypatch, "probe", run)

    code = main(["probe", "--n", "5"])

    assert code == 0
    assert seen["n"] == 5


def test_main_resolves_hyphenated_names_to_underscored_modules(monkeypatch, tmp_path):
    seen: dict[str, bool] = {}

    @figure(name="my-sheet", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        seen["ran"] = True

    _register(monkeypatch, "my_sheet", run)

    assert main(["my-sheet"]) == 0
    assert seen["ran"] is True


def test_main_surfaces_a_nonzero_exit_from_the_figure(monkeypatch, tmp_path):
    @figure(name="boom", figures_dir=tmp_path)
    def run(ctx: FigureContext) -> None:
        raise SystemExit(3)

    _register(monkeypatch, "boom", run)

    assert main(["boom"]) == 3


def test_main_surfaces_a_usage_error_from_the_figure(monkeypatch, tmp_path):
    @figure(name="typed", figures_dir=tmp_path)
    @click.option("--n", type=int, required=True)
    def run(ctx: FigureContext, n: int) -> None:
        pass

    _register(monkeypatch, "typed", run)

    assert main(["typed", "--n", "not-an-int"]) == 2


def test_main_reports_a_module_that_defines_no_command(monkeypatch, capsys):
    blank = types.ModuleType("figures.blank")
    monkeypatch.setitem(sys.modules, "figures.blank", blank)

    code = main(["blank"])

    assert code == 2
    assert "no click command" in capsys.readouterr().err
