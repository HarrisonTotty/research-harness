"""Tests for the experiment launcher."""

import sys
import types

import click
import pandas as pd
import pytest

from experiments import ExperimentContext, experiment
from experiments.__main__ import main


def _register(
    monkeypatch: pytest.MonkeyPatch, name: str, command: click.Command
) -> None:
    """Install ``command`` as the sole command of a fake ``experiments.<name>``."""
    module = types.ModuleType(f"experiments.{name}")
    setattr(module, "command", command)  # noqa: B010 - module has no static attr
    monkeypatch.setitem(sys.modules, f"experiments.{name}", module)


def test_main_returns_usage_code_when_no_experiment_named():
    assert main([]) == 2


def test_main_reports_unknown_experiment_with_its_name(capsys):
    code = main(["ghost", "--n", "10"])

    assert code == 2
    assert "ghost" in capsys.readouterr().err


def test_main_dispatches_to_the_named_experiment(monkeypatch, tmp_path):
    out = tmp_path / "probe.json"

    @experiment(name="probe", results_dir=tmp_path)
    @click.option("--n", type=int, default=1)
    def run(ctx: ExperimentContext, n: int) -> None:
        ctx.write_result(pd.DataFrame({"n": [n]}))

    _register(monkeypatch, "probe", run)

    code = main(["probe", "--n", "5", "--out", str(out)])

    assert code == 0
    assert pd.read_json(out)["n"].tolist() == [5]


def test_main_resolves_hyphenated_names_to_underscored_modules(monkeypatch):
    seen: dict[str, bool] = {}

    @experiment(name="my-run")
    def run(ctx: ExperimentContext) -> None:
        seen["ran"] = True

    _register(monkeypatch, "my_run", run)

    assert main(["my-run"]) == 0
    assert seen["ran"] is True


def test_main_surfaces_a_nonzero_exit_from_the_experiment(monkeypatch):
    @experiment(name="boom")
    def run(ctx: ExperimentContext) -> None:
        raise SystemExit(3)

    _register(monkeypatch, "boom", run)

    assert main(["boom"]) == 3


def test_main_surfaces_a_usage_error_from_the_experiment(monkeypatch):
    @experiment(name="typed")
    @click.option("--n", type=int, required=True)
    def run(ctx: ExperimentContext, n: int) -> None:
        pass

    _register(monkeypatch, "typed", run)

    assert main(["typed", "--n", "not-an-int"]) == 2


def test_main_reports_a_module_that_defines_no_command(monkeypatch, capsys):
    blank = types.ModuleType("experiments.blank")
    monkeypatch.setitem(sys.modules, "experiments.blank", blank)

    code = main(["blank"])

    assert code == 2
    assert "no click command" in capsys.readouterr().err
