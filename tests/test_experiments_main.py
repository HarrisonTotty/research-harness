"""Tests for the experiment launcher."""

import pytest

from experiments.__main__ import main


def test_main_returns_usage_code_when_no_experiment_named():
    assert main([]) == 2


def test_main_reports_unknown_experiment_with_its_name():
    with pytest.raises(NotImplementedError, match="'ghost'"):
        main(["ghost", "--n", "10"])
