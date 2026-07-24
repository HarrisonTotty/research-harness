"""Tests for experiment logging configuration."""

import io
import logging

import pytest

from experiments.logging import LogLevel, configure_logging


@pytest.mark.parametrize(
    ("level", "expected"),
    [
        (LogLevel.DEBUG, logging.DEBUG),
        (LogLevel.INFO, logging.INFO),
        (LogLevel.WARNING, logging.WARNING),
        (LogLevel.ERROR, logging.ERROR),
        (LogLevel.CRITICAL, logging.CRITICAL),
    ],
)
def test_numeric_maps_each_level_to_its_logging_threshold(level, expected):
    assert level.numeric == expected


def test_configure_logging_sets_the_requested_root_level():
    root = configure_logging(level=LogLevel.WARNING, stream=io.StringIO())
    assert root.level == logging.WARNING


def test_configure_logging_emits_records_to_the_given_stream():
    stream = io.StringIO()
    root = configure_logging(level=LogLevel.INFO, stream=stream)
    root.info("hello %s", "world")
    assert "hello world" in stream.getvalue()


def test_configure_logging_does_not_duplicate_its_handlers_across_calls():
    before = len(logging.getLogger().handlers)
    configure_logging(stream=io.StringIO())
    configure_logging(stream=io.StringIO())
    after = len(logging.getLogger().handlers)
    assert after == before + 1


def test_configure_logging_mirrors_records_to_the_log_file(tmp_path):
    log_file = tmp_path / "nested" / "run.log"
    root = configure_logging(
        level=LogLevel.INFO, log_file=log_file, stream=io.StringIO()
    )
    root.warning("to file")
    for handler in root.handlers:
        handler.flush()
    assert "to file" in log_file.read_text(encoding="utf-8")
