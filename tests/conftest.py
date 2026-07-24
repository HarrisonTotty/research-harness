"""Shared fixtures for the experiment harness tests."""

import logging

import pytest


@pytest.fixture(autouse=True)
def _restore_root_logger():
    """Snapshot the root logger and restore it after each test.

    Handlers installed during the test (in particular file handlers) are closed
    before the snapshot is restored, so an open file is never left for the
    garbage collector to finalize — which pytest would surface as an error.
    """
    root = logging.getLogger()
    saved_handlers = root.handlers[:]
    saved_level = root.level
    yield
    for handler in root.handlers:
        if handler not in saved_handlers:
            handler.close()
    root.handlers[:] = saved_handlers
    root.setLevel(saved_level)
