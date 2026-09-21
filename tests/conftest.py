"""Shared fixtures.

Removed LNG_* flags are rejected at runtime (see `_REMOVED_ENV`), so a stale
value exported in a developer's shell or a CI runner would turn every test that
builds a network into a spurious failure. Tests should exercise the code, not
the ambient environment, so clear those names for the whole suite. The guard
itself is tested explicitly in tests/test_removed_env_flags.py.
"""

import pytest

from src.logic_network_generator import _REMOVED_ENV


@pytest.fixture(autouse=True)
def _clear_removed_env_flags(monkeypatch):
    for name in _REMOVED_ENV:
        monkeypatch.delenv(name, raising=False)
    yield
