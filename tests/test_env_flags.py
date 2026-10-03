"""LNG_* switches take exactly 0 or 1, and an unknown LNG_* name is an error
(code review 2026-10-02: `!= "0"` switches turned on for "false", `== "1"`
switches turned off for "true", and a misspelt name was never read)."""
import pytest

from src.env_flags import BOOL_FLAGS, env_flag, validate_env
from src.logic_network_generator import _REMOVED_ENV, _reject_removed_env


@pytest.mark.parametrize("name", sorted(BOOL_FLAGS))
def test_switch_values(monkeypatch, name):
    monkeypatch.delenv(name, raising=False)
    assert env_flag(name) is BOOL_FLAGS[name]
    monkeypatch.setenv(name, "")
    assert env_flag(name) is BOOL_FLAGS[name]
    monkeypatch.setenv(name, "1")
    assert env_flag(name) is True
    monkeypatch.setenv(name, "0")
    assert env_flag(name) is False
    for bad in ("true", "false", "yes", "2", "on"):
        monkeypatch.setenv(name, bad)
        with pytest.raises(ValueError, match=name):
            env_flag(name)


def test_the_two_reported_inversions_are_errors(monkeypatch):
    monkeypatch.setenv("LNG_HANDOFF_EDGES", "false")
    with pytest.raises(ValueError):
        _reject_removed_env()
    monkeypatch.delenv("LNG_HANDOFF_EDGES")
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "true")
    with pytest.raises(ValueError):
        _reject_removed_env()


def test_a_misspelt_name_is_an_error(monkeypatch):
    monkeypatch.setenv("LNG_BOUNDRY_HIERARCHY", "0")
    with pytest.raises(ValueError, match="LNG_BOUNDRY_HIERARCHY"):
        _reject_removed_env()


def test_known_names_pass(monkeypatch):
    for n in BOOL_FLAGS:
        monkeypatch.setenv(n, "0")
    monkeypatch.setenv("LNG_MAX_VARIANTS", "512")
    monkeypatch.setenv("LNG_PYTHON", "/usr/bin/python3")
    validate_env(removed=_REMOVED_ENV)
