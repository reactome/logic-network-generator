"""The removed-flag guard, tested for real.

This file exists because the first version of the guard was tested only by
asserting that a string appeared in `inspect.getsource(...)`. That passes
whether or not the guard runs: mutating it to `if False and ...` left the whole
suite green. These tests call the guard and the entry point that uses it.
"""
import subprocess
import sys
from pathlib import Path

import pytest

from src.logic_network_generator import _REMOVED_ENV, _reject_removed_env

REPO = Path(__file__).resolve().parents[1]


def test_clean_environment_is_fine(monkeypatch):
    for name in _REMOVED_ENV:
        monkeypatch.delenv(name, raising=False)
    _reject_removed_env()          # must not raise


@pytest.mark.parametrize("name", sorted(_REMOVED_ENV))
@pytest.mark.parametrize("value", ["1", "0", "any", "", "downstream_free"])
def test_any_value_of_a_removed_flag_raises(monkeypatch, name, value):
    # Including "0" and "": the point is that the name is gone, so there is no
    # value that means "use the old behaviour" and none that means "off".
    monkeypatch.setenv(name, value)
    with pytest.raises(ValueError, match=f"{name} was removed"):
        _reject_removed_env()


@pytest.mark.parametrize("name", sorted(_REMOVED_ENV))
def test_the_message_names_the_spec_that_holds_the_evidence(monkeypatch, name):
    monkeypatch.setenv(name, "1")
    with pytest.raises(ValueError) as exc:
        _reject_removed_env()
    assert "specs/" in str(exc.value)


@pytest.mark.parametrize("name", sorted(_REMOVED_ENV))
def test_entry_point_exits_nonzero_before_doing_any_work(monkeypatch, name, tmp_path):
    # The generator's per-pathway guard is not enough on its own: the pathway
    # loop catches every exception and continues, so a stale flag would leave
    # the previous run's files on disk and still exit 0.
    env = {**dict(__import__("os").environ), name: "1"}
    out = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, "bin/create-pathways.py", "--pathway-id", "R-HSA-73894",
         "--output-dir", str(out)],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=180,
    )
    assert proc.returncode != 0, "a removed flag must stop the run"
    assert f"{name} was removed" in (proc.stderr + proc.stdout)
    assert not out.exists() or not any(out.iterdir()), "must fail before writing output"


def test_a_failing_pathway_makes_the_run_exit_nonzero(tmp_path):
    """A partial catalog must not look like a successful run.

    The pathway loop catches every exception and continues, so before this the
    command logged "N failed" and still exited 0. A benchmark that globs the
    output directory would then score whatever mix of new and stale pathways
    happened to be there. Any failure reason serves here -- a bogus stable id
    or an unreachable database both take the same path.
    """
    listing = tmp_path / "list.tsv"
    listing.write_text("id\tpathway_name\nR-HSA-NOT-A-REAL-PATHWAY\tbogus\n")
    out = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, "bin/create-pathways.py",
         "--pathway-list", str(listing), "--output-dir", str(out)],
        cwd=REPO, capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode != 0, (
        "a run in which every pathway failed must not exit 0\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
