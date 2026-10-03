"""A pathway whose bundle cannot be written must FAIL, not count as built.

Code review 2026-10-02: every export after logic_network.csv was wrapped in a
try/except that logged and carried on. The orphan guard in
export_node_reaction_context raises ValueError on purpose, and that error
skipped node_resolution.csv, cofactors.csv and containment.csv, while
generate_pathway_file returned normally. The pathway counted as successful,
the build exited 0 and catalog.sh moved `current`; the consumer then ran its
default model with the self-inhibitor rules inert. A regeneration into an
existing directory also left the previous run's files behind.

These drive generate_pathway_file itself with the database and network
construction mocked, so the try/except blocks are what is under test.
"""
import pathlib

import pandas as pd
import pytest

import src.diagram_connectivity as dc
import src.pathway_generator as pg
from src.logic_network_generator import PathwayResult

EXPORTS = ["export_uuid_to_reactome_mapping", "export_entity_reaction_proxy_mapping",
           "export_nodes", "export_node_reaction_context", "export_node_resolution",
           "export_cofactors", "export_containment", "export_containment_structure",
           "export_drugs", "export_pools"]


@pytest.fixture
def mocked(monkeypatch):
    monkeypatch.setattr(pg, "get_reactome_release", lambda: 97)
    monkeypatch.setattr(pg, "get_reaction_connections", lambda pid: pd.DataFrame(
        {"preceding_reaction_id": ["R-HSA-1"], "following_reaction_id": [None],
         "event_status": ["No Preceding Event"]}))
    monkeypatch.setattr(pg, "get_decomposed_uid_mapping",
                        lambda pid, rc: (pd.DataFrame({"uid": ["u"]}), [("a", "b", "R-HSA-1")]))
    monkeypatch.setattr(dc, "augment_reaction_connections", lambda pid, rc: rc)
    net = pd.DataFrame([{"source_id": "x", "target_id": "y", "pos_neg": "pos", "and_or": "and",
                         "edge_type": "input", "stoichiometry": 1}])
    monkeypatch.setattr(pg, "create_pathway_logic_network", lambda *a, **k: PathwayResult(
        net, {"x": "R-HSA-2"}, pd.DataFrame(), pd.DataFrame({"uid": [], "reactome_id": []}), {}))

    def writes(*a, **k):
        for p in (x for x in a if isinstance(x, str) and x.endswith(".csv")):
            pathlib.Path(p).write_text("x\n")
    for name in EXPORTS:
        monkeypatch.setattr(pg, name, writes)
    return monkeypatch


@pytest.mark.parametrize("failing", EXPORTS)
def test_a_failed_export_fails_the_pathway(tmp_path, mocked, failing):
    def boom(*a, **k):
        raise ValueError("12 of 40 context rows name nodes absent from the logic network")
    mocked.setattr(pg, failing, boom)
    with pytest.raises(Exception):
        pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))


def test_a_complete_run_writes_every_bundle_file(tmp_path, mocked):
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    written = {p.name for p in (tmp_path / "R-HSA-1").glob("*.csv")}
    assert written == set(pg.BUNDLE_FILES)


def test_a_regeneration_leaves_no_stale_file(tmp_path, mocked):
    d = tmp_path / "R-HSA-1"
    d.mkdir()
    (d / "containment.csv").write_text("stale\n")

    def boom(*a, **k):
        raise ValueError("refusing")
    mocked.setattr(pg, "export_node_reaction_context", boom)
    with pytest.raises(Exception):
        pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    assert not (d / "containment.csv").exists()
