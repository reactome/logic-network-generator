"""A failed Neo4j lookup must fail the build, not be replaced by an empty or
fallback answer that changes the network (code review 2026-10-02)."""
import pytest

import src.logic_network_generator as lng
import src.neo4j_connector as nc
import src.reaction_generator as rg


def test_modifier_set_failure_raises_and_is_not_cached(monkeypatch):
    monkeypatch.setattr(rg, "_modifier_set_cache", None)

    def down():
        raise ConnectionError("neo4j down")
    monkeypatch.setattr(nc, "get_modifier_isoform_entity_set_ids", down)
    with pytest.raises(ConnectionError):
        rg.modifier_isoform_set_ids()
    assert rg._modifier_set_cache is None          # the failure did not stick
    monkeypatch.setattr(nc, "get_modifier_isoform_entity_set_ids", lambda: {"R-HSA-X"})
    assert "R-HSA-X" in rg.modifier_isoform_set_ids()


def test_containment_lookup_failure_raises(monkeypatch, tmp_path):
    # export_containment imports these inside the function
    monkeypatch.setattr(nc, "get_reactome_release", lambda: 97)

    def down(stid):
        raise ConnectionError("neo4j down")
    monkeypatch.setattr(rg, "get_terminal_components", down)
    with pytest.raises(ConnectionError):
        lng.export_containment({"u1": "R-HSA-1"}, str(tmp_path / "containment.csv"))
    assert not (tmp_path / "containment.csv").exists()
