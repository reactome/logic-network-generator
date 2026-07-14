"""Unit tests for diagram-sourced reaction connectivity (no Neo4j required).

Builds a tiny synthetic diagram (two reactions sharing one entity glyph) and
checks that the shared-glyph product->substrate pair is extracted and unioned
into reaction_connections. See reactome/logic-network-generator#39.
"""
import json

import pandas as pd
import pytest

import src.diagram_connectivity as dc


def _write_diagram(dirpath, stid):
    """Two reactions: A (R-HSA-100) outputs glyph 1; B (R-HSA-200) inputs glyph 1."""
    layout = {
        "nodes": [{"id": 1, "reactomeId": 500, "displayName": "X"}],
        "edges": [
            {"id": 10, "reactomeId": 100, "inputs": [], "outputs": [{"id": 1}], "catalysts": []},
            {"id": 20, "reactomeId": 200, "inputs": [{"id": 1}], "outputs": [], "catalysts": []},
        ],
    }
    graph = {"edges": [{"dbId": 100, "stId": "R-HSA-100"},
                       {"dbId": 200, "stId": "R-HSA-200"}]}
    (dirpath / f"{stid}.json").write_text(json.dumps(layout))
    (dirpath / f"{stid}.graph.json").write_text(json.dumps(graph))


@pytest.fixture
def diagram_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("LNG_DIAGRAM_DIR", str(tmp_path))
    # isolate this pathway to its two reactions (avoid Neo4j)
    monkeypatch.setattr(dc, "_pathway_reaction_stids",
                        lambda pid: {"R-HSA-100", "R-HSA-200"})
    return tmp_path


def test_shared_glyph_pair_extracted(diagram_dir):
    _write_diagram(diagram_dir, "R-HSA-TEST")
    pairs = dc.diagram_shared_product_pairs("R-HSA-TEST")
    assert pairs == {("R-HSA-100", "R-HSA-200")}


def test_no_diagram_returns_empty(diagram_dir, monkeypatch):
    # no diagram file written; ancestor lookup also finds nothing
    monkeypatch.setattr(dc, "_covering_diagram_stid", lambda pid: "")
    assert dc.diagram_shared_product_pairs("R-HSA-MISSING") == set()


def test_augment_adds_missing_pair(diagram_dir):
    _write_diagram(diagram_dir, "R-HSA-TEST")
    rc = pd.DataFrame({"preceding_reaction_id": ["R-HSA-100"],
                       "following_reaction_id": ["R-HSA-999"],
                       "event_status": ["Has Preceding Event"]})
    out = dc.augment_reaction_connections("R-HSA-TEST", rc)
    added = out[out["event_status"] == "Diagram Shared Product"]
    assert list(zip(added["preceding_reaction_id"], added["following_reaction_id"])) \
        == [("R-HSA-100", "R-HSA-200")]


def test_augment_skips_when_already_present(diagram_dir):
    _write_diagram(diagram_dir, "R-HSA-TEST")
    rc = pd.DataFrame({"preceding_reaction_id": ["R-HSA-100"],
                       "following_reaction_id": ["R-HSA-200"],
                       "event_status": ["Has Preceding Event"]})
    out = dc.augment_reaction_connections("R-HSA-TEST", rc)
    assert len(out) == 1  # nothing added; pair already linked


def test_augment_disabled_by_env(diagram_dir, monkeypatch):
    _write_diagram(diagram_dir, "R-HSA-TEST")
    monkeypatch.setenv("LNG_DIAGRAM_CONNECTIVITY", "0")
    rc = pd.DataFrame({"preceding_reaction_id": ["R-HSA-1"],
                       "following_reaction_id": ["R-HSA-2"],
                       "event_status": ["x"]})
    assert dc.augment_reaction_connections("R-HSA-TEST", rc) is rc
