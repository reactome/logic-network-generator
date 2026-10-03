"""Edges derived from complex/set structure at root inputs and terminal outputs
go to boundary_edges.csv; logic_network.csv holds only what was curated
(deltasignal specs/044)."""
import pandas as pd

import src.logic_network_generator as lng
import src.pathway_generator as pg
from src.logic_network_generator import PathwayResult
from tests.test_export_failfast import mocked  # noqa: F401  (fixture)


def test_the_boundary_pass_tags_exactly_what_it_adds(monkeypatch):
    data = [{"source_id": "r", "target_id": "x", "edge_type": "input"}]

    def inner(d, _m):
        d.append({"source_id": "a", "target_id": "X", "edge_type": "assembly"})
        d.append({"source_id": "Y", "target_id": "s", "edge_type": "dissociation"})
    monkeypatch.setattr(lng, "_emit_boundary_decomposition_edges_inner", inner)
    lng._emit_boundary_decomposition_edges(data, {})
    assert [bool(e.get("_boundary")) for e in data] == [False, True, True]


def test_the_export_splits_by_the_mask(tmp_path, mocked):  # noqa: F811
    net = pd.DataFrame([
        {"source_id": "x", "target_id": "y", "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-1"},
        {"source_id": "a", "target_id": "x", "pos_neg": "pos", "and_or": "and",
         "edge_type": "assembly", "stoichiometry": 1, "edge_reaction_id": None},
    ])
    mocked.setattr(pg, "create_pathway_logic_network", lambda *a, **k: PathwayResult(
        net, {"x": "R-HSA-2"}, pd.DataFrame(), pd.DataFrame({"uid": [], "reactome_id": []}), {},
        boundary_mask=[False, True]))
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    d = tmp_path / "R-HSA-1"
    curated = pd.read_csv(d / "logic_network.csv")
    boundary = pd.read_csv(d / "boundary_edges.csv")
    assert list(curated["edge_type"]) == ["input"]
    assert list(boundary["edge_type"]) == ["assembly"]
    assert list(curated.columns) == list(boundary.columns)


def test_a_mask_of_the_wrong_length_fails(tmp_path, mocked):  # noqa: F811
    net = pd.DataFrame([{"source_id": "x", "target_id": "y", "pos_neg": "pos", "and_or": "and",
                         "edge_type": "input", "stoichiometry": 1}])
    mocked.setattr(pg, "create_pathway_logic_network", lambda *a, **k: PathwayResult(
        net, {"x": "R-HSA-2"}, pd.DataFrame(), pd.DataFrame({"uid": [], "reactome_id": []}), {},
        boundary_mask=[False, True]))
    try:
        pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    except Exception as e:
        assert "boundary mask" in str(e)
    else:
        raise AssertionError("a misaligned mask must fail the pathway")
