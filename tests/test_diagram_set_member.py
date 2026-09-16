"""The diagram's set-member links, and the derived cofactor set.

Both cover defects found by auditing the generated networks against Reactome:
a realisation relationship curators drew that the networks did not carry, and a
hand-written cofactor list with six of thirteen entries stale or mislabelled.
No Neo4j: the derivation is stubbed.
"""
import json

import src.logic_network_generator as m
from src import diagram_connectivity as dc
from src import neo4j_connector


def _layout(links):
    return {"nodes": [{"id": 1, "reactomeId": 100}, {"id": 2, "reactomeId": 200},
                      {"id": 3, "reactomeId": 300}],
            "links": links}


def _graph():
    return {"nodes": [{"dbId": 100, "stId": "R-HSA-100"},
                      {"dbId": 200, "stId": "R-HSA-200"},
                      {"dbId": 300, "stId": "R-HSA-300"}]}


def test_only_set_member_links_are_taken(tmp_path, monkeypatch):
    """Interaction and set-to-set links point at entities with no node here;
    consuming them means ADDING nodes, which is a different change (#41)."""
    links = [
        {"renderableClass": "EntitySetAndMemberLink",
         "inputs": [{"id": 1}], "outputs": [{"id": 2}]},
        {"renderableClass": "Interaction",
         "inputs": [{"id": 1}], "outputs": [{"id": 3}]},
        {"renderableClass": "EntitySetAndEntitySetLink",
         "inputs": [{"id": 2}], "outputs": [{"id": 3}]},
    ]
    (tmp_path / "R-HSA-1.json").write_text(json.dumps(_layout(links)))
    (tmp_path / "R-HSA-1.graph.json").write_text(json.dumps(_graph()))
    monkeypatch.setattr(dc, "_diagram_dir", lambda: tmp_path)
    monkeypatch.setattr(dc, "_covering_diagram_stid", lambda p: "R-HSA-1")

    # member first: the specific realises the generic.
    assert dc.diagram_set_member_pairs("R-HSA-1") == {("R-HSA-100", "R-HSA-200")}


def test_one_edge_per_pair_not_the_cartesian_product():
    """Positional decomposition gives an entity many uuids. All-pairs turned
    two curated relationships into 48 edges in EPH-Ephrin, which is the blow-up
    that made the all-pairs silo bridge unusable."""
    mapping = {"m1": "R-HSA-100", "m2": "R-HSA-100", "m3": "R-HSA-100",
               "s1": "R-HSA-200", "s2": "R-HSA-200"}
    edges = [{"source_id": "x", "target_id": "m2", "pos_neg": "pos",
              "and_or": "and", "edge_type": "input", "stoichiometry": 1},
             {"source_id": "s1", "target_id": "y", "pos_neg": "pos",
              "and_or": "and", "edge_type": "output", "stoichiometry": 1}]
    n = m._emit_diagram_set_member_edges(edges, mapping,
                                         {("R-HSA-100", "R-HSA-200")})
    assert n == 1, "3 members x 2 sets must not become 6 edges"
    new = [e for e in edges if e["edge_type"] == "diagram_set_member"][0]
    # the best-connected occurrence on each side
    assert new["source_id"] == "m2" and new["target_id"] == "s1"
    # never imposes AND-completeness on the generic, and asserts realisation
    assert new["and_or"] == "or" and new["pos_neg"] == "pos"


def test_no_edge_when_an_endpoint_has_no_node():
    """A link whose other end participates in no curated reaction here would
    need a node invented for it."""
    edges: list = []
    n = m._emit_diagram_set_member_edges(
        edges, {"m1": "R-HSA-100"}, {("R-HSA-100", "R-HSA-999")})
    assert n == 0 and edges == []


def test_cofactor_set_is_derived_and_seed_is_a_fallback(monkeypatch):
    """The list this replaced had 6 of 13 entries stale or mislabelled."""
    m._cofactor_stids_cache = None
    monkeypatch.setattr(neo4j_connector, "get_cofactor_species",
                        lambda: [{"stable_id": "R-ALL-1"}, {"stable_id": "R-ALL-2"}])
    assert m._cofactor_stids() == frozenset({"R-ALL-1", "R-ALL-2"})

    # Unreachable Neo4j must not silently yield an EMPTY set, which would stop
    # excluding cofactors from bridges and depletion edges everywhere.
    m._cofactor_stids_cache = None

    def boom():
        raise RuntimeError("no neo4j")

    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", boom)
    assert m._cofactor_stids() == m._COFACTOR_STIDS_SEED
    assert m._COFACTOR_STIDS_SEED, "the fallback must not be empty"
    m._cofactor_stids_cache = None


def test_the_six_bad_ids_are_not_in_the_seed():
    """Verified against Release97: three do not exist, three are mislabelled."""
    for bad in ("R-ALL-217093", "R-ALL-110114", "R-ALL-29986",   # absent
                "R-ALL-29390",   # PXLP, was commented "Pi variant"
                "R-ALL-29438",   # GTP, was commented "PPi"
                "R-ALL-29360"):  # NAD+, was commented "ADP variant"
        assert bad not in m._COFACTOR_STIDS_SEED, bad
