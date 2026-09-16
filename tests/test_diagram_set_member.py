"""The diagram's set-member links, and the derived cofactor set.

Both cover defects found by auditing the generated networks against Reactome:
a realisation relationship curators drew that the networks did not carry, and a
hand-written cofactor list with six of thirteen entries stale or mislabelled.
No Neo4j: the derivation is stubbed.
"""
import json

import pytest

import src.logic_network_generator as m
from src import diagram_connectivity as dc
from src import neo4j_connector


@pytest.fixture(autouse=True)
def _clear_module_caches():
    """Reset the caches these tests touch, on failure as well as success.

    An assertion failure used to leave a stub cofactor set latched in
    `_cofactor_stids_cache`, and `_handoff_leaf_cache` holds leaves computed
    under it, so every later test in the session saw the stub.
    """
    m._cofactor_stids_cache = None
    m._handoff_leaf_cache.clear()
    yield
    m._cofactor_stids_cache = None
    m._handoff_leaf_cache.clear()


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


def test_seed_drops_only_the_genuinely_wrong_entries():
    """A wrong COMMENT is not a wrong ENTRY, and conflating the two regressed
    this seed once already: GTP and NAD+ were dropped because their labels were
    wrong, which made the offline path exclude FEWER real cofactors than the
    list it replaced."""
    for absent in ("R-ALL-217093", "R-ALL-110114", "R-ALL-29986"):
        assert absent not in m._COFACTOR_STIDS_SEED, absent
    # PXLP is the one true false positive: not a cofactor at all.
    assert "R-ALL-29390" not in m._COFACTOR_STIDS_SEED
    # Mis-commented but genuine — both are in _COFACTOR_CHEBI.
    assert "R-ALL-29438" in m._COFACTOR_STIDS_SEED, "GTP is a cofactor"
    assert "R-ALL-29360" in m._COFACTOR_STIDS_SEED, "NAD+ is a cofactor"


def test_a_failed_derivation_is_never_cached(monkeypatch):
    """Memoising the failure would let one transient reset on pathway 1 build
    every later pathway with the seed, while export_cofactors queries
    separately, succeeds, and ships a mismatched cofactors.csv."""
    from src import neo4j_connector

    calls = {"n": 0}

    def flaky():
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("transient")
        return [{"stable_id": "R-ALL-1"}]

    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", flaky)
    assert m._cofactor_stids() == m._COFACTOR_STIDS_SEED   # first call degrades
    assert m._cofactor_stids() == frozenset({"R-ALL-1"})   # and RETRIES


def test_the_emitter_is_actually_wired_in(monkeypatch):
    """Every other test here exercises the helpers directly, so deleting the
    call in create_pathway_logic_network would leave them all green."""
    import inspect
    src = inspect.getsource(m.create_pathway_logic_network)
    assert "_emit_diagram_set_member_edges(" in src
    assert "diagram_set_member_pairs" in src

    from src import pathway_generator as pg
    gen = inspect.getsource(pg)
    assert "LNG_DIAGRAM_SET_MEMBER" in gen
    # off by default, and the documented kill switch still disables it
    assert '"LNG_DIAGRAM_SET_MEMBER", "0"' in gen
    assert "LNG_DIAGRAM_SET_MEMBER" in pg._FINGERPRINTED_ENV
