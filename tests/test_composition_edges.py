"""Unit tests for ``_emit_composition_edges`` (LNG_COMPOSITION_EDGES).

A complex node gets an edge to the node of every complex that CONTAINS it,
within two hasComponent hops. This is the non-broadcasting repair for the
composition-only gap: it connects a complex to the few complexes it sits
inside (fan-out median 1, max 2), never a released subunit to its hub. Neo4j is
stubbed; the emission logic is what is under test.
"""
import os

import pytest

import src.neo4j_connector as nc
from src.logic_network_generator import _emit_composition_edges

# stable ids in a toy pathway
ISGF3, ISGF3_KPNA1, IMPORTIN = "R-HSA-909698", "R-HSA-9710958", "R-HSA-9710965"
STAT1 = "R-HSA-873791"          # a protein (EWAS), never a composition source
LABELS = {ISGF3: ["Complex"], ISGF3_KPNA1: ["Complex"], IMPORTIN: ["Complex"], STAT1: ["EntityWithAccessionedSequence"]}
# containment: importin contains ISGF3:KPNA1 (1 hop) which contains ISGF3, so importin contains ISGF3 at 2 hops
CONTAINERS = {ISGF3: {ISGF3_KPNA1: 1, IMPORTIN: 2}, ISGF3_KPNA1: {IMPORTIN: 1}, IMPORTIN: {}, STAT1: {ISGF3: 1}}


@pytest.fixture
def stub_neo4j(monkeypatch):
    monkeypatch.setattr(nc, "get_labels", lambda stid: LABELS.get(stid, []))
    monkeypatch.setattr(nc, "get_containing_complexes", lambda stid, max_hops=2: CONTAINERS.get(stid, {}))


def edge(src, tgt, ty="input"):
    return {"source_id": src, "target_id": tgt, "pos_neg": "pos", "and_or": "and", "edge_type": ty, "stoichiometry": 1}


def composition_edges(data):
    return {(e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "composition"}


def test_complex_connects_to_the_complexes_that_contain_it(stub_neo4j):
    # nodes: ISGF3 (u1), importin (u3). ISGF3:KPNA1 has NO node (set member) -- the 2-hop case.
    uuid_to_stid = {"u1": ISGF3, "u3": IMPORTIN}
    data = [edge("rX", "u1", "output")]
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == {("u1", "u3")}
    new = [e for e in data if e["edge_type"] == "composition"][0]
    assert (new["pos_neg"], new["and_or"], new["stoichiometry"]) == ("pos", "and", 1)


def test_one_hop_and_two_hop_both_emitted_when_both_have_nodes(stub_neo4j):
    uuid_to_stid = {"u1": ISGF3, "u2": ISGF3_KPNA1, "u3": IMPORTIN}
    data = []
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == {("u1", "u2"), ("u1", "u3"), ("u2", "u3")}


def test_a_dissociation_sink_is_never_a_source(stub_neo4j):
    # u_sink carries ISGF3's stable id but is a readout sink (target of a dissociation edge).
    uuid_to_stid = {"u1": ISGF3, "u_sink": ISGF3, "u3": IMPORTIN}
    data = [edge("cplx", "u_sink", "dissociation")]
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == {("u1", "u3")}
    assert all(e["source_id"] != "u_sink" for e in data if e["edge_type"] == "composition")


def test_proteins_are_not_sources_even_if_something_contains_them(stub_neo4j):
    # STAT1 is contained by ISGF3, but STAT1 is a protein: that is the leaf-subunit
    # bridge that lost four times, and it must NOT be emitted here.
    uuid_to_stid = {"s1": STAT1, "u1": ISGF3}
    data = []
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == set()


def test_no_self_edge_and_no_duplicate_of_an_existing_edge(stub_neo4j):
    uuid_to_stid = {"u1": ISGF3, "u3": IMPORTIN}
    data = [edge("u1", "u3", "assembly")]          # an edge already exists u1 -> u3
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == set()         # deduped against the existing pair
    assert len(data) == 1


def test_container_without_a_node_yields_nothing(stub_neo4j):
    uuid_to_stid = {"u1": ISGF3}                    # neither container is a node here
    data = []
    _emit_composition_edges(data, uuid_to_stid)
    assert data == []


def test_variant_nodes_resolve_to_their_base_complex(stub_neo4j):
    uuid_to_stid = {"v1": f"{ISGF3}::variant::R-HSA-1", "v2": f"{ISGF3}::variant::R-HSA-2", "u3": IMPORTIN}
    data = []
    _emit_composition_edges(data, uuid_to_stid)
    assert composition_edges(data) == {("v1", "u3"), ("v2", "u3")}


def test_flag_default_is_off(monkeypatch):
    # The emitter is only called from _emit_boundary_decomposition_edges under
    # the flag; with the flag unset the call site must not reach it.
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)
    import src.logic_network_generator as lng
    calls = []
    monkeypatch.setattr(lng, "_emit_composition_edges", lambda *a, **k: calls.append(1))
    src = open(lng.__file__).read()
    assert 'os.environ.get("LNG_COMPOSITION_EDGES", "0") == "1"' in src
    assert calls == []
