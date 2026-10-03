"""A SET component of a root complex is one node, fed by its members over OR
`set_member` edges (LNG code review F7, deltasignal specs/044). It used to be
flattened into its members, each a required AND input of the complex, so a
knockout of one alternative removed the complex.

Shape: GPVI, "VAV1 Rho/Rac effectors:GDP" = GDP + CandidateSet{RHOA, RAC1, CDC42}."""
import pytest
import src.neo4j_connector as nc
import src.logic_network_generator as lng
import src.reaction_generator as rg
from src.logic_network_generator import _emit_boundary_decomposition_edges

CX, GDP, SET, RHOA, RAC1, CDC42 = "R-HSA-CX", "R-ALL-GDP", "R-HSA-SET", "R-HSA-RHOA", "R-HSA-RAC1", "R-HSA-CDC42"
INNER, SUB1, SUB2 = "R-HSA-INNER", "R-HSA-SUB1", "R-HSA-SUB2"      # a member that is itself a complex
LABELS = {CX: ["Complex"], SET: ["EntitySet", "CandidateSet"], INNER: ["Complex"],
          GDP: ["SimpleEntity"], **{x: ["EntityWithAccessionedSequence"] for x in (RHOA, RAC1, CDC42, SUB1, SUB2)}}
COMPONENTS = {CX: {GDP: 1, SET: 1}, INNER: {SUB1: 1, SUB2: 1}}
MEMBERS = {SET: {RHOA, RAC1, CDC42, INNER}}


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(nc, "get_labels", lambda s: LABELS.get(s, []))
    monkeypatch.setattr(nc, "get_complex_components", lambda s: COMPONENTS.get(s, {}))
    monkeypatch.setattr(nc, "get_set_members", lambda s: MEMBERS.get(s, set()))
    # as the real walk does: a set or complex flattens to all its leaves
    leaves = {SET: {RHOA, RAC1, CDC42, SUB1, SUB2}, INNER: {SUB1, SUB2}, CX: {GDP, RHOA, RAC1, CDC42, SUB1, SUB2}}
    monkeypatch.setattr(lng, "get_terminal_components", lambda s: leaves.get(s, {s}), raising=False)
    monkeypatch.setattr(rg, "_modifier_set_cache", set())
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)


def run():
    data = [{"source_id": "u_CX", "target_id": "r1", "pos_neg": "pos", "and_or": "and",
             "edge_type": "input", "stoichiometry": 1},
            {"source_id": "r1", "target_id": "u_Z", "pos_neg": "pos", "and_or": "or",
             "edge_type": "output", "stoichiometry": 1}]
    r2u = {"u_CX": CX, "u_Z": "R-HSA-Z"}
    _emit_boundary_decomposition_edges(data, r2u)
    return data, r2u


def test_the_set_is_one_node_feeding_the_complex(stub):
    data, r2u = run()
    into_cx = [(r2u[e["source_id"]], e["edge_type"], e["and_or"]) for e in data if e["target_id"] == "u_CX"]
    assert sorted(into_cx) == [(GDP, "assembly", "and"), (SET, "assembly", "and")]


def test_members_feed_the_set_node_as_alternatives(stub):
    data, r2u = run()
    set_uuid = next(u for u, s in r2u.items() if s == SET)
    into_set = {(r2u[e["source_id"]], e["edge_type"], e["and_or"]) for e in data if e["target_id"] == set_uuid}
    assert into_set == {(m, "set_member", "or") for m in (RHOA, RAC1, CDC42, INNER)}
    # no member is a required input of the complex any more
    assert not any(e["target_id"] == "u_CX" and r2u[e["source_id"]] in (RHOA, RAC1, CDC42) for e in data)


def test_a_complex_member_is_built_from_its_components(stub):
    data, r2u = run()
    inner = next(u for u, s in r2u.items() if s == INNER)
    assert {(r2u[e["source_id"]], e["and_or"]) for e in data if e["target_id"] == inner} == \
        {(SUB1, "and"), (SUB2, "and")}


def test_every_new_edge_is_in_the_boundary_layer(stub):
    data, _ = run()
    assert all(e.get("_boundary") for e in data[2:])
