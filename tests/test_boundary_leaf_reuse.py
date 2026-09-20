"""A root complex's subunit leaf must not reuse a node that a reaction PRODUCES.

Reusing it welds root complex -> reactions -> produced protein -> (assembly) -> root
complex, a cycle Neo4j never had (1,994 of 2,077 cycle-carrying assembly edges on the
v97 catalog; deltasignal specs/018)."""
import os
import pytest
import src.neo4j_connector as nc
import src.logic_network_generator as lng
from src.logic_network_generator import _emit_boundary_decomposition_edges

C, P, Q = "R-HSA-100", "R-HSA-200", "R-HSA-300"       # C: root complex of P and Q
LABELS = {C: ["Complex"], P: ["EntityWithAccessionedSequence"], Q: ["EntityWithAccessionedSequence"]}


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(nc, "get_labels", lambda stid: LABELS.get(stid, []))
    monkeypatch.setattr(lng, "get_labels", lambda stid: LABELS.get(stid, []), raising=False)
    monkeypatch.setattr(lng, "get_terminal_components", lambda stid: {P, Q} if stid == C else {stid}, raising=False)
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)


def network():
    # C (root complex) -> r1 -> P_out (a PRODUCED copy of P);  Q_free is a root protein with no producer
    data = [
        {"source_id": "u_C", "target_id": "r1", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r1", "target_id": "u_Pout", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "u_Qfree", "target_id": "r2", "pos_neg": "pos", "and_or": "and", "edge_type": "catalyst", "stoichiometry": 1},
        {"source_id": "r2", "target_id": "u_X", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    r2u = {"u_C": C, "u_Pout": P, "u_Qfree": Q, "u_X": "R-HSA-400"}
    return data, r2u


def assembly_sources(data):
    return {e["source_id"]: e for e in data if e["edge_type"] == "assembly"}


def test_produced_node_is_not_reused_as_a_boundary_leaf(stub, monkeypatch):
    monkeypatch.delenv("LNG_BOUNDARY_LEAF_REUSE", raising=False)
    data, r2u = network()
    _emit_boundary_decomposition_edges(data, r2u)
    asm = assembly_sources(data)
    targets = {e["target_id"] for e in asm.values()}
    assert targets == {"u_C"}
    # P's leaf is a FRESH node (u_Pout is produced by r1), Q's leaf reuses the unproduced free node
    assert "u_Pout" not in asm
    assert "u_Qfree" in asm
    fresh = [u for u in asm if u not in ("u_Qfree",)]
    assert len(fresh) == 1 and r2u[fresh[0]] == P
    # and therefore no cycle: nothing produced feeds back into the root complex
    sources_into_C = {e["source_id"] for e in data if e["target_id"] == "u_C"}
    produced = {e["target_id"] for e in data if e["edge_type"] == "output"}
    assert not (sources_into_C & produced)


def test_legacy_reuse_is_available_and_welds_the_cycle(stub, monkeypatch):
    monkeypatch.setenv("LNG_BOUNDARY_LEAF_REUSE", "any")
    data, r2u = network()
    _emit_boundary_decomposition_edges(data, r2u)
    asm = assembly_sources(data)
    assert "u_Pout" in asm          # the old behaviour: the produced copy is reused -> C -> r1 -> u_Pout -> C


def test_bad_mode_is_an_error(stub, monkeypatch):
    monkeypatch.setenv("LNG_BOUNDARY_LEAF_REUSE", "sometimes")
    data, r2u = network()
    with pytest.raises(ValueError):
        _emit_boundary_decomposition_edges(data, r2u)
