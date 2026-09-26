"""A root complex's subunit leaf must not reuse a node that is DOWNSTREAM of that complex.

Reusing it welds root complex -> reactions -> produced protein -> (assembly) -> root
complex, a cycle Neo4j never had (1,994 of 2,077 cycle-carrying assembly edges on the
v97 catalog; deltasignal specs/018). A produced node that is NOT downstream is a real
feed-forward link (hasComponent with no reaction) and must still be reused -- severing
those cost Mitotic G1 28 cases in the first version of this fix."""
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
    # One hasComponent level, derived from the same stub, so these reuse rules are
    # exercised under LNG_BOUNDARY_HIERARCHY=1, the default since specs/030.
    def _components(stid):
        leaves = lng.get_terminal_components(stid)
        return {} if leaves == {stid} else {x: 1 for x in leaves}
    monkeypatch.setattr(nc, "get_complex_components", _components)


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


def test_produced_node_is_not_reused_as_a_boundary_leaf(stub):
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


def test_a_copy_reached_only_through_a_bridge_is_downstream_too(stub):
    # P has two copies: u_Pout (reaction output, downstream of C) and u_Pin (an input copy whose
    # only incoming edge is a diagram bridge from u_Pout -- also downstream). Neither may be reused.
    data, r2u = network()
    data.append({"source_id": "u_Pout", "target_id": "u_Pin", "pos_neg": "pos", "and_or": "and", "edge_type": "diagram_bridge", "stoichiometry": 1})
    data.append({"source_id": "u_Pin", "target_id": "r3", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})
    data.append({"source_id": "r3", "target_id": "u_Y", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1})
    r2u.update({"u_Pin": P, "u_Y": "R-HSA-500"})
    _emit_boundary_decomposition_edges(data, r2u)
    asm = assembly_sources(data)
    assert "u_Pout" not in asm and "u_Pin" not in asm


def test_a_produced_copy_that_is_not_downstream_is_reused(stub):
    # The Mitotic G1 shape: P's only copy is produced by a reaction the root complex C does NOT
    # reach (r2, fed by Q). Reusing it adds no cycle and keeps the feed-forward link
    # P -> (assembly) -> C that Reactome joins only by hasComponent.
    data = [
        {"source_id": "u_C", "target_id": "r1", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r1", "target_id": "u_W", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "u_Qfree", "target_id": "r2", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r2", "target_id": "u_Pelse", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    r2u = {"u_C": C, "u_W": "R-HSA-700", "u_Qfree": Q, "u_Pelse": P}
    _emit_boundary_decomposition_edges(data, r2u)
    asm = assembly_sources(data)
    assert "u_Pelse" in asm and asm["u_Pelse"]["target_id"] == "u_C"      # reused: a real link, no cycle
    assert "u_Qfree" in asm
    assert sum(1 for e in data if e["edge_type"] == "assembly") == 2


def test_two_root_complexes_sharing_a_subunit_each_get_an_acyclic_leaf(stub, monkeypatch):
    # C reaches P's produced copy (u_Pout) -> C gets a fresh leaf. C2 does NOT reach it ->
    # C2 reuses u_Pout (a real feed-forward link). Neither assembly edge closes a cycle.
    data, r2u = network()
    data.append({"source_id": "u_C2", "target_id": "r4", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})
    data.append({"source_id": "r4", "target_id": "u_Z", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1})
    r2u.update({"u_C2": "R-HSA-101", "u_Z": "R-HSA-600"})
    LABELS["R-HSA-101"] = ["Complex"]
    monkeypatch.setattr(lng, "get_terminal_components", lambda stid: {P, Q} if stid in (C, "R-HSA-101") else {stid}, raising=False)
    _emit_boundary_decomposition_edges(data, r2u)
    into = {t_: {e["source_id"] for e in data if e["edge_type"] == "assembly" and e["target_id"] == t_} for t_ in ("u_C", "u_C2")}
    assert "u_Pout" not in into["u_C"] and "u_Pout" in into["u_C2"]
    assert any(r2u.get(u) == P for u in into["u_C"])
    # acyclicity: no assembly source into a complex is reachable from that complex
    succ = {}
    for e in data: succ.setdefault(e["source_id"], []).append(e["target_id"])
    def reach(u):
        seen, st = set(), [u]
        while st:
            x = st.pop()
            for v in succ.get(x, []):
                if v not in seen: seen.add(v); st.append(v)
        return seen
    for cpx, srcs in into.items():
        assert not (srcs & reach(cpx))


def test_removed_flag_is_an_error_not_a_noop(stub, monkeypatch):
    # The escape hatch was removed; a stale value must fail loudly rather than
    # letting a run silently measure the default while claiming the old mode.
    for val in ("any", "downstream_free", "unproduced"):
        monkeypatch.setenv("LNG_BOUNDARY_LEAF_REUSE", val)
        data, r2u = network()
        with pytest.raises(ValueError, match="was removed"):
            _emit_boundary_decomposition_edges(data, r2u)
