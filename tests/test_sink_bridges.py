"""Acyclic sink bridges (LNG_SINK_BRIDGES): a released subunit reconnects to its consuming copies only where no cycle is closed."""
import pytest
from src.logic_network_generator import _emit_sink_bridge_edges

P, Q = "R-HSA-200", "R-HSA-300"

def net():
    # T (terminal complex) --dissociation--> P_sink ;  P_in --input--> r2 --output--> X   (P_in is a consuming copy, NOT upstream of the sink)
    # Q_sink from T too; Q_in consumes into r3 whose output feeds ... T's producer r1 (Q_in IS upstream of the sink -> bridging Q would close a cycle)
    data = [
        {"source_id": "u_Qin", "target_id": "r3", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r3", "target_id": "u_A", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "u_A", "target_id": "r1", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r1", "target_id": "u_T", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "u_T", "target_id": "u_Psink", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "u_T", "target_id": "u_Qsink", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "u_Pin", "target_id": "r2", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r2", "target_id": "u_X", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    r2u = {"u_Psink": P, "u_Pin": P, "u_Qsink": Q, "u_Qin": Q, "u_T": "R-HSA-100", "u_A": "R-HSA-400", "u_X": "R-HSA-500"}
    return data, r2u

def bridges(data): return [(e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "sink_bridge"]

def test_sink_bridges_to_a_consuming_copy_that_is_not_upstream():
    data, r2u = net(); n = _emit_sink_bridge_edges(data, r2u)
    assert n == 1 and bridges(data) == [("u_Psink", "u_Pin")]
    b = [e for e in data if e["edge_type"] == "sink_bridge"][0]
    assert b["pos_neg"] == "pos" and b["and_or"] == "or"

def test_no_bridge_when_the_consumer_is_upstream_of_the_sink():
    data, r2u = net(); _emit_sink_bridge_edges(data, r2u)
    assert ("u_Qsink", "u_Qin") not in bridges(data)      # Q_in -> r3 -> A -> r1 -> T -> Q_sink: would close a cycle

def test_result_is_acyclic_where_input_was_acyclic():
    data, r2u = net(); _emit_sink_bridge_edges(data, r2u)
    succ = {}
    for e in data: succ.setdefault(e["source_id"], []).append(e["target_id"])
    def reach(u):
        seen, st = set(), [u]
        while st:
            x = st.pop()
            for v in succ.get(x, []):
                if v not in seen: seen.add(v); st.append(v)
        return seen
    for s, t in bridges(data):
        assert s not in reach(t)

def test_only_true_sinks_bridge_and_a_sink_with_no_consumer_copy_is_left_alone():
    data, r2u = net()
    data.append({"source_id": "u_Psink", "target_id": "r9", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})  # P_sink already consumed -> not a sink
    r2u["r9"] = "R-HSA-900"
    n = _emit_sink_bridge_edges(data, r2u)
    assert n == 0

def test_variant_copies_join_on_the_base_stable_id():
    data, r2u = net(); r2u["u_Pin"] = P + "::variant::R-HSA-1"
    n = _emit_sink_bridge_edges(data, r2u)
    assert n == 1 and bridges(data) == [("u_Psink", "u_Pin")]

def test_idempotent_and_order_free():
    d1, r1 = net(); _emit_sink_bridge_edges(d1, r1)
    d2, r2 = net(); d2.reverse(); _emit_sink_bridge_edges(d2, r2)
    assert sorted(bridges(d1)) == sorted(bridges(d2))
