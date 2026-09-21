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


def two_sinks():
    """Two sinks of ONE entity plus one real consumer: the cascade case.

    The consumer filter must be evaluated against the pre-emitter out-degree,
    or the first sink to receive a bridge becomes an eligible 'consumer' for
    the second (69% of the first catalog run's edges were exactly this)."""
    data = [
        {"source_id": "u_T", "target_id": "s_aaa", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "u_T", "target_id": "s_zzz", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "u_cons", "target_id": "r9", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r9", "target_id": "u_X", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    r2u = {"s_aaa": P, "s_zzz": P, "u_cons": P, "u_T": "R-HSA-100", "u_X": "R-HSA-500", "r9": "R-HSA-900"}
    return data, r2u


def test_a_sink_is_never_a_bridge_target():
    data, r2u = two_sinks()
    n = _emit_sink_bridge_edges(data, r2u)
    assert n == 2
    assert sorted(bridges(data)) == [("s_aaa", "u_cons"), ("s_zzz", "u_cons")]
    sinks = {e["target_id"] for e in data if e["edge_type"] == "dissociation"}
    assert not [b for b in bridges(data) if b[1] in sinks]


def test_sink_visit_order_does_not_change_the_emitted_set():
    # identical graph, sink labels swapped: uuid4 ordering used to pick a different EDGE
    d1, r1 = two_sinks()
    d2, r2 = two_sinks()
    swap = {"s_aaa": "s_zzz", "s_zzz": "s_aaa"}
    d2 = [{**e, "source_id": swap.get(e["source_id"], e["source_id"]), "target_id": swap.get(e["target_id"], e["target_id"])} for e in d2]
    _emit_sink_bridge_edges(d1, r1); _emit_sink_bridge_edges(d2, r2)
    assert {(swap.get(s, s), t) for s, t in bridges(d1)} == set(bridges(d2))


def accumulation():
    """Four nodes of ONE entity X: two sinks (S0, S2) and two consumers (c, c2), wired so
    that a bridge emitted for S0 changes what c2 reaches.

        c  -> rA -> T2p --dissoc--> S2        (c reaches S2)
        c2 -> rC -> T0  --dissoc--> S0        (c2 reaches S0)

    Processing S0 first emits S0 -> c. That makes c2 reach S2 *through* the new edge, so
    S2 -> c2 is now cycle-closing and must be skipped. With a stale reachability cache it
    is emitted and closes S2 -> c2 -> rC -> T0 -> S0 -> c -> rA -> T2p -> S2."""
    X = P
    data = [
        {"source_id": "u_T0", "target_id": "s0", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "u_T2p", "target_id": "s2", "pos_neg": "pos", "and_or": "and", "edge_type": "dissociation", "stoichiometry": 1},
        {"source_id": "c", "target_id": "rA", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "rA", "target_id": "u_T2p", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "c2", "target_id": "rC", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "rC", "target_id": "u_T0", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    r2u = {"s0": X, "s2": X, "c": X, "c2": X, "u_T0": "R-HSA-100", "u_T2p": "R-HSA-101", "rA": "R-HSA-900", "rC": "R-HSA-901"}
    return data, r2u


def acyclic(data):
    succ = {}
    for e in data: succ.setdefault(e["source_id"], []).append(e["target_id"])
    def reach(u):
        seen, st = set(), [u]
        while st:
            x = st.pop()
            for v in succ.get(x, []):
                if v not in seen: seen.add(v); st.append(v)
        return seen
    return all(s not in reach(t) for s, t in bridges(data))


def test_accumulation_a_bridge_makes_a_later_candidate_cycle_closing():
    data, r2u = accumulation()
    n = _emit_sink_bridge_edges(data, r2u)
    assert acyclic(data), f"cycle-closing bridge emitted: {bridges(data)}"
    assert n == 1 and bridges(data) == [("s0", "c")]


def test_sink_visit_order_is_a_function_of_the_data_not_the_uuids():
    # Same graph, sink labels swapped. Ordering by (stable id, first appearance) must give
    # the structurally same bridge; ordering by uuid gives a different one.
    d1, r1 = accumulation()
    _emit_sink_bridge_edges(d1, r1)
    d2, r2 = accumulation()
    swap = {"s0": "s2", "s2": "s0"}
    d2 = [{**e, "source_id": swap.get(e["source_id"], e["source_id"]), "target_id": swap.get(e["target_id"], e["target_id"])} for e in d2]
    r2 = {swap.get(k, k): v for k, v in r2.items()}
    _emit_sink_bridge_edges(d2, r2)
    assert acyclic(d2)
    assert {(swap.get(s, s), t) for s, t in bridges(d1)} == set(bridges(d2)), (bridges(d1), bridges(d2))


def test_negative_fanout_cap_is_an_error(monkeypatch):
    monkeypatch.setenv("LNG_SINK_BRIDGE_MAX_FANOUT", "-1")
    data, r2u = net()
    with pytest.raises(ValueError):
        _emit_sink_bridge_edges(data, r2u)


def test_fanout_cap_skips_broadcast_sinks(monkeypatch):
    data, r2u = net()
    # give P three more acyclic consuming copies -> fan-out 4
    for i in range(3):
        data.append({"source_id": f"u_Pin{i}", "target_id": f"r2{i}", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})
        data.append({"source_id": f"r2{i}", "target_id": f"u_X{i}", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1})
        r2u[f"u_Pin{i}"] = P; r2u[f"u_X{i}"] = f"R-HSA-50{i}"
    monkeypatch.setenv("LNG_SINK_BRIDGE_MAX_FANOUT", "3")
    assert _emit_sink_bridge_edges(data, r2u) == 0
    monkeypatch.setenv("LNG_SINK_BRIDGE_MAX_FANOUT", "4")
    data, r2u = net()
    for i in range(3):
        data.append({"source_id": f"u_Pin{i}", "target_id": f"r2{i}", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})
        data.append({"source_id": f"r2{i}", "target_id": f"u_X{i}", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1})
        r2u[f"u_Pin{i}"] = P; r2u[f"u_X{i}"] = f"R-HSA-50{i}"
    assert _emit_sink_bridge_edges(data, r2u) == 4
