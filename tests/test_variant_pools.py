"""specs/046: pool nodes for variant depleters and for capped-reaction pool
references, and the regulator's produced-copy preference (no database)."""
from src import logic_network_generator as lng


def edge(s, t, et, pn="pos", ao="and"):
    return {"source_id": s, "target_id": t, "pos_neg": pn, "and_or": ao, "edge_type": et, "stoichiometry": 1}


def test_variants_of_one_depleter_become_one_pool_edge():
    ids = {"c1": "C::variant::S=a", "c2": "C::variant::S=b", "d": "D", "tp53": "TP53"}
    data = [edge("x", "y", "input"),
            edge("c1", "tp53", "depletion", "neg"), edge("c2", "tp53", "depletion", "neg"),
            edge("d", "tp53", "depletion", "neg")]
    lng._pool_variant_depleters(data, 1, ids)
    dep = [e for e in data if e["edge_type"] == "depletion"]
    assert len(dep) == 2                              # D kept, C's two variants -> one pool
    pool = next(e["source_id"] for e in dep if e["source_id"] != "d")
    assert ids[pool] == "C::pool"
    assert sorted(e["source_id"] for e in data if e["edge_type"] == "variant_pool" and e["target_id"] == pool) == ["c1", "c2"]
    assert all(e["and_or"] == "or" for e in data if e["edge_type"] == "variant_pool")
    assert data[0] == edge("x", "y", "input")         # edges before the pass untouched


def test_a_single_depleter_is_left_alone():
    ids = {"c1": "C::variant::S=a", "t": "T"}
    data = [edge("c1", "t", "depletion", "neg")]
    lng._pool_variant_depleters(data, 0, ids)
    assert data == [edge("c1", "t", "depletion", "neg")]


def test_pool_reference_reads_produced_variants_else_roots():
    ids = {"p": "R::pool", "v1": "R::variant::S=a", "v1b": "R::variant::S=a",
           "v2": "R::variant::S=b", "plain": "R", "rx": "RX"}
    data = [edge("rx", "v1", "output"),             # v1 produced, v1b an unfed copy of the same key
            edge("p", "rx2", "regulator")]
    lng._wire_variant_pool_refs(data, ids)
    members = sorted(e["source_id"] for e in data if e["edge_type"] == "variant_pool")
    assert members == ["v1", "v2"]                   # produced copy of S=a; S=b only exists as a root


def test_regulator_prefers_the_copy_a_preceding_reaction_produces(monkeypatch):
    import pandas as pd
    monkeypatch.setenv("LNG_COMPLEX_AS_NODE", "1")
    key = "PTEN::variant::S=x"
    cat = pd.DataFrame([{"reaction_id": "R-DEPHOS", "entity_id": key, "edge_type": "catalyst",
                         "uuid": "u", "reaction_uuid": "rxcopy"}])
    empty = pd.DataFrame(columns=cat.columns)
    registry = {(key, "vr_in", "input"): "unfed", (key, "vr_a", "output"): "fed_other",
                (key, "vr_b", "output"): "fed_preceding"}
    data, ids = [], {}
    lng.append_regulators(cat, empty, empty, data, ids, entity_uuid_registry=registry,
                          produced_by_eid={key: [("fed_other", "R-OTHER"), ("fed_preceding", "R-TRANSLATE")]},
                          preceding_by_reaction={"R-DEPHOS": {"R-TRANSLATE"}})
    assert [e["source_id"] for e in data] == ["fed_preceding"]
    data = []
    lng.append_regulators(cat, empty, empty, data, ids, entity_uuid_registry=registry)
    assert [e["source_id"] for e in data] == ["unfed"]    # default: first registry entry (unchanged)
