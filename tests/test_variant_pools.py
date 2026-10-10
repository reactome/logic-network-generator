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


import pytest
from src import variant_keys as vk


@pytest.fixture
def world():
    vk._lookups.update(labels=lambda s: {"R": ["Complex"], "S": ["EntitySet"], "SET": ["EntitySet"]}.get(s, ["EWAS"]),
                       components=lambda s: {"R": ["S", "G"]}.get(s, []),
                       members=lambda s: {"S": ["a", "b", "c"], "SET": ["m1", "m2"]}.get(s, []),
                       atomic_sets=lambda: set(), max_variants=lambda: 512)
    vk.reset_caches()
    yield
    vk._lookups.clear(); vk.reset_caches()


def test_pool_reference_reads_produced_variants_else_roots(world):
    ids = {"p": "R::pool", "v1": "R::variant::S=a", "v1b": "R::variant::S=a",
           "v2": "R::variant::S=b", "plain": "R", "rx": "RX"}
    data = [edge("rx", "v1", "output"),             # v1 produced, v1b an unfed copy of the same key
            edge("p", "rx2", "regulator")]
    lng._wire_variant_pool_refs(data, ids)
    members = sorted(e["source_id"] for e in data if e["edge_type"] == "variant_pool")
    assert "v1" in members and "v2" in members and "v1b" not in members and "plain" not in members
    new = [m for m in members if m not in ("v1", "v2")]
    assert [ids[m] for m in new] == ["R::variant::S=c"]   # S=c existed nowhere: a new root


def test_a_pooled_bare_set_reads_its_members(world):
    # vn5 RPL10: a bare set's variant keys are its MEMBERS' keys
    ids = {"p": "SET::pool", "x": "m1"}
    data = [edge("p", "rx", "regulator")]
    lng._wire_variant_pool_refs(data, ids)
    members = sorted(ids[e["source_id"]] for e in data if e["edge_type"] == "variant_pool")
    assert members == ["m1", "m2"]                       # m1's existing node, m2 a new root


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


def redge(s, t, et, rx, pn="pos"):
    e = edge(s, t, et, pn); e["edge_reaction_id"] = rx; return e


def test_a_capped_plain_output_feeds_its_variants(world):
    # vn5 RAF: one copy writes plain R; consumers read R's variant keys
    ids = {"out": "R", "va": "R::variant::S=a", "vb": "R::variant::S=b", "vc": "R::variant::S=c"}
    data = [redge("rx", "out", "output", "P"),
            redge("va", "c1", "input", "F"), redge("vb", "c2", "input", "F"),
            redge("vc", "c3", "input", "UNRELATED")]           # not a curated follower
    lng._wire_capped_outputs(data, ids, {"P": {"F"}})
    split = sorted((e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "variant_split")
    assert split == [("out", "va"), ("out", "vb")]


def test_split_skips_a_variant_another_reaction_produces(world):
    ids = {"out": "R", "va": "R::variant::S=a", "vb": "R::variant::S=b"}
    data = [redge("rx", "out", "output", "P"), redge("rx9", "va", "output", "Q"),
            redge("va", "c1", "input", "F"), redge("vb", "c2", "input", "F")]
    lng._wire_capped_outputs(data, ids, {"P": {"F"}})
    split = sorted(e["target_id"] for e in data if e["edge_type"] == "variant_split")
    assert split == ["vb"]                              # va already has a producer


def test_an_unproduced_plain_node_is_not_split(world):
    ids = {"out": "R", "va": "R::variant::S=a"}
    data = [redge("va", "rx2", "input", "F")]
    lng._wire_capped_outputs(data, ids, {"P": {"F"}})
    assert not [e for e in data if e["edge_type"] == "variant_split"]


def test_members_of_one_set_catalyst_deplete_once(world):
    # review of vn6: six DUSP members of ONE set each depleted the MAPK3 dimer.
    # Their ids share no prefix; the curated participant (_origin) groups them.
    ids = {"d1": "DUSP1", "d2": "DUSP6", "o": "OTHER", "t": "MAPK3"}
    data = [dict(edge("d1", "t", "depletion", "neg"), _origin="DUSPSET"),
            dict(edge("d2", "t", "depletion", "neg"), _origin="DUSPSET"),
            dict(edge("o", "t", "depletion", "neg"), _origin="OTHER")]
    lng._pool_variant_depleters(data, 0, ids)
    dep = [e for e in data if e["edge_type"] == "depletion"]
    assert len(dep) == 2 and all("_origin" not in e for e in data)
    pool = next(e["source_id"] for e in dep if e["source_id"] != "o")
    assert ids[pool] == "DUSPSET::pool"


def test_an_over_cap_set_pool_reads_its_members(world, monkeypatch):
    # review of vn6: RAS GEFs (645 variants) was a pool with no input
    monkeypatch.setattr(vk, "_max_variants", lambda: 1)
    ids = {"p": "SET::pool", "x": "m1", "y": "m2::variant::Q=z", "other": "zz"}
    data = []
    lng._wire_variant_pool_refs(data, ids)
    members = sorted(ids[e["source_id"]] for e in data if e["edge_type"] == "variant_pool")
    assert members == ["m1", "m2::variant::Q=z"]
