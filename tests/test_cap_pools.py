"""LNG_CAP_POOLS (deltasignal specs/036): in a reaction the variant cap BUNDLED
(every alternative of every set merged into one variant), a set input becomes
one pool: each alternative is an AND unit of its leaves, the pool ORs the
alternatives, and the reaction reads the pool once.

Shape: RAF "MAP2Ks and MAPKs bind to the activated RAF complex". Its input set S
("RAF/MAPK scaffolds") has alternatives A (a single protein) and K (a complex,
leaves K1, K2); over the cap all three leaves were separate required inputs.

Review of #combined-037: the first version pooled ANY set resolving to several
nodes, which also fires in UNCAPPED reactions where the one chosen alternative
is a complex (miR-93 RISC in PIP3: its four subunits were averaged, so a miR-93
knockout read 0.75x instead of 0)."""
import pandas as pd
import pytest

import src.logic_network_generator as m
from src import neo4j_connector

S, C, A, K, K1, K2 = "R-HSA-S", "R-HSA-C", "R-HSA-A", "R-HSA-K", "R-HSA-K1", "R-HSA-K2"
LEAVES = {S: {A, K1, K2}, A: {A}, K: {K1, K2}, C: {C}}


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(neo4j_connector, "get_reaction_input_output_ids",
                        lambda rid, io: {S, C} if io == "input" else {"R-HSA-OUT"})
    monkeypatch.setattr(neo4j_connector, "get_reaction_io_stoichiometry",
                        lambda rid, io: {S: 2} if io == "input" else {})
    monkeypatch.setattr(neo4j_connector, "get_labels",
                        lambda e: ["CandidateSet"] if e == S else ["EntityWithAccessionedSequence"])
    monkeypatch.setattr(neo4j_connector, "get_set_members", lambda e: [A, K] if e == S else [])
    monkeypatch.setattr(m, "_matching_leaves", lambda e: LEAVES.get(e, {e}))
    monkeypatch.setattr(m, "modifier_isoform_set_ids", lambda: set())
    monkeypatch.setattr(m, "CAPPED_IDS", {"1"})


def resolve(present, rid=1):
    uid_index = {"in": ([], set(present) | {C}, {}), "out": ([], {"R-HSA-OUT"}, {})}
    rmap = pd.DataFrame({"uid": ["vr1"], "input_hash": ["in"], "output_hash": ["out"], "reactome_id": [rid]})
    return m._resolve_vr_entities(rmap, uid_index)


def test_a_capped_bundle_becomes_a_pool_of_alternatives(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    ve = resolve({A, K1, K2})
    assert m._vr_input_pools == {"vr1": {S: {"alts": {A: {A}, K: {K1, K2}}, "stoich": 2}}}
    assert set(ve["vr1"][0]) == {A, K1, K2, C}      # members stay VR inputs (Phase 2 joins)


def test_an_uncapped_reaction_is_never_pooled(stub, monkeypatch):
    # the miR-93 case: one chosen alternative that is a complex resolves to its
    # subunits, which ARE co-required
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({K1, K2}, rid=2)
    assert m._vr_input_pools == {}


def test_the_cap_gate_itself(stub, monkeypatch):
    # two alternatives present in a reaction the cap did NOT bundle: not pooled
    # (only a capped reaction merged its alternatives)
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({A, K1, K2}, rid=2)
    assert m._vr_input_pools == {}


def test_one_alternative_present_is_not_a_pool(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({K1, K2})
    assert m._vr_input_pools == {}


def test_a_leaf_another_input_also_maps_to_stays_direct(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    monkeypatch.setattr(m, "_matching_leaves", lambda e: {A, K1, K2} if e == S else ({K1, K2, C} if e == C else LEAVES.get(e, {e})))
    monkeypatch.setattr(neo4j_connector, "get_labels",
                        lambda e: ["CandidateSet"] if e in (S,) else (["Complex"] if e == C else ["EntityWithAccessionedSequence"]))
    monkeypatch.setattr(m, "_map_annotated_entity_to_nodes",
                        lambda e, mem: {A, K1, K2} & mem if e == S else ({K1} if e == C else {e}))
    resolve({A, K1, K2})
    alts = m._vr_input_pools["vr1"][S]["alts"]
    assert alts == {A: {A}, K: {K2}}                       # K1 is also C's node: not pooled


def test_off_and_no_leak_between_pathways(stub, monkeypatch):
    monkeypatch.delenv("LNG_CAP_POOLS", raising=False)
    resolve({A, K1, K2})
    assert m._vr_input_pools == {}
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({A, K1, K2})
    assert m._vr_input_pools
    resolve({A}, rid=2)
    assert m._vr_input_pools == {}


def test_emission_wires_alternatives_as_and_units_into_one_pool():
    reg = {(x, "vr1", "input"): f"u-{x}" for x in (A, K1, K2, C)}
    r2u, data = {}, []
    pools = {S: {"alts": {A: {A}, K: {K1, K2}}, "stoich": 2}}
    m._emit_vr_inputs("vr1", [A, K1, K2, C], {C: 1}, pools, reg, r2u, data, "R-HSA-RX")
    into_vr = [(e["source_id"], e["edge_type"], e["stoichiometry"]) for e in data if e["target_id"] == "vr1"]
    pool = next(u for u, s in r2u.items() if s == S)
    alt = next(u for u, s in r2u.items() if s == K)
    assert sorted(into_vr) == sorted([(pool, "input", 2), ("u-" + C, "input", 1)])  # one input per component
    assert {(e["source_id"], e["and_or"], e["edge_type"]) for e in data if e["target_id"] == pool} == \
        {("u-" + A, "or", "set_member"), (alt, "or", "set_member")}
    assert {(e["source_id"], e["and_or"], e["edge_type"]) for e in data if e["target_id"] == alt} == \
        {("u-" + K1, "and", "assembly"), ("u-" + K2, "and", "assembly")}
    assert not any(e["target_id"] == "vr1" and e["source_id"] in ("u-" + A, "u-" + K1, "u-" + K2) for e in data)


def test_without_pools_emission_is_the_plain_input_fan_in():
    reg = {(x, "vr1", "input"): f"u-{x}" for x in (A, C)}
    data = []
    m._emit_vr_inputs("vr1", [A, C], {}, {}, reg, {}, data, "R-HSA-RX")
    assert sorted(e["source_id"] for e in data) == ["u-" + A, "u-" + C]
    assert {e["edge_type"] for e in data} == {"input"}
