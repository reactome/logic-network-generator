"""LNG_CAP_POOLS (deltasignal specs/036): a set input that resolves to SEVERAL
nodes in one virtual reaction -- its alternatives were bundled by the variant
cap -- is recorded as a set pool, so Phase 3 wires members -> pool -> reaction
and the reaction has one input per curated component.

Shape: RAF "MAP2Ks and MAPKs bind to the activated RAF complex". Its input set S
("RAF/MAPK scaffolds") is a CandidateSet; over the cap every candidate landed in
one virtual reaction as a separate required input."""
import pandas as pd
import pytest

import src.logic_network_generator as m
from src import neo4j_connector

S, C, M1, M2 = "R-HSA-S", "R-HSA-C", "R-HSA-M1", "R-HSA-M2"


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(neo4j_connector, "get_reaction_input_output_ids",
                        lambda rid, io: {S, C} if io == "input" else {"R-HSA-OUT"})
    monkeypatch.setattr(neo4j_connector, "get_reaction_io_stoichiometry", lambda rid, io: {})
    monkeypatch.setattr(neo4j_connector, "get_labels",
                        lambda e: ["CandidateSet"] if e == S else ["EntityWithAccessionedSequence"])
    monkeypatch.setattr(m, "_matching_leaves", lambda e: {M1, M2} if e == S else {e})
    monkeypatch.setattr(m, "modifier_isoform_set_ids", lambda: set())


def resolve(present):
    uid_index = {"in": ([], set(present) | {C}, {}), "out": ([], {"R-HSA-OUT"}, {})}
    rmap = pd.DataFrame({"uid": ["vr1"], "input_hash": ["in"], "output_hash": ["out"], "reactome_id": [1]})
    return m._resolve_vr_entities(rmap, uid_index)


def test_bundled_alternatives_become_a_pool(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    ve = resolve({M1, M2})
    assert m._vr_input_pools == {"vr1": {S: {M1, M2}}}
    # members stay in the VR's inputs, so Phase 2 still joins them to producers
    assert set(ve["vr1"][0]) == {M1, M2, C}


def test_a_single_chosen_member_is_not_a_pool(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({M1})
    assert m._vr_input_pools == {}


def test_off_records_no_pools(stub, monkeypatch):
    monkeypatch.delenv("LNG_CAP_POOLS", raising=False)
    resolve({M1, M2})
    assert m._vr_input_pools == {}


def test_pools_do_not_leak_between_pathways(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    resolve({M1, M2})
    assert m._vr_input_pools
    resolve({M1})                     # the next pathway's resolution starts clean
    assert m._vr_input_pools == {}


def test_a_modifier_isoform_set_is_not_pooled(stub, monkeypatch):
    monkeypatch.setenv("LNG_CAP_POOLS", "1")
    monkeypatch.setattr(m, "modifier_isoform_set_ids", lambda: {S})
    monkeypatch.setattr(m, "_matching_leaves", lambda e: {e})
    resolve({M1, M2})
    assert m._vr_input_pools == {}
