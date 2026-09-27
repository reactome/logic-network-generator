"""LNG_SET_POOL (deltasignal specs/033): a bare EntitySet catalyst or regulator
becomes ONE pool node fed by its members, instead of every member wired onto
every reaction copy as a separate required (catalyst, positive regulator) or
blocking (negative regulator) term.

Shape: RAF/MAP kinase. "RAF activating kinases" (S) catalyses "Phosphorylation
of RAF", which exists as several reaction copies; each copy used to receive all
members as separate AND catalysts."""
from typing import Any, Dict, List
from unittest.mock import patch

import pandas as pd
import pytest

with patch("py2neo.Graph"):
    import src.logic_network_generator as lng
    from src.logic_network_generator import append_regulators

S, K = "R-HSA-S", "R-HSA-K"                  # S = a set; K = a complex (not a set)
MEMBERS = {S: [("R-HSA-M1", 1), ("R-HSA-M2", 1), ("R-HSA-M3", 1)],
           K: [("R-HSA-A", 1), ("R-HSA-B", 1)]}
LABELS = {S: ["DefinedSet"], K: ["Complex"]}


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(lng, "_decompose_regulator_entity",
                        lambda e, variant_decomposition=False, bundle_complex=False:
                        MEMBERS.get(e, [(e, 1)]))
    import src.neo4j_connector as nc
    monkeypatch.setattr(nc, "get_labels", lambda e: LABELS.get(e, []))
    monkeypatch.setattr(nc, "get_set_members", lambda e: [m for m, _ in MEMBERS.get(e, [])])
    monkeypatch.setattr(lng, "modifier_isoform_set_ids", lambda: set())
    monkeypatch.delenv("LNG_SET_MEMBERS_OR", raising=False)


def row(entity, rxn):
    return {"reaction_id": "R-HSA-RX", "entity_id": entity, "edge_type": "x",
            "uuid": "u", "reaction_uuid": rxn}


def run(cat=(), neg=(), pos=()):
    data: List[Dict[str, Any]] = []
    r2u: Dict[str, str] = {}
    append_regulators(pd.DataFrame([row(*x) for x in cat]), pd.DataFrame([row(*x) for x in neg]),
                      pd.DataFrame([row(*x) for x in pos]), data, r2u)
    return data, r2u


def reach(data, src):
    fwd: Dict[str, set] = {}
    for e in data:
        fwd.setdefault(e["source_id"], set()).add(e["target_id"])
    seen, stack = {src}, [src]
    while stack:
        for t in fwd.get(stack.pop(), ()):
            if t not in seen:
                seen.add(t)
                stack.append(t)
    return seen


def test_off_is_the_legacy_member_fan_out(stub, monkeypatch):
    monkeypatch.delenv("LNG_SET_POOL", raising=False)
    data, r2u = run(cat=[(S, "rx1"), (S, "rx2")])
    assert len(data) == 6                                    # 3 members x 2 copies
    assert {e["and_or"] for e in data} == {"and"} and {e["edge_type"] for e in data} == {"catalyst"}
    assert not any(e["edge_type"] == "set_member" for e in data)


def test_on_one_pool_node_feeds_every_copy_once(stub, monkeypatch):
    monkeypatch.setenv("LNG_SET_POOL", "1")
    data, r2u = run(cat=[(S, "rx1"), (S, "rx2")])
    pools = [u for u, s in r2u.items() if s == S]
    assert len(pools) == 1
    p = pools[0]
    members = [e for e in data if e["edge_type"] == "set_member"]
    assert len(members) == 3 and {e["target_id"] for e in members} == {p}
    assert {r2u[e["source_id"]] for e in members} == {"R-HSA-M1", "R-HSA-M2", "R-HSA-M3"}
    assert all(e["and_or"] == "or" and e["pos_neg"] == "pos" for e in members)
    into = [e for e in data if e["source_id"] == p]
    assert sorted(e["target_id"] for e in into) == ["rx1", "rx2"]
    assert all(e["edge_type"] == "catalyst" and e["and_or"] == "and" and e["pos_neg"] == "pos"
               for e in into)
    # reachability is what it was: every member still reaches every copy
    for e in members:
        assert {"rx1", "rx2"} <= reach(data, e["source_id"])


def test_negative_regulator_set_is_one_inhibitor_term(stub, monkeypatch):
    monkeypatch.setenv("LNG_SET_POOL", "1")
    data, r2u = run(neg=[(S, "rx1")])
    p = next(u for u, s in r2u.items() if s == S)
    into = [e for e in data if e["target_id"] == "rx1"]
    assert len(into) == 1 and into[0]["source_id"] == p
    assert into[0]["pos_neg"] == "neg" and into[0]["and_or"] == "or" and into[0]["edge_type"] == "regulator"


def test_one_pool_per_set_across_roles_and_reactions(stub, monkeypatch):
    monkeypatch.setenv("LNG_SET_POOL", "1")
    data, r2u = run(cat=[(S, "rx1")], neg=[(S, "rx2")], pos=[(S, "rx3")])
    assert len([u for u, s in r2u.items() if s == S]) == 1
    assert len([e for e in data if e["edge_type"] == "set_member"]) == 3   # emitted once


def test_a_complex_and_a_modifier_set_are_not_pooled(stub, monkeypatch):
    monkeypatch.setenv("LNG_SET_POOL", "1")
    data, r2u = run(cat=[(K, "rx1")])
    assert not any(e["edge_type"] == "set_member" for e in data)
    assert len(data) == 2                                    # complex members as before
    monkeypatch.setattr(lng, "modifier_isoform_set_ids", lambda: {S})
    data, r2u = run(cat=[(S, "rx1")])
    assert not any(e["edge_type"] == "set_member" for e in data)


def test_a_member_already_in_the_network_is_reused(stub, monkeypatch):
    monkeypatch.setenv("LNG_SET_POOL", "1")
    data: List[Dict[str, Any]] = []
    r2u: Dict[str, str] = {}
    registry = {("R-HSA-M2", "some-rxn", "input"): "u-m2-existing"}
    append_regulators(pd.DataFrame([row(S, "rx1")]), pd.DataFrame(), pd.DataFrame(), data, r2u, registry)
    assert any(e["source_id"] == "u-m2-existing" and e["edge_type"] == "set_member" for e in data)
