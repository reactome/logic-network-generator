"""Unit tests for the schema-backed provenance exports (nodes.csv,
node_reaction_context.csv). No Neo4j: entity labels / terminal components are
mocked. See schema/logic_network.linkml.yaml."""
import pandas as pd

import src.logic_network_generator as m
from src import neo4j_connector


def test_export_nodes_classifies_kinds(tmp_path, monkeypatch):
    # simple entity labels; variant set-derivation stubbed out (no Neo4j)
    monkeypatch.setattr(neo4j_connector, "get_labels",
                        lambda e: ["EntityWithAccessionedSequence"])
    monkeypatch.setattr(m, "_derive_sets_and_chosen",
                        lambda parent, members: (["R-HSA-75202"], ["R-HSA-68891"]))
    monkeypatch.setattr(m, "get_terminal_components", lambda s: {s})

    variant = "R-HSA-141608::variant::R-HSA-68365_R-HSA-68891"
    # Realistic UUIDs (export_nodes detects mapping direction by UUID shape).
    s1 = "aaaaaaaa-0000-0000-0000-000000000001"
    v1 = "aaaaaaaa-0000-0000-0000-000000000002"
    d1 = "aaaaaaaa-0000-0000-0000-000000000003"
    rxn1 = "aaaaaaaa-0000-0000-0000-0000000000r1"
    edges = pd.DataFrame([
        {"source_id": s1, "target_id": rxn1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": rxn1, "target_id": v1, "pos_neg": "pos", "and_or": None,
         "edge_type": "output", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": v1, "target_id": d1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "dissociation", "stoichiometry": 1, "edge_reaction_id": None},
    ])
    reaction_id_map = pd.DataFrame({"uid": [rxn1], "reactome_id": ["R-HSA-100"]})
    uuid_mapping = {s1: "R-HSA-999", v1: variant, d1: "R-HSA-888"}

    out = tmp_path / "nodes.csv"
    m.export_nodes(edges, reaction_id_map, uuid_mapping, str(out))
    rows = {r["uuid"]: r for r in pd.read_csv(out, dtype=str, keep_default_na=False)
            .to_dict("records")}
    assert rows[rxn1]["node_kind"] == "reaction"
    assert rows[rxn1]["diagram_entity_id"] == "R-HSA-100"
    assert rows[v1]["node_kind"] == "set_variant"
    assert rows[v1]["diagram_entity_id"] == "R-HSA-141608"
    assert "R-HSA-68891" in rows[v1]["member_leaves"]
    assert rows[s1]["node_kind"] == "simple_entity"
    assert rows[d1]["node_kind"] == "dissociation_sink"


def test_export_node_reaction_context(tmp_path):
    entity_uuid_registry = {
        ("R-HSA-999", "rxn1", "input"): "s1",
        ("R-HSA-141608::variant::x", "rxn1", "output"): "v1",
    }
    reaction_id_map = pd.DataFrame({"uid": ["rxn1"], "reactome_id": ["R-HSA-100"]})
    catreg = pd.DataFrame([{"reaction_id": "R-HSA-100", "entity_id": "R-HSA-7",
                            "edge_type": "catalyst", "uuid": "c1", "reaction_uuid": "rxn1"}])
    out = tmp_path / "ctx.csv"
    m.export_node_reaction_context(entity_uuid_registry, reaction_id_map, catreg, str(out))
    ctx = pd.read_csv(out).to_dict("records")
    triples = {(r["context_node"], r["reaction_id"], r["role"]) for r in ctx}
    assert ("s1", "R-HSA-100", "input") in triples
    assert ("v1", "R-HSA-100", "output") in triples
    assert ("c1", "R-HSA-100", "catalyst") in triples


def test_parse_variant_members():
    assert m._parse_variant_members("R-HSA-1::variant::R-HSA-2_R-HSA-3") == {"R-HSA-2", "R-HSA-3"}
    assert m._parse_variant_members("R-HSA-1") == set()
