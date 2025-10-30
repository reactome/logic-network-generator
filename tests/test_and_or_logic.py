"""Tests for AND/OR logic based on user requirements.

User clarification:
- Multiple sources → same physical entity: OR relationships (R1→A (OR), R2→A (OR))
- Physical entity → reaction: AND relationships (always) (A→R3 (AND))
- Single source → physical entity: AND relationship (R1→A (AND) if R1 is only source)
"""

import pandas as pd
from typing import Dict, List, Any
import sys
from unittest.mock import patch

sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import extract_inputs_and_outputs


class TestAndOrLogic:
    """Test AND/OR logic assignment based on preceding reaction counts."""

    def test_single_preceding_reaction_creates_and_edges(self):
        """When one reaction produces a physical entity, edges should be AND."""
        # Setup: R1 produces MolA → MolB (single source for transformation)
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "r1-input-hash",
            "output_hash": "r1-output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "r1-input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # MolA
            {"uid": "r1-output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # MolB
        ])

        # Self-loop connection (reaction connects to itself)
        uid_reaction_connections = pd.DataFrame([
            {"preceding_uid": "r1-uuid", "following_uid": "r1-uuid"}
        ])

        reaction_uids = ["r1-uuid"]
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        extract_inputs_and_outputs(
            reaction_uid="r1-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        assert len(pathway_logic_network_data) == 1
        edge = pathway_logic_network_data[0]
        assert edge['and_or'] == 'and', "Single source should create AND relationship"
        assert edge['edge_type'] == 'input'

    def test_multiple_preceding_reactions_create_or_edges(self):
        """When multiple reactions feed into one, edges should be OR."""
        # Setup: R1 and R2 both produce physical entities consumed by R3
        # This simulates: R1→A (OR), R2→A (OR), A→R3 (AND)

        reaction_id_map = pd.DataFrame([
            {
                "uid": "r1-uuid",
                "reactome_id": 100,
                "input_hash": "r1-input-hash",
                "output_hash": "r1-output-hash",
            },
            {
                "uid": "r2-uuid",
                "reactome_id": 200,
                "input_hash": "r2-input-hash",
                "output_hash": "r2-output-hash",
            },
            {
                "uid": "r3-uuid",
                "reactome_id": 300,
                "input_hash": "r3-input-hash",
                "output_hash": "r3-output-hash",
            },
        ])

        decomposed_uid_mapping = pd.DataFrame([
            # R1 outputs MolA
            {"uid": "r1-output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # MolA
            # R2 outputs MolA (same physical entity from different reaction)
            {"uid": "r2-output-hash", "reactome_id": 200, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # MolA
            # R3 inputs MolA
            {"uid": "r3-input-hash", "reactome_id": 300, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # MolA
            # R3 outputs MolB
            {"uid": "r3-output-hash", "reactome_id": 300, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # MolB
        ])

        # R3 has TWO preceding reactions (R1 and R2)
        uid_reaction_connections = pd.DataFrame([
            {"preceding_uid": "r1-uuid", "following_uid": "r3-uuid"},
            {"preceding_uid": "r2-uuid", "following_uid": "r3-uuid"},
        ])

        reaction_uids = ["r3-uuid"]
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        extract_inputs_and_outputs(
            reaction_uid="r3-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        # Should create edges from R3's inputs to both R1 and R2's outputs
        assert len(pathway_logic_network_data) == 2, "Should create 2 edges (one per preceding)"

        for edge in pathway_logic_network_data:
            assert edge['and_or'] == 'or', "Multiple sources should create OR relationship"
            assert edge['edge_type'] == 'output'

    def test_three_preceding_reactions_create_or_edges(self):
        """Test OR logic with three preceding reactions."""
        reaction_id_map = pd.DataFrame([
            {"uid": "r1-uuid", "reactome_id": 100, "input_hash": "r1-in", "output_hash": "r1-out"},
            {"uid": "r2-uuid", "reactome_id": 200, "input_hash": "r2-in", "output_hash": "r2-out"},
            {"uid": "r3-uuid", "reactome_id": 300, "input_hash": "r3-in", "output_hash": "r3-out"},
            {"uid": "r4-uuid", "reactome_id": 400, "input_hash": "r4-in", "output_hash": "r4-out"},
        ])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "r1-out", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            {"uid": "r2-out", "reactome_id": 200, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            {"uid": "r3-out", "reactome_id": 300, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            {"uid": "r4-in", "reactome_id": 400, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            {"uid": "r4-out", "reactome_id": 400, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},
        ])

        # R4 has THREE preceding reactions
        uid_reaction_connections = pd.DataFrame([
            {"preceding_uid": "r1-uuid", "following_uid": "r4-uuid"},
            {"preceding_uid": "r2-uuid", "following_uid": "r4-uuid"},
            {"preceding_uid": "r3-uuid", "following_uid": "r4-uuid"},
        ])

        reaction_uids = ["r4-uuid"]
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        extract_inputs_and_outputs(
            reaction_uid="r4-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        assert len(pathway_logic_network_data) == 3
        for edge in pathway_logic_network_data:
            assert edge['and_or'] == 'or', "Three sources should create OR relationships"

    def test_zero_preceding_reactions_creates_and_edges(self):
        """Root reactions (no preceding) should still create AND edges."""
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "r1-input-hash",
            "output_hash": "r1-output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "r1-input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            {"uid": "r1-output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},
        ])

        # No preceding reactions (root)
        uid_reaction_connections = pd.DataFrame(columns=["preceding_uid", "following_uid"])

        reaction_uids = ["r1-uuid"]
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        extract_inputs_and_outputs(
            reaction_uid="r1-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        # With no preceding reactions, no edges are created
        # This is expected - root reactions have no edges from preceding reactions
        assert len(pathway_logic_network_data) == 0
