"""Tests for transformation semantics.

Verify that edges correctly represent biochemical transformations:
- Edges connect inputs to outputs within reactions
- Multiple inputs × multiple outputs = cartesian product
- Transformations flow in the correct direction
"""

import pytest
import pandas as pd
from typing import Dict, List, Any
import sys
sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')
from src.logic_network_generator import extract_inputs_and_outputs


class TestTransformationSemantics:
    """Test that edges correctly represent biochemical transformations."""

    def test_single_input_single_output_creates_one_edge(self):
        """Reaction: A → B should create exactly one edge A→B."""
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "input-hash",
            "output_hash": "output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # Input: MolA
            {"uid": "output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # Output: MolB
        ])

        uid_reaction_connections = pd.DataFrame([
            {"preceding_uid": "r1-uuid", "following_uid": "r1-uuid"}  # Self-loop
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

        assert len(pathway_logic_network_data) == 1, "Should create exactly one edge"

        edge = pathway_logic_network_data[0]
        entity_a_uuid = reactome_id_to_uuid[1001]
        entity_b_uuid = reactome_id_to_uuid[1002]

        assert edge['source_id'] == entity_a_uuid, "Source should be input physical entity A"
        assert edge['target_id'] == entity_b_uuid, "Target should be output physical entity B"

    def test_two_inputs_one_output_creates_two_edges(self):
        """Reaction: A + B → C should create edges A→C and B→C."""
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "input-hash",
            "output_hash": "output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # Input: MolA
            {"uid": "input-hash", "reactome_id": 100, "component_id": 1,
             "component_id_or_reference_entity_id": 1, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # Input: MolB
            {"uid": "output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1003},  # Output: MolC
        ])

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

        assert len(pathway_logic_network_data) == 2, "Should create 2 edges (A→C, B→C)"

        entity_a_uuid = reactome_id_to_uuid[1001]
        entity_b_uuid = reactome_id_to_uuid[1002]
        entity_c_uuid = reactome_id_to_uuid[1003]

        sources = {edge['source_id'] for edge in pathway_logic_network_data}
        targets = {edge['target_id'] for edge in pathway_logic_network_data}

        assert sources == {entity_a_uuid, entity_b_uuid}, "Sources should be A and B"
        assert targets == {entity_c_uuid}, "All targets should be C"

    def test_one_input_two_outputs_creates_two_edges(self):
        """Reaction: A → B + C should create edges A→B and A→C."""
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "input-hash",
            "output_hash": "output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # Input: MolA
            {"uid": "output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # Output: MolB
            {"uid": "output-hash", "reactome_id": 100, "component_id": 1,
             "component_id_or_reference_entity_id": 1, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1003},  # Output: MolC
        ])

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

        assert len(pathway_logic_network_data) == 2, "Should create 2 edges (A→B, A→C)"

        entity_a_uuid = reactome_id_to_uuid[1001]
        entity_b_uuid = reactome_id_to_uuid[1002]
        entity_c_uuid = reactome_id_to_uuid[1003]

        sources = {edge['source_id'] for edge in pathway_logic_network_data}
        targets = {edge['target_id'] for edge in pathway_logic_network_data}

        assert sources == {entity_a_uuid}, "All sources should be A"
        assert targets == {entity_b_uuid, entity_c_uuid}, "Targets should be B and C"

    def test_two_inputs_two_outputs_cartesian_product(self):
        """Reaction: A + B → C + D should create 4 edges (cartesian product).

        Edges: A→C, A→D, B→C, B→D
        """
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "input-hash",
            "output_hash": "output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            # Inputs: A, B
            {"uid": "input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # MolA
            {"uid": "input-hash", "reactome_id": 100, "component_id": 1,
             "component_id_or_reference_entity_id": 1, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # MolB
            # Outputs: C, D
            {"uid": "output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1003},  # MolC
            {"uid": "output-hash", "reactome_id": 100, "component_id": 1,
             "component_id_or_reference_entity_id": 1, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1004},  # MolD
        ])

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

        assert len(pathway_logic_network_data) == 4, "Should create 4 edges (2×2 cartesian product)"

        entity_a_uuid = reactome_id_to_uuid[1001]
        entity_b_uuid = reactome_id_to_uuid[1002]
        entity_c_uuid = reactome_id_to_uuid[1003]
        entity_d_uuid = reactome_id_to_uuid[1004]

        # Check that all 4 combinations exist
        edge_pairs = {(edge['source_id'], edge['target_id']) for edge in pathway_logic_network_data}
        expected = {
            (entity_a_uuid, entity_c_uuid),  # A→C
            (entity_a_uuid, entity_d_uuid),  # A→D
            (entity_b_uuid, entity_c_uuid),  # B→C
            (entity_b_uuid, entity_d_uuid),  # B→D
        }

        assert edge_pairs == expected, f"Expected all 4 combinations, got {edge_pairs}"

    def test_transformation_direction_input_to_output(self):
        """Verify edges always flow from inputs to outputs (not backwards)."""
        reaction_id_map = pd.DataFrame([{
            "uid": "r1-uuid",
            "reactome_id": 100,
            "input_hash": "input-hash",
            "output_hash": "output-hash",
        }])

        decomposed_uid_mapping = pd.DataFrame([
            {"uid": "input-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},  # Input
            {"uid": "output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},  # Output
        ])

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

        edge = pathway_logic_network_data[0]
        input_uuid = reactome_id_to_uuid[1001]
        output_uuid = reactome_id_to_uuid[1002]

        # Critical assertion: verify direction
        assert edge['source_id'] == input_uuid, "Source must be INPUT physical entity (reactant)"
        assert edge['target_id'] == output_uuid, "Target must be OUTPUT physical entity (product)"
        assert edge['source_id'] != edge['target_id'], "Should not be a self-loop"
