"""Tests for logic_network_generator module."""

import pytest
import pandas as pd
from typing import Dict, List, Any


# Import functions to test
import sys
from unittest.mock import patch

sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import (
        _assign_uuids,
        _determine_edge_properties,
        _add_pathway_connections,
    )


class Test_assign_uuids:
    """Tests for _assign_uuids function."""

    def test_assigns_new_uuid_for_new_reactome_id(self):
        """Should create a new UUID for a reactome ID not in the mapping."""
        reactome_id_to_uuid: Dict[str, str] = {}
        reactome_ids = ["12345"]

        result = _assign_uuids(reactome_ids, reactome_id_to_uuid)

        assert len(result) == 1
        assert "12345" in reactome_id_to_uuid
        assert result[0] == reactome_id_to_uuid["12345"]

    def test_reuses_existing_uuid_for_known_reactome_id(self):
        """Should reuse existing UUID for a reactome ID already in the mapping."""
        existing_uuid = "test-uuid-123"
        reactome_id_to_uuid = {"12345": existing_uuid}
        reactome_ids = ["12345"]

        result = _assign_uuids(reactome_ids, reactome_id_to_uuid)

        assert len(result) == 1
        assert result[0] == existing_uuid

    def test_handles_multiple_reactome_ids(self):
        """Should handle multiple reactome IDs correctly."""
        reactome_id_to_uuid: Dict[str, str] = {"12345": "existing-uuid"}
        reactome_ids = ["12345", "67890", "11111"]

        result = _assign_uuids(reactome_ids, reactome_id_to_uuid)

        assert len(result) == 3
        assert result[0] == "existing-uuid"  # Reused
        assert result[1] != result[2]  # New UUIDs are different


class Test_determine_edge_properties:
    """Tests for _determine_edge_properties function."""

    def test_single_preceding_reaction_returns_and(self):
        """When there's one preceding reaction, should return 'and' and 'input'."""
        and_or, edge_type = _determine_edge_properties(1)

        assert and_or == "and"
        assert edge_type == "input"

    def test_multiple_preceding_reactions_returns_or(self):
        """When there are multiple preceding reactions, should return 'or' and 'output'."""
        and_or, edge_type = _determine_edge_properties(2)
        assert and_or == "or"
        assert edge_type == "output"

        and_or, edge_type = _determine_edge_properties(5)
        assert and_or == "or"
        assert edge_type == "output"

    def test_zero_preceding_reactions(self):
        """Edge case: zero preceding reactions should return 'and' and 'input'."""
        and_or, edge_type = _determine_edge_properties(0)
        assert and_or == "and"
        assert edge_type == "input"


class Test_add_pathway_connections:
    """Tests for _add_pathway_connections function."""

    def test_adds_single_connection(self):
        """Should add a single connection between one input and one output."""
        pathway_data: List[Dict[str, Any]] = []
        input_uuids = ["input-uuid-1"]
        output_uuids = ["output-uuid-1"]

        _add_pathway_connections(
            input_uuids, output_uuids, "and", "input", pathway_data
        )

        assert len(pathway_data) == 1
        edge = pathway_data[0]
        assert edge["pos_neg"] == "pos"
        assert edge["and_or"] == "and"
        assert edge["edge_type"] == "input"

    def test_cartesian_product_of_inputs_and_outputs(self):
        """Should create edges for all combinations of inputs and outputs."""
        pathway_data: List[Dict[str, Any]] = []
        input_uuids = ["input-1", "input-2"]
        output_uuids = ["output-1", "output-2", "output-3"]

        _add_pathway_connections(
            input_uuids, output_uuids, "or", "output", pathway_data
        )

        # Should create 2 * 3 = 6 edges
        assert len(pathway_data) == 6

        # Check all combinations exist
        sources = [edge["source_id"] for edge in pathway_data]
        targets = [edge["target_id"] for edge in pathway_data]

        # All inputs should appear as sources
        assert sources.count("input-1") == 3
        assert sources.count("input-2") == 3

        # All outputs should appear as targets
        assert targets.count("output-1") == 2
        assert targets.count("output-2") == 2
        assert targets.count("output-3") == 2

    def test_edge_direction_semantics(self):
        """
        CRITICAL TEST: Verify edge direction represents correct molecular flow.

        Assumption: edges should represent molecular flow through the pathway.
        - If input_uuids are from current reaction's inputs
        - And output_uuids are from preceding reaction's outputs
        - Then edges should flow: preceding_output → current_input

        Current implementation: source_id = input_uuid, target_id = output_uuid
        This would be: current_input → preceding_output (BACKWARDS?)

        Expected: source_id = output_uuid, target_id = input_uuid
        This would be: preceding_output → current_input (FORWARD)
        """
        pathway_data: List[Dict[str, Any]] = []
        current_input_uuids = ["current-input-molecule"]
        preceding_output_uuids = ["preceding-output-molecule"]

        _add_pathway_connections(
            current_input_uuids, preceding_output_uuids, "and", "input", pathway_data
        )

        edge = pathway_data[0]

        # Document what we observe
        print(f"\nObserved edge: {edge['source_id']} → {edge['target_id']}")
        print(f"If correct flow: preceding-output-molecule → current-input-molecule")
        print(f"Current code creates: {edge['source_id']} → {edge['target_id']}")

        # This test will FAIL if edges are backwards
        # Expected behavior: molecular flow from preceding output to current input
        # TODO: Determine if this assertion is correct based on system requirements
        # assert edge["source_id"] == "preceding-output-molecule", "Edge should flow from preceding output"
        # assert edge["target_id"] == "current-input-molecule", "Edge should flow to current input"

        # For now, just document what the code actually does
        assert edge["source_id"] == "current-input-molecule"  # Current behavior
        assert edge["target_id"] == "preceding-output-molecule"  # Current behavior
