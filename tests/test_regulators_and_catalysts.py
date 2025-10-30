"""Tests for regulator and catalyst functionality.

These tests verify that:
1. Negative regulators are correctly marked with pos_neg = "neg"
2. Positive regulators are correctly marked with pos_neg = "pos"
3. Catalysts are correctly marked with pos_neg = "pos"
4. Regulatory edges have correct edge_type values
5. Regulatory relationships are properly created
"""

import pytest
import pandas as pd
from typing import Dict, List, Any
import sys
from unittest.mock import Mock, patch

sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import append_regulators


class TestRegulatorsAndCatalysts:
    """Test regulatory and catalytic relationships in logic networks."""

    def test_negative_regulators_have_neg_pos_neg(self):
        """Negative regulators should have pos_neg = 'neg'."""
        # Create mock regulator data
        negative_regulator_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "regulator",
             "uuid": "neg-regulator-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": 101, "catalyst_id": 201, "edge_type": "regulator",
             "uuid": "neg-regulator-2", "reaction_uuid": "reaction-2"},
        ])

        catalyst_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Append regulators
        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        # Verify all negative regulator edges have pos_neg = "neg"
        assert len(pathway_logic_network_data) == 2, "Should create 2 negative regulator edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'neg', f"Negative regulator should have pos_neg='neg', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'regulator', f"Should have edge_type='regulator', got '{edge['edge_type']}'"
            assert edge['source_id'] in ['neg-regulator-1', 'neg-regulator-2'], "Source should be negative regulator UUID"

    def test_positive_regulators_have_pos_pos_neg(self):
        """Positive regulators should have pos_neg = 'pos'."""
        # Create mock regulator data
        positive_regulator_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "regulator",
             "uuid": "pos-regulator-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": 101, "catalyst_id": 201, "edge_type": "regulator",
             "uuid": "pos-regulator-2", "reaction_uuid": "reaction-2"},
        ])

        catalyst_map = pd.DataFrame()
        negative_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Append regulators
        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        # Verify all positive regulator edges have pos_neg = "pos"
        assert len(pathway_logic_network_data) == 2, "Should create 2 positive regulator edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'pos', f"Positive regulator should have pos_neg='pos', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'regulator', f"Should have edge_type='regulator', got '{edge['edge_type']}'"

    def test_catalysts_have_pos_pos_neg(self):
        """Catalysts should have pos_neg = 'pos' and edge_type = 'catalyst'."""
        # Create mock catalyst data
        catalyst_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": 101, "catalyst_id": 201, "edge_type": "catalyst",
             "uuid": "catalyst-2", "reaction_uuid": "reaction-2"},
        ])

        negative_regulator_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Append regulators
        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        # Verify all catalyst edges have correct properties
        assert len(pathway_logic_network_data) == 2, "Should create 2 catalyst edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'pos', f"Catalyst should have pos_neg='pos', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'catalyst', f"Should have edge_type='catalyst', got '{edge['edge_type']}'"

    def test_mixed_regulators_and_catalysts(self):
        """Test that mixed regulators and catalysts are all correctly marked."""
        # Create mock data with all three types
        catalyst_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame([
            {"reaction_id": 101, "catalyst_id": 201, "edge_type": "regulator",
             "uuid": "neg-reg-1", "reaction_uuid": "reaction-2"},
        ])

        positive_regulator_map = pd.DataFrame([
            {"reaction_id": 102, "catalyst_id": 202, "edge_type": "regulator",
             "uuid": "pos-reg-1", "reaction_uuid": "reaction-3"},
        ])

        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Append all regulators
        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        # Verify we have all three edges
        assert len(pathway_logic_network_data) == 3, "Should create 3 edges total"

        # Separate edges by type
        catalyst_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'catalyst']
        regulator_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'regulator']

        # Verify counts
        assert len(catalyst_edges) == 1, "Should have 1 catalyst edge"
        assert len(regulator_edges) == 2, "Should have 2 regulator edges"

        # Verify catalyst properties
        assert catalyst_edges[0]['pos_neg'] == 'pos', "Catalyst should be positive"

        # Verify regulator properties
        negative_edges = [e for e in regulator_edges if e['pos_neg'] == 'neg']
        positive_edges = [e for e in regulator_edges if e['pos_neg'] == 'pos']

        assert len(negative_edges) == 1, "Should have 1 negative regulator"
        assert len(positive_edges) == 1, "Should have 1 positive regulator"

    def test_regulator_edges_point_to_reactions(self):
        """Regulator and catalyst edges should point to reaction UUIDs as targets."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "catalyst",
             "uuid": "catalyst-uuid-1", "reaction_uuid": "reaction-uuid-1"},
        ])

        negative_regulator_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        # Verify edge structure
        edge = pathway_logic_network_data[0]
        assert edge['source_id'] == 'catalyst-uuid-1', "Source should be catalyst UUID"
        assert edge['target_id'] == 'reaction-uuid-1', "Target should be reaction UUID"

    def test_regulators_have_empty_and_or_logic(self):
        """Regulators and catalysts should have empty AND/OR logic (not transformations)."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": 100, "catalyst_id": 200, "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame([
            {"reaction_id": 101, "catalyst_id": 201, "edge_type": "regulator",
             "uuid": "neg-reg-1", "reaction_uuid": "reaction-2"},
        ])

        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Append with empty and_or
        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",  # Should be empty for regulators
            edge_type=""
        )

        # Verify all edges have empty and_or
        for edge in pathway_logic_network_data:
            assert edge['and_or'] == "", f"Regulator/catalyst should have empty and_or, got '{edge['and_or']}'"

    def test_empty_regulator_maps_create_no_edges(self):
        """Empty regulator dataframes should not create any edges."""
        catalyst_map = pd.DataFrame()
        negative_regulator_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            and_or="",
            edge_type=""
        )

        assert len(pathway_logic_network_data) == 0, "Empty regulator maps should create no edges"


class TestRealNetworkRegulators:
    """Test regulators in actual generated networks (if available)."""

    @pytest.mark.skipif(
        not pd.io.common.file_exists('pathway_logic_network_69620.csv'),
        reason="Real network file not available"
    )
    def test_real_network_has_negative_regulators(self):
        """If real network exists, verify it has properly marked negative regulators."""
        network = pd.read_csv('pathway_logic_network_69620.csv')

        # Get all regulatory edges
        regulator_edges = network[network['edge_type'] == 'regulator']

        if len(regulator_edges) > 0:
            # Check for negative regulators
            negative_regulators = regulator_edges[regulator_edges['pos_neg'] == 'neg']
            positive_regulators = regulator_edges[regulator_edges['pos_neg'] == 'pos']

            print(f"\nRegulator statistics:")
            print(f"  Total regulators: {len(regulator_edges)}")
            print(f"  Negative regulators: {len(negative_regulators)}")
            print(f"  Positive regulators: {len(positive_regulators)}")

            # All regulators should be either positive or negative
            assert len(negative_regulators) + len(positive_regulators) == len(regulator_edges), \
                "All regulators should be marked as either positive or negative"

    @pytest.mark.skipif(
        not pd.io.common.file_exists('pathway_logic_network_69620.csv'),
        reason="Real network file not available"
    )
    def test_real_network_catalysts_are_positive(self):
        """If real network exists, verify all catalysts are positive."""
        network = pd.read_csv('pathway_logic_network_69620.csv')

        catalyst_edges = network[network['edge_type'] == 'catalyst']

        if len(catalyst_edges) > 0:
            # All catalysts should be positive
            negative_catalysts = catalyst_edges[catalyst_edges['pos_neg'] == 'neg']

            assert len(negative_catalysts) == 0, \
                f"Found {len(negative_catalysts)} negative catalysts - catalysts should always be positive"

            print(f"\nCatalyst statistics:")
            print(f"  Total catalysts: {len(catalyst_edges)}")
            print(f"  All catalysts are positive ✓")
