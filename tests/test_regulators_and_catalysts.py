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
from pathlib import Path
from unittest.mock import patch

# Add project root to Python path dynamically
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import append_regulators


def _mock_decompose(entity_id):
    """Return entity as-is (no decomposition) for unit tests."""
    return [(entity_id, 1)]


class TestRegulatorsAndCatalysts:
    """Test regulatory and catalytic relationships in logic networks."""

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_negative_regulators_have_neg_pos_neg(self, mock_decompose):
        """Negative regulators should have pos_neg = 'neg'."""
        negative_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "regulator",
             "uuid": "neg-regulator-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": "R-HSA-101", "entity_id": "R-HSA-201", "edge_type": "regulator",
             "uuid": "neg-regulator-2", "reaction_uuid": "reaction-2"},
        ])

        catalyst_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        assert len(pathway_logic_network_data) == 2, "Should create 2 negative regulator edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'neg', f"Negative regulator should have pos_neg='neg', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'regulator', f"Should have edge_type='regulator', got '{edge['edge_type']}'"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_positive_regulators_have_pos_pos_neg(self, mock_decompose):
        """Positive regulators should have pos_neg = 'pos'."""
        positive_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "regulator",
             "uuid": "pos-regulator-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": "R-HSA-101", "entity_id": "R-HSA-201", "edge_type": "regulator",
             "uuid": "pos-regulator-2", "reaction_uuid": "reaction-2"},
        ])

        catalyst_map = pd.DataFrame()
        negative_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        assert len(pathway_logic_network_data) == 2, "Should create 2 positive regulator edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'pos', f"Positive regulator should have pos_neg='pos', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'regulator', f"Should have edge_type='regulator', got '{edge['edge_type']}'"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_catalysts_have_pos_pos_neg(self, mock_decompose):
        """Catalysts should have pos_neg = 'pos' and edge_type = 'catalyst'."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
            {"reaction_id": "R-HSA-101", "entity_id": "R-HSA-201", "edge_type": "catalyst",
             "uuid": "catalyst-2", "reaction_uuid": "reaction-2"},
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
        )

        assert len(pathway_logic_network_data) == 2, "Should create 2 catalyst edges"

        for edge in pathway_logic_network_data:
            assert edge['pos_neg'] == 'pos', f"Catalyst should have pos_neg='pos', got '{edge['pos_neg']}'"
            assert edge['edge_type'] == 'catalyst', f"Should have edge_type='catalyst', got '{edge['edge_type']}'"
            assert edge['and_or'] == 'and', f"Catalyst should have and_or='and', got '{edge['and_or']}'"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_mixed_regulators_and_catalysts(self, mock_decompose):
        """Test that mixed regulators and catalysts are all correctly marked."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-101", "entity_id": "R-HSA-201", "edge_type": "regulator",
             "uuid": "neg-reg-1", "reaction_uuid": "reaction-2"},
        ])

        positive_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-102", "entity_id": "R-HSA-202", "edge_type": "regulator",
             "uuid": "pos-reg-1", "reaction_uuid": "reaction-3"},
        ])

        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        assert len(pathway_logic_network_data) == 3, "Should create 3 edges total"

        catalyst_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'catalyst']
        regulator_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'regulator']

        assert len(catalyst_edges) == 1, "Should have 1 catalyst edge"
        assert len(regulator_edges) == 2, "Should have 2 regulator edges"

        assert catalyst_edges[0]['pos_neg'] == 'pos', "Catalyst should be positive"

        negative_edges = [e for e in regulator_edges if e['pos_neg'] == 'neg']
        positive_edges = [e for e in regulator_edges if e['pos_neg'] == 'pos']

        assert len(negative_edges) == 1, "Should have 1 negative regulator"
        assert len(positive_edges) == 1, "Should have 1 positive regulator"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_regulator_edges_point_to_reactions(self, mock_decompose):
        """Regulator and catalyst edges should point to reaction UUIDs as targets."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
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
        )

        edge = pathway_logic_network_data[0]
        assert edge['target_id'] == 'reaction-uuid-1', "Target should be reaction UUID"
        # source_id is now a new UUID (from decomposition), verify it maps back
        assert reactome_id_to_uuid[edge['source_id']] == 'R-HSA-200', \
            "Source UUID should map back to entity stId"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_and_or_logic_per_type(self, mock_decompose):
        """Catalysts and regulators should both propagate AND/OR from decomposition."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-101", "entity_id": "R-HSA-201", "edge_type": "regulator",
             "uuid": "neg-reg-1", "reaction_uuid": "reaction-2"},
        ])

        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        catalyst_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'catalyst']
        regulator_edges = [e for e in pathway_logic_network_data if e['edge_type'] == 'regulator']

        # Catalysts (pos) → and; negative regulators (neg) → or
        for edge in catalyst_edges:
            assert edge['and_or'] == "and", f"Catalyst should have and_or='and', got '{edge['and_or']}'"
        for edge in regulator_edges:
            assert edge['pos_neg'] == "neg"
            assert edge['and_or'] == "or", (
                f"Negative regulator should have and_or='or', got '{edge['and_or']}'"
            )

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_empty_regulator_maps_create_no_edges(self, mock_decompose):
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
        )

        assert len(pathway_logic_network_data) == 0, "Empty regulator maps should create no edges"

    @patch('src.logic_network_generator._decompose_regulator_entity')
    def test_complex_catalyst_decomposed_to_and_members(self, mock_decompose):
        """Complex catalysts should be decomposed into AND members."""
        mock_decompose.return_value = [
            ("R-HSA-301", 1),
            ("R-HSA-302", 1),
            ("R-HSA-303", 1),
        ]

        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-300", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
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
        )

        assert len(pathway_logic_network_data) == 3, "Complex with 3 components should create 3 edges"

        for edge in pathway_logic_network_data:
            assert edge['edge_type'] == 'catalyst'
            assert edge['pos_neg'] == 'pos'
            assert edge['and_or'] == 'and', "Complex members should have AND logic"
            assert edge['target_id'] == 'reaction-1'

        # Verify all decomposed members are in the UUID mapping
        mapped_stids = set(reactome_id_to_uuid.values())
        assert mapped_stids == {"R-HSA-301", "R-HSA-302", "R-HSA-303"}

    @patch('src.logic_network_generator._decompose_regulator_entity')
    def test_entityset_catalyst_emits_one_edge_per_member(self, mock_decompose):
        """EntitySet catalysts should emit one edge per decomposed member.

        and_or reflects reaction-level requirement (all catalysts are
        required → 'and'), not within-entity decomposition.
        """
        mock_decompose.return_value = [
            ("R-HSA-401", 1),
            ("R-HSA-402", 1),
        ]

        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-400", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
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
        )

        assert len(pathway_logic_network_data) == 2, "EntitySet with 2 members should create 2 edges"

        for edge in pathway_logic_network_data:
            assert edge['edge_type'] == 'catalyst'
            assert edge['pos_neg'] == 'pos'
            assert edge['and_or'] == 'and', "Catalyst edges are reaction-required → and"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_stoichiometry_defaults_to_one(self, mock_decompose):
        """Edges should have stoichiometry=1 by default."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
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
        )

        assert len(pathway_logic_network_data) == 1
        assert pathway_logic_network_data[0]['stoichiometry'] == 1

    @patch('src.logic_network_generator._decompose_regulator_entity')
    def test_nested_complex_stoichiometry_multiplication(self, mock_decompose):
        """Nested Complex with stoichiometry: Complex with 2x SubComplex that has 3x Protein -> stoichiometry 6."""
        mock_decompose.return_value = [
            ("R-HSA-PROTEIN", 6),  # 2 * 3 = 6
        ]

        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-OUTER-COMPLEX", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
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
        )

        assert len(pathway_logic_network_data) == 1
        edge = pathway_logic_network_data[0]
        assert edge['stoichiometry'] == 6, f"Expected stoichiometry 6 (2*3), got {edge['stoichiometry']}"
        assert edge['edge_type'] == 'catalyst'
        assert edge['and_or'] == 'and'

    @patch('src.logic_network_generator._decompose_regulator_entity')
    def test_complex_with_mixed_stoichiometry(self, mock_decompose):
        """Complex with components having different stoichiometries."""
        mock_decompose.return_value = [
            ("R-HSA-A", 2),
            ("R-HSA-B", 1),
            ("R-HSA-C", 3),
        ]

        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-COMPLEX", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
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
        )

        assert len(pathway_logic_network_data) == 3
        stoichs = [e['stoichiometry'] for e in pathway_logic_network_data]
        assert sorted(stoichs) == [1, 2, 3], f"Expected stoichiometries [1, 2, 3], got {sorted(stoichs)}"


class TestRegulatorUuidReuse:
    """Test that regulators reuse existing pathway UUIDs when available."""

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_regulator_reuses_pathway_uuid(self, mock_decompose):
        """When entity_uuid_registry contains the same stId, its UUID should be reused."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Simulate entity_uuid_registry with R-HSA-200 already registered
        existing_uuid = "existing-uuid-for-200"
        entity_uuid_registry = {
            ("R-HSA-200", "some-vr-uid", "input"): existing_uuid,
        }

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            entity_uuid_registry=entity_uuid_registry,
        )

        assert len(pathway_logic_network_data) == 1
        edge = pathway_logic_network_data[0]
        assert edge['source_id'] == existing_uuid, \
            f"Should reuse existing UUID '{existing_uuid}', got '{edge['source_id']}'"

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_regulator_creates_fresh_uuid_when_no_pathway_match(self, mock_decompose):
        """When entity_uuid_registry has no matching stId, a fresh UUID should be created."""
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-100", "entity_id": "R-HSA-200", "edge_type": "catalyst",
             "uuid": "catalyst-1", "reaction_uuid": "reaction-1"},
        ])

        negative_regulator_map = pd.DataFrame()
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        # Registry with a DIFFERENT entity - no match for R-HSA-200
        entity_uuid_registry = {
            ("R-HSA-999", "some-vr-uid", "input"): "uuid-for-999",
        }

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
            entity_uuid_registry=entity_uuid_registry,
        )

        assert len(pathway_logic_network_data) == 1
        edge = pathway_logic_network_data[0]
        assert edge['source_id'] != "uuid-for-999", \
            "Should NOT reuse UUID from a different entity"
        assert edge['source_id'] != "", "Should have a valid UUID"


class TestRegulatorSharedAcrossReactions:
    """Same regulator entity on multiple reactions must share one UUID.

    Without this, MDM2 regulating both R1 and R2 (and not appearing as a
    regular input/output of either) would get a fresh UUID per regulator-
    row, leaving the network with two disconnected MDM2 nodes that look
    like different proteins to a perturbation tool. The fix updates
    stid_to_existing_uuid as fresh UUIDs are minted, so subsequent
    emissions for the same stId reuse it.
    """

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_same_regulator_on_two_reactions_shares_uuid(self, mock_decompose):
        # MDM2 (R-HSA-MDM2) is a positive regulator of two different reactions.
        # No entity_uuid_registry — it's a pure regulator, not an input/output.
        positive_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-R1", "entity_id": "R-HSA-MDM2",
             "edge_type": "regulator", "uuid": "reg-link-1",
             "reaction_uuid": "vr-1"},
            {"reaction_id": "R-HSA-R2", "entity_id": "R-HSA-MDM2",
             "edge_type": "regulator", "uuid": "reg-link-2",
             "reaction_uuid": "vr-2"},
        ])
        catalyst_map = pd.DataFrame()
        negative_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        # Two regulator edges emitted — one per reaction
        assert len(pathway_logic_network_data) == 2
        # Both must share the same source UUID for MDM2
        sources = {e["source_id"] for e in pathway_logic_network_data}
        assert len(sources) == 1, (
            f"Same regulator on different reactions must share one UUID, "
            f"got {len(sources)} distinct UUIDs: {sources}"
        )
        # And the targets are different reactions
        targets = {e["target_id"] for e in pathway_logic_network_data}
        assert targets == {"vr-1", "vr-2"}

    @patch('src.logic_network_generator._decompose_regulator_entity', side_effect=_mock_decompose)
    def test_catalyst_and_regulator_for_same_protein_share_uuid(self, mock_decompose):
        # Same protein appears as catalyst of R1 and negative regulator of R2.
        # Distinct edges (different edge_type) but the same source node.
        catalyst_map = pd.DataFrame([
            {"reaction_id": "R-HSA-R1", "entity_id": "R-HSA-PROT",
             "edge_type": "catalyst", "uuid": "cat-1", "reaction_uuid": "vr-1"},
        ])
        negative_regulator_map = pd.DataFrame([
            {"reaction_id": "R-HSA-R2", "entity_id": "R-HSA-PROT",
             "edge_type": "regulator", "uuid": "reg-1", "reaction_uuid": "vr-2"},
        ])
        positive_regulator_map = pd.DataFrame()
        pathway_logic_network_data: List[Dict[str, Any]] = []
        reactome_id_to_uuid: Dict[str, str] = {}

        append_regulators(
            catalyst_map,
            negative_regulator_map,
            positive_regulator_map,
            pathway_logic_network_data,
            reactome_id_to_uuid,
        )

        assert len(pathway_logic_network_data) == 2
        sources = {e["source_id"] for e in pathway_logic_network_data}
        assert len(sources) == 1, (
            "Same protein in catalyst and regulator role must share one UUID"
        )


class TestRegulatorDecompositionConsistency:
    """Test that regulator decomposition is consistent with pathway decomposition."""

    @patch('src.neo4j_connector.get_set_members')
    @patch('src.neo4j_connector.get_complex_components')
    @patch('src.neo4j_connector.get_labels')
    @patch('src.logic_network_generator._complex_contains_entity_set')
    def test_simple_complex_regulator_kept_intact(
        self, mock_contains_set, mock_labels, mock_components, mock_members
    ):
        """Simple complexes (no EntitySets) should be kept intact, not decomposed."""
        from src.logic_network_generator import _decompose_regulator_entity

        mock_labels.return_value = ["Complex", "PhysicalEntity"]
        mock_contains_set.return_value = False
        mock_components.return_value = {"R-HSA-A": 1, "R-HSA-B": 1}

        result = _decompose_regulator_entity("R-HSA-SIMPLE-COMPLEX")

        assert len(result) == 1, f"Simple complex should return single entity, got {len(result)}"
        assert result[0] == ("R-HSA-SIMPLE-COMPLEX", 1)

    @patch('src.neo4j_connector.get_set_members')
    @patch('src.neo4j_connector.get_complex_components')
    @patch('src.neo4j_connector.get_labels')
    @patch('src.logic_network_generator._complex_contains_entity_set')
    def test_complex_with_entityset_regulator_decomposed(
        self, mock_contains_set, mock_labels, mock_components, mock_members
    ):
        """Complexes containing EntitySets should be fully decomposed."""
        from src.logic_network_generator import _decompose_regulator_entity

        # Return different labels based on entity_id
        def labels_side_effect(entity_id):
            if entity_id == "R-HSA-COMPLEX-WITH-SET":
                return ["Complex", "PhysicalEntity"]
            elif entity_id == "R-HSA-PROTEIN-A":
                return ["EntityWithAccessionedSequence", "PhysicalEntity"]
            elif entity_id == "R-HSA-PROTEIN-B":
                return ["EntityWithAccessionedSequence", "PhysicalEntity"]
            return ["PhysicalEntity"]

        mock_labels.side_effect = labels_side_effect
        mock_contains_set.return_value = True
        mock_components.return_value = {"R-HSA-PROTEIN-A": 2, "R-HSA-PROTEIN-B": 1}

        result = _decompose_regulator_entity("R-HSA-COMPLEX-WITH-SET")

        assert len(result) == 2, f"Complex with 2 components should return 2 members, got {len(result)}"
        member_ids = {r[0] for r in result}
        assert member_ids == {"R-HSA-PROTEIN-A", "R-HSA-PROTEIN-B"}
        # Check stoichiometry is preserved
        stoich_map = dict(result)
        assert stoich_map["R-HSA-PROTEIN-A"] == 2
        assert stoich_map["R-HSA-PROTEIN-B"] == 1


class TestRealNetworkRegulators:
    """Test regulators in actual generated networks (if available)."""

    @pytest.mark.skipif(
        not any(
            (d / "logic_network.csv").exists()
            for d in Path("output").iterdir()
            if d.is_dir()
        ) if Path("output").exists() else True,
        reason="No generated pathway directories found in output/"
    )
    def test_real_network_has_negative_regulators(self):
        """If real network exists, verify it has properly marked negative regulators."""
        network_path = next(
            d / "logic_network.csv"
            for d in sorted(Path("output").iterdir())
            if d.is_dir() and (d / "logic_network.csv").exists()
        )
        network = pd.read_csv(network_path)

        # Get all regulatory edges
        regulator_edges = network[network['edge_type'] == 'regulator']

        if len(regulator_edges) > 0:
            # Check for negative regulators
            negative_regulators = regulator_edges[regulator_edges['pos_neg'] == 'neg']
            positive_regulators = regulator_edges[regulator_edges['pos_neg'] == 'pos']

            print("\nRegulator statistics:")
            print(f"  Total regulators: {len(regulator_edges)}")
            print(f"  Negative regulators: {len(negative_regulators)}")
            print(f"  Positive regulators: {len(positive_regulators)}")

            # All regulators should be either positive or negative
            assert len(negative_regulators) + len(positive_regulators) == len(regulator_edges), \
                "All regulators should be marked as either positive or negative"

    @pytest.mark.skipif(
        not any(
            (d / "logic_network.csv").exists()
            for d in Path("output").iterdir()
            if d.is_dir()
        ) if Path("output").exists() else True,
        reason="No generated pathway directories found in output/"
    )
    def test_real_network_catalysts_are_positive(self):
        """If real network exists, verify all catalysts are positive."""
        network_path = next(
            d / "logic_network.csv"
            for d in sorted(Path("output").iterdir())
            if d.is_dir() and (d / "logic_network.csv").exists()
        )
        network = pd.read_csv(network_path)

        catalyst_edges = network[network['edge_type'] == 'catalyst']

        if len(catalyst_edges) > 0:
            # All catalysts should be positive
            negative_catalysts = catalyst_edges[catalyst_edges['pos_neg'] == 'neg']

            assert len(negative_catalysts) == 0, \
                f"Found {len(negative_catalysts)} negative catalysts - catalysts should always be positive"

            print("\nCatalyst statistics:")
            print(f"  Total catalysts: {len(catalyst_edges)}")
            print("  All catalysts are positive")
