"""Integration test for edge direction using synthetic pathway data.

This test creates a simple synthetic pathway to verify edge direction:

Pathway: MoleculeA → Reaction1 → MoleculeX → Reaction2 → MoleculeY

Expected edges in the logic network:
  1. MoleculeA → MoleculeX (A is consumed by R1, X is produced by R1)
  2. MoleculeX → MoleculeY (X is consumed by R2, Y is produced by R2)

This represents forward flow: root input → intermediate → terminal output
"""

import pytest
import pandas as pd
from typing import Dict, List, Any
import sys
from unittest.mock import patch

sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import extract_inputs_and_outputs


class TestEdgeDirectionIntegration:
    """Integration test for edge direction in pathway logic network."""

    def test_simple_two_reaction_pathway(self):
        """
        Test a simple pathway: R1 produces X, R2 consumes X.

        Reaction 1 (preceding):
          - No inputs (root)
          - Output: MoleculeX (Reactome ID: 1001)

        Reaction 2 (following):
          - Input: MoleculeX (Reactome ID: 1001)
          - Output: MoleculeY (Reactome ID: 1002)

        Expected edge: MoleculeX (from R1 output) → MoleculeX (to R2 input)
        Since it's the same physical entity, we expect UUID to be reused.
        Expected flow semantics: preceding_output → current_input
        """

        # Create synthetic reaction_id_map
        # Each reaction has a UUID, reactome_id, input_hash, and output_hash
        reaction_id_map = pd.DataFrame([
            {
                "uid": "reaction-1-uuid",
                "reactome_id": 100,
                "input_hash": "input-hash-r1",  # R1 has no terminal inputs (root)
                "output_hash": "output-hash-r1",  # R1 outputs MoleculeX
            },
            {
                "uid": "reaction-2-uuid",
                "reactome_id": 200,
                "input_hash": "input-hash-r2",  # R2 inputs MoleculeX
                "output_hash": "output-hash-r2",  # R2 outputs MoleculeY
            }
        ])

        # Create synthetic decomposed_uid_mapping
        # This maps hashes to their terminal reactome IDs
        decomposed_uid_mapping = pd.DataFrame([
            # Reaction 1 output: MoleculeX (ID: 1001)
            {
                "uid": "output-hash-r1",
                "reactome_id": 100,
                "component_id": 0,
                "component_id_or_reference_entity_id": 0,
                "input_or_output_uid": None,
                "input_or_output_reactome_id": 1001,  # MoleculeX
            },
            # Reaction 2 input: MoleculeX (ID: 1001)
            {
                "uid": "input-hash-r2",
                "reactome_id": 200,
                "component_id": 0,
                "component_id_or_reference_entity_id": 0,
                "input_or_output_uid": None,
                "input_or_output_reactome_id": 1001,  # MoleculeX
            },
            # Reaction 2 output: MoleculeY (ID: 1002)
            {
                "uid": "output-hash-r2",
                "reactome_id": 200,
                "component_id": 0,
                "component_id_or_reference_entity_id": 0,
                "input_or_output_uid": None,
                "input_or_output_reactome_id": 1002,  # MoleculeY
            },
        ])

        # Create uid_reaction_connections: R1 precedes R2
        uid_reaction_connections = pd.DataFrame([
            {
                "preceding_uid": "reaction-1-uuid",
                "following_uid": "reaction-2-uuid",
            }
        ])

        # Prepare data structures
        reaction_uids = ["reaction-2-uuid"]  # Process reaction 2
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        # Run the function
        extract_inputs_and_outputs(
            reaction_uid="reaction-2-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        # Verify results
        assert len(pathway_logic_network_data) == 1, "Should create exactly one edge"

        edge = pathway_logic_network_data[0]

        # Both source and target should have the same UUID (it's the same physical entity)
        molecule_x_uuid = reactome_id_to_uuid.get(1001) or reactome_id_to_uuid.get(1001.0)
        assert molecule_x_uuid is not None, "MoleculeX should have been assigned a UUID"

        print(f"\n=== Test Results ===")
        print(f"MoleculeX UUID: {molecule_x_uuid}")
        print(f"Edge created: {edge['source_id']} → {edge['target_id']}")
        print(f"AND/OR: {edge['and_or']}, Edge Type: {edge['edge_type']}")

        # CRITICAL VERIFICATION: Check edge direction
        # Scenario: R1 produces MoleculeX, R2 consumes MoleculeX
        # Expected: MoleculeX flows from R1's output to R2's input

        # The key question: what do source_id and target_id represent?
        # Option A (forward flow): source = R1's output X, target = R2's input X
        #   Both are the same molecule, so source_id == target_id == molecule_x_uuid
        # Option B (backward flow): source = R2's input X, target = R1's output X
        #   Both are the same molecule, so source_id == target_id == molecule_x_uuid

        # Since they're the same molecule, we can't distinguish forward from backward!
        # This is a self-loop edge, which reveals a problem with the test design.

        assert edge['source_id'] == molecule_x_uuid
        assert edge['target_id'] == molecule_x_uuid

        print("\n=== Issue Identified ===")
        print("When the same molecule appears as both output of R1 and input of R2,")
        print("we get a self-loop edge. This doesn't help us verify direction.")
        print("\nWe need a test with DIFFERENT molecules at each stage.")

    def test_three_reaction_pathway_with_distinct_molecules(self):
        """
        Test pathway with distinct molecules at each stage.

        Pathway structure:
          R1: produces MolA (1001)
          R2: consumes MolA, produces MolB (1002)
          R3: consumes MolB, produces MolC (1003)

        Expected edges for forward flow (output → input):
          R1_output(MolA) → R2_input(MolA) - but these are same molecule!
          R2_output(MolB) → R3_input(MolB) - but these are same molecule!

        The issue: we're creating molecule→molecule edges, not reaction→reaction edges.
        And molecules are identified by their Reactome ID, not by which reaction they belong to.

        So MolA from R1's output is THE SAME NODE as MolA in R2's input.

        This means we CANNOT have edges between them - they're the same node!

        The real edges must be connecting DIFFERENT molecules:
          MolA → MolB (representing the transformation through R2)
          MolB → MolC (representing the transformation through R3)

        But wait - that's not what the code does. Let me re-examine...

        The code connects:
          current reaction's INPUT molecules → preceding reaction's OUTPUT molecules

        For R2 (current), R1 (preceding):
          R2_inputs = [MolA]
          R1_outputs = [MolA]
          Creates edge: MolA → MolA (self-loop!)

        This seems wrong. Unless... the molecules have different representations?
        Or maybe the logic is different than I think?
        """

        # Actually, let me check what happens when inputs and outputs are DIFFERENT
        # R1: no inputs, output = MolA
        # R2: input = MolA, output = MolB

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
        ])

        decomposed_uid_mapping = pd.DataFrame([
            # R1 outputs MolA
            {"uid": "r1-output-hash", "reactome_id": 100, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            # R2 inputs MolA
            {"uid": "r2-input-hash", "reactome_id": 200, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1001},
            # R2 outputs MolB
            {"uid": "r2-output-hash", "reactome_id": 200, "component_id": 0,
             "component_id_or_reference_entity_id": 0, "input_or_output_uid": None,
             "input_or_output_reactome_id": 1002},
        ])

        uid_reaction_connections = pd.DataFrame([
            {"preceding_uid": "r1-uuid", "following_uid": "r2-uuid"}
        ])

        reaction_uids = ["r2-uuid"]
        reactome_id_to_uuid: Dict[str, str] = {}
        pathway_logic_network_data: List[Dict[str, Any]] = []

        extract_inputs_and_outputs(
            reaction_uid="r2-uuid",
            reaction_uids=reaction_uids,
            uid_reaction_connections=uid_reaction_connections,
            reaction_id_map=reaction_id_map,
            decomposed_uid_mapping=decomposed_uid_mapping,
            reactome_id_to_uuid=reactome_id_to_uuid,
            pathway_logic_network_data=pathway_logic_network_data,
        )

        print(f"\n=== Test Results for Distinct Molecules ===")
        print(f"Number of edges created: {len(pathway_logic_network_data)}")
        print(f"Reactome ID to UUID mapping: {reactome_id_to_uuid}")

        for i, edge in enumerate(pathway_logic_network_data):
            print(f"Edge {i}: {edge['source_id']} → {edge['target_id']}")
            # Find which physical entity this is
            for reactome_id, uuid in reactome_id_to_uuid.items():
                if uuid == edge['source_id']:
                    print(f"  Source is Physical Entity with Reactome ID {reactome_id}")
                if uuid == edge['target_id']:
                    print(f"  Target is Physical Entity with Reactome ID {reactome_id}")

        # Get UUIDs for our physical entities (keys might be int or float)
        entity_a_uuid = reactome_id_to_uuid.get(1001) or reactome_id_to_uuid.get(1001.0)
        entity_b_uuid = reactome_id_to_uuid.get(1002) or reactome_id_to_uuid.get(1002.0)

        assert len(pathway_logic_network_data) == 1
        edge = pathway_logic_network_data[0]

        print(f"\nEntityA UUID: {entity_a_uuid}")
        print(f"EntityB UUID: {entity_b_uuid}")
        print(f"Edge: {edge['source_id']} → {edge['target_id']}")

        # NOW we can test direction!
        # Current code: input_uuid → output_uuid
        # Where input_uuid = R2's input = EntityA
        # And output_uuid = R1's output = EntityA
        # So edge would be: EntityA → EntityA (self-loop again!)

        # Hmm, still a self-loop. The issue is that EntityA appears in both
        # R2's input list and R1's output list, and they get the SAME UUID.

        assert edge['source_id'] == entity_a_uuid, "Current code creates self-loop"
        assert edge['target_id'] == entity_a_uuid, "Both ends are the same physical entity"

        print("\n=== Conclusion ===")
        print("We're still getting self-loops because:")
        print("  R2's input (EntityA) and R1's output (EntityA) have the same UUID")
        print("\nThis suggests the edges DON'T represent physical entity flow between reactions.")
        print("Instead, they might represent something else entirely.")
        print("\nNeed to re-examine the actual pathway_logic_network_69620.csv data")
        print("to understand what non-self-loop edges actually represent.")
