#!/usr/bin/env python3
"""
Comprehensive validation script for generated logic networks.

This script validates that the logic network generation is working correctly by:
1. Checking the structure of the logic network
2. Validating UUID mappings
3. Reconstructing Reactome reactions from the logic network
4. Comparing with Neo4j to verify correctness
5. Validating regulator and catalyst propagation

Usage:
    python scripts/validate_logic_network.py --pathway-id 69620
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd
from py2neo import Graph
import os


class ValidationResult:
    """Container for validation results."""

    def __init__(self, test_name: str):
        self.test_name = test_name
        self.passed = True
        self.errors = []
        self.warnings = []
        self.info = []

    def fail(self, message: str):
        """Mark test as failed with error message."""
        self.passed = False
        self.errors.append(message)

    def warn(self, message: str):
        """Add warning message."""
        self.warnings.append(message)

    def add_info(self, message: str):
        """Add informational message."""
        self.info.append(message)

    def print_result(self):
        """Print the validation result."""
        status = "✅ PASS" if self.passed else "❌ FAIL"
        print(f"\n{status}: {self.test_name}")

        for info in self.info:
            print(f"  ℹ️  {info}")

        for warning in self.warnings:
            print(f"  ⚠️  {warning}")

        for error in self.errors:
            print(f"  ❌ {error}")


class LogicNetworkValidator:
    """Validates a generated logic network against Neo4j."""

    def __init__(self, pathway_id: int):
        self.pathway_id = pathway_id
        self.output_dir = Path("output")

        # Connect to Neo4j
        uri = os.getenv("NEO4J_URI", "bolt://localhost:7687")
        self.graph = Graph(uri, auth=("neo4j", "test"))

        # Load generated files
        self.logic_network = None
        self.uuid_to_reactome = None
        self.decomposed_uid_mapping = None

    def load_files(self) -> ValidationResult:
        """Load all required files."""
        result = ValidationResult("File Loading")

        try:
            # Load logic network
            logic_network_file = self.output_dir / f"pathway_logic_network_{self.pathway_id}.csv"
            if not logic_network_file.exists():
                result.fail(f"Logic network file not found: {logic_network_file}")
                return result

            self.logic_network = pd.read_csv(logic_network_file)
            result.add_info(f"Loaded logic network: {len(self.logic_network)} edges")

            # Load UUID to Reactome mapping
            uuid_to_reactome_file = self.output_dir / f"uuid_to_reactome_{self.pathway_id}.csv"
            if not uuid_to_reactome_file.exists():
                result.fail(f"UUID mapping file not found: {uuid_to_reactome_file}")
                return result

            self.uuid_to_reactome = pd.read_csv(uuid_to_reactome_file)
            result.add_info(f"Loaded UUID mappings: {len(self.uuid_to_reactome)} entries")

            # Load decomposed UID mapping
            decomposed_file = self.output_dir / f"decomposed_uid_mapping_{self.pathway_id}.csv"
            if not decomposed_file.exists():
                result.fail(f"Decomposed mapping file not found: {decomposed_file}")
                return result

            self.decomposed_uid_mapping = pd.read_csv(decomposed_file)
            result.add_info(f"Loaded decomposed mappings: {len(self.decomposed_uid_mapping)} entries")

        except Exception as e:
            result.fail(f"Error loading files: {str(e)}")

        return result

    def validate_structure(self) -> ValidationResult:
        """Validate the structure of the logic network."""
        result = ValidationResult("Logic Network Structure")

        # Check required columns
        required_cols = {'source_id', 'target_id', 'pos_neg', 'and_or', 'edge_type'}
        actual_cols = set(self.logic_network.columns)

        if not required_cols.issubset(actual_cols):
            missing = required_cols - actual_cols
            result.fail(f"Missing required columns: {missing}")
            return result

        result.add_info("All required columns present")

        # Check edge types
        edge_types = self.logic_network['edge_type'].unique()
        valid_edge_types = {'input', 'output', 'catalyst', 'regulator'}
        invalid_types = set(edge_types) - valid_edge_types

        if invalid_types:
            result.fail(f"Invalid edge types found: {invalid_types}")
        else:
            result.add_info(f"Valid edge types: {list(edge_types)}")

        # Check pos_neg values
        pos_neg_values = self.logic_network['pos_neg'].dropna().unique()
        valid_pos_neg = {'pos', 'neg'}
        invalid_pos_neg = set(pos_neg_values) - valid_pos_neg

        if invalid_pos_neg:
            result.fail(f"Invalid pos_neg values found: {invalid_pos_neg}")
        else:
            result.add_info(f"Valid pos_neg values: {list(pos_neg_values)}")

        # Check for null UUIDs
        null_sources = self.logic_network['source_id'].isna().sum()
        null_targets = self.logic_network['target_id'].isna().sum()

        if null_sources > 0 or null_targets > 0:
            result.fail(f"Found null UUIDs: {null_sources} sources, {null_targets} targets")

        # Print edge type distribution
        edge_dist = self.logic_network['edge_type'].value_counts()
        result.add_info(f"Edge distribution: {edge_dist.to_dict()}")

        return result

    def validate_uuid_mapping(self) -> ValidationResult:
        """Validate that all entity UUIDs can be mapped to Reactome IDs."""
        result = ValidationResult("UUID Mapping Completeness")

        # Get all UUIDs from logic network
        all_uuids_in_network = set(self.logic_network['source_id'].unique()) | \
                               set(self.logic_network['target_id'].unique())

        # Build UUID lookup from mapping file (only contains entity UUIDs, not reaction UUIDs)
        entity_uuids_in_mapping = set(self.uuid_to_reactome['uuid'].unique())

        # Identify reaction UUIDs (appear as targets of input edges or sources of output edges)
        input_edges = self.logic_network[self.logic_network['edge_type'] == 'input']
        output_edges = self.logic_network[self.logic_network['edge_type'] == 'output']
        reaction_uuids = set(input_edges['target_id'].unique()) | set(output_edges['source_id'].unique())

        # Entity UUIDs are all UUIDs minus reaction UUIDs
        entity_uuids_in_network = all_uuids_in_network - reaction_uuids

        result.add_info(f"Total UUIDs in logic network: {len(all_uuids_in_network)}")
        result.add_info(f"  Entity UUIDs: {len(entity_uuids_in_network)}")
        result.add_info(f"  Reaction UUIDs: {len(reaction_uuids)}")

        # Check if all entity UUIDs are in the mapping file
        unmappable_entities = entity_uuids_in_network - entity_uuids_in_mapping

        if unmappable_entities:
            result.fail(f"Found {len(unmappable_entities)} entity UUIDs not in mapping file")
            for uuid_val in list(unmappable_entities)[:5]:  # Show first 5
                result.fail(f"  Unmappable entity: {uuid_val}")
        else:
            result.add_info(f"All {len(entity_uuids_in_network)} entity UUIDs are in mapping file")

        # Check for empty mappings
        empty_mappings = 0
        for _, row in self.uuid_to_reactome.iterrows():
            entity_ids_str = row['entity_ids']
            if pd.isna(entity_ids_str) or not entity_ids_str or entity_ids_str.strip() == '':
                empty_mappings += 1

        if empty_mappings > 0:
            result.warn(f"{empty_mappings} UUIDs have empty entity_ids mappings")
        else:
            result.add_info("All entity UUIDs map to at least one Reactome entity ID")

        return result

    def validate_regulator_propagation(self) -> ValidationResult:
        """Validate that regulators are properly propagated from Neo4j."""
        result = ValidationResult("Regulator Propagation")

        # Query Neo4j for regulators
        positive_query = f"""
        MATCH (pathway:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(regulator:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_pos_count = self.graph.run(positive_query).data()[0]['count']

        negative_query = f"""
        MATCH (pathway:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(regulator:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_neg_count = self.graph.run(negative_query).data()[0]['count']

        catalyst_query = f"""
        MATCH (pathway:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:catalystActivity]->(ca:CatalystActivity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_catalyst_count = self.graph.run(catalyst_query).data()[0]['count']

        # Count in logic network
        regulator_edges = self.logic_network[self.logic_network['edge_type'] == 'regulator']
        logic_pos_reactions = len(regulator_edges[regulator_edges['pos_neg'] == 'pos']['target_id'].unique())
        logic_neg_reactions = len(regulator_edges[regulator_edges['pos_neg'] == 'neg']['target_id'].unique())

        catalyst_edges = self.logic_network[self.logic_network['edge_type'] == 'catalyst']
        logic_catalyst_reactions = len(catalyst_edges['target_id'].unique())

        result.add_info(f"Neo4j: {neo4j_pos_count} reactions with positive regulators")
        result.add_info(f"Logic network: {logic_pos_reactions} virtual reactions with positive regulators")

        result.add_info(f"Neo4j: {neo4j_neg_count} reactions with negative regulators")
        result.add_info(f"Logic network: {logic_neg_reactions} virtual reactions with negative regulators")

        result.add_info(f"Neo4j: {neo4j_catalyst_count} reactions with catalysts")
        result.add_info(f"Logic network: {logic_catalyst_reactions} virtual reactions with catalysts")

        # Note: Logic network may have more because of EntitySet decomposition
        if logic_pos_reactions >= neo4j_pos_count:
            result.add_info("Positive regulators: ✓ (may be duplicated for virtual reactions)")
        else:
            result.warn(f"Missing positive regulators: expected >={neo4j_pos_count}, got {logic_pos_reactions}")

        if logic_neg_reactions >= neo4j_neg_count:
            result.add_info("Negative regulators: ✓ (may be duplicated for virtual reactions)")
        else:
            result.warn(f"Missing negative regulators: expected >={neo4j_neg_count}, got {logic_neg_reactions}")

        if logic_catalyst_reactions >= neo4j_catalyst_count:
            result.add_info("Catalysts: ✓ (may be duplicated for virtual reactions)")
        else:
            result.warn(f"Missing catalysts: expected >={neo4j_catalyst_count}, got {logic_catalyst_reactions}")

        return result

    def validate_reconstruction(self) -> ValidationResult:
        """Validate that the logic network can reconstruct the original pathway."""
        result = ValidationResult("Pathway Reconstruction")

        # Build UUID lookup
        uuid_dict = {}
        for _, row in self.uuid_to_reactome.iterrows():
            uuid_val = row['uuid']
            entity_ids_str = row['entity_ids']
            if pd.notna(entity_ids_str) and entity_ids_str:
                entity_ids = set(int(eid) for eid in entity_ids_str.split('|') if eid)
                uuid_dict[uuid_val] = entity_ids

        # Get input and output edges
        input_edges = self.logic_network[self.logic_network['edge_type'] == 'input']
        output_edges = self.logic_network[self.logic_network['edge_type'] == 'output']

        # Find all virtual reactions (they appear as targets of input edges and sources of output edges)
        reaction_uuids = set(input_edges['target_id'].unique()) | set(output_edges['source_id'].unique())

        # For each virtual reaction, reconstruct its input→output pairs
        all_edges = []
        unconvertible_reactions = 0

        for reaction_uuid in reaction_uuids:
            # Get inputs to this reaction
            reaction_inputs = input_edges[input_edges['target_id'] == reaction_uuid]
            input_entity_uuids = set(reaction_inputs['source_id'].unique())

            # Get outputs from this reaction
            reaction_outputs = output_edges[output_edges['source_id'] == reaction_uuid]
            output_entity_uuids = set(reaction_outputs['target_id'].unique())

            # Convert to Reactome IDs
            input_reactome_ids = set()
            for uuid_val in input_entity_uuids:
                if uuid_val in uuid_dict:
                    input_reactome_ids.update(uuid_dict[uuid_val])

            output_reactome_ids = set()
            for uuid_val in output_entity_uuids:
                if uuid_val in uuid_dict:
                    output_reactome_ids.update(uuid_dict[uuid_val])

            if not input_reactome_ids or not output_reactome_ids:
                unconvertible_reactions += 1
                continue

            # Create all input×output pairs for this reaction
            for inp in input_reactome_ids:
                for outp in output_reactome_ids:
                    all_edges.append((inp, outp))

        # Deduplicate
        unique_edges = set(all_edges)

        result.add_info(f"Found {len(reaction_uuids)} virtual reactions in logic network")
        result.add_info(f"Reconstructed {len(all_edges)} Reactome input→output pairs")
        result.add_info(f"After deduplication: {len(unique_edges)} unique pairs")

        if unconvertible_reactions > 0:
            result.warn(f"{unconvertible_reactions} virtual reactions could not be fully converted")
        else:
            result.add_info("All virtual reactions successfully converted")

        # Get Neo4j reactions
        query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        OPTIONAL MATCH (r)-[:input]->(inp)
        OPTIONAL MATCH (r)-[:output]->(out)
        WITH r, collect(DISTINCT inp.dbId) AS inputs, collect(DISTINCT out.dbId) AS outputs
        RETURN r.dbId AS reaction_id,
               [x IN inputs WHERE x IS NOT NULL] AS inputs,
               [x IN outputs WHERE x IS NOT NULL] AS outputs
        """

        neo4j_reaction_pairs = set()
        reactions_data = self.graph.run(query).data()

        for row in reactions_data:
            inputs = row["inputs"]
            outputs = row["outputs"]
            for inp in inputs:
                for outp in outputs:
                    neo4j_reaction_pairs.add((inp, outp))

        result.add_info(f"Neo4j: {len(neo4j_reaction_pairs)} input→output pairs")

        # Compare
        missing = neo4j_reaction_pairs - unique_edges
        extra = unique_edges - neo4j_reaction_pairs
        matches = len(neo4j_reaction_pairs) - len(missing)
        accuracy = (matches / len(neo4j_reaction_pairs) * 100) if len(neo4j_reaction_pairs) > 0 else 0

        result.add_info(f"Matching: {matches}/{len(neo4j_reaction_pairs)} ({accuracy:.1f}%)")

        if accuracy == 100.0:
            result.add_info("🎉 Perfect reconstruction!")
        elif accuracy >= 90:
            result.add_info("Good reconstruction (>90%)")
        else:
            result.warn(f"Reconstruction accuracy below 90%: {accuracy:.1f}%")

        if missing:
            result.warn(f"{len(missing)} edges in Neo4j but not in logic network")

        if extra:
            result.warn(f"{len(extra)} edges in logic network but not in Neo4j")

        return result

    def validate_no_spurious_self_loops(self) -> ValidationResult:
        """Verify no inappropriate self-loops exist at UUID level."""
        result = ValidationResult("Self-Loop Detection")

        # Check each edge type for self-loops
        for edge_type in ['input', 'output', 'catalyst', 'regulator']:
            edges = self.logic_network[self.logic_network['edge_type'] == edge_type]
            self_loops = edges[edges['source_id'] == edges['target_id']]

            if len(self_loops) > 0:
                result.warn(f"{edge_type} has {len(self_loops)} self-loops at UUID level")
                # Show examples
                for _, edge in self_loops.head(3).iterrows():
                    result.warn(f"  Example: {edge['source_id']} → {edge['target_id']}")
            else:
                result.add_info(f"{edge_type}: No self-loops ✓")

        return result

    def validate_entity_coverage(self) -> ValidationResult:
        """Verify all Neo4j entities appear in logic network."""
        result = ValidationResult("Entity Coverage")

        # Get all entities from Neo4j (inputs and outputs)
        query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:input|output]->(entity:PhysicalEntity)
        RETURN COLLECT(DISTINCT entity.dbId) as entity_ids
        """
        neo4j_result = self.graph.run(query).data()
        neo4j_entities = set(neo4j_result[0]['entity_ids']) if neo4j_result else set()

        # Get all entities from logic network via uuid_to_reactome mapping
        ln_entities = set()
        for _, row in self.uuid_to_reactome.iterrows():
            entity_ids_str = row['entity_ids']
            if pd.notna(entity_ids_str):
                entity_ids = set(int(eid) for eid in str(entity_ids_str).split('|') if eid)
                ln_entities.update(entity_ids)

        missing_entities = neo4j_entities - ln_entities
        extra_entities = ln_entities - neo4j_entities

        result.add_info(f"Neo4j entities: {len(neo4j_entities)}")
        result.add_info(f"Logic network entities: {len(ln_entities)}")

        if missing_entities:
            result.fail(f"Missing {len(missing_entities)} entities from Neo4j")
            for entity_id in list(missing_entities)[:5]:
                result.fail(f"  Missing entity: {entity_id}")
        else:
            result.add_info("All Neo4j entities present ✓")

        if extra_entities:
            result.add_info(f"Logic network has {len(extra_entities)} extra entities (from catalysts/regulators)")

        return result

    def validate_catalyst_completeness(self) -> ValidationResult:
        """Verify all Neo4j catalysts are present in logic network."""
        result = ValidationResult("Catalyst Completeness")

        # Get catalysts from Neo4j
        query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:catalystActivity]->(ca)-[:physicalEntity]->(catalyst)
        RETURN COLLECT(DISTINCT catalyst.dbId) as catalyst_ids
        """
        neo4j_result = self.graph.run(query).data()
        neo4j_catalysts = set(neo4j_result[0]['catalyst_ids']) if neo4j_result else set()

        # Get catalysts from logic network
        catalyst_edges = self.logic_network[self.logic_network['edge_type'] == 'catalyst']
        ln_catalysts = set()

        for catalyst_uuid in catalyst_edges['source_id'].unique():
            # Look up in uuid_to_reactome
            mapping = self.uuid_to_reactome[self.uuid_to_reactome['uuid'] == catalyst_uuid]
            if not mapping.empty:
                entity_ids_str = mapping.iloc[0]['entity_ids']
                if pd.notna(entity_ids_str):
                    entity_ids = set(int(eid) for eid in str(entity_ids_str).split('|') if eid)
                    ln_catalysts.update(entity_ids)

        missing = neo4j_catalysts - ln_catalysts

        result.add_info(f"Neo4j catalysts: {len(neo4j_catalysts)}")
        result.add_info(f"Logic network catalysts: {len(ln_catalysts)}")

        if missing:
            result.fail(f"Missing {len(missing)} catalysts from Neo4j")
            for catalyst_id in list(missing)[:5]:
                result.fail(f"  Missing catalyst: {catalyst_id}")
        else:
            result.add_info("All catalysts present ✓")

        return result

    def validate_regulator_polarity(self) -> ValidationResult:
        """Verify regulator pos_neg values match Neo4j."""
        result = ValidationResult("Regulator Polarity")

        # Get positive regulators from Neo4j
        pos_query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:regulatedBy]->(reg:PositiveRegulation)-[:regulator]->(pe)
        RETURN COLLECT(DISTINCT pe.dbId) as regulator_ids
        """
        pos_result = self.graph.run(pos_query).data()
        neo4j_positive = set(pos_result[0]['regulator_ids']) if pos_result else set()

        # Get negative regulators from Neo4j
        neg_query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:regulatedBy]->(reg:NegativeRegulation)-[:regulator]->(pe)
        RETURN COLLECT(DISTINCT pe.dbId) as regulator_ids
        """
        neg_result = self.graph.run(neg_query).data()
        neo4j_negative = set(neg_result[0]['regulator_ids']) if neg_result else set()

        # Check logic network regulators
        regulator_edges = self.logic_network[self.logic_network['edge_type'] == 'regulator']

        pos_mismatches = []
        neg_mismatches = []
        checked_count = 0

        for _, edge in regulator_edges.iterrows():
            reg_uuid = edge['source_id']
            pos_neg = edge['pos_neg']

            # Look up Reactome ID
            mapping = self.uuid_to_reactome[self.uuid_to_reactome['uuid'] == reg_uuid]
            if mapping.empty:
                continue

            entity_ids_str = mapping.iloc[0]['entity_ids']
            if pd.notna(entity_ids_str):
                entity_id = int(str(entity_ids_str).split('|')[0])
                checked_count += 1

                # Check if polarity matches
                if entity_id in neo4j_positive and pos_neg != 'pos':
                    pos_mismatches.append(entity_id)
                if entity_id in neo4j_negative and pos_neg != 'neg':
                    neg_mismatches.append(entity_id)

        result.add_info(f"Checked {checked_count} regulator edges")
        result.add_info(f"Neo4j: {len(neo4j_positive)} positive, {len(neo4j_negative)} negative")

        if pos_mismatches:
            result.fail(f"Positive regulators with wrong polarity: {pos_mismatches}")
        if neg_mismatches:
            result.fail(f"Negative regulators with wrong polarity: {neg_mismatches}")

        if not pos_mismatches and not neg_mismatches:
            result.add_info("All regulator polarities correct ✓")

        return result

    def validate_reaction_coverage(self) -> ValidationResult:
        """Verify all Neo4j reactions are represented in logic network."""
        result = ValidationResult("Reaction Coverage")

        # Get all reactions from Neo4j
        query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        RETURN COUNT(DISTINCT r) as reaction_count
        """
        neo4j_result = self.graph.run(query).data()
        neo4j_reaction_count = neo4j_result[0]['reaction_count'] if neo4j_result else 0

        # Count reactions in logic network (reaction UUIDs are targets of input edges)
        input_edges = self.logic_network[self.logic_network['edge_type'] == 'input']
        ln_reaction_count = input_edges['target_id'].nunique()

        result.add_info(f"Neo4j reactions: {neo4j_reaction_count}")
        result.add_info(f"Logic network reactions: {ln_reaction_count}")

        if ln_reaction_count < neo4j_reaction_count:
            result.fail(f"Missing {neo4j_reaction_count - ln_reaction_count} reactions")
        elif ln_reaction_count > neo4j_reaction_count:
            extra = ln_reaction_count - neo4j_reaction_count
            result.add_info(f"Logic network has {extra} virtual reactions (from EntitySet expansion) ✓")
        else:
            result.add_info("All reactions present (no EntitySet expansion) ✓")

        return result

    def validate_edge_counts(self) -> ValidationResult:
        """Compare edge counts with Neo4j."""
        result = ValidationResult("Edge Count Verification")

        # Query Neo4j for unique entity counts per edge type
        query = f"""
        MATCH (p:Pathway {{dbId: {self.pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        OPTIONAL MATCH (r)-[:input]->(inp)
        OPTIONAL MATCH (r)-[:output]->(out)
        OPTIONAL MATCH (r)-[:catalystActivity]->(ca)-[:physicalEntity]->(cat)
        OPTIONAL MATCH (r)-[:regulatedBy]->(reg)-[:regulator]->(regulator)
        RETURN
            COUNT(DISTINCT inp) as input_count,
            COUNT(DISTINCT out) as output_count,
            COUNT(DISTINCT cat) as catalyst_count,
            COUNT(DISTINCT regulator) as regulator_count
        """

        neo4j_result = self.graph.run(query).data()
        neo4j_counts = neo4j_result[0] if neo4j_result else {}

        # Get logic network edge counts
        ln_inputs = len(self.logic_network[self.logic_network['edge_type'] == 'input'])
        ln_outputs = len(self.logic_network[self.logic_network['edge_type'] == 'output'])
        ln_catalysts = len(self.logic_network[self.logic_network['edge_type'] == 'catalyst'])
        ln_regulators = len(self.logic_network[self.logic_network['edge_type'] == 'regulator'])

        result.add_info(f"Input edges: Neo4j entities={neo4j_counts.get('input_count', 0)}, LN edges={ln_inputs}")
        result.add_info(f"Output edges: Neo4j entities={neo4j_counts.get('output_count', 0)}, LN edges={ln_outputs}")
        result.add_info(f"Catalyst edges: Neo4j entities={neo4j_counts.get('catalyst_count', 0)}, LN edges={ln_catalysts}")
        result.add_info(f"Regulator edges: Neo4j entities={neo4j_counts.get('regulator_count', 0)}, LN edges={ln_regulators}")

        # Note: Logic network can have MORE edges due to EntitySet expansion
        result.add_info("Note: Logic network may have more edges due to EntitySet expansion")

        return result

    def run_all_validations(self) -> bool:
        """Run all validations and return overall success."""
        print("=" * 80)
        print(f"LOGIC NETWORK VALIDATION - Pathway {self.pathway_id}")
        print("=" * 80)

        results = []

        # Load files
        load_result = self.load_files()
        load_result.print_result()
        results.append(load_result)

        if not load_result.passed:
            print("\n❌ Cannot continue validation - failed to load files")
            return False

        # Run validations
        results.append(self.validate_structure())
        results[-1].print_result()

        results.append(self.validate_uuid_mapping())
        results[-1].print_result()

        results.append(self.validate_no_spurious_self_loops())
        results[-1].print_result()

        results.append(self.validate_entity_coverage())
        results[-1].print_result()

        results.append(self.validate_catalyst_completeness())
        results[-1].print_result()

        results.append(self.validate_regulator_polarity())
        results[-1].print_result()

        results.append(self.validate_reaction_coverage())
        results[-1].print_result()

        results.append(self.validate_edge_counts())
        results[-1].print_result()

        results.append(self.validate_regulator_propagation())
        results[-1].print_result()

        results.append(self.validate_reconstruction())
        results[-1].print_result()

        # Print summary
        print("\n" + "=" * 80)
        print("VALIDATION SUMMARY")
        print("=" * 80)

        passed = sum(1 for r in results if r.passed)
        total = len(results)

        print(f"\nTests passed: {passed}/{total}")

        all_passed = all(r.passed for r in results)
        if all_passed:
            print("\n✅ ALL VALIDATIONS PASSED")
        else:
            print("\n❌ SOME VALIDATIONS FAILED")

        return all_passed


def main():
    parser = argparse.ArgumentParser(description="Validate generated logic network")
    parser.add_argument(
        "--pathway-id",
        type=int,
        required=True,
        help="Reactome pathway ID to validate"
    )

    args = parser.parse_args()

    validator = LogicNetworkValidator(args.pathway_id)
    success = validator.run_all_validations()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
