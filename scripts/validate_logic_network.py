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
from typing import Set

import pandas as pd
from py2neo import Graph
import os

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.credential_redaction import install as _install_credential_redaction

_install_credential_redaction()

# Depth limit when resolving an EntitySet to its leaf members. Sets nest a few
# levels at most; an unbounded traversal is slow and can revisit cycles.
MAX_SET_NESTING = 5


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

    def __init__(self, pathway_id, output_dir: Path = None):
        # Reactome stable ids are the primary identifier; a dbId is accepted
        # only because the generated output directories are named with one and
        # some callers still pass it. Both forms are resolved to the pair once,
        # and every query below matches on stId when one is available.
        self.pathway_id = pathway_id
        self.pathway_stid = None
        self.pathway_dbid = None
        self.output_dir = Path(output_dir) if output_dir else Path("output")

        # Connect to Neo4j. Every other module reads NEO4J_URL/USER/PASSWORD;
        # this script read NEO4J_URI and hardcoded auth=("neo4j","test"), so
        # pointing it at a secured or remote database silently validated
        # localhost instead. NEO4J_URI is still honoured as a fallback.
        uri = os.getenv("NEO4J_URL") or os.getenv("NEO4J_URI", "bolt://localhost:7687")
        user = os.getenv("NEO4J_USER", "neo4j")
        password = os.getenv("NEO4J_PASSWORD", "test")
        self.graph = Graph(uri, auth=(user, password))
        self._resolve_pathway_ids()

        # Load generated files
        self.logic_network = None
        self.uuid_to_reactome = None
        self.decomposed_uid_mapping = None

    def _resolve_pathway_ids(self):
        """Fill in both the stable id and the dbId from whichever was supplied."""
        raw = str(self.pathway_id).strip()
        if raw.upper().startswith("R-"):
            row = self.graph.run(
                "MATCH (p:Pathway {stId: $sid}) RETURN p.dbId AS db", sid=raw
            ).evaluate()
            self.pathway_stid = raw
            self.pathway_dbid = row
        else:
            try:
                dbid = int(raw)
            except ValueError:
                raise SystemExit(
                    f"--pathway-id must be a Reactome stable id (R-HSA-...) or a "
                    f"numeric dbId; got {raw!r}"
                )
            row = self.graph.run(
                "MATCH (p:Pathway {dbId: $db}) RETURN p.stId AS sid", db=dbid
            ).evaluate()
            self.pathway_dbid = dbid
            self.pathway_stid = row
        # Unreachable as originally written: both branches above assign one of
        # the two fields unconditionally, so the conjunction was never true and
        # a nonexistent id sailed through to produce the malformed glob
        # "*_R-HSA-R-HSA-99999999/". Test the field the LOOKUP was meant to
        # fill, which is the one that is None when the pathway does not exist.
        if self.pathway_dbid is None or self.pathway_stid is None:
            raise SystemExit(
                f"Pathway {raw!r} not found in Neo4j (Reactome release "
                f"{os.getenv('LNG_REACTOME_RELEASE', 'current')})."
            )
        # Directory names carry the stable id, Cypher prefers it too.
        self.pathway_id = self.pathway_dbid if self.pathway_dbid is not None else raw

    def _pathway_match(self, var: str = "p") -> str:
        """A Cypher pathway match clause bound to $pathway_key.

        stId when Reactome gives us one; dbId only as the fallback.
        """
        prop = "stId" if self.pathway_stid else "dbId"
        return f"({var}:Pathway {{{prop}: $pathway_key}})"

    @property
    def pathway_key(self):
        """The value bound to $pathway_key — the stable id where available."""
        return self.pathway_stid if self.pathway_stid else self.pathway_dbid

    def _resolve_layout(self):
        """Locate this pathway's files under either output layout.

        The generator emits output/<Name>_R-HSA-<id>/{logic_network.csv,
        stid_to_uuid_mapping.csv} plus cache/decomposed_uid_mapping.csv. This
        script only knew the legacy flat output/pathway_logic_network_{dbId}.csv
        naming, so load_files() failed before ANY check ran — every validation
        in here was unreachable against current output.
        """
        for pattern in (f"*_R-HSA-{self.pathway_id}", f"*_{self.pathway_id}"):
            for directory in sorted(self.output_dir.glob(pattern)):
                network = directory / "logic_network.csv"
                mapping = directory / "stid_to_uuid_mapping.csv"
                if network.exists() and mapping.exists():
                    return (
                        network,
                        mapping,
                        directory / "cache" / "decomposed_uid_mapping.csv",
                        f"per-pathway ({directory.name})",
                    )
        legacy_network = self.output_dir / f"pathway_logic_network_{self.pathway_id}.csv"
        if legacy_network.exists():
            return (
                legacy_network,
                self.output_dir / f"uuid_to_reactome_{self.pathway_id}.csv",
                self.output_dir / f"decomposed_uid_mapping_{self.pathway_id}.csv",
                "legacy flat",
            )
        return None, None, None, None

    def load_files(self) -> ValidationResult:
        """Load all required files."""
        result = ValidationResult("File Loading")

        try:
            network_file, mapping_file, decomposed_file, layout = self._resolve_layout()
            if network_file is None:
                result.fail(
                    f"No generated output found for pathway {self.pathway_id} under "
                    f"{self.output_dir}/ (looked for a per-pathway directory "
                    f"*_R-HSA-{self.pathway_id}/ and the legacy flat "
                    f"pathway_logic_network_{self.pathway_id}.csv)"
                )
                return result
            result.add_info(f"Layout: {layout}")

            self.logic_network = pd.read_csv(network_file)
            result.add_info(f"Loaded logic network: {len(self.logic_network)} edges")

            if not mapping_file.exists():
                result.fail(f"UUID mapping file not found: {mapping_file}")
                return result
            self.uuid_to_reactome = pd.read_csv(mapping_file)
            # The per-pathway layout names the column `stable_id` and carries one
            # stable id per uuid; the checks below expect the legacy
            # `entity_ids`/`reactome_id` names.
            if "stable_id" in self.uuid_to_reactome.columns:
                if "reactome_id" not in self.uuid_to_reactome.columns:
                    self.uuid_to_reactome["reactome_id"] = self.uuid_to_reactome["stable_id"]
                if "entity_ids" not in self.uuid_to_reactome.columns:
                    self.uuid_to_reactome["entity_ids"] = self.uuid_to_reactome["stable_id"]
            result.add_info(f"Loaded UUID mappings: {len(self.uuid_to_reactome)} entries")

            if decomposed_file.exists():
                self.decomposed_uid_mapping = pd.read_csv(decomposed_file)
                result.add_info(
                    f"Loaded decomposed mappings: {len(self.decomposed_uid_mapping)} entries"
                )
            else:
                # Optional: the cache dir may have been pruned. Checks that need
                # it should guard on None rather than the whole run failing.
                self.decomposed_uid_mapping = None
                result.warn(f"Decomposed mapping not found ({decomposed_file}); "
                            f"checks that need it will be skipped")

        except Exception as e:
            result.fail(f"Error loading files: {str(e)}")

        return result

    def _mapping_uses_stid(self) -> bool:
        """True when the uuid mapping carries Reactome stable ids, not dbIds.

        The per-pathway layout stores stIds (R-HSA-...); the legacy flat layout
        stored pipe-delimited numeric dbIds. Every comparison against Neo4j has
        to query the matching property, and parsing has to skip int().
        """
        values = [str(v) for v in self.uuid_to_reactome['entity_ids'].dropna()]
        return any(v.startswith('R-') for v in values)

    def _parse_entity_ids(self, entity_ids_str, uses_stid: bool) -> Set:
        """Split a pipe-delimited entity-id cell into ids of the right type."""
        parts = [eid for eid in str(entity_ids_str).split('|') if eid]
        return set(parts) if uses_stid else {int(eid) for eid in parts}

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
        # Kept in step with what the generator actually emits. This allowlist
        # had drifted to the original four while the generator gained
        # assembly/dissociation/depletion, so the check failed on a perfectly
        # healthy network — it never surfaced because load_files() aborted the
        # run before reaching it. Counts across all 92 catalog pathways:
        # input 181835, output 121669, dissociation 83058, assembly 61660,
        # catalyst 33531, regulator 6258, depletion 3548.
        # Every edge_type the generator can emit. Keep in sync with the
        # "edge_type" literals in src/logic_network_generator.py — a stale
        # allowlist here fails healthy networks, which is how `diagram_bridge`
        # (added with the additive diagram-connectivity work) first showed up.
        valid_edge_types = {
            'input', 'output', 'catalyst', 'regulator',
            'assembly', 'dissociation', 'depletion',
            'diagram_bridge', 'handoff',
        }
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
        pathway_match = self._pathway_match("pathway")
        positive_query = f"""
        MATCH {pathway_match}-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(regulator:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_pos_count = self.graph.run(positive_query, pathway_key=self.pathway_key).data()[0]['count']

        pathway_match = self._pathway_match("pathway")
        negative_query = f"""
        MATCH {pathway_match}-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(regulator:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_neg_count = self.graph.run(negative_query, pathway_key=self.pathway_key).data()[0]['count']

        pathway_match = self._pathway_match("pathway")
        catalyst_query = f"""
        MATCH {pathway_match}-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:catalystActivity]->(ca:CatalystActivity)
        RETURN COUNT(DISTINCT reaction) AS count
        """
        neo4j_catalyst_count = self.graph.run(catalyst_query, pathway_key=self.pathway_key).data()[0]['count']

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

        # Expansion means the network legitimately has MORE regulated virtual
        # reactions than Neo4j has regulated reactions, so a shortfall cannot
        # be an exact assertion. But a role Neo4j records must not vanish
        # entirely: that is a floor, and it is checkable.
        #
        # Every branch below used to be warn() only, so this check could not
        # fail — deleting every regulator edge in the network still gave
        # "PASS: Regulator Propagation" and 11/11. That is the same defect
        # fixed in validate_edge_counts, one method away, and missed here.
        for label, expected_count, emitted in (
            ("positive regulator", neo4j_pos_count, logic_pos_reactions),
            ("negative regulator", neo4j_neg_count, logic_neg_reactions),
            ("catalyst", neo4j_catalyst_count, logic_catalyst_reactions),
        ):
            if emitted >= expected_count:
                result.add_info(
                    f"{label.capitalize()}s: ✓ (may be duplicated for virtual reactions)"
                )
            elif emitted == 0 and expected_count > 0:
                result.fail(
                    f"Neo4j has {expected_count} reactions with a {label} but the "
                    f"logic network has none at all"
                )
            else:
                result.warn(
                    f"Missing {label}s: expected >={expected_count}, got {emitted}"
                )

        return result

    def validate_reconstruction(self) -> ValidationResult:
        """Validate that the logic network can reconstruct the original pathway."""
        result = ValidationResult("Pathway Reconstruction")

        # The uuid mapping may carry stIds (per-pathway layout) or dbIds
        # (legacy). Comparing across id spaces yields 0% matching by
        # construction — which the reconstruction floor would then report as a
        # hard failure on a perfectly healthy network.
        uses_stid = self._mapping_uses_stid()
        id_property = "stId" if uses_stid else "dbId"

        # Build UUID lookup
        uuid_dict = {}
        for _, row in self.uuid_to_reactome.iterrows():
            uuid_val = row['uuid']
            entity_ids_str = row['entity_ids']
            if pd.notna(entity_ids_str) and entity_ids_str:
                entity_ids = self._parse_entity_ids(entity_ids_str, self._mapping_uses_stid())
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
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
        OPTIONAL MATCH (r)-[:input]->(inp)
        OPTIONAL MATCH (r)-[:output]->(out)
        WITH r, collect(DISTINCT inp.{id_property}) AS inputs, collect(DISTINCT out.{id_property}) AS outputs
        RETURN r.dbId AS reaction_id,
               [x IN inputs WHERE x IS NOT NULL] AS inputs,
               [x IN outputs WHERE x IS NOT NULL] AS outputs
        """

        neo4j_reaction_pairs = set()
        reactions_data = self.graph.run(query, pathway_key=self.pathway_key).data()

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

        # A reconstruction accuracy of 0% used to report PASS, because this
        # branch only warned. Fail below a floor so catastrophic breakage is
        # caught, and keep warning in the band between the floor and 90% where a
        # shortfall may be a known structural gap rather than a defect.
        raw_min = os.environ.get("LNG_VALIDATE_MIN_RECONSTRUCTION", "50")
        try:
            min_accuracy = float(raw_min)
        except (TypeError, ValueError):
            result.warn(
                f"LNG_VALIDATE_MIN_RECONSTRUCTION={raw_min!r} is not a number; "
                f"using the default of 50."
            )
            min_accuracy = 50.0
        if accuracy == 100.0:
            result.add_info("🎉 Perfect reconstruction!")
        elif accuracy >= 90:
            result.add_info("Good reconstruction (>90%)")
        elif accuracy >= min_accuracy:
            result.warn(f"Reconstruction accuracy below 90%: {accuracy:.1f}%")
        else:
            result.fail(
                f"Reconstruction accuracy {accuracy:.1f}% is below the "
                f"{min_accuracy:.0f}% floor (LNG_VALIDATE_MIN_RECONSTRUCTION)"
            )

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
                # This check exists to reject self-loops, so finding them must
                # FAIL. Reporting them as a warning left `passed = True`, i.e.
                # "ALL VALIDATIONS PASSED" on a network full of the exact defect
                # the check is named after.
                result.fail(f"{edge_type} has {len(self_loops)} self-loops at UUID level")
                # Show examples
                for _, edge in self_loops.head(3).iterrows():
                    result.warn(f"  Example: {edge['source_id']} → {edge['target_id']}")
            else:
                result.add_info(f"{edge_type}: No self-loops ✓")

        return result

    def _set_leaf_members(self, ids: set, id_property: str) -> dict:
        """Map each EntitySet id to its non-set leaf members.

        The generator *splits* EntitySets: the set's own id never becomes a
        node, each member does. Comparing raw Neo4j ids against the network
        therefore reports every set as missing.
        """
        if not ids:
            return {}
        query = f"""
        MATCH (s:EntitySet) WHERE s.{id_property} IN $ids
        MATCH (s)-[:hasMember|hasCandidate*1..{MAX_SET_NESTING}]->(m:PhysicalEntity)
        WHERE NOT m:EntitySet
        RETURN s.{id_property} AS set_id,
               COLLECT(DISTINCT m.{id_property}) AS members
        """
        rows = self.graph.run(query, ids=list(ids)).data()
        return {row["set_id"]: set(row["members"]) for row in rows}

    def _decomposed_ids(self) -> set:
        """Ids the generator decomposed into components.

        A decomposed Complex contributes its *components* as nodes, so its own
        id is absent from stid_to_uuid_mapping.csv while still being fully
        represented. It is recorded here instead.
        """
        if self.decomposed_uid_mapping is None:
            return set()

        # Only the PARENT columns. `component_id` and
        # `input_or_output_reactome_id` hold the components themselves — every
        # participating entity in the pathway — so including them made
        # `expected - present` empty by construction and gutted the checks this
        # helper exists to support: with all 116 catalyst edges deleted from
        # R-HSA-69620, coverage reported 7 missing instead of 18. A decomposed
        # entity is one that has components, not one that IS a component.
        ids = set()
        for column in ("reactome_id", "source_entity_id"):
            if column not in self.decomposed_uid_mapping.columns:
                continue
            frame = self.decomposed_uid_mapping
            if "component_id" in frame.columns:
                # A row only attests decomposition when the component differs
                # from the parent; a self-referential row decomposes nothing.
                frame = frame[frame[column] != frame["component_id"]]
            for value in frame[column].dropna():
                text = str(value)
                if text not in ("", "None"):
                    ids.add(text)
        return ids

    def _uncovered(self, expected: set, present: set, id_property: str):
        """Expected ids not represented in the network, allowing for splitting.

        Returns ``(missing, partial)``. ``partial`` maps a set id to the
        members that did *not* survive: a set the generator split only halfway
        is a real defect, so it is reported separately rather than excused.
        """
        # ORDER MATTERS. The decomposed-Complex excuse must be applied AFTER
        # the EntitySet analysis, never before it.
        #
        # A split EntitySet is itself a decomposition parent — its stId sits in
        # decomposed_uid_mapping.reactome_id — so unioning _decomposed_ids()
        # into `present` up front excused every set before its members were
        # ever examined, took the early return, and never called
        # _set_leaf_members at all. The `partial` branch below was therefore
        # unreachable: deleting ALL seven members of R-HSA-1445138 still
        # reported "9 EntitySets represented by their members" and 11/11.
        #
        # That is the same over-broad excuse this helper was rewritten to
        # remove, one layer down: the earlier fix restored sensitivity for
        # plain entities and left sets fully masked.
        direct = expected - present
        if not direct:
            return set(), {}

        members = self._set_leaf_members(direct, id_property)
        decomposed = self._decomposed_ids()
        missing, partial = set(), {}
        for entity_id in direct:
            leaves = members.get(entity_id)
            if leaves:
                # A set is covered only when its members actually survived.
                # Being a decomposition parent does not excuse the SET — that
                # was the masking bug — but it does excuse a MEMBER: a member
                # that is itself a decomposed Complex is represented by its
                # components rather than by its own id, so it is present in
                # the only sense available to it.
                absent = leaves - present - decomposed
                if absent:
                    partial[entity_id] = absent
            elif entity_id in decomposed:
                # A decomposed Complex is represented by its components rather
                # than by its own id.
                continue
            else:
                missing.add(entity_id)
        return missing, partial

    def validate_entity_coverage(self) -> ValidationResult:
        """Verify all Neo4j entities appear in logic network."""
        result = ValidationResult("Entity Coverage")

        # Which identifier does the mapping carry? The per-pathway layout uses
        # Reactome stable ids (R-HSA-...); the legacy flat layout used
        # pipe-delimited numeric dbIds. Query whichever property matches, or
        # this compares stIds against dbIds and int() blows up on "R-HSA-...".
        raw_values = [str(v) for v in self.uuid_to_reactome['entity_ids'].dropna()]
        uses_stid = any(v.startswith('R-') for v in raw_values)
        id_property = 'stId' if uses_stid else 'dbId'

        # Parameterised: this query previously interpolated self.pathway_id and
        # was safe only because argparse coerces it with type=int.
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:input|output]->(entity:PhysicalEntity)
        RETURN COLLECT(DISTINCT entity.{id_property}) as entity_ids
        """
        neo4j_result = self.graph.run(query, pathway_key=self.pathway_key).data()
        neo4j_entities = set(neo4j_result[0]['entity_ids']) if neo4j_result else set()

        # Get all entities from logic network via uuid_to_reactome mapping
        ln_entities = set()
        for _, row in self.uuid_to_reactome.iterrows():
            entity_ids_str = row['entity_ids']
            if pd.notna(entity_ids_str):
                parts = [eid for eid in str(entity_ids_str).split('|') if eid]
                if uses_stid:
                    ln_entities.update(parts)
                else:
                    ln_entities.update(int(eid) for eid in parts)

        missing_entities, partial_sets = self._uncovered(
            neo4j_entities, ln_entities, id_property
        )
        extra_entities = ln_entities - neo4j_entities
        split_sets = (neo4j_entities - ln_entities) - missing_entities - set(partial_sets)

        result.add_info(f"Neo4j entities: {len(neo4j_entities)}")
        result.add_info(f"Logic network entities: {len(ln_entities)}")
        if split_sets:
            result.add_info(f"{len(split_sets)} EntitySets represented by their members ✓")

        if missing_entities:
            result.fail(f"Missing {len(missing_entities)} entities from Neo4j")
            for entity_id in list(missing_entities)[:5]:
                result.fail(f"  Missing entity: {entity_id}")
        elif not partial_sets:
            result.add_info("All Neo4j entities present ✓")

        if partial_sets:
            result.fail(f"{len(partial_sets)} EntitySets only partially split")
            for set_id, dropped in list(partial_sets.items())[:5]:
                result.fail(f"  {set_id}: {len(dropped)} member(s) absent, e.g. {sorted(dropped)[0]}")

        if extra_entities:
            result.add_info(f"Logic network has {len(extra_entities)} extra entities (from catalysts/regulators)")

        return result

    def validate_catalyst_completeness(self) -> ValidationResult:
        """Verify all Neo4j catalysts are present in logic network."""
        result = ValidationResult("Catalyst Completeness")

        # Get catalysts from Neo4j, using whichever identifier the uuid mapping
        # carries — comparing Neo4j dbIds against network stIds reported every
        # catalyst as "missing".
        uses_stid = self._mapping_uses_stid()
        id_property = "stId" if uses_stid else "dbId"
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:catalystActivity]->(ca)-[:physicalEntity]->(catalyst)
        RETURN COLLECT(DISTINCT catalyst.{id_property}) as catalyst_ids
        """
        neo4j_result = self.graph.run(query, pathway_key=self.pathway_key).data()
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
                    ln_catalysts.update(
                        self._parse_entity_ids(entity_ids_str, self._mapping_uses_stid())
                    )

        missing, partial_sets = self._uncovered(
            neo4j_catalysts, ln_catalysts, id_property
        )
        split_sets = (neo4j_catalysts - ln_catalysts) - missing - set(partial_sets)

        result.add_info(f"Neo4j catalysts: {len(neo4j_catalysts)}")
        result.add_info(f"Logic network catalysts: {len(ln_catalysts)}")
        if split_sets:
            result.add_info(f"{len(split_sets)} set-valued catalysts flattened to members ✓")

        if missing:
            result.fail(f"Missing {len(missing)} catalysts from Neo4j")
            for catalyst_id in list(missing)[:5]:
                result.fail(f"  Missing catalyst: {catalyst_id}")
        elif not partial_sets:
            result.add_info("All catalysts present ✓")

        if partial_sets:
            result.fail(f"{len(partial_sets)} set-valued catalysts only partially flattened")
            for set_id, dropped in list(partial_sets.items())[:5]:
                result.fail(f"  {set_id}: {len(dropped)} member(s) absent, e.g. {sorted(dropped)[0]}")

        return result

    def _entity_for_uuid(self, node_uuid, uses_stid):
        """First Reactome id a network node uuid maps to, or None."""
        mapping = self.uuid_to_reactome[self.uuid_to_reactome['uuid'] == node_uuid]
        if mapping.empty:
            return None
        parsed = self._parse_entity_ids(mapping.iloc[0]['entity_ids'], uses_stid)
        return next(iter(parsed), None)

    def validate_regulator_polarity(self) -> ValidationResult:
        """Verify regulator pos_neg values match Neo4j, per regulated reaction."""
        result = ValidationResult("Regulator Polarity")

        uses_stid = self._mapping_uses_stid()
        id_property = "stId" if uses_stid else "dbId"

        # Polarity is a property of a (regulator, reaction) pair, not of the
        # entity: Reactome legitimately curates the same entity as a positive
        # regulator of one reaction and a negative regulator of another (PI5P
        # in R-HSA-1257604). Comparing flat per-entity sets reported every such
        # dual-role regulator as wrong whatever the generator emitted.
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:regulatedBy]->(reg)-[:regulator]->(pe)
        WHERE reg:PositiveRegulation OR reg:NegativeRegulation
        RETURN r.{id_property} AS reaction, pe.{id_property} AS regulator,
               CASE WHEN reg:PositiveRegulation THEN 'pos' ELSE 'neg' END AS polarity
        """
        rows = self.graph.run(query, pathway_key=self.pathway_key).data()

        # A set-valued regulator is flattened onto the reaction as its members,
        # so accept any leaf member wherever Neo4j names the set itself.
        leaves = self._set_leaf_members({row["regulator"] for row in rows}, id_property)
        expected = {}
        for row in rows:
            for regulator_id in {row["regulator"]} | leaves.get(row["regulator"], set()):
                expected.setdefault((row["reaction"], regulator_id), set()).add(row["polarity"])

        dual_role = {r for (_, r), pol in expected.items() if len(pol) > 1}

        regulator_edges = self.logic_network[self.logic_network['edge_type'] == 'regulator']
        mismatches = []
        unattested = 0
        checked_count = 0

        for _, edge in regulator_edges.iterrows():
            entity_id = self._entity_for_uuid(edge['source_id'], uses_stid)
            reaction_id = self._entity_for_uuid(edge['target_id'], uses_stid)
            if entity_id is None or reaction_id is None:
                continue

            polarities = expected.get((reaction_id, entity_id))
            if polarities is None:
                # Neo4j records no regulation of this reaction by this entity —
                # a provenance gap, not a polarity error.
                unattested += 1
                continue

            checked_count += 1
            if edge['pos_neg'] not in polarities:
                mismatches.append(
                    f"{entity_id} on {reaction_id}: emitted {edge['pos_neg']}, "
                    f"Neo4j has {'/'.join(sorted(polarities))}"
                )

        result.add_info(f"Checked {checked_count} regulator edges against Neo4j")
        result.add_info(f"Neo4j: {len(rows)} regulator/reaction pairs")
        if dual_role:
            result.add_info(f"{len(dual_role)} entity/reaction pairs are dual-role (pos and neg)")

        if mismatches:
            result.fail(f"{len(mismatches)} regulator edges with wrong polarity")
            for detail in mismatches[:5]:
                result.fail(f"  {detail}")
        else:
            result.add_info("All regulator polarities correct \u2713")

        if unattested:
            result.warn(f"{unattested} regulator edges have no matching Neo4j regulation")

        return result

    def validate_reaction_coverage(self) -> ValidationResult:
        """Verify all Neo4j reactions are represented in logic network."""
        result = ValidationResult("Reaction Coverage")

        # Compare reaction IDENTITIES, not counts. Counting virtual reactions
        # against Neo4j reactions can only ever take the "extra VRs" branch,
        # because VR expansion always inflates the total — so a genuinely
        # dropped reaction was invisible. R-HSA-170834 loses 8 of its 100
        # reactions (the whole SMAD/receptor degradation arm, which have an
        # input but no output) and this check still passed.
        uses_stid = self._mapping_uses_stid()
        id_property = "stId" if uses_stid else "dbId"
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
        RETURN COLLECT(DISTINCT r.{id_property}) AS reaction_ids
        """
        rows = self.graph.run(query, pathway_key=self.pathway_key).data()
        neo4j_reactions = set(rows[0]["reaction_ids"]) if rows else set()

        # The reaction each edge came from, where the generator recorded it.
        ln_reactions = set()
        if "edge_reaction_id" in self.logic_network.columns:
            for value in self.logic_network["edge_reaction_id"].dropna():
                text = str(value)
                if text not in ("", "None"):
                    ln_reactions.add(text)

        # Fall back to the virtual-reaction nodes' own mapping when the column
        # is absent or unpopulated (it is NaN on regulator edges).
        if not ln_reactions:
            input_edges = self.logic_network[self.logic_network["edge_type"] == "input"]
            for uuid in input_edges["target_id"].unique():
                entity_id = self._entity_for_uuid(uuid, uses_stid)
                if entity_id is not None:
                    ln_reactions.add(entity_id)

        virtual_reactions = self.logic_network[
            self.logic_network["edge_type"] == "input"
        ]["target_id"].nunique()

        result.add_info(f"Neo4j reactions: {len(neo4j_reactions)}")
        result.add_info(f"Logic network reactions: {len(ln_reactions)}")
        result.add_info(f"Virtual reactions: {virtual_reactions}")

        missing = neo4j_reactions - ln_reactions
        if missing:
            result.fail(f"Missing {len(missing)} reactions from Neo4j")
            for reaction_id in sorted(missing)[:5]:
                result.fail(f"  Missing reaction: {reaction_id}")
        else:
            result.add_info("All Neo4j reactions represented \u2713")

        extra = ln_reactions - neo4j_reactions
        if extra:
            result.warn(f"{len(extra)} reactions in the network are not in this pathway")

        return result

    def validate_edge_counts(self) -> ValidationResult:
        """Compare edge counts with Neo4j."""
        result = ValidationResult("Edge Count Verification")

        # Query Neo4j for unique entity counts per edge type
        p_match = self._pathway_match("p")
        query = f"""
        MATCH {p_match}-[:hasEvent*]->(r:ReactionLikeEvent)
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

        neo4j_result = self.graph.run(query, pathway_key=self.pathway_key).data()
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

        # This check emitted only add_info, so it was structurally incapable of
        # failing. Expansion means the network legitimately has MORE edges than
        # Neo4j has distinct entities, so the only sound assertion is the floor:
        # a role Neo4j records must not vanish entirely.
        result.add_info("Logic network may have more edges than Neo4j entities (expansion)")

        for role, ln_count, neo4j_key in (
            ("input", ln_inputs, "input_count"),
            ("output", ln_outputs, "output_count"),
            ("catalyst", ln_catalysts, "catalyst_count"),
            ("regulator", ln_regulators, "regulator_count"),
        ):
            expected = neo4j_counts.get(neo4j_key, 0) or 0
            if expected > 0 and ln_count == 0:
                result.fail(
                    f"Neo4j has {expected} {role} entities but the logic network "
                    f"has no {role} edges at all"
                )

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
        required=True,
        help="Reactome pathway stable id (R-HSA-69620). A numeric dbId is "
             "accepted for back-compat but the stable id is preferred.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory holding the generated network (default: ./output)",
    )

    args = parser.parse_args()

    validator = LogicNetworkValidator(args.pathway_id, args.output_dir)
    success = validator.run_all_validations()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
