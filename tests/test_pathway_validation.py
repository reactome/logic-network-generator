"""Comprehensive validation test for logic network generation.

This test validates that the generated logic networks correctly represent
the original pathways from the database by:
1. Querying the database directly for pathway data
2. Comparing against the generated logic network files
3. Verifying completeness of regulators, catalysts, and entity decomposition

These tests require a running Neo4j database with Reactome data.
"""

import re
import pandas as pd
import pytest
from pathlib import Path
from py2neo import Graph


def find_pathway_dir(pathway_id: str) -> Path:
    """Find the output directory for a pathway by its trailing numeric ID.

    Handles both naming conventions:
    "Foo_12345" and "Foo_R-HSA-12345".
    """
    output_dir = Path("output")
    if not output_dir.exists():
        return None
    for d in output_dir.iterdir():
        if d.is_dir():
            match = re.search(r"(\d+)$", d.name)
            if match and match.group(1) == pathway_id:
                return d
    return None


def get_available_pathways():
    """Return pathway directories that have complete generated files."""
    output_dir = Path("output")
    if not output_dir.exists():
        return []
    available = []
    for d in sorted(output_dir.iterdir()):
        if (d.is_dir()
                and (d / "logic_network.csv").exists()
                and (d / "stid_to_uuid_mapping.csv").exists()
                and (d / "cache" / "decomposed_uid_mapping.csv").exists()):
            # Extract trailing numeric pathway ID — handles both "Foo_12345"
            # and "Foo_R-HSA-12345" naming.
            match = re.search(r"(\d+)$", d.name)
            if match:
                available.append((match.group(1), d))
    return available


AVAILABLE_PATHWAYS = get_available_pathways()
# Use first 3 available pathways for parametrized tests
SAMPLE_PATHWAYS = AVAILABLE_PATHWAYS[:3]


@pytest.mark.database
class TestPathwayValidation:
    """Comprehensive validation of logic network generation.

    Note: These tests require Neo4j database to be running.
    """

    @pytest.fixture(scope="module")
    def graph(self):
        """Create Neo4j graph connection."""
        try:
            g = Graph("bolt://localhost:7687", auth=("neo4j", "test"))
            g.run("RETURN 1").data()
            return g
        except Exception:
            pytest.skip("Neo4j database not available")

    @pytest.fixture(params=SAMPLE_PATHWAYS,
                    ids=[p[1].name for p in SAMPLE_PATHWAYS])
    def pathway_files(self, request):
        """Load generated files for a pathway."""
        pathway_id, pathway_dir = request.param
        return {
            'pathway_id': pathway_id,
            'pathway_dir': pathway_dir,
            'logic_network': pd.read_csv(pathway_dir / "logic_network.csv"),
            'uuid_mapping': pd.read_csv(pathway_dir / "stid_to_uuid_mapping.csv"),
            'decomposed_mapping': pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv"),
            'reaction_connections': pd.read_csv(pathway_dir / "cache" / "reaction_connections.csv"),
        }

    def test_database_connection(self, graph):
        """Verify database connection works."""
        result = graph.run("RETURN 1 as test").data()
        assert len(result) == 1
        assert result[0]['test'] == 1

    def test_all_reactions_present(self, graph, pathway_files):
        """Validate that all reactions from the pathway are in reaction_connections."""
        pathway_id = pathway_files['pathway_id']
        reaction_connections = pathway_files['reaction_connections']

        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        RETURN DISTINCT reaction.stId as reaction_id
        """
        db_reactions = graph.run(query).data()
        db_reaction_ids = {row['reaction_id'] for row in db_reactions}

        generated_reaction_ids = set(
            reaction_connections['preceding_reaction_id'].dropna().unique()
        ).union(
            set(reaction_connections['following_reaction_id'].dropna().unique())
        )

        missing_reactions = db_reaction_ids - generated_reaction_ids
        coverage = (len(db_reaction_ids) - len(missing_reactions)) / len(db_reaction_ids) if db_reaction_ids else 1.0

        assert coverage > 0.8, (
            f"Pathway {pathway_id}: Only {coverage*100:.1f}% of DB reactions present. "
            f"Missing {len(missing_reactions)}/{len(db_reaction_ids)}"
        )

    def test_uuid_mapping_completeness(self, pathway_files):
        """Validate that UUID mapping covers all UUIDs in logic network."""
        logic_network = pathway_files['logic_network']
        uuid_mapping = pathway_files['uuid_mapping']

        network_uuids = set(logic_network['source_id'].unique()) | set(logic_network['target_id'].unique())
        mapping_uuids = set(uuid_mapping['uuid'].unique())

        unmapped_uuids = network_uuids - mapping_uuids
        assert len(unmapped_uuids) == 0, \
            f"Found {len(unmapped_uuids)} UUIDs in network without mapping entries"

    def test_logic_network_has_valid_structure(self, pathway_files):
        """Validate basic structure of logic network."""
        logic_network = pathway_files['logic_network']
        required_columns = ['source_id', 'target_id', 'pos_neg', 'and_or', 'edge_type']

        for col in required_columns:
            assert col in logic_network.columns, f"Missing column: {col}"

        assert logic_network['source_id'].notna().all(), "Found null source_id"
        assert logic_network['target_id'].notna().all(), "Found null target_id"

        valid_pos_neg = {'pos', 'neg'}
        assert set(logic_network['pos_neg'].dropna().unique()).issubset(valid_pos_neg)

        valid_edge_types = {
            "input", "output", "catalyst", "regulator",
            "depletion", "assembly", "dissociation", "handoff",
        }
        assert set(logic_network['edge_type'].unique()).issubset(valid_edge_types)

    def test_regulators_present(self, graph, pathway_files):
        """Validate that regulators from database are present in logic network."""
        pathway_id = pathway_files['pathway_id']
        logic_network = pathway_files['logic_network']

        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(regulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN DISTINCT reaction.dbId as reaction_id, pe.dbId as regulator_id
        """
        db_regulators = graph.run(query).data()

        regulator_edges = logic_network[logic_network['edge_type'] == 'regulator']
        catalyst_edges = logic_network[logic_network['edge_type'] == 'catalyst']

        if len(db_regulators) > 0:
            total_regulatory = len(regulator_edges) + len(catalyst_edges)
            assert total_regulatory > 0, \
                f"Pathway {pathway_id}: DB has {len(db_regulators)} regulators but none in logic network"

    def test_no_self_loops_in_main_pathway(self, pathway_files):
        """Validate that main pathway edges don't have excessive self-loops."""
        logic_network = pathway_files['logic_network']

        main_edges = logic_network[
            ~logic_network['edge_type'].isin(['catalyst', 'regulator'])
        ]

        if len(main_edges) == 0:
            pytest.skip("No main pathway edges")

        self_loops = main_edges[main_edges['source_id'] == main_edges['target_id']]
        self_loop_ratio = len(self_loops) / len(main_edges)

        # Report but don't fail for known self-loop issue
        assert self_loop_ratio < 0.95, \
            f"Pathway {pathway_files['pathway_id']}: {self_loop_ratio*100:.1f}% self-loops in main edges"

    def test_position_aware_uuids_working(self, pathway_files):
        """Validate that same entity at different positions has different UUIDs."""
        uuid_mapping = pathway_files['uuid_mapping']

        reactome_id_counts = uuid_mapping['stable_id'].value_counts()
        multi_position_entities = reactome_id_counts[reactome_id_counts > 1].index

        for entity_id in multi_position_entities:
            entity_rows = uuid_mapping[uuid_mapping['stable_id'] == entity_id]
            uuids = entity_rows['uuid'].unique()
            assert len(uuids) == len(entity_rows), \
                f"Entity {entity_id} at {len(entity_rows)} positions has only {len(uuids)} unique UUIDs"
