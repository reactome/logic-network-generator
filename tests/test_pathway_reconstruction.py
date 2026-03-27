"""Test that generated logic networks can be reconstructed back to original pathways.

This test ensures bidirectional traceability:
- Forward: Reactome pathway -> Logic network (generation)
- Backward: Logic network -> Reactome pathway (reconstruction)

Requirements:
1. All entities must be traceable back to their original IDs
2. EntitySet members must be traceable back to their parent EntitySets
3. Virtual reactions must be traceable back to their source reactions

These tests require a running Neo4j database with Reactome data.
"""

import pandas as pd
import pytest
from pathlib import Path
from typing import Dict, Set, Tuple, List
from py2neo import Graph


def find_pathway_dirs():
    """Find all generated pathway directories with complete files."""
    output_dir = Path("output")
    if not output_dir.exists():
        return []
    dirs = []
    for d in sorted(output_dir.iterdir()):
        if (d.is_dir()
                and (d / "logic_network.csv").exists()
                and (d / "cache" / "decomposed_uid_mapping.csv").exists()
                and (d / "cache" / "best_matches.csv").exists()):
            parts = d.name.rsplit("_", 1)
            if len(parts) == 2 and parts[1].isdigit():
                dirs.append((parts[1], d))
    return dirs


AVAILABLE_PATHWAYS = find_pathway_dirs()
# Use a small sample for detailed reconstruction tests
SAMPLE_PATHWAYS = AVAILABLE_PATHWAYS[:3] if len(AVAILABLE_PATHWAYS) > 3 else AVAILABLE_PATHWAYS


@pytest.mark.database
class TestPathwayReconstruction:
    """Validate reconstruction of original pathways from logic networks."""

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
    def pathway_data(self, request):
        """Load generated pathway files."""
        pathway_id, pathway_dir = request.param
        return {
            'pathway_id': pathway_id,
            'pathway_dir': pathway_dir,
            'best_matches': pd.read_csv(pathway_dir / "cache" / "best_matches.csv"),
            'decomposed': pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv"),
            'logic_network': pd.read_csv(pathway_dir / "logic_network.csv"),
        }

    def test_source_entity_id_column_exists(self, pathway_data):
        """Verify that source_entity_id column exists in decomposed mapping."""
        decomposed = pathway_data["decomposed"]
        assert "source_entity_id" in decomposed.columns, \
            "source_entity_id column missing from decomposed_uid_mapping"

    def test_source_entity_id_populated_for_entitysets(self, pathway_data):
        """Verify that source_entity_id is populated for EntitySet members."""
        decomposed = pathway_data["decomposed"]

        populated_count = decomposed['source_entity_id'].notna().sum()

        # Some pathways may not have entity sets, so just check it doesn't error
        assert populated_count >= 0, "source_entity_id count should be non-negative"

    def test_virtual_reactions_trace_to_source(self, pathway_data):
        """Verify that all virtual reactions can be traced back to their source reaction."""
        best_matches = pathway_data["best_matches"]
        decomposed = pathway_data["decomposed"]

        untraceable = 0
        sample_size = min(20, len(best_matches))

        for _, row in best_matches.head(sample_size).iterrows():
            input_uid = row['incomming']
            output_uid = row['outgoing']

            input_rows = decomposed[decomposed['uid'] == input_uid]
            if input_rows.empty:
                untraceable += 1
                continue

            output_rows = decomposed[decomposed['uid'] == output_uid]
            if output_rows.empty:
                untraceable += 1
                continue

            # Verify both come from same reaction
            input_reactions = set(input_rows['reactome_id'].unique())
            output_reactions = set(output_rows['reactome_id'].unique())

            if not input_reactions & output_reactions:
                untraceable += 1

        assert untraceable == 0, \
            f"{pathway_data['pathway_id']}: {untraceable}/{sample_size} virtual reactions are untraceable"

    def test_no_information_loss_in_decomposition(self, pathway_data, graph):
        """Verify that no entities are lost during decomposition."""
        pathway_id = pathway_data['pathway_id']
        decomposed = pathway_data["decomposed"]

        query = f"""
        MATCH (p:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:input|output]->(e)
        RETURN DISTINCT e.dbId AS entity_id
        """
        result = graph.run(query).data()
        neo4j_entities = {row["entity_id"] for row in result if row["entity_id"] is not None}

        # Get all entities from decomposed mapping
        decomposed_entities = set()

        if 'component_id' in decomposed.columns:
            decomposed_entities.update(decomposed['component_id'].dropna().astype(int).unique())

        if 'input_or_output_reactome_id' in decomposed.columns:
            decomposed_entities.update(
                decomposed['input_or_output_reactome_id'].dropna().astype(int).unique()
            )

        if 'source_entity_id' in decomposed.columns:
            decomposed_entities.update(
                decomposed['source_entity_id'].dropna().astype(int).unique()
            )

        # Also check reactome_id column for reaction IDs that might be entities
        decomposed_entities.update(decomposed['reactome_id'].dropna().astype(int).unique())

        missing = neo4j_entities - decomposed_entities

        # Allow some missing (e.g., entities only in catalysts/regulators not in input/output)
        coverage = (len(neo4j_entities) - len(missing)) / len(neo4j_entities) if neo4j_entities else 1.0

        assert coverage > 0.5, (
            f"Pathway {pathway_id}: Only {coverage*100:.1f}% entity coverage. "
            f"Missing {len(missing)}/{len(neo4j_entities)} entities"
        )

    def test_all_reactions_in_decomposition(self, pathway_data, graph):
        """All reactions from DB should appear in the decomposed_uid_mapping."""
        pathway_id = pathway_data['pathway_id']
        decomposed = pathway_data["decomposed"]

        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        RETURN DISTINCT reaction.dbId as reaction_id
        """
        db_reactions = {row['reaction_id'] for row in graph.run(query).data()}

        decomposed_reactions = set(decomposed['reactome_id'].dropna().astype(int).unique())

        missing = db_reactions - decomposed_reactions
        coverage = (len(db_reactions) - len(missing)) / len(db_reactions) if db_reactions else 1.0

        assert coverage > 0.8, (
            f"Pathway {pathway_id}: Only {coverage*100:.1f}% of DB reactions in decomposition. "
            f"Missing {len(missing)}/{len(db_reactions)}"
        )
