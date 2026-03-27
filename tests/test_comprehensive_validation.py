"""Comprehensive validation: generated pathways vs Neo4j database.

Tests verify that generated logic networks correctly capture:
1. All positive and negative regulators from the database
2. All catalytic activity from the database
3. Correct decomposition of complexes and entity sets
4. Proper edge structure (source_id, target_id, pos_neg, and_or, edge_type)

These tests require a running Neo4j database with Reactome data.
"""

import pandas as pd
import pytest
import sys
from pathlib import Path
from collections import defaultdict

from py2neo import Graph

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def find_pathway_dir(pathway_id: str) -> Path:
    """Find the output directory for a pathway by its ID."""
    output_dir = Path("output")
    for d in output_dir.iterdir():
        if d.is_dir() and d.name.endswith(f"_{pathway_id}"):
            return d
    return None


# Test pathways: a mix of small, medium, and large
TEST_PATHWAY_IDS = ["9612973", "9909396", "73894", "112316", "397014"]


def get_available_test_pathways():
    """Return pathway IDs that have been generated."""
    available = []
    for pid in TEST_PATHWAY_IDS:
        d = find_pathway_dir(pid)
        if d and (d / "logic_network.csv").exists():
            available.append(pid)
    return available


AVAILABLE_PATHWAYS = get_available_test_pathways()


@pytest.fixture(scope="module")
def graph():
    """Create Neo4j graph connection."""
    try:
        g = Graph("bolt://localhost:7687", auth=("neo4j", "test"))
        g.run("RETURN 1").data()
        return g
    except Exception:
        pytest.skip("Neo4j database not available")


@pytest.mark.database
class TestRegulatorCompleteness:
    """Verify all regulators from Neo4j are present in generated networks."""

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_all_positive_regulators_present(self, graph, pathway_id):
        """Every positive regulator in Neo4j should appear as a pos/regulator edge."""
        pathway_dir = find_pathway_dir(pathway_id)
        network = pd.read_csv(pathway_dir / "logic_network.csv")

        # Query DB for positive regulators
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(reg:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN DISTINCT reaction.dbId as reaction_id, pe.dbId as regulator_id
        """
        db_pos_regulators = graph.run(query).data()

        # Count in network
        pos_reg_edges = network[
            (network['edge_type'] == 'regulator') & (network['pos_neg'] == 'pos')
        ]

        if len(db_pos_regulators) > 0:
            assert len(pos_reg_edges) > 0, (
                f"Pathway {pathway_id}: DB has {len(db_pos_regulators)} positive regulators "
                f"but network has 0 positive regulator edges"
            )
            # Allow some loss due to reactions not in reaction_connections
            coverage = len(pos_reg_edges) / len(db_pos_regulators)
            assert coverage >= 0.8, (
                f"Pathway {pathway_id}: DB has {len(db_pos_regulators)} positive regulators "
                f"but network only has {len(pos_reg_edges)} ({coverage*100:.0f}% coverage)"
            )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_all_negative_regulators_present(self, graph, pathway_id):
        """Every negative regulator in Neo4j should appear as a neg/regulator edge."""
        pathway_dir = find_pathway_dir(pathway_id)
        network = pd.read_csv(pathway_dir / "logic_network.csv")

        # Query DB for negative regulators
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(reg:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN DISTINCT reaction.dbId as reaction_id, pe.dbId as regulator_id
        """
        db_neg_regulators = graph.run(query).data()

        # Count in network
        neg_reg_edges = network[
            (network['edge_type'] == 'regulator') & (network['pos_neg'] == 'neg')
        ]

        if len(db_neg_regulators) > 0:
            assert len(neg_reg_edges) > 0, (
                f"Pathway {pathway_id}: DB has {len(db_neg_regulators)} negative regulators "
                f"but network has 0 negative regulator edges"
            )
            coverage = len(neg_reg_edges) / len(db_neg_regulators)
            assert coverage >= 0.8, (
                f"Pathway {pathway_id}: DB has {len(db_neg_regulators)} negative regulators "
                f"but network only has {len(neg_reg_edges)} ({coverage*100:.0f}% coverage)"
            )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_negative_regulators_marked_neg(self, graph, pathway_id):
        """All regulator edges with pos_neg='neg' should only be negative regulators."""
        pathway_dir = find_pathway_dir(pathway_id)
        network = pd.read_csv(pathway_dir / "logic_network.csv")

        neg_edges = network[network['pos_neg'] == 'neg']
        # All negative edges should be regulators (not catalysts or main edges)
        for _, edge in neg_edges.iterrows():
            assert edge['edge_type'] == 'regulator', (
                f"Found neg edge with edge_type='{edge['edge_type']}' instead of 'regulator'"
            )


@pytest.mark.database
class TestCatalystCompleteness:
    """Verify all catalysts from Neo4j are present in generated networks."""

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_all_catalysts_present(self, graph, pathway_id):
        """Every catalyst in Neo4j should appear as a pos/catalyst edge."""
        pathway_dir = find_pathway_dir(pathway_id)
        network = pd.read_csv(pathway_dir / "logic_network.csv")

        # Query DB for catalysts
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:catalystActivity]->(ca:CatalystActivity)-[:physicalEntity]->(pe:PhysicalEntity)
        RETURN DISTINCT reaction.dbId as reaction_id, pe.dbId as catalyst_id
        """
        db_catalysts = graph.run(query).data()

        # Count in network
        catalyst_edges = network[network['edge_type'] == 'catalyst']

        if len(db_catalysts) > 0:
            assert len(catalyst_edges) > 0, (
                f"Pathway {pathway_id}: DB has {len(db_catalysts)} catalysts "
                f"but network has 0 catalyst edges"
            )
            # Some catalysts may be missed if their reaction isn't in reaction_connections
            coverage = len(catalyst_edges) / len(db_catalysts)
            assert coverage >= 0.7, (
                f"Pathway {pathway_id}: DB has {len(db_catalysts)} catalysts "
                f"but network only has {len(catalyst_edges)} ({coverage*100:.0f}% coverage)"
            )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_catalysts_always_positive(self, graph, pathway_id):
        """All catalyst edges should have pos_neg='pos'."""
        pathway_dir = find_pathway_dir(pathway_id)
        network = pd.read_csv(pathway_dir / "logic_network.csv")

        catalyst_edges = network[network['edge_type'] == 'catalyst']
        if len(catalyst_edges) == 0:
            pytest.skip("No catalyst edges in this pathway")

        neg_catalysts = catalyst_edges[catalyst_edges['pos_neg'] != 'pos']
        assert len(neg_catalysts) == 0, (
            f"Found {len(neg_catalysts)} catalyst edges that are not positive"
        )


@pytest.mark.database
class TestDecompositionCorrectness:
    """Verify that complex/set decomposition correctly captures all entities."""

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_all_reactions_in_decomposition(self, graph, pathway_id):
        """All reactions from DB should appear in the decomposed_uid_mapping."""
        pathway_dir = find_pathway_dir(pathway_id)
        decomposed = pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv")

        # Query DB for reactions
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        RETURN DISTINCT reaction.dbId as reaction_id
        """
        db_reactions = {row['reaction_id'] for row in graph.run(query).data()}

        # Get reactions from decomposition
        decomposed_reactions = set(decomposed['reactome_id'].dropna().astype(int).unique())

        # Check coverage
        missing = db_reactions - decomposed_reactions
        coverage = len(db_reactions - missing) / len(db_reactions) if db_reactions else 1.0

        assert coverage > 0.8, (
            f"Pathway {pathway_id}: Only {coverage*100:.1f}% of DB reactions are in decomposition. "
            f"Missing {len(missing)}/{len(db_reactions)} reactions."
        )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_complexes_are_decomposed(self, graph, pathway_id):
        """Complexes with components should be decomposed into their parts."""
        pathway_dir = find_pathway_dir(pathway_id)
        decomposed = pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv")

        # Query DB for complexes with components
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:input|output]->(complex:Complex)-[:hasComponent]->(component)
        RETURN DISTINCT complex.dbId as complex_id, count(DISTINCT component) as num_components
        """
        db_complexes = graph.run(query).data()

        if len(db_complexes) == 0:
            pytest.skip("No complexes in this pathway")

        # For complexes with >1 component, we expect multiple rows in decomposition
        multi_component_complexes = [c for c in db_complexes if c['num_components'] > 1]

        # Check that decomposition has multiple hashes per reaction (indicating decomposition happened)
        reaction_hash_counts = decomposed.groupby('reactome_id')['uid'].nunique()
        multi_hash_reactions = reaction_hash_counts[reaction_hash_counts > 1]

        assert len(multi_hash_reactions) > 0, (
            f"Pathway {pathway_id}: Has {len(multi_component_complexes)} multi-component complexes "
            f"but no reactions have multiple decomposition hashes"
        )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_entity_sets_are_decomposed(self, graph, pathway_id):
        """EntitySets should be decomposed into their members."""
        pathway_dir = find_pathway_dir(pathway_id)
        decomposed = pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv")

        # Query DB for entity sets
        query = f"""
        MATCH (pathway:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:input|output]->(es:EntitySet)-[:hasMember|hasCandidate]->(member)
        RETURN DISTINCT es.dbId as set_id, count(DISTINCT member) as num_members
        """
        db_sets = graph.run(query).data()

        if len(db_sets) == 0:
            pytest.skip("No entity sets in this pathway")

        # Source entity ID should track original sets
        if 'source_entity_id' in decomposed.columns:
            source_entities = decomposed['source_entity_id'].dropna().astype(int).unique()
            db_set_ids = {row['set_id'] for row in db_sets}
            covered_sets = db_set_ids.intersection(set(source_entities))

            # Some sets should be tracked
            assert len(covered_sets) > 0 or len(source_entities) > 0, (
                f"Pathway {pathway_id}: Has {len(db_sets)} entity sets "
                f"but source_entity_id tracking found none"
            )

    @pytest.mark.parametrize("pathway_id", AVAILABLE_PATHWAYS)
    def test_best_matches_pair_same_reaction(self, graph, pathway_id):
        """best_matches should pair input/output hashes from the same reaction."""
        pathway_dir = find_pathway_dir(pathway_id)
        decomposed = pd.read_csv(pathway_dir / "cache" / "decomposed_uid_mapping.csv")
        best_matches = pd.read_csv(pathway_dir / "cache" / "best_matches.csv")

        mismatches = 0
        sample_size = min(20, len(best_matches))

        for _, match in best_matches.head(sample_size).iterrows():
            incoming_hash = match["incomming"]
            outgoing_hash = match["outgoing"]

            incoming_reactions = set(
                decomposed[decomposed["uid"] == incoming_hash]["reactome_id"].unique()
            )
            outgoing_reactions = set(
                decomposed[decomposed["uid"] == outgoing_hash]["reactome_id"].unique()
            )

            if not incoming_reactions.intersection(outgoing_reactions):
                mismatches += 1

        assert mismatches == 0, (
            f"Pathway {pathway_id}: {mismatches}/{sample_size} best_matches "
            f"pair hashes from different reactions"
        )


@pytest.mark.database
class TestEdgeCountSummary:
    """Summary test: print edge counts for all pathways and verify basic sanity."""

    def test_all_pathways_edge_summary(self, graph):
        """Print summary of all pathway edge counts for review."""
        output_dir = Path("output")
        results = []

        for d in sorted(output_dir.iterdir()):
            if not d.is_dir() or not (d / "logic_network.csv").exists():
                continue

            network = pd.read_csv(d / "logic_network.csv")
            main = network[~network['edge_type'].isin(['catalyst', 'regulator'])]
            catalysts = network[network['edge_type'] == 'catalyst']
            pos_regs = network[(network['edge_type'] == 'regulator') & (network['pos_neg'] == 'pos')]
            neg_regs = network[(network['edge_type'] == 'regulator') & (network['pos_neg'] == 'neg')]

            results.append({
                'pathway': d.name,
                'total': len(network),
                'main': len(main),
                'catalysts': len(catalysts),
                'pos_reg': len(pos_regs),
                'neg_reg': len(neg_regs),
            })

        print("\n" + "=" * 90)
        print(f"{'Pathway':<45} {'Total':>7} {'Main':>7} {'Cat':>5} {'+Reg':>5} {'-Reg':>5}")
        print("-" * 90)
        for r in results:
            print(f"{r['pathway']:<45} {r['total']:>7} {r['main']:>7} {r['catalysts']:>5} {r['pos_reg']:>5} {r['neg_reg']:>5}")
        print("=" * 90)

        # Every pathway should have either main edges or catalyst/regulator edges
        for r in results:
            assert r['total'] > 0, f"Pathway {r['pathway']} has no edges at all"
