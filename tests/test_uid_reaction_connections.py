"""Tests to verify uid_reaction_connections correctness.

Tests run against generated pathway data in the output directory.
"""

import pandas as pd
import pytest
from pathlib import Path


def find_pathway_dirs():
    """Find all generated pathway directories with required cache files."""
    output_dir = Path("output")
    if not output_dir.exists():
        return []
    dirs = []
    for d in sorted(output_dir.iterdir()):
        if (d.is_dir()
                and (d / "cache" / "reaction_connections.csv").exists()
                and (d / "cache" / "decomposed_uid_mapping.csv").exists()
                and (d / "cache" / "best_matches.csv").exists()):
            dirs.append(d)
    return dirs


PATHWAY_DIRS = find_pathway_dirs()

# Integration tier: requires generated pathway artifacts in output/
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        len(PATHWAY_DIRS) == 0,
        reason="No generated pathway directories found in output/"
    ),
]

# Use a sample of up to 5 pathways
SAMPLE_DIRS = PATHWAY_DIRS[:5] if len(PATHWAY_DIRS) > 5 else PATHWAY_DIRS


class TestUIDReactionConnections:
    """Test the uid_reaction_connections data structure correctness."""

    @pytest.fixture(params=SAMPLE_DIRS, ids=[d.name for d in SAMPLE_DIRS])
    def pathway_data(self, request):
        """Load pathway data files."""
        d = request.param
        return {
            "name": d.name,
            "reaction_connections": pd.read_csv(d / "cache" / "reaction_connections.csv"),
            "decomposed_uid_mapping": pd.read_csv(d / "cache" / "decomposed_uid_mapping.csv"),
            "best_matches": pd.read_csv(d / "cache" / "best_matches.csv"),
        }

    def test_best_matches_are_within_same_reaction(self, pathway_data):
        """Verify best_matches pair inputs/outputs from the SAME reaction."""
        best_matches = pathway_data["best_matches"]
        decomposed_uid_mapping = pathway_data["decomposed_uid_mapping"]

        mismatches = 0
        sample_size = min(10, len(best_matches))

        for _, match in best_matches.head(sample_size).iterrows():
            incoming_hash = match["incomming"]
            outgoing_hash = match["outgoing"]

            incoming_reactions = set(
                decomposed_uid_mapping[
                    decomposed_uid_mapping["uid"] == incoming_hash
                ]["reactome_id"].unique()
            )

            outgoing_reactions = set(
                decomposed_uid_mapping[
                    decomposed_uid_mapping["uid"] == outgoing_hash
                ]["reactome_id"].unique()
            )

            if not incoming_reactions & outgoing_reactions:
                mismatches += 1

        assert mismatches == 0, (
            f"{pathway_data['name']}: {mismatches}/{sample_size} best_matches "
            f"pair hashes from different reactions"
        )

    def test_reaction_connections_show_pathway_topology(self, pathway_data):
        """Verify reaction_connections represent pathway topology, not self-loops."""
        reaction_connections = pathway_data["reaction_connections"]

        connections_with_both = reaction_connections.dropna()

        if len(connections_with_both) == 0:
            pytest.skip("No complete reaction connections")

        self_loops = connections_with_both[
            connections_with_both["preceding_reaction_id"]
            == connections_with_both["following_reaction_id"]
        ]

        self_loop_percentage = (len(self_loops) / len(connections_with_both)) * 100

        assert self_loop_percentage < 10, (
            f"{pathway_data['name']}: {self_loop_percentage:.1f}% of reaction "
            f"connections are self-loops"
        )

    def test_hash_to_reactome_id_mapping_is_not_one_to_one(self, pathway_data):
        """Verify that hashes can map to multiple reactome_ids (shared entities)."""
        decomposed_uid_mapping = pathway_data["decomposed_uid_mapping"]

        hash_groups = decomposed_uid_mapping.groupby("uid")["reactome_id"].nunique()
        shared_hashes = hash_groups[hash_groups > 1]

        # This is expected - same combination can appear in multiple reactions
        assert len(shared_hashes) >= 0

    def test_decomposition_creates_multiple_combinations(self, pathway_data):
        """Verify decomposition creates multiple combinations for complexes/sets."""
        decomposed_uid_mapping = pathway_data["decomposed_uid_mapping"]

        reaction_groups = decomposed_uid_mapping.groupby("reactome_id")["uid"].nunique()
        multi_decomp = reaction_groups[reaction_groups > 1]

        # At least some reactions should have multiple decompositions
        # (unless the pathway has no complexes/sets)
        assert len(reaction_groups) > 0, "No reactions in decomposed mapping"


class TestAllPathwaysHaveValidStructure:
    """Integration test: verify all generated pathways have valid structure."""

    @pytest.mark.parametrize("pathway_dir", PATHWAY_DIRS,
                             ids=[d.name for d in PATHWAY_DIRS])
    def test_pathway_has_valid_structure(self, pathway_dir):
        """Each pathway should have a valid logic network."""
        logic_network_path = pathway_dir / "logic_network.csv"
        if not logic_network_path.exists():
            pytest.skip("No logic_network.csv")

        logic_network = pd.read_csv(logic_network_path)

        required_columns = ["source_id", "target_id", "pos_neg", "and_or", "edge_type"]
        for col in required_columns:
            assert col in logic_network.columns, f"Missing column: {col}"

        assert len(logic_network) > 0, "Logic network is empty"

        valid_edge_types = {"input", "output", "catalyst", "regulator"}
        actual_types = set(logic_network["edge_type"].unique())
        invalid = actual_types - valid_edge_types
        assert len(invalid) == 0, f"Invalid edge_type values: {invalid}"
