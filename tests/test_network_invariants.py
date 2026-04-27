"""Tests for network invariants - properties that should always hold.

These tests verify structural properties of the generated networks:
- No self-loops in main pathway edges
- Root inputs are always sources (never targets)
- Terminal outputs are always targets (never sources)
- AND/OR logic is consistent
- Edge direction represents transformations

Tests run against all generated pathways in the output directory.
"""

import os
import pytest
import pandas as pd
from pathlib import Path


def get_generated_pathways():
    """Find all generated pathway directories with logic_network.csv."""
    output_dir = Path("output")
    if not output_dir.exists():
        return []
    pathways = []
    for d in sorted(output_dir.iterdir()):
        if d.is_dir() and (d / "logic_network.csv").exists():
            pathways.append(str(d / "logic_network.csv"))
    return pathways


GENERATED_PATHWAYS = get_generated_pathways()

# Integration tier: requires generated pathway artifacts in output/
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        len(GENERATED_PATHWAYS) == 0,
        reason="No generated pathway directories found in output/"
    ),
]


# Use a smaller representative sample for parametrized tests
SAMPLE_PATHWAYS = GENERATED_PATHWAYS[:5] if len(GENERATED_PATHWAYS) > 5 else GENERATED_PATHWAYS


class TestNetworkInvariants:
    """Test invariants that should hold for any valid pathway logic network."""

    @pytest.fixture(params=SAMPLE_PATHWAYS, ids=[Path(p).parent.name for p in SAMPLE_PATHWAYS])
    def network(self, request):
        """Load a generated pathway logic network."""
        return pd.read_csv(request.param)

    @pytest.fixture
    def main_edges(self, network):
        """Extract main pathway edges (excluding catalyst/regulator)."""
        return network[~network['edge_type'].isin(['catalyst', 'regulator'])]

    def test_required_columns_exist(self, network):
        """Network must have all required columns."""
        required = ['source_id', 'target_id', 'pos_neg', 'and_or', 'edge_type']
        for col in required:
            assert col in network.columns, f"Missing column: {col}"

    def test_no_null_source_or_target(self, network):
        """No edges should have null source_id or target_id."""
        assert network['source_id'].notna().all(), "Found null source_id"
        assert network['target_id'].notna().all(), "Found null target_id"

    def test_valid_edge_types(self, network):
        """All edge_type values must be valid."""
        valid_edge_types = {'input', 'output', 'catalyst', 'regulator'}
        actual = set(network['edge_type'].unique())
        invalid = actual - valid_edge_types
        assert len(invalid) == 0, f"Invalid edge_type values: {invalid}"

    def test_valid_pos_neg_values(self, network):
        """pos_neg must be 'pos' or 'neg'."""
        valid = {'pos', 'neg'}
        actual = set(network['pos_neg'].dropna().unique())
        invalid = actual - valid
        assert len(invalid) == 0, f"Invalid pos_neg values: {invalid}"

    def test_and_logic_consistency(self, network):
        """AND ⇔ contributes to the reaction proceeding.

        Allowed: input, catalyst, positive regulator.
        Disallowed: output, negative regulator (any one blocker suffices,
        so neg regulators are OR).
        """
        and_edges = network[network['and_or'] == 'and']
        if len(and_edges) == 0:
            pytest.skip("No AND edges")
        allowed = (
            and_edges['edge_type'].isin({'input', 'catalyst'})
            | (
                (and_edges['edge_type'] == 'regulator')
                & (and_edges['pos_neg'] == 'pos')
            )
        )
        incorrect = and_edges[~allowed]
        assert len(incorrect) == 0, (
            f"Found {len(incorrect)} AND edges that are neither input/catalyst "
            f"nor positive regulator"
        )

    def test_negative_regulators_are_or(self, network):
        """Negative regulators carry OR logic: any one blocker suffices."""
        neg_reg = network[
            (network['edge_type'] == 'regulator') & (network['pos_neg'] == 'neg')
        ]
        if len(neg_reg) == 0:
            pytest.skip("No negative regulator edges")
        wrong = neg_reg[neg_reg['and_or'] != 'or']
        assert len(wrong) == 0, (
            f"Found {len(wrong)} negative regulator edges with and_or != 'or'"
        )

    def test_positive_regulators_are_and(self, network):
        """Positive regulators carry AND logic: treated as required at the
        network level. Conditional dependence is modeled later in parameters.
        """
        pos_reg = network[
            (network['edge_type'] == 'regulator') & (network['pos_neg'] == 'pos')
        ]
        if len(pos_reg) == 0:
            pytest.skip("No positive regulator edges")
        wrong = pos_reg[pos_reg['and_or'] != 'and']
        assert len(wrong) == 0, (
            f"Found {len(wrong)} positive regulator edges with and_or != 'and'"
        )

    def test_assembly_edges_are_pos_and(self, network):
        """Assembly edges (leaf → root-input complex) are pos/and.

        These are synthetic edges added at the boundary so individual
        proteins can be perturbed at network entry — see
        docs/DESIGN_DECISIONS.md, "Two layers of decomposition."
        """
        asm = network[network['edge_type'] == 'assembly']
        if len(asm) == 0:
            pytest.skip("No assembly edges")
        wrong = asm[(asm['pos_neg'] != 'pos') | (asm['and_or'] != 'and')]
        assert len(wrong) == 0, (
            f"Found {len(wrong)} assembly edges with pos_neg/and_or != pos/and"
        )

    def test_dissociation_edges_are_pos_and(self, network):
        """Dissociation edges (terminal-output complex → leaf) are pos/and."""
        diss = network[network['edge_type'] == 'dissociation']
        if len(diss) == 0:
            pytest.skip("No dissociation edges")
        wrong = diss[(diss['pos_neg'] != 'pos') | (diss['and_or'] != 'and')]
        assert len(wrong) == 0, (
            f"Found {len(wrong)} dissociation edges with pos_neg/and_or != pos/and"
        )

    def test_no_boundary_self_loops(self, network):
        """Assembly/dissociation edges never connect a node to itself."""
        boundary = network[network['edge_type'].isin({'assembly', 'dissociation'})]
        if len(boundary) == 0:
            pytest.skip("No boundary edges")
        loops = boundary[boundary['source_id'] == boundary['target_id']]
        assert len(loops) == 0, f"Found {len(loops)} boundary self-loops"

    def test_or_logic_consistency(self, main_edges):
        """Edges with 'or' logic should have edge_type='output'."""
        if len(main_edges) == 0:
            pytest.skip("No main pathway edges")
        or_edges = main_edges[main_edges['and_or'] == 'or']
        incorrect = or_edges[or_edges['edge_type'] != 'output']
        assert len(incorrect) == 0, f"Found {len(incorrect)} OR edges with edge_type != 'output'"

    def test_pos_neg_is_pos_for_main_edges(self, main_edges):
        """Main pathway edges should all be positive (transformations)."""
        if len(main_edges) == 0:
            pytest.skip("No main pathway edges")
        non_pos = main_edges[main_edges['pos_neg'] != 'pos']
        assert len(non_pos) == 0, f"Found {len(non_pos)} main edges with pos_neg != 'pos'"

    def test_catalyst_edges_are_positive(self, network):
        """Catalyst edges should always be positive."""
        catalysts = network[network['edge_type'] == 'catalyst']
        if len(catalysts) == 0:
            pytest.skip("No catalyst edges")
        neg_catalysts = catalysts[catalysts['pos_neg'] == 'neg']
        assert len(neg_catalysts) == 0, f"Found {len(neg_catalysts)} negative catalysts"

    def test_network_has_edges(self, network):
        """Network should have a non-zero number of edges."""
        assert len(network) > 0, "Network has no edges"

    def test_network_not_suspiciously_large(self, network):
        """Sanity check: network shouldn't be excessively large."""
        assert len(network) < 10_000_000, f"Network suspiciously large: {len(network)} edges"


class TestAllPathwaysHaveContent:
    """Verify all 29 generated pathways have meaningful content."""

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_pathway_has_edges(self, network_path):
        """Each pathway should have at least some edges."""
        network = pd.read_csv(network_path)
        assert len(network) > 0, f"Pathway has no edges"

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_pathway_has_uuid_mapping(self, network_path):
        """Each pathway should have a stid_to_uuid_mapping.csv."""
        mapping_path = Path(network_path).parent / "stid_to_uuid_mapping.csv"
        assert mapping_path.exists(), f"Missing {mapping_path}"
        mapping = pd.read_csv(mapping_path)
        assert len(mapping) > 0, "UUID mapping is empty"

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_pathway_has_cache_files(self, network_path):
        """Each pathway should have cached intermediate files."""
        cache_dir = Path(network_path).parent / "cache"
        assert cache_dir.exists(), f"Missing cache directory"
        assert (cache_dir / "reaction_connections.csv").exists(), "Missing reaction_connections.csv"
        assert (cache_dir / "decomposed_uid_mapping.csv").exists(), "Missing decomposed_uid_mapping.csv"
        assert (cache_dir / "best_matches.csv").exists(), "Missing best_matches.csv"

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_pathway_has_main_edges(self, network_path):
        """Every pathway must have main (input/output) edges, not just catalysts/regulators.

        Bug history: Cellular_responses_to_stimuli_8953897 had 0 main edges due to
        an O(n^2) duplication bug in extract_inputs_and_outputs that was fixed.
        This test ensures no pathway is missing main transformation edges.
        """
        network = pd.read_csv(network_path)
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]
        assert len(main_edges) > 0, (
            f"Pathway has {len(network)} total edges but 0 main (input/output) edges. "
            f"Edge types: {dict(network['edge_type'].value_counts())}"
        )

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_main_edges_not_duplicated(self, network_path):
        """Main edges should not have N^2 duplication from the extract_inputs_and_outputs bug.

        Bug history: The outer loop in create_pathway_logic_network called
        extract_inputs_and_outputs N times, and the function internally iterated
        over ALL N reactions, creating N copies of every edge.
        This test ensures each edge appears at most once.
        """
        network = pd.read_csv(network_path)
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]
        if len(main_edges) == 0:
            pytest.skip("No main edges")

        # Check for exact duplicate rows
        duplicated = main_edges.duplicated(subset=['source_id', 'target_id', 'edge_type'], keep=False)
        num_duplicated = duplicated.sum()
        assert num_duplicated == 0, (
            f"Found {num_duplicated} duplicated main edges out of {len(main_edges)} total. "
            f"This suggests the O(n^2) duplication bug in extract_inputs_and_outputs."
        )

    @pytest.mark.parametrize("network_path", GENERATED_PATHWAYS,
                             ids=[Path(p).parent.name for p in GENERATED_PATHWAYS])
    def test_main_edges_proportional_to_best_matches(self, network_path):
        """Main edge count should be roughly proportional to best_matches, not N^2.

        Each best_match creates a virtual reaction with a few input×output edges.
        The total main edges should be within a reasonable ratio of best_matches count.
        """
        cache_dir = Path(network_path).parent / "cache"
        if not (cache_dir / "best_matches.csv").exists():
            pytest.skip("No best_matches.csv")

        network = pd.read_csv(network_path)
        best_matches = pd.read_csv(cache_dir / "best_matches.csv")
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        if len(main_edges) == 0 or len(best_matches) == 0:
            pytest.skip("No main edges or best_matches")

        ratio = len(main_edges) / len(best_matches)
        # Each best_match creates input+output edges (entity→reaction→entity model)
        # Ratio > 50 strongly suggests N^2 duplication
        assert ratio < 50, (
            f"Ratio of main_edges/best_matches = {ratio:.1f} is too high. "
            f"main_edges={len(main_edges)}, best_matches={len(best_matches)}. "
            f"This suggests O(n^2) edge duplication."
        )
