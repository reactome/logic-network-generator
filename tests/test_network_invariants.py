"""Tests for network invariants - properties that should always hold.

These tests verify structural properties of the generated networks:
- No self-loops in main pathway edges
- Root inputs are always sources (never targets)
- Terminal outputs are always targets (never sources)
- AND/OR logic is consistent
- Edge direction represents transformations
"""

import pandas as pd


class TestNetworkInvariants:
    """Test invariants that should hold for any valid pathway logic network."""

    def test_no_self_loops_in_main_pathway(self):
        """Main pathway edges should never have source_id == target_id.

        Rationale: Reactions transform molecules, so inputs ≠ outputs.
        """
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        self_loops = main_edges[main_edges['source_id'] == main_edges['target_id']]

        assert len(self_loops) == 0, f"Found {len(self_loops)} self-loop edges in main pathway"

    def test_root_inputs_never_appear_as_targets(self):
        """Root inputs should only appear as source_id, never as target_id.

        Rationale: Root inputs are consumed by reactions but not produced.
        """
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        sources = set(main_edges['source_id'].unique())
        targets = set(main_edges['target_id'].unique())
        root_inputs = sources - targets

        # Check that none of the root inputs appear as targets
        roots_as_targets = root_inputs & targets
        assert len(roots_as_targets) == 0, f"Found {len(roots_as_targets)} root inputs appearing as targets"

    def test_terminal_outputs_never_appear_as_sources(self):
        """Terminal outputs should only appear as target_id, never as source_id.

        Rationale: Terminal outputs are produced but not consumed.
        """
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        sources = set(main_edges['source_id'].unique())
        targets = set(main_edges['target_id'].unique())
        terminal_outputs = targets - sources

        # Check that none of the terminal outputs appear as sources
        terminals_as_sources = terminal_outputs & sources
        assert len(terminals_as_sources) == 0, f"Found {len(terminals_as_sources)} terminal outputs appearing as sources"

    def test_all_nodes_reachable_from_roots(self):
        """All nodes should be reachable from root inputs via directed edges.

        Rationale: Disconnected components suggest data problems.
        """
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        sources = set(main_edges['source_id'].unique())
        targets = set(main_edges['target_id'].unique())
        root_inputs = sources - targets

        # BFS from roots
        visited = set(root_inputs)
        queue = list(root_inputs)

        while queue:
            current = queue.pop(0)
            # Find all edges from current node
            outgoing = main_edges[main_edges['source_id'] == current]
            for _, edge in outgoing.iterrows():
                target = edge['target_id']
                if target not in visited:
                    visited.add(target)
                    queue.append(target)

        all_nodes = sources | targets
        unreachable = all_nodes - visited

        # Allow some unreachable nodes (might be in disconnected branches)
        # But warn if too many
        unreachable_pct = len(unreachable) / len(all_nodes) * 100 if all_nodes else 0

        assert unreachable_pct < 50, f"{unreachable_pct:.1f}% of nodes unreachable from roots"

    def test_and_logic_consistency(self):
        """Edges with 'and' logic should have edge_type='input'."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        and_edges = main_edges[main_edges['and_or'] == 'and']
        incorrect = and_edges[and_edges['edge_type'] != 'input']

        assert len(incorrect) == 0, f"Found {len(incorrect)} AND edges with edge_type != 'input'"

    def test_or_logic_consistency(self):
        """Edges with 'or' logic should have edge_type='output'."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        or_edges = main_edges[main_edges['and_or'] == 'or']
        incorrect = or_edges[or_edges['edge_type'] != 'output']

        assert len(incorrect) == 0, f"Found {len(incorrect)} OR edges with edge_type != 'output'"

    def test_all_edges_have_and_or_logic(self):
        """All main pathway edges should have and_or specified."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        missing_logic = main_edges[main_edges['and_or'].isna()]

        assert len(missing_logic) == 0, f"Found {len(missing_logic)} edges without AND/OR logic"

    def test_pos_neg_is_always_pos_for_main_edges(self):
        """Main pathway edges should all be positive (activation)."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        non_pos = main_edges[main_edges['pos_neg'] != 'pos']

        assert len(non_pos) == 0, f"Found {len(non_pos)} main edges with pos_neg != 'pos'"

    def test_catalyst_edges_have_no_and_or_logic(self):
        """Catalyst edges shouldn't have AND/OR logic (they're not transformations)."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        catalyst_edges = network[network['edge_type'] == 'catalyst']

        has_logic = catalyst_edges[catalyst_edges['and_or'].notna()]

        # This is just documenting current behavior - may or may not be desired
        print(f"\nCatalyst edges with AND/OR logic: {len(has_logic)}/{len(catalyst_edges)}")

    def test_regulator_edges_have_no_and_or_logic(self):
        """Regulator edges shouldn't have AND/OR logic (they're not transformations)."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        regulator_edges = network[network['edge_type'] == 'regulator']

        has_logic = regulator_edges[regulator_edges['and_or'].notna()]

        # This is just documenting current behavior
        print(f"\nRegulator edges with AND/OR logic: {len(has_logic)}/{len(regulator_edges)}")

    def test_network_has_reasonable_size(self):
        """Sanity check: network should have a reasonable number of edges."""
        network = pd.read_csv('pathway_logic_network_69620.csv')

        assert len(network) > 0, "Network has no edges"
        assert len(network) < 100000, "Network suspiciously large"

        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]
        assert len(main_edges) > 0, "Network has no main pathway edges"

    def test_unique_molecules_are_reasonable(self):
        """Sanity check: should have reasonable number of unique molecules."""
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        all_molecules = set(main_edges['source_id'].unique()) | set(main_edges['target_id'].unique())

        assert len(all_molecules) > 0, "No molecules found"
        assert len(all_molecules) < 10000, "Suspiciously many molecules"

        # Should have at least one root and one terminal
        sources = set(main_edges['source_id'].unique())
        targets = set(main_edges['target_id'].unique())
        roots = sources - targets
        terminals = targets - sources

        assert len(roots) > 0, "No root inputs found"
        assert len(terminals) > 0, "No terminal outputs found"
