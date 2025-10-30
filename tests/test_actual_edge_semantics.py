"""Test to understand what edges actually represent by examining real data."""

import pandas as pd


class TestActualEdgeSemantics:
    """Examine real pathway data to understand edge semantics."""

    def test_examine_real_non_self_loop_edges(self):
        """
        Load the real pathway data and examine non-self-loop edges
        to understand what they actually represent.
        """
        # Load the real data
        network = pd.read_csv('pathway_logic_network_69620.csv')
        main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

        # Find non-self-loop edges
        non_self_loops = main_edges[main_edges['source_id'] != main_edges['target_id']]

        print("\n=== Real Pathway Data Analysis ===")
        print(f"Total main pathway edges: {len(main_edges)}")
        print(f"Self-loop edges: {len(main_edges) - len(non_self_loops)}")
        print(f"Non-self-loop edges: {len(non_self_loops)}")

        if len(non_self_loops) > 0:
            print("\nSample non-self-loop edges:")
            for idx, edge in non_self_loops.head(5).iterrows():
                print(f"  {edge['source_id']} → {edge['target_id']}")
                print(f"    AND/OR: {edge['and_or']}, Edge Type: {edge['edge_type']}")

            # Get the unique physical entities involved
            all_sources = set(non_self_loops['source_id'].unique())
            all_targets = set(non_self_loops['target_id'].unique())
            all_entities = all_sources | all_targets

            print(f"\nUnique physical entities in non-self-loop edges: {len(all_entities)}")

            # Check if these entities also appear in self-loop edges
            self_loop_entities = set(main_edges[main_edges['source_id'] == main_edges['target_id']]['source_id'].unique())
            overlap = all_entities & self_loop_entities

            print(f"Physical entities that appear in BOTH self-loops and non-self-loops: {len(overlap)}")

            # This tells us if the same entities can have both types of edges
            if len(overlap) > 0:
                print("\nThis suggests physical entities can have edges to themselves AND to other entities")
                print("Which means edges might represent different types of relationships")
            else:
                print("\nPhysical entities either have self-loop edges OR non-self-loop edges, not both")
                print("This suggests different categories of physical entities")

        # NOW the key question: what do these different entities represent?
        # Are they from different reactions? Different stages of decomposition?

        # Let's also check: do source and target entities cluster?
        sources_only = set(non_self_loops['source_id'].unique()) - set(non_self_loops['target_id'].unique())
        targets_only = set(non_self_loops['target_id'].unique()) - set(non_self_loops['source_id'].unique())
        both = set(non_self_loops['source_id'].unique()) & set(non_self_loops['target_id'].unique())

        print("\n=== Node Role Analysis ===")
        print(f"Physical entities that are ONLY sources: {len(sources_only)}")
        print(f"Physical entities that are ONLY targets: {len(targets_only)}")
        print(f"Physical entities that are BOTH: {len(both)}")

        # If we have clear sources and targets, that suggests directed flow
        # If most are "both", that suggests a more interconnected structure

    def test_hypothesis_multiple_reactions_same_entity(self):
        """
        Hypothesis: Non-self-loop edges occur when multiple reactions
        produce or consume variations of the same physical entity.

        For example:
          - R1 outputs Complex(A,B)
          - R2 outputs Complex(A,C)
          - R3 inputs Complex(A,B) and Complex(A,C)

        After decomposition, both complexes might share component A,
        leading to edges between different complex representations.
        """
        print("\n=== Hypothesis Testing ===")
        print("This hypothesis requires examining the decomposed_uid_mapping")
        print("to see if different complexes share components.")
        print("\nFor now, this is a placeholder for future investigation.")

        # TODO: Load decomposed_uid_mapping and check if physical entities
        # that have non-self-loop edges represent decomposed components
        # from different parent entities
