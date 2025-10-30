"""Example: Generate and analyze a pathway logic network.

This script demonstrates how to:
1. Generate a logic network for a specific Reactome pathway
2. Analyze network properties (root inputs, terminal outputs, edge counts)
3. Export the network for further analysis

Prerequisites:
- Neo4j database with Reactome data running at localhost:7687
- Poetry environment with dependencies installed

Usage:
    poetry run python examples/generate_pathway_example.py
"""

import sys
sys.path.insert(0, '.')

import pandas as pd
from src.pathway_generator import generate_pathway_file
from src.logic_network_generator import find_root_inputs, find_terminal_outputs


def main():
    """Generate and analyze a pathway logic network."""

    # Example pathway: Cell Cycle (Reactome ID: 69620)
    # This is a well-studied pathway with moderate complexity
    pathway_id = "69620"
    pathway_name = "Cell Cycle, Mitotic"
    taxon_id = "9606"  # Homo sapiens

    print("="*70)
    print("Logic Network Generator - Example Usage")
    print("="*70)
    print(f"\nGenerating logic network for pathway: {pathway_name}")
    print(f"Pathway ID: {pathway_id}")
    print(f"Taxon ID: {taxon_id}\n")

    try:
        # Generate the pathway logic network
        # This will create several CSV files:
        # - reaction_connections_{pathway_id}.csv
        # - decomposed_uid_mapping_{pathway_id}.csv
        # - best_matches_{pathway_id}.csv
        # - pathway_logic_network_{pathway_id}.csv (the final output)
        print("Step 1: Fetching reactions from Neo4j...")
        print("Step 2: Decomposing complexes and entity sets...")
        print("Step 3: Matching inputs and outputs...")
        print("Step 4: Creating logic network...\n")

        generate_pathway_file(
            pathway_id=pathway_id,
            taxon_id=taxon_id,
            pathway_name=pathway_name,
            decompose=False
        )

        print("\n" + "="*70)
        print("Generation Complete!")
        print("="*70)

        # Load the generated network for analysis
        network_file = f"pathway_logic_network_{pathway_id}.csv"
        network = pd.read_csv(network_file)

        # Analyze network properties
        print(f"\nNetwork Analysis:")
        print(f"  Total edges: {len(network)}")

        # Count edge types
        edge_types = network['edge_type'].value_counts()
        print(f"\n  Edge types:")
        for edge_type, count in edge_types.items():
            print(f"    - {edge_type}: {count}")

        # Count AND/OR relationships
        print(f"\n  Logic relationships:")
        and_edges = len(network[network['and_or'] == 'and'])
        or_edges = len(network[network['and_or'] == 'or'])
        print(f"    - AND edges (required): {and_edges}")
        print(f"    - OR edges (alternatives): {or_edges}")

        # Find root inputs and terminal outputs
        root_inputs = find_root_inputs(network)
        terminal_outputs = find_terminal_outputs(network)
        print(f"\n  Network structure:")
        print(f"    - Root inputs (starting points): {len(root_inputs)}")
        print(f"    - Terminal outputs (endpoints): {len(terminal_outputs)}")

        # Unique physical entities
        unique_sources = network['source_id'].nunique()
        unique_targets = network['target_id'].nunique()
        all_entities = set(network['source_id'].unique()) | set(network['target_id'].unique())
        print(f"    - Unique physical entities: {len(all_entities)}")

        # Sample edges
        print(f"\n  Sample edges (first 5):")
        sample_edges = network.head(5)
        for idx, edge in sample_edges.iterrows():
            print(f"    {edge['source_id'][:8]}... → {edge['target_id'][:8]}... "
                  f"({edge['and_or'].upper()}, {edge['edge_type']})")

        print("\n" + "="*70)
        print("Output Files:")
        print("="*70)
        print(f"  Main output: {network_file}")
        print(f"  Cached files:")
        print(f"    - reaction_connections_{pathway_id}.csv")
        print(f"    - decomposed_uid_mapping_{pathway_id}.csv")
        print(f"    - best_matches_{pathway_id}.csv")

        print("\n" + "="*70)
        print("Next Steps:")
        print("="*70)
        print("  1. Load the network in your analysis tool (Cytoscape, NetworkX, etc.)")
        print("  2. Run perturbation experiments by removing root inputs")
        print("  3. Analyze pathway flow from roots to terminals")
        print("  4. Identify key intermediate nodes")
        print("\nFor more pathways, see: https://reactome.org/PathwayBrowser/\n")

    except ConnectionError as e:
        print(f"\n❌ Connection Error: {e}")
        print("\nTroubleshooting:")
        print("  1. Ensure Neo4j is running: docker ps")
        print("  2. Start Neo4j if needed:")
        print("     docker run -p 7474:7474 -p 7687:7687 \\")
        print("       -e NEO4J_dbms_memory_heap_maxSize=8g \\")
        print("       public.ecr.aws/reactome/graphdb:Release94")
        sys.exit(1)

    except ValueError as e:
        print(f"\n❌ Validation Error: {e}")
        print("\nTroubleshooting:")
        print("  1. Verify the pathway ID is correct")
        print("  2. Check that the pathway exists in Reactome database")
        print("  3. Try a different pathway ID (e.g., 69620, 68875)")
        sys.exit(1)

    except Exception as e:
        print(f"\n❌ Unexpected Error: {e}")
        print("\nPlease report this issue at:")
        print("  https://github.com/reactome/logic-network-generator/issues")
        sys.exit(1)


if __name__ == "__main__":
    main()
