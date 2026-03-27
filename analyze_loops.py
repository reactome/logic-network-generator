#!/usr/bin/env python3
"""
Analyze biological loops (cycles) in Reactome database vs generated logic network.

A biological loop occurs when a molecule/reaction participates in a pathway
that eventually produces itself.
"""

import pandas as pd
from pathlib import Path
from py2neo import Graph
from typing import Set, List, Dict
import networkx as nx


def find_loops_in_reactome(graph: Graph, pathway_id: int) -> List[List[int]]:
    """Find loops in Reactome database for a pathway.

    A loop exists when reaction R1 has an output that is eventually an input to R1
    through a chain of reactions.
    """
    # Get all reactions in pathway
    query = f'''
    MATCH (p:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
    RETURN DISTINCT r.dbId AS reaction_id
    '''
    reactions = [row['reaction_id'] for row in graph.run(query).data()]

    # Build reaction connectivity graph
    print(f"Found {len(reactions)} reactions in pathway {pathway_id}")

    # Get all precedingEvent relationships
    query = f'''
    MATCH (p:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(r1:ReactionLikeEvent)
    MATCH (r1)-[:precedingEvent]->(r2:ReactionLikeEvent)
    RETURN DISTINCT r1.dbId AS from_reaction, r2.dbId AS to_reaction
    '''

    edges = graph.run(query).data()
    print(f"Found {len(edges)} precedingEvent edges in pathway")

    # Build directed graph
    G = nx.DiGraph()
    for edge in edges:
        G.add_edge(edge['from_reaction'], edge['to_reaction'])

    # Find all cycles
    try:
        cycles = list(nx.simple_cycles(G))
        return cycles
    except:
        return []


def find_loops_in_generated_network(network_path: Path) -> List[List[str]]:
    """Find loops in generated logic network.

    A loop exists when entity A has a path back to itself through the network.
    """
    network = pd.read_csv(network_path)

    # Only use main pathway edges (not catalyst/regulator)
    main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

    print(f"\nGenerated network has {len(main_edges)} main pathway edges")

    # Build directed graph
    G = nx.DiGraph()
    for _, edge in main_edges.iterrows():
        G.add_edge(edge['source_id'], edge['target_id'])

    print(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

    # Find all cycles
    try:
        cycles = list(nx.simple_cycles(G))
        return cycles
    except:
        return []


def analyze_entity_level_loops_in_reactome(graph: Graph, pathway_id: int) -> List[List[int]]:
    """Find loops at the entity level (not reaction level) in Reactome.

    A loop exists when entity E is consumed by a reaction that produces E
    (directly or through a chain).
    """
    # Build entity-level network from Reactome
    query = f'''
    MATCH (p:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
    MATCH (r)-[:input]->(inp)
    MATCH (r)-[:output]->(out)
    WHERE inp.dbId IS NOT NULL AND out.dbId IS NOT NULL
    RETURN DISTINCT inp.dbId AS input_entity, out.dbId AS output_entity
    '''

    edges = graph.run(query).data()

    # Build directed graph at entity level
    G = nx.DiGraph()
    for edge in edges:
        G.add_edge(edge['input_entity'], edge['output_entity'])

    print(f"\nReactome entity-level network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Find all cycles
    try:
        cycles = list(nx.simple_cycles(G))
        return cycles
    except:
        return []


def main():
    """Compare loops between Reactome and generated network."""

    print("=" * 80)
    print("LOOP ANALYSIS: Reactome Database vs Generated Logic Network")
    print("=" * 80)

    pathway_id = 69620
    output_dir = Path('output')
    network_path = output_dir / 'pathway_logic_network_69620.csv'

    # Connect to Reactome
    graph = Graph('bolt://localhost:7687', auth=('neo4j', 'test'))

    # 1. Reaction-level loops in Reactome
    print("\n" + "=" * 80)
    print("1. REACTION-LEVEL LOOPS IN REACTOME")
    print("=" * 80)
    reactome_reaction_loops = find_loops_in_reactome(graph, pathway_id)
    print(f"\n✓ Found {len(reactome_reaction_loops)} reaction-level loops in Reactome")

    if reactome_reaction_loops:
        print("\nReaction loops:")
        for i, cycle in enumerate(reactome_reaction_loops[:5], 1):
            print(f"  {i}. Cycle of length {len(cycle)}: {' → '.join(map(str, cycle))} → {cycle[0]}")
        if len(reactome_reaction_loops) > 5:
            print(f"  ... and {len(reactome_reaction_loops) - 5} more")

    # 2. Entity-level loops in Reactome
    print("\n" + "=" * 80)
    print("2. ENTITY-LEVEL LOOPS IN REACTOME")
    print("=" * 80)
    reactome_entity_loops = analyze_entity_level_loops_in_reactome(graph, pathway_id)
    print(f"\n✓ Found {len(reactome_entity_loops)} entity-level loops in Reactome")

    if reactome_entity_loops:
        print("\nEntity loops (top 10):")
        # Sort by cycle length for readability
        sorted_loops = sorted(reactome_entity_loops, key=len)
        for i, cycle in enumerate(sorted_loops[:10], 1):
            print(f"  {i}. Cycle of length {len(cycle)}: {' → '.join(map(str, cycle[:5]))}{'...' if len(cycle) > 5 else ''}")
        if len(reactome_entity_loops) > 10:
            print(f"  ... and {len(reactome_entity_loops) - 10} more")

    # 3. Entity-level loops in generated network
    print("\n" + "=" * 80)
    print("3. ENTITY-LEVEL LOOPS IN GENERATED LOGIC NETWORK")
    print("=" * 80)
    generated_loops = find_loops_in_generated_network(network_path)
    print(f"\n✓ Found {len(generated_loops)} entity-level loops in generated network")

    if generated_loops:
        print("\nGenerated network loops (top 10):")
        sorted_loops = sorted(generated_loops, key=len)
        for i, cycle in enumerate(sorted_loops[:10], 1):
            # Show first 80 chars of each UUID
            cycle_str = ' → '.join([str(node)[:8] + '...' for node in cycle[:3]])
            if len(cycle) > 3:
                cycle_str += '...'
            print(f"  {i}. Cycle of length {len(cycle)}: {cycle_str}")
        if len(generated_loops) > 10:
            print(f"  ... and {len(generated_loops) - 10} more")

    # 4. Summary comparison
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nReactome Database:")
    print(f"  - Reaction-level loops: {len(reactome_reaction_loops)}")
    print(f"  - Entity-level loops: {len(reactome_entity_loops)}")

    print(f"\nGenerated Logic Network:")
    print(f"  - Entity-level loops: {len(generated_loops)}")

    print("\n" + "=" * 80)

    # Analysis
    if len(reactome_entity_loops) == 0 and len(generated_loops) == 0:
        print("✅ PERFECT MATCH: Neither Reactome nor generated network have loops")
    elif len(reactome_entity_loops) > 0 and len(generated_loops) > 0:
        print(f"✅ BOTH HAVE LOOPS: Reactome has {len(reactome_entity_loops)}, Generated has {len(generated_loops)}")
        print("   This is expected for pathways with feedback mechanisms.")
    elif len(reactome_entity_loops) > 0 and len(generated_loops) == 0:
        print(f"⚠️  MISMATCH: Reactome has {len(reactome_entity_loops)} loops, but generated network has 0")
        print("   The generated network may be missing feedback loops.")
    else:
        print(f"⚠️  MISMATCH: Reactome has 0 loops, but generated network has {len(generated_loops)}")
        print("   The generated network may have spurious cycles.")

    print("=" * 80)


if __name__ == "__main__":
    main()
