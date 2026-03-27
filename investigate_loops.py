#!/usr/bin/env python3
"""
Investigate the specific loops found in Reactome vs generated network.
"""

import pandas as pd
from pathlib import Path
from py2neo import Graph
import networkx as nx


def get_entity_name(graph: Graph, entity_id: int) -> str:
    """Get display name for an entity."""
    query = f'''
    MATCH (e {{dbId: {entity_id}}})
    RETURN e.displayName AS name, labels(e) AS labels
    '''
    result = graph.run(query).data()
    if result:
        return f"{result[0]['name']} ({result[0]['labels'][0]})"
    return str(entity_id)


def analyze_reactome_loops(graph: Graph, pathway_id: int):
    """Analyze the 5 loops found in Reactome."""
    print("=" * 80)
    print("REACTOME LOOPS - DETAILED ANALYSIS")
    print("=" * 80)

    # Build entity network
    query = f'''
    MATCH (p:Pathway {{dbId: {pathway_id}}})-[:hasEvent*]->(r:ReactionLikeEvent)
    MATCH (r)-[:input]->(inp)
    MATCH (r)-[:output]->(out)
    WHERE inp.dbId IS NOT NULL AND out.dbId IS NOT NULL
    RETURN DISTINCT inp.dbId AS input_entity, out.dbId AS output_entity,
           r.dbId AS reaction_id, r.displayName AS reaction_name
    '''

    edges = graph.run(query).data()

    # Build graph
    G = nx.DiGraph()
    edge_details = {}
    for edge in edges:
        inp = edge['input_entity']
        out = edge['output_entity']
        G.add_edge(inp, out)
        if (inp, out) not in edge_details:
            edge_details[(inp, out)] = []
        edge_details[(inp, out)].append({
            'reaction_id': edge['reaction_id'],
            'reaction_name': edge['reaction_name']
        })

    cycles = list(nx.simple_cycles(G))
    print(f"\nFound {len(cycles)} loops:")

    for i, cycle in enumerate(cycles, 1):
        print(f"\n{'='*80}")
        print(f"Loop {i}: Length {len(cycle)}")
        print('='*80)

        # Print cycle with entity names
        for j, entity_id in enumerate(cycle):
            entity_name = get_entity_name(graph, entity_id)
            next_entity_id = cycle[(j + 1) % len(cycle)]
            next_entity_name = get_entity_name(graph, next_entity_id)

            print(f"\n{entity_id}: {entity_name}")
            print(f"  ↓")

            # Show reactions connecting these entities
            if (entity_id, next_entity_id) in edge_details:
                for reaction in edge_details[(entity_id, next_entity_id)]:
                    print(f"  via Reaction {reaction['reaction_id']}: {reaction['reaction_name']}")

        print(f"\n  ↓ (back to {cycle[0]})")

        # Check if entities in this loop appear in decomposed network
        print(f"\n🔍 Checking if loop entities appear in generated network...")
        check_entities_in_generated_network(cycle, pathway_id)


def check_entities_in_generated_network(entity_ids: list, pathway_id: int):
    """Check if entities from a Reactome loop appear in the generated network."""
    decomposed = pd.read_csv(f'output/decomposed_uid_mapping_{pathway_id}.csv')

    for entity_id in entity_ids:
        # Check if this entity appears in decomposition
        matches = decomposed[decomposed['component_id_or_reference_entity_id'] == entity_id]

        if len(matches) > 0:
            uuids = matches['uid'].unique()
            print(f"  - Entity {entity_id}: Found in {len(matches)} decomposed rows, {len(uuids)} unique UUIDs")
        else:
            print(f"  - Entity {entity_id}: NOT FOUND in decomposed network")


def analyze_generated_loop(pathway_id: int):
    """Analyze the 1 loop found in generated network."""
    print("\n" + "=" * 80)
    print("GENERATED NETWORK LOOP - DETAILED ANALYSIS")
    print("=" * 80)

    network = pd.read_csv(f'output/pathway_logic_network_{pathway_id}.csv')
    main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]

    # Build graph
    G = nx.DiGraph()
    for _, edge in main_edges.iterrows():
        G.add_edge(edge['source_id'], edge['target_id'])

    cycles = list(nx.simple_cycles(G))

    if cycles:
        cycle = cycles[0]
        print(f"\nLoop of length {len(cycle)}:")

        # Load UUID mapping to get entity info
        uuid_mapping = pd.read_csv(f'output/uuid_mapping_{pathway_id}.csv')
        decomposed = pd.read_csv(f'output/decomposed_uid_mapping_{pathway_id}.csv')

        for i, uuid in enumerate(cycle):
            next_uuid = cycle[(i + 1) % len(cycle)]

            # Get entity info
            uuid_info = uuid_mapping[uuid_mapping['uuid'] == uuid]
            if len(uuid_info) > 0:
                entity_name = uuid_info.iloc[0]['entity_name']
                position = uuid_info.iloc[0]['position']
                print(f"\nUUID: {uuid[:16]}...")
                print(f"  Entity: {entity_name}")
                print(f"  Position: {position}")
            else:
                print(f"\nUUID: {uuid[:16]}... (no name found)")

            # Get component details
            components = decomposed[decomposed['uid'] == uuid]
            if len(components) > 0:
                comp_ids = components['component_id_or_reference_entity_id'].unique()
                print(f"  Components: {list(comp_ids)}")

            print(f"  ↓ connects to {next_uuid[:16]}...")


def main():
    pathway_id = 69620
    graph = Graph('bolt://localhost:7687', auth=('neo4j', 'test'))

    analyze_reactome_loops(graph, pathway_id)
    analyze_generated_loop(pathway_id)

    print("\n" + "=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print("\nReactome has 5 loops, generated network has 1.")
    print("This difference may occur because:")
    print("  1. Decomposition breaks complexes into components")
    print("  2. Some loops at the complex level don't exist at component level")
    print("  3. Position-aware UUIDs distinguish same entity at different positions")
    print("=" * 80)


if __name__ == "__main__":
    main()
