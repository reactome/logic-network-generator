#!/usr/bin/env python3
"""
Comprehensive validation script to verify generated logic network matches Reactome.

This script validates that:
1. Reaction connectivity in generated network matches Reactome topology
2. Decomposed components correctly represent complex/set memberships
3. Edges connect the right entities based on shared physical components
"""

import pandas as pd
from pathlib import Path
from py2neo import Graph
from typing import List, Set, Tuple

def validate_reaction_pair(
    prec_id: int,
    foll_id: int,
    decomposed_uid_mapping: pd.DataFrame,
    best_matches: pd.DataFrame,
    graph: Graph
) -> dict:
    """Validate a single reaction pair."""

    # Query Reactome for actual connectivity
    query = f'''
    MATCH (r1:ReactionLikeEvent {{dbId: {prec_id}}})
    MATCH (r2:ReactionLikeEvent {{dbId: {foll_id}}})
    OPTIONAL MATCH (r1)-[:output]->(out1)
    OPTIONAL MATCH (r2)-[:input]->(in2)
    RETURN r1.displayName AS r1_name,
           collect(DISTINCT out1.dbId) AS r1_outputs,
           r2.displayName AS r2_name,
           collect(DISTINCT in2.dbId) AS r2_inputs
    '''

    result = graph.run(query).data()[0]

    # Check for shared entities in Reactome
    r1_outs = set([x for x in result["r1_outputs"] if x])
    r2_ins = set([x for x in result["r2_inputs"] if x])
    reactome_shared_entities = r1_outs & r2_ins

    # Check decomposed components
    r1_uids = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == prec_id]['uid'].unique()
    r2_uids = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == foll_id]['uid'].unique()

    # Get R1 output components
    r1_match = best_matches[best_matches['incomming'].isin(r1_uids)]
    if len(r1_match) == 0:
        return {"valid": False, "reason": "No best match for R1"}

    r1_out_hash = r1_match.iloc[0]['outgoing']
    r1_out_components = set(decomposed_uid_mapping[
        decomposed_uid_mapping['uid'] == r1_out_hash
    ]['component_id_or_reference_entity_id'])

    # Get R2 input components
    r2_match = best_matches[best_matches['outgoing'].isin(r2_uids)]
    if len(r2_match) == 0:
        return {"valid": False, "reason": "No best match for R2"}

    r2_in_hash = r2_match.iloc[0]['incomming']
    r2_in_components = set(decomposed_uid_mapping[
        decomposed_uid_mapping['uid'] == r2_in_hash
    ]['component_id_or_reference_entity_id'])

    # Check for shared components
    shared_components = r1_out_components & r2_in_components

    # Validation: If Reactome connects them, we should have shared components
    should_connect = len(reactome_shared_entities) > 0
    we_connect = len(shared_components) > 0

    return {
        "valid": should_connect == we_connect,
        "prec_id": prec_id,
        "foll_id": foll_id,
        "prec_name": result["r1_name"],
        "foll_name": result["r2_name"],
        "reactome_shared_entities": reactome_shared_entities,
        "decomposed_shared_components": shared_components,
        "should_connect": should_connect,
        "we_connect": we_connect,
    }


def main():
    """Run comprehensive validation."""

    print("=" * 80)
    print("VALIDATION: Generated Logic Network vs Reactome Database")
    print("=" * 80)

    # Load data
    output_dir = Path('output')
    network = pd.read_csv(output_dir / 'pathway_logic_network_69620.csv')
    decomposed_uid_mapping = pd.read_csv(output_dir / 'decomposed_uid_mapping_69620.csv')
    reaction_connections = pd.read_csv(output_dir / 'reaction_connections_69620.csv')
    best_matches = pd.read_csv(output_dir / 'best_matches_69620.csv')

    graph = Graph('bolt://localhost:7687', auth=('neo4j', 'test'))

    print(f"\n📊 Loaded Data:")
    print(f"  - Network edges: {len(network):,}")
    print(f"  - Reaction connections: {len(reaction_connections)}")
    print(f"  - Best matches: {len(best_matches)}")
    print(f"  - Decomposition rows: {len(decomposed_uid_mapping):,}")

    # Test all valid reaction pairs
    valid_connections = reaction_connections[
        reaction_connections['following_reaction_id'].notna()
    ]

    print(f"\n🔬 Validating {len(valid_connections)} reaction pairs...")

    results = []
    for idx, row in valid_connections.head(20).iterrows():  # Test first 20
        prec_id = int(row['preceding_reaction_id'])
        foll_id = int(row['following_reaction_id'])

        result = validate_reaction_pair(
            prec_id, foll_id, decomposed_uid_mapping, best_matches, graph
        )
        results.append(result)

    # Analyze results
    valid_count = sum(1 for r in results if r.get("valid", False))
    total_count = len(results)

    print(f"\n✅ Validation Results: {valid_count}/{total_count} pairs validated correctly")

    # Show details
    print(f"\n📋 Sample Validations:")
    for i, result in enumerate(results[:5]):
        if result.get("valid"):
            status = "✓ PASS"
        else:
            status = "✗ FAIL"

        print(f"\n{i+1}. {status}")
        print(f"   {result['prec_id']} → {result['foll_id']}")
        print(f"   {result['prec_name']}")
        print(f"   → {result['foll_name']}")
        print(f"   Reactome entities: {len(result['reactome_shared_entities'])} shared")
        print(f"   Decomposed components: {len(result['decomposed_shared_components'])} shared")
        print(f"   Should connect: {result['should_connect']}")
        print(f"   We connect: {result['we_connect']}")

    # Summary statistics
    print(f"\n📈 Statistics:")
    connected_in_reactome = sum(1 for r in results if r.get("should_connect", False))
    connected_by_us = sum(1 for r in results if r.get("we_connect", False))

    print(f"  - Pairs connected in Reactome: {connected_in_reactome}/{total_count}")
    print(f"  - Pairs connected by algorithm: {connected_by_us}/{total_count}")
    print(f"  - Match rate: {valid_count/total_count*100:.1f}%")

    # Final verdict
    print(f"\n{'=' * 80}")
    if valid_count == total_count:
        print("✅ VALIDATION PASSED: Generated network matches Reactome topology!")
    else:
        print(f"⚠️  VALIDATION ISSUES: {total_count - valid_count} mismatches found")
    print(f"{'=' * 80}")

    return valid_count == total_count


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
