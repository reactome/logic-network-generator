#!/usr/bin/env python3
"""Quick test of position-aware UUID implementation."""

import pandas as pd
from src.logic_network_generator import create_pathway_logic_network
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types

# Use pathway 1227986 which has cached files
pathway_id = "1227986"

print(f"Testing position-aware UUIDs with pathway {pathway_id}")
print("=" * 80)

# Load cached data
print("\n1. Loading cached data...")
reaction_connections = pd.read_csv(f"output/reaction_connections_{pathway_id}.csv")
decomposed_uid_mapping = pd.read_csv(
    f"output/decomposed_uid_mapping_{pathway_id}.csv",
    dtype=decomposed_uid_mapping_column_types
)
best_matches = pd.read_csv(f"output/best_matches_{pathway_id}.csv")

print(f"   - Reaction connections: {len(reaction_connections)} rows")
print(f"   - Decomposed UID mapping: {len(decomposed_uid_mapping)} rows")
print(f"   - Best matches: {len(best_matches)} rows")

# Generate logic network
print("\n2. Generating logic network...")
try:
    result = create_pathway_logic_network(
        decomposed_uid_mapping, reaction_connections, best_matches
    )
    print(f"   ✓ Success! Generated {len(result.logic_network)} edges")
except Exception as e:
    print(f"   ✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Analyze UUID mapping
print("\n3. Analyzing UUID mapping...")
print(f"   - Total unique UUIDs: {len(result.uuid_mapping)}")

# Count how many entities appear at multiple positions
from collections import Counter
entity_positions = Counter(result.uuid_mapping.values())
multi_position = {entity: count for entity, count in entity_positions.items() if count > 1}

print(f"   - Entities at single position: {len(entity_positions) - len(multi_position)}")
print(f"   - Entities at multiple positions: {len(multi_position)}")

if multi_position:
    max_positions = max(multi_position.values())
    example_entity = [e for e, c in multi_position.items() if c == max_positions][0]
    print(f"   - Max positions for one entity: {max_positions} (dbId: {example_entity})")

# Check for position-aware behavior
print("\n4. Checking position-aware behavior...")
# Find an entity that appears multiple times
if len(multi_position) > 0:
    # Look for this entity in the logic network
    example_entity_uuids = [uuid for uuid, dbId in result.uuid_mapping.items() if dbId == example_entity]
    print(f"   - Entity {example_entity} has {len(example_entity_uuids)} UUIDs:")
    for i, uuid in enumerate(example_entity_uuids[:3]):  # Show first 3
        # Find where this UUID appears in logic network
        as_source = result.logic_network[result.logic_network['source_id'] == uuid]
        as_target = result.logic_network[result.logic_network['target_id'] == uuid]
        print(f"     UUID {i+1} ({uuid[:8]}...): {len(as_source)} as source, {len(as_target)} as target")

    if len(example_entity_uuids) > 1:
        print(f"   ✓ Position-aware: same entity has different UUIDs at different positions!")
    else:
        print(f"   ✗ Warning: expected multiple UUIDs but found only one")
else:
    print("   - No multi-position entities found (pathway might be too simple)")

print("\n5. Checking for self-loops...")
self_loops = result.logic_network[result.logic_network['source_id'] == result.logic_network['target_id']]
self_loop_ratio = len(self_loops) / len(result.logic_network) if len(result.logic_network) > 0 else 0
print(f"   - Self-loops: {len(self_loops)} / {len(result.logic_network)} ({self_loop_ratio*100:.2f}%)")

if self_loop_ratio < 0.05:
    print(f"   ✓ Self-loop ratio is low (< 5%)")
else:
    print(f"   ✗ Warning: high self-loop ratio")

print("\n" + "=" * 80)
print("Test complete!")
print("=" * 80)
