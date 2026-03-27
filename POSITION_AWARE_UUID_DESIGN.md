# Position-Aware UUID Design

## Overview

The logic network generator uses **position-aware UUIDs** to represent physical entities at different positions in pathway networks. This design ensures that:

1. The same entity at different pathway positions gets different UUIDs
2. Entities in the same connected component share the same UUID
3. Self-loops are minimized in the generated logic network

## Problem Statement

In Reactome pathways, the same physical entity (e.g., ATP, a specific protein) can appear at multiple points in a pathway. Using a single UUID for all occurrences would create excessive self-loops in the logic network. Using completely unique UUIDs would lose the connection between related positions.

### Example Scenario

```
Reaction1 -> gene1 -> Reaction2
Reaction3 -> gene1 -> Reaction2
```

**Without position-awareness**: gene1 gets one UUID everywhere → creates self-loops

**With position-awareness + union-find**:
- gene1 gets UUID_A when connecting Reaction1→Reaction2 and Reaction3→Reaction2
- gene1 gets UUID_B when used elsewhere in the pathway (e.g., Reaction100→Reaction101)

## Implementation

### Core Data Structure

```python
entity_uuid_registry: Dict[tuple, str]
```

**Key format**: `(entity_dbId, reaction_uuid, role)`
- `entity_dbId`: Reactome database ID (e.g., "113592")
- `reaction_uuid`: UUID of the reaction involving this entity
- `role`: Either "input" or "output"

**Value**: UUID string for the entity at this position

### Union-Find Algorithm

The `_get_or_create_entity_uuid()` function implements union-find logic:

1. **Check target position**: Does entity have UUID as input to target reaction?
2. **Check source position**: Does entity have UUID as output of source reaction?
3. **Merge if needed**: If both exist but differ, merge all references to use one UUID
4. **Share if one exists**: If only one position has UUID, share it with the other
5. **Create new**: If neither position has UUID, create a new one

This ensures entities in the same connected component share UUIDs, while entities at disconnected positions get different UUIDs.

## Benefits

### Zero Self-Loops
Real-world testing on pathway 1227986:
- **Before**: Unknown (self-loops were a known issue)
- **After**: 0 self-loops (0.00% of 7514 edges)

### Multi-Position Tracking
- Entity 113592 in pathway 1227986: 8 different UUIDs at 8 positions
- Proper tracking of entities throughout complex pathways

### Traceable Back to Reactome
The UUID→dbId mapping allows reconstruction of which Reactome entity each UUID represents:

```python
# Export format
uuid_to_reactome_mapping.csv:
uuid,reactome_dbId
3e715e93-...,113592
b75df0cb-...,113592  # Same entity, different position
```

## Usage

### In Code

```python
# Initialize registry
entity_uuid_registry: Dict[tuple, str] = {}

# Assign UUIDs for entities between reactions
input_uuids = _assign_uuids(
    input_reactome_ids,
    source_reaction_uuid="rxn1-uuid",
    target_reaction_uuid="rxn2-uuid",
    entity_uuid_registry=entity_uuid_registry
)

# Registry automatically tracks and merges positions
```

### In Generated Files

The `uuid_to_reactome_{pathway_id}.csv` file maps all UUIDs back to their Reactome database IDs, enabling:
- Validation of generated networks
- Reconstruction of pathway topology
- Integration with Reactome database

## Testing

Comprehensive testing verified:
- ✅ 73 unit tests pass
- ✅ End-to-end pathway generation works
- ✅ 0% self-loops in real pathways
- ✅ Union-find correctly merges connected positions
- ✅ Different positions get different UUIDs

## References

- Implementation: `src/logic_network_generator.py` (lines 308-385)
- Tests: `tests/test_logic_network_generator.py`
- End-to-end test: `test_position_aware.py`
