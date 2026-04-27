# Position-Aware UUIDs

## Why

In a Reactome pathway, the same physical entity (a protein, a small molecule like ATP) often appears at many positions. Two naive choices both fail:

- **One UUID per entity, everywhere** → every reuse becomes a self-loop in the logic network.
- **A fresh UUID at every position** → the connection between adjacent reactions that share an entity is lost.

Position-aware UUIDs sit between the two: the same entity gets the *same* UUID across positions in one connected component, and *different* UUIDs across disconnected positions.

## How it works

The assignment is keyed on the entity's role in a specific edge:

```
key   = (entity_dbId, reaction_uuid, role)   # role ∈ {"input", "output"}
value = uuid_string
```

When assigning a UUID for an entity flowing from `source_reaction` to `target_reaction`, `_get_or_create_entity_uuid` does a small union-find:

1. Look up the entity as input to `target_reaction`.
2. Look up the entity as output of `source_reaction`.
3. If both exist and differ — merge: rewrite all keys pointing at one UUID to point at the other.
4. If only one exists — share it with the other position.
5. If neither exists — mint a new UUID for both.

The result: entities reachable through shared reactions collapse to one UUID; entities at independent positions stay distinct.

## What gets exported

`uuid_to_reactome_{pathway_id}.csv` is the inverse map — every UUID alongside its Reactome `dbId`. Multiple UUIDs sharing the same `dbId` is normal: it means the entity appears at multiple disconnected positions.

```
uuid,reactome_dbId
3e715e93-...,113592
b75df0cb-...,113592
```

This file is what makes the network round-trippable back to Reactome.

## Where to look

- Implementation: `src/logic_network_generator.py` (`_get_or_create_entity_uuid`, `_assign_uuids`)
- Tests: `tests/test_uuid_position_bug.py`, `tests/test_logic_network_generator.py`
