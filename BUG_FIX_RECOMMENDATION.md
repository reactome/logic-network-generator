# Bug Fix Recommendation: create_uid_reaction_connections

## Problem Statement

**Current Behavior**: Pathway 69620 generates ZERO main pathway edges (only 37 catalysts + 8 regulators)

**Expected Behavior**: Should generate input→output transformation edges representing the biochemical reactions

## Root Cause Analysis

### The Fundamental Misunderstanding

The current code confuses two different concepts:

1. **Input/Output pairing WITHIN reactions** (`best_matches`)
   - Pairs decomposed inputs with decomposed outputs for the SAME reaction
   - Example: Reaction 141429 has input_hash `ae0ebb...` → output_hash `33a1d5...`
   - Both hashes have `reactome_id = 141429`

2. **Pathway connectivity BETWEEN reactions** (what `create_uid_reaction_connections` should do)
   - Should connect reactions based on shared physical entities
   - Example: If Reaction A outputs Entity X, and Reaction B inputs Entity X, then A→B

### The Bug (lines 109-144 in src/logic_network_generator.py)

```python
def create_uid_reaction_connections(
    reaction_id_map: pd.DataFrame,
    best_matches: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame
) -> pd.DataFrame:
    # BUG: This loses 27% of virtual reactions (74 → 54)
    reactome_id_to_uid_mapping = dict(
        zip(reaction_id_map["reactome_id"], reaction_id_map["uid"])
    )

    uid_reaction_connections_data = []

    for _, match in best_matches.iterrows():
        incomming_hash = match["incomming"]
        outgoing_hash = match["outgoing"]

        # BUG: These are ALWAYS equal (both from same reaction!)
        preceding_reaction_id = _get_reactome_id_from_hash(decomposed_uid_mapping, incomming_hash)
        following_reaction_id = _get_reactome_id_from_hash(decomposed_uid_mapping, outgoing_hash)

        # BUG: Maps same reactome_id to same UID → self-loop!
        preceding_uid = reactome_id_to_uid_mapping.get(preceding_reaction_id)
        following_uid = reactome_id_to_uid_mapping.get(following_reaction_id)

        # Creates self-loop 87% of the time
        uid_reaction_connections_data.append({
            "preceding_uid": preceding_uid,
            "following_uid": following_uid
        })
```

**Empirical Evidence**:
- 62 connections created
- 54 are self-loops (87%)
- Only 8 valid connections (13%)
- Result: extract_inputs_and_outputs() finds almost no preceding reactions → no edges created

## Recommended Fix

### Option A: Use Original reaction_connections (RECOMMENDED)

The correct pathway topology already exists in `reaction_connections` (from Neo4j `precedingEvent` relationships). Just map it to virtual reactions:

```python
def create_uid_reaction_connections_FIXED(
    reaction_id_map: pd.DataFrame,
    reaction_connections: pd.DataFrame,  # Add this parameter!
    decomposed_uid_mapping: pd.DataFrame,
    best_matches: pd.DataFrame
) -> pd.DataFrame:
    """Create connections between virtual reactions based on pathway topology."""

    uid_reaction_connections_data = []

    # Iterate over ORIGINAL pathway connections
    for _, conn in reaction_connections.iterrows():
        preceding_reactome_id = conn["preceding_reaction_id"]
        following_reactome_id = conn["following_reaction_id"]

        # Skip rows with no preceding event
        if pd.isna(preceding_reactome_id) or pd.isna(following_reactome_id):
            continue

        # Get all virtual reactions for these reactome_ids
        preceding_virtual_reactions = reaction_id_map[
            reaction_id_map["reactome_id"] == preceding_reactome_id
        ]
        following_virtual_reactions = reaction_id_map[
            reaction_id_map["reactome_id"] == following_reactome_id
        ]

        # Connect virtual reactions based on shared physical entities
        for _, prec_vr in preceding_virtual_reactions.iterrows():
            prec_output_hash = prec_vr["output_hash"]
            prec_output_entities = decomposed_uid_mapping[
                decomposed_uid_mapping["uid"] == prec_output_hash
            ]["component_id_or_reference_entity_id"].tolist()

            for _, foll_vr in following_virtual_reactions.iterrows():
                foll_input_hash = foll_vr["input_hash"]
                foll_input_entities = decomposed_uid_mapping[
                    decomposed_uid_mapping["uid"] == foll_input_hash
                ]["component_id_or_reference_entity_id"].tolist()

                # Check for shared entities
                shared = set(prec_output_entities) & set(foll_input_entities)

                if len(shared) > 0:
                    # Create connection
                    uid_reaction_connections_data.append({
                        "preceding_uid": prec_vr["uid"],
                        "following_uid": foll_vr["uid"],
                        "shared_entities": len(shared)
                    })

    return pd.DataFrame(uid_reaction_connections_data)
```

### Option B: Infer from Shared Physical Entities

If `reaction_connections` isn't available, infer connectivity from shared physical entities:

```python
def create_uid_reaction_connections_from_entities(
    reaction_id_map: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame
) -> pd.DataFrame:
    """Infer virtual reaction connections from shared physical entities."""

    uid_reaction_connections_data = []

    # For each virtual reaction
    for idx1, vr1 in reaction_id_map.iterrows():
        vr1_output_hash = vr1["output_hash"]
        vr1_outputs = decomposed_uid_mapping[
            decomposed_uid_mapping["uid"] == vr1_output_hash
        ]["component_id_or_reference_entity_id"].tolist()

        # Find virtual reactions whose inputs match vr1's outputs
        for idx2, vr2 in reaction_id_map.iterrows():
            if idx1 == idx2:
                continue  # Skip self

            vr2_input_hash = vr2["input_hash"]
            vr2_inputs = decomposed_uid_mapping[
                decomposed_uid_mapping["uid"] == vr2_input_hash
            ]["component_id_or_reference_entity_id"].tolist()

            # Check for shared entities
            shared = set(vr1_outputs) & set(vr2_inputs)

            if len(shared) > 0:
                uid_reaction_connections_data.append({
                    "preceding_uid": vr1["uid"],
                    "following_uid": vr2["uid"],
                    "shared_entities": len(shared)
                })

    return pd.DataFrame(uid_reaction_connections_data)
```

## Implementation Steps

1. **Backup current code**
   ```bash
   cp src/logic_network_generator.py src/logic_network_generator.py.backup
   ```

2. **Implement Option A** (recommended - uses existing Reactome topology)
   - Modify `create_uid_reaction_connections` signature to accept `reaction_connections`
   - Implement the fixed logic
   - Update call site in `create_pathway_logic_network` (line 674)

3. **Add test for correctness**
   ```python
   def test_uid_reaction_connections_no_self_loops():
       """Verify uid_reaction_connections doesn't create excessive self-loops."""
       # Generate pathway 69620
       # Load uid_reaction_connections
       # Assert: self-loops < 10% of connections
       # Assert: len(uid_reaction_connections) > 50
   ```

4. **Regenerate pathway 69620**
   ```bash
   rm output/pathway_logic_network_69620.csv
   poetry run python bin/create-pathways.py --pathway-id 69620
   ```

5. **Verify results**
   - Check that main pathway edges exist
   - Verify edge count is reasonable (should be >> 45)
   - Run full test suite

## Expected Outcomes After Fix

### Before Fix:
- **Total edges**: 45
  - Main pathway edges: 0 ❌
  - Catalyst edges: 37
  - Regulator edges: 8
- **uid_reaction_connections**: 87% self-loops

### After Fix (Expected):
- **Total edges**: 500-2000 (estimated)
  - Main pathway edges: 400-1900 ✅
  - Catalyst edges: 37
  - Regulator edges: 8
- **uid_reaction_connections**: < 10% self-loops

## Testing Strategy

1. **Unit test for the fix**
   - Mock data with 2-3 reactions
   - Verify correct connections created
   - Verify no self-loops

2. **Integration test with pathway 69620**
   - Regenerate network
   - Verify main edges exist
   - Compare against manual Reactome query

3. **Regression test with multiple pathways**
   - Test 5-10 different pathways
   - Ensure all generate reasonable edge counts
   - Verify no pathway has 0 main edges

## Alternative: Is This By Design?

**Question**: Could pathway 69620 be a special case where no main edges is correct?

**Answer**: NO. Evidence:
1. Reactome shows reaction 141429 has inputs (141412, 141447) and output (141408)
2. These entities should create transformation edges
3. The 87% self-loop rate is clearly a bug, not a feature
4. Catalysts/regulators working suggests Neo4j queries are fine, so the issue is specific to main edge logic

## Priority

**CRITICAL** - This prevents the system from generating the core functionality (transformation edges). All generated networks are missing their primary content.

---

## Additional Notes

- The cartesian product edge creation in `extract_inputs_and_outputs` is fine
- The Hungarian algorithm best matching is working correctly
- The decomposition algorithm is sound
- Only this specific function needs fixing

**Estimated Effort**: 4-8 hours (implementation + testing)
