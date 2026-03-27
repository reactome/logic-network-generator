# EntitySet Parent Tracking Fix

## Problem

When we decompose EntitySets into their members, we lose track of which EntitySet they came from. This makes it impossible to accurately reconstruct the original pathway.

### Example

**Reaction 69598:** Ubiquitination of phosphorylated CDC25A
- **Neo4j Input:** EntitySet `9943734` (p-S82-CDC25A)
- **Generated:** Members `[9943706, 9943732]` (the alternatives)

**Current state:** We have the members but don't know they came from EntitySet `9943734`
**Needed:** Track that `9943706` and `9943732` both came from parent EntitySet `9943734`

## Current Data Structure

`decomposed_uid_mapping` has columns:
```
- uid: The virtual complex UID
- reactome_id: The REACTION ID (not entity!)
- component_id: The component ID
- component_id_or_reference_entity_id: Resolved reference
- input_or_output_uid: If component is a nested UID
- input_or_output_reactome_id: If component is a simple entity
```

## Proposed Solution

Add a new column `parent_entity_set_id` to track EntitySet lineage:

```python
{
    "uid": "abc123...",
    "reactome_id": 69598,  # reaction ID
    "component_id": 9943706,
    "component_id_or_reference_entity_id": 9943706,
    "input_or_output_uid": None,
    "input_or_output_reactome_id": 9943706,
    "parent_entity_set_id": 9943734  # NEW: which EntitySet this came from
}
```

## Implementation Plan

### 1. Update DataFrame Schema

**File:** `src/reaction_generator.py`
**Line:** ~34

```python
decomposed_uid_mapping = pd.DataFrame(
    columns=[
        "uid",
        "reactome_id",
        "component_id",
        "component_id_or_reference_entity_id",
        "input_or_output_uid",
        "input_or_output_reactome_id",
        "parent_entity_set_id",  # NEW COLUMN
    ]
)
```

### 2. Modify `break_apart_entity` Function

Need to pass parent EntitySet ID through the recursion:

```python
def break_apart_entity(entity_id: int, parent_set_id: Optional[int] = None) -> Set[str]:
    """Break apart entity, tracking which EntitySet (if any) it came from."""

    if "EntitySet" in labels:
        # When decomposing an EntitySet, pass its ID as the parent
        for member_id in member_ids:
            members = break_apart_entity(member_id, parent_set_id=entity_id)  # Pass EntitySet ID
            ...
```

### 3. Update Row Creation

**Locations:**
- `get_broken_apart_ids()` - Lines 116-138
- `get_uids_for_iterproduct_components()` - Lines 166-187

Add `parent_entity_set_id` to every row dict:

```python
row = {
    "uid": uid,
    "component_id": member,
    "reactome_id": reactome_id,
    "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(member),
    "input_or_output_uid": None,
    "input_or_output_reactome_id": member,
    "parent_entity_set_id": parent_set_id  # NEW
}
```

### 4. Update All Call Sites

Every call to `break_apart_entity` needs to handle the new return structure or pass parent info:
- `get_reaction_inputs()` - Line ~358
- `get_reaction_outputs()` - Line ~375
- Complex decomposition - Line ~291

### 5. Update Reconstruction Logic

With this information, reconstruction becomes:

```python
# Get components from generated data
components = [9943706, 9943732]

# Check if they share a parent EntitySet
parent_sets = decomposed[decomposed['component_id'].isin(components)]['parent_entity_set_id'].unique()

if len(parent_sets) == 1 and pd.notna(parent_sets[0]):
    # These came from an EntitySet, use the parent ID
    original_entity_id = int(parent_sets[0])  # 9943734
else:
    # These are independent entities
    original_entity_ids = components
```

## Files to Modify

1. **src/reaction_generator.py**
   - Line 34: Add column to DataFrame schema
   - Line 233: Modify `break_apart_entity()` signature
   - Line 268: Pass parent when decomposing EntitySets
   - Lines 116-138, 166-187: Add field to row dicts

2. **Tests** (update expected DataFrames):
   - tests/test_uuid_mapping_export.py
   - tests/test_and_or_logic.py
   - tests/test_transformation_semantics.py
   - tests/test_uuid_position_bug.py

## Expected Results

After this fix:
- **Perfect reconstruction:** Should go from 50% → ~90%+
- **EntitySet tracking:** Full traceability from member → parent EntitySet
- **Backward compatible:** Cells without EntitySet parents have NULL/NaN

## Testing Strategy

1. Unit tests: Verify `parent_entity_set_id` is populated correctly
2. Integration test: Reconstruct pathway 69620, expect 90%+ match rate
3. Regression test: Existing functionality unchanged (simple entities, complexes)
