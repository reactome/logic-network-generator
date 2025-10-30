# Test-Based Analysis of Edge Direction

## Test Suite Created

1. **Unit tests** (`test_logic_network_generator.py`): ✅ All 9 tests pass
   - `_assign_uuids`: Correctly creates/reuses UUIDs for Reactome IDs
   - `_determine_edge_properties`: Correctly returns AND/OR based on preceding reaction count
   - `_add_pathway_connections`: Creates cartesian product of input×output edges

2. **Integration tests** (`test_edge_direction_integration.py`): ✅ Tests pass
   - Synthetic pathway test: R1 → R2 with shared molecule
   - **Result**: Creates self-loop edges (MolA → MolA)
   - **Conclusion**: When the same molecule appears in connected reactions, we get self-loops

3. **Real data analysis** (`test_actual_edge_semantics.py`): ✅ Test passes
   - Analyzed actual pathway_logic_network_69620.csv
   - **Critical Finding**: **ZERO self-loop edges** in real data!

## Key Discoveries

### Discovery 1: Real Data Has No Self-Loops

```
Total main pathway edges: 4,995
Self-loop edges: 0
Non-self-loop edges: 4,995
```

**All edges connect DIFFERENT molecules.**

### Discovery 2: Clear Directional Flow

```
Node Analysis:
- Sources only (never targets): 9 molecules
- Targets only (never sources): 11 molecules
- Both source and target: 2 molecules
```

This pattern strongly suggests **correct forward flow**: `roots → intermediates → terminals`

### Discovery 3: Contradiction with Synthetic Test

**Synthetic test** (R1 outputs MolA, R2 inputs MolA):
- Result: Self-loop (MolA → MolA)

**Real pathway data**:
- Result: No self-loops at all

**Implication**: The synthetic test doesn't accurately model real pathway structure.

## Why No Self-Loops in Real Data?

### Hypothesis 1: Different Molecules at Each Stage
Real reactions might transform molecules such that:
- R1 consumes A, produces B
- R2 consumes C, produces D
- Edges: A→B, C→D (no shared molecules)

But this doesn't explain pathway connectivity...

### Hypothesis 2: Decomposition Creates Distinct Representations
When complexes are decomposed:
- Complex1(A,B) → components A and B (with UIDs tied to Complex1)
- Complex2(A,C) → components A and C (with UIDs tied to Complex2)
- Even though both contain "A", they get different UUIDs because they're from different complexes

**This is more likely!** The decomposition process might create molecule representations that are context-dependent.

### Hypothesis 3: UUID Assignment Strategy
The `reactome_id_to_uuid` mapping might be more complex than assumed. Perhaps:
- Same Reactome ID in different contexts gets different UUIDs?
- Or the "input_or_output_reactome_id" values are already unique per context?

## Current Understanding: Edge Direction

Given the real data shows:
- **9 root inputs** (source only)
- **11 terminal outputs** (target only)
- **Clear forward flow pattern**

### Tentative Conclusion

**The edges appear to flow in the CORRECT direction** for biological pathway flow:
```
source_id (roots) → target_id (terminals)
```

However, we still don't fully understand:
1. Why synthetic test creates self-loops but real data doesn't
2. What causes edges between different molecules in real data
3. Whether the current code at line 281-282 (`source_id: input_uuid, target_id: output_uuid`) is semantically correct or backwards

## Recommended Next Steps

1. **Examine decomposed_uid_mapping structure** to understand how molecules get unique representations
2. **Trace through ONE real reaction pair** to see exactly which molecules get connected and why they're different
3. **Create better synthetic test** that matches real data structure (no self-loops)
4. **Add comprehensive documentation** explaining the data flow and edge semantics

## Test Files Created

- `tests/__init__.py`
- `tests/test_logic_network_generator.py` - Unit tests for helper functions
- `tests/test_edge_direction_integration.py` - Integration test with synthetic data
- `tests/test_actual_edge_semantics.py` - Analysis of real pathway data

All tests pass: `poetry run pytest tests/ -v`
