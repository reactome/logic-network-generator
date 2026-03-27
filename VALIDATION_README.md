# Pathway Logic Network Validation System

## Overview

Comprehensive validation system that verifies the correctness of generated logic networks by comparing them against the source Neo4j database.

## What It Validates

### 1. **Completeness Checks**
- ✅ All reactions from pathway are present
- ✅ All physical entities are accounted for
- ✅ All reaction connections are preserved
- ✅ All regulators and catalysts are included

### 2. **Correctness Checks**
- ✅ UUID mapping covers all UUIDs in logic network
- ✅ No orphaned UUIDs (unused mappings)
- ✅ Logic network has valid structure (columns, data types)
- ✅ Position-aware UUIDs working (same entity at different positions has different UUIDs)

### 3. **Integrity Checks**
- ✅ No excessive self-loops in main pathway (with position-aware UUIDs)
- ✅ Decomposition preserves information
- ✅ Reaction connections match database

### 4. **Statistics**
- 📊 Comprehensive summary comparing DB vs generated files
- 📊 Position-aware UUID effectiveness metrics
- 📊 Coverage percentages for all validations

## Usage

### Quick Validation (Recommended)
Run validation on the default pathway (69620):

```bash
poetry run python validate_pathway.py
```

### Validate Specific Pathway
```bash
poetry run python validate_pathway.py <pathway_id>
```

Example:
```bash
poetry run python validate_pathway.py 1257604
```

### Run Individual Tests
```bash
# Run all validation tests
poetry run pytest tests/test_pathway_validation.py -v -s

# Run specific validation
poetry run pytest tests/test_pathway_validation.py::TestPathwayValidation::test_all_reactions_present -v -s

# Run with summary statistics
poetry run pytest tests/test_pathway_validation.py::TestPathwayValidation::test_summary_statistics -v -s
```

## What Gets Validated

### Input: Database Pathway
- Queries Neo4j database for pathway structure
- Extracts reactions, entities, connections, regulators

### Generated Files (in `output/` directory)
- `output/pathway_logic_network_<pathway_id>.csv` - Main logic network
- `output/uuid_mapping_<pathway_id>.csv` - UUID to Reactome ID mapping
- `output/decomposed_uid_mapping_<pathway_id>.csv` - Decomposition details
- `output/reaction_connections_<pathway_id>.csv` - Reaction connectivity

### Validation Tests

#### Test 1: `test_all_reactions_present`
Verifies all reactions from the database pathway are in the generated reaction_connections file.

**What it checks:**
- Queries DB for all reactions in pathway
- Compares with reactions in generated files
- Reports missing or extra reactions

**Expected:** All DB reactions should be present (100% coverage)

#### Test 2: `test_all_physical_entities_have_uuids`
Verifies all physical entities from reactions have UUID mappings.

**What it checks:**
- Extracts entities from DB
- Checks if they appear in UUID mapping or decomposed mapping
- Accounts for decomposition (sets/complexes)

**Expected:** All entities should be accounted for

#### Test 3: `test_reaction_connections_are_complete`
Verifies reaction connections match database relationships.

**What it checks:**
- Queries DB for reaction→entity→reaction connections
- Compares with generated reaction_connections
- Calculates coverage percentage

**Expected:** >70% coverage (some differences due to decomposition/matching)

#### Test 4: `test_uuid_mapping_completeness`
Verifies UUID mapping covers all UUIDs used in logic network.

**What it checks:**
- Extracts all UUIDs from logic network edges
- Checks if all are in UUID mapping file
- Reports any unmapped UUIDs

**Expected:** 100% coverage - no unmapped UUIDs

#### Test 5: `test_no_orphaned_uuids_in_mapping`
Checks for UUIDs in mapping that aren't used in logic network.

**What it checks:**
- Finds UUIDs in mapping not used in network
- Calculates usage rate
- Reports orphaned UUIDs

**Expected:** High usage rate (>80%), some orphans are OK (terminal entities)

#### Test 6: `test_logic_network_has_valid_structure`
Validates basic structure and data integrity.

**What it checks:**
- All required columns present
- No null values in critical columns
- Valid values for categorical columns (pos_neg, and_or, edge_type)

**Expected:** All structural checks pass

#### Test 7: `test_position_aware_uuids_working`
Validates the UUID position bug fix is working.

**What it checks:**
- Finds entities appearing at multiple positions
- Verifies each position has a unique UUID
- Reports multi-position entities

**Expected:** Each position has unique UUID (this validates the fix!)

#### Test 8: `test_regulators_present`
Verifies regulators from database are in logic network.

**What it checks:**
- Queries DB for all regulators
- Counts regulator/catalyst edges in logic network
- Ensures regulatory edges exist if DB has regulators

**Expected:** Regulator edges present if DB has regulators

#### Test 9: `test_no_self_loops_in_main_pathway`
Validates position-aware UUIDs eliminated most self-loops.

**What it checks:**
- Counts self-loops in main pathway edges
- Calculates self-loop ratio
- Verifies it's very low (<5%)

**Expected:** Very few self-loops with position-aware UUIDs

#### Test 10: `test_decomposition_preserves_information`
Validates complexes and sets are properly decomposed.

**What it checks:**
- Queries DB for all complexes and entity sets
- Checks if they appear in decomposed_mapping
- Calculates decomposition coverage

**Expected:** >50% coverage (some may not be in active connections)

#### Test 11: `test_summary_statistics`
Comprehensive summary comparing DB vs generated files.

**What it reports:**
- Pathway name and ID
- DB statistics (reactions, entities)
- Generated file statistics (edges, UUIDs, mappings)
- Position-aware UUID statistics
- Multi-position entity counts

**Expected:** Produces comprehensive summary for analysis

## Expected Runtime

- **Small pathways** (<50 reactions): 30-60 seconds
- **Medium pathways** (50-200 reactions): 1-3 minutes
- **Large pathways** (>200 reactions): 3-10 minutes

Runtime includes:
- Database queries
- Logic network generation
- File I/O
- Validation checks

## Interpreting Results

### ✅ All Tests Pass
Logic network is valid and correctly represents the pathway!

### ⚠️ Coverage Warnings
- **Reaction connections <70%:** May indicate complex matching issues
- **Entity coverage <100%:** Check for missing decomposition
- **UUID usage <80%:** May indicate disconnected entities (could be OK)

### ❌ Test Failures
- **Missing reactions:** Critical - investigate database query or filters
- **Unmapped UUIDs:** Critical - UUID assignment bug
- **Self-loop ratio >5%:** Position-aware UUIDs may not be working
- **Invalid structure:** Critical - data corruption or generation bug

## Example Output

```
=================================================================================
PATHWAY VALIDATION SUMMARY - Pathway 69620
=================================================================================

Pathway: Pathway Name

Database Statistics:
  Reactions: 150
  Physical Entities: 300

Generated Files Statistics:
  Reaction Connections: 145
  Logic Network Edges: 500
  - Main pathway edges: 400
  - Catalyst edges: 75
  - Regulator edges: 25
  UUID Mappings: 320
  Unique UUIDs in network: 315

Position-Aware UUID Statistics:
  Entities at multiple positions: 45
  Total position instances: 120
  Average positions per multi-position entity: 2.7

=================================================================================
```

## Troubleshooting

### Database Connection Errors
```bash
# Check database is running
poetry run python -c "from py2neo import Graph; g = Graph('bolt://localhost:7687', auth=('neo4j', 'test')); print(g.run('RETURN 1').data())"
```

### Test Timeouts
- Increase pytest timeout: `pytest --timeout=300`
- Or run individual tests separately

### File Not Found Errors
- Ensure you're running from project root
- Check that pathway files were generated successfully

### Low Coverage Warnings
- Check pathway complexity (highly interconnected pathways may have complex matching)
- Verify decomposition settings
- Review database query results

## Files

- `tests/test_pathway_validation.py` - Main validation test suite
- `validate_pathway.py` - Convenience script for running validation
- `VALIDATION_README.md` - This file

## Benefits

1. **Confidence:** Know your logic networks are correct
2. **Bug Detection:** Catch issues early
3. **Regression Testing:** Ensure changes don't break correctness
4. **Documentation:** Understand pathway complexity
5. **Quality Metrics:** Track coverage and accuracy

## Future Enhancements

Potential additions:
- Validate edge directionality semantically
- Check for biological validity (e.g., impossible reactions)
- Compare multiple pathways for consistency
- Generate validation reports in HTML/PDF
- Automated regression testing in CI/CD

---

**Created:** 2025-11-11
**Purpose:** Validate logic network generation correctness
**Status:** Production Ready ✅
