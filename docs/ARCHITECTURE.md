# Architecture

## Overview

The Logic Network Generator transforms Reactome pathway data into directed logic networks suitable for perturbation analysis and pathway flow studies. The system decomposes complex biochemical structures (complexes and entity sets) into individual components and creates a network where edges represent biochemical transformations.

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                          Reactome Neo4j Database                     │
│                       (Biological Pathway Data)                      │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ Neo4j Queries
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                    reaction_connections_{pathway_id}.csv             │
│    (Connections between reactions: preceding → following)            │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ Decomposition
                                    │ (Break complexes/sets into components)
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                 decomposed_uid_mapping_{pathway_id}.csv              │
│  (Maps hashes to individual physical entities - proteins, etc.)     │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ Hungarian Algorithm
                                    │ (Optimal input/output pairing)
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                    best_matches_{pathway_id}.csv                     │
│        (Pairs of input/output combinations within reactions)         │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ Logic Network Generation
                                    │ (Create transformation edges)
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                    pathway_logic_network.csv                         │
│  (source_id → target_id edges with AND/OR logic annotations)        │
└─────────────────────────────────────────────────────────────────────┘
```

## Key Concepts

### 1. Physical Entities

In Reactome, a `:PhysicalEntity` represents any biological molecule or complex:
- Simple molecules (ATP, water)
- Proteins (individual gene products)
- Complexes (protein complexes like Complex(A,B,C))
- Entity sets (alternative molecules like EntitySet(IsoformA, IsoformB))

### 2. Decomposition

Complex structures are broken down into individual components:

```
Input: Complex(ProteinA, ProteinB, EntitySet(ATP, GTP))
                    ↓ decomposition
Output:
  - Combination 1: ProteinA, ProteinB, ATP
  - Combination 2: ProteinA, ProteinB, GTP
```

This creates all possible molecular combinations through cartesian product, preserving biological alternatives.

### 3. Virtual Reactions

A single biological reaction in Reactome may represent multiple transformations after decomposition:

```
Biological Reaction (Reactome ID: 12345):
  Inputs: Complex(A,B), ATP
  Outputs: Complex(A,B,P), ADP

After decomposition and best matching:
  Virtual Reaction 1 (UID: uuid-1, Reactome ID: 12345):
    input_hash: "hash-of-[A,B,ATP]"
    output_hash: "hash-of-[A,B,P,ADP]"

  Virtual Reaction 2 (UID: uuid-2, Reactome ID: 12345):
    input_hash: "hash-of-[A,B,ATP]"
    output_hash: "hash-of-[A,P,B,ADP]"
  ...
```

Each virtual reaction gets a unique UID (UUID v4) while preserving the link to the original Reactome reaction ID.

### 4. Edge Semantics

**CRITICAL**: Edges represent transformations WITHIN reactions, not connections BETWEEN reactions.

```
Reaction: ATP + Water → ADP + Phosphate

Creates 4 edges (cartesian product):
  ATP       → ADP
  ATP       → Phosphate
  Water     → ADP
  Water     → Phosphate
```

Reactions connect **implicitly** through shared physical entities:

```
Reaction 1: A → B (creates edge where B is target)
Reaction 2: B → C (creates edge where B is source)

Result: Pathway flow A → B → C (B connects the reactions)
```

**No self-loops** exist because reactions transform molecules (inputs ≠ outputs).

### 5. AND/OR Logic

The logic network assigns AND/OR relationships based on how many reactions produce the same physical entity:

**OR Relationship** (Multiple sources):
```
R1: Glycolysis → ATP
R2: Oxidative Phosphorylation → ATP
R3: ATP → Energy

For R3: ATP can come from R1 OR R2
Edges: R1→ATP (OR), R2→ATP (OR)
Then:  ATP→R3 (AND - ATP is required)
```

**AND Relationship** (Single source):
```
R1: Glucose → Glucose-6-Phosphate
R2: Glucose-6-Phosphate → ...

Only one source produces Glucose-6-Phosphate
Edge: R1→G6P (AND - required)
```

**Rule**:
- Multiple preceding reactions → OR (alternatives)
- Single preceding reaction → AND (required)
- All inputs to reactions are AND (required)

## Component Architecture

### Core Components

#### 1. `src/neo4j_connector.py`
**Purpose**: Query Reactome Neo4j database

**Key Functions**:
- `get_reaction_connections()`: Get preceding/following reaction pairs
- `get_catalysts_for_reaction()`: Get catalyst relationships
- `get_positive/negative_regulators_for_reaction()`: Get regulatory relationships

**Output**: Raw Reactome data as DataFrames

#### 2. `src/reaction_generator.py`
**Purpose**: Decompose complexes and sets into components

**Key Functions**:
- `get_decomposed_uid_mapping()`: Main decomposition orchestrator
- Handles complexes (using `itertools.product` for combinations)
- Handles entity sets (using `itertools.product` for alternatives)
- Recursively decomposes nested structures

**Output**: `decomposed_uid_mapping` with all molecular combinations

#### 3. `src/best_reaction_match.py`
**Purpose**: Pair input/output combinations optimally

**Algorithm**: Hungarian algorithm (optimal assignment)

**Input**: Input combinations and output combinations from same reaction

**Output**: `best_matches` DataFrame with optimal pairings

#### 4. `src/logic_network_generator.py`
**Purpose**: Generate the final logic network

**Key Functions**:
- `create_pathway_logic_network()`: Main orchestrator
- `create_reaction_id_map()`: Create virtual reactions from best_matches
- `extract_inputs_and_outputs()`: Create transformation edges
- `_determine_edge_properties()`: Assign AND/OR logic
- `_add_pathway_connections()`: Add edges with cartesian product
- `append_regulators()`: Add catalyst/regulator edges

**Output**: Logic network DataFrame with edges and logic annotations

### Bin Scripts

#### `bin/create-pathways.py`
**Purpose**: Command-line interface for generating pathways

**Usage**:
```bash
# Single pathway
poetry run python bin/create-pathways.py --pathway-id 69620

# Multiple pathways
poetry run python bin/create-pathways.py --pathway-list pathways.tsv
```

#### `bin/create-db-id-name-mapping-file.py`
**Purpose**: Create human-readable mapping of database IDs to names

## Network Properties

### Node Types
- **Root Inputs**: Physical entities that only appear as sources (pathway starting points)
- **Intermediate Entities**: Appear as both sources and targets (connect reactions)
- **Terminal Outputs**: Physical entities that only appear as targets (pathway endpoints)

### Edge Types
- **Main edges**: Transformation edges within reactions
  - `edge_type`: "input" (single source, AND) or "output" (multiple sources, OR)
  - `pos_neg`: "pos" (positive transformation)
  - `and_or`: "and" (required) or "or" (alternative)

- **Regulatory edges**: Catalysts and regulators
  - `edge_type`: "catalyst" or "regulator"
  - `pos_neg`: "pos" (positive regulation) or "neg" (negative regulation)
  - `and_or`: Empty (not applicable to regulation)

### Network Structure
- **Directed**: Edges have direction (source → target)
- **Acyclic**: No cycles in main transformation edges
- **Bipartite-like**: Entities and reactions connect through transformations
- **No self-loops**: Reactions always transform inputs to different outputs

## Testing Strategy

### Test Categories

1. **Unit Tests** (`tests/test_logic_network_generator.py`)
   - Individual helper functions
   - UUID assignment
   - Edge property determination

2. **Integration Tests** (`tests/test_edge_direction_integration.py`)
   - Multi-reaction pathways
   - End-to-end data flow

3. **Semantic Tests** (`tests/test_transformation_semantics.py`)
   - Cartesian product correctness
   - Edge direction validation
   - Transformation logic

4. **Invariant Tests** (`tests/test_network_invariants.py`)
   - No self-loops
   - Root inputs only as sources
   - Terminal outputs only as targets
   - AND/OR logic consistency

5. **Logic Tests** (`tests/test_and_or_logic.py`)
   - Multiple sources → OR
   - Single source → AND
   - User requirement validation

6. **Validation Tests** (`tests/test_input_validation.py`)
   - Empty DataFrame handling
   - Missing column detection
   - Error message clarity

### Test Coverage
- **43 tests** total (100% passing)
- Covers core functionality, edge semantics, and network properties
- See `TEST_SUITE_SUMMARY.md` for detailed breakdown

## Design Decisions

### Why Virtual Reactions?
- **Problem**: A biological reaction may have multiple input/output combinations after decomposition
- **Solution**: Create multiple "virtual reactions" representing each combination
- **Benefit**: Clean mapping from combinations to transformations

### Why Cartesian Product for Edges?
- **Problem**: How to represent transformation within a reaction with multiple inputs/outputs?
- **Solution**: Every input connects to every output (cartesian product)
- **Rationale**: Biochemically accurate - all reactants contribute to all products

### Why Implicit Reaction Connections?
- **Problem**: How do reactions connect in the network?
- **Solution**: Through shared physical entities (molecule appears as target in R1, source in R2)
- **Benefit**: Natural representation - pathways flow through molecules, not abstract connections

### Why AND/OR Based on Preceding Count?
- **User Requirement**: Multiple sources should be OR, inputs to reactions should be AND
- **Implementation**: Count preceding reactions - if >1 then OR, otherwise AND
- **Rationale**: Matches biological intuition (alternatives vs requirements)

## Performance Considerations

### Caching
- Files are cached: `reaction_connections_{id}.csv`, `decomposed_uid_mapping_{id}.csv`, `best_matches_{id}.csv`
- Subsequent runs reuse cached data
- UUID assignments cached in `reactome_id_to_uuid` dictionary

### Scalability
- Decomposition uses itertools.product (efficient for combinatorics)
- Hungarian algorithm is O(n³) but pathways are typically small (<1000 reactions)
- Pandas operations are vectorized where possible

### Typical Performance
- Small pathway (10-20 reactions): <1 second
- Medium pathway (100-200 reactions): 1-5 seconds
- Large pathway (500+ reactions): 5-30 seconds

## Future Improvements

See `IMPROVEMENT_RECOMMENDATIONS.md` for comprehensive list. Key areas:

1. **Remove global database connection** - Use dependency injection
2. **Add more comprehensive tests** - Decomposition logic, Neo4j queries
3. **Performance benchmarks** - Track generation time across versions
4. **Better error handling** - Graceful handling of edge cases

## References

- **Reactome Database**: https://reactome.org/
- **Test Suite Documentation**: `TEST_SUITE_SUMMARY.md`
- **Test Findings**: `TEST_FINDINGS.md`
- **Complete Understanding**: `COMPLETE_UNDERSTANDING.md`
- **Improvement Recommendations**: `IMPROVEMENT_RECOMMENDATIONS.md`
