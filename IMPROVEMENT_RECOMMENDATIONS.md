# Repository Improvement Recommendations

## Priority 1: Critical for Quality 🔴

### 1. Clean Up Debug Code

**Issue**: Production code contains debug logging and print statements from investigation.

**Location**: `src/logic_network_generator.py` lines 300-357

```python
# Current (verbose debug logging):
logger.debug("\n" + "="*80)
logger.debug("INSTRUMENTATION: Starting extract_inputs_and_outputs")
logger.debug(f"Processing {len(reaction_uids)} reaction UIDs")
print("row")
print(row)
```

**Recommendation**:
- Remove or gate debug logging behind a flag
- Remove all `print()` statements
- Use proper logging levels (DEBUG, INFO, WARNING, ERROR)

**Impact**: Professional code, easier to read, better performance

---

### 2. Remove Global State

**Issue**: Global database connection creates testing/maintenance problems.

**Location**: `src/logic_network_generator.py` lines 9-10

```python
# Current (global):
uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))
```

**Recommendation**:
```python
# Better: Dependency injection
class PathwayGenerator:
    def __init__(self, graph: Graph):
        self.graph = graph

    def create_pathway_logic_network(self, ...):
        # Use self.graph instead of global
```

**Benefits**:
- Testable (can inject mock database)
- Configurable (different databases for dev/prod)
- Thread-safe
- Follows best practices

---

### 3. Add Input Validation

**Issue**: No validation of inputs - can crash with confusing errors.

**Recommendation**:
```python
def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: Any,
) -> pd.DataFrame:
    """Create a pathway logic network from decomposed UID mappings."""

    # Validate inputs
    if decomposed_uid_mapping.empty:
        raise ValueError("decomposed_uid_mapping cannot be empty")

    required_cols = ['uid', 'reactome_id', 'input_or_output_reactome_id']
    missing = set(required_cols) - set(decomposed_uid_mapping.columns)
    if missing:
        raise ValueError(f"decomposed_uid_mapping missing columns: {missing}")

    # ... rest of function
```

**Impact**: Better error messages, easier debugging, prevents silent failures

---

### 4. Fix Confusing Variable Names

**Issue**: `input_uuid` and `output_uuid` suggest inter-reaction flow but actually represent intra-reaction transformations.

**Location**: `src/logic_network_generator.py` lines 270-286, 340-354

**Recommendation**:
```python
# Current (confusing):
def _add_pathway_connections(
    input_uuids: List[str],    # Unclear
    output_uuids: List[str],   # Unclear
    ...
):
    for input_uuid in input_uuids:
        for output_uuid in output_uuids:
            pathway_logic_network_data.append({
                "source_id": input_uuid,
                "target_id": output_uuid,
                ...
            })

# Better (clear):
def _add_transformation_edges(
    reactant_molecule_uuids: List[str],  # What goes in
    product_molecule_uuids: List[str],    # What comes out
    and_or: str,
    edge_type: str,
    pathway_logic_network_data: List[Dict[str, Any]]
) -> None:
    """Add edges representing biochemical transformations.

    Creates directed edges from reactant molecules to product molecules,
    representing the transformation that occurs within a reaction.

    Args:
        reactant_molecule_uuids: Molecules consumed (inputs to reaction)
        product_molecule_uuids: Molecules produced (outputs from reaction)
        ...
    """
    for reactant_uuid in reactant_molecule_uuids:
        for product_uuid in product_molecule_uuids:
            pathway_logic_network_data.append({
                "source_id": reactant_uuid,  # Reactant (consumed)
                "target_id": product_uuid,    # Product (produced)
                "pos_neg": "pos",
                "and_or": and_or,
                "edge_type": edge_type,
            })
```

**Impact**: Code is self-documenting, easier to understand

---

## Priority 2: Important for Maintainability 🟡

### 5. Add Type Hints Everywhere

**Issue**: Many functions lack type hints, making code harder to understand.

**Current Coverage**: ~40% (estimated)
**Target**: 100%

**Example**:
```python
# Before:
def _get_reactome_id_from_hash(decomposed_uid_mapping, hash_value):
    return decomposed_uid_mapping.loc[
        decomposed_uid_mapping["uid"] == hash_value, "reactome_id"
    ].values[0]

# After:
def _get_reactome_id_from_hash(
    decomposed_uid_mapping: pd.DataFrame,
    hash_value: str
) -> int:
    """Extract reactome_id for a given hash from decomposed_uid_mapping.

    Args:
        decomposed_uid_mapping: DataFrame containing uid to reactome_id mappings
        hash_value: Hash string to look up

    Returns:
        Reactome ID as integer

    Raises:
        IndexError: If hash_value not found in mapping
    """
    result = decomposed_uid_mapping.loc[
        decomposed_uid_mapping["uid"] == hash_value, "reactome_id"
    ].values

    if len(result) == 0:
        raise ValueError(f"Hash not found in mapping: {hash_value}")

    return int(result[0])
```

**Benefits**:
- IDE autocomplete works better
- Catch bugs earlier (with mypy)
- Self-documenting code

---

### 6. Break Down Large Functions

**Issue**: Some functions do too much (50+ lines).

**Example**: `extract_inputs_and_outputs` (80+ lines) does:
1. Iterates through reactions
2. Extracts input/output information
3. Processes preceding reactions
4. Determines edge properties
5. Adds connections
6. Logs everything

**Recommendation**:
```python
# Split into focused functions:

def _process_reaction_pair(
    current_reaction_uid: str,
    preceding_reaction_uid: str,
    reaction_id_map: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
) -> List[Dict[str, Any]]:
    """Process a single pair of connected reactions.

    Returns edges representing the transformation.
    """
    # Extract molecules
    input_molecules = _extract_terminal_molecules(...)
    output_molecules = _extract_terminal_molecules(...)

    # Determine logic
    and_or, edge_type = _determine_edge_properties(...)

    # Create edges
    return _create_transformation_edges(
        input_molecules, output_molecules, and_or, edge_type
    )

def extract_inputs_and_outputs(...):
    """Main orchestration - delegates to helper functions."""
    for reaction_uid in reaction_uids:
        preceding_uids = _get_preceding_reactions(...)

        for preceding_uid in preceding_uids:
            edges = _process_reaction_pair(
                reaction_uid, preceding_uid, ...
            )
            pathway_logic_network_data.extend(edges)
```

**Benefits**:
- Easier to test (test individual pieces)
- Easier to understand (clear responsibilities)
- Easier to modify (change one piece without affecting others)

---

### 7. Add Comprehensive Docstrings

**Issue**: Many functions lack docstrings explaining their purpose and data structures.

**Recommendation**: Use numpy/Google style docstrings:

```python
def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: pd.DataFrame,
) -> pd.DataFrame:
    """Create a pathway logic network from Reactome data.

    This function generates a directed graph representing biochemical pathways
    where:
    - Nodes are molecules (identified by UUIDs)
    - Edges are transformations within reactions (input → output)
    - AND/OR logic indicates whether multiple sources are alternatives

    The network is suitable for perturbation analysis and pathway flow studies.

    Args:
        decomposed_uid_mapping: DataFrame with columns:
            - uid: Hash of molecule combination
            - reactome_id: Biological reaction ID
            - input_or_output_reactome_id: Terminal molecule ID
        reaction_connections: DataFrame with columns:
            - preceding_reaction_id: Upstream reaction
            - following_reaction_id: Downstream reaction
        best_matches: DataFrame with columns:
            - incomming: Input hash (within reaction)
            - outgoing: Output hash (within reaction)

    Returns:
        DataFrame representing the logic network with columns:
            - source_id: UUID of input molecule (reactant)
            - target_id: UUID of output molecule (product)
            - and_or: Logic type ('and' or 'or')
            - edge_type: Edge category ('input', 'output', 'catalyst', etc.)
            - pos_neg: Positive or negative regulation

    Raises:
        ValueError: If input DataFrames are empty or missing required columns

    Examples:
        >>> mapping = pd.read_csv('decomposed_uid_mapping.csv')
        >>> connections = pd.read_csv('reaction_connections.csv')
        >>> matches = pd.read_csv('best_matches.csv')
        >>> network = create_pathway_logic_network(mapping, connections, matches)
        >>> print(f"Created network with {len(network)} edges")

    Notes:
        - Edges represent transformations within reactions, not connections
          between reactions
        - Reactions connect implicitly through shared molecules
        - No self-loops in the network (reactions transform molecules)
        - Root inputs appear only as sources, terminal outputs only as targets
    """
    # ... implementation
```

**Impact**: Self-documenting code, easier onboarding for new developers

---

### 8. Set Up CI/CD Pipeline

**Issue**: No automated testing on commits/PRs.

**Recommendation**: Create `.github/workflows/test.yml`:

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.9", "3.10", "3.11", "3.12"]

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}

    - name: Install Poetry
      run: pip install poetry

    - name: Install dependencies
      run: poetry install

    - name: Run tests
      run: poetry run pytest tests/ -v --cov=src --cov-report=xml

    - name: Upload coverage
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml

    - name: Run type checking
      run: poetry run mypy src/

    - name: Run linting
      run: poetry run ruff check src/
```

**Benefits**:
- Catch bugs before they're merged
- Ensure tests pass on all Python versions
- Track code coverage over time
- Enforce code quality standards

---

### 9. Add Code Coverage Reporting

**Current**: Unknown coverage
**Target**: >80%

**Setup**:
```bash
poetry add --group dev pytest-cov
poetry run pytest tests/ --cov=src --cov-report=html
```

**Add to CI** (see #8 above)

**Benefits**:
- Identify untested code
- Track coverage trends
- Ensure new code is tested

---

## Priority 3: Nice to Have 🟢

### 10. Add More Comprehensive Tests

**Current Coverage Gaps**:
- Decomposition logic (`src/reaction_generator.py`)
- Best matching algorithm (`src/best_reaction_match.py`)
- Neo4j query functions (`src/neo4j_connector.py`)
- Catalyst/regulator logic
- Edge cases (empty inputs, malformed data, etc.)

**Recommendation**:
```python
# tests/test_decomposition.py
class TestSetDecomposition:
    def test_simple_set_breaks_into_components(self):
        """EntitySet(A,B,C) should decompose into [A, B, C]."""
        # ...

    def test_nested_set_recursive_decomposition(self):
        """EntitySet(A, EntitySet(B,C)) should fully decompose."""
        # ...

    def test_complex_with_sets_combinatorial(self):
        """Complex(EntitySet(A,B), C) should create combinations."""
        # ...

# tests/test_neo4j_queries.py (with mock database)
class TestNeo4jQueries:
    def test_get_reaction_connections_returns_expected_structure(self):
        # ...

    def test_handles_reactions_with_no_preceding(self):
        # ...
```

**Target**: 80%+ code coverage

---

### 11. Add Performance Benchmarks

**Issue**: No baseline for performance monitoring.

**Recommendation**:
```python
# tests/test_performance.py
import pytest
import time

class TestPerformance:
    def test_pathway_generation_time(self):
        """Pathway 69620 should generate in <5 seconds."""
        start = time.time()

        # Generate pathway
        result = create_pathway_logic_network(...)

        elapsed = time.time() - start
        assert elapsed < 5.0, f"Took {elapsed:.2f}s (expected <5s)"

    @pytest.mark.parametrize("pathway_id", [69620, 68875, ...])
    def test_multiple_pathways(self, pathway_id):
        """All pathways should generate without errors."""
        result = create_pathway_logic_network(...)
        assert len(result) > 0
```

**Benefits**:
- Detect performance regressions
- Optimize slow code
- Set SLAs for generation time

---

### 12. Add Architecture Documentation

**Create**: `docs/ARCHITECTURE.md`

```markdown
# Architecture

## Overview

The logic network generator transforms Reactome pathway data into
logic networks suitable for perturbation analysis.

## Data Flow

```
Reactome DB (Neo4j)
    ↓ (query)
reaction_connections.csv
    ↓ (decompose)
decomposed_uid_mapping.csv
    ↓ (match)
best_matches.csv
    ↓ (generate)
pathway_logic_network.csv
```

## Components

### 1. Neo4j Connector (`neo4j_connector.py`)
- Queries Reactome database
- Extracts reaction connections
- Gets entity components

### 2. Reaction Generator (`reaction_generator.py`)
- Decomposes complexes and sets
- Creates combinatorial expansions
- Generates hash-based UIDs

### 3. Best Match Algorithm (`best_reaction_match.py`)
- Pairs input/output combinations
- Uses Hungarian algorithm
- Maximizes molecule overlap

### 4. Logic Network Generator (`logic_network_generator.py`)
- Creates molecule-to-molecule edges
- Assigns AND/OR logic
- Adds catalysts and regulators

## Key Concepts

### Transformations Within Reactions
Edges represent transformations WITHIN reactions, not connections
BETWEEN reactions. See COMPLETE_UNDERSTANDING.md for details.

### AND/OR Logic
- Single source → AND (required)
- Multiple sources → OR (alternatives)

### No Self-Loops
Reactions transform molecules, so inputs ≠ outputs, therefore
no self-loops in the network.
```

---

### 13. Improve Error Handling

**Issue**: Limited error handling and recovery.

**Recommendation**:
```python
# Custom exceptions
class LogicNetworkError(Exception):
    """Base exception for logic network generation."""
    pass

class InvalidMappingError(LogicNetworkError):
    """Raised when decomposed_uid_mapping is invalid."""
    pass

class DatabaseConnectionError(LogicNetworkError):
    """Raised when cannot connect to Neo4j."""
    pass

# Use in code
def create_pathway_logic_network(...):
    try:
        # Validate inputs
        _validate_inputs(decomposed_uid_mapping, ...)

        # Generate network
        result = _generate_network(...)

        return result

    except pd.errors.EmptyDataError as e:
        raise InvalidMappingError(
            "decomposed_uid_mapping is empty or malformed"
        ) from e
    except Exception as e:
        logger.error(f"Failed to generate pathway: {e}")
        raise LogicNetworkError(
            f"Network generation failed: {e}"
        ) from e
```

**Benefits**:
- Better error messages
- Easier debugging
- Graceful failure modes

---

### 14. Add Configuration Management

**Issue**: Hard-coded values scattered through code.

**Recommendation**: Create `config.py`:

```python
from dataclasses import dataclass
from typing import Optional
import os

@dataclass
class Config:
    """Configuration for logic network generator."""

    # Neo4j connection
    neo4j_uri: str = "bolt://localhost:7687"
    neo4j_user: str = "neo4j"
    neo4j_password: str = "test"

    # Generation settings
    max_decomposition_depth: int = 10
    cache_intermediate_results: bool = True
    output_directory: str = "output"

    # Logging
    log_level: str = "INFO"
    debug_instrumentation: bool = False

    @classmethod
    def from_env(cls) -> 'Config':
        """Load configuration from environment variables."""
        return cls(
            neo4j_uri=os.getenv("NEO4J_URI", cls.neo4j_uri),
            neo4j_user=os.getenv("NEO4J_USER", cls.neo4j_user),
            neo4j_password=os.getenv("NEO4J_PASSWORD", cls.neo4j_password),
            log_level=os.getenv("LOG_LEVEL", cls.log_level),
            debug_instrumentation=os.getenv("DEBUG", "false").lower() == "true",
        )

# Usage
config = Config.from_env()
graph = Graph(config.neo4j_uri, auth=(config.neo4j_user, config.neo4j_password))
```

**Benefits**:
- Easy to configure for different environments
- No hard-coded values
- Environment variable support

---

### 15. Add Examples and Tutorials

**Create**: `examples/` directory

```python
# examples/basic_usage.py
"""
Basic usage example for logic network generator.

This example shows how to generate a logic network for a single pathway.
"""

from src.logic_network_generator import create_pathway_logic_network
from src.pathway_generator import generate_pathway_file
import pandas as pd

# Generate pathway 69620 (Jak-STAT signaling)
print("Generating pathway 69620...")
generate_pathway_file(
    pathway_id="69620",
    taxon_id="9606",  # Homo sapiens
    pathway_name="Jak-STAT signaling pathway"
)

# Load the generated data
decomposed = pd.read_csv("decomposed_uid_mapping_69620.csv")
connections = pd.read_csv("reaction_connections_69620.csv")
matches = pd.read_csv("best_matches_69620.csv")

# Create logic network
network = create_pathway_logic_network(decomposed, connections, matches)

# Analyze results
print(f"\nGenerated network with {len(network)} edges")

main_edges = network[~network['edge_type'].isin(['catalyst', 'regulator'])]
print(f"Main pathway edges: {len(main_edges)}")

sources = set(main_edges['source_id'].unique())
targets = set(main_edges['target_id'].unique())
roots = sources - targets
terminals = targets - sources

print(f"Root inputs: {len(roots)}")
print(f"Terminal outputs: {len(terminals)}")
print(f"Intermediate molecules: {len(sources & targets)}")

# Save network
network.to_csv("pathway_logic_network_69620.csv", index=False)
print("\nNetwork saved to pathway_logic_network_69620.csv")
```

---

## Implementation Priority

### Phase 1 (Week 1): Critical Cleanup
1. Remove debug code
2. Fix confusing variable names
3. Add input validation
4. Clean up print statements

### Phase 2 (Week 2): Infrastructure
5. Set up CI/CD
6. Add code coverage
7. Remove global state
8. Add configuration management

### Phase 3 (Week 3): Documentation
9. Add comprehensive docstrings
10. Create architecture documentation
11. Add examples and tutorials

### Phase 4 (Ongoing): Testing & Quality
12. Add missing tests (target 80%+ coverage)
13. Add performance benchmarks
14. Improve error handling
15. Add type hints everywhere

---

## Metrics to Track

**Code Quality:**
- [ ] Type hint coverage: 100%
- [ ] Test coverage: >80%
- [ ] Docstring coverage: 100% of public functions
- [ ] No print statements in production code
- [ ] No global state

**Performance:**
- [ ] Pathway generation: <5s for typical pathway
- [ ] Memory usage: <2GB for large pathways
- [ ] Test suite: <10s total runtime

**Maintainability:**
- [ ] Average function length: <30 lines
- [ ] Cyclomatic complexity: <10
- [ ] Code duplication: <5%

---

## Quick Wins (Can Do Today)

1. **Remove print statements** (5 minutes)
   ```bash
   # Find all print statements
   grep -r "print(" src/
   # Remove them
   ```

2. **Add type hints to main functions** (30 minutes)
   - Start with `create_pathway_logic_network`
   - Add to `extract_inputs_and_outputs`

3. **Set up basic CI** (30 minutes)
   - Copy GitHub Actions workflow above
   - Commit and push

4. **Add input validation** (15 minutes)
   - Add to `create_pathway_logic_network`
   - Check for empty DataFrames

5. **Update README with test instructions** (10 minutes)
   ```markdown
   ## Testing

   Run tests:
   ```bash
   poetry run pytest tests/ -v
   ```

   With coverage:
   ```bash
   poetry run pytest tests/ --cov=src
   ```
   ```

**Total Time**: ~90 minutes for significant quality improvement!

---

## Long-Term Vision

**Goal**: Production-ready, maintainable, well-documented codebase

**Success Criteria:**
- ✅ 80%+ test coverage
- ✅ CI/CD pipeline running
- ✅ Comprehensive documentation
- ✅ No confusing variable names
- ✅ Type hints everywhere
- ✅ Easy for new developers to understand
- ✅ Performance benchmarks established
- ✅ Error handling is robust

**Benefits:**
- Faster development (less debugging)
- Easier collaboration (clear code)
- Fewer bugs (better testing)
- Better performance (benchmarks)
- Professional quality (CI/CD)
