# Quick Wins: Improvements You Can Make Today

These are simple, high-impact improvements that take <2 hours total.

## 1. Remove Debug Print Statements (5 minutes)

### Find them:
```bash
grep -n "print(" src/logic_network_generator.py
```

### Remove these lines:
- Line 48: `print("row")`
- Line 49: `print(row)`
- Line 34: `print("Checking best_matches contents:")`

### Why: Professional code shouldn't have print statements

---

## 2. Update README with Test Instructions (5 minutes)

Add this section to `README.md`:

```markdown
## Testing

Run the test suite:
```bash
poetry run pytest tests/ -v
```

Run with coverage report:
```bash
poetry run pytest tests/ --cov=src --cov-report=html
open htmlcov/index.html
```

Run specific test file:
```bash
poetry run pytest tests/test_and_or_logic.py -v
```

### Test Suite

- **34 tests** covering core functionality
- Tests for AND/OR logic, transformations, network invariants
- See `TEST_SUITE_SUMMARY.md` for details
```

### Why: Makes it easy for others to run tests

---

## 3. Add GitHub Actions CI (15 minutes)

Create `.github/workflows/test.yml`:

```yaml
name: Tests

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.12'

    - name: Install Poetry
      run: pip install poetry

    - name: Install dependencies
      run: poetry install

    - name: Run tests
      run: poetry run pytest tests/ -v

    - name: Run type checking
      run: poetry run mypy --ignore-missing-imports src/
      continue-on-error: true  # Don't fail build yet
```

### Why: Automatically runs tests on every commit

---

## 4. Add Type Hints to Main Function (20 minutes)

Edit `src/logic_network_generator.py`:

```python
# Before (line 418):
def create_pathway_logic_network(
    decomposed_uid_mapping,
    reaction_connections,
    best_matches,
):

# After:
from typing import Any
import pandas as pd

def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: pd.DataFrame,
) -> pd.DataFrame:
    """Create a pathway logic network from decomposed UID mappings.

    Args:
        decomposed_uid_mapping: Mapping from hashes to molecules
        reaction_connections: Connections between reactions
        best_matches: Pairings of input/output hashes

    Returns:
        DataFrame representing the logic network

    Raises:
        ValueError: If input DataFrames are empty or invalid
    """
```

### Why: Better IDE support, catches bugs earlier

---

## 5. Add Input Validation (15 minutes)

Add to `create_pathway_logic_network` at the start:

```python
def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: pd.DataFrame,
) -> pd.DataFrame:
    """..."""

    # Validate inputs
    if decomposed_uid_mapping.empty:
        raise ValueError("decomposed_uid_mapping cannot be empty")

    required_cols = {'uid', 'reactome_id', 'input_or_output_reactome_id'}
    missing = required_cols - set(decomposed_uid_mapping.columns)
    if missing:
        raise ValueError(
            f"decomposed_uid_mapping missing required columns: {missing}"
        )

    if best_matches.empty:
        raise ValueError("best_matches cannot be empty")

    # Continue with rest of function...
```

### Why: Better error messages, catch problems early

---

## 6. Rename Confusing Variables (30 minutes)

In `_add_pathway_connections` (line 270):

```python
# Before:
def _add_pathway_connections(
    input_uuids: List[str],
    output_uuids: List[str],
    ...
):
    for input_uuid in input_uuids:
        for output_uuid in output_uuids:
            pathway_logic_network_data.append({
                "source_id": input_uuid,
                "target_id": output_uuid,
                ...
            })

# After:
def _add_pathway_connections(
    reactant_molecule_uuids: List[str],  # Clearer: molecules consumed
    product_molecule_uuids: List[str],    # Clearer: molecules produced
    and_or: str,
    edge_type: str,
    pathway_logic_network_data: List[Dict[str, Any]]
) -> None:
    """Add edges representing biochemical transformations.

    Creates edges from reactant molecules to product molecules,
    representing transformations within reactions.
    """
    for reactant_uuid in reactant_molecule_uuids:
        for product_uuid in product_molecule_uuids:
            pathway_logic_network_data.append({
                "source_id": reactant_uuid,   # Reactant (consumed)
                "target_id": product_uuid,     # Product (produced)
                "pos_neg": "pos",
                "and_or": and_or,
                "edge_type": edge_type,
            })
```

**Also update the call site** (line 353):

```python
# Before:
_add_pathway_connections(
    input_uuids, output_uuids, and_or, edge_type, pathway_logic_network_data
)

# After:
_add_pathway_connections(
    reactant_molecule_uuids=input_uuids,  # Current reaction's inputs
    product_molecule_uuids=output_uuids,   # Preceding reaction's outputs
    and_or=and_or,
    edge_type=edge_type,
    pathway_logic_network_data=pathway_logic_network_data
)
```

### Why: Self-documenting code, matches terminology in papers/docs

---

## 7. Add .gitignore Entries (2 minutes)

Add to `.gitignore`:

```
# Test artifacts
.pytest_cache/
.coverage
htmlcov/
*.coverage

# IDE
.vscode/
.idea/
*.swp

# Python
__pycache__/
*.pyc
*.pyo
*.pyd
.Python
*.egg-info/

# OS
.DS_Store
Thumbs.db

# Temporary files
*.tmp
*.bak
debug_log.txt
```

### Why: Keeps repo clean

---

## 8. Add Coverage Configuration (5 minutes)

Add to `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
addopts = [
    "--verbose",
    "--strict-markers",
]

[tool.coverage.run]
source = ["src"]
omit = [
    "*/tests/*",
    "*/test_*.py",
]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
    "if __name__ == .__main__.:",
    "if TYPE_CHECKING:",
]
```

### Why: Better test configuration, coverage reporting

---

## 9. Document Key Functions (20 minutes)

Add docstrings to these functions:

### `_determine_edge_properties` (line 249):

```python
def _determine_edge_properties(num_preceding_reactions: int) -> tuple:
    """Determine AND/OR logic and edge type.

    Logic:
    - Single source (num_preceding == 1) → AND relationship (required)
    - Multiple sources (num_preceding > 1) → OR relationship (alternatives)

    This implements the user requirement:
    - R1→A (OR), R2→A (OR) when multiple sources feed same molecule
    - A→R3 (AND) for any molecule going into reaction

    Args:
        num_preceding_reactions: Number of reactions feeding into current one

    Returns:
        Tuple of (and_or, edge_type):
        - ('and', 'input') for single source
        - ('or', 'output') for multiple sources
    """
```

### `extract_inputs_and_outputs` (line 289):

```python
def extract_inputs_and_outputs(
    reaction_uid: str,
    reaction_uids: List[str],
    uid_reaction_connections: pd.DataFrame,
    reaction_id_map: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
    pathway_logic_network_data: List[Dict[str, Any]],
) -> None:
    """Extract inputs and outputs for reactions and create transformation edges.

    This function creates edges representing biochemical transformations
    WITHIN each reaction (not connections BETWEEN reactions).

    For each reaction:
    1. Get terminal molecules from inputs (reactants)
    2. Get terminal molecules from outputs (products)
    3. Create edges: reactants → products
    4. Assign AND/OR logic based on number of preceding reactions

    Reactions connect IMPLICITLY through shared molecules:
    - Molecule X is output from Reaction 1 (appears as target)
    - Molecule X is input to Reaction 2 (appears as source)
    - Result: X connects R1 and R2

    Args:
        reaction_uid: Current reaction being processed
        reaction_uids: List of all reactions to process
        uid_reaction_connections: Connections between reactions
        reaction_id_map: Mapping of reaction UIDs to hashes
        decomposed_uid_mapping: Mapping of hashes to molecules
        reactome_id_to_uuid: Cache of molecule UUIDs
        pathway_logic_network_data: Output list (modified in-place)
    """
```

### Why: Code is self-documenting, easier to understand

---

## Total Time: ~2 hours

These 9 improvements will significantly increase code quality with minimal effort:

- ✅ Remove debug code
- ✅ Add test documentation
- ✅ Set up CI
- ✅ Add type hints
- ✅ Add validation
- ✅ Rename confusing variables
- ✅ Clean up .gitignore
- ✅ Configure coverage
- ✅ Document key functions

## After These Changes

Your code will:
- ✅ Run tests automatically on every commit (CI)
- ✅ Have better error messages (validation)
- ✅ Be easier to understand (clear names, docstrings)
- ✅ Be more professional (no debug prints)
- ✅ Have IDE support (type hints)

## Next Steps

After these quick wins, see `IMPROVEMENT_RECOMMENDATIONS.md` for:
- Comprehensive refactoring
- Additional testing
- Architecture documentation
- Performance optimization
