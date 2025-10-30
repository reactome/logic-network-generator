# Logic Network Generator

[![Tests](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml/badge.svg)](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml)

Generate logic networks from Reactome pathways by decomposing sets and complexes into their individual components.

## Setup

### Prerequisites

- [Python 3](https://www.python.org/downloads/)
- [Poetry](https://python-poetry.org/)
- [Docker](https://www.docker.com/) (for Neo4j database)

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/reactome/logic-network-generator.git
   cd logic-network-generator
   ```

2. Install dependencies:

   ```bash
   poetry install
   ```

3. Start the Neo4j Reactome database:

   ```bash
   docker run -p 7474:7474 -p 7687:7687 \
     -e NEO4J_dbms_memory_heap_maxSize=8g \
     public.ecr.aws/reactome/graphdb:Release94
   ```

   **Note:** Replace `Release94` with the desired Reactome version.

   The database will be accessible at:
   - Neo4j Browser: http://localhost:7474
   - Bolt protocol: bolt://localhost:7687

## Usage

### Generate Pathway Logic Networks

Generate logic networks for pathways using a pathway ID:

```bash
poetry run python bin/create-pathways.py --pathway-id 69620
```

Or generate for multiple pathways using a pathway list file:

```bash
poetry run python bin/create-pathways.py --pathway-list pathway_list.tsv
```

The pathway list file should be tab-separated with columns: `id` and `pathway_name`.

### Create Database ID to Name Mapping

```bash
poetry run python bin/create-db-id-name-mapping-file.py
```

## Examples

The `examples/` directory contains complete working examples:

### Generate and Analyze a Pathway

```bash
poetry run python examples/generate_pathway_example.py
```

This example demonstrates:
- Generating a logic network for the Cell Cycle pathway
- Analyzing network properties (edges, nodes, logic relationships)
- Finding root inputs and terminal outputs
- Error handling and troubleshooting

See **[examples/README.md](examples/README.md)** for:
- Additional usage patterns
- Example pathways to try
- Cytoscape export
- Troubleshooting guide

## Testing

The project has a comprehensive test suite with 52 tests covering core functionality, AND/OR logic, transformation semantics, network invariants, and regulatory relationships.

### Run All Tests

```bash
poetry run pytest tests/ -v
```

### Run Tests with Coverage

```bash
poetry run pytest tests/ --cov=src --cov-report=html
```

View the coverage report:
```bash
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

### Run Specific Test Files

```bash
# Test AND/OR logic
poetry run pytest tests/test_and_or_logic.py -v

# Test input validation
poetry run pytest tests/test_input_validation.py -v

# Test network invariants
poetry run pytest tests/test_network_invariants.py -v

# Test transformation semantics
poetry run pytest tests/test_transformation_semantics.py -v
```

### Test Suite Overview

- **52 tests** total (100% passing)
- **Unit tests**: Core helper functions
- **Integration tests**: End-to-end pathway generation
- **Validation tests**: Input validation and error handling
- **Invariant tests**: Network structural properties
- **Semantics tests**: Transformation logic and edge direction
- **Regulatory tests**: Negative regulators, positive regulators, and catalysts

For detailed test documentation, see `TEST_SUITE_SUMMARY.md`.

## Development

### Run Type Checking

```bash
poetry run mypy --ignore-missing-imports .
```

### Run Linting

```bash
poetry run flake8 .
```

## Documentation

### Architecture
- **[Architecture Overview](docs/ARCHITECTURE.md)** - Complete system architecture, data flow, and key concepts
  - Data flow from Neo4j to logic network
  - Virtual reactions and edge semantics
  - AND/OR logic rules
  - Design decisions and rationale

### Test Documentation
- **[Test Suite Summary](TEST_SUITE_SUMMARY.md)** - Overview of all 52 tests
- **[Test Findings](TEST_FINDINGS.md)** - Investigation results from edge direction analysis
- **[Complete Understanding](COMPLETE_UNDERSTANDING.md)** - Definitive explanation of edge semantics

### Improvement Documentation
- **[Improvement Recommendations](IMPROVEMENT_RECOMMENDATIONS.md)** - Prioritized list of 15 improvements
- **[Quick Wins](QUICK_WINS.md)** - 9 quick improvements (~2 hours total)
- **[Changelog](CHANGELOG.md)** - Detailed history of all changes
