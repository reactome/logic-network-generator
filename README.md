# Logic Network Generator

[![Tests](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml/badge.svg)](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml)
[![Code Style](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Generate logic networks from Reactome pathways by decomposing sets and complexes into their individual components.

## Features

- ✅ **Position-Aware UUIDs** - Same entity at different positions gets unique identifiers
- ✅ **Comprehensive Validation** - 100% validated against source database
- ✅ **Identifier Resolution** - Find entities by UniProt, gene symbol, or Reactome ID
- ✅ **Batch Processing** - Generate multiple pathways from a list
- ✅ **Production Ready** - Full test coverage, error handling, and logging

## Quick Start

### Prerequisites

- [Python 3.9+](https://www.python.org/downloads/)
- [Poetry](https://python-poetry.org/)
- [Docker](https://www.docker.com/) (for Neo4j database)

### Installation

```bash
# Clone and install
git clone https://github.com/reactome/logic-network-generator.git
cd logic-network-generator
poetry install

# Start Neo4j Reactome database (easiest method)
docker-compose up -d

# Or using plain docker
docker run -p 7474:7474 -p 7687:7687 \
  -e NEO4J_dbms_memory_heap_maxSize=8g \
  public.ecr.aws/reactome/graphdb:Release94
```

### Generate a Pathway

```bash
# Single pathway
poetry run python bin/create-pathways.py --pathway-id 69620

# Multiple pathways
poetry run python bin/create-pathways.py --pathway-list pathways.tsv
```

## Output Files

All generated files are saved to the `output/` directory:

- **`pathway_logic_network_{id}.csv`** - Main logic network with edges
- **`uuid_mapping_{id}.csv`** - UUID to Reactome ID mapping with position info
- **`decomposed_uid_mapping_{id}.csv`** - Complex/set decomposition details
- **`reaction_connections_{id}.csv`** - Reaction connectivity graph
- **`best_matches_{id}.csv`** - Input/output matching for reactions

## Logic Network Format

The generated logic network CSV has these columns:

| Column | Description |
|--------|-------------|
| `source_id` | UUID of source entity |
| `target_id` | UUID of target entity |
| `pos_neg` | `pos` (activation) or `neg` (inhibition) |
| `and_or` | `and` (all inputs required) or `or` (any input sufficient) |
| `edge_type` | `input`, `output`, `catalyst`, or `regulator` |

## Utilities

### Create Database ID Mapping

Generate a mapping file from Reactome database IDs to human-readable names:

```bash
# Basic usage (human entities only)
poetry run python bin/create-db-id-name-mapping-file.py

# All species
poetry run python bin/create-db-id-name-mapping-file.py --all-species

# Custom output location
poetry run python bin/create-db-id-name-mapping-file.py --output my_mapping.tsv
```

Output columns: `database_identifier`, `node_type`, `display_name`, `reference_entity_name`, `reference_entity_identifier`, `instance_class`

## Validation

Comprehensive validation ensures generated networks match the source database:

```bash
# Run all validation tests
poetry run pytest tests/test_pathway_validation.py -v

# Run comprehensive validation (includes loop analysis, regulator matching, identifier resolution)
poetry run pytest tests/test_comprehensive_validation.py -v

# Quick validation script
poetry run python validate_pathway.py 69620
```

See [VALIDATION_README.md](VALIDATION_README.md) for details.

## Testing

```bash
# Run unit tests (no database required - fast)
poetry run pytest tests/ -v -m "not database"

# Run all tests including database tests (requires Neo4j)
poetry run pytest tests/ -v

# Run only database/integration tests
poetry run pytest tests/ -v -m "database"

# Run with coverage
poetry run pytest tests/ --cov=src --cov-report=html -m "not database"
open htmlcov/index.html

# Run specific test categories
poetry run pytest tests/test_and_or_logic.py -v
poetry run pytest tests/test_regulators_and_catalysts.py -v
poetry run pytest tests/test_network_invariants.py -v
```

**Test Suite**: 82 tests total
- **62 unit tests** - Core functionality, AND/OR logic, regulators, invariants (no database required)
- **20 integration tests** - Comprehensive validation against Neo4j database (requires database)

## Examples

Complete working examples in the `examples/` directory:

```bash
poetry run python examples/generate_pathway_example.py
```

See [examples/README.md](examples/README.md) for more usage patterns and example pathways.

## Documentation

- **[Architecture](docs/ARCHITECTURE.md)** - System architecture, data flow, and design decisions
- **[Position-Aware UUIDs](POSITION_AWARE_UUID_DESIGN.md)** - Design and implementation of position-aware UUID system
- **[Validation](VALIDATION_README.md)** - Comprehensive validation system documentation
- **[Examples](examples/README.md)** - Usage examples and patterns
- **[Changelog](CHANGELOG.md)** - Version history and notable changes

## Development

```bash
# Start Neo4j database
docker-compose up -d

# Stop Neo4j database
docker-compose down

# Type checking
poetry run mypy --ignore-missing-imports src/

# Linting
poetry run ruff check src/

# Formatting
poetry run ruff format src/

# Pre-commit hooks
poetry run pre-commit install
poetry run pre-commit run --all-files
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed development guidelines.

## License

Apache 2.0 - See [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
Logic Network Generator - Reactome Pathway Logic Network Generation Tool
https://github.com/reactome/logic-network-generator
```
