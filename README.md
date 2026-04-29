# Logic Network Generator

[![Tests](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml/badge.svg)](https://github.com/reactome/logic-network-generator/actions/workflows/test.yml)
[![Code Style](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Generate logic networks from Reactome pathways by decomposing complexes and EntitySets into their components, expanding alternatives, and emitting a graph of input/output/catalyst/regulator/assembly/dissociation edges suitable for perturbation modeling.

## Features

- **Position-aware UUIDs** — the same entity at different pathway positions gets distinct identifiers, so perturbing a protein in one location doesn't unintentionally perturb it elsewhere.
- **Full EntitySet expansion with provenance** — every alternative input becomes its own virtual reaction; `decomposed_uid_mapping.csv` records `source_entity_id` so leaves can be traced back to their parent set.
- **Boundary decomposition** — root-input and terminal-output complexes get synthetic assembly / dissociation edges so individual proteins are perturbable at the network's edges; intermediate complexes stay intact (they're real species in the pathway).
- **Reaction-level AND/OR semantics** — input/catalyst/positive-regulator edges are AND, negative regulators are OR (any one blocks). See `docs/DESIGN_DECISIONS.md`.
- **Bulk Cypher pre-fetch** — pathway generation pulls all entity/reaction data in five queries up front, making per-reaction processing cache-only. Cell_Cycle's 474 reactions go from hours to minutes.
- **Validated against curator predictions** — 70.55% end-to-end agreement with the MP-BioPath curator-prediction test set across 12,895 valid cases; 98.3% on cases where the network is the deciding factor. See `validation_results/`.

## Quick Start

### Prerequisites

- [Python 3.10+](https://www.python.org/downloads/)
- [Poetry](https://python-poetry.org/)
- [Docker](https://www.docker.com/) (for the Reactome Neo4j database)

### Installation

```bash
git clone https://github.com/reactome/logic-network-generator.git
cd logic-network-generator
poetry install

# Start a local Reactome Neo4j database
docker-compose up -d
```

By default the connection points at `bolt://localhost:7687` with user `neo4j` / password `test`. Override via env vars `NEO4J_URL`, `NEO4J_USER`, `NEO4J_PASSWORD`.

### Generate a Pathway

```bash
# Single pathway (use the R-HSA-prefixed stable ID)
poetry run python bin/create-pathways.py --pathway-id R-HSA-69620

# Batch from a TSV with `id` and `pathway_name` columns
poetry run python bin/create-pathways.py --pathway-list pathways.tsv

# Every Homo sapiens top-level pathway
poetry run python bin/create-pathways.py --top-level-pathways
```

## Output

Each pathway generates a directory under `output/`:

```
output/<Pathway_Name>_R-HSA-<id>/
├── logic_network.csv          # Main output: edges of the perturbation graph
├── stid_to_uuid_mapping.csv   # UUID → Reactome stable ID
└── cache/
    ├── reaction_connections.csv
    ├── decomposed_uid_mapping.csv
    └── best_matches.csv
```

## Logic Network Format

`logic_network.csv` columns:

| Column | Description |
|---|---|
| `source_id` | UUID of source node (entity or virtual reaction) |
| `target_id` | UUID of target node |
| `pos_neg` | `pos` (activates / produces) or `neg` (negative regulator) |
| `and_or` | `and` (required), `or` (alternative source), or empty (single producer) |
| `edge_type` | `input`, `output`, `catalyst`, `regulator`, `assembly`, or `dissociation` |
| `stoichiometry` | Stoichiometric coefficient from Reactome |

`assembly` and `dissociation` edges only appear at boundaries: a leaf protein assembles into a root-input complex, or a terminal-output complex dissociates into its components.

## Validation

The generated networks have been benchmarked against the MP-BioPath curator-prediction test set ([Sundararaman et al., 2017](https://reactome.org/community/publications)).

- **70.55% end-to-end accuracy** on 12,895 valid test cases (vs ~75% published for MP-BioPath on its 10-pathway empirical subset)
- **98.3% network correctness** on cases where the network's connectivity is the deciding factor (excludes propagator limitations and v86→v96 test-set drift)
- 1.2% of valid tests flagged as `bug_candidate`; on per-case investigation all examined cases were structural limitations of directed-flow Boolean propagation (substrate-consumption mass-action effects), not network-generation bugs

Full methodology, per-pathway reports, and failure-category analysis: `validation_results/README.md`.

## Utilities

```bash
# Generate a database-ID-to-name mapping file
poetry run python bin/create-db-id-name-mapping-file.py

# Validate a generated set of networks against the MP-BioPath test set
poetry run python bin/validate-against-mpbiopath.py

# For each "no_path" failure, ask Neo4j whether the original pathway has a path
poetry run python bin/check-no-path-cases-in-neo4j.py
```

## Testing

Three tiers, distinguished by what they need to run:

```bash
# Unit tier (CI runs this — no database, no generated artifacts)
poetry run pytest -m "not database and not integration"

# Integration tier (needs output/ artifacts from a prior run)
poetry run pytest -m integration

# Database tier (needs a running Reactome Neo4j)
poetry run pytest -m database
```

See `tests/README.md` for details on each tier.

### Tracking new Reactome releases

The database-tier tests are version-agnostic — they discover pathways
from the loaded graph rather than hard-coding stable IDs — so they
should survive a Reactome bump. The workflow when a new Reactome
release ships:

1. Update the image tag in `docker-compose.yml`
   (`public.ecr.aws/reactome/graphdb:Release<N>`).
2. `docker-compose pull && docker-compose up -d`.
3. `poetry run pytest -m database -v` — `test_reactome_version.py`
   records which release the run actually used; other tests will
   surface anything Reactome changed schematically.
4. If accuracy numbers in `validation_results/` matter for that
   release, re-run `poetry run python bin/validate-against-mpbiopath.py`
   and update the README headline numbers.

## Documentation

- [Architecture](docs/ARCHITECTURE.md) — system architecture and data flow
- [Design Decisions](docs/DESIGN_DECISIONS.md) — behaviors that look surprising but are intentional (Complex vs EntitySet semantics, the two-layer decomposition model, surplus input/output fan-out)
- [Position-Aware UUIDs](docs/UUID_DESIGN.md) — why and how UUIDs are assigned per pathway position
- [Examples](examples/README.md) — usage examples and patterns
- [Validation Results](validation_results/README.md) — full benchmark methodology and per-pathway numbers

## Development

```bash
# Start the Neo4j database
docker-compose up -d

# Lint and format
poetry run ruff check src/
poetry run ruff format src/

# Pre-commit hooks
poetry run pre-commit install
poetry run pre-commit run --all-files
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

## License

Apache 2.0 — see [LICENSE](LICENSE).
