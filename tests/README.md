# Tests

Three tiers, distinguished by what they need to run.

## Tiers

| Tier | Marker | Needs | Runs in CI? |
|---|---|---|---|
| **Unit** | _(none)_ | Just the source tree | Yes |
| **Integration** | `integration` | Generated pathway artifacts in `output/` from a prior local run | No |
| **Database** | `database` | A live Neo4j on `bolt://localhost:7687` (auth `neo4j` / `test`) | No |

The CI workflow runs `pytest -m "not database and not integration"` — i.e., the unit tier only.

## Running locally

```bash
# Unit only (what CI runs)
poetry run pytest -m "not database and not integration"

# Integration tier (requires output/ artifacts from a prior pathway generation)
poetry run pytest -m integration

# Database tier (requires Neo4j running locally)
poetry run pytest -m database

# Everything
poetry run pytest
```

To populate `output/` for the integration tier, run a pathway generation first:

```bash
poetry run python bin/create-pathways.py --pathway-id 69620
```

## Adding a new test

- **No DB, no `output/`** — write it as a normal unit test, no marker needed.
- **Reads files from `output/`** — add `pytest.mark.integration` (module-level `pytestmark`).
- **Hits Neo4j** — add `pytest.mark.database` (class- or module-level).

If a test needs both, mark it `database` — the assumption is that you generated `output/` from the same Neo4j you're testing against.
