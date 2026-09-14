# CLAUDE.md

Guidance for Claude Code working in this repository.

This repo is the **upstream** of DeltaSignal: it reads Reactome from Neo4j and
emits the logic-network TSV/CSV that `~/gitroot/deltasignal` consumes. The
division of labour is deliberate — **the generator represents pathways the way
curators authored them; DeltaSignal decides how to process that
representation.** Do not push modelling choices upstream into the generator.

## Where design decisions live

This repo uses spec-kit. **`specs/NNN-name/` is the record of why** — what was
measured, what was rejected, and the numbers. This file carries only current
state; per-feature results belong in that feature's `research.md`, not here.
Before changing generation behaviour, read the spec for the last feature that
touched it rather than re-deriving from the code.

- `specs/001-validator-fidelity/` — the validator must be able to fail.
- `specs/002-break-apart-memo/` — the memo returned components where the
  function returned combinations, so a repeated Complex degraded into "any
  one of my subunits".
- `specs/003-bridge-cofactor-guard/` — cofactor guard on diagram bridges.
- `.specify/memory/constitution.md` — principles the specs are checked
  against.
- `docs/DESIGN_DECISIONS.md`, `docs/UUID_DESIGN.md`,
  `docs/ARCHITECTURE.md` — durable representational semantics (Complex vs
  EntitySet, the two decomposition layers, uuid scheme). These are current
  state, not history, and stay where they are.
- Behaviour that is not yet a spec is tracked as GitHub issues on
  `reactome/logic-network-generator`.

## Running

Configuration is `.env` (see `.env.example`); `NEO4J_URL`, `NEO4J_USER`,
`NEO4J_PASSWORD` reach the connector via `load_dotenv`.

```bash
poetry run python bin/create-pathways.py --pathway-id R-HSA-69620
poetry run python bin/create-pathways.py --pathway-list pathways.tsv
poetry run python bin/create-pathways.py --top-level-pathways
```

Tests are tiered; the default tier needs no database:

```bash
poetry run pytest -m "not database and not integration"
poetry run pytest -m integration   # needs output/ from a prior run
poetry run pytest -m database      # needs a running Reactome Neo4j
```

See `README.md` for the full CLI and `PIPELINE_SETTINGS.md` for the
`LNG_*` switches.

## Traps that have cost real time

- **Regenerating into a populated directory can be a silent no-op.**
  Pre-fingerprint caches have no fingerprint, so `_cache_is_reusable` adopts
  them and reproduces the old networks while reporting success. `rm -rf
  <output-dir>/*/cache` before any regeneration you intend to be real. This
  has invalidated an A/B before.
- **`PYTHONHASHSEED` must be set before the interpreter starts** or output is
  not reproducible; `bin/create-pathways.py` re-execs itself to enforce it.
  `LNG_ALLOW_NONDETERMINISM=1` opts out.
- **`output/` is gitignored.** Tests that read it skip in CI, so a local
  failure there is not necessarily a CI failure — and a CI pass is not
  evidence those tests ran.
- **Version skew is the first thing to check, not the last.** More than one
  "the generator drops X" investigation has ended in the artifact being from
  a different Reactome release than the database being queried.
- Credentials can reach logs through py2neo's own exception reprs;
  `src/credential_redaction.py` scrubs at the output boundary. Do not
  reintroduce per-call-site redaction — it was proven insufficient.

## Validating changes against DeltaSignal

Generation changes are only meaningful if they move accuracy, and the
benchmark is noisy in a specific way: uuid4 node ids are minted fresh per
build, dict order follows them, and that sets solve order inside a cyclic
component — so **a non-converged solve returns different values for a
structurally identical network.** Any A/B must (1) build both arms from the
same catalog directory or report the network-hash overlap, (2) give a
per-pathway net breakdown, and (3) count changed predictions where both arms
converged. See `~/gitroot/deltasignal/specs/002-upregulation-propagation/quickstart.md`.
