# Pipeline Settings

How to configure the **logic-network-generator (LNG)** + **DeltaSignal** pipeline
when running it. This is aimed at someone standing the pipeline up, not tuning
the model.

## TL;DR

**Do not set the tuning flags.** The `LNG_*` (generator) and `DS_*` (solver)
defaults are the validated winning configuration — they are the correct values
baked into code, and overriding them is only for benchmarking experiments.

The only things you configure for a normal run are **environment settings**:
where Neo4j is, and which Reactome release you're on. Everything else is a
default.

## 1. Settings you MUST set (environment)

### Neo4j connection (LNG)
The generator reads the Reactome graph from Neo4j. Defaults are for a local dev
instance and are almost certainly not yours:

| Variable | Default | Set to |
|---|---|---|
| `NEO4J_URL` | `bolt://localhost:7687` | your Reactome Neo4j bolt URL |
| `NEO4J_USER` | `neo4j` | your user |
| `NEO4J_PASSWORD` | `test` | your password |

### Reactome release alignment (the #1 gotcha)
Everything must be the **same Reactome release**. Current canonical release is
**Release97**. LNG already defaults to 97, so this lines up automatically —
you just must not mix versions:

- The **Neo4j graph** must be the Release97 database.
- The **diagram set** on disk must be the Release97 diagrams (see
  `LNG_DIAGRAM_DIR` below). Diagram connectivity is on by default and reads
  these.
- Any **ground truth** you score against must also be Release97.

Mixing releases (e.g. a v96 diagram set against a v97 graph) does not error —
it silently produces networks that are unscoreable or subtly wrong. If in doubt,
verify all three are Release97 before trusting results.

### Diagram directory (LNG)
| Variable | Default | Set it if |
|---|---|---|
| `LNG_DIAGRAM_DIR` | `~/reactome-diagrams/97` | your Release97 diagram JSON lives elsewhere |

## 2. Settings you MIGHT set (situational)

| Variable | Repo | Default | Set it if |
|---|---|---|---|
| `DS_PATHWAY_CATALOG` | DeltaSignal | `/app/pathway_catalog` | using the `GET /api/pathways` catalog flow — point it at LNG's output directory |
| `DS_REACTOME_ENRICH` | DeltaSignal | `1` (on) | set `0` when running offline or when you don't need Reactome display names — skips reactome.org calls (degrades gracefully either way; this just avoids the latency) |
| `LNG_MAX_VARIANTS` | LNG | `512` | leave it; only change if a specific very large complex still blows up generation |

## 3. Tuning flags — leave at defaults

These encode the validated model/generator behavior. **Do not change them for a
production run** — they are listed here only so you know what you're getting and
what NOT to touch. Overriding is for A/B benchmarking only.

**Generator (LNG):**

| Variable | Default | Meaning |
|---|---|---|
| `LNG_COMPLEX_AS_NODE` | `1` | a Complex is one node (the key connectivity fix) |
| `LNG_SET_EXPAND` | `1` | expand EntitySets to members at matching granularity |
| `LNG_DIAGRAM_CONNECTIVITY` | `1` | union diagram-drawn reaction pairs into connectivity |
| `LNG_DIAGRAM_BRIDGE` | `1` | add (don't merge) diagram bridge edges |
| `LNG_CATALYST_BUNDLE` | `0` | legacy; subsumed by `LNG_COMPLEX_AS_NODE` |
| `LNG_HANDOFF_EDGES` | `0` | precedingEvent hand-off bridges (net-negative, off) |

**Solver (DeltaSignal):**

| Variable | Default | Meaning |
|---|---|---|
| `DS_INHIBITION_MODE` | `divide` | de-repression-capable inhibition |
| `DS_AND_MODE` | `hill_log` | continuous log-fold AND aggregation |
| `DS_OR_MODE` | `mean` | OR alternatives each carry weight |
| `DS_ASSEMBLY_LIMITING` | `1` | a complex can't exceed its scarcest subunit |
| `DS_HILL_LOG_ZMAX` | `10.0` | AND saturation point |
| `DS_SCC_SOLVE` | `1` | SCC-condensation solve (handles feedback loops) |

(There are more `DS_*` switches for other inhibition/AND modes, loop handling,
etc. — all default to the validated config. If you didn't set it, it's right.)

## 4. Running the solver

If you drive DeltaSignal over its HTTP API, two things matter:

- **Activity scale is asymmetric.** Send `observations` on the **0–100** scale
  (`0` = none, `1` = normal baseline, `100` = 100× normal). The returned
  `node_activities` are on the **0–1** internal scale — multiply by `100` to put
  them back on the display scale.
- **Parse once, solve many.** For many perturbations against the same pathway,
  call `POST /api/parse` once, keep the returned `network_id`, and pass it to
  each `POST /api/solve` instead of re-shipping the whole network. `network_id`
  is an in-memory, per-process cache, so re-parse if the server restarts.

Full request/response schemas: `deltasignal/docs/API.md`.
