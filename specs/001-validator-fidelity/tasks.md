# Tasks: Validator fidelity to the generated representation

**Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md)

## Phase 1 — Make it run (FR1, FR2)

- [x] **T001** `_resolve_layout()`: try `output/*_R-HSA-{id}` and `output/*_{id}`,
      fall back to the legacy flat `pathway_logic_network_{dbId}.csv`.
- [x] **T002** `load_files()`: normalise `stable_id` to `reactome_id` +
      `entity_ids`; treat `cache/decomposed_uid_mapping.csv` as optional with a
      warning rather than aborting.
- [x] **T003** Connect via `NEO4J_URL` (fallback `NEO4J_URI`) +
      `NEO4J_USER`/`NEO4J_PASSWORD`. Was hardcoded `auth=("neo4j","test")`, so
      pointing it at a secured database silently validated localhost.
- [x] **T004** `_mapping_uses_stid()` / `_parse_entity_ids()`; select `stId` vs
      `dbId` per check. This alone moved reconstruction from a spurious 0.0% to
      a real 88.2% on R-HSA-69620.
- [x] **T005** Widen `valid_edge_types` to `{input, output, catalyst,
      regulator, assembly, dissociation, depletion}` — the original four failed
      on healthy networks.
- [x] **T006** Parameterise all Cypher on `$pathway_id` (FR6). Previously
      interpolated, safe only because `argparse type=int` coerced it.

## Phase 2 — Make it right (FR3, FR4, FR5)

- [x] **T007** `_set_leaf_members()` + `_uncovered()`; wire into entity coverage
      and catalyst completeness. Verified first: all 9 "missing entities" on
      R-HSA-69620 are `EntitySet`s with 28/28 members present.
- [x] **T008** Report partially-split sets as a distinct FAIL naming the absent
      members (FR4) — do not let T007 mask a half-dropped set.
- [x] **T009** `_decomposed_ids()`: count a decomposed Complex as covered.
      Verified first: the 14 "partially split" members are all present in
      `decomposed_uid_mapping.csv`.
- [x] **T010** Rewrite `validate_regulator_polarity()` per
      `(regulator, reaction)` pair, expanding set-valued regulators to leaves.
      Verified first: PI5P really is curated both ways in R-HSA-1257604.
- [x] **T011** Count regulator edges with no matching Neo4j regulation
      separately as a warning, not as a polarity failure.

## Phase 3 — Make it trustworthy (FR7)

- [x] **T012** Report comparison counts on every check (`Checked 271 regulator
      edges`, `N EntitySets represented by their members`).
- [x] **T013** Add `--output-dir` so a corrupted copy can be validated.
- [x] **T014** Corruption control: 5 flipped `pos_neg` values produce exactly 5
      failures with correct attribution.
- [x] **T015** Harden the `LNG_VALIDATE_MIN_RECONSTRUCTION` parse — degrade to
      the default with a warning instead of raising.
- [x] **T016** `ruff` clean (also removed 3 pre-existing unused `typing`
      imports).

## Phase 4 — Open

- [ ] **T017** Investigate finding **F1**: `R-HSA-9980233` ("PIP3 activates
      mTORC2") is absent from every generated artifact for R-HSA-1257604
      (88/89 reactions represented). Not a self-loop, not disease, not
      inferred, output complex decomposes, `precedingEvent` present both ways.
      Determine the drop point in the generator. The pathway has no diagram of
      its own, so diagram-driven filtering is a candidate — unconfirmed.
      Per the architecture principle, any fix needs a DeltaSignal before/after,
      not just an edge-count delta.
- [ ] **T018** Triage the reconstruction gap on R-HSA-69620: 16 Neo4j
      input→output pairs absent, 53 extra pairs emitted. Currently warned at
      88.2%; decide whether the floor should rise once explained.
- [ ] **T019** Run the validator across the wider curator pathway set, not just
      the 3 benchmark-eligible ones, and record the F1-class rate.
- [ ] **T020** Wire the validator into CI once T017/T018 are resolved and the
      three pathways are green.

## Notes

The pattern across this whole review — retracted "122 non-converged", the
dead-code cache fingerprint, and now three classes of validator false
failure — is that a check or a claim must be verified against ground truth
before it is reported. Every relaxation in Phase 2 was confirmed with a Neo4j
query first; F1 is reported as real only because the same queries ruled out
every structural explanation tried.
