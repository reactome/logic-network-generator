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
- [x] **T005** Widen `valid_edge_types` to every type the generator emits —
      the original four failed on healthy networks. See also T021.
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

- [x] **T017** ~~Investigate finding F1~~ **RETRACTED — not a defect.**
      `R-HSA-9980233` is a Reactome **v97** addition; the checked-in output was
      generated 2026-07-16 against v96, before the v97 bump, and carries no
      cache fingerprint. Traced the generator stage by stage: connections
      (89/89), decomposition (1 input + 1 output combination) and
      `find_best_reaction_match` (emits a VR) all handle it correctly.
      Regenerated against Release97 -> present in every artifact, coverage
      passes. Validating a stale artifact against a live database is version
      skew, not a generator gap.
- [x] **T021** Add `diagram_bridge` and `handoff` to the edge-type allowlist —
      surfaced by the fresh regeneration; the stale July artifact predates the
      additive diagram-connectivity work so it never exercised them.
- [x] **T022** `_cache_is_reusable`: return False for a cache directory with no
      CSVs. It was warning "reusing it, provenance UNVERIFIED" about an empty
      directory and stamping `provenance="adopted"` onto the cache the run was
      about to *generate*, mislabelling fresh output as unverified.
- [ ] **T018** Triage the reconstruction gap on R-HSA-69620: 16 Neo4j
      input→output pairs absent, 53 extra pairs emitted. Currently warned at
      88.2%; decide whether the floor should rise once explained.
- [ ] **T019** Run the validator across the wider curator pathway set, not just
      the 3 benchmark-eligible ones.
- [ ] **T023** Regenerate the checked-in `output/` tree against Release97. It
      is v96-era and unfingerprinted; anything validated or benchmarked against
      it inherits the skew that produced the retracted F1.
- [ ] **T020** Wire the validator into CI once T018 is resolved. All three
      pathways are green (11/11) on Release97-fresh output; CI must regenerate
      rather than validate checked-in artifacts, or it will re-run T017.

## Notes

The pattern across this whole review — the retracted "122 non-converged"
figure, the dead-code cache fingerprint, three classes of validator false
failure, and now the retracted F1 — is that a claim must be verified against
ground truth before it is reported.

F1 is the sharpest version: every *structural* explanation was correctly ruled
out with Neo4j queries, and it was still wrong, because the thing never
questioned was the **provenance of the artifact being validated**. Ruling out
explanations is not the same as establishing a cause. When a check disagrees
with an artifact, the artifact's origin — release, code version, cache
fingerprint — is part of the evidence, not background.
