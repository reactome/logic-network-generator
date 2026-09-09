# Tasks: `break_apart_entity` memo must agree with the function

**Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md) | **Issue**: #58

## Phase 1 — Establish the contract

- [x] **T001** Read both writers and both return statements; determine which
      column each entity kind's return value lives in.
- [x] **T002** Confirm the memo is correct for EntitySets and wrong only for
      Complexes, so the fix does not break the working half.
- [x] **T003** Confirm no id can collide between the two writers
      (`reactome_id` is the set id for provenance rows and the complex id for
      iterproduct rows; reaction rows use the reaction id).

## Phase 2 — Fix

- [x] **T004** Branch the memo on `"EntitySet" not in labels`; return
      `{r["uid"]}` for a Complex, the `input_or_output_*` union otherwise.
- [x] **T005** Verify the `MAX_VARIANTS` cap path still round-trips (FR4): the
      cap emits one bundled uid and the memo returns that uid.

## Phase 3 — Verify

- [x] **T006** Idempotence probe on R-HSA-69242: 11 of 48 entities disagreed
      before, 0 after. `R-HSA-187501` now returns the same 2 combination uids
      on both calls.
- [x] **T007** Lockdown test in `tests/test_decomposition_semantics.py`,
      extending `TestCrossCallStability` — which asserted stability for
      EntitySets only, the case the memo got right, which is why the bug
      survived. Confirmed to FAIL with the fix reverted.
- [x] **T008** Full suite: 924 passed vs 923 on main; the same 7 pre-existing
      failures (stale checked-in `output/` tree), no new ones. `ruff` clean.
- [x] **T009** Validator: 11/11 on all three benchmark-eligible pathways.
- [x] **T010** DeltaSignal A/B, ten evaluation pathways, experimental ground
      truth. NEUTRAL: 452/627, macro-F1 0.664946 in both arms, 0 of 847
      predictions changed.
- [x] **T011** Confirm the A/B was a real comparison: only 3 of 43 network
      hashes shared between arms, and 221 raw values / 389 gene-uuid counts /
      22 convergence flags differ. A neutral result from an accidentally
      identical input would look the same in the summary.

## Phase 4 — Open

- [ ] **T012** Re-run the catalog-wide regeneration with the memo fixed and
      re-measure #61: bridges should shrink by roughly the VR-inflation factor
      before any cofactor guard is added, which may change how that fix is
      scoped.
- [ ] **T013** Re-check the `MAX_VARIANTS` cap hits catalog-wide. Several were
      triggered by inflated combination counts and should now disappear; a cap
      hit silently bundles components into one opaque node, so fewer is a
      correctness win worth quantifying.
- [ ] **T014** Consider whether `get_reaction_connections` should carry an
      explicit `ORDER BY`. The memo bug made decomposition order-dependent
      (reversing the reaction list changed R-HSA-69242 from 134 to 116 VRs);
      that is fixed, but row order is still relied upon implicitly.
