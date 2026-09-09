# Implementation Plan: `break_apart_entity` memo must agree with the function

**Spec**: [spec.md](./spec.md) | **Branch**: `fix/break-apart-memo` | **Issue**: #58

## Approach

The fix is three lines of dispatch, but the risk is entirely in getting the
*semantics* right, so the work is mostly establishing what each branch owes its
caller before changing anything.

1. **Establish what the uncached branches actually return.** Read the two
   writers (`get_uids_for_iterproduct_components`, `_emit_entityset_provenance_rows`)
   and the two return statements, and confirm which column holds what. The memo
   is correct for EntitySets and wrong for Complexes; a blanket rewrite would
   have broken the working half.
2. **Make the memo branch on label**, returning `uid` for a Complex and the
   `input_or_output_*` union for an EntitySet.
3. **Prove idempotence empirically**, not by reading — re-decompose every
   entity in a pathway and count disagreements.
4. **Measure through DeltaSignal** (FR3), with a control that differs only in
   the memo.

## Key decisions

**Branch on the absence of `EntitySet`, not the presence of `Complex`.** The
guard is `if "EntitySet" not in labels`. Reactome labels are a set and the
outer condition already established that one of the two is present, so this
routes anything set-like down the member path and only a pure Complex down the
combination path. Getting this backwards would silently reintroduce the bug for
any entity carrying both labels.

**Do not "fix" the row schema.** It would be tidier for the Complex writer to
record its combination in a column named for it, but changing the emitted
`decomposed_uid_mapping` schema would ripple into the validator, the benchmark
harness, and every cached artifact. The memo is the defect; the schema is
merely awkward.

**No feature flag.** Project convention is default-OFF until measured, but that
convention exists for *modelling choices* whose direction is unknown. This is a
function disagreeing with its own memo, and the A/B shows zero prediction
changes, so there is nothing to hedge. `LNG_SET_MEMBERS_OR` is the contrasting
case: faithful, directionally uncertain, measurably costly, and correctly
flagged.

## Verification

- Unit: the new lockdown test must FAIL with the fix reverted (confirmed) and
  pass with it.
- Idempotence: 0 entities disagreeing on re-decomposition across R-HSA-69242.
- Suite: no new failures. 7 pre-existing failures on `main` (they validate
  against the stale checked-in v96 `output/` tree) remain exactly 7.
- Validator: all three benchmark-eligible pathways 11/11 on the regenerated
  networks.
- DeltaSignal A/B: both arms on DeltaSignal `main`, Release97, ten evaluation
  pathways, harness defaults, differing only in the memo. Verify the catalogs
  genuinely differ before believing a neutral result.
