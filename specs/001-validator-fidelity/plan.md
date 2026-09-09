# Implementation Plan: Validator fidelity to the generated representation

**Spec**: [spec.md](./spec.md) | **Branch**: `fix/adversarial-review`

## Approach

Three layers, applied in order, because each one only becomes observable once
the previous is fixed:

1. **Make it run at all.** Layout resolution + schema normalisation. Until
   `load_files()` succeeds, no check has ever executed, so no check has ever
   been validated.
2. **Make it right.** Teach each check the generator's representation, one
   structural difference at a time, verifying against Neo4j that each excused
   difference is genuinely structural before excusing it.
3. **Make it trustworthy.** Prove non-vacuity: report comparison counts, and
   corrupt a network on purpose to confirm the check still fails.

Step 3 is the one that is normally skipped and is the reason the original
script's `validate_reconstruction` could pass at 0%.

## Key decisions

**Excuse by structure, never by symptom.** Each relaxation is tied to a
specific Reactome construct and checked against Neo4j at validation time
(`_set_leaf_members` resolves real `hasMember`/`hasCandidate` edges; the
decomposed mapping is read from the generator's own output). No relaxation is
a blanket "ignore missing ids".

**Keep the sharp edge.** FR4 exists because FR3 is dangerous on its own: if
"a set counts as present when its members are present" were the whole rule, a
generator that dropped half a set would pass. `_uncovered` therefore returns
`(missing, partial)` and reports partial splits as failures with the specific
absent members named.

**Per-pair, not per-entity, for polarity.** `edge_reaction_id` is `NaN` on all
271 regulator edges (it is only populated for input/output edges), so the
reaction is recovered from the target virtual-reaction node's mapping, which
resolves to the reaction stId. This is exact, not approximate.

**Unattested ≠ wrong polarity.** A regulator edge with no matching Neo4j
regulation is counted and warned separately, not failed as a polarity error —
those are different defects with different causes.

## Structure

```
scripts/validate_logic_network.py
  MAX_SET_NESTING              module constant, bounds set traversal
  _resolve_layout()            per-pathway dir, else legacy flat
  load_files()                 schema normalisation; decomposed mapping optional
  _mapping_uses_stid()         stId vs dbId discrimination
  _parse_entity_ids()          pipe-delimited ids, typed per layout
  _set_leaf_members()          EntitySet -> non-set leaves (FR3)
  _decomposed_ids()            ids represented by their components (FR3)
  _uncovered()                 -> (missing, partial)                (FR3, FR4)
  _entity_for_uuid()           node uuid -> Reactome id
  validate_entity_coverage()   FR3, FR4
  validate_catalyst_completeness()  FR3, FR4
  validate_regulator_polarity()     FR5, per (regulator, reaction)
  validate_reconstruction()    id-space aware; graded floor
  main()                       --output-dir  (needed for FR7 controls)
```

## Verification

Run all three benchmark-eligible pathways against Release97; every excused
difference must have been individually confirmed structural via a Neo4j query
recorded in the spec's Measurements section; and the FR7 corruption control
must fail with correct attribution.
