# Feature Specification: `break_apart_entity` memo must agree with the function

**Feature Branch**: `fix/break-apart-memo`
**Created**: 2026-09-09
**Status**: Measured — see Results
**Issue**: [#58](https://github.com/reactome/logic-network-generator/issues/58)

## Problem

`break_apart_entity` memoizes its result by re-reading the rows it wrote to
`decomposed_uid_mapping`. For a Complex, the memo returns **different values
from the function it memoizes**, and the caller cannot tell.

The two entity kinds return different things, and the memo only handles one:

| kind | uncached branch returns | rows written by | memo read |
|---|---|---|---|
| `EntitySet` | its **members** (alternatives) | `_emit_entityset_provenance_rows` — members land in `input_or_output_uid` / `input_or_output_reactome_id` | union of those two columns — **correct** |
| `Complex` | `get_broken_apart_ids(...)` = the **combination uids**, one per set-variant, each standing for the whole complex | `get_uids_for_iterproduct_components` — combination in `uid`, individual **components** in `input_or_output_*` | union of `input_or_output_*` — **the components, not the combinations** |

The caller feeds the result into a per-reaction `itertools.product` as
*alternatives*. So on its second and every later appearance, a Complex is
modelled as "any one of my subunits".

## Why this is a correctness bug, not a modelling choice

`R-HSA-187501 = CCNA:CDK2` is used by two reactions. `R-HSA-174273` is
processed first and correctly yields 2 variants; `R-HSA-174054` — the reaction
that *forms* the complex — hits the memo. From the generated `S_Phase` cache:

```
6af74b4209 -> component_id R-HSA-157440   (CDK2 alone)
e18b12b9e3 -> component_id R-HSA-157444   (CCNA2 alone)
3fda7d321c -> component_id R-HSA-157446   (CCNA1 alone)
```

`CCNA + CDK2 → CCNA:CDK2` is emitted with **"CDK2 alone" as a complete
output**. The module docstring forbids exactly this ("Complex = a bound
species… treated atomically"). No curator asserts that a free subunit is a
product of complex formation, so this is not a representation choice the
constitution's "represent curator intent" principle protects — it is the
generator stating something false.

## Requirements

- **FR1** — The memo MUST return what the uncached branch returns, for both
  entity kinds. Complexes return combination uids; EntitySets return members.
- **FR2** — Re-decomposing any entity within a pathway MUST be idempotent.
- **FR3** — The fix MUST be validated end-to-end through DeltaSignal, not by
  edge counts. Connectivity metrics have repeatedly been the wrong metric in
  this project.
- **FR4** — Behaviour under `MAX_VARIANTS` truncation must be preserved: the
  cap path returns a single bundled uid, and the memo must return that uid.

## Measurements

### Idempotence (FR2), R-HSA-69242

| | entities disagreeing on re-decomposition |
|---|---|
| before | 11 of 48 |
| after | **0** |

`R-HSA-187501` before: 2 flat leaf stIds. After: the same 2 combination uids
the uncached call returns.

### Network size, benchmark-eligible pathways (Release97, `PYTHONHASHSEED=0`)

| pathway | edges before → after | VRs before → after |
|---|---|---|
| R-HSA-1257604 PIP3 | 4,595 → 4,560 | 603 → 593 |
| R-HSA-453279 Mitotic G1 | 3,419 → 3,108 | 580 → 493 |
| R-HSA-69620 Cell Cycle Checkpoints | 3,836 → 3,750 | 798 → 784 |

Catalog-wide the effect is far larger on the worst pathways (NER 2,175 → 205
VRs; catalog total 53,782 → 39,848), which is why this is the root cause of the
`diagram_bridge` blow-up (#61) and of spurious `MAX_VARIANTS` truncation.

### DeltaSignal A/B — see `plan.md` for the protocol

Both arms: DeltaSignal `main` (unchanged), Reactome Release97, the same ten
MP-BioPath evaluation pathways, `--ground-truth experimental`, harness
defaults. The only difference between arms is the memo.

| arm | scored | correct | accuracy | macro-F1 | balanced acc |
|---|---|---|---|---|---|
| control (memo as-is) | 627 | 452 | 0.7209 | 0.664946 | 0.6959 |
| memo agrees | 627 | 452 | 0.7209 | 0.664946 | 0.6959 |

**NEUTRAL, and verified to be a real comparison rather than a no-op.** The two
catalogs genuinely differ — only 3 of 43 recorded network file hashes are
shared — and the networks change substantially:

| quantity | cases differing (of 847) |
|---|---|
| **prediction (the classification)** | **0** |
| `predicted_ui` (the raw value) | 221 |
| `gene_uuid_count` | 389 |
| `output_uuid_count` | 42 |
| `converged` | 22 |
| `iterations` | 323 |

So the fix moves a quarter of the raw activity values and reduces the number of
UUIDs a gene perturbation pins in nearly half the cases, without flipping a
single classification. Convergence improves slightly (non-converged 155 → 149;
14 cases gain convergence, 8 lose it) at a negligible cost of +0.8% total solver
iterations.

### Interpretation

This is the outcome the architecture principle predicts for a pure correctness
fix. The generator stops asserting that a free subunit is a complete product of
complex formation; DeltaSignal's classifications are unchanged, so nothing in
the solver was depending on the false structure. Unlike `LNG_SET_MEMBERS_OR`
— a faithful change that cost 51 of 223 cases and had to stay default-OFF
pending a solver capability — this one has no benchmark cost and needs no flag.

The benchmark cannot detect the improvement, because the duplicated VRs were
byte-identical copies computing the same value: removing redundancy changes
magnitudes and iteration counts but not answers. The case for landing it is
correctness plus the downstream consequences (#61 bridge scale, spurious
`MAX_VARIANTS` truncation, catalyst-edge multiplication), not accuracy.

Control differs slightly from the committed GSoC evaluation (452 vs 456 of 627)
because both arms here use a fresh Release97 regeneration on current LNG `main`,
rather than that report's frozen catalog pinned at LNG `7aca90d`. Control and
treatment differ only in the memo, which is what the comparison requires.

## Out of scope

- The `diagram_bridge` all-pairs/cofactor work (#61), which this shrinks but
  does not fix.
- Dropped zero-input/zero-output reactions (#59).
- AND-flattened assembly edges (#60).
