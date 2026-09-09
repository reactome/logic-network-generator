# Feature Specification: cofactor guard on diagram bridges

**Feature Branch**: `fix/bridge-cofactor-guard` (on top of `fix/break-apart-memo`)
**Created**: 2026-09-09
**Status**: Measured — neutral, landing as a correctness fix
**Issue**: [#61](https://github.com/reactome/logic-network-generator/issues/61)

## Problem

The additive diagram-bridge pass carries the comment "Layout-filtered pairs
only (cofactor hubs already excluded by the caller)". Nothing excluded them.
`_node_leaves` (the handoff path) *does* subtract
`_COFACTOR_STIDS | _UBIQUITIN_STIDS`; the bridge path did not.

A bridge asserts "this producer feeds that consumer". For a shared cofactor or
free ubiquitin that is meaningless — every reaction in a pathway touches ATP —
so bridging on one couples reactions with no causal relationship. Curators
sometimes draw a single shared Ub glyph, and VR multiplicity turns two glyphs
into hundreds of edges.

## Requirement

- **FR1** — No `diagram_bridge` may be emitted on a species in
  `_COFACTOR_STIDS | _UBIQUITIN_STIDS`, matching `_node_leaves`.

## Measurements

Measured **on top of the memo fix** (#58 / PR #64), because bridge counts are
quadratic in virtual-reaction count and the memo fix already removes most of
them. Measuring against unfixed `main` would have attributed the memo fix's
effect to this change.

### Bridges removed, ten evaluation pathways

| pathway | bridges before guard | after |
|---|---|---|
| R-HSA-69242 S Phase | 77 | **12** |
| R-HSA-453279 Mitotic G1 | 105 | **41** |
| R-HSA-69620 Cell Cycle Checkpoints | 983 | 946 |
| R-HSA-1227986 ERBB2 | 214 | 213 |
| R-HSA-5673001 RAF/MAPK | 732 | 731 |
| (5 others) | unchanged | unchanged |
| **total** | **5,526** | **5,358** |

Exactly the 168 cofactor/ubiquitin bridges predicted, and nothing else.

For scale, the memo fix alone had already taken S Phase 158 → 77, Mitotic G1
174 → 105, Mitotic Prophase 3,567 → 12 and HDR 20,552 → 3,004.

### DeltaSignal A/B — NEUTRAL

Both arms on DeltaSignal `main`, Release97, ten evaluation pathways,
experimental ground truth, differing only in the guard.

| arm | scored | correct | accuracy | macro-F1 |
|---|---|---|---|---|
| memo fix (base) | 627 | 452 | 0.7209 | 0.664946 |
| + cofactor guard | 627 | 455 | 0.7257 | 0.668600 |

**The apparent +3 is NOT the guard.** All 12 changed predictions are in
R-HSA-3700989 (TP53) — a pathway with **zero** diagram bridges in both arms,
whose canonical structure is byte-identical between them. Every changed case is
a non-converged solve:

```
R-HSA-3700989   AKT1  dir=0  conv True->False   pred 0->1
R-HSA-3700989   AKT1  dir=2  conv False->False  pred 2->0
R-HSA-3700989   BRCA1 dir=0  conv False->False  pred 2->0
... (12 total, all TP53, all non-converged)
```

In the five pathways the guard actually changed, **no prediction moved**.
Case-level McNemar on the raw 7 gains / 4 losses gives p = 0.55; clustered by
perturbation solve it is p = 1.00 with 0 net gains.

## Methodological finding (affects every LNG A/B in this project)

`uuid4` node ids are minted fresh on every regeneration — 0 of TP53's 2,307
uuids are shared between the two arms. Dict iteration order follows those ids,
which sets Gauss-Seidel sweep order inside an SCC, so **non-converged solves
return different values for a structurally identical network**. This matches
the independent finding that edge-order deviations appear only in solves that
honestly report `converged=false`.

Consequence: any LNG change measured by regenerating carries this noise, and it
can manufacture a flattering delta. An A/B here must either restrict to
converged cases or hold uuids stable, and a per-pathway breakdown is the cheap
way to catch it — a gain in a pathway the change did not touch is the tell.

This compounds the related problem in the DeltaSignal evaluation
(reactome/deltasignal#14): half the diagram-on/off "gains" there also come from
pairs with a non-converged arm.

## Decision

Land it. It is prediction-neutral, structurally correct, mirrors the guard
`_node_leaves` already applies, and removes 168 edges that assert causal
coupling Reactome does not curate — concentrated in two evaluation pathways
(S Phase was 84% cofactor bridges, Mitotic G1 61%). The case is correctness,
not accuracy; the benchmark cannot see it.
