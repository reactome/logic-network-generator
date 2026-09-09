# Feature Specification: Validator fidelity to the generated representation

**Feature Branch**: `fix/adversarial-review`
**Created**: 2026-09-09
**Status**: Implemented (one finding open)
**Input**: Adversarial review of logic-network-generator, 2026-09-08/09

## Problem

`scripts/validate_logic_network.py` is the only automated check that the
generated logic network faithfully represents Reactome. It was **completely
inert**: `load_files()` expected the legacy flat layout
`output/pathway_logic_network_{dbId}.csv`, while the generator has since moved
to `output/<Name>_R-HSA-<id>/logic_network.csv` plus
`stid_to_uuid_mapping.csv` and `cache/decomposed_uid_mapping.csv`. Loading
aborted before any check ran, so every validation — including the
warn-to-FAIL hardening added earlier in this same review — never executed.

Once the loader was fixed, the checks failed on *healthy* networks. The cause
was a single mistaken assumption running through the whole script: that a
correct logic network is an **id-for-id mirror** of Reactome. It is not. The
generator deliberately re-represents Reactome:

| Reactome construct | Generated representation | Naive comparison reports |
|---|---|---|
| `EntitySet` | *split* — the set id is never a node, each member is | set "missing" |
| decomposed `Complex` | represented by its components | complex "missing" |
| dual-role regulator | one `pos` edge and one `neg` edge, on different reactions | polarity "wrong" both ways |
| set-valued catalyst | flattened onto the target VR as N member edges | catalyst "missing" |
| dbId vs stId | per-pathway layout carries stIds | reconstruction 0.0% |

A validator that fails on correct output is worse than no validator: it trains
the reader to ignore it, and it masks the real defects underneath.

## Architecture principle (Adam, 2026-09)

> We want to make the LNGs in a way that represents them based on how curators
> intended to design them. And DeltaSignal would figure out how best to
> process them.

Validation follows from this: the validator must check **faithfulness to
curator intent under the generator's chosen representation**, not
representational identity with the Neo4j graph. Every excused difference must
be excused for a stated structural reason, and each such reason must still
leave a way for the underlying defect to be caught.

## Requirements

- **FR1** — The validator MUST locate the generated files under the current
  per-pathway layout, and continue to support the legacy flat layout.
- **FR2** — The validator MUST query whichever Reactome identifier the uuid
  mapping carries (`stId` or `dbId`), never assume one.
- **FR3** — Entity and catalyst coverage MUST treat an `EntitySet` as covered
  when its non-set leaf members are present, and a decomposed `Complex` as
  covered when it is recorded in the decomposed mapping.
- **FR4** — A set that was only *partially* split MUST still FAIL, naming the
  members that did not survive. Excusing splitting must not excuse dropping.
- **FR5** — Regulator polarity MUST be checked per `(regulator, reaction)`
  pair. Reactome legitimately curates one entity as a positive regulator of
  one reaction and a negative regulator of another.
- **FR6** — All Cypher MUST be parameterised, not string-interpolated.
- **FR7** — Every check MUST be demonstrably non-vacuous: it reports how many
  items it actually compared, and a deliberately corrupted network must fail it.

## Measurements

Against Neo4j Reactome Release97, on the three benchmark-eligible pathways.

| pathway | before port | after loader fix | after FR3/FR5 |
|---|---|---|---|
| R-HSA-69620 | could not load | 9/11 | **11/11** |
| R-HSA-453279 | could not load | 10/11 | **11/11** |
| R-HSA-1257604 | could not load | 9/11 | **10/11** (1 real finding) |

Reconstruction accuracy on R-HSA-69620 went from a spurious **0.0%** — an
id-space artifact that would have tripped the new FAIL floor on a healthy
network — to a measured **88.2%** (120/136).

Triage of the failures that the port surfaced, all verified against Neo4j
rather than assumed:

- **9 "missing entities" (R-HSA-69620): artifact.** All nine are `EntitySet`s;
  28 of 28 leaf members are present in the network. Confirms FR3.
- **11 + 3 "partially split sets": artifact.** The absent members are
  decomposed complexes, recorded in `decomposed_uid_mapping.csv`. Confirms the
  decomposed-mapping half of FR3.
- **PI5P (R-ALL-1806240) polarity: artifact.** Genuinely both a
  `PositiveRegulation` regulator (1 reaction) and a `NegativeRegulation`
  regulator (1 reaction) in R-HSA-1257604. Confirms FR5.
- **mTORC2 (R-HSA-198626): REAL.** See below.

Non-vacuity evidence (FR7): the polarity check compares 271 of 271 regulator
edges on R-HSA-1257604 with 0 unattested; injecting 5 flipped `pos_neg` values
into a copy of the network produces exactly 5 failures with correct
attribution.

## Open finding

**F1 — one reaction is silently dropped from R-HSA-1257604.**
`R-HSA-9980233` ("PIP3 activates mTORC2") appears in **no** generated artifact:
88 of the pathway's 89 reactions are represented, this one is not. It is an
ordinary non-disease *Homo sapiens* `Reaction`; it is not a self-loop
(inputs PIP3 `R-ALL-179838` + mTORC2 `R-HSA-198626`, output
mTORC2:PIP3 `R-HSA-9980228`); its output complex decomposes normally; and it
carries `precedingEvent` links in both directions. mTORC2's five components
are all present in the network via other reactions, so only this reaction and
the mTORC2 complex node itself are lost.

Cause not yet established. The pathway has no diagram of its own
(`representedPathway` returns nothing), so diagram-driven filtering is a
candidate but unconfirmed. Deliberately **not** chased inside the validator
port — it is a generator question, tracked as T007.

## Out of scope

- Fixing F1 (tracked separately; needs its own before/after DeltaSignal
  measurement per the architecture principle).
- The reconstruction gap below 90% (88.2%): 16 Neo4j edges absent, 53 extra.
  Warned, not failed, pending its own triage.
