# Design Decisions

Behaviors that look surprising at first but are intentional. Read this before assuming something is a bug.

## Complex vs EntitySet — they look alike, they aren't

Both are written `{A, B}` in informal notation but they mean opposite things in biology and the pipeline treats them oppositely. Conflating them is the single most common way to misread this codebase.

| | **Complex** | **EntitySet** |
|---|---|---|
| Semantic | A bound species (`A:B` dimer is a distinct molecule) | A role marker — *any one of* `{A, B}` plays this role |
| Cross-reaction matching | **Atomic.** R1 producing complex `A:B` and R2 consuming free `A` are *not* linked. They are different species. To get free A from `A:B` you need a dissociation reaction. | **See-through.** R1 producing set `{A, B}` and R2 consuming free `A` *are* linked — the set's A-alternative *is* free A. |
| `break_apart_entity` returns | `{A:B}` (itself, intact) when it's a simple complex; cartesian-product UIDs when it contains an EntitySet | A flat set of alternatives `{A, B}` |

If you find yourself thinking "the matcher needs to see components inside a complex," stop. That's the bug, not what you're fixing — it would create biologically false links between unrelated species.

## Two layers of decomposition — don't collapse them

There are two distinct artifacts. Different rules.

- **`decomposed_uid_mapping.csv` — plumbing.** Used by `find_best_reaction_match` to compute which input combinations of a reaction line up with which output combinations (Hungarian assignment over component overlap). EntitySets are decomposed here so cross-reaction linking through alternatives works. Simple Complexes are *not* decomposed here, by design (see above).
- **`logic_network.csv` — output.** What downstream simulation consumes. Boundary complexes (root inputs, terminal outputs) have synthetic assembly/dissociation edges to their leaf components so individual proteins can be perturbed and read; intermediate complexes are kept as single nodes (they're real species in the pathway).

Adding rows to `decomposed_uid_mapping.csv` does not "expose" an entity to the final output, and removing nodes from `logic_network.csv` doesn't change the matcher's behavior. They are independent.

## EntitySet expansion produces multiple virtual reactions

Reactome's `EntitySet` represents biological alternatives — "any one of {A, B, C}" — and `Complex` represents a structured combination — "A bound to B". When a reaction's input is an EntitySet (or a Complex containing one), the logic network expands it: one virtual reaction per concrete combination of members.

**Example.** Reaction 69598 (Ubiquitination of phosphorylated CDC25A) in Neo4j has inputs `[68524, 9943734]`, where 68524 is an EntitySet of 14 ubiquitin alternatives and 9943734 is an EntitySet of 2 CDC25A alternatives. The logic network produces 14 × 2 = 28 virtual reactions, one per combination.

**Consequence.** Direct 1:1 reaction matching against Neo4j is ~40% on pathways with many EntitySets — the rest don't fail, they expand. Use `decomposed_uid_mapping.source_entity_id` to trace any expanded component back to the EntitySet it came from.

## Source entity tracking

`decomposed_uid_mapping` carries two columns that exist purely for traceability:

- `source_entity_id`: the parent EntitySet or Complex this row was decomposed *from*. `None` for entities that weren't decomposed.
- `source_reaction_id`: reserved — will hold the original Reactome reaction ID once virtual-reaction creation populates it.

Reconstruction logic: if a set of `component_id`s share a single non-null `source_entity_id`, they came from one decomposed parent and the original input was that parent. Otherwise they were independent entities.

## Loop count differs from Reactome

Reactome's pathway 69620 contains 5 entity-level cycles. The generated logic network for the same pathway contains 1. This is correct.

Reactome's loops live at the **complex level** — e.g., the MAD2-kinetochore cycle is `Kinetochore:Mad1:MAD2*` → `Mad1:kinetochore` → `Kinetochore:Mad1:MAD2`, all complexes. After decomposition those complexes don't exist as nodes; their components do, and the components don't form the same cycle. Empirically, 12 of 14 entities involved in Reactome's 5 loops on 69620 don't appear as nodes in the decomposed network at all — they were complexes that got broken apart.

The single remaining loop in the generated network represents true component-level feedback (a component that genuinely cycles back to itself, e.g., free ubiquitin recycling).

If complex-level loops matter for some downstream analysis, the answer is to keep complexes intact (don't decompose them) — not to add edges to the current network.

## Reactions without inputs or outputs are dropped

Pathway 69620 has 63 reactions in Reactome; 50 (79.4%) appear in the logic network. The 13 missing ones have no inputs or no outputs in Neo4j — typically pure regulatory reactions or polymerizations. A logic network is built from input → output edges, so these can't be represented and are skipped.
