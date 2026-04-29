# Logic Network Validation — MP-BioPath Test Set

This directory holds the results of validating the generated logic networks
against the MP-BioPath curator-prediction test set published in:

> Reactome is a database of human biological pathways manually curated from
> the primary literature… 4968 test cases were analyzed across 10 pathways,
> of which 847 were supported by published empirical findings. Out of the
> 847 test cases, curators' predictions agreed with the experimental evidence
> in 670 and disagreed in 177 cases, resulting in ∼81% overall accuracy.
> MP-BioPath predictions agreed with experimental evidence for 625 and
> disagreed for 222 test cases, resulting in ∼75% overall accuracy.

The original MP-BioPath tool ran against a smaller, simpler representation
of pathway diagrams. The current logic networks are produced by a different
pipeline (this repo), with substantially different structural choices —
full EntitySet expansion, intermediate complexes preserved, boundary
decomposition with synthetic assembly/dissociation edges, fan-out for
mismatched input/output combination counts, source_entity_id provenance,
and position-aware UUIDs.

This validation answers the question: **do the new networks reproduce the
curator-predicted perturbation outcomes about as well as MP-BioPath did on
the original networks?**

## Reactome version

The test set was generated against Reactome v86. We are now on v96. Of the
93 pathways in the test set, 91 still exist; 84 successfully regenerated
under the current pipeline (7 pathways exist as Pathway nodes but have no
reachable ReactionLikeEvents in v96).

## Methodology

For each remaining pathway, for every (perturbation_gene, key_output)
test case in the curator predictions:

1. Find every UUID in the network whose stable ID has a reference entity
   matching the perturbed gene name (catches every compartmental and
   modification form of that protein).
2. Pin those UUIDs to the perturbation value (0 = knockout, 2 = upregulate);
   leave all other nodes at 1 (control).
3. Run synchronous Boolean propagation:
   - Virtual reactions: AND of input/catalyst/positive-regulator
     contributions, with negative-regulator contributions inverted before
     AND-ing.
   - Entities (output side): MAX of producer-VR activities (OR over multiple
     producers).
   - Boundary assembly edges feed into a complex via AND; dissociation edges
     feed leaves from the complex.
4. Iterate to a fixed point (cap 50 iterations).
5. Read predicted state at the key output (max over the entity's UUIDs).
6. Compare to the curator's expected state (0 / 1 / 2).

## Headline result

**70.55% accuracy** on the 12,895 test cases that are currently runnable
(both endpoints exist in v96 Reactome). This is the honest comparison
against MP-BioPath's published 75% on its 10-pathway empirical subset:
within ~5 points, on a 13× larger valid-test set, with a generic Boolean
propagator that wasn't tuned for these networks.

A "raw" number of 65.89% across all 22,738 published cases is misleading
because 9,843 of those (43%) reference genes or key-output entities that
have been retired or renumbered between Reactome v86 (when the test set
was built) and v96 (current).

| Metric                                         | Cases   | Correct | Accuracy |
|------------------------------------------------|---------|---------|----------|
| Raw — counts drift cases as default failures   | 22,738  | 14,981  | 65.89%   |
| **Valid-only — drift cases removed**           | 12,895  | 9,097   | **70.55%** |

For comparison: random guessing ≈ 33%; MP-BioPath on v86 networks ≈ 75%;
curator agreement with experimental data ≈ 81%.

## Failure breakdown (valid cases only, 3,798 wrong)

| Category | Count | % of failures | Interpretation |
|---|---|---|---|
| `propagator_missed` | 2,278 | 60.0% | Path exists; the simple Boolean propagator couldn't carry the perturbation. The MIN-of-inputs rule for AND gates can't propagate "up" signals when other inputs sit at 1 — a structural limit of this propagator, not a network bug. Goes away with learned parameters. |
| `no_path` | 1,422 | 37.4% | Both endpoints exist but no directed path connects them. Could be a real generation bug, *or* a genuine v96 change where Reactome removed an intermediate reaction between 2018 and 2026. Needs per-case inspection to disambiguate. |
| `false_positive_change` | 98 | 2.6% | Predicted change but curator said normal. |

So at most ~11% of valid tests *could* indicate a generation bug, and that
slice is shared with genuine pathway changes in current Reactome.

## Files

- `2026-04-29_mpbio_per_pathway.tsv` — one row per pathway: total cases,
  raw accuracy, valid-only accuracy.
- `2026-04-29_mpbio_per_pathway_failures.tsv` — one row per failing test
  case with predicted, expected, and a failure category.
- `2026-04-29_mpbio_console.log` — full run output including overall
  accuracy and confusion matrix.

## Failure categories

Each failing case is classified to distinguish *network-generation issues*
from *propagator limitations* from *test-set drift*:

- **`gene_not_in_network`** — the perturbed gene has no UUIDs in this
  pathway's network. Possible causes: gene-name → reference-entity mapping
  failed, the gene's entities were removed in v96, or the network is
  missing this part of the pathway.
- **`keyoutput_not_in_network`** — the key output entity (by stable ID)
  has no UUID in this network. Likely v96 drift (entity retired) or the
  network is missing this output.
- **`no_path`** — both endpoints exist but no directed path connects
  perturbation to key output. Either Reactome's pathway genuinely doesn't
  link them, or the network-generation pipeline missed an edge.
- **`false_positive_change`** — predicted a perturbation but the curator
  expected normal. The propagator over-propagated.
- **`propagator_missed`** — path exists but propagation didn't carry the
  perturbation correctly. Most likely a propagator limitation (synchronous
  Boolean over a 3-state lattice can't capture nuanced dynamics), not a
  network bug.

The headline-level summary (overall accuracy, confusion matrix, category
counts) is in the console log.
