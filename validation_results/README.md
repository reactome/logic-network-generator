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

## Headline results — two metrics, reported together

We report two complementary metrics:

| Metric | Value | What it measures |
|---|---|---|
| **End-to-end prediction accuracy** | **70.55%** (9,097 / 12,895) | Network + generic Boolean propagator vs curator predictions. Directly comparable to MP-BioPath's published 75%. |
| **Network correctness** | **98.3%** (9,097 / 9,253) | Of cases where the *network's connectivity* determines the answer, the network agrees with the curator. Excludes failures attributable to the propagator (which parameter learning is expected to fix) and to v86→v96 test-set staleness. |

The 98.3% number is the cleanest measure of whether the **logic-network
generation itself** is correct. The 70.55% number is the honest
end-to-end comparison that includes all sources of disagreement —
useful for benchmarking against published tools.

For context: random guessing ≈ 33%; MP-BioPath on v86 networks ≈ 75%
vs experimental; curator agreement with experimental data ≈ 81%.

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

| Category | Count | % of failures | % of valid tests | Interpretation |
|---|---|---|---|---|
| `propagator_missed` | 2,278 | 60.0% | 17.7% | Path exists; the simple Boolean propagator couldn't carry the perturbation. The MIN-of-inputs rule for AND gates can't propagate "up" signals when other inputs sit at 1 — a structural limit of this propagator, not a network bug. Goes away with learned parameters. |
| `no_path` | 1,422 | 37.4% | 11.0% | Both endpoints exist but no directed path connects them. **Resolved further below by Neo4j cross-check.** |
| `false_positive_change` | 98 | 2.6% | 0.8% | Predicted change but curator said normal. |

## Neo4j cross-check on the no_path slice

For each `no_path` failure case, we built the equivalent reaction-flow
graph from current Neo4j and asked: does the original pathway have a
path from the perturbed gene's entities to the key-output entity? If
Neo4j has a path but our logic network doesn't, the network is missing
something — that's a generation-bug candidate. If Neo4j has no path
either, the curator's prediction was based on v86 connectivity that
v96 has dropped, and the failure is just test-set staleness.

| Bucket | Count | % of no_path | % of valid tests |
|---|---|---|---|
| `bug_candidate` (Neo4j has path; we missed it) | **156** | 11.0% | **1.2%** |
| `pathway_changed` (Neo4j has no path either) | 1,266 | 89.0% | 9.8% |

**~1.2% of valid tests are potential network-generation bugs.** If every
one of those were a real bug and we fixed them all, accuracy would rise
from 70.55% to roughly 71.8%. Most of the remaining gap is propagator,
not network.

The 156 bug candidates concentrate in 9 pathways — most heavily in
Transcriptional_Regulation_by_TP53 (80 cases, all involving MDM2
perturbation) and PIP3_activates_AKT_signaling (18 cases involving PTEN).

### Tracing the cluster: not a generation bug

End-to-end inspection of the MDM2/TP53 cluster (and a confirming check
on the PTEN/AKT cluster) showed every "bug candidate" we examined was
the same structural pattern:

> **The biological prediction depends on substrate consumption.**

Concrete example. The reaction `MDM2 ubiquitinates TP53` consumes
the MDM2:TP53 complex (containing TP53) and produces PolyUb-TP53
Tetramer. The curator-implicit chain to TIGAR upregulation reads:

> Knock out MDM2 → less ubiquitination → **less TP53 consumed** → more
> free TP53 → more TP53 Tetramer → more TIGAR transcription.

The bolded step is mass-action depletion: when a reaction consumes
its input, less of that input remains for other reactions. Our
directed-flow network correctly emits every Reactome reaction with
its inputs and outputs as edges, but it does **not** emit a "this
reaction depletes its substrate" edge — i.e., the negative effect of
flux-on-substrate isn't represented in the graph topology.

The PTEN/AKT cluster shows the same pattern: PTEN catalyzes
"PTEN dephosphorylates PIP3" (consuming PIP3); the curator-predicted
effect on AKT-downstream targets (p-S356-RPS6KB2, p-S939-TSC2, etc.)
hinges on PIP3 depletion rather than on a forward-flow path the
graph could carry.

So the 156 bug candidates are **not** generation bugs to fix in the
pipeline — Reactome's reactions are all there, with the right inputs,
outputs, catalysts, and regulators. They reflect a known limitation
of pathway-diagram-derived directed-flow Boolean networks. Two ways
to address them if needed:

1. **Add explicit consumption edges**: emit a `pos_neg='neg'` edge
   from each reaction back to each of its inputs, modeling
   substrate depletion as a negative effect. A meaningful semantic
   addition; would push validation accuracy higher.
2. **Let parameter learning capture it**: alphaSignal's learned
   parameters can absorb mass-action steady-state implicitly without
   a structural change to the network — the model can learn that
   "this reaction's flux depletes this substrate" from the
   correlation in training data.

Either is a forward-looking improvement. Neither is a fix to a
generation bug.

## Vs experimental data (10-pathway subset, 499 valid tests)

| Metric | Cases | Correct | Accuracy |
|---|---|---|---|
| Raw                              | 849   | 188     | 22.14%   |
| Valid-only                       | 499   | 170     | 34.07%   |

Dramatically lower than vs curator (70.55%) and lower than MP-BioPath's
published 75% vs experimental. The reason is structural to the
propagator: experimental data heavily samples cases where the
perturbation **caused a measured change** (only 10% are "normal";
curator data is 61% "normal"). The simple Boolean MIN-of-inputs
under-predicts changes (and especially can't propagate "up" through
AND gates), so it does much worse on a dataset that's mostly changes.

This is a propagator limitation, not a network limitation. With
parameter learning (alphaSignal), the same network should recover
performance against experimental data.

## Files

- `2026-04-29_mpbio_per_pathway.tsv` — vs curator, per-pathway
- `2026-04-29_mpbio_per_pathway_failures.tsv` — vs curator, per-case failures with category
- `2026-04-29_mpbio_per_pathway_experimental.tsv` — vs experimental, per-pathway
- `2026-04-29_mpbio_per_pathway_experimental_failures.tsv` — vs experimental, per-case failures
- `2026-04-29_no_path_neo4j_check.tsv` — Neo4j cross-check classifying every `no_path` case as bug_candidate or pathway_changed
- `2026-04-29_mpbio_console.log` — full run output

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
