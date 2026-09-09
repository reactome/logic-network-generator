# logic-network-generator Constitution

The generator turns curated Reactome pathways into logic networks that
DeltaSignal solves. These principles govern what belongs here and what does
not.

## Core Principles

### I. Represent curator intent (NON-NEGOTIABLE)

Build the logic network to reflect **how curators intended the pathway to be
designed**. DeltaSignal decides how best to process it.

The corollary is the one that actually bites: when a *faithful* representation
regresses the benchmark, the fix goes in DeltaSignal's processing, not in an
unfaithful encoding here. `LNG_SET_MEMBERS_OR` is the worked example — marking
set-valued catalysts as OR is biologically correct and costs 51 of 223 cases,
because DeltaSignal has no aggregator for "one of N redundant alternatives was
lost". The change stays default-OFF until the solver can process it. Do not
resolve that tension by making the network lie.

### II. Set-splitting already encodes the OR

Because EntitySets are *split* — each alternative becomes its own virtual
reaction — a member genuinely IS required within its own VR, so `and` is
correct for inputs and outputs. OR is almost never needed.

The only exception is **sets we explicitly do not split on**: catalysts and
positive regulators, which `append_regulators` flattens onto the reaction's
single existing VR. Only that path may be given OR semantics. Never extend it
to input/output edges.

### III. Validate mapping changes through DeltaSignal

Any change to how regulators, catalysts, or connections are mapped must be
validated end-to-end against DeltaSignal before it is called an improvement.
An edge-count delta, a connectivity percentage, or a reachability figure is
not evidence of a better model — several such "wins" in this repo's history
were accuracy-neutral or negative. Report the benchmark before/after.

### IV. Verify before reporting

A claim about the pipeline must be checked against ground truth — Neo4j, the
generated artifacts, or a benchmark run — before it is stated as fact. This
principle exists because it has been violated: a cache-fingerprint feature was
reported as landed while having zero call sites, and a solver convergence
figure was reported that turned out to be the measuring bug's own artifact.

When a check disagrees with the output, establish which one is wrong before
"fixing" either.

### V. Checks must be able to fail

A validator that cannot fail is worse than none: it trains readers to ignore
it and it hides the defects underneath. Every check states how much it
actually compared, and every relaxation is tied to a named structural reason
with the underlying defect still reachable. Confirm sensitivity by corrupting
the input on purpose.

### VI. Determinism

Generation must be reproducible: same inputs, same Reactome release, same
environment, same output. Generation-affecting settings belong in the cache
fingerprint, and a cache is reused only when that fingerprint matches — never
on mere file existence.

## Additional Constraints

- **Reactome release is part of the input.** Record it; a network is only
  comparable to another at the same release. Version skew has previously been
  misread as a connectivity bug.
- **Never log or echo credentials.** `NEO4J_URL` may carry a password; redact
  before it reaches a log, an exception message, or an artifact.
- **Parameterise all Cypher.** Never interpolate, even where an argument
  parser happens to constrain the type.

## Development Workflow

- Generation-affecting changes ship **default-OFF** behind an `LNG_*` flag
  until measured, and the flag joins the cache fingerprint.
- A negative or neutral result is recorded in the spec with its numbers, not
  discarded. The record is what stops the next attempt from repeating it.
- Findings that are out of scope for the change at hand become tracked tasks
  rather than being chased inline or dropped.

## Governance

This constitution supersedes convention and habit. Specs and plans under
`specs/` are checked against it; a plan that conflicts with a principle must
either change or state the justification explicitly. Amendments require a
stated rationale and a version bump.

**Version**: 1.0.0 | **Ratified**: 2026-09-09 | **Last Amended**: 2026-09-09
