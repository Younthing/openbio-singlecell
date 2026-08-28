# Cassiopeia EffectivePlasticity: module design review

## Decision and current interface

**Retain and deeply correct `OpenBioSingleCellCassiopeiaPlasticity` as one fixed-tree/one-categorical-state descriptive analysis.** It annotates a defensive AnnData copy and returns a transparent cell table. It does not reconstruct the tree, recluster cells, compare Conditions or mutate the shared tree artifact.

Inputs, in order:

1. `adata`;
2. `tree` (`OPENBIO_CASSIOPEIA_TREE`);
3. exact `annotation_key`;
4. `annotation_status`, enum `unknown|provisional|curated`;
5. `analysis_mode`, enum `exploratory|report_grade`;
6. `minimum_state_fraction=0.025`;
7. advanced `output_key="sc_effective_plasticity"`;
8. advanced `overwrite_existing=False`;
9. advanced `max_working_gib=4` covering the AnnData copy and topology × state dynamic program.

Outputs are `adata`, `table`, `summary`, and `code`.

Report-grade mode checks `adata.uns['openbio_singlecell']['annotations'][annotation_key]` provenance against the caller-declared status. Missing, malformed, or mismatched provenance no longer blocks an otherwise computable expert analysis: it produces a prominent warning and remains a report limitation. Curated remains a caller/workflow assertion.

## Exact implementation and ownership

Validate the current tree artifact, exact unique string identities and a conservative working-memory estimate before copying. Determine state frequencies over original tree leaves, retain states with `frequency >= minimum_state_fraction`, prune excluded/missing-annotation leaves on a defensive topology, collapse unifurcations, and insert deterministic state-group nodes for terminal exhausted polytomies as in the released study code. The inserted node IDs are collision-safe and canonical. AnnData omissions and null annotations are disclosed exclusions. At least two retained leaves and one edge are computationally required; a one-state result is allowed, reported as zero plasticity, and warned as lacking observed state diversity.

Compute every rooted-subtree Fitch-Hartigan minimum in one bounded dynamic program and divide by that subtree's edge count. For each retained leaf average the defined scores for non-leaf ancestors on the root-to-leaf path. Independently call the official backend for the complete corrected tree under a locked/restored RNG state and require its global parsimony to match. All scores must be finite in `[0,1]`.

The output AnnData preserves axes and source annotation, adds the score plus a fixed `<output_key>_status` column and bounded OpenBio provenance, and honors collision policy. The table reports every AnnData cell with cell ID, lineage, exact state, original state count/fraction, inclusion/status, path-subtree count and score. Excluded/non-tree cells have null score. Tree contents and AnnData input remain unchanged.

The strict report includes report-ready methods/results, global and single-cell key results, annotation and topology provenance, versions, citations, warnings and limitations without NaN/Infinity. Standalone `code` returns `(annotated_adata, cell_table, summary_dict)` with identical structure/numerics.

## Legacy/current migration matrix

Legacy inputs were `adata, tree, annotation_key, output_key, summary_key`; current inputs replace mutable `summary_key` with `annotation_status, analysis_mode, minimum_state_fraction, overwrite_existing`. Legacy output was only `adata`; current outputs append `table, summary, code`.

Legacy `scPlasticity` values are methodologically invalid because they divided global parsimony by nodes, subtree parsimony by leaves, omitted state filtering/pruning/polytomy preprocessing, and wrote zero-like semantics. They cannot be renamed or reused. Migration must mark the old value stale and require a rerun. `annotation_key` maps directly; `output_key` maps only after an explicit collision review. `summary_key` is removed. Without explicit annotation status/mode, migration defaults to `unknown/exploratory`, never Curated/report-grade. No automatic migration may overwrite an existing output column. All graph changes must be atomic/idempotent and preserve link/scope encodings. Migration preflight also validates the exact STRING/BOOLEAN/FLOAT type of every connected widget port rather than classifying a node from port names alone.

## Verification gate

Tests use hand-checkable rooted trees to distinguish edge/node/leaf denominators; verify 2.5% equality, rare-state pruning, root-preserving unifurcation collapse, exhausted terminal polytomies, zero-edge exclusion and path averaging; compare the dynamic program with the official backend; cover null/numeric/non-categorical states, duplicate/missing IDs, AnnData extra cells, exact category identity, provenance/mode rules, collisions, defensive copies, RNG restoration, malicious backend disagreement, strict JSON, compiled code and parity. A real `cassiopeia-mt==2.1.3` smoke covers official parsimony.
