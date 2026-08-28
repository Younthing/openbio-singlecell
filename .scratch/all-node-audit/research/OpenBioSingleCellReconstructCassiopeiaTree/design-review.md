# Cassiopeia VanillaGreedy reconstruction: module design review

## Decision

**Retain and deepen `OpenBioSingleCellReconstructCassiopeiaTree` as one computable character payload to one VanillaGreedy topology.** The small interface selects one lineage and the two topology-affecting solver policies; the implementation hides dependency/signature gates, defensive ownership, mutation guards and graph validation. Parsing and character conversion stay upstream; expansion and plasticity stay downstream.

Supporting Neighbor Joining in the same dropdown was rejected: it has a distinct distance function, optional synthetic root, observed-outgroup semantics, fast external implementation and thread interface. ILP/Hybrid additionally introduce Gurobi license, seed, time, logfile and resource concerns. Those solvers require separate nodes and their own two pre-change documents.

## Current interface and artifact

Inputs, in order:

1. `characters` (`OPENBIO_CASSIOPEIA_CHARACTERS`);
2. exact canonical `lineage_id`, required to have a computable payload; upstream QC warning status is disclosed but is not an execution gate;
3. `prior_transformation`, enum `negative_log|inverse|square_root_inverse`, default `negative_log`;
4. `collapse_mutationless_edges=False`.

Outputs are `tree`, `summary`, and `code`.

The node never reads a file. It validates the current character artifact fingerprints, uses defensive matrix/prior copies, fixes `missing_data_classifier` to the official `assign_missing_average` callable, and calls one public VanillaGreedy solve. It serializes solver execution behind a process lock, snapshots/restores global NumPy state, and rejects input/backend mutation. Upstream QC warnings are carried into the report; only fewer than two matrix cells or another backend/structural invariant blocks computation.

`OPENBIO_CASSIOPEIA_TREE` owns a defensive solved-tree copy plus immutable canonical payload: lineage ID, exact ordered leaf/character axes, missing indicator, character/prior/upstream fingerprints, solver class/configuration, solver-generated root, canonical sorted directed edges and branch lengths, topology fingerprint, dependency identity and diagnostics. `copy_tree()` validates current contents and returns a copy; consumers never share the stored backend object. Validation catches forged class/schema/producer values and current-content tampering.

The strict `summary` contains report-ready methods/results, key topology diagnostics, exact upstream linkage, versions, citations, warnings and limitations without NaN/Infinity. Standalone `code` accepts the typed character artifact and returns `(tree_artifact, summary_dict)` with structural/numerical parity.

## Legacy/current migration matrix

The legacy node accepted eight file widgets and emitted only `tree`; the current node accepts four typed/scientific widgets and emits `tree, summary, code`. A safe migration must atomically insert or reuse a current QC/prepare node for the exact same file and wire its `characters` output, then map legacy `tumor` to `lineage_id` only after canonical lineage identity is proven. Legacy `allele_representation_threshold` belongs on the inserted QC node. The file-column widgets likewise move upstream.

Automatic migration is blocked because the legacy `mutation_family_column` was used as an unproven empirical-prior grouping, its cut-site schema/missing allele/conflict policy were absent, and no QC artifact/root/version fingerprint existed. The only executable migration is one carrying an explicit reviewed preparation payload satisfying the QC migration matrix. The old solver behavior maps to `prior_transformation=negative_log` and `collapse_mutationless_edges=False`; the root is recorded as `solver_generated`, never relabeled synthetic. Existing output links may remain on tree output 0; summary/code are appended. Graph mutation must be atomic/idempotent and preserve all link encodings and scopes. Connected widget ports are schema evidence too: their concrete STRING/BOOLEAN/FLOAT types must match the reviewed legacy or current widget contract before any mutation.

## Verification gate

Tests cover artifact type/schema/producer/tamper checks, exact lineage selection, failed QC, all-missing/uninformative rejection upstream, callable and transformation wiring, collapse policy, canonical input non-mutation, RNG restoration, backend mutation/malicious cyclic/disconnected/leaf-mismatch results, duplicate profiles, exact topology hash, defensive copies, strict JSON, compiled standalone code and parity. The real smoke solves a small informative matrix with `cassiopeia-mt==2.1.3` and repeats it to confirm the canonical topology.
