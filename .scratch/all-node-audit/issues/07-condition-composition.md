# Sample-level Condition contrast and composition

Type: task
Status: resolved

Blocked by: 06

## Initial audit focus

Typed pseudobulk counts, sample-to-condition integrity, estimable designs, replicate disclosure, formal versus exploratory inference, and compositional-model diagnostics.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: The required `official-usage.md` and `design-review.md` records landed for Pseudobulk,
  Pseudobulk EdgeR, legacy Decoupler Pseudobulk Contrast, and Pseudobulk PyDESeq2; the primary agent read all
  eight documents before any implementation. The reviewed seam introduces a typed Sample-by-population count
  artifact, keeps EdgeR QL and PyDESeq2 as separate one-population engines with a shared pure design preflight,
  and retires the removed Decoupler 1.x contrast executor behind a non-executing migration shim. Implementation
  remains blocked on completion of Batch 06.
- 2026-08-28: The scVI Differential Expression research records landed and were read in full. The node is retained
  only as one two-population `scVI Model DE Evidence` operation with explicit shared-Technical-batch support,
  posterior-FDR language, and mandatory `exploratory_model_evidence` status; it cannot replace Sample-level
  pseudobulk Condition inference.
- 2026-08-28: Sample Composition Summary and Differential Composition Test each have both required pre-change
  records, read in full. The summary is retained as a pure descriptive complete Sample-by-annotation zero grid.
  The under-specified Kruskal-Wallis/Mann-Whitney production node is slated for fail-closed retirement rather than
  being relabeled as formal compositional inference; no automatic mapping to an annotation-reference model is safe.
- 2026-08-28: Milo Differential Abundance now has its two required pre-change records, both read in full by the
  primary agent. The retained deep module is one explicit pairwise Sample-level Condition contrast: it subsets before
  graph construction, fixes Pertpy 1.3 to `solver="edger"`, validates an estimable design with at least three
  independent Samples per arm, controls both Scanpy and neighborhood RNG, treats `SpatialFDR` as primary, preserves
  majority annotation separately from a derived Mixed label, and emits the full `table, summary, code` contract.
- 2026-08-28: Schist Nested Model has exactly two pre-change records and both were read in full. It is retained under
  its stable node ID but moved from differential abundance to graph clustering: the reviewed operation is one exact
  Schist 0.10.0 unweighted undirected nested-SBM hierarchy fit on a validated named graph, not a Condition test.
  The interface returns every hierarchy level in copied AnnData plus `summary, code`, fixes one-thread/nonzero-seed
  behavior, validates marginals/nesting/stats transactionally, fails closed on package-name/version ambiguity, and
  stays out of packaged examples until the exact graph-tool environment has a real-backend smoke test.
- 2026-08-28: scCODA and tascCODA each have exactly two pre-change records, all four read in full. Both remain
  separate Sample-level Condition modules behind a shared private aggregation/design/NUTS seam. scCODA becomes one
  explicit flat two-group contrast with honest posterior expected-FDR threshold/realized quantity disclosure.
  tascCODA becomes one explicit hierarchy-aware contrast with a validated root-to-leaf manifest, signed aggregation
  bias, fixed Pertpy 1.3.0 `lambda_0=50`, `lambda_1=5`, `theta=0.5`, direct hierarchy-node versus derived-leaf scopes,
  and no inert `estimated_fdr` input or claim. Both report the pinned one-chain diagnostic limitation and emit
  `table, summary, code`; their ambiguous legacy formula-based schemas require atomic migration rejection.
- 2026-08-28: Sample Composition Summary is now a pure descriptive `table, summary, code` module. It rejects any
  missing/ambiguous metadata, preserves categorical annotation order, emits the complete Sample-by-annotation zero
  grid with explicit denominators and strict closure/count conservation, bounds both the primary grid and report
  preview, and has runtime/generated full-summary parity without reading expression values. The unsafe legacy
  Differential Composition Test now fails closed with an actionable distinction between Condition baseline and
  compositional reference annotation. Focused tests pass 8/8 under `-W error`; migration/example work remains.
- 2026-08-28: scVI Model DE Evidence now has a deep process-local model seam: the artifact fingerprints the exact
  registered count source (including a declared count layer), fitted feature/observation axes, setup metadata, and
  hashable fitted weights. One two-population contrast defaults to supported shared-Technical-batch decoding,
  `change` mode with a two-sided practical-effect event, and posterior expected-FDR disclosure. Backend results must
  contain exactly one row per fitted feature and pass a pinned 23-column schema, probability/complement, realized
  pseudocount, effect-range, Boolean tag, group-label, and finite-value postconditions. The strict JSON explicitly
  labels the result exploratory model evidence rather than Sample-level Condition inference, while generated code
  consumes the same live fitted model and never retrains. Twenty-seven focused fake-backend/legacy tests pass with
  warnings as errors, and a real scvi-tools 1.5.0.post1 one-epoch smoke test produced the expected full table using
  the ordered shared batches `b1,b2`.
- 2026-08-28: Independent frozen-core follow-up closed the remaining backend-boundary findings. EdgeR and PyDESeq2
  now receive private count/design/contrast copies, snapshot all axes, pre-existing metadata, exact count content,
  design values and contrast, and enforce postconditions against both the passed objects and any adapter-owned
  copies. Runtime and generated code reject gene/count (including fractional count)/design/contrast mutation
  identically. EdgeR preflight now verifies and reports rpy2, R, edgeR, BiocParallel and RhpcBLASctl before fitting.
  Pseudobulk role metadata rejects non-scalar values before Decoupler and the AnnData citation is corrected to
  JOSS 2021;6:4371.
- 2026-08-28: PyDESeq2 now captures backend warnings across construction, fit and test, deterministically aggregates
  repeats, and separately reports requested versus realized dispersion trends. A real Pertpy 1.3.0/PyDESeq2 0.5.4
  smoke returned all four fitted genes, disclosed `parametric` requested / `mean` realized, and normalized 95 emitted
  warnings into six records. The local EdgeR preflight stopped before fitting with the exact missing-rpy2 action and
  named all required R components.
- 2026-08-28: scVI cleanup was verified against the signed 1.5.0.post1 manager implementation. Transferred managers
  enter the per-instance store first, so the adapter retrieves the exact manager through `get_anndata_manager`,
  publishes that same object through `register_manager`, and finally calls only
  `deregister_manager(analysis_adata)`; it never invokes the broad no-argument cleanup or touches private stores.
  Success and exception tests preserve unrelated managers. Generated scVI DE source embeds its complete validator,
  canonicalizer and summary builder without importing OpenBio helpers, and has full table/summary parity. Reporting
  now includes positive/negative tagged evidence with Bayes factors, unused declared categories, registered count
  state/evidence, Gayoso et al., audited scvi-tools version, Torch/CUDA/cuDNN/MPS runtime, and immutable references.
  A real one-epoch CPU scvi-tools 1.5.0.post1 smoke returned the complete six-gene schema and left the fitted model
  attached to its registered AnnData.
- 2026-08-28: Final frozen-core verification passed 218 focused tests with warnings treated as errors (one optional
  environment-dependent test skipped), plus Ruff, `py_compile`, and `git diff --check`. Obsolete abundance-module
  helper implementations and imports were removed without changing the reviewed composition modules. Workflow
  migration/example changes remain under the separate migration owner and are not claimed by this core pass.
- 2026-08-28: Closure rerun on the integrated worktree passed 241 focused tests with `-W error` (two optional
  environment-specific tests skipped) across pseudobulk, sample composition, scVI evidence, Milo, Schist, scCODA,
  and tascCODA, plus all seven registration checks. The full workflow-migration suite passed 64/64 after the
  migration owner added exact current-schema no-ops, a safe additive Sample Composition translation, and atomic
  fail-closed handling for retired or scientifically underdetermined legacy interfaces, including nested graphs,
  object/array links, exposures, and global identifiers. Batch07 is resolved; final repository-wide checks will
  rerun these frozen contracts together with later batches.
