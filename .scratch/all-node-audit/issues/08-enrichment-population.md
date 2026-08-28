# Scoring, enrichment, drug analysis, and population prioritization

Type: task
Status: resolved

Blocked by: 07

## Initial audit focus

Delete/replace cell-level pathway t-test pseudoreplication, preserve the tested gene universe, disclose exploratory marker-derived enrichment, initialize/version external resources, and isolate cNMF/Augur artifacts.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: Superseding Augur boundary clarification after the user's open-expert-tool requirement: explicit Raw
  selection is sufficient and mutable analysis history is not a gate. Raw feature completeness is user-declared.
  Expert-selected X/layer and nonnegative finite non-integer expression are preserved and disclosed; biological
  Sample support, repeated-measure mappings, and Condition/Technical-batch confounding are advisory audits rather
  than execution gates. Only the two-arm/cell-support/backend/numeric/alignment/safety invariants remain hard.
- 2026-08-28: Pre-change research completed for all 13 Batch08 node IDs; each now has exactly
  `official-usage.md` and `design-review.md`. Key freezes: decoupler nodes target the 2.2 `mt.*` interfaces;
  observation-level scores are descriptive only; `PathwayScoreTTest` is replaced by one-population Sample-mean
  Welch inference; ranked/selected enrichment requires direct complete evidence plus its exact universe; legacy
  DGIdb AnnData annotation is retired in favor of an explicit version/license/hash resource artifact; drug score,
  ORA, and GSEA consume that artifact without internal marker ranking; Augur consumes Raw snapshot counts and returns
  a typed result artifact while remaining explicitly cell-CV exploratory; Augur Results becomes a pure artifact view;
  population correlation emits a canonical table plus plot. No runtime, test, dependency, workflow, or ledger file
  was changed in this research task.
- 2026-08-28: `OpenBioSingleCellCellTypeCorrelation` was retained under its stable node ID and renamed
  “Population Centroid Correlation”. The implementation is read-only, requires one explicit finite dense
  representation, independently verifies Scanpy's centroid correlation matrix and dendrogram metadata, rejects
  undefined/complex inputs, and emits a deterministic pair table, PNG, strict report, and equivalent function.
  Sample and Technical-batch limitations are explicit. Fifteen focused tests plus seven registration checks pass
  with warnings treated as errors; Ruff and bytecode compilation are clean.
- 2026-08-28: Historical implementation checkpoint, superseded by the final open-expert-boundary entries above and
  below: `OpenBioSingleCellAugur` and `OpenBioSingleCellAugurResults` were enhanced at a typed artifact seam. At this
  checkpoint Augur required a proven post-QC full-gene Raw count snapshot, explicit Sample/population/Condition arms,
  classifier-only Pertpy 1.3.0, two Samples per arm and per eligible population, and both global and population-level
  Technical-batch audits; perfect confounding is rejected and partial imbalance is disclosed. It returns immutable,
  tamper-evident canonical priorities, cell-level cross-validation, and feature-importance tables with strict JSON,
  dynamic versions/references, AUC dispersion, fingerprints, and equivalent code. Augur Results only validates and
  selects one defensive table copy. Thirty-seven focused tests (including a real Pertpy 1.3.0 CPU smoke) plus seven
  registration tests pass with warnings as errors; Ruff and bytecode compilation are clean.
- 2026-08-28: The four DGIdb nodes were rebuilt around a local typed `OPENBIO_DGIDB_RESOURCE` seam. “DGIdb
  Annotation” is now a strict CSV/TSV loader with exact ten-field release/license metadata, streamed SHA-256,
  same-size/mtime tamper checks, canonical content/accounting fingerprints, duplicate provenance aggregation, and
  generated portable code. Single-drug scoring fixes Pertpy 1.3 mean scoring and independently verifies the matched
  target arithmetic on a private copy. Drug Hypergeometric and Drug GSEA consume approved selected/complete-rank
  evidence plus its exact paired universe, reuse the hardened full-family decoupler 2.2 ORA/GSEA cores, and disclose
  target coverage, direction limits, adjusted-p semantics, and non-clinical interpretation. Twenty focused tests,
  including real Pertpy and decoupler CPU smokes and runtime/generated parity, plus a 72-test enrichment/registration
  cross-suite pass with warnings as errors; Ruff, bytecode compilation, and diff checks are clean.
- 2026-08-28: Release-status audit corrected the distinction between frozen scientific cores and workflow
  compatibility. All 13 Batch08 design reviews require either an explicit safe translation or an atomic actionable
  rejection for their breaking legacy schemas; no Batch08 rule exists yet in `workflow_migrations.mjs`. Therefore
  all 13 ledger entries remain provisional until the central migration owner covers current-schema idempotence,
  legacy object/array links, nested graphs/exposures/global IDs, mixed-family atomicity, and the documented
  no-invention rules for resource metadata, tested universes, Sample roles, and typed artifacts.
- 2026-08-28: Final open-expert-boundary pass supersedes earlier statements that described provenance confidence,
  Raw/count-like use, recommended permutations/replication, license review, or batch balance as execution gates.
  AUCell, Gaussian/empirical GSVA, panel and drug scores now preserve explicit Raw/count-like choices with warnings;
  only Poisson GSVA retains its finite nonnegative integer kernel precondition. Pathway Welch requires the true
  two-observation mathematical minimum per arm while two Samples, cross-batch Samples, and perfect confounding run
  with `formal_interpretation_invalid`. Ranked/Drug GSEA require only two permutations; low values and seed zero are
  explicit cautions. Resource/evidence producer, approval, Sample, batch, curation and license-review claims are
  caller attestations, while structure, axes, exact tested universe, current-content fingerprints, local file safety,
  and backend arithmetic remain hard. Augur likewise uses explicit source selection and advisory Sample/batch audits;
  its generated code and Population Centroid Correlation code are now self-contained.
- 2026-08-28: Central migration now covers all 13 current schemas, safe AUCell/GSVA translations, typed Augur Results
  rewiring, population-correlation slot changes, and atomic rejection where legacy resources/universes/roles cannot be
  inferred. Main-thread verification passed 143 Batch08 + registration tests under `-W error`, Ruff, pycompile and
  diff checks; the full central migration suite passed 98/98. Batch08 is frozen and resolved.
