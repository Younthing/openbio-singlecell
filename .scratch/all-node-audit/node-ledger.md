# Node audit ledger

Last updated: 2026-08-29

Legend:

- `[x]` inventory and initial code audit complete.
- `[R]` required official-usage and design-review documents landed; implementation and tests complete.
- `[ ]` pre-change research/design documents and implementation still pending.
- `adapter` nodes are audited for correctness and interface depth but do not receive forced scientific `summary/code` outputs unless they perform an analysis.

Current-only repository invariant: exactly 94 unique registered node IDs across 23 node modules. Reporting scope is
85 processing/analysis/report/view nodes; 5 input/config adapters and 4 output adapters are audit-only. No deprecated,
migration-only, or compatibility-shim node is registered.

## Batch 01 — quality control

- [R] `OpenBioSingleCellCalculateQC`
- [R] `OpenBioSingleCellFilterCells`
- [R] `OpenBioSingleCellFilterGenes`
- [R] `OpenBioSingleCellQCPlots`

## Batch 02 — QC decisions and artifact correction

- [R] `OpenBioSingleCellMarkMADOutliers` — corrected SciPy scaling; per-Sample thresholds disclosed.
- [R] `OpenBioSingleCellScrublet` — explicit count source and per-Sample/capture execution.
- [R] `OpenBioSingleCellFilterDoublets` — strict boolean prediction contract.

## Batch 03 — input, study, and expression state

- [R] `OpenBioSingleCellLoadH5AD` — strict cell identity and portable embedded provenance.
- [R] `OpenBioSingleCellLoad10xMTX` — strict count/identity contract; variable-only name repair.
- [R] `OpenBioSingleCellLoad10xStudy` — strict complete Sample inventory and feature-alignment disclosure.
- [R] `OpenBioSingleCellLoad10xH5` — official reader with strict count/identity and filter disclosure.
- [R] `OpenBioSingleCellCoreStudyParameters` — retained unchanged as a configuration adapter.
- [R] `OpenBioSingleCellAnnDataSummary` — bounded strict structural report and equivalent code.
- [R] `OpenBioSingleCellSnapshotExpression` — sole canonical Raw snapshot and `counts` owner.
- [R] `OpenBioSingleCellSubsetObservations` — exact selection with explicit missing-value semantics.
- [R] `OpenBioSingleCellMergeObservationAnnotations` — explicit conflict, provenance, and dtype policy.
- [R] `OpenBioSingleCellMapGeneIdsFromGTF` — stable-ID-preserving symbol annotation without feature filtering.

## Batch 04 — normalization and feature modeling

- [R] `OpenBioSingleCellNormalizeTotal` — explicit count source; no implicit snapshot or stale log marker.
- [R] `OpenBioSingleCellLog1p` — atomic active-X transform; never creates log-normalized `raw`.
- [R] `OpenBioSingleCellNormalizeToLayer` — one declared count source to one protected derived layer.
- [R] `OpenBioSingleCellPearsonResidualsToLayer` — strict count model with explicit clipping and dense-memory guard.
- [R] `OpenBioSingleCellHighlyVariableGenes` — flavor-specific state, grouping, forced-gene, and Raw-subset contracts.
- [R] `OpenBioSingleCellScale` — layer-preserving centered or sparse variance scaling with memory preflight.
- [R] `OpenBioSingleCellCNMFRankSurvey` — method-author cNMF 1.7.1 managed-workspace survey with immutable typed run artifact and explicit expert-selected count source.
- [R] `OpenBioSingleCellCNMF` — method-author cNMF 1.7.1 consensus with continuous usages, complete file/axis validation, strict report, and equivalent source.

## Batch 05 — representation, Technical batch integration, graph, and clustering

- [R] `OpenBioSingleCellPCA`
- [R] `OpenBioSingleCellHarmonyIntegration`
- [R] `OpenBioSingleCellSCVIIntegration`
- [R] `OpenBioSingleCellNeighbors` — representation dimensions are preflighted against the selected latent space.
- [R] `OpenBioSingleCellUMAP`
- [R] `OpenBioSingleCellTSNE`
- [R] `OpenBioSingleCellForceDirectedGraph` — initialization and layout-engine contracts are explicit; no silent fallback.
- [R] `OpenBioSingleCellLeiden`
- [R] `OpenBioSingleCellLeidenResolutionSweep`
- [R] `OpenBioSingleCellPCAMetadataAssociations` — Sample-level SciPy diagnostics replace removed decoupler APIs.

## Batch 06 — Cluster marker evidence and annotation

- [R] `OpenBioSingleCellMarkerGenes` — complete-family marker evidence plus exact tested universe.
- [R] `OpenBioSingleCellFilterMarkerGenes` — direct paired-universe filtering with tamper checks.
- [R] `OpenBioSingleCellMarkerExpressionPlot` — read-only selected-data plotting with isolated jitter RNG.
- [R] `OpenBioSingleCellUMAPPlot` — read-only embedding evidence with explicit color semantics.
- [R] `OpenBioSingleCellCellTypistAnnotation` — provisional labels, matched confidence, model SHA, atomic ownership.
- [R] `OpenBioSingleCellMarkerORAEvidence` — explicit-set decoupler 2.x evidence; no label commitment.
- [R] `OpenBioSingleCellMapClusterAnnotations` — reviewed mapping/commitment with complete disclosure.

## Batch 07 — Sample-level contrasts and composition

- [R] `OpenBioSingleCellPseudobulk` — immutable Sample-by-population count artifact with complete design and count-state provenance.
- [R] `OpenBioSingleCellPseudobulkEdgeR` — one-population Sample-level edgeR quasi-likelihood contrast with explicit R runtime preflight.
- [R] `OpenBioSingleCellPseudobulkDESeq2` — one-population Sample-level PyDESeq2 contrast with realized-dispersion disclosure.
- [R] `OpenBioSingleCellSCVIDifferentialExpression` — provenance-bound posterior model evidence, explicitly not Condition inference.
- [R] `OpenBioSingleCellSampleCompositionSummary` — pure complete Sample-by-annotation descriptive grid with structural-zero accounting.
- [R] `OpenBioSingleCellSchistNestedModel` — exact Schist nested-SBM hierarchy fit on one validated named graph.
- [R] `OpenBioSingleCellMiloDifferentialAbundance` — Pertpy 1.3 Sample-level neighborhood contrast with SpatialFDR and batch audit.
- [R] `OpenBioSingleCellSccodaDifferentialComposition` — explicit two-Condition Sample-level flat compositional model with posterior diagnostics.
- [R] `OpenBioSingleCellTasccodaDifferentialComposition` — hierarchy-bound Sample-level compositional model with signed aggregation validation.

## Batch 08 — scoring, enrichment, drugs, and population prioritization

- [R] `OpenBioSingleCellAUCellScores`
- [R] `OpenBioSingleCellGSVAScores`
- [R] `OpenBioSingleCellGenePanelScores`
- [R] `OpenBioSingleCellPathwayScoreTTest` — replaced by one-population Sample-mean Welch evidence.
- [R] `OpenBioSingleCellRankedGSEA`
- [R] `OpenBioSingleCellGeneSetOverrepresentation`
- [R] `OpenBioSingleCellDGIdbAnnotation` — typed pinned local resource loader.
- [R] `OpenBioSingleCellDrugScores`
- [R] `OpenBioSingleCellDrugHypergeometric`
- [R] `OpenBioSingleCellDrugGSEA`
- [R] `OpenBioSingleCellAugur`
- [R] `OpenBioSingleCellAugurResults` — pure typed-artifact view.
- [R] `OpenBioSingleCellCellTypeCorrelation` — Population Centroid Correlation.

## Batch 09 — regulation and communication

- [R] `OpenBioSingleCellCollecTRIULM` — typed decoupler 2.2.0 core with an explicit resource and complex-policy contract.
- [R] `OpenBioSingleCellRankTFActivities` — complete-family typed-artifact core.
- [R] `OpenBioSingleCellImportPySCENICResults` — strict manifest/SHA-bound typed importer.
- [R] `OpenBioSingleCellSCENICRegulonSpecificity` — exact RSS typed-artifact core.
- [R] `OpenBioSingleCellSCENICActivityBinarization` — deterministic typed binary core.
- [R] `OpenBioSingleCellSCENICTFModules` — final motif-pruned membership view.
- [R] `OpenBioSingleCellLianaCommunication` — typed by-Sample LIANA 1.9.0 core.
- [R] `OpenBioSingleCellLianaDotPlot` — typed-result-only pure view.

## Batch 10 — trajectory, velocity, CNV, and lineage

- [R] `OpenBioSingleCellCellCycleScore`
- [R] `OpenBioSingleCellDiffusionMap` — canonical consumed graph and independently verified transition eigenpairs.
- [R] `OpenBioSingleCellPAGA` — scratch-isolated categorical v1.2 abstraction with independently recomputed weights/MST and bounded report.
- [R] `OpenBioSingleCellDPT` — graph-bound verified DiffMap, exact diffusion distance, and strict typed/stable root.
- [R] `OpenBioSingleCellVelocityFilterAndNormalize`
- [R] `OpenBioSingleCellVelocityMoments`
- [R] `OpenBioSingleCellEstimateVelocity`
- [R] `OpenBioSingleCellVelocityGraph`
- [R] `OpenBioSingleCellRecoverDynamics`
- [R] `OpenBioSingleCellVelocityGeneRanking`
- [R] `OpenBioSingleCellVelocityStreamPlot`
- [R] `OpenBioSingleCellInferCNV` — explicit finite current-axis source/reference/Sample/assembly, unverified feature-completeness disclosure, and typed state.
- [R] `OpenBioSingleCellCNVPCA` — verified ARPACK/uncentered typed-state reduction with bounded output.
- [R] `OpenBioSingleCellCNVScore` — independently verified group score over matched typed state plus downstream annotations.
- [R] `OpenBioSingleCellCassiopeiaLineageQC` — exact-character/QC core with a defensive typed artifact.
- [R] `OpenBioSingleCellReconstructCassiopeiaTree` — typed VanillaGreedy reconstruction core.
- [R] `OpenBioSingleCellCassiopeiaExpansionTest` — independently verified within-tree statistic with complete-family BH.
- [R] `OpenBioSingleCellCassiopeiaPlasticity` — published edge-denominator/preprocessing implementation.

## Batch 11 — output adapters and release audit

- [R] `OpenBioSingleCellPreviewResult` — typed bounded UI preview with stable workflow/node identity and verified atomic PNG replacement.
- [R] `OpenBioSingleCellSaveH5AD` — verified failure-atomic native H5AD persistence with portable compression control.
- [R] `OpenBioSingleCellExportCSV` — failure-atomic UTF-8 export preserving explicit, collision-safe row-index identity.
- [R] `OpenBioSingleCellSavePNG` — full-decode validation and exact-byte failure-atomic durable PNG publication.

## Removed after audit — not registered

- [R] `OpenBioSingleCellUseExpressionLayer` — source selection belongs to each consuming analysis; the pass-through state switch was deleted.
- [R] `OpenBioSingleCellNormalizeGeneNames` — an unsafe mechanical identity editor; evidence-backed GTF annotation remains.
- [R] `OpenBioSingleCellMarkerORAAnnotation` — split into marker ORA evidence and explicit reviewed annotation commitment.
- [R] `OpenBioSingleCellDecouplerPseudobulkContrast` — obsolete decoupler 1.x/underspecified contrast; current Sample-level methods remain.
- [R] `OpenBioSingleCellDifferentialCompositionTest` — invalid cell-level rank-test design; current Sample-level composition methods remain.
- [R] `OpenBioSingleCellRunPySCENIC` — obsolete in-process environment contract; the strict external-result importer remains.
- [R] `OpenBioSingleCellLianaResults` — redundant hidden-state extractor; the producing node returns the typed result directly.
- [R] `OpenBioSingleCellCNVStructure` — non-atomic composite; its five independently meaningful current nodes remain.

## Completion checks

- [x] Exactly two pre-change documents exist for every modified node.
- [x] Every processing/analysis node exposes primary result, `summary`, and `code` (or an explicitly documented merge/delete decision).
- [x] All `summary` values satisfy the strict report schema and are JSON-serializable without NaN/Infinity.
- [x] All `code` values compile and method-equivalence tests cover each implementation family.
- [x] Public workflows, README, release manifest, frontend payload, and output adapters match the final interface.
- [x] Full tests and lint pass in the configured ComfyUI environment.
