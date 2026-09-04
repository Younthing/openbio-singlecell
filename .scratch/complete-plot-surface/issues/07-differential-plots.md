# Add differential-expression and differential-abundance plots

Type: task
Status: resolved
Blocked by: 01, 03

## Scope

- Pseudobulk QC, edgeR/PyDESeq2 contrast, scVI DE Evidence, Milo, scCODA, and tascCODA Plot adapters.

## Acceptance criteria

- Sample-level Condition contrasts remain distinct from cell-level exploratory evidence.
- Posterior, frequentist, neighborhood, and compositional quantities retain method-specific names and limitations.

## Comments

- 2026-09-04: Claimed with the abundance evidence contract as a blocker.
- 2026-09-04: The Milo and composition Diagnostic-evidence contracts were handed off with exact typed readers,
  axes, fingerprints, method metadata, posterior draws, and sampler statistics.
- 2026-09-04: All six companion Plot adapters, every closed view, standalone-code parity, typed one-shot Worker
  round trips, input immutability, global-state restoration, and method-specific failure paths are GREEN. Domain
  regressions pass (180 producer/codec/worker tests plus 87 composition/sample/plot tests, with one optional skip).
- 2026-09-04: Final review removed renderer-side log-CPM/SVD/correlation from Pseudobulk QC, leaving only retained
  profile QC evidence, and added producer-bound content fingerprints plus valid-range tamper tests for edgeR,
  PyDESeq2, and scVI DE tables. Genuine tascCODA scope-specific nullable fields now pass all four plot views.

## Answer

Added separate read-only Plot adapters for Pseudobulk QC, frequentist Pseudobulk Condition contrasts, cell/model-
conditional scVI population evidence, Milo neighborhood differential abundance, scCODA relative composition, and
tascCODA hierarchy-aware relative composition. Each consumes only its exact stored result, validates producer,
schema, axes and available fingerprints before rendering, emits `plot`, strict `summary`, and standalone equivalent
`code`, and keeps independent-Sample inference, posterior evidence, overlapping-neighborhood evidence, and
compositional reference effects explicitly distinct.
