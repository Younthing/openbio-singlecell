# Retain Augur predictions and complete its Plot adapter

Type: task
Status: resolved

## Scope

- Extend AugurResult with fold-level labels/scores or deterministic ROC points.
- Add priorities, cross-validation, feature-importance, and ROC views.

## Acceptance criteria

- Prediction evidence is tied to population/subsample/fold and validates against existing AUC evidence.
- Existing typed artifact identity and immutability remain intact.

## Comments

- 2026-09-04: Claimed as an upstream-contract prerequisite and its direct consumer.
- 2026-09-04: Resolved through artifact, fitted-estimator extraction, codec, and four plot-view red-green slices;
  60 Augur tests plus 38 related artifact/population tests and Ruff passed.

## Answer

- AugurResult/codec v2 retains canonical fold-level predictions keyed by population, subsample, fold, and
  observation, with binary labels and finite probability scores. Validation recomputes each fold AUC against the
  retained cross-validation evidence and rejects coverage, range, identity, and fingerprint tampering.
- The Pertpy adapter reads predictions from each already-fitted fold estimator without retraining. The portable
  artifact identity remains `OPENBIO_AUGUR_RESULT`.
- Augur Plot provides priorities, cross-validation, feature-importance, and ROC views through the existing
  one-shot PlotResult path, with read-only input handling, bounded static rendering, strict summary, and equivalent
  standalone code. The generated code embeds the portable validator and renderer, imports no OpenBio package, and
  preserves artifact-tamper and resource-limit checks.
