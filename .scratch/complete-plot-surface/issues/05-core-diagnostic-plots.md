# Add correction, factorization, clustering, and diagnostic plots

Type: task
Status: resolved
Blocked by: 01

## Scope

- MAD Outlier, Scrublet, cNMF Rank/Programs, PCA Loadings, Neighbor Graph, Leiden Sweep, Schist Hierarchy, and PCA Metadata Association plots.

## Acceptance criteria

- Every Plot consumes stored evidence, defaults to the primary diagnostic view, and does not select a scientific optimum.

## Comments

- 2026-09-04: Claimed; implementation starts after the shared convention is frozen.
- 2026-09-04: MAD/Scrublet, cNMF Rank/Programs, PCA Loadings, Neighbor Graph, Leiden Sweep, and Schist Hierarchy
  are implemented and focused regressions pass. PCA Metadata Association remains in progress under the exclusive
  results-module owner.
- 2026-09-04: PCA Metadata Associations Plot is implemented under the exclusive results-module owner with closed
  association-heatmap/effect-size views, canonical producer/schema/content-fingerprint validation, standalone
  code parity, one-shot Worker input preservation, and 9 focused tests GREEN. Ticket-wide resolution remains with
  the coordinating owner/integration pass.
- 2026-09-04: Representative heatmap QA added an axis-external significance legend bound to the retained producer
  alpha (`BH-adjusted p <= alpha`); focused runtime/standalone parity remains 9 passed and Ruff is GREEN. The
  re-rendered figure passed visual QA with no overlap or clipping.
- 2026-09-04: Resolved after aggregate producer/codec/plot regression, resource-bound audit, and representative
  visual QA. The completed first-wave cross-domain run passed 245 tests with one optional-backend skip.

## Answer

- MAD Outlier and Scrublet Diagnostics retain and validate their exact threshold, simulated-score, Sample/group,
  and observation-axis evidence, with closed distribution/rate views.
- cNMF Rank and Programs plots read retained restart/stability/error, program-usage, and top-gene evidence without
  selecting K or changing the consensus result.
- PCA Loadings, Neighbor Graph Diagnostics, Leiden Resolution Sweep, Schist Hierarchy, and PCA Metadata
  Associations render their stored axes and metrics without rerunning PCA, graphs, clustering, or association
  testing. Static item/pixel limits and sparse-memory preflight bound rendering work.
- Every adapter has `plot`, strict `summary`, equivalent `code`, a registered one-shot Worker operation, immutable
  input tests, and explicit limitations that leave component, resolution, and hierarchy choices to the analyst.
