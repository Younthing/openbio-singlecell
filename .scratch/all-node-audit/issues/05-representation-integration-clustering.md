# Representation, Technical batch integration, graph, and clustering

Type: task
Status: resolved

Blocked by: 04

## Initial audit focus

Representation dimensionality, actual neighbor parameters, optional dependencies, integration diagnostics without overclaiming, graph prerequisites, and clustering result disclosure.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: Pre-change research is complete for all ten nodes. Each node has exactly
  `official-usage.md` and `design-review.md` under `.scratch/all-node-audit/research/<node-id>/`, and the
  primary agent read all twenty documents before implementation began.
- 2026-08-28: The reviewed seams keep PCA, Harmony, scVI, Neighbors, UMAP, t-SNE, force-directed graph,
  single-resolution Leiden, resolution sweep, and PCA metadata diagnostics as separate atomic analyses.
  Shared internals are limited to representation/graph validation and the fixed Leiden adapter.
- 2026-08-28: Confirmed release-blocking defaults: Neighbors/t-SNE can request more dimensions than a
  short learned latent representation supplies; the Best Practice scVI graph requests 30 dimensions from
  the default 10-dimensional model; Force-Directed Graph defaults to unavailable PAGA initialization; and
  its requested ForceAtlas2 layout silently falls back to Fruchterman-Reingold when `fa2-modified` is absent.
- 2026-08-28: PCA Metadata Associations will replace the removed decoupler 1.x API with an explicit
  Sample-level SciPy implementation so cells are not treated as independent biological replicates.
- 2026-08-28: All ten nodes were refactored with atomic result, strict JSON `summary`, and equivalent
  `code` contracts. Final hardening added global/atomic legacy-workflow preflight, collision-safe Leiden
  output keys, and an unambiguous Harmony orientation contract, including square embeddings.
- 2026-08-28: Release verification passed 294 targeted Python tests under `-W error` (the combined
  Batch05/Batch06 suite) and all 40 workflow-migration tests. Independent review findings were resolved
  before this issue was closed.
- 2026-08-28: A final release-dependency audit found that the implementation had combined the
  `harmonypy` 2.0.0 call signature with the historical 0.x `Z_corr` transpose convention. The two
  Harmony research documents were corrected before code changes; runtime and generated code now require
  exact `harmonypy==2.0.0` and accept only its documented cell-by-component output orientation. The
  focused integration suite passes 36 tests under `-W error`, and Ruff reports no findings.

## Answer

The ten representation, integration, graph, clustering, and metadata-diagnostic nodes remain separate
because each represents one scientific or computational decision. Their shared code is limited to strict
representation/graph validation, report construction, and deterministic adapters. Parameters that change
the scientific result are exposed; implementation mechanics and validated defaults remain internal.

Every modified node has the required `official-usage.md` and `design-review.md` written and read before
implementation. Runtime reports disclose method assumptions, package/software versions, key quantitative
results, warnings, and references without claiming biological replication from cells. Generated source is
tested against the runtime behavior, and legacy workflow bundles are migrated only after a complete,
zero-mutation preflight.
