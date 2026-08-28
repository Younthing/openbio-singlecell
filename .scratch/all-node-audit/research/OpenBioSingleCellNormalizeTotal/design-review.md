# Normalize Total: module design review

## Current module

The module copies `AnnData`, silently creates `layers["counts"]` from the pre-normalized `X` when that key is absent, and then applies `scanpy.pp.normalize_total` to `X`. Its visible interface contains only `target_sum`, so callers cannot declare which expression state is being treated as counts.

## Correctness and interface problems

- Implicit count-layer creation gives this module ownership of a Raw snapshot decision that belongs to `OpenBioSingleCellSnapshotExpression`.
- A normalized, logged, scaled, negative, non-finite, all-zero, or partly zero-total matrix can enter without a local scientific invariant.
- The interface cannot select an existing count layer without an upstream state-switch node, making the scientific source depend on workflow ordering.
- `target_sum` is constrained by the UI but direct execution is not validated for positivity or finiteness.
- No report exposes input depth, achieved output totals, non-integer input, limitations, citations, or software versions.
- No equivalent source is returned, and callers can currently rely unknowingly on the undocumented `counts` side effect.

## Decision: keep and enhance as an atomic normalization module

Keep the module because total-count scaling is a coherent, independently useful transformation. Do not merge it with `Log1p`: users can validly inspect or consume linear depth-normalized values, and keeping the numerical steps distinct makes repeated or omitted transforms visible. Do not merge it with `SnapshotExpression`: a Raw snapshot is an auditable expression-state checkpoint, while normalization is a reproducible derived representation.

The deletion test supports retention. Removing the module would reproduce source resolution, count validation, sparse-safe totals, Scanpy parameter handling, reporting, and code generation in every caller that needs a normalized `X`. Those behaviors provide leverage behind a small interface.

## Target interface and parameter policy

- Inputs: `adata`; explicit expression `source` restricted to aligned `X` or a named layer, default `X`; visible positive `target_sum` defaulting to 10,000.
- Optional advanced scientific controls: `exclude_highly_expressed=False` and its conditional `max_fraction=0.05`, if the implementation elects to support the documented Scanpy alternative. They must be exposed together, never inferred from the data.
- Outputs: copied `adata`, structured `summary`, equivalent Python `code`, in that order.
- Fixed hidden policy: normalize into output `X`; copy the input; do not mutate the selected source layer; do not create/overwrite `raw` or any layer; use natural storage dtype behavior from Scanpy; do not expose `inplace`, `copy`, `key_added`, `obsm`, threads, or chunk controls.

The expression source is visible because it determines the scientific state used as counts. `target_sum` is visible because it defines the numeric result. Copy/storage controls remain hidden because varying them increases coupling without changing the scientific question. `raw` is deliberately not a selectable source: a Raw snapshot can retain a wider variable axis than current `X`, so assigning its normalized matrix to `X` would not have a stable aligned-output contract.

For source compatibility, append the new source/advanced inputs rather than reordering the existing `adata, target_sum` prefix. Existing saved workflows that depended on this node to manufacture `layers["counts"]` require an explicit Snapshot node before normalization; this semantic migration should be stated in release notes rather than emulated.

## Invariants and failure modes

- Require an in-memory, nonempty AnnData and an aligned selected matrix.
- Require a finite numeric aligned matrix. Treat non-integer input, negative values/totals, no positive values, and zero-total cells as expert-visible warnings; Scanpy leaves zero-total rows unchanged rather than making execution impossible.
- Accept non-integer corrected counts only with a warning and explicit limitation.
- Require a finite positive target; validate the high-expression threshold when enabled.
- Preserve cell/feature identity and order, annotations, existing layers, and existing `raw`; do not mutate the caller.
- Proven logged/scaled state and unknown provenance both proceed with differentiated warnings. Provenance is evidence for interpretation, not a numeric execution gate, and neither case may be reported as proven raw UMI input.

## `summary` and `code` contract

The report should disclose the source, dimensions, sparse/dense type and dtype, integer-like status, input per-cell total distribution, requested target, exclusion policy, achieved output-total distribution, count of zero-total cells, and whether existing canonical count/raw states were preserved. Its report-ready result must distinguish library-depth scaling from Technical batch integration or Condition inference. Capture Scanpy/AnnData/NumPy/SciPy/Pandas/Python/openbio-singlecell versions dynamically and include method/software references.

Generated code must implement the same source resolution, validation, copy, normalization, and preservation contract. It should be executable with public scientific packages and reproduce the primary AnnData state except plugin-only analysis history.

## Cohesion, coupling, and seam placement

The external seam is the node interface; Scanpy and AnnData are in-process dependencies, so no adapter is justified. Source resolution, count validation, sparse-safe summaries, normalization, and provenance stay local to the implementation. A shared private normalization helper may be used by `NormalizeToLayer`, but it is an internal seam and must not become another user-facing state-switch module.

This shape has high cohesion: one call maps one declared count representation to one normalized active representation. Removing count snapshot creation also lowers coupling because only `SnapshotExpression` owns the canonical `counts`/`raw` convention.

## Verification plan

- Dense and sparse `X`/layer inputs reach the requested total without mutating the caller or selected source layer.
- Existing `layers["counts"]` and `raw` are preserved byte/value-equivalently; absent ones remain absent.
- Empty axes, missing/misaligned/non-numeric/non-finite matrices, invalid targets, and backed objects fail clearly. Negative/all-zero/zero-total or proven transformed sources execute with warnings, unchanged zero rows, exact disclosure, and finite-output postvalidation.
- Non-integer corrected counts produce a report warning without being mislabeled as raw UMI counts.
- High-expression exclusion, if exposed, matches Scanpy and reports its different total interpretation.
- Strict-JSON `summary` and compiling generated `code` reproduce the normalized `X` and all scientific input/error invariants.
