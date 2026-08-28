# Normalize to Layer: module design review

## Current module

The module already selects `X` or a named layer, copies AnnData, normalizes a private working object, optionally applies log1p or square root, and stores the result under a user-selected layer without replacing `X`. This is the normalization path used by all packaged production workflows.

## Correctness and interface problems

- Count suitability is not validated: negative, non-finite, transformed, all-zero, or zero-total-cell sources can enter.
- `target_sum=0` is silently converted to Scanpy's dataset-dependent `None`; direct NaN/negative input is not guarded.
- A direct-call transform outside `none/log1p/sqrt` silently behaves like `none` instead of failing.
- The output key may be blank only after a string call, may equal the selected source, may be `counts`, and silently overwrites any existing layer.
- The implementation copies full `obs` and `var` into a temporary AnnData even though the numerical operation needs only the aligned matrix.
- No report distinguishes intermediate normalized totals from post-transform layer values or discloses expression state, collisions, references, limitations, and versions.
- No equivalent source is returned.

## Decision: keep, deepen, and do not merge into consumers

Keep this module as the canonical constructor of a reusable derived normalized layer. The deletion test is decisive: packaged workflows route the same `log1p_norm` layer to HVG selection, PCA, marker evidence, scoring, and visualization. Merging normalization into each consumer would duplicate count validation and computation, allow branch-specific preprocessing drift, and make it harder to report exactly which representation was shared.

Although the implementation performs depth scaling and an optional elementary transform internally, its external outcome is atomic: materialize one declared derived expression representation while preserving all source states. That gives callers substantial leverage behind one interface. Do not merge it with `SnapshotExpression`, whose independently auditable responsibility is to establish the canonical Raw snapshot. Do not make downstream consumers responsible for creating layers; they should only declare which representation they read.

Retain the separate `NormalizeTotal` and `Log1p` modules for workflows that intentionally transform active `X`. Their output contract differs from this layer-producing module. Centralize shared numerical validation and rendering behind an internal seam rather than exposing another public adapter.

## Target interface and parameter policy

- Inputs: `adata`; visible dynamic `source` restricted to aligned `X` or named layer, default `X`; visible positive `target_sum=10000`; visible `transform` in `none/log1p/sqrt`; advanced `output_layer="log1p_norm"`; advanced `overwrite_existing=False`.
- Outputs: copied `adata`, structured `summary`, equivalent Python `code`, in that order.
- Fixed hidden policy: `exclude_highly_expressed=False`, natural-log base, copy-on-write, canonical pseudocount one, no `raw` source, and no mutation of `X`, `raw`, `counts`, or unrelated layers. `inplace`, `copy`, `key_added`, `obsm`, chunking, and thread controls remain implementation details.

Source, target, and transform are visible because they define the scientific representation. The destination key is advanced because it is workflow storage plumbing. Overwrite is advanced and false by default because replacement is a provenance/safety decision, not routine tuning. The rare high-expression exclusion branch is fixed to the documented Scanpy default here to keep this common recipe interface small; users who need it can use the atomic Normalize Total path and must report that choice there.

Preserve the existing positional prefix `adata, source, target_sum, transform, output_layer` and append overwrite control for saved-workflow compatibility. Existing calls with `target_sum=0` should receive a migration error, not retain ambiguous sentinel behavior.

## Invariants and failure modes

- Require an in-memory, nonempty AnnData and a selected aligned `X`/layer matrix.
- Require a finite numeric aligned source. Warn on non-integer, signed, all-zero/zero-total, or provenance-inconsistent input; then hard-check the selected transform's numeric domain (`sqrt`: intermediate values >= 0; `log1p`: intermediate values > -1).
- Require a finite positive target and an explicitly supported transform.
- Require a nonempty destination different from the selected source; reserve `counts` for `SnapshotExpression`.
- Reject existing destination layers unless overwrite is explicit; report replacements.
- Preserve current `X`, existing `raw`, canonical counts, annotations, identity/order, and all unrelated layers; never mutate the caller.
- Proven normalized/logged/scaled state and unknown state both proceed with differentiated warnings; provenance is evidence, not an execution veto. Never describe either as proven raw counts.

## `summary` and `code` contract

The report must state the source and destination, cells/features, sparse/dense storage and dtype, integer-like status, input depth distribution, target and fixed high-expression policy, intermediate normalized-depth distribution, selected transform/formula, final value distribution, overwrite behavior, and preservation of `X`/Raw snapshot/counts. The core result sentence should say that a named derived layer was created for a stated number of cells and genes; it must not call that layer counts or claim Technical batch correction. Include relevant Scanpy, AnnData, and transformation references plus dynamic runtime versions.

Generated code should resolve the source, validate count state and collisions, copy AnnData, compute normalization and the selected transform on a private working matrix/object, and assign only the destination layer. It must reproduce errors and warnings that affect scientific interpretation; plugin-only history is excluded.

## Cohesion, coupling, and seam placement

The external seam hides source adaptation, sparse-safe validation, two-step numerical composition, collision safety, and provenance/reporting. Scanpy, AnnData, NumPy, and SciPy are in-process dependencies, so no adapter is warranted. An internal pure helper shared with the atomic nodes improves locality without widening any node interface.

The module is highly cohesive when described as “derive one reusable expression representation.” Keeping this computation before the consumer fan-out reduces coupling: downstream modules name the layer they require but do not know how it was built. Reserving `counts` and `raw` for one Snapshot module eliminates competing ownership of canonical state.

## Verification plan

- Dense and sparse `X`/layer sources produce Scanpy-equivalent `none` and `log1p` results and NumPy-equivalent square-root results.
- Input `X`, source layer, canonical counts, existing `raw`, and unrelated layers remain unchanged; the caller is not mutated.
- Missing source, empty axes, non-numeric/non-finite input, invalid transform domain/target, backed objects, blank/reserved/same destination, and unapproved collisions fail clearly. Signed/all-zero/per-cell-zero and transformed-state inputs execute when the selected transform remains finite, with exact warnings and disclosure.
- Existing destination fails by default and is reproducibly replaced only with explicit overwrite and report disclosure.
- Reports distinguish intermediate normalized totals from final transformed sums and mark the output as non-count derived expression.
- Strict-JSON `summary` and compiling generated `code` reproduce the layer and all scientific error invariants.

## 2026-08-29 generated-code provenance repair

Equivalent source still creates exactly one requested derived layer and preserves all other states, but it no longer
creates or extends OpenBio `analysis_history`. The `summary` remains the reproducibility/provenance artifact for the
node run.
