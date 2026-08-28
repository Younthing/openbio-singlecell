# cNMF Consensus Programs: module design review

## Decision

Retain `OpenBioSingleCellCNMF` as the second atomic staged node, preserve its current schema and typed `OPENBIO_CNMF_RUN` input, and migrate its backend to the method authors' file-backed `cnmf==1.7.1`. Do not merge it back into Rank Survey and do not restore hard cell classification.

Inputs remain `run`, `selected_k`, `density_threshold`, `local_neighborhood_size`, `n_top_genes`, and `overwrite_existing`. Outputs remain `adata`, `summary`, and `code`.

## Cohesive responsibility

The node performs one decision-dependent analysis: final consensus and aligned gene-program annotation for one surveyed K/threshold. It does not validate arbitrary count sources, prepare or factorize again, choose K automatically, persist a user path, plot, or assign hard cell labels.

This retains leverage: the expensive survey state can support multiple downstream attempts, while the consensus node's interface exposes only analyst choices that can change the scientific result.

## Exact execution design

Under `run.locked_backend()` the node:

1. verifies concrete-or-explicit-portable artifact identity, class-owned marker/schema, construction-bound complete metadata/metric snapshot, liveness, private-root identity, selected-K membership, backend-required parameter domains, and overwrite preflight;
2. imports and guards exact `cnmf==1.7.1` and the public consensus/load signatures;
3. validates the selected K merged-spectrum and normalized-count files through the contained-path adapter;
4. removes only the exact contained K density cache so neighborhood changes cannot reuse stale upstream state;
5. independently computes local densities and KMeans-label postconditions from the merged spectrum;
6. calls official final `consensus` with all controls explicit and `build_ref=False`;
7. calls official `load_results`, unpacks its exact four-tuple, and validates identifiers, shapes, values, and program family;
8. validates official result files at their expected contained paths;
9. builds the complete strict report, including direct source/source-feature/current-feature context, and applies all canonical writes to one source-aligned AnnData copy;
10. returns without closing the reusable run.

No OmicVerse in-memory fields are part of the adapter.

## Density/cache seam

In 1.7.1 the density cache key does not encode `local_neighborhood_size`. Cache handling belongs inside the module rather than becoming a public `clear_cache` parameter. The adapter may unlink only the exact expected density-cache file after containment, symlink, and regular-file checks. It never recursively deletes a backend-selected path.

The independent density/label implementation is intentionally narrow and pinned to the audited algorithm. It validates that the requested fraction gives at least one neighbor, retained components are sufficient, labels cover `1..K`, and all diagnostics are finite/aligned. Its version and role are disclosed in summary and generated code.

## Transactional output contract

All collision checks occur before the external final call where possible. All AnnData mutation occurs on a fresh copy only after every external result passes validation. An `X` or layer run uses the original current AnnData feature axis. A Raw run uses a Survey-owned materialization of the complete `adata.raw` axis and records the Survey-time current-feature count separately; Consensus therefore returns a Raw-axis AnnData without requiring Raw/current equality or claiming Raw history/integrity. Canonical storage is limited to:

- `obsm["X_cnmf_usage"]`;
- `varm["cnmf_gep_scores"]`;
- `varm["cnmf_gep_tpm"]`;
- `uns["cnmf"]`.

The module never sprays program columns into `obs`/`var`, never writes `cNMF_cluster`, and never copies broad backend namespaces. Program names and axes are explicit and identical across matrices. Existing keys fail closed unless overwrite is explicit.

## Run ownership and failure behavior

The run owns its private directory until explicit close or finalization. Consensus serializes cache mutation with the run lock but does not assume ownership of cleanup. Metadata and metrics are getter-only views bound to a private complete snapshot; the local concrete run or an explicitly class-marked generated run must validate that snapshot before use. This portable seam exists only because generated Survey and Consensus source may run in separate namespaces; it is not generic duck typing.

Equivalent Python exposes `close()`/context-manager cleanup; in a graph, finalization after cache release is the reachable cleanup mechanism because automatic close would invalidate sibling Consensus consumers. A cached run may retain its directory until eviction or process exit and this limitation is reported. A cleanup exception leaves the object closed while retaining the temporary owner and attached finalizer for fallback cleanup. If cleanup is secondary to an existing construction, context-body, or node-report failure, the primary exception remains outward and the cleanup failure is attached as a note. Explicit `close()` cleanup failure remains outward and retryable. All backend paths are revalidated at use time, providing replacement/substitution evidence between Survey and Consensus, not cryptographic authenticity against a hostile writer. Path escape, symlink/reparse substitution, nonregular files, replacement during read, missing restart family, metadata/metric snapshot mismatch, or closed state fail before output mutation.

Known upstream warnings are matched narrowly by exact category/message/source. Unknown warnings are re-emitted. An external exception does not close a still-valid upstream run and does not partially annotate input AnnData.

## Summary and generated-code parity

The report inherits the entire Survey context—source state/evidence/advisories, candidate metrics, fixed policy, expected/completed jobs, input fingerprint, and resource envelope—then adds the one final decision, filter accounting, result distributions, leading genes, canonical keys, versions, references, and limitations. Its key results directly expose `source`, `source_features`, and `current_features`: `X`/layer output is current-axis aligned, while Raw output is aligned to the complete Raw source axis and may have a different feature count from the current AnnData. The typed run snapshot carries both counts so Consensus does not reconstruct current context from mutable input or indirect evidence. This makes the staged artifact self-describing without coupling the node to mutable backend defaults.

Standalone code uses the same owned-run class and public official calls. Runtime/generated parity covers all three matrices, source-aligned output axes, direct source/source-feature/current-feature reporting (including Raw Survey through Consensus), canonical scientific `uns["cnmf"]`, density/labels/filter accounting, strict JSON, and lifecycle behavior. Per-execution managed-file SHA manifests and absolute-layout fingerprints are independently verified before/after within each run but are not expected to be byte-identical between two separately created temporary artifacts; plugin-only execution history is also outside parity.

## Migration matrix

The existing migration remains correct and must not change for this backend-only replacement:

- legacy `adata` plus source/K-range/restarts/HVG/seed become Rank Survey inputs;
- the inserted typed `run` edge feeds this node;
- legacy `selected_k`, `density_threshold`, `local_neighborhood_size`, `n_top_genes`, and overwrite intent remain mapped here;
- removed working-directory/name/worker/GPU/classifier controls remain removed or fail closed when connected/exposed;
- unsupported legacy raw selection remains visibly flagged;
- downstream legacy AnnData consumers remain connected to this node's AnnData output;
- already-current staged schemas are migration no-ops.

Public schemas and typed semantics are unchanged, so no workflow migration edit is required.

## Verification plan

- fake official file backend asserting exact final calls and tuple return orientation;
- real tiny exact-1.7.1 Survey-to-Consensus smoke;
- selected-K, threshold, neighborhood, top-gene, overwrite, and closed-run boundaries;
- repeated consensus calls with different neighborhood fractions proving stale cache isolation;
- independent density/label parity and adversarial missing/nonfinite/misaligned labels/results;
- all result-axis permutations and identifier mismatch/duplicate rejection;
- positive Raw Survey-to-Consensus runtime/generated parity on the complete Raw feature axis, with direct `source`, `source_features`, and `current_features` summary assertions;
- canonical storage, input immutability, transactional failure, and no hard clusters;
- exact concrete/portable artifact acceptance, fake duck/instance-marker rejection, and metadata/metric rebinding or low-level mutation rejection;
- live directory during use, explicit cleanup, finalizer cleanup, retryable cleanup failure, primary-exception preservation, and path replacement rejection;
- exact known-warning isolation and unknown-warning propagation under `-W error`;
- strict summary and standalone generated full-state parity;
- registration plus unchanged legacy migration tests.

## Final assessment

The node should be retained. It is a deep, high-cohesion module around one final consensus decision and hides a comparatively complex file/cache/backend boundary. Merging it with Survey would lose reuse and couple analyst choice to expensive computation; splitting storage or tuple loading into public nodes would expose accidental backend detail.
