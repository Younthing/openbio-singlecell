# Diffusion Pseudotime: module design review

## Decision

**Retain and deeply refactor `OpenBioSingleCellDPT`.** It performs exactly one operation: orient an already validated diffusion geometry from one explicit biological root. Keep Diffusion Map upstream and PAGA separate.

## Target interface

Inputs:

- `adata`;
- `neighbors_key="neighbors"`, visible;
- `root_mode`, bounded combo `cell_id|group_medoid`, default `cell_id`;
- `root_cell_id`, used by `cell_id` and matched exactly to unique `obs_names`;
- `root_column` and `root_value`, used by `group_medoid`, with exact typed label matching;
- `n_dcs=10`, visible;
- `overwrite_existing=False`, advanced.

For `group_medoid`, select the nominated population's cell nearest its centroid over informative diffusion components `1:n_dcs`, breaking exact ties by lexical observation ID. This is deterministic and reportable; it does not pretend that the algorithm discovered the biological root. Reject an empty root population. The default `cell_id` forces the most scientifically important decision to be explicit.

Outputs are `adata`, `summary`, and `code`. Fix `n_branchings=0`; branch detection is PAGA's responsibility. Use canonical `dpt_pseudotime`/`iroot` keys because Scanpy 1.12.3 does not expose an output namespace, and guard collisions.

## State and postconditions

Validate the complete named graph and require one connected component for this supported single-root contract. Require a valid Diffusion Map bundle whose recorded graph and observation fingerprints match the requested graph/current observations. Validate `n_dcs` against the stored basis and finite values. Never generate Diffusion Map implicitly.

After Scanpy, require one finite pseudotime per cell, values within `[0, 1]` up to tolerance, selected root value zero up to tolerance, no input mutation, and unchanged graph/diffusion bundle. Store root ID/mode/rationale and graph/diffusion fingerprints in provenance.

`summary` reports exact root ID and population/medoid rule, graph and diffusion identity, n_dcs, pseudotime min/quantiles/max, counts near start/end, and bounded summaries by an optional annotation only if explicitly selected. It states that pseudotime is relative/exploratory and no branch, measured-time, `Sample`, or `Condition` inference occurred.

`code` defines a standalone function reproducing exact root selection, prerequisites, collision checks, `uns["iroot"]`, public `sc.tl.dpt`, and postconditions; it returns `(output_adata, summary_dict)`.

## Migration and tests

Old `root_column/root_value` workflows migrate to `group_medoid`, never to old first-row behavior, and receive a result-changing migration notice. Insert Diffusion Map immediately upstream when absent and bind it to the same `neighbors_key`. If graph connectivity is unknown, migration is allowed but execution remains fail-closed.

Tests cover missing/stale/malformed diffusion state, graph-fingerprint mismatch, exact ID roots, typed group labels and conversion collisions, row-order-invariant medoid selection/ties, disconnected graphs, component bounds, collisions, known pseudotime properties, input preservation, strict JSON, compiled generated code, and runtime/code equivalence.

## 2026-08-28 adversarial implementation audit

The first implementation was not release-ready. Its inherited graph fingerprint could accept a stale Diffusion Map computed from a different actual matrix, and range/root checks allowed arbitrary signature-compatible pseudotime to receive trusted provenance. Although the private core distinguished numeric and string population labels, the public node exposed only a string `root_value`; numeric categoricals were therefore unreachable. The medoid centroid and exact-float tie check also depended on observation row order in cancellation cases.

Release closure requires independent recomputation of the Scanpy 1.12.3 DPT distance row and normalized ordering, an explicit bounded label-type selector with strict parsing, cell-ID-sorted stable centroid accumulation and tolerance-aware lexical tie breaking, a DPT-only standalone dependency closure, and adversarial tests through the public node seam. Legacy roots must migrate to `group_medoid` with an explicit result-changing note and an inserted or proven same-graph Diffusion Map; unsafe graph-boundary cases fail closed.

## 2026-08-28 implemented closure

DPT now verifies the upstream diffusion eigenpairs independently and reconstructs Scanpy 1.12.3's root distance exactly, including the `0.9994` stationary-state branch, `lambda/(1-lambda)` scaling, square-root distance, and normalization by the finite maximum. Backend pseudotime must equal that independent vector within disclosed float32 tolerances. The report fingerprints a stable `(pseudotime, lexical cell ID)` ordering without serializing the entire cell axis.

Before accepting that vector, DPT also requires the post-backend stored graph matrices to equal the exact canonical matrices consumed, rather than merely reproducing their tolerance-averaged fingerprint.

The public `group_medoid` branch now requires `root_value_type=string|integer|number`; text is parsed strictly and finite numeric categoricals are reachable without string coercion. Candidate rows are sorted by cell ID, centroids and squared distances use stable `math.fsum`, and values within a disclosed relative/absolute tolerance tie are resolved lexically. NumPy state is restored and generated code contains only the DPT dependency closure.

Whole-graph migration rewrites the old string root to explicit string-typed `group_medoid`, preserves and renames connected/direct/proxied root controls, and either verifies a directly upstream same-key Diffusion Map with enough components or inserts one immediately on the same AnnData edge. It uses one collision-free allocator across root/subgraphs and object/array links, preserves every DPT slot-zero consumer, adds a result-changing review note, and returns zero on a second pass. Unconnected, boundary/reroute, graph-key/component conflicts, broken links, and partial/mixed serialization reject the complete workflow before mutation.

An inserted prerequisite uses `overwrite_existing=true` deliberately on its output copy so a bare or indirectly inherited legacy diffusion basis cannot collide with or masquerade as the new graph-bound provenance. Direct current/legacy Diffusion Map producers with the same key and enough components are reused instead. DPT-to-DPT chains are rejected during preflight because their already oriented state cannot be a valid new-basis prerequisite.
