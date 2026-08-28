# PAGA: module design review

## Decision

**Retain and deepen `OpenBioSingleCellPAGA` as one undirected partition-graph abstraction.** Do not merge it with DPT, PAGA plotting, force-layout initialization, or RNA-velocity PAGA. Those operations have different prerequisites and claims.

Target inputs: `adata`; visible `groupby="leiden"`; visible `neighbors_key="neighbors"`; advanced `overwrite_existing=False`. Fix `model="v1.2"`, `use_rna_velocity=False`, and `copy=False` internally. Target outputs: `adata`, `summary`, `code`.

## Deep implementation

Use the shared named-graph contract. Require a complete categorical `groupby`, no blank/string-collision labels, at least two represented categories, and deterministic removal of globally unused categories on the working copy with disclosure. Reject collisions for both `uns["paga"]` and the official `uns[f"{groupby}_sizes"]` sidecar before calling Scanpy. Verify the returned `uns["paga"]` group key, the sidecar counts, and both sparse matrices. The full connectivity matrix must be square, finite, nonnegative, symmetric, zero-diagonal, and category-aligned. Scanpy 1.12.3 stores the tree with one orientation per edge rather than symmetrically: validate it as a finite nonnegative zero-diagonal subgraph with exact parent-edge weights, no reciprocal duplicates, and an undirected spanning forest after explicit symmetrization.

Store category labels/order, group cell counts, graph/obs fingerprints and output-matrix fingerprints in portable provenance. This named graph artifact remains inside AnnData because existing UMAP/draw-graph consumers use Scanpy's canonical bundle; the explicit provenance is the interface, not a bare `uns` name.

`summary` includes group counts, graph diagnostics, all nonzero group-pair edges in a guarded/bounded table, top edges, matrix shapes/nnz, component sizes, model, references and versions. Limitations explicitly say the weights are connectivity evidence rather than p-values or trajectories.

`code` performs equivalent validation and returns `(output_adata, summary_dict)` after an explicit public `sc.tl.paga` call. It does not plot or set PAGA positions.

## Migration and tests

The old inputs map directly; add non-overwrite behavior. A future directed velocity PAGA requires a separate node and a typed velocity-graph artifact rather than another mode here.

Tests cover categorical order and unused categories, missing/blank/colliding labels, one group, custom/malformed graphs, output collision, known small graph edge weights, matrix/category postconditions, disconnected graph disclosure, input preservation, strict JSON, generated-code compilation, and runtime/code equivalence.

## 2026-08-28 adversarial implementation audit

The first implementation was not release-ready. Structural output checks never used the validated distance graph, so a signature-compatible backend could return an unrelated or all-zero abstraction and receive trusted provenance. Scanpy 1.12.3 also sanitizes every string `obs`/`var` column during `tl.paga`; calling it on the output object changed unrelated state. Changing `groupby` with overwrite left the former `<groupby>_sizes` sidecar, and the permitted 256-group dense summary exceeded five megabytes.

Release closure requires independent model-v1.2 recomputation from the exact stored distance sparsity and categorical codes, public-API execution on a minimal scratch AnnData followed by transfer of only the audited PAGA bundle, complete old-sidecar cleanup, an explicit summary edge budget/truncation disclosure, and a PAGA-only standalone dependency closure. Malicious-backend, unrelated-state, overwrite, report-budget, and migration tests are mandatory.

## 2026-08-28 implemented closure

PAGA now independently reconstructs the exact Scanpy 1.12.3 model-v1.2 calculation from the stored coordinates of the named distance matrix: directed within/between-category edge counts, the observed-to-random-null scaling capped at one, and SciPy's inverse-weight minimum spanning forest. Both backend matrices must equal this recomputation, so a structurally valid zero or arbitrary bundle cannot acquire provenance.

The public `scanpy.tl.paga` call runs on a minimal zero-variable scratch AnnData containing only the selected categorical observation column, audited neighbor metadata, and canonical graph matrices. Only the verified PAGA bundle and current `{groupby}_sizes` are transferred; unrelated `obs`/`var` strings cannot be sanitized, and overwrite removes the former grouping's sizes sidecar. Partition hashing streams the categorical codes instead of materializing per-cell JSON records. Reports retain the exact total edge count but cap serialized edge rows at 2,000 with an explicit truncation flag, limit, and warning. Standalone code is PAGA-only, preserves NumPy process state, and has runtime/code parity. Migration atomically adds non-overwrite/report controls while preserving links/exposures in root and subgraphs.

The scratch graph is also compared byte-canonically after the backend call, closing equal-and-opposite sparse-matrix tampering that would survive an averaged symmetry fingerprint.
