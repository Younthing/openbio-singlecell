# Scale: module design review

## Current module

The node copies `AnnData` and calls `scanpy.pp.scale(output, max_value=...)` on `X`. It always uses Scanpy's default `zero_center=True`, translates `max_value <= 0` to no clipping, and returns the transformed copy.

## Correctness and interface problems

1. The source is implicitly `X`, and the scaled values replace the primary expression state in the returned object. A caller cannot scale a named normalized layer while preserving it.
2. `zero_center` is hidden at `True`. Sparse input therefore densifies by default, even though uncentered scaling is the official sparse-preserving option and Scanpy PCA can center sparse input implicitly.
3. `max_value=0` is overloaded as “no clipping,” although the official function treats zero as a valid clipping threshold. The current interface cannot faithfully represent Scanpy's parameter.
4. There is no dense-allocation estimate, ceiling, output-density disclosure, or early failure for matrices that exceed available memory.
5. Generic `.var["mean"]` and `.var["std"]` are overwritten, with no source/destination key tying them to a particular scaling run.
6. Non-finite input, fewer than two observations (`ddof=1`), backed data, duplicate axes, destination conflicts, and unsupported matrix types are not validated explicitly.
7. The node emits neither a structured `summary` nor equivalent `code`.

## Decision: keep, redefine as layer-preserving, and enhance

Keep Scale because centered z-scaling and sparse-preserving variance scaling are independently useful transformations. Redefine the module to read one explicit current-axis expression source and write one explicit output layer without replacing the source. Do not merge it with normalization, regression, HVG selection, or PCA: each is a distinct scientific decision, and PCA can often avoid the densification that explicit centering causes.

The default clustering workflow is correct to omit Scale when PCA can consume normalized/logarithmized sparse data directly. The node should be selected only when a downstream method or visualization genuinely needs a materialized scaled matrix.

## Target interface

- Visible inputs: `adata`, explicit `X`/named-layer expression source, `zero_center`, and clipping mode (`none` or `custom`).
- Conditional scientific input: non-negative custom `max_value`.
- Advanced technical inputs: output layer, overwrite-existing opt-in, and dense-allocation ceiling in GiB.
- Outputs: copied `adata` with a scaled destination layer, structured `summary`, equivalent Python `code`.

`zero_center` stays visible because it changes both the statistic and the storage complexity. Clipping remains exposed because it changes extremes used by downstream analyses. `copy=False`, `mask_obs=None`, public layer routing, `ddof=1`, finite checks, memory estimation, and result summarization are hidden implementation policy. Do not expose `obsm` or `mask_obs` on this gene-expression node; observation-mask fitting is a different cohort-definition question, and arbitrary `obsm` scaling would broaden the interface beyond its cohesive purpose.

The source and output layer names are storage controls, but source remains visible because it identifies the expression state. Output naming and overwrite are advanced. `raw` is excluded because it can have a different variable axis after feature selection and should remain the full-gene Raw snapshot.

## Validation and state contract

- Require an in-memory `AnnData` with at least two observations, at least one feature, unique axes, and a finite numeric selected matrix.
- Accept dense and SciPy CSR/CSC sources. Reject backed/lazy sources with an actionable `adata.to_memory()` instruction unless a supported chunked implementation is deliberately added later.
- Treat clipping “none” as `max_value=None`; require custom `max_value >= 0` and preserve zero as a real threshold.
- Compute/report feature means and `ddof=1` standard deviations from the selected source. Retain constant features following the pinned Scanpy behavior and report their count.
- For sparse input with `zero_center=True`, estimate one dense array conservatively at eight bytes per entry, warn that peak use is higher, and fail before copying/allocation when the configured ceiling is exceeded. For `zero_center=False`, preserve sparse storage.
- Reject blank output names, a destination identical to the selected source layer, and existing destinations unless overwrite is explicitly enabled. Never modify the selected source, `X`, other layers, Raw snapshot, observations, or unrelated feature annotations.
- Store scale parameters and feature-statistic provenance keyed by destination rather than silently reusing ambiguous generic columns. If compatibility requires Scanpy's generic columns, the report must identify and overwrite them only under explicit policy.

## Report contract

`summary` includes:

- methods text that distinguishes centered z-scaling from uncentered variance scaling and states `ddof=1` and clipping bounds;
- results text suitable for a preprocessing report without describing scaling as batch correction or inference;
- key results for dimensions/source, input/output sparse status and dtype, feature mean/std summaries, constant features, effective clipping interval, clipped entries, output range, and dense-memory lower bound;
- all resolved parameters, warnings for densification/uncentered clipping/expression-state limitations, and the pinned constant-feature behavior;
- Scanpy and AnnData references plus dynamic Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy versions.

`code` is a self-contained function that copies the input, validates the same shape/finite/destination/memory invariants, creates an isolated working AnnData or matrix from the selected source, invokes public `scanpy.pp.scale` with explicit `zero_center` and `max_value`, stores only the scaled result and namespaced provenance in the destination, and returns an equivalent `AnnData`. It must reproduce sparse preservation for `zero_center=False` and dense output for centered sparse input without plugin history.

## Cohesion and coupling assessment

The redesigned node is a deep in-process module: its small interface hides source routing, Scanpy's `ddof=1` and constant-feature behavior, clipping semantics, sparse/dense transitions, resource safeguards, namespaced provenance, reporting, and version capture. Merging it into PCA would reduce leverage because explicit scaling is optional and often less memory-efficient than PCA's implicit centering; merging it with normalization would conflate depth normalization/variance stabilization with feature standardization.

## Verification plan

- Dense, CSR, and CSC sources from `X` and a named layer match direct Scanpy results while leaving the input/source unchanged.
- With clipping disabled, `zero_center=True` produces zero-centered unit-variance nonconstant features; `zero_center=False` preserves sparse storage and does not claim zero means.
- No clipping, custom zero, and custom positive clipping reproduce Scanpy's centered symmetric and uncentered upper-only behavior.
- Constant genes, integer input casts, small two-cell data, non-finite values, one-cell data, empty features, duplicate axes, backed data, and unsupported types follow the contract.
- Sparse centering reports the allocation and fails under a deliberately low ceiling before allocating; uncentered scaling does not trigger that failure.
- Blank/same-source/existing destinations obey overwrite policy, source matrices and Raw remain intact, and scale provenance is destination-specific.
- `summary` is strict JSON and `code` compiles and reproduces primary matrices, storage type, annotations, and error invariants.

## Final audit amendment

The normalized-log layer suggestion must be the selected interface default. In the layer-preserving preprocessing design, leaving `X` selected by default can scale retained counts while merely warning after the scientific choice has already been made; selecting `log1p_norm` makes the safe common path explicit while preserving `X` as an available deliberate source. Legacy direct Python calls that omit the dynamic source may fall back to `X` only when the suggested layer is absent; normal UI calls carry their explicit selection.

The open-expert boundary audit retains only finite/aligned matrix checks, unique axes, the one-observation `ddof=1` backend prerequisite, destination integrity, and the explicit dense-memory budget as hard gates. Expression-state provenance, common clipping values, and whether scaling is advisable are interpretation warnings and summary limitations; no further runtime change is required.
