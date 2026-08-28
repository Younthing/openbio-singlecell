# Neighbors: module design review

## Decision

**Keep and deepen `OpenBioSingleCellNeighbors` as the single representation-to-observation-graph module.** Do not merge it with PCA/scVI/Harmony, UMAP, force-directed drawing, or Leiden. Those modules produce or consume the graph through a real AnnData graph seam; deleting Neighbors would redistribute representation selection, key conventions, collision handling, and graph validation across every consumer.

The atomic responsibility is:

```text
one declared cell representation -> one named distance/connectivity graph
```

It must not calculate a missing PCA implicitly, integrate a `Technical batch`, create an embedding, choose clusters, or claim biological neighborhoods.

## Interface and the P0 fix

The current `n_pcs=50` interface leaks PCA terminology into arbitrary `.obsm` matrices and breaks 10D scVI. Replace it with a representation-neutral control:

- `adata`, required;
- `use_rep`, visible, explicit `X` or `.obsm` key;
- `n_dimensions`, visible, default `0` meaning all available dimensions;
- `n_neighbors`, visible, default `15`;
- `metric`, visible bounded combo, default chosen and documented consistently (prefer official Euclidean unless migration policy preserves cosine);
- `method`, advanced combo `umap|gauss|jaccard`, default `umap`;
- `key_added`, advanced, default `neighbors`;
- `overwrite_existing`, advanced, default `False`;
- `random_seed`, advanced, default `0`.

Internally, `n_dimensions=0` becomes `n_pcs=None` for `.obsm`, using the complete latent representation. A positive value slices only after validation and is reported as “dimensions,” never universally as PCs. For old saved workflows, migrate the old `n_pcs` value, but detect values larger than the selected width and convert them to all dimensions with an explicit migration record rather than leaving a runtime failure.

Keep `knn=True`, `transformer=None`, `metric_kwds={}`, `distances=None`, and `copy=False` hidden. `knn=False` changes the graph/memory contract and has no demonstrated workflow need; transformer objects and arbitrary metric kwargs are shallow backend exposure. If a second concrete scalable backend becomes required, add a bounded advanced transformer choice after parity tests rather than an object seam with one adapter.

## Deep implementation behind the seam

The implementation should:

1. validate in-memory nonempty AnnData and unique observations;
2. resolve the representation without Scanpy's implicit PCA fallback;
3. validate axis alignment, numeric finiteness, available dimensions, and metric-specific degeneracy;
4. require `n_neighbors >= 2`; allow values at or above `n_obs`, then disclose Scanpy's effective resolved neighborhood size;
5. resolve output keys and reject collisions unless overwrite is explicit;
6. call `scanpy.pp.neighbors` with every supported/hidden parameter explicit;
7. validate the stored uns/obsp contract after the call;
8. compute graph diagnostics and emit primary AnnData, `summary`, and `code`.

This earns depth: callers learn a small representation/graph interface while backend selection, Scanpy key indirection, sparse validation, collision handling, and diagnostics remain local.

## Summary contract

`summary` must be strict JSON and include:

- methods and references for Scanpy, the selected connectivity kernel, and any non-default metric rationale;
- declared/resolved representation, source shape, available/used dimensions;
- requested/effective `n_neighbors`, metric, method, hidden `knn=True`, transformer policy, and seed;
- uns/distances/connectivities keys and overwrite status;
- sparse matrix shapes/dtypes/nnz, connected components, component-size distribution, isolated-cell count, degree and weighted-degree summaries, and symmetry diagnostics;
- warnings for multiple components, tiny components, unusually large/small neighborhood choice, and approximate-search reproducibility;
- exact Scanpy, AnnData, NumPy, SciPy, scikit-learn, umap-learn, and OpenBio versions;
- limitations that the graph is exploratory and neither `Technical batch` correction nor `Sample`-level `Condition` inference.

The result prose may say that a graph was constructed; it must not say clusters or biological populations were discovered.

## Equivalent code

`code` should define a standalone function that takes AnnData, resolves the explicit representation, validates dimensions/collisions, calls `sc.pp.neighbors` with explicit parameters, validates the generated keys, and returns an independent AnnData. It must reproduce warnings and overwrite policy and must not import ComfyUI or rely on OpenBio history.

Tests should cross only this interface and cover X/PCA/scVI/Harmony widths, the 10D-scVI regression, all-dimensions behavior, custom keys, collisions, every supported metric/method, invalid/constant/nonfinite representations, oversized neighborhoods, sparse graph invariants, seed repeatability, strict JSON, and generated-code equivalence.

The open-tool implementation preserves oversized expert requests and reports Scanpy's effective `n_neighbors`. Runtime and generated source share the same representation-only preconditions, including support for zero-variable AnnData when an aligned `.obsm` representation is selected.
