# Force-Directed Graph: module design review

## Decision

**Keep and deepen `OpenBioSingleCellForceDirectedGraph` as an atomic named-graph-to-2D-layout module.** It shares graph-resolution helpers with UMAP and Leiden internally, but should not be merged with them: each layout has different scientific meaning, dependencies, storage, and failure modes. It should not be merged with PAGA because PAGA is a coarse partition graph analysis whose output is only one optional initializer.

The seam is:

```text
one named validated observation graph + one explicit layout policy -> one named 2D graph drawing
```

## Proposed interface and P0 correction

- `adata`, required;
- `neighbors_key`, visible, default `neighbors`;
- `layout`, visible bounded combo, default `fr` so the base installation works without silently changing method;
- `init_mode`, advanced bounded combo `random|existing|paga`, default `random`;
- `init_key`, advanced and used only for `existing`;
- `key_suffix`, advanced, default empty (resolved from actual layout);
- `overwrite_existing`, advanced, default `False`;
- `random_seed`, advanced, default `0`.

`init_pos="paga"` must no longer be the default. It remains useful only as an explicit choice after strict validation of PAGA positions, grouping, categories and coarse connectivities.

For `layout="fa"`, fail with an actionable optional-dependency error when fa2-modified is absent. Silent fallback violates the interface because the caller requested one method and receives another. Defaulting to `fr` provides a functioning igraph adapter in the installed environment. Include official `grid_fr` and `rt_circular` only if their behavior is covered; tree layouts may need a separate advanced `root` control or should be omitted until a real workflow requires them.

Keep arbitrary `**kwds`, direct adjacency/obsp sockets, `n_jobs`, and backend objects hidden. If iteration control is needed, expose one validated advanced value translated per supported layout; an untyped kwargs seam is shallow and untestable.

## Deep implementation and locality

The implementation should share a private graph resolver with UMAP/Leiden while retaining its own public interface. It should validate AnnData, graph metadata/matrix, components, layout dependency, initialization, and both coordinate/uns collisions before calling Scanpy.

Wrap Scanpy's legacy `numpy.random` and Python `random` global-state mutation in a module-level lock with snapshot/restore. This internal seam creates locality: one fix protects all calls without making RNG internals part of the node interface. Capture warnings to identify any unexpected backend fallback and fail if executed layout differs from the resolved contract.

After execution, verify the actual key and finite `n_obs × 2` coordinates. Preserve X, raw, layers, graph matrices, annotations, and unrelated embeddings.

## Summary contract

Strict JSON `summary` must include:

- graph-drawing/Scanpy reference plus the selected layout's method reference;
- resolved neighbors/connectivities keys, shape, nnz, degree summaries, connected-component sizes and isolated cells;
- requested and actual layout, backend (`fa2-modified` or python-igraph), dependency version/status, initialization mode/key/PAGA provenance, seed, and fixed iteration policy;
- resolved coordinate key, global `draw_graph` parameter-key overwrite status, coordinate shape/finite status and per-axis range/median/IQR;
- warnings for disconnected/tiny components, tree layout on unsuitable graphs, optional backend constraints, and stochastic/version sensitivity;
- Scanpy, igraph, fa2-modified when used, AnnData, NumPy, SciPy, and OpenBio versions;
- limitations that coordinates are exploratory, not lineage, cluster truth, `Technical batch` integration, or `Sample`-level `Condition` inference.

The results sentence should say that a graph layout was computed. It must name the actual layout and must not say ForceAtlas2 if Scanpy ran Fruchterman-Reingold.

## Equivalent code and tests

`code` should define a standalone function with graph/initialization/collision validation, explicit dependency checks, RNG snapshot/restore, an explicit `sc.tl.draw_graph` call, and postcondition validation. It must reproduce the selected backend and fail rather than fall back. No ComfyUI/OpenBio history dependency is allowed.

Tests should cover the no-PAGA default regression, explicit valid/invalid PAGA, existing-coordinate initialization, missing/malformed/custom graphs, each supported igraph layout, missing/present fa2-modified behavior, fallback detection, disconnected/isolated graphs, coordinate and `uns["draw_graph"]` collisions, global RNG restoration, seed repeatability, strict JSON, input preservation, and generated-code equivalence.

The undirected policy is selected explicitly at the shared graph seam rather than imposed on UMAP. Runtime and generated source both accept zero-variable AnnData when the graph and initialization axes are valid.
