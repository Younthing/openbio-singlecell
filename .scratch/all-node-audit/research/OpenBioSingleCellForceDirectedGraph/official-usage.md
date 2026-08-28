# Force-Directed Graph: official usage and scientific practice

## Audited runtime and primary sources

The installed runtime is Scanpy 1.12.3 with python-igraph 1.0.0. Neither `fa2-modified` nor the legacy `fa2` package is installed. The installed Scanpy signature is:

```python
scanpy.tl.draw_graph(
    adata, layout="fa", *, init_pos=None, root=None, random_state=0,
    n_jobs=None, adjacency=None, key_added_ext=None,
    neighbors_key=None, obsp=None, copy=False, **kwds,
)
```

Primary references:

- [Scanpy `tl.draw_graph` documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.draw_graph.html)
- [Scanpy 1.12.3 implementation and fa2-modified adapter](https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_draw_graph.py)
- [python-igraph layout documentation](https://python.igraph.org/en/main/api/igraph.layout.html)
- [python-igraph visualization/layout guide](https://python.igraph.org/en/latest/visualisation.html)
- [fa2-modified source](https://github.com/AminAlam/fa2_modified)
- [Jacomy et al., ForceAtlas2, DOI 10.1371/journal.pone.0098679](https://doi.org/10.1371/journal.pone.0098679)
- [Fruchterman and Reingold, DOI 10.1002/spe.4380211102](https://doi.org/10.1002/spe.4380211102)
- [Scanpy paper, DOI 10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

## Official operation and graph prerequisite

`draw_graph` lays out an existing sparse observation graph. Scanpy explicitly requires Neighbors first unless a direct `adjacency`/`obsp` matrix is supplied. With `neighbors_key`, it resolves the graph through `uns[neighbors_key]["connectivities_key"]`.

Supported Scanpy layouts include `fa`, `fr`, `grid_fr`, `kk`, `drl`, `lgl`, `rt`, and `rt_circular`. `fa` uses fa2-modified when installed. In Scanpy 1.12, a missing fa2-modified package causes a warning and silently coerces `layout="fa"` to `"fr"`; that changes the executed method and output key. A scientific wrapper must fail clearly or report the resolved layout/backend, not claim ForceAtlas2 when Fruchterman-Reingold ran.

For igraph layouts, Scanpy creates an undirected graph from connectivities. Fruchterman-Reingold, Kamada-Kawai, and DrL accept initialization coordinates; tree layouts have a distinct root control. Algorithm kwargs such as iteration counts vary by backend and should not be passed as an untyped workflow dictionary.

## The `init_pos="paga"` P0

Official `init_pos` options are no initialization/random (`None`/`False`), PAGA (`"paga"`/`True`), or an existing two-dimensional `.obsm` key. Scanpy's own default is `None`, which initializes randomly.

The pre-refactor OpenBio node changed the default to `"paga"`. Scanpy's implementation immediately calls `get_init_pos_from_paga`, which requires at least:

- a valid neighbors graph;
- `adata.uns["paga"]["pos"]`;
- `adata.uns["paga"]["groups"]` naming a categorical obs column;
- coarse PAGA connectivities aligned to those categories.

Neither this node nor the standard upstream graph flow creates PAGA. The refactored node therefore defaults to PAGA-independent random initialization; PAGA remains an explicitly validated advanced choice.

## Output and randomness

Scanpy writes coordinates to `obsm[f"X_draw_graph_{key_added_ext or resolved_layout}"]` and always writes parameters to `uns["draw_graph"]`. Thus `key_added_ext` changes only the coordinate suffix; it is not a full independent namespace for the parameter record.

The implementation seeds NumPy's legacy global RNG for random/PAGA initialization and Python's module-global `random` generator for igraph layouts. Without snapshot/restore and a lock, execution can contaminate unrelated stochastic nodes and race with concurrent calls. A wrapper should localize this upstream side effect.

Force-directed layouts are iterative/stochastic. Fixed seed, resolved backend/layout, graph, initialization and version are all required for reproducibility.

## Preconditions and failure modes

- Resolve and validate the named connectivities graph: square `n_obs × n_obs`, finite, nonnegative, meaningful edges, and compatible uns metadata.
- Require unique observation identifiers.
- Validate layout against the supported set and its optional dependency.
- For `fa`, either require fa2-modified or explicitly choose a different layout; never silently misreport fallback.
- For PAGA initialization, validate all PAGA fields and category alignment before Scanpy.
- For an existing coordinate initializer, require `n_obs × 2`, finite numeric coordinates.
- Tree layouts require a meaningful root policy; disconnected/non-tree graphs can make the result hard to interpret.
- Disconnected graph components, isolated cells, and tiny components materially affect layout and must be disclosed.
- Reject output-coordinate collisions and the global `uns["draw_graph"]` collision unless overwrite is explicit.
- `key_suffix` must be nonblank after normalization and must not create reserved/conflicting keys.

## Scientific interpretation and reporting

A force-directed coordinate map is an exploratory graph drawing. Attractive/repulsive forces and initialization determine visual geometry; Euclidean distances, area, density, and separation in the map are not inferential quantities. It is not trajectory inference even when a tree-like layout is used, and PAGA initialization does not itself establish lineage.

Report the resolved graph key, layout and actual backend, optional-dependency status, initialization and its provenance, seed, output key, graph components/degree diagnostics, coordinate ranges, and software versions. State explicitly that no `Technical batch` correction, clustering, `Sample` aggregation, trajectory inference, or `Condition` test was performed.

## Open expert-tool boundary

The reviewed drawing contract remains a finite, nonnegative, symmetric, zero-diagonal undirected observation graph because the selected layouts are deliberately adapted as undirected graph drawings. Missing optional backends, invalid initialization state, axis mismatch, and non-finite output remain hard. Disconnected components, isolates, small components, and graph-only AnnData with no gene variables are valid and disclosed rather than rejected.
