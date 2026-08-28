# OpenBioSingleCellLeiden — official usage research

## Scope and current environment

This record covers one Leiden partition of one declared observation-neighbor graph. It does not cover graph construction, representation learning, marker ranking, cell-type annotation, or Condition inference.

Repository environment inspected on 2026-08-28:

- `scanpy 1.12.3` (repository requirement `scanpy[leiden,scrublet]>=1.12.3,<1.13`);
- `python-igraph 1.0.0`;
- `leidenalg 0.12.0`.

The refactored node calls `scanpy.tl.leiden` with visible `resolution=1.0`, fixed `flavor="igraph"` and `directed=False`, explicit iteration/repeat controls, and the selected `neighbors_key`; it stores one categorical `obs` column and returns strict `summary` plus equivalent `code`.

## Primary sources

- Scanpy Leiden reference: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html
- Scanpy 1.12.3 implementation: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_leiden.py
- Scanpy neighbor-graph reference: https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html
- python-igraph `Graph.community_leiden`: https://igraph.org/python/versions/latest/api/igraph.Graph.html#community_leiden
- leidenalg `find_partition` and partition reference: https://leidenalg.readthedocs.io/en/stable/reference.html#leidenalg.find_partition
- leidenalg method-author repository: https://github.com/vtraag/leidenalg
- Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*. 2019;9:5233. https://doi.org/10.1038/s41598-019-41695-z; open full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC6435756/
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0

The Leiden paper establishes the algorithmic motivation and connected-community guarantees relative to Louvain. It does not establish that a graph partition is a biological cell type, that one resolution is universally correct, or that clustering cells supplies replicate-aware evidence about a Condition.

## Scanpy 1.12.3 contract

The inspected signature is:

```python
scanpy.tl.leiden(
    adata,
    resolution=1,
    random_state=0,
    key_added="leiden",
    adjacency=None,
    directed=None,
    use_weights=True,
    n_iterations=-1,
    partition_type=None,
    neighbors_key=None,
    obsp=None,
    copy=False,
    flavor=None,
    **clustering_args,
)
```

Scanpy requires a graph produced by `scanpy.pp.neighbors`, BBKNN, or an equivalent aligned adjacency. With `neighbors_key`, Scanpy resolves the connectivities matrix through the neighbor metadata and corresponding `obsp` key. `scanpy.pp.neighbors(key_added=None)` conventionally stores metadata in `uns["neighbors"]` and matrices in `obsp["distances"]` and `obsp["connectivities"]`; a custom key stores custom-key-prefixed matrices.

The function writes categorical string memberships to `obs[key_added]` and parameters plus a backend-reported modularity value to `uns[key_added]`. Cluster identifiers are nominal labels (`"0"`, `"1"`, …), not ordered quantities.

## Flavor, directedness, objective, weights, and iterations

These controls are coupled and must not be treated as interchangeable toggles:

- In Scanpy 1.12.3, `flavor="igraph"` requires an undirected graph. Scanpy rejects `directed=True`, constructs an undirected igraph graph, uses edge weights by default, and sets `objective_function="modularity"` unless overridden.
- `igraph.Graph.community_leiden` itself defaults to CPM, `resolution=1.0`, `beta=0.01`, and two iterations; Scanpy's igraph adapter deliberately changes the objective default to modularity.
- In the `leidenalg` flavor, Scanpy defaults to `RBConfigurationVertexPartition`, uses a linear `resolution_parameter`, defaults `directed` to true when unspecified, and passes the random seed to `leidenalg.find_partition`.
- The two flavors therefore differ in adapter behavior, directedness, default partition/objective semantics, and performance. Their resolution values and modularity/quality results should not be pooled as if they were the same method.
- `n_iterations=2` is the underlying-package fast default and the future igraph-oriented behavior described by Scanpy 1.12.3. python-igraph documents any negative value as iteration until a stable iteration, zero as returning the initial singleton partition without optimization, and positive values as a fixed iteration count. The exact integer must be disclosed because it can change the partition.

For a cohesive, reproducible baseline, keep `flavor="igraph"`, `directed=False`, `use_weights=True`, and modularity objective fixed and report them. Expose `n_iterations` as an advanced scientific integer with default 2. Supporting `leidenalg`, directed graphs, CPM, custom partitions, or unweighted graphs would require a separately documented mode because those choices change the scientific contract rather than merely performance.

## Resolution and random initialization

`resolution` must be finite and non-negative. Zero is a valid zero-penalty modularity setting that commonly collapses the graph to few communities and is therefore warned, not rejected; negative values do not satisfy this node's finite-modularity result contract. Higher values generally produce more, smaller communities. Resolution has no universal biological calibration and is not directly comparable across different neighbor graphs, representations, weights, objectives, or backend flavors. `resolution=1.0` is the official Scanpy/igraph default under the resolved modularity path, not a claim that it is optimal for single-cell annotation.

`random_state` controls stochastic optimization. A fixed integer makes the chosen run reproducible under a fixed graph and software stack, but one seed does not demonstrate robustness. The stored partition should use the declared base seed. Optimization stability can be assessed without changing the primary labels by repeating the same graph/resolution/settings under deterministically derived seeds and reporting pairwise adjusted Rand index (ARI). ARI is label-invariant and therefore appropriate when arbitrary cluster IDs differ. It measures optimization reproducibility, not biological validity.

- scikit-learn adjusted Rand reference: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html
- Hubert L, Arabie P. Comparing partitions. *Journal of Classification*. 1985;2:193–218. https://doi.org/10.1007/BF01908075

A practical default is five total starts (base start plus four deterministic alternatives), with mean, minimum, and maximum pairwise ARI disclosed. Any positive repeat count is accepted without an arbitrary cap; one start means stability was not assessed and yields null ARI summaries, while large workloads are warned.

## Graph preconditions and diagnostics

Before clustering, resolve the exact connectivity key and validate the matrix rather than relying on a later library error:

- in-memory, nonempty AnnData with at least two observations and unique `obs_names`;
- `uns[neighbors_key]` exists and identifies an existing connectivity matrix in `obsp`;
- exact shape `n_obs × n_obs`, numeric finite values, no negative weights, and at least one positive off-diagonal edge;
- zero diagonal and symmetry within a documented tolerance for the fixed undirected igraph mode;
- canonical sparse form with duplicates summed before graph statistics;
- disclose positive edge count, graph density, connected-component count, and isolated-cell count.

An empty graph, wrong shape, negative/non-finite weight, asymmetric graph, or missing pointer is an error. Disconnected graphs are not necessarily invalid, but Leiden cannot merge observations across components; report the component count. Isolated cells can become singleton communities; retain them but issue a prominent warning and disclose their identifiers only in a bounded form.

The report should inherit relevant upstream neighbor parameters from `uns[neighbors_key]["params"]`—for example representation, PC count, neighbor count, metric, transformer/kernel, and random state—without claiming values that are absent. Clusters depend materially on these upstream choices and on whether the representation modeled a Technical batch.

## Output and reporting requirements

The primary output is a copy of AnnData with:

- categorical labels in `obs[key_added]`;
- Scanpy parameters/modularity plus OpenBio diagnostic metadata under `uns[key_added]`;
- input expression matrices, Raw snapshot, representations, graph, and unrelated annotations unchanged.

Fail on collisions in both `obs[key_added]` and `uns[key_added]` unless explicit advanced overwrite is enabled. One collision is never overwritable: after resolving the named graph, reject `key_added == neighbors_key` for both overwrite settings because `uns[neighbors_key]` is the graph metadata/pointer bundle required by downstream graph consumers. Never silently overwrite a curated annotation or graph provenance.

The strict JSON `summary` should include:

- methods: exact graph key, Scanpy/igraph Leiden, modularity objective, weighted undirected graph, resolution, iteration policy, base seed, and stability-repeat definition;
- results: number of clusters, complete cluster-size/fraction mapping, minimum/median/maximum size, singleton count, largest-cluster fraction, modularity, graph components/isolates, and pairwise-ARI stability summary;
- parameters: all visible and hidden resolved controls plus actual connectivity key and inherited neighbor settings;
- warnings/limitations: resolution and graph dependence, isolated/disconnected graph effects, seed stability versus biological validity, exploratory status, and software-version dependence;
- references: Leiden paper, Scanpy, igraph, ARI when used, and AnnData aligned-storage documentation;
- software versions: Python, openbio-singlecell, Scanpy, AnnData, python-igraph, NumPy, pandas, SciPy, and scikit-learn. Do not list leidenalg as the active backend when fixed igraph flavor is used.

Report-ready prose may say that Leiden partitioned the declared neighbor graph into K communities and report their sizes/stability. It must not call the labels true cell types, curated annotations, or Condition differences. Under this repository's domain language, later marker rankings are Cluster marker evidence and remain exploratory; Condition inference requires Sample-level replication.

The `code` endpoint should contain an ordinary Python function that resolves and validates the same graph, copies AnnData, calls `scanpy.tl.leiden` with every resolved parameter explicit, computes the same graph/cluster/stability diagnostics, and returns the annotated AnnData. It must not depend on ComfyUI or plugin history bookkeeping.

## Open expert-tool boundary

Finite `resolution=0` is an executable zero-penalty expert choice and is accepted with a one-community/low-resolution warning; negative resolution remains unsupported by the finite-modularity report contract. `n_iterations` accepts every integer: negative values request iteration to stability, zero preserves the initial partition without optimization, and positive values are fixed iteration counts. Repeat count is positive but has no arbitrary scientific maximum. Complete categorical memberships may use any backend labels; unused categories are removed and labels are normalized deterministically without changing the partition. The gene axis is irrelevant to graph partitioning.
