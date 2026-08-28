# OpenBioSingleCellLeidenResolutionSweep — official usage research

## Scope and current state

This node compares Leiden partitions over several non-negative resolution values on one unchanged neighbor graph. The refactored implementation parses an ordered comma-separated grid, reruns fixed `flavor="igraph"` and `directed=False` with declared iteration/repeat controls, writes one `obs` column per resolution, and returns AnnData plus structured metrics, strict `summary`, and equivalent `code`.

Current inspected environment is Scanpy 1.12.3, python-igraph 1.0.0, and leidenalg 0.12.0. The official behavior and primary sources are the same as recorded for `OpenBioSingleCellLeiden`:

- Scanpy Leiden: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html
- Scanpy 1.12.3 implementation: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_leiden.py
- Scanpy neighbor graphs: https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html
- python-igraph Leiden: https://igraph.org/python/versions/latest/api/igraph.Graph.html#community_leiden
- leidenalg reference: https://leidenalg.readthedocs.io/en/stable/reference.html#leidenalg.find_partition
- Traag VA, Waltman L, van Eck NJ. *Scientific Reports*. 2019;9:5233. https://doi.org/10.1038/s41598-019-41695-z

## What a resolution sweep can and cannot establish

Under Scanpy's fixed igraph/modularity path, higher resolution generally yields more, smaller communities. The mapping is graph- and objective-dependent; it is not a calibrated cell-type scale. A sweep is decision-support evidence for inspecting granularity, size pathologies, and robustness. It must not automatically select a “best” resolution from maximum modularity or maximum cluster count.

Partitions at different resolutions have nominal labels that cannot be matched by label string. Label-invariant adjusted Rand index (ARI) can quantify agreement between adjacent resolutions. Call this `adjacent_resolution_ari`, not stochastic stability: it measures how much the partition changes when resolution changes.

Stochastic stability is a separate within-resolution calculation. For each resolution, rerun the same graph/objective/iteration policy under the analyst-declared deterministic seed sequence and report mean/min/max pairwise ARI. This measures optimization reproducibility only. A stable but biologically unhelpful partition is possible, and an unstable result warrants inspection rather than automatic deletion.

- scikit-learn ARI documentation: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html
- Hubert L, Arabie P. Comparing partitions. *Journal of Classification*. 1985;2:193–218. https://doi.org/10.1007/BF01908075

## Shared graph and algorithm contract

Every resolution must use the exact same resolved connectivity matrix, weights, backend, objective, directedness, iteration count, and base/repeat seed policy. Validate once before the loop:

- in-memory AnnData, at least two observations, unique observation names;
- declared neighbor metadata and connectivity `obsp` key exist;
- matrix shape exactly `n_obs × n_obs`, numeric finite non-negative weights, zero diagonal, symmetric for the fixed undirected flavor, and at least one positive edge;
- report edge count/density, connected components, and isolated observations.

Fix and disclose `flavor="igraph"`, `objective_function="modularity"`, `directed=False`, and `use_weights=True`. Expose `n_iterations` as an advanced integer (`2` by default; negative until stable, zero no optimization, positive fixed iterations). Do not mix igraph and leidenalg results in one table: Scanpy's adapters use different default objective/partition semantics and directedness.

## Resolution-list contract

The current comma-separated interface can be retained for workflow migration, but parsing must be strict:

- input is a string with one or more nonempty tokens; large workloads are disclosed rather than rejected by an arbitrary token cap;
- every token parses to a finite non-negative float;
- duplicate numeric resolutions are rejected rather than silently removed;
- values are sorted ascending for the report/table while preserving exact resolved values;
- generated result keys are unique and collision-checked in both `obs` and `uns`;
- no generated result key may equal the resolved `neighbors_key`, regardless of overwrite policy, because that `uns` slot stores the graph metadata/pointer bundle;
- a collision-proof canonical decimal token is preferred; if compact formatting is used, reject any key collision before the first Leiden run.

The default `0.25,0.5,1.0,2.0` is a broad exploratory grid, not an official biological recommendation. The exact list belongs in methods/results.

## Required structured output

Multiple `obs` columns alone are not a usable resolution analysis. Keep aligned categorical memberships in AnnData because they are the natural cell-level storage, but also return a structured table with one row per resolution and stable columns:

- `resolution` and generated `key`;
- `n_clusters`;
- `min_cluster_size`, `median_cluster_size`, `max_cluster_size`;
- `smallest_cluster_fraction`, `largest_cluster_fraction`;
- `singleton_clusters`;
- backend-reported `modularity`;
- `stability_repeats`, `stability_mean_ari`, `stability_min_ari`, `stability_max_ari`;
- `adjacent_previous_resolution` and `adjacent_resolution_ari` (`None` for the first row in JSON form).

The table must be numeric/JSON-safe and sorted by resolution. It must contain no `best_resolution`, automatic recommendation, cell-type name, or Condition claim. The complete cell memberships remain categorical `obs` columns; repeat-run memberships stay internal and must not create dozens of extra columns.

## Cluster-size and stability disclosure

For each resolution, cluster sizes must sum exactly to `n_obs`. Report all cluster sizes in the JSON summary using their resolution key, plus bounded table statistics. Flag singleton clusters, very small clusters, one-cluster results, or partitions dominated by one cluster as review warnings, not proof of error. Thresholds used for warnings must be disclosed and described as heuristics.

Report graph connected components and isolates once at sweep level. A partition cannot merge disconnected graph components, so a high minimum cluster count may be imposed by the graph rather than resolution. Stability summaries use identical runs except for seed; adjacent-resolution ARI uses the stored base-seed memberships.

Modularity is diagnostic, not a universal selection score. It depends on graph construction, weights, objective, and resolution, and maximizing it across a resolution grid can favor an analystically inappropriate granularity.

## Outputs, summary, and code

Outputs should be:

- annotated `adata` containing one categorical base-seed partition per resolution;
- structured `k_metrics`-style table renamed for clustering, e.g. `resolution_metrics`;
- strict JSON `summary`;
- equivalent Python source `code`.

The `summary` includes exact graph provenance and diagnostics, ordered resolution records, complete cluster-size mappings, within-resolution ARI stability, adjacent-resolution ARI, hidden fixed policy, seed derivation, warnings/limitations, method/software references, and versions for Python, openbio-singlecell, Scanpy, AnnData, python-igraph, NumPy, pandas, SciPy, and scikit-learn.

Report-ready prose should state that several granularities of the same graph were surveyed and summarize cluster-count/size/stability changes. It must not say that a resolution discovered the true number of cell types. Cluster marker evidence remains exploratory, and formal Condition inference remains Sample-level under the repository ADR.

`code` should implement one ordinary Python function that validates the graph and resolution list, copies AnnData, executes every base and repeat run with all arguments explicit, returns annotated AnnData plus a pandas metrics table, and reproduces runtime failures. It must not depend on ComfyUI or plugin history bookkeeping.

## Failure conditions

Reject before the first run:

- all graph failures listed above;
- empty, nonnumeric, NaN, infinite, negative, duplicate, or key-colliding resolutions; zero is valid and large grids are workload warnings;
- an empty prefix;
- invalid `n_iterations`, repeat count, or random seed;
- any target `obs` or `uns` collision unless overwrite is explicit, and unconditionally any target key equal to the resolved graph's `neighbors_key`;
- a resolved workload that cannot be executed by the backend or host; large but executable workloads are reported and warned rather than rejected by a fixed wrapper ceiling.

After each run, reject missing/noncategorical memberships, missing cells, non-finite modularity, invalid cluster-size totals, or a backend result that cannot be aligned exactly to the original observation order.

## Open expert-tool boundary

Resolution zero is accepted and disclosed; negative/non-finite values remain invalid. One start is sufficient to survey resolutions, with all within-resolution stability fields represented as JSON `null` and a warning that stability was not assessed. Fixed ceilings of 20 starts, 256 resolutions, or 256 total runs are not scientific or safety invariants and must not reject expert choices; the resolved workload is reported before execution. Any nonblank key prefix is allowed because output keys are value data, not paths or executable source. Output-key collisions remain hard and are preflighted atomically.
