# OpenBioSingleCellSchistNestedModel — official usage research

Research date: 2026-08-28

## Scope and scientific identity

Schist fits stochastic block models to a cell-neighborhood graph. In nested mode it returns a hierarchy of graph partitions, with cells at the finest level recursively assigned to coarser blocks. It is a cell-group discovery method.

It is **not** a differential-abundance or differential-composition method. It has no Sample, Condition, subject, Technical-batch, design-matrix, effect-size, or hypothesis-test input. The word “inference” in the stochastic-block-model literature means inference of a generative graph model; it does not turn pooled cells into independent biological replicates and does not estimate a Condition contrast.

The single-cell paper describes Schist as a Scanpy-compatible wrapper around graph-tool that identifies cell groups from a k-nearest-neighbor graph, returns every hierarchy level, and can quantify assignment uncertainty. The paper used Schist 0.7.6 and the then-public `schist.inference.nested_model()` entry point. That function name is historical, not the supported current call.

Primary sources:

- Morelli L, Giansanti V, Cittaro D. Nested Stochastic Block Models applied to the analysis of single cell data. *BMC Bioinformatics*. 2021;22:576. [doi:10.1186/s12859-021-04489-7](https://doi.org/10.1186/s12859-021-04489-7)
- Peixoto TP. Hierarchical Block Structures and High-Resolution Model Selection in Large Networks. *Physical Review X*. 2014;4:011047. [doi:10.1103/PhysRevX.4.011047](https://doi.org/10.1103/PhysRevX.4.011047)
- Peixoto TP. Nonparametric Bayesian inference of the microcanonical stochastic block model. *Physical Review E*. 2017;95:012317. [doi:10.1103/PhysRevE.95.012317](https://doi.org/10.1103/PhysRevE.95.012317)

## Supported official release and installation facts

The latest signed upstream release at the research date is [Schist v0.10.0](https://github.com/dawe/schist/releases/tag/v0.10.0), released 2026-08-04 with graph-tool 3 support. The exact released source, rather than an unversioned documentation build or `master`, is the reproducible contract for this review.

Official Schist installation documentation recommends:

~~~text
conda install -c conda-forge schist
~~~

and explains that this installs `graph-tool-base`; full graph-tool is needed only for graphical capabilities. It also documents installing Schist from its official GitHub source after graph-tool is available. Sources: [Schist README](https://schist.readthedocs.io/en/latest/readme.html) and [official repository](https://github.com/dawe/schist).

There are two packaging hazards that must be treated as correctness issues:

- The conda-forge Schist page still listed 0.9.4 when checked, while the upstream repository had released 0.10.0. Version 0.9.4 has materially different analysis parameters and consensus behavior.
- The package currently named `schist` on PyPI is a different knowledge-graph project, not the dawe/schist single-cell package. `pip install schist` must not be suggested or added as a Python extra for this node. Registry sources: [conda-forge Schist package](https://anaconda.org/conda-forge/schist) and [PyPI `schist`](https://pypi.org/project/schist/).

Graph-tool is a compiled C++/Python dependency. Its official installation page recommends conda-forge on GNU/Linux and macOS, and WSL or the official Docker image on Windows; native Windows installation is not offered as the ordinary path. Source: [graph-tool installation instructions](https://graph-tool.skewed.de/installation.html).

A production Adapter therefore needs an explicitly tested environment and must verify package identity, Schist version, graph-tool import, and public-call signature. A bare successful `import schist` is not sufficient evidence that the correct package is installed.

## Current v0.10.0 public call

The exact v0.10.0 entry point is `schist.inference.fit_model`; `schist.inference.__init__` exports only `fit_model` and `fit_model_multi`. A direct nested-model call with the reviewed policy is:

~~~python
import schist as scs

scs.inference.fit_model(
    adata,
    nested=True,
    assortative=False,
    collect_marginals=True,
    n_samples=100,
    key_added="nsbm",
    adjacency=validated_csr_connectivities,
    neighbors_key="neighbors",
    constraint_key=None,
    deg_corr=True,
    directed=False,
    use_weights=False,
    bisection=True,
    simple_init=False,
    n_jobs=1,
    n_iter=10,
    beta=1.0,
    save_model=None,
    copy=False,
    random_seed=123,
)
~~~

Exact sources:

- [v0.10.0 `fit_model` implementation](https://github.com/dawe/schist/blob/v0.10.0/schist/inference/_model.py)
- [v0.10.0 inference exports](https://github.com/dawe/schist/blob/v0.10.0/schist/inference/__init__.py)
- [current Schist inference reference](https://schist.readthedocs.io/en/latest/autoapi/schist/inference/index.html)

The released v0.10.0 signature contains `nested`, `assortative`, `collect_marginals`, `n_samples`, `key_added`, `adjacency`, `neighbors_key`, `constraint_key`, `deg_corr`, `directed`, `use_weights`, `bisection`, `simple_init`, `n_jobs`, `n_iter`, `beta`, `save_model`, `copy`, and `random_seed`. It does not contain the old `n_init`, `tolerance`, `max_iter`, `refine_model`, or `dispatch_backend` controls.

## Version-specific semantic break

Schist 0.9 unified the older inference functions under `fit_model`, but 0.9.4 and 0.10.0 do not merely rename the same arguments:

| Behavior | v0.9.4 | v0.10.0 |
| --- | --- | --- |
| Nested-model selector | `model="nsbm"` | `nested=True, assortative=False` |
| Main robustness mechanism | `n_init` independently minimized states, then consensus | one description-length minimization; when marginals are enabled, posterior partition sampling followed by consensus |
| Repetition control | `n_init=100` | `n_samples=100` and `n_iter=10` |
| Optimization controls | `tolerance`, `n_sweep`, `beta`, `max_iter` | `bisection`, `simple_init`; backend minimizer owns convergence |
| Parallel backend | joblib `dispatch_backend`, `n_jobs` | graph-tool/OpenMP `n_jobs` |
| Refinement | `refine_model`, `refine_iter` | absent from public v0.10.0 call |

Source: [v0.9.4 `fit_model`](https://github.com/dawe/schist/blob/v0.9.4/schist/inference/_model.py), [v0.10.0 `fit_model`](https://github.com/dawe/schist/blob/v0.10.0/schist/inference/_model.py), and [Schist 0.9 release note](https://schist.readthedocs.io/en/latest/release_notes.html).

The older PBMC tutorial still shows `n_init=100` and explains multiple independent minimizations. It is valid evidence for the older interface and the paper's practice, but it is not the v0.10.0 signature. Source: [official PBMC tutorial](https://schist.readthedocs.io/en/latest/clustering_pbmc.html).

An Adapter that conditionally passes whichever parameter happens to exist would make an unchanged workflow mean different analyses. The initial production compatibility matrix should therefore contain only the exact audited v0.10.0 release. Future releases need a new audit branch and explicit regression evidence.

## Input prerequisite: one named cell-neighborhood graph

Schist requires a neighborhood graph built before fitting, ordinarily with `scanpy.pp.neighbors()` or BBKNN. UMAP coordinates are not required. Expression values are not directly read by the model once the graph is supplied, although the upstream representation, dimensions, metric, and neighbor count used to construct that graph materially determine the result.

Official examples:

~~~python
import scanpy as sc
import schist as scs

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, key_added="neighbors")
scs.inference.fit_model(adata, neighbors_key="neighbors")
~~~

The official spatial tutorial demonstrates a non-default named graph:

~~~python
sc.pp.neighbors(adata, key_added="spectral_neighbors", use_rep="X_spectral")
scs.inference.fit_model(adata, neighbors_key="spectral_neighbors")
~~~

Sources: [Schist PBMC tutorial](https://schist.readthedocs.io/en/latest/clustering_pbmc.html), [Schist spatial ATAC tutorial](https://schist.readthedocs.io/en/latest/Spatial_ATAC/spatial_atac.html), and [Scanpy neighbors reference](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html).

For the repository's Scanpy 1.12.3 baseline, `adata.uns[neighbors_key]["connectivities_key"]` must be a nonempty string pointing to an observation-aligned matrix in `adata.obsp`. Schist v0.10.0 also contains a legacy Scanpy <=1.4.6 branch in which connectivities live inside `uns`; the refactored node should not enlarge its Interface to support that obsolete storage layout.

The v0.10.0 graph conversion source has additional consequences:

- it handles `scipy.sparse.csr_matrix` specially, so the Adapter should normalize the resolved matrix to canonical CSR;
- for `directed=False`, it builds the graph from the strict upper triangle;
- for `use_weights=False`, connectivity values define edge presence but are not used as model weights;
- it constructs graph-tool vertices from observed edge endpoints, so isolated trailing observations can produce fewer graph vertices than AnnData observations.

Source: [v0.10.0 graph conversion implementation](https://github.com/dawe/schist/blob/v0.10.0/schist/_utils/_gt_utils.py).

Safe prerequisites are therefore: in-memory AnnData, at least two observations, unique observation identifiers, a sparse square numeric graph aligned to `n_obs`, finite non-negative values, zero diagonal, symmetry within a declared tolerance followed by canonical symmetrization, at least one edge, and no zero-degree observations. The Adapter should pass the validated CSR explicitly through `adjacency` while still recording the resolved `neighbors_key` and connectivity key.

## Actual model, seed, and backend semantics

With `nested=True` and `assortative=False`, v0.10.0 calls graph-tool's `minimize_nested_blockmodel_dl()` and constructs a `NestedBlockState`. The general SBM is not restricted to assortative communities; it can represent other block-to-block connectivity patterns. This differs from the planted-partition model and from modularity optimization.

The nested SBM recursively models block edge-count graphs. Level 0 is the finest cell partition; increasing level numbers are coarser. The implementation prunes redundant one-group levels and appends a single one-group root, so the number of output levels is data-dependent, not a fixed “up to ten” contract. The hierarchy is the model output; a level is not automatically a validated biological cell type.

Graph-tool documents description-length minimization as a stochastic, generally NP-hard optimization that can return different answers because of similar-probability partitions or difficult optima. Its documentation recommends repeated fits and comparing description length in demanding use cases. It gives `minimize_nested_blockmodel_dl()` a complexity of `O(E log^2 V)` and notes optional OpenMP parallel execution. Sources: [graph-tool inference guide](https://graph-tool.skewed.de/static/docs/stable/demos/inference/inference.html) and [`minimize_nested_blockmodel_dl` reference](https://graph-tool.skewed.de/static/docs/stable/autosummary/graph_tool.inference.minimize_nested_blockmodel_dl.html).

Schist v0.10.0 has several seed details that must be disclosed:

- it tests `if random_seed`, so `random_seed=0` means “seed not set” rather than seed zero;
- a nonzero seed is passed to `graph_tool.seed_rng`;
- the source forces description-length minimization to one OpenMP thread when a nonzero seed is present because thread allocation is nondeterministic;
- posterior sampling runs inside an OpenMP context using `n_jobs`.

For a reproducible baseline, require `random_seed >= 1` and fix `n_jobs=1` for both minimization and posterior sampling. Even then, promise repeatability only for the same Schist, graph-tool, platform/build, graph bytes, and parameters; do not promise cross-version or cross-platform bitwise identity.

In v0.10.0, `n_iter` is the number of multiflip MCMC iterations per posterior sample. It is not a maximum optimizer iteration or an observed convergence statistic. `n_samples` is the number of posterior partitions collected when marginals are enabled. It is not the old number of independent initializations. Schist warns when marginal collection uses fewer than 100 posterior samples, so a report-ready node should enforce at least 100 rather than silently proceeding with weaker evidence.

The public call exposes no R-hat, effective sample size, trace, accepted-move summary, or binary convergence flag. Completion of the call must be reported as completed minimization and fixed posterior sampling, not as proven MCMC convergence.

## Parameter exposure for an nSBM-only node

Scientifically meaningful public controls are:

- named graph key;
- output namespace (`key_added`);
- posterior partition count (`n_samples`);
- degree correction;
- a nonzero graph-tool seed;
- overwrite and memory-budget policy.

For this atomic nested-model node, the following settings should be fixed and disclosed rather than turned into a generic model-builder Interface:

- `nested=True`;
- `assortative=False`;
- `collect_marginals=True`;
- `constraint_key=None`;
- `directed=False`;
- `use_weights=False`;
- `bisection=True`;
- `simple_init=False`;
- `n_jobs=1`;
- `n_iter=10`;
- `beta=1.0`;
- `save_model=None`;
- `copy=False` on the one caller-owned output copy.

Directed, weighted, constrained, flat-SBM, and planted-partition modes change the model rather than merely tuning the same operation. Weighted behavior also changed between Schist releases. They require separately motivated Interfaces and audits, not boolean accumulation here.

## Exact AnnData reads and writes in v0.10.0

### Reads

- `adata.uns[neighbors_key]` for graph metadata when an explicit adjacency is not passed;
- the pointer-selected `adata.obsp[connectivities_key]` under modern Scanpy;
- `adata.obs[constraint_key]` only if a constraint is requested; the proposed node fixes this to `None`.

### Writes for nested mode

- `adata.obs[f"{key_added}_level_0"]`, `..._level_1`, and consecutive coarser categorical columns through the one-group root;
- when marginals are collected, `adata.obsm[f"CM_{key_added}_level_{level}"]` for every output level;
- `adata.uns["schist"][key_added]["stats"]`, containing total entropy, modularity by level, and level entropy;
- `adata.uns["schist"][key_added]["blocks"]`, a dictionary of backend block arrays;
- `adata.uns["schist"][key_added]["params"]`, containing a subset of parameters.

The exact v0.10.0 implementation does **not** add a new matrix to `obsp`. The input graph remains there because the node works on an AnnData copy. It also does not store a live graph-tool state in `uns`, despite the API prose mentioning a state; `save_model` instead writes a pickle to an external path. The production node should fix `save_model=None` to avoid an unrelated filesystem side effect and an opaque pickle artifact.

Marginal matrices hold cell-to-block probabilities. Schist's source explicitly notes that some rows may not sum to one because sampled partitions can contain different numbers of groups. Reports must not silently renormalize those arrays or call `1 - max_probability` a calibrated error rate.

Schist silently removes every existing `obs` column beginning with the chosen level prefix and overwrites `uns["schist"][key_added]`; it can also replace marginal matrices with the same names while leaving stale higher-level matrices from an older run. A production Adapter must preflight the entire output family and fail by default, or remove the complete family on the output copy only when overwrite is explicit.

## Resource, memory, and large-data limits

The single-cell paper reports that nSBM was roughly 6–30 times slower than Leiden in its benchmarks, used 100 historical model initializations, and benefited from concurrent processes and graph-tool CPU parallelism. It also states that runtime is driven by graph edges and can be substantial. These figures are context, not a runtime guarantee for v0.10.0.

The official large-sample tutorial says execution time is difficult to predict from dataset size alone. It experiments with subsampling and label projection, but concludes that a subsampled nSBM can reconstruct an inconsistent hierarchy and discourages that approximation without substantial refinement. The node must not silently sketch, subsample, or transfer labels to escape a resource limit. Source: [official large-dataset tutorial](https://schist.readthedocs.io/en/stable/Large_Samples/Large_Samples.html).

Memory is not limited to the input CSR graph:

- caller immutability requires an AnnData copy, including expression matrices and aligned annotations;
- Schist creates upper-triangle coordinate arrays and a graph-tool graph;
- posterior sampling keeps multiple hierarchy partitions;
- `PartitionModeState` and the nested state hold model data;
- marginal collection creates dense cell-by-block matrices at every hierarchy level;
- output categoricals and `uns` block arrays add aligned state.

The finest block count is unknown before fitting, so exact peak memory cannot be predicted. A conservative preflight can use the worst-case `n_obs × n_obs` dense marginal upper bound, a posterior-hierarchy bound proportional to `n_samples × n_obs × hierarchy_depth`, the measurable owned bytes of the AnnData copy, CSR/COO conversion buffers, and a disclosed graph-tool safety factor. It must run before copying AnnData or constructing coordinate arrays. Such an estimate is a guard, not a hard process-RSS cap, because graph-tool's C++ allocations are opaque to Python.

Backed AnnData should fail rather than be implicitly materialized. A failed budget check should recommend an intentionally smaller upstream analysis branch or an appropriately provisioned supported environment; it should not change model mode, drop marginals, enable a fast heuristic, or subsample without a new analyst decision.

## Current implementation deviations

The repository's current node:

- is categorized under `differential-abundance`, although it performs graph clustering;
- calls `schist.inference.nested_model`, which is not exported by Schist 0.9.4 or 0.10.0;
- passes only `random_seed`, leaving graph identity and every model/backend choice implicit;
- permits seed zero even though v0.10.0 interprets zero as unset;
- does not verify that the imported `schist` is the single-cell package;
- has no Schist/graph-tool version contract;
- does not validate or canonicalize the graph and cannot choose a non-default graph explicitly;
- relies on the vendor's destructive output-prefix behavior without collision policy;
- copies full AnnData without a resource preflight;
- does not expose or validate hierarchy levels, nesting, marginals, entropy, or modularity;
- provides no summary, equivalent code, references, versions, limitations, or convergence disclosure;
- has no focused tests or packaged example.

## Reporting references

At minimum the node's `summary.references` should include:

- Morelli L, Giansanti V, Cittaro D. 2021. [doi:10.1186/s12859-021-04489-7](https://doi.org/10.1186/s12859-021-04489-7)
- Peixoto TP. 2014. [doi:10.1103/PhysRevX.4.011047](https://doi.org/10.1103/PhysRevX.4.011047)
- Peixoto TP. 2017. [doi:10.1103/PhysRevE.95.012317](https://doi.org/10.1103/PhysRevE.95.012317)
- Peixoto TP. The graph-tool python library. 2014. [doi:10.6084/m9.figshare.1164194](https://doi.org/10.6084/m9.figshare.1164194)
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. *JOSS*. 2021;6:4371. [doi:10.21105/joss.04371](https://doi.org/10.21105/joss.04371)

Runtime software versions must be collected dynamically for openbio-singlecell, Schist, graph-tool, Scanpy, AnnData, NumPy, pandas, SciPy, Python, and the active platform/build information needed to interpret graph-tool reproducibility.
