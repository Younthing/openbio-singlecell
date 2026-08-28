# OpenBioSingleCellSchistNestedModel — module design review

Research date: 2026-08-28

## Decision

**Keep the node ID, move the node from differential abundance to graph clustering, and rebuild it as a version-gated nested-SBM hierarchy Module. Do not merge it with Leiden or Leiden Resolution Sweep. Do not advertise or add it to an example until the exact optional environment has a reproducible installation and a real-backend smoke test.**

The retained atomic question is:

> Given one already constructed, named, undirected cell-neighborhood graph, what hierarchy of cell blocks does the reviewed Schist nested stochastic block model infer, and what uncertainty and description-length diagnostics accompany that hierarchy?

This is cohesive. Graph validation, one nSBM fit, hierarchy validation, marginal diagnostics, and reporting all serve that question. Sample composition, Condition inference, annotation curation, level selection, plotting, graph construction, flat SBM, PPBM, and multimodal Schist answer different questions and remain outside the Module.

Move the public class out of `nodes_abundance.py` into a clustering-owned module and register it under `openbio/single-cell/clustering`. Preserve `OpenBioSingleCellSchistNestedModel` so saved workflow identity is stable. A display such as “Schist Nested-SBM Hierarchy” is more precise than “Schist Nested Model.”

## Why retain rather than delete

The deletion test favors retention. Without this Module, a caller must understand modern Scanpy graph indirection, CSR normalization, Schist's version-specific signature, graph-tool seed/thread behavior, destructive prefix writes, hierarchy ordering, marginal storage, and backend validation. That complexity would be recreated in notebooks or workflow glue.

The method is also not redundant with Leiden:

- Leiden optimizes a chosen graph-partition objective at an analyst-selected resolution and returns a flat partition.
- nSBM fits a nonparametric generative graph model and returns an intrinsic hierarchy of nested partitions.
- the general SBM is not restricted to assortative connectivity, while the repository's Leiden Module fixes a weighted undirected modularity objective;
- Schist has a heavy graph-tool optional dependency, posterior sampling, different stochastic behavior, different reports, and substantially different resource costs.

Do not merge it with Resolution Sweep. A sweep compares independent Leiden resolutions and explicitly does not assert that they form a hierarchy. Schist returns one model-defined hierarchy, and no resolution grid is involved. Combining them would create a broad clustering workbench with low Locality and incompatible parameter semantics.

## Current Interface and correctness failures

The current Interface has `adata` and `random_seed`, returns only `adata`, and calls:

~~~python
schist.inference.nested_model(output, random_seed=random_seed)
~~~

This fails against current official releases because `nested_model` is no longer exported. It also hides the graph key, model family, marginal/consensus policy, degree correction, directedness, weight policy, optimizer policy, backend parallelism, output namespace, overwrite behavior, and optional-package version.

The current category makes the more serious scientific error: it places pooled-cell graph clustering beside Sample-aware differential-abundance Modules. Nothing in the Implementation compares Conditions or uses Samples. A downstream reader could mistake model “inference” for formal biological inference.

Other shallow behavior includes:

- seed zero is accepted even though Schist v0.10.0 treats it as no seed;
- any importable package named `schist` is accepted, despite the current PyPI name collision;
- the default graph is implicit and unvalidated;
- AnnData is copied before any graph/resource failure;
- existing Schist outputs can be silently deleted or overwritten by upstream code;
- no level order, nesting invariant, root, marginals, entropy, or modularity is checked;
- no convergence limitation, graph provenance, references, versions, summary, or equivalent code is returned.

## Deep Module and correct Seam placement

Use the repository's named-graph logic as the shared internal Seam with Leiden and Resolution Sweep. That Seam should resolve:

~~~text
neighbors_key
  -> adata.uns[neighbors_key]["connectivities_key"]
  -> adata.obsp[connectivities_key]
  -> validated canonical symmetric CSR graph + diagnostics + fingerprint
~~~

Schist-specific adaptation belongs behind a second private Seam. It verifies package identity/version/signature, maps the stable node policy to exact v0.10.0 keywords, invokes the backend, discovers the returned hierarchy, validates all backend state, and returns normalized diagnostics. It must not expose graph-tool states, `PartitionModeState`, OpenMP contexts, old joblib backend names, vendor cache/state dictionaries, or pickle paths.

This gives the public Module Depth: the Interface contains the few scientific and workflow choices, while version drift, graph representation, storage families, backend invariants, and report construction remain in the Implementation. The Schist Adapter is a local-substitutable dependency in tests; graph-tool and Schist remain true external dependencies.

Do not merge neighbor construction into this Module. The upstream representation, number of dimensions, metric, neighbor count, and integration choice are scientific decisions with their own report and fingerprint. The Schist summary consumes and discloses that graph provenance rather than recreating it.

## Proposed atomic Interface

### Inputs

Visible:

- `adata`.

Advanced but explicit:

- `neighbors_key="neighbors"` — the exact Scanpy named graph;
- `key_added="nsbm"` — owned output namespace;
- `posterior_samples=100` — mapped to v0.10.0 `n_samples`, minimum 100;
- `degree_correction=True` — mapped to `deg_corr` and always reported;
- `overwrite_existing=False` — all-or-nothing output-family replacement on the copy;
- `max_working_memory_gib=8.0` — conservative preflight, not a process-RSS promise;
- `random_seed=123` — valid range `1..2**31-1`; zero is rejected.

These controls are exposed because they change the graph analyzed, output ownership, posterior evidence, model family, or execution safety. Keeping them advanced avoids making ordinary use noisy without hiding them from a reproducible workflow.

### Fixed and disclosed policy

- exact initial compatibility target: official Schist v0.10.0 release only;
- `nested=True`, `assortative=False`;
- `collect_marginals=True`;
- `constraint_key=None`;
- `directed=False`, `use_weights=False`;
- `bisection=True`, `simple_init=False`;
- `n_jobs=1`;
- `n_iter=10`, `beta=1.0`;
- `save_model=None`;
- pass one validated CSR through `adjacency`;
- make exactly one AnnData output copy, then call Schist with `copy=False`;
- return all levels including the one-group root;
- no automatic preferred-level selection.

The fixed policy deliberately omits a model dropdown, directed/weighted switches, constraint annotations, speed heuristics, state-file path, and backend selector. Each would widen the scientific operation or leak an Implementation detail. If a real workflow later needs weighted nSBM, constrained partitions, PPBM, flat SBM, multimodal fitting, or approximate large-data fitting, it should receive its own research review and Interface.

### Outputs

- `adata` — copied input with validated Schist hierarchy, marginal arrays, and metadata;
- `summary` — strict JSON report with graph/model/results/references/versions;
- `code` — source text defining one equivalent public-package function.

No separate table is needed. Per-cell memberships and marginals are observation-aligned AnnData state, while the bounded per-level aggregate records fit naturally in `summary`. A table output would duplicate the same level diagnostics without serving a registered consumer.

## Version gate and optional-dependency Adapter

The initial Adapter must fail closed unless all of the following hold:

1. the imported module exposes the single-cell Schist authorship/repository identity or another robust package-identity signal;
2. its canonical parsed version is exactly the audited `0.10.0` release;
3. `schist.inference.fit_model` exists;
4. its public signature contains the exact required v0.10.0 keyword set;
5. `graph_tool.all` imports and exposes version/build information.

The error should name the imported module path and detected version and point to the official conda/source and graph-tool installation documentation. It must not advise `pip install schist`, because that name currently resolves to a different project on PyPI.

Do not add a `schist` Python extra to `pyproject.toml`. A reproducible optional environment should be documented/tested with conda-forge graph-tool and the exact official Schist tag. Native Windows should fail with guidance to WSL or a supported container/environment rather than an OSError-shaped generic message.

The v0.9.4 Adapter should not be retained as an invisible fallback. Its `n_init`/joblib consensus means a different analysis. If supporting 0.9.4 is operationally necessary, create a named compatibility branch with its own parameter/report contract and parity tests; do not make one workflow change meaning based on installed version.

## Graph validation and canonicalization

Before importing the expensive backend, copying AnnData, or materializing graph coordinates:

- reject backed AnnData and fewer than two observations;
- require unique observation identifiers;
- require a nonempty string `neighbors_key`;
- require mapping-like `adata.uns[neighbors_key]` with a nonempty `connectivities_key` pointer;
- require the pointer target in `adata.obsp`;
- require a SciPy sparse matrix with exact `(n_obs, n_obs)` shape;
- normalize to CSR, sum duplicates, eliminate explicit zeros, and reject nonnumeric data;
- reject non-finite or negative entries;
- require a zero diagonal within the repository's fixed tolerance;
- require symmetry within a disclosed tolerance, then canonicalize by averaging with the transpose;
- reject a graph with no positive edges or any zero-degree observation;
- compute graph diagnostics and a content fingerprint from observation identity, CSR structure, and numeric values.

Multiple connected components may be analyzed, but their count and sizes are reported and a warning is emitted. Automatically connecting components would alter the graph. Graph weights are fingerprinted and reported as upstream connectivity values but are ignored by the fixed unweighted nSBM policy; the report must say so explicitly.

Passing the canonical CSR through `adjacency` avoids Schist's CSR/other-sparse branching and legacy Scanpy storage logic. Passing the resolved `neighbors_key` too keeps vendor parameters and the OpenBio report traceable to the upstream graph.

## Output namespace and transaction policy

For a proposed `key_added=K`, the owned family is:

- all `obs` columns matching `K_level_<nonnegative integer>`;
- all `obsm` matrices matching `CM_K_level_<nonnegative integer>`;
- `uns["schist"][K]`;
- OpenBio history/report metadata for this operation.

Reject a blank/unsafe namespace. Before copying or fitting, inspect the complete family. If any member exists and `overwrite_existing=False`, fail with the colliding keys. If overwrite is true, delete the complete family on the output copy immediately before fitting; never modify the caller. Do not let Schist's broad `startswith` deletion define ownership, and reject confusing noncanonical prefix collisions that Schist itself would delete.

No partial output is published if validation, fitting, or postconditions fail. Work only on the private copy and return it after all invariants pass.

## Backend call and output postconditions

Pass every reviewed Schist keyword explicitly. Do not rely on upstream defaults, because the default for `collect_marginals` has already differed between the v0.10.0 tag and later documentation/source builds.

After the call, require:

- consecutive categorical `K_level_0..K_level_L` columns and no other owned level keys;
- exact observation index/order, no missing labels, canonical used categories, and positive cluster sizes summing to `n_obs` at each level;
- level 0 finest, nonincreasing cluster count with level, and a one-cluster final root;
- every fine cluster maps to exactly one cluster at the next coarser level;
- `uns["schist"][K]` with mapping-like `stats`, `params`, and `blocks`;
- finite scalar entropy plus finite modularity and level-entropy arrays aligned to the output levels;
- recorded vendor parameters consistent with the requested model, graph, seed, marginal, directedness, weight, and degree-correction policy;
- one `CM_K_level_i` matrix per level because marginals are fixed on;
- marginal matrices aligned to observations, with column count equal to the used categories at that level, finite values in `[0, 1]`, and no empty-probability observation;
- no new or changed `obsp` key except preservation of the copied named graph;
- all untouched input axes, X/layers, Raw snapshot, annotations, and graph metadata preserved;
- successful H5AD round-trip without a live graph-tool object or external state reference.

Do not force marginal rows to sum to one: upstream documents that differing sampled partition sizes can leave incomplete mass in the represented consensus groups. Report row-sum coverage and maximum-membership distributions instead.

## Resource preflight

Normalize and validate `max_working_memory_gib` before any observation-scale Python list, AnnData copy, COO coordinate array, or backend state is created. Compute a conservative estimate containing:

- measurable owned buffers that `adata.copy()` will duplicate, including X, layers, aligned matrices, DataFrames, and sparse graph storage;
- canonical CSR plus strict-upper-triangle index buffers used to construct graph-tool;
- a disclosed graph-tool vertex/edge safety factor;
- posterior hierarchy samples bounded by `posterior_samples`, `n_obs`, integer width, and a conservative maximum depth;
- the worst-case dense finest-level marginal matrix `n_obs × n_obs`, plus coarser matrices;
- output categorical codes, block arrays, and bounded report structures.

Use Python integer arithmetic to avoid overflow. Reject before the copy if the estimate exceeds the budget and include the estimated components in the error. The summary records the same estimate, budget, safety factors, graph size, and that opaque C++ allocations mean this is not a guaranteed RSS ceiling.

Do not implement an automatic fallback to `collect_marginals=False`, `simple_init=True`, disabled bisection, fewer posterior samples, a flat model, PPBM, subsampling, or label transfer. Each changes the analysis. A user can intentionally make a smaller upstream subset or provision a larger supported environment.

## Hierarchy reporting without automatic level selection

The Module returns every hierarchy level. It does not add `selected_level`, `best_level`, a cell-type alias, or an automatic recommendation.

For every level, `summary.key_results.hierarchy_levels` contains an ordered strict-JSON record with:

- level index and exact `obs`/`obsm` keys;
- `is_finest`, `is_root`, and number of clusters;
- minimum, median, maximum, singleton count, and largest-cluster fraction;
- modularity and level entropy from Schist;
- marginal matrix shape;
- marginal row-sum minimum/median/maximum;
- maximum-membership-probability minimum/median/maximum;
- whether the nesting invariant to the next level passed.

Cluster-size maps can be bounded deterministically for JSON size while recording the total and truncation policy; the full memberships remain in AnnData. The results prose may say that Schist inferred a hierarchy with a stated number of levels and block counts. It may not call a particular level optimal, identify blocks as biological cell types, or claim differential abundance.

Level selection is a downstream analyst decision informed by markers, stability, Sample composition, known biology, and potentially method-specific diagnostics. The Schist paper discusses modularity and random-matrix strategies, but automatically implementing one here would merge a second analysis decision into the fitting Module.

## `summary` contract

The strict JSON output contains no NaN or Infinity and includes:

- node ID, title, operation, schema version, and input dimensions;
- exact named-graph key, connectivity pointer, graph fingerprint, edge count, density, degrees, isolates, connected components, and available upstream graph parameters/representation provenance;
- exact public and fixed model parameters, including the fact that connectivity weights were ignored;
- Schist/graph-tool backend identity, seed, one-thread policy, posterior sample count, MCMC sweeps per sample, and beta;
- hierarchy records described above, total entropy, informative-level count, root level, and output storage keys;
- a completion field that says the backend returned and passed postconditions;
- `convergence_diagnostic_available=false` and an explicit explanation that no ESS, R-hat, trace, or formal convergence status is exposed;
- the conservative memory estimate/budget and runtime duration;
- warnings for disconnected graphs, weak marginal mass, dominant/singleton-heavy levels, or output overwriting;
- limitations that memberships are exploratory graph-model evidence, depend on the upstream graph, may reflect Sample/Technical batch, and are not Curated annotation or Sample-aware Condition inference;
- method/package references from `official-usage.md`;
- dynamic versions for openbio-singlecell, Schist, graph-tool, Scanpy, AnnData, NumPy, pandas, SciPy, Python, OS/platform, and graph-tool build/OpenMP information.

References are a JSON list of structured citation records, not one prose blob. Versions must come from the imported runtime, not hard-coded documentation values.

## Equivalent `code` contract

The source output defines one callable such as `run_schist_nsbm(adata)` and uses only public packages. It must:

- verify the exact official Schist version/package identity and graph-tool availability;
- validate settings, graph, namespace, and memory before copying;
- canonicalize and fingerprint the named CSR graph;
- copy once and apply the complete overwrite transaction on that copy;
- call `schist.inference.fit_model` with every fixed/public argument explicit;
- run the same level, nesting, stats, parameter, marginal, graph-preservation, and alignment postconditions;
- return the annotated AnnData.

It may omit Comfy result wrapping, plugin execution history, and timing. It may not call `nested_model`, rely on Schist defaults, choose another version branch, write a pickle, omit the marginal/consensus mode, skip resource validation, mutate the caller, or silently select one hierarchy level.

Because graph-tool output can be build-sensitive, generated-code parity tests should run with an injected deterministic fake Adapter for exact state/failure equivalence and with the exact optional real environment for a small seeded smoke comparison.

## Migration

1. Preserve the node ID but move its class and registry ownership out of `nodes_abundance.py` into clustering; update the category to `openbio/single-cell/clustering`.
2. Change the display to “Schist Nested-SBM Hierarchy” while keeping the old ID for saved graphs.
3. Migrate the legacy widget vector so its former second value remains `random_seed`; insert new advanced inputs at stable defaults without shifting that value into `neighbors_key` or `posterior_samples`.
4. Replace the legacy `nested_model` call with the exact v0.10.0 Adapter. Do not auto-convert an installed 0.9.x execution.
5. Add `summary` and `code` after `adata`; existing links from output zero remain valid.
6. Remove every differential-abundance description or grouping implication. Documentation must direct Condition questions to separately reviewed Sample-aware abundance/composition Modules.
7. Add explicit optional-environment guidance. Do not add the unrelated PyPI package as an extra.
8. Keep the node out of packaged examples until Linux/WSL/container installation and an exact real-backend test pass. No current packaged example contains the node, so this does not remove demonstrated behavior.
9. If the exact dependency cannot be made reproducible for the supported product environments, keep the ID as a fail-closed migration notice rather than registering a node that imports an arbitrary package and fails deep inside execution.

## Required tests

### Schema, registration, and migration

- node remains registered exactly once under the clustering category and is absent from abundance ownership;
- output order is `adata`, `summary`, `code`;
- saved legacy seed widget migrates correctly and unrelated graph links/widgets are preserved;
- README/install text never recommends `pip install schist`;
- examples, registration, and optional-dependency claims remain mutually consistent.

### Dependency/version behavior

- missing Schist, missing graph-tool, OSError, wrong PyPI package, missing `inference`, missing `fit_model`, v0.9.4, development-master signature drift, and unsupported future versions all fail actionably;
- exact v0.10.0 fake/real signature passes;
- the error reports detected module path/version and official installation route;
- no backend import happens after an already-determined graph/settings/memory failure if import ordering is not required for the check.

### Graph and transaction validation

- default and non-default named graph keys;
- missing/malformed `uns` metadata or pointer, missing `obsp` matrix, wrong shape, dense input, nonnumeric/negative/non-finite values, diagonal, asymmetry, no edges, isolates, duplicate observation names, empty/backed AnnData;
- CSR, CSC, and sparse-array inputs normalize to the same canonical CSR and fingerprint where mathematically identical;
- disconnected graphs succeed with a warning and exact component diagnostics;
- every collision in the owned `obs`/`obsm`/`uns` family fails before copying/backend work by default;
- overwrite clears stale higher-level family members on the output copy only;
- backend failure leaves caller state byte-for-byte unchanged.

### Exact backend call

- capture and assert every v0.10.0 keyword and value;
- nonzero seed only, fixed `n_jobs=1`, `collect_marginals=True`, unweighted undirected nested model, no constraint, no state file, and `copy=False`;
- the passed adjacency is the validated canonical CSR and the named graph metadata is reported;
- `posterior_samples < 100`, invalid degree-correction boolean, invalid namespace, invalid budget, or seed zero fail before fitting.

### Output postconditions

- consecutive finest-to-root categorical columns, all categories used, no missing cells, and cluster sizes conserve `n_obs`;
- nonincreasing cluster counts, exact parent mapping, and one-cluster root;
- malformed gaps, noncategorical labels, unused categories, wrong index/order, missing root, split parents, or extra owned keys fail;
- stats/params/blocks presence and finite aligned entropy/modularity diagnostics;
- one finite bounded marginal matrix per level with exact shape; malformed, missing, negative, greater-than-one, all-zero, or misaligned matrices fail;
- output graph is preserved, no new `obsp` result is invented, input X/layers/Raw/annotations remain unchanged, and H5AD round-trip succeeds.

### Resource, report, and code

- memory preflight runs before AnnData copy, coordinate materialization, or backend call;
- estimate covers AnnData copy, graph conversion, posterior partitions, worst-case marginals, hierarchy outputs, and reports without integer overflow;
- strict JSON serialization with `allow_nan=False`;
- ordered per-level records, graph fingerprint/diagnostics, exact parameters, completion-versus-convergence distinction, references, dynamic versions, platform/build, and limitations;
- no p-value, Condition effect, differential-abundance wording, selected/best level, or biological-truth claim;
- generated source compiles, uses `fit_model`, contains every explicit policy, preserves input, and matches runtime scientific state and failure categories with a fake backend;
- a small exact-v0.10.0 real Schist/graph-tool test runs only in the declared supported optional environment; otherwise it skips with the precise environment reason.

## Final cohesion and coupling assessment

The current node is shallow: it exposes one seed while making callers absorb a broken vendor call, implicit graph choice, package ambiguity, hierarchy semantics, destructive storage, and absent reporting. Moving it beside graph clustering corrects the domain model. The shared named-graph Seam provides Leverage with Leiden without merging distinct algorithms, while the private Schist Adapter localizes volatile version/backend details.

The proposed Module remains intentionally coupled to one reviewed Schist release and one unweighted undirected nSBM policy. That narrow compatibility is a strength: an unchanged workflow retains one scientific meaning. Future model variants or package upgrades deepen or extend the Adapter only after explicit evidence; they do not silently widen the public Interface.
