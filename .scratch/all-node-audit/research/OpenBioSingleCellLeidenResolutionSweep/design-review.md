# OpenBioSingleCellLeidenResolutionSweep — design review

## Decision

**Keep and substantially enhance Resolution Sweep as a separate atomic diagnostic module. Do not merge it into single-resolution Leiden, and do not delete it.**

Its atomic question is: “How do graph partitions, sizes, and optimization robustness change over this declared resolution grid?” Single-resolution Leiden answers a later question: “Store one chosen partition.” They occur at different analyst decision points and have different outputs. Combining them would make one interface carry both exploratory comparison and final annotation semantics.

The current node should no longer be considered complete when it only sprays multiple columns into `obs`. Retain those aligned columns because downstream plots need memberships, but make the structured resolution table and report first-class outputs.

## Current interface problems

- comma parsing silently removes duplicates and accepts non-finite floats until the backend fails;
- compact resolution key formatting can collide;
- no preflight collision check, so partial mutation is possible before a later resolution fails;
- graph key and matrix invariants are delegated to Scanpy and not reported;
- flavor, objective, directedness, weighting, and iterations are hidden;
- one seed provides no stochastic-stability evidence;
- no cluster sizes, modularity, graph components, isolates, or adjacent-resolution agreement;
- no structured table, `summary`, or `code`;
- implementation is physically separated from single Leiden despite identical graph/backend logic.

## Proposed interface

Visible inputs:

- `adata`;
- `resolutions="0.25,0.5,1.0,2.0"` as a migration-friendly strict list.

Advanced inputs:

- `key_prefix="leiden"`;
- `neighbors_key="neighbors"`;
- `n_iterations=2`, accepting every integer under python-igraph semantics;
- `stability_repeats=5`;
- `random_seed=0`;
- `overwrite_existing=False`.

Hidden and disclosed fixed policy:

- igraph flavor;
- undirected weighted graph;
- modularity objective;
- deterministic repeat-seed derivation;
- no automatic selection or recommendation.

Outputs:

- annotated `adata`;
- `resolution_metrics` structured table;
- strict JSON `summary`;
- equivalent source `code`.

## Deep module and shared internal seam

The enhanced module hides strict list/key planning, graph resolution and validation, repeated stochastic optimization, cluster-size aggregation, ARI calculation, categorical storage, table creation, and reporting behind a focused interface. That is leverage: callers specify a graph and resolutions without learning Scanpy private storage or igraph objects.

Move the single Leiden and sweep implementations behind one internal clustering seam in the clustering module. The internal interface takes a validated graph contract plus resolved algorithm policy and returns immutable partition diagnostics. Both public nodes use the same adapter, making this a real internal seam with two consumers and concentrating backend/version changes for locality.

Do not merge neighbor construction into this module. The upstream representation and graph parameters are distinct scientific decisions and must remain inspectable. Do not expose a generic Python partition object or Scanpy private `NeighborsView`; those would leak implementation into the workflow interface.

## Structured table design

One row per sorted resolution is the stable comparison interface. Cell memberships remain in AnnData, while table fields cover cluster count/sizes, modularity, repeat-seed ARI, and adjacent-resolution ARI. Use `None` for the first adjacent comparison in JSON; avoid NaN/Infinity in the report.

No “best” column is produced. The module may warn about singleton-heavy, dominant-cluster, unstable, or one-cluster solutions, but it does not rank resolutions. Biological interpretability, marker evidence, Sample composition, and sensitivity to upstream representation require downstream review.

`stability_mean_ari` and `stability_min_ari` are explicitly within-resolution repeat-start diagnostics. `adjacent_resolution_ari` is an across-resolution sensitivity diagnostic. Keeping these names distinct prevents a shallow report that conflates two different questions.

## Cohesion and coupling

Resolution parsing, repeat runs, and diagnostics all contribute to one comparative analysis, so they are cohesive. The base partition for each resolution is stored; repeat partitions are hidden implementation state because exposing them would enlarge the interface without supporting a current consumer.

The module is intentionally coupled to the fixed graph and single backend policy. Every resolution reuses exactly those inputs. It is not coupled to expression X, Raw snapshot, marker testing, or annotation truth. Pooled-cell graph partitions may reflect a Sample or Technical batch and are not Condition inference.

The deletion test supports retention: without this module, every workflow would recreate resolution parsing, N repeated Leiden runs, metrics, collisions, and comparison tables. By contrast, deleting automatic best-resolution selection removes unsupported complexity and does not move any required behavior elsewhere.

## Summary and code contract

`summary` carries the ordered table records, complete per-resolution cluster-size mappings, graph/neighbor provenance, fixed backend policy, seed schedule, warnings, limitations, references, and runtime versions. Results prose is suitable for a methods/results report but calls resolutions analyst-reviewed granularities, not ground truth.

`code` returns `(adata, resolution_metrics)` using public packages. It validates before mutation, calls Scanpy with all fixed/advanced arguments explicit, computes the same ARIs and table, and matches runtime membership/diagnostic failures. Graph resolution precedes output collision checks, and every planned output key is rejected if it equals the resolved `neighbors_key`, even under overwrite, because that `uns` slot is reserved for graph metadata. Plugin history and timing are the only acceptable omissions.

## Failure and verification plan

- Shared graph-malformation tests with single Leiden: missing pointer/key, wrong shape, negative/non-finite weights, diagonal, asymmetry, no edges, duplicate observation names.
- Resolution parser tests: whitespace, order, duplicates, NaN/Infinity, negative values, zero, large grids with workload warning, and generated-key collision.
- Preflight all `obs`/`uns` collisions before any run; verify overwrite behavior is all-or-nothing. Test that runtime and generated code reject a planned key equal to `neighbors_key` under both overwrite settings without damaging graph metadata, and that the input/result graph remains resolvable downstream.
- Assert every partition is categorical, aligned, complete, and its cluster sizes sum to `n_obs`.
- Verify deterministic base memberships for fixed seed and unchanged graph.
- Compare pairwise stability and adjacent-resolution ARI with direct scikit-learn calculations; ensure label renaming does not change ARI.
- Verify table row order/schema, strict JSON summary, H5AD round-trip, and absence of any best-resolution field.
- Verify generated code produces the same annotated columns/table and fails identically on malformed graphs/resolutions.
- Workload disclosure: report the resolved total run count and warn for a large request without imposing an arbitrary fixed maximum.

For the open expert tool, remove the fixed run-count ceiling and expose workload only through summary/warnings; resource cancellation belongs to the execution host. The metrics table remains stable when `stability_repeats=1`: its stability mean/min/max columns are null, while adjacent-resolution ARI is still computed. Runtime and generated code share this behavior.
