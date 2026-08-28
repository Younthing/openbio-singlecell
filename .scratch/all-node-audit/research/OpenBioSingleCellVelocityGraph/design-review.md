# Velocity Graph: module design review

## Decision

**Retain and deepen as one velocity-vector-to-directed-correlation-graph module.** Its output remains a staged velocity artifact, not a CellRank transition kernel or estimator.

Inputs: `velocity_state` at stage `velocity_estimated`; advanced `vkey="velocity"`; advanced `xkey="Ms"` fixed/bounded to a verified artifact layer; advanced `mode_neighbors="distances"` combo `distances|connectivities`; advanced `n_jobs=1` (portable/reproducible default); advanced `overwrite_existing=False`. Keep approximate projection, recursion, random subsampling, uncertainty and arbitrary gene subsets hidden until separately researched. Fix `n_recurse_neighbors=1` so “named graph support” means direct edges rather than the 0.3.4 distance-mode default's implicit two-hop expansion, and fix `sqrt_transform=False` so policy does not change silently by velocity model.

Outputs: new `OPENBIO_VELOCITY_STATE` at stage `velocity_graph`, `summary`, `code`. Validate the upstream model/layers/mask/named neighbor graph and output collisions. Call scVelo explicitly, capture warnings, and verify two observation-aligned finite sparse matrices with zero diagonal, disjoint signs, direct-edge support contained in the selected named matrix, and exact parameter metadata. Record fingerprints for both directed matrices and all upstream state.

`summary` reports positive/negative nnz and weight quantiles, in/out degree and zero-degree counts, source graph/model/gene counts, compute settings, methods/references, versions and limitations. `code` returns portable `(graph_adata, summary_dict)` with the same checks.

## CellRank seam and migration

Do not add a `cellrank=True` mode. A future `Build CellRank Kernel` node should consume this artifact and emit a distinct immutable `OPENBIO_TRANSITION_MATRIX` containing a row-stochastic matrix, direction, kernel expression/weights, obs fingerprint and CellRank version. A future GPCCA node should consume that artifact and emit a separate estimator artifact. CellRank 2's official split between kernels and estimators is a real seam with multiple adapters, not optional internals of scVelo.

Old workflows map `adata` to the staged artifact chain and retain `n_jobs`; missing stages receive inserted preparation/moment/velocity nodes or block. Tests cover stage/fingerprint corruption, collision handling, sparse matrix postconditions, zero-degree disclosure, n_jobs, no input mutation, strict JSON, 0.3.4 smoke, source compilation and equivalence.
