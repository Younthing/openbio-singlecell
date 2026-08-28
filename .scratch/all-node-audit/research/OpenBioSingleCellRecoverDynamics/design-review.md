# Recover Velocity Dynamics: module design review

## Decision

**Retain and deepen as one explicit dynamical-parameter fit.** It returns reusable fitted state; it does not estimate final dynamical velocities, build a graph, compute latent time, rank by cell group, or run CellRank.

Inputs: `velocity_state` at stage `moments` or a verified stochastic/deterministic `velocity_estimated` stage; visible `gene_selection` combo `velocity_genes|all`, default `velocity_genes`; advanced `n_top_genes=0` (0 means all selected); advanced `max_iter=10`; advanced `n_jobs=1`; advanced `max_dense_gib=2.0`; advanced `overwrite_existing=False`. `velocity_genes` requires a complete upstream mask and fails otherwise. For a positive `n_top_genes`, select deterministically by descending summed `Ms` with original feature order as the tie-break, matching the audited 0.3.4 backend before passing an explicit list. Keep assignment mode, basal/steady-state/scaling/time flags, add key, backend and arbitrary kwargs fixed to the audited policy.

Outputs: new `OPENBIO_VELOCITY_STATE` stage `dynamics_recovered`, `summary`, `code`.

When recovery consumes a steady-state estimate, first copy `r2` into the recovered-fit namespace when needed for
the later dynamical model, then remove the entire stale upstream `vkey` result family (`offset`, `offset2`, `beta`,
`gamma`, `qreg_ratio`, `r2`, `genes`, velocity/variance layers and parameter record). The recovered artifact owns
the fit and remains collision-free for a subsequent dynamical estimate with the same `vkey`.

Preflight selected genes, artifact/graph/layer fingerprints, estimated dense output and collisions. After execution verify the complete `fit_*` schema, aligned shapes, finite-or-documented-failed genes, finite nonnegative (not probability-bounded) `fit_likelihood`, selected/failed counts, and `uns["recover_dynamics"]` parameters. Fail if scVelo internally creates a selection/model that conflicts with the request.

The locked Pandas 3.0.5 runtime exposes two upstream 0.3.4 incompatibilities: `make_unique_list` forwards a Python list to `pandas.unique`, and `_read_pars` returns read-only `Series.values` arrays that `align_dynamics` mutates. Keep the public scVelo call and exact gene family, but surround it with a lock-protected compatibility context that verifies and temporarily replaces only those two callable globals with stable-unique object-array and writable-copy semantics, restoring both unconditionally. Record the two compatibility actions in parameters/summary and test restoration on success and failure; do not patch Pandas globally, alter scientific values, or silently choose another gene set.

`summary` contains exact selected gene identity/count/fingerprint, convergence/failure counts, fit-likelihood and kinetic-parameter quantiles, iterations/jobs/memory, warnings, references and versions. Never serialize cell-by-gene matrices into JSON. `code` returns portable `(fitted_adata, summary_dict)` with identical validation/call/postconditions.

Migration maps the old node to `gene_selection="velocity_genes"`; workflows without a verified mask must insert an upstream stochastic estimate or explicitly choose `all`. Add a separate future `OpenBioSingleCellVelocityLatentTime` only after its own two research documents; it must consume a dynamical velocity-graph state and require explicit root/end evidence.

Tests cover stage/selection combinations, absent/corrupt mask, all/top-N selection, memory/collision guards, partial fit failures, backend fallback, deterministic ordering, defensive artifact ownership, strict JSON, real 0.3.4 smoke, compilation and equivalence.
