# Velocity Moments: module design review

## Decision

**Retain and deepen as a single prepared-state plus named-graph to moment-state transformation.** Keep graph construction in the existing Neighbors node and keep kinetic fitting downstream.

Inputs: `velocity_state` at stage `prepared`; `adata_with_graph` only if the artifact does not already own the selected graph; visible `neighbors_key="neighbors"`; advanced `mode="connectivities"` as a bounded combo `connectivities|distances`; advanced `max_dense_gib=2.0`; advanced `overwrite_existing=False`. Remove `n_pcs` and `n_neighbors` from this node because changing them would authorize hidden graph reconstruction. Prefer storing the validated graph in the velocity artifact during preparation or binding it explicitly here by exact obs fingerprint.

Outputs: a new `OPENBIO_VELOCITY_STATE` at stage `moments`, `summary`, and `code`.

Before scVelo, validate artifact stage/fingerprints, canonical normalized layers, graph matrices/metadata and obs alignment; reject existing `Ms`/`Mu` unless overwrite. Alias the requested graph on a private copy and call `scv.pp.moments(..., n_neighbors=None, n_pcs=None, use_rep=None, mode=..., copy=False)` so no graph is recomputed. Capture warnings and fail if scVelo announces normalization or neighbor calculation.

The artifact records graph key/fingerprint, mode, `Ms`/`Mu` fingerprints and dense byte count. `summary` reports graph diagnostics, moment shapes/dtypes/ranges/quantiles, dense-memory estimate, warnings, methods/references and versions. `code` takes portable prepared AnnData plus named graph, reproduces aliases/validation/call, and returns `(moment_adata, summary_dict)`.

Migration maps old zero/default dimension values to explicit graph reuse. Nonzero old values require inserting/reconfiguring a Neighbors node and receive a migration warning; they must not remain inert controls.

Tests include stage/type rejection, graph-name/fingerprint/obs mismatch, no-hidden-recompute assertion, mode handling, memory guard, collisions, malformed outputs, input/artifact immutability, fake/real 0.3.4 calls, strict JSON, compiled code and equivalence.
