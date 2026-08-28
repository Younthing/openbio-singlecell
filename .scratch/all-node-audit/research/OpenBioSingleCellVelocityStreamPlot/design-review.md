# Velocity Stream Plot: module design review

## Decision

**Retain as a read-only plot adapter over a verified velocity-graph artifact.** It does not belong in the computational graph node because figures have separate styling, mutation and output semantics.

Inputs: `velocity_state` stage `velocity_graph`; visible `basis="umap"`; visible `color_key="leiden"` (blank means no color); advanced `density=2.0`; advanced `smooth=0.5`; advanced `min_mass=1.0`; fixed vkey from the artifact. Outputs: `plot`, `summary`, `code`.

Validate basis/color/artifact fingerprints and work on a defensive AnnData copy. Use an explicitly owned Matplotlib figure/axis, `show=False`, and verify PNG bytes. Any temporary `velocity_<basis>` cache remains on the private copy. Do not mutate global scVelo/Matplotlib settings.

`summary` includes basis/color semantics, graph/model identity, vector/embedding diagnostics, grid parameters, warnings, references and package versions. Results prose says a visualization was rendered, not that a lineage/fate was discovered. `code` returns a portable `(figure_or_png, summary_dict)` and reproduces the exact styling and validation.

Migration maps old `groupby` to `color_key`. Tests cover missing/malformed basis/graph/color, blank/numeric/categorical color, private cache behavior, no input/global-setting mutation, deterministic rendering tolerance, strict JSON, code compilation and equivalent figure metadata.
