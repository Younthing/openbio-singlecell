# Estimate RNA Velocity: module design review

## Decision

**Retain and deepen as one explicit kinetic-model estimation.** It consumes a validated staged artifact and returns a new one; it does not compute moments, recover dynamics, build transition graphs, latent time or CellRank fates.

Target inputs: `velocity_state`; visible `mode` combo `stochastic|deterministic|dynamical`, default `stochastic`; advanced `vkey="velocity"`; advanced `min_r2=0.01`; advanced nonnegative and unbounded `min_likelihood=0.001` used only for dynamical (likelihood density is not a probability); advanced `overwrite_existing=False`. Keep offsets, grouping, raw use, latent-time regularization, gene filtering and arbitrary kwargs hidden at audited values.

Require stage `moments` for deterministic/stochastic and stage `dynamics_recovered` for dynamical. Validate all model-specific fields before the call. Capture warnings and fail if requested mode differs from `uns[f"{vkey}_params"]["mode"]` or if fallback occurs. Verify an aligned velocity layer whose selected velocity-gene columns are finite, a complete boolean velocity-gene mask, finite parameter columns, and a minimum evidence floor (report/fail when fewer than 10 genes are selected, matching scVelo's own severe warning). Versioned 0.3.4 source initializes non-fitted dynamical-gene columns to `NaN`; those non-selected values may remain missing only when counted and disclosed, rather than being falsely reported as a wholly finite matrix.

Treat every exact 0.3.4 steady-state parameter output (`offset`, `offset2`, `beta`, `gamma`, `qreg_ratio`,
`r2`, `genes`) plus velocity/variance/auxiliary layers and the parameter record as one reserved `vkey` family for
collision, overwrite and artifact fingerprints. Snapshot and restore unrelated `obs` and `var` frames around the
backend's unconditional categorical conversion, then reattach only the generated `vkey` columns. In dynamical
mode, additionally require the recovered-fit fingerprint to remain identical after estimation.

For dynamical mode under locked Pandas 3.0.5, scope a lock-protected `compute_divergence` compatibility wrapper that copies only read-only ndarray inputs before the audited function mutates them. Verify the exact helper signature, restore the original even on failure, record the action, and leave deterministic/stochastic paths untouched.

Outputs: new `OPENBIO_VELOCITY_STATE` at stage `velocity_estimated`, `summary`, `code`. The artifact records mode, fixed settings, vkey, upstream fingerprints, selected-gene fingerprint and output-layer fingerprint.

`summary` includes model, prerequisites, selected/failed gene counts, R2/likelihood distributions as applicable, velocity magnitude quantiles, warnings, references and versions. It prominently labels the result exploratory/model-dependent. `code` returns `(velocity_adata, summary_dict)` with identical public scVelo call and fail-closed fallback detection.

Migration preserves explicitly selected old mode; old default values become `deterministic` to preserve intent and receive a notice that new workflows default to stochastic. Tests cover every stage/mode combination, missing/corrupt fit fields, fallback interception, full-family collisions, too-few genes, finite output checks, unrelated annotation preservation, input ownership, strict JSON, real 0.3.4 smoke, compilation and equivalence.
