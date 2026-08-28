# Merge Observation Annotations: module design review

## Current module

The node transfers one `obs` column from a second AnnData by shared `obs_names`, copies the target object, and leaves its observation and variable sets unchanged. Requiring unique indices and using label-based lookup are correct foundations.

## Correctness and interface problems

- `overwrite=False` silently ignores conflicting non-missing target values; the caller cannot distinguish conflicts from identical values.
- Source observations absent from the target are silently discarded even though `subset_adata` is expected to derive from the target lineage.
- Exact barcode equality is treated as sufficient provenance.
- `prefix` performs a second string-transformation concern and coerces non-string annotations.
- The implementation converts the target to object dtype, potentially losing categorical semantics without disclosure.
- It does not report overlap, conflicts, assignments, or final category counts.

## Decision: keep and enhance; do not merge with subsetting

Keep this node as the single annotation-transfer module. `SubsetObservations` changes which observations exist; this node preserves all target observations and changes one annotation column. Combining them would mix axis selection with index-aligned data transfer and obscure both error modes.

## Target interface

- Visible: target `adata`, source `subset_adata`, `source_column`, `target_column`.
- Advanced: `conflict_policy = error | keep_target | overwrite`, default `error`.
- Fixed hidden policy: source missing values are skipped; source-only observation identifiers are rejected; target-only observations remain unchanged.
- Outputs: annotated target copy, `summary`, `code`.

The current `overwrite` boolean should be replaced because it conflates “fill missing” with “ignore conflicts.” For saved-workflow compatibility, a temporary shim may map legacy `overwrite=False` to `keep_target` and `True` to `overwrite`, while new workflows receive `error` by default. The legacy `prefix` input should be removed after migration. If retained for one cycle, require a string/categorical source, apply it before conflict classification, and report the type conversion.

## Invariants and conflict classification

- Require unique `obs_names` in both inputs; never synthesize uniqueness.
- Require nonempty source/target column names and an existing source column.
- Require at least one shared observation and reject every source-only identifier by default. A caller needing intersection behavior should first use an explicit subset step.
- If both inputs carry compatible OpenBio source provenance, verify it. If provenance is absent or insufficient, proceed with a limitation warning; if it demonstrably conflicts, fail.
- Do not require equal variable axes because full-gene and HVG branches are an intended use case.
- Align by `obs_names`, independent of source order.
- Classify each shared, non-missing source value as: fills missing target; identical to target; or conflicts with target. Apply the selected policy only to the conflict class.
- Preserve the source dtype when creating a new target column. When merging categorical columns, explicitly union compatible categories or report a justified dtype conversion.
- Preserve target `X`, layers, raw, axes, and every unrelated annotation exactly.

## Report contract

`summary` includes both input shapes; source/target/shared/source-only/target-only counts; source missing, fill, identical, conflict, overwrite, and keep-target counts; source/target columns; conflict policy; output non-missing and category counts; dtype before/after; provenance evidence and limitations; AnnData/pandas references; and dynamic software versions. The methods/results text must say “transferred an annotation by exact observation identity,” not “validated cell types.”

`code` implements the same index validation, source-only check, alignment, conflict classification, dtype handling, policy, and copy semantics using public AnnData/pandas operations.

## Parameter policy and depth

Column names and the two data inputs define the transfer and stay visible. Conflict policy is advanced but result-defining and must be explicit. Prefix formatting, identifier repair, automatic intersection, and variable-axis reconciliation are excluded. This concentrates complex alignment and audit logic behind one small interface and keeps observation selection in its own module.

## Verification plan

- Shuffled source order proves label alignment rather than positional assignment.
- Duplicate identifiers on either side, no overlap, and source-only identifiers fail with actionable errors.
- Source missing values, target missing fills, identical values, and all three conflict policies produce separately counted outcomes.
- Categorical inputs preserve/union categories or disclose a deliberate dtype conversion; legacy prefix compatibility rejects non-string misuse.
- Compatible, absent, and demonstrably conflicting provenance follow the stated verify/warn/fail behavior.
- Target `X`, layers, `raw`, axes, and unrelated annotations are unchanged; inputs are immutable.
- Summary is strict JSON and generated code compiles and returns the equivalent annotation column.

## Open expert-boundary decision (2026-08-28)

Relax exact source-provenance equality and subset-only lineage assumptions to audit statuses. The implementation
aligns the intersection by unique `obs_names`, ignores source-only rows, and records zero-overlap/all-missing no-op
outcomes. This concentrates deterministic alignment at the module seam while leaving biological-source confidence to
the report. Runtime and generated code share the same behavior; no provenance-tamper gate is tested.
