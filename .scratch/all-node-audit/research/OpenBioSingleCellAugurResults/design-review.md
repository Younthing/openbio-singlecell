# OpenBioSingleCellAugurResults — design review

## Decision

**Retain as an explicit, intentionally narrow result adapter.** Replace the AnnData/`uns` interface with the typed
Augur artifact. Its shallowness is appropriate because it is an adapter at a real type seam, not a scientific module.
Do not merge it into Augur unless the UI can expose three first-class table outputs without losing a reusable result
artifact.

## Adapter interface

Inputs:

1. typed `result` from `OpenBioSingleCellAugur`;
2. `view`: `priorities`, `cross_validation`, or `feature_importance`.

Outputs:

```text
table, summary, code
```

`table` is a `TableResult` containing an immutable copy of the selected canonical DataFrame and inherited input/
producer/fingerprint metadata. `summary` preserves the complete Augur summary and appends only selected table name,
row/column counts, and content fingerprint. `code` defines a tiny equivalent selector returning
`(DataFrame, summary_dict)`. It performs no model call, ranking, filtering, aggregation, or scientific reinterpretation.

## Integrity and failure contract

Require exact artifact class/schema/producer/version, recompute all canonical table content fingerprints, validate
view-specific columns/dtypes/unique keys/ranges, and ensure the artifact's summary is strict JSON. Reject mutable
dicts, AnnData `uns`, unknown tables, malformed or tampered frames, and output beyond a configured preview/export
guard. Preserve full rows; this adapter selects a table, not a top-N subset.

## Sample and Technical batch semantics

All Sample, Condition, Technical batch, expression-state, annotation-status, and cell-level-CV limitations are
inherited unchanged. The adapter must not weaken them or describe a selected table as a new inference.

## Migration and tests

The migrator rewires only a direct old Augur-to-AugurResults link to the new artifact port and maps
`summary_metrics -> priorities`, `full_results -> cross_validation`, and `feature_importances -> feature_importance`.
Any arbitrary AnnData producer must stop with an actionable error.

Tests cover each exact view/schema, direct producer requirement, parent/table fingerprint tampering, immutable copies,
unknown/malformed views, row guards, parent summary/reference/limitation propagation, strict JSON, compiled selector
code, and exact runtime/generated table+summary equality. No Pertpy import should occur in adapter tests.

## Cohesion assessment

The adapter has one responsibility: convert one validated table inside AugurResult into the generic table interface.
Deleting it would push artifact knowledge into generic preview/export nodes, so the seam earns locality despite its
small implementation.

