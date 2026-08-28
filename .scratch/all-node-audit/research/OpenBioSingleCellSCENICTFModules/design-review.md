# OpenBioSingleCellSCENICTFModules — design review

## Decision

**Enhance and redefine as an atomic final-regulon membership report over `SCENICResultArtifact`.** Remove the raw
adjacency/expression recomputation. The old behavior is pre-cisTarget candidate-module generation and is not the
final SCENIC result users expect. Do not merge this report into Import: the importer validates and materializes the
complete artifact, while this node supplies a focused, filterable scientific table and report.

## Atomic interface

Inputs, in order:

1. immutable `scenic_result` (`SCENICResultArtifact`);
2. optional exact `transcription_factor`, blank for all TFs;
3. `max_output_rows` guard.

Outputs, in order:

```text
table, summary, code
```

The canonical table columns are
`regulon,transcription_factor,regulation,context,target,target_weight,motif_evidence_count,motif_ids,target_rank`.
Rows come only from motif-pruned final regulons imported from `ctx`, not from recomputed candidate modules. Ordering
is complete artifact regulon order then descending target weight and target ID. A TF filter changes the returned
view and is disclosed; it never renumbers module/regulon identities.

No AnnData, expression source, raw adjacency, file path, arbitrary module ID, or pySCENIC runtime is part of the
interface. External resource/run provenance is carried by the artifact.

## Invariants and reporting

- Artifact edge fields are strict nonblank scalars; weights are finite; regulon-to-TF/sign/context identity is
  consistent; target duplicates and motif evidence were normalized at import.
- An explicit TF must exist or return an empty, schema-stable table with an actionable warning according to the
  repository empty-result policy; the complete-family counts remain in the summary.
- Output rows are a defensively owned exact subset of artifact membership. No resource read, expression read,
  correlation, motif test, or mutation occurs.
- `summary` follows the strict report schema, gives bounded report-ready regulon/target results and complete/query
  accounting, preserves external resource/container provenance, distinguishes inferred membership from binding,
  and emits dynamic OpenBio/pandas versions. JSON rejects NaN/Infinity.
- `code` returns `(table, summary_dict)` and reproduces artifact validation, exact filtering, ordering, and reporting.

## Scientific boundary

This is a reporting/projection node, not a new network inference or hypothesis test. Target weights and motif support
are method evidence, not causal effect estimates. It performs no Sample/Condition/Technical batch analysis.

## Migration and tests

Saved workflows remove `adata`, `network`, and expression `source`; connect `scenic_result`; keep the exact TF query
if one existed; and add `summary`/`code`. Old cached module tables cannot be migrated as final regulons because motif
pruning information was absent. The semantic change must be explicit in workflow migration notes.

Tests must cover all-TF and one-TF views, multiple regulons per TF, activation/repression context, target weights and
ties, motif evidence, absent TF behavior, schema-stable empty results, row guards, defensive ownership, no
pySCENIC/file/expression access, strict JSON, provenance/limitation wording, and exact runtime/generated equivalence.

## Cohesion assessment

Although projection is intentionally simple, the node earns its interface by presenting stable final-regulon
semantics, query behavior, provenance, ordering, and report-ready interpretation. Import remains the security and
normalization seam; this module remains the scientific reporting seam.

