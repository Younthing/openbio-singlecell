# OpenBioSingleCellSCENICRegulonSpecificity — design review

## Decision

**Enhance and retain as one atomic descriptive RSS module.** Fix the reversed output semantics, consume the verified
SCENIC result artifact, and use an exact local implementation of pySCENIC 0.12.1's formula compatible with NumPy 2.
Do not merge with binarization: RSS compares continuous activity distributions with annotation membership, whereas
binarization estimates per-regulon on/off thresholds.

## Atomic interface

Inputs, in order:

1. `adata`, used only for aligned labels;
2. immutable `scenic_result` (`SCENICResultArtifact`);
3. `annotation_key`, default `cell_type`;
4. `annotation_status` (`unknown`, `provisional`, or `curated`), default `unknown`;
5. `max_output_rows` guard.

Outputs, in order:

```text
table, summary, code
```

The complete table has canonical columns
`group,regulon,rss,rank_within_group,group_cells,total_cells` and one row for every observed group-regulon pair.
Group order follows collision-safe observed categorical order; regulon order follows the artifact. Sorting for
display is group order, descending RSS, and regulon ID. Top-N is a bounded summary concern, not a table filter.

The implementation fixes the pySCENIC 0.12.1 formula and SciPy Jensen-Shannon distance definition. The obsolete
`pyscenic.rss` module is not a runtime dependency; pySCENIC remains a method/reference and an externally declared
version in provenance.

## Invariants and reporting

- Artifact observation IDs exactly equal AnnData observation IDs; activity rows are reordered only after set
  equality, duplicate, and nonblank checks.
- AUC values and computed RSS are finite. AUC must be nonnegative, and every regulon must have positive total
  activity; degenerate regulons fail with their exact names.
- Labels are complete, nonblank scalars after collision-safe normalization; at least one observed group is required.
  A single observed group is returned with a prominent warning that it has no between-group specificity
  interpretation; unused categorical levels are excluded and disclosed.
- The output grid, ranks, ranges, and orientation are independently verified. No silent clipping or replacement of
  undefined values is allowed.
- `summary` follows the strict report schema, provides report-ready top specificity patterns and complete accounting,
  calls the output descriptive cell-identity evidence, and includes dynamic OpenBio, NumPy, pandas, and SciPy
  versions plus external SCENIC provenance. JSON rejects NaN/Infinity.
- `code` returns `(table, summary_dict)` and reproduces exact alignment, formula, validation, sorting, and reporting.

## Scientific boundary

RSS describes concentration of per-cell regulon activity across annotations. It performs no permutation or
hypothesis test, does not account for Sample/Condition/Technical batch, and must not be worded as differential
regulation. Provisional annotations remain explicitly provisional.

## Migration and tests

Saved workflows connect `scenic_result`, rename `groupby` to `annotation_key`, add annotation status, remove
`activity_key`, and add `summary`/`code`. Existing cached tables are invalid because group and regulon meanings were
reversed.

Tests must cover hand-calculated RSS, exact orientation, multiple groups/regulons, categorical order and unused
levels, missing/colliding labels, zero-sum/negative/non-finite AUC, one-group warning, ID misalignment, full grid/ranks,
row guards, input immutability, no pySCENIC import, strict JSON, external provenance wording, and exact
runtime/generated equivalence.

## Cohesion assessment

Alignment, a complete mathematical definition, degenerate-case policy, ordering, provenance, and reporting sit
behind three scientific inputs. The module has one reason to change: the RSS definition and result contract.
