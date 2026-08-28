# Core Study Parameters: authoritative usage research

Researched: 2026-08-28

## Authoritative repository contract

This is an OpenBio configuration adapter rather than a wrapper around an external scientific method. Its official interface is the versioned `CoreStudyParameters` schema in `openbio_singlecell/study_inputs.py`, interpreted using the domain vocabulary in `CONTEXT.md`.

The canonical payload is one strict JSON object:

```json
{
  "schema_version": 1,
  "sample_column": "sample",
  "condition_column": "group",
  "batch_column": "",
  "annotation_column": "cell_type",
  "reference": "control",
  "comparison": "treated"
}
```

The parser requires exactly these seven fields, integer `schema_version: 1`, string values, and already-trimmed text. Missing or unexpected fields, a non-object top level, non-string values, non-finite JSON constants, and leading/trailing whitespace are errors. It emits the six values after `schema_version` as reusable strings.

## Scientific meaning

- `sample_column`: observation column containing **Sample**, the independent biological specimen and replicate unit.
- `condition_column`: observation column containing the biological **Condition**.
- `batch_column`: optional observation column containing **Technical batch**, which must not be substituted for Sample or Condition.
- `annotation_column`: observation column containing the population annotation selected for downstream analysis; provisional and curated annotations remain distinct.
- `reference` and `comparison`: levels of the Condition column defining an oriented Condition contrast, not column names.

The node validates configuration syntax only because it receives no dataset. Downstream consumers must validate that configured columns exist, contain no missing values where required, and contain the requested levels.

## Reporting reference

Füllgrabe A et al. Guidelines for reporting single-cell RNA-seq experiments. *Nature Biotechnology*. 2020;38:1384-1386. https://doi.org/10.1038/s41587-020-00744-z

The reporting guideline motivates explicit disclosure of experimental design, Samples, biological conditions, technical factors, annotations, and comparisons. It does not define this repository-specific JSON schema.

## Open expert-boundary re-review (2026-08-28)

The versioned JSON shape and string types remain the adapter's serialization contract. String contents are user-owned
column names and levels, however, and pandas can legitimately contain leading or trailing whitespace. The parser will
therefore preserve every string exactly rather than rejecting whitespace. Dataset consumers remain responsible for
checking whether the exact declared label exists.
