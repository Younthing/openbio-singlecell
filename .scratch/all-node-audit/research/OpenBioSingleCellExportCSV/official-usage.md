# Export CSV: official usage research

Researched: 2026-08-28

## Official pandas usage

Pandas exports a DataFrame with `DataFrame.to_csv`. Its documented default is `index=True`,
because the index is an axis of the table rather than disposable display chrome. The API
also exposes `index_label`, `encoding`, and `lineterminator` for an explicit file contract.

```python
frame.to_csv(
    "results.csv",
    index=True,
    index_label="__openbio_row_id__",
    encoding="utf-8",
    lineterminator="\n",
)
```

- Official `DataFrame.to_csv` API:
  https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html

CSV has no standard representation for pandas dtypes, categories, index metadata, or the
difference between every possible missing/string value. It is therefore a human-readable
exchange/report format, not an equivalent replacement for H5AD or Parquet.

Python documents same-filesystem `os.replace` as atomic on success, and secure temporary
files can be created in the destination directory with `NamedTemporaryFile`.

- `os.replace`: https://docs.python.org/3/library/os.html#os.replace
- `NamedTemporaryFile`:
  https://docs.python.org/3/library/tempfile.html#tempfile.NamedTemporaryFile

## OpenBio CSV contract

- Preserve every row and column in order and always serialize the DataFrame index.
- Require a single-level column axis. Pandas can emit multi-row CSV headers for MultiIndex
  columns, but that representation is not an unambiguous one-header exchange table; callers
  must flatten and name those columns explicitly before export.
- Use the existing non-empty index name when it is unambiguous. Give unnamed index levels
  deterministic reserved labels. Reject labels that collide with data columns rather than
  emit duplicate headers.
- Use UTF-8, comma delimiter, LF line endings, pandas' documented quoting behavior, and a
  fixed `NA` missing-value representation. CSV remains lossy because a literal string `NA`
  cannot be distinguished from a missing value without an external schema.
- Stage, flush/fsync, and atomically commit. A failed export must not truncate an existing CSV.
- `overwrite=False` must never replace a file created between name selection and commit.

No statistical method is performed by this adapter, so no method citation applies.
