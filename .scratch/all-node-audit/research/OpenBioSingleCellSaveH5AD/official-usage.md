# Save H5AD: official usage research

Researched: 2026-08-28

## Official AnnData usage

AnnData's native writer is `AnnData.write_h5ad(filename, compression=..., ...)`. It stores the
annotated matrix and its axis annotations, layers, embeddings, graphs, raw snapshot, and
unstructured metadata in the H5AD format. The official API supports no compression, `gzip`,
and `lzf`; optional third-party HDF5 filters can make files unreadable without the same
plugin and are inappropriate as a portable workflow default.

```python
adata.write_h5ad("study.h5ad", compression="gzip")
```

- Official writer API:
  https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.write.html
- AnnData on-disk format:
  https://anndata.readthedocs.io/en/stable/fileformat-prose.html

The corresponding official reader can open a written file in backed read-only mode. This is
sufficient for a bounded post-write check of shape and axis identity without loading the
whole expression matrix again.

- Official reader API:
  https://anndata.readthedocs.io/en/stable/generated/anndata.io.read_h5ad.html

Python documents `os.replace` as atomic on success when source and destination are on the
same filesystem. `NamedTemporaryFile` provides a securely created staging name; on Windows
the handle must be closed before AnnData/HDF5 reopens it.

- `os.replace`: https://docs.python.org/3/library/os.html#os.replace
- `NamedTemporaryFile`:
  https://docs.python.org/3/library/tempfile.html#tempfile.NamedTemporaryFile

## Method/software reference

Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data
matrices. *Journal of Open Source Software*. 2024;9(101):4371.
https://doi.org/10.21105/joss.04371

## Reproducible export contract

- Write with the official AnnData API; never reconstruct only selected slots.
- Stage in the final directory, close the HDF5 writer, flush the staged file to disk, reopen
  it read-only, and verify shape plus exact observation/variable identifiers before commit.
- Commit with atomic replacement. If validation or writing fails, retain an existing target
  byte-for-byte and remove the staging file.
- `overwrite=False` must never replace a concurrently created file.
- Expose the standard compression choice as an advanced storage parameter. Use portable
  built-in codecs only; scientific contents must not depend on compression.
