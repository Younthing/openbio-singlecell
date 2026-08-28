# Load H5AD: official usage research

Researched: 2026-08-28

## Official AnnData usage

AnnData's native reader is `anndata.io.read_h5ad(filename, backed=None, *, as_sparse=(), as_sparse_fmt=csr_matrix, chunk_size=6000)`. The ordinary call loads the complete object into memory:

```python
import anndata as ad

adata = ad.read_h5ad("study.h5ad")
```

`backed="r"` is an optional file-backed mode, but the official documentation notes that backed mutation support is limited: only changes to `X` are written in place, while changes to slots such as `obs` must be saved to a new file. This workflow expects ordinary mutable `AnnData`, so the loader should deliberately use memory mode rather than exposing a partly compatible object.

- Official reader: https://anndata.readthedocs.io/en/latest/generated/anndata.io.read_h5ad.html
- Official `AnnData` structure: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
- On-disk format specification: https://anndata.readthedocs.io/en/stable/fileformat-prose.html

The official `var_names_make_unique()` operation appends numeric suffixes to duplicated variable names while leaving the first occurrence unchanged. It is a naming repair, not gene-identity inference.

## Method/software reference

Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Scientific contract

- H5AD is a general annotated-data container, not proof that `X` contains raw counts. A loader must preserve the stored expression state and must not relabel normalized, logged, scaled, imputed, or count data.
- Load all supported AnnData slots in memory and require at least one observation and one variable.
- Observation identifiers represent cells and must be non-empty and unique. Ambiguous cell identity is an error; never repair it silently.
- Variable identifiers should be unique for reliable downstream indexing. The optional repair may suffix duplicate variable names, but the action and duplicate count must be recorded in OpenBio source metadata. Disabling repair preserves the file exactly and must disclose any duplicates as a warning.
- Preserve `raw`, layers, axis annotations, embeddings, graphs, and unstructured metadata exactly except for the explicitly requested variable-index repair and OpenBio provenance metadata.
- Record the input-relative path, stable file fingerprint, original shape, index-uniqueness state, load mode, and naming policy in provenance.
- Do not apply count validation: the format can legitimately contain non-count analysis states.

## Open expert-boundary re-review (2026-08-28)

The official AnnData reader accepts any structurally valid H5AD, including objects with empty axes or non-unique
indexes. Those states can limit later name-based analyses, but they are not H5AD decoding failures. The loader will
therefore load them and record bounded axis advisories instead of rejecting them. `make_var_names_unique=True`
remains an explicit repair choice; observation names are never rewritten automatically.

Embedded `uns["openbio_singlecell"]` is file content, not trusted execution history. A malformed or older embedded
metadata block must not prevent a valid H5AD from loading. Compatible metadata is preserved; incompatible metadata
is quarantined into sanitized source provenance and replaced by a fresh local metadata envelope with a warning.
