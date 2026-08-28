# Snapshot Expression: official usage research

Researched: 2026-08-28

## Official AnnData usage

AnnData defines `raw` as a stored version of `X` and `var`. The official initialization is:

```python
adata.raw = adata.copy()
```

Observation slicing also slices `raw`, while variable slicing leaves `raw` unchanged. A feature removed from current `X` can consequently still be recovered from `raw`. This behavior fits the repository's domain term **Raw snapshot**: the post-QC, full-gene count state retained before highly variable gene selection.

- AnnData raw API and slicing behavior: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html
- AnnData object and aligned slicing: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html

AnnData layers hold matrices on the same current observation and variable axes as `X`. A canonical `layers["counts"]` therefore remains aligned to the current object and is automatically narrowed by later feature subsetting; `raw` retains the full feature axis that existed when the snapshot was created.

- AnnData layers API: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.layers.html

When counts are already in a layer rather than `X`, the equivalent official-storage pattern uses a temporary copy so that the stored `raw.X` contains the selected counts while the returned object's current `X` remains unchanged:

```python
output = adata.copy()
counts = output.layers["input_counts"].copy()
output.layers["counts"] = counts.copy()
raw_source = output.copy()
raw_source.X = counts.copy()
output.raw = raw_source
```

`raw` and layers are storage representations, not proof that the matrix is raw UMI count data or that all genes are present. Those scientific preconditions require validation, workflow provenance, and a limitation in the report.

## Scientific practice

The repository's intended sequence is: load and establish stable feature identity, perform QC and accepted cell/gene filtering, create the Raw snapshot, then normalize and/or subset to highly variable genes. Saving a log-normalized or scaled matrix as `raw` would make later count-dependent analysis incorrect. Saving after HVG subsetting would defeat `raw`'s full-gene recovery role.

Integer UMI counts are the common case, but some scientifically legitimate corrected-count matrices are non-integer. Finite, non-negative, non-zero input is required. Non-integer values should be retained with a warning and explicit matrix-type limitation rather than universally rejected as invalid.

## Method and software reference

This is a provenance/storage checkpoint rather than a statistical method. Cite:

- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

The report should state the selected source; cell/feature dimensions; total and nonzero counts; sparse/dense storage and dtype; whether values are integer-like; canonical layer and Raw snapshot creation; overwrite status; and whether repository provenance confirms creation before HVG subsetting. It must explicitly disclose that an arbitrary external H5AD cannot prove “full-gene” status from matrix values alone.

## Open expert-boundary re-review (2026-08-28)

AnnData's public storage operation accepts any aligned matrix, including empty, negative, non-finite, or non-integer
values and non-unique axes. These properties affect whether the snapshot should be interpreted as UMI counts, but the
expert's explicit source selection is authoritative for this storage checkpoint. The node will copy the selected `X`
or layer and report a machine-readable expression audit; it will not inspect analysis history, prove full-gene state,
or bind current content to an earlier Raw state. Existing destinations still require explicit overwrite to avoid
irreversible replacement.
