# Load 10x H5: official usage research

Researched: 2026-08-28

## Official Scanpy usage

Scanpy documents `scanpy.read_10x_h5(filename, *, genome=None, gex_only=True, backup_url=None)`. It returns an `AnnData` with the count matrix in `X`, barcodes in `obs_names`, feature names in `var_names`, gene IDs and feature types in `var`, and any additional `/matrix/features` metadata.

```python
import scanpy as sc

adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5", gex_only=True)
adata.var_names_make_unique()
```

The official Scanpy multi-Sample tutorial uses this pattern before concatenation. For legacy 10x H5 files containing more than one genome, the API requires a `genome` selection.

- Official API: https://scanpy.readthedocs.io/en/stable/api/scanpy.read_10x_h5.html
- Official tutorial: https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering.html

## Official 10x format

10x documents the HDF5 feature-barcode matrix as a top-level `matrix` group containing matrix values/index arrays, dimensions, barcodes, and feature metadata. H5 is the binary counterpart of MEX and represents feature-by-barcode UMI counts.

- 10x HDF5 format: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/outputs/cr-outputs-h5-matrices

## Method/software references

- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I et al. anndata. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Scientific contract

- Accept only a local 10x feature-barcode `.h5`/`.hdf5` file and call the pinned Scanpy reader without network fallback.
- Require a non-empty matrix with finite, non-negative, integer-like UMI counts.
- Require non-empty unique cell barcodes and never repair them. Cell identity errors must be visible.
- Preserve the original gene IDs, feature types, genomes, and other feature metadata returned by Scanpy. Duplicate displayed feature names may be suffixed only under the explicit variable-name policy.
- `gex_only=True` deliberately excludes non-Gene-Expression modalities. `genome` is a legacy-format selector, not a substitute for gene-annotation mapping.
- Record file fingerprint, original and retained dimensions, selected genome, feature filtering, count validation, and variable-index repair in source metadata.

## Open expert-boundary re-review (2026-08-28)

The official Scanpy reader returns the file's AnnData representation and exposes `gex_only` as feature selection.
After successful decoding, empty retained axes, all-zero matrices, and non-unique/blank axis labels are usability
advisories rather than reasons to reject a representable object. Finite non-negative integer-like values remain hard
because this adapter explicitly declares the payload as a 10x count matrix. Cell identifiers are never repaired;
variable-index suffixing occurs only when the expert selects `make_unique=True`.
