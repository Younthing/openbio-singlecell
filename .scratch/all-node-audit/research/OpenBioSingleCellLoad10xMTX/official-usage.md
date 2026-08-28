# Load 10x MTX: official usage research

Researched: 2026-08-28

## Official Scanpy usage

Scanpy documents `scanpy.read_10x_mtx(path, *, var_names="gene_symbols", make_unique=True, cache=False, gex_only=True, prefix=None, compressed=True)`. Its PBMC example reads a feature-barcode directory and then makes gene-symbol variable names unique:

```python
import scanpy as sc

adata = sc.read_10x_mtx(
    "filtered_feature_bc_matrix/",
    var_names="gene_symbols",
    cache=False,
    gex_only=True,
)
adata.var_names_make_unique()
```

The official `make_unique` parameter applies to the variable index. It does not authorize changing cell barcodes.

- Official API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_mtx.html
- Official Scanpy PBMC example: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html

## Official 10x format

Cell Ranger's MEX output contains `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, and `features.tsv[.gz]`. Features are matrix rows, barcodes are columns, and each matrix entry is a UMI count. Cell Ranger v3+ feature rows provide gene/feature ID, display name, and feature type. A filtered matrix contains cell-associated barcodes; a raw matrix also contains background/non-cell barcodes.

- 10x MEX specification: https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/outputs/cr-outputs-mex-matrices

## Method/software references

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I et al. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Scientific contract

- Require exactly one compatible matrix, barcode table, and feature table, with dimensions consistent with the feature and barcode row counts.
- Construct `AnnData.X` as observations-by-variables, so the 10x feature-by-barcode matrix is transposed once into a sparse cell-by-feature matrix.
- Require finite, non-negative, integer-like UMI counts. This loader defines an explicit count-state source.
- Require non-empty, unique cell barcodes and never suffix or otherwise repair them. A duplicate barcode is ambiguous cell identity.
- Require stable feature IDs to be non-empty and unique. Gene symbols may repeat; when gene symbols are selected as the variable index, the optional uniqueness policy may suffix only that displayed index while retaining original IDs and symbols in `var`.
- `gex_only=True` keeps only `Gene Expression` features and excludes modalities such as Antibody Capture and CRISPR Guide Capture. If the feature-type column is absent in a legacy file, disclose that filtering could not be verified.
- Record the three input files and fingerprints, compression/layout variant, selected variable-name policy, feature-type filtering, original and retained dimensions, count validation, and naming repairs in source metadata.

## Open expert-boundary re-review (2026-08-28)

Matrix/feature/barcode dimensional agreement and finite non-negative integer-like 10x values remain format-semantic
requirements. In contrast, an all-zero or empty retained matrix, blank display labels, and duplicate barcode/feature
identifiers can still be represented by AnnData. They will load with explicit warnings and axis/count audit fields.
The adapter does not silently repair cell barcodes; `make_unique` continues to affect only the selected variable
index, matching the official Scanpy interface.
