# Map Gene IDs from GTF: official usage research

Researched: 2026-08-28

## Official 10x and Scanpy feature identity

The 10x Cell Ranger MEX format stores feature ID and feature name separately in the first and second columns of `features.tsv.gz`. For Gene Expression features, those values come from GTF `gene_id` and `gene_name`, respectively. Stable identity and human-readable symbol are therefore distinct fields.

- 10x Cell Ranger feature-barcode MEX format: https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/outputs/cr-outputs-mex-matrices
- 10x Cell Ranger HDF5 matrix format: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/outputs/cr-outputs-h5-matrices

`scanpy.read_10x_mtx` defaults to `var_names="gene_symbols"`; the alternate index is `gene_ids`. The official PBMC example shows a symbol-indexed AnnData with `var: 'gene_ids'`. A GTF mapper that assumes `var_names` contains Ensembl IDs will therefore fail on the normal default loader output.

- Scanpy `read_10x_mtx`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_mtx.html
- Scanpy PBMC example showing `gene_ids` in `var`: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html

## Official GENCODE and Ensembl identity semantics

GENCODE GTF gene records carry `gene_id` and `gene_name`. Gene/transcript IDs include numeric versions from GENCODE v7 onward. IDs for chrY pseudoautosomal-region copies may also carry `_PAR_Y`, for example the shape `ENSG....7_PAR_Y`. Removing everything after the first period incorrectly discards `_PAR_Y` and can collapse distinct identities.

- GENCODE GTF data format: https://www.gencodegenes.org/pages/data_format.html

Ensembl distinguishes stable IDs from gene names, which may change as annotation improves. Stable IDs should remain the principal feature identity; GTF `gene_name` is an annotation and can legitimately repeat.

- Ensembl stable IDs: https://mart.ensembl.org/info/genome/stable_ids/index.html

The safe reconciliation order is:

1. match the declared input gene ID exactly to GTF `gene_id`;
2. only after exact failure, remove a terminal numeric version while retaining `_PAR_Y`, e.g. `re.sub(r"\.\d+(?=_PAR_Y$|$)", "", identifier)`;
3. count version-normalized matches separately because they indicate a possible matrix/GTF release mismatch;
4. never resolve conflicting normalized IDs with first-wins behavior.

## Scientific practice and references

Mapping should annotate identity, not silently filter expression features. Unmapped features remain present with a `gtf_mapped` flag; feature filtering is a separate transformation. Multiple stable IDs may share a gene symbol, so symbol duplication is reportable but not an error when stable IDs remain unique.

Relevant references:

- Frankish A et al. GENCODE 2025. *Nucleic Acids Research*. 2025. https://doi.org/10.1093/nar/gkae1078
- Martin FJ et al. Ensembl 2025. *Nucleic Acids Research*. 2025. https://doi.org/10.1093/nar/gkae1071
- Virshup I et al. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

Report input features; exact, version-normalized, unmapped, and ambiguous counts/rates; duplicate gene-symbol count while stable IDs remain unique; GTF path fingerprint/provider/release/genome-build when discoverable; malformed/missing-attribute lines; existing gene-name agreement/fill/conflict counts; and a limitation when annotation release/build metadata is missing.

## Open expert-boundary re-review (2026-08-28)

GTF parsing ambiguity is a true mapping error, but match coverage and metadata confidence are not. Zero matches,
duplicate/version-colliding declared IDs, existing gene-name disagreements, empty axes, and an already-present Raw
snapshot can all be represented without changing expression alignment. The node will annotate each feature
positionally, preserve existing conflicting names, report identity/coverage advisories, and leave Raw untouched. It
does not assert that current and Raw feature namespaces must match. Missing/blank declared IDs remain hard because no
declared stable identity exists for those positions.
