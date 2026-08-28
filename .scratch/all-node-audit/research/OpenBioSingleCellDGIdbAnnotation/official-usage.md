# OpenBioSingleCellDGIdbAnnotation — official usage research

## Audited baseline and P0 dependency behavior

The current node calls `pertpy.md.Drug().annotate(adata, source="dgidb")`. Pertpy 1.3.0 lazily downloads
`https://exampledata.scverse.org/pertpy/dgidb.tsv`, groups `drug_claim_name` to `gene_claim_name`, and writes a
pipe-separated `adata.var["compounds"]` column. The node does not report the DGIdb release, derived-file creation
date, hash, contributing sources/licenses, organism/identifier namespace, unmatched genes, or the network actually
used downstream. It performs hidden network I/O when the cache is cold.

Pertpy is not declared in the repository dependency metadata even though the runtime currently has 1.3.0 installed.
That makes installation behavior non-reproducible. The sibling Drug Scores/Hypergeometric/GSEA nodes construct a new
`Drug` object and access `drug.dgidb.dictionary` before calling `set()`/`dict()`, so the attribute does not exist and
all three fail on a fresh object. Annotation happens to initialize only its own short-lived object and does not fix
those nodes.

## Official Pertpy 1.3.0 interface

Primary sources:

- https://pertpy.readthedocs.io/en/stable/api/metadata/pertpy.metadata.Drug.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/metadata/_drug.py
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/metadata/_metadata.py

The documented annotation call is:

```python
drug = pertpy.md.Drug()
drug.annotate(adata, source="dgidb", copy=False)
targets = drug.dgidb.dict()  # lazy-loads before returning the dictionary
```

`DrugDataBase.dict()` and `.df()` are the safe lazy accessors. Direct `.dictionary` access is invalid before
initialization. Pertpy documents that genes must be HGNC-compatible. Its cached derived TSV has no scientific
metadata contract, so it is unsuitable as the sole reproducibility record.

## Official DGIdb resource practice

- Versioned downloads: https://dgidb.org/downloads
- Data accessibility/licensing: https://dgidb.org/about/overview/data-accessibility
- DGIdb v5 publication: https://doi.org/10.1093/nar/gkad1040

DGIdb provides raw versioned data dumps. Its downloads page warns that some imported sources have redistribution
restrictions and directs users to source-specific license information. The MIT license applies to DGIdb software,
not automatically to every aggregated data source. A scientific adapter must therefore require a locally prepared
snapshot plus explicit release/source-license metadata and must never silently download “latest”. A string such as
`reviewed` is only the caller's attestation and cannot prove a legal audit; pending, unknown, or incomplete review
remains loadable for expert inspection but must trigger a prominent redistribution/license warning.

For the downstream OpenBio methods, the canonical resource is a long table of exact nonblank drug and HGNC gene
identifiers. Exact duplicate pairs may collapse, while source/evidence fields are retained for disclosure when
present. No drug-name normalization, gene alias mapping, therapeutic direction inference, or clinical ranking is
performed.

## References to emit

- Cannon M, et al. DGIdb 5.0: rebuilding the drug-gene interaction database for precision medicine and drug
  discovery platforms. *Nucleic Acids Research*. 2024;52:D1227-D1235. https://doi.org/10.1093/nar/gkad1040
- Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025.
  https://doi.org/10.1038/s41592-025-02909-7

## Report/code implications

The report must identify source/release/date/download URL, organism and gene namespace, caller-declared license and
review status without presenting either as programmatically verified,
raw and canonical row counts, drug/gene/source counts, invalid/duplicate handling, SHA-256, and the research-only
DGIdb disclaimer. Equivalent generated code must load only the pinned local file, validate the same metadata and
schema, and return `(resource_dataframe, summary_dict)`; it must not fetch or annotate AnnData.
