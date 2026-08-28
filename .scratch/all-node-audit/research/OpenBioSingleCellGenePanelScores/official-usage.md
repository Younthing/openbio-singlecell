# OpenBioSingleCellGenePanelScores — official usage research

## Final open-expert boundary

Any finite explicitly selected X/Raw/layer source remains executable. Count/scaled/residual/unknown evidence is
reported as a limitation, and Raw never has to match current X or a mutable history record.

## Audited baseline

The current node parses a multi-panel file, loops over every panel, and calls Scanpy `tl.score_genes`, writing an
unbounded set of dynamically named `obs` columns. It does not carry resource metadata/hash, permits Raw snapshot
selection even when it contains counts, silently skips panels with no overlap, and fixes `ctrl_size=50` irrespective
of panel size. Its `use_raw=False` override is correct, but the public operation is too broad and the scientific
inputs are under-specified.

The reviewed runtime is Scanpy 1.12.3 with AnnData 0.13.2.

## Official Scanpy 1.12.3 usage

Primary sources:

- https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html
- https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_score_genes.py
- https://scanpy.readthedocs.io/en/stable/how-to/cell-cycle.html

The installed 1.12.3 signature is:

```text
score_genes(adata, gene_list, *, ctrl_as_ref=True, ctrl_size=50,
            gene_pool=None, n_bins=25, score_name="score",
            random_state=0, copy=False, use_raw=None, layer=None)
```

The score is mean panel expression minus mean expression of a matched reference set sampled from expression bins.
Official documentation notes that `ctrl_as_ref` will change to `False` in Scanpy 2.0 and suggests setting
`ctrl_size=len(gene_list)` when the list is not very small. The adapter must set `use_raw=False` and the selected
layer explicitly rather than inherit Scanpy's `use_raw=None` behavior. For stable semantics across Scanpy versions,
`ctrl_as_ref=False` should be fixed and reported.

```python
sc.tl.score_genes(
    work,
    gene_list=matched_genes,
    gene_pool=work.var_names,
    ctrl_as_ref=False,
    ctrl_size=resolved_ctrl_size,
    n_bins=n_bins,
    score_name=output_key,
    random_state=random_seed,
    use_raw=False,
)
```

Scanpy's how-to applies library-size normalization and log1p before scoring. This repository's Raw snapshot normally
contains counts, so normalized expression remains the recommended default. The public `score_genes` calculation is
still arithmetically defined for an explicitly selected Raw/count-like matrix; expert use is therefore allowed with a
prominent warning rather than blocked solely because it departs from the recommended practice.

## Scientific-use findings

One score is descriptive for one predeclared panel in each cell. It is not a formal gene-set enrichment p-value,
does not prove pathway activation, and is sensitive to the gene pool and expression cohort used to construct
controls. Cells are not independent biological replicates. Any Condition comparison must aggregate/model the score
at the Sample level within a population and account for Technical batch as nuisance where appropriate.

The panel resource needs exact organism, identifier namespace, version/date, scope, license, citation, and file
SHA-256. Target genes are matched exactly; no case conversion, aliases, or online translation. Missing panel genes,
matched fraction, resolved controls, bins, and all warnings belong in the report.

## References to emit

- Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.
  *Science*. 2016;352:189-196. https://doi.org/10.1126/science.aad0501
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I, et al. anndata: Annotated data. *JOSS*. 2024;9:4371. https://doi.org/10.21105/joss.04371

## Report/code implications

Report the exact matched/requested genes (bounded preview plus counts), panel/resource identity and SHA, expression
state, gene pool, resolved control size, bins, seed, score summary, excluded targets, lack of inference, references,
and dynamic versions. Equivalent generated code must validate one named panel and return
`(output_adata, summary_dict)` with exact report equality and no download.
