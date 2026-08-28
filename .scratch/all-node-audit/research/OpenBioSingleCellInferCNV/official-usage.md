# Infer CNV: official usage and scientific practice

## Versioned official API

The reviewed optional extra now locks and real-smokes infercnvpy 0.6.1. Its public interface is:

```python
infercnvpy.tl.infercnv(
    adata, *, reference_key=None, reference_cat=None, reference=None,
    lfc_clip=3, window_size=100, step=10, dynamic_threshold=1.5,
    exclude_chromosomes=("chrX", "chrY"), chunksize=5000,
    n_jobs=None, inplace=True, layer=None, key_added="cnv",
    calculate_gene_values=False,
)
```

It stores the smoothed CNV representation in `obsm["X_<key_added>"]` and genomic-window metadata under `uns[key_added]`.

The 0.6.1 implementation indexes its convolution output into `var_names` with a two-dimensional NumPy index. Under
Pandas 3, an Arrow-backed string index rejects that operation. The reviewed private backend copy therefore uses an
equivalent object-string Index, while the public input/output axes and fingerprints remain exact and unchanged.

Primary sources:

- [infercnvpy `tl.infercnv`](https://infercnvpy.readthedocs.io/en/latest/generated/infercnvpy.tl.infercnv.html)
- [The inferCNV method and input preparation](https://infercnvpy.readthedocs.io/en/latest/infercnv.html)
- [Official example with normal reference cells](https://infercnvpy.readthedocs.io/en/latest/notebooks/reproduce_infercnv.html)
- [GTF coordinate annotation](https://infercnvpy.readthedocs.io/en/latest/generated/infercnvpy.io.genomic_position_from_gtf.html)
- [infercnvpy 0.6.1 release](https://github.com/icbi-lab/infercnvpy/releases/tag/v0.6.1)
- Tirosh et al., *Science* 2016, DOI [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
- Patel et al., *Science* 2014, DOI [10.1126/science.1254257](https://doi.org/10.1126/science.1254257)

## Expression, coordinate and reference contracts

Official infercnvpy recommends a post-QC, **normalized and log-transformed**, preferably full-gene expression matrix. It subtracts reference expression in the supplied value space, smooths genes in genomic order, centers each cell, thresholds noise and median-filters the result. OpenBio therefore defaults to the named `log1p_norm` layer and reports any expression-state evidence it can discover, but it does not use mutable OpenBio history as an execution credential. An explicit expert selection is sufficient: unknown, count-like, scaled, or otherwise nonstandard finite input remains executable and is disclosed with a warning. The implementation cannot prove that a current feature axis is complete merely from AnnData contents.

`adata.raw` is not a hidden prerequisite or an axis-binding authority for this node. The selected current-axis `X` or layer, its current feature identifiers, and its genomic-coordinate table define the calculation. No history-deletion or Raw/current equality test is part of the runtime contract. Raw may still be selected in other nodes that expose it, but this node intentionally exposes current-axis `X`/layer inputs because the coordinate annotations consumed by infercnvpy live on the current `var` axis.

Genomic coordinates must match the genome annotation used for quantification. The canonical columns are `chromosome`, `start`, and `end`; merely checking column names is insufficient. Require unique feature IDs, complete chromosome labels, finite integer coordinates, `start < end`, supported chromosome naming, and enough ordered genes per retained chromosome/window. infercnvpy 0.6.1's implementation only processes labels beginning with `chr` and always omits `chrM`, independently of the public `exclude_chromosomes` default; OpenBio must make that effective exclusion explicit. Gene-ID/version mapping and the user-declared GTF assembly belong in provenance.

For best results official examples use known non-malignant cells as the background. Although upstream infercnvpy permits all-cell average reference, doing that silently changes the estimand and can attenuate widespread tumor signal. The OpenBio node currently has an empty default `reference_categories`, so its default path fails before calling infercnvpy. A safe scientific interface requires explicit nonempty reference categories and verifies exact typed membership, number of cells and `Sample` distribution. It must not string-coerce categories or silently fall back to all cells.

## Limits and disclosure

infercnvpy's own project warns that the package is experimental and its results have not been formally validated beyond qualitative similarity to inferCNV. Expression-derived CNV is indirect evidence, sensitive to cell type, technical effects, reference choice, gene coverage and smoothing parameters; it does not replace DNA CNV measurement.

Report the exact selected source, any non-authoritative state evidence, a clear `full_gene_completeness_verified=false` disclosure, reference cells/categories and declared biological-Sample coverage, coordinate/assembly provenance, genes/chromosomes/windows used/excluded, algorithm parameters, CNV amplitude summaries, output fingerprints, warnings, references and versions. No `Condition` inference is performed.
