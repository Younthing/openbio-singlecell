# OpenBioSingleCellCellTypeCorrelation — official usage research

## Audited baseline

The current node explicitly calls Scanpy `tl.dendrogram` on a copied AnnData and then renders
`pl.correlation_matrix`, which is broadly consistent with official usage. However, it calls a matrix a “cell type
correlation” without disclosing that Scanpy first averages representation coordinates within each label and then
correlates those centroid vectors. It outputs only a PNG, allows implicit representation behavior through a blank
value, does not expose/record linkage settings or centroid/group support, and does not provide summary/code.

## Official Scanpy 1.12.3 interface

Primary sources:

- https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.dendrogram.html
- https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.correlation_matrix.html
- https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_dendrogram.py

Reviewed usage:

```python
info = sc.tl.dendrogram(
    work,
    groupby=population_key,
    use_rep=representation_key,
    cor_method=correlation_method,
    linkage_method=linkage_method,
    optimal_ordering=True,
    inplace=False,
)
matrix = info["correlation_matrix"]
order = info["categories_ordered"]
```

The result dict also contains linkage, group labels, representation, methods, and dendrogram coordinates. Official
documentation states that the correlation matrix is computed from average values of the selected representation or
gene profiles for predefined groups. `pl.correlation_matrix` is only a renderer of that result.

Pearson, Spearman, and Kendall are supported. A representation must have at least two informative dimensions and at
least two nonempty labels; non-finite/constant centroid vectors can yield undefined correlations and must fail rather
than appear as blank pixels. Using a named `obsm` representation avoids Scanpy's automatic PCA/X selection and
surprise computation.

## Scientific interpretation

This is descriptive group-centroid similarity and hierarchical ordering. Correlation across latent dimensions is
not invariant to arbitrary rotations/rescalings of an embedding and must not be interpreted as lineage, interaction,
causality, or formal population difference. Integrated representations can reduce Technical batch effects but can
also change biological geometry. Sample is not an inference unit here because there is no hypothesis test; cells
contribute to centroids and large Samples can dominate unless upstream sampling is balanced. Both facts belong in
limitations.

Population labels must disclose Provisional versus Curated status. The exact representation producer, dimensions,
cells per label, category order, correlation/linkage choices, and input observation fingerprint are reporting inputs.

## References to emit

- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15.
  https://doi.org/10.1186/s13059-017-1382-0
- Virtanen P, et al. SciPy 1.0. *Nature Methods*. 2020;17:261-272.
  https://doi.org/10.1038/s41592-019-0686-2
- Virshup I, et al. anndata: Annotated data. *JOSS*. 2024;9:4371.
  https://doi.org/10.21105/joss.04371

## Report/code implications

Report centroid construction, representation/dimensions/provenance, group cell counts/order, correlation matrix
extrema/ties, linkage, all limitations, references, and exact software versions. Generated code returns the
canonical pair table/matrix plus PNG-equivalent plotting inputs and `summary_dict`; it must not rely on a preexisting
dendrogram in caller AnnData.

