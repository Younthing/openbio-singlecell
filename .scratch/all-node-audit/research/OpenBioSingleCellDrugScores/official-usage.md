# OpenBioSingleCellDrugScores — official usage research

## Final open-expert boundary

The mean target score is numerically defined for any finite explicitly selected X/Raw/layer source. Count, scaled,
residual, transformed, and unknown-state evidence changes interpretation warnings only; no Raw/current binding or
history credential is required.

## Audited baseline and P0

The current node creates `pertpy.md.Drug()` and immediately reads `drug.dgidb.dictionary`. In Pertpy 1.3.0 that
attribute is created only by `DrugDataBase.set()`; a new object is unloaded, so execution raises `AttributeError`.
Even if initialized, the node would use an unversioned cached remote resource, score every DGIdb drug at once, leave
Pertpy's heterogeneous arrays in `uns`, and report neither target coverage nor clinical limitations.

## Official Pertpy 1.3.0 scoring interface

Primary sources:

- https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Enrichment.html
- https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/enrichment.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_enrichment.py
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/metadata/_drug.py

The public scoring call is:

```python
enrichment = pertpy.tl.Enrichment()
enrichment.score(
    work,
    layer=selected_layer,
    targets={drug_name: target_genes},
    nested=False,
    method="mean",
    key_added="openbio_drug",
)
```

Pertpy recommends log-normalized expression. With `method="mean"`, the score is the mean expression of matched
targets and is deterministic for any finite aligned numeric matrix; it does not itself enforce normalization.
Experts may explicitly select Raw/count-like input or retain all-zero observations, with the exact state and zero
support disclosed because score scale and interpretation differ from the recommendation. `method="seurat"` samples
controls with `numpy.random.default_rng()` without a seed in
Pertpy 1.3.0, so it is not reproducible through this interface and must not be offered. Pertpy writes score matrix,
variables, matched genes, and all genes to `uns` keys derived from `key_added`; the adapter must extract and validate
those values rather than expose backend storage as a public contract.

## Scientific-use findings

One DGIdb drug-target score is a descriptive expression summary, not predicted sensitivity, response, efficacy,
mechanism direction, dose, safety, or treatment recommendation. DGIdb interaction claims can include heterogeneous
evidence and do not specify that high target expression makes a drug beneficial. The exact resource release/license,
drug identifier, matched HGNC targets, expression source, and score definition must be visible.

Cells are observations, not biological replicates. A Condition comparison of scores must be a separate Sample-level
analysis within a defined population; Technical batch is nuisance downstream. This node does no statistical test.

## References to emit

- Cannon M, et al. DGIdb 5.0. *Nucleic Acids Research*. 2024;52:D1227-D1235.
  https://doi.org/10.1093/nar/gkad1040
- Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025.
  https://doi.org/10.1038/s41592-025-02909-7

## Report/code implications

Report one selected drug, exact requested/matched targets, caller-selected expression state and all-zero-observation
count, cells/features, score distribution,
resource release/hash/license/evidence sources, method/no-test status, warnings/limitations, references, and exact
versions. Generated code consumes the explicit resource DataFrame, performs no download, and returns
`(output_adata, summary_dict)` with identical results.
