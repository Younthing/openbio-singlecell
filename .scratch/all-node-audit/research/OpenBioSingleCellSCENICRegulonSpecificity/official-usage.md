# OpenBioSingleCellSCENICRegulonSpecificity — official usage research

## Audited baseline

pySCENIC 0.12.1 defines `pyscenic.rss.regulon_specificity_scores(auc_mtx, cell_type_series)`. Its return orientation
is **annotation groups by regulons**. The current node names the return index `regulon` and melts columns into
`group`, reversing both biological fields. Even in a compatible runtime this yields a semantically wrong table.

In the configured NumPy 2.4 environment, the official 0.12.1 implementation also fails because it allocates with the
removed alias `numpy.float`. The method is small and fully specified, so the supported node should implement and
test an exact local equivalent rather than import an obsolete runtime module.

## Official method and formula

Official source and protocol:

- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/rss.py
- https://pyscenic.readthedocs.io/en/stable/faq.html
- https://doi.org/10.1038/s41596-020-0336-2
- https://doi.org/10.1016/j.celrep.2018.10.045

The official call is:

```python
from pyscenic.rss import regulon_specificity_scores

rss = regulon_specificity_scores(
    auc_mtx=auc_matrix,                 # cells x regulons
    cell_type_series=annotation_labels # indexed by the same cells
)
# rss.index is annotation group; rss.columns are regulons
```

For one regulon and one group, pySCENIC normalizes the nonnegative AUC vector across all cells and compares it with
the normalized binary group-membership vector. It returns:

```text
RSS = 1 - Jensen-Shannon distance(normalized regulon AUC, normalized group indicator)
```

SciPy's `jensenshannon` returns the square root of Jensen-Shannon divergence. RSS is a descriptive similarity score,
not a hypothesis test, effect size between Conditions, probability, or FDR. No p-value is produced. The function's
cell-level distribution makes large groups influential and does not account for Sample replication or Technical
batch.

Zero-sum regulon AUC makes normalization undefined; non-finite/negative values, missing labels, and index
misalignment likewise require a hard preflight. A one-observed-group input remains numerically defined but has no
between-group specificity interpretation, so an open expert tool should return it with a prominent warning rather
than reject it. Categorical levels with zero observed cells should be excluded and disclosed. Results must retain
the complete group-by-regulon grid before any top-hit summary.

## Scientific interpretation

RSS can prioritize regulons whose AUCell activity pattern is concentrated in a curated or provisional annotation.
It does not prove TF binding or direct regulation, and a high score does not imply differential activity between
Conditions. Cells from one Sample remain correlated observations; Condition inference requires a separate
Sample-level design. Annotation status must therefore be visible in the report.

## References to emit

- Suo S, et al. Revealing the critical regulators of cell identity in the mouse cell atlas. *Cell Reports*.
  2018;25:1436-1445.e3. https://doi.org/10.1016/j.celrep.2018.10.045
- Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature
  Protocols*. 2020;15:2247-2276. https://doi.org/10.1038/s41596-020-0336-2
- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463

## Report and generated-code implications

The report must state the exact formula and non-inferential status; disclose annotation status/order and group
sizes, total cells/regulons, full grid size, zero-sum checks, bounded top regulons per group, SCENIC artifact/resource
provenance, a one-group interpretability warning when applicable, limitations, references, and dynamic versions.
Equivalent code must compute the audited NumPy/SciPy
formula without importing `pyscenic.rss`, preserve orientation, validate the complete grid, and return
`(table, summary_dict)`.
