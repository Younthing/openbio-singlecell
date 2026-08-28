# CNV Score: official usage and scientific practice

## Data-flow clarification

The official function reads both `adata.obsm[f"X_{use_rep}"]` and `adata.obs[groupby]`. In the atomic OpenBio graph,
the immutable `OPENBIO_CNV_STATE` proves the inferred matrix while a separate downstream AnnData carries annotations
created by PCA/Neighbors/Leiden or another audited partitioning path. The reviewed node therefore accepts both,
verifies exact axes and CNV-state fingerprints, and only then calls the official API on the downstream copy.

## Audited API and estimand

The replacement is audited against infercnvpy 0.6.1:

```python
infercnvpy.tl.cnv_score(
    adata,
    groupby="cnv_leiden",
    *,
    use_rep="cnv",
    key_added="cnv_score",
    inplace=True,
    obs_key=None,
)
```

For every observed group, the implementation calculates one scalar: the mean absolute value over **all cells and
all inferred-CNV windows in that group**. It then assigns the same group scalar to each member cell. The deprecated
`obs_key` alias must remain `None`.

Primary sources:

- [infercnvpy `tl.cnv_score`](https://infercnvpy.readthedocs.io/en/stable/generated/infercnvpy.tl.cnv_score.html)
- [infercnvpy 0.6.1 source](https://github.com/icbi-lab/infercnvpy/blob/v0.6.1/src/infercnvpy/tl/_scores.py)
- [infercnvpy analysis tutorial](https://infercnvpy.readthedocs.io/en/latest/notebooks/reproduce_infercnv.html)
- Tirosh et al., *Science* 2016, DOI [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)

## Interpretation and validation

The score is group-dependent descriptive evidence. It is not absolute copy number, a per-cell independent
measurement, a tumor probability, a calibrated threshold, or a hypothesis test. Changing the grouping changes the
score even when `X_cnv` is unchanged. Tumor/normal interpretation requires domain knowledge and preferably DNA
validation; the node must never create a universal label or cutoff.

Require a complete categorical grouping aligned to a validated `OPENBIO_CNV_STATE`; reject missing labels, implicit
string conversion, unused/category-display collisions, one empty family, or a bare CNV matrix. Verify finite real
CNV values and output collisions. Independently recompute every group score and per-cell assignment from the exact
matrix, require agreement with the public backend, and report group size, score, absolute-value distribution,
reference/Sample support, grouping status and limitations. No `Condition` comparison is performed.
