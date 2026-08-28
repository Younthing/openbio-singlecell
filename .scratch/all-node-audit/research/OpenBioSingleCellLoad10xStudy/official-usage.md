# Load 10x Study: official usage research

Researched: 2026-08-28

## Official Scanpy/AnnData pattern

Scanpy's current clustering tutorial reads each 10x H5 Sample separately, stores the objects in a mapping keyed by Sample ID, and concatenates them with an explicit Sample label:

```python
import anndata as ad
import scanpy as sc

samples = {
    "donor_1": sc.read_10x_h5("donor_1/filtered_feature_bc_matrix.h5"),
    "donor_2": sc.read_10x_h5("donor_2/filtered_feature_bc_matrix.h5"),
}
adata = ad.concat(samples, label="sample", index_unique="-")
```

- Scanpy tutorial: https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering.html
- Official concatenation API and examples: https://anndata.readthedocs.io/en/stable/generated/anndata.concat.html

`anndata.concat` uses mapping keys as dataset keys. `label` writes those keys to an observation column, and `index_unique` appends the key to original indices. `join="inner"` takes the feature intersection; `join="outer"` takes the union. The official documentation warns that outer joins pad missing sparse variables with zeros, which can be confused with measured biological zero unless disclosed.

## Reporting and domain references

- Füllgrabe A et al. Guidelines for reporting single-cell RNA-seq experiments. *Nature Biotechnology*. 2020;38:1384-1386. https://doi.org/10.1038/s41587-020-00744-z
- Virshup I et al. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

Repository domain rule: a **Sample** is an independent biological specimen and the replicate unit for quality assessment and Condition-level inference. A **Technical batch** is a separate source of unwanted variation. Directory names used by this loader therefore assert Sample identity; they must not be described as mere cell batches.

## Scientific contract

- The selected root has a strict inventory: every first-level directory is declared to be one Sample and must contain one valid 10x MTX triplet. Never silently omit an invalid child, because that changes the Study and may remove a biological replicate.
- Use each directory name as the Sample label; require non-empty, distinct labels and preserve the complete Sample inventory in provenance.
- Apply the same count, barcode, feature-ID, variable-name, and `gex_only` contract as `Load 10x MTX` to every Sample.
- `join="inner"` retains only features shared by all Samples and must fail if the intersection is empty. `join="outer"` retains the union and must warn that absent sparse features are zero-filled by AnnData.
- Preserve each cell's original barcode while creating a globally unique observation index with the Sample key. The Sample column remains the explicit biological identity.
- Report per-Sample input/retained dimensions, common/union feature counts, naming policy, join semantics, and all source fingerprints in OpenBio provenance.

## Open expert-boundary re-review (2026-08-28)

`anndata.concat(join="inner")` is defined even when the feature intersection is empty. An empty intersection is thus
a reportable Study outcome rather than a concatenation error. Outer-join zero filling remains disclosed. Invalid
Sample directories still fail because there is no unambiguous payload to concatenate and silently omitting a declared
Sample would change the Study. Per-Sample axis-identity concerns inherited from the single-Sample loader are warnings;
an actual duplicate feature index that AnnData cannot align remains a computation precondition.
