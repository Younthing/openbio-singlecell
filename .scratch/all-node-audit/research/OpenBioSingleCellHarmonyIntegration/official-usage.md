# Harmony Integration: official usage research

Researched: 2026-08-28

## Scope and primary sources

This record covers `OpenBioSingleCellHarmonyIntegration`: correction of a previously calculated PCA embedding for user-declared **Technical batch** covariates. It does not cover PCA, neighbor-graph construction, clustering, biological annotation, or formal **Condition** inference.

Primary sources:

- Scanpy `scanpy.external.pp.harmony_integrate` reference and example: https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html
- harmonypy method-author repository and current `run_harmony` reference implementation: https://github.com/slowkow/harmonypy and https://github.com/slowkow/harmonypy/blob/master/harmonypy/harmony.py
- Harmony method site and multi-covariate tutorial: https://portals.broadinstitute.org/harmony/index.html and https://portals.broadinstitute.org/harmony/articles/quickstart.html
- Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*. 2019;16:1289-1296. https://doi.org/10.1038/s41592-019-0619-0
- Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods*. 2022;19:41-50. https://doi.org/10.1038/s41592-021-01336-8

The repository pins Scanpy `>=1.12.3,<1.13`. This audit targets the method-author **harmonypy 2.0.0** release
(`v2.0.0`, released 2026-04-26), whose C++ rewrite and result orientation are breaking changes from the 0.x
series. The release extra therefore remains pinned to the reviewed 2.0.0 baseline. Runtime may accept a 2.0.x patch
only when module/distribution identities agree, the required public interface is present, and the returned
shape/orientation is validated rather than guessed; rejecting a compatible patch by number alone is not a
computational safeguard. The 0.x family and unversioned future major/minor interfaces remain unsupported.

## Canonical workflow

Scanpy documents Harmony as an embedding correction after PCA and before graph construction:

```python
import scanpy as sc
import scanpy.external as sce

sc.pp.pca(adata)
sce.pp.harmony_integrate(
    adata,
    key="technical_batch",
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
)
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
```

`key` may be a string or a sequence of observation-column names. Harmony's own tutorial likewise supports several covariates through `vars_use=[...]`. The corrected result is another `n_obs x n_components` embedding; Harmony does not produce corrected gene expression and the original PCA coordinates should remain available.

The declared covariates are effects that the analyst intends to remove. In this repository they must be described as **Technical batch** covariates. A **Sample** is the biological replicate and a **Condition** is the biological cohort/treatment of interest; neither may be silently renamed to “batch.” Correcting a Sample or Condition can remove biological signal, especially when it is confounded with laboratory, lane, chemistry, or acquisition effects.

## Audited harmonypy interface

The current method-author implementation exposes:

```python
harmonypy.run_harmony(
    data_mat,
    meta_data,
    vars_use,
    theta=None,
    lamb=None,
    sigma=0.1,
    nclust=None,
    tau=0,
    block_size=0.05,
    max_iter_harmony=10,
    max_iter_kmeans=4,
    epsilon_cluster=1e-3,
    epsilon_harmony=1e-2,
    alpha=0.2,
    batch_prop_cutoff=1e-5,
    verbose=True,
    random_state=0,
    ncores=0,
)
```

The returned Harmony object exposes corrected coordinates (`Z_corr`), original coordinates, soft cluster assignments, centroids, Harmony objective values, k-means objective values, and k-means rounds. Version 2.0.0 changed several historical defaults to match the newer R implementation; notably `max_iter_kmeans` is 4, not the legacy node's 100. All scientific and convergence defaults used by the node must be passed explicitly and recorded so an upstream release cannot silently change the result.

Important parameter meanings from the official implementation and Harmony reference are:

- `theta`: diversity penalty; larger values encourage greater Technical-batch diversity within soft clusters;
- `lamb`: ridge penalty; current harmonypy uses automatic estimation for `None`, while a fixed value must be positive;
- `sigma`: soft-clustering bandwidth;
- `nclust`: number of soft clusters, with an automatic size-dependent default;
- `tau`: protection against over-correction of small batches;
- `max_iter_harmony` and `max_iter_kmeans`: iteration ceilings, not evidence that convergence occurred;
- `random_state`: k-means initialization seed;
- `ncores`: BLAS thread count; `1` is the explicit single-thread setting.

The harmonypy 2.0.0 README, quick start, source contract, and release test all define `Z_corr` as **cells × PCs**
with the same shape and orientation as a cells × PCs input. The official 2.0.0 Scanpy example assigns
`harmony_out.Z_corr` directly to `adata.obsm[adjusted_basis]`. The historical Scanpy wrapper currently stores
`harmony_out.Z_corr.T`; that adapter was written for harmonypy 0.x, whose result used PCs × cells, and is not the
2.0.0 contract. Therefore the reviewed adapter must require the audited 2.0.x interface family, require
`result.Z_corr.shape == input_basis.shape == (n_obs, n_components)`, and store it without transposition. It must
not support both orientations by guessing from shape: a square embedding makes such guessing undecidable and
would silently accept an unsupported package contract.

Primary version/orientation evidence:

- harmonypy v2.0.0 release: https://github.com/slowkow/harmonypy/releases/tag/v2.0.0
- method-author quick start and Scanpy example: https://github.com/slowkow/harmonypy/tree/v2.0.0
- v2.0.0 implementation and release test: https://github.com/slowkow/harmonypy/blob/v2.0.0/harmonypy/harmony.py and https://github.com/slowkow/harmonypy/blob/v2.0.0/tests/test_harmony.py

## Input contract

The selected `obsm` basis must be a numeric, finite, two-dimensional matrix with exactly `n_obs` rows and at least one retained component. One-dimensional correction and zero-variance components are weak or degenerate scientific choices, but they are not wrapper-level type/axis failures: warn and let the reviewed backend either produce a finite exact-shape result or fail. Observation names must be unique. Because Harmony consumes `.obsm` and `.obs`, an empty or duplicate current **variable** axis is irrelevant and must not veto an otherwise aligned embedding. Every declared Technical-batch column must:

- exist in `adata.obs`;
- contain no null or blank value;
- disclose when it contains only one level, because that covariate contributes no correction information;
- disclose singleton and other small levels as unstable rather than imposing an undocumented sample-size gate;
- be distinct from any explicitly protected Sample or Condition column.

The node cannot establish from a column name that an effect is genuinely technical. The report must state that the keys were user-declared and warn that a confounded design cannot distinguish unwanted and biological variation. It must never claim that Sample or Condition effects were safely preserved merely because Harmony completed.

Reject an empty observation axis, backed inputs, malformed or non-finite embeddings, duplicate Technical-batch keys, and unsupported runtime signatures. The destination must be nonblank and distinct from the input basis: `overwrite_existing` authorizes replacement of a prior destination, not destruction of the source embedding. Preserve the original PCA basis and every gene-expression representation.

## Audited cluster-count boundary

harmonypy 2.0.0's Python layer resolves automatic `nclust` with `int(min(round(N / 30), 100))`; for very small `N` that expression can yield zero even though the C++ implementation allocates `K` clusters. Its without-replacement centroid initialization also cannot terminate when `K > N`. The adapter therefore resolves automatic mode to `max(1, min(round(N / 30), 100))` and hard-requires `1 <= nclust <= n_obs`. These are backend computability/safety bounds, not a recommendation that two or more clusters are scientifically superior.

The 2.0.0 Python layer leaves a scalar `sigma` as a Python float when `nclust=1` and later calls `.astype` on it. To support the otherwise valid one-cluster expert choice, the adapter supplies an explicit length-`K` NumPy sigma vector for every resolved `K`. Both runtime and generated code use that public argument form. Large but executable `nclust`, iteration counts, and parameter values are disclosed rather than blocked by arbitrary UI ceilings.

## Output and diagnostic semantics

The scientific output is an AnnData copy with the corrected embedding stored at the declared `adjusted_basis`. Basic execution diagnostics may report:

- input, raw `Z_corr`, and stored output embedding shapes;
- ordered Technical-batch keys, level counts, and cell counts per level;
- resolved model parameters and seed;
- number of Harmony rounds and k-means rounds exposed by the backend;
- initial/final objective values when available and finite;
- per-cell displacement norms between input and corrected coordinates;
- finite-value and shape checks on the corrected embedding.

These diagnostics show what the algorithm did, not whether integration is scientifically successful. The output must not say “batch effects removed,” “conditions preserved,” “clusters improved,” or “integration succeeded” based only on completion, objective decrease, mixing, or displacement. The Luecken et al. benchmark evaluates integration as a trade-off between Technical-batch removal and biological-variation conservation using multiple metrics. Those comparisons belong to separate downstream diagnostic nodes and require analyst-declared biological labels or trajectories.

## References and software reporting

The strict JSON summary should cite the Harmony paper, Scanpy wrapper, harmonypy implementation, and the integration-benchmark paper. Runtime versions should include openbio-singlecell, Scanpy, harmonypy, AnnData, NumPy, pandas, and the numerical backend where discoverable. CPU execution and the resolved BLAS thread policy must be disclosed.

`original_basis_preserved` is a verified structural fact only when input and destination keys are distinct and the source coordinates remain value-identical on the returned copy. It must never be hard-coded true after an in-place key replacement.
