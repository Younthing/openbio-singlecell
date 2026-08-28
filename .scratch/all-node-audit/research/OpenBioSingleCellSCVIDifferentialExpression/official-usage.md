# OpenBioSingleCellSCVIDifferentialExpression: official usage research

## Scope and verification snapshot

This pre-change record covers one operation: reuse one already-fitted `scvi.model.SCVI` model to compare exactly two cell populations with `SCVI.differential_expression`. It does not cover model fitting, integration-quality assessment, multi-group marker ranking, or formal Sample-level Condition inference.

The public API and implementation were checked on 2026-08-28 against:

- the current scvi-tools documentation;
- the current released scvi-tools `1.5.0.post1` source at signed release commit `56520c7`;
- the scVI and effect-size-aware differential-expression method papers; and
- primary studies of pseudoreplication in single-cell condition comparisons.

The repository does not pin or directly declare scvi-tools. Any implementation must therefore verify the installed public signature at runtime, fail on an incompatible API, and report the installed version rather than assuming `1.5.0.post1`.

## Primary and official sources

### Software documentation and source

- scvi-tools `SCVI` API, including `setup_anndata`, `train`, and `differential_expression`: https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html
- scvi-tools introductory workflow and current differential-expression examples: https://docs.scvi-tools.org/en/latest/tutorials/notebooks/quick_start/api_overview.html
- scvi-tools differential-expression background: https://docs.scvi-tools.org/en/latest/user_guide/background/differential_expression.html
- scvi-tools preprocessing guidance: https://docs.scvi-tools.org/en/latest/tutorials/notebooks/use_cases/preprocessing.html
- scvi model guide, including the interpretation of decoded expression: https://docs.scvi-tools.org/en/latest/user_guide/models/scvi.html
- signed scvi-tools releases (`1.5.0.post1`, commit `56520c7`): https://github.com/scverse/scvi-tools/releases
- released `RNASeqMixin.differential_expression` implementation: https://github.com/scverse/scvi-tools/blob/56520c7/src/scvi/model/base/_rnamixin.py#L330-L454
- released DE selection, sorting, and FDR-tag implementation: https://github.com/scverse/scvi-tools/blob/56520c7/src/scvi/model/base/_de_core.py#L97-L190
- released posterior comparison and batch-conditioning implementation: https://github.com/scverse/scvi-tools/blob/56520c7/src/scvi/model/base/_differential.py#L63-L456
- released raw-count summary implementation: https://github.com/scverse/scvi-tools/blob/56520c7/src/scvi/model/_utils.py#L159-L218
- released model-instance/class manager registration and targeted deregistration implementation: https://github.com/scverse/scvi-tools/blob/56520c7/src/scvi/model/base/_base_model.py#L299-L378

No community example was needed.

### Method and reporting papers

- Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. *Deep generative modeling for single-cell transcriptomics*. Nature Methods 15, 1053-1058 (2018). https://doi.org/10.1038/s41592-018-0229-2
- Boyeau P, Lopez R, Regier J, Gayoso A, Jordan MI, Yosef N. *Deep Generative Models for Detecting Differential Expression in Single Cells*. MLCB 2019 proceedings. https://mlcb.github.io/mlcb2019_proceedings/papers/paper_32.pdf ; preprint DOI: https://doi.org/10.1101/794289
- Gayoso A, Lopez R, Xing G, et al. *A Python library for probabilistic analysis of single-cell omics data*. Nature Biotechnology 40, 163-166 (2022). https://doi.org/10.1038/s41587-021-01206-w
- Squair JW, Gautier M, Kathe C, et al. *Confronting false discoveries in single-cell differential expression*. Nature Communications 12, 5692 (2021). https://doi.org/10.1038/s41467-021-25960-2
- Zimmerman KD, Espeland MA, Langefeld CD. *A practical solution to pseudoreplication bias in single-cell studies*. Nature Communications 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1

## Required upstream model state

The official workflow is `SCVI.setup_anndata(...)`, construct `SCVI(...)`, `model.train(...)`, then call a downstream method. `setup_anndata` registers where the model's data live; its `layer` argument denotes raw count data and its `batch_key` denotes the primary registered batch. The current tutorial states that scvi-tools models normally require raw, non-log-library-size-normalized counts and recommends roughly 1,000-10,000 HVGs for scVI depending on context.

The DE call is conditional on the fitted model and its registered data. Its `adata` must have structure equivalent to the model's initial AnnData. The model consumes the registered count representation and learned parameters; it does not reinterpret an arbitrary downstream `adata.X` as a new expression source.

For this repository, the relevant producer is `OpenBioSingleCellSCVIIntegration`. Its current `SCVIModel` artifact:

- contains a live, trained, process-local model rather than serialized fitted weights;
- owns a private registered AnnData copy;
- binds the model to exact ordered `var_names` and to the set of training `obs_names`;
- allows a downstream observation subset but rejects unknown observations or a changed variable order;
- copies requested downstream observation metadata onto the registered count view;
- serializes temporary AnnData-manager use and deregisters temporary managers;
- snapshots training parameters and training diagnostics; and
- does not retrain during DE.

Implementation verification on the installed `1.5.0.post1` release exposed an important public-manager lifecycle
detail. A transferred analysis AnnData is initially registered only in the fitted model's per-instance manager store,
whereas `deregister_manager(adata)` resolves its target through the class store. Therefore safe targeted cleanup is:
retrieve the exact transferred manager with `get_anndata_manager(analysis_adata, required=False)`, publish that same
manager with `register_manager`, and then call `deregister_manager(analysis_adata)`. This removes only the temporary
analysis manager. Calling `deregister_manager()` without an AnnData is not acceptable because it also clears unrelated
class and instance managers. A real one-epoch CPU smoke test confirmed this exact sequence on `1.5.0.post1`.

The current integration artifact was designed for a representation-learning view, normally the repository's 5,000-HVG view. Consequently, a downstream scVI DE table is limited to the genes fitted by that model. It is not a full-gene analysis merely because the downstream AnnData still has other representations elsewhere.

Necessary model evidence for every DE report is therefore:

- model class and exact fitted-model fingerprint when the artifact supports one;
- process-local/serialization status;
- ordered training observation and feature identities or their fingerprints;
- registered raw-count source, declared expression state, and feature count;
- primary Technical batch key plus other categorical/continuous nuisance covariates;
- constructor parameters, train/validation split, requested and actual epochs, early-stopping policy, batch size, random seed, accelerator/device, and finite exposed training metrics;
- installed scvi-tools, PyTorch, Lightning, AnnData, NumPy, pandas, and openbio-singlecell versions; and
- an explicit statement that training completion and finite loss do not prove adequate fit, batch removal, or preservation of biology.

## Current public API

The released `1.5.0.post1` signature is:

```python
SCVI.differential_expression(
    adata=None,
    groupby=None,
    group1=None,
    group2=None,
    idx1=None,
    idx2=None,
    mode="vanilla",
    delta=0.25,
    batch_size=None,
    all_stats=True,
    batch_correction=False,
    batchid1=None,
    batchid2=None,
    fdr_target=0.05,
    silent=False,
    weights="uniform",
    filter_outlier_cells=False,
    importance_weighting_kwargs=None,
    dataloader=None,
    **kwargs,
)
```

`**kwargs` are passed to the public `DifferentialComputation.get_bayes_factors` method. Defaults relevant to this node in release `1.5.0.post1` are `n_samples_overall=5000`, `use_permutation=False`, `m_permutation=10000`, `pseudocounts=None`, `threshold_counts=0.01`, `test_mode="three"`, and no credible intervals.

The repository adapter must check the installed signature for every keyword it uses. It must not silently catch `TypeError` and retry with a scientifically different subset of parameters.

## Population-selection semantics

- `groupby` names an `adata.obs` column.
- The annotated type of `group1` is a list of group labels. Official tutorials also pass a scalar string, and the released source explicitly wraps a scalar into a one-element list.
- `group2` is a scalar label. If it is `None`, each `group1` is compared with the union of all other groups.
- `idx1`/`idx2` are an alternative interface accepting masks, indices, or `pandas.DataFrame.query` strings. `idx1` overrides `group1`/`group2`.
- A list of multiple `group1` labels produces multiple result blocks. This is a convenience batch operation, not one atomic contrast.
- With explicit `group1` and `group2`, observations in other labels and observations with missing group labels are not members of either comparison population.

For an atomic workflow node, one scalar `group1` and one distinct scalar `group2` should be required. The broad `idx1`/`idx2` query interface, implicit one-versus-rest behavior, and multiple-`group1` loop should remain hidden.

## What the method compares

The scvi-tools background describes the population-specific normalized expression distribution by aggregating cell-level approximate posteriors. For each population, the implementation repeatedly samples cells/latent states and decodes normalized expression through the trained nonlinear model. The comparison is therefore conditional on:

- the selected cells;
- the fitted weights and approximate posterior;
- the model's registered count view and fitted features;
- the batch/covariate conditioning policy; and
- Monte Carlo sampling.

It is not a test on the downstream AnnData's normalized `X`. It is also not a test whose independent observations are Samples.

## `mode`, `delta`, and direction

The two modes are scientifically different.

| Mode | Released computation | Interpretation and reporting |
|---|---|---|
| `vanilla` | Per gene, estimate `P(scale1 > scale2)` and its complement, then return their log posterior odds as `bayes_factor`. | Directional evidence for group 1 exceeding group 2. The public background calls this the point-null formulation; the released source is more precise about the computed directional comparison. It has no practical-effect-size threshold and should not be the new default. |
| `change` | Form posterior log2 fold-change samples and compare them with `delta`. | Effect-size-aware model evidence. The official background recommends this mode to avoid prioritizing tiny effects. Positive `lfc_mean` means decoded expression is higher in `group1` than `group2`. |

For `change`, if `pseudocounts=None`, the released code estimates a positive offset from genes with low observed mean counts. The posterior hypothesis uses that offset. The reported LFC distribution summaries are computed with `1e-3 * pseudocounts`, a released implementation detail that must be recorded rather than reinterpreted.

The lower-level `test_mode` also changes `proba_de`:

- `test_mode="two"` adds posterior probabilities for `LFC >= delta` and `LFC < -delta`, matching the composite event `abs(LFC) >= delta` up to the boundary convention;
- the released default `test_mode="three"` takes the maximum of the two directional probabilities, treating the two directions as separate alternatives.

An adapter must pass and disclose this choice explicitly. Merely reporting `mode="change"` and `delta` is insufficient.

## Posterior FDR tag

For `change` mode, scvi-tools sorts genes by `proba_de` and adds a dynamic column named `is_de_fdr_<fdr_target>`. It selects the largest prefix whose mean posterior non-DE probability is at most `fdr_target`.

This is posterior expected FDR under the fitted model and selected cells. It is not a p-value, not Benjamini-Hochberg correction, and not evidence that Type-I error is controlled across independent biological Samples. Reports should say “posterior-FDR-tagged” rather than “statistically significant.” `vanilla` mode does not receive this Boolean FDR column from `_de_core`.

## Batch-conditioning semantics

`batch_correction` controls decoder conditioning; it does not retroactively assess whether integration succeeded.

| Public arguments | Released behavior |
|---|---|
| `batch_correction=False` | `_de_core` sets `use_observed_batches=True`. Each selected cell is decoded conditional on its observed primary registered batch. Different batch composition between groups therefore remains part of the population comparison. `batchid1`/`batchid2` must not be supplied. |
| `batch_correction=True`, both batch-ID lists omitted | Both groups are decoded over all categories of the primary registered `batch_key`; per-batch decoded values are averaged and same-batch comparisons are formed. This can include group-by-batch combinations with no observed support. |
| `batch_correction=True`, identical `batchid1` and `batchid2` sets | Both groups are decoded on the same requested batches and posterior pairs are formed within batch. This is the clearest standardized comparison. |
| `batch_correction=True`, disjoint ID sets | Each group is decoded on its own disjoint set. The implementation permits it, but the result mixes biological comparison with different decoder batch targets. |
| Partially overlapping, non-identical ID sets | The released source emits a warning that correction is not trustworthy and then falls back to an unmatched comparison. A strict workflow adapter must reject this case. |

The API documentation requires `batchid1` and `batchid2` to be either exactly equal or disjoint. Batch identifiers refer only to the primary `batch_key` registered by `setup_anndata`. Additional categorical and continuous nuisance covariates remain conditioned according to the model/data path; `batchid1`/`batchid2` do not standardize their population composition. A report must therefore list those covariates and must not claim that `batch_correction=True` removed every confounder.

Counterfactual decoder conditioning also is not a substitute for overlap. If group and Technical batch are completely confounded, decoding unobserved combinations is extrapolation. A conservative adapter should use only Technical-batch categories observed in both populations, pass the same ordered set for both groups, and fail when there is no shared support.

## Released output fields

The returned object is a pandas DataFrame indexed by the model's feature names and sorted descending by `bayes_factor` for `vanilla` or `proba_de` for `change`.

Core `vanilla` fields are:

- `proba_m1`, `proba_m2`, `bayes_factor`;
- `scale1`, `scale2`; and
- when `all_stats=True`, the raw-count fields below.

Core `change` fields are:

- `proba_de`, `proba_not_de`, `bayes_factor`;
- `scale1`, `scale2`, `pseudocounts`, `delta`;
- `lfc_mean`, `lfc_median`, `lfc_std`, `lfc_min`, `lfc_max`;
- optional credible-interval fields only if requested; and
- `is_de_fdr_<fdr_target>`.

With `all_stats=True`, both modes append:

- `raw_mean1`, `raw_mean2` from the registered count matrix;
- `non_zeros_proportion1`, `non_zeros_proportion2`; and
- `raw_normalized_mean1`, `raw_normalized_mean2`, which are mean per-cell count values after linear counts-per-10,000 scaling, not log1p values.

Group-selection calls also append `comparison`, `group1`, and `group2`. An OpenBio adapter should materialize the index as one canonical `gene` column, normalize the dynamic FDR column to stable `is_de_fdr` while retaining the target in metadata, and validate rather than silently accepting version-dependent malformed output.

## Narrow, scientifically explicit reference call

The following is the intended official call shape after a compatible analysis AnnData has been constructed. It deliberately fixes lower-level policy; the proposed node interface is specified in `design-review.md`.

```python
result = trained_model.differential_expression(
    adata=analysis_adata,
    groupby=groupby,
    group1=[group1],
    group2=group2,
    mode="change",
    delta=delta,
    all_stats=True,
    batch_correction=True,
    batchid1=shared_technical_batches,
    batchid2=shared_technical_batches,
    fdr_target=fdr_target,
    silent=True,
    weights="uniform",
    filter_outlier_cells=False,
    n_samples_overall=n_samples_overall,
    test_mode="two",
)
```

This call uses `test_mode="two"` intentionally so `proba_de` represents the two-sided composite practical-change event. That is an OpenBio adapter policy, not the released scvi-tools default. The realized policy must appear in `summary` and generated `code`.

For an explicitly selected observed-batch comparison, the call instead uses `batch_correction=False` and omits both batch-ID arguments. That path must carry a warning that observed batch/covariate composition may contribute to the result.

## Why this is exploratory model evidence, not Sample-level Condition inference

The scVI DE implementation forms aggregate posterior mixtures from selected cells and samples normalized expression through a model shared across cells. It does not receive a `sample_key`, construct one independent response per Sample, estimate between-Sample variance for a Condition coefficient, or model the nesting of cells within Samples.

The distinction remains even when a donor/Sample column was registered as a Technical batch:

- treating Sample as a nuisance may remove Sample-associated biology rather than estimate a Condition effect;
- leaving Sample unregistered does not make cells from the same Sample independent;
- decoding both groups on shared batches standardizes a decoder covariate but does not create biological replication; and
- posterior expected FDR remains conditional on the fitted cell-level generative model.

Squair et al. show that cell-level analyses that ignore biological-replicate variation can produce false discoveries and that adding more cells does not solve this problem. Zimmerman et al. likewise identify cells within one individual as correlated subsamples and emphasize inference at the experimental-unit level.

Accordingly, even if `groupby` contains Condition labels, this node may only report exploratory, model-dependent evidence. The repository's formal Condition contrast remains full-gene count aggregation by Sample within one defined population, followed by a replicate-aware method such as edgeR or DESeq2. The scVI result can generate hypotheses, characterize clusters/states, or provide sensitivity evidence; it cannot replace that Condition analysis.

Additional reasons for cautious language are:

- the result is usually restricted to an HVG representation-learning view;
- the same cells fitted the model and may also have generated the clustering used as `groupby`;
- model fit, nuisance declarations, and counterfactual support are assumptions rather than verified facts;
- the operation is Monte Carlo and device/version dependent; and
- the current workflow artifact is process-local and cannot be reconstructed from emitted code alone.

## Minimum scientific disclosure

Every report must state:

1. that the result is `exploratory_model_evidence`;
2. the exact two groups, optional population restriction, cell counts, missing/other-label exclusions, and ordered feature scope;
3. fitted-model identity/fingerprint, registered count source/state, Technical batch and other nuisance covariates, training parameters, diagnostics, random seed, and device;
4. `mode`, exact hypothesis, `delta`, `test_mode`, realized pseudocount, posterior sample count, weighting, outlier policy, raw-stat policy, and posterior FDR target;
5. observed versus counterfactual batch behavior, the exact Technical batches used, per-group support in each, and whether any extrapolation was allowed;
6. number of tested model features, number posterior-FDR-tagged, signed top effects with posterior probabilities and raw-count summaries, or an explicit empty result;
7. all method/package references and software versions; and
8. limitations: no Sample-level Condition inference, no p-values/BH claim, no full-gene claim unless the model actually used all genes, no integration-quality claim, and conditional reproducibility on the same fitted model.

## Current-node deviations requiring refactor

The current repository node passes only `adata`, `groupby`, two scalar labels, `mode`, and `delta`. It leaves all other scvi-tools behavior implicit. It also:

- defaults to `vanilla` despite the official effect-size recommendation;
- does not expose or disclose batch conditioning, posterior FDR, posterior sample count, `test_mode`, weighting, pseudocount policy, or randomness;
- does not use the integration artifact's training diagnostics in its report;
- does not validate the backend table beyond renaming an `index` column;
- does not emit references, software versions, `summary`, or `code`;
- uses a generic “differential expression” title without the required exploratory qualification; and
- can be misread as Condition inference although its statistical unit is the selected cell/posterior mixture.

Those are adapter deficiencies, not reasons to alter scvi-tools itself.
