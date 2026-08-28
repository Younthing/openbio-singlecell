# OpenBioSingleCellAugur — official usage research

## Audited baseline

The current node is close to Pertpy's call sequence but selects `log1p_norm` by default even though Pertpy 1.3.0
`Augur.load` explicitly checks for raw counts and warns when they are absent. It does not require/audit biological
Sample or Technical batch, exposes a categorical comparison while also offering an inappropriate regressor, hides
important subsampling/fold settings, and returns a backend-produced AnnData whose `uns` contains heterogeneous Python
results. A second node later guesses tables from that mutable `uns` dictionary.

The reviewed implementation and release-extra pin target Pertpy 1.3.0. Runtime accepts the 1.3.x patch family only
when the audited public constructor/load/predict signatures are still exact; the actual executing patch version is
reported. Rejecting a signature-compatible patch solely because it is not 1.3.0 would be an unnecessary expert-tool
restriction.

## Official Pertpy 1.3.0 interface

Primary sources:

- https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Augur.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_augur.py
- https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/augur.html

Official binary-classification usage is:

```python
augur = pertpy.tl.Augur("random_forest_classifier", random_state=random_seed)
loaded = augur.load(
    count_adata,
    label_col=condition_key,
    cell_type_col=population_key,
    condition_label=control,
    treatment_label=treatment,
)
augur_adata, results = augur.predict(
    loaded,
    n_subsamples=50,
    subsample_size=20,
    folds=3,
    feature_perc=0.5,
    var_quantile=0.5,
    span=0.75,
    filter_negative_residuals=False,
    n_threads=n_threads,
    augur_mode="default",
    select_variance_features=True,
    key_added="openbio_augur",
    random_state=random_seed,
    zero_division=0,
)
```

The result contains `summary_metrics`, `full_results`, and `feature_importances`. For classifiers the summary metrics
include mean Augur/AUC, accuracy, precision, F1, and recall. `full_results` records AUC by population, subsample, and
fold. Feature importances/standardized coefficients are model-dependent.

Versioned source review identified two reporting details. First, `load` only filters `condition_label` and
`treatment_label` in the categorical-label path, so OpenBio must explicitly subset and categorize exactly two
Conditions before calling it. Second, Pertpy seeds cell subsampling with `subsample_idx` rather than the supplied
global `random_state`; the supplied seed controls estimator/fold behavior but does not alter that deterministic
subsample sequence. This 1.3.0 limitation must be disclosed, not hidden behind a claim that one seed controls all RNG.

## Scientific interpretation

Augur ranks populations by how well a classifier separates perturbation labels from cell profiles after repeated
balanced cell subsampling. It is exploratory prioritization, not a replicate-aware Condition hypothesis test and not
a p-value. Random cross-validation folds cells, so cells from the same biological Sample can occur in training and
test folds; Sample/donor or Technical batch effects can drive separability. OpenBio uses the requested Sample and
optional Technical-batch columns to audit mappings, support, and confounding, but these are reporting diagnostics:
they must not prevent an expert from running the descriptive cell-level method when its documented cell/subsampling
requirements are met.
Sample identities are caller-declared. If `sample_key` deliberately reuses the population, Condition, or Technical
batch column, the calculation remains available but the Sample audit is marked `role_reused_unverified`; its unique
label counts must not be interpreted as biological replicate counts.

Raw snapshot counts are the recommended expression state, while the neutral UI default remains the current `X`
representation so OpenBio does not silently switch the user's active matrix. Selecting `raw` explicitly is sufficient for
that data contract; Augur validates the selected Raw axes and matrix but does not require or trust mutable
analysis-history entries. The existing explicit `X`/layer selector remains usable for expert-supplied nonnegative
finite expression. A non-integer selected source is disclosed as user-selected rather than rejected or silently cast
to integers.
Neither an AnnData object nor its history can prove that the user-declared snapshot contains every assayed post-QC
feature. Reports therefore use `user_declared_not_programmatically_verifiable` rather than upgrading the source to
“validated full-gene”.
The comparison is restricted to two explicit Conditions and population labels must disclose whether they are
Provisional or Curated. Two biological Samples per arm, stable Sample mappings, and non-confounded Technical batches
are recommended for reviewability, but Augur does not use them as its computation unit. They are therefore warnings,
not execution gates. Formal Condition claims still require a Sample-level method.

## References to emit

- Skinnider MA, et al. Cell type prioritization in single-cell data. *Nature Biotechnology*. 2021;39:30-34.
  https://doi.org/10.1038/s41587-020-0605-1
- Squair JW, et al. Prioritization of cell types responsive to biological perturbations with Augur.
  *Nature Protocols*. 2021;16:3836-3873. https://doi.org/10.1038/s41596-021-00561-x
- Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025.
  https://doi.org/10.1038/s41592-025-02909-7

## Report/code implications

Report the exact explicitly selected expression state, integer-like diagnostic, user-declared feature-completeness
status, Conditions, population/annotation status, Sample and Technical-batch audits, model and all fixed/exposed
Augur settings, eligible/skipped populations/reasons,
full AUC distributions and priority ties, feature-importance caveats, RNG limitation, no-test/no-replicate-aware
status, references, and dynamic versions.
Generated code returns `(results_dict_of_dataframes, summary_dict)` and never stores backend dictionaries in caller
AnnData.
