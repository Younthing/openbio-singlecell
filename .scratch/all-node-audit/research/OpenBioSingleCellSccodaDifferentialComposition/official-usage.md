# OpenBioSingleCellSccodaDifferentialComposition: official usage and primary-source audit

## Audit scope

- Audited: 2026-08-28.
- Repository implementation: `OpenBioSingleCellSccodaDifferentialComposition` in `openbio_singlecell/nodes_abundance.py`.
- Upstream contract checked against Pertpy 1.3.0, released 2026-08-24. PyPI's trusted-publishing provenance identifies source commit [`e2976e13944bbdb3e30db25bb4b231f730e83ea5`](https://github.com/scverse/pertpy/tree/e2976e13944bbdb3e30db25bb4b231f730e83ea5). Facts that depend on implementation details below refer to that immutable source, not an unversioned `latest` branch.
- Sources are restricted to official Pertpy/Patsy/NumPyro/ArviZ documentation and source, the Pertpy paper, and the original scCODA/NUTS/compositional-data papers.

The repository's domain contract is a **Sample-level Condition contrast**. Cells contribute cell-type counts within a biological replicate; they are not independent experimental units. scCODA estimates relative compositional effects, not absolute cell abundance and not cell-level differential expression.

## Current node behavior

The current node exposes `adata`, `sample_key`, `annotation_key`, comma-separated `covariate_keys`, an arbitrary Patsy `formula`, `reference_cell_type`, `estimated_fdr`, NUTS sample/warm-up counts, and `random_seed`. It calls:

1. `Sccoda.load(..., type="cell_level", generate_sample_level=True)` on `adata.copy()`;
2. `prepare(..., formula=..., reference_cell_type=...)`;
3. `run_nuts(...)`;
4. `set_fdr(...)` and `summary_prepare(...)`;
5. concatenation of upstream intercept/effect tables into one generic table.

Only `table` is returned. There is no explicit focal comparison, Sample-level uniqueness/rank check, sampler diagnostic report, strict JSON summary, software provenance, or equivalent source-code output.

This is close to Pertpy's mechanical call sequence but not yet a safe scientific interface. In particular, Pertpy silently skips a requested covariate that is not constant within Sample, and Patsy may choose categorical bases from data order. A formal node must reject those ambiguities before fitting.

## Official data-loading contract

The official [`pertpy.tools.Sccoda`](https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Sccoda.html) API accepts either Sample-level counts or cell-level AnnData. The cell-level route used here is documented in the [official scCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda.html) and implemented in Pertpy's pinned [`_base_coda.py`](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py).

For cell-level input, `load`/`from_scanpy` aggregates observations by Sample identifier and cell-type identifier. The generated MuData contains:

- an `"rna"` modality representing the cell-level input;
- a `"coda"` modality whose `.X` is a Sample-by-cell-type count matrix;
- Sample covariates in `mdata["coda"].obs`;
- cell-type count totals in `mdata["coda"].var["n_cells"]`.

The statistical row is therefore one Sample. A valid `sample_key` must identify the biological experimental unit whose cells were jointly collected. A subject identifier alone is not a valid Sample key when that subject has multiple tissues, time points, treatments, or libraries that constitute distinct observations. Conversely, technical libraries from the same biological unit must not be presented as independent Samples merely to increase replication.

Pertpy allows a column or list of columns as the Sample identifier. The pinned implementation first wraps even a scalar column name in a list, then writes a working-copy `scCODA_sample_id` by applying `"-".join` to the selected values. This assumes joinable string values and can collapse distinct composite tuples to the same text. OpenBio should require one upstream-defined canonical `sample_key`, validate its typed values, construct its own collision-free internal codes on a temporary copy, retain the original values for disclosure, and verify that each Sample has exactly one value for every modeled covariate.

The pinned loader checks whether covariates are unique within each Sample. When they are not, it logs and skips the covariate rather than failing. That is unsuitable for an auditable node: silently fitting a different design is a correctness failure. OpenBio must preflight and reject any Sample that maps to multiple Condition or adjustment values.

Counts, rather than proportions, are required. During `prepare`, Pertpy converts the Sample count matrix to `float64`, warns if rows resemble proportions, replaces zeros with the fixed pseudocount `0.5`, and uses the resulting row totals. The report must distinguish observed integer cell counts from the internal pseudocount-adjusted values.

## Formula, covariates, and comparison semantics

Pertpy builds a design matrix with [`patsy.dmatrix`](https://patsy.readthedocs.io/en/latest/formulas.html). It removes the intercept column and records the remaining covariate names. Categorical variables use Patsy's treatment coding; an explicit base can be expressed with [`C(variable, Treatment(reference=...))`](https://patsy.readthedocs.io/en/latest/categorical-coding.html). The pinned Pertpy docstring also notes that an implicit categorical reference can follow the first Sample/category ordering.

An arbitrary formula is useful in a programming API but too broad for this workflow node:

- it can encode multiple Condition effects, interactions, splines, or transformations, so one node run no longer answers one atomic comparison;
- the current `covariate_keys` and `formula` can disagree;
- category order can silently reverse the reported direction;
- formula names and levels can be syntactically unsafe;
- the concatenated result does not identify which row is the intended writing contrast.

For a formal two-condition node, the primary contract should instead name `condition_key`, `reference_condition`, and `comparison_condition`, with optional declared adjustment covariates. The implementation should build a safe internal design with the exact comparison coefficient, then verify that the design is full rank. Other Condition levels must be rejected with an instruction to subset upstream rather than silently included in a multi-level analysis. The summary must state direction as “comparison relative to reference”.

Adjustment covariates are Sample-level nuisance variables, not additional focal results. Numeric covariates must be finite. Categorical adjustments must have a deterministic, disclosed base; an ordered categorical declaration is preferable to order-of-first-observation behavior. Perfect confounding or rank deficiency is a preflight error. scCODA does not add random effects, so paired/repeated-measures data require an explicit supported design policy and must not be described as independent replicates merely because a subject column was added.

## Cell-type reference

scCODA is identifiable relative to one cell type. The selected reference's covariate effect is fixed to zero, and all other effects are interpreted relative to it. An explicit, biologically defensible reference is preferable.

`reference_cell_type="automatic"` is official functionality. Pertpy's pinned source selects, among cell types with a zero fraction below `automatic_reference_absence_threshold=0.05`, the one with the smallest dispersion in relative abundance. The implementation threshold therefore means presence in more than 95% of Samples. One Pertpy docstring still says 90%; the executable 0.05 rule and original scCODA/tascCODA method description are the auditable contract.

Automatic selection is data-dependent. If retained, the node must disclose the resolved reference, zero fraction, relative-abundance dispersion, candidate set, and threshold. No eligible reference is a hard error. The report must not imply that an automatically stable reference is known to be biologically invariant.

## Model and inference

The original scCODA paper models Sample-by-cell-type counts with a Dirichlet-multinomial likelihood and a log-linear covariate model. A spike-and-slab prior supports sparse cell-type effects, while a fixed reference cell type makes the compositional parameterization identifiable. Effects are relative to the reference and to the remaining composition. Increasing one component necessarily changes the normalized shares of others.

The official Pertpy workflow recommends NUTS, the No-U-Turn Sampler introduced by Hoffman and Gelman. The pinned signature defaults to:

```text
num_samples = 10000
num_warmup = 1000
rng_key = 0
```

`rng_key` seeds both NumPy initialization and JAX sampling in Pertpy. Reproducibility still requires identical input ordering, dependency versions, numerical backend/device, precision policy, and model settings; a seed alone is not a cross-platform bitwise guarantee.

Pertpy forwards additional keyword arguments to NumPyro's NUTS kernel. The released wrapper constructs NumPyro `MCMC` without setting `num_chains`, so the NumPyro default is one chain. Consequently:

- R-hat cannot establish between-chain convergence;
- a successful return must not be summarized as `converged=true`;
- within-chain effective sample size and Monte Carlo standard error can be reported, but are not a substitute for independent chains;
- OpenBio should report a stable status such as `limited_single_chain`, with stronger warnings when acceptance or within-chain diagnostics are poor.

The pinned Pertpy runner stores potential energy, number of integration steps, adapted step size, acceptance probability, and mean acceptance probability. It does not request NumPyro's `diverging` field. A divergence count is therefore unavailable in Pertpy 1.3.0 and must be reported as unavailable, never inferred as zero. Pertpy warns when mean acceptance is below 0.6 or above 0.95; those official bounds should be preserved in the node's diagnostic policy.

Although HMC exists in the upstream modeling layer, a sampler switch would enlarge this node into multiple inference methods with different tuning contracts. The atomic node should use version-pinned NUTS and disclose the fixed kernel defaults. If another sampler is later required, it deserves a separately reviewed interface.

Official references:

- [Pertpy scCODA API](https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Sccoda.html)
- [Official extended scCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda_extended.html)
- [Pinned Pertpy 1.3.0 CODA source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py)
- [NumPyro MCMC documentation](https://num.pyro.ai/en/stable/mcmc.html)
- Hoffman MD, Gelman A. The No-U-Turn Sampler. *JMLR*. 2014;15:1593-1623. [Primary paper](https://jmlr.org/papers/v15/hoffman14a.html).

## Credible effects and estimated FDR

scCODA does not return classical p-values and does not apply a Benjamini-Hochberg correction. In Pertpy 1.3.0, the classic spike-and-slab path computes each coefficient's posterior inclusion probability as the fraction of sampled `beta` values whose absolute magnitude exceeds `1e-3`. For each candidate cutoff it evaluates the mean of `1 - inclusion_probability` among included effects, takes the first cutoff strictly below the target, floors that cutoff to three decimal places, and falls back to `1.0` when no candidate qualifies. It then sets excluded final parameters to zero. Because flooring happens after the candidate check, OpenBio should recompute and report the realized posterior expected-FDR quantity at the applied threshold instead of assuming it exactly equals the requested target.

`estimated_fdr` is therefore a target for posterior expected false discoveries, not a frequentist guarantee. Results should be called **credible effects under the selected posterior expected-FDR threshold**, not “statistically significant genes/cell types”. The report must include:

- requested `estimated_fdr`;
- realized inclusion-probability threshold;
- realized posterior expected-FDR quantity after applying the rounded threshold;
- posterior inclusion probability per cell type;
- the exact rule/version used;
- the number and identities of credible effects;
- reference cell type and comparison direction.

Changing `estimated_fdr` changes selection and is a visible result-defining parameter. It should be constrained strictly inside `(0, 1)` under OpenBio policy; endpoint values are technically accepted upstream but do not represent useful reporting choices.

## Official result fields and interpretation

`summary_prepare` returns an intercept table and an effect table. In Pertpy 1.3.0 the displayed intercept fields include:

- `Final Parameter`;
- HDI lower/upper bounds;
- posterior `SD`;
- `Expected Sample`.

Effect rows are indexed by covariate and cell type and include:

- `Final Parameter`;
- HDI lower/upper bounds, calculated conditionally for included spike-and-slab samples;
- posterior `SD`;
- `Inclusion probability`;
- `Expected Sample`;
- compositional `log2-fold change`.

The compositional log2-fold change is calculated after exponentiating and renormalizing the composition to the same average total. It can differ in sign or magnitude from the raw model coefficient, and a cell type whose selected coefficient is zero can still receive a nonzero compositional fold change because other components changed. The canonical OpenBio table must label both quantities and must not treat expected cell counts as absolute abundance estimates.

The current node concatenates intercept and effect sections with heterogeneous columns. For one atomic Condition contrast, the primary table should instead contain only the focal effect rows in a stable schema; intercept and nuisance-covariate summaries belong in structured model details, not in the writing-ready result table.

## Dependency and version contract

Pertpy 1.3.0 requires Python 3.12 or newer and publishes a `jax` extra. Its tagged project metadata includes AnnData, Scanpy, MuData, ArviZ `>=1.2.0`, and an optional JAX stack containing NumPyro, JAX-related optimization packages, and their transitive dependencies. Installation for scCODA is documented conceptually as `pertpy[jax]`.

OpenBio currently targets Python 3.12 and pins AnnData `>=0.13.2,<0.14`, but does not declare Pertpy, JAX, JAXlib, NumPyro, MuData, or ArviZ as project dependencies. A production implementation must therefore add and test an explicit optional dependency set rather than rely on an incidental environment. For this adapter, pinning and testing `pertpy==1.3.0` is preferable to an open-ended lower bound because output fields and diagnostics are source-sensitive.

Every summary must report actual imported versions of at least:

- OpenBio SingleCell and Python;
- Pertpy and its source/release identifier where available;
- NumPyro, JAX, and JAXlib;
- ArviZ, MuData, AnnData, Scanpy, Patsy, pandas, NumPy;
- backend/device and whether JAX x64 was enabled.

Normative version sources:

- [Pertpy 1.3.0 on PyPI](https://pypi.org/project/pertpy/1.3.0/)
- [Pertpy 1.3.0 tagged project metadata](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/pyproject.toml)
- Heumos L, et al. Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025. [doi:10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7).

## Auditable official-use sequence

The upstream sequence, adapted to make all relevant arguments explicit, is:

```python
import pertpy as pt

model = pt.tl.Sccoda()
mdata = model.load(
    adata.copy(),
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier=annotation_key,
    sample_identifier=sample_key,
    covariate_obs=sample_covariates,
)
mdata = model.prepare(
    mdata,
    modality_key="coda",
    formula=explicit_formula,
    reference_cell_type=reference_cell_type,
)
model.run_nuts(
    mdata,
    modality_key="coda",
    num_samples=num_samples,
    num_warmup=num_warmup,
    rng_key=random_seed,
)
model.set_fdr(mdata, estimated_fdr, modality_key="coda")
intercept_df, effect_df = model.summary_prepare(
    mdata["coda"], est_fdr=estimated_fdr
)
```

This call sequence is necessary but not sufficient. OpenBio must add Sample semantics, design/reference validation, output normalization, diagnostic extraction, strict reporting, and version gating around it.

## Minimum scientific disclosure

A report suitable for writing must disclose:

- biological Sample definition and per-condition Sample counts;
- cell-type annotation field and whether the annotation is curated or provisional;
- Sample-by-cell-type table dimensions, observed zeros, total cells per Sample, and zero pseudocount policy;
- exact Condition comparison direction and every adjustment covariate with its coding/base;
- design-matrix dimensions, rank, and any collinearity warning;
- explicit or automatically resolved reference cell type and its selection diagnostics;
- scCODA/Pertpy version, NUTS draws, warm-up, seed, backend, and single-chain limitation;
- mean acceptance and available ESS/MCSE/energy/step diagnostics; divergences as unavailable under this pinned adapter;
- posterior expected-FDR target and realized inclusion-probability threshold;
- writing-ready credible effects with inclusion probabilities, relative coefficients, HDIs, and compositional fold changes;
- the limitations that effects are relative, cells are not replicates, low Sample replication limits inference, annotation error propagates into counts, and the model does not prove absolute abundance changes.

## Primary references

- Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for compositional single-cell data analysis. *Nature Communications*. 2021;12:6876. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8616929/).
- Aitchison J. The statistical analysis of compositional data. *Journal of the Royal Statistical Society: Series B*. 1982;44:139-160. [doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x).
- [Official Pertpy scCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda.html).
- [Official Pertpy tutorial source repository](https://github.com/scverse/pertpy-tutorials).
- [Pinned Pertpy `Sccoda` source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_sccoda.py).
- [Pinned Pertpy shared CODA source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py).
