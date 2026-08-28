# OpenBioSingleCellTasccodaDifferentialComposition: official usage and primary-source audit

## Audit scope

- Audited: 2026-08-28.
- Repository implementation: `OpenBioSingleCellTasccodaDifferentialComposition` in `openbio_singlecell/nodes_abundance.py`.
- Upstream contract checked against Pertpy 1.3.0, released 2026-08-24, with trusted-publishing provenance at source commit [`e2976e13944bbdb3e30db25bb4b231f730e83ea5`](https://github.com/scverse/pertpy/tree/e2976e13944bbdb3e30db25bb4b231f730e83ea5).
- Sources are restricted to official Pertpy/Patsy/NumPyro/ArviZ documentation and source, the Pertpy paper, and the original tascCODA/scCODA/NUTS/compositional-data papers.

tascCODA is a Sample-level hierarchical compositional model. Its input tree is not presentation metadata: it changes the parameterization and sparse selection problem. Results describe relative effects on hierarchy nodes and their propagation to leaves, not absolute cell abundance.

## Current node behavior and critical discrepancies

The current node adds `hierarchy_levels`, `phi`, `lambda_1`, and `estimated_fdr` to the same broad `covariate_keys`/arbitrary-`formula` interface used by the scCODA wrapper. It calls `Tasccoda.load(..., type="cell_level", levels_orig=levels, add_level_name=True)`, prepares the model with `pen_args={"phi": phi, "lambda_1": lambda_1}`, runs NUTS, and concatenates intercept, leaf-effect, and hierarchy-node frames. It returns only `table`.

Four discrepancies are result-defining:

1. The UI constrains `phi >= 0`, but the original tascCODA paper defines meaningful negative values that favor effects on larger subtrees nearer the root.
2. The node defaults `lambda_1=3.5`; the paper and Pertpy 1.3.0 default to `lambda_1=5`. There is no local rationale or disclosure for the changed prior.
3. The node exposes `estimated_fdr`, but Pertpy 1.3.0 explicitly documents that `est_fdr` has no effect for the spike-and-slab LASSO selection used by tascCODA and that `set_fdr` does not work for it.
4. `hierarchy_levels` is passed through without proving root-to-leaf ordering, one parent per child, leaf equality with `annotation_key`, name uniqueness, or reference-path validity.

The current table adds a generic `summary_section` label but still mixes three statistical scopes in one heterogeneous schema without normalizing their distinct selection semantics for writing. It also omits hierarchy fingerprint, penalty settings, selection threshold (`Delta`), diagnostics, versions, references, and code.

## Official cell-level loading and MuData contract

The official [`pertpy.tools.Tasccoda`](https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Tasccoda.html) API supports cell-level or already aggregated Sample-level AnnData. In the pinned source, cell-level `load`:

- aggregates cells to a Sample-by-cell-type count matrix using the shared scCODA loader;
- returns MuData with cell-level `"rna"` and aggregated `"coda"` modalities;
- copies requested Sample covariates into `mdata["coda"].obs`;
- imports a hierarchy under `mdata["coda"].uns[key_added]`, defaulting to `"tree"`.

The API states that `levels_orig` for cell-level data must begin at the root level and end at the leaf level. `add_level_name=True` prefixes internal node names with their level name to reduce collisions. See the pinned [`Tasccoda.load` source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_tasccoda.py) and the [official tascCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/tasccoda.html).

As with scCODA, the experimental row is one biological Sample, not one cell or one technical library. Each Sample must have exactly one Condition and one value per adjustment covariate. The pinned shared loader wraps even one Sample column in a list and creates `scCODA_sample_id` with `"-".join`; nonstring values and ambiguous composite concatenations are therefore unsafe. OpenBio should build collision-free internal Sample codes on a copy while preserving original typed identifiers. The loader also logs and skips covariates that are not constant within Sample, which is unsafe for a formal node; OpenBio must reject them before model fitting.

The input to the statistical model is an `n Samples × p leaves` count matrix, not proportions. Pertpy converts counts to `float64`, replaces zeros with the fixed `0.5` pseudocount during preparation, and uses the adjusted row totals. The original tascCODA paper emphasizes that row totals reflect sampling/sequencing scale and that only relative compositional information is available.

## Hierarchy contract

The original method accepts a bifurcating or multifurcating rooted tree whose `p` leaves exactly correspond to modeled cell types. An ancestor matrix maps every hierarchy-node effect to its descendant leaves. Effects placed on an internal node propagate equally to its descendant leaves; without sparse regularization, identical leaf effects and one ancestral effect can represent the same fitted composition.

For cell-level Pertpy loading, every cell must carry a complete root-to-leaf path through the columns supplied as `levels_orig`. A production wrapper must validate more than column existence:

- the sequence is ordered root to leaf;
- the final leaf values equal the `annotation_key` categories exactly;
- every leaf maps to one and only one path;
- every child at each level maps to exactly one parent;
- no path component is missing or empty;
- all modeled leaves are reachable from a single root and no unused category creates a phantom node;
- prefixed node names remain unique and cannot collide with leaf names or reserved suffixes;
- deterministic level/category ordering yields a reproducible tree/ancestor-matrix fingerprint;
- the chosen reference leaf and all of its ancestors exist and are identified.

To remove current ambiguity, a public OpenBio input should define `hierarchy_keys` as ancestor columns from root through immediate parent; the implementation appends `annotation_key` as the leaf. Requiring users to repeat the leaf key in a comma-separated hierarchy list creates an easy off-by-one semantic error.

The paper describes the method as most suitable for moderate feature dimensions, typically fewer than 100 leaves. This is not a universal hard theorem, but it warrants a resource/scope warning and explicit reporting when the modeled leaf count approaches or exceeds that regime.

## Formula, Condition contrast, and reference leaf

Pertpy delegates design construction to Patsy, removes the intercept column, and accepts a broad R-style formula. Treatment bases can be expressed explicitly with `C(..., Treatment(reference=...))`; otherwise categorical reference behavior can follow category/Sample ordering. The same atomicity concerns as scCODA apply: an arbitrary formula can create multiple covariates, interactions, or implicit directions.

For the workflow node, one run should represent one explicitly directed two-condition contrast with optional Sample-level adjustments. `condition_key`, `reference_condition`, and `comparison_condition` should replace the arbitrary formula as the focal interface. The internal design must be finite and full rank, and exactly one result coefficient must map back to the intended comparison.

tascCODA also requires a reference leaf for compositional identifiability. The chosen leaf's effect and the effects of all its ancestors are fixed to zero. Consequently, changing the reference can change the numerical parameterization and interpretation of every other node/leaf effect. The summary must report the exact reference leaf and complete reference path.

`reference_cell_type="automatic"` uses the same pinned Pertpy rule as scCODA: candidates have zero fraction below `0.05` (present in more than 95% of Samples), then the least-dispersed relative-abundance candidate is chosen. It is data-selected and must be disclosed with candidate diagnostics; it is not proof of biological invariance.

## Tree-adaptive prior and `phi`

tascCODA places a spike-and-slab LASSO prior on hierarchy-node effects. Pertpy 1.3.0 documents and implements these defaults:

```text
lambda_0 = 50
lambda_1 = 5
phi = 0
theta = 0.5
```

`lambda_0` strongly shrinks spike effects. `lambda_1` controls slab penalization. Pertpy scales `lambda_1` for each node using its descendant-leaf count and `phi`:

- `phi = 0`: all nodes receive equal slab penalization;
- `phi < 0`: larger subtrees nearer the root are penalized less, favoring more aggregated effects;
- `phi > 0`: detailed leaf-level effects are favored.

The original paper calls `phi` an aggregation-bias trade-off between generalization/parsimonious higher-level effects and detailed predictive effects. Simulations show that a misspecified value can increase false discoveries; the authors explicitly recommend cross-validation over different `phi` values in practice. Therefore `phi` is not a cosmetic advanced knob. A single-fit node may expose a signed `aggregation_bias`, but its report must say whether the value was pre-specified and recommend a separate sensitivity/CV workflow before confirmatory writing. The node must not secretly select `phi` on the same full data and then report the fit as if it were pre-specified.

The current `lambda_1=3.5` must not survive as a hidden compatibility default. Unless OpenBio implements and validates a separately documented prior-sensitivity interface, keep `lambda_0=50`, `lambda_1=5`, and `theta=0.5` fixed to the pinned Pertpy contract and disclose them. Exposing several correlated prior knobs makes the node a tuning workbench rather than one atomic analysis.

The original paper assigns the global mixture weight `theta` a Beta prior. Pertpy release 1.1.1 fixed two tascCODA defects: theta collapse causing no credible effects and results always being marked insignificant. In 1.3.0, `theta` is instead held fixed at `0.5`; the pinned source explains that its posterior was weakly identified under HMC/NUTS and collapsed toward zero, inflating the selection threshold. This is a material, release-specific departure from the published generative specification and must be reported. Versions earlier than 1.1.1 are scientifically unsafe for this adapter, and pinning/testing 1.3.0 is the stronger policy.

Primary sources:

- Ostner J, Carcy S, Müller CL. tascCODA: Bayesian Tree-Aggregated Analysis of Compositional Amplicon and Single-Cell Data. *Frontiers in Genetics*. 2021;12:766405. [doi:10.3389/fgene.2021.766405](https://doi.org/10.3389/fgene.2021.766405); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8689185/).
- [Pinned Pertpy tascCODA source and prior defaults](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_tasccoda.py).
- [Pertpy 1.1.1 official release notes](https://github.com/scverse/pertpy/releases/tag/1.1.1).

## NUTS, seed, and diagnostic limitations

tascCODA uses the same Pertpy NUTS runner as scCODA. The pinned defaults are 10,000 posterior draws, 1,000 warm-up draws, and `rng_key=0`. The seed is used for NumPy initialization and JAX sampling. Identical results also depend on input ordering, dependency versions, backend/device, and precision.

The Pertpy wrapper constructs NumPyro `MCMC` without setting `num_chains`, so the official NumPyro default is one chain. R-hat is unavailable as a convergence assessment. The runner records potential energy, number of integration steps, adapted step size, acceptance probability, and mean acceptance probability, but it does not request the `diverging` field. Reports must therefore state:

- chain count `1` and a single-chain limitation;
- no assertion of multi-chain convergence;
- mean acceptance, with Pertpy warnings outside `[0.6, 0.95]`;
- within-chain ESS/MCSE if computed safely through ArviZ;
- available energy/step summaries;
- `divergences_available=false`, not a fabricated zero count.

HMC/NUTS behavior is especially important in tascCODA because the release history contains a prior-identification failure. The node should keep NUTS as the fixed, versioned inference method rather than expose a sampler switch. Nonfinite posteriors, missing node effects, or malformed shapes are hard failures; poor available diagnostics produce a prominent warning and suppress unqualified writing claims.

Official sources:

- [Pinned Pertpy shared NUTS source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py)
- [NumPyro MCMC documentation](https://num.pyro.ai/en/stable/mcmc.html)
- [ArviZ summary API](https://python.arviz.org/en/stable/api/generated/arviz.summary.html)
- Hoffman MD, Gelman A. The No-U-Turn Sampler. *JMLR*. 2014;15:1593-1623. [Primary paper](https://jmlr.org/papers/v15/hoffman14a.html).

## Credible effects and why `estimated_fdr` must be removed

tascCODA does not use scCODA's inclusion-probability expected-FDR rule. Pertpy identifies node effects under the spike-and-slab LASSO using a node-specific `Delta` threshold. A node's selected final parameter is its posterior median when the absolute median exceeds `Delta`; otherwise it is zero. `Is credible` reflects that rule.

The pinned shared source says both:

- `set_fdr` does not work for the spike-and-slab LASSO selection method;
- the `est_fdr` parameter has no effect for that selection method.

Therefore the current tascCODA `estimated_fdr` input is misleading and must be removed, not hidden. A report must not claim FDR control at the user-entered value and must not use the scCODA term “inclusion probability threshold” for tascCODA node selection.

The original paper reports empirical FDR in simulations, but that benchmark metric is not an adjustable per-analysis FDR control. Result language should be “credible hierarchy-node effects under the pinned tascCODA spike-and-slab LASSO/Delta rule,” not “FDR-adjusted significant effects.”

Normative source: [pinned Pertpy `set_fdr` and `credible_effects` implementation](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py).

## Official output fields and two statistical scopes

Pertpy 1.3.0 returns three frames:

1. Intercepts by leaf: `Final Parameter`, HDI bounds, `SD`, `Expected Sample`.
2. Leaf effects by covariate/cell type: `Effect`, posterior `Median`, HDI bounds, `SD`, `Expected Sample`, compositional `log2-fold change`.
3. Node effects by covariate/node: `Final Parameter`, `Median`, HDI bounds, `SD`, `Delta`, `Is credible`.

The node table is the direct sparse-selection evidence. Leaf `Effect` values are obtained by multiplying selected node effects through the ancestor matrix; they are sums propagated from credible nodes, not independently selected leaf discoveries. A canonical output must preserve this distinction with an explicit scope field such as `hierarchy_node` versus `derived_leaf`.

Compositional leaf log2-fold changes are calculated after renormalization to a common average total. They are relative changes and may differ in sign/magnitude from raw node or leaf coefficients. Expected counts are model-derived values under the documented covariate profile, not proof of absolute tissue abundance.

The writing report should foreground credible node effects, name their level and descendant leaves, and then describe leaf-level consequences as derived. Calling every affected descendant a separate credible discovery would inflate and misstate the selected evidence.

## Dependency and version contract

Pertpy 1.3.0 requires Python `>=3.12`. Its tagged metadata defines:

- core dependencies including AnnData, Scanpy, MuData, and ArviZ `>=1.2.0`;
- a `jax` extra with NumPyro and JAX ecosystem packages;
- a `tcoda` extra that includes `pertpy[jax]`, `toytree>=3.0`, `ete4`, and `pyqt6`.

OpenBio currently does not declare Pertpy, JAX/JAXlib, NumPyro, MuData, ArviZ, toytree, or ete4. tascCODA must therefore have an explicit optional dependency installation and a version gate. Pin and test `pertpy==1.3.0`; at absolute minimum, reject Pertpy earlier than 1.1.1 because of the documented credibility defects.

Report actual versions of OpenBio, Python, Pertpy, NumPyro, JAX/JAXlib, ArviZ, MuData, AnnData, Scanpy, Patsy, pandas, NumPy, toytree, and ete4, plus backend/device/x64 state. PyQt is an installation dependency of Pertpy's `tcoda` extra but is not a scientific result dependency; report it only if the actual tree path imports/uses it.

Normative version sources:

- [Pertpy 1.3.0 on PyPI](https://pypi.org/project/pertpy/1.3.0/)
- [Pertpy 1.3.0 tagged project metadata](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/pyproject.toml)
- Heumos L, et al. Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025. [doi:10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7).

## Auditable official-use sequence

An explicit upstream call sequence for cell-level input is:

```python
import pertpy as pt

model = pt.tl.Tasccoda()
mdata = model.load(
    adata.copy(),
    type="cell_level",
    cell_type_identifier=annotation_key,
    sample_identifier=sample_key,
    covariate_obs=sample_covariates,
    levels_orig=[*ancestor_keys, annotation_key],
    add_level_name=True,
    key_added="tree",
)
mdata = model.prepare(
    mdata,
    modality_key="coda",
    formula=explicit_formula,
    reference_cell_type=reference_cell_type,
    tree_key="tree",
    pen_args={"phi": aggregation_bias, "lambda_0": 50, "lambda_1": 5, "theta": 0.5},
)
model.run_nuts(
    mdata,
    modality_key="coda",
    num_samples=num_samples,
    num_warmup=num_warmup,
    rng_key=random_seed,
)
intercept_df, leaf_df, node_df = model.summary_prepare(mdata["coda"])
```

`summary_prepare` should not be passed a user-visible FDR value for tascCODA. OpenBio must wrap this official sequence with Sample/hierarchy/design validation, diagnostic extraction, stable node-versus-leaf output normalization, strict reporting, and dependency checks.

## Minimum scientific disclosure

A writing-ready tascCODA result must state:

- biological Sample definition, comparison direction, per-group Sample counts, and adjustment coding;
- annotation field/status and Sample-by-leaf count-table dimensions/totals/zeros/pseudocount;
- hierarchy columns in root-to-leaf order, tree fingerprint, leaf/internal-node counts, multifurcation summary, and validation result;
- reference leaf and every reference ancestor;
- signed `phi`/aggregation bias, whether it was pre-specified, any CV/sensitivity workflow, and fixed `lambda_0=50`, `lambda_1=5`, `theta=0.5`;
- Pertpy/tascCODA version, NUTS draws/warm-up/seed/backend, single-chain limitation, and available diagnostics;
- credible node effects with level, descendants, median, HDI, SD, Delta, and direction;
- derived leaf effects/LFCs clearly labeled as propagation from selected nodes;
- no claim that `estimated_fdr` controlled this fit;
- limitations: relative composition only, cells are not replicates, low Sample counts, hierarchy/annotation uncertainty, `phi` sensitivity and possible misspecification, moderate-feature regime, no random effects, and lack of demonstrated multi-chain convergence.

## Primary references

- Ostner J, Carcy S, Müller CL. tascCODA: Bayesian Tree-Aggregated Analysis of Compositional Amplicon and Single-Cell Data. *Frontiers in Genetics*. 2021;12:766405. [doi:10.3389/fgene.2021.766405](https://doi.org/10.3389/fgene.2021.766405); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8689185/).
- Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for compositional single-cell data analysis. *Nature Communications*. 2021;12:6876. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6).
- Aitchison J. The statistical analysis of compositional data. *JRSS B*. 1982;44:139-160. [doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x).
- [Official Pertpy tascCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/tasccoda.html).
- [Official Pertpy tutorial source repository](https://github.com/scverse/pertpy-tutorials).
- [Pinned Pertpy `Tasccoda` source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_tasccoda.py).
- [Pinned Pertpy shared CODA source](https://github.com/scverse/pertpy/blob/e2976e13944bbdb3e30db25bb4b231f730e83ea5/src/pertpy/tools/_coda/_base_coda.py).
