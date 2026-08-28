# Pseudobulk aggregation: official usage research

Researched: 2026-08-28

## Open expert boundary correction (normative)

The analyst's explicit `Raw`, `X`, or named-layer selection is the source action. OpenBio validates the selected matrix as a finite, non-negative, integer-like count matrix aligned to the selected observation/feature axes, but it does **not** require OpenBio Snapshot history, equality to `adata.raw.X`, a current/Raw full-gene-axis match, or software provenance that claims to prove biological completeness. Those properties cannot be established from numeric equality or workflow history. The report records the selected source and feature axis as caller declarations and warns that raw-count identity, assay completeness, QC suitability, population annotation, and biological Sample identity remain analyst responsibilities.

Likewise, Curated annotation is recommended for formal population interpretation but is not required to aggregate counts. `formal` versus `exploratory`, annotation status, and role keys are disclosure fields rather than authorization gates. Suspicious role aliases and weak profile/replicate strata are retained when aggregation is well-defined, with `formal_interpretation_invalid` and actionable warnings where appropriate. Hard failures are reserved for computation-critical invariants: invalid/unaligned matrix axes, non-finite/negative/non-integer counts, overflow, invalid labels, a Sample with internally contradictory Sample-level metadata, malformed backend output, or no retained profile.

For integer-sum safety, overflow must be decided from each profile's **exact non-negative Python-integer total**, not from a conservative `maximum_value * number_of_cells` bound. Because all entries have already been validated as non-negative integers, a profile total at or below signed-int64 maximum also proves every per-gene sum is representable. Values such as one large count plus zeros are valid and must not be rejected merely because a worst-case product could overflow; an exact profile total above signed-int64 maximum remains a hard computational failure before NumPy or Decoupler aggregation.

Decoupler 2.2 allocates its aggregate matrix and `psbulk_counts` QC as float64, so integer agreement above `2^53` cannot be proven from its returned values. OpenBio's pre-backend Python-integer totals and subsequent int64 per-gene sums are therefore authoritative for the emitted artifact and profile filtering. The adapter cross-checks Decoupler exactly through float64's exact-integer range and requires rounding-equivalent values above it, while reporting every precision-limited profile; it must never describe those high-count backend values as an exact independent check.

## Scientific scope

Pseudobulk differential expression aggregates raw gene counts from cells belonging to the same independent biological **Sample** and the same defined population, then fits a bulk RNA-seq count model to those Sample-level profiles. The inferential unit is the Sample, never the cell. Summing cells within Sample prevents the cell-level pseudoreplication that otherwise makes nominal sample size and significance anticonservative.

The repository domain vocabulary defines report meaning, while the caller remains responsible for the biological declarations:

- `Sample` is an independent biological replicate and the unit of Condition inference;
- `Technical batch` should be declared as a separate sample-level nuisance role, never silently used as replicate identity or a surrogate Condition; intentional role aliases remain executable but are reported as interpretation-invalid;
- Curated population annotation is preferred for formal population-specific interpretation; provisional or unknown status remains executable but must be disclosed as interpretation-limiting;
- a Condition contrast is tested within exactly one defined population downstream.

The following studies provide the general-practice basis:

- Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*. 2021;12:738. https://doi.org/10.1038/s41467-021-21038-1
- Crowell HL, Soneson C, Germain P-L, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications*. 2020;11:6077. https://doi.org/10.1038/s41467-020-19894-4

## Current official Decoupler API

The current stable Decoupler release researched was 2.2.0, released 2026-07-11. Its public aggregation signature is:

```python
decoupler.pp.pseudobulk(
    adata,
    sample_col,
    groups_col,
    layer=None,
    raw=False,
    empty=False,
    mode="sum",
    skip_checks=False,
    bsize=250_000,
    verbose=False,
) -> anndata.AnnData
```

Official usage for count-model differential expression is:

```python
import decoupler as dc

pdata = dc.pp.pseudobulk(
    adata,
    sample_col="sample",
    groups_col="cell_type",
    layer="counts",
    mode="sum",
)
dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000, inplace=True)
```

The API expects raw integer counts and sums by Sample and group. Decoupler 2.2.0 emits Cartesian group-by-Sample candidates, including synthetic zero rows for empty combinations; OpenBio must identify and remove those empty candidates before constructing the represented-profile artifact. It also emits total-count QC, cell-count QC, and the per-gene fraction of contributing cells with nonzero values. Current 2.2.0 source has an upstream naming inconsistency worth insulating callers from: the docstring says `obs["psbulk_n_cells"]`, but the implementation and `filter_samples` use `obs["psbulk_cells"]`; both use `obs["psbulk_counts"]` and `layers["psbulk_props"]`. OpenBio must expose stable canonical artifact fields rather than make its public contract depend on either Decoupler spelling. The node defaults to the explicit named `counts` layer; experts may select Raw or X directly, and selecting Raw adds no history or equality gate.

`decoupler.pp.filter_samples(adata, min_cells=10, min_counts=1000, inplace=True)` removes low-information Sample-by-population profiles. The official vignette calls these thresholds dataset-dependent and arbitrary; 10 cells/1,000 counts is a rule of thumb, not a universal validity boundary. Report both requested thresholds and the profiles removed.

Feature filtering is a downstream, population-specific model concern. Current Decoupler exposes:

```python
decoupler.pp.filter_by_expr(
    adata,
    group=None,
    lib_size=None,
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
    inplace=True,
)
```

This is adapted from edgeR `filterByExpr`. It should be applied after selecting one population and the modeled samples, because the design and library sizes define whether a gene has sufficient expression. Aggregation must keep the full aligned gene axis and must not silently run a global HVG or detection-proportion filter.

Primary sources:

- Decoupler pseudobulk API: https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.pseudobulk.html
- Decoupler sample filter API: https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.filter_samples.html
- Decoupler expression filter API: https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.filter_by_expr.html
- Official single-cell pseudobulk vignette: https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_psbk.html
- Decoupler 2.2.0 source captured at tag `v2.2.0`, commit `7e8e957cbbb2230079cdd3c13a0ac114a67b9655`: https://github.com/scverse/decoupler/tree/v2.2.0
- Badia-I-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances*. 2022;2(1):vbac016. https://doi.org/10.1093/bioadv/vbac016

## Current Pertpy API and why it is not the preferred aggregation seam

Pertpy 1.3.0, released 2026-08-24, still exposes:

```python
pertpy.tools.PseudobulkSpace.compute(
    adata,
    *,
    target_col="perturbation",
    groups_col=None,
    layer_key=None,
    embedding_key=None,
    mode="sum",  # also count_nonzero, mean, var, median
) -> anndata.AnnData
```

It delegates to `scanpy.get.aggregate` and preserves only observation columns that are constant within each aggregate. The class belongs to Pertpy's perturbation-space utilities and its generic modes include non-count summaries. By contrast, Pertpy's own current differential-expression tutorial explicitly uses `decoupler.pp.pseudobulk(..., mode="sum")` followed by `decoupler.pp.filter_samples` before EdgeR or PyDESeq2. That is the official family workflow the node should follow.

- Pertpy differential-expression tutorial: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html
- Pertpy PseudobulkSpace API: https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.PseudobulkSpace.html
- Pertpy 1.3.0 source captured at tag `1.3.0`, commit `e2976e13944bbdb3e30db25bb4b231f730e83ea5`: https://github.com/scverse/pertpy/tree/1.3.0
- Heumos L, Schaar AC, Lance C, et al. Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025. https://doi.org/10.1038/s41592-025-02909-7

## Count-state and metadata assumptions

The aggregation source is the analyst-declared Raw, `X`, or named layer and must satisfy the computation-critical conditions below before aggregation:

- nonempty cells and genes with unique observation and variable identifiers;
- numeric, finite, non-negative, integer-like values, with at least one positive count;
- exact observation/feature alignment for the selected source; Raw may legitimately carry a different feature axis from current `var_names`, and that Raw axis becomes the aggregate axis without positional truncation;
- finite, non-negative, integer-like values. Software cannot prove from values alone that the analyst did not select normalized, imputed, or otherwise inappropriate integer data, so the source declaration and this limitation are reported;
- non-missing, non-blank Sample and population labels without string-coercion collisions;
- one non-missing Condition value per Sample;
- when supplied, one non-missing Technical batch value per Sample.

Condition and Technical batch are metadata, not aggregation keys. If cells assigned to one Sample disagree on Condition or Technical batch, the dataset does not satisfy the simple independent-sample contract. Fail with the conflicting Sample IDs; do not split that Sample into pseudo-replicates or choose a majority label. A repeated-measures/technical-replicate design would require a separate explicitly modeled node.

After aggregation, require exactly one row per represented Sample-by-population pair, unique profile identities, finite non-negative integer sums, no zero-cell or zero-library retained profile, the selected-source gene order, and metadata values that are constant and correctly aligned. `min_cells`/`min_counts` filtering may remove a complete stratum, but that loss is recorded and warned; a downstream request for an absent Condition remains a hard design error.

## Statistical unit, downstream design, and degenerate inputs

Aggregation itself fits no model and constructs no design matrix. It creates the count substrate on which a downstream engine selects exactly one population and constructs a Sample-level design such as:

```text
~ TechnicalBatch + categorical_covariates + continuous_covariates + Condition
```

The aggregation report must therefore disclose the number of unique Samples, profiles and contributing cells per population and Condition, not present cell counts as replication. Profiles missing because a Sample had no cells in a population are absent observations, not zero-count substitute replicates. Do not synthesize empty Sample-by-population rows.

Reject computation-breaking degeneracies and report interpretation-limiting ones separately:

- no surviving profiles or genes;
- only one Sample overall or sparse Condition/population strata are report warnings unless they leave a downstream requested design non-estimable;
- a Sample appearing under multiple Condition or Technical batch values;
- all-zero aggregated libraries, overflow, non-finite sums, or loss/reordering of genes;
- population labels that are provisional/unknown when formal mode is requested produce `formal_interpretation_invalid` plus a warning, not an aggregation failure;
- filters that leave few independent Samples are disclosed. Three or more per Condition is a practice recommendation, not an aggregation or universal engine-validity boundary.

## Required outputs and disclosure

The scientific output should be a typed pseudobulk-count artifact, not a generic AnnData that can be confused with cell-level data. Its report-ready `summary` JSON must include:

- `methods`: analyst-declared count source/feature axis, sum aggregation, Sample and population key declarations, profile QC policy, and statement that Sample is the declared inferential unit;
- `results`: input cell/gene/Sample counts, output profile count, profiles and contributing cells/counts by population and Condition, retained/removed profile records, and library-size/cell-count distributions;
- `key_results`: concise prose suitable for methods/report writing without claiming differential expression;
- `parameters`: visible thresholds and every resolved fixed policy, including `mode="sum"`, no synthetic empty profiles, annotation status, Condition key, and optional Technical batch key;
- `warnings`: source/completeness/Sample-identity declarations, role aliases, sparse strata, low replication, provisional/unknown labels, removed profiles, and dataset-specific threshold caveat;
- `references`: the Decoupler API/source, pseudobulk practice papers, AnnData, and any package actually used;
- `software_versions`: openbio-singlecell, Decoupler, AnnData, NumPy, pandas, SciPy, and Scanpy only if executed.

The `code` output must be executable equivalent function source. It must select and validate the same count matrix, call the same public aggregation API with explicit arguments, normalize unstable upstream QC names into the stable artifact contract, enforce Sample/Condition/Technical-batch integrity, apply the same profile filters, and return the same primary artifact. It may omit only plugin-history/timing bookkeeping; it must not omit scientific validation, hide defaults, or depend on ComfyUI internals.

AnnData references:

- AnnData aligned storage contract: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371
