# OpenBioSingleCellSampleCompositionSummary — design review

Research date: 2026-08-28

## Decision

**Keep and enhance this node as a pure descriptive Sample-composition module. Do not merge it with differential-composition inference, a compositional model, plotting, filtering, or annotation curation.**

Its atomic analysis is one validated aggregation from cell-level metadata to a complete Sample-by-annotation count/proportion grid. It performs no hypothesis test. This separation preserves the **Sample** as the potential future inference unit while ensuring the descriptive result remains useful for quality review, plotting, export, and model preparation.

The deletion test favors retention: without this module, every plot and inferential adapter would duplicate Sample-to-Condition validation, zero-grid construction, denominator semantics, ordering, and composition diagnostics. Those rules provide real Leverage and Locality behind a small Interface.

## Current Interface problems

- group_key and output column group conflict with the domain term **Condition**.
- Missing Sample, Condition, or annotation values are silently excluded, changing denominators and possibly removing Samples.
- Conversion to strings can collapse typed identifiers.
- The denominator is implicit because sample_total_cells is not in the table.
- There is no explicit annotation-status disclosure.
- No output-size preflight protects the Sample-by-annotation product.
- The result lacks report-ready methods/results, references, software versions, limitations, strict JSON, and equivalent source code.
- Tests exercise one four-Sample demonstration only and do not establish the scientific contract.

## Proposed atomic Interface

Visible inputs:

- adata;
- sample_key, default "sample";
- condition_key, default "condition";
- annotation_key, default "cell_type";
- annotation_status, one of unknown, provisional, or curated, default unknown.

Advanced input:

- max_output_rows, a positive integer guard checked against n_samples × n_annotations before constructing the complete grid.

Outputs, in this order:

- table;
- summary;
- code.

annotation_status affects interpretation and report wording, not counts. A curated value is a caller assertion unless exact upstream provenance can verify it. Keeping that assertion visible is preferable to inferring curation from a column name.

Fixed hidden policy:

- every input cell belongs to the denominator; missing or blank required metadata fails;
- the annotation universe is all globally observed nonmissing categories;
- every Sample is crossed with that universe and missing combinations receive count zero;
- globally unused declared categorical levels are excluded from the grid but disclosed;
- no pseudocount, smoothing, imputation, sample filtering, minimum-cell threshold, or condition averaging;
- Sample and Condition order follows first appearance; categorical annotation order is respected after restricting to observed levels, otherwise first appearance is used;
- input AnnData is read-only;
- report previews are bounded deterministically, while the primary table is complete.

Do not expose pandas observed/dropna switches, denominator choice, arbitrary aggregation functions, pseudocounts, or a missing-value policy. Those switches would turn one clear descriptive operation into several incompatible meanings.

## Primary table contract

The table has one row per complete Sample-annotation pair and these stable columns:

| Column | Meaning |
| --- | --- |
| sample | Canonical Sample label |
| condition | The unique Condition assigned to the Sample |
| annotation | Population annotation label |
| cell_count | Nonnegative integer count of retained input cells |
| sample_total_cells | Positive integer denominator shared by every row for that Sample |
| proportion | cell_count divided by sample_total_cells |

Required invariants:

- row count equals n_samples × n_annotations;
- each (sample, annotation) key is unique;
- each Sample maps to one Condition;
- cell_count is integer and nonnegative;
- sample_total_cells is constant within Sample and positive;
- counts summed across the whole table equal adata.n_obs;
- counts summed within Sample equal sample_total_cells;
- proportions are finite and in [0, 1];
- proportions sum to one within each Sample to a documented floating tolerance;
- cell_count zero implies proportion zero;
- the input obs, axes, expression matrices, and unstructured metadata are unchanged.

One Sample and one annotation are valid for description. They are not evidence for a Condition contrast; the summary must say so. Zero genes are also acceptable because this module reads only observation metadata.

## Deep Module and Seam placement

Create a private composition Implementation that owns:

- direct-call parameter and identifier validation;
- collision-safe label normalization;
- Sample-to-Condition mapping validation;
- deterministic annotation-universe resolution;
- complete zero-grid aggregation;
- invariant checks and diagnostics;
- bounded summary construction;
- equivalent-code rendering.

Its internal result should be a canonical table plus diagnostics, not a partially processed AnnData. The same zero-grid constructor is a justified Seam for future reviewed compositional-model Adapters because both description and inference require identical Sample/count semantics. Model formula parsing, reference-population choice, posterior inference, p-values, and plotting remain outside this module.

A public generic frequency-table Module would be shallower: it would make every caller reconstruct Sample, Condition, denominator, and closure semantics. The domain-specific Interface earns its depth by hiding those rules.

## Summary contract

summary is a SummaryResult whose payload round-trips through strict JSON with allow_nan=False. It contains at least:

~~~json
{
  "methods": [],
  "results": [],
  "key_results": {
    "input_cells": 0,
    "samples": 0,
    "conditions": 0,
    "annotations": 0,
    "complete_grid_rows": 0,
    "zero_count_rows": 0,
    "zero_count_fraction": 0.0,
    "samples_per_condition": {},
    "cells_per_sample": {
      "min": 0,
      "q1": 0.0,
      "median": 0.0,
      "q3": 0.0,
      "max": 0
    },
    "denominator": "all input cells within each Sample",
    "maximum_proportion_sum_error": 0.0,
    "annotation_status": "unknown",
    "condition_descriptives": []
  },
  "parameters": {},
  "references": [],
  "software_versions": {},
  "warnings": [],
  "limitations": []
}
~~~

methods states that cell counts were aggregated by Sample over the complete observed annotation universe and normalized within Sample. results provides report-ready prose naming the number of independent Samples, Condition balance, cell-depth range, leading annotation proportions, and zero-grid prevalence. condition_descriptives reports Sample-level n, mean, median, quartiles, min, max, and zero-Sample count for each Condition-annotation pair; if bounded for JSON size, the truncation count and deterministic ordering are explicit.

references includes the AnnData and pandas software references plus Aitchison and the single-cell composition references used to justify interpretation. software_versions dynamically records Python, openbio-singlecell, anndata, pandas, and NumPy when used. It must not list SciPy, statsmodels, pertpy, or Scanpy merely because they exist in the environment.

Mandatory limitations:

- proportions are closed relative captured-cell composition, not absolute abundance;
- dissociation, capture, QC, and sampling depth can affect composition;
- a zero captured count does not prove biological absence;
- no statistical Condition contrast was performed;
- Provisional or unknown annotation status does not become Curated annotation.

## Code contract

code is compilable, self-contained Python defining one function such as summarize_sample_composition(adata). It returns the same six-column pandas DataFrame, in the same row order and dtypes, and reproduces:

- distinct/present/nonblank column validation;
- missing-value and typed-label-collision failures;
- unique cell identifiers;
- one-Condition-per-Sample validation;
- annotation-universe and ordering policy;
- max_output_rows preflight before grid materialization;
- complete zero counts and denominator/proportion calculation;
- every table invariant and input immutability.

The code may omit plugin timing and Result wrapping. It may not silently drop invalid cells, rely on a pandas default that omits zero combinations, or return only observed nonzero pairs. The literal parameters and annotation-status assertion must be embedded.

## Failure and boundary contract

Reject:

- non-AnnData-like input or missing obs;
- zero observations;
- blank, duplicated, or absent column names;
- missing, blank, non-scalar, or display-colliding labels;
- duplicated observation identifiers;
- a Sample mapped to zero or multiple Conditions;
- no observed annotation;
- boolean/noninteger/nonpositive max_output_rows;
- n_samples × n_annotations above the row budget;
- any post-aggregation count, denominator, range, key, or closure invariant failure.

Do not reject backed AnnData solely because it is backed: this operation reads obs and does not require the expression matrix. Do not require genes or a minimum number of Samples because the module is descriptive.

## Required tests

1. Schema pins input order/defaults and outputs table, summary, code.
2. A category absent from one or several Samples receives explicit zero rows.
3. Row count, total-cell conservation, denominators, and closure hold with unequal cells per Sample.
4. Missing Sample, Condition, and annotation values fail without dropping cells.
5. Whitespace labels and typed string collisions fail.
6. A Sample assigned to two Conditions fails even if one conflicting cell has missing annotation.
7. Categorical order is preserved; globally unused levels are excluded and reported.
8. One Sample/one annotation is valid but reported as descriptive only.
9. The grid budget fails before MultiIndex/table materialization.
10. Caller obs, X, layers, raw, axes, and uns remain unchanged on success and failure.
11. Summary contains only finite JSON values and actual dependency versions/references.
12. Generated code compiles and reproduces dense table values, dtypes, ordering, zeros, and representative failures.

## Migration

- Rename group_key to condition_key and the table column group to condition in the node, generator, examples, tests, and downstream consumers together.
- Append summary and code after table so legacy output-zero links remain meaningful.
- Migrate the Sample Composition Comparison workflow to the descriptive node only unless a separately audited formal compositional model is added.
- Mark its annotation_status explicitly. The end-to-end workflow currently summarizes CellTypist output and therefore must say provisional, not curated.
- Update CSV names and README language to "Sample composition summary"; do not imply that the descriptor compares Conditions statistically.
- Do not keep simultaneous public aliases group/condition indefinitely; a graph migration should rewrite the old widget name once.
