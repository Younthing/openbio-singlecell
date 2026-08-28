# OpenBioSingleCellFilterMarkerGenes — official usage research

## Scope and current environment

This record covers deterministic threshold selection from an already computed canonical **Cluster marker evidence** table. It does not recompute statistics, modify AnnData, validate a cell type, or perform a Condition contrast.

Environment inspected on 2026-08-28:

- `scanpy 1.12.3`;
- `pandas 3.0.5`;
- `numpy 2.4.4`.

The current node takes any `TableResult` with four named columns, coerces those columns to numeric, retains rows satisfying inclusive fold-change/prevalence/adjusted-p thresholds, and returns one filtered `TableResult`. It does not verify marker-result provenance, group/gene keys, finite values, or duplicates, and has no `summary` or `code` output.

## Primary sources

- Scanpy `filter_rank_genes_groups`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.filter_rank_genes_groups.html
- Scanpy 1.12.3 filter implementation: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_rank_genes_groups.py
- Scanpy `rank_genes_groups`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
- Scanpy `rank_genes_groups_df`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.rank_genes_groups_df.html
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Repository ADR: `docs/adr/0001-separate-marker-evidence-from-condition-inference.md`.

## Official Scanpy filter contract

The installed signature is:

```python
scanpy.tl.filter_rank_genes_groups(
    adata,
    *,
    key=None,
    groupby=None,
    use_raw=None,
    key_added="rank_genes_groups_filtered",
    min_in_group_fraction=0.25,
    min_fold_change=1,
    max_out_group_fraction=0.5,
    compare_abs=False,
)
```

Scanpy filters an AnnData `rank_genes_groups` bundle using log fold change, fraction expressed within the group, and fraction expressed outside the group. It stores a second structured result in `adata.uns[key_added]`; rejected gene names are replaced by NaN so the original structured shape is retained. It returns `None` and has no adjusted-p-value threshold.

The 1.12.3 implementation uses strict comparisons:

```python
(fraction_in > min_in_group_fraction)
& (fraction_out < max_out_group_fraction)
& (fold_change > min_fold_change)
```

`compare_abs=True` applies the fold-change threshold to absolute magnitude; that is unsuitable for this node's positive cluster-characterizing marker contract and should remain fixed false. When compatible `logfoldchanges`, `pts`, and `pts_rest` are absent, Scanpy may recompute values from the selected AnnData expression. The OpenBio table filter has no AnnData and must never pretend to recompute missing statistics.

The current OpenBio node is therefore not a wrapper around `scanpy.tl.filter_rank_genes_groups`. It is a useful table-native selection with one additional adjusted-p-value criterion and inclusive threshold semantics. Its report and generated code must state that distinction. It should not cite Scanpy defaults as if they were universal biological cutoffs.

## Proposed canonical table and threshold semantics

Require the canonical MarkerGenes columns:

- identifiers: `group`, `gene`, `rank`;
- statistics: `score`, `log2_fold_change_approx`, `p_value`, `p_adjusted`;
- prevalence: `fraction_in_group`, `fraction_reference`.

Require provenance showing that the input operation is **directly** `marker_genes` and uses a non-logistic method. Reject `filter_marker_genes` input so repeated filtering cannot silently turn several analyst choices into one indistinguishable selection. Reject missing, non-numeric, non-finite, or out-of-range promised statistics instead of coercing them to NaN and silently dropping rows. Require unique `(group, gene)` pairs, positive ranks, fractions and p-values in `[0,1]`, tolerance-aware `p_adjusted >= p_value`, whitespace-unambiguous identifiers, and consistent group labels. Preserve the source table's within-group rank ordering.

The input marker table and its tested-gene universe must carry the same analysis, complete-ranking, and ordered-universe fingerprints. The consumer recomputes each artifact's current-content fingerprint before thresholding. This detects numeric mutation and rejects a universe from a same-axis/same-parameter run whose underlying expression produced a different complete ranking.

The proposed filter remains inclusive, matching current OpenBio behavior:

```python
log2_fold_change_approx >= min_log2_fold_change
fraction_in_group >= min_fraction_in_group
fraction_reference <= max_fraction_reference
p_adjusted <= max_p_adjusted
```

Inclusive versus strict equality must be explicit in methods/code because it differs from Scanpy's helper. All four thresholds must be finite; fraction and adjusted-p thresholds must lie in `[0,1]`. A fold-change threshold can be any finite number, although a positive default expresses the positive-marker contract.

`p_adjusted` is the upstream per-group BH result. Filtering it at 0.05 does not create a new correction, correct for trying several thresholds, or turn pooled-cell evidence into Sample-level inference. Thresholds are analyst choices, not universal marker definitions. Prevalence cutoffs are especially sensitive to sequencing depth, normalization, dropout, and group size.

## Diagnostics and output contract

Primary outputs are a filtered `TableResult` with exactly the canonical input columns in their original order and the original tested-gene-universe `TableResult` unchanged. The filtered table preserves the upstream complete-ranking fingerprint, records the direct MarkerGenes table content fingerprint, and carries its own recomputable current-content fingerprint. The inputs remain unchanged.

`summary` should report:

- upstream marker method, group/reference policy, expression source, and original ranking parameters available from provenance;
- all four thresholds and inclusive comparison operators;
- input and retained row counts and fractions, globally and per group;
- rows failing each individual criterion, overlap-aware total rejection, groups with zero retained markers, and minimum/median/maximum retained counts;
- numeric ranges/missing counts checked for each criterion;
- explicit warnings that the thresholds are user choices and that marker evidence is exploratory;
- Scanpy filtering semantics as related documentation, BH, Scanpy, and the replicate-awareness caveat;
- Python, openbio-singlecell, pandas, NumPy, Scanpy, and statsmodels versions. Scanpy and statsmodels are listed because their upstream ranking/correction provenance and related semantics are disclosed even though the filter mask itself is pandas logic.

`code` should be a self-contained function accepting `(marker_table, universe)` pandas DataFrames, validating both exact schemas and gene membership, checking the input content identities embedded for this source operation, applying the same inclusive mask, preserving row order, and returning `(filtered_table, universe_copy)`. It should not require AnnData, ComfyUI, OpenBio history, or Scanpy.

## Scientific caveats

- Filtering selects from the finite ranking exported upstream. If MarkerGenes was truncated, a gene absent from the table cannot be recovered even if it would satisfy thresholds; summary must disclose upstream truncation.
- Thresholding after inspecting results is exploratory and can introduce selection bias. Repeated threshold tuning is not represented by the upstream adjusted p-values.
- A gene passing every threshold is stronger descriptive Cluster marker evidence, not a validated Curated annotation and not proof of exclusivity.
- Cell-level adjusted p-values remain unsuitable for formal Condition inference because cells within a Sample are not independent replicates. Condition contrasts require Sample-level methods.
- Empty output is a valid possible result only when the validated input genuinely has no passing rows; it must be reported prominently, not caused silently by missing logistic-regression statistics.
- A filtered result is terminal for this node. To revise thresholds, rerun FilterMarkerGenes from the direct MarkerGenes table and universe rather than chaining FilterMarkerGenes outputs.
