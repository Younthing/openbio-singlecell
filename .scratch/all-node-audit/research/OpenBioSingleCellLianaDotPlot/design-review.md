# OpenBioSingleCellLianaDotPlot — design review

## Decision

**Enhance and retain as a read-only Sample-resolved LIANA result renderer.** Consume the typed direct communication
result, derive method semantics from its verified metadata, preselect rows deterministically, call
`dotplot_by_sample`, and output `plot`, `summary`, and `code`. Keep it separate from computation so visual selection
cannot alter the scientific result.

## Atomic interface

Inputs, in order:

1. typed `result` from `OpenBioSingleCellLianaCommunication`;
2. optional exact sender `source_labels` and receiver `target_labels`;
3. method-specific `selection` dynamic choice:
   - rank aggregate: `max_specificity_rank` and ascending `magnitude_rank`; or
   - CellPhoneDB: `max_cellphone_pvalue` and descending `lr_means`;
4. `top_n`, with an explicit fixed scope such as per Sample/sender/receiver pair;
5. `figure_width` and `figure_height`;
6. advanced `max_plot_rows` and image-size guard.

Outputs, in order:

```text
plot, summary, code
```

There is no `adata`, `result_key`, or separate `method` input. The typed result owns the method, Sample key, native
score fields/directions, resource provenance, and Condition mapping. The node rejects a schema/metadata mismatch
rather than guessing. A generic `significance_threshold` is removed.

The implementation applies sender/receiver/cutoff/top-N selection to a defensive table copy with deterministic ties,
then calls `li.pl.dotplot_by_sample(liana_res=selected, ...)`. Colour/size inversions are fixed per method only for
visual readability and are explicitly reported. The source table and metadata remain byte-for-byte unchanged.

## Invariants and reporting

- Required identity and method-native score columns exist, are finite/in-range, and match typed metadata. Sample,
  source, target, ligand-complex, and receptor-complex identities are nonblank and unique at the declared result
  grain; CellPhoneDB's additional subunit fields are validated without fabricating them for rank aggregation.
- Requested labels must exist; selection directions and threshold meanings are method-specific; top-N ties have a
  stable secondary key. An empty selection returns a clear schema-stable warning or controlled error, never a
  misleading blank figure.
- The selection plan is computed before rendering and independently checked against plotted rows. Figure dimensions,
  point count, facet count, and pixel/memory bounds are enforced.
- `summary` follows the strict report schema, records the full visual selection audit trail, bounded plotted key
  interactions, Sample/Condition/resource provenance, warnings/limitations, references, and dynamic OpenBio, LIANA,
  plotnine, Matplotlib, pandas, and NumPy versions. JSON rejects NaN/Infinity.
- `code` defines an equivalent function returning `(png_bytes, summary_dict)`, repeats validation/selection on a
  copy, uses the exact LIANA 1.9.0 `dotplot_by_sample(..., return_fig=True)` call, draws the returned plotnine object,
  and returns `(png_bytes, summary_dict)`.

## Scientific boundary

The plot displays within-Sample candidate interactions. Rank-aggregate specificity is not a p-value;
`cellphone_pvals` is unadjusted cell-label-permutation specificity, not a Condition test. No visual size, colour, or
facet establishes secretion, binding, spatial contact, causality, or differential communication.

## Migration and tests

The exact legacy-to-current matrix is:

| Legacy contract | Current contract | Safe translation |
| --- | --- | --- |
| input 0 `adata` plus input 6 `result_key` | input 0 `OPENBIO_LIANA_RESULT` | Rewire only to output 0 of a freshly rebuilt Communication producer; pooled hidden tables fail closed. |
| input 1 `method` | removed | Validate it against the rebuilt producer, then derive the method exclusively from the typed artifact. A mismatch blocks migration. |
| inputs 2–3 `source_labels`,`target_labels` | inputs 1–2 with the same names | Preserve only schema-valid strings; empty items or duplicate labels that the old parser silently discarded block migration. |
| input 4 generic `significance_threshold` | input 3 method branch | For a proven rank producer map to `rank_aggregate.max_specificity_rank`; for a proven CellPhoneDB producer map to `cellphonedb.max_cellphone_pvalue`. Unknown/mismatched methods block. |
| input 5 `top_n` | input 4 `top_n` | Preserve exact valid positive integers. |
| inputs 7–8 `figure_width`,`figure_height` | inputs 5–6 | Preserve only finite values within the new 1–30 bounds; otherwise require review. |
| absent | inputs 7–8 `max_plot_rows`,`max_image_pixels` | Append current safe defaults. |
| output 0 `plot` | output 0 `plot` | Existing output-zero links remain type-compatible after the producer/result rewrite. |
| absent | outputs 1 `summary`, 2 `code` | Trailing additions; do not fabricate legacy history. |

Migration must preflight and rewrite the Communication→Results→DotPlot chain atomically. Standalone plots, arbitrary
AnnData producers, missing/pooled tables, partial or mixed schemas, method disagreement, and incompatible exposed or
connected widgets reject before mutation. Current-schema nodes are no-ops. Implementation belongs to the shared
workflow-migration owner. That final owner pass is complete: the whole-graph rejection cases, exact current no-ops,
atomicity, and idempotence are covered by the frozen migration suite.

Tests must cover both method profiles and official 1.9.0 call signatures, explicit rejection of 1.10.x, sample facets, exact label filtering,
method-specific cutoffs/directions/inversions, deterministic top-N ties, empty selections, missing/malformed columns,
metadata mismatches, input immutability despite LIANA inversion, figure guards, renderer return types, strict JSON,
scientific wording, dynamic versions, and exact runtime/generated equivalence.

## Cohesion assessment

The module owns one visual transformation. It hides method-specific selection, mutation-prone plotting behavior,
rendering, guards, and report construction behind an explicit typed-result interface while leaving inference and
resource handling upstream.
