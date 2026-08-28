# OpenBioSingleCellCellTypeCorrelation — design review

## Decision

**Enhance and rename the display to “Population Centroid Correlation”.** Retain one cohesive descriptive analysis
that emits both machine-readable correlations and a plot. Do not split computation and rendering until a second
renderer exists; the canonical table/matrix is the test surface and the plot is a deterministic view of it.

## Atomic interface

Inputs, in order:

1. `adata`;
2. `population_key`, default `cell_type`;
3. required explicit `representation_key`, default `X_pca`;
4. `correlation_method` (`pearson`, `spearman`, `kendall`);
5. advanced `linkage_method`, default `complete`;
6. advanced `n_dimensions`, where `0` means all validated representation dimensions;
7. `annotation_status` (`unknown`, `provisional`, `curated`);
8. visual-only `color_map` and `show_numbers`;
9. maximum groups/table rows.

Outputs:

```text
table, plot, summary, code
```

Require a named `obsm` matrix; no blank/automatic X/PCA selection. The table contains one deterministic upper-triangle
row per distinct population pair: population A/B, correlation, dendrogram order indices, cell counts A/B, method,
representation, and dimension count. The plot renders the full symmetric matrix in dendrogram order. The input is
read-only and no dendrogram state is written back.

## Validation/report contract

Validate unique observation IDs, nonblank scalar labels, categorical observed-level order, finite 2-D representation,
dimension bounds, at least two groups and two informative dimensions, finite symmetric unit-diagonal correlation,
valid linkage/order, and plotting palette. Exclude unused categorical levels and disclose them. Fail undefined
constant-centroid correlations.

`summary` follows the strict report schema with methods, centroid/group support, representation provenance, bounded
strongest positive/negative pairs and ties, category/dendrogram order, warnings/limitations, references, and dynamic
OpenBio/Scanpy/AnnData/NumPy/pandas/SciPy/matplotlib versions. `code` defines an equivalent function returning
`(pair_dataframe, figure_or_png, summary_dict)` and uses the same canonical computation rather than reading `uns`.

## Sample and Technical batch semantics

No p-value or Condition claim is produced. Sample is not treated as a replicate; unequal cells per Sample can weight
centroids, and the summary reports cell counts but cannot infer Sample support without changing the method. Technical
batch is not adjusted here; representation provenance must say whether/how upstream integration modeled it. The
result is exploratory similarity, not lineage or cell-type truth.

## Migration and tests

Legacy ports add `table`, `summary`, and `code`; `plot` remains available. Blank/implicit representation values need
an explicit user choice. Rename `groupby`/`use_rep` to the clearer population/representation terms and preserve known
links where exact.

Tests: hand-computed centroid/correlation matrices for all methods; Scanpy 1.12.3 result parity; group order and unused
categories; dimension subsetting; non-finite/constant/one-dimension/one-group failures; representation producer and
Technical-batch disclosure; Provisional/Curated status; symmetric table/plot consistency; no AnnData mutation or
global matplotlib state; palette/number rendering; row guard; strict JSON; code compilation; and runtime/generated
table+summary equivalence.

## Cohesion assessment

The module owns one descriptive centroid-correlation result and its direct view. Representation construction and
statistical inference remain outside; Scanpy result normalization, validation, reporting, and rendering stay local.

## 2026-08-29 repair record

- New category: `openbio/single-cell/diagnostics`.
- Classification rationale: the atomic analysis computes descriptive population-centroid similarities and emits a canonical table plus its deterministic view; returning a plot does not make the underlying operation merely `visualization`.
- Merge/delete decision: retain the cohesive analysis node. Do not merge it with representation construction or inferential models, and do not split or delete it until a second renderer or independent table consumer creates a real seam.
