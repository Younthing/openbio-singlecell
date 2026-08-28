# OpenBioSingleCellMarkerExpressionPlot — design review

Research date: 2026-08-28

## Recommendation

**Keep and enhance. Do not merge with Marker Genes and do not delete.**

Marker Genes performs statistical ranking and returns a stable table. Marker Expression Plot renders a user-declared panel from a selected expression representation. Combining them would couple exploratory feature selection to visualization and would blur the repository ADR separating **Cluster marker evidence** from **Condition contrast**. The plot node is a useful atomic analysis when it only renders and reports a panel; it must not rank genes, change annotations, or imply differential evidence.

The four render modes can remain internal adapters behind one module because they share genuine validation, expression-source, grouping, reporting, and PNG-lifecycle behavior. The external interface should use a dynamic plot-mode choice so callers do not see parameters that a selected adapter ignores.

## Current interface and behavior

Current inputs:

| Input | Current behavior | Finding |
| --- | --- | --- |
| `adata` | Passed directly to Scanpy plotters | Usually read-only, but implicit dendrogram computation can mutate `uns` |
| `genes` | Comma-separated, trimmed, deduplicated | Good compact interface; missing strict string/type/size validation |
| `groupby` | Required `obs` key | Numeric values can be silently binned by some Scanpy plotters; missing/blank labels are not audited |
| `plot_type` | dotplot/matrixplot/tracksplot/violin | Direct Python calls with any other string silently fall through to violin |
| `source` | X/raw/layer | Correct seam; raw feature count is currently reported incorrectly as `adata.n_vars` |
| `standard_scale` | Always visible, passed only to dot/matrix | Shallow interface: ignored by tracks/violin but still reported as if used |
| `dendrogram` | Always visible, passed to dot/matrix/tracks | Ignored by violin; can mutate caller and may use the wrong representation |
| `log` | Shared across all modes | Ambiguous/inconsistent scientific behavior across plot families |

Current output is only `plot: PlotResult`. Its `parameters`, prose `description`, warnings, elapsed time, and input dimensions are useful transport metadata, but it does not satisfy the required JSON scientific summary or equivalent source-code output.

Known consumers:

- `OpenBioSingleCellPreviewResult`
- `OpenBioSingleCellSavePNG`
- the Single-Cell Best Practice workflow, using dotplot on `layer:log1p_norm`, `groupby=leiden_scvi`, no dendrogram

There are no direct behavior tests for Marker Expression Plot. The example test only checks two saved widget values.

## Correctness findings

### Critical: hidden mutation

With dotplot, matrixplot, or tracksplot and `dendrogram=True`, Scanpy computes and stores `uns['dendrogram_<groupby>']` when it is missing. The node passes caller-owned `adata`, so its read-only description is incorrect. This violates both locality and the plot node's atomic contract.

Required design: plot only from a private copy. Compute any requested dendrogram explicitly from the selected marker genes and selected expression representation under a private key. Assert that the caller's `X`, layers, raw, obs, obsm, and uns are unchanged in tests.

### Critical: dendrogram can represent different data from the panel

Implicit Scanpy defaults can use PCA or `.X`, not the named layer/raw marker panel. The tree may therefore reorder groups using unrelated dimensions. The node should compute a group-mean dendrogram from exactly the selected genes and source, or disable dendrogram when that cannot be guaranteed.

### High: backend depends on external global state

The implementation assumes legacy Matplotlib return types. Scanpy's preview preset selects HoloViews adapters. Use a local `sc.settings.override(preset=sc.Preset.ScanpyV1, autoshow=False)` context and validate the public return object. Never change the process-global preset permanently.

### High: parameters are accepted or reported when unused

`standard_scale` is irrelevant to tracksplot and violin; dendrogram is irrelevant to violin. The current summary metadata still lists them. A dynamic plot configuration should expose only meaningful choices and the report must contain only effective parameters.

### High: `log` mixes transformation and presentation

Remove the common `log` parameter. Expression transformation belongs upstream in the explicit source layer. If needed, a violin-only `y_scale=linear|log` is a presentation choice. Generated code and the summary must name the actual operation, never merely `log=True`.

### Medium: source and grouping validation is incomplete

Add validation for:

- in-memory, nonempty AnnData and unique observation identifiers;
- selected source shape, numeric type, and finite values for selected genes;
- unique selected-source feature names;
- complete nonblank group labels;
- explicit categorical grouping with preserved category order;
- at least one observed category and counts per category;
- strict `plot_type`, scale, and boolean types for Python callers;
- bounded gene count or a resource estimate before materializing a large dense panel.

One group can remain valid but should warn that cross-group comparison is unavailable. A numeric `groupby` should be rejected unless the node deliberately exposes a separate binning policy; hidden Scanpy binning is not acceptable.

### Medium: figure lifecycle is implicit

Render to PNG in a `try/finally` and close/release the Matplotlib figure after bytes are captured. Multiple workflow executions must not accumulate pyplot-managed figures.

### Medium: scientific claims are too thin

The existing description only says that Scanpy rendered the image. It should explain what color, dot size, violin width, tracks, scaling, and dendrogram mean for the selected mode. It must state that the result is descriptive **Cluster marker evidence**, not a replicate-aware **Condition contrast** or cell-type truth.

## Proposed atomic interface

```text
inputs
  adata
  genes                       comma-separated exact feature identifiers
  groupby                     categorical obs column
  plot                        DynamicCombo
    dotplot
      standard_scale          none | var | group; advanced
      expression_cutoff       default 0; advanced
      mean_only_expressed     default false; advanced
    matrixplot
      standard_scale          none | var | group; advanced
    tracksplot
      # no standard_scale
    violin
      density_norm            width | area | count; advanced
      show_cells              default false; advanced
      y_scale                 linear | log; advanced
  source                      X | raw | layer
  group_order                 observed | dendrogram; advanced
  figure_width                automatic by default; advanced only if needed
  figure_height               automatic by default; advanced only if needed

outputs
  plot                         PlotResult
  summary                      SummaryResult (JSON-compatible payload)
  code                         STRING
```

`group_order=dendrogram` can internally use fixed, disclosed Pearson correlation plus complete linkage to keep the interface narrow. If researchers need alternatives in real workflows, add correlation/linkage as nested advanced inputs under the dendrogram option; do not expose them for observed order.

Hidden implementation choices that should be fixed and reported:

- legacy Matplotlib Scanpy preset in a local override;
- `show=False` and no Scanpy save side effect;
- dot expression cutoff `0.0` and mean over all cells unless explicitly exposed;
- PNG DPI and bounding-box policy;
- a collision-safe private dendrogram key;
- figure cleanup;
- stable group/category ordering rules.

Do not expose arbitrary `**kwargs`, Scanpy save paths, global theme settings, or Matplotlib axes objects. Those widen the interface without serving the workflow contract.

## Summary contract

The `summary` payload should be JSON serializable with finite values only and should use the repository-wide report envelope. At minimum:

```json
{
  "methods": [
    "Rendered a Scanpy dot plot from layer:log1p_norm grouped by leiden_scvi; color is group mean and size is the fraction above 0."
  ],
  "results": [
    "Displayed 12 genes across 8 groups and 4231 cells as descriptive Cluster marker evidence."
  ],
  "key_results": {
    "plot_type": "dotplot",
    "expression_source": "layer:log1p_norm",
    "expression_state": "logged",
    "genes": ["ACTA2", "TAGLN"],
    "groupby": "leiden_scvi",
    "group_order": ["0", "1"],
    "group_cell_counts": {"0": 123, "1": 98},
    "missing_group_cells": 0,
    "display_semantics": {
      "color": "mean expression in group",
      "size": "fraction of cells with expression > 0",
      "standard_scale": "var"
    },
    "dendrogram": {
      "enabled": false,
      "source_matches_panel": true,
      "correlation": null,
      "linkage": null
    },
    "plotted_statistics": [
      {"group": "0", "gene": "ACTA2", "n_cells": 123, "mean": 1.23, "median": 0.8, "fraction_above_cutoff": 0.72}
    ]
  },
  "parameters": {},
  "references": [],
  "software_versions": {},
  "warnings": [],
  "limitations": []
}
```

`plotted_statistics` should match the plot's unscaled descriptive inputs. For dotplot it must include mean and expressing fraction; for matrixplot mean; for violin useful quantiles; for tracksplot per-group/gene descriptive ranges. If size becomes excessive, impose and document a gene-panel limit rather than silently truncating the scientific summary.

Required warnings/limitations include, as applicable:

- expression provenance unknown or source appears to be raw/untransformed when the chosen plot normally expects normalized values;
- one group, small groups, missing labels, or omitted cells;
- standard scaling destroys absolute magnitude comparisons along the scaled dimension;
- dendrogram order is descriptive and derived from group means of the selected panel;
- the visualization is Cluster marker evidence, not Condition inference or a Curated annotation.

References should include the Scanpy paper, the exact public plot function documentation used, `tl.dendrogram` when enabled, and the Matplotlib paper. Software versions should include Scanpy, Matplotlib, Seaborn, anndata, NumPy, pandas, and SciPy actually used in the chosen path.

## Code contract

`code` must be compilable, self-contained Python containing one callable such as:

```python
def plot_marker_expression(adata):
    ...
    return figure_or_png
```

It must reproduce the scientific and rendering state:

- exact gene order and group order;
- exact expression source selection;
- validation of source features and grouping;
- plot-specific effective parameters only;
- explicit private-copy dendrogram computation from the same panel/source;
- local ScanpyV1/autoshow override;
- fixed figure size/DPI and no input mutation;
- explicit figure cleanup if returning PNG bytes.

The code should not import OpenBio helpers and should not depend on hidden workflow history to run. It can omit OpenBio report/history construction, but it cannot silently fall back to another source, another group order, another backend, or Scanpy defaults that differ from runtime.

## Failure contract and tests

### Validation failures before plotting

- backed or empty AnnData where safe rendering is unsupported;
- nonunique observation or selected-source feature identifiers;
- empty/non-string gene specification;
- selected source missing or misaligned;
- missing, duplicate, nonnumeric, or nonfinite selected gene values;
- missing groupby, blank/missing groups, numeric groupby without explicit binning;
- invalid plot variant or variant-specific option;
- dendrogram requested with fewer than three groups: warn and retain observed order, or reject consistently;
- dendrogram cannot be computed because group means are constant/nonfinite;
- estimated plot/report size exceeds a documented resource limit.

### Backend failures

- installed Scanpy lacks required public parameters/return contract;
- expected Matplotlib figure cannot be obtained;
- dendrogram result is incompatible with current categories;
- PNG is empty or lacks the PNG signature.

### Required tests

1. Schema pins `plot`, `summary`, `code` outputs and dynamic option visibility.
2. Dot/matrix/tracks/violin fake-adapter tests pin effective public kwargs.
3. Input is deep-unchanged for every mode, especially dendrogram cache miss.
4. Dendrogram uses the selected genes and source, not caller PCA/X, and reports the resulting order.
5. Existing caller dendrogram metadata is neither overwritten nor trusted implicitly.
6. Global Scanpy preset/autoshow state is restored after success and failure.
7. X/raw/layer feature-axis and `input_genes` semantics are correct.
8. Numeric/missing/blank groups and nonfinite expression fail clearly.
9. Plot-specific summaries contain the exact descriptive statistics used.
10. Summary is strict JSON; code compiles and reproduces the image semantics.
11. No figures remain registered after repeated success or backend failure.
12. A real Scanpy smoke test is optional/skipped when dependencies are unavailable.

## Migration

- Append outputs in order `plot, summary, code` so existing first-output Preview/SavePNG links remain meaningful; migrate serialized node output metadata.
- Replace `plot_type`, global `standard_scale`, `dendrogram`, and `log` widgets with dynamic plot/group-order options. Named widgets require explicit mapping; positional legacy workflows require an ordered migration.
- Map existing `log=False` to no expression transformation. For `log=True`, do not guess equivalence across modes: migrate to the closest mode-specific option with a visible warning, or require review.
- Existing best-practice dotplot maps directly to dotplot + `standard_scale=var` + observed order and `layer:log1p_norm`.
- Preserve node ID to avoid needless graph replacement.

## Decision against merging

- Do not merge into Marker Genes: ranking and plotting are distinct atomic analyses with different outputs and claims.
- Do not merge into UMAP Plot: the only common behavior is PNG/report plumbing, which belongs in private reporting/render helpers rather than a broad public plotting interface.
- If another marker-panel renderer later needs identical source/group/statistics behavior, extract a private deep module then; there is already enough variation among the four current Scanpy adapters to justify an internal seam, not four duplicated public nodes.

## Final-review resolution

The frozen public input order is:

```text
adata, genes, groupby, plot, source, group_order, random_seed
```

`random_seed` is an advanced integer, default `0`, range `0..2^31-1`. It is a rendering seed for violin individual-cell
jitter, not an inferential seed. The implementation applies it in an isolated NumPy RNG context for every branch so
the interface and generated function have one reproducibility rule; the complete caller RNG state is restored after
success and backend failure. Runtime and generated violin PNG bytes must agree for the same environment/input/seed.

The private module no longer copies caller-owned AnnData. It validates and slices the explicit X/Raw/layer source,
budgets the complete estimated plotting working set, and constructs a minimal AnnData containing only selected panel
X, one groupby obs column, and selected gene identifiers. Scanpy receives `use_raw=False, layer=None` because the
adapter's X is already the exact bound representation. This is not a source fallback: the report and generated code
retain the original source identity, and dense, sparse, and Raw tests pin the selected values. Dendrogram input is a
separate minimal copy of this same panel and never sees caller PCA, unrelated X, layers, Raw, obsm, or cached uns.

Raw state reporting is evidence-based: an OpenBio Raw snapshot with matching history is `counts`; an external Raw
panel equal to current X with proven logged/normalized/etc. state inherits that evidence; every other Raw source is
`unknown` with an explicit no-assumption warning. No branch silently transforms the expression representation.
SciPy is always versioned/referenced because sparse-source handling is common. Seaborn is versioned/referenced only
for violin, where the supported Scanpy adapter actually invokes it.
