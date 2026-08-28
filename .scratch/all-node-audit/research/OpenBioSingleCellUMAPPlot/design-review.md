# OpenBioSingleCellUMAPPlot — design review

Research date: 2026-08-28

## Recommendation

**Keep and enhance. Do not merge with UMAP computation and do not delete.**

UMAP computation mutates a copied AnnData by deriving coordinates from a named graph; UMAP Plot consumes already-computed coordinates and returns a plot. Keeping them separate preserves an atomic compute-versus-render distinction, allows the same embedding to be colored several ways without recomputation, and keeps plot styling out of the scientific embedding interface.

Do not generalize the public node to every embedding yet. There is currently one public plot adapter, so a generic `EmbeddingPlot` seam would be hypothetical. The implementation can use a private reusable scatter renderer. If PCA/t-SNE/force-directed plot nodes later share the same validated interface, a public generalization will then have real adapters and leverage.

## Current interface and behavior

Current inputs:

| Input | Current behavior | Finding |
| --- | --- | --- |
| `adata` | Read-only manual Matplotlib plotting | Good atomic behavior |
| `color` | Empty, numeric obs, or inferred categorical obs | Ambiguous for boolean/integer-coded categories; missing values are mishandled |
| `point_size` | Passed to `Axes.scatter(s=...)` | Needs strict finite positive validation; units are points squared |
| `color_map` | Used only for numeric color | Reported even when categorical/empty; no early registry validation |

The embedding key is hard-coded to `obsm['X_umap']`. This is coupled to only the default output of the upstream UMAP node, even though that node and official Scanpy support a custom `key_added`.

Current output is only `plot: PlotResult`. Consumers are Preview Result and Save PNG. Example workflows use the default UMAP key and color by several observation columns. Existing tests cover categorical/numeric PNG output, shallow copy-on-write checks, missing default coordinates, fewer than two columns, and missing color columns. They do not cover custom keys, row alignment, nonnumeric/nonfinite coordinates, missing color values, palette collisions, parameter types, category order, scientific summary, or equivalent code.

## Correctness findings

### High: upstream key mismatch

Add `embedding_key`, default `X_umap`, and use that exact `obsm` key. For provenance, pair `X_umap` with `uns['umap']`; a custom upstream UMAP key is paired with the same `uns` key under Scanpy's current contract. If matching OpenBio history/metadata is absent, allow imported AnnData but warn that UMAP provenance could not be verified.

### High: coordinate validation is incomplete

Before figure creation require:

- in-memory AnnData with at least one observation;
- unique observation identifiers;
- `embedding_key` present;
- dense/materializable numeric two-dimensional coordinates;
- row count exactly `n_obs`;
- at least two columns;
- finite displayed coordinates.

Fail rather than allowing Matplotlib to silently omit invalid positions. Report both total observations and plotted observations; for valid coordinates they must agree.

### High: color inference and missingness are ambiguous

Use `color_mode=auto|categorical|continuous`:

- auto: pandas categorical, string, object, and boolean are categorical; other numeric values are continuous;
- categorical: preserve declared category order, otherwise deterministic first appearance;
- continuous: strict numeric conversion, reject infinities, and explicitly render/report missing values.

The current `pd.Categorical(values.astype(str))` converts missing values to a visible `"nan"` label and can discard the declared categorical order. It also uses Matplotlib's default cycle, which repeats for more categories than available distinct colors. Replace it with an explicit palette mapping and report that mapping.

### Medium: continuous display policy is implicit

Current continuous points remain in observation order. Dense high-valued points may be hidden. Fix `sort_order=True` to match the documented Scanpy behavior or expose it as an advanced boolean. Validate the colormap using the public Matplotlib registry. Report effective linear `vmin`/`vmax`; warn when the color vector is constant. Robust percentile clipping should only be added as an explicit advanced policy, never silently.

### Medium: categorical legend omission loses meaning

Omitting the legend beyond 20 categories is reasonable for readability, but the output must warn and the summary must retain the category/color mapping and counts. Consider a lower visual threshold if labels are long; the threshold should be fixed and reported rather than hidden in prose only.

### Medium: parameter and figure lifecycle checks

- Trim and validate `embedding_key`, `color`, and colormap strings.
- Require `point_size` to be numeric, non-boolean, finite, and positive.
- Validate `color_mode`, booleans, and any palette name for direct Python calls.
- Fix and report alpha, marker, linewidth, figure size, DPI, and bounding-box policy.
- Render with object-oriented Matplotlib/Agg and release the figure after PNG bytes are captured.

### Medium: scientific reporting is absent

The plot result description says it did not modify AnnData but does not disclose the embedding source, coordinate diagnostics, color statistics/mapping, missing values, UMAP parameters, or interpretation limits. A writer cannot use the output alone to describe what was plotted.

## Proposed atomic interface

```text
inputs
  adata
  embedding_key                default X_umap
  color                        default leiden; empty means uncolored
  color_mode                   auto | categorical | continuous; advanced
  point_size                   default automatic or 10 points^2
  continuous_color_map         default viridis; advanced
  categorical_palette          default tab20/automatic; advanced
  sort_order                   default true; advanced, continuous only
  missing_color                default lightgray; advanced
  legend_policy                automatic | show | hide; advanced, categorical only

outputs
  plot                          PlotResult
  summary                       SummaryResult (JSON-compatible payload)
  code                          STRING
```

To keep the module deep, prefer a DynamicCombo color input:

```text
coloring
  none
  obs
    key
    mode auto|categorical|continuous
```

Then place continuous colormap/sort settings and categorical palette/legend settings only under their applicable mode if the UI supports nested dynamic inputs. Avoid arbitrary Matplotlib kwargs, axes injection, theme mutation, or file paths; Preview/Save PNG already own display/export.

`point_size=automatic` is scientifically safer over wide cell-count ranges. A fixed documented rule such as `min(max_constant / n_obs, upper_bound)` may be hidden and disclosed in the summary, with a fixed-size option for publication tuning. If automatic size would complicate migration, retain explicit size but keep strict validation.

## Summary contract

The `summary` payload should use the repository-wide JSON report envelope and contain finite JSON values only. At minimum:

```json
{
  "methods": [
    "Rendered dimensions 1 and 2 from obsm['X_umap'] with Matplotlib; cells were colored categorically by leiden."
  ],
  "results": [
    "Plotted 4231 cells across 8 displayed categories; this is exploratory visualization, not inferential evidence."
  ],
  "key_results": {
    "embedding_key": "X_umap",
    "dimensions": [1, 2],
    "cells_total": 4231,
    "cells_plotted": 4231,
    "coordinate_ranges": {
      "UMAP1": {"min": -8.2, "max": 11.4},
      "UMAP2": {"min": -7.1, "max": 9.6}
    },
    "coordinates_finite": true,
    "embedding_provenance": {
      "verified": true,
      "neighbors_key": "neighbors",
      "min_dist": 0.5,
      "spread": 1.0,
      "random_seed": 0
    },
    "color": {
      "key": "leiden",
      "mode": "categorical",
      "missing": 0,
      "category_order": ["0", "1"],
      "category_counts": {"0": 123, "1": 98},
      "palette": {"0": "#1f77b4", "1": "#ff7f0e"},
      "legend_shown": true
    },
    "rendering": {
      "point_size_points_squared": 10.0,
      "alpha": 0.8,
      "sort_order": false,
      "dpi": 120
    }
  },
  "parameters": {},
  "references": [],
  "software_versions": {},
  "warnings": [],
  "limitations": []
}
```

For continuous color, replace category fields with count, missing count, min, max, median/quantiles, constant flag, colormap, normalization, effective `vmin`/`vmax`, and sort order. For no color, report the fixed point color.

Required warnings/limitations include, as applicable:

- UMAP provenance unavailable or incomplete;
- missing colors rendered with a special color;
- constant continuous color;
- categorical legend omitted;
- many categories reuse/approximate colors only if the palette policy permits it (prefer rejection or a generated nonrepeating mapping);
- UMAP axes, apparent global distance, islands, and mixing are not inferential measurements;
- cell-level Condition coloring is not a replicate-aware Condition contrast;
- Technical-batch visual mixing alone does not establish successful integration.

References should include the UMAP method/software papers, Scanpy computation/embedding documentation, and Matplotlib paper. Software versions should include Matplotlib, NumPy, pandas, anndata, and Scanpy only if its metadata/functions are actually used.

## Code contract

`code` must be compilable, self-contained Python with one callable such as:

```python
def plot_umap(adata):
    ...
    return png_bytes
```

It must reproduce:

- exact `embedding_key` and displayed dimensions;
- coordinate shape/type/finite validation;
- exact color-mode inference or explicit mode;
- missing-value policy;
- category order and literal color mapping, or continuous sort/normalization/cmap;
- point size, alpha, marker edge, figure size, labels, title, legend policy, DPI, and tight bounding box;
- input immutability and figure cleanup.

For categorical reproducibility, embed the resolved category-to-hex-color mapping rather than asking Matplotlib to choose it again. For continuous values, embed effective normalization limits. The function can omit OpenBio summary construction but cannot omit scientific failure checks or silently fall back from a missing custom key to `X_umap`.

## Failure contract and tests

### Validation failures

- backed/empty AnnData if unsupported;
- duplicate observation names;
- blank/non-string or absent embedding key;
- coordinate array not numeric, not 2D, row-misaligned, fewer than two dimensions, or nonfinite;
- blank color treated inconsistently with explicit none mode;
- missing color key;
- invalid color mode or a continuous request for nonnumeric data;
- infinite color values;
- invalid colormap/palette/missing color;
- point size boolean, nonnumeric, nonfinite, zero, or negative;
- too many categorical levels for a nonrepeating supported palette, unless a deterministic expansion policy exists.

### Required tests

1. Schema pins outputs `plot, summary, code` and input defaults/order.
2. Default and custom UMAP keys both render and report paired provenance.
3. Imported `X_umap` without provenance renders with a warning rather than false attribution.
4. Coordinate row mismatch, string/object dtype, NaN, and infinity fail before figure creation.
5. Categorical dtype order and boolean auto-mode are preserved.
6. Integer-coded metadata can be forced categorical and continuous.
7. Missing categorical and continuous values follow the disclosed missing-color policy; infinities fail.
8. Numeric sort order and effective normalization are pinned.
9. Category/color mapping contains no unintended repeats within the supported category limit.
10. Legend omission warning and summary mapping are present.
11. Caller `X`, obs, uns, and obsm remain byte/value equivalent.
12. Summary serializes with strict JSON and reports actual plotted cell count.
13. Generated code compiles and reproduces category mapping/continuous normalization and failures.
14. Repeated rendering and failure leave no registered/open figures.

## Migration

- Append outputs as `plot, summary, code`; existing output-zero Preview/SavePNG links remain meaningful, but serialized schemas need migration.
- Add `embedding_key=X_umap` without changing existing scientific behavior.
- Rename `color_map` to `continuous_color_map` if a separate categorical palette is added; migrate the old value only to the continuous branch.
- Map current `color=""` to explicit no-color mode.
- Map current nonempty colors to `mode=auto`; warn on legacy numeric integer columns where old behavior was continuous and new users may prefer categorical.
- Preserve the node ID and current default title/axes labels for default `X_umap` unless the reporting contract intentionally changes them.

## Atomicity and seam placement

- UMAP computation owns graph choice, initialization, optimization, and random seed.
- UMAP Plot owns coordinate/color validation, visual mapping, PNG rendering, and descriptive reporting.
- Preview Result and Save PNG own UI preview and filesystem output.
- A private scatter-rendering module is the correct seam for shared validation/palette/figure behavior. It should return bytes plus resolved rendering diagnostics, giving tests the same interface callers use.
- Do not make UMAP Plot write palettes or provenance into AnnData merely to reuse Scanpy plotting defaults; returning results rather than producing hidden side effects keeps the module deep and testable.

## Final-review resolution

The renderer accepts `(n_obs, n_dimensions)` for every `n_dimensions >= 2` and always displays the first two
columns. The summary fields `available_coordinate_dimensions` and `used_coordinate_dimensions` make this projection
explicit, and generated code enforces the same rule. This replaces the overly narrow exact-two-column validation;
three-dimensional UMAP output is an official upstream result, not malformed input.

The no-color branch is fully specified as `#246bfe`, alpha `0.8`, zero linewidth, and the same point-area policy as
colored branches. Missing categorical and continuous observations are rendered with the validated `missing_color`,
counted, and warned about explicitly. SciPy is always listed because sparse coordinate validation/materialization is
inside the implementation even when a particular input happens to be dense.
