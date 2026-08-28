# OpenBioSingleCellMarkerExpressionPlot — official usage research

Research date: 2026-08-28

## Scope

This review covers the public Scanpy plotting functions currently selected by the node:

- `scanpy.pl.dotplot`
- `scanpy.pl.matrixplot`
- `scanpy.pl.tracksplot`
- `scanpy.pl.violin`
- `scanpy.tl.dendrogram`, because three plot modes can invoke it implicitly

The installed test environment was Scanpy 1.12.3. The stable/latest public documentation was also checked because current Scanpy exposes a preview plotting backend whose return types differ from the legacy Matplotlib backend.

## Official semantics

### Dot plot

Official usage:

```python
import scanpy as sc

markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
plot = sc.pl.dotplot(
    adata,
    markers,
    groupby="cluster",
    use_raw=False,
    layer="log1p_norm",
    expression_cutoff=0.0,
    mean_only_expressed=False,
    standard_scale="var",
    dendrogram=False,
    show=False,
    return_fig=True,
)
plot.make_figure()
figure = plot.fig
```

For each gene/group pair, dot color represents mean expression and dot size represents the fraction of cells above `expression_cutoff`; the default cutoff is zero. `standard_scale="var"` or `"group"` rescales the selected dimension to `[0, 1]` after subtracting its minimum. `return_fig=True` returns a `DotPlot` object rather than a Matplotlib axes mapping.

Source: [Scanpy `pl.dotplot` reference](https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.dotplot.html).

### Matrix plot

Official usage:

```python
plot = sc.pl.matrixplot(
    adata,
    markers,
    groupby="cluster",
    use_raw=False,
    layer="log1p_norm",
    standard_scale="var",
    dendrogram=False,
    show=False,
    return_fig=True,
)
plot.make_figure()
figure = plot.fig
```

The matrix color encodes mean expression for each gene/group pair. It supports the same `use_raw`, `layer`, category ordering, standard scaling, and dendrogram concepts as the legacy dot plot.

Source: [Scanpy `pl.matrixplot` reference](https://scanpy.readthedocs.io/en/latest/api/scanpy.pl.matrixplot.html).

### Tracks plot

Official usage:

```python
axes = sc.pl.tracksplot(
    adata,
    markers,
    groupby="cluster",
    use_raw=False,
    layer="counts",
    dendrogram=False,
    log=False,
    show=False,
)
figure = next(iter(axes.values())).figure
```

The legacy Matplotlib implementation draws one filled expression track per gene across cells ordered by the categorical grouping. The official documentation says that its best results are obtained from raw counts that are not logged. It requires `groupby` to order observations and documents a categorical grouping as the intended input. Unlike dot and matrix plots, the public function has no `standard_scale` parameter.

Source: [Scanpy `pl.tracksplot` reference](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.tracksplot.html).

### Violin plot

Official usage:

```python
axes_or_grid = sc.pl.violin(
    adata,
    markers,
    groupby="cluster",
    use_raw=False,
    layer="log1p_norm",
    stripplot=False,
    density_norm="width",
    multi_panel=len(markers) > 1,
    show=False,
)
```

The legacy implementation wraps Seaborn's violin plot. With `density_norm="width"`, every violin has the same maximum width; width therefore does not encode group sample size. `density_norm="count"` is the documented option whose width reflects observation count. `stripplot=False` suppresses individual points. The function has no dendrogram or `standard_scale` parameter.

Source: [Scanpy `pl.violin` reference](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.violin.html).

## Expression-source precedence

The plotting functions expose `use_raw` and `layer`; their documentation states that a named layer takes precedence over `use_raw`. If neither is controlled, several Scanpy plotters default to `adata.raw` when it is present. The node is correct to resolve the source explicitly and pass both values rather than accepting Scanpy's context-dependent default.

Gene names must belong to the selected expression representation. `raw` can have a different feature axis from `adata.var_names`, so validation and `input_genes` reporting must use `adata.raw.var_names`/`adata.raw.n_vars` when raw is selected.

## The `log` parameter is not a safe shared abstraction

The public references describe `log=True` as plotting on a logarithmic axis. The current Scanpy 1.12.3 legacy group-plot implementation prepares dot/matrix/tracks data through `_prepare_dataframe`, where `log=True` applies `log1p` to expression values before aggregation, while the violin implementation uses its own plotting path. Consequently one shared node parameter named `log` does not guarantee the same scientific transformation across the four modes.

OpenBio should make the expression representation explicit through the source selector. A plotting node should not silently become another normalization/transformation node. If a violin logarithmic y-axis remains desirable, it should be a violin-specific display option named `y_scale`, not a shared expression-processing switch.

Primary technical sources: the public function references above and their linked official Scanpy source.

## Dendrogram behavior and mutation

The public dendrogram function is:

```python
sc.tl.dendrogram(
    adata,
    groupby="cluster",
    var_names=markers,
    use_raw=False,
    cor_method="pearson",
    linkage_method="complete",
    optimal_ordering=False,
    key_added="_openbio_marker_dendrogram",
    inplace=True,
)
```

Important official behavior:

- Hierarchical clustering is performed on group-level mean vectors, not individual cells.
- With no `var_names`, Scanpy selects PCA by default unless `.X` has fewer than 50 variables; it can even compute PCA if needed.
- With `var_names`, the selected genes are used and `use_rep`/`n_pcs` are ignored.
- `use_raw` is honored only when `var_names` is supplied.
- There is no public `layer` parameter in `tl.dendrogram`.
- The defaults are Pearson correlation, complete linkage, and `optimal_ordering=False`.
- With `inplace=True`, results are stored in `adata.uns[key_added]`, by default `dendrogram_{groupby}`.

Source: [Scanpy `tl.dendrogram` reference](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dendrogram.html).

When `dendrogram=True` is passed to dotplot, matrixplot, or tracksplot and the expected key is absent, the plotting implementation automatically calls `tl.dendrogram` with default parameters and writes the result to `adata.uns`. Therefore the current node's claim that it renders "without modifying AnnData" is false in that path.

The implicit default is also scientifically misaligned with the node's expression-source interface: it may order groups from PCA or `.X`, while the visible marker panel is drawn from `raw` or a named layer. The tree then does not necessarily summarize the displayed genes or displayed expression state.

The safe wrapper pattern is:

1. Never pass caller-owned `adata` to an implicit dendrogram path.
2. Create a private plotting copy.
3. Construct the dendrogram from the same selected genes and the same selected expression matrix. For a named layer, use a private minimal AnnData whose `.X` is that selected layer, because `tl.dendrogram` has no `layer` argument.
4. Request `inplace=False` or store the result under a collision-safe private key on the private copy.
5. Pass that explicit key to the plotter.
6. Report the correlation method, linkage method, selected genes, and resulting category order.

## Plotting backend stability

Current Scanpy documentation states that the active plotting backend is selected by `scanpy.settings.preset`. `ScanpyV1` uses the legacy Matplotlib interface assumed by the node; `ScanpyV2Preview` can return HoloViews objects with different parameters and return types. Scanpy documents `settings.override()` as the local context manager for settings.

For a PNG-producing node, use a local context such as:

```python
with sc.settings.override(preset=sc.Preset.ScanpyV1, autoshow=False):
    ...
```

Do not mutate global Scanpy settings without restoring them, and do not rely on the caller's current preset.

Sources: [Scanpy settings](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.settings.html), [Scanpy preset](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.Preset.html).

## Reproducibility

These renderers do not expose a random seed in their documented legacy interfaces when `stripplot=False`. Reproducibility depends on:

- the exact expression matrix and group/category order;
- whether group ordering is observed or dendrogram-derived;
- explicit aggregation, cutoff, scaling, and density settings;
- Scanpy, Matplotlib, Seaborn, NumPy, pandas, and SciPy versions;
- Matplotlib backend, font availability, figure size, DPI, and save options.

A dendrogram is deterministic for fixed inputs and method parameters, but ties and dependency versions can affect ordering. PNG bytes should not be promised bit-for-bit stable across software/font environments.

## Scientific interpretation

The output is **Cluster marker evidence** under the repository vocabulary. It helps annotation review; it is not a **Condition contrast**, does not establish a **Curated annotation**, and must not describe a cluster label or marker panel as biological truth.

Dot/matrix means and cell fractions are descriptive cell-level summaries. Violin and track shapes are descriptive distributions. None are replicate-aware inferential results, and cell count is not biological replicate count. Standard scaling further removes absolute magnitude comparability along the scaled dimension.

## Primary references

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0).
- Scanpy development team. [`scanpy.pl.dotplot`](https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.dotplot.html).
- Scanpy development team. [`scanpy.pl.matrixplot`](https://scanpy.readthedocs.io/en/latest/api/scanpy.pl.matrixplot.html).
- Scanpy development team. [`scanpy.pl.tracksplot`](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.tracksplot.html).
- Scanpy development team. [`scanpy.pl.violin`](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.violin.html).
- Scanpy development team. [`scanpy.tl.dendrogram`](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dendrogram.html).
- Hunter JD. Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*. 2007;9(3):90-95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55).

## Final-review clarification

### Raw is a container, not an expression-state guarantee

AnnData `raw` preserves an expression matrix and feature annotations; the official interface does not declare that
the matrix contains integer counts. External Scanpy objects commonly use `raw` for log-normalized full-feature
expression. The node may report Raw as `counts` only when OpenBio Snapshot Expression history binds the Raw
destination to the count snapshot. It may report another state when selected Raw values equal a current-X panel with
independently verifiable state provenance. Otherwise the state is `unknown`, values are plotted unchanged, and the
summary warns that it did not assume counts, normalization, or logarithmization.

### Individual-cell jitter consumes NumPy global RNG

`sc.pl.violin(..., stripplot=True, jitter=...)` delegates individual-cell point placement to the legacy
Scanpy/Seaborn Matplotlib path. Empirical contract testing against the supported versions confirms that jitter draws
from NumPy's legacy global RNG. Therefore `show_cells=True` needs an explicit rendering seed. The safe wrapper saves
the caller's complete `np.random` state, seeds locally, performs plotting, and restores the saved state in `finally`
on success or failure. The same seed must be embedded in generated code. Seaborn is an actually used/reported
dependency only for the violin branch; SciPy is used for source sparse/dense handling in every branch.

### Resource-safe adapter shape

Passing `adata.copy()` to Scanpy copies unrelated X, Raw, layers, embeddings, graphs, and metadata, so a budget based
only on the selected panel is false. The correct adapter is a minimal private AnnData whose X is exactly the selected
source-by-gene panel, whose obs contains only the validated grouping column, and whose var index contains only the
selected genes. Dendrogram computation receives a separate minimal panel copy. A pre-materialization resource
estimate covers the selected slice, dense panel, plotting/dendrogram matrix copies, group aggregates/correlation,
summary rows, and figure/PNG working buffers.
