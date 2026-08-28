# OpenBioSingleCellUMAPPlot — official usage research

Research date: 2026-08-28

## Scope

This node does not compute UMAP; it renders two stored coordinates. The relevant official interfaces are:

- `scanpy.tl.umap`, which defines where default and custom embeddings/parameters are stored;
- `scanpy.pl.umap` and the more general `scanpy.pl.embedding`;
- Matplotlib `Axes.scatter`, colormap registry, and normalization behavior;
- the original UMAP method/software papers for interpretation limits.

## Upstream UMAP storage contract

Official computation example:

```python
import scanpy as sc

sc.pp.neighbors(adata, key_added="neighbors")
sc.tl.umap(
    adata,
    neighbors_key="neighbors",
    min_dist=0.5,
    spread=1.0,
    n_components=2,
    init_pos="spectral",
    random_state=0,
)
```

With `key_added=None`, Scanpy stores coordinates in `adata.obsm['X_umap']` and parameters in `adata.uns['umap']`. With a custom `key_added`, both coordinates and parameters are stored under that exact key in `obsm` and `uns`. The OpenBio UMAP analysis node exposes this custom key, so a plotting node fixed to `X_umap` cannot consume all valid outputs of its upstream node.

Source: [Scanpy `tl.umap` reference](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html).

## Official Scanpy plotting usage

Default UMAP:

```python
figure = sc.pl.umap(
    adata,
    color="leiden",
    size=10.0,
    palette="tab20",
    color_map="viridis",
    sort_order=True,
    na_color="lightgray",
    show=False,
    return_fig=True,
)
```

Generic stored embedding:

```python
figure = sc.pl.embedding(
    adata,
    basis="custom_umap",
    color="leiden",
    size=10.0,
    palette="tab20",
    color_map="viridis",
    sort_order=True,
    na_color="lightgray",
    show=False,
    return_fig=True,
)
```

Official semantics:

- `basis` identifies the embedding in `obsm`.
- `sort_order=True` draws larger continuous color values above smaller values to reduce occlusion.
- `palette` is for categorical annotations and can be a qualitative colormap name, mapping, cycler, or sequence.
- `color_map`/`cmap` is for continuous values.
- `na_color` controls null or masked observations; `na_in_legend` controls their categorical legend entry.
- `size` is marker area.
- `color` can reference observations or genes, but the current OpenBio node intentionally limits it to observation metadata.

Sources: [Scanpy `pl.umap` reference](https://scanpy.readthedocs.io/en/1.10.x/api/generated/scanpy.pl.umap.html), [Scanpy `pl.embedding` reference](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html).

## Mutation behavior

Scanpy's categorical embedding plotting can read an existing `adata.uns['<color>_colors']` palette and can set that key when a palette is supplied. Passing caller-owned AnnData through this path therefore does not provide a robust read-only guarantee.

The current OpenBio node uses Matplotlib directly and does not mutate AnnData in its normal paths. That is a defensible implementation choice. If it is changed to `sc.pl.embedding`, it must plot a private copy and test caller state. Manual rendering still must copy/normalize arrays before sorting so that it never reorders or writes caller data.

## Matplotlib scatter semantics

Official usage:

```python
from matplotlib.figure import Figure

figure = Figure(figsize=(8, 7), constrained_layout=True)
axis = figure.subplots()
points = axis.scatter(
    x,
    y,
    c=continuous_values,
    s=10.0,
    cmap="viridis",
    alpha=0.8,
    linewidths=0,
)
figure.colorbar(points, ax=axis, label="score")
```

Matplotlib documents `s` as marker area in points squared, not radius. Marker edge linewidth can visibly enlarge small markers; a fixed `linewidths=0` or `edgecolors='none'` should be disclosed for predictable sizing.

For continuous values, Matplotlib normally applies linear normalization from the data range to `[0, 1]` and then maps through the colormap. A named colormap should be validated through the public `matplotlib.colormaps` registry before rendering.

For categorical values, use a qualitative palette rather than sampling a sequential metric colormap. Matplotlib's official colormap guidance distinguishes qualitative maps (`tab10`, `tab20`, `Set2`, and others) from perceptual maps intended for numeric magnitude.

Sources: [Matplotlib `Axes.scatter`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.scatter.html), [colormap normalization](https://matplotlib.org/stable/users/explain/colors/colormapnorms.html), [colormap registry](https://matplotlib.org/stable/api/cm_api.html), [choosing colormaps](https://matplotlib.org/stable/tutorials/colors/colormaps.html).

## Correct handling of color data

### Continuous

- Require a numeric, one-dimensional series aligned to `n_obs`.
- Distinguish finite values, missing values, and infinities. Infinities should fail; missing values need an explicit policy/color and count.
- Report the normalization type and effective `vmin`/`vmax`.
- Drawing higher values last, as Scanpy's `sort_order=True` does, is a meaningful display rule and must be fixed or exposed.
- A constant numeric column is valid but uninformative; warn and avoid a misleading gradient range.

### Categorical

- Preserve declared pandas categorical order and its `ordered` flag. For plain strings, choose and report a deterministic order such as first appearance.
- Boolean values should normally be categorical, not continuous 0/1.
- Integer-coded clusters are ambiguous; an `auto|categorical|continuous` color-mode control avoids silently assigning metric meaning.
- Create and report an explicit category-to-color mapping. The default Matplotlib property cycle has limited distinct colors and can repeat; silent reuse makes categories visually ambiguous.
- Missing values should be a separate disclosed state, not the literal category string `"nan"`.
- If the legend is omitted because of category count, retain the mapping in the JSON summary and warn that the PNG alone is not self-decoding.

## Coordinate validation

The plot requires a dense or safely materialized numeric array with:

- shape `(n_obs, n_dimensions)` with `n_dimensions >= 2`;
- the first two dimensions selected explicitly for this two-axis renderer;
- all displayed coordinates finite;
- observation order aligned to `adata.obs_names` by the AnnData `obsm` contract.

The current node checks only dimensionality and column count. It does not check row count, numeric dtype, or finite values. Matplotlib can silently omit nonfinite points, making the reported cell count differ from the displayed count.

## Reproducibility

The rendering step has no scientific RNG. Given fixed coordinate/color arrays and fixed ordering, it should be deterministic in one environment. PNG bytes can still vary with:

- Matplotlib, NumPy, pandas, and Pillow versions;
- renderer/backend and fonts;
- figure size, DPI, bounding-box layout, antialiasing, and legend layout;
- category ordering, palette assignment, and continuous draw order.

The plot should not expose a `random_seed`. It should report the upstream UMAP seed and parameters only when verifiable from OpenBio analysis history or the paired Scanpy `uns` metadata, labeling that value as embedding provenance rather than plot randomness.

## Scientific interpretation

UMAP is a nonlinear low-dimensional representation optimized from a neighborhood graph. A scatter plot is exploratory visualization, not an inferential measurement. Axis orientation, sign, origin, and absolute scale have no independent biological meaning. Dense islands, apparent separation, or visual mixing do not by themselves establish a cell type, a successful Technical-batch correction, a Condition effect, or preservation of biology.

Coloring by Sample, Technical batch, Condition, cluster, or annotation is descriptive. In particular, plotting individual cells by Condition does not create a replicate-aware **Condition contrast**; the **Sample** remains the inference unit.

## Primary references

- McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426. [doi:10.48550/arXiv.1802.03426](https://doi.org/10.48550/arXiv.1802.03426).
- McInnes L, Healy J, Saul N, Großberger L. UMAP: Uniform Manifold Approximation and Projection. *Journal of Open Source Software*. 2018;3(29):861. [doi:10.21105/joss.00861](https://doi.org/10.21105/joss.00861).
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0).
- Scanpy development team. [`scanpy.tl.umap`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html).
- Scanpy development team. [`scanpy.pl.embedding`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html).
- Hunter JD. Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*. 2007;9(3):90-95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55).

## Final-review clarification

Scanpy/AnnData embeddings are not restricted to two stored dimensions. A two-axis UMAP renderer must accept any
numeric `obsm` value shaped `(n_obs, n_dimensions)` with at least two columns, display columns one and two, and
report both the available dimension count and the displayed dimension identities. Extra dimensions are retained in
the caller object and are not copied into the renderer. Only displayed coordinates must be finite; an unused third
coordinate does not justify silently changing or rejecting the displayed two-dimensional view.

The read-only implementation uses SciPy's sparse interface when validating/materializing a sparse stored embedding,
so SciPy is an actually used software dependency and belongs in the report versions/references. When `color` is
empty, the fixed point color is `#246bfe` and must be disclosed. When categorical or continuous color contains
missing observations, the summary must report their count and explicit special color and the warnings must say that
the special missing-value color was used; the PNG alone is otherwise ambiguous.
