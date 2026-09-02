# Core AnnData plots

Fixed point: `24a018fca9621a484afede7f842e8bc491bb4744`

## Goal

Complete the general AnnData plotting path without adding downstream-analysis-specific plots or a generic plotting DSL.

## Required behavior

1. Upgrade the existing UMAP plot node in place to a general **Embedding Plot** while preserving its node ID. It reads any validated `obsm` coordinates, exposes a selectable one-based dimension pair, and renders one target per run: no color, one `obs` column, or one gene from an explicit X/Raw/layer Expression source. Existing categorical and continuous metadata behavior remains available. Generated code, summary, PNG output, and input immutability remain equivalent.
2. Add **Sample Composition Plot** as a read-only consumer of the canonical table emitted by Sample Composition Summary. It renders Sample-level stacked bars with `proportion` as the default and `cell_count` as the alternate value. Samples are grouped by Condition; the plot performs no Condition aggregation or inference.
3. Add **HVG Selection Plot** as a read-only AnnData plot. It uses the stored flavor-appropriate mean and variability columns, distinguishes algorithm-selected, forced, and unselected genes, and fails clearly when the complete feature-selection evidence needed by the plot is unavailable.
4. Add **PCA Variance Plot** as a read-only AnnData plot over stored PCA variance ratios. It renders per-PC explained variance and cumulative explained variance without recomputing PCA. `n_pcs=0` means all available components.
5. Extend **QC Plots** with `overview` and `grouped` views. The grouped view requires one categorical `obs` key and renders the same source-coherent total-expression, detected-gene, and optional mitochondrial distributions by group. The default overview remains unchanged.
6. Rename only the display name of Marker Expression Plot to **Grouped Gene Expression Plot**; retain its node ID and four existing modes.

## Public test seams

- Node schemas: exact visible inputs, defaults, dynamic branches, display names, categories, and `plot, summary, code` outputs.
- Owned scientific operations: valid PNG, strict JSON summary, compilable/equivalent generated code, exact expression/table/AnnData source selection, realistic validation failures, and caller input immutability.
- Worker operations: typed input artifacts remain unchanged and plot artifacts round-trip through the existing PNG codec.

## Constraints

- Reuse the existing PlotResult, worker, artifact, Preview, Save PNG, and Persist paths. Do not add a dependency, frontend protocol, plot artifact type, or interactive renderer.
- Preserve Sample as the descriptive composition unit and keep Cluster marker evidence separate from Condition inference.
- Every expression-derived quantity in one plot uses one explicit Expression source.
- Do not add Highest Expressed Genes, generic x-y scatter, marker-ranking plots, PCA loadings, additional expression modes, or downstream-specific plots in this change.
