# openbio-singlecell

`openbio-singlecell` is an independent ComfyUI custom-node pack for in-memory single-cell RNA-seq analysis. It works with the standard ComfyUI frontend and with the paired OpenBio frontend build. The pack does not modify ComfyUI's executor, cache, queue, WebSocket protocol, or built-in nodes.

## Supported baseline

- ComfyUI commit `924743af083c151296cc16f925aeab113b6484e8`
- ComfyUI version 0.33.0 or newer
- ComfyUI_frontend commit `795292d90777efdfb92fc5d3f928dcf13a7cf8a8` / version 1.52.3
- Python 3.12 or newer
- Node.js 25 and pnpm 11.13.1 for building the paired frontend
- `scanpy[leiden]>=1.12.3,<1.13`
- `anndata>=0.13.2,<0.14`

Notebook-derived optional nodes import their own libraries only when executed. Install the relevant packages in the ComfyUI Python environment for Harmony, scVI, CellTypist, decoupler, pertpy, LIANA, scVelo, infercnvpy, Schist, OmicVerse cNMF, Cassiopeia, or pySCENIC/loompy analyses.

The exact source pairing is recorded in `release_manifest.json`.

## Install with ComfyUI Manager

ComfyUI Manager/Registry is the production installation path for this pack. After the first Registry release is published, search for `openbio-singlecell` in Manager, install it, and restart ComfyUI. Manager owns the installed copy under `ComfyUI/custom_nodes`; the plugin source does not belong in the ComfyUI source repository. Manager installs `requirements.txt` and then the root `install.py` generates the local demonstration H5AD. That entry point does not install packages or make network requests.

This checkout is not yet claiming a public Registry listing. Formal publication still requires two real values that must not be guessed:

- the Publisher ID created for the project in the Comfy Registry;
- the public Git repository URL for this independent repository.

Before publishing, add the actual URL as `[project.urls].Repository` in `pyproject.toml` and add the exact `PublisherId` under `[tool.comfy]`. Keep Registry API keys in the publisher's secret store, never in this repository. Once those values exist, validate and inspect the package before publishing:

```sh
comfy node validate
comfy node pack
comfy node publish
```

Registry archives are built from Git-tracked files. The included `.comfyignore` keeps tests out of the installed node pack while retaining the example workflows, local demo generator, launch scripts, licenses, and release manifest.

## Local source development

Keep the backend, frontend, and plugin as three sibling repositories:

```text
parent/
  ComfyUI/
  ComfyUI_frontend/
  openbio-singlecell/
```

ComfyUI discovers custom nodes from its custom-node roots, so expose only this plugin directory at its normal load location. On Windows, create an exact NTFS junction from the parent directory:

```powershell
New-Item -ItemType Junction `
  -Path .\ComfyUI\custom_nodes\openbio-singlecell `
  -Target .\openbio-singlecell
```

On Linux or macOS, create the equivalent symbolic link:

```sh
ln -s "$PWD/openbio-singlecell" "$PWD/ComfyUI/custom_nodes/openbio-singlecell"
```

Do not also expose the same checkout through `extra_model_paths.yaml`; loading one plugin through two custom-node roots can register it twice. A Manager installation is a real directory under `custom_nodes` and does not need this development link.

For a sibling source checkout, install its Python dependencies and create the local demo AnnData from the plugin root. Use the same Python interpreter that runs ComfyUI.

Windows PowerShell:

```powershell
cd .\openbio-singlecell
powershell -ExecutionPolicy Bypass -File .\scripts\install.ps1 `
  -ComfyRoot ..\ComfyUI `
  -Python D:\path\to\python.exe
```

Linux or macOS:

```sh
cd openbio-singlecell
OPENBIO_COMFYUI_ROOT=../ComfyUI \
OPENBIO_PYTHON=/path/to/comfyui/python \
sh scripts/install.sh
```

When the pack is installed directly below `ComfyUI/custom_nodes`, the installer can infer the ComfyUI root. The installer only installs this pack's `requirements.txt`, checks imports, and generates the local demo AnnData. It uses `uv pip` when `uv` is available, so a normal `uv venv` does not need an embedded `pip`; otherwise it falls back to `python -m pip`. Use `-ForceDemo` on Windows or `--force` on Linux to replace that AnnData file.

The generated file is:

```text
ComfyUI/input/openbio-singlecell/openbio_singlecell_demo.h5ad
```

It contains 600 cells, 500 genes, three known cell groups, mitochondrial genes, `sample`, `condition`, and `batch`
metadata, known marker genes, and CSR sparse integer counts. A valid existing file is left unchanged; replacing an
invalid file requires the explicit force option.

## Example workflows

Open one of the five production starting-point templates in `example_workflows`:

- `Quality Control and Clean Counts.json` loads raw counts, calculates and filters QC metrics, previews the
  retained data, and writes a clean H5AD only through the explicit save node.
- `Cell Clustering and Marker Discovery.json` reads raw counts from `adata.X`, leaves `X` and existing layers
  unchanged, and creates a normalized `log1p_norm` layer for highly variable genes, PCA, clustering, marker ranking,
  plots, and explicit outputs. It deliberately omits Scale so the production path does not densify the expression
  matrix merely to reach PCA.
- `Sample Composition Comparison.json` loads AnnData that already contains sample metadata and fans reusable
  sample, condition, annotation, reference, and comparison strings into two composition nodes whose concrete
  table outputs can be previewed or exported.
- `scVI Batch Integration and Contrast.json` trains scVI from raw counts in `adata.X`, constructs neighbors
  from `X_scVI`, visualizes the integrated embedding, and passes the concrete in-memory model to scVI differential
  expression instead of retraining it.
- `Single-Cell Best Practice.json` is the production-oriented end-to-end example for the supplied multi-sample
  study. It re-audits QC, removes predicted doublets, applies sample-wise MAD rules, preserves full-gene counts in
  both `layers["counts"]` and `raw`, trains scVI on a 5,000-HVG count view, joins its Leiden labels back to the
  full-gene object, compares two exploratory marker rankings, adds provisional CellTypist labels, and runs one
  explicit PyDESeq2 pseudobulk contrast.

Every workflow JSON has a same-stem 768 x 768 JPEG cover in the same directory. The JSON and cover lists are
explicit release artifacts, so a renamed, missing, extra, malformed, or mismatched template asset fails the release
contract tests. `example_workflows` remains the only packaged workflow source; there is no mirrored copy under
`web` that can drift from the reviewed files.

The bundled values are reviewed starting points for the generated demonstration data, not universal scientific
thresholds. Before using a template in production, replace the input path and review its visible QC thresholds,
feature count, grouping columns, labels, comparison groups, and output names. In particular:

- QC assumes gene symbols use the `MT-` mitochondrial prefix. Its starting cell thresholds are 90 minimum genes
  and 20% maximum mitochondrial counts, and its starting gene threshold is 3 minimum cells.
- composition requires `adata.obs["sample"]`, `adata.obs["condition"]`, and `adata.obs["cell_type"]`. The demonstration
  compares the `control` and `treated` biological conditions, balanced across `adata.obs["batch"]`. Its Kruskal-Wallis
  and Mann-Whitney results are an exploratory sample-level comparison; use a covariate-aware compositional model when
  the study design requires it.
- QC and the complete scVI template require non-negative integer counts in `adata.X` at entry because QC and
  filtering operate on `X`. If counts only exist in a layer, first use Use Expression Layer to restore that layer to
  `X`, then keep the scVI template's source set to `X`. The clustering template can instead read a counts layer by
  changing Normalize to Layer's `source` and `source_layer`; it does not overwrite `X` or existing layers.
- scVI also requires `adata.obs["batch"]`; its downstream contrast requires the grouping column and labels shown on
  that node. Install the optional `scvi-tools` dependency in the same Python environment that runs ComfyUI before
  executing this template.
- `Single-Cell Best Practice` expects the external file
  `ComfyUI/input/openbio-singlecell/anndata_qc.h5ad`; it is not bundled with the plugin or generated by the
  installer. The reviewed contract uses raw integer counts in `X` and the `sample`, `batch`, and `group` observation
  columns. Install `scvi-tools`, `celltypist`, and `pertpy` with its PyDESeq2 requirements before executing the
  corresponding branches.

The best-practice workflow keeps two deliberate output objects. The full-gene H5AD retains active full-gene counts,
`raw`, `log1p_norm`, HVG flags, provisional CellTypist labels, and the joined `leiden_scvi` column, but it does not
contain the scVI embedding or graph. The HVG-scVI H5AD retains the 5,000-gene active matrices, `X_scVI`, neighbors,
UMAP, and Leiden results; AnnData's full-gene `raw` snapshot remains available for recovery after variable
subsetting. Neither H5AD stores the trained scVI model weights.

Within this existing-node implementation, sample-wise MAD decisions are exact but the before/after QC plots are
global rather than sample-faceted. The HVG and full-gene top-marker tables are Scanpy cluster rankings used as
annotation evidence; they are not scVI all-cluster DE or replicate-aware condition tests. The reference dotplot is a
fixed canonical marker panel because marker tables cannot currently drive a gene-list input dynamically. CellTypist
uses `Adult_Human_Vascular.pkl` as a provisional reference. The pseudobulk branch selects its provisional
`smc_pc_intermediate` label before fitting `~group` for `nonDM_ED` versus `Normal`; this is a real model label but
is only a runnable starting choice. Replace it with a populated provisional label after reviewing CellTypist output,
or connect the same branch to the later curated `cell_type` column. Duplicate the explicit population-selection
branch for additional cell types.

The scVI `model` output is an ephemeral Python object for downstream nodes in the current execution. Saving the
workflow records the connection and parameters, not trained model weights; after restarting ComfyUI, rerun scVI
Integration before executing its differential-expression consumer.

## Study parameters

`Core Study Parameters` is an optional, interactive parameter source. It outputs six ordinary `STRING` values:
`sample_column`, `condition_column`, `batch_column`, `annotation_column`, `reference`, and `comparison`. Connect only
the values an analysis needs, and fan out the same output to several nodes when that avoids repeated typing. The node
does not load or modify AnnData, bundle the selections into an opaque design object, or make analysis nodes depend on
it; every compatible string input remains independently editable when it is not connected. Its custom parameter
editor is implemented entirely inside this plugin and does not patch ComfyUI or ComfyUI_frontend.

## Run modes

### Standard ComfyUI frontend

Start ComfyUI normally after installation. All `openbio/single-cell/*` nodes and their node-local result previews remain available without the paired OpenBio frontend build.

### OpenBio frontend suite

From the three-repository parent layout, build the paired frontend and generate its license inventory:

```sh
cd ComfyUI_frontend
corepack pnpm install --frozen-lockfile
corepack pnpm build:openbio
corepack pnpm exec node ../openbio-singlecell/scripts/generate_frontend_notices.mjs
```

The notices script defaults to the sibling `ComfyUI_frontend` checkout. It copies the frontend GPL license and generates a deduplicated production-dependency notice at `dist/LICENSE` and `dist/THIRD_PARTY_NOTICES.md`. Run it after every frontend build because Vite replaces `dist`. Verify the existing files without rewriting them with:

```sh
corepack pnpm exec node ../openbio-singlecell/scripts/generate_frontend_notices.mjs --check
```

Then launch the suite from the plugin repository.

Windows PowerShell:

```powershell
cd ..\openbio-singlecell
powershell -ExecutionPolicy Bypass -File .\scripts\start.ps1 `
  -ComfyRoot ..\ComfyUI `
  -Python D:\path\to\python.exe
```

Linux or macOS:

```sh
cd ../openbio-singlecell
OPENBIO_COMFYUI_ROOT=../ComfyUI \
OPENBIO_PYTHON=/path/to/comfyui/python \
sh scripts/start.sh
```

The start scripts require `ComfyUI_frontend/dist/index.html`, `dist/LICENSE`, and `dist/THIRD_PARTY_NOTICES.md`, then launch ComfyUI with:

```text
--disable-api-nodes --enable-assets --front-end-root <sibling frontend dist>
```

Additional arguments are forwarded to ComfyUI. ComfyUI's own requirements must already be installed. In particular, `--enable-assets` requires ComfyUI's database dependencies and a working local database connection.

## Data and output behavior

- AnnData values use `adata` ports and the `OPENBIO_ANNDATA` wire type.
- Core Study Parameters exposes its six selections as standard ComfyUI `STRING` wires, so they can connect to any
  compatible string input.
- Tables, plots, and structured summaries use distinct `table`, `plot`, and `summary` ports with the
  plugin-owned `OPENBIO_SINGLE_CELL_TABLE`, `OPENBIO_SINGLE_CELL_PLOT`, and
  `OPENBIO_SINGLE_CELL_SUMMARY` wire types.
- Modifying nodes copy their input before changing it; read-only nodes do not copy it.
- AnnData transforms output `adata`; analysis artifacts output their concrete table, plot, or summary type;
  preview and save nodes are terminal and have no data output.
- Preview accepts all three artifact types, while CSV and PNG outputs accept only tables and plots respectively,
  so incompatible links are rejected before execution.
- Frequently tuned analysis choices, expression sources, and result-defining thresholds stay visible. Core Study
  Parameters is an optional source for reusing sample, condition, annotation, and primary-contrast strings; analysis
  nodes keep their ordinary explicit inputs and do not consume a combined design object. Random seeds, internal
  storage keys, output column names, and iteration limits are advanced inputs.
- Dataset-specific condition, tumor, cluster, and cell-type values are never supplied as defaults; nodes that need them require an explicit value.
- Load H5AD and Load 10x H5 accept a browser file selection or a file dropped directly onto the node; uploads are
  stored under `ComfyUI/input/openbio-singlecell` and the node keeps the resulting relative path.
- Inputs are limited to relative paths under the ComfyUI `input` directory and are checked again at the read boundary.
- Temporary plots are written below `temp/openbio-singlecell`; ComfyUI clears its temp directory at startup.
- CSV, PNG, and H5AD files are only made permanent by explicit output nodes and are written below `output/openbio-singlecell`.
- Node UI payloads are bounded; complete tables are exported as CSV.
- DataFrames remain inside table artifacts and plots carry rendered PNG data. Library-specific models, trees,
  figures, and other opaque Python objects remain node implementation details unless a workflow has a concrete
  reusable downstream contract for them.
- Trained scVI models, pySCENIC regulatory networks, and solved Cassiopeia trees have explicit
  `OPENBIO_SCVI_MODEL`, `OPENBIO_SCENIC_NETWORK`, and `OPENBIO_CASSIOPEIA_TREE` contracts because each has a
  real downstream consumer. No generic Python-object or generic model wire is exposed.
- A transformation keeps one `adata` output when it has no reusable secondary product. The deliberate exceptions
  are scVI Integration (`adata`, `model`) and Run pySCENIC (`adata`, `network`); Cassiopeia reconstruction produces
  one reusable `tree` for its expansion and plasticity consumers.

## Scope and limitations

The current in-memory workflow includes 10x study loading, QC, expression snapshots, layer-aware preprocessing, Harmony/scVI integration, clustering, annotation, pseudobulk differential analysis, enrichment, LIANA communication, compositional abundance testing, OmicVerse cNMF, pySCENIC, Cassiopeia lineage analysis, PAGA/DPT, RNA velocity, and inferCNV. The pack does not declare a maximum AnnData size and does not yet provide AnnData backed mode, out-of-core processing, automatic disk caching, ATAC, or spatial analysis. Notebook stages implemented only in R, destructive file organization, publication-only plots, the empty Geneformer notebook, and the import-only spatial notebook are intentionally not represented as nodes. Scale and GSVA may densify data as their underlying libraries normally do. A running scientific operation may finish before a stop request takes effect; stopping prevents later nodes from starting.

## Offline behavior

The plugin contains no telemetry, analytics, update checks, or remote configuration. The OpenBio launcher disables Comfy API nodes. Optional scientific libraries may download a selected model or knowledge resource, such as CellTypist, CollecTRI, or DGIdb data, when it is not already cached; prepare those resources in advance for an offline run. The installer may contact the Python package index configured for `pip`.

## License and attribution

`openbio-singlecell` is licensed under GPL-3.0-or-later. The complete license text is included in `LICENSE`. It is based on ComfyUI and its paired ComfyUI_frontend, both distributed under GPLv3 terms. Keep upstream copyright notices, source attribution, and the corresponding GPL license text with redistributed source or binaries. See `THIRD_PARTY_NOTICES.md` for source-package notices; paired prebuilt frontend releases also carry their generated production-dependency notice inside `dist`.
