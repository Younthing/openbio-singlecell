# openbio-singlecell

`openbio-singlecell` is an independent ComfyUI custom-node pack for file-artifact single-cell RNA-seq analysis. It works with the standard ComfyUI frontend and with the paired OpenBio frontend build. The pack does not modify ComfyUI's executor, cache, queue, WebSocket protocol, or built-in nodes.

## Supported baseline

- ComfyUI commit `924743af083c151296cc16f925aeab113b6484e8`
- ComfyUI version 0.33.0 or newer
- ComfyUI_frontend commit `795292d90777efdfb92fc5d3f928dcf13a7cf8a8` / version 1.52.3
- Python 3.12 or newer
- Node.js 25 and pnpm 11.13.1 for building the paired frontend
- `scanpy[leiden,scrublet]>=1.12.3,<1.13`
- `anndata>=0.13.2,<0.14`

Notebook-derived optional nodes import their own scientific backends only when executed. The core
`requirements.txt` deliberately does not install those large or platform-specific stacks.

### Optional backend installation

Install an extra from the plugin root with the same Python interpreter that runs ComfyUI, for example
`python -m pip install ".[velocity]"`. The reviewed extras are:

| Extra | Reviewed backend | Node family |
| --- | --- | --- |
| `celltypist` | `celltypist>=1.7,<2` | CellTypist annotation |
| `decoupler` | `decoupler==2.2.0` | scoring, enrichment, CollecTRI and pseudobulk preparation |
| `harmony` | `harmonypy==2.0.0` | Harmony integration |
| `scvi` | `scvi-tools>=1.5,<1.6` | scVI integration and model differential evidence |
| `pertpy` | `pertpy==1.3.0` | Augur and drug perturbation nodes |
| `composition` | `pertpy[tcoda]==1.3.0` | scCODA and tascCODA, including the audited JAX/tree stack |
| `milo` | `pertpy[milo-edger]==1.3.0` | Milo plus its Python-to-R bridge |
| `pseudobulk` | `decoupler==2.2.0`, `pertpy[de]==1.3.0`, `pydeseq2>=0.5,<0.6` | pseudobulk preparation, edgeR and PyDESeq2 |
| `liana` | `liana==1.9.0`, `pandas<3` | LIANA communication and its result plot |
| `velocity` | `scvelo==0.3.4` | RNA-velocity processing and analysis |
| `cnv` | `infercnvpy==0.6.1` | inferCNV, CNV PCA and group scores |
| `cnmf` | `cnmf==1.7.1` | method-author cNMF rank survey and consensus programs |
| `lineage` | `cassiopeia-mt==2.1.3; platform_system == 'Linux'` | Cassiopeia lineage preparation, reconstruction and analysis on the verified Linux/WSL release boundary |
| `forceatlas2` | `fa2-modified==0.4` | explicit ForceAtlas2 graph layout |

The `milo` extra installs `rpy2`, but pip cannot install R or Bioconductor packages. Milo additionally requires an
external R runtime with `edgeR`, `limma`, and `statmod`. The edgeR pseudobulk engine requires R with `edgeR`,
`BiocParallel`, and `RhpcBLASctl`. Install those packages in the R installation used by `rpy2`; neither engine
silently falls back to PyDESeq2. The `pseudobulk` extra is sufficient for the audited PyDESeq2 path.

`harmonypy==2.0.0` publishes official wheels for Linux and macOS, not native Windows. On Windows, run the
plugin under WSL with a supported Linux wheel or provide and validate a supported source build; the extra must not
be expected to produce a native Windows wheel.

Schist is intentionally not a pip extra. The package named `schist` on PyPI is an unrelated knowledge-graph
project. Install the audited single-cell Schist and graph-tool into the same environment that runs ComfyUI using
conda-forge (or build the official sources and graph-tool yourself):

```sh
conda install -c conda-forge "schist=0.10.0" graph-tool
```

Graph-tool's supported ordinary paths are Linux and macOS; use WSL on Windows. The Schist nodes verify the official
package identity, exact version, graph-tool backend, and public API before analysis.

pySCENIC is an external workflow, not an in-process extra. pySCENIC 0.12.1 is incompatible with this plugin's
Python 3.12+/NumPy 2 baseline. Run the official `aertslab/pyscenic:0.12.1` image or a dedicated Python 3.10
environment pinned to `pyscenic==0.12.1`, `ctxcore==0.2.0`, `arboreto==0.1.6`, `loompy==3.0.8`,
`numpy>=1.21,<1.24`, and `pandas>=1.3.5,<2`, then import the complete
`openbio-singlecell/pyscenic-external-run/v1` bundle. OpenBio does not patch removed NumPy aliases or start a
container engine on the user's behalf.

The Cassiopeia extra is intentionally marked for Linux only. This is the verified release boundary for
`cassiopeia-mt==2.1.3`; use WSL/Linux on Windows. Native Windows and macOS execution are not claimed by this release.

OmicVerse is not declared as an extra. Its 2.3.1 distribution requires `anndata<0.12.0`, which cannot coexist with
this release's `anndata>=0.13.2,<0.14` core contract. The cNMF extra uses the method authors' `cnmf==1.7.1`
distribution; installing OmicVerse does not satisfy that adapter.

All extras remain opt-in and are not installed by ComfyUI Manager's core `requirements.txt` path.

CellTypist models are also opt-in. Before running the node, either provide a trusted local `.pkl` path or deliberately
download one with CellTypist's official API (replace the example model as appropriate):

```python
from celltypist import models

models.download_models(model="Adult_Human_Vascular.pkl")
```

That explicit command uses the network; node execution itself never downloads or enumerates the remote catalog.
Python pickle model files can execute code when loaded, so use only an official or otherwise trusted source. A
reported SHA-256 identifies the bytes used but does not establish that the file is safe or biologically appropriate.

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
  Sample, Condition, and annotation strings into one descriptive Sample Composition Summary node. Its complete
  Sample-by-population table and report can be previewed, and the table can be exported.
- `scVI Batch Integration and Contrast.json` trains scVI from raw counts in `adata.X`, constructs neighbors
  from `X_scVI`, visualizes the integrated embedding, and passes its native session-only model artifact to scVI
  differential expression instead of retraining it.
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
  records the declared `control` and `treated` Condition mapping while constructing the complete descriptive
  Sample-by-cell-type count/proportion grid, including structural zeros. It performs no Condition hypothesis test;
  choose a reviewed Sample-level compositional model when inferential comparison is required.
- Every consumer declares its expression source. Count-model methods enforce the numeric domain their algorithms
  require; QC, filtering, and normalization nodes that support broader finite expert inputs disclose signed or
  non-count-like values as warnings rather than imposing a universal count gate. For conventional practice, select
  a non-negative count representation. When counts live in `layers["counts"]`, select that layer directly on each
  count-dependent node. Each consumer uses the source selected on that node. The clustering template reads counts
  into `log1p_norm` without overwriting `X` or existing layers.
- The two bundled scVI templates explicitly set `technical_batch_key="batch"`, so their input data require
  `adata.obs["batch"]`; another workflow may declare a different Technical batch key. The downstream contrast
  requires the grouping column and labels shown on that node. Install the optional `scvi-tools` dependency in the
  same Python environment that runs ComfyUI before executing this template.
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

The scVI `model` output is a native session-only artifact for downstream nodes in the current Artifact Runtime
session, never a live Python object in ComfyUI. Saving the workflow records the connection and parameters, not the
native model directory; after restarting ComfyUI, rerun scVI Integration before executing its
differential-expression consumer.

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

The start scripts require a built `ComfyUI_frontend/dist/index.html`, then launch ComfyUI with:

```text
--disable-api-nodes --enable-assets --front-end-root <sibling frontend dist> --cache-classic
```

Additional arguments are forwarded to ComfyUI. The scripts add `--cache-classic` only when no cache mode was
supplied; `--cache-classic`, `--cache-none`, and `--cache-lru N` are supported. OpenBio rejects ComfyUI's default
RAM-pressure cache because `ArtifactTicket` values do not measure their retained disk results. ComfyUI's own
requirements must already be installed. In particular, `--enable-assets` requires ComfyUI's database dependencies
and a working local database connection.

## Data and output behavior

- AnnData values use `adata` ports and the `OPENBIO_ANNDATA` wire type. Large typed ports carry immutable
  `ArtifactTicket` values; the ComfyUI process does not load their AnnData, tables, plots, models, or other
  scientific state.
- Core Study Parameters exposes its six selections as standard ComfyUI `STRING` wires, so they can connect to any
  compatible string input.
- Tables, plots, and structured summaries use distinct `table`, `plot`, and `summary` ports with the
  plugin-owned `OPENBIO_SINGLE_CELL_TABLE`, `OPENBIO_SINGLE_CELL_PLOT`, and
  `OPENBIO_SINGLE_CELL_SUMMARY` wire types. Tables and plots are File artifacts; summaries, code, scalar values,
  and small metrics remain strict JSON-compatible values or strings.
- Refactored analysis and transformation nodes expose their primary result first, followed by a structured
  `summary` and a standard `STRING` `code` port. The summary wire carries an OpenBio result artifact whose
  `.summary` payload is strict JSON-compatible data with report-ready
  Methods and Results text, key results, resolved parameters, warnings and limitations, references, and runtime
  software versions; the preview/API payload adapter serializes that inner payload as JSON. The code port contains
  an equivalent Python function with the resolved scientific choices and scientific input checks; it intentionally
  omits plugin-specific `analysis_history` bookkeeping because the adjacent summary carries that provenance.
- Every scientific node invocation starts a fresh one-shot Worker, reads its immutable inputs into private state,
  computes in place, publishes new File artifacts, and exits. Its optional `worker` input accepts an
  `OPENBIO_WORKER` value from the Python Worker node, which selects one exact local Python executable; when
  unconnected it uses the ComfyUI interpreter. A selected interpreter is probed for the OpenBio protocol, Python
  3.12, and required dependencies. The path is trusted local workflow input, not a malicious-code sandbox.
- AnnData artifacts are complete uncompressed H5AD files. A one-shot Worker never opens an input artifact for
  writing and does not make a defensive full-object copy merely to protect upstream state. Algorithm-required
  workspaces, subsets, and the Raw snapshot may still allocate memory.
- AnnData transforms use `adata` as their primary output; analysis artifacts use their concrete table, plot, or
  summary type as the primary output; preview and save nodes are terminal and have no data output.
- Preview accepts all three artifact types, while CSV and PNG outputs accept only tables and plots respectively,
  so incompatible links are rejected before execution. CSV export always retains the row index under an explicit,
  collision-checked header; H5AD export streams the canonical uncompressed artifact without loading or recompressing it.
- Frequently tuned analysis choices, expression sources, and result-defining thresholds stay visible. Core Study
  Parameters is an optional source for reusing sample, condition, annotation, and primary-contrast strings; analysis
  nodes keep their ordinary explicit inputs and do not consume a combined design object. Random seeds, internal
  storage keys, output column names, and iteration limits are advanced inputs.
- PCA, Harmony, scVI, Neighbors, embedding, and clustering nodes keep representations and named graph bundles
  explicit. `n_dimensions=0` means all available dimensions after validation; this avoids requesting nonexistent PCs
  from short learned latent spaces. Harmony/scVI covariates are Technical batch or nuisance variables, not biological
  Sample replicates or Condition labels. Existing representation, graph, embedding, or clustering keys are protected
  by default and require an explicit overwrite choice.
- Force-Directed Graph defaults to Fruchterman-Reingold. PAGA and existing-coordinate initialization are explicit and
  validated; Fruchterman-Reingold and Kamada-Kawai use python-igraph, while ForceAtlas2 requires `fa2-modified` and
  never falls back silently. Every result reports the actual method and backend.
- PCA Metadata Associations aggregates scores within biological Samples before association testing. Its results are
  Sample-level diagnostics with one Benjamini-Hochberg adjustment across all valid PC-by-metadata hypotheses, not
  cell-level independent-replicate tests or proof that a Technical batch correction is appropriate.
- Leiden Resolution Sweep reports size, modularity, repeat-start stability, and adjacent-resolution agreement for
  every declared resolution. It does not select or label a biologically "best" resolution automatically.
- Marker Genes ranks every selected-source gene for each categorical cluster, recomputes within-cluster
  Benjamini-Hochberg adjustment over that complete tested family, and emits both a canonical marker table and the
  exact ordered tested-gene universe. Filter Marker Genes accepts only that direct pair; rerun filtering from the
  original ranking instead of chaining filters. These are exploratory Cluster marker results, not Sample-level
  Condition inference.
- UMAP Plot and Marker Expression Plot are read-only renderers. They report the actual embedding dimensions,
  expression source, missing-value handling, package versions, and plot-specific method; marker plots copy only
  selected genes and grouping metadata. Seeded violin jitter is isolated from NumPy's global random state.
- CellTypist Annotation produces provisional per-cell model labels and score evidence from an explicitly verified
  count or CP10K/log1p source. A selected model must already exist locally; execution does not enumerate or download
  the model catalog. Marker ORA Evidence consumes a direct filtered-marker table plus its pinned universe and an
  explicitly licensed resource CSV, but does not assign cell types. Map Cluster Annotations is the separate reviewed
  commitment step and records whether labels are provisional or curated.
- Dataset-specific condition, tumor, cluster, and cell-type values are never supplied as defaults; nodes that need them require an explicit value.
- Load H5AD and Load 10x H5 accept a browser file selection or a file dropped directly onto the node; uploads are
  stored under `ComfyUI/input/openbio-singlecell` and the node keeps the resulting relative path.
- Inputs are limited to relative paths under the ComfyUI `input` directory and are checked again at the read boundary.
- Temporary File artifacts are stored below `temp/openbio-singlecell`. Each node run publishes all of its File
  artifacts atomically under one `RunLease`; the shared run directory remains available while an `ArtifactTicket`
  is held by ComfyUI's cache or a running consumer and becomes eligible for cleanup after the final reference is
  released. UUIDs provide unique artifact identity; they are not content hashes or workflow-cache keys.
- Classic, None, and LRU cache modes retain or recompute node results according to ComfyUI's normal workflow
  signatures. A cache hit reuses the cached `ArtifactTicket` without starting a one-shot Worker. None retains no
  reusable entries, LRU releases evicted entries, and Classic may keep current intermediates until invalidation,
  cache reset, or process exit. If immediate cleanup is prevented, a later startup scavenges a stale session only
  after confirming that its owning process is dead.
- Independent DAG branches may read the same parent artifact concurrently and publish separate results. OpenBio
  uses ComfyUI's native concurrency without a memory scheduler, so concurrent one-shot Workers' AnnData and
  algorithm workspace memory use adds together.
- CSV, PNG, and H5AD files are only made permanent by explicit output nodes and are written below
  `output/openbio-singlecell`. They are staged and validated in the destination directory before an atomic commit;
  a failed write leaves an existing destination unchanged. PNG previews use the same validation/commit discipline
  below the temporary directory.
- Persist Artifact copies a portable File artifact to `output/openbio-singlecell/artifacts` as a Persisted artifact
  using staged validation, a hash manifest, and atomic publication. It neither moves the temporary source nor
  returns a workflow value. There is no generic Restore node; an existing H5AD can re-enter only through the normal
  H5AD input path after it is placed under ComfyUI's input directory.
- Node UI payloads are bounded; complete tables are exported as CSV.
- Tables use JSONL plus an explicit schema, plots use PNG, and specialized typed results use validated file codecs.
  Arbitrary Python objects are never carried by `ArtifactTicket` values. The registry exposes only the
  artifact-backed scientific nodes; there is no live-object compatibility mode.
- Completed cNMF rank surveys, trained scVI models, validated external pySCENIC result bundles, and solved
  Cassiopeia trees have explicit typed contracts because each has a real downstream consumer. No generic
  Python-object or generic model wire is exposed. cNMF and scVI results that use their libraries' native directories
  are native session-only artifacts: they are bound to the originating Artifact Runtime session and the producer's
  probed Python interpreter identity. Each consumer still starts a fresh one-shot Worker with that identity. These
  artifacts cannot become Persisted artifacts, survive restart, be restored, or be consumed by a different
  interpreter. Imported pySCENIC evidence is immutable and bound to its external-run manifest.
- A transformation keeps one `adata` output when it has no reusable secondary product. The deliberate exceptions
  include scVI Integration (`adata`, `model`) and Import pySCENIC Results (`adata`, `scenic_result`); cNMF Rank Survey
  produces reusable `run` state plus a K-metrics table, and Cassiopeia reconstruction produces one reusable `tree`
  for its expansion and plasticity consumers.

## Scope and limitations

The current workflow includes 10x study loading, QC, expression snapshots, layer-aware preprocessing, Harmony/scVI
integration, clustering, annotation, pseudobulk differential analysis, enrichment, LIANA communication,
compositional abundance testing, method-author cNMF, Schist, external pySCENIC result import, Cassiopeia lineage
analysis, PAGA/DPT, RNA velocity, and inferCNV. The pack does not declare a maximum AnnData size and does not provide
general out-of-core computation, ATAC, or spatial analysis. File artifacts bound retained ComfyUI-process memory,
but each one-shot Worker still materializes one private scientific state plus any algorithm workspace.
Notebook stages implemented only in R, destructive file organization, publication-only plots, the empty Geneformer
notebook, and the import-only spatial notebook are intentionally not represented as nodes. Scale and GSVA may
densify data as their underlying libraries normally do. A running scientific operation may finish before a stop
request takes effect; stopping prevents later nodes from starting.

## Offline behavior

The plugin contains no telemetry, analytics, update checks, or remote configuration. The OpenBio launcher disables
Comfy API nodes. Node execution never performs an implicit model or knowledge-resource download. CellTypist requires
an explicitly preinstalled local model. The DGIdb loader accepts only a user-provided, versioned local CSV/TSV
resource and has no network path. CollecTRI defaults to a user-provided local network; it accesses decoupler's
official remote resource only when the user selects the official-resource mode and separately enables its visible
`allow_network_access` control. Prepare every selected resource in advance and keep that control disabled for an
offline run. The installer may contact the Python package index configured for `pip`.

## License and attribution

`openbio-singlecell` is licensed under GPL-3.0-or-later. The complete license text is included in `LICENSE`. It is based on ComfyUI and its paired ComfyUI_frontend, both distributed under GPLv3 terms. Keep upstream copyright notices, source attribution, and the corresponding GPL license text with redistributed source or binaries. See `THIRD_PARTY_NOTICES.md` for source-package notices; paired prebuilt frontend releases also carry their generated production-dependency notice inside `dist`.
