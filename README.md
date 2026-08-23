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

It contains 600 cells, 500 genes, three known cell groups, mitochondrial genes, `sample` and `batch` metadata, known marker genes, and CSR sparse integer counts. A valid existing file is left unchanged; replacing an invalid file requires the explicit force option.

## Example workflows

Open either workflow from the repository's `example_workflows` directory:

- `openbio_singlecell_basic_qc.json` loads the demo AnnData, calculates and filters QC metrics, previews a structural summary and QC plots, and writes H5AD only through the explicit save node.
- `openbio_singlecell_full_analysis.json` adds normalization, log transformation, highly variable genes, scaling, PCA, neighbors, UMAP, Leiden clustering, marker genes, previews, and explicit CSV, PNG, and H5AD outputs.

`example_workflows` is the only packaged workflow source. There is no mirrored copy under `web`, so examples cannot drift from the files reviewed and released here.

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
- Single-cell result values use `result` ports and the plugin-owned `OPENBIO_SINGLE_CELL_RESULT` wire type.
- Modifying nodes copy their input before changing it; read-only nodes do not copy it.
- AnnData transforms output only `adata`; tables and plots output only `result`; preview and save nodes are terminal and have no data output.
- Frequently tuned analysis choices stay visible. Random seeds, internal storage keys, output column names, and iteration limits are advanced inputs.
- Inputs are limited to relative paths under the ComfyUI `input` directory and are checked again at the read boundary.
- Temporary plots are written below `temp/openbio-singlecell`; ComfyUI clears its temp directory at startup.
- CSV, PNG, and H5AD files are only made permanent by explicit output nodes and are written below `output/openbio-singlecell`.
- Node UI payloads are bounded; complete tables are exported as CSV.

## Scope and limitations

The current in-memory workflow includes layer-aware preprocessing, batch-aware HVG selection, Harmony integration, representation-aware neighbor graphs, annotation, and PAGA/DPT trajectory tools. Harmony, CellTypist, and marker ORA use the corresponding optional packages from the active ComfyUI Python environment. The pack does not declare a maximum AnnData size and does not yet provide AnnData backed mode, out-of-core processing, automatic disk caching, scVI, CellChat, RNA velocity, ATAC, or spatial analysis. Scale may densify data as Scanpy normally does. A running Scanpy/igraph operation may finish before a stop request takes effect; stopping prevents later nodes from starting.

## Offline behavior

Analysis and demo generation are local-only. The plugin contains no telemetry, analytics, update checks, model downloads, remote configuration, or other runtime network request path. The OpenBio launcher disables Comfy API nodes. The installer may contact the Python package index configured for `pip`; install dependencies before disconnecting for a fully offline run.

## License and attribution

`openbio-singlecell` is licensed under GPL-3.0-or-later. The complete license text is included in `LICENSE`. It is based on ComfyUI and its paired ComfyUI_frontend, both distributed under GPLv3 terms. Keep upstream copyright notices, source attribution, and the corresponding GPL license text with redistributed source or binaries. See `THIRD_PARTY_NOTICES.md` for source-package notices; paired prebuilt frontend releases also carry their generated production-dependency notice inside `dist`.
