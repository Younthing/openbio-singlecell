# Native R H5AD interoperability probe — throwaway evidence

Question: can the existing project H5AD artifact move through a native R whole-object read/write without losing unrelated scientific state?

**Verdict: not as an unrestricted whole-object transport in the tested stack.** Native R matrix I/O worked on actual files, but successful rewrites lost `raw`, altered dtypes, and encountered metadata encoding and project codec incompatibilities. This probe makes no production serializer changes.

Tested Python: 3.13.5, anndata 0.13.2, pandas 3.0.5, NumPy 2.4.4. Tested R: 4.5.3, anndataR 1.0.2 (Bioconductor 3.22), rhdf5 2.54.1. Reticulate was an imported dependency only: Python initialization remained `FALSE` before and after both successful native read/write operations. The successful runs used reticulate 1.46.0; a later conda replacement with 1.47.0 did not resolve the Windows startup failure described below.

The input uses the project's `write_anndata` codec: 4 cells × 3 genes, float32 CSR X, int32 CSR counts, ordered categorical metadata with an unused level and a missing value, boolean/int64/nullable Int64 metadata, raw with a separate 5-feature axis, a float32 embedding and CSR graph, float64 loadings and CSC feature graph, and nested JSON `uns` with mapping-order preservation.

| State | pandas 3 default strings | Explicit object-dtype strings |
|---|---|---|
| Matrix numbers, shape, CSR/CSC orientation | Preserved in X, counts, obsm, obsp, varm, varp | Preserved |
| float32 X/embedding/graph and int32 counts | Rewritten as float64 | Rewritten as float64 |
| Cell/gene ID values and metadata columns | Lost: reader warns on `nullable-string-array`, returns empty obs/var, and writes numeric IDs | IDs, column values and dtypes preserved |
| Ordered categories, unused level, missing entry | Lost with the affected frame | Preserved |
| Named indexes (`cell_id`, `gene_id`) | Lost | Lost; ID values still preserved |
| `raw` and its independent feature axis | Completely absent from rewritten H5AD | Completely absent |
| Existing nested `uns`, decoded by project codec | Preserved when R changes no keys | Preserved when R changes no keys |
| Add `uns$r_probe` in R, retain opaque codec manifest | Native Python reads it; project codec rejects stale mapping-order keys | Same rejection |

The metadata loss is attributable to anndataR 1.0.2's unsupported nullable string encoding; it is not a claim that all R H5AD readers lose metadata. The exact warning was `No function for reading H5AD encoding <nullable-string-array> for element "obs/cell_id"`. Matrix numeric comparisons are exact for this small fixture; this does not certify all possible numeric magnitudes or object types.

The object-string case was an additional successful controlled run: the R console showed a four-column obs data.frame with an ordered factor, logical column, integer column and nullable integer column. Its later rerun files were cleared by the probe before the Windows process failed, so the permanently preserved binary evidence is the first default-string run. The controlled-case result above records that observed run; it is not presented as a successful latest rerun.

The project manifest `__openbio_json_uns_v1__` stores both JSON payload paths and mapping key order. An unchanged R rewrite preserves its opaque strings, allowing `read_anndata` to reconstruct nested values. Adding a new `uns` key without regenerating the manifest causes `ValueError: AnnData nested-uns mapping-order keys do not match the stored mapping.` Plain `anndata.read_h5ad` understands the H5AD structure but does not decode this project-specific contract.

## Runtime result and evidence

Initial native-R runs succeeded, but repeated fresh processes intermittently failed while loading the mandatory reticulate dependency with `Mingw-w64 runtime failure: 32 bit pseudo relocation ... out of range`. The latest attempts failed with process code `3221226505` before analysis. Installing the conda-native reticulate build did not cure this. rhdf5 and Matrix loaded independently. Setup was stopped without changing production dependencies, base R, Matrix, or global configuration.

- `h5ad_input/data.h5ad`, `h5ad_roundtrip/data.h5ad`, `h5ad_r_modified/data.h5ad`: preserved input and two real outputs from the successful default-string run.
- `h5ad_run.log`: successful native R output, metadata warnings, no raw binding, and Python-uninitialized checks.
- `h5ad_initial_report.json`: regenerated exact comparisons of those preserved files through both native Python and project codec readers.
- `h5ad_report.json` and `h5ad_*_run.log`: latest rerun results, including actual startup failures; these must not be cited as passes.
- `h5ad_setup.R`, `h5ad_setup.log`, `h5ad_conda_reticulate.log`: isolated installation attempt and version evidence.

Read the preserved successful outputs without starting R:

```powershell
& 'D:\learn\ComfyUI\.venv\Scripts\python.exe' '.scratch\r-worker-prototype\h5ad_probe.py' --compare-archived
```

Attempt both fresh native R cases (returns nonzero if either R process fails):

```powershell
& 'D:\learn\ComfyUI\.venv\Scripts\python.exe' '.scratch\r-worker-prototype\h5ad_probe.py'
```

For the first R analysis nodes, the evidence supports Python retaining the full AnnData artifact and applying only declared R outputs returned through explicit matrix/table/JSON files. That preserves untouched `raw`, metadata dtypes, provenance, and the project codec manifest. A native whole-object R H5AD path needs a pinned, reliable runtime and per-node encoding/slot fidelity evidence before it can claim to preserve a whole artifact. No Seurat conversion or large-data scaling was tested.

Primary references: [anndataR native read interface](https://anndatar.scverse.org/reference/read_h5ad.html), [native write interface](https://anndatar.scverse.org/reference/write_h5ad.html), [anndataR development status](https://anndatar.scverse.org/articles/development_status.html) (raw support remains blank in the displayed 1.3.1 matrix), [tested release DESCRIPTION](https://raw.githubusercontent.com/scverse/anndataR/RELEASE_3_22/DESCRIPTION). The current website documents a newer development version than the tested Bioconductor/R-compatible binary; measured losses above are scoped to the tested versions.
