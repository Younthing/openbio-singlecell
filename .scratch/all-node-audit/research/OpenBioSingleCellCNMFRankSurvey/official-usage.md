# cNMF Rank Survey: official usage research

Researched: 2026-08-28

## Scope and authoritative implementation

`OpenBioSingleCellCNMFRankSurvey` implements only the expensive upstream half of a staged cNMF analysis: validate a declared non-negative expression matrix from current `X`, one named layer, or an explicitly selected `adata.raw`; prepare the official file-backed run; execute every stochastic restart for every candidate rank; combine the complete restart family; and report one stability/reconstruction-error record per rank. It does not select a rank automatically, apply the final density filter, or interpret programs. The report separately states whether local provenance supports the method's usual UMI-count interpretation; numeric values alone never prove that provenance.

The supported backend is the method authors' package, exactly `cnmf==1.7.1`, imported as `cnmf`. OmicVerse is not an implementation dependency for this node.

Primary sources:

- method-author repository and stepwise guide: https://github.com/dylkot/cNMF and https://github.com/dylkot/cNMF/blob/master/Stepwise_Guide.md
- versioned PyPI release: https://pypi.org/project/cnmf/1.7.1/
- cNMF method paper: Kotliar D, Veres A, Nagy MA, et al. *eLife*. 2019;8:e43803. https://doi.org/10.7554/eLife.43803
- AnnData file format/API: https://anndata.readthedocs.io/

PyPI records 1.7.1 as released on 2026-04-11 under the MIT license. The pure-Python wheel is `cnmf-1.7.1-py3-none-any.whl`, SHA-256 `a4b237169626f632fee04a31441365669b12a02c72a93b079686b0ce47a29ec9`. PyPI reports that upload was not made through Trusted Publishing, and the distribution metadata does not establish a source-commit mapping. The adapter therefore pins and reports the package version and distribution hash evidence; it must not claim a commit pin.

## Scientific interpretation

cNMF fits non-negative matrix factorization repeatedly at each K and clusters the resulting components. The rank survey exposes:

- `stability`: the silhouette of replicate spectra assigned to K consensus clusters;
- `prediction_error`: the Frobenius reconstruction residual on the normalized count matrix.

These are stochastic NMF restarts on one prepared matrix, not bootstrap or biological replicates. K remains an analyst decision supported by the stability/error trade-off and biological interpretability. Prediction error is scale- and dataset-dependent and is only comparable within the same run. Official examples commonly use many restarts and the node defaults to 100. The adapter requires only two restarts: cNMF's reported silhouette statistic requires at least two clusters and fewer clusters than replicate-component observations, which fails for the one-restart family. Fewer than 100 restarts is disclosed as weak stability evidence rather than blocked.

## Exact public API used

For `cnmf==1.7.1`, the audited public signatures are:

```python
cNMF(output_dir=".", name=None)

prepare(
    counts_fn,
    components,
    n_iter=100,
    densify=False,
    tpm_fn=None,
    seed=None,
    beta_loss="frobenius",
    num_highvar_genes=2000,
    genes_file=None,
    alpha_usage=0.0,
    alpha_spectra=0.0,
    init="random",
    max_NMF_iter=1000,
)
factorize(worker_i=0, total_workers=1, skip_completed_runs=False)
combine(components=None, skip_missing_files=False)
consensus(
    k,
    density_threshold=0.5,
    local_neighborhood_size=0.30,
    show_clustering=True,
    skip_density_and_return_after_stats=False,
    close_clustergram_fig=False,
    refit_usage=True,
    normalize_tpm_spectra=False,
    norm_counts=None,
    build_ref=True,
)
```

The node guards the installed distribution version and the named parameters before execution. Accepting a callable with only `**kwargs` is insufficient: the underlying function must expose the audited controls.

## Canonical file-backed chain

The equivalent official chain is:

```python
from cnmf import cNMF

runner = cNMF(output_dir=private_root, name=private_name)
runner.prepare(
    counts_fn=private_counts_h5ad,
    components=components,
    n_iter=n_iter,
    densify=False,
    tpm_fn=None,
    seed=random_seed,
    beta_loss="frobenius",
    num_highvar_genes=num_highvar_genes,
    genes_file=None,
    alpha_usage=0.0,
    alpha_spectra=0.0,
    init="random",
    max_NMF_iter=1000,
)
runner.factorize(worker_i=0, total_workers=1, skip_completed_runs=False)
runner.combine(components=components, skip_missing_files=False)
```

For every K, statistics mode is then called with the normalized counts loaded from the run file:

```python
stats = runner.consensus(
    k=K,
    density_threshold=2.0,
    local_neighborhood_size=0.30,
    show_clustering=False,
    skip_density_and_return_after_stats=True,
    close_clustergram_fig=True,
    refit_usage=True,
    normalize_tpm_spectra=False,
    norm_counts=normalized_counts,
    build_ref=False,
)
```

The returned `stats` DataFrame must have the four exact rows `k`, `local_density_threshold`, `silhouette`, and `prediction_error`. The node emits one ordered table row per surveyed K and checks that K, silhouette, error, and threshold are finite and in their scientific domains. Statistics mode does not create final usages or gene spectra.

## Managed files and completeness evidence

The official package is file-backed. The node creates one unique private `TemporaryDirectory`, a private run name, and a private H5AD containing only the selected count matrix plus exact observation/variable identifiers. It reopens that H5AD and verifies shape, axes, and a matrix fingerprint before preparing cNMF.

The owned run directory contains the official artifacts addressed through `runner.paths`, including:

- normalized counts H5AD and TPM/full-gene files;
- replicate parameters/seeds and run parameters;
- one spectrum NPZ for every `(K, restart)` job;
- one merged-spectrum NPZ for every K;
- density caches and final consensus files created later.

After factorization, every expected restart file is present, regular, contained within the private root, and has the expected K-by-HVG axes. After combination, each merged spectrum has exactly `n_iter * K` component rows and the same unique HVG columns. The normalized-count file is reopened and checked against the exact input cell axis. Missing jobs are never tolerated: `skip_completed_runs=False` and `skip_missing_files=False` are fixed.

Paths obtained from the external backend are treated as untrusted. Before every read or targeted cleanup, the adapter resolves the path, requires strict containment under the owned root, rejects `..`, absolute escape, symlink/reparse traversal, non-regular files, and file replacement during a read. No working-directory or name input is exposed to users.

## Count and resource policy

Official UMI practice recommends filtered, untransformed integer counts. Runtime hard gates are limited to backend-computability and alignment: the selected matrix must be numeric, finite, non-negative, have nonempty unique axes, and contain no zero-total cell. A zero-total gene is disclosed but remains calculable—the official backend can exclude it from HVGs and retain a zero full-gene result. A non-integer or provenance-confirmed transformed-but-still-non-negative source likewise remains available to an expert, with a prominent summary warning that the default UMI interpretation may not apply. Negative scaled/residual matrices still fail because non-negative matrix factorization cannot calculate them. Sparse input remains sparse (`densify=False`).

Raw is never substituted implicitly, but `source=raw` is a supported expert choice. The adapter uses `adata.raw.X`, `adata.raw.obs_names`, and the complete `adata.raw.var_names` exactly as that source's own analysis axes. Raw may have more or different features than current `adata.var_names`; no axis intersection, Raw/current equality, value comparison, processing-history proof, or historical-integrity claim is required. Raw must exist and have nonempty unique identifiers. For a Raw run, the owned base/result AnnData is materialized on Raw's own complete observation and feature axes so final usages and gene-program matrices remain aligned to the selected source rather than being forced into current axes.

Candidate K is an inclusive ascending range. Each K is at least two because the reported silhouette requires at least two clusters. `n_iter` is at least two because a one-restart family has exactly K replicate components for K clusters and cannot satisfy the silhouette sample/label domain. The adapter imposes no arbitrary K, restart, or HVG upper bound, does not require K to be no greater than the cell count, and does not require K to be no greater than the realized HVG or positive-variance feature counts: the official random-initialized factorization accepts overcomplete shapes. Those regimes are scientifically weakly identified and potentially expensive, so they are disclosed for expert judgment. A requested HVG count above the available genes is likewise disclosed and the official backend may retain at most eligible genes. Resource preflight reports the candidate family, expected restart jobs, requested and realized matrix/spectra estimates, safety factor, and a 2 GiB reference envelope; exceeding that conservative estimate produces an expert-facing warning rather than a runtime gate.

The fixed scientific policy is Frobenius loss, random initialization, no regularization, no supplied TPM or gene file, `max_NMF_iter=1000`, and one complete worker. Worker partitioning and filesystem naming are hidden implementation details.

## Typed run lifetime

The returned `OPENBIO_CNMF_RUN` owns its `TemporaryDirectory`; the directory remains present while the run is live and may be reused for multiple downstream threshold/K decisions. Survey metadata and K metrics are exposed as read-only properties. At construction the run records a canonical snapshot fingerprint over their complete typed contents; liveness checks and Consensus recompute that fingerprint, so rebinding or even low-level mutation of source provenance, advisories, metrics, resource accounting, or artifact hashes is rejected before scientific use.

Equivalent Survey and Consensus code may be executed in separate Python namespaces, so generated Consensus accepts a portable cNMF run only through an explicit class-owned artifact marker/schema plus the same complete snapshot, concrete field-type, lifecycle-owner, root, and backend checks. A class with merely similarly named methods, or an instance-provided marker, is not an artifact. This is a strict process-local type/interface check, not an authenticity claim against arbitrary hostile Python code.

In equivalent Python, `run.close()` and context-manager exit are explicit, idempotent cleanup paths and make all subsequent access fail. In a node graph there is deliberately no close node or auto-close in Consensus, because either would invalidate sibling Consensus consumers of the staged result; cleanup occurs when the engine releases its cached run object, at which point the run's finalizer removes the directory. The summary discloses that a cached graph can therefore retain the private directory until cache eviction or process exit. Preparation or postcondition failure cleans the directory immediately rather than returning a partial artifact. Cleanup is secondary to an already-active analysis, context-body, or report exception: a cleanup failure is attached to that primary exception as a note instead of replacing it. If explicit `close()` itself is the primary operation and cleanup raises, the cleanup error remains outward; the run is closed but retains its cleanup owner and live finalizer so `close()` or finalization can retry. The run is process-local and intentionally non-serializable.

## Known upstream warnings and compatibility boundary

On Python 3.12 with NumPy 2, AnnData 0.13, Scanpy 1.12, and pandas 3, real 1.7.1 prepare/factorize/combine/statistics succeeds. The package currently emits two known upstream warnings:

- `ResourceWarning` for the unclosed `*.nmf_idvrun_params.yaml` reader in `cnmf/cnmf.py` during factorize/combine;
- `pandas.errors.Pandas4Warning`: `Starting with pandas version 4.0 all arguments of sum will be keyword-only.` from the statistics path in `cnmf/cnmf.py`.

The adapter captures only these exact category/message/source combinations, records them in the report, and re-emits every unknown warning. It does not install a global warning filter. The ResourceWarning is an upstream file-handle defect; the adapter does not broaden it into permission to ignore arbitrary resource leaks.

## Required report and generated source

The strict-JSON summary reports the exact declared source/fingerprint, `source_features` (the selected source's feature-axis size), current feature count for context, source state/evidence, a structured verbatim `input_advisories` list, dimensions, total counts, candidate Ks, expected/completed jobs, seed, requested/realized HVGs, all K metrics, fixed policy, resource envelope, owned-file mode, package versions, known warnings, references, and limitations. For Raw it reports `source="raw"` and never implies Raw/current feature or value equality. Unknown provenance is described without asserting that the values are integer-valued or UMI counts; pooled-cell programs can reflect Sample or Technical batch structure and are not replicate-aware Condition inference.

The `code` output is standalone source using public `cnmf==1.7.1`, a private `TemporaryDirectory`, the same count/path/completeness checks, exact explicit parameters, and a returned owned run object with `close()`/finalizer semantics. It neither imports `openbio_singlecell` nor uses a shared ComfyUI path.
