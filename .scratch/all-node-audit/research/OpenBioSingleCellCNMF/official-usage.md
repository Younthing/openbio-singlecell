# cNMF Consensus Programs: official usage research

Researched: 2026-08-28

## Scope

`OpenBioSingleCellCNMF` is the downstream atomic operation in the staged workflow. It consumes one live, integrity-checked `OPENBIO_CNMF_RUN`, applies one analyst-selected surveyed K and local-density threshold, computes final consensus programs through the method authors' exact `cnmf==1.7.1` package, loads all official result matrices, and annotates a copy of the original AnnData. It never repeats preparation/factorization and never accepts arbitrary AnnData or a filesystem path in place of the typed run. The SHA/path checks detect replacement relative to the Survey snapshot; they are not a cryptographic authenticity or hostile-writer guarantee.

Authoritative sources:

- method-author repository/stepwise guide: https://github.com/dylkot/cNMF and https://github.com/dylkot/cNMF/blob/master/Stepwise_Guide.md
- versioned release: https://pypi.org/project/cnmf/1.7.1/
- Kotliar D, Veres A, Nagy MA, et al. *eLife*. 2019;8:e43803. https://doi.org/10.7554/eLife.43803
- AnnData: https://anndata.readthedocs.io/

The package is MIT licensed. Version 1.7.1 is pinned; its PyPI distribution is not presented as a source-commit pin or Trusted Publishing attestation.

## Official final-consensus calls

The supported final call is explicit:

```python
runner.consensus(
    k=selected_k,
    density_threshold=density_threshold,
    local_neighborhood_size=local_neighborhood_size,
    show_clustering=False,
    skip_density_and_return_after_stats=False,
    close_clustergram_fig=True,
    refit_usage=True,
    normalize_tpm_spectra=False,
    norm_counts=normalized_counts,
    build_ref=False,
)

usage, spectra_scores, spectra_tpm, top_genes = runner.load_results(
    K=selected_k,
    density_threshold=density_threshold,
    n_top_genes=n_top_genes,
    norm_usage=True,
)
```

In 1.7.1, `load_results` returns exactly the tuple `(usage, spectra_scores, spectra_tpm, top_genes)`. Code must not assume an OmicVerse `consensus_results` dictionary or any in-memory `iter_spectra_dict`, `merged_spectra_dict`, or `norm_counts` attribute.

`build_ref=False` avoids the optional reference-model branch, which is outside this node's scientific contract. `show_clustering=False` and `close_clustergram_fig=True` avoid an unmanaged plot side effect. Program usages remain continuous mixtures; no RFC/hard cluster is produced.

## Density filtering and independent verification

The upstream merged spectrum contains `n_iter * K` restart components. cNMF L2-normalizes component spectra, computes Euclidean distances to a neighborhood of size `int(local_neighborhood_size * n_iter)`, defines each component's local density from those neighbors, retains components whose density is below the selected threshold, and clusters the retained components into K groups with scikit-learn KMeans (`n_init=10`, `random_state=1`).

The node independently reconstructs this narrow diagnostic from the official merged-spectrum file using the audited 1.7.1 algorithm. It verifies:

- local density aligns exactly to the complete merged-component index and is finite/non-negative;
- the retained index derived from the requested threshold matches the components eligible for clustering;
- enough components remain for K nonempty groups;
- independently derived labels are integer-like and cover exactly `1..K`;
- official result artifacts exist only at canonical contained paths and correspond to the same K/threshold.

This verifier is a postcondition, not a replacement backend. The official public `consensus` and `load_results` remain the source of final usages and spectra. Before a repeated call on one run, the adapter removes only the exact contained density cache for that K so a different `local_neighborhood_size` cannot reuse the upstream cache whose filename does not encode that fraction.

## Result axes and AnnData contract

The four loaded results are validated by identifiers rather than position.  The
feature axis is the source axis captured by Survey, not an assumed alias for the
current AnnData axis:

- `usage`: cells by K programs, exact original `obs_names`, finite, non-negative, each positive row normalized approximately to one;
- `spectra_scores`: source genes by K programs after the official return orientation, exact Survey source feature identifiers, finite;
- `spectra_tpm`: source genes by K programs, exact Survey source feature identifiers, finite and non-negative;
- `top_genes`: exactly K unique program lists, identifiers drawn from the Survey source feature axis, with bounded unique entries.

For an `X` or named-layer Survey, the source feature axis is the input AnnData's
current `var_names`; Consensus therefore annotates a copy on that same current
axis.  For an explicit `Raw` Survey, Survey materializes `adata.raw` as its own
AnnData and the run retains both the complete Raw source-feature count and the
original current-feature count.  Consensus returns the Raw-materialized AnnData
on the complete Raw `var_names` axis; it does not require Raw and current axes or
values to be equal, and it makes no claim about Raw history or cryptographic
integrity.

Program names are canonicalized to `cNMF_1` through `cNMF_K` consistently across all matrices. On a transactional copy of that source-aligned AnnData the node stores:

- `obsm["X_cnmf_usage"]`: aligned cells-by-program DataFrame;
- `varm["cnmf_gep_scores"]`: aligned genes-by-program DataFrame;
- `varm["cnmf_gep_tpm"]`: aligned genes-by-program DataFrame;
- `uns["cnmf"]`: full survey decision context, chosen consensus controls/diagnostics, top genes, storage keys, versions, fingerprints, warnings, and limitations.

Input X, layers, Raw snapshot, identifiers, and annotations are not mutated. Existing canonical result keys cause failure unless `overwrite_existing=True`; failure leaves the input and prior result state untouched.

## Run lifetime and concurrency

Consensus requires either this module's concrete live `CNMFRun` or the explicitly marked portable `CNMFRun` produced by equivalent Survey code. It does not accept a generic object merely because it has `metadata`, `metrics`, and similarly named methods. Before copying data or consulting a metric, it checks the class-owned artifact marker/schema, complete construction snapshot, concrete metadata/metric field types and values, lifecycle owner, root identity, and exact backend. Rebinding or low-level mutation of source state/evidence/advisories, stability/error values, or artifact manifests therefore fails before output mutation. These checks establish process-local domain equivalence; they are not hostile-code authentication.

Consensus executes under the artifact lock because upstream density/result caches are mutable. It verifies the root and every official path before and after reads/calls. In equivalent Python, `run.close()` or context-manager exit remains the caller's explicit cleanup path after the last consensus attempt. In the node graph, the engine releases the cached typed run after it is no longer referenced and its finalizer performs cleanup; Consensus cannot auto-close without breaking sibling K/threshold consumers. The summary therefore discloses that the private directory can remain until graph-cache release or process exit. Closed, finalized, substituted, escaped, symlinked, or integrity-mismatched runs are rejected.

The run is process-local. Saving only the annotated AnnData is supported; serializing or reconnecting the live run across process restarts is not.

## Reporting requirements

The strict-JSON summary includes:

- the complete ordered K-metric family and all Survey parameters/job accounting/resource estimates;
- direct structured `source`, `source_features`, and `current_features` fields, so an `X`/layer result is visibly current-axis aligned and a Raw result is visibly Raw-axis aligned without inferring this from nested provenance;
- the Survey source state, provenance evidence, and every expert input advisory without loss at the staged boundary;
- selected K explicitly identified as an analyst choice;
- threshold/neighborhood size, total/retained/filtered component counts and fraction;
- selected-K stability and prediction error;
- per-program usage distribution and bounded leading genes;
- exact storage keys, overwrite behavior, fingerprints, package versions, known upstream warnings, references, and limitations.

Report prose must state that programs require biological annotation and external validation; pooled-cell patterns may reflect Sample or Technical batch effects; continuous usages are exploratory program evidence, not replicate-aware Condition inference; non-negativity cannot directly represent repression; and low nonzero usages can reflect overfitting.

## Known upstream boundary

The same exact 1.7.1 warning adapter applies. Real final-consensus smoke additionally established an unclosed-reader `ResourceWarning` for `*.overdispersed_genes.txt` at `cnmf/cnmf.py:964` and an AnnData `ImplicitModificationWarning` for upstream view materialization at line 969. The surveyed YAML-reader `ResourceWarning` and pandas keyword-only `Pandas4Warning` may also be collected when delayed upstream objects are finalized. Only the exact audited category, stable message, source file, and line combinations are recorded and consumed. Unknown warnings are re-emitted and remain fatal under `-W error`.

## Equivalent source

The `code` output is standalone, imports public `cnmf==1.7.1`, consumes the owned run object returned by the generated Survey function, validates the same official files and source axis, performs the same cache isolation/independent density-label postcondition, calls official consensus/load_results, and writes the same canonical AnnData state. The generated typed run carries the Survey-time current-feature count in addition to the selected source-feature count, so final `source`, `source_features`, and `current_features` disclosure is runtime/generated equivalent for `X`, layer, and Raw. It neither imports `openbio_singlecell` nor changes temporary ownership. Runtime and generated outputs have full scientific-state parity.
