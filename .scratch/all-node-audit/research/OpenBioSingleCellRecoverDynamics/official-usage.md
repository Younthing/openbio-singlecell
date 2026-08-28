# Recover Velocity Dynamics: official usage and scientific practice

## Official API and outputs

For scVelo 0.3.4:

```python
scvelo.tl.recover_dynamics(
    data, var_names="velocity_genes", n_top_genes=None, max_iter=10,
    assignment_mode="projection", t_max=None, fit_time=True,
    fit_scaling=True, fit_steady_states=True, fit_connected_states=None,
    fit_basal_transcription=None, use_raw=False, load_pars=None,
    return_model=None, plot_results=False, steady_state_prior=None,
    add_key="fit", copy=False, n_jobs=None, backend="loky",
    show_progress_bar=True, **kwargs,
)
```

The EM model recovers transcription, splicing and degradation parameters plus gene/cell-specific dynamical times and states. It writes many `fit_*` fields/layers and `uns["recover_dynamics"]`. Official scVelo workflow runs this before `velocity(..., mode="dynamical")`.

Sources:

- [scVelo `recover_dynamics`](https://scvelo.readthedocs.io/en/stable/scvelo.tl.recover_dynamics.html)
- [scVelo dynamical modeling tutorial](https://scvelo.readthedocs.io/en/latest/DynamicalModeling.html)
- [scVelo 0.3.4 release](https://github.com/theislab/scvelo/releases/tag/v0.3.4)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)

## Selection, prerequisites and latent time

The default `var_names="velocity_genes"` can trigger an internal steady-state velocity fit when no velocity-gene mask exists. That is a second analysis hidden inside recovery. A production wrapper must either require an upstream verified mask or explicitly request all genes/a bounded subset; it must not silently invent the selection.

The function needs canonical spliced/unspliced abundances, moments/connectivities and enough informative genes. It is expensive and writes dense cell-by-gene fitted layers, so gene count, iteration count, jobs and memory must be bounded and reported. Existing fit bundles must be collision-guarded. In 0.3.4, `n_top_genes` orders candidates by descending summed `Ms` (or raw spliced data when `use_raw=True`); a wrapper that promises an exact selected-gene identity should reproduce that stable selection before the call and pass the resulting explicit gene list with backend `n_top_genes=None`.

The audited public implementation has two incompatibilities with the repository's resolved Pandas 3.0.5. It converts explicit `var_names` to a Python list and passes that list from `make_unique_list` to `pandas.unique`, which no longer accepts bare lists. Its `_read_pars` also takes `Series.values` views and `align_dynamics` writes into them, while Pandas 3 exposes read-only arrays. These are upstream software-compatibility defects rather than scientific reasons to change genes or parameters. A wrapper may install narrowly scoped, lock-protected replacements for those two callable globals: convert only list/tuple inputs to an object ndarray before stable `pandas.unique`, and return writable copies of parameter arrays. It must verify both original helper signatures and restore them in `finally`. The exact public `scv.tl.recover_dynamics` signature and explicit gene family remain unchanged, and the summary must disclose both actions plus scVelo/Pandas versions.

`recover_dynamics` does **not** compute a gene-shared `latent_time`. Official `scv.tl.latent_time` is a later operation after dynamical velocity and its graph; it may infer root/end states if none are supplied. The repository currently has no latent-time node. That capability must be a separately researched node with explicit root/end policy, not a claim added to this node's summary.

Recovered parameters are model fits, not directly measured kinetic rates or fate probabilities. The 0.3.4 model computes `fit_likelihood` as the exponential of a normal-noise log-likelihood/density score. Successful fits must have finite nonnegative values, but the score is not a probability and is not mathematically bounded above by 1. It ranks fit quality under the model and is not a differential-expression p-value.
