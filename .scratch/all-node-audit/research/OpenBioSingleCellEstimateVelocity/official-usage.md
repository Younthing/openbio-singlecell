# Estimate RNA Velocity: official usage and scientific practice

## Official API and model prerequisites

For scVelo 0.3.4:

```python
scvelo.tl.velocity(
    data, vkey="velocity", mode="stochastic", fit_offset=False,
    fit_offset2=False, filter_genes=False, groups=None, groupby=None,
    groups_for_fit=None, constrain_ratio=None, use_raw=False,
    use_latent_time=None, perc=None, min_r2=0.01,
    min_likelihood=0.001, r2_adjusted=None,
    use_highly_variable=True, diff_kinetics=None, copy=False, **kwargs,
)
```

Deterministic and stochastic modes require prepared `Ms`/`Mu`; stochastic additionally uses second-order information. Dynamical mode requires `tl.recover_dynamics` first. When dynamical parameters are absent, scVelo can warn and fall back to stochastic. A scientific wrapper must reject the request instead of reporting a dynamical result that did not run.

The exact 0.3.4 steady-state implementation writes `layers[vkey]`, `layers[f"variance_{vkey}"]`,
`uns[f"{vkey}_params"]`, and the complete `var` family `f"{vkey}_{offset,offset2,beta,gamma,qreg_ratio,r2,genes}"`;
dynamical mode writes the velocity layer, `f"{vkey}_u"`, mask and parameter record. Its unconditional
`strings_to_categoricals(adata)` call can also change unrelated `obs`/`var` string-column dtypes. Those backend
side effects are not scientific outputs: an atomic wrapper must restore unrelated annotations while retaining and
fingerprinting every mode-specific output, and its collision policy must cover the entire key family.

Official sources and method references:

- [scVelo `tl.velocity`](https://scvelo.readthedocs.io/en/stable/scvelo.tl.velocity.html)
- [scVelo getting started](https://scvelo.readthedocs.io/en/stable/getting_started.html)
- [scVelo citation guide](https://scvelo.readthedocs.io/en/stable/)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
- La Manno et al. 2018, DOI [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6)
- Bergen et al., *Molecular Systems Biology* 2021, DOI [10.15252/msb.202110282](https://doi.org/10.15252/msb.202110282)

## Current interface findings

The current default is deterministic whereas official 0.3.4 defaults to stochastic. Either can be supported, but the executed model, prerequisites and fixed thresholds must be explicit. The current node validates none of them and does not verify the output layer/selected genes or detect backend fallback. It also leaves the complete steady-state parameter family and unrelated annotation-dtype mutation outside any collision or preservation contract.

In the locked Pandas 3.0.5 runtime, dynamical mode exposes another upstream compatibility defect: 0.3.4 obtains read-only arrays from fitted `adata.var` columns and `compute_divergence` mutates `std_u` in place. A wrapper may use a lock-protected compatibility context that verifies the exact `compute_divergence` signature, temporarily replaces it with a pass-through that copies only read-only ndarray inputs, and restores it in `finally`. This does not change the public `scv.tl.velocity` call or numeric values and must be disclosed with dynamic package versions.

Velocity is a per-gene vector field estimate conditional on a kinetic model. It is not fate probability, lineage proof, measured time, or a replicate-aware differential analysis. Poor spliced/unspliced signal, non-steady biology, transcriptional bursts, batch effects and model misspecification can reverse or destabilize directions. Report model-specific fit diagnostics and number/fraction of usable genes.

`fit_likelihood` is an exponentiated normal-density log likelihood, not a probability; finite values and
`min_likelihood` must be nonnegative but are not bounded above by 1.
