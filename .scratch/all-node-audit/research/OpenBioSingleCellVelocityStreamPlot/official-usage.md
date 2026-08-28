# Velocity Stream Plot: official usage and scientific practice

## Official operation

scVelo 0.3.4 provides:

```python
scvelo.pl.velocity_embedding_stream(
    adata, basis=None, vkey="velocity", density=2,
    smooth=None, min_mass=None, cutoff_perc=None,
    arrow_color=None, arrow_style="-|>", arrow_size=1,
    max_length=4, integration_direction="both", ...
)
```

The public 0.3.4 defaults for `smooth` and `min_mass` are `None` and are resolved internally from grid density. The reviewed OpenBio adapter deliberately exposes and passes explicit reproducible values `smooth=0.5` and `min_mass=1.0`; those are wrapper policy, not claimed scVelo signature defaults.

It projects a previously computed velocity graph/vector field onto an existing embedding and draws a smoothed grid stream field. If the embedded velocities are absent, plotting may compute/cache them. A read-only plot adapter must operate on a private copy and verify that the original scientific object is unchanged.

Sources:

- [scVelo getting started and embedding plots](https://scvelo.readthedocs.io/en/stable/getting_started.html)
- [scVelo plotting API](https://scvelo.readthedocs.io/en/stable/api.html)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
- La Manno et al. 2018, DOI [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6)

## Preconditions and interpretation

Require a verified velocity-graph state, a finite numeric `n_obs × 2` embedding (`obsm["X_<basis>"]`), and an optional complete categorical/numeric color column. Validate exact vkey/graph/obs fingerprints. Do not accept an arbitrary AnnData merely because it has similarly named keys.

Streamlines are interpolation/visualization, not individual-cell paths, transition probabilities, branch confidence, fate probabilities or measured time. Density, smoothing, grid support and the embedding geometry influence appearance. Report these parameters and the fraction of cells/locations with usable vectors; never infer a biological direction solely from the picture.
