# Diffusion Pseudotime: official usage and scientific practice

## Audited API and prerequisites

The locked Scanpy 1.12.3 interface is:

```python
scanpy.tl.dpt(
    adata, n_dcs=10, *, n_branchings=0, min_group_size=0.01,
    allow_kendall_tau_shift=True, neighbors_key=None, copy=False
)
```

Official documentation requires, in order, a neighbor graph, `scanpy.tl.diffmap`, and an explicitly annotated root cell in `adata.uns["iroot"]`. With `n_branchings=0` it writes `obs["dpt_pseudotime"]`. Scanpy recommends PAGA, rather than DPT's historical recursive mode, for branch detection.

Primary sources:

- [Scanpy `tl.dpt`](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dpt.html)
- [Scanpy `tl.diffmap`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.diffmap.html)
- Haghverdi et al., *Nature Methods* 2016, DOI [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971)
- Wolf et al., *Genome Biology* 2019, DOI [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x)

## P0 findings in the current node

The current implementation checks only the neighbor `uns` key. It neither requires nor computes `X_diffmap`/`diffmap_evals`, so an ordinary documented invocation fails or can accidentally consume stale diffusion state. This is release-blocking.

It selects `flatnonzero(obs[root_column].astype(str) == root_value)[0]`. The root is therefore the first row in a group after lossy string conversion. Row order is not a biological criterion; reordering identical data can reverse/change pseudotime. It can also merge distinct typed labels. A root must be an exact cell ID or be chosen by a deterministic, disclosed rule inside a biologically nominated root population.

The requested `n_dcs` must not exceed the available, graph-bound diffusion basis. Existing pseudotime/`iroot` values must not be overwritten silently. A disconnected graph yields incomparable/unreachable regions relative to one root and should normally be subset before a single ordering is reported.

## Interpretation

DPT is a relative ordering from a chosen root through one graph; it is not measured time, a lineage-tracing result, a causal transition direction, or a `Condition` comparison. Orientation and numeric values depend on root, graph, diffusion basis and preprocessing. The report must name the exact selected root cell and selection rationale, root population if used, graph/diffusion fingerprints, pseudotime distribution, and limitations.
