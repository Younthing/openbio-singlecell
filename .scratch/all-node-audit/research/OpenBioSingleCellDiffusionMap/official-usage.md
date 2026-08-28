# Diffusion Map: official usage and scientific practice

## Audited API and references

The locked Scanpy 1.12.3 interface is:

```python
scanpy.tl.diffmap(
    adata, n_comps=15, *, neighbors_key=None, random_state=0, copy=False
)
```

It requires a previously computed observation-neighbor graph. With `neighbors_key`, Scanpy follows `uns[neighbors_key]` to the named connectivity and distance matrices. It writes `obsm["X_diffmap"]` and `uns["diffmap_evals"]`. The zeroth column is the stationary solution; the first informative diffusion component is column 1.

Primary sources:

- [Scanpy `tl.diffmap`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.diffmap.html)
- [Scanpy `tl.dpt`](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dpt.html)
- Haghverdi et al., *Bioinformatics* 2015, DOI [10.1093/bioinformatics/btv325](https://doi.org/10.1093/bioinformatics/btv325)
- Haghverdi et al., *Nature Methods* 2016, DOI [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971)
- Coifman et al., *PNAS* 2005, DOI [10.1073/pnas.0500334102](https://doi.org/10.1073/pnas.0500334102)
- Wolf et al., *Genome Biology* 2018, DOI [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

## Graph and output contract

The diffusion operator is determined by the exact graph, including its upstream representation, distance metric, neighbor count and connectivity method. Checking only that `uns[neighbors_key]` exists, as the current node does, is insufficient. Both matrices and metadata pointers must be present, observation-aligned, finite, nonnegative and internally consistent.

The output namespace is fixed in Scanpy 1.12.3. Existing `X_diffmap` or `diffmap_evals` values must not be overwritten silently. The number of requested components must be feasible for the number of observations and the eigensolver; returned coordinates/eigenvalues must be finite and have matching dimensions.

Disconnected components and isolated observations materially change diffusion geometry. They must be reported; for a downstream single-root pseudotime analysis, a disconnected graph is normally a reason to subset to the biological lineage rather than compare incomparable components.

## Scientific interpretation

Diffusion components are an exploratory manifold representation. Component signs are arbitrary, the stationary component is non-informative, and the geometry depends on preprocessing and graph construction. Diffusion map alone does not choose a biological root, orient time, detect branches, prove lineage relationships, or test a `Condition`.

Report named graph identity/provenance, component count, eigenvalues, graph components/isolates, coordinate diagnostics, seed, output keys, warnings, method references and versions. If the result will feed DPT, record a fingerprint that binds the diffusion basis to the exact graph and observation order.
