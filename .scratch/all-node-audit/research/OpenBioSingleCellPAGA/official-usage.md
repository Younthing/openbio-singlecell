# PAGA: official usage and scientific practice

## Audited API and references

The locked Scanpy 1.12.3 interface is:

```python
scanpy.tl.paga(
    adata, groups=None, *, use_rna_velocity=False,
    model="v1.2", neighbors_key=None, copy=False
)
```

It coarse-grains a named single-cell neighbor graph over a categorical partition. Scanpy stores the grouping key and sparse group-level `connectivities` and `connectivities_tree` under `uns["paga"]`, and also writes the represented cell counts to `uns[f"{groups}_sizes"]`. In Scanpy 1.12.3, `connectivities` is symmetric, while `connectivities_tree` is the one-orientation sparse encoding returned from the minimum-spanning-tree construction; the latter is not a symmetric adjacency matrix until it is explicitly symmetrized for graph diagnostics. PAGA connectivity is a ratio of observed to expected inter-partition connections under its null graph model; the authors explicitly do not present it as a p-value.

Primary sources:

- [Scanpy `tl.paga`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.paga.html)
- [Scanpy `pl.paga`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.paga.html)
- Wolf et al., *Genome Biology* 2019, DOI [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x)
- Wolf et al., *Genome Biology* 2018, DOI [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

## Required contract

`groups` must name a complete categorical observation column aligned to the graph. Category order defines the rows and columns of the PAGA matrices, so implicit string conversion, unused categories, missing labels, or changed category order can corrupt interpretation. At least two represented groups are required.

`neighbors_key` identifies the exact upstream graph; existence of an `uns` key alone is not enough. Graph shape, pointer keys, finite nonnegative weights, observation alignment and provenance must be validated. Model `v1.2` uses the named distance graph's sparsity pattern to count edges, so the distance matrix is part of the scientific input even though connectivity diagnostics are also reported. Existing `uns["paga"]` or `uns[f"{groups}_sizes"]` must not be overwritten without explicit authorization.

`use_rna_velocity=True` switches to a directed velocity-derived graph and different scientific claim. It additionally requires `uns["velocity_graph"]` and remains marked change-prone in official documentation. It does not belong behind an unnoticed boolean in the undirected cluster-connectivity node.

Plot layout/edge thresholding is separate from computing PAGA. `pl.paga` may add `uns["paga"]["pos"]`; this analysis node should not create or claim a particular drawing.

## Interpretation and reporting

PAGA is exploratory coarse topology conditional on the selected partition and graph. Edge weights are not probabilities, p-values, lineage proof, or branch confidence intervals. Report represented category labels/counts, graph identity, both matrices, edge-weight distribution and strongest edges, components, fixed model, and versions. State that no direction/root/pseudotime or `Sample`-level `Condition` inference was performed.
