# Cassiopeia VanillaGreedy tree reconstruction: official usage

## Versioned public workflow

This node uses the exact build documented by the upstream preparation node: `cassiopeia-mt==2.1.3`, tag `v2.1.3-mt`, commit `09868dd04c073d81dc60611f6f91cdbab17783f1`, sdist SHA-256 `5934043af6506a6439576816be49e170bc41dbe3766eaf49d14d05c3577a3e24`. Runtime must check distribution version and callable signatures. The verified execution boundary is Linux/WSL Ubuntu 24.04 with Python 3.12.3; native Windows/Python 3.13 installation fails, and macOS/Python 3.13 are unverified. See the preparation-node audit for the exact Python/NumPy/pandas metadata constraints and release marker.

The audited reconstruction is deliberately one solver:

```python
tree = cassiopeia.data.CassiopeiaTree(
    character_matrix=character_matrix,
    missing_state_indicator=-1,
    priors=priors,
)

solver = cassiopeia.solver.VanillaGreedySolver(
    missing_data_classifier=(
        cassiopeia.solver.missing_data_methods.assign_missing_average
    ),
    prior_transformation="negative_log",
)
solver.solve(tree, collapse_mutationless_edges=False)
```

`VanillaGreedySolver` is a top-down approximate parsimony method. Its public constructor takes the classifier callable, not the text `"average"`; passing the text is not equivalent. `prior_transformation` supports the documented `negative_log`, `inverse`, and `square_root_inverse` transformations. `collapse_mutationless_edges=True` reconstructs ancestral characters and changes topology, so it is an explicit scientific parameter.

VanillaGreedy contains no advertised random seed and the audited 2.1.3 source does not sample RNG state. Nevertheless, tied mutations are resolved by deterministic encounter order. The upstream character artifact must therefore provide canonical cell/character/state order, and execution must prove that Python/NumPy random state and input objects were not changed. A repeated canonical input should yield one topology fingerprint.

The root created by VanillaGreedy is a **solver-generated internal root**. It is not an observed outgroup and it is not the implicit all-zero sample used by `NeighborJoiningSolver(add_root=True)`. The previous design's “synthetic uncut root for Greedy” statement was incorrect and is retired. Neighbor joining, Hybrid and ILP are different solver/resource/rooting interfaces and are not hidden behind this node.

Sources:

- [Official reconstruction tutorial](https://cassiopeia-lineage.readthedocs.io/en/latest/notebooks/reconstruct.html)
- [Official `CassiopeiaTree` API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.data.CassiopeiaTree.html)
- [Official `VanillaGreedySolver` API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.solver.VanillaGreedySolver.html)
- [Pinned compatible source](https://github.com/andrecossa5/Cassiopeia/tree/v2.1.3-mt)
- Jones et al., *Genome Biology* 2020, DOI [10.1186/s13059-020-02000-8](https://doi.org/10.1186/s13059-020-02000-8)

## Postflight and interpretation

A valid result is one rooted directed tree with exactly the selected computable character-matrix cells as leaves, one root, no cycle, no duplicate node ID, `edges = nodes - 1`, every non-root node having exactly one parent, finite nonnegative branch lengths, and a connected underlying topology. The matrix, priors, missing indicator and upstream fingerprints must still match after solving. Duplicate character profiles may legitimately be placed as sister leaves and remain distinct biological cell identities. Upstream QC warnings are disclosed but do not prevent expert reconstruction.

The report must disclose the solver-generated root, classifier, prior transformation, mutationless-edge policy, topology/axis/prior fingerprints, multifurcation and depth diagnostics, and any duplicate profiles. A reconstructed tree is an inferred topology conditional on one character/QC/solver policy; it is not ground truth and does not quantify topology uncertainty.
