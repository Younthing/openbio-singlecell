# Cassiopeia clade expansion testing: official usage

## Audited API and estimand

The pinned implementation is `cassiopeia-mt==2.1.3` (`v2.1.3-mt`, commit `09868dd04c073d81dc60611f6f91cdbab17783f1`). Its verified execution boundary is Linux/WSL Ubuntu 24.04 with Python 3.12.3; native Windows/Python 3.13 installation fails, and macOS/Python 3.13 are unverified. The public call is:

```python
tested_tree = cassiopeia.tl.compute_expansion_pvalues(
    solved_tree,
    min_clade_size=10,
    min_depth=1,
    copy=True,
)
```

The routine initializes every node's `expansion_pvalue` to 1, then tests a child clade only when its leaf count is at least `min_clade_size` and its node-count depth from the root is at least `min_depth`. For an eligible child with `b` leaves, whose parent subtree has `n` leaves and `k` children, the audited source computes the one-sided neutral-coalescent probability

```text
C(n - b, k - 1) / C(n - 1, k - 1).
```

Consequently a stored value of 1 does not prove that a node was tested; eligibility must be independently enumerated from topology. Use `copy=True`, verify the exact public signature, and compare every eligible backend value with the formula before reporting it.

Sources:

- [Official expansion-p-value API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.tl.compute_expansion_pvalues.html)
- [Pinned compatible topology source](https://github.com/andrecossa5/Cassiopeia/blob/v2.1.3-mt/cassiopeia/tools/topology.py)
- Yang, Jones et al., *Cell* 2022, DOI [10.1016/j.cell.2022.04.015](https://doi.org/10.1016/j.cell.2022.04.015)
- Griffiths and Tavaré, *Stochastic Models* 1998, DOI [10.1080/15326349808807471](https://doi.org/10.1080/15326349808807471)
- Benjamini and Hochberg, *JRSS B* 1995, DOI [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## Statistical family and reporting

One run defines one hypothesis family: all topology nodes satisfying the declared size/depth rules in one fixed inferred lineage tree. Apply Benjamini-Hochberg once to the complete finite eligible family, including raw p-values of 1, with stable tie handling. Never prefilter raw significance before correction. Nodes outside the family receive null p/q values and an exclusion reason; they must not be labeled tested.

The unit here is a **clade within one reconstructed lineage**, not a biological Sample or independent clone replicate. Dependence among nested clades is expected; BH controls the stated within-tree family under its assumptions but does not control a collection of trees run separately. Cross-lineage or Condition claims require a separately defined replicate-aware model and multiplicity family.

The full table reports lineage, node/parent IDs, depth, leaf count/fraction, eligibility/reason, raw p, BH q, tested and significant status, plus topology fingerprint. Report the strongest adjusted expansions and complete family size. The neutral coalescent null and inferred topology are assumptions; a call does not establish selection, fitness, mechanism, causality or timing.
