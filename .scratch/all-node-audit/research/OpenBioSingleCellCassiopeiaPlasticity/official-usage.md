# Cassiopeia EffectivePlasticity: official method usage

## Public primitive and published implementation

The pinned backend is `cassiopeia-mt==2.1.3` (`v2.1.3-mt`, commit `09868dd04c073d81dc60611f6f91cdbab17783f1`). Its verified execution boundary is Linux/WSL Ubuntu 24.04 with Python 3.12.3; native Windows/Python 3.13 installation fails, and macOS/Python 3.13 are unverified. Cassiopeia exposes the Fitch-Hartigan primitive:

```python
parsimony = cassiopeia.tl.score_small_parsimony(
    solved_tree,
    meta_item="cell_state",
    root=None,
)
```

The function copies internally but its top-down tie resolution calls global `numpy.random.choice`; the minimum parsimony count is invariant to the selected optimal labeling, while global RNG mutation is not acceptable. Execution must isolate/restore NumPy state and independently verify the score.

Yang, Jones et al. define:

```text
EffectivePlasticity(subtree) = Fitch-Hartigan parsimony / number of subtree edges

scEffectivePlasticity(cell) = mean(
    EffectivePlasticity(subtree rooted at ancestor)
    for each non-leaf node on the root-to-cell path
)
```

Their released reference code additionally does the following before scoring: (1) removes cells in states represented below 2.5% of the original tree leaves, (2) prunes the tree and collapses unifurcations, and (3) resolves an “exhausted” terminal polytomy with at least three leaf children and multiple states by inserting one state-specific internal node per observed state. It then excludes the zero-edge leaf subtree from each path average. These operations affect both numerator and denominator and must be reproduced or explicitly declared as a deviation. This node reproduces them, using equality-preserving `frequency >= threshold` to match the paper's “less than” wording rather than the reference script's integer-floor artifact.

Sources:

- [Official `score_small_parsimony` API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.tl.score_small_parsimony.html)
- [Pinned compatible small-parsimony source](https://github.com/andrecossa5/Cassiopeia/blob/v2.1.3-mt/cassiopeia/tools/small_parsimony.py)
- [Published KP-Tracer reference implementation](https://github.com/mattjones315/KPTracer-release/blob/master/reproducibility/Figure4_S4/scripts/compute_plasticity_indices.py)
- [Archived study code, DOI 10.5281/zenodo.6354596](https://doi.org/10.5281/zenodo.6354596)
- Yang, Jones et al., *Cell* 2022, DOI [10.1016/j.cell.2022.04.015](https://doi.org/10.1016/j.cell.2022.04.015)
- Fitch, *Systematic Zoology* 1971, DOI [10.2307/2412116](https://doi.org/10.2307/2412116)
- Hartigan, *Biometrics* 1973, DOI [10.2307/2529686](https://doi.org/10.2307/2529686)

## Identity, annotation and inference scope

Tree and AnnData axes must be unique exact strings. Non-categorical state columns remain a type error because the estimand is categorical. Missing tree leaves, null states, and rare-state-filtered leaves are excluded with explicit null score/status and report warnings rather than silently becoming zero. At least two annotated retained leaves and one edge are computationally required; one retained state is allowed and yields zero plasticity with a no-diversity warning.

Both the defensive AnnData copy and the topology × categorical-state dynamic program require an explicit finite working-memory preflight. An input exceeding the declared bound fails before copying or calling Cassiopeia.

The phenotype is a categorical annotation. Current OpenBio provenance marking the exact column as caller-declared Curated is the preferred report-grade evidence. An expert may still execute report-grade mode with missing or mismatched provenance, but the summary must prominently state that limitation; unknown or Provisional status remains visibly labeled. No string coercion may merge distinct category identities.

EffectivePlasticity is descriptive evidence on one inferred lineage topology, not an independent-replicate Condition test. Cell rows are correlated descendants, and the node performs no p-value or population comparison. The statistic is sensitive to state definition, sampling, rare-state cutoff, polytomy handling, root and topology uncertainty; it does not prove reversible transitions, transition direction, timing or causal plasticity.

Report the retained/excluded state counts, pre/post topology sizes, inserted exhausted-polytomy groups, global parsimony/edge denominator, cell-score distribution, annotation/provenance status, exact tree/solver fingerprints, versions, citations and limitations.
