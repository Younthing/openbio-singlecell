# OpenBioSingleCellSCENICActivityBinarization — official usage research

## Audited baseline

pySCENIC 0.12.1 exposes `pyscenic.binarization.binarize`. The current node calls it without a seed or worker count,
so stochastic threshold derivation is not reproducible and the report records neither backend settings nor per-
regulon on/off results. It mutates arbitrary AnnData keys and emits no `summary` or `code`.

The upstream implementation has a subtle seed contract: `derive_threshold` calls `numpy.random.seed(seed)` only
when `seed` is truthy. A user-facing seed of `0` therefore behaves like no seed. This repository generally uses zero
as a default seed, so blindly exposing that convention would falsely claim reproducibility.

## Official pySCENIC 0.12.1 usage

Versioned source and documentation:

- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/binarization.py
- https://pyscenic.readthedocs.io/en/stable/faq.html
- https://doi.org/10.1038/s41596-020-0336-2

The public call is:

```python
from pyscenic.binarization import binarize

binary, thresholds = binarize(
    auc_mtx=auc_matrix,
    threshold_overides=None,  # upstream spelling
    seed=1,                   # nonzero is required by 0.12.1
    num_workers=1,
)
```

For each regulon, `derive_threshold` first asks whether its cell-level AUC distribution is bimodal. The default
`method="hdt"` uses Hartigan's dip test. A distribution not called bimodal receives `mean + 2*standard deviation`.
A called-bimodal distribution receives a two-component Gaussian-mixture fit and a kernel-density trough between the
component means. A cell is active only when `AUC > threshold`; equality is off. The function returns a cell-by-
regulon integer matrix plus one threshold per regulon.

The procedure is a data-dependent dichotomization heuristic. In 0.12.1, the bundled dip-test implementation returns
no p-value for a constant distribution or at most four distinct values, so that regulon follows the deterministic
`mean + 2*standard deviation` branch; this is defined but commonly yields all cells off and must be disclosed. The
dip test, mixture fit, and KDE can still fail or become unstable for non-finite or other degenerate distributions.
The thresholds are not
p-values, and binary activity is not evidence of direct binding or a Condition effect. Running with multiple
processes complicates RNG and platform behavior; one worker and a nonzero seed are the auditable default.

The 0.12.1 implementation mutates NumPy's legacy global RNG inside threshold derivation and always constructs a
`multiprocessing.Pool`, even for one worker. A local exact implementation must isolate/save/restore RNG state and
avoid the unnecessary process boundary while retaining the same seeded algorithm. Manual overrides change the method and need an
explicit, typed, fingerprinted input if supported; hidden overrides are unacceptable.

The bundled dip-test wrapper fixes `numt=1000` Monte Carlo simulations and materializes a dense
`1000 * number_of_cells` uniform array for every threshold that reaches the HDT branch. A memory preflight that only
counts the public cell-by-regulon matrices is therefore unsafe. The local single-worker implementation must include
that transient simulation allocation whenever at least one regulon lacks a manual threshold; an all-regulon override
does not execute or allocate the dip-test simulation.

The tagged 0.12.1 code itself uses removed `numpy.msort` and `numpy.float` aliases in its bundled dip test under the
repository's NumPy 2 runtime. A compatible implementation may replace only those aliases with `numpy.sort` and
builtin `float`; this is a compatibility repair, not a change to the statistical procedure, and must be checked
against the real 0.12.1 backend in a compatible smoke fixture.

## References to emit

- Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature
  Protocols*. 2020;15:2247-2276. https://doi.org/10.1038/s41596-020-0336-2
- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463
- Hartigan JA, Hartigan PM. The dip test of unimodality. *The Annals of Statistics*. 1985;13:70-84.
  https://doi.org/10.1214/aos/1176346577

## Report and generated-code implications

The report must disclose the exact threshold algorithm, strict greater-than rule, nonzero seed, one-worker policy,
cell/regulon counts, threshold range, active-cell counts/proportions for every regulon (bounded key-result preview),
manual overrides if any, degenerate checks, SCENIC artifact provenance, limitations, references, and dynamic
versions. Equivalent code must isolate global RNG, call the audited 0.12.1 function or an explicitly verified exact
implementation, validate the complete binary/threshold outputs, and return them with the same strict summary.
