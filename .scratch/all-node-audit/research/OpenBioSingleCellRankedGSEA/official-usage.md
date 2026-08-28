# OpenBioSingleCellRankedGSEA — official usage research

## Audited baseline

The current node calls removed decoupler 1.x `run_gsea` and expects a three-value return. It takes a generic marker
table without validating producer, completeness, tested-gene universe, ranking fingerprint, direction, or
truncation; drops duplicate genes silently; sizes sets before intersecting with the tested universe; and emits no
adjusted p-value. It is non-executable with supported decoupler 2.x and can silently analyze an invalid partial rank.

## Official decoupler 2.2 interface

Primary sources:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.gsea.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_gsea.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_run.py
- https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_psbk.html

For one complete ranked feature vector:

```python
scores, p_adjusted = dc.mt.gsea(
    data=ranked_frame,  # one row, unique genes as columns
    net=network,
    tmin=min_targets,
    raw=False,
    empty=False,
    verbose=False,
    times=n_permutations,
    seed=random_seed,
)
```

With `times > 1`, decoupler 2.2 stores/returns normalized enrichment score (NES) and BH-adjusted empirical p-values.
The public interface does not expose raw p-values or leading-edge members. The report/table must not relabel the
adjusted values as raw p-values.

Versioned source review found a reproducibility hazard: `_ridx` only shuffles when `seed` is truthy, so `seed=0`
leaves every permutation unshuffled. The backend accepts zero, so an expert may execute it, but the report must warn
that the resulting empirical null is not a valid shuffled-permutation analysis. `times` defaults to 1000 upstream;
`times > 1` is the computational minimum for the returned adjusted empirical p-values. Values below 1000 remain
exploratory/low-resolution and require a strong warning rather than an unnecessary hard failure. Exact ties and the
ordered gene vector must be pinned because rank order affects results.

## Scientific-use findings

GSEA consumes a complete continuous ranking, not a thresholded marker list. For formal Condition interpretation the
ranking should be a Sample-level differential statistic within one defined population. A cluster-marker one-vs-rest
ranking can be accepted only as exploratory Cluster marker evidence and must never be presented as Condition
inference. Sample remains the replicate upstream; Technical batch belongs in the differential model as nuisance.

Those producer, scope, Sample, and Technical-batch fields are scientific declarations supplied with a generic table;
they are not facts that this downstream node can prove. The adapter may hard-validate the table/universe schema,
current-content fingerprints, axis identity, rank completeness, and resource bytes. Missing, unfamiliar, or
non-replicate-aware inference declarations must instead be preserved in `summary` with a caution that formal
Condition interpretation is unsupported; they must not prevent an expert from running the mathematically defined
enrichment calculation.

The exact tested-gene universe is the complete ranking. Truncation, arbitrary duplicate selection, missing/non-finite
scores, gene translation, or changing the resource after ranking invalidates the result. Gene-set sizes are computed
after exact resource/universe intersection. Resource organism, namespace, scope, version/date, license, citation,
and byte hash must be reported.

## References to emit

- Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide
  expression profiles. *PNAS*. 2005;102:15545-15550. https://doi.org/10.1073/pnas.0506580102
- Badia-i-Mompel P, et al. decoupleR. *Bioinformatics Advances*. 2022;2:vbac016.
  https://doi.org/10.1093/bioadv/vbac016
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995;57:289-300.
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## Report/code implications

Disclose evidence scope/contrast/direction as caller/producer declarations, upstream inference unit/design (including
missing or internally inconsistent declarations), universe and ranking fingerprints,
resource/hash/intersection, permutations/positive seed, NES and adjusted-p family, set exclusions, unresolved ties,
lack of raw p/leading edge in this backend, references, and exact versions. Equivalent code must return
`(evidence_dataframe, summary_dict)` using the same complete-rank/resource validation and public 2.2 call.
