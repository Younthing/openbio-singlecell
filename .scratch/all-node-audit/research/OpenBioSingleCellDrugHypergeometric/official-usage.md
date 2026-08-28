# OpenBioSingleCellDrugHypergeometric — official usage research

## Audited baseline and P0

The current node has two independent defects. It accesses uninitialized `drug.dgidb.dictionary`, and it internally
runs `scanpy.tl.rank_genes_groups(..., method="wilcoxon")` on cells before Pertpy hypergeometric enrichment. Thus it
both fails at runtime and, if minimally patched, would couple exploratory cell-level marker ranking, selection,
resource acquisition, and ORA without an exact tested universe or upstream provenance.

## Official Pertpy behavior and reviewed replacement

Pertpy 1.3.0 sources:

- https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Enrichment.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_enrichment.py
- https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/enrichment.html

`Enrichment.hypergeometric` expects an AnnData with existing `rank_genes_groups`, takes a target dictionary,
selects markers using adjusted p-values/direction, intersects targets with `adata.var_names`, computes one-sided
hypergeometric survival p-values, and corrects within each cluster. That official workflow is suitable for an
interactive exploratory object, but it does not meet OpenBio's direct-artifact/universe integrity requirement.

The hardened replacement uses the mathematically equivalent explicit-set seam, decoupler 2.2
`mt.query_set(features, net, alternative="greater", n_bg=len(universe), ha_corr=0.5, tmin=...)`, followed by
independent contingency verification and BH. Official API:
https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.query_set.html.

## Scientific-use findings

Drug ORA requires three already defined objects: an explicit selected gene set, its exact tested-gene universe, and
a pinned DGIdb drug-target resource. Marker ranking/filtering and resource loading must remain upstream. For Cluster
marker evidence the result is exploratory population-annotation evidence; for Condition results the upstream model
must use Sample as inference unit and Technical batch as nuisance. The node cannot repair invalid upstream inference.

For a generic evidence artifact, producer identity/approval and those inference fields are user declarations rather
than programmatically verified facts. The node hard-checks the selected genes, exact universe, paired fingerprints,
and typed DGIdb content needed for the contingency calculation. Missing, unfamiliar, cell-level, or inconsistent
scope/Sample/batch declarations remain runnable for experts but must be disclosed with a warning that formal
Condition interpretation is not validated.

An overlap indicates that selected genes are overrepresented among a drug's known target claims. It does not indicate
benefit, reversal of a disease signature, target direction, dose, safety, or efficacy. Resource source/license and
gene namespace are central limitations.

## References to emit

- Fisher RA. On the interpretation of chi-square from contingency tables. *JRSS*. 1922;85:87-94.
  https://doi.org/10.2307/2340521
- Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies.
  https://doi.org/10.1111/j.1469-1809.1955.tb01285.x
- Benjamini Y, Hochberg Y. Controlling the false discovery rate.
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Cannon M, et al. DGIdb 5.0. https://doi.org/10.1093/nar/gkad1040
- Heumos L, et al. pertpy. https://doi.org/10.1038/s41592-025-02909-7

## Report/code implications

Report upstream evidence scope/selection/inference unit as caller/producer declarations, exact universe and
fingerprints, all drug-set contingencies,
overlap genes, raw/BH p, alternatives/ties, DGIdb release/hash/license/source accounting, non-clinical limitations,
references, and exact versions. Generated code returns `(evidence_dataframe, summary_dict)` without AnnData ranking
or network access.
