# OpenBioSingleCellDrugGSEA — official usage research

## Audited baseline and P0

The current node reads uninitialized `drug.dgidb.dictionary`, internally computes per-cell Wilcoxon marker rankings,
and passes the resulting AnnData state to `pertpy.tl.Enrichment().gsea`. It therefore fails on a fresh resource and
otherwise conflates ranking, resource download, and enrichment without proving rank completeness or the tested-gene
universe.

## Official Pertpy and replacement interfaces

Pertpy 1.3.0 `Enrichment.gsea`:

- https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Enrichment.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_enrichment.py
- https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/enrichment.html

It expects `adata.uns["rank_genes_groups"]`, constructs a gene/score frame per cluster, calls BlitzGSEA, and returns a
dictionary of result DataFrames. That is the documented Pertpy workflow, but fabricating Scanpy private storage from
an external rank would create a brittle adapter. OpenBio should instead use the same public decoupler 2.2
`mt.gsea(data, net, times, seed)` ranked-vector seam as generic Ranked GSEA:
https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.gsea.html.

For `times>1`, decoupler returns NES and BH-adjusted empirical p-values. It does not expose raw p or leading edge via
the public interface. Versioned 2.2 source does not shuffle when seed is zero: the backend accepts that choice, but a
zero-seed result must be labeled as lacking a valid shuffled empirical null. Fewer than the recommended/default 1000
permutations is likewise executable for exploration, with a low-resolution warning; the hard computational minimum
for this interface is two.

## Scientific-use findings

Drug GSEA needs one complete continuous gene ranking and a pinned DGIdb target resource. It must never create the
ranking itself. A formal Condition rank must derive from a Sample-level model within one population, with Technical
batch modeled as nuisance where applicable. A complete Cluster marker rank is allowed only as exploratory Cluster
marker evidence.

Generic artifact producer/approval, scope, direction, Sample-level status, and Technical-batch handling are
caller-supplied declarations, not downstream proof. Complete rank/universe structure and fingerprints remain hard
requirements because they define GSEA. Missing, unfamiliar, cell-level, or inconsistent inference declarations
remain runnable for experts but must force cautious summary language and prohibit claiming validated Condition
inference.

Positive/negative enrichment means drug target claims concentrate toward one end of the supplied statistic. It is
not a connectivity-map signature-reversal test and conveys no treatment direction, efficacy, dose, safety, or
recommendation. Exact DGIdb release/source licenses and HGNC matching must be disclosed.

## References to emit

- Subramanian A, et al. GSEA. *PNAS*. 2005;102:15545-15550.
  https://doi.org/10.1073/pnas.0506580102
- Badia-i-Mompel P, et al. decoupleR. https://doi.org/10.1093/bioadv/vbac016
- Cannon M, et al. DGIdb 5.0. https://doi.org/10.1093/nar/gkad1040
- Heumos L, et al. pertpy. https://doi.org/10.1038/s41592-025-02909-7
- Benjamini Y, Hochberg Y. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## Report/code implications

Report rank scope/direction and upstream design as caller/producer declarations, exact universe/ranking fingerprints, DGIdb release/hash/license and
target coverage, permutations/seed, positive/negative NES and adjusted p, ties/set exclusions, lack of raw p/leading
edge, non-clinical limitation, references, and exact versions. Generated code returns
`(evidence_dataframe, summary_dict)` without Scanpy ranking or downloads.
