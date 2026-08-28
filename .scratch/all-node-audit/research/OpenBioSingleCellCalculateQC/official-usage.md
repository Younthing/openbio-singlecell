# Calculate QC Metrics: official usage research

Researched: 2026-08-28

## Official interface and example

Scanpy documents `scanpy.pp.calculate_qc_metrics(adata, *, expr_type="counts", var_type="genes", qc_vars=(), percent_top=(50, 100, 200, 500), layer=None, use_raw=False, inplace=False, log1p=True)`. Its example annotates mitochondrial genes in `adata.var`, calls the function with `qc_vars=["mito"]` and `inplace=True`, then visualizes total counts, detected genes, and mitochondrial percentage.

- API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html
- PBMC workflow: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html
- Pearson-residual preprocessing tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html

The official API says `percent_top` is 1-indexed, accepts `None`, and measures cumulative expression among the top-ranked genes. It supports dense, CSR, and CSC matrices, and can read `X`, a named layer, or `raw`.

## Scientific practice

The metrics are descriptive QC evidence. Mitochondrial, ribosomal, or hemoglobin fractions can be tissue-, species-, assay-, and Sample-dependent; matching zero genes must be disclosed rather than silently interpreted as zero biological signal. Fixed thresholds are not inferred by this node.

## References

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- McCarthy DJ, Campbell KR, Lun ATL, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. *Bioinformatics*. 2017;33(8):1179-1186. https://doi.org/10.1093/bioinformatics/btw777
- Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control. *Bioinformatics*. 2021;37(7):963-967. https://doi.org/10.1093/bioinformatics/btaa751

## Implications for implementation

- Expose a dynamic expression source with `X`, `raw`, and named layer options; default to `X` for compatibility.
- Pass the resolved layer/use-raw choice to Scanpy instead of copying the matrix into `X`.
- Report matched gene counts, requested and effective `percent_top`, distribution summaries, warnings, and software versions.
- Keep species-specific gene-identification patterns user-configurable; treat them as advanced technical inputs.

## Open expert boundary clarification

The official function accepts an explicitly selected `raw` representation and does not require analysis-history evidence. The wrapper trusts an explicit `X`/`raw`/layer choice and does not require Raw to be a current-axis provenance binding. When Raw has a broader feature axis, QC is calculated in that selected feature space; observation metrics return to the current object and variable-level metrics are mapped only where feature identity is unambiguous, with coverage disclosed.

Finite signed or non-integer values are not proof of UMI counts, but count-likeness is a scientific suitability judgment rather than a universal Python precondition. Such inputs proceed with a prominent machine-readable warning whenever the requested Scanpy calculation remains mathematically defined. Non-finite values and an actually invalid requested transform domain remain hard failures.

When a requested `percent_top` rank exceeds the selected feature count, the wrapper reports the all-feature fraction explicitly. For any non-zero signed expression sum that fraction is 100%; a zero denominator is undefined rather than evidence of biological zero and must be disclosed as missing in the derived percentage instead of silently encoded as 0%.
