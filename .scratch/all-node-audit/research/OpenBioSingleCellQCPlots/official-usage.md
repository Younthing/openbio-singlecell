# QC Plots: official usage research

Researched: 2026-08-28

## Official examples

Scanpy's PBMC and Pearson-residual tutorials visualize detected genes, total counts, and mitochondrial percentage before selecting thresholds. Examples use violin/histogram views and scatter plots of total counts against detected genes and mitochondrial percentage.

- QC metrics API/example: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html
- PBMC workflow: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html
- Pearson-residual preprocessing tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html

## Scientific practice

QC plots support threshold selection; they do not themselves certify cell quality. Missing mitochondrial annotations must be visible, not silently plotted as a valid zero-valued distribution. Summary statistics should accompany the image so a methods/results report can disclose the inspected distributions.

## References

- Wolf FA, Angerer P, Theis FJ. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- McCarthy DJ, Campbell KR, Lun ATL, Wills QF. *Bioinformatics*. 2017;33(8):1179-1186. https://doi.org/10.1093/bioinformatics/btw777

## Implications for implementation

- Prefer calculated `obs` metrics when present; otherwise derive counts/detected genes from a disclosed expression source.
- If mitochondrial information is unavailable, omit/mark that panel rather than encode missing values as biological zero.
- Return the plot, a structured report with distribution summaries, and equivalent plotting code.

## Open expert boundary clarification

The plotting adapter trusts the user's explicit fallback source and does not inspect analysis history or require Raw to match current `var`. Raw feature annotations are used when Raw is selected and available; otherwise mitochondrial evidence is omitted with a warning. Finite signed/fractional fallback values remain plottable and are disclosed as unverified count-like evidence. Only absence of any finite plottable primary values is a hard plotting failure.
