# Filter Cells: official usage research

Researched: 2026-08-28

## Official interface and examples

Scanpy documents `scanpy.pp.filter_cells(data, *, min_counts=None, min_genes=None, max_counts=None, max_genes=None, inplace=True, copy=False)`. The official function accepts only one optional threshold per call. The PBMC tutorial first filters low-feature cells, calculates QC metrics, inspects plots, and then slices `AnnData` with a conjunction of detected-gene and mitochondrial-percentage criteria.

- API: https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.filter_cells.html
- PBMC workflow: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html
- Pearson-residual preprocessing tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html

## Scientific practice

Filtering thresholds are dataset decisions, not universal defaults. Count depth, detected genes, and mitochondrial percentage should be inspected jointly and, for multi-Sample studies, checked by Sample because technical depth and tissue composition can differ. The node must report exactly what was removed; it must not claim that retained cells are universally “high quality.”

## References

- Wolf FA, Angerer P, Theis FJ. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- McCarthy DJ, Campbell KR, Lun ATL, Wills QF. *Bioinformatics*. 2017;33(8):1179-1186. https://doi.org/10.1093/bioinformatics/btw777
- Osorio D, Cai JJ. *Bioinformatics*. 2021;37(7):963-967. https://doi.org/10.1093/bioinformatics/btaa751

## Implications for implementation

- Preserve conjunction semantics with one boolean mask; do not pass multiple thresholds to one Scanpy call or imply that Scanpy supports that.
- Expose the expression source because counts and detected genes must come from the intended count matrix.
- Report input, retained, and removed cell counts; retention percentage; each active threshold; zero-count cells; and mitochondrial-column availability.
- Emit a warning when no criterion is active rather than silently presenting the unchanged object as filtered.

## Open expert boundary clarification

Scanpy returns a boolean subset mask; an all-false mask is still a well-defined filtering result. The wrapper allows a zero-cell AnnData and reports that every cell was removed instead of enforcing a universal retention rule. Likewise, finite signed or fractional values can still be summed and tested against an expert-specified threshold. They are disclosed as non-count-like expression sums, not rejected solely for violating recommended count practice. Non-finite inputs and ambiguous/non-finite mitochondrial criteria remain hard failures.

## Threshold activation correction

Scanpy uses `None` to distinguish an absent threshold from a numeric threshold. The node therefore must not overload numeric zero as the only disabled sentinel. Total-expression bounds have explicit enable controls and accept any finite signed or fractional value, including zero. The mitochondrial upper bound also has a separate enable control so 0% is expressible. Detected-gene bounds remain non-negative integers because they count detected features; their historical zero-disabled behavior is retained.

Legacy serialized values require an explicit migration: each old non-zero total-expression or mitochondrial threshold becomes the same numeric value with its enable control set to true; each old zero sentinel becomes a disabled control. Merely defaulting all new controls to false would silently disable old non-zero filters.
