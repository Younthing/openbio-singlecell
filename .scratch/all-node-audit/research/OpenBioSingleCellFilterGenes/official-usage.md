# Filter Genes: official usage research

Researched: 2026-08-28

## Official interface and examples

Scanpy documents `scanpy.pp.filter_genes(data, *, min_counts=None, min_cells=None, max_counts=None, max_cells=None, inplace=True, copy=False)`. Only one optional threshold is accepted per call. Official tutorials commonly remove genes detected in fewer than a small number of cells before downstream modeling, but the chosen threshold remains a disclosed analysis decision.

- API: https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.filter_genes.html
- PBMC workflow: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html
- Pearson-residual preprocessing tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html

## Scientific practice

Low-detection filtering reduces uninformative features and computation, but can remove rare-population markers. Filtering should use the intended count matrix and report the number and percentage of retained genes. It is distinct from highly variable gene selection: the latter creates a modeling view, while a post-QC full-gene Raw snapshot should remain recoverable.

## References

- Wolf FA, Angerer P, Theis FJ. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- McCarthy DJ, Campbell KR, Lun ATL, Wills QF. *Bioinformatics*. 2017;33(8):1179-1186. https://doi.org/10.1093/bioinformatics/btw777

## Implications for implementation

- Keep conjunction semantics in a single mask.
- Add an explicit expression source.
- Report retained/removed genes and active criteria, and distinguish this operation from highly variable gene selection.
- Reject contradictory bounds; permit an empty result with explicit disclosure.

## Open expert boundary clarification

An all-false feature mask is a defined outcome and is not blocked merely because downstream tools commonly require genes. Finite signed or fractional expression can also drive an expert-declared numerical filter, but the report stops calling its totals verified counts. Explicit Raw selection is sufficient: no history/full-gene proof is required. When Raw has more features than current `X`, statistics are matched to current features by unambiguous feature identity; unmatched current features are retained and disclosed rather than silently misaligned.

## Threshold activation correction

Total-expression bounds use explicit enable controls rather than treating zero as the disabled state. An enabled bound accepts any finite signed or fractional value, including zero; detected-cell thresholds remain non-negative integers because they count observations. Legacy workflow migration maps each old non-zero total-expression threshold to enabled and each old zero sentinel to disabled while preserving the numeric values.
