# Mark MAD Outliers: official and method research

Researched: 2026-08-28

## Official API semantics

SciPy defines `scipy.stats.median_abs_deviation(x, scale=...)` as the median absolute deviation divided by `scale`. A numeric `scale=1.4826` therefore makes the result smaller; it does **not** apply the common normal-consistency multiplier. The supported `scale="normal"` divides by approximately `0.67449` and is the correct SciPy spelling for a normal-consistent MAD estimate.

- SciPy API and examples: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html

SciPy also provides `nan_policy="omit"`. Missing observations must not silently poison an entire Sample threshold, but their absence must be counted and disclosed.

## Community single-cell usage and method literature

The single-cell best-practices QC workflow defines an MAD-based helper and applies it to QC metrics, including Sample-wise calculations. It also emphasizes that cutoffs are dataset-specific and should be inspected rather than treated as universal constants.

- Community workflow: https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
- Leys C, Ley C, Klein O, Bernard P, Licata L. Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology*. 2013;49(4):764-766. https://doi.org/10.1016/j.jesp.2013.03.013
- McCarthy DJ, Campbell KR, Lun ATL, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. *Bioinformatics*. 2017;33(8):1179-1186. https://doi.org/10.1093/bioinformatics/btw777

## Scientific implications

- Use one direction per node invocation. Metrics that require different biological directions (for example upper-only mitochondrial percentage versus two-sided library complexity) should use separate nodes.
- When grouping is enabled, thresholds are estimated independently within each biological `Sample` or capture, not across a merged study.
- Multiple metrics in one invocation form an explicit union of the same-direction rule; per-metric flags and the final union must be reported.
- `nmads` and the scaling convention are user choices, not validated biological constants.
- Report medians, MADs, effective lower/upper thresholds, missing values, zero-MAD groups, per-metric flags, per-Sample flags, and the union count.

SciPy defines MAD for a one-observation finite slice (the result is zero). A minimum finite group size is therefore an expert stability policy, not a mathematical requirement. The interface permits `minimum_group_size=1`; singleton/zero-MAD thresholds run and are reported with a warning rather than being prohibited.

Grouped-label validation exists only to keep group-specific thresholds and report entries unambiguous. Its diagnostic must describe MAD reporting and must not claim that Scrublet provenance is involved.
