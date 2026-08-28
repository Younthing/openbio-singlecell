# OpenBioSingleCellSampleCompositionSummary — official usage research

Research date: 2026-08-28

## Scope

This node is not a Scanpy statistical test. Its scientific operation is a descriptive aggregation of cell-level observation metadata into one complete Sample-by-annotation count and proportion table.

The relevant official interfaces are:

- AnnData obs, which stores one-dimensional annotations aligned to cells;
- pandas cross-tabulation/reindexing for a complete frequency grid;
- compositional-data literature for interpreting proportions whose rows sum to one.

Under this repository's domain language, an AnnData observation is a cell, while a **Sample** is an independent biological specimen. A Technical batch or a collection of cells must not be substituted for a Sample merely because both appear in obs.

Sources:

- [AnnData obs reference](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.obs.html)
- [AnnData object reference](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html)
- [pandas crosstab reference](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.crosstab.html)
- [pandas reindex reference](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.reindex.html)

## Canonical descriptive operation

For each Sample s and annotation category k:

- cell_count(s, k) is the number of retained input cells assigned to k;
- sample_total_cells(s) is the number of retained input cells in s across the complete annotation universe;
- proportion(s, k) = cell_count(s, k) / sample_total_cells(s).

The canonical table contains every Sample crossed with every annotation category observed anywhere in the input. A category absent from one Sample is a measured zero, not a missing row:

~~~python
import pandas as pd

sample_order = list(dict.fromkeys(obs["sample"]))
annotation_order = list(dict.fromkeys(obs["cell_type"]))
full_index = pd.MultiIndex.from_product(
    [sample_order, annotation_order],
    names=["sample", "annotation"],
)
counts = (
    obs.groupby(["sample", "cell_type"], observed=True)
    .size()
    .reindex(full_index, fill_value=0)
    .rename("cell_count")
    .reset_index()
)
counts["sample_total_cells"] = counts.groupby("sample")["cell_count"].transform("sum")
counts["proportion"] = counts["cell_count"] / counts["sample_total_cells"]
~~~

pandas documents the same zero-preserving behavior for categorical cross-tabulation with dropna=False. An implementation may use crosstab or explicit product/reindex, but it must validate the same complete grid and must not depend on pandas defaults that can omit all-zero categories.

The proposed fixed annotation universe is the set of nonmissing categories observed in at least one input cell. For a pandas Categorical column, preserve declared category order after restricting it to observed categories and report globally unused declared levels separately. This avoids silently treating stale unused categories as measured populations while retaining deterministic ordering. A future explicit expected-category input would be a different interface and is not needed for this refactor.

## Sample and Condition metadata

Each Sample must map to exactly one **Condition**. This is a sample-design invariant, not an optional cleanup:

~~~python
condition_counts = obs.groupby("sample", observed=True)["condition"].nunique(dropna=False)
if (condition_counts != 1).any():
    raise ValueError("Each Sample must map to exactly one Condition.")
~~~

The three metadata columns must be distinct, present, nonmissing, and nonblank. Values that collide only after display-string conversion, such as integer 1 and string "1", must fail rather than merge biological specimens or annotations. Sample identifiers must be unique as labels even when multiple cells share the same Sample, and cell identifiers should be unique so accidental duplicated cells are not silently counted twice.

Missing metadata must not be silently excluded. Dropping an unannotated cell changes that Sample's denominator; dropping all invalid cells from one Sample can remove an independent replicate entirely. The safe default is to fail and require an explicit upstream decision to correct metadata or assign an intentional "Unassigned" annotation.

## Proportions and compositional closure

For every valid Sample, the reported proportions satisfy:

~~~text
sum over annotations of proportion(sample, annotation) = 1
~~~

This closure means the columns are not independent. Increasing one annotation's captured proportion necessarily decreases the sum of the others, even if their absolute abundance in the specimen did not change. The total number of captured cells also depends on acquisition, dissociation, filtering, and sampling depth. Consequently:

- the table describes relative captured-cell composition, not absolute tissue abundance;
- a zero count means no cell of that annotation was captured and retained in that Sample, not proof of biological absence;
- per-Sample proportions are suitable for plots and descriptive reporting, but the table alone is not a **Condition contrast**;
- no pseudocount is needed or appropriate for the descriptive count/proportion table.

The foundational compositional-data reference is Aitchison's analysis of the simplex. The scCODA paper applies the same closure problem directly to single-cell population counts and shows why independent per-population interpretations can be misleading.

Sources:

- Aitchison J. The Statistical Analysis of Compositional Data. JRSS B. 1982;44:139–177. [doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)
- Büttner M et al. scCODA is a Bayesian model for compositional single-cell data analysis. Nature Communications. 2021;12:6876. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)

## Biological replication and annotation status

The propeller method documentation calls its sample input the biological replicate and constructs a cell-type-by-biological-replicate count/proportion matrix. The official pertpy scCODA tutorial likewise aggregates one row per statistical sample and requires every sample-level covariate to be unique within that sample.

Sources:

- [speckle propeller manual](https://bioconductor.org/packages/release/bioc/manuals/speckle/man/speckle.pdf)
- [pertpy scCODA tutorial](https://pertpy.readthedocs.io/en/1.0.1/tutorials/notebooks/sccoda.html)
- Phipson B et al. propeller: testing for differences in cell type proportions in single cell data. Bioinformatics. 2022;38:4720–4726. [doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582)

A descriptive table can summarize a **Provisional annotation**, but its report must name that status and must not present it as a **Curated annotation**. Formal population interpretation should use an expert-reviewed Curated annotation. If annotation provenance cannot be verified, the status is "unknown" and that limitation belongs in the summary.

## Current implementation audit

The current helper in openbio_singlecell/nodes_abundance.py has two important correct behaviors:

- it checks that Sample, grouping, and annotation columns are distinct and present;
- it constructs a complete product of observed Samples and annotations and fills absent Sample-annotation combinations with count zero;
- it checks a one-group-per-Sample mapping;
- it computes proportions within Sample and does not mutate AnnData.

The remaining deviations are material:

- the public and output language says group rather than the domain term Condition;
- all cells missing any of the three metadata values are silently dropped before the Sample-to-Condition mapping check;
- string conversion can merge typed identifiers that render to the same string and can hide blank values;
- a Sample whose cells are all dropped disappears without a study-design failure;
- globally unused categorical levels and ordering policy are not reported;
- the table omits sample_total_cells, so the proportion denominator is not explicit on each row;
- no invariant proves row count equals Samples times annotations, counts sum to all input cells, or each Sample's proportions sum to one;
- no cell-identifier uniqueness or direct-call parameter validation is present;
- the result contains only a table and lacks the required strict-JSON summary and equivalent code output.

The bundled Sample Composition Comparison workflow has four Samples, two per Condition, and three annotation categories. The descriptive node therefore returns the expected 4 × 3 = 12 rows. Existing tests pin that row count and workflow wiring, but do not cover zero cells, missing metadata, identifier collisions, unused categories, denominator invariants, input immutability, strict JSON, or generated-code equivalence.

## Report-ready disclosure

A general scientific report should state:

- the Sample, Condition, and annotation columns used;
- whether the annotation was Curated, Provisional, or of unknown provenance;
- the number of independent Samples and Samples per Condition;
- total input cells and cells per Sample, including min, quartiles, median, and max;
- the observed annotation universe and any unused declared levels;
- the complete grid size and number/fraction of zero-count combinations;
- the denominator definition and a checked closure invariant;
- per-Condition descriptive summaries of each annotation's Sample-level proportion;
- that no hypothesis test was performed;
- that proportions describe captured relative composition and not absolute abundance.

The report must not use cell count as inferential sample size, call a Technical batch a Sample, claim a population is biologically absent from a zero captured count, or turn a Provisional annotation into a Curated annotation.

## Software and method references

- Virshup I et al. anndata: Annotated data. Journal of Open Source Software. 2021;6:4371. [doi:10.21105/joss.04371](https://doi.org/10.21105/joss.04371)
- McKinney W. Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference. 2010. [doi:10.25080/Majora-92bf1922-00](https://doi.org/10.25080/Majora-92bf1922-00)
- Aitchison J. The Statistical Analysis of Compositional Data. JRSS B. 1982;44:139–177. [doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)
- Büttner M et al. scCODA is a Bayesian model for compositional single-cell data analysis. Nature Communications. 2021;12:6876. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)
- Phipson B et al. propeller: testing for differences in cell type proportions in single cell data. Bioinformatics. 2022;38:4720–4726. [doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582)
