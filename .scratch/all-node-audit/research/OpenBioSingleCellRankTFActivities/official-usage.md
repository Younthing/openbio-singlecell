# OpenBioSingleCellRankTFActivities — official usage research

## Audited baseline

The node calls the removed decoupler 1.x helpers `get_acts` and `rank_sources_groups`, then expects plural legacy
columns `names` and `pvals`. Decoupler 2.2.0 moved these operations to `pp.get_obsm` and `tl.rankby_group`; its result
uses singular `name`, `pval`, and `padj`. The current implementation therefore fails before producing a table.

## Official decoupler 2.2 usage

Official interfaces and examples:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.get_obsm.html
- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.tl.rankby_group.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/tl/_rankby_group.py
- https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html

The official pattern after ULM is:

```python
import decoupler as dc

scores = dc.pp.get_obsm(adata, key="score_ulm")
ranked = dc.tl.rankby_group(
    adata=scores,
    groupby="cell_type",
    reference="rest",
    method="t-test_overestim_var",
)
```

`rankby_group` accepts `wilcoxon`, Welch `t-test`, or `t-test_overestim_var`. In 2.2.0, `wilcoxon` calls SciPy
`ranksums`, a two-sided large-sample rank-sum approximation without tie correction. Both t-test branches are
two-sided Welch tests. For every group and regulator it returns
`group`, `reference`, `name`, `stat`, `meanchange`, `pval`, and `padj`. `padj` is Benjamini-Hochberg correction within
the complete regulator family for that group. `t-test_overestim_var` deliberately substitutes the tested group's
size for the reference size while retaining the reference mean and standard deviation; it is conservative only when
that replacement decreases the reference size. Positive statistics/mean changes identify
activities higher in the group; negative results are equally part of the tested family.

### Versioned-source corrections discovered during implementation

`dc.pp.get_obsm` constructs a new AnnData but passes the source `obs`, `uns`, and `obsm` objects into its constructor;
it is not a documented defensive-copy seam. The adapter therefore calls it only on a private carrier and audits the
carrier before and after both backend calls.

For a list reference, decoupler 2.2 builds the reference mask as the union of the listed labels but still emits
backend rows for labels that are themselves members of that list; those rows compare a group with a reference pool
that includes the same cells. The canonical adapter excludes every explicit reference-member label from the tested
groups, while verifying that the remaining rows use the declared union. For a single reference string, 2.2 skips
that reference group directly. `rest` tests every observed group against its complement.

The implementation visits annotation labels in first-observed order and excludes globally unused categorical
levels. It does not validate missing, duplicate, or unknown list entries and renders a list reference with
`", ".join(...)`; the adapter consequently performs exact, collision-safe parsing and requires a nonempty, distinct,
all-known list before calling the backend.

The 2.2 implementation replaces a non-finite backend `pval` with 1 before BH, but can leave `stat` non-finite (for
example under constant values or undersized t-test groups). Such rows are not valid report evidence. The adapter
preflights method-specific group sizes and requires every canonical statistic, mean change, p-value, and adjusted
p-value to be finite; it independently recomputes BH within each complete group-by-regulator family.

The returned `meanchange` is the arithmetic difference of activity means. Its sign can disagree with the rank-sum
statistic, so canonical `direction` is defined by the statistic while `mean_change` is preserved and sign
discordance is counted in the report. Canonical ordering uses `(p_adjusted, p_value, -statistic)` as the scientific
key; identical keys share the minimum/competition rank, and regulator ID is only a stable display tie-breaker.

The official implementation iterates all observed annotation groups and all score features. Filtering to p-values
before returning top rows is not part of the method and loses the complete tested family. A report should preserve
all rows, then mark or summarize results using `padj`, not the raw `pval`.

## Statistical-unit and interpretation findings

This procedure treats observations as independent. In a single-cell object, cells from the same Sample violate that
assumption for Condition inference. It is acceptable only as exploratory annotation-associated TF activity evidence,
analogous to Cluster marker evidence. The public interface should therefore name an annotation key and annotation
status, not a generic group key that invites a Condition column. Formal Condition inference requires Sample-level
aggregation/modeling within one defined population.

Group labels must be complete and scalar, and each tested/reference group must have enough finite observations for
the selected test. Category order, group sizes, reference membership, and the full multiple-testing family must be
reported. No result establishes that a TF binds a target or causes the observed phenotype.

## References to emit

- Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics
  data. *Bioinformatics Advances*. 2022;2:vbac016. https://doi.org/10.1093/bioadv/vbac016
- Müller-Dott S, et al. Expanding the coverage of regulons from high-confidence prior knowledge for accurate
  estimation of transcription factor activities. *Nucleic Acids Research*. 2023;51:10934-10949.
  https://doi.org/10.1093/nar/gkad841
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple
  testing. *Journal of the Royal Statistical Society B*. 1995;57:289-300.
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## Report and generated-code implications

The report must disclose that the unit is a cell and the result is exploratory; record annotation status, group and
reference sizes, method, complete tested-family size, positive/negative counts, adjusted-p-value threshold used only
for summaries, bounded leading regulators, TF-activity artifact/network provenance, limitations, references, and
dynamic versions. SciPy and the selected Welch or rank-sum method reference are reported dynamically. Equivalent
code must reconstruct the score AnnData through `dc.pp.get_obsm`, call
`dc.tl.rankby_group`, validate the complete canonical result, and return `(table, summary_dict)`.
