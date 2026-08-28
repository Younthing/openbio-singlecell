# OpenBioSingleCellLianaCommunication — official usage research

## Audited baseline

LIANA is not currently pinned as an optional dependency. LIANA 1.10.0 became the current PyPI release on 2026-08-27;
this implementation deliberately targets and exact-locks the separately audited LIANA 1.9.0 wheel released on
2026-08-20. Version 1.10.x must fail closed until its evolving method/resource/plotting contracts receive a separate
audit. The current node runs
one pooled cell-level analysis with `liana.method.<method>`, stores a hidden DataFrame in `adata.uns`, and records no
Sample, Condition, organism, resource content, permutations, seed, identity status, or result-column semantics.

Pooling cells across biological specimens treats cell abundance as replication and can let Samples with more cells
dominate scores. It cannot support a Condition claim. LIANA's documented `.by_sample` interface is the appropriate
descriptive seam when multiple Samples are present.

## Official LIANA 1.9 usage

Official interfaces and examples:

- https://liana-py.readthedocs.io/en/stable/generated/liana.method.AggregateClass.html
- https://liana-py.readthedocs.io/en/stable/generated/liana.method.rank_aggregate.__call__.html
- https://liana-py.readthedocs.io/en/stable/notebooks/basic_usage.html
- https://liana-py.readthedocs.io/en/stable/notebooks/mofatalk.html
- https://liana-py.readthedocs.io/en/stable/generated/liana.resource.select_resource.html
- https://liana-py.readthedocs.io/en/stable/notebooks/prior_knowledge.html

`AggregateClass.by_sample` exposes `adata`, `sample_key`, `key_added`, `inplace`, and `verbose`, then forwards method
arguments through `**kwargs` to the underlying audited `__call__`. For a multi-Sample study, the audited pattern is:

```python
import liana as li

resource = li.rs.select_resource("consensus")  # bundled HUMAN-symbol resource
work = adata.copy()  # mandatory: LIANA 1.9 by_sample mutates its argument even with inplace=False
table = li.mt.rank_aggregate.by_sample(
    work,
    sample_key="sample",
    groupby="cell_type",
    resource=resource,
    expr_prop=0.1,
    min_cells=5,
    return_all_lrs=False,
    use_raw=False,
    layer="log1p_norm",
    n_perms=1000,
    seed=1337,
    n_jobs=1,
    inplace=False,
    verbose=False,
)
```

`cellphonedb.by_sample` has the same wrapper interface. `.by_sample` splits by the declared Sample, runs the method
inside each Sample, and returns a long DataFrame with the exact Sample-key column added. It does **not** test
Conditions. Each Sample remains an independent context for a later replicate-aware comparison. A LIANA 1.9
implementation detail converts the Sample column to categorical and always creates/fills `adata.uns[key_added]`
before returning, even when `inplace=False`; OpenBio must therefore call it only on a private copy and independently
verify that the caller-owned AnnData is unchanged.

LIANA recommends full-gene library-size-normalized, log1p expression. OpenBio passes `use_raw=False` and the
explicitly selected current-axis layer/X so that the backend does not silently switch representations. Normalized
log1p input remains the default and report recommendation, but expression history and a Raw/current binding are not
runtime credentials. A finite expert-selected nonstandard representation is executed and its scale/completeness
limitations are disclosed rather than prohibited by provenance code.
`expr_prop` requires each ligand/receptor (including complex subunits) to be expressed in a minimum fraction of the
corresponding sender/receiver identity. `min_cells` is the within-Sample minimum for an identity. `n_perms` and seed
affect permutation-based specificity measures.

## Method-specific result semantics

`rank_aggregate` combines ranks from several methods. Lower `magnitude_rank` and lower `specificity_rank` are more
highly prioritized. These values are aggregate rank probabilities/scores and must not be labeled generic adjusted
p-values or “significant.” The full table also carries component method scores whose directions differ.

`cellphonedb` returns `lr_means` (larger is stronger expression magnitude) and `cellphone_pvals` (smaller is more
specific under cell-label permutations). The latter is not a Sample-level Condition p-value and LIANA does not add a
multiple-testing-adjusted column here. A generic `significance_threshold` across both methods is scientifically
wrong.

Both methods share sender `source`, receiver `target`, `ligand_complex`, and `receptor_complex`. The audited 1.9
CellPhoneDB result additionally carries the individual `ligand`/`receptor` subunits and their group means/
proportions; rank aggregation drops those subunit columns and retains the complex identifiers plus its complete
component-score family. The node must preserve the exact method-native fields, document directionality, and never
manufacture missing subunit fields or a universal score.

LIANA 1.9 declares `pandas<3`. An isolated smoke with pandas 2.3.3 confirms both public methods and the plotting
contract; pandas 3.x is outside this audited runtime and must fail before backend execution rather than surface an
opaque internal error.

## Resource, organism, and replicate findings

LIANA's `consensus` and most bundled resources use **human gene symbols**; 1.9 also ships a distinct
`mouseconsensus`, so “all bundled resources are human” is false. This node's `bundled_human` branch excludes
`mouseconsensus`; non-human analysis uses only an explicit local resource and documented mapping/orthology policy.
`li.rs.select_resource` materializes the packaged canonical `ligand,receptor` table, while custom resources are
passed through `resource=`. The public selector returns no release, license, or source-license metadata, so both
bundled and local modes require a strict caller-reviewed metadata object and record exact canonical resource SHA-256.
Silent capitalization, implicit aliasing, and online homology translation are prohibited.

Sample must map to exactly one Condition. Identity labels should be stable across Samples, with provisional versus
curated status disclosed. Insufficient Sample-by-identity strata must be reported rather than silently disappearing.
The result remains descriptive candidate communication: dissociated scRNA-seq provides expression compatibility,
not spatial contact, protein abundance, ligand secretion, receptor activation, direction of causal signaling, or
Condition inference.

## References to emit

- Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference from single-cell
  RNA-Seq data. *Nature Communications*. 2022;13:3224. https://doi.org/10.1038/s41467-022-30755-0
- Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication inference. *Nature Cell
  Biology*. 2024;26:1613-1622. https://doi.org/10.1038/s41556-024-01469-w
- Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit
  ligand-receptor complexes. *Nature Protocols*. 2020;15:1484-1506.
  https://doi.org/10.1038/s41596-020-0292-x
- Kolde R, et al. Robust rank aggregation for gene list integration and meta-analysis. *Bioinformatics*.
  2012;28:573-580. https://doi.org/10.1093/bioinformatics/btr709

## Report and generated-code implications

The report must disclose Sample as context, Sample-to-Condition mapping, identity status and eligible/excluded
Sample-by-identity strata, method and score directions, expression state, resource bytes/metadata/organism, tested
and returned interaction counts, permutations/seed, bounded per-Sample key results, absence of Condition inference,
limitations, references, and dynamic versions. Equivalent code must call the 1.9 `.by_sample` interface with an
explicit resource and expression state, validate native columns and Sample identity, and return `(table,
summary_dict)` from a private work copy without downloads or caller-owned AnnData mutation. The portable result is
the canonical DataFrame plus strict metadata; the process-local typed wrapper is not required outside OpenBio.
