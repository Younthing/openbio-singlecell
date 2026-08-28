# MarkerORAAnnotation: official usage and primary-source audit

## Audit scope

- Audited: 2026-08-28.
- Repository surface inspected: `OpenBioSingleCellMarkerORAAnnotation` in
  `openbio_singlecell/nodes_annotation.py`, its registered schema, the Cluster marker evidence domain definition,
  the annotation ADR, `OpenBioSingleCellMarkerGenes`, `OpenBioSingleCellFilterMarkerGenes`, and the test suite.
- External sources are limited to official decoupler API documentation, versioned official source, the decoupler
  paper, and original statistical/resource papers. No community example is needed because official examples exist.
- decoupler is an optional, unconstrained dependency in this repository and is not installed in the inspected
  environment. Therefore every completed run must report the version actually imported; this audit must not be
  mistaken for a promise that an arbitrary installed version is compatible.

Within this repository's domain, the scientifically appropriate product is **Cluster marker evidence** supporting a
**Provisional annotation** decision. It is not a Curated annotation, biological truth, or a Sample-level Condition
contrast.

## Current implementation and test coverage

The public schema is:

```text
adata
resource_csv="openbio-singlecell/markers.csv"
groupby="leiden"
use_raw=True
source_column="cell_type"
target_column="genesymbol"
output_column="ora_annotation"
min_n=3
seed=123
    -> adata
```

Execution reads the resource, silently drops rows missing either resource identifier, converts both identifiers to
strings, and calls the decoupler 1.x sequence:

```python
decoupler.run_ora(..., min_n=min_n, seed=seed, use_raw=use_raw)
acts = decoupler.get_acts(output, obsm_key="ora_estimate")
ranked = decoupler.rank_sources_groups(
    acts,
    groupby=groupby,
    reference="rest",
    method="t-test_overestim_var",
)
```

It then takes the first returned source for every cluster and assigns it to all cells in that cluster. There is no
direct behavior test for this node; registration is the only test coverage. No repository copy of the default
`openbio-singlecell/markers.csv` was found, so a fresh checkout does not provide the advertised default resource.

The implementation has material correctness gaps:

1. It uses a removed decoupler 1.x API without constraining the optional dependency.
2. It forces the copied Raw snapshot to float64/CSR and reassigns `.raw` merely to satisfy a backend detail. That
   changes the returned object's audit snapshot instead of treating it as immutable provenance.
3. `use_raw` chooses a storage slot, not an expression-state contract. The ORA input normalization and gene universe
   are not declared or reported.
4. It does not expose `n_up` or `n_background`; decoupler 1.9.2 therefore silently selects the top 5% of each cell's
   values against a nominal background of 20,000.
5. It discards per-cell ORA p-values, performs a second cell-level group-versus-rest test, discards that table, and
   writes one top label even when the statistic is negative, adjusted p-value is nonsignificant, or candidates tie.
6. It silently discards missing resource rows, does not remove duplicate `(source, target)` pairs, does not validate
   organism/gene-identifier namespace, and does not report feature overlap or sources removed by `min_n`.
7. It emits neither an evidence table nor the required strict JSON `summary` and equivalent `code`.

## The decoupler version boundary is semantic, not cosmetic

The current node exactly targets decoupler 1.9.2. The
[official 1.9.2 `run_ora` API](https://decoupler.readthedocs.io/en/v1.9.2/generated/decoupler.run_ora.html)
documents `run_ora`, `min_n`, `seed`, `use_raw`, `.obsm["ora_estimate"]`, and `.obsm["ora_pvals"]`. In that version:

- the selected set is the top 5% of values by default;
- feature order is shuffled by the supplied seed before ordinal ranking to break ties;
- the test is one-tailed Fisher exact;
- the score is `-log10(p)` rather than an odds ratio;
- the default numeric background is 20,000.

The [official 1.9.2 source](https://github.com/scverse/decoupler/blob/v1.9.2/decoupler/method_ora.py)
confirms those operations. Its `break_ties` helper randomizes feature order, so the node's `seed` affects which tied
features cross the top-set boundary; in sparse single-cell rows, the large zero tie makes that a real scientific
choice rather than generic bookkeeping.

decoupler 2.0 deliberately broke that interface. The
[official changelog](https://decoupler.readthedocs.io/en/stable/changelog.html#id14) says 1.x calls are expected to
fail under 2.x and records the relevant changes:

- `run_ora` became `decoupler.mt.ora`;
- `min_n` became `tmin`;
- `get_acts` became `decoupler.pp.get_obsm`;
- `rank_sources_groups` became `decoupler.tl.rankby_group`;
- ORA now reports an odds-ratio statistic, uses a two-sided Fisher exact test, and stores BH-adjusted p-values.

The current 2.x output keys are `.obsm["score_ora"]` and `.obsm["padj_ora"]`, not the 1.x keys. A version branch
that merely renames functions would still change the statistic, tail, adjustment, output meaning, and tie policy.
Silent compatibility dispatch is therefore unacceptable. The supported major version and exact runtime version must
be explicit provenance, and old and new results must never be presented as directly equivalent.

As of this audit, [PyPI lists decoupler 2.2.0 as the latest release](https://pypi.org/project/decoupler/). There is an
additional release-specific reason not to build the refactor around observation-level `mt.ora`: the
[v2.2.0 implementation](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_ora.py#L253-L277)
uses `row > n_up` after ascending ranks, although `n_up` is documented as the number of top features. The
[current main source](https://github.com/scverse/decoupler/blob/main/src/decoupler/mt/_ora.py#L253-L284)
has corrected this to `row > (nvar - n_up)`. Any future use of `mt.ora` needs a release containing that correction
and a regression test with a known ranked row. The preferred explicit-set design below uses `mt.query_set` and is
not exposed to this boundary defect.

## Official observation-level ORA semantics

The current [official `decoupler.mt.ora` API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.ora.html)
accepts AnnData, a DataFrame, or a matrix/axis tuple and a long-format network with `source` and `target` columns.
It documents these material parameters:

- `tmin`: minimum targets per source after pruning to the data features;
- `n_up`: number of top-ranked features selected per observation, defaulting to the greater of top 5% and two;
- `n_bm`: optionally also select bottom-ranked features;
- `n_bg`: numeric background size, default 20,000; `None` requests a feature-specific background;
- `raw` or `layer`: expression storage selection;
- `ha_corr`: Haldane-Anscombe correction, default 0.5.

The documentation recommends normalized values for observational counts and warns against raw integer counts.
Normalization is irrelevant only when the rows are already contrast-level statistics such as log fold changes or
Wald statistics. Thus `use_raw=True` is not a generally correct default for observation scoring; `.raw` says nothing
about whether its values are counts, normalized/logged data, or a biologically appropriate feature universe.

The [official ORA source](https://github.com/scverse/decoupler/blob/main/src/decoupler/mt/_ora.py) builds, for selected
feature set `S`, resource set `F`, and universe `U`:

```text
a = |S intersection F|
b = |F minus S|
c = |S minus F|
d = |U minus (S union F)|
```

Its current score is the Haldane-Anscombe-corrected log odds ratio:

```text
log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
```

`mt.ora` uses a two-tailed Fisher exact test. The shared
[official method runner](https://github.com/scverse/decoupler/blob/main/src/decoupler/mt/_run.py#L109-L118)
then applies Benjamini-Hochberg across sources within each observation and stores the adjusted values, not raw
p-values, under `padj_ora`.

The runner calls `extract(..., shuffle=True)`. The versioned
[official extraction source](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/data.py#L118-L126)
uses a fixed seed of zero to permute features before ordinal ranking. Ties are therefore deterministic for a fixed
version and original feature ordering, but their membership at `n_up` remains an implementation-derived choice.
With many equal zeros or a tied marker-ranking cutoff, reports must disclose that fact rather than imply a unique
biological ordering.

## The official explicit-set API is the better atomic primitive

decoupler now provides [official `decoupler.mt.query_set`](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.query_set.html):

```python
dc.mt.query_set(
    features=selected_genes,
    net=network,
    alternative="greater",
    n_bg=len(universe),
    ha_corr=0.5,
    tmin=min_targets,
)
```

It tests one explicit feature set against every set in a long-format network and returns `source`, `stat`, `pval`,
and `padj`. Its [official source](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_query_set.py)
shows that it:

- converts the queried features and each source's targets to sets;
- builds the same four contingency counts;
- computes the corrected log odds ratio;
- delegates the exact test to `scipy.stats.fisher_exact` with an explicit alternative;
- applies Benjamini-Hochberg across the retained sources in that call.

For positive cell-type marker enrichment, `alternative="greater"` is the defined scientific question. It should be
fixed and disclosed, not changed to a two-sided test merely because observation-level `mt.ora` uses one. A test for
under-representation would be a different capability.

`n_bg=None` in `query_set` does **not** mean the measured-gene universe. The source constructs `d` from the union of
network targets after subtracting queried/network members. This can omit measured genes absent from the resource.
For marker evidence the implementation must instead pre-intersect the marker network and selected marker sets with
the explicit tested-gene universe and pass `n_bg=len(universe)`. This makes all four cells subsets of the same `U`
and guarantees a nonnegative `d`.

## Required marker, ranking, and universe contract

The existing `OpenBioSingleCellMarkerGenes` table has stable columns:

```text
group, gene, score, logFC, p, p_adj, rank, pct_in_group, pct_rest
```

`OpenBioSingleCellFilterMarkerGenes` selects that evidence by minimum log fold change and within-group prevalence,
maximum reference prevalence, and adjusted p-value. The ORA node should consume this **already selected positive
Cluster marker evidence**, not compute another hidden marker/ranking analysis from expression.

The atomic ORA contract is set-based:

- one nonmissing, nonblank `group` and gene identifier per row;
- one unique `(group, gene)` member after explicit validation;
- the table must identify the upstream marker method, reference, direction, filtering criteria, and feature
  identifier namespace;
- all selected genes must belong to an exact, unique tested-gene universe from that same marker analysis;
- the universe must contain every gene that could have passed the marker filter, not only the reported top genes;
- ranking statistics remain provenance and report context; ORA membership is not weighted by rank, score, logFC, or
  marker p-value.

The present MarkerGenes result stores only `input_genes`, a count, and emits only the top `n_genes` per group. It
does not carry the universe's gene identities. Consequently it is not yet sufficient input for a scientifically
defined ORA. The upstream contract must be extended to emit the exact universe (or a typed artifact containing the
ordered identifiers and fingerprint). `adata.var_names`, Raw snapshot variables, the resource target union, and a
hard-coded genome size are not interchangeable substitutes.

An unfiltered top-N ranking is also not a valid implicit selected set: it can contain negative or nonsignificant
genes and truncates at a boundary that may split tied statistics. The safe workflow is:

```text
MarkerGenes (full tested universe recorded)
  -> FilterMarkerGenes (explicit positive marker criteria)
  -> Marker ORA evidence (set test only)
  -> human evidence review
  -> MapClusterAnnotations (provisional/curated decision)
```

If a future UI offers `top_n` selection, it must either include all genes tied at the boundary or require an explicit
deterministic tie policy and report the expanded/truncated set. This parameter belongs in marker selection, not in
ORA.

## Resource contract, duplicates, missing values, and identifier overlap

The resource must resolve to a regular file and have the declared source/target columns. The node should normalize
those columns to decoupler's canonical `source` and `target` names only after validation. Required behavior:

- reject null, blank, non-scalar, or non-stringable biological identifiers before converting values; never turn
  missing values into literal `"nan"`, `"None"`, or `"<NA>"` genes;
- trim surrounding whitespace but do not uppercase, lowercase, remove Ensembl versions, or translate species
  implicitly; those operations can merge distinct identifiers and belong in an explicit identifier-mapping node;
- require unique gene identifiers in the tested universe;
- remove exact duplicate `(source, target)` resource pairs deterministically **before** applying the minimum target
  size, report how many were removed, and preserve the original file fingerprint;
- reject contradictory metadata or identifiers that collapse after the declared normalization;
- allow the same gene in different source sets and the same gene in different cluster selected sets;
- require every selected marker to occur in the tested universe; a violation is an upstream contract mismatch;
- prune resource targets outside the tested universe and report retained/lost counts per source;
- keep selected markers that are absent from all resource sets: they correctly contribute to `c` in the
  contingency table and should be reported as unmatched selected genes;
- fail globally when no resource target overlaps the universe or no source survives pruning, because this usually
  indicates wrong organism or identifier namespace;
- a single cluster with zero overlap may remain in the evidence output as unresolved, but an all-cluster zero-overlap
  run must fail.

The current official decoupler network validator itself rejects duplicate `(source, target)` rows; see
[`decoupler.pp.prune` and `_validate_net`](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/net.py#L49-L125).
The official single-cell tutorial likewise
[removes duplicate cell-type/gene pairs and formats `source`/`target`](https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html#cell-type-scoring-based-on-marker-genes).
Relying on a backend assertion would produce a late implementation error and lose the data-quality counts, so the
node should enforce the contract first.

The resource provenance must include canonical path, byte size, SHA-256, resource/display name, resource version or
retrieval date, organism/taxon, gene identifier namespace, tissue/organ scope when applicable, license, and the
resource's own citation. A custom CSV with only two columns cannot prove these facts; the caller must provide them or
the file must be paired with trusted metadata. The file hash identifies the exact bytes but does not establish their
biological suitability.

## `tmin`/minimum-set-size semantics

`tmin` is the minimum number of resource targets for a source **after intersection with the eligible feature
universe**. It is not:

- the number of selected markers in a cluster;
- the number of overlapping markers (`a`);
- a marker effect-size threshold;
- an annotation confidence threshold.

The old name `min_n` obscures that distinction. Rename it to `min_targets` or `min_set_size`. The existing default of
three is defensible for compatibility because both old and current official single-cell examples use `tmin/min_n=3`
for filtered cell-marker resources, while the generic current API defaults to five. In either case the exact value,
pre/post-universe sizes, and sources removed must be reported. A separate `min_overlap` rule, if used to flag
candidate labels, must be named and reported separately.

## P-values, effect size, and multiple testing

For each cluster/source pair the reportable evidence is not a single opaque score. It must contain at least:

- selected marker count `|S|` and universe size `|U|`;
- resource set size before and after universe intersection;
- `a`, `b`, `c`, and `d`;
- overlap count and exact overlap gene identifiers;
- corrected log odds ratio (`stat`);
- raw Fisher exact p-value and the stated alternative;
- BH-adjusted p-value with a named adjustment family.

Calling `query_set` once per cluster produces a BH family consisting of retained resource sources within that
cluster. That is the natural decision family for choosing among candidate labels for one cluster, but a report that
interprets all cluster/source pairs should also provide a global BH adjustment across the full table or clearly say
why it does not. Both columns should have unambiguous names such as `p_adj_within_group` and `p_adj_global`.

Marker sets overlap and the tests are not independent; an adjusted p-value is evidence under the declared test
family, not posterior probability that a cell type is correct. The log odds ratio is an effect-size direction and
magnitude, not a calibrated confidence score. Small resource sets and a one-gene overlap can yield fragile ranks;
overlap count/genes and resource size must remain visible.

Primary statistical references:

- Fisher RA. On the Interpretation of χ² from Contingency Tables, and the Calculation of P. *Journal of the Royal
  Statistical Society*. 1922;85:87-94. [doi:10.2307/2340521](https://doi.org/10.2307/2340521).
- Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies. *Annals of Human
  Genetics*. 1956;20:309-311.
  [doi:10.1111/j.1469-1809.1955.tb01285.x](https://doi.org/10.1111/j.1469-1809.1955.tb01285.x).
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple
  Testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).

## Cluster ranking and annotation semantics

The old official annotation example performs per-cell scoring, extracts an activity AnnData, and ranks sources by
cluster. The [current official single-cell tutorial](https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html)
uses `tl.rankby_group(..., method="t-test_overestim_var")`, explicitly filters to `stat > 0`, inspects several top
candidate cell types, and describes them cautiously. It does not justify assigning the first row unconditionally.

The [official `rankby_group` source](https://github.com/scverse/decoupler/blob/main/src/decoupler/tl/_rankby_group.py)
shows that `t-test_overestim_var` compares cells in a group with cells in the rest, uses the group sample size for
both sides of a Welch-style t-test, applies BH across sources within each group, and sorts by adjusted p-value, raw
p-value, then statistic. It is cell-level descriptive ranking, not replicate-aware Condition inference.

The current node violates the minimum interpretation contract by discarding `stat`, `meanchange`, `pval`, and
`padj` and selecting the first row even if all candidates are negative/nonsignificant. It also lacks a deterministic
scientific tie rule: resource/source order may decide equal rows.

The refactored evidence node should not write one cell-type label. It should return all tested candidate rows and a
bounded top-candidate summary. Exact ties on the declared ranking tuple must remain multiple candidates and mark the
cluster ambiguous; alphabetical order may stabilize display only and must never manufacture a winner. Clusters with
no positive, adequately overlapping, adjusted-significant candidate remain unresolved. A reviewer may then supply a
mapping to `OpenBioSingleCellMapClusterAnnotations`, which records whether the resulting annotation is provisional or
caller-declared curated.

The original PanglaoDB paper provides a useful primary-source precedent: it used a one-sided Fisher/hypergeometric
test and BH correction and assigned `Unknown` when the top adjusted p-value exceeded 0.05 rather than forcing a
known type. Resource-specific use must cite the actual supplied resource; PanglaoDB is cited only when it is used.

- Franzén O, Gan L-M, Björkegren JLM. PanglaoDB: a web server for exploration of mouse and human single-cell RNA
  sequencing data. *Database*. 2019;2019:baz046.
  [doi:10.1093/database/baz046](https://doi.org/10.1093/database/baz046).

## Background choice and scientific limitations

The background must be all genes eligible to be selected by the upstream cluster marker test under the same
organism and identifier namespace. A fixed 20,000, the whole theoretical genome, the resource target union, only
expressed markers, only HVGs, or only returned top-N genes changes the null model. The choice must be exact and
reportable.

Timmons and colleagues demonstrated that biological, technological, and detection biases can confound functional
enrichment and explicitly argued for appropriate background expression universes. That result applies directly here:
marker-resource enrichment can support annotation review, but it cannot validate the upstream marker analysis or
prove a biological identity.

- Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound functional enrichment analysis of global
  -omics data. *Genome Biology*. 2015;16:186.
  [doi:10.1186/s13059-015-0761-7](https://doi.org/10.1186/s13059-015-0761-7).

Additional report limitations must include resource incompleteness and context specificity, overlapping/nonindependent
marker sets, identifier/species/tissue mismatch risk, dependence on upstream clustering and marker criteria,
within-cluster heterogeneity, and the fact that cells are not independent biological replicates.

## Minimum stable outputs and provenance

The atomic primary output should be a table, followed by `summary` and `code`; it should not mutate AnnData. The
summary must be strict JSON without NaN/Infinity and contain:

- method name, exact alternative, Haldane correction, contingency definition, adjustment families, and all resolved
  thresholds;
- upstream marker operation/method/reference/filter parameters and exact tested-universe identity/fingerprint;
- resource path, SHA-256, name/version/date, organism, identifier namespace, scope, license, citation, row/duplicate/
  missing counts, overlap losses, and sources before/after `min_targets`;
- groups analyzed, selected-marker counts, unmatched markers, complete tested-row count, positive/significant counts,
  unresolved/ambiguous groups, and bounded top candidates with effect size, p-values, overlap count/genes, and
  alternatives;
- limitations and warnings stated in provisional-evidence language;
- dynamic Python, openbio-singlecell, decoupler, SciPy, pandas, NumPy, anndata, and relevant upstream Scanpy versions.

The equivalent source must define one importable function using the resolved inputs. It must validate/deduplicate the
resource, validate the selected marker and exact universe sets, explicitly intersect the network, call
`dc.mt.query_set(..., alternative="greater", n_bg=len(universe), ha_corr=0.5, tmin=...)` per group, reconstruct the
contingency/overlap columns, perform the documented adjustment families and tie/unresolved logic, and return the same
evidence table and report data without ComfyUI-only wrappers.

## Primary software and method references

- [Current official ORA API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.ora.html)
- [Current official explicit-set API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.query_set.html)
- [Current official cluster-ranking API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.tl.rankby_group.html)
- [Current official ORA source](https://github.com/scverse/decoupler/blob/main/src/decoupler/mt/_ora.py)
- [Versioned 2.2.0 explicit-set source](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_query_set.py)
- [Official decoupler 2.x changelog](https://decoupler.readthedocs.io/en/stable/changelog.html)
- [Official current single-cell enrichment tutorial](https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html)
- [Official 1.9.2 cell-annotation example](https://decoupler.readthedocs.io/en/v1.9.2/notebooks/cell_annotation.html)
- Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer
  biological activities from omics data. *Bioinformatics Advances*. 2022;2:vbac016.
  [doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016).

## Replacement generated-report contract

The active Marker ORA Evidence replacement freezes its generated API as
`run_marker_ora_evidence(marker_table, universe) -> (evidence_dataframe, summary_dict)`. `summary_dict` is exactly the
runtime scientific summary and includes the full methods/results/parameters/references/software/warnings/limitations
payload plus dynamic resource citation, byte hash, and all target accounting. It is not a reduced diagnostics object.
The deprecated legacy annotation node remains a migration-only shim and does not produce either output.

## Public compatibility status

Because the audited implementation is deliberately non-executing, its stable ID is retained only to make old
workflows inspectable. The public schema must advertise deprecation explicitly; it must not appear to be an active
annotation method. No official backend invocation, scientific result, `summary`, or `code` is attributable to this
shim. The active, cited method boundary is Marker ORA Evidence followed by reviewed annotation mapping.
