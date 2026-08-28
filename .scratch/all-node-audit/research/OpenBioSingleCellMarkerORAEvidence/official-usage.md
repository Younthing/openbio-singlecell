# OpenBioSingleCellMarkerORAEvidence — official usage research

## Scope

This node is the reviewed replacement capability for the statistical part of the legacy
`OpenBioSingleCellMarkerORAAnnotation` node. It tests explicit Cluster marker sets against an explicit marker-set
resource and returns evidence. It does not score cells, mutate AnnData, or assign labels.

The full pre-change investigation of the defective legacy compound implementation is recorded in
`../MarkerORAAnnotation/official-usage.md`. This record narrows the official API and method evidence for the new
node identity before implementation begins.

## Official decoupler interface

The supported API is decoupler 2.x `decoupler.mt.query_set`:

- API documentation: https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.query_set.html
- versioned 2.2.0 source: https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_query_set.py
- network validation/pruning source: https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/net.py
- official changelog documenting the 1.x to 2.x break:
  https://decoupler.readthedocs.io/en/stable/changelog.html

Verified against the versioned 2.2.0 source above, the public signature is
`query_set(features, net, alternative="greater", n_bg=20000, ha_corr=0.5, tmin=5, verbose=False)`.
The official function accepts one query feature set through `features`, a long-form network whose canonical columns
are `source` and `target`, an explicit background size (`n_bg`), minimum target count (`tmin`), one-sided alternative,
and Haldane correction. It returns one enrichment result per retained source. The correct OpenBio adapter therefore
calls it once per Cluster group:

```python
result = decoupler.mt.query_set(
    features=selected_genes,
    net=network,
    n_bg=len(tested_gene_universe),
    tmin=min_targets,
    alternative="greater",
    ha_corr=0.5,
    verbose=False,
)
```

The adapter must require a compatible decoupler 2.x installation and must not fall back to the removed 1.x
`run_ora`, `get_acts`, or `rank_sources_groups` sequence.

## OpenBio provenance adapter (schema v2)

The following SHA-256 fields are OpenBio adapter provenance, not arguments or outputs of
`decoupler.mt.query_set`:

- `analysis_fingerprint` identifies the ranking analysis axes, parameters, ordered observations, and group
  membership;
- `ranking_fingerprint` identifies the complete, non-truncated all-group by all-tested-gene ranking, including its
  numeric evidence;
- `universe_fingerprint` identifies the exact ordered tested-gene universe;
- each artifact's `content_fingerprint` identifies its current canonical DataFrame rows and columns;
- a Filter Marker Genes artifact additionally carries `upstream_content_fingerprint`, identifying the direct Marker
  Genes table that was filtered.

Schema-v2 current-content fingerprints serialize rows in their existing order, strings exactly, and numeric values
as `float(value).hex()` before canonical JSON SHA-256. Marker ORA Evidence accepts only a direct Filter Marker Genes
table paired with a direct Marker Genes universe. The pair validator recomputes both current-content fingerprints,
recomputes the ordered universe identity, and cross-checks the shared analysis, full-ranking, and universe identities.
The ORA adapter then pins the three shared identities, both current-content identities, and the Filter upstream
binding in its result and summary; generated code pins the available identities and recomputes the two supplied
tables' current content and the universe identity before calling decoupler. The complete ranking is not
reconstructible from a filtered table, so its fingerprint is cross-checked across the two direct artifacts and
against the runtime-pinned validator result rather than falsely recomputed from the subset.

## Reviewed scientific usage

For each group, let the exact tested-gene universe be `U`, selected positive markers be `S`, and one resource set be
`R`, all after exact identifier validation. The one-sided over-representation table is:

```text
a = |S ∩ R|
b = |R \ S| within U
c = |S \ R| within U
d = |U \ (S ∪ R)|
```

The node must:

1. require a direct Filter Marker Genes table plus direct Marker Genes universe and validate their shared analysis,
   complete-ranking, and universe fingerprints, their recomputed current-content fingerprints, and the Filter
   artifact's direct-upstream content fingerprint;
2. require every selected gene to be in the exact universe;
3. validate and fingerprint a caller-supplied resource plus its organism, identifier namespace, scope, version,
   date, license, and citation metadata;
4. collapse exact duplicate `(source, target)` rows, intersect targets with the universe, then apply `min_targets`;
5. call `query_set` once per group with `alternative="greater"`, `n_bg=len(U)`, `ha_corr=0.5`, and the declared
   `min_targets`;
6. independently reconstruct `a:b:c:d` and overlap genes and verify the backend result;
7. retain raw p-values, BH-adjust within each group, and also report a global BH column across all tested pairs;
8. retain all evidence rows and use thresholds only for descriptive candidate flags.

The exact gene identities, not only their count, define the null model. A same-size but different universe is a hard
contract failure. Identifier case-folding, online gene translation, fixed 20,000-gene backgrounds, expression
normalization, and per-cell top-expression sets are outside this node.

## Reporting practice

Results are exploratory Cluster marker evidence for annotation review. A significant positive enrichment is not a
probability of cell identity, Curated annotation, or Sample-level Condition inference. Reports must retain overlap
genes, alternatives and ties, unresolved groups, selection/resource provenance, both multiple-testing families,
software versions, and the resource citation/license.

No default marker CSV is shipped by this repository. A resource must be supplied explicitly and must carry metadata
sufficient for scientific and licensing review.

## Method and software references

- Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics
  data. *Bioinformatics Advances*. 2022;2:vbac016. https://doi.org/10.1093/bioadv/vbac016
- Fisher RA. On the interpretation of chi-square from contingency tables, and the calculation of P.
  *Journal of the Royal Statistical Society*. 1922;85:87-94. https://doi.org/10.2307/2340521
- Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies.
  *Annals of Human Genetics*. 1956;20:309-311. https://doi.org/10.1111/j.1469-1809.1955.tb01285.x
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995;57:289-300.
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound functional enrichment analysis.
  *Genome Biology*. 2015;16:186. https://doi.org/10.1186/s13059-015-0761-7

## Generated-report contract (implementation freeze)

The generated function is frozen as
`run_marker_ora_evidence(marker_table, universe) -> (evidence_dataframe, summary_dict)`. The returned summary is the
complete strict-JSON scientific report, not only diagnostics or `key_results`: it contains the same methods, results,
parameters, warnings, limitations, references, software versions, resource metadata/reference/SHA-256, target
accounting, and bounded group results as the runtime `summary` output. Runtime and generated execution call the same
pure summary builder after the same validated statistical core. ComfyUI timing/history wrappers are deliberately not
part of this scientific payload and therefore do not create a reproducibility exception.
