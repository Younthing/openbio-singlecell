# OpenBioSingleCellGeneSetOverrepresentation — official usage research

## Audited baseline

The current node combines generic-table filtering, optional group selection, background construction from unrelated
AnnData `X`/`raw` feature names, resource parsing, hypergeometric testing, and reporting. It cannot prove that the
AnnData features are the genes actually tested upstream, permits selected genes outside that test family, and treats
an arbitrary `TableResult` as differential evidence. This can produce a wrong null universe even when universe size
looks plausible.

## Official explicit-set interface

For a preselected feature set and an exact tested universe, decoupler 2.2 provides `mt.query_set`:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.query_set.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_query_set.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/net.py

The versioned public signature is
`query_set(features, net, alternative="greater", n_bg=20000, ha_corr=0.5, tmin=5, verbose=False)`. Correct usage is:

```python
result = dc.mt.query_set(
    features=selected_genes,
    net=network,
    alternative="greater",
    n_bg=len(tested_universe),
    ha_corr=0.5,
    tmin=min_targets,
    verbose=False,
)
```

The OpenBio adapter must independently reconstruct the in-universe `a,b,c,d` contingency table, overlap genes,
one-sided Fisher p-value, and Haldane-corrected log odds ratio and verify backend results. BH is applied across the
complete tested resource family for the one selected set. Thresholds flag candidates but never delete rows.

## Scientific-use findings

The selected set should be produced upstream by an explicit filter; this node must not decide what “significant”
means. For generic tables the filter identity and approval flag are producer/caller declarations, not programmatic
proof. Their absence or inconsistency is reported as a caution rather than blocking an expert calculation. The exact
gene identities in the tested universe, not AnnData's current features and not a fixed 20,000 background, define the
null. Resource targets are intersected with that universe before `min_targets`. Changing same-sized universe genes
changes the computation and must fail fingerprint checks.

If the selected genes derive from Cluster marker evidence, results are exploratory annotation evidence. If they
derive from a formal Condition contrast, the upstream model must use Sample as inference unit and Technical batch as
nuisance. ORA cannot repair pseudoreplication or confounding upstream. Identifier translation, case folding, resource
downloads, and label assignment are outside this node.

Scope, inference-unit, replicate-awareness, and batch-handling metadata on a generic artifact remain unverified
declarations. Missing, unfamiliar, cell-level, or table/universe-inconsistent declarations must be surfaced in the
strict summary as limiting formal interpretation, not treated as evidence that the ORA arithmetic itself is invalid.

## References to emit

- Fisher RA. On the interpretation of chi-square from contingency tables. *JRSS*. 1922;85:87-94.
  https://doi.org/10.2307/2340521
- Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies.
  *Annals of Human Genetics*. 1956;20:309-311. https://doi.org/10.1111/j.1469-1809.1955.tb01285.x
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995;57:289-300.
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Timmons JA, et al. Multiple sources of bias confound functional enrichment analysis. *Genome Biology*.
  2015;16:186. https://doi.org/10.1186/s13059-015-0761-7
- Badia-i-Mompel P, et al. decoupleR. *Bioinformatics Advances*. 2022;2:vbac016.
  https://doi.org/10.1093/bioadv/vbac016

## Report/code implications

Report upstream evidence scope/design and selection provenance explicitly as caller/producer declarations, exact
universe identity, all contingency counts,
overlaps/effects/raw and adjusted p, tested family, resource metadata/hash/license, alternatives/ties, limitations,
references, and exact versions. Generated code returns `(evidence_dataframe, summary_dict)` and repeats all artifact,
resource, contingency, backend, and strict-JSON validation.
