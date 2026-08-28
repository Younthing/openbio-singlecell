# OpenBioSingleCellCollecTRIULM — official usage research

## Audited baseline

The optional dependency contract is currently `decoupler>=2.1.4,<3`; this audit targets decoupler 2.2.0, the
current release on 2026-08-28. The node still uses the removed decoupler 1.x calls `get_collectri` and `run_ulm`,
passes the removed `split_complexes` argument, and expects `obsm["ulm_estimate"]`/`obsm["ulm_pvals"]`. It is
non-executable with the supported 2.x package family.

Decoupler 2.2 exposes CollecTRI through `decoupler.op.collectri` and ULM through `decoupler.mt.ulm`. The official
loader downloads the publication snapshot from Zenodo and returns a long table with `source`, `target`, `weight`,
supporting resources, and references. The affiliation-dependent `license` argument is scientifically and legally
material and cannot remain hidden.

### Versioned-source corrections discovered during implementation

The 2.2.0 API page documents `license={"academic","commercial","nonprofit"}`, but the versioned
`decoupler.op.collectri` function body does not read that argument after binding it: all three values reach the same
fixed Zenodo URL and no affiliation-specific filtering is applied. The adapter must still validate, forward, and
report the caller's affiliation policy for forward compatibility, but it must not claim that decoupler 2.2 enforced
or transformed the network under that policy. The Zenodo record also does not expose that affiliation choice as the
dataset file's license.

The 2.2.0 `_download` implementation performs a direct HTTP download on every loader call and has no package cache
contract. Consequently official-loader execution with `allow_network_access=false` must fail before calling the
loader; an allowed call is reported as network access, never as a cache hit. A live 2.2.0 human load on the audit
date returned 42,990 unique `source,target` edges with columns `source,target,weight,resources,references,
sign_decision`, weights in `{-1,+1}`, and retained `AP1`/`NFKB` complexes. The fixed download contained 43,175 raw
rows; the loader's unqualified `dropna()` removed 185 rows with missing reference fields before the complex policy.
The downloaded file was 4,345,649 bytes with SHA-256
`4473c9189dd53dacc80297709ad1452dda1086a1cc2185f9a56146c261668701` (Zenodo reports MD5
`cee4a3943c059e6dd8796ced6dec44f6`). Removing AP1 and NFKB removed 1,395 more edges. These observed counts and hashes
are snapshot smoke-test facts, not API guarantees; the runtime canonical content hash remains authoritative.

## Official decoupler 2.2 usage

Versioned official interfaces:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.op.collectri.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/op/_collectri.py
- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.ulm.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_ulm.py
- https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html

The audited public call is:

```python
import decoupler as dc

network = dc.op.collectri(
    organism="human",
    remove_complexes=False,
    license="academic",
    verbose=False,
)
dc.mt.ulm(
    data=work,
    net=network,
    tmin=5,
    raw=False,
    empty=False,
    bsize=250_000,
    verbose=False,
    tval=True,
)
scores = work.obsm["score_ulm"]
adjusted_pvalues = work.obsm["padj_ulm"]
```

`remove_complexes=False` retains the AP1 and NFKB complexes; it is not equivalent to the legacy
`split_complexes` flag. ULM fits, for every observation and regulator, a univariate linear model across measured
features using signed CollecTRI weights. The regression uses the complete selected feature axis: network targets
carry their signed weights and every non-target feature carries weight zero, so the residual degrees of freedom are
`n_features - 2` and background-feature selection can change both scores and p-values. With `tval=True`, the score
is the t statistic for the fitted slope: positive means activity consistent with the regulon, negative means
repression/inactivity, and values near zero are inconclusive. The second matrix is Benjamini-Hochberg-adjusted
p-values (`padj_ulm`), not raw p-values. Adjustment is performed separately for each observation across that
observation's complete retained-regulator family. Reproducibility therefore includes a fingerprint of the complete
ordered expression-feature axis, not only the intersecting network targets.

Decoupler documents normalized continuous expression (for example library-size normalization followed by log1p)
for observational counts. Its 2.2 API also explicitly permits feature scaling after normalization, while warning
that scaling changes feature importance and makes results depend on the included observations. OpenBio keeps
normalized expression as the default and reports the consequence of count-like, residual, scaled, or unknown input,
but finite aligned data are not rejected merely because mutable history assigns a nonrecommended state. The selected
source is the expert user's declaration; history is descriptive evidence, not an authorization credential.

## Resource and scientific-use findings

`dc.op.collectri` retrieves the publication resource from Zenodo record 8192729, whose file-level license is
CC BY 4.0; the CollecTRI project separately directs users to the original licenses of incorporated resources. The
loader removes every row containing any missing field, normalizes resource/reference strings, and drops duplicate
source-target pairs. For non-human organisms it also makes a second uncached HCOP request with hidden translation
defaults, translates both source and target, and keeps the first translated duplicate. Official-loader mode in this
node is therefore restricted to the directly published human resource; a non-human run must provide a separately
licensed, fingerprinted, pre-translated local snapshot. A local snapshot must declare whether its bytes represent
the raw Zenodo table or a post-decoupler table. Reproducibility therefore requires the
runtime network's canonical content SHA-256, row/source/target counts, organism, complex policy, license selection,
and decoupler version; a package version or URL alone is insufficient. Here `license selection` is explicitly a
caller affiliation declaration, not proof of 2.2 filtering or a dataset-license assertion. Execution must not
silently download a resource in an offline workflow. A local, pre-materialized network with strict metadata is the
reproducible path; an official-loader mode must disclose that network access occurred and that decoupler 2.2 has no
loader cache contract.

Scores are per cell and exploratory. Cells from the same Sample are not biological replicates. Per-cell adjusted
p-values only test the within-cell feature-weight regression and do not provide a Condition contrast, correct a
Technical batch, or establish TF causality. Any Condition inference must be a separate Sample-level analysis within
a defined population.

## References to emit

- Müller-Dott S, et al. Expanding the coverage of regulons from high-confidence prior knowledge for accurate
  estimation of transcription factor activities. *Nucleic Acids Research*. 2023;51:10934-10949.
  https://doi.org/10.1093/nar/gkad841
- Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics
  data. *Bioinformatics Advances*. 2022;2:vbac016. https://doi.org/10.1093/bioadv/vbac016
- Virshup I, et al. anndata: Annotated data. *Journal of Open Source Software*. 2024;9:4371.
  https://doi.org/10.21105/joss.04371

## Report and generated-code implications

The report must disclose ULM's score definition and adjusted-p-value family, expression state, input dimensions,
retained regulators/targets/edges after exact feature intersection, `tmin`, complex and license policies, resource
SHA-256/metadata, score and adjusted-p-value ranges, bounded strongest positive/negative activities, warnings,
references, and dynamic software versions. Software disclosure also includes Numba because 2.2 performs its BH
step through Numba (with float32 output), plus Requests when official loading is used. Equivalent code must call the
2.2 interfaces above, verify the complete feature-axis hash and canonical
`score_ulm`/`padj_ulm` frames, perform no download when using a local resource, and return the same result and strict
JSON summary.
