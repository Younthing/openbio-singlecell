# OpenBioSingleCellGSVAScores — official usage research

## Final open-expert boundary

Gaussian and empirical kernels accept any finite explicitly selected source with scale warnings. Poisson keeps only
its arithmetic domain—finite nonnegative integer values—as a hard requirement; X/Raw/layer location and history are
not credentials. No Raw/current binding is required or tested.

## Audited baseline

The current node calls removed decoupler 1.x `run_gsva`, forces the entire selected matrix dense, defaults
`abs_rnk=True`, exposes a seed that GSVA 2.x does not use, and expects `obsm["gsva_estimate"]`. With the supported
decoupler 2.x range it is non-executable and its default differs from the official method default.

The implementation target reviewed here is decoupler 2.2.0; the base environment locks Scanpy 1.12.3 and AnnData
0.13.2. Decoupler remains optional and its exact installed version must be checked and reported.

## Official decoupler 2.2 usage

Primary sources:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.gsva.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_gsva.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/_Method.py
- https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html

```python
dc.mt.gsva(
    data=work,
    net=network.rename(columns={source_column: "source", target_column: "target"}),
    tmin=min_targets,
    layer=None,
    raw=False,
    empty=False,
    bsize=batch_size,
    verbose=False,
    kcdf="gaussian",
    maxdiff=True,
    absrnk=False,
    tau=1.0,
)
scores = work.obsm["score_gsva"]
```

In 2.2.0 `kcdf` is `"gaussian"`, `"poisson"`, or `None`; the documented default is Gaussian. The implementation
parameter is named `maxdiff`, while the prose/API page also describes it as `mx_diff`. `absrnk=False` is the official
default. GSVA returns an observation-by-set score frame and does not perform a hypothesis test or return p-values.

The implementation converts sparse batches to dense before density transformation. Gaussian GSVA is recommended for
normalized continuous expression and Poisson requires finite nonnegative integer values. The public interface does
not itself prohibit an expert from explicitly applying Gaussian or empirical kernels to count-like input, so that
choice remains executable but is prominently disclosed as nonrecommended. A Poisson request on incompatible values
remains a hard method-definition error. Because single-cell matrices can be very large, a preflight dense-byte
estimate and a hard caller-visible memory limit are required.

## Scientific-use findings

GSVA yields relative, cohort-dependent observation-level gene-set scores: the density transformation uses the
observations included in the run. Subsetting or adding cells/Samples can change scores. Scores are descriptive and
are not p-values. Cells from one Sample cannot be treated as independent biological replicates in a Condition test.
Downstream inference must use Sample-level summaries/models within a defined population, with Technical batch
handled as a nuisance variable where supported.

The exact expression feature universe and gene-set resource matter. No identifier case folding, synonym expansion,
organism inference, or online resource fetch belongs in this node. Resource name/version/date/organism/namespace/
scope/license/citation and SHA-256 must accompany every run. Intersected target counts and filtered sets must be
reported.

## References to emit

- Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data.
  *BMC Bioinformatics*. 2013;14:7. https://doi.org/10.1186/1471-2105-14-7
- Badia-i-Mompel P, et al. decoupleR. *Bioinformatics Advances*. 2022;2:vbac016.
  https://doi.org/10.1093/bioadv/vbac016
- Virshup I, et al. anndata: Annotated data. *JOSS*. 2024;9:4371. https://doi.org/10.21105/joss.04371

## Report/code implications

The report must disclose kernel, `maxdiff`, `absrnk`, `tau`, input state, observation cohort, dense memory estimate,
feature/set/edge accounting, per-set score distributions, lack of statistical testing, resource metadata/hash,
warnings, limitations, references, and exact versions. Generated code must use the same 2.2 public call, perform
the same memory preflight and postconditions, avoid downloads, and return `(output_adata, summary_dict)`.
