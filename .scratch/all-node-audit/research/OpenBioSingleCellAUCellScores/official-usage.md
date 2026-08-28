# OpenBioSingleCellAUCellScores — official usage research

## Final open-expert boundary

The explicit X/Raw/layer selection is sufficient. Finite count, scaled, residual, transformed, and unknown-state
inputs remain executable; history-derived state changes warnings and report wording only. No Raw/current binding is
required or tested.

## Audited baseline

The current node calls the removed decoupler 1.x `run_aucell` interface and expects
`obsm["aucell_estimate"]`. The optional dependency range is `decoupler>=2.1.4,<3`, and the implementation target
validated for this audit is decoupler 2.2.0. In 2.2.0 the public method is `decoupler.mt.aucell`; an AnnData input is
modified at `obsm["score_aucell"]`. The old call is therefore non-executable with the supported package family.

The base environment locks Scanpy 1.12.3 and AnnData 0.13.2. Decoupler is optional and is not present in the base
lock, so execution must check the installed version and public interface and report the exact runtime version.

## Official decoupler 2.2 usage

Official API and versioned implementation:

- https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.aucell.html
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_aucell.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_run.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/data.py
- https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/_Method.py
- https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html

The public call shape is a `Method` interface:

```python
import decoupler as dc

dc.mt.aucell(
    data=work,
    net=network.rename(columns={source_column: "source", target_column: "target"}),
    tmin=min_targets,
    raw=False,
    empty=False,
    bsize=batch_size,
    verbose=False,
    n_up=n_top_features,
)
scores = work.obsm["score_aucell"]
```

`network` is long-form with unique `source,target` pairs. `tmin` is applied after overlap with the measured feature
universe. `n_up=None` selects the top 5% of measured features, clipped to the valid range. AUCell ranks each
observation and measures how early a set's targets appear. It returns a score in `[0,1]`; it performs no hypothesis
test and returns no p-value.

The complete public call path matters for tie handling. `Method.__call__` enters `mt._run._run`, which calls
`pp.data.extract(..., shuffle=True)`; `_break_ties` applies a deterministic feature permutation from
`numpy.random.default_rng(seed=0)` before the inner AUCell function uses
`scipy.stats.rankdata(..., method="ordinal")`. Thus ties follow the fixed seed-0 permuted feature order. The seed is
an internal Decoupler policy, not a user parameter, and the node must disclose it without exposing an inert
`random_seed` control. The exact measured feature order and universe remain reproducibility inputs because the
fixed permutation is applied to their incoming positions.

## Data-state and scientific-use findings

Decoupler's method documentation recommends normalized continuous data (for example library-size normalization
followed by log1p), but AUCell itself ranks features within each observation and the public 2.2 interface accepts a
numeric matrix rather than enforcing one expression state. Experts may therefore explicitly select `X`, a named
layer, or the post-QC Raw snapshot. Raw/count-like input is disclosed prominently because sparsity and ties can change
rank recovery; it is not rejected solely as a practice preference. The adapter still rejects non-finite values,
duplicate feature identifiers, a zero-feature or zero-cell input, all-zero observations whose complete tie rank would
be arbitrary, and any output that changes observation order.

AUCell scores are per observation and descriptive. Cells from the same Sample are not biological replicates, and
the score node must not perform or imply Condition inference. A downstream Condition contrast must aggregate or
model scores at the Sample level within a defined population. Technical batch is a nuisance variable for that
downstream analysis; AUCell does not correct it.

The gene-set resource determines biological meaning. It must be caller-supplied, organism/identifier compatible,
versioned, licensed, cited, byte-fingerprinted, and intersected exactly with the expression feature universe. A
bundled unversioned default is not adequate for a scientific report.

## References to emit

- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463
- Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics
  data. *Bioinformatics Advances*. 2022;2:vbac016. https://doi.org/10.1093/bioadv/vbac016
- Virshup I, et al. anndata: Annotated data. *Journal of Open Source Software*. 2024;9:4371.
  https://doi.org/10.21105/joss.04371

## Report and generated-code implications

The report must name AUCell, state explicitly that there is no statistical test, identify expression state, exact
feature and observation counts, `n_up`, retained/excluded set counts, per-set in-universe target counts, score range
and bounded top/bottom summaries, resource metadata/SHA-256, tie policy, warnings, references, and dynamic package
versions. Equivalent generated code must validate the same resource and state, call `dc.mt.aucell`, verify the
canonical score frame, and return `(output_adata, summary_dict)` without downloading data.
