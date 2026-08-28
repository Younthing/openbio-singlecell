# Normalize Total: official usage research

Researched: 2026-08-28

## Official Scanpy interface and example

Scanpy documents `scanpy.pp.normalize_total(adata, *, target_sum=None, exclude_highly_expressed=False, max_fraction=0.05, key_added=None, layer=None, obsm=None, inplace=True, copy=False)`. The operation rescales every observation by its total expression so that, when all genes contribute to the size factor, each nonempty cell has the requested total. `target_sum=1e6` is CPM; `target_sum=None` uses a dataset-derived median total.

The official example obtains a returned normalized matrix without mutating its input:

```python
X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)["X"]
```

The stable Scanpy clustering tutorial shows the common RNA workflow as a deliberate sequence: save counts, run total-count normalization, then run `log1p`. In this repository, the first action belongs to `OpenBioSingleCellSnapshotExpression`; normalization must not reproduce that storage side effect.

- Scanpy `normalize_total` API and parameter semantics: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
- Scanpy official preprocessing tutorial: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html#normalization
- Scanpy preprocessing index: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html

`exclude_highly_expressed=True` excludes a gene from size-factor calculation when it exceeds `max_fraction` of a cell's original total in at least one cell. The excluded genes are still rescaled by the resulting factor, so the sum over *all* genes can exceed `target_sum`; only the not-excluded genes are constrained to that target. This choice therefore changes the scientific result and must never be inferred silently. `key_added`, `inplace`, `copy`, `layer`, and `obsm` are storage/execution controls in the upstream interface, not all of which belong in the node interface.

## Scientific input state and practice

The usual input is a declared count-like expression state, normally post-QC UMI counts in `X` or the canonical `layers["counts"]`. That is a scientific recommendation, not something the numeric values or plugin history can prove. Scanpy can execute on finite signed or non-integer matrices, and it leaves zero-total rows unchanged while warning. An expert-facing wrapper therefore hard-validates only the aligned finite numeric matrix and the explicitly positive target; it warns and discloses non-integer values, negative values/totals, zero-total rows, all-zero input, and provenance inconsistent with counts. The report must never upgrade such input to “raw UMI counts.”

Total-count normalization corrects per-cell library depth by a single scalar. It is not Technical batch integration, variance stabilization, or a Sample-level inferential model. Its output is a derived expression representation and must not be supplied to count-dependent scVI, pseudobulk, or count-model differential-expression interfaces. Those consumers should continue to select the Raw snapshot/count layer explicitly.

The target total is an analysis choice. A positive explicit default such as 10,000 is reproducible and avoids using `0` as a sentinel for Scanpy's dataset-dependent `None` behavior. If median-depth normalization is ever supported, it should be represented by an explicit mode rather than an overloaded numeric value, and the resolved median must be reported.

## Edge conditions and failure behavior

- Reject zero axes and non-finite/non-numeric or misaligned matrices. Negative values and zero/non-positive cell totals are warnings because the public backend still defines a result; zero-total rows remain zero and cannot achieve the requested target, which must be explicit in `summary`.
- Require a finite `target_sum > 0`; do not translate `0`, falsey values, or NaN to `None`.
- If high-expression exclusion is exposed, require `0 < max_fraction <= 1` and report the genes excluded from size-factor estimation.
- Reject a missing/blank selected layer and backed AnnData with an actionable `to_memory()` instruction; the node's copy-on-write contract is an in-memory operation.
- Preserve observation/variable order, annotations, `raw`, and every layer. In particular, do not create or overwrite `layers["counts"]`.
- Accept dense and supported sparse matrices without densifying merely to normalize or summarize them.
- Warn when the selected matrix is non-integer-like, contains negative values/non-positive totals, is all zero, or when provenance is absent or indicates a transformed state.

## Methods and software references

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Zheng GXY et al. Massively parallel digital transcriptional profiling of single cells. *Nature Communications*. 2017;8:14049. https://doi.org/10.1038/ncomms14049
- Weinreb C, Wolock S, Klein AM. SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. *Bioinformatics*. 2018;34(7):1246-1248. https://doi.org/10.1093/bioinformatics/btx792
- Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. *Nature Methods*. 2023;20:665-672. https://doi.org/10.1038/s41592-023-01814-1

The Scanpy software paper is the required software citation. The Zheng and Weinreb papers are methods/examples named by the official Scanpy documentation; the transformation comparison supports the limitation that depth scaling is one preprocessing representation rather than a universally sufficient expression model.

## Required `summary` and `code` disclosure

`summary` should identify the selected count source, input dimensions and storage type, integer-like status, per-cell total distribution before normalization, target total, high-expression-exclusion policy, output per-cell total distribution, zero-total-cell count, preservation of existing `raw`/layers, and all warnings or limitations. The methods sentence must state that one per-cell scalar was used; the results sentence should report the number of cells normalized and the achieved totals without implying Technical batch correction. References must include Scanpy and the relevant method/practice sources. Runtime versions should include Python, openbio-singlecell, Scanpy, AnnData, NumPy, and SciPy (plus Pandas when used in reporting).

`code` should define a standalone function that resolves the same `X`/layer source, validates the same invariants, copies the AnnData, applies the resolved Scanpy parameters, writes only the normalized matrix to output `X`, and preserves all pre-existing layers and `raw`. Plugin-only history bookkeeping is intentionally omitted.
