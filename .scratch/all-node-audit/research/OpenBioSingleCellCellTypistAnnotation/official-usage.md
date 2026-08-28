# CellTypistAnnotation: official usage and primary-source audit

## Audit scope

- Audited: 2026-08-28.
- Repository implementation: `OpenBioSingleCellCellTypistAnnotation` in `openbio_singlecell/nodes_annotation.py`, its schema registration/example, and existing tests.
- Inspected environment: CellTypist 1.7.1 is installed. The repository does not pin CellTypist as a core dependency, so an execution report must disclose the actually imported version rather than assume this audit version.
- Sources are restricted to CellTypist's official documentation, official source repository/source rendering, and the original method paper. API facts below describe the documented 1.7.1 behavior and must be rechecked when the dependency changes.

The domain meaning of this node is **Provisional annotation**: a model-derived label and its audit evidence. It is neither a curated annotation nor biological ground truth, and it must not be conflated with cluster marker evidence or a sample-level Condition contrast.

## Current node behavior

The current public inputs are `adata`, `model`, `use_raw=True`, `majority_voting=True`, `label_column`, and `confidence_column`; the only output is `adata`. Execution copies either `adata.raw` or `adata`, always applies total-count normalization to 10,000 followed by `log1p`, calls `celltypist.annotate`, and inserts one selected label and `conf_score` into `.obs`.

This departs from, or incompletely represents, the official contract in four important ways:

1. In-memory AnnData input is documented as already log1p-normalized to 10,000 counts per cell. Unconditional normalization/logging can double-transform an already valid matrix.
2. The official `annotate` default is `majority_voting=False`; the node defaults to `True` and leaves automatic over-clustering implicit.
3. `AnnotationResult.to_adata()` defaults to confidence insertion by `predicted_labels`. The node then selects `majority_voting` as its label when voting was requested. The written label and written confidence can therefore refer to different predictions.
4. The node drops the probability/decision matrices, model provenance, feature-overlap diagnostics, warnings, software versions, and reproducible code.

## Official inference contract

The documented signature is:

```python
celltypist.annotate(
    filename,
    model=None,
    transpose_input=False,
    gene_file=None,
    cell_file=None,
    mode="best match",
    p_thres=0.5,
    majority_voting=False,
    over_clustering=None,
    use_GPU=False,
    min_prop=0,
)
```

`model` may be a loaded `celltypist.models.Model`, a local model path, or a model name resolvable by CellTypist. In `best match` mode each cell receives the class with the greatest decision score. `prob match` is a separate multi-label mode controlled by `p_thres`; it must not be silently substituted for single-label provisional annotation. The official API and source are the normative references for these details:

- [Official `celltypist.annotate` API](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html)
- [Official `annotate` implementation](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/annotate.html)

### Expression input requirements

For in-memory AnnData, CellTypist documents `.X` (or alternatively `.raw.X`) as log1p-transformed expression normalized to 10,000 counts per cell. Features should be gene symbols, and retaining non-expressed genes is preferred because the classifier intersects query genes with model features. The official classifier source performs only heuristic scale checks on AnnData and can fall back from `.X` to `.raw.X`; it does not establish the upstream provenance of the matrix.

Consequences for this node:

- `use_raw` identifies a container slot, not whether values are raw counts. AnnData `.raw` commonly contains normalized/logged values, despite its name.
- A declared raw-count source may be normalized once with the fixed CellTypist target sum of 10,000 and `log1p`. A source declared already CellTypist-compatible must pass through unchanged. Ambiguous or contradictory state must fail rather than guess.
- Validate a two-dimensional, nonempty, numeric, finite, nonnegative matrix; unique and exactly aligned observation/feature axes; nonzero cell totals when counts are to be normalized; gene-symbol-like feature identifiers; and query/model feature overlap.
- Report the number and fraction of model features matched. Zero overlap is an official hard failure, but merely nonzero overlap can still be scientifically inadequate; a project threshold or prominent low-overlap warning must be documented rather than invented silently.
- Do not silently make duplicate gene identifiers unique, because suffixing can obscure identifier corruption and change feature matching.

Primary official references:

- [CellTypist tutorial: input preparation and model inspection](https://celltypist.readthedocs.io/en/stable/notebook/celltypist_tutorial.html)
- [Official classifier and `AnnotationResult` source](https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html)

### Prediction tables and confidence alignment

An `AnnotationResult` contains:

- `predicted_labels`, whose rows correspond to cells;
- `decision_matrix`, with cells by model classes;
- `probability_matrix`, with cells by model classes;
- optional `majority_voting` labels after over-cluster refinement.

The official `to_adata(insert_conf_by=...)` source defines two distinct confidence mappings:

- `insert_conf_by="predicted_labels"`: the row maximum of the probability matrix;
- `insert_conf_by="majority_voting"`: probability assigned to that voted class, falling back to the row maximum if the voted label is not a model class.

The current node calls `to_adata()` without `insert_conf_by`, so it gets confidence for the individual `predicted_labels`, even when it later writes `majority_voting` as the selected label. This is a correctness defect, not merely missing metadata.

Furthermore, majority voting can emit `Heterogeneous` when `min_prop` is not met. `Heterogeneous` is not a classifier class, so a row maximum is not the probability of that label. The node should expose the meanings separately:

- selected provisional label;
- selected-label probability only when the label is a model class, otherwise null;
- individual best-match label and maximum class probability;
- majority-vote label, when enabled;
- majority-vote support fraction within the over-cluster, which is not a classifier probability.

CellTypist obtains per-class probabilities by applying a sigmoid independently to decision scores. They lie in `[0, 1]` but are not a multiclass softmax distribution and need not sum to one. They are model scores, not calibrated biological certainty. See the official sources:

- [Official `AnnotationResult` API](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html)
- [Official classifier and result conversion source](https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html)
- [Official model prediction source](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/models.html)

For auditability, store a validated probability matrix in a named `.obsm` slot and its exact ordered class labels in `.uns`, rather than creating an unbounded number of `.obs` columns. The decision matrix may be an advanced opt-in because it duplicates a potentially large cell-by-class array. Both matrices must have the exact observation index and class ordering returned by the loaded model, the expected shape, and finite values; probability values must be within `[0, 1]`.

## Majority voting and over-clustering

Majority voting is a refinement over groups of cells, not another independent classifier. Official source behavior is material:

- If no over-clustering is supplied, CellTypist may use existing connectivities or construct a neighbor graph and run over-clustering heuristics. The heuristic resolution varies with cell count and its implementation depends on Scanpy/Leiden versions.
- If the data contain 50 or fewer cells, `annotate(..., majority_voting=True)` warns and returns without the majority-voting step. The current node then expects a missing `majority_voting` column and fails late.
- With `min_prop > 0`, an over-cluster whose dominant individual label does not reach that proportion receives `Heterogeneous`.

For a reproducible scientific node, an explicit categorical over-clustering key is preferred when majority voting is enabled. It keeps graph construction and clustering outside this annotation node and makes the grouping inspectable. Validate exact observation alignment, no missing group memberships, and at least one cell per group. If automatic over-clustering remains supported, disclose that it occurred, the graph source, actual resolution/backend, relevant software versions, and the lack of a CellTypist-level random seed control; never present it as equivalent to a caller-supplied partition.

`majority_voting=False` is the safer public default and agrees with the official API. If explicit voting is requested on 50 or fewer cells, fail before inference with an actionable message rather than silently changing the requested method.

## Determinism

- Base logistic-regression inference is deterministic for an identical matrix, model artifact, model/package implementation, and class ordering.
- Majority aggregation over an explicitly supplied over-clustering vector is deterministic for identical inputs.
- Automatic over-clustering is version/backend/graph dependent and has no `random_seed` argument in the CellTypist `annotate` API. It should not be described as fully deterministic.
- GPU inference is a different execution path and should remain hidden/off unless separately supported, tested, and disclosed.

Reproducibility therefore requires recording the normalized matrix source/transform, resolved model identity and hash, query-feature overlap, package versions, majority-voting mode, grouping provenance, and actual output class order—not just the user-entered model string.

## Model provenance and trust

Official `Model.load` accepts a model object, a path, or a known name and resolves names through CellTypist's local model cache. Model metadata include description fields such as date, details, source, and version, plus cell types and features. A name can resolve to a locally updated artifact, so the report should include:

- requested model identifier;
- resolved local path or stable artifact identifier;
- SHA-256 of the exact model file when path-backed;
- official description/date/details/source/version fields;
- model class count and feature count;
- matched query feature count and fraction;
- CellTypist, scikit-learn, Scanpy, anndata, NumPy, and Python versions used.

Resolve and load the model once, then pass the loaded object to annotation so validation, reported provenance, and inference refer to the same artifact. Execution should not silently download or update model files. A missing model must fail with instructions for an explicit installation step. Because models are Python pickle artifacts, loading one is a code-execution trust boundary: only trusted official or user-authorized model files are acceptable.

References:

- [Official CellTypist model source](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/models.html)
- [Official CellTypist repository](https://github.com/Teichlab/celltypist)

## Required preconditions and failure semantics

All inexpensive checks should complete before model inference, and failure must leave the caller's AnnData unchanged.

Preflight failures should include:

- absent or malformed AnnData, empty axes, duplicate observation or feature identifiers;
- missing selected expression slot, ambiguous expression state, nonnumeric/nonfinite/negative values, or zero-total cells for declared counts;
- invalid or untrusted model artifact, missing model, malformed metadata/classes/features, or inadequate feature overlap under the documented policy;
- invalid output names, conflicting output/metadata keys, or any collision when `overwrite_existing=False`;
- invalid majority parameters, missing/misaligned/noncategorical over-clustering membership, or explicit majority voting on too few cells;
- unsupported mode combinations.

Postconditions should require:

- result rows exactly equal `adata.obs_names`—not a `reindex` that silently inserts missing values;
- complete, nonempty categorical individual predictions from the model class set;
- complete majority labels when voting completed, allowing only model classes plus the documented `Heterogeneous` sentinel;
- probability/decision matrices with exact row index, exact ordered model classes, expected shapes, and valid numeric ranges;
- selected-label confidence derived explicitly from the selected label semantics;
- complete group support metrics when majority voting is used;
- no partial columns, `.obsm`, or `.uns` state after any error.

Errors should name the operation, selected expression source, requested/resolved model where safe, installed CellTypist version, expected invariant, and observed value. Upstream dependency warnings that change semantics—especially the `<=50` majority-voting fallback—must be converted to a preflight error or an explicit structured warning, not lost.

## General scientific practice and reportable caveats

The original CellTypist study trained logistic-regression models on curated, harmonized immune-cell references spanning multiple studies and tissues and reported strong performance within its evaluated settings. It also found that representation of a cell type in training data was a major determinant of prediction accuracy. These results support model-assisted annotation, not universal validity for every tissue, species, protocol, disease state, or model version.

The node's report must state:

- predictions are provisional and require marker/context review before curation;
- the class repertoire is closed, so unseen states, doublets, low-quality cells, and out-of-domain cells may be forced to a known best match;
- model probability is not biological truth or guaranteed calibration;
- majority voting smooths within the supplied/derived groups and can erase rare, transitional, or mixed states;
- batch robustness reported in the paper is context-specific. Tissue, species, protocol, disease, Technical batch, and training imbalance can shift results;
- label frequencies and low-score rates may be summarized descriptively by Sample or Technical batch, but cell-level labels must not be used as independent replicates for Condition inference. Formal Condition contrast remains Sample-level.

Primary method paper:

- Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*. 2022;376:eabl5197. [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197); [PubMed 35549406](https://pubmed.ncbi.nlm.nih.gov/35549406/).
- [Official CellTypist project portal and publication list](https://www.celltypist.org/)

## Minimum report disclosure

An analysis summary should be strict JSON and contain at least:

- method name, exact public and resolved parameters, expression source and transformation performed;
- model provenance/fingerprint and query-feature overlap;
- individual/selected/majority label semantics and whether majority voting actually completed;
- per-label cell counts, missing/heterogeneous counts, individual-versus-majority agreement, selected-score distribution, and majority-support distribution when applicable;
- probability/decision storage keys, shapes, and exact class ordering;
- descriptive Sample/Technical-batch diagnostics when such metadata were explicitly supplied;
- warnings and the scientific limitations above;
- software versions and the official API/source/method-paper references listed in this document.

The equivalent source-code output must reproduce the same model resolution and fingerprint check, expression validation/one-time transformation, explicit annotation parameters, result alignment checks, confidence mapping, canonical `.obs`/`.obsm`/`.uns` writes, summary calculations, and failures.

## Final implementation clarifications

- The generated `run_celltypist_annotation(adata)` function returns `(output_adata, summary_dict)`. Runtime and generated execution call the same pure summary builder, so the complete scientific payload is equal: methods, results, key results, resolved parameters, warnings, limitations, references, and software versions. Wall-clock timing and plugin history records are wrapper metadata and are intentionally outside this strict summary payload.
- Dependency versions are collected once by the standalone core, stored in canonical `uns[metadata_key]["software_versions"]`, carried in diagnostics, and consumed unchanged by the shared summary builder. The AnnData provenance and returned summary therefore cannot disagree about the software environment for one execution.
- Validation of `verified_cp10k_log1p` follows the CellTypist 1.7.1 scale boundary per cell: every inverse-`log1p` row total must be within an absolute 1 count of 10,000 (`rtol=0`, `atol=1`). A total such as 9,950 is not accepted as CP10K. Genuine float32 CP10K matrices remain valid within that absolute numerical tolerance.
- Individual and selected label count maps cover every ordered model class with integer counts, including zero. `Heterogeneous` is added only when it is actually returned by majority voting.
- `overwrite_existing=True` is ownership-aware. Existing outputs may be replaced only when `uns[metadata_key]` is a schema-version-1 `celltypist_annotation` record whose `columns` and `matrices` explicitly own the colliding artifacts. On a verified rerun, prior majority columns, decision matrices, renamed outputs, and a renamed selected-label annotation entry are removed atomically when the new run no longer produces them. Missing, malformed, or conflicting ownership fails without deleting caller data.
- Cache invalidation fingerprints the resolved local model using canonical path, byte size, nanosecond modification time, and SHA-256. Explicit paths and official locally resolved model names use the same artifact fingerprint. Fingerprinting and execution use only an existing path or the official `get_model_path`; they never enumerate, download, or update models.
- A missing model error names the requested identifier and instructs the user to either provide a trusted local `.pkl` or explicitly run the network-enabled official command `celltypist.models.download_models(model='MODEL.pkl')` with their actual model identifier. The node itself never performs that command.
- CellTypist model files are Python pickle artifacts and can execute code while loading. Reports and load errors therefore instruct users to accept only official or explicitly trusted sources. SHA-256 identifies exact bytes for reproducibility; it does not establish that a model is safe.
