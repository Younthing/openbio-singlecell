# CellTypistAnnotation: pre-change design review

## Decision

**Retain and enhance this node as a separate atomic analysis module. Do not merge it with marker testing, enrichment, curated-label mapping, clustering, or Condition inference. Do not delete it.**

Its single scientific question is: *given one validated expression source and one pinned CellTypist model, what provisional per-cell annotations and auditable model evidence does that model produce?* Marker evidence and expert curation may assess those predictions later, but they are different domain operations.

The deletion test supports retention: removing the module would force workflows to redistribute nontrivial expression-state validation, trusted model resolution, feature matching, label/probability alignment, majority-vote semantics, provenance, and output validation across callers. The module has useful depth once those concerns are hidden behind a small interface.

## Current depth and locality problems

The present interface is small, but its implementation omits rather than hides complexity:

- `use_raw` couples a storage slot to an assumed expression state, then unconditional normalization/logging can double-transform data.
- the default enables a compound hidden analysis—automatic graph/over-clustering and majority voting—inside an annotation node;
- the selected majority-vote label can be paired with the confidence for the individual prediction;
- result matrices, model identity, feature overlap, dependency versions, and scientific limitations disappear;
- late missing-column/index failures leak CellTypist internals instead of enforcing module postconditions;
- the caller receives no structured report or equivalent code.

This gives low information hiding and poor locality: users must know CellTypist source behavior to interpret two apparently simple `.obs` columns.

## Proposed atomic seam

Make the node a deep module whose public interface expresses scientific choices while an internal CellTypist adapter owns dependency mechanics.

Conceptual execution seam:

```text
validate names/collisions and expression declaration
  -> resolve and fingerprint trusted model
  -> validate/prepare one expression matrix exactly once
  -> run CellTypist through a narrow adapter
  -> validate cell/class axes and numeric postconditions
  -> derive explicitly named individual/majority evidence
  -> commit canonical AnnData outputs atomically
  -> emit strict summary JSON and equivalent source code
```

The adapter is justified because CellTypist is a true external dependency. A production adapter and a deterministic fake make source/result edge cases testable without turning the adapter into a public extension socket. Matrix storage and reporting stay in-process; adding pluggable sinks would create leverage debt without a present external boundary.

## Public interface

Recommended outputs:

```text
adata, summary, code
```

Recommended inputs, expressed conceptually rather than as a commitment to exact schema spelling:

| Input | Visibility | Design |
|---|---|---|
| `adata` | primary | Required AnnData; input must remain unchanged on failure. |
| expression source | primary | Explicit `X`, Raw snapshot, or named layer. `use_raw` remains only as a migration alias. |
| expression state | primary | Required choice between verified raw counts and verified CellTypist-compatible CP10K/log1p; ambiguous state fails. |
| `model` | primary | Required trusted model name/path/artifact. There is no scientifically universal default model. |
| `majority_voting` | primary | Default `False`, matching the official API; changes the scientific output. |
| `over_clustering_key` | advanced | Preferred explicit categorical grouping when voting is enabled. Automatic over-clustering, if retained, must be an explicit mode. |
| `min_prop` | advanced | Expose only with majority voting; validate `[0, 1]`; produces `Heterogeneous` below the threshold. |
| selected label/confidence column names | advanced | Keep for downstream compatibility and validate names/collisions. |
| probability storage key | advanced | Canonical default; enable by default for audit evidence. |
| store decision matrix | advanced | Default off because it duplicates a cell-by-class matrix. |
| output metadata key | advanced | Canonical default such as `celltypist`; validate collisions. |
| `overwrite_existing` | advanced | Default `False`; preflight every prospective `.obs`, `.obsm`, and `.uns` collision. |

The fixed normalization target of 10,000 and `log1p`, result-axis rules, trusted pickle policy, and matrix postconditions are hidden policy rather than exposed tuning knobs. `use_GPU`, transpose/file parsing, and `prob match` should remain hidden for this atomic single-label node; supporting multi-label annotation would be a separate public capability with different output semantics.

The current `model=""` effectively requires an explicit model. Preserve that good behavior rather than silently adopting the library's immune-model default.

## Canonical primary output

The returned AnnData should be a copy with a canonical, versioned record in `.uns[metadata_key]`. Suggested evidence layout:

- selected provisional label: categorical `.obs[label_column]`;
- selected-label probability: numeric nullable `.obs[confidence_column]`, null when the selected label is not a model class;
- individual best-match label and row-maximum class probability: canonical `.obs` fields;
- when enabled, majority-vote label and majority support fraction: distinct canonical `.obs` fields;
- probability matrix: `.obsm[probability_key]`, cells by exact ordered model classes;
- optional decision matrix in another `.obsm` slot;
- exact ordered classes, model provenance/hash, feature overlap, expression provenance, parameters, dependency versions, and warning/limitation flags in `.uns[metadata_key]`.

Do not call the majority support fraction a confidence score. Do not treat the row maximum as probability of `Heterogeneous`. Avoid one `.obs` column per class, which makes interface size scale with model vocabulary.

The primary result must remain usable by existing downstream nodes that accept an annotation key, while its metadata identifies the field as Provisional annotation rather than Curated annotation.

## Structured `summary` contract

`summary` is strict, standards-compliant JSON with no NaN/Infinity, pandas/NumPy scalars, or non-string object keys. It should include:

- `schema_version`, method name, run status, resolved parameters, and actual expression transform;
- references: official API, official implementation/model sources, repository, and original CellTypist paper;
- model: requested/resolved identifier, SHA-256, description/date/details/source/version, class/feature counts, and feature overlap;
- results: cell count, per-label complete count map, missing/heterogeneous count, score quantiles/range, and probability matrix shape/class order;
- majority diagnostics when applicable: grouping source, whether voting completed, group count/sizes, majority-support distribution, and individual-versus-majority agreement;
- optional descriptive label/low-score distributions by explicitly selected Sample or Technical-batch metadata;
- warnings/limitations that say “provisional”, avoid “ground truth”/“validated”, and state closed-vocabulary and batch/domain-shift caveats;
- software: Python, CellTypist, scikit-learn, Scanpy, anndata, NumPy, and relevant backend versions.

The summary must not invent a best biological label, assert calibration, perform cell-level Condition inference, or call a model prediction curated truth.

## Equivalent `code` contract

The text output should define one importable function with explicit arguments and return the same scientifically meaningful AnnData plus the same summary object. Equivalence requires more than repeating `celltypist.annotate`:

1. validate axes, selected expression source/state, values, collisions, majority settings, and trusted model before inference;
2. resolve/load one model, verify its fingerprint and metadata, and calculate query/model feature overlap;
3. transform verified counts exactly once or preserve verified CP10K/log1p values;
4. call CellTypist with explicit `mode`, voting, grouping, and `min_prop` parameters;
5. enforce identical observation and class axes, label completeness, matrix shape/finite/range invariants, and selected confidence semantics;
6. commit the same canonical `.obs`, `.obsm`, and `.uns` structure atomically;
7. construct strict JSON with the same result diagnostics, provenance, versions, references, and caveats.

Runtime and generated code must reject the same malformed backend results. Truly incidental plugin history or wall-clock timing may differ; scientific output and disclosed resolved policy may not.

## Diagnostics and failure design

### Before invoking CellTypist

- Validate the full prospective output key set and reject all collisions unless overwrite was explicitly authorized.
- Validate nonempty unique observation and feature identifiers and exact source-axis alignment.
- Validate numeric, finite, nonnegative expression and state-specific conditions; reject ambiguous/double-transformed inputs.
- Load the exact trusted model and validate its classes/features/metadata; calculate feature overlap under a documented warning/error policy.
- Validate majority-voting combinations, `min_prop`, categorical over-clustering membership and exact observation alignment. Explicit voting on 50 or fewer cells should fail early because official CellTypist otherwise returns without voting.

No backend work should occur after any preflight error.

### After invoking CellTypist

- Require prediction rows to equal `adata.obs_names` exactly. Do not silently repair with `reindex`.
- Require complete individual labels from the model class set; majority labels must be model classes or the documented `Heterogeneous` sentinel.
- Require exact cell-by-class probability and decision matrices, exact class order, finite decision values, and probabilities in `[0, 1]`.
- Derive selected-label probability by explicit lookup, not by positional maximum when the selected label differs from the individual prediction.
- Require complete majority support/group metrics when voting ran.
- Build all changes on a copy and publish only after every invariant passes.

Every error should identify the operation, invariant, observed value, relevant key/model, and installed CellTypist version. Warnings that alter requested semantics are structured report data or fatal preconditions, not console-only messages.

## Scientific boundary and batch diagnostics

The node may summarize output label distribution, score distribution, feature overlap, and disagreement. If the caller explicitly identifies Sample or Technical-batch metadata, it may report descriptive differences in label and low-score frequencies to expose possible domain shift. It must not:

- infer that model labels are validated cell identities;
- merge technical groups or correct expression;
- perform marker validation or curated mapping;
- test Condition differences using cells as independent replicates;
- choose labels because they make downstream biology appear cleaner.

The report must note that model coverage, tissue/species/protocol/disease mismatch, class imbalance, low-quality cells, doublets, and novel/transitional states can affect predictions. Majority voting can suppress rare or mixed states and depends on its grouping.

## Retain/merge/delete assessment

- **Retain:** yes. A pinned pretrained-classifier application is a coherent atomic analysis and has enough hidden complexity to warrant a deep module.
- **Enhance:** yes. Fix expression-state handling and label/confidence alignment; add probability evidence, model/input provenance, postconditions, summary, and code.
- **Merge:** no. Marker enrichment, marker evidence, curated annotation, clustering, and sample-level contrasts have different inputs, uncertainties, and failure semantics. Coupling them would erase the Provisional/Curated boundary.
- **Delete:** no. The deletion test shows the adapter and validation/provenance policy would otherwise be duplicated in workflows.

## Migration plan

1. Preserve the node identifier and existing selected output column defaults (`celltypist_cell_type`, `celltypist_confidence`) so downstream annotation-key workflows continue to resolve.
2. Add `summary` and `code` outputs; retain `adata` as the primary output.
3. Introduce explicit expression source/state. Interpret persisted `use_raw` only as a deprecated source alias; do not infer state from it. Saved workflows lacking state should receive a precise migration error or an explicitly selected legacy compatibility mode with a prominent warning—never a silent double transform.
4. Change new-node majority voting default to the official `False`. Persisted workflows that explicitly stored `True` keep that intent, but must supply valid grouping or deliberately opt into disclosed automatic over-clustering and satisfy the cell-count precondition.
5. Correct confidence semantics even if numeric values change: existing majority-vote workflows currently expose a misleading pairing. Record the old/new meaning in metadata rather than preserve the defect.
6. Add canonical individual/majority evidence and probability matrix without changing the user-configured selected label key. Downstream marker review/curation can migrate incrementally.
7. Version the `.uns` and summary schemas. Generated code and runtime share the same resolved fixed policy and postconditions.

## Test obligations before implementation is accepted

- Raw-count and pre-normalized inputs each undergo exactly the intended transform; ambiguous and double-transformed cases fail.
- Observation/feature duplicates, nonfinite/negative values, zero count rows, missing source, insufficient feature overlap, invalid model, and every collision fail before the adapter is called.
- Majority off/on, explicit grouping, `Heterogeneous`, `<=50` cells, invalid grouping/minimum proportion, and automatic-mode disclosure are covered.
- A fake result proves selected labels align with explicitly looked-up selected confidence, including disagreement between individual and majority labels.
- Probability/decision rows and ordered classes are checked; malformed index/shape/range/classes fail without partial writes.
- Model name/path resolution, SHA-256, metadata, query-feature overlap, and installed versions appear in `.uns` and strict JSON.
- Input AnnData is unchanged; outputs are categorical/numeric with complete count maps and no NaN/Infinity in JSON.
- Generated code produces equivalent scientific `.obs`, `.obsm`, `.uns`, and summary and rejects the same malformed fake backend results.
- Reports always use provisional language and include closed-vocabulary, score-calibration, majority-smoothing, and batch/domain-shift limitations.

## Final implementation resolution

The accepted implementation deepens two internal seams without expanding the node's public scientific inputs. The annotation core returns canonical diagnostics, while one pure, standalone-embeddable summary builder converts those complete diagnostics into the strict report dictionary. Both the node runtime and generated source call that same builder. The generated function interface is frozen as `run_celltypist_annotation(adata) -> (output_adata, summary_dict)`; timing and plugin history remain explicitly non-scientific wrapper metadata.

Software versions are collected once inside the standalone execution core, persisted in `uns[metadata_key]`, and passed through diagnostics to the pure summary builder. This single-source design prevents dependency provenance in AnnData from drifting from the `summary` payload.

The existing `adata, summary, code` node outputs and input order remain unchanged. The class additionally implements `fingerprint_inputs(model, **kwargs)` so cache identity follows the resolved local model artifact rather than only the model string. The fingerprint contains canonical path, size, `mtime_ns`, and SHA-256, resolves names only with official `get_model_path`, and never invokes model enumeration or download APIs.

Overwrite is treated as a provenance-owned replacement operation, not permission to delete arbitrary fields. A colliding prior artifact must expose schema version 1, operation `celltypist_annotation`, exact column/matrix ownership, and consistent stored matrix shapes. A changed selected-label key also requires the corresponding OpenBio annotation entry to equal the verified prior provenance before removal. All stale owned columns, matrices, and annotation entries are removed only on the private output copy after backend validation; ownership ambiguity fails before inference and leaves the input unchanged.

The expression adapter now enforces an absolute per-row CP10K tolerance of one count, retains zero-count entries for every model class in both individual and selected label maps, and adds `Heterogeneous` only when observed. Runtime/generated equivalence tests cover verified counts and float32 CP10K/log1p, majority voting, optional decision storage, model path/SHA/software versions/warnings, stale-artifact replacement, and identical malformed-backend failures.

Model provenance language deliberately avoids claiming that a resolved artifact is trusted. Missing-model guidance shows the explicit network-enabled official download command using the actual requested identifier, while execution remains offline. Summary warnings and load errors disclose that `.pkl` loading can execute code and that a SHA-256 fingerprint proves byte identity, not safety.

## Primary references informing this design

- [CellTypist official annotate API](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html)
- [CellTypist official annotate implementation](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/annotate.html)
- [CellTypist official classifier and AnnotationResult source](https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html)
- [CellTypist official model source](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/models.html)
- [CellTypist official tutorial](https://celltypist.readthedocs.io/en/stable/notebook/celltypist_tutorial.html)
- [CellTypist official repository](https://github.com/Teichlab/celltypist)
- Domínguez Conde C, Xu C, Jarvis LB, et al. *Science*. 2022;376:eabl5197. [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197); [PubMed](https://pubmed.ncbi.nlm.nih.gov/35549406/).
