# OpenBioSingleCellAugur — design review

## Decision

**Enhance and retain as an atomic exploratory population-prioritization module with a typed Augur result artifact.**
It must not return a backend work AnnData. Augur Results remains a pure artifact-to-table adapter. This creates a real
seam: the analysis owns computation/integrity, while views do not know Pertpy's mutable `uns` layout.

## Atomic interface

Inputs, in order:

1. `adata` with an explicit expression source; Raw snapshot counts are recommended but explicitly selected, while
   the UI defaults to the current `X` representation; feature
   completeness remains user-declared rather than programmatically proven, and mutable analysis history is not an
   execution gate. Expert-selected `X`/layer values are accepted when numeric, finite, nonnegative, aligned, and
   nonempty; non-integer state is disclosed instead of silently truncated;
2. `sample_key`, default `sample`;
3. `population_key`, default `cell_type`;
4. `condition_key`, default `condition`;
5. required `control` and `treatment`;
6. classifier (`random_forest_classifier` or `logistic_regression_classifier` only);
7. explicit expression `source`, default `X`; selecting Raw is sufficient and adds no history/binding gate;
8. `annotation_status` (`unknown`, `provisional`, `curated`);
9. optional advanced `technical_batch_key` for audit only;
10. advanced `n_subsamples=50`, `subsample_size=20`, `folds=3`, `n_threads=1`;
11. nonzero `random_seed`;
12. result-row/memory guards.

Outputs:

```text
result, summary, code
```

Every downstream Augur Results view reports dimensions from `selected_source_features`; Raw-specific diagnostics are
nullable provenance fields and must never be used as a required view dimension when the expert selected X or a layer.

Fix original Augur variance feature selection (`feature_perc=0.5`, `var_quantile=0.5`, `span=0.75`, negative-residual
filter off), default mode, zero-division behavior, and estimator hyperparameters to Pertpy 1.3 defaults. They are
reported but hidden to keep the interface deep. Remove regressor, normalized-expression default, arbitrary result
key, and mutable AnnData output.

## Typed result artifact

The immutable process-local artifact contains canonical DataFrames:

- `priorities`: one row/population, with mean Augur/AUC/accuracy/precision/F1/recall and deterministic rank;
- `cross_validation`: one row/population/subsample/fold with AUC;
- `feature_importance`: population/subsample/fold/gene/importance (possibly large and guarded);
- validated parameters, selected-expression/observation fingerprints, advisory Sample/Condition/Technical batch
  audit, skipped
  populations and reasons, dynamic versions, strict scientific summary, artifact-content fingerprints, and producer
  schema.

Transpose/canonicalize Pertpy `summary_metrics`, validate exact eligible population families and finite metric ranges,
and reject the “no population worked” warning as a hard error. Do not mutate input AnnData. `summary` is the same
payload stored in the artifact. `code` defines an equivalent function returning portable DataFrames plus
`summary_dict`; it may omit the wrapper class but not validation or report content.

## Scientific boundary

Sample is audited but Pertpy 1.3 cross-validation remains cell-level; the report must not call it Sample-level
Condition inference. Low biological-Sample support, repeated-measure Sample mappings, and perfect or partial
Condition/Technical-batch confounding are warnings and limitations rather than hard errors because they affect
interpretation, not whether Pertpy's cell-level calculation is defined.
Reusing another role as `sample_key` similarly remains possible for an expert, but the summary must mark Sample
identity unverified and must not upgrade those labels to proven biological replicates.
Population labels and Provisional/Curated status are explicit. Augur AUC measures separability/responsiveness priority,
not effect direction or statistical significance.

## Migration and tests

Migration changes the primary type from AnnData to AugurResult, renames `cell_type_key` to `population_key`, defaults
to Raw counts, removes regressor/result key/span/select-variance surface controls, and adds Sample/audit/summary/code.
Old downstream Augur Results links can be rewired by type, but arbitrary AnnData consumers cannot be guessed.

Tests: Pertpy 1.3.x family plus exact audited public signatures; explicit two-Condition subset/categorical conversion; selected-source numeric/
axis validation and non-integer disclosure; classifier-only choices; advisory Sample mapping/replicate support and
Technical batch confounding; eligible/skipped populations;
metrics, transpose, rank/ties and range checks; feature-importance schemas for RF/logistic; Pertpy RNG-scope
disclosure; deterministic rerun; malformed backend/no-success; guards; input immutability; artifact tampering; strict
JSON; generated code compilation; and runtime/generated canonical table+summary equality.

## Cohesion assessment

The analysis module owns everything necessary to obtain and validate one Augur prioritization. Pertpy storage,
feature selection, cross-validation, audits, canonicalization, and reporting are hidden. Table selection is delegated
to a narrow adapter rather than coupled to AnnData state.
