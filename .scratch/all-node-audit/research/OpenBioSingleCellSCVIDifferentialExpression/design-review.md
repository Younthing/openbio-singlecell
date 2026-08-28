# OpenBioSingleCellSCVIDifferentialExpression: design review before refactor

## Decision

Retain the node ID, deepen the module, and change its user-facing meaning to **scVI Model DE Evidence**. Do not merge it into `OpenBioSingleCellSCVIIntegration`, do not merge it with Sample-level pseudobulk contrasts, and do not delete it.

The retained atomic responsibility is:

> Given one compatible observation subset and one already-fitted `SCVIModel`, compare exactly two cell populations once, validate the released scvi-tools result, and report explicitly exploratory model evidence.

The node must never train or load a model, perform clustering, run multiple contrasts, choose markers for annotation, or claim Sample-level Condition inference.

## Repository evidence reviewed

The review covered:

- `openbio_singlecell/nodes_differential.py`, including neighboring pseudobulk nodes;
- `openbio_singlecell/scvi_model.py`;
- the refactored `OpenBioSingleCellSCVIIntegration` implementation and its `SCVIModel` output contract;
- `tests/test_scvi_model.py`, `tests/test_integration_nodes.py`, registration tests, and workflow/example assertions;
- `.scratch/all-node-audit/research/OpenBioSingleCellSCVIIntegration/{official-usage.md,design-review.md}`;
- `CONTEXT.md`; and
- `docs/adr/0001-separate-marker-evidence-from-condition-inference.md`.

The domain boundary is already clear: a Sample is the independent biological specimen, a Technical batch is a nuisance source, a Condition is the biological cohort/treatment, and a formal Condition contrast is replicate-aware inference within a defined population. This node does not implement that final contract.

## Current interface and behavior

The current input order is:

1. `adata`
2. `model`
3. `groupby="group"`
4. `group1=""`
5. `group2=""`
6. `subset_column=""`
7. `subset_value=""`
8. `mode="vanilla"`
9. `delta=0.25`

It returns only `table`.

The node delegates most state handling to `SCVIModel.differential_expression`. That wrapper already supplies a valuable seam: it rejects detached/untrained models, requires exact feature order and known observation identities, builds a private registered-count AnnData in downstream observation order, copies only required observation metadata, serializes temporary-manager access, preserves a primary backend error, and does not retrain.

The node then calls `.reset_index()`, conditionally renames `index` to `gene`, wraps whatever columns the backend returned, and reports only a small parameter dictionary. Existing fake-backend tests establish lifecycle and alignment but intentionally return only `lfc_mean` and `proba_de`, so they do not test the real scvi-tools output contract.

## Correctness findings

1. **Scientific status is unstated.** The current title and description imply generic differential-expression inference. The method is conditional cell/population model evidence and is not replicate-aware Condition inference.
2. **`vanilla` is the wrong new default.** The official guide recommends `change` to avoid prioritizing negligible effects. Released `vanilla` is also directional (`scale1 > scale2`), which is easy to misread as a two-sided point-null test.
3. **`delta` is always visible but irrelevant in `vanilla`.** The interface permits meaningless parameter combinations.
4. **Batch behavior is implicit.** The current call inherits `batch_correction=False`, so cells are decoded in their observed primary batches. Group-specific batch composition can drive the comparison.
5. **No overlap policy exists.** Simply enabling all-batch correction would allow unsupported group-by-batch counterfactuals. Complete Condition/Sample confounding must fail, not be hidden by decoding.
6. **Other nuisance covariates are undisclosed.** `batchid1`/`batchid2` standardize only the primary registered batch, not every categorical/continuous covariate supplied during training.
7. **The hypothesis is under-specified.** `test_mode`, pseudocount policy, posterior sample count, weighting, outlier filtering, and posterior-FDR target are implicit lower-level defaults.
8. **Monte Carlo state is unmanaged.** Posterior expression sampling and NumPy sampling make the result stochastic. The node has no DE seed, no scoped RNG restoration, and no explicit serialized RNG section.
9. **Input expression semantics are opaque.** The downstream `adata` supplies identities and grouping metadata, while expression and raw-stat fields come from the model's private registered count view. Users can otherwise assume current `adata.X` was analyzed.
10. **Model evidence is incomplete.** Only model class and training parameters reach table metadata. Training diagnostics, fitted-weight identity, registered-count identity, primary Technical batch, device, and process-local limitation are absent.
11. **Backend output is trusted.** Missing genes, duplicate features, a colliding `gene` column, wrong probability ranges, non-finite values, inconsistent complements, unexpected FDR columns, wrong group labels, and malformed booleans pass through.
12. **No stable table schema exists.** The FDR column contains a floating-point target in its name, and real `vanilla`/`change` schemas differ substantially.
13. **No result narrative exists.** Counts, batch support, tagged-gene count, signed leading effects, empty findings, and interpretation for writing are absent.
14. **No references, versions, summary, or code are emitted.** This violates the all-node audit output contract.
15. **The live model cannot be reconstructed from code alone.** Generated code that silently retrained a replacement model would not be scientifically equivalent.

## Module boundary and depth

The trained-model artifact is the correct seam. The public node should be thin, while the `SCVIModel` module remains deep:

- The **node interface** expresses one scientific comparison and report policy.
- The **model adapter** owns fitted-model compatibility, registered data, manager lifecycle, locking, exact backend invocation, and immutable provenance.
- A narrow **result adapter** owns mode-specific schema validation and canonicalization.
- The shared reporting layer owns strict JSON, references, software versions, and code text.

This division keeps backend complexity local. Callers should not need to understand scvi-tools AnnData managers, temporary registration, global RNG state, dynamic FDR column names, or posterior batch pairing.

The model artifact should be deepened before or with this refactor to expose immutable evidence rather than forcing the node to reach into backend internals. Its public evidence should include an artifact schema version, model class, fitted-state fingerprint, registered-data/axis fingerprint, training parameters, training diagnostics, ordered observation/feature identities, and process-local status. Fingerprints are OpenBio adapter provenance, not scvi-tools fields.

## Why not merge or delete

### Do not merge with scVI Integration

Training and model-DE have different rates of change and different scientific responsibilities. One model is expensive to fit and should support several independently configured evidence contrasts without retraining. Merging would broaden the interface, entangle training failures with contrast failures, and make each additional contrast repeat stateful work.

### Do not merge with pseudobulk Condition contrast

The products have different statistical units and meanings:

- scVI model evidence aggregates cell-level approximate posteriors and decoded expression;
- pseudobulk Condition inference aggregates full-gene counts within Sample and estimates between-Sample variation.

A single “DE” node with a method switch would hide this epistemic boundary and invite invalid substitutions.

### Do not delete

Deletion would remove a useful validated adapter for exploratory cluster/state comparisons and force users to manipulate a live scvi-tools model and manager lifecycle manually. The node earns its place if it adds strict alignment, supported batch conditioning, stable results, and cautious reporting. Without those protections, deletion would be preferable to the current ambiguous wrapper; with them, retention is justified.

## Proposed atomic interface

Keep the class/node ID `OpenBioSingleCellSCVIDifferentialExpression` for serialized-workflow migration, but change the display name to `scVI Model DE Evidence` and the operation name to `scvi_model_de_evidence`.

### Inputs in exact order

1. `adata: OPENBIO_ANNDATA`
2. `model: OPENBIO_SCVI_MODEL`
3. `groupby: STRING = "group"`
4. `group1: STRING = ""`
5. `group2: STRING = ""`
6. `population_scope: DYNAMIC_COMBO = all`
   - `all`: no nested inputs;
   - `obs_value`: `subset_column: STRING = ""`, `subset_value: STRING = ""`.
7. `mode: DYNAMIC_COMBO = change`
   - `change`: `delta: FLOAT = 0.25` and advanced `fdr_target: FLOAT = 0.05`;
   - `vanilla`: no nested inputs and no FDR claim.
8. `batch_handling: COMBO = shared_technical_batches`
   - `shared_technical_batches`;
   - `observed_technical_batches`.
9. `n_samples_overall: INT = 5000` (advanced)
10. `random_seed: INT = 0` (advanced)

### Outputs in exact order

1. `table: OPENBIO_TABLE`
2. `summary: OPENBIO_SUMMARY`
3. `code: STRING`

### Exposed versus fixed policy

Expose parameters that alter the scientific question or Monte Carlo precision:

- exact population labels and optional restriction;
- practical-change versus legacy directional mode;
- `delta` and posterior FDR target for `change`;
- observed versus supported shared-batch conditioning;
- posterior sample count; and
- seed.

Fix and disclose the following:

- one contrast only; explicit `group2`; no one-versus-rest;
- `group1=[group1]` at the official backend boundary;
- `all_stats=True`;
- `silent=True`;
- `weights="uniform"`;
- `filter_outlier_cells=False`;
- `dataloader=None` and backend default batch size;
- `pseudocounts=None`, allowing released data-driven offset estimation;
- `change_fn=None`, `m1_domain_fn=None`, and no credible intervals;
- `test_mode="two"` for `change`, so `proba_de` represents the practical two-sided event rather than the released default maximum of signed alternatives;
- `use_permutation=False`; and
- no arbitrary `**kwargs`, custom query strings, importance weighting, custom batch-ID lists, or custom callable hypotheses.

`test_mode="two"` is an intentional OpenBio adapter decision and must be stated in the summary and code. It aligns the probability, Bayes factor, and posterior-FDR tag with the reported hypothesis `abs(LFC) >= delta`. This decision must have a direct fake-backend assertion because it differs from the lower-level released default.

## Batch policy

### `shared_technical_batches` (new default)

1. Read the primary `technical_batch_key` from immutable model training metadata.
2. Within the selected population scope, compute observed Technical-batch levels separately for `group1` and `group2` without string-coercion collisions.
3. Preserve registered categorical order and take the ordered intersection.
4. Require at least one shared level. Complete confounding is an actionable failure directing Condition questions to Sample-level inference.
5. Pass the same list as `batchid1` and `batchid2` with `batch_correction=True`.
6. Do not include a batch level unless each population actually has cells in it. This avoids unsupported group-by-batch combinations at the population level.
7. Report all registered levels, observed levels by group, shared levels used, excluded non-shared levels, and the group-by-batch cell-count table.

This is counterfactual standardization on supported primary Technical batches. It is not evidence that all nuisance variation was removed, and it does not standardize additional categorical/continuous covariates.

### `observed_technical_batches`

Pass `batch_correction=False` and omit both batch-ID arguments. Report the actual group-by-batch composition and warn that it contributes to the estimand. This option preserves the legacy backend policy for workflow migration and for carefully justified descriptive comparisons.

Do not expose the released disjoint-batch behavior or “all registered batches” as ordinary UI options. Both can encourage extrapolation or compare different decoder targets. A future expert node may add them only with a separate research record and overlap diagnostics.

## Validation and boundary conditions

### Model and axes

- Require exactly the supported `SCVIModel` artifact schema and a currently trained, attached backend.
- Recompute/verify available artifact fingerprints; never trust copied metadata alone.
- Require nonempty, unique ordered `obs_names` and `var_names`; identifiers must be strings, nonempty, finite where applicable, and free of surrounding whitespace.
- Require every downstream observation to be in the fitted-model observation set and exact equality of ordered downstream `var_names` with model `var_names`.
- State explicitly that downstream expression matrices are ignored; the model-owned registered counts supply decoding inputs and `all_stats` fields.
- Require a fitted-state/weight fingerprint for a fully hardened contract. Until the artifact provides one, disclose that exact fitted weights cannot be independently identified from the result.

### Labels and scope

- Require canonical nonempty column names and label strings; do not strip into a different identifier or coerce numbers, booleans, objects, or missing values to strings.
- Require `groupby` to exist and not equal any declared Technical/nuisance covariate used during training.
- Require distinct, observed `group1` and `group2` values after optional restriction.
- Require at least two cells per population; disclose and warn for small populations or small group-by-batch cells rather than presenting them as stable.
- Preserve categorical order and report unused declared categories, missing-label cells, other-label cells, and cells excluded by the optional restriction.
- `obs_value` requires both nested strings, exact type-compatible matching, and a nonempty result. `all` rejects stray legacy subset values during migration.
- Detect reserved/colliding columns used by scvi-tools temporary DE and raw-normalization caching; fail or isolate them on a private copy without overwriting user data.

### Numeric parameters

- `mode` must resolve to exactly `change` or `vanilla`.
- `delta` must be a finite positive float for `change`; it must be absent for `vanilla`.
- `fdr_target` must be finite and strictly between zero and one for `change`; it must be absent for `vanilla`.
- `n_samples_overall` must be a non-boolean positive integer. For shared-batch decoding it must yield at least two posterior samples per selected batch; warn when below a documented stability threshold.
- `random_seed` must be a non-boolean integer in `[0, 2**31 - 1]`.

### Model-training evidence

- Require the producer's declared count source/state and primary Technical batch key.
- Reject a `groupby` that the model was told to treat as Technical batch or nuisance; asking the decoder to remove a variable and then treating it as biology is incoherent.
- Carry the integration training parameters and diagnostics into summary by immutable copy.
- Report feature scope as the model's fitted feature set/HVG view, never as all genes unless verified.
- Training completion, a finite latent representation, and finite loss are necessary execution evidence, not scientific validation.

## Execution pipeline

1. Resolve the dynamic inputs without accepting orphaned legacy values.
2. Validate model artifact, axes, identifiers, grouping columns, and parameters before backend work.
3. Build the population mask and exact group masks; compute all exclusion and batch-support counts.
4. Resolve the strict batch policy. Fail before inference if shared support is absent.
5. Acquire the model's manager lock and a process-wide scoped scVI-DE RNG lock.
6. Snapshot Python, NumPy, Torch CPU, available CUDA, available MPS, and scvi global seed state; set the requested seed; execute outside an inherited inference-mode restriction where required by the backend.
7. Construct a private registered-count AnnData in selected observation order and copy only validated observation metadata.
8. Call the public `differential_expression` method once with every fixed and exposed keyword explicit.
9. In `finally`, retrieve and publish the exact transferred analysis manager through the released public
   `get_anndata_manager`/`register_manager` seam, then call `deregister_manager(analysis_adata)`; never call the
   broad no-argument cleanup. Preserve the primary exception and restore every RNG/global state in `finally`.
10. Validate and canonicalize the backend DataFrame.
11. Build the table wrapper, strict JSON summary, and equivalent code string only after all postconditions succeed.

The operation is copy-atomic: neither input AnnData, model-owned registered AnnData, downstream categorical metadata, global RNG settings, nor backend manager stores may be observably altered after success or failure.

## Table contract and postconditions

The canonical leading columns are `gene`, `comparison`, `group1`, and `group2`. Keep released scientific field names after those columns, except rename the one exact dynamic `is_de_fdr_<target>` column to stable `is_de_fdr`; store its numeric target in table parameters and summary.

Require:

- a pandas DataFrame with one row per fitted model feature;
- exact feature-set equality with the model's ordered features, no duplicates, no missing/unknown features, and one unambiguous feature index;
- exact group labels and comparison text in every row;
- exact mode-specific required columns and no colliding `gene` or canonical FDR column;
- numeric columns with numeric dtypes, finite values, and no booleans masquerading as numbers;
- probabilities and nonzero proportions in `[0, 1]`;
- complement pairs summing to one within a declared floating tolerance;
- nonnegative scales/raw means, positive realized pseudocount, constant expected `delta`, and a Boolean FDR tag;
- stable sorting: descending `proba_de` for `change` or `bayes_factor` for `vanilla`, with fitted feature order as a deterministic tie-breaker; and
- JSON-safe metadata with no NaN or infinity.

Mode-specific required fields are the released fields documented in `official-usage.md`. Unexpected additional numeric credible-interval fields should be rejected because this interface never requests them; silently passing version-dependent columns would defeat the stable contract.

For `change`, positive `lfc_mean` is reported as higher decoded expression in `group1`; negative is higher in `group2`. For `vanilla`, do not invent an FDR tag or two-sided effect claim.

## `summary` JSON contract

The `summary` endpoint is strict JSON with top-level schema identifier `openbio-singlecell/scvi-model-de-evidence/v1` and status `exploratory_model_evidence`. It must contain at least:

```json
{
  "schema_version": "openbio-singlecell/scvi-model-de-evidence/v1",
  "status": "exploratory_model_evidence",
  "node_id": "OpenBioSingleCellSCVIDifferentialExpression",
  "method": {},
  "comparison": {},
  "model_evidence": {},
  "analysis_summary": {},
  "parameters": {},
  "references": [],
  "software_versions": {},
  "warnings": [],
  "limitations": []
}
```

Required content is:

- `method`: scVI model-DE method, exact mode/hypothesis, decoded-expression interpretation, `test_mode`, pseudocount policy and realized value, posterior sampling, weighting, outlier, raw-stat, sorting, and posterior-FDR semantics;
- `comparison`: exact groups, optional scope, group cell counts, missing/other/excluded cells, feature count/order fingerprint, Technical-batch policy, registered/observed/shared/excluded batch levels, and group-by-batch counts;
- `model_evidence`: model/artifact schema, model class, fitted-state fingerprint if available, registered-data/axis fingerprints, process-local status, count source/state, fitted feature scope, Technical/nuisance covariates, training parameters, training diagnostics, seed, and device;
- `analysis_summary`: tested feature count, posterior-FDR-tagged count or `null` for `vanilla`, tagged counts split by positive/negative `lfc_mean`, and up to ten leading tagged positive and negative effects with gene, `lfc_mean`, `proba_de`, Bayes factor, tag, and raw-count summaries. If no feature is tagged, both lists remain empty and a separate, clearly non-selected `top_model_evidence` list may show the highest posterior probabilities without calling them discoveries;
- `references`: structured citations for Lopez et al., Boyeau et al., Gayoso et al., scvi-tools API/source, and the replicate-inference limitation papers;
- `software_versions`: openbio-singlecell, scvi-tools, PyTorch, Lightning, AnnData, NumPy, pandas, Python, and relevant accelerator runtime;
- `warnings`: data-specific instability, observed-batch composition, small support, missing fitted fingerprint, or runtime-version qualification; and
- `limitations`: always state no Sample-level Condition inference, no p-value/BH interpretation, model/HVG scope, conditionality on fitted weights and nuisance choices, possible cluster/model double use, process-local model, and lack of integration-quality proof.

The writing-oriented core must use cautious phrases such as “posterior-FDR-tagged model evidence.” It must never call genes “significant” or describe a Condition effect as established.

## `code` contract

`code` is a UTF-8 source string defining an equivalent function such as:

```python
def run_scvi_model_de_evidence(adata, trained_model):
    ...
    return table, evidence
```

The source must:

- accept the same already-fitted official model (or a clearly documented compatible artifact) as an argument;
- contain the resolved node parameters as literals or explicit function defaults;
- perform the same axis/label/scope/batch-support validation;
- construct the same private analysis AnnData without mutating inputs;
- scope and restore the same RNG/global states;
- make one public `trained_model.differential_expression` call with all arguments explicit;
- preserve cleanup/primary errors;
- apply the same strict output validation, canonicalization, and deterministic order; and
- return the canonical DataFrame plus the evidence dictionary used to form the report.

It must not train a replacement model, download data, save/load hidden checkpoints, use private scvi-tools APIs, omit the batch policy, or claim it can reproduce fitted weights. The header and summary must say that numerical equivalence is conditional on the same fitted model weights, registered count view, software/device, and input observation metadata. This limitation is part of scientific equivalence, not a reason to generate misleading self-contained retraining code.

## Test plan for implementation

### Interface and migration

- Assert the exact new input order, dynamic options/defaults, output order/types, display name, and operation/status strings.
- Assert legacy workflows explicitly migrate `vanilla` plus `observed_technical_batches`; new instances receive `change` plus `shared_technical_batches`.
- Reject orphaned `delta`, FDR, or subset values created by malformed migrations.

### Fake backend and call fidelity

- Assert one call, no retraining, scalar-to-one-element-list handling, and every explicit public keyword.
- Assert `test_mode="two"`, uniform weights, no outlier filter, all stats, FDR, posterior sample count, silence, and mode-specific omission of irrelevant parameters.
- Assert observed-batch mode supplies no batch IDs; shared mode supplies identical ordered IDs.
- Make the fake mutate the temporary AnnData and manager store; prove all user/model-owned inputs remain unchanged and cleanup occurs once.

### Batch and population support

- Cover balanced, partially shared, no-shared, unused categorical, missing group, missing batch, other-label, and restricted-population cases.
- Reject complete group/batch confounding under the default rather than decoding unsupported combinations.
- Prove categorical order is retained and per-group support is disclosed.
- Prove `groupby` cannot reuse a training nuisance role.

### Malformed results

- Missing/duplicate/unknown genes; wrong feature count; colliding `gene`; missing/extra required columns.
- String/object numeric columns; boolean numerics; NaN/infinity; negative scales/raw means; probabilities/proportions out of range; complement mismatch.
- Wrong group/comparison values; wrong or multiple dynamic FDR columns; non-Boolean tag; wrong realized delta/pseudocount.
- Backend row reordering and ties, verifying canonical stable order.
- Separate complete real-schema fixtures for `change` and `vanilla`; the current two-column fake is insufficient.

### Randomness, failures, and code

- Same seed/same fake model gives the same result; different seed reaches the backend as a distinct scoped state.
- Python/NumPy/Torch/CUDA/MPS/scvi states are restored on success and failure.
- Concurrent consumers remain serialized; no temporary manager leaks.
- A primary DE error wins over cleanup failure; cleanup-only failure is raised.
- Generated code compiles, contains no training/download/filesystem action, calls the same public API, and produces table/evidence parity with runtime for both batch modes.

### Optional real-scvi smoke test

Use a tiny, version-gated, CPU-only raw-count fixture with explicit skip reason when scvi-tools is absent. Verify the released mode-specific fields and manager cleanup, but keep fake-backend tests authoritative for deterministic edge cases.

## Migration risks

1. Output arity changes from one to three. Registration, generator metadata, migrations, examples, and downstream tests must update together after implementation is approved.
2. The new default changes from `vanilla`/observed batch to `change`/shared supported batches. Existing workflows must preserve their old semantics explicitly; silent scientific reinterpretation is unacceptable.
3. `subset_column`/`subset_value`, `mode`, and `delta` move into dynamic inputs. Positional widget migrations require exact mapping.
4. Stable `is_de_fdr` and stricter mode-specific columns may break consumers that depended on raw backend column names or tolerated partial fake results.
5. The `SCVIModel` public adapter must grow batch/FDR/sampling/test-mode support and immutable model provenance without exposing the raw backend.
6. Adding fitted-state/data fingerprints may make older in-memory artifacts incompatible. Since the artifact is already process-local, an actionable rerun of SCVI Integration is safer than silently accepting unverifiable state.
7. scvi-tools is optional and unpinned. Compatibility guards and version-qualified tests are required; a backend signature/output change must fail explicitly.
8. Saved workflow code cannot embed a live fitted model. The `code` contract requires the same model argument, and this limitation must remain visible.
9. The conservative shared-batch default will fail comparisons previously allowed under complete confounding. That failure is intentional and directs formal Condition questions to the correct Sample-level path.

## Acceptance criteria

The refactor is acceptable only when:

- one execution represents exactly one two-population model-evidence contrast;
- the official call and batch semantics are explicit and tested;
- no unsupported counterfactual or malformed backend result is silently accepted;
- input/model/RNG/manager state is atomic on success and failure;
- the table is stable and scientifically labeled;
- `summary` contains references, writing-ready key results, versions, model/training evidence, and mandatory limitations;
- `code` is equivalent conditional on the same fitted model and says so;
- no output or text can reasonably be read as replicate-aware Sample-level Condition inference; and
- scVI Integration and pseudobulk Condition inference remain separate modules with narrow, high-leverage seams.
