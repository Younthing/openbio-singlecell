# scVI Integration: module design review

## Current module, outputs, and consumers

At the start of the original audit, the node performed four actions in one call:

1. resolves `X` or a layer and registers that location with `SCVI.setup_anndata`;
2. constructs and trains `scvi.model.SCVI`;
3. always writes the latent mean and optionally writes posterior parameters and MDE coordinates;
4. returns both mutated-copy AnnData and a live process-local `SCVIModel` wrapper.

The concrete `OPENBIO_SCVI_MODEL` output has one registered consumer: `OpenBioSingleCellSCVIDifferentialExpression`. The wrapper binds a trained backend to exact `var_names` and a set of training observations, rebuilds aligned temporary AnnData for downstream subsets, calls differential expression without retraining, and deregisters temporary managers. `tests/test_scvi_model.py` verifies this lifecycle, inference-mode escape, source routing, optional MDE/distribution outputs, exact-axis failures, and wrapper immutability of parameter dictionaries. Registration and example tests assert the model connection and `X_scVI` use by Neighbors.

The refactored node and model artifact now cover source routing, covariates, scoped RNG state, output collision, training diagnostics, strict reporting, generated code, and concurrent consumers. This follow-up corrects the remaining overvalidation of provenance, fractional values, zero-total genes, architecture size, train fraction, and constant latent dimensions.

## Correctness and cohesion problems

1. The selected source is never validated as raw counts. A log-normalized matrix can be trained silently even though official guidance requires raw counts.
2. `batch_key` defaults to empty, so the Integration node can perform no integration while retaining that name.
3. Batch/covariate controls do not distinguish Technical batch from Sample or Condition and do not protect biologically relevant columns.
4. Missing categorical values, non-finite/constant continuous covariates, invalid linear size factors, duplicate roles, zero-total cells/genes, duplicate axes, and backed inputs are not checked.
5. Current architecture defaults differ from official scVI defaults without rationale; `gene-label` is offered without registering labels, and `gene-cell` resource risk is unbounded.
6. `max_epochs=0` is an undocumented UI sentinel for omission. Early stopping is enabled without exposing/reporting validation split, monitor, patience, actual epochs, or stopping outcome.
7. Device selection is implicit `auto`; actual CPU/GPU device, resource controls, and nondeterminism limitations are absent.
8. Setting global `scvi.settings.seed` is not scoped/restored. Concurrent executions and later scvi operations can influence each other.
9. MDE and posterior-parameter booleans broaden the node beyond its core fitted-model/latent result and multiply storage-key collision modes.
10. The wrapper publicly exposes mutable registered AnnData, has no lock for concurrent consumers, and cannot persist across process restart.
11. The output has no `summary` or `code`, and generic history cannot support a scientific report.

## Decision: retain one fit-plus-canonical-latent node; remove auxiliary modes

Do not split training and canonical latent extraction in this batch. They share the same required registered AnnData, latent extraction is the direct fitted-model result, and most workflows need both the annotated AnnData and reusable model. A separate mandatory embedding node would add a fragile model/AnnData pairing and an extra interface without removing meaningful complexity. The current model artifact already earns its seam through the differential-expression consumer.

However, narrow the node to one atomic analysis: fit one scVI model and materialize its canonical latent mean. Remove `compute_mde`, `store_latent_distribution`, `qzm_key`, `qzv_key`, and `mde_key` from the core interface through an explicit saved-workflow migration. If posterior parameters or MDE remain desired, create separately researched model-consuming nodes; MDE is visualization, not training.

Do not merge upstream HVG selection, downstream Neighbors, clustering, integration-quality assessment, or differential expression. These have independent scientific parameters and claims. In particular, integration-quality metrics must compare Technical-batch removal with Condition/biology conservation and should not be hidden in training.

## Target interface

- Visible scientific inputs: `adata`, explicit count `source`, required `technical_batch_key`, optional additional technical categorical/continuous covariates, `n_latent`, and `gene_likelihood`.
- Advanced model inputs: `n_layers`, `dropout_rate`, and supported dispersion (`gene` or `gene-batch`; omit misleading `gene-label`, bound or omit `gene-cell`).
- Advanced training inputs: automatic/custom epoch mode, early-stopping mode plus explicit validation/monitor policy, minibatch size, and seed.
- Advanced compute/state inputs: accelerator (`auto`, `cpu`, `gpu`, `mps` as supported), single device selection, `output_key`, and `overwrite_existing`.
- Hidden policy: copy input, no filesystem writes/checkpoints, one device, public scvi-tools calls only, finite non-negative selected-expression/covariate/latent checks, scoped seed mutation, and ComfyUI inference-mode escape. Provenance and integer likeness are reported evidence, not claims that the input is verified raw UMI.
- Outputs: copied `adata`, concrete process-local `model`, strict JSON `summary`, equivalent Python `code`.

`size_factor_key` may remain advanced only with strict positive-linear validation. The simple/default path uses library size. The UI should represent automatic epochs as a named mode rather than magic integer zero.

## Deep model-artifact seam

Keep `OPENBIO_SCVI_MODEL`, but deepen `SCVIModel`:

- own a private registered training copy and do not expose it for mutation;
- freeze observation/variable identities, source and training metadata, runtime versions, actual device, and training diagnostics;
- serialize backend inference/manager registration with a lock;
- provide narrow methods needed by real consumers rather than exposing the backend;
- state process-local/nonpersistent lifetime explicitly;
- keep filesystem save/load outside this object and outside Integration.

The existing exact-variable-order and observation-subset checks are appropriate. Consumer-specific observation metadata may be joined onto an aligned private copy, but Condition comparisons remain subject to the repository's replicate-aware inference policy. The existence of a scVI DE method does not convert pooled cells into independent Samples.

## Validation and scientific claim contract

- Require in-memory AnnData, nonempty unique axes, and complete nonblank Technical-batch labels. Treat one-level and singleton strata as noninformative/unstable warnings rather than hidden sample-size gates.
- Validate the exact selected source as finite and non-negative with nonzero cell totals. Integer-like UMI counts are recommended; fractional values, all-zero genes, unknown provenance, and provenance indicating log/scaled/residual/other transformed expression run with strong machine-readable warnings because the public backend remains executable. Never report these states as verified raw UMI.
- The explicit source selector is the analyst's opt-in to train on that representation. Do not add a second Raw/provenance binding gate; report corrected/pseudocount or transformed evidence exactly as observed.
- Require finite continuous nuisance covariates and complete categorical nuisance covariates; disclose constant continuous covariates as noninformative. Require positive finite linear size factors and no duplicate/incompatible observation-column roles.
- Reject explicit overlap between Technical-batch/covariate keys and supplied Sample/Condition keys. Always warn about untestable confounding.
- Require positive integer `n_latent`, `n_layers`, batch size, and fixed epochs, plus the backend's exact `0 < train_size <= 1` domain and PyTorch dropout domain `0 <= dropout_rate <= 1`; remove arbitrary UI maxima, the 0.5/0.99 train-size narrowing, and the artificial `dropout_rate=0.99` ceiling. Warn for overcomplete latent width, unusually deep models, small training fractions, no validation holdout, and the executable-but-degenerate `dropout_rate=1` endpoint. Resolve `check_val_every_n_epoch=None` when `train_size=1.0`; reject only the impossible combination of full training allocation with validation-based early stopping. Reject output-key collisions unless overwrite is true.
- Verify backend trained state, exact latent shape, and finite values. Count/warn constant latent dimensions instead of rejecting the trained model. Capture actual epochs and finite history values.
- Never infer integration quality from training completion, loss, cluster appearance, or Technical-batch mixing alone.

## Summary contract

`summary` includes:

- methods: scVI generative model on the explicitly selected count-intended source, its resolved state/evidence and numeric audit, registered Technical-batch/nuisance design, architecture, training and latent-mean extraction; do not claim verified raw UMI without such evidence;
- results: trained-state and latent structural diagnostics, actual epochs, final/best loss values when available, without claiming successful batch removal;
- key results: input dimensions/totals, feature count, Technical-batch/covariate level sizes, train/validation cell counts, latent dimension/variance diagnostics, actual device, and stopping information;
- parameters: every source/model/training/device/output value plus resolved automatic values;
- warnings: process-local model, stochastic/cross-device reproducibility limits, possible confounding and over-correction, HVG-view scope, integration-quality assessment still required, and no Sample-level Condition inference;
- references: scVI method, scvi-tools software, official API/tutorial, integration benchmark, and AnnData semantics as applicable;
- software versions: openbio-singlecell, scvi-tools, PyTorch, Lightning, AnnData, NumPy, pandas, and accelerator libraries where available.

All values must be JSON-native and finite; training histories are summarized with bounded scalar/list previews rather than embedded DataFrames.

## Equivalent-code feasibility

Generated code is scientifically equivalent and feasible. It should define one self-contained function that:

1. copies and validates AnnData/non-negative count-intended source/covariate state while reproducing advisory provenance, fractional, and zero-gene warnings;
2. snapshots and temporarily sets scvi/NumPy/Torch seed state where public facilities allow;
3. runs `SCVI.setup_anndata`, constructor, and `train` with every resolved parameter explicit inside `torch.inference_mode(False)`;
4. validates/stores the canonical latent mean;
5. restores global settings in `finally` and returns `(output, model)`.

It does not reproduce the plugin socket wrapper or history bookkeeping and does not contain trained weights. The text must explain that executing it retrains the model. Runtime and code must agree on accelerator/device, epoch/early-stopping policy, source validation, key collisions, and latent state, including advisory warnings and constant-dimension accounting. If deterministic algorithms are not forced by the runtime, the code must not promise bitwise equality.

## Migration implications

Saved workflows should:

- migrate old `batch_key` to `technical_batch_key` with a visible warning that semantics must be reviewed;
- preserve source, architecture, likelihood, epochs, early stopping, output key, and seed;
- remove optional MDE/posterior widgets and record that these outputs are no longer generated by Integration;
- map `max_epochs=0` to explicit automatic mode;
- require the user to resolve any removed `gene-label`/`gene-cell` choice rather than silently replacing it;
- preserve `adata` and model links and append new summary/code outputs without shifting the primary adata/model slots if canvas compatibility requires stable indices.

## Verification plan

- Dense/CSR/CSC sources from `X` and named layers; reject negative, non-finite, zero-total-cell, backed, empty, and duplicate-axis inputs. Verify fractional, provenance-marked transformed, and zero-total-gene inputs run with exact warnings and are not called verified raw UMI.
- Technical-batch/covariate tests for missing/blank/null values, warning-only single-level/tiny/constant choices, continuous non-finite values, key overlap, and protected Sample/Condition conflicts.
- Fake scvi backend asserts every setup/constructor/train/latent argument, automatic/custom epochs, early-stopping policy, CPU/GPU selection, scoped seed restoration, and no filesystem writes.
- Preserve the existing inference-mode gradient test and exact model/AnnData identity tests; add lock/concurrent-consumer and private-state mutation tests.
- Validate actual training-history extraction, finite latent shape/variance, output collision handling, strict summary/citations/versions, and cautious wording.
- Generated code compiles and matches runtime primary AnnData state and fake-backend calls; plugin history, elapsed timing, and wrapper type are the only allowed differences.
- Optional real-scvi smoke tests should be version-gated, tiny, device-explicit, and skipped with an actionable reason when the optional dependency is absent.
- Cover `n_latent` wider than cells/features, `n_layers>20`, `train_size<0.5` and `train_size=1`, and a fake constant latent dimension as warning-only cases; keep only positive scalar/backend-domain failures hard. The audited 1.5.0.post1 real smoke demonstrates zero-total-gene, overcomplete latent, and full-training-set execution.
