# scVI Integration: official usage research

Researched: 2026-08-28

## Scope and primary sources

This record covers `OpenBioSingleCellSCVIIntegration`: registration of an explicitly selected count-intended representation, fitting one scVI model, and materializing its canonical latent representation for graph-based exploratory analysis. It also audits the trained-model output because the current downstream `OpenBioSingleCellSCVIDifferentialExpression` consumer depends on it. The selection is an analyst action, not proof that the values are raw UMI counts.

Primary sources:

- scvi-tools introduction/tutorial: https://docs.scvi-tools.org/en/latest/tutorials/notebooks/quick_start/api_overview.html
- scvi-tools `SCVI` API: https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html
- scvi-tools data-loading guidance: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/data_loading.html
- scvi-tools FAQ on batch/covariates and count input: https://docs.scvi-tools.org/en/stable/faq.html
- scvi-tools training configuration and callbacks: https://docs.scvi-tools.org/en/stable/user_guide/use_case/training_configuration.html and https://docs.scvi-tools.org/en/latest/user_guide/use_case/using_callbacks.html
- scvi-tools model save/load guide: https://docs.scvi-tools.org/en/latest/user_guide/use_case/saving_and_loading_models.html
- Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. *Nature Methods*. 2018;15:1053-1058. https://doi.org/10.1038/s41592-018-0229-2
- Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*. 2022;40:163-166. https://doi.org/10.1038/s41587-021-01206-w
- Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods*. 2022;19:41-50. https://doi.org/10.1038/s41592-021-01336-8

scvi-tools is optional and unpinned in this repository. The node must runtime-check the installed public methods/keywords and report the installed scvi-tools, PyTorch, Lightning, AnnData, NumPy, pandas, and openbio-singlecell versions.

## Canonical official workflow

The official tutorial uses a raw-count layer, registers study covariates, constructs the model, trains it, and stores its latent mean:

```python
import scvi

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="technical_batch",
    categorical_covariate_keys=["technical_protocol"],
    continuous_covariate_keys=["percent_mito"],
)
model = scvi.model.SCVI(
    adata,
    n_latent=10,
    n_layers=1,
    dropout_rate=0.1,
    dispersion="gene",
    gene_likelihood="zinb",
)
model.train(
    max_epochs=None,
    accelerator="auto",
    devices="auto",
    early_stopping=False,
)
adata.obsm["X_scVI"] = model.get_latent_representation()
```

The scvi-tools introduction recommends roughly 1,000-10,000 HVGs for scVI, depending on context. The repository ADR fixes the production workflow to a 5,000-HVG representation-learning view and explicitly separates it from full-gene, replicate-aware Condition inference.

## Count and feature requirements

Official guidance says scVI expects raw counts, not library-size-normalized or log-transformed values. A named `counts` layer is the standard way to preserve them. The tutorial notes that models can numerically run on non-negative real values, but those values should only be intentional corrected/pseudocounts; normalized data with altered variance/covariance is not an acceptable substitute.

For this repository's standard UMI path, the selected source should therefore be finite, non-negative, integer-like raw counts. That is the recommended interpretation, not a property that provenance metadata can prove. The public scvi-tools backend can train on explicitly selected finite non-negative fractional/corrected values and emits an advisory count-data warning; a real 1.5.0.post1 smoke completed on such input. The open expert node therefore warns prominently for fractional values or provenance indicating normalized/logged/scaled/residual expression, records the exact state evidence, and never relabels the source as verified raw UMI. It still rejects negative/non-finite values and zero-total cells, which make the count likelihood/library-size path invalid. Empty/duplicate axes remain alignment failures. All-zero genes are disclosed and warned rather than rejected: a real 1.5.0.post1 one-epoch smoke retained one such gene and produced a finite latent representation.

Feature selection is upstream. The node does not calculate HVGs or silently subset genes. It should disclose the training feature count and warn outside the usual 1,000-10,000 guidance. `n_latent` is an architecture choice, not a matrix-rank request: scvi-tools 1.5.0.post1 successfully trained with `n_latent > n_vars`. Require a positive integer, warn when it exceeds the cell or feature dimension, and let the backend own resource/training failures.

## Technical batch, Sample, Condition, and covariates

The official `setup_anndata` documentation is explicit: `batch_key`, categorical covariates, and continuous covariates are nuisance factors whose effects the model tries to minimize in latent space; biologically relevant factors that should be retained must not be passed there. The FAQ recommends `batch_key` for the main technical effect (lab, lane, chemistry, dataset origin) and categorical covariates for additional nuisance effects.

In this repository:

- **Technical batch** is the unwanted technical source supplied to `batch_key` or additional technical covariates;
- **Sample** is the biological replicate and must remain available for quality assessment and Condition-level inference;
- **Condition** is the biological cohort/treatment of interest and must not be silently regressed out.

The node cannot prove these roles from names. It must label the controls `technical_batch_key`, `technical_categorical_covariates`, and `technical_continuous_covariates`, reject overlap with explicitly supplied Sample/Condition keys, and warn that confounding prevents guaranteed separation of technical and biological variation. A one-level batch means no cross-batch integration, but the scVI representation remains computable; this is an explicit warning and status downgrade rather than a runtime gate.

Categorical columns require complete, nonblank values. One-level or singleton strata are disclosed as noninformative/unstable but remain expert-selectable. Continuous covariates require finite numeric values; a constant covariate is likewise noninformative and warned rather than rejected. `size_factor_key`, when supplied, must be finite, strictly positive, and on the linear scale; the official default is to use library size, so an external size factor should remain advanced and exceptional. Duplicate use of one observation column in incompatible roles must fail.

## Model and training controls

The current official constructor defaults are `n_latent=10`, `n_layers=1`, `dropout_rate=0.1`, `dispersion="gene"`, and `gene_likelihood="zinb"`. The pre-refactor node overrode three of these (`n_layers=2`, `dropout_rate=0.3`, `dispersion="gene-batch"`) without reporting why; the refactored defaults now match the public constructor. PyTorch's actual dropout domain is the closed interval `[0, 1]`, and scvi-tools 1.5.0.post1 constructs an SCVI model at `dropout_rate=1.0`; that endpoint is computationally valid but drops every affected activation and therefore warrants a degeneracy warning rather than an artificial `0.99` UI ceiling. `gene-batch` is meaningful only with a declared batch. `gene-label` should not be exposed because this node does not register labels; a one-category implicit label makes the option misleading. `gene-cell` has a much larger parameter/resource footprint and remains outside this cohesive interface.

`SCVI.train` accepts an automatic epoch heuristic when `max_epochs=None`, explicit CPU/GPU/MPS/auto acceleration and devices, train/validation fractions, minibatch size, and early stopping. Official `SCVI.train` defaults early stopping to false, although enabling it is supported when validation is configured. When enabled, the monitored validation metric, patience, minimum delta, and validation frequency are part of the result and must be explicit. The actual number of epochs, train/validation cell counts, stopping reason where available, and final/best finite loss values should be reported from `model.history`/trainer state.

The audited 1.5.0.post1 `validate_data_split` contract is exactly `0 < train_size <= 1`; a real CPU smoke with `train_size=1.0` completes when the adapter does not force a validation loop. Values below 0.5 are small-training-set choices and `1.0` leaves no validation holdout, so both warrant interpretation/training warnings but not wrapper rejection. At `1.0`, `check_val_every_n_epoch` must be `None`; explicitly requesting early stopping is a real configuration error because the declared validation monitor has no validation data. `n_layers` is a positive architecture depth with no official upper bound; 21 layers construct successfully in the audited backend. Fixed UI maxima for latent width, layers, epochs, or batch size are not backend contracts and must be removed.

All scvi-tools models generally run faster on GPU and use a detected GPU under `accelerator="auto"`. A seed controls NumPy/Torch initialization through `scvi.settings.seed`, but fixed seed does not establish bitwise equality across devices, CUDA kernels, library versions, or hardware. The selected accelerator/device, actual model device, precision if configurable, and reproducibility limitation must be disclosed. Multi-device distributed training should remain unsupported until its model and summary parity are tested.

## Latent and auxiliary outputs

`get_latent_representation()` returns an `n_obs x n_latent` representation, normally the posterior mean, suitable for a downstream neighbor graph. Exact shape/alignment and finite values are hard output requirements. A constant latent dimension is degenerate for downstream geometry but is still a typed finite backend result; count and warn rather than discard the entire trained model. It is not corrected gene expression and does not itself prove that Technical batch effects were removed or biology preserved.

Posterior distribution parameters (`qzm`, `qzv`) are model diagnostics/advanced inference artifacts, and MDE is a two-dimensional visualization transform. Neither is required to fit the model or expose its canonical integration embedding. Combining them behind booleans in the training node adds unrelated result modes and extra key collisions. They should be removed from the core integration interface and, if retained, implemented as separate model-consuming analysis/visualization nodes.

The Luecken et al. benchmark assesses integration through both batch-removal and biological-conservation metrics. A training loss or visually mixed embedding is not enough. The summary may report training and structural diagnostics but must not state “batch removed,” “conditions preserved,” “integrated successfully,” or “biological populations recovered.”

## Model lifecycle and persistence

scvi-tools officially saves model weights/configuration with `model.save()` and reconstructs with `SCVI.load()`, optionally with AnnData. The current custom `SCVIModel` wire instead holds a live trained model and its registered AnnData in memory. It is process-local, non-serializable through a saved workflow, tied to exact variable order and registered observation identities, and vulnerable to mutation of its publicly exposed registered AnnData.

If the live object remains, its interface must hide the mutable registered object, own a private aligned training copy, serialize all model calls through a lock, validate trained/attached state, and expose immutable training metadata. Saving/loading model weights is a separate explicit filesystem node, not a hidden side effect of Integration.

## Summary and code implications

The strict JSON summary should include the explicitly selected **count-intended** source, resolved expression state/evidence, integer-like status, zero-total gene count, dimensions/totals, technical-covariate roles and level sizes, architecture/likelihood/dispersion, train/validation policy, actual epochs and loss diagnostics, seed, requested and actual device, latent shape/finite/constant-dimension checks, warnings/limitations, references, and versions. Methods text must not say “declared raw UMI” unless that state is actually proven; numeric non-negativity alone is not proof. The report should explicitly say the model/embedding support exploratory representation and are not Sample-level Condition inference.

Equivalent source can reproduce the scientific operation: copy AnnData, validate the selected finite non-negative count-intended values and covariates, disclose their resolved state/evidence and numeric audit, register with `setup_anndata`, construct/train with every resolved value explicit, store the latent representation, and return `(adata, model)`. The source does not serialize the already-trained weights and need not recreate the plugin-specific wrapper/history. Re-executing it retrains a stochastic model; same seed is reproducibility metadata, not a bitwise guarantee.
