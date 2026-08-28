# Run Scrublet: module design review

## Current module

The node wraps Scanpy's Scrublet call and annotates `AnnData`, which is a coherent atomic analysis. It currently copies the input.

## Correctness and interface problems

- It implicitly uses `X`; normalized/log/scaled `X` can therefore be supplied while described as Scrublet input counts.
- It defaults to technical `batch` rather than biological `sample`/capture and does not validate the grouping column or missing labels.
- It silently makes observation names unique, changing cell identifiers rather than rejecting an ambiguous input.
- Expected doublet rate, threshold, simulation ratio, and PCA dimension are hidden entirely.
- It does not surface the effective automatic thresholds, score distributions, group results, resolved defaults, limitations, or code.

## Decision: enhance, do not merge

Keep Scrublet detection separate from `Filter Predicted Doublets`. Detection is stochastic model evidence that should be inspected and possibly re-thresholded; filtering is a reversible-looking but destructive dataset decision. Keeping the seam also permits prediction columns from other tools.

## Target interface

- Inputs: `adata`; explicit `X`/`raw`/layer count source; Sample/capture grouping key; expected doublet rate; explicit automatic/manual threshold mode and manual value; simulation ratio; PCA components; optional neighbor override; random seed.
- Outputs: annotated `adata`, structured `summary`, equivalent `code`.
- Visible scientific choices: grouping key and expected rate. Manual threshold is visible because it directly defines calls.
- Advanced controls: simulation ratio, PCA components, seed, and layer name selected through the dynamic source.
- Hidden Scanpy defaults are fixed explicitly in generated code and disclosed in summary: rate SD `0.02`, UMI subsampling `1.0`, Euclidean metric, variance normalization, no log transform, mean centering, automatic neighbor count, and default approximate-neighbor selection.

## Invariants

Require non-empty data, present/non-missing Sample labels when grouped, finite non-negative input with positive signal, more selected features than requested PCs, and at least `n_prin_comps + 1` cells in every Sample/capture. Fractional values and duplicate observation identifiers are advisory states: run the backend, preserve positional identity exactly, and disclose that the input is not verified raw UMI counts / that Cell IDs are non-unique. Permit one requested principal component when the backend dimensions support it. Scanpy internally filters cells, genes, and variable genes, so a remaining PCA-dimension failure is translated into a domain-specific instruction to reduce PCs or provide a larger capture. Build an isolated working `AnnData` from the selected matrix, run Scanpy there, and transfer only Scrublet annotations/provenance onto a copy of the original object. This prevents the wrapper's preprocessing from changing the caller's chosen expression state.

No analysis-history or Raw-binding proof is consulted; explicit source selection is the user's declaration.

## Freeze-review validation correction

- Validate `threshold` only when `threshold_mode="manual"`, require finiteness rather than an undocumented `[0, 1]` range, and emit a machine-readable advisory when it lies outside the conventional score interval.
- Require a finite positive simulated-doublet ratio.
- Keep generated-code validation order equivalent to runtime: non-finite expression first, then negative values, then absence of positive signal. This preserves the module's documented error modes.
- Automatic mode fixes the generated/backend threshold to `None` and ignores the inactive manual value.

## Report contract

Report input dimensions/source, grouping sizes, score distribution, predicted count/rate overall and per group, requested/effective threshold(s), all effective model parameters, warnings/limitations, Wolock and Scanpy references, and plugin/Scanpy/AnnData/NumPy/Pandas/SciPy versions.
