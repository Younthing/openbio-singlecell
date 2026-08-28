# OpenBioSingleCellAugurResults — official usage research

## Audited baseline

The current node accepts any AnnData, looks up `uns[result_key]`, assumes a dict, and converts one of three values to a
DataFrame. It does not validate that Augur produced the data, the selected table schema, parameter/input provenance,
row family, content integrity, or software version. `full_results` is documented by Pertpy as a DataFrame in 1.3.0
but historically described as a dict; `DataFrame(stored[value]).reset_index()` can silently reshape malformed input.

This node performs no Augur computation. Its legitimate role is an adapter over a typed, validated result artifact,
not a second analysis module and not an AnnData `uns` parser.

## Official Pertpy result shapes

Pertpy 1.3.0 sources:

- https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Augur.html
- https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_augur.py

`Augur.predict` returns `(updated_adata, results)`. Versioned source constructs:

- `summary_metrics`: DataFrame whose columns are population names and whose rows are mean metrics;
- `full_results`: DataFrame with `idx`, `augur_score`, `folds`, and `cell_type`;
- `feature_importances`: DataFrame with gene, importance, subsample, fold, and population fields;
- per-population internal lists of cross-validation dicts, which are implementation detail and should not be exposed.

The analysis adapter must canonicalize and validate these once. A results-view node should merely select a canonical
table from the resulting artifact, preserving analysis summary/provenance and content fingerprint.

## Scientific reporting practice

The priorities table ranks exploratory population responsiveness by classifier AUC; the cross-validation table shows
dispersion across cell subsamples/folds; feature importance is model-dependent and not causal. None is a p-value or
replicate-aware Condition result. Sample/Technical batch audits and annotation status belong to the parent artifact
and must remain attached to every selected view.

## References to propagate

The view propagates, rather than reinvents, the parent references:

- Skinnider MA, et al. https://doi.org/10.1038/s41587-020-0605-1
- Squair JW, et al. https://doi.org/10.1038/s41596-021-00561-x
- Heumos L, et al. https://doi.org/10.1038/s41592-025-02909-7

## Report/code implications

The output summary must be the complete parent scientific summary plus a small deterministic `selected_view`
section. Generated code is a portable selector/validator over the canonical result mapping, returning
`(selected_dataframe, summary_dict)`; it must not reconstruct or rerun Augur.

