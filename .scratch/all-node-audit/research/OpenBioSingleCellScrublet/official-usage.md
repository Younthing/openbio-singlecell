# Run Scrublet: official usage research

Researched: 2026-08-28

## Official Scanpy usage

`scanpy.pp.scrublet` predicts doublets from observed transcriptomes and simulated doublets. The official API says the input should be a raw, unnormalized count matrix and works best for a single Sample or similar Samples from the same experiment. `batch_key` partitions the calculation. The standard output is:

- `obs['doublet_score']`: continuous score;
- `obs['predicted_doublet']`: boolean automatic/manual-threshold call;
- `uns['scrublet']`: simulated scores, parents, effective parameters, and threshold information (nested by batch when `batch_key` is used).

Official API: https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pp.scrublet.html

The official Scanpy clustering tutorial calls `sc.pp.scrublet(adata, batch_key="sample")` before normalization and clustering:

- https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering.html

The repository pins Scanpy `>=1.12.3,<1.13`. The installed 1.12.3 signature uses `random_state`; current development documentation names the forward RNG parameter `rng`. Implementation must follow the supported pinned signature and disclose the runtime version.

Automatic thresholding uses scikit-image. The distributable dependency must therefore request Scanpy's `scrublet` extra; relying on an incidental package in a developer environment would make the default node fail after installation.

## Original method and best practices

The Scrublet project documents raw UMI counts and recommends running each Sample separately because doublets arise through random co-encapsulation within a capture. It also recommends checking the score threshold and the location of predictions in an embedding instead of treating automatic calls as unquestionable truth.

- Official project and example links: https://github.com/swolock/scrublet
- Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. *Cell Systems*. 2019;8(4):281-291.e9. https://doi.org/10.1016/j.cels.2018.11.005
- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0

## Scientific implications

- The expression source is a result-defining input and must be explicit, non-negative, finite, and contain usable signal. Fractional values are disclosed as not verified raw UMI counts but are not rejected when Scanpy can run them.
- Default grouping is `sample`; a blank key is an explicit global run and should carry a warning for multi-Sample data.
- Preserve cell identifiers; do not silently call `obs_names_make_unique()` because that changes Sample/cell identity.
- Expose expected doublet rate and an explicit automatic/manual threshold mode. Keep simulation ratio, principal components, neighbor override, and seed as advanced controls.
- Treat a missing effective threshold in `uns['scrublet']` (including any missing batched threshold) as analysis failure. Scanpy can otherwise leave an all-false prediction column after automatic threshold selection fails.
- Report predicted count/rate, score distribution, per-Sample counts/rates, automatic or manual effective thresholds, input source, resolved hidden defaults, references, and software versions.
- Predictions are QC evidence, not ground truth. Filtering remains a separate node.

The official wording says Scrublet “works best” on raw unnormalized counts; it does not establish integer dtype or unique observation labels as universal backend preconditions. Duplicate Cell IDs therefore proceed with an identity warning because row order is preserved, and finite non-negative fractional input proceeds with a count-suitability warning. PCA dimensions, non-negative finite signal, complete grouping assignment, and documented output fields remain true computational contracts.

## Manual-threshold correction

The official Scanpy/Scrublet interface types `threshold` as `float | None` and defines the call by comparing observed scores with that threshold; it does not specify a `[0, 1]` validation range. Automatic mode uses `None`, so an otherwise unused manual widget value must not block it. In manual mode the wrapper requires only a finite threshold. Values outside the observed or conventional score range are allowed and disclosed because they intentionally produce an extreme all/none call. The simulated-doublet ratio remains a true backend parameter and must be finite and greater than zero.
