# Changelog

## Unreleased

- Add browser file selection and node-level drag-and-drop uploads to the H5AD and 10x H5 input nodes.
- Add an interactive Core Study Parameters node whose six selections are ordinary `STRING` outputs that can be
  connected and fanned out wherever useful, without coupling analysis nodes to an opaque combined design object.
- Replace the two development examples with five production starting-point templates for QC, layer-aware
  clustering and marker discovery, sample composition, scVI batch integration with a reusable model contrast, and
  an end-to-end multi-sample best-practice workflow built entirely from the existing node catalog.
- Add the `Single-Cell Best Practice` template with hard and sample-wise robust QC, Scrublet, full-gene raw/counts
  snapshots, 5,000-HVG scVI clustering, parallel marker evidence, provisional CellTypist annotation, and a
  `nonDM_ED` versus `Normal` PyDESeq2 pseudobulk branch for one explicitly selected population.
- Keep the existing pseudobulk and PyDESeq2 nodes compatible with the current Pertpy API while preserving their
  public node inputs and applying aggregate thresholds within the pseudobulk node.
- Add a generated 768 x 768 JPEG cover for every workflow and release contracts that keep template names, covers,
  node schemas, reviewed parameters, and required metadata aligned.
- Document that template values must be reviewed for each study and that the optional scVI model output is an
  ephemeral in-memory workflow object rather than persisted model weights.

## 0.2.0 - 2026-08-24

- Add notebook-derived AnnData preparation, QC correction, Harmony/scVI, clustering, annotation, pseudobulk, enrichment, communication, abundance, regulatory, trajectory, RNA velocity, and inferCNV nodes.
- Complete the notebook-backed Python workflows for pySCENIC, Cassiopeia plasticity, cNMF assignment, Augur state, scVI posterior/MDE, velocity gene ranking, sample composition, and gene-set analysis.
- Keep frequently tuned analysis choices visible while moving storage keys, output column names, seeds, and iteration limits to advanced inputs.
- Require dataset-specific labels explicitly and keep each node's natural domain output explicit.
- Replace the catch-all single-cell result wire with concrete table, plot, and structured-summary contracts so
  incompatible CSV/PNG/table links are rejected on the canvas.
- Preserve trained scVI models and pySCENIC networks as concrete reusable outputs instead of retraining or asking
  users to pass temporary intermediate files.
- Make scVI's counts source and DE mode explicit, keep optional latent distributions/MDE advanced and off by
  default, and preserve pySCENIC's optional `HVG ∪ TF` inference-gene selection.
- Reuse solved Cassiopeia trees through a dedicated reconstruction node rather than rebuilding the same tree in
  each lineage analysis node.

## 0.1.0 - 2026-08-23

- Add the first in-memory single-cell analysis node set for ComfyUI.
- Add bounded node result previews and explicit H5AD, CSV, and PNG outputs.
- Add a deterministic 600-cell demonstration AnnData generator.
- Use explicit `OPENBIO_ANNDATA` and `OPENBIO_SINGLE_CELL_RESULT` workflow contracts.
- Package validated example workflows from the single `example_workflows` source.
- Add standard-plugin and paired OpenBio frontend launch instructions.
- Add reproducible frontend production-license notices and release artifact checks.
- Add Windows and Linux install and start scripts.
- Add a Manager install entry point that generates the local demo without runtime dependency installation.
- Support pip-less `uv venv` environments in the source installation scripts.
