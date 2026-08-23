# Changelog

## Unreleased

## 0.2.0 - 2026-08-24

- Add notebook-derived AnnData preparation, QC correction, Harmony/scVI, clustering, annotation, pseudobulk, enrichment, communication, abundance, regulatory, trajectory, RNA velocity, and inferCNV nodes.
- Complete the notebook-backed Python workflows for pySCENIC, Cassiopeia plasticity, cNMF assignment, Augur state, scVI posterior/MDE, velocity gene ranking, sample composition, and gene-set analysis.
- Keep frequently tuned analysis choices visible while moving storage keys, output column names, seeds, and iteration limits to advanced inputs.
- Require dataset-specific labels explicitly and keep each node's natural domain output explicit.
- Replace the catch-all single-cell result wire with concrete table, plot, and structured-summary contracts so
  incompatible CSV/PNG/table links are rejected on the canvas.

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
