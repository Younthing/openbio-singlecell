# Add trajectory, velocity, and CNV plots

Type: task
Status: resolved
Blocked by: 01

## Scope

- PAGA, Diffusion Spectrum, DPT Gene Trend, Velocity Dynamics/Ranking, CNV Heatmap, and CNV Score Plot adapters.

## Acceptance criteria

- Graph, pseudotime, kinetic, genome-window, and score axes are validated and no new fate/tumor/causal claim is introduced.

## Answer

- Added seven domain-specific, read-only companion adapters in their existing scientific categories:
  `OpenBioSingleCellPAGAPlot`, `OpenBioSingleCellDiffusionSpectrumPlot`,
  `OpenBioSingleCellDPTGeneTrendPlot`, `OpenBioSingleCellVelocityDynamicsPlot`,
  `OpenBioSingleCellVelocityGeneRankingPlot`, `OpenBioSingleCellCNVHeatmapPlot`, and
  `OpenBioSingleCellCNVScorePlot`.
- Every adapter emits `plot`, strict `summary`, and standalone equivalent `code`, and is registered as a
  one-shot Worker operation. Plot Workers leave AnnData, VelocityState, CNVState, and table artifact files unchanged.
- Validation binds plots to the exact stored graph/partition/diffusion/pseudotime, recovered-fit, genomic-window,
  score, observation, feature, and content fingerprints. The velocity-ranking and CNV-score producers now retain
  the table/state evidence needed for fail-closed downstream plotting.
- Rendering remains descriptive: no PAGA direction, DPT measured-time/causal interpretation, velocity fate/driver
  inference, absolute-copy-number/tumor cutoff, or Sample-level Condition test is introduced.

## Comments

- 2026-09-04: Claimed; implementation starts after the shared convention is frozen.
- 2026-09-04: Completed red-green slices for all seven schemas, owned operations, standalone-code parity,
  tamper rejection, and Worker artifact round-trips. Focused trajectory/velocity/CNV regression set: 94 passed,
  1 skipped; targeted Ruff and `git diff --check` pass.
- 2026-09-04: Representative Velocity phase portrait, CNV group-mean heatmap, PAGA graph, and Diffusion Spectrum
  PNGs were visually inspected. Labels, chromosome boundary, color semantics, and PAGA/Diffusion legends passed.
