# Changelog

## Unreleased

- Add 45 domain-specific companion Plot nodes across correction, integration, clustering, marker evidence,
  annotation, factorization, differential testing, enrichment, regulation, trajectory, velocity, copy number,
  lineage, and cell prioritization. Companion nodes consume retained evidence read-only and consistently return a
  PNG, strict summary, and equivalent code without rerunning the scientific analysis.
- Retain complete Harmony/scVI training diagnostics, Milo neighborhood topology and coordinates,
  scCODA/tascCODA posterior diagnostics, and Augur fold predictions in validated typed artifacts so convergence,
  neighborhood, posterior/trace, and ROC plots remain auditable downstream.
- Complete the core AnnData plotting path with Sample-level composition bars, flavor-aware HVG diagnostics, PCA
  explained and cumulative variance, grouped QC distributions, and a general stored-embedding renderer that colors
  cells by observation metadata or one explicitly selected expression-source gene. Grouped gene-expression plots
  retain their existing dot, matrix, track, and violin modes under a clearer display name.
- Make the public surface current-only: remove eight retired or compatibility-only nodes, all hidden direct-call
  parameter aliases, the `ora` installation-extra alias, old analysis-history readers, incomplete scVI artifact
  provenance defaults, and automatic browser-side workflow migration. The registry now contains exactly 144
  non-deprecated nodes; saved graphs and process-local artifacts must match their current schemas.
- Make result previews and H5AD/CSV/PNG output adapters failure-atomic. PNGs are fully decoded before publication,
  H5AD files are reopened to verify shape and axis identity, CSV files preserve explicit row-index identity, and
  non-overwrite publication fails closed under a late filename collision.
- Refactor cluster-marker analysis into a canonical ranked marker table plus its exact tested-gene universe. Marker
  statistics retain the complete hypothesis family for per-group Benjamini-Hochberg correction, preserve true
  nonzero fractions for arbitrarily small signals, disclose ambiguous or non-count-like expert inputs in warnings,
  preflight bounded working memory, and use streaming content fingerprints. Filter Marker Genes accepts only a direct Marker Genes pair so
  threshold revisions remain independently auditable.
- Make UMAP and marker-expression plots read-only, expression-source-aware report artifacts with strict `summary`
  and equivalent `code` outputs. UMAP uses the first two coordinates of any validated >=2-D embedding; marker
  dot/matrix/track/violin plots construct only the selected private matrix, disclose missing values and scale
  evidence, and isolate seeded violin jitter without advancing NumPy's global random state.
- Separate annotation into provisional CellTypist model evidence, explicit-resource Marker ORA evidence, and
  reviewed cluster-label mapping. CellTypist requires a preinstalled local model and an explicit verified expression
  state; Marker ORA uses the decoupler 2.x `mt.query_set` API, retains the full group-by-source evidence grid, and
  never commits labels. The former all-in-one Marker ORA Annotation node is no longer registered.
- Refactor PCA, Harmony, scVI, Neighbors, UMAP, t-SNE, force-directed layouts, Leiden, resolution sweeps, and PCA
  metadata diagnostics into separate validated analyses with strict `summary` and equivalent `code` outputs. Named
  graph bundles and representation dimensionality are checked end to end; Technical batch variables stay distinct
  from biological Samples, clustering stability is disclosed, and PCA metadata associations use Sample-level
  aggregation.
- Use `n_dimensions` consistently; explicit `X` and non-PCA learned representations use `0` for
  all available dimensions, preserving Scanpy's prior representation semantics. Force-directed layouts now default
  to Fruchterman-Reingold with random initialization, require an
  explicit validated PAGA prerequisite, and fail instead of silently replacing unavailable ForceAtlas2.
- Refactor normalization and feature-modeling nodes around explicit expression states, copy-on-write layer outputs,
  flavor-specific validation, strict JSON summaries, and equivalent Python source. Raw snapshots are no longer
  created implicitly by normalization or log transformation.
- Split the former all-in-one cNMF node into cNMF Rank Survey and cNMF Consensus Programs, connected by the
  process-local `OPENBIO_CNMF_RUN` contract. Rank stability/error is now disclosed before K selection, continuous
  usages and both GEP matrices are retained, and the unrelated RFC/argmax hard-label behavior plus shared
  worker/directory controls are removed. Survey resource use is preflighted against a conservative fixed budget,
  while consensus density and neighborhood choices remain explicit. The adapter uses the method authors'
  `cnmf==1.7.1` file-backed API inside an artifact-owned managed temporary workspace; incompatible OmicVerse is not
  installed. Current workflows connect Rank Survey explicitly before Consensus Programs.
- Remove in-process pySCENIC execution. pySCENIC 0.12.1 runs only in its pinned external Python 3.10/container
  environment; OpenBio imports a complete versioned, hash-bound expression/adjacency/regulon/AUCell bundle into an
  immutable `OPENBIO_SCENIC_RESULT` artifact. Downstream specificity, binarization, and membership nodes consume
  that typed evidence, and no node patches NumPy aliases, downloads ranking resources, or starts a container engine.
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
- Add the initial notebook-derived regulatory, lineage, cNMF, perturbation, scVI, velocity, composition, and
  gene-set analysis adapters; their current audited contracts are described under Unreleased above.
- Keep frequently tuned analysis choices visible while moving storage keys, output column names, seeds, and iteration limits to advanced inputs.
- Require dataset-specific labels explicitly and keep each node's natural domain output explicit.
- Replace the catch-all single-cell result wire with concrete table, plot, and structured-summary contracts so
  incompatible CSV/PNG/table links are rejected on the canvas.
- Preserve trained scVI models as concrete reusable outputs; pySCENIC was subsequently moved to a validated external
  bundle-import seam because its supported runtime is incompatible with the host Python/NumPy baseline.
- Make scVI's counts source and DE mode explicit and keep optional latent distributions/MDE advanced and off by
  default. pySCENIC gene selection and resource choices are now disclosed by the imported external-run manifest.
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
