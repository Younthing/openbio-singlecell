# Complete domain-specific plotting surface

Fixed point: `cce793af406f5be5cb311769fe80ecf534a3630e`

## Goal

Give every analysis result that materially benefits from visualization a discoverable, domain-specific companion Plot node, while keeping transformations, inference, and rendering separate.

## Uniform public Interface

- Plot nodes live in the existing scientific category of the result they interpret. Only genuinely cross-domain AnnData renderers remain in `visualization`.
- One scientific result contract maps to one Plot node. Different views of that same result use a closed `DynamicCombo`; different result semantics use different nodes.
- Inputs are one or more exact AnnData, canonical `TableResult`, or existing/new typed artifacts. Table consumers require an allowed `source.operation`, exact schema, and relevant provenance/axis fingerprints; they never infer column roles.
- Plot nodes are read-only and never rerun preprocessing, model fitting, statistical testing, or resource download.
- Outputs are always `plot`, strict `summary`, and standalone equivalent `code` through the existing PlotResult/one-shot Worker path.
- Rendering uses current in-process dependencies, object-oriented Matplotlib/Agg or the installed method-author renderer when it is the narrowest correct Adapter. No new plotting dependency, wildcard artifact socket, arbitrary x/y/mark/facet DSL, default plot output on every producer, compatibility shim, or silent fallback is allowed.
- Every renderer owns deterministic ordering, missing-value policy, resource/pixel limits, global-state cleanup, warnings, interpretation limitations, and actual plotted-count disclosure.

## Diagnostic evidence contract changes

1. Harmony retains the complete validated objective/convergence sequence with the adjusted representation provenance. Its Plot node renders the stored sequence and never reruns Harmony.
2. scVI's native model artifact retains complete supported training/validation metric sequences with epoch identity. Its Plot node consumes that session-only artifact and preserves Worker affinity.
3. Milo publishes a portable typed result containing the canonical result table plus neighborhood membership, neighborhood graph, representative coordinates, observation/neighborhood axes, and fingerprints. Its Plot node supports both table evidence and neighborhood-on-embedding views without rerunning Milo.
4. scCODA/tascCODA publish a typed composition-model result containing their canonical table, model metadata/hierarchy, and validated ArviZ `InferenceData`/posterior diagnostics sufficient for effect, posterior, and trace views.
5. Augur's existing typed result additionally retains fold-level true labels, prediction scores or deterministic ROC points tied to population/subsample/fold. Its Plot node supports priority, CV, feature-importance, and ROC views.

## Required companion Plot coverage

### correction

- MAD Outlier Plot: metric distributions/threshold evidence and Sample-level marked fraction.
- Scrublet Diagnostics Plot: observed/simulated score distributions, thresholds, and Sample-level predicted-doublet fraction.

### preprocessing, dimension reduction, integration, clustering, diagnostics

- Retain existing QC, HVG, PCA variance, Embedding, and grouped-expression plots.
- PCA Loadings Plot; Neighbor Graph Diagnostics Plot; Harmony Convergence Plot; scVI Training Plot.
- Leiden Resolution Sweep Plot; Schist Hierarchy Plot; PCA Metadata Association Plot.

### marker evidence and annotation

- Marker Evidence Plot consumes AnnData plus the direct Marker table/universe, selects upstream-ranked top genes, and reuses grouped-expression semantics without reranking.
- CellTypist Diagnostics Plot and Cell Cycle Score Plot cover evidence not already clear from a single Embedding Plot.

### factorization

- cNMF Rank Plot renders stored stability/prediction-error evidence without selecting K.
- cNMF Programs Plot covers program usage and top-gene evidence already stored in the output AnnData.

### differential expression and abundance

- Pseudobulk QC Plot; frequentist edgeR/PyDESeq2 contrast Plot; separately worded scVI DE Evidence Plot.
- Milo Differential Abundance Plot; scCODA Differential Composition Plot; tascCODA Differential Composition Plot.
- Retain Sample Composition Plot and move it into `differential-abundance`.

### enrichment and regulatory

- Score/Activity Plot for stored AUCell, GSVA, gene-panel, drug, TF, and SCENIC activity matrices/columns, with exact producer provenance.
- Pathway Score Contrast Plot; Ranked GSEA Plot; ORA Evidence Plot; Drug Enrichment Plot.
- TF Activity Ranking Plot; SCENIC Regulon Specificity Plot; SCENIC Binarization Plot; SCENIC Regulon Membership Plot.

### trajectory, velocity, copy number, lineage, prioritization

- PAGA Plot; Diffusion Spectrum Plot; DPT Gene Trend Plot.
- Retain Velocity Stream Plot; add Velocity Dynamics Plot and Velocity Gene Ranking Plot.
- CNV Heatmap Plot and CNV Score Plot; CNV PCA remains covered by Embedding Plot.
- Cassiopeia Lineage QC Plot, Tree Plot, Expansion Plot, and Plasticity Plot.
- Augur Plot over all retained Diagnostic evidence views.

## Explicit non-goals

- No dedicated plot for runtime/config/input/output nodes; expression snapshot, subset, annotation merge, GTF mapping; Normalize/Log1p/Scale/Pearson transform; Filter Cells/Genes/Doublets; single Leiden; cluster annotation mapping; DGIdb resource; velocity preparation/moments/estimate/graph intermediate steps.
- No new LIANA network/chord view; the existing dot plot is the basic companion.
- No automatic best K/resolution/threshold, annotation curation, causal claim, Condition inference from cells, or Technical batch/Condition conflation.
- No backwards-compatibility aliases. Current workflows and examples are migrated to current schemas together.

## Public test seams

- Node schemas: IDs, categories, typed inputs, active view parameters, output order, and defaults.
- Producer contracts: complete Diagnostic evidence survives codec/runtime round-trip, validates axes/fingerprints, remains immutable, and rejects tampering.
- Owned Plot operations: scientific meaning, deterministic order, realistic failures, PNG decoding, strict summary, standalone code equivalence, global state restoration, input immutability, and bounded resources.
- Worker operations: typed artifacts and source files are unchanged; wrong producer/schema/axis/Worker identity fails before rendering.
- Example workflows and documentation expose representative defaults without requiring every view in every template.

## Completion

- All required plots and evidence contracts are registered and documented.
- Representative plots are rendered and visually inspected.
- Narrow tests, the complete Python suite, Ruff, workflow generation, release checks, and `git diff --check` pass.
- Standards and Spec reviews of the complete staged/unstaged/untracked change report no unresolved findings.
