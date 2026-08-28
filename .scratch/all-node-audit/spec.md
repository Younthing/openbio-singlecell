# All-node scientific audit and refactor

## Goal

Audit every registered OpenBio single-cell node, grouping only nodes with a genuinely shared scientific operation or implementation seam. Before changing a node, create two node-local records:

1. `research/<node-id>/official-usage.md`: official examples when available, otherwise a community example plus method literature.
2. `research/<node-id>/design-review.md`: scientific correctness, interface depth, parameter policy, cohesion/coupling, and a keep/merge/delete/enhance decision.

## Required analysis-node outputs

Every analysis or transformation node must expose its primary result plus:

- `summary`: an `OPENBIO_SINGLE_CELL_SUMMARY` artifact whose `.summary` payload (and preview/API representation) is strict JSON-serializable data containing methods text, results text, key results, parameters, warnings, package/method references, and runtime software versions.
- `code`: a `STRING` containing an equivalent Python function with the resolved parameters represented explicitly.

Equivalent source reproduces the scientific primary result and scientific input/error invariants. It does not recreate plugin-only `openbio_singlecell.analysis_history` bookkeeping; the adjacent `summary` artifact is the portable provenance record.

Nodes whose primary result is already the structured `summary` artifact expose that artifact once, enriched to the same reporting schema, plus `code`. Pure input adapters, output adapters, and parameter/UI nodes are audited but receive these ports only if they perform an analysis.

## Scientific invariants

- `Sample` is the biological replicate and condition-level inference unit.
- Cluster marker evidence remains exploratory annotation evidence; it is not a condition contrast.
- Provisional annotation remains distinct from curated annotation.
- Count-dependent operations disclose and validate their expression source.
- Dataset-specific filtering thresholds are reported as user choices, not universal biological constants.
- Transforming nodes copy inputs; report/plot nodes are read-only.

## Implementation policy

- Centralize report construction, JSON normalization, runtime version capture, and code rendering behind one in-process module seam.
- Keep scientific result extraction local to each node so key results remain method-specific.
- Frequently tuned, result-defining parameters remain visible. Storage keys, output names, seeds, and performance controls are advanced unless they define the scientific question.
- Tests assert through node interfaces: primary result, `summary`, and `code`.

## Batches

1. Quality control: calculate metrics, filter cells, filter genes, QC plots.
2. Data ingestion and expression-state management.
3. Normalization, feature selection, scaling, and artifact correction.
4. Representation, integration, neighborhood graph, clustering, and embeddings.
5. Annotation, cluster marker evidence, and population summaries.
6. Replicate-aware condition contrasts and differential abundance.
7. Enrichment, factorization, regulation, communication, CNV, trajectory, lineage, and velocity.
8. Result and output adapters, examples, release documentation, and final registry audit.

## Verification

- Baseline on 2026-08-28: `D:\learn\ComfyUI\.venv\Scripts\python.exe -m pytest -q` → `138 passed`.
- Each batch adds focused tests, runs Ruff, then runs the full suite in the ComfyUI environment.
