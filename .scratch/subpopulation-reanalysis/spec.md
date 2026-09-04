# Subpopulation reanalysis

## Objective

Make the existing single-cell nodes support a complete standard-RNA subpopulation workflow without adding a
general-purpose AnnData query language or a monolithic analysis node.

## Confirmed workflow

1. Select one or more parent populations with `Subset Observations`.
2. Materialize the selected observations' Raw snapshot as a fresh active AnnData.
3. Recreate the explicit expression snapshots/layers needed by existing preprocessing nodes.
4. Recompute HVG, representation, graph, embedding, clustering, marker evidence, and annotation within the
   subpopulation.
5. Transfer the reviewed subtype annotation back to the full parent AnnData.

## Requirements

### Raw Snapshot to AnnData

- Add one scientific node with one required `adata` input and standard `adata`, `summary`, and `code` outputs.
- Its matrix restoration is based on `adata.raw.to_adata()` after observation subsetting.
- Preserve the selected observation axis and `obs` metadata exactly, and restore `raw.X` plus `raw.var` as the
  active matrix and variable metadata.
- Return a fresh analysis object without inherited layers, Raw, embeddings, pairwise observation graphs, variable
  embeddings, pairwise variable graphs, or stale analysis results.
- Clear non-OpenBio `uns` analysis results, retain the existing `openbio_singlecell` provenance metadata (including
  source, annotations, warnings, and history), then append this operation's normal OpenBio history.
- Do not infer or claim that the Raw snapshot contains counts.
- Reject missing Raw, zero selected observations, or zero Raw variables with an actionable error.
- Do not mutate the input; emitted equivalent code must produce the same primary AnnData behavior while following
  the repository convention of omitting plugin-specific history bookkeeping.

### Annotation round trip

- Keep the existing `Merge Observation Annotations` interface.
- When its source column carries OpenBio annotation provenance, transfer that provenance to the target column.
- Reject different nonempty OpenBio source provenance. External inputs without sufficient provenance remain usable
  with a disclosed limitation.
- Continue transferring only the requested observation annotation, never subset-only expression, embeddings, or
  graphs.

### Workflow template

- Add one reviewed subpopulation-reclustering template built from ordinary nodes.
- Include the required diagnostic/marker plots and explicit saves for both the subpopulation object and the updated
  parent object.
- Keep full-gene marker evidence by using `Highly Variable Genes` with `subset=false` on the standard path.
- Do not add automatic all-population loops, generic selection algebra, or arbitrary AnnData concatenation.

## Verification seams

- Scientific operation seam: execute the Raw conversion and annotation merge operations through the same public
  operation registry used by nodes.
- Node seam: schema/registration exposes the new node with the standard output contract.
- Workflow seam: generated packaged workflow contains and correctly connects the reviewed node sequence.
- Repository seam: focused tests, complete pytest suite, Ruff, workflow generation checks, and two-axis code review.

## Fixed point

`ca47278930ebf2055349e1e2cd746f52a2bac78f`
