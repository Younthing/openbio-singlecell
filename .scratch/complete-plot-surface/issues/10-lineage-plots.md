# Add Cassiopeia lineage plots

Type: task
Status: resolved
Blocked by: 01

## Scope

- Lineage QC, Tree, Clade Expansion, and Effective Plasticity Plot adapters.

## Acceptance criteria

- Typed character/tree identity, leaf/cell alignment, topology, annotations, and table provenance are validated without reconstructing the tree.

## Answer

- Added `OpenBioSingleCellCassiopeiaLineageQCPlot`, `OpenBioSingleCellCassiopeiaTreePlot`, `OpenBioSingleCellCassiopeiaExpansionPlot`, and `OpenBioSingleCellCassiopeiaPlasticityPlot` in the lineage category.
- Added matching one-shot Worker operations `openbio.node.cassiopeialineageqcplot`, `openbio.node.cassiopeiatreeplot`, `openbio.node.cassiopeiaexpansionplot`, and `openbio.node.cassiopeiaplasticityplot`.
- QC, expansion, and plasticity tables now carry content and upstream artifact/topology fingerprints; nullable table fields and plasticity scores round-trip through the portable codecs.
- Renderers validate exact artifact types, rooted topology, leaf/AnnData alignment, categorical annotation, canonical table schemas and producers, per-row evidence, and immutable inputs. Generated code reproduces each PNG without rerunning reconstruction or inference.
- Narrow verification: 35 lineage tests passed; Ruff and `git diff --check` passed. Representative Cassiopeia Tree PNG visual QA passed with readable topology depth, leaf labels, and legend and no clipping/overlap.

## Comments

- 2026-09-04: Claimed; implementation starts after the shared convention is frozen.
- 2026-09-04: Resolved with four domain-specific companion plots and portable, fingerprint-bound table evidence.
