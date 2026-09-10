# DDRTree, rooting and evidence nodes

Type: task
Status: resolved

Implement Monocle 2 stages, native trajectory/gene plots, pseudotime differential tests, optional BEAM and precise
AnnData attachment. Use real native calls and equivalent-code verification, preserving expert controls.

## Answer

All ten registered nodes are implemented. Native execution and standalone generated-code comparisons passed for
trajectory construction, orientation, gene evidence, plots and AnnData attachment. Typed metadata and NA are retained;
new DDRTree graphs invalidate old native ordering. See ../verification.md.
