# Cassiopeia clade expansion testing: module design review

## Decision and current interface

**Retain and deepen `OpenBioSingleCellCassiopeiaExpansionTest` as one fixed-tree, within-lineage hypothesis family.** It consumes a validated tree copy, invokes and independently checks the public statistic, performs the one required multiple-testing correction, and returns evidence. It never reconstructs or mutates the tree and never aggregates across lineages.

Inputs, in order:

1. `tree` (`OPENBIO_CASSIOPEIA_TREE`);
2. positive integer `minimum_clade_size=10`; singleton families are permitted with an interpretation warning;
3. integer `minimum_depth=1`;
4. `fdr_threshold=0.05`.

Outputs are `table`, `summary`, and `code`. A fractional clade cutoff is intentionally absent: the official interface and null-family boundary use an integer count, while a fraction adds a second competing policy and rounding rule.

The node validates artifact/provenance/topology fingerprints, enumerates every node and eligibility before the backend call, calls `compute_expansion_pvalues(copy=True)` on a defensive copy, proves source-tree and global-state non-mutation, independently recomputes each eligible probability, and rejects missing/non-finite/out-of-range/disagreeing values. It applies a local stable Benjamini-Hochberg implementation to the whole eligible family. Output order is canonical depth-first order with exact node identity.

The strict report names the lineage/clade inference unit, family size, FDR method/threshold, number and strongest adjusted results, solver/root/topology provenance, exact versions, warnings and limitations. `code` is standalone and returns `(table, summary_dict)` with the same eligibility, formula verification and BH results.

## Legacy/current migration matrix

Legacy inputs were `tree, minimum_clade_fraction, minimum_depth, expansion_pvalue_threshold`; current inputs are `tree, minimum_clade_size, minimum_depth, fdr_threshold`. Legacy outputs had only `table`; current outputs append `summary, code`.

`minimum_depth` maps directly. Neither `minimum_clade_fraction` nor the raw `expansion_pvalue_threshold` has a semantics-preserving mapping: the former depended on runtime tree size and passed a float to an integer contract, while the latter labeled unadjusted p-values and cannot become an FDR threshold. Migration must therefore stop for explicit integer cutoff and FDR review, without mutating the graph. Once supplied, links to table output 0 remain valid and new outputs append. The migration must be atomic/idempotent across array/object links and all workflow scopes. A connected `minimum_depth` remains an INT widget port (and the clade/FDR controls retain their respective INT/FLOAT types); a same-named STRING connection is an incompatible schema, not a current no-op.

## Verification gate

Tests cover exact topology eligibility, root/leaves/shallow nodes, boundary sizes/depths, source copy semantics, independent combinatorial probabilities, backend omission/NaN/out-of-range/malicious disagreement, zero/one/tied hypotheses, trusted BH fixtures, raw-significant but FDR-nonsignificant clades, stable ordering/hash, tree tampering, strict JSON, compiled code and parity. A real `cassiopeia-mt==2.1.3` smoke verifies the public API on a solved tree.
