# OpenBioSingleCellSCENICTFModules — official usage research

## Audited baseline

The current node calls `pyscenic.utils.modules_from_adjacencies` on the raw GRN adjacency and a newly selected
expression matrix. It then reports only `transcription_factor,module,gene`, discarding target weights, module
context, activation/repression status, thresholds, and motif evidence. The resulting objects are **pre-cisTarget
co-expression modules**, not the motif-pruned final SCENIC regulons scored by AUCell. Calling them generic “SCENIC TF
modules” next to imported final results is scientifically ambiguous.

The node also hides all method-defining pySCENIC defaults and may recalculate correlations from an expression state
different from the external `ctx` run. Its `ScenicNetwork` input comes from the retired in-process Run node and
contains only the GRN stage, so it cannot establish final regulon membership.

## Official pySCENIC module and regulon semantics

Official sources and tutorial:

- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/utils.py
- https://pyscenic.readthedocs.io/en/stable/tutorial.html
- https://pyscenic.readthedocs.io/en/stable/installation.html
- https://doi.org/10.1038/s41596-020-0336-2

The official pre-pruning call is:

```python
from pyscenic.utils import modules_from_adjacencies

candidate_modules = modules_from_adjacencies(
    adjacencies,
    expression,
    thresholds=(0.75, 0.90),
    top_n_targets=(50,),
    top_n_regulators=(5, 10, 50),
    min_genes=20,
    absolute_thresholds=False,
    rho_dichotomize=True,
    keep_only_activating=True,
    rho_threshold=0.03,
    rho_mask_dropouts=False,
)
```

These candidate modules are built by weight percentiles, top targets, and top regulators. Pearson correlations
split activating/repressing candidates, the TF is added, and modules below `min_genes` are removed. The subsequent
`ctx` stage tests cis-regulatory motif enrichment and prunes targets. Only those motif-pruned objects are final
regulons, and only those are inputs to AUCell.

For a result/reporting node downstream of external pySCENIC, recomputing candidate modules is neither equivalent nor
necessary. The authoritative membership is the safely imported, motif-pruned regulon artifact. A final membership
table should preserve regulon name, TF, sign/context, target, target weight, and motif-support count/identifiers when
available. Multiple final regulons may share a TF; TF identity alone is not a unique module key.

## Scientific interpretation

Final regulons combine expression-derived co-expression and motif prior knowledge. Membership remains inferred and
does not prove physical binding or causal regulation. The resource organism, genome build, motif/ranking database,
input gene universe, and external run parameters determine which targets can appear. A per-TF query is a view of the
complete imported family, not a new statistical test.

## References to emit

- Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature
  Protocols*. 2020;15:2247-2276. https://doi.org/10.1038/s41596-020-0336-2
- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463
- Moerman T, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks.
  *Bioinformatics*. 2019;35:2159-2161. https://doi.org/10.1093/bioinformatics/bty916

## Report and generated-code implications

The report must name the table as final motif-pruned regulon membership, disclose complete/query counts, TF and
target counts, shared/multi-regulon TFs, sign/context and motif evidence, target-weight ranges, external resource/run
provenance, limitations, references, and dynamic versions. Equivalent code should project the verified immutable
artifact without importing pySCENIC or rereading expression/resources and return `(table, summary_dict)`.

