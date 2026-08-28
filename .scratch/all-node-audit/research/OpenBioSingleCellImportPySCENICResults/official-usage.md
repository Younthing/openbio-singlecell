# OpenBioSingleCellImportPySCENICResults — official usage research

## Audited baseline

The legacy node reads `RegulonsAUC` and `CellID` from a final loom, then calls `pyscenic.utils.load_motifs`,
`pyscenic.transform.df2regulons`, and `pyscenic.export.add_scenic_metadata`. With pySCENIC 0.12.1 in the configured
NumPy 2.4 environment, importing `pyscenic.transform` fails at the removed `numpy.object` alias. More seriously,
`pyscenic.utils.load_motifs` parses two CSV columns with Python `eval`; using it on a user-supplied file is not an
acceptable trust boundary.

Current validation only checks that every AnnData observation occurs somewhere in the loom. It does not reject
duplicate or extra loom cells, validate regulon/AUC identity, verify a completed three-stage run, hash file content,
or bind outputs to the input expression/resource state. `loompy.connect(..., validate=False)` further suppresses
loom validation.

## Official result structures

Official pySCENIC 0.12.1 documentation:

- https://pyscenic.readthedocs.io/en/stable/installation.html
- https://pyscenic.readthedocs.io/en/stable/tutorial.html
- https://pyscenic.readthedocs.io/en/stable/faq.html
- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/utils.py
- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/transform.py
- https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/export.py

The `ctx` output is a multi-index motif-enrichment CSV from which final motif-pruned `Regulon` objects are formed.
Target-gene fields encode Python-literal-like collections of `(gene, weight)` pairs. The `aucell` CLI can write a
cell-by-regulon CSV directly; without `--transpose`, its first column is the cell index and its remaining columns are
regulons. A final loom commonly stores cell IDs in column attributes and AUC in `ca.RegulonsAUC`, but pySCENIC itself
opens these looms with `validate=False` and its append path also invokes the fragile binarization helper. The official
helpers are convenient inside the trusted 0.12.1 environment, but their use of
`eval` and removed NumPy aliases means the plugin importer must implement a strict, safe equivalent parser.

The supported import boundary therefore uses the official CSV AUCell output and the `ctx` multi-index CSV. This
avoids treating a broad expression-containing loom as a trusted result container. Pickled/DAT regulon objects must
never be accepted because they are code-execution artifacts; YAML is also outside this strict parser boundary.

The tagged CLI expression loader delegates CSV values to `pandas.read_csv` and does not require integer or
nonnegative counts. Raw post-QC counts remain the conventional/recommended input, but normalized or otherwise
transformed finite numeric matrices are backend-valid expert choices. Their explicit `expression_state` must be
preserved and warned rather than rejected. Likewise, the CLI has defaults for method, worker count, mode, AUC
threshold, and optional random seeds; omission can reduce reproducibility but does not make an already materialized,
content-bound result invalid. Both expression and AUCell CSV files use their first column as the index; pandas permits
an arbitrary index-column header, so only its row identities—not a preferred header spelling—are semantic. The
manifest must declare pySCENIC, ctxcore, arboreto, and Python, but additional package-version entries (for example
NumPy, pandas, Dask, and loompy) are useful provenance and must not be rejected.

## Required validation semantics

A safe import must:

- require a versioned external-run bundle manifest, resolve only canonical relative member paths beneath its
  directory, and stream-verify the SHA-256 of the expression export, all three stage outputs, and every resource
  before parsing results;
- parse JSON with finite-value enforcement, target-pair lists with `ast.literal_eval`, and Context using an explicit
  `frozenset(<set-literal>)` grammar plus `ast.literal_eval` only for the inner literal; the official CSV emits a
  `frozenset(...)` constructor expression that plain `ast.literal_eval` cannot accept;
- reject duplicate/nonblank violations in cells, genes, final regulons, TF/target identities, and adjacency edges;
  repeated motif evidence and repeated TF-target support across motif rows are legitimate, but exact duplicate motif
  records fail and TF-target weights are consolidated with pySCENIC's documented maximum-weight union rule;
- validate the official `TF,target,importance` adjacency schema, unique TF list, required motif-annotation columns,
  and Feather-v2 ranking-database framing under explicit file and edge budgets;
- require the AUCell CSV cell order to equal the AnnData observation order exactly; set equality followed by silent
  reordering is not sufficient provenance;
- require AUC columns to equal the imported final-regulon family exactly, with finite values in valid AUCell range;
- validate the manifest's input cell/gene/value hashes, organism, namespace, explicitly named numeric expression
  state, resource bundle, CLI argument arrays, exact pySCENIC 0.12.1 declaration, and other declared versions. A
  container digest and explicit `grn`/`aucell` seeds are strongly preferred but may be absent with prominent
  reproducibility warnings; only files declared and actually consumed are content-bound, while unrelated companion
  files in the directory are ignored;
- materialize immutable activity and regulon membership before returning; no live file handle survives import;
- avoid `add_scenic_metadata`, which couples the importer to broad, version-sensitive AnnData mutations.

SCENIC results are regulatory-network evidence inferred from expression and prior motif databases. They are not
proof of direct binding, causal regulation, or a Sample-level Condition effect.

Official GRN importance and final target weights are nonnegative. A zero value remains a well-typed backend output
(for example from an expert top-N selection), so it is retained and counted in warnings; only negative or non-finite
weights violate this official-flow contract.

## References to emit

- Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature
  Protocols*. 2020;15:2247-2276. https://doi.org/10.1038/s41596-020-0336-2
- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463
- Moerman T, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks.
  *Bioinformatics*. 2019;35:2159-2161. https://doi.org/10.1093/bioinformatics/bty916
- Virshup I, et al. anndata: Annotated data. *Journal of Open Source Software*. 2024;9:4371.
  https://doi.org/10.21105/joss.04371

## Report and generated-code implications

The report must identify the external workflow/environment, exact cell/gene/regulon/edge/motif counts, expression
state, resource bundle and hashes, CLI parameters/seeds or omitted-default warnings, imported score range, any excluded/merged motif records,
manifest/member hashes, limitations, references, and dynamic importer/runtime versions. Equivalent code must use the
same safe parser and hash/identity postconditions, perform no network access or code evaluation, and return the same
AnnData/artifact/strict summary.
