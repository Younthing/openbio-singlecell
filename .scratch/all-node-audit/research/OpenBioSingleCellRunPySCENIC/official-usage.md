# OpenBioSingleCellRunPySCENIC — official usage research

## Audited baseline

pySCENIC's latest official release is 0.12.1 (21 November 2022). Its official installation documentation uses a dedicated
Python 3.10 environment or the versioned `aertslab/pyscenic:0.12.1` container. The repository requires Python 3.12
and currently resolves NumPy 2.4.4. In the configured ComfyUI environment, importing the pySCENIC CLI fails because
0.12.1 still references removed aliases `numpy.object`; `pyscenic.rss` also uses removed `numpy.float`. The in-process
node cannot be made correct by changing one call.

The current node invokes all three CLI stages in a temporary directory, deletes every reusable output except a
pre-cisTarget adjacency object, and does not emit `summary` or `code`. It defaults to a log-normalized layer although
the protocol starts from post-QC full-gene counts. It also mixes external workflow orchestration, resource handling,
result import, AnnData mutation, and network construction in one shallow interface.

## Official pySCENIC 0.12.1 workflow

Official references:

- https://pyscenic.readthedocs.io/en/stable/installation.html
- https://pyscenic.readthedocs.io/en/stable/tutorial.html
- https://pyscenic.readthedocs.io/en/stable/faq.html
- https://github.com/aertslab/pySCENIC/tree/0.12.1
- https://doi.org/10.1038/s41596-020-0336-2

The documented container workflow is three explicit stages:

```bash
docker run --rm -v /data:/data aertslab/pyscenic:0.12.1 \
  pyscenic grn --method grnboost2 --num_workers 6 --seed 777 \
  -o /data/expression.adjacencies.tsv \
  /data/expression.csv /data/allTFs_hg38.txt

docker run --rm -v /data:/data aertslab/pyscenic:0.12.1 \
  pyscenic ctx /data/expression.adjacencies.tsv \
  /data/database1.rankings.feather /data/database2.rankings.feather \
  --annotations_fname /data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /data/expression.csv \
  --mode custom_multiprocessing --num_workers 6 \
  --output /data/regulons.csv

docker run --rm -v /data:/data aertslab/pyscenic:0.12.1 \
  pyscenic aucell --num_workers 6 --seed 777 --auc_threshold 0.05 \
  -o /data/aucell.csv /data/expression.csv /data/regulons.csv
```

The SCENIC protocol conventionally begins with a post-QC count matrix and full gene space. The 0.12.1 CSV loader also
accepts other finite numeric expression matrices, so an expert transformed state is backend-valid when named and
disclosed rather than silently relabeled as counts. The workflow first infers co-expression modules
with GRNBoost2 or GENIE3, prunes indirect targets using species/genome-build-matched cisTarget ranking databases and
motif annotations, then scores the motif-pruned regulons per cell with AUCell. Restricting the input to HVGs changes
the regulatory search universe and is not the protocol default.

Since pySCENIC 0.12.0, `ctxcore>=0.2` Feather-v2 ranking databases are the supported ranking format. The TF list,
ranking databases, motif annotations, organism, genome build, gene-identifier namespace, and database
release are a coherent resource bundle. Mixing hg19/hg38/mm9/mm10 resources or identifier namespaces invalidates
the motif-pruning stage. Image tags and file paths alone are insufficient provenance; content SHA-256 hashes are
required for consumed inputs/outputs, while an immutable container digest is strongly preferred and may be absent
for a fully version-recorded non-container environment.

## External-run artifact requirements

An auditable external run must retain at least:

- the exact input cell and gene identifiers plus a canonical expression-state fingerprint;
- the exported expression input, adjacency TSV/CSV, motif-enrichment/regulon CSV, and AUCell CSV matrix;
- pySCENIC/Python/dependency versions and, when used, the container image plus immutable digest;
- complete CLI argument arrays; explicit integer seeds for `grn` and `aucell`, worker counts, dropout-mask and
  AUCell-threshold settings are strongly preferred, while official omitted defaults must be warned;
- TF list, every ranking database, and motif annotation path, release, organism/build, license, and SHA-256;
- output file SHA-256 values and a versioned machine-readable run manifest.

The computation is stochastic at GRN inference and AUCell's per-cell gene-ranking tie break. Both CLI stages expose
`--seed`; an explicit integer, including zero, is passed through as a real seed. This differs from the separate
pySCENIC 0.12.1 `binarize` helper, whose truthiness check makes zero behave as unspecified. An explicit seed is the
reproducible choice; an omitted seed is still backend-valid but must be recorded as omitted and warned. Parallel tree
fitting may still require environment/thread provenance and empirical reproducibility checks.

## References to emit

- Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature
  Protocols*. 2020;15:2247-2276. https://doi.org/10.1038/s41596-020-0336-2
- Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*.
  2017;14:1083-1086. https://doi.org/10.1038/nmeth.4463
- Moerman T, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks.
  *Bioinformatics*. 2019;35:2159-2161. https://doi.org/10.1093/bioinformatics/bty916
- Huynh-Thu VA, et al. Inferring regulatory networks from expression data using tree-based methods. *PLoS ONE*.
  2010;5:e12776. https://doi.org/10.1371/journal.pone.0012776

## Consequence for this node

The official, reproducible use is an external container/workflow followed by a strict import node. Automatically
starting Docker/Podman/Singularity from a scientific ComfyUI node would add host privileges, filesystem mounts,
network access, long-lived processes, and platform-specific failure modes to an already incompatible runtime. This
node should not claim that it can execute pySCENIC inside the plugin environment.

## Public compatibility status

The retained ID is a non-executing compatibility surface for old graphs, not an official pySCENIC runner in this
runtime. Its schema must be explicitly deprecated and direct users to the pinned external workflow plus strict
import seam. Since no computation occurs, the shim has no generated method citation, runtime software record,
`summary`, or equivalent `code`; those belong to the external manifest and active importer.
