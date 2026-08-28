# OpenBioSingleCellLeiden — design review

## Decision

**Keep and enhance the node as the atomic one-resolution clustering module. Do not merge it with neighbor construction or the resolution sweep, and do not delete it.**

The atomic result is one declared partition of one already-built graph. Neighbor construction is a different scientific operation whose representation, metric, and neighborhood size must be reviewable before clustering. The sweep is a diagnostic decision-support operation across several resolutions. Merging all three would enlarge the interface, obscure ordering constraints, and prevent a reviewer from distinguishing graph selection from partition selection.

## Current interface problems

The current implementation is a shallow adapter around one Scanpy call:

- it does not validate or disclose the actual connectivity matrix selected by `neighbors_key`;
- it silently hides flavor, directedness, weights, objective, and iteration count;
- it fixes two iterations without allowing an analyst to request convergence;
- it does not assess stochastic stability;
- it silently overwrites existing `obs`/`uns` keys;
- it returns no cluster-size, graph, modularity, or robustness diagnostics;
- it has no required `summary` and `code` outputs;
- its physical implementation is separated from the resolution-sweep node, encouraging duplicated behavior.

## Proposed interface

Visible inputs:

- `adata`;
- `resolution=1.0`, finite and non-negative; zero is executable and warned.

Advanced inputs:

- `neighbors_key="neighbors"`;
- `key_added="leiden"`;
- `n_iterations=2`, accepting every integer (negative means until stable, zero means no optimization);
- `stability_repeats=5`, any positive integer including the base start; one start yields null stability statistics and large workloads are warned;
- `random_seed=0`;
- `overwrite_existing=False`.

Hidden, fixed, and fully disclosed policy:

- `flavor="igraph"`;
- `directed=False`;
- `use_weights=True`;
- `objective_function="modularity"`;
- no custom partition, initial membership, beta, node weights, or maximum community size.

Those fixed values are not mere implementation details: they define the partition objective. Hiding them is appropriate only if summary/code make them visible to report readers. Exposing all igraph/leidenalg arguments on one node would create a shallow interface whose modes do not share one scientific meaning.

Outputs:

- annotated `adata`;
- strict JSON `summary`;
- equivalent Python source `code`.

## Deep module and seam placement

The node becomes a deep module by hiding graph-key resolution, sparse validation, fixed-backend adaptation, deterministic repeat-seed generation, ARI calculation, aligned storage, and report construction behind one small clustering interface. The AnnData/graph seam is the external seam: callers supply an aligned graph and receive one partition plus auditable evidence.

Create an internal clustering implementation shared with Resolution Sweep. It should accept a validated graph contract and resolved run settings, then return memberships, modularity, cluster sizes, and stability diagnostics. This is a real internal seam because the same adapter has two consumers. Do not expose Scanpy private utilities or igraph objects on node sockets.

The interface is also the test surface: production tests can construct small valid/malformed AnnData graphs and exercise the same graph resolver and clustering path used by callers. Tests should not reach into Scanpy internal graph objects.

## Cohesion and coupling

Every public input contributes to the single partition or its reproducibility evidence. Stability repeats are cohesive because they test stochastic robustness of that exact partitioning problem; they do not create additional public annotations. Upstream neighbor parameters remain coupled by scientific necessity but are inherited read-only rather than duplicated as Leiden inputs.

The deletion test supports retention: deleting this module would force graph validation, backend policy, collision handling, stability calculation, reporting, and canonical storage into every workflow that needs one chosen resolution. Deleting a generic flavor/partition argument, by contrast, removes complexity without moving required behavior elsewhere, so those controls stay out of the first refactored interface.

## Summary and code contract

`summary` contains methods, report-ready results, key results, parameters, warnings, limitations, references, and software versions. Key results include graph/component diagnostics, complete cluster sizes/fractions, modularity, and mean/min/max pairwise ARI. The prose explicitly calls the result a graph partition, not a cell-type truth or Condition result.

`code` reproduces the scientific output and diagnostics with public Scanpy/scikit-learn functions and every hidden policy value explicit. Plugin timing/history may be absent; membership, categorical labels, Scanpy `uns` result, stability diagnostics, and failures must be equivalent.

## Failure and verification plan

- Reject backed/empty/single-observation inputs and duplicate observation identifiers.
- Reject missing neighbor metadata, missing connectivity key, wrong shape, dense/object/non-numeric ambiguity, NaN/Infinity, negative weights, self-loops, asymmetry, and edge-free graphs.
- Report but do not automatically reject disconnected components or isolated cells.
- Reject non-finite/negative resolution, non-integer iterations, non-positive repeat count, empty result key, and collisions without overwrite. After resolving the graph, reject `key_added == graph.neighbors_key` unconditionally—even when overwrite is requested—so a clustering result cannot delete the graph's `uns` pointer bundle.
- Assert input immutability and preservation of X, layers, Raw snapshot, `obsm`, and graph.
- Assert categorical complete membership, stable natural labels, finite modularity, exact cluster-size totals, and deterministic output for a fixed seed/software stack.
- Verify repeated-start ARI against direct scikit-learn calculation and distinguish it from biological validity.
- Compare generated-code membership, diagnostics, graph-key reservation, and malformed-input failures with runtime; verify rejected and successful runs leave the named graph resolvable for downstream embeddings.
- Require strict JSON and H5AD round-trip of the annotated output.

The open-tool implementation keeps type/axis/finite/undirected graph failures hard, but turns zero resolution, zero iterations, unusual repeat counts, singleton/dominant communities, and weak repeat-start agreement into disclosed expert diagnostics. Runtime and generated source cross the same graph seam and accept graph-only zero-variable AnnData.

## 2026-08-29 repair record

- New category: `openbio/single-cell/clustering`.
- Classification rationale: the atomic analysis partitions one existing named neighbor graph; it does not construct or reduce a feature representation, so `dimension-reduction` misstates the operation.
- Merge/delete decision: retain the one-resolution node. Do not merge it with neighbor construction or the multi-resolution sweep and do not delete it; those nodes answer distinct graph-construction and resolution-selection questions.
