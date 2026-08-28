# Batch09 regulatory nodes — exact legacy/current migration and release matrix

Date: 2026-08-28. This is the input contract for the central workflow-migration owner. The pySCENIC owner did not
edit `web/workflow_migrations.mjs`. Every legacy detection below must be completed as a graph-wide preflight before
any node, link, array entry, subgraph, exposure, or global ID is mutated. A blocked workflow remains byte-for-byte
unchanged and allocates no IDs.

## Global migration rule

- Traverse root nodes, nested arrays, object-valued node maps, subgraphs, and exposed inputs/outputs before applying
  any edit. Match both node type and the exact legacy input/output signature; an unrecognized hybrid is blocked.
- Resolve links by stable node ID and port identity, not array position alone. Preflight every incoming/outgoing link
  and every exposed port. Never drop an edge or exposure to make a migration pass.
- Treat nested array/object node containers and nested/object subgraph-definition containers as part of the global
  inventory. If the central mutator cannot preserve their global-ID semantics, detect the Batch09 node and fail
  closed before any allocator or unrelated migration runs; never silently skip the container.
- Exact interface evidence includes runtime scalar types, converted-widget `widget.name`, link types, both origin and
  target endpoint types, subgraph-boundary types, and optional direct/proxy exposure type metadata. Names or equal
  widget counts alone do not prove a current schema.
- Current signatures below are idempotent no-ops. Mixed legacy/current graphs are not partially upgraded.
- All ten legacy signatures require an explicit scientific choice that cannot be reconstructed. The safe automatic
  action is therefore a diagnostic block, not a best-effort rewrite. Future interactive migration may use the
  conditional paths below only after the user supplies the missing values.

## Final adversarial container review

Independent release review found ten discovery/closure counterexample families that the earlier heuristic inventory
missed:

1. an object-valued `definitions.subgraphs` map with its own reserved `id` key was misclassified as a graph and hid
   an inner legacy Batch09 node;
2. an object-valued node map with its own `type`/`id` keys was misclassified as one node and hid an inner blocker;
3. a type-only partial object such as `{type: "OpenBioSingleCellRunPySCENIC"}` escaped node classification entirely;
4. malformed converted-widget metadata (`widget={}`, null name, array, or string) was accepted;
5. non-array `graph.widgets` was silently interpreted as no exposure metadata;
6. non-array `graph.links` reached the global allocator and raised a native `TypeError` instead of a stable
   migration diagnostic;
7. malformed or ambiguous serialized identity/structure still let an unrelated migration commit: non-scalar or
   non-safe-integer node/subgraph IDs, one Batch09 object owned by more than one graph position, duplicate Batch09
   IDs across graph definitions, non-node entries in a flat `nodes` array, incomplete link records, incomplete
   `proxyWidgets` tuples, a type-only subgraph host that hid an incompatible proxy exposure, and a direct node whose
   ID collided with the reserved `-10`/`-20` subgraph boundary sentinels;
8. canonical graph discovery stopped after seeing one graph-like object and did not inspect its nonstandard nested
   values. A hidden graph could therefore conceal an incompatible subgraph proxy, a larger node ID, or a complete
   link whose ID was outside the allocator; an unrelated cNMF migration then allocated that same hidden node/link ID
   and committed a real global-ID collision;
9. dangling node input/output link references were not represented in `graph.links`, so the allocator reused their
   IDs and accidentally bound old unrelated ports to a newly inserted cNMF link. All direct-node and subgraph-boundary
   port references therefore require bidirectional closure against exactly one complete graph-link record.
10. unsafe serialized high-water counters were silently ignored while later synchronization preserved them, and
    conservative capacity estimates rejected valid same-node migrations at the safe-integer boundary. Present
    `state.lastNodeId`/`lastLinkId` and legacy `last_node_id`/`last_link_id` values must be nonnegative safe integers
    (numbers or canonical unsigned decimal strings); a present non-null `state` must be a non-array object so synchronized counters
    survive JSON serialization. Capacity reservation must count only exact insertions proven by completed plans.

The final design uses one deep Batch09 **presence scan** over the original workflow object, independent of graph/node
shape heuristics. Every object whose `type` is a Batch09 node ID must occur as a direct member of a proven flat
`graph.nodes` array; otherwise the entire workflow is rejected before any allocator or mutation. When any occurrence
exists, every discovered graph must expose flat array `nodes`, array-or-absent `links`, and array-or-absent `widgets`.
An independent cycle-safe deep graph-candidate scan inventories every object with an own `nodes` or `links` field,
matching canonical graph discovery. The only excluded `links` owners are object identities proven to be output-port
records of direct nodes in canonical graphs; their arrays contain scalar link references rather than graph link
records. This catches marker-free hidden graphs with array, object-map, or otherwise malformed `links` containers.
Each candidate must already belong to the canonical root/`definitions.subgraphs` graph inventory; nonstandard
hidden graph containers are rejected before global ID allocation rather than guessed into scope.
Every direct `widgets` entry must itself be a non-array object with a scalar string/number node ID and string name;
malformed entries cannot be treated as evidence that no exposure exists.
Subgraph discovery no longer treats `id` alone as graph evidence, so reserved keys in object maps remain containers.
Converted-widget metadata, when present, must be a non-array object with a string `name` exactly matching the input.
The trusted inventory requires each direct node to be a non-array object with a nonempty string or finite-safe-integer ID
that is not the reserved `-10`/`-20` boundary sentinel. A normal registered node type is a nonempty string; a numeric
type is accepted only when it exactly identifies one canonical subgraph definition, matching ComfyUI's serialized
subgraph-host form. Each Batch09 object must have exactly one direct graph position, and Batch09 IDs must be globally
unique across discovered graphs. A subgraph containing Batch09 must likewise have a stable definition ID.
Every serialized link must be a complete array/object record with stable string/finite-safe-integer
link/origin/target IDs, nonnegative integer slots, and a string type. Numeric-looking string IDs outside JavaScript's
safe-integer range are rejected instead of being rounded into allocator counters. Every node input `link`, node
output `links`, subgraph input/output `linkIds`, and graph-link endpoint must resolve bidirectionally and uniquely.
Every host `proxyWidgets` entry must be a tuple with a stable node ID and nonempty string widget name, and host
discovery uses the same trusted flat-node inventory rather than a permissive node heuristic.
These are serialization/atomicity boundaries only; they do not inspect analysis history, bind Raw to current X, or
constrain expert scientific parameter choices.

Allocator counter metadata is validated independently of graph IDs. Missing counters mean zero; a present counter
may be a nonnegative safe-integer number or canonical unsigned decimal string, but negative, fractional, nonnumeric, or unsafe values
fail before mutation. Missing/null `state` is treated as absent; every other present state must be a non-array object,
and allocation plus synchronization use that same shape rule. Legacy `last_node_id` and `last_link_id` are accepted
and synchronized independently; an absent sibling is not synthesized, preserving the graph's link-serialization mode
while preventing a present counter from becoming stale. Reservation is exact: only non-reused legacy DPT plans, CNV Structure expansion, active staged
cNMF, newly inserted Cassiopeia QC producers, Cassiopeia reconstruction links, and proven marker-universe links
consume IDs. Same-node schema upgrades and current no-ops consume no capacity.

## Exact node matrix

| Node | Legacy inputs → outputs | Current inputs → outputs | Safe decision and reason |
|---|---|---|---|
| `OpenBioSingleCellCollecTRIULM` | `adata, organism(human|mouse), split_complexes=false, source(X|raw|layer), activity_key, pvalue_key` → `adata` | `adata, resource(local_network{network_csv,resource_metadata_json}|official_collectri{affiliation_license,allow_network_access}), source(X|layer), complex_policy(retain|remove), min_targets, batch_size, max_output_rows, max_working_memory_gib, overwrite_existing` → `adata, activities, summary, code` | **Block.** `split_complexes` cannot be reinterpreted as either resource identity or `complex_policy`; `true` means split complex members and has no current equivalent, while even `false` does not establish a licensed/hash-bound resource. Mouse official mode is deliberately absent. Raw input, arbitrary keys, and network provenance also cannot be reconstructed. Interactive replacement must choose a resource/license, expression state, and complex policy. |
| `OpenBioSingleCellRankTFActivities` | `adata, groupby, reference, method(wilcoxon|t-test...), top_n, max_pvalue(raw), activity_key` → `table` | `adata, activities, annotation_key, annotation_status, reference, method(t-test_overestim_var|t-test|wilcoxon), report_p_adjusted, max_output_rows` → `table, summary, code` | **Block.** `groupby→annotation_key` and method/reference text are mechanically translatable, but no typed activity artifact exists. `max_pvalue` was a raw-p filter and cannot become the adjusted-p reporting flag; `top_n` truncation was removed and annotation status is unknown. It can only be replaced in the same reviewed operation as a new CollecTRI producer and explicit resource choice. |
| `OpenBioSingleCellRunPySCENIC` | `adata, source, use_highly_variable, highly_variable_key, tf_list_file, ranking_database_files, motif_annotations_file, grn_method, mask_dropouts, auc_threshold, num_workers, random_seed, activity_key` → `adata, network` | same retained legacy widgets → no outputs; `is_output_node=true` retired shim | **Block and instruct external rerun.** The old implementation deleted its temporary adjacency/regulon/AUCell files and captured neither a complete manifest nor immutable container/resource provenance. Neither output can be recreated. Even an isolated instance should require acknowledgement before becoming the nonexecuting shim; any output link/exposure is an unconditional block. |
| `OpenBioSingleCellImportPySCENICResults` | `adata, final_loom_file, regulons_csv, activity_key` → `adata` | `adata, run_manifest_json, overwrite, max_file_bytes, max_adjacency_edges, max_regulon_edges, max_dense_bytes` → `adata, scenic_result, summary, code` | **Block.** A loom path plus regulon CSV cannot supply the expression export, adjacency, licensed resources, commands, exact versions, expression-state declaration, or file hashes required by the manifest. A container digest is recommended but optional and does not remove those missing scientific/provenance inputs. Do not invent metadata from filenames or hash only the two surviving files. |
| `OpenBioSingleCellSCENICRegulonSpecificity` | `adata, groupby, activity_key` → `table` | `adata, scenic_result, annotation_key, annotation_status, max_output_rows` → `table, summary, code` | **Block.** `groupby→annotation_key` is safe, but the activity key is untyped hidden AnnData state and cannot establish imported artifact provenance or complete regulon membership. A new typed Import connection and annotation-status choice are mandatory. |
| `OpenBioSingleCellSCENICActivityBinarization` | `adata, activity_key, binary_key, threshold_key` → `adata` | `scenic_result, random_seed(1..2^31-1), threshold_overrides?(typed table), max_dense_bytes` → `binary, thresholds, summary, code` | **Block.** The legacy call had no recorded seed, mutated arbitrary AnnData keys, and exposed neither threshold family nor a typed source. Seed zero/unspecified cannot be relabeled deterministic. Downstream AnnData links are not semantically replaceable by a binary artifact link without an explicit adapter. |
| `OpenBioSingleCellSCENICTFModules` | `adata, network(OPENBIO_SCENIC_NETWORK), transcription_factor, source` → `table` | `scenic_result(OPENBIO_SCENIC_RESULT), transcription_factor, max_output_rows` → `table, summary, code` | **Block unconditionally.** Legacy output recomputed pre-cisTarget co-expression modules from raw adjacency and expression. Current output is final motif-pruned regulon membership. They answer different questions; changing the wire type would falsely relabel evidence even if a new Import artifact were available. |
| `OpenBioSingleCellLianaCommunication` | `adata, groupby, method, resource_name, source(X\|raw\|layer), result_key` → `adata` | `adata, sample_key, condition_key, identity_key, annotation_status, organism, method, resource(bundled_human\|local_resource with metadata), source(X\|layer), expression_proportion, min_cells_per_identity_sample, permutations, random_seed, jobs, max_output_rows, max_working_memory_gib` → `result, summary, code` | **Block.** The pooled legacy run lacks biological Sample/Condition roles, annotation status, reviewed resource metadata, deterministic execution details, and a typed direct result. Preserve neither its hidden `uns` table nor its Raw branch as evidence. A reviewed replacement must choose the explicit current source/resource/roles and rerun by Sample. Runtime source selection is an expert choice and does not require history or a Raw/current binding. |
| `OpenBioSingleCellLianaResults` | `adata, result_key` → generic `table` | registered fail-before-read retirement shim; no active result | **Block.** A hidden `uns` key cannot recover LIANA method, resource, Sample/Condition roles, score directions, or complete-family provenance. Remove it only in the same reviewed operation that reruns Communication and rewires typed consumers; arbitrary table consumers are incompatible. |
| `OpenBioSingleCellLianaDotPlot` | `adata, method, source_labels, target_labels, significance_threshold, top_n, result_key, figure_width, figure_height` → `plot` | `result(OPENBIO_LIANA_RESULT), source_labels, target_labels, selection(rank_aggregate.max_specificity_rank\|cellphonedb.max_cellphone_pvalue), top_n, figure_width, figure_height, max_plot_rows, max_image_pixels` → `plot, summary, code` | **Block.** The old generic threshold and hidden table do not prove a method-specific typed result. Rewire only after a fresh reviewed Communication run; the threshold branch may be selected from that proven producer method, but standalone/mismatched/pooled plots remain untranslatable. |

## Conditional interactive replacement order

1. Obtain or rerun a pySCENIC 0.12.1 external bundle, choose a licensed CollecTRI resource, or declare the LIANA
   Sample/Condition/identity roles and reviewed resource metadata; do not reuse hidden legacy AnnData keys as
   provenance.
2. Replace the producer first and connect its typed artifact output.
3. Replace each consumer with explicit annotation status, row/memory guard, and (for binarization) nonzero seed.
4. Review all old table/AnnData consumers and exposures because primary output meanings changed. Only after all links
   validate should one atomic graph transaction commit; preserve existing node IDs where the type is unchanged and
   allocate deterministic new IDs only for explicitly requested adapters.

## External pySCENIC dependency recommendation

Do not add pySCENIC, ctxcore, or loompy to the plugin's Python >=3.12 core. Prefer an immutable external Python 3.10
container; a fully version-recorded non-container environment remains importable with a reproducibility warning. A reproducible legacy-compatible test environment should pin at least
`pyscenic==0.12.1`, `ctxcore==0.2.0`, `arboreto==0.1.6`, `loompy==3.0.8`, `numpy>=1.21,<1.24`, and
`pandas>=1.3.5,<2`. The published pySCENIC metadata leaves NumPy unbounded, but tagged 0.12.1 code uses removed
`np.object`/`np.float`; its dip implementation also uses removed `np.msort`. Therefore `numpy<2` is not a sufficient
release constraint. The plugin importer itself has no pySCENIC/ctxcore/loompy runtime dependency and accepts strict
CSV plus Feather-v2 framing only.

The CollecTRI optional dependency must be `decoupler==2.2.0`, not `>=2.1.4,<3`, because the audited implementation
guards that exact public `op.collectri`/`mt.ulm` version, signatures, and output keys. The current release manifest
already reflects the exact pin; retain it.

## Release artifact types

- `OPENBIO_SCENIC_RESULT` → internal artifact `openbio-singlecell/scenic-result`, schema 1, sole producer
  `OpenBioSingleCellImportPySCENICResults`; consumers RSS, binarization, and final-regulon membership.
- `OPENBIO_SCENIC_BINARY` → internal artifact `openbio-singlecell/scenic-binary`, schema 1, sole producer
  `OpenBioSingleCellSCENICActivityBinarization`.
- `OPENBIO_SCENIC_NETWORK` has no current producer or consumer and is removed from the Python/public release
  contract. The central migration preflight must recognize that literal legacy wire name directly; recognition does
  not require retaining a dead connectable type.
