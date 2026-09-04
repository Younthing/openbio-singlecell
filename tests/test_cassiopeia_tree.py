from __future__ import annotations

import copy
import hashlib
import itertools
import json
import math
import random
import uuid
from collections import Counter
from dataclasses import replace
from types import SimpleNamespace

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import openbio_singlecell.cassiopeia_tree as lineage
from openbio_singlecell.artifact_codecs import read_plot, write_anndata, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.cassiopeia_codec import write_characters, write_tree
from openbio_singlecell.cassiopeia_tree import CassiopeiaCharacters, CassiopeiaTree
from openbio_singlecell.node_types import (
    AnnDataType,
    CassiopeiaCharactersType,
    CassiopeiaTreeType,
    PlotResultType,
    SummaryResultType,
    TableResultType,
)
from openbio_singlecell.nodes_lineage import (
    LINEAGE_NODE_CLASSES,
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaLineageQC,
    OpenBioSingleCellCassiopeiaPlasticity,
    OpenBioSingleCellReconstructCassiopeiaTree,
)
from openbio_singlecell.operations_lineage import (
    cassiopeia_expansion_plot,
    cassiopeia_expansion_plot_owned,
    cassiopeia_lineage_qc_plot,
    cassiopeia_plasticity_plot,
    cassiopeia_plasticity_plot_owned,
    cassiopeia_tree_plot,
    cassiopeia_tree_plot_owned,
    expansion_owned,
    lineage_qc_owned,
    lineage_qc_plot_owned,
    plasticity_owned,
)
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def assign_missing_average(*args, **kwargs):
    del args, kwargs
    return None


class FakeCassiopeiaTree:
    copy_count = 0

    def __init__(
        self,
        character_matrix=None,
        missing_state_indicator=-1,
        cell_meta=None,
        character_meta=None,
        priors=None,
        tree=None,
        dissimilarity_map=None,
        parameters=None,
        root_sample_name=None,
    ):
        del character_meta, dissimilarity_map, parameters, root_sample_name
        self.character_matrix = None if character_matrix is None else character_matrix.copy(deep=True)
        self.missing_state_indicator = missing_state_indicator
        self.cell_meta = None if cell_meta is None else cell_meta.copy(deep=True)
        self.priors = copy.deepcopy(priors) if priors is not None else {}
        self._graph = None if tree is None else tree.copy()

    def copy(self):
        type(self).copy_count += 1
        return copy.deepcopy(self)

    @property
    def root(self):
        roots = [node for node in self._graph if self._graph.in_degree(node) == 0]
        return roots[0]

    @property
    def leaves(self):
        return [node for node in self._graph if self._graph.out_degree(node) == 0]

    def get_tree_topology(self):
        return self._graph.copy()

    def set_attribute(self, node, attribute_name, value):
        self._graph.nodes[node][attribute_name] = value

    def get_attribute(self, node, attribute_name):
        return self._graph.nodes[node][attribute_name]


def _fake_priors(allele_table, grouping_variables=None, cut_sites=None):
    if grouping_variables is None:
        grouping_variables = ["intBC"]
    if cut_sites is None:
        raise AssertionError("Audited caller must provide exact cut sites.")
    groups = list(allele_table.groupby(grouping_variables, observed=True, sort=True, dropna=False))
    counts: Counter[str] = Counter()
    for _, group in groups:
        alleles = {
            str(value)
            for cut_site in cut_sites
            for value in group[cut_site].tolist()
            if not pd.isna(value) and "none" not in str(value).lower()
        }
        counts.update(alleles)
    if not counts:
        return pd.DataFrame(columns=["count", "freq"])
    return pd.DataFrame({"count": counts, "freq": {allele: count / len(groups) for allele, count in counts.items()}})


def _fake_converter(
    alleletable,
    ignore_intbcs=None,
    allele_rep_thresh=1.0,
    missing_data_allele=None,
    missing_data_state=-1,
    mutation_priors=None,
    cut_sites=None,
    collapse_duplicates=True,
):
    del ignore_intbcs, missing_data_allele, collapse_duplicates
    if cut_sites is None:
        raise AssertionError("Audited caller must provide exact cut sites.")
    cells = list(dict.fromkeys(alleletable["cellBC"].tolist()))
    keys: list[tuple[str, str]] = []
    observations: dict[tuple[str, str], list[object]] = {}
    by_cell: dict[str, dict[tuple[str, str], object]] = {cell: {} for cell in cells}
    for row in alleletable.itertuples(index=False):
        record = row._asdict()
        cell = str(record["cellBC"])
        intbc = str(record["intBC"])
        for cut_site in cut_sites:
            key = (intbc, cut_site)
            if key not in observations:
                keys.append(key)
                observations[key] = []
            value = record[cut_site]
            by_cell[cell][key] = value
    retained = []
    for key in keys:
        values = [by_cell[cell].get(key, np.nan) for cell in cells]
        normalized = ["<missing>" if pd.isna(value) else str(value) for value in values]
        counts = Counter(normalized)
        if not any(count / len(normalized) > allele_rep_thresh for count in counts.values()):
            retained.append(key)
    data: list[list[int]] = []
    state_maps: dict[int, dict[int, str]] = {position: {} for position in range(len(retained))}
    priors: dict[int, dict[int, float]] = {position: {} for position in range(len(retained))}
    allele_states: dict[tuple[str, str], dict[str, int]] = {key: {} for key in retained}
    for cell in cells:
        row_values = []
        for position, key in enumerate(retained):
            value = by_cell[cell].get(key, np.nan)
            if pd.isna(value):
                row_values.append(missing_data_state)
            elif str(value) == "NONE":
                row_values.append(0)
            else:
                allele = str(value)
                state = allele_states[key].setdefault(allele, len(allele_states[key]) + 1)
                state_maps[position][state] = allele
                priors[position][state] = float(mutation_priors.loc[allele, "freq"])
                row_values.append(state)
        data.append(row_values)
    frame = pd.DataFrame(data, index=cells, columns=[f"r{index + 1}" for index in range(len(retained))])
    return frame, priors, state_maps


class FakeVanillaGreedySolver:
    mode = "balanced"
    solve_count = 0
    last_classifier = None
    last_transformation = None
    last_collapse = None

    def __init__(
        self,
        missing_data_classifier=assign_missing_average,
        prior_transformation="negative_log",
    ):
        type(self).last_classifier = missing_data_classifier
        type(self).last_transformation = prior_transformation

    def solve(self, cassiopeia_tree, layer=None, collapse_mutationless_edges=False, logfile="stdout.log"):
        del layer, logfile
        np.random.random()
        random.random()
        type(self).solve_count += 1
        type(self).last_collapse = collapse_mutationless_edges
        cells = list(cassiopeia_tree.character_matrix.index)
        root = f"volatile-root-{type(self).solve_count}"
        graph = nx.DiGraph()
        if type(self).mode == "star" or len(cells) < 4:
            graph.add_edges_from((root, cell, {"length": 1.0}) for cell in cells)
        else:
            left = f"volatile-left-{type(self).solve_count}"
            right = f"volatile-right-{type(self).solve_count}"
            graph.add_edge(root, left, length=1.0)
            graph.add_edge(root, right, length=1.0)
            midpoint = len(cells) // 2
            graph.add_edges_from((left, cell, {"length": 1.0}) for cell in cells[:midpoint])
            graph.add_edges_from((right, cell, {"length": 1.0}) for cell in cells[midpoint:])
        cassiopeia_tree._graph = graph


def _leaf_counts(graph, node):
    children = list(graph.successors(node))
    if not children:
        return 1
    return sum(_leaf_counts(graph, child) for child in children)


def _fake_expansion(tree, min_clade_size=10, min_depth=1, copy=False):
    result = tree.copy() if copy else tree
    graph = result._graph
    root = result.root
    depths = {root: 0}
    for parent, child in nx.dfs_edges(graph, root):
        depths[child] = depths[parent] + 1
    for node in graph:
        result.set_attribute(node, "expansion_pvalue", 1.0)
    for parent in nx.dfs_preorder_nodes(graph, root):
        children = list(graph.successors(parent))
        n = _leaf_counts(graph, parent)
        k = len(children)
        for child in children:
            b = _leaf_counts(graph, child)
            if b >= min_clade_size and depths[child] >= min_depth:
                result.set_attribute(child, "expansion_pvalue", math.comb(n - b, k - 1) / math.comb(n - 1, k - 1))
    return result if copy else None


def _fake_small_parsimony(
    cassiopeia_tree,
    meta_item,
    root=None,
    infer_ancestral_states=True,
    label_key="label",
):
    del infer_ancestral_states, label_key
    graph = cassiopeia_tree._graph
    root = cassiopeia_tree.root if root is None else root
    reachable = set(nx.descendants(graph, root)) | {root}
    leaves = [node for node in reachable if graph.out_degree(node) == 0]
    states = sorted(str(cassiopeia_tree.cell_meta.loc[leaf, meta_item]) for leaf in leaves)
    states = list(dict.fromkeys(states))
    internals = [node for node in reachable if node not in leaves]
    best = math.inf
    fixed = {leaf: str(cassiopeia_tree.cell_meta.loc[leaf, meta_item]) for leaf in leaves}
    for assignments in itertools.product(states, repeat=len(internals)):
        labels = {**fixed, **dict(zip(internals, assignments, strict=True))}
        score = sum(labels[parent] != labels[child] for parent, child in graph.edges if parent in reachable)
        best = min(best, score)
    return int(best)


def make_fake_cassiopeia():
    return SimpleNamespace(
        __version__="2.1.3",
        pp=SimpleNamespace(
            compute_empirical_indel_priors=_fake_priors,
            convert_alleletable_to_character_matrix=_fake_converter,
        ),
        data=SimpleNamespace(CassiopeiaTree=FakeCassiopeiaTree),
        solver=SimpleNamespace(VanillaGreedySolver=FakeVanillaGreedySolver),
        tl=SimpleNamespace(
            compute_expansion_pvalues=_fake_expansion,
            score_small_parsimony=_fake_small_parsimony,
        ),
    )


@pytest.fixture
def fake_backend(monkeypatch):
    fake = make_fake_cassiopeia()
    FakeVanillaGreedySolver.mode = "balanced"
    FakeVanillaGreedySolver.solve_count = 0
    FakeCassiopeiaTree.copy_count = 0
    monkeypatch.setattr(lineage, "_require_cassiopeia", lambda: fake)
    return fake


def _rows():
    return [
        {"Tumor": "T1", "cellBC": "c1", "intBC": "i1", "r1": "A", "r2": "NONE"},
        {"Tumor": "T1", "cellBC": "c2", "intBC": "i1", "r1": "A", "r2": "B"},
        {"Tumor": "T1", "cellBC": "c3", "intBC": "i1", "r1": "C", "r2": "B"},
        {"Tumor": "T1", "cellBC": "c4", "intBC": "i1", "r1": "C", "r2": "NONE"},
    ]


def _prepare(tmp_path, rows=None, module=lineage):
    path = tmp_path / "alleles.tsv"
    pd.DataFrame(_rows() if rows is None else rows).to_csv(path, sep="\t", index=False)
    return module.prepare_cassiopeia_characters(
        str(path),
        first_column_as_index=False,
        cut_site_columns="r1,r2",
        prior_grouping_columns="Tumor,intBC",
        allele_representation_threshold=0.98,
        minimum_cells=2,
        maximum_missing_fraction=1.0,
        maximum_uncut_fraction=1.0,
        minimum_unique_fraction=0.0,
        minimum_informative_character_fraction=0.0,
        openbio_version="test",
    )


def _qc_owned(tmp_path):
    path = tmp_path / "alleles.tsv"
    pd.DataFrame(_rows()).to_csv(path, sep="\t", index=False)
    return lineage_qc_owned(
        path,
        first_column_as_index=False,
        lineage_column="Tumor",
        cell_barcode_column="cellBC",
        integration_barcode_column="intBC",
        cut_site_columns="r1,r2",
        prior_grouping_columns="Tumor,intBC",
        missing_data_allele="",
        allele_representation_threshold=0.98,
        minimum_cells=2,
        maximum_missing_fraction=1.0,
        maximum_uncut_fraction=1.0,
        minimum_unique_fraction=0.0,
        minimum_informative_character_fraction=0.0,
        max_file_mib=512,
        max_matrix_gib=2.0,
    )


def _assert_rng_equal(before, after):
    assert before[0] == after[0]
    assert np.array_equal(before[1], after[1])
    assert before[2:] == after[2:]


def _directory_fingerprint(path):
    digest = hashlib.sha256()
    for member in sorted(path.rglob("*")):
        if member.is_file():
            digest.update(member.relative_to(path).as_posix().encode("utf-8"))
            digest.update(member.read_bytes())
    return digest.hexdigest()


def test_schemas_are_atomic_typed_and_report_code():
    schemas = {node.define_schema().node_id: node.define_schema() for node in LINEAGE_NODE_CLASSES}
    assert [item.id for item in schemas["OpenBioSingleCellCassiopeiaLineageQC"].inputs] == [
        "allele_table_file",
        "first_column_as_index",
        "lineage_column",
        "cell_barcode_column",
        "integration_barcode_column",
        "cut_site_columns",
        "prior_grouping_columns",
        "missing_data_allele",
        "allele_representation_threshold",
        "minimum_cells",
        "maximum_missing_fraction",
        "maximum_uncut_fraction",
        "minimum_unique_fraction",
        "minimum_informative_character_fraction",
        "max_file_mib",
        "max_matrix_gib",
    ]
    assert [(item.display_name, item.io_type) for item in schemas["OpenBioSingleCellCassiopeiaLineageQC"].outputs] == [
        ("characters", CassiopeiaCharactersType.io_type),
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert schemas["OpenBioSingleCellReconstructCassiopeiaTree"].inputs[0].io_type == CassiopeiaCharactersType.io_type
    assert schemas["OpenBioSingleCellCassiopeiaExpansionTest"].inputs[0].io_type == CassiopeiaTreeType.io_type
    assert schemas["OpenBioSingleCellCassiopeiaPlasticity"].inputs[0].io_type == AnnDataType.io_type
    assert schemas["OpenBioSingleCellCassiopeiaPlasticity"].inputs[1].io_type == CassiopeiaTreeType.io_type
    assert "minimum_clade_fraction" not in {
        item.id for item in schemas["OpenBioSingleCellCassiopeiaExpansionTest"].inputs
    }
    assert "summary_key" not in {item.id for item in schemas["OpenBioSingleCellCassiopeiaPlasticity"].inputs}
    assert [schema.node_id for schema in schemas.values()] == [
        node.define_schema().node_id for node in LINEAGE_NODE_CLASSES
    ]
    assert OpenBioSingleCellCassiopeiaLineageQC in LINEAGE_NODE_CLASSES
    assert OpenBioSingleCellReconstructCassiopeiaTree in LINEAGE_NODE_CLASSES
    assert OpenBioSingleCellCassiopeiaExpansionTest in LINEAGE_NODE_CLASSES
    assert OpenBioSingleCellCassiopeiaPlasticity in LINEAGE_NODE_CLASSES


def test_lineage_plot_schemas_are_domain_specific_and_typed():
    expected = {
        "OpenBioSingleCellCassiopeiaLineageQCPlot": [
            ("characters", CassiopeiaCharactersType.io_type),
            ("table", TableResultType.io_type),
        ],
        "OpenBioSingleCellCassiopeiaTreePlot": [
            ("tree", CassiopeiaTreeType.io_type),
            ("adata", AnnDataType.io_type),
            ("annotation_key", "STRING"),
        ],
        "OpenBioSingleCellCassiopeiaExpansionPlot": [
            ("tree", CassiopeiaTreeType.io_type),
            ("table", TableResultType.io_type),
        ],
        "OpenBioSingleCellCassiopeiaPlasticityPlot": [
            ("adata", AnnDataType.io_type),
            ("tree", CassiopeiaTreeType.io_type),
            ("table", TableResultType.io_type),
            ("output_key", "STRING"),
        ],
    }
    schemas = {node.define_schema().node_id: node.define_schema() for node in LINEAGE_NODE_CLASSES}
    for node_id, inputs in expected.items():
        schema = schemas[node_id]
        assert schema.category == "openbio/single-cell/lineage"
        assert [(item.id, item.io_type) for item in schema.inputs] == inputs
        assert [(item.display_name, item.io_type) for item in schema.outputs] == [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ]


def test_qc_strict_summary_denominators_and_defensive_artifact(tmp_path, fake_backend):
    artifact, table, summary, code = _prepare(tmp_path)
    assert isinstance(artifact, CassiopeiaCharacters)
    assert artifact.lineage_ids == ("T1",)
    assert table.loc[0, "status"] == "pass"
    assert table.loc[0, "mean_missing_fraction"] == 0.0
    assert table.loc[0, "mean_uncut_fraction_observed"] == pytest.approx(0.25)
    assert summary["key_results"]["lineages_passed"] == 1
    json.dumps(summary, allow_nan=False)
    compile(code, "<cassiopeia-code>", "exec")
    matrix, priors, state_maps = artifact.copy_lineage("T1")
    matrix.iloc[0, 0] = 99
    priors[0][1] = 0.1
    state_maps[0][1] = "tampered"
    pristine, pristine_priors, pristine_maps = artifact.copy_lineage("T1")
    assert pristine.iloc[0, 0] != 99
    assert pristine_priors[0][1] != 0.1
    assert pristine_maps[0][1] != "tampered"


def test_lineage_qc_plot_uses_exact_artifact_bound_table_without_mutation(tmp_path, fake_backend):
    characters, table, _, _ = _qc_owned(tmp_path)
    before_fingerprint = characters.fingerprint
    before_table = table.table.copy(deep=True)

    plotted, report, code = lineage_qc_plot_owned(characters, table)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["lineages_plotted"] == 1
    assert details["lineages_available"] == 1
    assert details["status_counts"] == {"pass": 1, "warning": 0, "unavailable": 0}
    assert table.parameters["characters_fingerprint"] == characters.fingerprint
    assert table.parameters["qc_table_sha256"] == characters.metadata["qc_table_sha256"]
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_cassiopeia_lineage_qc"](characters, table.table) == plotted.png
    assert characters.fingerprint == before_fingerprint
    assert_frame_equal(table.table, before_table)


def test_lineage_qc_plot_rejects_non_character_artifacts(tmp_path, fake_backend):
    _, table, _, _ = _qc_owned(tmp_path)

    with pytest.raises(TypeError, match="character artifact"):
        lineage_qc_plot_owned(object(), table)


def test_lineage_qc_plot_rejects_internally_inconsistent_table_metadata(tmp_path, fake_backend):
    characters, table, _, _ = _qc_owned(tmp_path)

    with pytest.raises(ValueError, match="table provenance"):
        lineage_qc_plot_owned(characters, replace(table, input_cells=table.input_cells + 1))


def test_lineage_qc_plot_worker_round_trips_artifacts_without_rewriting_inputs(tmp_path, fake_backend):
    characters, table, _, _ = _qc_owned(tmp_path)
    characters_root = tmp_path / "characters"
    characters_root.mkdir()
    write_characters(characters_root, characters)
    table_root = tmp_path / "qc-table"
    table_root.mkdir()
    write_table(
        table_root,
        table.table,
        result_metadata(table),
        json_list_columns=("qc_warnings", "unavailability_reasons"),
    )
    before = (_directory_fingerprint(characters_root), _directory_fingerprint(table_root))
    staging = tmp_path / "qc-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cassiopeia_lineage_qc_plot(
        context,
        {
            "characters": {
                "type": "artifact",
                "path": str(characters_root.resolve()),
                "kind": "OPENBIO_CASSIOPEIA_CHARACTERS",
                "codec": "openbio-cassiopeia-characters",
            },
            "table": {
                "type": "artifact",
                "path": str(table_root.resolve()),
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
            },
        },
        {},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _ = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert (_directory_fingerprint(characters_root), _directory_fingerprint(table_root)) == before


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (
            _rows() + [{"Tumor": "T1", "cellBC": "c1", "intBC": "i1", "r1": "C", "r2": "NONE"}],
            "Conflicting repeated",
        ),
        (
            _rows() + [{"Tumor": "T2", "cellBC": "c1", "intBC": "i2", "r1": "A", "r2": "B"}],
            "not globally unique",
        ),
        (
            [dict(_rows()[0], r1="DeletionNoneLike"), *_rows()[1:]],
            "reserved text",
        ),
    ],
)
def test_qc_rejects_conflicts_cross_lineage_identity_and_reserved_tokens(tmp_path, fake_backend, rows, message):
    with pytest.raises(ValueError, match=message):
        _prepare(tmp_path, rows)


def test_qc_reports_all_missing_and_all_uncut_lineages(tmp_path, fake_backend):
    rows = _rows()
    rows.extend({"Tumor": "M", "cellBC": cell, "intBC": "i2", "r1": "", "r2": ""} for cell in ("m1", "m2"))
    rows.extend({"Tumor": "U", "cellBC": cell, "intBC": "i3", "r1": "NONE", "r2": "NONE"} for cell in ("u1", "u2"))
    artifact, table, summary, _ = _prepare(tmp_path, rows)
    records = table.set_index("lineage_id")
    assert records.loc["M", "status"] == "unavailable"
    assert records.loc["U", "status"] == "unavailable"
    assert "T1" in artifact.lineage_ids
    assert summary["key_results"]["lineages_failed"] == 2


def test_qc_threshold_and_prior_policy_warnings_do_not_block_expert_tree(tmp_path, fake_backend):
    path = tmp_path / "alleles.tsv"
    pd.DataFrame(_rows()).to_csv(path, sep="\t", index=False)
    artifact, table, summary, _ = lineage.prepare_cassiopeia_characters(
        str(path),
        first_column_as_index=False,
        cut_site_columns="r1,r2",
        prior_grouping_columns="intBC",
        minimum_cells=10,
        maximum_missing_fraction=1.0,
        maximum_uncut_fraction=1.0,
        minimum_unique_fraction=0.0,
        minimum_informative_character_fraction=0.0,
        openbio_version="test",
    )
    assert table.loc[0, "status"] == "warning"
    assert "T1" in artifact.lineage_ids
    assert summary["key_results"]["lineages_with_qc_warnings"] == 1
    assert any("independence" in warning for warning in summary["warnings"])
    _, tree_summary, _ = lineage.reconstruct_cassiopeia_tree(artifact, "T1", openbio_version="test")
    assert any("upstream QC thresholds" in warning for warning in tree_summary["warnings"])


def test_character_artifact_detects_internal_tampering(tmp_path, fake_backend):
    artifact, _, _, _ = _prepare(tmp_path)
    artifact._CassiopeiaCharacters__matrices["T1"].iloc[0, 0] = 777
    with pytest.raises(ValueError, match="tampered"):
        artifact.copy_lineage("T1")


def test_tree_is_canonical_rng_safe_and_defensive(tmp_path, fake_backend):
    characters, _, _, _ = _prepare(tmp_path)
    numpy_before = np.random.get_state()
    python_before = random.getstate()
    first, summary, code = lineage.reconstruct_cassiopeia_tree(
        characters,
        "T1",
        prior_transformation="inverse",
        collapse_mutationless_edges=True,
        openbio_version="test",
    )
    numpy_after = np.random.get_state()
    python_after = random.getstate()
    _assert_rng_equal(numpy_before, numpy_after)
    assert python_before == python_after
    second, _, _ = lineage.reconstruct_cassiopeia_tree(
        characters,
        "T1",
        prior_transformation="inverse",
        collapse_mutationless_edges=True,
        openbio_version="test",
    )
    assert isinstance(first, CassiopeiaTree)
    assert first.topology_fingerprint == second.topology_fingerprint
    assert first.root.startswith("__openbio_internal__")
    assert FakeVanillaGreedySolver.last_classifier is assign_missing_average
    assert FakeVanillaGreedySolver.last_transformation == "inverse"
    assert FakeVanillaGreedySolver.last_collapse is True
    assert summary["key_results"]["root_semantics"] == "solver_generated"
    json.dumps(summary, allow_nan=False)
    compile(code, "<cassiopeia-tree-code>", "exec")
    external = first.copy_tree()
    external._graph.add_node("evil")
    assert "evil" not in first.copy_tree()._graph


def test_worker_reconstruction_reuses_the_codec_owned_character_matrix(tmp_path, fake_backend):
    from openbio_singlecell.operations_lineage import reconstruct_owned

    created = []

    class TrackingTree(FakeCassiopeiaTree):
        def __init__(self, *args, character_matrix=None, priors=None, **kwargs):
            super().__init__(*args, character_matrix=character_matrix, priors=priors, **kwargs)
            self.character_matrix = character_matrix
            self.priors = {} if priors is None else priors
            created.append(self)

    fake_backend.data.CassiopeiaTree = TrackingTree
    characters, _, _, _ = _prepare(tmp_path)
    source_matrix, _, _ = characters.copy_lineage("T1", _owned=True)

    tree, _, _ = reconstruct_owned(
        characters,
        lineage_id="T1",
        prior_transformation="negative_log",
        collapse_mutationless_edges=False,
    )

    assert created[0].character_matrix is source_matrix
    assert created[1].character_matrix is source_matrix
    assert tree.copy_tree(_owned=True) is created[1]


def test_tree_artifact_detects_topology_tampering(tmp_path, fake_backend):
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    tree._CassiopeiaTree__tree._graph.add_node("evil")
    with pytest.raises(ValueError, match="topology"):
        tree.copy_tree()


def test_tree_plot_uses_stored_topology_and_exact_leaf_annotations(tmp_path, fake_backend):
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata()
    before_tree = tree.fingerprint
    before_obs = adata.obs.copy(deep=True)
    solves_before = FakeVanillaGreedySolver.solve_count

    plotted, report, code = cassiopeia_tree_plot_owned(tree, adata, annotation_key="cell_type")

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["lineage_id"] == "T1"
    assert details["leaf_count"] == 4
    assert details["annotation_categories"] == ["A", "B"]
    assert details["ann_data_cells_not_in_tree"] == 1
    assert details["topology_sha256"] == tree.topology_fingerprint
    assert FakeVanillaGreedySolver.solve_count == solves_before
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_cassiopeia_tree"](tree, adata) == plotted.png
    assert tree.fingerprint == before_tree
    assert_frame_equal(adata.obs, before_obs)


def test_tree_plot_rejects_non_tree_artifacts():
    with pytest.raises(TypeError, match="tree artifact"):
        cassiopeia_tree_plot_owned(object(), _plasticity_adata(), annotation_key="cell_type")


def test_tree_plot_worker_reads_tree_and_anndata_without_rewriting_them(tmp_path, fake_backend, monkeypatch):
    monkeypatch.setattr("openbio_singlecell.cassiopeia_codec._require_cassiopeia", lambda: fake_backend)
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata()
    tree_root = tmp_path / "tree-input"
    tree_root.mkdir()
    write_tree(tree_root, tree)
    adata_root = tmp_path / "annotation-input"
    adata_root.mkdir()
    write_anndata(adata_root, adata)
    before = (_directory_fingerprint(tree_root), _directory_fingerprint(adata_root))
    staging = tmp_path / "tree-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cassiopeia_tree_plot(
        context,
        {
            "tree": {
                "type": "artifact",
                "path": str(tree_root.resolve()),
                "kind": "OPENBIO_CASSIOPEIA_TREE",
                "codec": "openbio-cassiopeia-tree",
            },
            "adata": {
                "type": "artifact",
                "path": str(adata_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
        },
        {"annotation_key": "cell_type"},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _ = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert (_directory_fingerprint(tree_root), _directory_fingerprint(adata_root)) == before


def test_expansion_independently_checks_formula_and_applies_bh(tmp_path, fake_backend):
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    table, summary, code = lineage.compute_cassiopeia_expansions(
        tree, minimum_clade_size=2, minimum_depth=1, fdr_threshold=0.5, openbio_version="test"
    )
    eligible = table.loc[table["eligible"]]
    assert eligible.shape[0] == 2
    assert eligible["p_value"].tolist() == pytest.approx([2 / 3, 2 / 3])
    assert eligible["q_value_bh"].tolist() == pytest.approx([2 / 3, 2 / 3])
    assert not bool(eligible["significant_fdr"].any())
    assert summary["key_results"]["eligible_hypotheses"] == 2
    json.dumps(summary, allow_nan=False)
    compile(code, "<cassiopeia-expansion-code>", "exec")


def test_expansion_plot_is_bound_to_the_tested_topology_and_full_table(tmp_path, fake_backend):
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    table, _, _ = expansion_owned(tree, minimum_clade_size=2, minimum_depth=1, fdr_threshold=0.5)
    before_tree = tree.fingerprint
    before_table = table.table.copy(deep=True)
    solves_before = FakeVanillaGreedySolver.solve_count

    plotted, report, code = cassiopeia_expansion_plot_owned(tree, table)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["lineage_id"] == "T1"
    assert details["topology_nodes_plotted"] == 7
    assert details["eligible_clades"] == 2
    assert details["significant_clades"] == 0
    assert details["topology_sha256"] == tree.topology_fingerprint
    assert table.parameters["tree_fingerprint"] == tree.fingerprint
    assert table.parameters["topology_sha256"] == tree.topology_fingerprint
    assert FakeVanillaGreedySolver.solve_count == solves_before
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_cassiopeia_expansion"](tree, table.table) == plotted.png
    assert tree.fingerprint == before_tree
    assert_frame_equal(table.table, before_table)


def test_expansion_plot_rejects_internally_inconsistent_table_metadata(tmp_path, fake_backend):
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    table, _, _ = expansion_owned(tree, minimum_clade_size=2, minimum_depth=1, fdr_threshold=0.5)

    with pytest.raises(ValueError, match="table provenance"):
        cassiopeia_expansion_plot_owned(tree, replace(table, input_cells=table.input_cells + 1))


def test_expansion_plot_worker_preserves_tree_and_table_artifacts(tmp_path, fake_backend, monkeypatch):
    monkeypatch.setattr("openbio_singlecell.cassiopeia_codec._require_cassiopeia", lambda: fake_backend)
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    table, _, _ = expansion_owned(tree, minimum_clade_size=2, minimum_depth=1, fdr_threshold=0.5)
    tree_root = tmp_path / "expansion-tree"
    tree_root.mkdir()
    write_tree(tree_root, tree)
    table_root = tmp_path / "expansion-table"
    table_root.mkdir()
    write_table(table_root, table.table, result_metadata(table))
    before = (_directory_fingerprint(tree_root), _directory_fingerprint(table_root))
    staging = tmp_path / "expansion-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cassiopeia_expansion_plot(
        context,
        {
            "tree": {
                "type": "artifact",
                "path": str(tree_root.resolve()),
                "kind": "OPENBIO_CASSIOPEIA_TREE",
                "codec": "openbio-cassiopeia-tree",
            },
            "table": {
                "type": "artifact",
                "path": str(table_root.resolve()),
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
            },
        },
        {},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _ = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert (_directory_fingerprint(tree_root), _directory_fingerprint(table_root)) == before


def test_worker_expansion_keeps_only_the_backend_required_tree_copy(tmp_path, fake_backend):
    from openbio_singlecell.operations_lineage import expansion_owned

    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    FakeCassiopeiaTree.copy_count = 0

    expansion_owned(tree, minimum_clade_size=2, minimum_depth=1, fdr_threshold=0.5)

    assert FakeCassiopeiaTree.copy_count == 1


def test_expansion_rejects_malicious_backend_disagreement(tmp_path, fake_backend, monkeypatch):
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")

    def wrong(tree, min_clade_size=10, min_depth=1, copy=False):
        result = _fake_expansion(tree, min_clade_size, min_depth, copy)
        for node in result._graph:
            if result.get_attribute(node, "expansion_pvalue") != 1.0:
                result.set_attribute(node, "expansion_pvalue", 0.123)
        return result

    fake_backend.tl.compute_expansion_pvalues = wrong
    with pytest.raises(RuntimeError, match="disagrees"):
        lineage.compute_cassiopeia_expansions(tree, minimum_clade_size=2, openbio_version="test")


def _plasticity_adata(states=("A", "A", "B", "B", "C")):
    adata = ad.AnnData(np.zeros((5, 1), dtype=float))
    adata.obs_names = ["c1", "c2", "c3", "c4", "extra"]
    adata.var_names = ["g1"]
    adata.obs["cell_type"] = pd.Categorical(states)
    return adata


def test_plasticity_uses_edge_denominator_preserves_inputs_and_reports_all_cells(tmp_path, fake_backend):
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata()
    before = adata.copy()
    numpy_before = np.random.get_state()
    output, table, summary, code = lineage.add_cassiopeia_plasticity(
        adata,
        tree,
        annotation_key="cell_type",
        annotation_status="unknown",
        analysis_mode="exploratory",
        minimum_state_fraction=0.025,
        openbio_version="test",
    )
    _assert_rng_equal(numpy_before, np.random.get_state())
    assert list(output.obs_names) == list(adata.obs_names)
    assert "sc_effective_plasticity" not in adata.obs
    assert before.obs.equals(adata.obs)
    included = output.obs.loc[["c1", "c2", "c3", "c4"], "sc_effective_plasticity"].astype(float)
    assert included.tolist() == pytest.approx([1 / 12] * 4)
    assert pd.isna(output.obs.loc["extra", "sc_effective_plasticity"])
    assert table.shape[0] == 5
    assert table.loc[table["cell_id"] == "extra", "status"].item() == "not_in_tree"
    assert summary["key_results"]["global_small_parsimony"] == 1
    assert summary["parameters"]["subtree_denominator"] == "directed edge count"
    json.dumps(summary, allow_nan=False)
    compile(code, "<cassiopeia-plasticity-code>", "exec")


def test_plasticity_plot_validates_scores_annotations_table_and_tree(tmp_path, fake_backend):
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    output, table, _, _ = plasticity_owned(
        _plasticity_adata(),
        tree,
        annotation_key="cell_type",
        annotation_status="unknown",
        analysis_mode="exploratory",
        minimum_state_fraction=0.025,
        output_key="sc_effective_plasticity",
        overwrite_existing=False,
        max_working_gib=4.0,
    )
    before_tree = tree.fingerprint
    before_obs = output.obs.copy(deep=True)
    before_uns = copy.deepcopy(output.uns["openbio_cassiopeia_plasticity"])
    solves_before = FakeVanillaGreedySolver.solve_count

    plotted, report, code = cassiopeia_plasticity_plot_owned(
        output,
        tree,
        table,
        output_key="sc_effective_plasticity",
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["lineage_id"] == "T1"
    assert details["included_cells"] == 4
    assert details["annotation_categories"] == ["A", "B"]
    assert details["status_counts"]["not_in_tree"] == 1
    assert details["topology_sha256"] == tree.topology_fingerprint
    assert table.parameters["tree_fingerprint"] == tree.fingerprint
    assert FakeVanillaGreedySolver.solve_count == solves_before
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_cassiopeia_plasticity"](output, tree, table.table) == plotted.png
    assert tree.fingerprint == before_tree
    assert_frame_equal(output.obs, before_obs)
    assert output.uns["openbio_cassiopeia_plasticity"] == before_uns


def test_plasticity_plot_rejects_internally_inconsistent_table_metadata(tmp_path, fake_backend):
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    output, table, _, _ = plasticity_owned(_plasticity_adata(), tree, annotation_key="cell_type")

    with pytest.raises(ValueError, match="table provenance"):
        cassiopeia_plasticity_plot_owned(
            output,
            tree,
            replace(table, input_genes=table.input_genes + 1),
            output_key="sc_effective_plasticity",
        )


def test_plasticity_plot_worker_preserves_all_three_input_artifacts(tmp_path, fake_backend, monkeypatch):
    monkeypatch.setattr("openbio_singlecell.cassiopeia_codec._require_cassiopeia", lambda: fake_backend)
    characters, _, _, _ = _qc_owned(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    output, table, _, _ = plasticity_owned(
        _plasticity_adata(),
        tree,
        annotation_key="cell_type",
        annotation_status="unknown",
        analysis_mode="exploratory",
        minimum_state_fraction=0.025,
        output_key="sc_effective_plasticity",
        overwrite_existing=False,
        max_working_gib=4.0,
    )
    adata_root = tmp_path / "plasticity-adata"
    adata_root.mkdir()
    write_anndata(adata_root, output)
    tree_root = tmp_path / "plasticity-tree"
    tree_root.mkdir()
    write_tree(tree_root, tree)
    table_root = tmp_path / "plasticity-table"
    table_root.mkdir()
    write_table(table_root, table.table, result_metadata(table))
    before = tuple(_directory_fingerprint(path) for path in (adata_root, tree_root, table_root))
    staging = tmp_path / "plasticity-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cassiopeia_plasticity_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(adata_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
            "tree": {
                "type": "artifact",
                "path": str(tree_root.resolve()),
                "kind": "OPENBIO_CASSIOPEIA_TREE",
                "codec": "openbio-cassiopeia-tree",
            },
            "table": {
                "type": "artifact",
                "path": str(table_root.resolve()),
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
            },
        },
        {"output_key": "sc_effective_plasticity"},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _ = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert tuple(_directory_fingerprint(path) for path in (adata_root, tree_root, table_root)) == before


def test_worker_owned_plasticity_updates_its_private_anndata_in_place(tmp_path, fake_backend):
    from openbio_singlecell.operations_lineage import plasticity_owned

    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata()
    FakeCassiopeiaTree.copy_count = 0

    output, _, _, _ = plasticity_owned(adata, tree, annotation_key="cell_type")

    assert output is adata
    assert "sc_effective_plasticity" in adata.obs
    assert FakeCassiopeiaTree.copy_count == 0


def test_plasticity_frequency_equality_exhausted_polytomy_and_report_grade_warning(tmp_path, fake_backend):
    FakeVanillaGreedySolver.mode = "star"
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata(("A", "A", "A", "B", "C"))
    output, _, summary, _ = lineage.add_cassiopeia_plasticity(
        adata,
        tree,
        annotation_key="cell_type",
        minimum_state_fraction=0.25,
        openbio_version="test",
    )
    assert summary["key_results"]["retained_states"] == ["A", "B"]
    assert summary["key_results"]["exhausted_polytomies_split"] == 1
    assert output.obs.loc["c4", "sc_effective_plasticity"] is not pd.NA
    _, _, unverified_report, _ = lineage.add_cassiopeia_plasticity(
        adata,
        tree,
        annotation_key="cell_type",
        annotation_status="curated",
        analysis_mode="report_grade",
        openbio_version="test",
    )
    assert any("without matching registered Curated" in warning for warning in unverified_report["warnings"])
    adata.uns["openbio_singlecell"] = {
        "annotations": {"cell_type": {"annotation_status": "curated", "method": "reviewed"}}
    }
    _, _, report_grade, _ = lineage.add_cassiopeia_plasticity(
        adata,
        tree,
        annotation_key="cell_type",
        annotation_status="curated",
        analysis_mode="report_grade",
        minimum_state_fraction=0.25,
        openbio_version="test",
    )
    assert report_grade["parameters"]["analysis_mode"] == "report_grade"


def test_plasticity_excludes_missing_identity_and_annotation_with_one_state_warning(tmp_path, fake_backend):
    FakeVanillaGreedySolver.mode = "star"
    characters, _, _, _ = _prepare(tmp_path)
    tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    adata = _plasticity_adata(("A", "A", "B", "B", "C"))[["c1", "c2", "c4", "extra"]].copy()
    adata.obs.loc["c4", "cell_type"] = pd.NA
    output, table, summary, _ = lineage.add_cassiopeia_plasticity(
        adata, tree, annotation_key="cell_type", openbio_version="test"
    )
    assert float(output.obs.loc["c1", "sc_effective_plasticity"]) == 0.0
    assert pd.isna(output.obs.loc["c4", "sc_effective_plasticity"])
    assert table.loc[table["cell_id"] == "c4", "status"].item() == "missing_annotation"
    assert summary["key_results"]["tree_leaves_missing_from_anndata"] == 1
    assert summary["key_results"]["tree_leaves_with_missing_annotation"] == 1
    assert summary["key_results"]["global_small_parsimony"] == 0
    assert any("Only one annotation state" in warning for warning in summary["warnings"])


def test_generated_source_compiles_executes_and_matches_all_four_operations(tmp_path, fake_backend):
    characters, qc, _, code = _prepare(tmp_path)
    namespace = {"__name__": "generated_cassiopeia"}
    exec(compile(code, "<generated-cassiopeia>", "exec"), namespace)
    namespace["_require_cassiopeia"] = lambda: fake_backend
    path = tmp_path / "alleles.tsv"
    generated_characters, generated_qc, _, _ = namespace["prepare_cassiopeia_characters"](
        str(path),
        first_column_as_index=False,
        cut_site_columns="r1,r2",
        prior_grouping_columns="Tumor,intBC",
        allele_representation_threshold=0.98,
        minimum_cells=2,
        maximum_missing_fraction=1.0,
        maximum_uncut_fraction=1.0,
        minimum_unique_fraction=0.0,
        minimum_informative_character_fraction=0.0,
        openbio_version="test",
    )
    assert_frame_equal(qc, generated_qc)
    original_tree, _, _ = lineage.reconstruct_cassiopeia_tree(characters, "T1", openbio_version="test")
    generated_tree, _, _ = namespace["reconstruct_cassiopeia_tree"](generated_characters, "T1", openbio_version="test")
    assert original_tree.topology_fingerprint == generated_tree.topology_fingerprint
    original_expansion, _, _ = lineage.compute_cassiopeia_expansions(
        original_tree, minimum_clade_size=2, openbio_version="test"
    )
    generated_expansion, _, _ = namespace["compute_cassiopeia_expansions"](
        generated_tree, minimum_clade_size=2, openbio_version="test"
    )
    assert_frame_equal(original_expansion, generated_expansion)
    adata = _plasticity_adata()
    original_output, original_table, _, _ = lineage.add_cassiopeia_plasticity(
        adata, original_tree, annotation_key="cell_type", openbio_version="test"
    )
    generated_output, generated_table, _, _ = namespace["add_cassiopeia_plasticity"](
        adata, generated_tree, annotation_key="cell_type", openbio_version="test"
    )
    assert_frame_equal(original_table, generated_table)
    assert original_output.obs["sc_effective_plasticity"].equals(generated_output.obs["sc_effective_plasticity"])
