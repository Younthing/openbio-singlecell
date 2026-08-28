from __future__ import annotations

import gzip
import json

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_data import (
    OpenBioSingleCellMapGeneIdsFromGTF,
    OpenBioSingleCellMergeObservationAnnotations,
    OpenBioSingleCellSnapshotExpression,
    OpenBioSingleCellSubsetObservations,
)


def _dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _assert_report(report, code, node_id):
    assert report.kind == "summary"
    assert report.summary["schema_version"] == 1
    assert report.summary["node_id"] == node_id
    assert report.summary["methods"]
    assert report.summary["results"]
    assert report.summary["key_results"]
    assert report.summary["references"]
    assert report.summary["software_versions"]["openbio-singlecell"] == "0.2.0"
    json.dumps(report.summary, allow_nan=False)
    assert code.endswith("\n")
    compile(code, f"<{node_id}-code>", "exec")


@pytest.fixture
def data_adata(science):
    counts = science.np.asarray(
        [
            [1.0, 0.0, 2.0, 0.0],
            [0.0, 3.0, 0.0, 1.0],
            [4.0, 0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0, 2.0],
        ]
    )
    obs = science.pd.DataFrame(
        {"sample": ["a", "b", None, "a"]},
        index=["cell_0", "cell_1", "cell_2", "cell_3"],
    )
    var = science.pd.DataFrame(index=["gene.0", "gene.1", "gene.2", "gene.3"])
    value = science.ad.AnnData(science.sparse.csr_matrix(counts), obs=obs, var=var)
    value.layers["candidate"] = science.sparse.csc_matrix(counts + 10.0)
    value.obsm["coordinates"] = science.np.arange(8, dtype=float).reshape(4, 2)
    value.obsp["links"] = science.sparse.eye(4, format="csr")
    ensure_metadata(value, display_name="synthetic", source={"kind": "test", "study": "same"})
    return value


def test_data_node_schemas_are_current_and_reportable():
    expected_inputs = {
        OpenBioSingleCellSnapshotExpression: ["adata", "source", "overwrite_existing"],
        OpenBioSingleCellSubsetObservations: ["adata", "column", "values", "invert", "missing_policy"],
        OpenBioSingleCellMergeObservationAnnotations: [
            "adata",
            "subset_adata",
            "source_column",
            "target_column",
            "conflict_policy",
        ],
        OpenBioSingleCellMapGeneIdsFromGTF: [
            "adata",
            "gtf_path",
            "id_source",
            "id_column",
            "gene_name_column",
        ],
    }

    for node_class, input_ids in expected_inputs.items():
        schema = node_class.GET_SCHEMA()
        assert [input_.id for input_ in schema.inputs] == input_ids
        assert all(input_.optional is False for input_ in schema.inputs)
        assert [output.display_name for output in schema.outputs] == ["adata", "summary", "code"]


def test_storage_and_selection_nodes_return_explicit_empty_results(science):
    empty = science.ad.AnnData(
        science.np.empty((0, 0)),
        obs=science.pd.DataFrame(
            {"group": science.pd.Series(dtype="string")},
            index=science.pd.Index([], dtype="string"),
        ),
        var=science.pd.DataFrame(index=science.pd.Index([], dtype="string")),
    )
    empty.layers["empty"] = science.np.empty((0, 0))
    ensure_metadata(empty)

    snapshotted, snapshot_report, snapshot_code = OpenBioSingleCellSnapshotExpression.execute(
        empty
    ).result
    subset, subset_report, subset_code = OpenBioSingleCellSubsetObservations.execute(
        empty,
        "group",
        "selected",
    ).result

    for output, report, code, node_id in (
        (snapshotted, snapshot_report, snapshot_code, "OpenBioSingleCellSnapshotExpression"),
        (subset, subset_report, subset_code, "OpenBioSingleCellSubsetObservations"),
    ):
        assert output.shape == (0, 0)
        assert report.summary["warnings"]
        _assert_report(report, code, node_id)

    namespaces = [{}, {}]
    for code, namespace in zip(
        (snapshot_code, subset_code),
        namespaces,
        strict=True,
    ):
        exec(code, namespace)
    assert namespaces[0]["create_raw_snapshot"](empty).shape == (0, 0)
    assert namespaces[1]["subset_observations"](empty).shape == (0, 0)


def test_snapshot_expression_creates_canonical_states_from_declared_layer(data_adata, science):
    original_x = data_adata.X.copy()
    expected = data_adata.layers["candidate"].copy()
    output, report, code = OpenBioSingleCellSnapshotExpression.execute(
        data_adata,
        source={"source": "layer", "source_layer": "candidate"},
        overwrite_existing=False,
    ).result

    science.np.testing.assert_array_equal(_dense(output.X, science), _dense(original_x, science))
    science.np.testing.assert_array_equal(
        _dense(output.layers["counts"], science),
        _dense(expected, science),
    )
    assert output.raw is not None
    science.np.testing.assert_array_equal(_dense(output.raw.X, science), _dense(expected, science))
    assert "counts" not in data_adata.layers
    assert report.summary["key_results"]["source"] == "layer:candidate"
    assert report.summary["key_results"]["x_preserved"] is True
    _assert_report(report, code, "OpenBioSingleCellSnapshotExpression")

    namespace = {}
    exec(code, namespace)
    equivalent = namespace["create_raw_snapshot"](data_adata)
    science.np.testing.assert_array_equal(_dense(equivalent.X, science), _dense(output.X, science))
    science.np.testing.assert_array_equal(
        _dense(equivalent.layers["counts"], science),
        _dense(output.layers["counts"], science),
    )
    science.np.testing.assert_array_equal(_dense(equivalent.raw.X, science), _dense(output.raw.X, science))

    sliced = output[["cell_3", "cell_0"], ["gene.1", "gene.0"]]
    assert sliced.obs_names.tolist() == ["cell_3", "cell_0"]
    assert sliced.var_names.tolist() == ["gene.1", "gene.0"]
    assert sliced.raw is not None
    assert sliced.raw.obs_names.tolist() == ["cell_3", "cell_0"]
    assert sliced.raw.var_names.tolist() == data_adata.var_names.tolist()


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (lambda value, science: setattr(value, "X", science.np.zeros(value.shape)), "no positive"),
        (
            lambda value, science: setattr(value, "X", science.np.full(value.shape, -1.0)),
            "negative expression",
        ),
        (
            lambda value, science: setattr(value, "X", science.np.full(value.shape, science.np.nan)),
            "non-finite",
        ),
        (lambda value, science: setattr(value, "obs_names", ["c", "c", "d", "e"]), "unique"),
        (lambda value, science: setattr(value, "var_names", ["g", "g", "h", "i"]), "unique"),
    ],
)
def test_snapshot_expression_discloses_selected_source_advisories(
    data_adata,
    science,
    mutator,
    message,
):
    mutator(data_adata, science)
    output, report, code = OpenBioSingleCellSnapshotExpression.execute(data_adata).result
    assert output.raw is not None
    assert any(message in warning for warning in report.summary["warnings"])
    assert report.summary["key_results"]["history_used"] is False
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["create_raw_snapshot"](data_adata)
    assert equivalent.shape == output.shape


def test_snapshot_expression_overwrite_and_noninteger_disclosure(data_adata, science):
    first, _, _ = OpenBioSingleCellSnapshotExpression.execute(data_adata).result
    with pytest.raises(ValueError, match="already exist"):
        OpenBioSingleCellSnapshotExpression.execute(first)

    replaced, report, _ = OpenBioSingleCellSnapshotExpression.execute(
        first,
        overwrite_existing=True,
    ).result
    assert report.summary["key_results"]["overwritten_destinations"] == [
        "raw",
        "layers['counts']",
    ]
    assert replaced.raw is not None

    noninteger = data_adata.copy()
    noninteger.X = science.np.asarray(_dense(noninteger.X, science), dtype=float) + 0.25
    _, noninteger_report, _ = OpenBioSingleCellSnapshotExpression.execute(noninteger).result
    assert noninteger_report.summary["key_results"]["integer_like"] is False
    assert any("non-integer" in warning for warning in noninteger_report.summary["warnings"])

def test_snapshot_generated_code_quotes_arbitrary_layer_names(data_adata):
    layer_name = 'candidate"with-quote'
    data_adata.layers[layer_name] = data_adata.X.copy()
    _, _, code = OpenBioSingleCellSnapshotExpression.execute(
        data_adata,
        source={"source": "layer", "source_layer": layer_name},
    ).result
    namespace = {}
    exec(code, namespace)
    output = namespace["create_raw_snapshot"](data_adata)
    assert output.raw is not None


def test_snapshot_rejects_backed_anndata_with_actionable_error(data_adata, science, tmp_path):
    path = tmp_path / "backed.h5ad"
    data_adata.write_h5ad(path)
    backed = science.ad.read_h5ad(path, backed="r")
    try:
        with pytest.raises(ValueError, match=r"to_memory\(\)"):
            OpenBioSingleCellSnapshotExpression.execute(backed)
        _, _, code = OpenBioSingleCellSnapshotExpression.execute(data_adata).result
        namespace = {}
        exec(code, namespace)
        with pytest.raises(ValueError, match=r"to_memory\(\)"):
            namespace["create_raw_snapshot"](backed)
    finally:
        backed.file.close()


def test_subset_observations_has_explicit_missing_semantics_and_equivalent_code(data_adata):
    selected, report, code = OpenBioSingleCellSubsetObservations.execute(
        data_adata,
        "sample",
        "a, absent",
        False,
        "exclude",
    ).result
    assert selected.obs_names.tolist() == ["cell_0", "cell_3"]
    assert report.summary["key_results"]["unmatched_requested_values"] == ["absent"]
    assert any("not present" in warning for warning in report.summary["warnings"])
    _assert_report(report, code, "OpenBioSingleCellSubsetObservations")
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["subset_observations"](data_adata)
    assert equivalent.obs_names.equals(selected.obs_names)

    included, _, _ = OpenBioSingleCellSubsetObservations.execute(
        data_adata,
        "sample",
        "a",
        True,
        "include",
    ).result
    assert included.obs_names.tolist() == ["cell_1", "cell_2"]

    excluded, _, _ = OpenBioSingleCellSubsetObservations.execute(
        data_adata,
        "sample",
        "a",
        True,
        "exclude",
    ).result
    assert excluded.obs_names.tolist() == ["cell_1"]

    with pytest.raises(ValueError, match="missing value"):
        OpenBioSingleCellSubsetObservations.execute(data_adata, "sample", "a", False, "error")
    empty, empty_report, empty_code = OpenBioSingleCellSubsetObservations.execute(
        data_adata,
        "sample",
        "absent",
    ).result
    assert empty.shape == (0, data_adata.n_vars)
    assert empty_report.summary["key_results"]["empty_output"] is True
    assert any("zero cells" in warning for warning in empty_report.summary["warnings"])
    empty_namespace = {}
    exec(empty_code, empty_namespace)
    assert empty_namespace["subset_observations"](data_adata).shape == empty.shape


def test_subset_observations_preserves_aligned_slots_order_and_raw(data_adata, science):
    data_adata.raw = data_adata
    output, report, _ = OpenBioSingleCellSubsetObservations.execute(
        data_adata,
        "sample",
        "b,a",
    ).result

    assert output.obs_names.tolist() == ["cell_0", "cell_1", "cell_3"]
    science.np.testing.assert_array_equal(
        output.obsm["coordinates"],
        data_adata.obsm["coordinates"][[0, 1, 3]],
    )
    science.np.testing.assert_array_equal(
        _dense(output.obsp["links"], science),
        _dense(data_adata.obsp["links"], science)[science.np.ix_([0, 1, 3], [0, 1, 3])],
    )
    assert output.raw is not None
    assert output.raw.obs_names.equals(output.obs_names)
    assert report.summary["key_results"]["observation_order_preserved"] is True

    colliding = data_adata.copy()
    colliding.obs["mixed"] = science.pd.Series(
        [1, "1", 2, 3],
        index=colliding.obs_names,
        dtype="object",
    )
    with pytest.raises(ValueError, match="collide"):
        OpenBioSingleCellSubsetObservations.execute(colliding, "mixed", "1")

    semantically_equal = data_adata.copy()
    semantically_equal.obs["mixed"] = science.pd.Series(
        [1, science.np.int64(1), 2, 3],
        index=semantically_equal.obs_names,
        dtype="object",
    )
    selected_equal, _, equal_code = OpenBioSingleCellSubsetObservations.execute(
        semantically_equal,
        "mixed",
        "1",
    ).result
    assert selected_equal.obs_names.tolist() == ["cell_0", "cell_1"]
    equal_namespace = {}
    exec(equal_code, equal_namespace)
    equivalent_equal = equal_namespace["subset_observations"](semantically_equal)
    assert equivalent_equal.obs_names.equals(selected_equal.obs_names)

    literal_missing = data_adata.copy()
    literal_missing.obs["label"] = ["<missing>", None, science.np.nan, science.pd.NA]
    _, literal_report, _ = OpenBioSingleCellSubsetObservations.execute(
        literal_missing,
        "label",
        "<missing>",
    ).result
    value_counts = literal_report.summary["key_results"]["input_value_counts"]
    counts = value_counts["items"]
    assert any(item["value"] == "<missing>" and item["is_missing"] is False for item in counts)
    missing_items = [item for item in counts if item["is_missing"] is True]
    assert missing_items == [{"value": None, "value_type": None, "is_missing": True, "count": 3}]
    assert value_counts["total_categories"] == 2


def _annotation_inputs(science):
    target = science.ad.AnnData(
        science.sparse.csr_matrix(science.np.eye(4)),
        obs=science.pd.DataFrame(
            {"cell_type": [None, "B", "old", None], "unchanged": [1, 2, 3, 4]},
            index=["c0", "c1", "c2", "c3"],
        ),
        var=science.pd.DataFrame(index=["g0", "g1", "g2", "g3"]),
    )
    source = target[["c2", "c0", "c1"], :].copy()
    source.obs["label"] = science.pd.Categorical(
        ["new", "T", "B"],
        categories=["B", "T", "new"],
    )
    ensure_metadata(target, source={"kind": "test", "study": "shared"})
    ensure_metadata(source, source={"kind": "test", "study": "shared"})
    return target, source


def test_merge_annotations_classifies_conflicts_and_generated_code_matches(science):
    target, source = _annotation_inputs(science)
    with pytest.raises(ValueError, match="conflicting value"):
        OpenBioSingleCellMergeObservationAnnotations.execute(
            target,
            source,
            source_column="label",
            target_column="cell_type",
            conflict_policy="error",
        )

    kept, kept_report, _ = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        source_column="label",
        target_column="cell_type",
        conflict_policy="keep_target",
    ).result
    assert kept.obs["cell_type"].astype("string").fillna("<missing>").tolist() == [
        "T",
        "B",
        "old",
        "<missing>",
    ]
    kept_results = kept_report.summary["key_results"]
    assert (kept_results["filled"], kept_results["identical"], kept_results["conflicts"]) == (1, 1, 1)
    assert kept_results["kept_target"] == 1

    output, report, code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        source_column="label",
        target_column="cell_type",
        conflict_policy="overwrite",
    ).result
    assert output.obs["cell_type"].astype("string").fillna("<missing>").tolist() == [
        "T",
        "B",
        "new",
        "<missing>",
    ]
    assert output.obs["unchanged"].equals(target.obs["unchanged"])
    assert report.summary["key_results"]["overwritten"] == 1
    _assert_report(report, code, "OpenBioSingleCellMergeObservationAnnotations")
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["merge_observation_annotations"](target, source)
    assert equivalent.obs["cell_type"].astype("string").equals(output.obs["cell_type"].astype("string"))


def test_merge_annotations_preserves_new_categorical_dtype_and_validates_identity(science):
    target, source = _annotation_inputs(science)
    output, report, _ = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        "label",
        "new_label",
        conflict_policy="error",
    ).result
    assert isinstance(output.obs["new_label"].dtype, science.pd.CategoricalDtype)
    assert output.obs["new_label"].astype("string").fillna("<missing>").tolist() == [
        "T",
        "B",
        "new",
        "<missing>",
    ]
    assert report.summary["key_results"]["target_dtype_after"] == "category"

    source_only = source.copy()
    source_only.obs_names = ["c2", "c0", "outside"]
    _, source_only_report, source_only_code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source_only,
        "label",
        conflict_policy="keep_target",
    ).result
    assert source_only_report.summary["key_results"]["source_only_cells"] == 1
    assert any("Ignored 1 source" in warning for warning in source_only_report.summary["warnings"])
    source_only_namespace = {}
    exec(source_only_code, source_only_namespace)
    source_only_equivalent = source_only_namespace["merge_observation_annotations"](
        target,
        source_only,
    )
    assert source_only_equivalent.obs_names.equals(target.obs_names)

    conflicting_source = source.copy()
    ensure_metadata(conflicting_source, source={"kind": "test", "study": "different"})
    _, conflicting_report, _ = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        conflicting_source,
        "label",
        conflict_policy="keep_target",
    ).result
    assert conflicting_report.summary["key_results"]["provenance_status"] == "conflicting_advisory"
    assert any("different OpenBio source" in warning for warning in conflicting_report.summary["warnings"])


def test_merge_annotations_preserves_native_ordered_categories_and_rejects_noop_sources(science):
    target, source = _annotation_inputs(science)
    source.obs["integer_label"] = science.pd.Categorical(
        [2, 1, 2],
        categories=[1, 2],
        ordered=True,
    )
    output, _, code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        "integer_label",
        "integer_label",
        conflict_policy="error",
    ).result
    assert isinstance(output.obs["integer_label"].dtype, science.pd.CategoricalDtype)
    assert output.obs["integer_label"].cat.categories.tolist() == [1, 2]
    assert output.obs["integer_label"].cat.ordered is True
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["merge_observation_annotations"](target, source)
    assert equivalent.obs["integer_label"].dtype == output.obs["integer_label"].dtype

    missing_source = source.copy()
    missing_source.obs["missing"] = science.pd.NA
    missing_output, missing_report, missing_code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        missing_source,
        "missing",
    ).result
    assert missing_output.obs["cell_type"].equals(target.obs["cell_type"])
    assert any("no non-missing" in warning for warning in missing_report.summary["warnings"])
    missing_namespace = {}
    exec(missing_code, missing_namespace)
    missing_equivalent = missing_namespace["merge_observation_annotations"](target, missing_source)
    assert missing_equivalent.obs["cell_type"].equals(missing_output.obs["cell_type"])

    categorical_target = target.copy()
    categorical_target.obs["ordered"] = science.pd.Categorical(
        ["A", "A", "A", "A"],
        categories=["A"],
        ordered=True,
    )
    categorical_source = source.copy()
    categorical_source.obs["ordered_source"] = science.pd.Categorical(
        ["B", "B", "B"],
        categories=["B"],
        ordered=True,
    )
    unchanged, unchanged_report, _ = OpenBioSingleCellMergeObservationAnnotations.execute(
        categorical_target,
        categorical_source,
        "ordered_source",
        "ordered",
        conflict_policy="keep_target",
    ).result
    assert unchanged.obs["ordered"].cat.categories.tolist() == ["A"]
    assert unchanged.obs["ordered"].cat.ordered is True
    assert unchanged_report.summary["key_results"]["kept_target"] == 3


def test_merge_annotations_preserves_nullable_dtype_and_reports_incompatible_fallback(science):
    target, source = _annotation_inputs(science)
    target.obs["nullable"] = science.pd.Series(
        [1, science.pd.NA, 3, science.pd.NA],
        index=target.obs_names,
        dtype="Int64",
    )
    source.obs["nullable_source"] = science.pd.Series(
        [3, 1, 2],
        index=source.obs_names,
        dtype="Int64",
    )
    output, report, code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        "nullable_source",
        "nullable",
        conflict_policy="error",
    ).result
    assert str(output.obs["nullable"].dtype) == "Int64"
    assert output.obs["nullable"].tolist() == [1, 2, 3, science.pd.NA]
    assert report.summary["key_results"]["dtype_fallback_to_object"] is False
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["merge_observation_annotations"](target, source)
    assert equivalent.obs["nullable"].dtype == output.obs["nullable"].dtype
    assert equivalent.obs["nullable"].equals(output.obs["nullable"])

    incompatible = source.copy()
    incompatible.obs["nullable_source"] = science.pd.Series(
        [3, 1, "B"],
        index=incompatible.obs_names,
        dtype="object",
    )
    fallback, fallback_report, fallback_code = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        incompatible,
        "nullable_source",
        "nullable",
        conflict_policy="error",
    ).result
    assert str(fallback.obs["nullable"].dtype) == "object"
    assert fallback.obs["nullable"].tolist() == [1, "B", 3, science.pd.NA]
    assert fallback_report.summary["key_results"]["dtype_fallback_to_object"] is True
    assert any("object dtype" in warning for warning in fallback_report.summary["warnings"])
    fallback_namespace = {}
    exec(fallback_code, fallback_namespace)
    fallback_equivalent = fallback_namespace["merge_observation_annotations"](target, incompatible)
    assert fallback_equivalent.obs["nullable"].equals(fallback.obs["nullable"])


def test_merge_annotation_value_counts_keep_distinct_scalar_types(science):
    target, source = _annotation_inputs(science)
    source.obs["mixed"] = science.pd.Series(
        [1, "1", "x"],
        index=source.obs_names,
        dtype="object",
    )
    _, report, _ = OpenBioSingleCellMergeObservationAnnotations.execute(
        target,
        source,
        "mixed",
        "mixed",
        conflict_policy="error",
    ).result
    counts = report.summary["key_results"]["final_value_counts"]
    assert counts["total_categories"] == 4
    assert any(item["value"] == 1 and item["value_type"] == "builtins.int" for item in counts["items"])
    assert any(
        item["value"] == "1" and item["value_type"] == "builtins.str" for item in counts["items"]
    )


def _gtf_line(gene_id, gene_name, feature="gene"):
    return (
        f"chr1\ttest\t{feature}\t1\t2\t.\t+\t.\t"
        f'gene_id "{gene_id}"; gene_name "{gene_name}";\n'
    )


def test_gtf_annotation_keeps_stable_identity_and_generated_code_matches(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    gtf_path = input_dir / "genes.gtf"
    gtf_path.write_text(
        "#!provider GENCODE\n"
        "#!release 49\n"
        "#!genome-build GRCh38\n"
        + _gtf_line("ENSG000001.1", "A")
        + _gtf_line("ENSG000002.7_PAR_Y", "B")
        + _gtf_line("ENSG000003", "A"),
        encoding="utf-8",
    )
    value = science.ad.AnnData(
        science.sparse.csr_matrix(science.np.eye(3)),
        var=science.pd.DataFrame(
            {"gene_ids": ["ENSG000001.1", "ENSG000002.9_PAR_Y", "ENSG000004"]},
            index=["display-A", "display-B", "display-C"],
        ),
    )

    output, report, code = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        value,
        "genes.gtf",
        id_source="var_column",
        id_column=" gene_ids ",
        gene_name_column=" gene_symbols ",
    ).result
    assert output.var_names.tolist() == ["ENSG000001.1", "ENSG000002.9_PAR_Y", "ENSG000004"]
    assert output.var["feature_names_before_gtf"].tolist() == ["display-A", "display-B", "display-C"]
    assert output.var["gene_symbols"].astype("string").fillna("<missing>").tolist() == [
        "A",
        "B",
        "<missing>",
    ]
    assert output.var["gtf_mapped"].tolist() == [True, True, False]
    assert output.var["gtf_match_type"].tolist() == ["exact", "version_normalized", "unmapped"]
    results = report.summary["key_results"]
    assert (results["exact_matches"], results["version_normalized_matches"], results["unmapped_features"]) == (
        1,
        1,
        1,
    )
    assert results["duplicate_gene_symbol_occurrences"] == 0
    assert results["gtf"]["provider"] == "GENCODE"
    assert results["gtf"]["release"] == "49"
    assert results["gtf"]["genome_build"] == "GRCh38"
    assert results["gtf"]["path"] == "genes.gtf"
    assert report.summary["parameters"]["id_column"] == "gene_ids"
    assert report.summary["parameters"]["gene_name_column"] == "gene_symbols"
    _assert_report(report, code, "OpenBioSingleCellMapGeneIdsFromGTF")

    namespace = {}
    exec(code, namespace)
    equivalent = namespace["annotate_gene_ids_from_gtf"](value, gtf_path)
    assert equivalent.var_names.equals(output.var_names)
    assert equivalent.var["gene_symbols"].astype("string").equals(
        output.var["gene_symbols"].astype("string")
    )
    assert equivalent.var["gtf_mapped"].equals(output.var["gtf_mapped"])


def test_gtf_annotation_supports_gzip_fallback_and_rejects_ambiguity(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    fallback_path = input_dir / "fallback.gtf.gz"
    with gzip.open(fallback_path, "wt", encoding="utf-8") as handle:
        handle.write("#!genome-version GRCh38\n")
        handle.write(_gtf_line("ENSG1.1_PAR_Y", "A", feature="transcript"))
        handle.write(_gtf_line("ENSG2", "A", feature="exon"))
    value = science.ad.AnnData(
        science.np.ones((2, 2)),
        var=science.pd.DataFrame(index=["ENSG1.9_PAR_Y", "ENSG2"]),
    )
    output, report, _ = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        value,
        "fallback.gtf.gz",
        id_source="var_names",
    ).result
    assert output.var_names.equals(value.var_names)
    assert output.var["gene_symbols"].tolist() == ["A", "A"]
    assert report.summary["key_results"]["duplicate_gene_symbol_occurrences"] == 1
    assert report.summary["key_results"]["gtf"]["used_fallback_records"] is True
    assert report.summary["key_results"]["gtf"]["release"] is None
    assert report.summary["key_results"]["gtf"]["genome_build"] == "GRCh38"

    ambiguous_path = input_dir / "ambiguous.gtf"
    ambiguous_path.write_text(
        _gtf_line("ENSG1.1", "A") + _gtf_line("ENSG1.2", "B"),
        encoding="utf-8",
    )
    ambiguous = science.ad.AnnData(science.np.ones((1, 1)))
    ambiguous.var_names = ["ENSG1"]
    with pytest.raises(ValueError, match="ambiguous"):
        OpenBioSingleCellMapGeneIdsFromGTF.execute(
            ambiguous,
            "ambiguous.gtf",
            id_source="var_names",
        )


@pytest.mark.parametrize(
    "reserved",
    ["feature_names_before_gtf", "gtf_mapped", "gtf_match_type"],
)
def test_gtf_annotation_rejects_reserved_symbol_columns(comfy_directories, science, reserved):
    input_dir, _, _ = comfy_directories
    (input_dir / "genes.gtf").write_text(_gtf_line("ENSG1", "A"), encoding="utf-8")
    value = science.ad.AnnData(science.np.ones((1, 1)))
    value.var_names = ["ENSG1"]
    with pytest.raises(ValueError, match="reserved"):
        OpenBioSingleCellMapGeneIdsFromGTF.execute(
            value,
            "genes.gtf",
            id_source="var_names",
            gene_name_column=reserved,
        )


def test_gtf_annotation_discloses_duplicate_ids_and_existing_symbol_conflicts(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    (input_dir / "genes.gtf").write_text(
        _gtf_line("ENSG1.1", "A") + _gtf_line("ENSG2", "B"),
        encoding="utf-8",
    )

    duplicate = science.ad.AnnData(science.np.ones((1, 2)))
    duplicate.var_names = ["ENSG1.1", "ENSG1.2"]
    duplicate_output, duplicate_report, duplicate_code = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        duplicate,
        "genes.gtf",
        id_source="var_names",
    ).result
    assert duplicate_output.var["gene_symbols"].tolist() == ["A", "A"]
    assert duplicate_report.summary["key_results"]["version_normalized_id_collisions"] == 1
    assert any("version removal" in warning for warning in duplicate_report.summary["warnings"])
    duplicate_namespace = {}
    exec(duplicate_code, duplicate_namespace)
    duplicate_equivalent = duplicate_namespace["annotate_gene_ids_from_gtf"](
        duplicate,
        input_dir / "genes.gtf",
    )
    assert duplicate_equivalent.var["gene_symbols"].tolist() == ["A", "A"]

    conflict = science.ad.AnnData(
        science.np.ones((1, 2)),
        var=science.pd.DataFrame(
            {"gene_symbols": [None, "wrong"]},
            index=["ENSG1.1", "ENSG2"],
        ),
    )
    conflict_output, conflict_report, conflict_code = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        conflict,
        "genes.gtf",
        id_source="var_names",
    ).result
    assert conflict_output.var["gene_symbols"].tolist() == ["A", "wrong"]
    assert conflict_report.summary["key_results"]["existing_name_conflicts"] == 1
    assert any("disagreed" in warning for warning in conflict_report.summary["warnings"])
    conflict_namespace = {}
    exec(conflict_code, conflict_namespace)
    conflict_equivalent = conflict_namespace["annotate_gene_ids_from_gtf"](
        conflict,
        input_dir / "genes.gtf",
    )
    assert conflict_equivalent.var["gene_symbols"].tolist() == ["A", "wrong"]

    unmatched = science.ad.AnnData(science.np.ones((1, 1)))
    unmatched.var_names = ["UNMATCHED"]
    unmatched_output, unmatched_report, unmatched_code = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        unmatched,
        "genes.gtf",
        id_source="var_names",
    ).result
    assert unmatched_output.var["gtf_mapped"].tolist() == [False]
    assert unmatched_report.summary["key_results"]["mapped_features"] == 0
    assert any("None of the declared" in warning for warning in unmatched_report.summary["warnings"])
    unmatched_namespace = {}
    exec(unmatched_code, unmatched_namespace)
    unmatched_equivalent = unmatched_namespace["annotate_gene_ids_from_gtf"](
        unmatched,
        input_dir / "genes.gtf",
    )
    assert unmatched_equivalent.var["gtf_mapped"].tolist() == [False]


def test_gtf_annotation_discloses_non_string_symbols_and_counts_final_duplicates(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    (input_dir / "genes.gtf").write_text(_gtf_line("ENSG1", "A"), encoding="utf-8")
    invalid = science.ad.AnnData(
        science.np.ones((1, 1)),
        var=science.pd.DataFrame({"gene_symbols": [123]}, index=["ENSG1"]),
    )
    invalid_output, invalid_report, invalid_code = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        invalid,
        "genes.gtf",
        id_source="var_names",
    ).result
    assert invalid_output.var["gene_symbols"].tolist() == [123]
    assert invalid_report.summary["key_results"]["existing_non_string_gene_names"] == 1
    assert any("non-string" in warning for warning in invalid_report.summary["warnings"])
    invalid_namespace = {}
    exec(invalid_code, invalid_namespace)
    assert invalid_namespace["annotate_gene_ids_from_gtf"](
        invalid,
        input_dir / "genes.gtf",
    ).var["gene_symbols"].tolist() == [123]

    duplicate_final = science.ad.AnnData(
        science.np.ones((1, 2)),
        var=science.pd.DataFrame(
            {"gene_symbols": [None, "A"]},
            index=["ENSG1", "UNMAPPED"],
        ),
    )
    output, report, _ = OpenBioSingleCellMapGeneIdsFromGTF.execute(
        duplicate_final,
        "genes.gtf",
        id_source="var_names",
    ).result
    assert output.var["gene_symbols"].tolist() == ["A", "A"]
    assert report.summary["key_results"]["duplicate_gene_symbol_occurrences"] == 1
