from __future__ import annotations

import json

import pytest

from openbio_singlecell.nodes_preprocess import (
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeToLayer,
    OpenBioSingleCellNormalizeTotal,
)
from openbio_singlecell.operations_preprocess import log1p, normalize_to_layer, normalize_total
from tests.artifact_operation_harness import run_anndata_operation


def _normalize_total(adata, target_sum=10000.0, source=None):
    return run_anndata_operation(
        normalize_total,
        adata,
        {"target_sum": target_sum, "source": source or {"source": "X"}},
    )


def _log1p(adata):
    return run_anndata_operation(log1p, adata, {})


def _normalize_to_layer(
    adata,
    source=None,
    target_sum=10000.0,
    transform="log1p",
    output_layer="log1p_norm",
    overwrite_existing=False,
):
    return run_anndata_operation(
        normalize_to_layer,
        adata,
        {
            "source": source or {"source": "X"},
            "target_sum": target_sum,
            "transform": transform,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
        },
    )


def _dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _named_layer_keys(adata):
    return {key for key in adata.layers.keys() if key is not None}


def _adata(science, *, sparse: bool = True, include_states: bool = True):
    counts = science.np.asarray(
        [
            [1.0, 2.0, 0.0, 3.0],
            [4.0, 0.0, 2.0, 2.0],
            [0.0, 3.0, 6.0, 1.0],
        ]
    )
    matrix = science.sparse.csr_matrix(counts) if sparse else counts
    adata = science.ad.AnnData(
        X=matrix.copy(),
        obs=science.pd.DataFrame(index=["cell_a", "cell_b", "cell_c"]),
        var=science.pd.DataFrame(index=["gene_a", "gene_b", "gene_c", "gene_d"]),
    )
    if include_states:
        adata.layers["counts"] = matrix.copy()
        adata.layers["alternate"] = matrix.copy()
        adata.obsm["X_test"] = science.np.arange(6, dtype=float).reshape(3, 2)
        adata.varm["test_loadings"] = science.np.arange(8, dtype=float).reshape(4, 2)
        adata.obsp["test_graph"] = science.sparse.eye(3, format="csr")
        adata.varp["test_graph"] = science.sparse.eye(4, format="csr")
        adata.uns["test_metadata"] = {"label": "preserve", "values": science.np.asarray([1, 2])}
        adata.raw = adata.copy()
    return adata


def _assert_matrix_equal(actual, expected, science):
    science.np.testing.assert_allclose(_dense(actual, science), _dense(expected, science), rtol=1e-12, atol=1e-12)


def _assert_uns_value_equal(actual, expected, science):
    if isinstance(expected, dict):
        assert isinstance(actual, dict)
        assert set(actual) == set(expected)
        for key in expected:
            _assert_uns_value_equal(actual[key], expected[key], science)
    elif isinstance(expected, science.pd.DataFrame):
        science.pd.testing.assert_frame_equal(actual, expected)
    elif isinstance(expected, science.pd.Series):
        science.pd.testing.assert_series_equal(actual, expected)
    elif hasattr(expected, "shape"):
        science.np.testing.assert_array_equal(actual, expected)
    else:
        assert actual == expected


def _assert_axis_mapping_equal(actual, expected, science):
    assert set(actual.keys()) == set(expected.keys())
    for key in expected.keys():
        if isinstance(expected[key], science.pd.DataFrame):
            science.pd.testing.assert_frame_equal(actual[key], expected[key])
        else:
            _assert_matrix_equal(actual[key], expected[key], science)


def _assert_expression_storage_equal(actual, expected, science):
    _assert_matrix_equal(actual.X, expected.X, science)
    science.pd.testing.assert_frame_equal(actual.obs, expected.obs)
    science.pd.testing.assert_frame_equal(actual.var, expected.var)
    assert _named_layer_keys(actual) == _named_layer_keys(expected)
    for key in _named_layer_keys(expected):
        _assert_matrix_equal(actual.layers[key], expected.layers[key], science)
    assert (actual.raw is None) is (expected.raw is None)
    if expected.raw is not None:
        _assert_matrix_equal(actual.raw.X, expected.raw.X, science)
        science.pd.testing.assert_frame_equal(actual.raw.var, expected.raw.var)
    actual_uns = {key: value for key, value in actual.uns.items() if key != "openbio_singlecell"}
    expected_uns = {key: value for key, value in expected.uns.items() if key != "openbio_singlecell"}
    _assert_uns_value_equal(actual_uns, expected_uns, science)
    for actual_mapping, expected_mapping in (
        (actual.obsm, expected.obsm),
        (actual.varm, expected.varm),
        (actual.obsp, expected.obsp),
        (actual.varp, expected.varp),
    ):
        _assert_axis_mapping_equal(actual_mapping, expected_mapping, science)


def _run_code(code: str, function_name: str, adata):
    namespace: dict[str, object] = {}
    exec(code, namespace)
    return namespace[function_name](adata)


def _assert_report(node_id: str, report, code: str):
    assert report.kind == "summary"
    assert report.summary["schema_version"] == 1
    assert report.summary["node_id"] == node_id
    assert report.summary["methods"]
    assert report.summary["results"]
    assert report.summary["key_results"]
    assert report.summary["references"]
    assert report.summary["software_versions"]["scanpy"] != "not-installed"
    json.dumps(report.summary, allow_nan=False)
    assert code.endswith("\n")
    compile(code, f"<{node_id}-code>", "exec")


def test_basic_preprocessing_schemas_expose_primary_summary_and_code():
    normalize_schema = OpenBioSingleCellNormalizeTotal.define_schema()
    log_schema = OpenBioSingleCellLog1p.define_schema()
    layer_schema = OpenBioSingleCellNormalizeToLayer.define_schema()

    assert [item.id for item in normalize_schema.inputs] == ["adata", "target_sum", "source"]
    assert normalize_schema.inputs[-1].optional is False
    assert [item.id for item in log_schema.inputs] == ["adata"]
    assert [item.id for item in layer_schema.inputs] == [
        "adata",
        "source",
        "target_sum",
        "transform",
        "output_layer",
        "overwrite_existing",
    ]
    assert layer_schema.inputs[-1].optional is False
    assert layer_schema.inputs[-1].advanced is True
    for schema in (normalize_schema, log_schema, layer_schema):
        assert [(item.display_name, item.io_type) for item in schema.outputs] == [
            ("adata", "OPENBIO_ANNDATA"),
            ("summary", "OPENBIO_SINGLE_CELL_SUMMARY"),
            ("code", "STRING"),
        ]


@pytest.mark.parametrize(
    ("operation", "function_name"),
    [
        ("normalize_total", "normalize_total_expression"),
        ("log1p", "log1p_expression"),
        ("normalize_to_layer", "normalize_expression_to_layer"),
    ],
)
def test_generated_basic_preprocessing_code_does_not_write_openbio_history(science, operation, function_name):
    adata = _adata(science, include_states=False)
    if operation == "normalize_total":
        _, _, code = _normalize_total(adata, 100.0).result
    elif operation == "log1p":
        _, _, code = _log1p(adata).result
    else:
        _, _, code = _normalize_to_layer(
            adata,
            {"source": "X"},
            100.0,
            "none",
            "normalized",
        ).result

    generated = _run_code(code, function_name, adata)
    assert "openbio_singlecell" not in generated.uns


@pytest.mark.parametrize("sparse", [False, True])
def test_normalize_total_uses_declared_source_and_preserves_all_stored_states(science, sparse):
    adata = _adata(science, sparse=sparse)
    adata.layers["counts"] = adata.layers["counts"] * 2
    original = adata.copy()

    output, report, code = _normalize_total(
        adata,
        100.0,
        {"source": "layer", "source_layer": "counts"},
    ).result

    science.np.testing.assert_allclose(_dense(output.X, science).sum(axis=1), 100.0)
    for key in _named_layer_keys(original):
        _assert_matrix_equal(output.layers[key], original.layers[key], science)
        _assert_matrix_equal(adata.layers[key], original.layers[key], science)
    _assert_matrix_equal(output.raw.X, original.raw.X, science)
    _assert_matrix_equal(adata.X, original.X, science)
    assert report.summary["key_results"]["source"] == "layer:counts"
    assert report.summary["key_results"]["input_integer_like"] is True
    assert report.summary["parameters"]["exclude_highly_expressed"] is False
    assert report.summary["key_results"]["raw_preserved"] is True
    _assert_report("OpenBioSingleCellNormalizeTotal", report, code)

    generated = _run_code(code, "normalize_total_expression", adata)
    _assert_expression_storage_equal(generated, output, science)
    with pytest.raises(ValueError, match="source layer not found"):
        _run_code(code, "normalize_total_expression", _adata(science, include_states=False))


def test_normalize_total_never_creates_counts_or_raw(science):
    adata = _adata(science, include_states=False)

    output, report, _ = _normalize_total(adata, 10_000.0).result

    assert _named_layer_keys(adata) == set()
    assert _named_layer_keys(output) == set()
    assert adata.raw is None
    assert output.raw is None
    assert report.summary["warnings"] == []


@pytest.mark.parametrize("target_sum", [0.0, -1.0, float("nan"), float("inf")])
def test_normalize_total_rejects_invalid_target(science, target_sum):
    with pytest.raises(ValueError, match="finite and greater than zero"):
        _normalize_total(_adata(science), target_sum)


def test_normalize_total_warns_for_zero_library_and_allows_repeated_explicit_normalization(science):
    zero_library = _adata(science, sparse=False)
    zero_library.X = zero_library.X.copy()
    zero_library.X[1, :] = 0
    zero_output, zero_report, _ = _normalize_total(zero_library, 100.0).result
    science.np.testing.assert_array_equal(zero_output.X[1], science.np.zeros(zero_library.n_vars))
    assert zero_report.summary["key_results"]["zero_total_cells"] == 1
    assert any("leaves those rows unchanged" in warning for warning in zero_report.summary["warnings"])

    normalized = _normalize_total(_adata(science), 100.0).result[0]
    repeated, report, _ = _normalize_total(normalized, 100.0).result
    assert repeated.shape == normalized.shape
    assert any("non-integer values" in warning for warning in report.summary["warnings"])
    assert not any("provenance" in warning or "state" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize(("bad_value", "message"), [(float("nan"), "non-finite")])
def test_count_normalization_nodes_reject_invalid_expression(science, bad_value, message):
    for execute in (
        lambda value: _normalize_total(value, 100.0),
        lambda value: _normalize_to_layer(
            value,
            {"source": "X"},
            100.0,
            "none",
            "derived",
        ),
    ):
        adata = _adata(science, sparse=False)
        adata.X[0, 0] = bad_value
        with pytest.raises(ValueError, match=message):
            execute(adata)


def test_count_normalization_nodes_allow_signed_finite_expert_input_with_warning(science):
    for execute in (
        lambda value: _normalize_total(value, 100.0),
        lambda value: _normalize_to_layer(
            value,
            {"source": "X"},
            100.0,
            "none",
            "derived",
        ),
    ):
        adata = _adata(science, sparse=False, include_states=False)
        adata.X[0, 0] = -0.5
        output, report, code = execute(adata).result
        assert any("negative expression values" in warning for warning in report.summary["warnings"])
        function_name = (
            "normalize_total_expression"
            if "def normalize_total_expression" in code
            else "normalize_expression_to_layer"
        )
        with pytest.warns(UserWarning) as caught:
            generated = _run_code(code, function_name, adata)
        assert any("negative expression values" in str(item.message) for item in caught)
        target = generated.X if function_name == "normalize_total_expression" else generated.layers["derived"]
        assert science.np.isfinite(_dense(target, science)).all()
        assert output.shape == generated.shape


@pytest.mark.parametrize(
    ("transform", "message"),
    [("log1p", "greater than -1"), ("sqrt", "non-negative normalized values")],
)
def test_normalize_to_layer_retains_transform_math_domain_hard_gates(science, transform, message):
    adata = _adata(science, sparse=False, include_states=False)
    adata.X[0, 0] = -0.5
    with pytest.raises(ValueError, match=message):
        _normalize_to_layer(
            adata,
            {"source": "X"},
            100.0,
            transform,
            "derived",
        )


def test_log1p_uses_the_explicit_source_without_inferred_state(science):
    normalized = _normalize_total(_adata(science), 100.0).result[0]
    before = normalized.copy()

    output, report, code = _log1p(normalized).result

    science.np.testing.assert_allclose(_dense(output.X, science), science.np.log1p(_dense(before.X, science)))
    for key in _named_layer_keys(before):
        _assert_matrix_equal(output.layers[key], before.layers[key], science)
    _assert_matrix_equal(output.raw.X, before.raw.X, science)
    _assert_matrix_equal(normalized.X, before.X, science)
    assert report.summary["warnings"] == []
    assert report.summary["key_results"]["base"] == "natural"
    assert "source_state" not in report.summary["key_results"]
    assert "source_state_evidence" not in report.summary["key_results"]
    _assert_report("OpenBioSingleCellLog1p", report, code)

    generated = _run_code(code, "log1p_expression", normalized)
    _assert_expression_storage_equal(generated, output, science)
    before_second = generated.X.copy()
    generated_twice = _run_code(code, "log1p_expression", generated)
    science.np.testing.assert_allclose(
        _dense(generated_twice.X, science),
        science.np.log1p(_dense(before_second, science)),
    )


def test_log1p_does_not_create_raw_or_layers(science):
    without_raw = _adata(science, include_states=False)
    output = _log1p(without_raw).result[0]
    assert output.raw is None
    assert _named_layer_keys(output) == set()

def test_log1p_rejects_invalid_domain_and_nonfinite_but_allows_repeated_explicit_transform(science):
    negative = _adata(science, sparse=False, include_states=False)
    negative.X[0, 0] = -1
    with pytest.raises(ValueError, match="finite real log1p domain"):
        _log1p(negative)

    nonfinite = _adata(science, sparse=False, include_states=False)
    nonfinite.X[0, 0] = science.np.nan
    with pytest.raises(ValueError, match="non-finite"):
        _log1p(nonfinite)

    logged = _log1p(_adata(science, include_states=False)).result[0]
    twice, report, _ = _log1p(logged).result
    assert science.np.isfinite(_dense(twice.X, science)).all()
    assert not any("already log-transformed" in warning for warning in report.summary["warnings"])

    signed = _adata(science, sparse=False, include_states=False)
    signed.X[0, 0] = -0.5
    transformed, signed_report, _ = _log1p(signed).result
    assert science.np.isfinite(transformed.X).all()
    assert any("values in (-1, 0)" in warning for warning in signed_report.summary["warnings"])


def test_normalize_total_clears_stale_log_marker_before_followup_log1p(science):
    adata = _adata(science, sparse=False)
    adata.X = science.np.log1p(adata.X)
    adata.uns["log1p"] = {"base": None}

    normalized, report, normalize_code = _normalize_total(
        adata,
        100.0,
        {"source": "layer", "source_layer": "counts"},
    ).result

    assert "log1p" not in normalized.uns
    assert report.summary["key_results"]["stale_log1p_marker_removed"] is True
    logged, log_report, log_code = _log1p(normalized).result
    assert log_report.summary["warnings"] == []

    generated_normalized = _run_code(normalize_code, "normalize_total_expression", adata)
    assert "log1p" not in generated_normalized.uns
    generated_logged = _run_code(log_code, "log1p_expression", generated_normalized)
    _assert_expression_storage_equal(generated_logged, logged, science)


@pytest.mark.parametrize("transform", ["none", "log1p", "sqrt"])
def test_normalize_to_layer_builds_one_derived_layer_and_generated_code_matches(science, transform):
    adata = _adata(science)
    original = adata.copy()

    output, report, code = _normalize_to_layer(
        adata,
        {"source": "layer", "source_layer": "counts"},
        100.0,
        transform,
        f"derived_{transform}",
        False,
    ).result

    work = science.ad.AnnData(X=adata.layers["counts"].copy())
    science.sc.pp.normalize_total(work, target_sum=100.0)
    if transform == "log1p":
        science.sc.pp.log1p(work)
    elif transform == "sqrt":
        work.X = work.X.sqrt() if science.sparse.issparse(work.X) else science.np.sqrt(work.X)
    _assert_matrix_equal(output.layers[f"derived_{transform}"], work.X, science)
    _assert_matrix_equal(output.X, original.X, science)
    _assert_matrix_equal(output.layers["counts"], original.layers["counts"], science)
    _assert_matrix_equal(output.raw.X, original.raw.X, science)
    assert report.summary["key_results"]["transform"] == transform
    assert report.summary["key_results"]["input_integer_like"] is True
    assert report.summary["parameters"]["exclude_highly_expressed"] is False
    assert report.summary["key_results"]["replaced_existing_layer"] is False
    _assert_report("OpenBioSingleCellNormalizeToLayer", report, code)

    generated = _run_code(code, "normalize_expression_to_layer", adata)
    _assert_expression_storage_equal(generated, output, science)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"target_sum": 0.0}, "finite and greater than zero"),
        ({"target_sum": float("nan")}, "finite and greater than zero"),
        ({"transform": "bogus"}, "Unsupported"),
        ({"output_layer": "  "}, "cannot be empty"),
        ({"output_layer": "counts"}, "already exists"),
    ],
)
def test_normalize_to_layer_rejects_invalid_configuration(science, kwargs, message):
    values = {
        "source": {"source": "layer", "source_layer": "counts"},
        "target_sum": 100.0,
        "transform": "log1p",
        "output_layer": "derived",
        "overwrite_existing": False,
        **kwargs,
    }
    with pytest.raises(ValueError, match=message):
        _normalize_to_layer(_adata(science), **values)


def test_normalize_to_layer_rejects_unauthorized_output_collisions(science):
    adata = _adata(science)
    with pytest.raises(ValueError, match="already exists"):
        _normalize_to_layer(
            adata,
            {"source": "layer", "source_layer": "alternate"},
            100.0,
            "none",
            "alternate",
        )

    adata.layers["derived"] = adata.X.copy()
    with pytest.raises(ValueError, match="already exists"):
        _normalize_to_layer(
            adata,
            {"source": "layer", "source_layer": "counts"},
            100.0,
            "none",
            "derived",
        )
    output, report, _ = _normalize_to_layer(
        adata,
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "none",
        "derived",
        True,
    ).result
    science.np.testing.assert_allclose(_dense(output.layers["derived"], science).sum(axis=1), 100.0)
    assert report.summary["key_results"]["replaced_existing_layer"] is True


def test_normalize_to_layer_warns_for_zero_library_and_allows_repeated_explicit_source(science):
    zero_library = _adata(science, sparse=False)
    zero_library.layers["counts"] = zero_library.layers["counts"].copy()
    zero_library.layers["counts"][0, :] = 0
    zero_output, zero_report, _ = _normalize_to_layer(
        zero_library,
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "none",
        "derived",
    ).result
    science.np.testing.assert_array_equal(zero_output.layers["derived"][0], science.np.zeros(zero_library.n_vars))
    assert zero_report.summary["key_results"]["zero_total_cells"] == 1

    derived = _normalize_to_layer(
        _adata(science),
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "log1p",
        "derived",
    ).result[0]
    second, report, _ = _normalize_to_layer(
        derived,
        {"source": "layer", "source_layer": "derived"},
        100.0,
        "none",
        "second",
    ).result
    assert second.layers["second"].shape == derived.shape
    assert any("non-integer values" in warning for warning in report.summary["warnings"])
    assert not any("provenance" in warning or "state" in warning for warning in report.summary["warnings"])


def test_normalize_to_layer_reports_noninteger_external_count_state(science):
    adata = _adata(science, sparse=False)
    adata.layers["counts"][0, 0] = 1.5

    output, report, code = _normalize_to_layer(
        adata,
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "none",
        "derived",
    ).result

    assert report.summary["key_results"]["input_integer_like"] is False
    assert any("non-integer" in warning for warning in report.summary["warnings"])
    with pytest.warns(UserWarning) as generated_warnings:
        generated = _run_code(code, "normalize_expression_to_layer", adata)
    assert [str(item.message) for item in generated_warnings] == report.summary["warnings"]
    _assert_expression_storage_equal(generated, output, science)


def test_generated_normalization_code_reproduces_zero_library_warning_and_output(science):
    adata = _adata(science, sparse=False)
    _, _, code = _normalize_to_layer(
        adata,
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "none",
        "derived",
    ).result
    adata.layers["counts"] = adata.layers["counts"].copy()
    adata.layers["counts"][2, :] = 0

    with pytest.warns(UserWarning) as caught:
        generated = _run_code(code, "normalize_expression_to_layer", adata)
    assert any("zero total expression" in str(item.message) for item in caught)
    science.np.testing.assert_array_equal(generated.layers["derived"][2], science.np.zeros(adata.n_vars))

    without_layer = _adata(science, sparse=False, include_states=False)
    with pytest.raises(ValueError, match="source layer not found"):
        _run_code(code, "normalize_expression_to_layer", without_layer)


def test_three_basic_preprocessing_nodes_reject_empty_and_accept_backed_at_file_boundary(science, tmp_path):
    empty = science.ad.AnnData(X=science.np.empty((0, 2)))
    with pytest.raises(ValueError, match="at least one cell and one feature"):
        _normalize_total(empty)
    with pytest.raises(ValueError, match="at least one cell and one feature"):
        _log1p(empty)
    with pytest.raises(ValueError, match="at least one cell and one feature"):
        _normalize_to_layer(empty)

    path = tmp_path / "backed.h5ad"
    _adata(science).write_h5ad(path)
    for operation in (
        _normalize_total,
        _log1p,
        _normalize_to_layer,
    ):
        backed = science.ad.read_h5ad(path, backed="r")
        try:
            assert operation(backed).result[0].isbacked is False
        finally:
            backed.file.close()
