from __future__ import annotations

import json

import pytest

from openbio_singlecell.nodes_preprocess import (
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellPearsonResidualsToLayer,
    OpenBioSingleCellScale,
)
from openbio_singlecell.operations_preprocess import (
    highly_variable_genes,
    normalize_to_layer,
    pearson_residuals_to_layer,
    scale,
)
from tests.artifact_operation_harness import run_anndata_operation


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
            "source": {"source": "X"} if source is None else source,
            "target_sum": target_sum,
            "transform": transform,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
        },
    )


def _pearson(
    adata,
    source=None,
    theta=100.0,
    clipping_mode="sqrt_n_obs",
    custom_clip=10.0,
    output_layer="analytic_pearson_residuals",
    overwrite_existing=False,
    max_dense_gib=2.0,
):
    return run_anndata_operation(
        pearson_residuals_to_layer,
        adata,
        {
            "source": {"source": "layer", "source_layer": "counts"} if source is None else source,
            "theta": theta,
            "clipping_mode": clipping_mode,
            "custom_clip": custom_clip,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
            "max_dense_gib": max_dense_gib,
        },
    )


def _hvg(
    adata,
    n_top_genes=2000,
    flavor="seurat",
    source=None,
    batch_key="",
    always_keep_genes="",
    subset=False,
    theta=100.0,
    clipping_mode="sqrt_n_obs",
    custom_clip=10.0,
    chunksize=1000,
    span=0.3,
    n_bins=20,
    overwrite_existing=False,
):
    return run_anndata_operation(
        highly_variable_genes,
        adata,
        {
            "n_top_genes": n_top_genes,
            "flavor": flavor,
            "source": {"source": "layer", "layer_name": "log1p_norm"} if source is None else source,
            "batch_key": batch_key,
            "always_keep_genes": always_keep_genes,
            "subset": subset,
            "theta": theta,
            "clipping_mode": clipping_mode,
            "custom_clip": custom_clip,
            "chunksize": chunksize,
            "span": span,
            "n_bins": n_bins,
            "overwrite_existing": overwrite_existing,
        },
    )


def _scale(
    adata,
    source=None,
    zero_center=True,
    clipping_mode="custom",
    custom_max_value=10.0,
    output_layer="scaled",
    overwrite_existing=False,
    max_dense_gib=2.0,
):
    return run_anndata_operation(
        scale,
        adata,
        {
            "source": {"source": "layer", "layer_name": "log1p_norm"} if source is None else source,
            "zero_center": zero_center,
            "clipping_mode": clipping_mode,
            "custom_max_value": custom_max_value,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
            "max_dense_gib": max_dense_gib,
        },
    )


def _dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _adata(science, *, sparse: bool = True, raw: bool = True):
    counts = science.np.asarray(
        [
            [7, 1, 0, 2, 4, 1],
            [8, 0, 1, 2, 5, 2],
            [6, 2, 1, 3, 4, 1],
            [9, 1, 2, 1, 3, 2],
            [1, 7, 2, 4, 1, 3],
            [2, 8, 1, 5, 2, 4],
            [1, 6, 3, 4, 1, 5],
            [2, 9, 2, 3, 2, 4],
        ],
        dtype=float,
    )
    matrix = science.sparse.csr_matrix(counts) if sparse else counts
    obs = science.pd.DataFrame(
        {"batch": ["a"] * 4 + ["b"] * 4},
        index=[f"cell_{index}" for index in range(counts.shape[0])],
    )
    var = science.pd.DataFrame(index=[f"gene_{index}" for index in range(counts.shape[1])])
    adata = science.ad.AnnData(X=matrix.copy(), obs=obs, var=var)
    adata.layers["counts"] = matrix.copy()
    if raw:
        adata.raw = adata.copy()
    return adata


def _logged(science, *, sparse: bool = True, raw: bool = True):
    return _normalize_to_layer(
        _adata(science, sparse=sparse, raw=raw),
        {"source": "layer", "source_layer": "counts"},
        100.0,
        "log1p",
        "log1p_norm",
    ).result[0]


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
    assert report.summary["references"]
    assert report.summary["software_versions"]["scanpy"] != "not-installed"
    json.dumps(report.summary, allow_nan=False)
    compile(code, f"<{node_id}-code>", "exec")


def _assert_preserved(actual, original, science):
    science.np.testing.assert_allclose(_dense(actual.X, science), _dense(original.X, science))
    science.np.testing.assert_allclose(
        _dense(actual.layers["counts"], science), _dense(original.layers["counts"], science)
    )
    assert (actual.raw is None) is (original.raw is None)
    if original.raw is not None:
        science.np.testing.assert_allclose(_dense(actual.raw.X, science), _dense(original.raw.X, science))


def test_feature_modeling_schemas_expose_primary_summary_and_code():
    pearson = OpenBioSingleCellPearsonResidualsToLayer.define_schema()
    hvg = OpenBioSingleCellHighlyVariableGenes.define_schema()
    scale = OpenBioSingleCellScale.define_schema()

    assert [item.id for item in pearson.inputs] == [
        "adata",
        "source",
        "theta",
        "clipping_mode",
        "custom_clip",
        "output_layer",
        "overwrite_existing",
        "max_dense_gib",
    ]
    assert [item.id for item in hvg.inputs][-7:] == [
        "theta",
        "clipping_mode",
        "custom_clip",
        "chunksize",
        "span",
        "n_bins",
        "overwrite_existing",
    ]
    assert [item.id for item in scale.inputs] == [
        "adata",
        "source",
        "zero_center",
        "clipping_mode",
        "custom_max_value",
        "output_layer",
        "overwrite_existing",
        "max_dense_gib",
    ]
    for schema in (pearson, hvg, scale):
        assert [(item.display_name, item.io_type) for item in schema.outputs] == [
            ("adata", "OPENBIO_ANNDATA"),
            ("summary", "OPENBIO_SINGLE_CELL_SUMMARY"),
            ("code", "STRING"),
        ]
        source = next(item for item in schema.inputs if item.id == "source")
        assert source.options[0].key == "layer"


@pytest.mark.parametrize(
    ("mode", "custom", "scanpy_clip"), [("sqrt_n_obs", 4.0, None), ("custom", 0.0, 0.0), ("none", 4.0, float("inf"))]
)
def test_pearson_residuals_matches_scanpy_preserves_source_and_code(science, mode, custom, scanpy_clip):
    adata = _adata(science)
    original = adata.copy()

    with pytest.warns(UserWarning, match="dense float64"):
        output, report, code = _pearson(
            adata,
            source={"source": "layer", "source_layer": "counts"},
            theta=75.0,
            output_layer="residuals",
            clipping_mode=mode,
            custom_clip=custom,
            max_dense_gib=1.0,
        ).result

    expected = science.sc.experimental.pp.normalize_pearson_residuals(
        adata,
        theta=75.0,
        clip=scanpy_clip,
        check_values=True,
        layer="counts",
        inplace=False,
    )["X"]
    science.np.testing.assert_allclose(output.layers["residuals"], expected)
    _assert_preserved(output, original, science)
    assert report.summary["parameters"]["clipping_mode"] == mode
    assert report.summary["key_results"]["sparse_to_dense"] is True
    _assert_report("OpenBioSingleCellPearsonResidualsToLayer", report, code)

    generated = _run_code(code, "pearson_residuals_to_layer", adata)
    science.np.testing.assert_allclose(generated.layers["residuals"], output.layers["residuals"])
    _assert_preserved(generated, original, science)


def test_pearson_residuals_warns_for_fractional_counts_and_rejects_undefined_inputs(science):
    noninteger = _adata(science, sparse=False)
    noninteger.layers["counts"][0, 0] = 1.5
    output, report, code = _pearson(
        noninteger,
        source={"source": "layer", "source_layer": "counts"},
    ).result
    assert science.np.isfinite(output.layers["analytic_pearson_residuals"]).all()
    assert report.summary["key_results"]["input_integer_like"] is False
    assert any("fractional non-negative" in warning for warning in report.summary["warnings"])
    with pytest.warns(UserWarning, match="fractional non-negative"):
        generated = _run_code(code, "pearson_residuals_to_layer", noninteger)
    science.np.testing.assert_allclose(
        generated.layers["analytic_pearson_residuals"],
        output.layers["analytic_pearson_residuals"],
    )

    zero_cell = _adata(science, sparse=False)
    zero_cell.layers["counts"][0, :] = 0
    with pytest.raises(ValueError, match="cells with zero total counts"):
        _pearson(
            zero_cell,
            source={"source": "layer", "source_layer": "counts"},
        )

    zero_gene = _adata(science, sparse=False)
    zero_gene.layers["counts"][:, 0] = 0
    with pytest.raises(ValueError, match="genes with zero total counts"):
        _pearson(
            zero_gene,
            source={"source": "layer", "source_layer": "counts"},
        )

    collision = _adata(science)
    collision.layers["residuals"] = collision.X.copy()
    with pytest.raises(ValueError, match="already exists"):
        _pearson(
            collision,
            source={"source": "layer", "source_layer": "counts"},
            output_layer="residuals",
        )
    with pytest.raises(ValueError, match="exceeding max_dense_gib"):
        _pearson(
            _adata(science),
            source={"source": "layer", "source_layer": "counts"},
            max_dense_gib=1e-12,
        )


def test_pearson_residuals_worker_result_is_materialized_and_generated_code_operates_in_place(science):
    adata = _adata(science, sparse=False)
    view = adata[:, :]
    assert view.is_view

    output, _, code = _pearson(
        view,
        source={"source": "layer", "source_layer": "counts"},
        output_layer="residuals",
    ).result

    assert view.is_view
    assert not output.is_view
    generated = _run_code(code, "pearson_residuals_to_layer", view)
    assert generated is view
    assert not view.is_view
    science.np.testing.assert_allclose(generated.layers["residuals"], output.layers["residuals"])


def test_pearson_residuals_allows_singleton_axes_with_disclosure(science):
    singleton_cell = _adata(science, sparse=False)[:1, :].copy()
    singleton_cell.layers["counts"] = singleton_cell.layers["counts"] + 1.0
    singleton_feature = _adata(science, sparse=False)[:, :1].copy()
    for subset in (singleton_cell, singleton_feature):
        output, report, code = _pearson(
            subset,
            source={"source": "layer", "source_layer": "counts"},
        ).result
        assert science.np.isfinite(output.layers["analytic_pearson_residuals"]).all()
        if subset.n_obs == 1:
            assert report.summary["key_results"]["residual_variances_ddof1"]["missing"] == subset.n_vars
        generated = _run_code(code, "pearson_residuals_to_layer", subset)
        science.np.testing.assert_allclose(
            generated.layers["analytic_pearson_residuals"],
            output.layers["analytic_pearson_residuals"],
        )


def test_hvg_log_flavor_attributes_forced_genes_and_code(science):
    logged = _logged(science)
    baseline = _hvg(
        logged,
        n_top_genes=3,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
        n_bins=4,
    ).result[0]
    algorithm_count = int(baseline.var["highly_variable"].sum())
    forced_gene = str(baseline.var_names[~baseline.var["highly_variable"].astype(bool)][0])

    output, report, code = _hvg(
        logged,
        n_top_genes=3,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
        always_keep_genes=f"{forced_gene},missing_gene,{forced_gene}",
        n_bins=4,
    ).result

    assert int(output.var["highly_variable_algorithm"].sum()) == algorithm_count
    assert bool(output.var.loc[forced_gene, "highly_variable_forced"])
    assert bool(output.var.loc[forced_gene, "highly_variable"])
    assert report.summary["key_results"]["forced_added_count"] == 1
    assert report.summary["key_results"]["final_selected_count"] == algorithm_count + 1
    assert report.summary["key_results"]["missing_forced_genes"] == ["missing_gene"]
    _assert_report("OpenBioSingleCellHighlyVariableGenes", report, code)

    generated = _run_code(code, "select_highly_variable_genes", logged)
    for column in ("highly_variable", "highly_variable_algorithm", "highly_variable_forced"):
        science.np.testing.assert_array_equal(generated.var[column], output.var[column])


def test_hvg_ignores_external_scanpy_log_marker(science):
    external = _adata(science, sparse=False)
    science.sc.pp.normalize_total(external, target_sum=100.0)
    science.sc.pp.log1p(external)

    output, report, code = _hvg(
        external,
        n_top_genes=3,
        flavor="seurat",
        source={"source": "X"},
        n_bins=4,
    ).result

    assert "source_state" not in report.summary["key_results"]
    assert "source_state_evidence" not in report.summary["key_results"]
    generated = _run_code(code, "select_highly_variable_genes", external)
    science.np.testing.assert_array_equal(generated.var["highly_variable"], output.var["highly_variable"])


def test_hvg_honors_explicit_sources_and_enforces_structural_contracts(science):
    counts = _adata(science)
    counts_output, _, _ = _hvg(
        counts,
        n_top_genes=3,
        flavor="seurat",
        source={"source": "layer", "layer_name": "counts"},
    ).result
    assert "highly_variable" in counts_output.var

    logged = _logged(science, raw=False)
    pearson_output, pearson_report, pearson_code = _hvg(
        logged,
        n_top_genes=3,
        flavor="pearson_residuals",
        source={"source": "layer", "layer_name": "log1p_norm"},
    ).result
    assert "highly_variable" in pearson_output.var
    assert any("fractional non-negative" in warning for warning in pearson_report.summary["warnings"])
    with pytest.warns(UserWarning, match="fractional non-negative"):
        generated_pearson = _run_code(pearson_code, "select_highly_variable_genes", logged.copy())
    science.np.testing.assert_array_equal(
        generated_pearson.var["highly_variable"],
        pearson_output.var["highly_variable"],
    )

    subset_output, subset_report, subset_code = _hvg(
        logged,
        n_top_genes=3,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
        subset=True,
    ).result
    assert subset_output.raw is None
    assert not any("Raw snapshot" in warning for warning in subset_report.summary["warnings"])
    generated_subset = _run_code(subset_code, "select_highly_variable_genes", logged)
    science.np.testing.assert_array_equal(generated_subset.var_names, subset_output.var_names)

    invalid_batch = _logged(science)
    invalid_batch.obs.loc[invalid_batch.obs_names[0], "batch"] = None
    with pytest.raises(ValueError, match="missing labels"):
        _hvg(
            invalid_batch,
            n_top_genes=3,
            flavor="seurat",
            source={"source": "layer", "layer_name": "log1p_norm"},
            batch_key="batch",
        )

    first = _hvg(
        _logged(science),
        n_top_genes=3,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
    ).result[0]
    with pytest.raises(ValueError, match="annotations already exist"):
        _hvg(
            first,
            n_top_genes=3,
            flavor="seurat",
            source={"source": "layer", "layer_name": "log1p_norm"},
        )
    replaced, replaced_report, _ = _hvg(
        first,
        n_top_genes=2,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
        overwrite_existing=True,
        subset=True,
    ).result
    assert replaced.n_vars == replaced_report.summary["key_results"]["final_selected_count"]
    assert replaced.raw is not None
    assert replaced.raw.n_vars == counts.n_vars


def test_hvg_pearson_flavor_resolves_custom_zero_and_reports_batches(science):
    output, report, code = _hvg(
        _adata(science),
        n_top_genes=3,
        flavor="pearson_residuals",
        source={"source": "layer", "layer_name": "counts"},
        batch_key="batch",
        theta=50.0,
        clipping_mode="custom",
        custom_clip=0.0,
        chunksize=2,
    ).result

    assert int(output.var["highly_variable_algorithm"].sum()) == 3
    assert report.summary["parameters"]["resolved_clip"] == 0.0
    assert report.summary["key_results"]["batch_group_sizes"] == [
        {"label": "a", "count": 4},
        {"label": "b", "count": 4},
    ]
    expected_intersection = int(output.var["highly_variable_intersection"].sum())
    assert report.summary["key_results"]["highly_variable_intersection_count"] == expected_intersection
    assert report.summary["key_results"]["highly_variable_intersection_rate"] == expected_intersection / output.n_vars
    _assert_report("OpenBioSingleCellHighlyVariableGenes", report, code)


def test_hvg_grouped_pearson_auto_clip_is_explicit_and_top_features_are_algorithmic(science):
    adata = _adata(science)
    baseline = _hvg(
        adata,
        n_top_genes=2,
        flavor="pearson_residuals",
        source={"source": "layer", "layer_name": "counts"},
        batch_key="batch",
        clipping_mode="sqrt_n_obs",
        chunksize=2,
    ).result[0]
    forced_gene = str(baseline.var_names[~baseline.var["highly_variable"].astype(bool)][0])

    output, report, code = _hvg(
        adata,
        n_top_genes=2,
        flavor="pearson_residuals",
        source={"source": "layer", "layer_name": "counts"},
        batch_key="batch",
        always_keep_genes=forced_gene,
        clipping_mode="sqrt_n_obs",
        chunksize=2,
    ).result

    resolved_clip = float(science.np.sqrt(adata.n_obs))
    assert report.summary["parameters"]["resolved_clip"] == resolved_clip
    expected = adata.copy()
    science.sc.experimental.pp.highly_variable_genes(
        expected,
        theta=100.0,
        clip=resolved_clip,
        n_top_genes=2,
        batch_key="batch",
        chunksize=2,
        flavor="pearson_residuals",
        check_values=True,
        layer="counts",
        subset=False,
        inplace=True,
    )
    science.np.testing.assert_allclose(output.var["residual_variances"], expected.var["residual_variances"])
    science.np.testing.assert_array_equal(output.var["highly_variable_algorithm"], expected.var["highly_variable"])
    algorithm_names = set(output.var_names[output.var["highly_variable_algorithm"]])
    assert set(report.summary["key_results"]["top_algorithm_features"]) <= algorithm_names
    assert forced_gene not in report.summary["key_results"]["top_algorithm_features"]

    generated = _run_code(code, "select_highly_variable_genes", adata)
    science.np.testing.assert_allclose(generated.var["residual_variances"], output.var["residual_variances"])
    science.np.testing.assert_array_equal(generated.var["highly_variable"], output.var["highly_variable"])


def test_hvg_allows_oversized_request_and_pearson_singleton_group(science):
    logged = _logged(science)
    oversized, report, code = _hvg(
        logged,
        n_top_genes=logged.n_vars + 5,
        flavor="seurat",
        source={"source": "layer", "layer_name": "log1p_norm"},
        n_bins=4,
    ).result
    assert report.summary["key_results"]["requested_algorithm_count"] == logged.n_vars + 5
    assert report.summary["key_results"]["algorithm_selected_count"] <= logged.n_vars
    assert any("return at most" in warning for warning in report.summary["warnings"])
    with pytest.warns(UserWarning):
        generated = _run_code(code, "select_highly_variable_genes", logged)
    science.np.testing.assert_array_equal(generated.var["highly_variable"], oversized.var["highly_variable"])

    singleton_batch = _adata(science, sparse=False)
    singleton_batch.obs["batch"] = ["singleton", *(["other"] * (singleton_batch.n_obs - 1))]
    output, singleton_report, singleton_code = _hvg(
        singleton_batch,
        n_top_genes=3,
        flavor="pearson_residuals",
        source={"source": "layer", "layer_name": "counts"},
        batch_key="batch",
    ).result
    assert "highly_variable" in output.var
    assert any("singleton group" in warning for warning in singleton_report.summary["warnings"])
    with pytest.warns(UserWarning, match="singleton group"):
        generated_singleton = _run_code(singleton_code, "select_highly_variable_genes", singleton_batch)
    science.np.testing.assert_array_equal(
        generated_singleton.var["highly_variable"],
        output.var["highly_variable"],
    )


@pytest.mark.parametrize("zero_center", [False, True])
def test_scale_matches_scanpy_preserves_source_storage_and_code(science, zero_center):
    adata = _logged(science, sparse=True)
    original = adata.copy()
    kwargs = {
        "source": {"source": "layer", "layer_name": "log1p_norm"},
        "zero_center": zero_center,
        "clipping_mode": "custom",
        "custom_max_value": 0.0,
        "output_layer": "scaled_test",
        "max_dense_gib": 1.0,
    }
    if zero_center:
        with pytest.warns(UserWarning, match="densifies"):
            output, report, code = _scale(adata, **kwargs).result
    else:
        output, report, code = _scale(adata, **kwargs).result

    expected = science.ad.AnnData(X=adata.layers["log1p_norm"].copy())
    if zero_center:
        with pytest.warns(UserWarning, match="densifies"):
            science.sc.pp.scale(expected, zero_center=zero_center, max_value=0.0)
    else:
        science.sc.pp.scale(expected, zero_center=zero_center, max_value=0.0)
    science.np.testing.assert_allclose(_dense(output.layers["scaled_test"], science), _dense(expected.X, science))
    assert science.sparse.issparse(output.layers["scaled_test"]) is (not zero_center)
    _assert_preserved(output, original, science)
    assert report.summary["parameters"]["resolved_max_value"] == 0.0
    _assert_report("OpenBioSingleCellScale", report, code)

    if zero_center:
        with pytest.warns(UserWarning, match="densifies"):
            generated = _run_code(code, "scale_expression_to_layer", adata)
    else:
        generated = _run_code(code, "scale_expression_to_layer", adata)
    science.np.testing.assert_allclose(
        _dense(generated.layers["scaled_test"], science),
        _dense(output.layers["scaled_test"], science),
    )


def test_scale_none_collision_and_sparse_center_memory_contract(science):
    adata = _logged(science)
    output, report, _ = _scale(
        adata,
        source={"source": "layer", "layer_name": "log1p_norm"},
        zero_center=False,
        clipping_mode="none",
        output_layer="unclipped",
        max_dense_gib=1e-12,
    ).result
    assert science.sparse.issparse(output.layers["unclipped"])
    assert report.summary["parameters"]["resolved_max_value"] is None

    with pytest.raises(ValueError, match="already exists"):
        _scale(
            output,
            source={"source": "layer", "layer_name": "log1p_norm"},
            output_layer="unclipped",
        )
    with pytest.raises(ValueError, match="exceeding max_dense_gib"):
        _scale(
            adata,
            source={"source": "layer", "layer_name": "log1p_norm"},
            zero_center=True,
            clipping_mode="none",
            max_dense_gib=1e-12,
        )
