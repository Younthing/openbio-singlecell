from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_results import (
    OpenBioSingleCellMarkerExpressionPlot,
    OpenBioSingleCellUMAPPlot,
)
from openbio_singlecell.operations_results import marker_expression_plot_owned, umap_plot_owned
from openbio_singlecell.result_plotting import _marker_expression_plot_impl, _plot_umap_impl

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@pytest.fixture
def plot_adata(science):
    obs_names = [f"cell_{index}" for index in range(12)]
    genes = ["G1", "G2", "G3", "G4"]
    groups = science.pd.Categorical(
        ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
        categories=["C", "A", "B"],
        ordered=True,
    )
    counts = science.np.asarray(
        [
            [9, 1, 4, 0],
            [8, 1, 3, 1],
            [7, 2, 4, 0],
            [9, 2, 5, 1],
            [1, 9, 3, 0],
            [2, 8, 4, 1],
            [1, 7, 5, 0],
            [2, 9, 4, 1],
            [2, 3, 9, 0],
            [1, 4, 8, 1],
            [2, 3, 7, 0],
            [1, 4, 9, 1],
        ],
        dtype=float,
    )
    obs = science.pd.DataFrame({"cluster": groups}, index=obs_names)
    var = science.pd.DataFrame(index=genes)
    adata = science.ad.AnnData(science.sparse.csr_matrix(counts + 100.0), obs=obs, var=var)
    adata.layers["log1p_norm"] = science.sparse.csr_matrix(science.np.log1p(counts))
    raw = science.ad.AnnData(
        science.sparse.csr_matrix(science.np.column_stack([counts, science.np.arange(12, dtype=float)])),
        obs=obs.copy(),
        var=science.pd.DataFrame(index=[*genes, "RAW_ONLY"]),
    )
    adata.raw = raw
    adata.obsm["custom_umap"] = science.pd.DataFrame(
        science.np.column_stack(
            [science.np.linspace(-2.0, 2.0, 12), science.np.linspace(3.0, -3.0, 12)]
        ),
        index=obs_names,
    )
    ensure_metadata(adata, display_name="plot fixture", source={"kind": "test"})
    adata.uns["openbio_singlecell"]["analysis_history"]["000000"] = {
        "operation": "normalize_to_layer",
        "parameters": {"output_layer": "log1p_norm", "transform": "log1p"},
    }
    return adata


def _dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _assert_adata_unchanged(actual, expected, science):
    science.np.testing.assert_array_equal(_dense(actual.X, science), _dense(expected.X, science))
    science.pd.testing.assert_frame_equal(actual.obs, expected.obs)
    science.pd.testing.assert_frame_equal(actual.var, expected.var)
    assert list(actual.layers) == list(expected.layers)
    for key in actual.layers:
        science.np.testing.assert_array_equal(_dense(actual.layers[key], science), _dense(expected.layers[key], science))
    assert list(actual.obsm) == list(expected.obsm)
    for key in actual.obsm:
        left, right = actual.obsm[key], expected.obsm[key]
        if isinstance(left, science.pd.DataFrame):
            science.pd.testing.assert_frame_equal(left, right)
        else:
            science.np.testing.assert_array_equal(left, right)
    assert actual.uns == expected.uns
    if actual.raw is not None:
        science.np.testing.assert_array_equal(_dense(actual.raw.X, science), _dense(expected.raw.X, science))
        science.pd.testing.assert_frame_equal(actual.raw.var, expected.raw.var)


def _assert_rng_state_equal(left, right, science):
    assert left[0] == right[0]
    science.np.testing.assert_array_equal(left[1], right[1])
    assert left[2:] == right[2:]


def test_plot_schemas_are_frozen_and_hide_inactive_marker_options():
    umap = OpenBioSingleCellUMAPPlot.define_schema()
    assert [item.id for item in umap.inputs] == [
        "adata",
        "embedding_key",
        "color",
        "color_mode",
        "point_size",
        "continuous_color_map",
        "categorical_palette",
        "sort_order",
        "missing_color",
        "legend_policy",
    ]
    assert [item.display_name for item in umap.outputs] == ["plot", "summary", "code"]
    assert umap.inputs[1].default == "X_umap"
    assert umap.inputs[2].default == "leiden"
    assert umap.inputs[3].default == "auto"
    assert umap.inputs[4].default == 10.0
    assert umap.inputs[5].default == "viridis"
    assert umap.inputs[6].default == "tab20"
    assert umap.inputs[7].default is True
    assert umap.inputs[8].default == "lightgray"
    assert umap.inputs[9].default == "automatic"

    marker = OpenBioSingleCellMarkerExpressionPlot.define_schema()
    assert [item.id for item in marker.inputs] == [
        "adata",
        "genes",
        "groupby",
        "plot",
        "source",
        "group_order",
        "random_seed",
    ]
    assert [item.display_name for item in marker.outputs] == ["plot", "summary", "code"]
    dynamic = marker.inputs[3]
    assert [option.key for option in dynamic.options] == ["dotplot", "matrixplot", "tracksplot", "violin"]
    assert [[item.id for item in option.inputs] for option in dynamic.options] == [
        ["standard_scale", "expression_cutoff", "mean_only_expressed"],
        ["standard_scale"],
        [],
        ["density_norm", "show_cells", "y_scale"],
    ]
    source = marker.inputs[4]
    assert [option.key for option in source.options] == ["layer", "X", "raw"]
    assert source.options[0].inputs[0].default == "log1p_norm"
    assert marker.inputs[5].default == "observed"
    assert marker.inputs[6].default == 0
    assert marker.inputs[6].advanced is True


def test_umap_categorical_missing_order_summary_code_and_immutability(plot_adata, science):
    plot_adata.obs["label"] = science.pd.Categorical(
        ["beta", "alpha", None, "beta", "alpha", "beta", "alpha", "beta", "alpha", "beta", "alpha", "beta"],
        categories=["beta", "alpha"],
        ordered=True,
    )
    snapshot = plot_adata.copy()

    plotted, summary, code = umap_plot_owned(
        plot_adata,
        embedding_key="custom_umap",
        color="label",
        color_mode="auto",
        point_size=7.5,
        continuous_color_map="viridis",
        categorical_palette="Set2",
        sort_order=True,
        missing_color="#cccccc",
        legend_policy="show",
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert summary.summary["key_results"]["color"]["category_order"] == ["beta", "alpha"]
    assert summary.summary["key_results"]["color"]["missing_count"] == 1
    assert summary.summary["key_results"]["color"]["legend_shown"] is True
    assert list(summary.summary["key_results"]["color"]["category_colors"]) == ["beta", "alpha"]
    assert any("special color" in warning for warning in summary.summary["warnings"])
    assert "scipy" in summary.summary["software_versions"]
    json.dumps(summary.summary, allow_nan=False)
    compile(code, "<umap-plot-code>", "exec")
    namespace = {}
    exec(code, namespace)
    generated = namespace["plot_umap"](plot_adata)
    assert generated == plotted.png
    _assert_adata_unchanged(plot_adata, snapshot, science)


def test_umap_continuous_missing_sorting_range_and_constant_disclosure(plot_adata, science):
    plot_adata.obs["score"] = [3.0, 1.0, science.np.nan, 2.0, 4.0, 5.0, 8.0, 7.0, 6.0, 0.0, 9.0, 10.0]
    _, summary, _ = umap_plot_owned(
        plot_adata,
        embedding_key="custom_umap",
        color="score",
        color_mode="continuous",
        sort_order=True,
    )
    color = summary.summary["key_results"]["color"]
    assert color["missing_count"] == 1
    assert color["finite_count"] == 11
    assert color["observed_range"] == [0.0, 10.0]
    assert "stable ascending" in color["draw_order"]
    assert any("special color" in warning for warning in summary.summary["warnings"])

    plot_adata.obs["score"] = 2.0
    _, constant, _ = umap_plot_owned(
        plot_adata,
        embedding_key="custom_umap",
        color="score",
        color_mode="continuous",
    )
    normalization = constant.summary["key_results"]["color"]["normalization"]
    assert normalization["vmin"] < 2.0 < normalization["vmax"]
    assert any("constant" in warning for warning in constant.summary["warnings"])


@pytest.mark.parametrize(
    ("coordinates", "message"),
    [
        (lambda np, n: np.ones((n - 1, 2)), "n_dimensions >= 2"),
        (lambda np, n: np.ones((n, 1)), "n_dimensions >= 2"),
        (lambda np, n: np.column_stack([np.ones(n), np.full(n, np.nan)]), "non-finite"),
        (lambda np, n: np.ones((n, 2), dtype=bool), "real numeric"),
    ],
)
def test_umap_rejects_malformed_coordinates(plot_adata, science, coordinates, message):
    proxy = SimpleNamespace(
        n_obs=plot_adata.n_obs,
        obs=plot_adata.obs,
        obs_names=plot_adata.obs_names,
        obsm={"broken": coordinates(science.np, plot_adata.n_obs)},
        uns={},
        isbacked=False,
    )
    with pytest.raises((TypeError, ValueError), match=message):
        _plot_umap_impl(proxy, embedding_key="broken", color="")


def test_umap_rejects_coordinate_index_misalignment_and_color_identity_collisions(plot_adata, science):
    proxy = SimpleNamespace(
        n_obs=plot_adata.n_obs,
        obs=plot_adata.obs,
        obs_names=plot_adata.obs_names,
        obsm={
            "misaligned": science.pd.DataFrame(
                science.np.ones((plot_adata.n_obs, 2)), index=list(reversed(plot_adata.obs_names))
            )
        },
        uns={},
        isbacked=False,
    )
    with pytest.raises(ValueError, match="not exactly aligned"):
        _plot_umap_impl(proxy, embedding_key="misaligned", color="")

    plot_adata.obs["collision"] = science.pd.Series(
        [1, "1"] * 6,
        index=plot_adata.obs_names,
        dtype=object,
    )
    with pytest.raises(ValueError, match="collapse after string conversion"):
        _plot_umap_impl(plot_adata, embedding_key="custom_umap", color="collision", color_mode="categorical")


def test_umap_uses_first_two_of_three_dimensions_and_discloses_uncolored_style(plot_adata, science):
    plot_adata.obsm["X_umap_3d"] = science.np.column_stack(
        [
            science.np.linspace(-1.0, 1.0, plot_adata.n_obs),
            science.np.linspace(2.0, -2.0, plot_adata.n_obs),
            science.np.linspace(100.0, 200.0, plot_adata.n_obs),
        ]
    )
    plotted, summary, code = umap_plot_owned(
        plot_adata,
        embedding_key="X_umap_3d",
        color="",
    )
    details = summary.summary["key_results"]
    assert details["available_coordinate_dimensions"] == 3
    assert details["used_coordinate_dimensions"] == [1, 2]
    assert details["color"]["fixed_color"] == "#246bfe"
    assert details["rendering"]["uncolored_point_color"] == "#246bfe"
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_umap"](plot_adata) == plotted.png


class _FakeSettings:
    def __init__(self):
        self.autoshow = True


class _FakeBasePlot:
    def __init__(self, figure_class):
        self._figure_class = figure_class
        self.fig = None

    def make_figure(self):
        self.fig = self._figure_class()
        self.fig.subplots()


class _FakeScanpy:
    def __init__(self, pyplot, *, fail=False, malformed=False, advance_rng=False):
        from matplotlib.figure import Figure

        self.settings = _FakeSettings()
        self.calls = []
        self.dendrogram_calls = []
        self._pyplot = pyplot
        self._figure_class = Figure
        self._fail = fail
        self._malformed = malformed
        self._advance_rng = advance_rng
        self.pl = SimpleNamespace(
            dotplot=self._plotter("dotplot"),
            matrixplot=self._plotter("matrixplot"),
            tracksplot=self._plotter("tracksplot"),
            violin=self._plotter("violin"),
        )
        self.tl = SimpleNamespace(dendrogram=self._dendrogram)

    def _plotter(self, name):
        def call(adata, **kwargs):
            adata.uns["private_plot_mutation"] = name
            self.calls.append((name, adata, kwargs))
            if self._advance_rng:
                import numpy as np

                np.random.random()
            if self._fail:
                raise RuntimeError("synthetic backend failure")
            if self._malformed:
                return object()
            if name in {"dotplot", "matrixplot"}:
                return _FakeBasePlot(self._figure_class)
            figure = self._figure_class()
            axis = figure.subplots()
            return {"main": axis} if name == "tracksplot" else axis

        return call

    def _dendrogram(self, adata, groupby, **kwargs):
        matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
        self.dendrogram_calls.append((adata, groupby, matrix, kwargs))
        categories = list(adata.obs[groupby].cat.categories)
        return {"categories_ordered": list(reversed(categories))}


def _fake_science(science, fake_scanpy):
    return SimpleNamespace(
        ad=science.ad,
        np=science.np,
        pd=science.pd,
        sc=fake_scanpy,
        sparse=science.sparse,
    )


@pytest.mark.parametrize(
    ("plot", "name", "required", "forbidden"),
    [
        (
            {"plot": "dotplot", "standard_scale": "group", "expression_cutoff": 1.0, "mean_only_expressed": True},
            "dotplot",
            {"standard_scale", "expression_cutoff", "mean_only_expressed", "categories_order", "return_fig"},
            {"density_norm", "stripplot", "order"},
        ),
        (
            {"plot": "matrixplot", "standard_scale": "none"},
            "matrixplot",
            {"standard_scale", "categories_order", "return_fig"},
            {"expression_cutoff", "mean_only_expressed", "density_norm"},
        ),
        (
            {"plot": "tracksplot"},
            "tracksplot",
            {"dendrogram", "use_raw", "layer"},
            {"standard_scale", "categories_order", "density_norm"},
        ),
        (
            {"plot": "violin", "density_norm": "area", "show_cells": True, "y_scale": "linear"},
            "violin",
            {"density_norm", "stripplot", "jitter", "order", "multi_panel"},
            {"standard_scale", "dendrogram", "categories_order"},
        ),
    ],
)
def test_marker_plot_fake_backend_receives_only_effective_branch_parameters(
    plot_adata, science, plot, name, required, forbidden
):
    import matplotlib.pyplot as plt

    fake = _FakeScanpy(plt)
    snapshot = plot_adata.copy()
    before_figures = set(plt.get_fignums())
    png, details = _marker_expression_plot_impl(
        plot_adata,
        genes="G1,G2",
        groupby="cluster",
        plot=plot,
        source_kind="layer",
        layer_name="log1p_norm",
        group_order="observed",
        _science=_fake_science(science, fake),
    )
    assert png.startswith(PNG_SIGNATURE)
    assert details["plot"] == plot
    assert fake.calls[0][0] == name
    kwargs = fake.calls[0][2]
    assert required <= set(kwargs)
    assert not (forbidden & set(kwargs))
    assert kwargs["layer"] is None
    assert kwargs["use_raw"] is False
    assert kwargs["show"] is False and kwargs["save"] is False
    assert list(fake.calls[0][1].obs["cluster"].cat.categories) == ["C", "A", "B"]
    assert fake.calls[0][1].shape == (plot_adata.n_obs, 2)
    assert list(fake.calls[0][1].var_names) == ["G1", "G2"]
    assert not [key for key in fake.calls[0][1].layers if key is not None]
    assert not list(fake.calls[0][1].obsm)
    assert fake.calls[0][1].raw is None
    assert details["resource_estimate"]["estimated_peak_bytes"] <= details["resource_estimate"]["limit_bytes"]
    assert fake.settings.autoshow is True
    assert set(plt.get_fignums()) == before_figures
    _assert_adata_unchanged(plot_adata, snapshot, science)


def test_marker_dendrogram_uses_exact_selected_layer_panel_and_resolved_order(plot_adata, science):
    import matplotlib.pyplot as plt

    fake = _FakeScanpy(plt)
    png, details = _marker_expression_plot_impl(
        plot_adata,
        genes="G1,G3",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source_kind="layer",
        layer_name="log1p_norm",
        group_order="dendrogram",
        _science=_fake_science(science, fake),
    )
    assert png.startswith(PNG_SIGNATURE)
    assert len(fake.dendrogram_calls) == 1
    _, _, dendrogram_panel, kwargs = fake.dendrogram_calls[0]
    expected = _dense(plot_adata.layers["log1p_norm"], science)[:, [0, 2]]
    science.np.testing.assert_allclose(dendrogram_panel, expected)
    assert kwargs["var_names"] == ["G1", "G3"]
    assert kwargs["use_raw"] is False
    assert kwargs["cor_method"] == "pearson"
    assert kwargs["linkage_method"] == "complete"
    assert kwargs["optimal_ordering"] is False
    assert details["resolved_group_order"] == ["B", "A", "C"]
    assert fake.calls[0][2]["dendrogram"] is False
    assert fake.calls[0][2]["categories_order"] == ["B", "A", "C"]


def test_marker_backend_failure_restores_scanpy_state_closes_figures_and_preserves_input(plot_adata, science):
    import matplotlib.pyplot as plt

    fake = _FakeScanpy(plt, fail=True, advance_rng=True)
    snapshot = plot_adata.copy()
    before_figures = set(plt.get_fignums())
    before_rng = science.np.random.get_state()
    with pytest.raises(RuntimeError, match="synthetic backend failure"):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1,G2",
            groupby="cluster",
            plot={"plot": "tracksplot"},
            source_kind="layer",
            layer_name="log1p_norm",
            _science=_fake_science(science, fake),
        )
    assert fake.settings.autoshow is True
    assert set(plt.get_fignums()) == before_figures
    _assert_rng_state_equal(science.np.random.get_state(), before_rng, science)
    _assert_adata_unchanged(plot_adata, snapshot, science)


def test_marker_malformed_backend_and_inactive_parameters_fail_strictly(plot_adata, science):
    import matplotlib.pyplot as plt

    fake = _FakeScanpy(plt, malformed=True)
    with pytest.raises(RuntimeError, match="legacy Matplotlib plot contract"):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1",
            groupby="cluster",
            plot={"plot": "dotplot"},
            source_kind="layer",
            layer_name="log1p_norm",
            _science=_fake_science(science, fake),
        )

    with pytest.raises(TypeError, match="random_seed must be an integer"):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1",
            groupby="cluster",
            plot={"plot": "tracksplot"},
            source_kind="layer",
            layer_name="log1p_norm",
            random_seed=True,
            _science=_fake_science(science, fake),
        )
    with pytest.raises(ValueError, match=r"between 0 and 2\^31 - 1"):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1",
            groupby="cluster",
            plot={"plot": "tracksplot"},
            source_kind="layer",
            layer_name="log1p_norm",
            random_seed=-1,
            _science=_fake_science(science, fake),
        )
    with pytest.raises(ValueError, match="inactive or unknown"):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1",
            groupby="cluster",
            plot={"plot": "tracksplot", "standard_scale": "var"},
            source_kind="layer",
            layer_name="log1p_norm",
            _science=_fake_science(science, fake),
        )


@pytest.mark.parametrize("source_kind", ["X", "layer", "raw"])
def test_marker_plot_builds_only_a_minimal_selected_source_anndata(
    plot_adata, science, monkeypatch, source_kind
):
    import matplotlib.pyplot as plt

    if source_kind == "X":
        plot_adata.X = _dense(plot_adata.X, science)
        layer_name = None
        source_matrix = plot_adata.X
        source_names = plot_adata.var_names
    elif source_kind == "layer":
        layer_name = "log1p_norm"
        source_matrix = plot_adata.layers[layer_name]
        source_names = plot_adata.var_names
    else:
        layer_name = None
        source_matrix = plot_adata.raw.X
        source_names = plot_adata.raw.var_names
    expected = _dense(source_matrix, science)[:, [source_names.get_loc("G1"), source_names.get_loc("G2")]]
    plot_adata.layers["unrelated_sentinel"] = science.sparse.csr_matrix(science.np.ones(plot_adata.shape))
    plot_adata.obsm["unrelated_sentinel"] = science.np.ones((plot_adata.n_obs, 64))
    plot_adata.uns["unrelated_sentinel"] = {"must_not_copy": True}
    original_copy = science.ad.AnnData.copy

    def guarded_copy(self, *args, **kwargs):
        if self is plot_adata:
            raise AssertionError("caller-owned AnnData was deep-copied")
        return original_copy(self, *args, **kwargs)

    monkeypatch.setattr(science.ad.AnnData, "copy", guarded_copy)
    fake = _FakeScanpy(plt)
    _, details = _marker_expression_plot_impl(
        plot_adata,
        genes="G1,G2",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source_kind=source_kind,
        layer_name=layer_name,
        _science=_fake_science(science, fake),
    )
    work = fake.calls[0][1]
    science.np.testing.assert_allclose(_dense(work.X, science), expected)
    assert list(work.obs) == ["cluster"]
    assert list(work.var_names) == ["G1", "G2"]
    assert not [key for key in work.layers if key is not None]
    assert not list(work.obsm) and work.raw is None
    assert work.uns == {"private_plot_mutation": "matrixplot"}
    resource = details["resource_estimate"]
    assert resource["dense_panel_bytes"] > 0
    assert resource["plotting_matrix_copy_bytes"] == resource["dense_panel_bytes"]
    assert resource["estimated_peak_bytes"] <= resource["limit_bytes"]


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda adata, science: setattr(
                adata, "obs", adata.obs.assign(cluster=science.np.arange(adata.n_obs))
            ),
            "categorical",
        ),
        (lambda adata, science: setattr(adata, "obs", adata.obs.assign(cluster=["A"] * 11 + [None])), "missing"),
        (lambda adata, science: adata.layers.__setitem__("log1p_norm", science.np.full(adata.shape, science.np.nan)), "non-finite"),
    ],
)
def test_marker_strict_group_and_expression_validation(plot_adata, science, mutate, message):
    import matplotlib.pyplot as plt

    mutate(plot_adata, science)
    fake = _FakeScanpy(plt)
    with pytest.raises((TypeError, ValueError), match=message):
        _marker_expression_plot_impl(
            plot_adata,
            genes="G1,G2",
            groupby="cluster",
            plot={"plot": "matrixplot"},
            source_kind="layer",
            layer_name="log1p_norm",
            _science=_fake_science(science, fake),
        )


def test_marker_dotplot_statistics_use_strict_cutoff_and_null_undefined_means(plot_adata, science):
    import matplotlib.pyplot as plt

    plot_adata.layers["log1p_norm"] = science.sparse.csr_matrix(
        science.np.column_stack(
            [science.np.asarray([1.0, 2.0, 2.0, 3.0] * 3), science.np.zeros(12), science.np.ones((12, 2))]
        )
    )
    fake = _FakeScanpy(plt)
    _, details = _marker_expression_plot_impl(
        plot_adata,
        genes="G1,G2",
        groupby="cluster",
        plot={"plot": "dotplot", "standard_scale": "none", "expression_cutoff": 2.0, "mean_only_expressed": True},
        source_kind="layer",
        layer_name="log1p_norm",
        _science=_fake_science(science, fake),
    )
    rows = {(row["group"], row["gene"]): row for row in details["plotted_statistics"]}
    assert rows[("C", "G1")]["fraction_above_cutoff"] == 0.25
    assert rows[("C", "G1")]["mean_expression"] == 3.0
    assert rows[("C", "G2")]["mean_expression"] is None
    assert any("reported as null" in warning for warning in details["warnings"])


def test_marker_raw_expression_state_requires_evidence(plot_adata, science):
    _, unknown, _ = marker_expression_plot_owned(
        plot_adata,
        genes="G1,G2",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source={"source": "raw"},
    )
    assert unknown.summary["key_results"]["expression_state"] == "unknown"
    assert unknown.summary["key_results"]["expression_evidence"] is None
    assert any("without assuming counts" in warning for warning in unknown.summary["warnings"])

    snapshot = plot_adata.copy()
    snapshot.layers["counts"] = snapshot.raw.X[:, : snapshot.n_vars].copy()
    snapshot.uns["openbio_singlecell"]["analysis_history"]["000001"] = {
        "operation": "snapshot_expression",
        "parameters": {"destination": "layer_and_raw", "layer_name": "counts", "source": "X"},
    }
    _, counts, counts_code = marker_expression_plot_owned(
        snapshot,
        genes="G1,RAW_ONLY",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source={"source": "raw"},
    )
    assert counts.summary["key_results"]["expression_state"] == "counts"
    assert "Snapshot Expression" in counts.summary["key_results"]["expression_evidence"]
    assert any("proven count-scale" in warning for warning in counts.summary["warnings"])
    namespace = {}
    exec(counts_code, namespace)
    assert namespace["plot_marker_expression"](snapshot).startswith(PNG_SIGNATURE)

    external_logged = plot_adata.copy()
    logged = science.np.log1p(_dense(external_logged.raw.X, science)[:, : external_logged.n_vars])
    external_logged.X = science.sparse.csr_matrix(logged)
    external_logged.raw = science.ad.AnnData(
        science.sparse.csr_matrix(logged),
        obs=external_logged.obs.copy(),
        var=external_logged.var.copy(),
    )
    external_logged.uns["log1p"] = {"base": None}
    _, logged_summary, _ = marker_expression_plot_owned(
        external_logged,
        genes="G1,G2",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source={"source": "raw"},
    )
    assert logged_summary.summary["key_results"]["expression_state"] == "logged"
    assert "Raw values equal proven current X" in logged_summary.summary["key_results"]["expression_evidence"]
    assert not any("count-scale" in warning for warning in logged_summary.summary["warnings"])


def test_marker_violin_seed_is_deterministic_and_restores_numpy_rng(plot_adata, science):
    outer_state = science.np.random.get_state()
    try:
        science.np.random.seed(20260828)
        before = science.np.random.get_state()
        first_plot, first_summary, first_code = marker_expression_plot_owned(
            plot_adata,
            genes="G1,G2",
            groupby="cluster",
            plot={"plot": "violin", "density_norm": "width", "show_cells": True, "y_scale": "linear"},
            source={"source": "layer", "layer_name": "log1p_norm"},
            group_order="observed",
            random_seed=17,
        )
        _assert_rng_state_equal(science.np.random.get_state(), before, science)
        second_plot, second_summary, second_code = marker_expression_plot_owned(
            plot_adata,
            genes="G1,G2",
            groupby="cluster",
            plot={"plot": "violin", "density_norm": "width", "show_cells": True, "y_scale": "linear"},
            source={"source": "layer", "layer_name": "log1p_norm"},
            group_order="observed",
            random_seed=17,
        )
        _assert_rng_state_equal(science.np.random.get_state(), before, science)
        assert first_plot.png == second_plot.png
        assert first_code == second_code
        assert first_summary.summary["parameters"]["random_seed"] == 17
        assert first_summary.summary["key_results"]["rendering"]["numpy_global_rng_state_restored"] is True
        assert "scipy" in first_summary.summary["software_versions"]
        assert "seaborn" in first_summary.summary["software_versions"]
        assert any("seaborn" in reference["citation"] for reference in first_summary.summary["references"])
        namespace = {}
        exec(first_code, namespace)
        generated = namespace["plot_marker_expression"](plot_adata)
        _assert_rng_state_equal(science.np.random.get_state(), before, science)
        assert generated == first_plot.png
        assert second_summary.summary["key_results"]["rendering"]["random_seed"] == 17
    finally:
        science.np.random.set_state(outer_state)


def test_marker_real_runtime_summary_code_and_immutability(plot_adata, science):
    snapshot = plot_adata.copy()
    plotted, summary, code = marker_expression_plot_owned(
        plot_adata,
        genes="G1,G2,G1",
        groupby="cluster",
        plot={"plot": "matrixplot", "standard_scale": "none"},
        source={"source": "layer", "layer_name": "log1p_norm"},
        group_order="observed",
    )
    assert plotted.png.startswith(PNG_SIGNATURE)
    assert summary.summary["key_results"]["genes"] == ["G1", "G2"]
    assert summary.summary["key_results"]["resolved_group_order"] == ["C", "A", "B"]
    assert summary.summary["key_results"]["selected_panel_range"][0] >= 0.0
    assert len(summary.summary["key_results"]["plotted_statistics"]) == 6
    assert "scipy" in summary.summary["software_versions"]
    assert "seaborn" not in summary.summary["software_versions"]
    assert not any("seaborn" in reference["citation"] for reference in summary.summary["references"])
    json.dumps(summary.summary, allow_nan=False)
    compile(code, "<marker-expression-code>", "exec")
    namespace = {}
    exec(code, namespace)
    generated = namespace["plot_marker_expression"](plot_adata)
    assert generated.startswith(PNG_SIGNATURE)
    assert generated == plotted.png
    _assert_adata_unchanged(plot_adata, snapshot, science)
