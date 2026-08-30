from __future__ import annotations

import json
from importlib.metadata import version
from pathlib import Path

import pytest

from openbio_singlecell.nodes_correction import (
    OpenBioSingleCellFilterDoublets,
    OpenBioSingleCellMarkMADOutliers,
    OpenBioSingleCellScrublet,
)
from openbio_singlecell.operations_correction import (
    filter_doublets_owned,
    mark_mad_outliers_owned,
    scrublet_owned,
)

_MAD_SCHEMA = OpenBioSingleCellMarkMADOutliers
_SCRUBLET_SCHEMA = OpenBioSingleCellScrublet
_FILTER_SCHEMA = OpenBioSingleCellFilterDoublets


def _adata_with_obs(science, values, *, samples=None):
    values = list(values)
    obs = science.pd.DataFrame(
        {"metric": values},
        index=[f"cell_{index}" for index in range(len(values))],
    )
    if samples is not None:
        obs["sample"] = list(samples)
    return science.ad.AnnData(
        science.sparse.csr_matrix(science.np.ones((len(values), 2), dtype=float)),
        obs=obs,
        var=science.pd.DataFrame(index=["gene_0", "gene_1"]),
    )


def _scrublet_adata(science):
    counts = science.np.asarray(
        [
            [2, 0, 1, 0],
            [3, 1, 0, 1],
            [2, 1, 1, 0],
            [0, 3, 1, 1],
            [1, 2, 0, 2],
            [0, 2, 2, 1],
        ],
        dtype=float,
    )
    obs = science.pd.DataFrame(
        {"sample": ["sample_a"] * 3 + ["sample_b"] * 3},
        index=[f"cell_{index}" for index in range(6)],
    )
    value = science.ad.AnnData(
        science.np.log1p(counts),
        obs=obs,
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(4)]),
    )
    value.layers["counts"] = science.sparse.csr_matrix(counts)
    return value


def _install_fake_scrublet(monkeypatch, science, *, include_threshold=True):
    calls = []

    def fake_scrublet(work, **parameters):
        calls.append((work.copy(), dict(parameters)))
        scores = science.np.asarray([0.1, 0.4, 0.2, 0.3, 0.8, 0.1], dtype=float)
        threshold = parameters["threshold"] if parameters["threshold"] is not None else 0.35
        work.obs["doublet_score"] = scores
        work.obs["predicted_doublet"] = scores > threshold
        group_key = parameters["batch_key"]
        if group_key:
            batches = {}
            for label in science.pd.unique(work.obs[group_key]):
                payload = {
                    "parameters": {
                        "expected_doublet_rate": parameters["expected_doublet_rate"],
                        "sim_doublet_ratio": parameters["sim_doublet_ratio"],
                        "n_neighbors": parameters["n_neighbors"] or 2,
                        "random_state": parameters["random_state"],
                    }
                }
                if include_threshold:
                    payload["threshold"] = threshold
                batches[label] = payload
            work.uns["scrublet"] = {"batches": batches, "batched_by": group_key}
        else:
            payload = {
                "parameters": {
                    "expected_doublet_rate": parameters["expected_doublet_rate"],
                    "sim_doublet_ratio": parameters["sim_doublet_ratio"],
                    "n_neighbors": parameters["n_neighbors"] or 2,
                    "random_state": parameters["random_state"],
                }
            }
            if include_threshold:
                payload["threshold"] = threshold
            work.uns["scrublet"] = payload

    monkeypatch.setattr(science.sc.pp, "scrublet", fake_scrublet)
    return calls


def test_mad_normal_scaling_uses_scipy_normal_convention_and_code_matches(science):
    value = _adata_with_obs(science, [0.0, 1.0, 2.0])

    raw, raw_report, _ = mark_mad_outliers_owned(
        value,
        metrics="metric",
        batch_key="",
        nmads=0.8,
        direction="upper",
        output_column="raw_outlier",
        scale_mad=False,
    )
    scaled, scaled_report, scaled_code = mark_mad_outliers_owned(
        value,
        metrics="metric",
        batch_key="",
        nmads=0.8,
        direction="upper",
        output_column="scaled_outlier",
        scale_mad=True,
    )

    assert raw.obs["raw_outlier"].tolist() == [False, False, True]
    assert scaled.obs["scaled_outlier"].tolist() == [False, False, False]
    raw_stats = raw_report.summary["key_results"]["metric_statistics"][0]["groups"][0]
    scaled_stats = scaled_report.summary["key_results"]["metric_statistics"][0]["groups"][0]
    assert raw_stats["mad"] == pytest.approx(1.0)
    assert scaled_stats["mad"] == pytest.approx(1.4826022185)
    assert scaled_stats["upper_threshold"] > raw_stats["upper_threshold"]
    namespace = {}
    exec(scaled_code, namespace)
    equivalent = namespace["mark_mad_outliers"](value)
    assert equivalent.obs["scaled_outlier"].equals(scaled.obs["scaled_outlier"])


def test_mad_is_grouped_by_sample_and_reports_thresholds(science):
    value = _adata_with_obs(
        science,
        [0.0, 1.0, 10.0, 100.0, 101.0, 110.0],
        samples=["a", "a", "a", "b", "b", "b"],
    )

    output, report, code = mark_mad_outliers_owned(
        value,
        metrics="metric",
        batch_key="sample",
        nmads=3.0,
        direction="upper",
    )

    assert output.obs["outlier"].tolist() == [False, False, True, False, False, True]
    groups = report.summary["key_results"]["metric_statistics"][0]["groups"]
    assert [(group["group"], group["median"], group["flagged"]) for group in groups] == [
        ("a", 1.0, 1),
        ("b", 101.0, 1),
    ]
    assert report.summary["key_results"]["flagged_cells"] == 2
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<mad-code>", "exec")


def test_mad_validates_group_labels_metrics_and_small_groups(science):
    finite = _adata_with_obs(science, [1.0, 2.0, 3.0], samples=["a", "a", "a"])
    with pytest.raises(ValueError, match="greater than zero"):
        mark_mad_outliers_owned(finite, metrics="metric", nmads=science.np.nan)

    missing_group = _adata_with_obs(science, [1.0, 2.0, 3.0], samples=["a", None, "a"])
    with pytest.raises(ValueError, match="missing labels"):
        mark_mad_outliers_owned(missing_group, metrics="metric")

    colliding_groups = _adata_with_obs(science, range(6), samples=[1, 1, 1, "1", "1", "1"])
    with pytest.raises(ValueError, match="grouped MAD reports"):
        mark_mad_outliers_owned(colliding_groups, metrics="metric")
    _, _, code = mark_mad_outliers_owned(finite, metrics="metric")
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="grouped MAD reports"):
        namespace["mark_mad_outliers"](colliding_groups)

    nonnumeric = _adata_with_obs(science, [1.0, 2.0, 3.0], samples=["a", "a", "a"])
    nonnumeric.obs["metric"] = ["low", "middle", "high"]
    with pytest.raises(ValueError, match="must be numeric"):
        mark_mad_outliers_owned(nonnumeric, metrics="metric")

    small = _adata_with_obs(science, [1.0, science.np.nan, 9.0], samples=["a", "a", "a"])
    output, report, _ = mark_mad_outliers_owned(
        small,
        metrics="metric",
        minimum_group_size=3,
    )
    assert not output.obs["outlier"].any()
    assert any("skipped" in warning for warning in report.warnings)
    assert report.summary["key_results"]["metric_statistics"][0]["missing"] == 1

    singleton = _adata_with_obs(science, [1.0, 2.0, 3.0], samples=["a", "b", "c"])
    singleton_output, singleton_report, singleton_code = mark_mad_outliers_owned(
        singleton,
        metrics="metric",
        minimum_group_size=1,
    )
    assert not singleton_output.obs["outlier"].any()
    assert any("zero MAD" in warning for warning in singleton_report.warnings)
    namespace = {}
    exec(singleton_code, namespace)
    assert namespace["mark_mad_outliers"](singleton).obs["outlier"].equals(
        singleton_output.obs["outlier"]
    )


def test_scrublet_uses_explicit_count_source_preserves_cell_ids_and_reports_groups(monkeypatch, science):
    value = _scrublet_adata(science)
    original_x = science.np.asarray(value.X).copy()
    original_counts = value.layers["counts"].copy()
    original_names = value.obs_names.copy()
    calls = _install_fake_scrublet(monkeypatch, science)

    output, report, code = scrublet_owned(
        value,
        batch_key="sample",
        random_seed=7,
        expected_doublet_rate=0.08,
        threshold_mode="manual",
        threshold=0.25,
        sim_doublet_ratio=3.0,
        n_prin_comps=2,
        n_neighbors=2,
        source={"source": "layer", "source_layer": "counts"},
    )

    science.np.testing.assert_array_equal(calls[0][0].X.toarray(), original_counts.toarray())
    assert calls[0][1]["threshold"] == 0.25
    assert calls[0][1]["batch_key"] == "sample"
    assert output.obs_names.equals(original_names)
    science.np.testing.assert_array_equal(value.X, original_x)
    science.np.testing.assert_array_equal(value.layers["counts"].toarray(), original_counts.toarray())
    assert output.obs["predicted_doublet"].dtype == bool
    assert output.obs["predicted_doublet"].tolist() == [False, True, False, True, True, False]
    assert report.summary["key_results"]["predicted_doublets"] == 3
    assert report.summary["key_results"]["effective_thresholds"] == {"sample_a": 0.25, "sample_b": 0.25}
    assert report.summary["parameters"]["source"] == "layer"
    assert report.summary["software_versions"]["scikit-image"] != "not-installed"
    json.dumps(report.summary, allow_nan=False)

    namespace = {}
    exec(code, namespace)
    equivalent = namespace["run_scrublet"](value)
    science.np.testing.assert_allclose(equivalent.obs["doublet_score"], output.obs["doublet_score"])
    assert equivalent.obs["predicted_doublet"].equals(output.obs["predicted_doublet"])
    assert equivalent.obs_names.equals(original_names)
    assert len(calls) == 2


def test_scrublet_reports_overwritten_result_fields(monkeypatch, science):
    value = _scrublet_adata(science)
    value.obs["doublet_score"] = 0.0
    value.obs["predicted_doublet"] = False
    value.uns["scrublet"] = {"old": True}
    _install_fake_scrublet(monkeypatch, science)

    _, report, _ = scrublet_owned(
        value,
        threshold_mode="manual",
        threshold=0.25,
        n_prin_comps=2,
        source={"source": "layer", "source_layer": "counts"},
    )

    assert any("were overwritten" in warning for warning in report.warnings)


def test_scrublet_warns_on_advisory_inputs_and_rejects_computational_contracts(monkeypatch, science):
    value = _scrublet_adata(science)
    _install_fake_scrublet(monkeypatch, science)

    duplicate = value.copy()
    duplicate.obs_names = ["duplicate", "duplicate", "c2", "c3", "c4", "c5"]
    duplicate_output, duplicate_report, duplicate_code = scrublet_owned(
        duplicate,
        n_prin_comps=2,
        source={"source": "layer", "source_layer": "counts"},
    )
    assert duplicate_output.obs_names.tolist() == duplicate.obs_names.tolist()
    assert duplicate_report.summary["key_results"]["duplicate_observation_identifiers"] is True
    assert any("not unique" in warning for warning in duplicate_report.warnings)
    namespace = {}
    exec(duplicate_code, namespace)
    assert namespace["run_scrublet"](duplicate).obs_names.tolist() == duplicate.obs_names.tolist()

    missing_sample = value.copy()
    missing_sample.obs.loc[missing_sample.obs_names[0], "sample"] = None
    with pytest.raises(ValueError, match="missing labels"):
        scrublet_owned(
            missing_sample,
            source={"source": "layer", "source_layer": "counts"},
        )

    noninteger = value.copy()
    noninteger.layers["counts"] = noninteger.layers["counts"].astype(float) * 0.5
    _, noninteger_report, _ = scrublet_owned(
        noninteger,
        n_prin_comps=2,
        source={"source": "layer", "source_layer": "counts"},
    )
    assert any("non-integer" in warning for warning in noninteger_report.warnings)

    one_component, _, _ = scrublet_owned(
        value,
        n_prin_comps=1,
        source={"source": "layer", "source_layer": "counts"},
    )
    assert one_component.n_obs == value.n_obs

    automatic, automatic_report, _ = scrublet_owned(
        value,
        threshold_mode="automatic",
        threshold=science.np.nan,
        n_prin_comps=2,
        source={"source": "layer", "source_layer": "counts"},
    )
    assert automatic.n_obs == value.n_obs
    assert automatic_report.summary["parameters"]["threshold"] is None
    json.dumps(automatic_report.summary, allow_nan=False)

    extreme, extreme_report, extreme_code = scrublet_owned(
        value,
        threshold_mode="manual",
        threshold=1.5,
        n_prin_comps=2,
        source={"source": "layer", "source_layer": "counts"},
    )
    assert not extreme.obs["predicted_doublet"].any()
    assert any("outside the conventional" in warning for warning in extreme_report.warnings)
    namespace = {}
    exec(extreme_code, namespace)
    assert namespace["run_scrublet"](value).obs["predicted_doublet"].equals(
        extreme.obs["predicted_doublet"]
    )

    with pytest.raises(ValueError, match="finite and greater than zero"):
        scrublet_owned(
            value,
            sim_doublet_ratio=science.np.nan,
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )

    negative = value.copy()
    negative.layers["counts"] = negative.layers["counts"].copy()
    negative.layers["counts"].data[0] = -1
    with pytest.raises(ValueError, match="negative expression"):
        scrublet_owned(
            negative,
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )

    nonfinite = value.copy()
    nonfinite.layers["counts"] = nonfinite.layers["counts"].copy()
    nonfinite.layers["counts"].data[:] = science.np.nan
    with pytest.raises(ValueError, match="non-finite expression"):
        scrublet_owned(
            nonfinite,
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )
    namespace = {}
    exec(extreme_code, namespace)
    with pytest.raises(ValueError, match="non-finite expression"):
        namespace["run_scrublet"](nonfinite)

    all_negative = value.copy()
    all_negative.layers["counts"] = all_negative.layers["counts"].copy()
    all_negative.layers["counts"].data[:] = -science.np.abs(all_negative.layers["counts"].data)
    with pytest.raises(ValueError, match="negative expression"):
        scrublet_owned(
            all_negative,
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )
    with pytest.raises(ValueError, match="negative expression"):
        namespace["run_scrublet"](all_negative)

    with pytest.raises(ValueError, match="requires at least 4 cells"):
        scrublet_owned(
            value,
            n_prin_comps=3,
            source={"source": "layer", "source_layer": "counts"},
        )

    colliding_groups = value.copy()
    colliding_groups.obs["sample"] = [1, 1, 1, "1", "1", "1"]
    with pytest.raises(ValueError, match="collide"):
        scrublet_owned(
            colliding_groups,
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )
    with pytest.raises(ValueError, match="collide"):
        namespace["run_scrublet"](colliding_groups)


def test_scrublet_rejects_missing_automatic_threshold(monkeypatch, science):
    value = _scrublet_adata(science)
    _install_fake_scrublet(monkeypatch, science, include_threshold=False)

    with pytest.raises(RuntimeError, match="did not determine a usable threshold"):
        scrublet_owned(
            value,
            batch_key="sample",
            threshold_mode="automatic",
            n_prin_comps=2,
            source={"source": "layer", "source_layer": "counts"},
        )


def test_scrublet_runtime_dependency_is_declared():
    root = Path(__file__).resolve().parents[1]

    assert version("scikit-image")
    assert "scanpy[leiden,scrublet]" in (root / "pyproject.toml").read_text(encoding="utf-8")
    assert "scanpy[leiden,scrublet]" in (root / "requirements.txt").read_text(encoding="utf-8")


def test_scrublet_scanpy_112_smoke(science):
    rng = science.np.random.default_rng(42)
    counts = rng.poisson(1.3, size=(60, 50)).astype(float)
    value = science.ad.AnnData(
        science.sparse.csr_matrix(counts),
        obs=science.pd.DataFrame(
            {"sample": ["sample_a"] * 30 + ["sample_b"] * 30},
            index=[f"cell_{index}" for index in range(60)],
        ),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(50)]),
    )

    output, report, code = scrublet_owned(
        value,
        batch_key="sample",
        random_seed=0,
        threshold_mode="manual",
        threshold=0.1,
        n_prin_comps=5,
    )

    assert output.obs["doublet_score"].dtype.kind == "f"
    assert output.obs["predicted_doublet"].dtype == bool
    assert report.summary["key_results"]["effective_thresholds"] == {"sample_a": 0.1, "sample_b": 0.1}
    assert report.summary["software_versions"]["scanpy"] == "1.12.3"
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["run_scrublet"](value)
    science.np.testing.assert_allclose(equivalent.obs["doublet_score"], output.obs["doublet_score"])
    assert equivalent.obs["predicted_doublet"].equals(output.obs["predicted_doublet"])


def test_filter_doublets_requires_boolean_predictions_and_code_matches(science):
    value = _adata_with_obs(science, [1.0, 2.0, 3.0, 4.0])
    value.obs["predicted_doublet"] = science.pd.Series(
        [False, True, False, True],
        index=value.obs_names,
        dtype="boolean",
    )
    snapshot = value.copy()

    output, report, code = filter_doublets_owned(value)

    assert output.obs_names.tolist() == ["cell_0", "cell_2"]
    assert value.obs.equals(snapshot.obs)
    assert report.summary["key_results"]["predicted_doublets_removed"] == 2
    assert report.summary["key_results"]["retained_singlet_candidates"] == 2
    assert report.summary["references"][0]["url"].startswith("https://anndata.readthedocs.io/")
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    equivalent = namespace["filter_predicted_doublets"](value)
    assert equivalent.obs_names.equals(output.obs_names)


@pytest.mark.parametrize(
    "predictions",
    [
        ["False", "True", "False"],
        [0, 1, 0],
    ],
)
def test_filter_doublets_rejects_non_boolean_columns(science, predictions):
    value = _adata_with_obs(science, [1.0, 2.0, 3.0])
    value.obs["predicted_doublet"] = predictions

    with pytest.raises(ValueError, match="boolean dtype"):
        filter_doublets_owned(value)


def test_filter_doublets_rejects_missing_and_reports_all_true_predictions(science):
    value = _adata_with_obs(science, [1.0, 2.0, 3.0])
    value.obs["predicted_doublet"] = science.pd.Series(
        [False, science.pd.NA, True],
        index=value.obs_names,
        dtype="boolean",
    )
    with pytest.raises(ValueError, match="missing values"):
        filter_doublets_owned(value)

    value.obs["predicted_doublet"] = True
    empty, empty_report, empty_code = filter_doublets_owned(value)
    assert empty.n_obs == 0
    assert any("Every input cell" in warning for warning in empty_report.warnings)
    namespace = {}
    exec(empty_code, namespace)
    assert namespace["filter_predicted_doublets"](value).n_obs == 0

    value.obs["predicted_doublet"] = False
    output, report, code = filter_doublets_owned(value)
    assert output.n_obs == value.n_obs
    assert any("retains every input cell" in warning for warning in report.warnings)
    namespace = {}
    exec(code, namespace)
    assert namespace["filter_predicted_doublets"](value).n_obs == value.n_obs


def test_correction_node_schema_prefixes_are_current():
    assert [item.id for item in _MAD_SCHEMA.define_schema().inputs[:7]] == [
        "adata",
        "metrics",
        "batch_key",
        "nmads",
        "direction",
        "output_column",
        "scale_mad",
    ]
    assert [item.id for item in _SCRUBLET_SCHEMA.define_schema().inputs[:3]] == [
        "adata",
        "batch_key",
        "random_seed",
    ]
    assert [item.id for item in _FILTER_SCHEMA.define_schema().inputs] == [
        "adata",
        "prediction_column",
    ]
