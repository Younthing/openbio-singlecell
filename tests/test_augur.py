from __future__ import annotations

import json
import uuid
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

from openbio_singlecell import augur as augur_core
from openbio_singlecell.artifact_codecs import read_plot, read_table, write_anndata
from openbio_singlecell.augur import (
    AUGUR_ARTIFACT_TYPE,
    AUGUR_CLASSIFIERS,
    AUGUR_VIEWS,
    CROSS_VALIDATION_COLUMNS,
    FEATURE_IMPORTANCE_COLUMNS,
    PREDICTION_COLUMNS,
    PRIORITY_COLUMNS,
    AugurResult,
    run_augur_analysis,
    run_augur_artifact,
    select_augur_view,
    validate_augur_result,
)
from openbio_singlecell.augur_codec import AUGUR_CODEC, read_augur, write_augur
from openbio_singlecell.node_types import AugurResultType, PlotResultType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_population import (
    OpenBioSingleCellAugur,
    OpenBioSingleCellAugurPlot,
    OpenBioSingleCellAugurResults,
)
from openbio_singlecell.operations_population import augur as augur_operation
from openbio_singlecell.operations_population import (
    augur_legacy_owned,
    augur_plot,
    augur_plot_owned,
    augur_results,
    augur_results_legacy_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse


def _make_augur_adata(
    science,
    *,
    n_genes: int = 12,
    cells_per_sample_population: int = 3,
    include_weak_population: bool = True,
    include_other_condition: bool = True,
):
    sample_design = (
        ("control_1", "control", "batch_1"),
        ("control_2", "control", "batch_2"),
        ("treatment_1", "treatment", "batch_1"),
        ("treatment_2", "treatment", "batch_2"),
    )
    records = []
    for population in ("population_A", "population_B"):
        for sample, condition, batch in sample_design:
            records.extend(
                {
                    "sample": sample,
                    "population": population,
                    "condition": condition,
                    "technical_batch": batch,
                }
                for _ in range(cells_per_sample_population)
            )
    if include_weak_population:
        records.extend(
            {
                "sample": sample,
                "population": "population_weak",
                "condition": condition,
                "technical_batch": batch,
            }
            for sample, condition, batch in sample_design
        )
    if include_other_condition:
        records.extend(
            {
                "sample": "other_1",
                "population": "population_A",
                "condition": "other",
                "technical_batch": "batch_1",
            }
            for _ in range(2)
        )

    obs_names = [f"cell_{index:04d}" for index in range(len(records))]
    obs = science.pd.DataFrame(records, index=obs_names)
    var = science.pd.DataFrame(index=[f"gene_{index:03d}" for index in range(n_genes)])
    rng = science.np.random.default_rng(20260828)
    counts = rng.poisson(2.0, size=(len(records), n_genes)).astype(science.np.int64)
    counts[0, 0] = max(1, int(counts[0, 0]))
    full = science.ad.AnnData(X=counts, obs=obs, var=var)
    full.layers["counts"] = counts.copy()
    full.raw = full.copy()
    full.uns["openbio_singlecell"] = {
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {"destination": "layer_and_raw", "layer_name": "counts"},
            }
        }
    }
    current_features = max(1, n_genes // 2)
    return full[:, :current_features].copy()


@pytest.fixture
def augur_adata(science):
    return _make_augur_adata(science)


def _fake_pertpy(science, *, backend_mutator=None, native_evidence=True):
    calls = []

    class FakeAugur:
        def __init__(
            self,
            estimator,
            *,
            n_estimators=100,
            max_depth=None,
            max_features=2,
            penalty="l2",
            random_state=None,
        ):
            self.estimator = estimator
            calls.append(
                {
                    "operation": "constructor",
                    "estimator": estimator,
                    "n_estimators": n_estimators,
                    "max_depth": max_depth,
                    "max_features": max_features,
                    "penalty": penalty,
                    "random_state": random_state,
                }
            )

        def load(
            self,
            input,
            layer=None,
            meta=None,
            label_col="label",
            cell_type_col="cell_type",
            condition_label=None,
            treatment_label=None,
        ):
            calls.append(
                {
                    "operation": "load",
                    "layer": layer,
                    "meta": meta,
                    "label_col": label_col,
                    "cell_type_col": cell_type_col,
                    "condition_label": condition_label,
                    "treatment_label": treatment_label,
                    "conditions": input.obs[label_col].tolist(),
                    "samples": (
                        input.obs["sample"].tolist()
                        if "sample" in input.obs
                        else [None] * int(input.n_obs)
                    ),
                    "n_obs": int(input.n_obs),
                    "n_vars": int(input.n_vars),
                }
            )
            loaded = input.copy()
            loaded.obs["label"] = input.obs[label_col].tolist()
            loaded.obs["cell_type"] = input.obs[cell_type_col].tolist()
            encoding = {value: index for index, value in enumerate(sorted({condition_label, treatment_label}))}
            loaded.obs["y_"] = [encoding[value] for value in loaded.obs["label"]]
            return loaded

        def run_cross_validation(
            self,
            subsample,
            *,
            subsample_idx,
            folds,
            random_state,
            zero_division,
        ):
            class ScoreEstimator:
                classes_ = science.np.asarray([0, 1])

                @staticmethod
                def predict_proba(values):
                    scores = science.np.asarray(values)[:, 0]
                    return science.np.column_stack([1.0 - scores, scores])

            return {
                "estimator": [ScoreEstimator() for _ in range(folds)],
                "test_augur_score": [float(subsample.uns["expected_auc"])] * folds,
            }

        def predict(
            self,
            adata,
            n_subsamples=50,
            subsample_size=20,
            folds=3,
            min_cells=None,
            feature_perc=0.5,
            var_quantile=0.5,
            span=0.75,
            filter_negative_residuals=False,
            n_threads=1,
            augur_mode="default",
            select_variance_features=True,
            key_added="augurpy_results",
            random_state=None,
            zero_division=0,
        ):
            calls.append(
                {
                    "operation": "predict",
                    "n_subsamples": n_subsamples,
                    "subsample_size": subsample_size,
                    "folds": folds,
                    "min_cells": min_cells,
                    "feature_perc": feature_perc,
                    "var_quantile": var_quantile,
                    "span": span,
                    "filter_negative_residuals": filter_negative_residuals,
                    "n_threads": n_threads,
                    "augur_mode": augur_mode,
                    "select_variance_features": select_variance_features,
                    "key_added": key_added,
                    "random_state": random_state,
                    "zero_division": zero_division,
                }
            )
            populations = list(dict.fromkeys(adata.obs["cell_type"].tolist()))
            metrics = (
                "mean_augur_score",
                "mean_auc",
                "mean_accuracy",
                "mean_precision",
                "mean_f1",
                "mean_recall",
            )
            scores = {population: 1.0 - 0.5 * index for index, population in enumerate(populations)}
            summary_metrics = science.pd.DataFrame(
                {
                    population: [
                        scores[population],
                        scores[population],
                        scores[population] - 0.02,
                        scores[population] - 0.03,
                        scores[population] - 0.04,
                        scores[population] - 0.05,
                    ]
                    for population in populations
                },
                index=list(metrics),
            )
            full_rows = []
            feature_rows = []
            population_runs = {population: [] for population in populations}
            for population in populations:
                for subsample in range(n_subsamples):
                    prediction_rows = []
                    for fold in range(folds):
                        full_rows.append(
                            {
                                "idx": subsample,
                                "augur_score": scores[population],
                                "folds": fold,
                                "cell_type": population,
                            }
                        )
                        importances = (
                            (0.6, 0.4)
                            if self.estimator == "random_forest_classifier"
                            else (0.6, -0.4)
                        )
                        for gene, importance in zip(("gene_000", "gene_001"), importances, strict=True):
                            feature_rows.append(
                                {
                                    "genes": gene,
                                    "feature_importances": importance,
                                    "subsample_idx": subsample,
                                    "fold": fold,
                                    "cell_type": population,
                                }
                            )
                        labels = [0, 0, 1] if fold == 0 else [0, 1, 1]
                        fold_scores = (
                            [0.1, 0.2, 0.9]
                            if population == "population_A"
                            else ([0.1, 0.9, 0.5] if fold == 0 else [0.5, 0.1, 0.9])
                        )
                        prediction_rows.extend(
                            {
                                "fold": fold,
                                "observation": f"{population}_{subsample}_{fold}_{index}",
                                "true_label": label,
                                "prediction_score": score,
                            }
                            for index, (label, score) in enumerate(
                                zip(labels, fold_scores, strict=True)
                            )
                        )
                    run = {"openbio_prediction_evidence": prediction_rows}
                    if not native_evidence:
                        run = self.run_cross_validation(
                            science.ad.AnnData(
                                X=science.np.asarray(
                                    [
                                        [0.1],
                                        [0.2],
                                        [0.3],
                                        [0.7],
                                        [0.8],
                                        [0.9],
                                    ]
                                    if population == "population_A"
                                    else [[0.5]] * 6
                                ),
                                obs=science.pd.DataFrame(
                                    {
                                        "y_": [0, 0, 0, 1, 1, 1],
                                        "label": ["control"] * 3 + ["treatment"] * 3,
                                    },
                                    index=[f"{population}_{subsample}_{index}" for index in range(6)],
                                ),
                                uns={"expected_auc": scores[population]},
                            ),
                            subsample_idx=subsample,
                            folds=folds,
                            random_state=random_state,
                            zero_division=zero_division,
                        )
                    population_runs[population].append(run)
            results = {
                "summary_metrics": summary_metrics,
                "full_results": science.pd.DataFrame(full_rows),
                "feature_importances": science.pd.DataFrame(feature_rows),
                **population_runs,
            }
            if backend_mutator is not None:
                backend_mutator(results)
            return adata.copy(), results

    return SimpleNamespace(__version__="1.3.0", tl=SimpleNamespace(Augur=FakeAugur), calls=calls)


def _run(adata, science, *, pertpy_module=None, _runner=run_augur_analysis, **overrides):
    parameters = {
        "sample_key": "sample",
        "population_key": "population",
        "condition_key": "condition",
        "control": "control",
        "treatment": "treatment",
        "classifier": "random_forest_classifier",
        "source_kind": "raw",
        "layer_name": None,
        "annotation_status": "provisional",
        "technical_batch_key": "technical_batch",
        "n_subsamples": 2,
        "subsample_size": 3,
        "folds": 2,
        "n_threads": 1,
        "random_seed": 17,
        "max_result_rows": 10_000,
        "max_result_mib": 10.0,
        "pertpy_module": pertpy_module or _fake_pertpy(science),
    }
    parameters.update(overrides)
    return _runner(adata, **parameters)


def _replace_obs_column(adata, column, values):
    obs = adata.obs.copy()
    obs[column] = values
    adata.obs = obs


def test_augur_portable_codec_round_trips_the_validated_table_family(tmp_path: Path, augur_adata, science):
    result = _run(augur_adata, science)
    tables, summary, metadata = validate_augur_result(result)
    root = tmp_path / "augur"
    root.mkdir()

    write_augur(root, tables, summary, metadata)
    restored_tables, restored_summary, restored_metadata = read_augur(root)

    assert AUGUR_CODEC == "augur-jsonl-v2"
    assert restored_summary == summary
    assert restored_metadata == metadata
    for view in AUGUR_VIEWS:
        science.pd.testing.assert_frame_equal(restored_tables[view], tables[view])


def test_augur_codec_rejects_tampered_prediction_evidence(tmp_path: Path, augur_adata, science):
    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))
    root = tmp_path / "augur"
    root.mkdir()
    write_augur(root, tables, summary, metadata)
    data_path = root / "tables" / "predictions" / "data.jsonl"
    rows = data_path.read_text(encoding="utf-8").splitlines()
    first = json.loads(rows[0])
    first[6] = 1.0
    rows[0] = json.dumps(first, separators=(",", ":"))
    data_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="AUC|fingerprint"):
        read_augur(root)


def test_augur_artifact_path_does_not_construct_the_legacy_live_result(monkeypatch, augur_adata, science):
    def reject_live_result(*_args, **_kwargs):
        raise AssertionError("worker artifact encoding must not construct AugurResult")

    monkeypatch.setattr(AugurResult, "__init__", reject_live_result)
    tables, summary, metadata = _run(augur_adata, science, _runner=run_augur_artifact)

    assert list(tables) == list(AUGUR_VIEWS)
    assert summary["node_id"] == "OpenBioSingleCellAugur"
    assert metadata["artifact_type"] == AUGUR_ARTIFACT_TYPE


def test_augur_results_operation_reads_portable_artifact_without_a_live_result(
    tmp_path: Path, augur_adata, science
):
    live = _run(augur_adata, science)
    tables, summary, metadata = validate_augur_result(live)
    input_root = tmp_path / "augur-input"
    input_root.mkdir()
    write_augur(input_root, tables, summary, metadata)
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = augur_results(
        context,
        {
            "result": {
                "type": "artifact",
                "kind": AUGUR_ARTIFACT_TYPE,
                "codec": AUGUR_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"view": "priorities"},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["table", "summary", "code"]
    restored, _ = read_table(staging / records[0]["payload"])
    science.pd.testing.assert_frame_equal(restored, tables["priorities"])


def test_augur_plot_operation_renders_priorities_without_changing_the_artifact(
    tmp_path: Path, augur_adata, science
):
    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))
    input_root = tmp_path / "augur-input"
    input_root.mkdir()
    write_augur(input_root, tables, summary, metadata)
    before = {path.relative_to(input_root): path.read_bytes() for path in input_root.rglob("*") if path.is_file()}
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = augur_plot(
        context,
        {
            "result": {
                "type": "artifact",
                "kind": AUGUR_ARTIFACT_TYPE,
                "codec": AUGUR_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"view": "priorities", "top_n": 20, "max_plot_rows": 100_000, "max_image_pixels": 40_000_000},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, plot_metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert plot_metadata["source"]["operation"] == "augur_plot"
    assert before == {
        path.relative_to(input_root): path.read_bytes() for path in input_root.rglob("*") if path.is_file()
    }


def test_augur_plot_renders_cross_validation_distribution(augur_adata, science):
    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))

    plotted, report, code = augur_plot_owned(
        tables,
        summary,
        metadata,
        view="cross_validation",
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["view"] == "cross_validation"
    assert "plot_augur_result" in code


def test_augur_plot_renders_fold_aggregated_feature_importance(augur_adata, science):
    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))

    plotted, report, _code = augur_plot_owned(
        tables,
        summary,
        metadata,
        view="feature_importance",
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["view"] == "feature_importance"


def test_augur_plot_renders_roc_from_retained_fold_predictions(augur_adata, science):
    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))

    plotted, report, _code = augur_plot_owned(
        tables,
        summary,
        metadata,
        view="roc",
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["view"] == "roc"
    assert report.summary["key_results"]["source_rows"] == 24
    assert report.summary["key_results"]["roc_curve_auc"] == pytest.approx(
        {"population_A": 1.0, "population_B": 0.5}
    )


@pytest.mark.parametrize("view", ["priorities", "cross_validation", "feature_importance", "roc"])
def test_augur_plot_generated_code_has_runtime_parity(augur_adata, science, view):
    result = _run(augur_adata, science)
    tables, summary, metadata = validate_augur_result(result)
    plotted, report, code = augur_plot_owned(
        tables,
        summary,
        metadata,
        view=view,
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )
    assert "openbio_singlecell" not in code
    namespace = {}
    exec(code, namespace)

    reproduced_png, reproduced_summary = namespace["plot_augur_result"](result)

    assert reproduced_png == plotted.png
    assert reproduced_summary == report.summary


def test_augur_plot_generated_code_rejects_tampered_portable_evidence(augur_adata, science):
    result = _run(augur_adata, science)
    code = augur_core.augur_plot_code(
        view="priorities",
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )
    namespace = {}
    exec(code, namespace)
    object.__getattribute__(result, "_tables")["priorities"].loc[0, "mean_accuracy"] = 0.25

    with pytest.raises(ValueError, match="fingerprint"):
        namespace["plot_augur_result"](result)


def test_augur_plot_enforces_closed_view_payload_and_resource_guards(augur_adata, science):
    assert OpenBioSingleCellAugurPlot.prepare_worker_arguments(
        {
            "result": "ticket",
            "view": {"view": "roc", "top_n": 7},
            "max_plot_rows": 100,
            "max_image_pixels": 1_000_000,
        }
    ) == {
        "result": "ticket",
        "view": "roc",
        "top_n": 7,
        "max_plot_rows": 100,
        "max_image_pixels": 1_000_000,
    }
    with pytest.raises(ValueError, match="inactive"):
        OpenBioSingleCellAugurPlot.prepare_worker_arguments(
            {"result": "ticket", "view": {"view": "roc", "top_n": 7, "population": "extra"}}
        )

    tables, summary, metadata = validate_augur_result(_run(augur_adata, science))
    with pytest.raises(ValueError, match="max_plot_rows"):
        augur_plot_owned(
            tables,
            summary,
            metadata,
            view="roc",
            top_n=20,
            max_plot_rows=23,
            max_image_pixels=40_000_000,
        )
    with pytest.raises(MemoryError, match="max_image_pixels"):
        augur_plot_owned(
            tables,
            summary,
            metadata,
            view="priorities",
            top_n=20,
            max_plot_rows=100_000,
            max_image_pixels=1,
        )

    result = _run(augur_adata, science)
    for code, exception, message in (
        (
            augur_core.augur_plot_code(
                view="roc", top_n=20, max_plot_rows=23, max_image_pixels=40_000_000
            ),
            ValueError,
            "max_plot_rows",
        ),
        (
            augur_core.augur_plot_code(
                view="priorities", top_n=20, max_plot_rows=100_000, max_image_pixels=1
            ),
            MemoryError,
            "max_image_pixels",
        ),
    ):
        namespace = {}
        exec(code, namespace)
        with pytest.raises(exception, match=message):
            namespace["plot_augur_result"](result)


def test_augur_operation_publishes_the_portable_codec(
    tmp_path: Path, monkeypatch, augur_adata, science
):
    live = _run(augur_adata, science)
    tables, summary, metadata = validate_augur_result(live)
    monkeypatch.setattr(
        "openbio_singlecell.operations_population.run_augur_artifact",
        lambda _adata, **_parameters: (tables, summary, metadata),
    )
    input_root = tmp_path / "adata-input"
    input_root.mkdir()
    write_anndata(input_root, augur_adata)
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = augur_operation(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {
            "sample_key": "sample",
            "population_key": "population",
            "condition_key": "condition",
            "control": "control",
            "treatment": "treatment",
            "classifier": "random_forest_classifier",
            "source": {"source": "raw"},
            "annotation_status": "provisional",
            "technical_batch_key": "technical_batch",
            "n_subsamples": 2,
            "subsample_size": 3,
            "folds": 2,
            "n_threads": 1,
            "random_seed": 17,
            "max_result_rows": 10_000,
            "max_result_mib": 10.0,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["result", "summary", "code"]
    assert records[0]["codec"] == AUGUR_CODEC
    restored_tables, restored_summary, restored_metadata = read_augur(staging / records[0]["payload"])
    assert restored_summary == summary
    assert restored_metadata == metadata
    for view in AUGUR_VIEWS:
        science.pd.testing.assert_frame_equal(restored_tables[view], tables[view])


def _set_raw_value(adata, value):
    raw = adata.raw.to_adata()
    raw.X = raw.X.astype(float)
    raw.X[0, 0] = value
    adata.raw = raw


def _collapse_control_samples(adata):
    obs = adata.obs.copy()
    selected = obs["condition"] == "control"
    obs.loc[selected, "sample"] = "control_1"
    obs.loc[selected, "technical_batch"] = "batch_1"
    adata.obs = obs


def _make_partially_imbalanced_batches(adata):
    mapping = {
        "control_1": "batch_1",
        "control_2": "batch_1",
        "treatment_1": "batch_1",
        "treatment_2": "batch_2",
    }
    obs = adata.obs.copy()
    obs["technical_batch"] = [mapping.get(sample, batch) for sample, batch in zip(
        obs["sample"], obs["technical_batch"], strict=True
    )]
    adata.obs = obs


def _make_population_specific_perfect_confounding(adata):
    obs = adata.obs.copy()
    obs.loc[obs["sample"] == "control_2", "technical_batch"] = "batch_1"
    treatment_2_rows = list(obs.index[obs["sample"] == "treatment_2"])
    for position, cell in enumerate(treatment_2_rows):
        obs.loc[cell, "sample"] = f"treatment_2_{position % 2 + 1}"
    population_a_treatment_1 = (
        (obs["population"] == "population_A") & (obs["sample"] == "treatment_1")
    )
    obs.loc[population_a_treatment_1, "population"] = "population_B"
    adata.obs = obs


def test_augur_schemas_are_atomic_and_typed():
    schema = OpenBioSingleCellAugur.define_schema()
    assert [item.id for item in schema.inputs] == [
        "adata",
        "sample_key",
        "population_key",
        "condition_key",
        "control",
        "treatment",
        "classifier",
        "source",
        "annotation_status",
        "technical_batch_key",
        "n_subsamples",
        "subsample_size",
        "folds",
        "n_threads",
        "random_seed",
        "max_result_rows",
        "max_result_mib",
    ]
    inputs = {item.id: item for item in schema.inputs}
    assert inputs["classifier"].options == list(AUGUR_CLASSIFIERS)
    assert [option.key for option in inputs["source"].options] == ["X", "raw", "layer"]
    assert all(inputs[name].advanced for name in (
        "technical_batch_key",
        "n_subsamples",
        "subsample_size",
        "folds",
        "n_threads",
        "random_seed",
        "max_result_rows",
        "max_result_mib",
    ))
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("result", AugurResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]

    result_schema = OpenBioSingleCellAugurResults.define_schema()
    assert [item.id for item in result_schema.inputs] == ["result", "view"]
    assert result_schema.inputs[0].io_type == AUGUR_ARTIFACT_TYPE
    assert result_schema.inputs[1].options == list(AUGUR_VIEWS)
    assert [(item.display_name, item.io_type) for item in result_schema.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]

    plot_schema = OpenBioSingleCellAugurPlot.define_schema()
    assert [item.id for item in plot_schema.inputs] == [
        "result",
        "view",
        "max_plot_rows",
        "max_image_pixels",
    ]
    view = plot_schema.inputs[1]
    assert [(option.key, [item.id for item in option.inputs]) for option in view.options] == [
        ("priorities", ["top_n"]),
        ("cross_validation", ["top_n"]),
        ("feature_importance", ["top_n"]),
        ("roc", ["top_n"]),
    ]
    assert [(item.display_name, item.io_type) for item in plot_schema.outputs] == [
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_augur_exact_backend_calls_summary_and_input_immutability(augur_adata, science):
    before = augur_adata.copy()
    fake = _fake_pertpy(science)
    result = _run(augur_adata, science, pertpy_module=fake)
    tables, summary, metadata = validate_augur_result(result)

    assert isinstance(result, AugurResult)
    assert result.artifact_type == AUGUR_ARTIFACT_TYPE
    assert len(result.fingerprint) == 64
    assert set(tables) == set(AUGUR_VIEWS)
    assert tuple(tables["priorities"].columns) == PRIORITY_COLUMNS
    assert tuple(tables["cross_validation"].columns) == CROSS_VALIDATION_COLUMNS
    assert tuple(tables["feature_importance"].columns) == FEATURE_IMPORTANCE_COLUMNS
    assert tuple(tables["predictions"].columns) == PREDICTION_COLUMNS
    assert tables["priorities"]["population"].tolist() == ["population_A", "population_B"]
    assert tables["priorities"]["rank"].tolist() == [1, 2]
    assert tables["cross_validation"].shape == (8, 4)
    assert tables["feature_importance"].shape == (16, 5)
    assert tables["predictions"].shape == (24, 6)
    assert summary["status"] == "exploratory_cell_level_cross_validation"
    assert summary["software_versions"]["pertpy"] == "1.3.0"
    assert len(summary["references"]) >= 3
    assert summary["key_results"]["raw_snapshot_features"] == augur_adata.raw.n_vars
    assert summary["key_results"]["input_current_features"] == augur_adata.n_vars
    assert summary["key_results"]["count_source"] == {
        "source": "raw.X",
        "selection": "explicit",
        "integer_like": True,
    }
    assert summary["key_results"]["prediction_label_encoding"] == {
        "negative_label": 0,
        "negative_condition": "control",
        "positive_label": 1,
        "positive_condition": "treatment",
    }
    assert summary["key_results"]["sample_audit"]["minimum_two_samples_per_arm"] is True
    assert summary["key_results"]["technical_batch_audit"]["perfect_condition_confounding"] is False
    auc_distribution = summary["key_results"]["cross_validation_auc_distribution"]
    assert [record["population"] for record in auc_distribution] == ["population_A", "population_B"]
    assert [record["observations"] for record in auc_distribution] == [4, 4]
    assert auc_distribution[0]["mean"] == pytest.approx(1.0)
    assert summary["key_results"]["skipped_populations"] == [
        {
            "population": "population_weak",
            "reasons": ["control_below_subsample_size", "treatment_below_subsample_size"],
            "cells_control": 2,
            "cells_treatment": 2,
            "samples_control": 2,
            "samples_treatment": 2,
        }
    ]
    assert any("splits cells" in warning for warning in summary["warnings"])
    assert "artifact_fingerprint_sha256" in metadata
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    constructor, load, predict = fake.calls
    assert constructor == {
        "operation": "constructor",
        "estimator": "random_forest_classifier",
        "n_estimators": 100,
        "max_depth": None,
        "max_features": 2,
        "penalty": "l2",
        "random_state": 17,
    }
    assert load["label_col"] == "condition"
    assert load["cell_type_col"] == "population"
    assert load["condition_label"] == "control"
    assert load["treatment_label"] == "treatment"
    assert set(load["conditions"]) == {"control", "treatment"}
    assert "other_1" not in load["samples"]
    assert load["n_vars"] == augur_adata.raw.n_vars
    assert predict == {
        "operation": "predict",
        "n_subsamples": 2,
        "subsample_size": 3,
        "folds": 2,
        "min_cells": None,
        "feature_perc": 0.5,
        "var_quantile": 0.5,
        "span": 0.75,
        "filter_negative_residuals": False,
        "n_threads": 1,
        "augur_mode": "default",
        "select_variance_features": True,
        "key_added": "openbio_augur",
        "random_state": 17,
        "zero_division": 0,
    }

    science.np.testing.assert_array_equal(augur_adata.X, before.X)
    science.np.testing.assert_array_equal(augur_adata.raw.X, before.raw.X)
    science.pd.testing.assert_frame_equal(augur_adata.obs, before.obs)
    science.pd.testing.assert_frame_equal(augur_adata.var, before.var)
    assert augur_adata.uns == before.uns


@pytest.mark.parametrize("classifier", list(AUGUR_CLASSIFIERS))
def test_augur_retains_predictions_from_fitted_fold_estimators(augur_adata, science, classifier):
    result = _run(
        augur_adata,
        science,
        pertpy_module=_fake_pertpy(science, native_evidence=False),
        classifier=classifier,
    )

    predictions = result.table("predictions")
    assert predictions.shape == (24, len(PREDICTION_COLUMNS))
    assert predictions.groupby(["population", "subsample"], observed=True).size().tolist() == [6, 6, 6, 6]


def test_augur_records_the_backend_positive_condition(augur_adata, science):
    _replace_obs_column(
        augur_adata,
        "condition",
        ["z_control" if value == "control" else "a_treatment" if value == "treatment" else value for value in augur_adata.obs["condition"]],
    )

    result = _run(augur_adata, science, control="z_control", treatment="a_treatment")

    assert result.summary["key_results"]["prediction_label_encoding"] == {
        "negative_label": 0,
        "negative_condition": "a_treatment",
        "positive_label": 1,
        "positive_condition": "z_control",
    }


def test_generated_augur_code_has_runtime_parity(augur_adata, science):
    fake = _fake_pertpy(science)
    result = _run(augur_adata, science, pertpy_module=fake)
    parameters = result.summary["parameters"]
    code = augur_core.augur_code(**parameters)
    assert "from openbio_singlecell" not in code
    compile(code, "<augur-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced_tables, reproduced_summary = namespace["run_augur_prioritization"](
        augur_adata, pertpy_module=_fake_pertpy(science)
    )
    for view in AUGUR_VIEWS:
        science.pd.testing.assert_frame_equal(reproduced_tables[view], result.table(view))
    assert reproduced_summary == result.summary


def test_generated_augur_code_handles_absent_technical_batch(augur_adata, science):
    result = _run(augur_adata, science, technical_batch_key="")
    assert result.summary["parameters"]["technical_batch_key"] is None
    code = augur_core.augur_code(**result.summary["parameters"])
    namespace = {}
    exec(code, namespace)
    reproduced_tables, reproduced_summary = namespace["run_augur_prioritization"](
        augur_adata, pertpy_module=_fake_pertpy(science)
    )
    for view in AUGUR_VIEWS:
        science.pd.testing.assert_frame_equal(reproduced_tables[view], result.table(view))
    assert reproduced_summary == result.summary


def test_augur_nodes_return_reports_code_and_pure_views(augur_adata, science, monkeypatch):
    monkeypatch.setattr(augur_core, "_require_pertpy", lambda: _fake_pertpy(science))
    result, report, code = augur_legacy_owned(
        augur_adata,
        sample_key="sample",
        population_key="population",
        condition_key="condition",
        control="control",
        treatment="treatment",
        source=None,
        annotation_status="curated",
        technical_batch_key="technical_batch",
        n_subsamples=2,
        subsample_size=3,
        folds=2,
        random_seed=17,
        max_result_rows=10_000,
        max_result_mib=10.0,
    )
    assert isinstance(result, AugurResult)
    assert report.summary == result.summary
    assert report.source["operation"] == "augur"
    compile(code, "<augur-node-code>", "exec")

    for view in AUGUR_VIEWS:
        table_result, view_report, view_code = augur_results_legacy_owned(
            result, view=view
        )
        assert table_result.source["operation"] == "augur_results"
        science.pd.testing.assert_frame_equal(table_result.table, result.table(view))
        assert view_report.summary["selected_view"]["view"] == view
        assert view_report.summary["selected_view"]["operation"] == "pure_artifact_view"
        namespace = {}
        exec(view_code, namespace)
        reproduced_table, reproduced_summary = namespace["select_augur_results"](result)
        science.pd.testing.assert_frame_equal(reproduced_table, table_result.table)
        assert reproduced_summary == view_report.summary


def test_augur_artifact_is_defensive_and_tamper_evident(augur_adata, science):
    result = _run(augur_adata, science)
    with pytest.raises(AttributeError, match="immutable"):
        result.extra = "forbidden"

    returned = result.table("priorities")
    returned.loc[0, "mean_auc"] = 0.01
    assert result.table("priorities").loc[0, "mean_auc"] != 0.01
    returned_summary = result.summary
    returned_summary["status"] = "changed"
    assert result.summary["status"] == "exploratory_cell_level_cross_validation"
    returned_predictions = result.table("predictions")
    returned_predictions.loc[0, "prediction_score"] = 0.0
    assert result.table("predictions").loc[0, "prediction_score"] != 0.0

    private_tables = object.__getattribute__(result, "_tables")
    private_tables["priorities"].loc[0, "mean_auc"] = 0.01
    with pytest.raises(ValueError, match="fingerprint"):
        validate_augur_result(result)
    with pytest.raises(ValueError, match="fingerprint"):
        select_augur_view(result, view="priorities")


def _reuse_prediction_observation_across_populations(table):
    table.loc[12, "observation"] = table.loc[0, "observation"]


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda table: table.__setitem__("prediction_score", [float("nan"), *table["prediction_score"].iloc[1:]]), "non-finite"),
        (lambda table: table.__setitem__("prediction_score", [2.0, *table["prediction_score"].iloc[1:]]), r"\[0, 1\]"),
        (lambda table: table.__setitem__("true_label", [2, *table["true_label"].iloc[1:]]), "binary"),
        (lambda table: table.__setitem__("observation", [*table["observation"].iloc[:3], table["observation"].iloc[0], *table["observation"].iloc[4:]]), "exactly once"),
        (_reuse_prediction_observation_across_populations, "more than one population"),
        (lambda table: table.__setitem__("prediction_score", [1.0, *table["prediction_score"].iloc[1:]]), "recomputed AUC"),
    ],
)
def test_augur_rejects_malformed_prediction_evidence(augur_adata, science, mutate, message):
    result = _run(augur_adata, science)
    predictions = object.__getattribute__(result, "_tables")["predictions"]
    mutate(predictions)

    with pytest.raises((TypeError, ValueError), match=message):
        validate_augur_result(result)


@pytest.mark.parametrize(
    ("mutate", "source_kind", "message"),
    [
        (lambda adata: setattr(adata, "raw", None), "raw", "Raw snapshot"),
        (lambda adata: adata.raw.X.__setitem__((0, 0), -1), "raw", "negative"),
        (lambda adata: _set_raw_value(adata, float("nan")), "raw", "non-finite"),
    ],
)
def test_augur_rejects_invalid_raw_snapshot_and_expression_source(
    augur_adata, science, mutate, source_kind, message
):
    mutate(augur_adata)
    with pytest.raises((TypeError, ValueError), match=message):
        _run(augur_adata, science, source_kind=source_kind)


@pytest.mark.parametrize(("source_kind", "layer_name"), [("raw", None), ("X", None), ("layer", "counts")])
def test_augur_preserves_expert_selected_noninteger_sources(
    augur_adata, science, source_kind, layer_name
):
    if source_kind == "raw":
        _set_raw_value(augur_adata, 0.5)
    elif source_kind == "X":
        augur_adata.X = augur_adata.X.astype(float)
        augur_adata.X[0, 0] = 0.5
        augur_adata.raw = None
    else:
        augur_adata.layers["counts"] = augur_adata.layers["counts"].astype(float)
        augur_adata.layers["counts"][0, 0] = 0.5
        augur_adata.raw = None

    summary = _run(
        augur_adata,
        science,
        source_kind=source_kind,
        layer_name=layer_name,
    ).summary

    source = summary["key_results"]["count_source"]
    assert source["integer_like"] is False
    assert "state" not in source
    assert "full_gene_status" not in source
    assert any("not integer-like" in warning for warning in summary["warnings"])
    if source_kind != "raw":
        assert summary["key_results"]["raw_snapshot_features"] is None
        assert summary["key_results"]["raw_snapshot_fingerprint_sha256"] is None


@pytest.mark.parametrize(
    ("mutate", "audit_path", "expected", "warning_text"),
    [
        (
            lambda adata: _replace_obs_column(
                adata,
                "sample",
                ["control_1" if value == "treatment_1" else value for value in adata.obs["sample"]],
            ),
            ("sample_audit", "sample_condition_mapping_valid"),
            False,
            "more than one Condition",
        ),
        (
            _collapse_control_samples,
            ("sample_audit", "minimum_two_samples_per_arm"),
            False,
            "fewer than two biological Samples",
        ),
        (
            lambda adata: _replace_obs_column(
                adata,
                "technical_batch",
                ["control_batch" if value == "control" else "treatment_batch"
                 for value in adata.obs["condition"]],
            ),
            ("technical_batch_audit", "perfect_condition_confounding"),
            True,
            "perfectly confounded",
        ),
        (
            lambda adata: _replace_obs_column(
                adata,
                "technical_batch",
                ["changed_batch" if index == 0 else value
                 for index, value in enumerate(adata.obs["technical_batch"])],
            ),
            ("technical_batch_audit", "sample_batch_mapping_valid"),
            False,
            "spans multiple Technical batches",
        ),
    ],
)
def test_augur_discloses_advisory_sample_and_batch_design(
    augur_adata, science, mutate, audit_path, expected, warning_text
):
    mutate(augur_adata)
    summary = _run(augur_adata, science).summary
    assert summary["key_results"][audit_path[0]][audit_path[1]] is expected
    assert any(warning_text in warning for warning in summary["warnings"])


def test_augur_rejects_noncanonical_population_labels(augur_adata, science):
    _replace_obs_column(
        augur_adata,
        "population",
        [1 if index == 0 else value for index, value in enumerate(augur_adata.obs["population"])],
    )
    with pytest.raises((TypeError, ValueError), match="canonical nonempty string"):
        _run(augur_adata, science)


def test_augur_discloses_partial_technical_batch_imbalance(augur_adata, science):
    _make_partially_imbalanced_batches(augur_adata)
    summary = _run(augur_adata, science).summary
    audit = summary["key_results"]["technical_batch_audit"]
    assert audit["perfect_condition_confounding"] is False
    assert audit["condition_batch_imbalanced"] is True
    assert audit["condition_batch_sample_counts"] == {
        "control": {"batch_1": 2, "batch_2": 0},
        "treatment": {"batch_1": 1, "batch_2": 1},
    }
    assert any("partial batch imbalance" in warning for warning in summary["warnings"])


def test_augur_discloses_population_specific_perfect_batch_confounding(augur_adata, science):
    _make_population_specific_perfect_confounding(augur_adata)
    summary = _run(augur_adata, science).summary
    records = summary["key_results"]["technical_batch_audit"]["population_batch_records"]
    record = next(item for item in records if item["population"] == "population_A")
    assert record["perfect_condition_confounding"] is True
    assert any("eligible population has perfect" in warning for warning in summary["warnings"])


@pytest.mark.parametrize(
    ("overrides", "exception", "message"),
    [
        ({"classifier": "random_forest_regressor"}, ValueError, "classifier-only"),
        ({"random_seed": 0}, ValueError, "random_seed"),
        ({"folds": 4, "subsample_size": 3}, ValueError, "cannot exceed"),
        ({"max_result_rows": 5}, ValueError, "estimated result size"),
        ({"max_result_mib": 0.0001}, MemoryError, "estimated canonical result memory"),
        ({"n_subsamples": 50}, ValueError, "no eligible population"),
    ],
)
def test_augur_rejects_invalid_parameters_and_guards(
    augur_adata, science, overrides, exception, message
):
    with pytest.raises(exception, match=message):
        _run(augur_adata, science, **overrides)


def test_augur_allows_advisory_sample_role_to_reuse_condition_column(augur_adata, science):
    summary = _run(augur_adata, science, sample_key="condition").summary
    audit = summary["key_results"]["sample_audit"]
    assert audit["minimum_two_samples_per_arm"] is False
    assert audit["identity_status"] == "role_reused_unverified"
    assert audit["reused_roles"] == ["Condition"]
    assert any("fewer than two biological Samples" in warning for warning in summary["warnings"])
    assert any("Sample identity is unverified" in warning for warning in summary["warnings"])


def _drop_backend_fold(results):
    results["full_results"] = results["full_results"].iloc[1:].reset_index(drop=True)


def _bad_backend_auc(results):
    results["full_results"].loc[0, "augur_score"] = 1.5


def _bad_backend_fractional_fold(results):
    frame = results["full_results"].copy()
    frame["folds"] = frame["folds"].astype(float)
    frame.at[0, "folds"] = 0.5
    results["full_results"] = frame


def _bad_backend_summary_mean(results):
    results["summary_metrics"].loc["mean_auc", "population_A"] = 0.1


def _bad_backend_identifier(results):
    frame = results["full_results"].copy()
    frame["cell_type"] = frame["cell_type"].astype(object)
    frame.at[0, "cell_type"] = 7
    results["full_results"] = frame


def _bad_backend_feature_schema(results):
    results["feature_importances"] = results["feature_importances"].drop(columns="genes")


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (_drop_backend_fold, "complete population/subsample/fold family"),
        (_bad_backend_auc, "AUC must lie"),
        (_bad_backend_fractional_fold, "exact integers"),
        (_bad_backend_summary_mean, "does not equal"),
        (_bad_backend_identifier, "canonical nonempty string"),
        (_bad_backend_feature_schema, "feature_importances schema changed"),
    ],
)
def test_augur_rejects_malformed_pertpy_results(augur_adata, science, mutator, message):
    fake = _fake_pertpy(science, backend_mutator=mutator)
    with pytest.raises((TypeError, ValueError, RuntimeError), match=message):
        _run(augur_adata, science, pertpy_module=fake)


def test_augur_accepts_capability_checked_pertpy_1_3_patch_family(augur_adata, science):
    compatible_patch = _fake_pertpy(science)
    compatible_patch.__version__ = "1.3.1"
    result = _run(augur_adata, science, pertpy_module=compatible_patch)
    assert result.summary["software_versions"]["pertpy"] == "1.3.1"
    validate_augur_result(result)

    wrong_version = _fake_pertpy(science)
    wrong_version.__version__ = "1.2.0"
    with pytest.raises(RuntimeError, match=r"compatible Pertpy 1\.3\.x"):
        _run(augur_adata, science, pertpy_module=wrong_version)

    class BadAugur:
        def __init__(self, estimator):
            self.estimator = estimator

        def load(self, input):
            return input

        def predict(self, adata):
            return adata, {}

    wrong_interface = SimpleNamespace(__version__="1.3.0", tl=SimpleNamespace(Augur=BadAugur))
    with pytest.raises(RuntimeError, match="interface is incompatible"):
        _run(augur_adata, science, pertpy_module=wrong_interface)


@pytest.mark.parametrize("classifier", list(AUGUR_CLASSIFIERS))
def test_augur_is_deterministic_for_each_classifier(augur_adata, science, classifier):
    first = _run(augur_adata, science, classifier=classifier)
    second = _run(augur_adata, science, classifier=classifier)
    assert first.fingerprint == second.fingerprint
    assert first.summary == second.summary
    for view in AUGUR_VIEWS:
        science.pd.testing.assert_frame_equal(first.table(view), second.table(view))


def test_augur_rejects_artifact_numeric_strings(augur_adata, science):
    result = _run(augur_adata, science)
    private_tables = object.__getattribute__(result, "_tables")
    private_tables["cross_validation"]["fold"] = private_tables["cross_validation"]["fold"].astype(str)
    with pytest.raises(TypeError, match="integer dtype"):
        validate_augur_result(result)


def test_real_pertpy_1_3_cpu_smoke(science):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            import pertpy
        except (ImportError, OSError) as exc:
            pytest.skip(f"Pertpy is unavailable: {exc}")
    if augur_core._pertpy_version(pertpy) != "1.3.0":
        pytest.skip("Real smoke requires exact Pertpy 1.3.0.")
    adata = _make_augur_adata(
        science,
        n_genes=60,
        cells_per_sample_population=6,
        include_weak_population=False,
        include_other_condition=False,
    )
    result = _run(
        adata,
        science,
        pertpy_module=pertpy,
        subsample_size=4,
        n_threads=2,
        max_result_rows=100_000,
        max_result_mib=100.0,
    )
    tables, summary, _ = validate_augur_result(result)
    assert tables["priorities"].shape == (2, len(PRIORITY_COLUMNS))
    assert tables["cross_validation"].shape == (8, len(CROSS_VALIDATION_COLUMNS))
    assert not tables["feature_importance"].empty
    assert tables["predictions"].shape == (32, len(PREDICTION_COLUMNS))
    assert summary["software_versions"]["pertpy"] == "1.3.0"
    plotted, report, _code = augur_plot_owned(
        tables,
        summary,
        result.metadata,
        view="roc",
        top_n=20,
        max_plot_rows=100_000,
        max_image_pixels=40_000_000,
    )
    assert plotted.png.startswith(b"\x89PNG")
    assert report.summary["key_results"]["roc_curve_auc"] == pytest.approx(
        dict(zip(tables["priorities"]["population"], tables["priorities"]["mean_auc"], strict=True))
    )
