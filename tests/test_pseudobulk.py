from __future__ import annotations

import json
import sys
import types
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_table, write_anndata
from openbio_singlecell.node_types import PseudobulkType, SummaryResultType
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellPseudobulk as _PseudobulkSchema,
)
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellPseudobulkDESeq2 as _PseudobulkDESeq2Schema,
)
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellPseudobulkEdgeR as _PseudobulkEdgeRSchema,
)
from openbio_singlecell.operations_differential import (
    pseudobulk as pseudobulk_operation,
)
from openbio_singlecell.operations_differential import (
    pseudobulk_deseq2 as pseudobulk_deseq2_operation,
)
from openbio_singlecell.operations_differential import (
    pseudobulk_edger as pseudobulk_edger_operation,
)
from openbio_singlecell.operations_differential import (
    run_pseudobulk_deseq2_owned,
    run_pseudobulk_edger_owned,
    run_pseudobulk_owned,
)
from openbio_singlecell.pseudobulk import PseudobulkArtifact, validate_pseudobulk_artifact
from openbio_singlecell.pseudobulk_artifact_codec import (
    PSEUDOBULK_METADATA,
    read_pseudobulk,
    write_pseudobulk,
)
from openbio_singlecell.sample_design import EDGER_COLUMNS, PYDESEQ2_COLUMNS, prepare_pseudobulk_design
from openbio_singlecell.worker_protocol import OperationContext


def _matrix_values(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _fixture_adata(science, *, sparse_source=False, annotation_status="curated"):
    rows = []
    counts = []
    sample_design = {
        "s1": ("control", "b1", 20.0),
        "s2": ("control", "b2", 34.0),
        "s3": ("control", "b1", 41.0),
        "s4": ("treated", "b2", 27.0),
        "s5": ("treated", "b1", 53.0),
        "s6": ("treated", "b2", 62.0),
    }
    index = 0
    for sample, (condition, batch, age) in sample_design.items():
        for replicate in range(2):
            rows.append(
                {
                    "sample": sample,
                    "cell_type": "T",
                    "condition": condition,
                    "batch": batch,
                    "age": age,
                }
            )
            base = 2 if condition == "control" else 5
            counts.append([base + replicate, 3 + replicate, 1, 2])
            index += 1
    for sample, repeats in (("s1", 1), ("s2", 2), ("s4", 1), ("s5", 2)):
        condition, batch, age = sample_design[sample]
        for replicate in range(repeats):
            rows.append(
                {
                    "sample": sample,
                    "cell_type": "B",
                    "condition": condition,
                    "batch": batch,
                    "age": age,
                }
            )
            counts.append([1 + replicate, 2, 2, 1])
            index += 1
    obs = science.pd.DataFrame(rows, index=[f"cell_{idx}" for idx in range(index)])
    var = science.pd.DataFrame(index=["G1", "G2", "G3", "G4"])
    matrix = science.np.asarray(counts, dtype=science.np.int64)
    if sparse_source == "csc":
        stored = science.sparse.csc_matrix(matrix)
    elif sparse_source:
        stored = science.sparse.csr_matrix(matrix)
    else:
        stored = matrix.copy()
    adata = science.ad.AnnData(stored.copy(), obs=obs, var=var)
    adata.layers["counts"] = stored.copy()
    raw_source = adata.copy()
    raw_source.X = stored.copy()
    adata.raw = raw_source
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "version": "0.2.0",
        "display_name": "fixture",
        "source": {},
        "random_seed": 0,
        "warnings": [],
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {
                    "destination": "layer_and_raw",
                    "layer_name": "counts",
                    "source": "X",
                },
            }
        },
        "annotations": {
            "cell_type": {
                "schema_version": 1,
                "operation": "map_cluster_annotations",
                "annotation_status": annotation_status,
                "curation_assertion": "caller_declared" if annotation_status == "curated" else None,
                "output_column": "cell_type",
            }
        },
    }
    return adata


def _fake_decoupler(science, *, qc_name="psbulk_cells", malformed=None):
    module = types.ModuleType("decoupler")
    module.__version__ = "2.2.0"
    module.calls = []

    def pseudobulk(
        adata,
        sample_col,
        groups_col,
        layer=None,
        raw=False,
        empty=False,
        mode="sum",
        skip_checks=False,
        bsize=250_000,
        verbose=False,
    ):
        module.calls.append(
            {
                "sample_col": sample_col,
                "groups_col": groups_col,
                "layer": layer,
                "raw": raw,
                "empty": empty,
                "mode": mode,
                "skip_checks": skip_checks,
                "bsize": bsize,
                "verbose": verbose,
                "input": adata.copy(),
            }
        )
        matrix = _matrix_values(adata.X, science)
        samples = sorted(set(adata.obs[sample_col].tolist()))
        populations = sorted(set(adata.obs[groups_col].tolist()))
        aggregate_rows = []
        aggregate_obs = []
        aggregate_names = []
        for population in populations:
            for sample in samples:
                mask = (adata.obs[sample_col] == sample).to_numpy() & (
                    adata.obs[groups_col] == population
                ).to_numpy()
                selected = matrix[mask]
                summed = selected.sum(axis=0) if selected.size else science.np.zeros(adata.n_vars, dtype=int)
                aggregate_rows.append(summed)
                aggregate_obs.append(
                    {
                        sample_col: sample,
                        groups_col: population,
                        qc_name: int(mask.sum()),
                        "psbulk_counts": float(summed.sum()),
                    }
                )
                aggregate_names.append(f"{sample}_{population}")
        result = science.ad.AnnData(
            science.np.asarray(aggregate_rows, dtype=float),
            obs=science.pd.DataFrame(aggregate_obs, index=aggregate_names),
            var=adata.var.copy(),
        )
        if malformed == "counts":
            result.X[0, 0] += 1
        elif malformed == "identity":
            result.obs.iloc[0, result.obs.columns.get_loc(sample_col)] = "wrong"
        elif malformed == "qc":
            result.obs.iloc[0, result.obs.columns.get_loc(qc_name)] += 1
        elif malformed == "genes":
            result.var_names = ["WRONG", *result.var_names[1:]]
        return result

    def filter_samples(adata, min_cells=10, min_counts=1000, inplace=True):
        module.calls.append(
            {
                "filter_samples": True,
                "min_cells": min_cells,
                "min_counts": min_counts,
                "inplace": inplace,
            }
        )
        keep = (adata.obs["psbulk_cells"] >= min_cells) & (adata.obs["psbulk_counts"] >= min_counts)
        names = adata.obs_names[keep].to_numpy(dtype=str)
        if malformed == "filter" and names.size:
            names = names[:-1]
        if inplace:
            adata._inplace_subset_obs(names)
            return None
        return names

    module.pp = types.SimpleNamespace(pseudobulk=pseudobulk, filter_samples=filter_samples)
    return module


def _run_aggregation(
    science,
    monkeypatch,
    *,
    adata=None,
    sparse_source=False,
    qc_name="psbulk_cells",
    source=None,
    min_cells=2,
    min_counts=1,
    inference_mode="formal",
    technical_batch_key="batch",
    categorical_covariate_keys="",
    continuous_covariate_keys="age",
):
    adata = adata if adata is not None else _fixture_adata(science, sparse_source=sparse_source)
    backend = _fake_decoupler(science, qc_name=qc_name)
    monkeypatch.setitem(sys.modules, "decoupler", backend)
    artifact, report, code = run_pseudobulk_owned(
        adata,
        sample_key="sample",
        population_key="cell_type",
        condition_key="condition",
        technical_batch_key=technical_batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        source=source or {"source": "raw"},
        inference_mode=inference_mode,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    return adata, backend, artifact, report, code


def test_pseudobulk_worker_operation_publishes_portable_state_without_mutating_input(
    science, monkeypatch, tmp_path
):
    adata = _fixture_adata(science)
    monkeypatch.setitem(sys.modules, "decoupler", _fake_decoupler(science))
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    input_bytes = input_path.read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = pseudobulk_operation(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            }
        },
        {
            "sample_key": "sample",
            "population_key": "cell_type",
            "condition_key": "condition",
            "technical_batch_key": "batch",
            "categorical_covariate_keys": "",
            "continuous_covariate_keys": "age",
            "source": {"source": "raw"},
            "inference_mode": "formal",
            "min_cells": 2,
            "min_counts": 1,
        },
    )

    assert input_path.read_bytes() == input_bytes
    assert [record["name"] for record in records] == ["pseudobulk", "summary", "code"]
    assert records[0] == {
        "type": "artifact",
        "name": "pseudobulk",
        "kind": "OPENBIO_SINGLE_CELL_PSEUDOBULK",
        "codec": "pseudobulk-h5ad-v1",
        "payload": "outputs/pseudobulk",
    }
    output_root = staging / records[0]["payload"]
    assert output_root != input_root
    assert (output_root / ANNDATA_PAYLOAD).is_file()
    assert json.loads((output_root / PSEUDOBULK_METADATA).read_text(encoding="utf-8"))
    restored = read_pseudobulk(output_root)
    restored_adata, restored_metadata = validate_pseudobulk_artifact(restored)
    assert restored_adata.shape == (8, 4)
    assert restored_metadata["role_keys"] == {
        "sample": "sample",
        "population": "cell_type",
        "condition": "condition",
        "technical_batch": "batch",
    }
    assert records[1]["type"] == "summary"
    assert isinstance(records[1]["value"], dict)
    assert records[2]["type"] == "string"
    assert "def run_pseudobulk" in records[2]["value"]


def _fake_pertpy(science, *, engine, malformed=None):
    module = types.ModuleType("pertpy")
    module.__version__ = "1.3.0"
    module.calls = []
    if engine == "edger":
        module._openbio_backend_versions = {
            "pertpy": "1.3.0",
            "rpy2": "3.6.3",
            "R": "R version 4.6.0",
            "edgeR": "4.10.3",
            "BiocParallel": "1.42.1",
            "RhpcBLASctl": "0.23-42",
        }

        class EdgeR:
            def __init__(self, adata, design, *, mask=None, layer=None, **kwargs):
                self.adata = adata.copy()
                self.design = design.copy()
                self.normalization_factors = science.np.linspace(0.8, 1.2, adata.n_obs)
                module.calls.append(("init", self.adata.copy(), self.design.copy(), mask, layer, kwargs))

            def fit(self, **kwargs):
                module.calls.append(("fit", dict(kwargs)))
                if malformed == "backend_gene_mutation":
                    self.adata.var_names = [f"MUT{index}" for index in range(self.adata.n_vars)]
                elif malformed == "backend_count_mutation":
                    self.adata.X[0, 0] += 1
                elif malformed == "backend_fractional_count_mutation":
                    values = science.np.asarray(self.adata.X, dtype=float).copy()
                    values[0, 0] += 0.5
                    self.adata.X = values
                elif malformed == "backend_design_mutation":
                    self.design.iloc[0, 0] += 1

            def test_contrasts(self, contrasts, **kwargs):
                module.calls.append(("test", science.np.asarray(contrasts).copy(), dict(kwargs)))
                if malformed == "backend_contrast_mutation":
                    contrasts[0] += 1
                genes = self.adata.var_names.tolist()
                frame = science.pd.DataFrame(
                    {
                        "variable": genes,
                        "log_fc": [1.5, -1.2, 0.3, -0.1][: len(genes)],
                        "logCPM": [5.0, 4.0, 3.0, 2.0][: len(genes)],
                        "F": [12.0, 9.0, 1.0, 0.2][: len(genes)],
                        "p_value": [0.001, 0.01, 0.3, 0.8][: len(genes)],
                        "adj_p_value": [0.004, 0.02, 0.4, 0.8][: len(genes)],
                        "contrast": None,
                    }
                )
                if malformed == "missing_gene":
                    frame = frame.iloc[:-1].copy()
                elif malformed == "nonfinite":
                    frame.loc[0, "F"] = science.np.nan
                elif malformed == "p_range":
                    frame.loc[0, "adj_p_value"] = 1.2
                elif malformed == "wrong_bh":
                    frame.loc[0, "adj_p_value"] = 0.9
                elif malformed == "schema":
                    frame["unexpected"] = 1
                return frame

        module.tl = types.SimpleNamespace(EdgeR=EdgeR)
    else:
        module._openbio_backend_versions = {
            "pertpy": "1.3.0",
            "pydeseq2": "0.5.4",
        }

        class PyDESeq2:
            def __init__(self, adata, design, *, mask=None, layer=None, **kwargs):
                self.adata = adata.copy()
                self.design = design.copy()
                dds_obs = self.adata.obs.copy(deep=True)
                dds_var = self.adata.var.copy(deep=True)
                dds_var["dispersions"] = science.np.linspace(0.1, 0.4, adata.n_vars)
                self.dds = types.SimpleNamespace(
                    X=self.adata.X.copy(),
                    obs_names=dds_obs.index.copy(),
                    var_names=dds_var.index.copy(),
                    obsm={
                        "design_matrix": self.design.copy(),
                        "size_factors": science.np.linspace(0.9, 1.1, adata.n_obs),
                    },
                    obs=dds_obs,
                    var=dds_var,
                    uns={"disp_function_type": "parametric"},
                    fit_type="parametric",
                    size_factors_fit_type="ratio",
                    control_genes=None,
                    min_mu=0.5,
                    min_disp=1e-8,
                    max_disp=max(10.0, float(adata.n_obs)),
                    refit_cooks=True,
                    min_replicates=7,
                    beta_tol=1e-8,
                    quiet=False,
                    low_memory=False,
                )
                module.calls.append(("init", self.adata.copy(), self.design.copy(), mask, layer, kwargs))

            def fit(self, **kwargs):
                module.calls.append(("fit", dict(kwargs)))
                dds_values = _matrix_values(self.dds.X, science)
                if bool((dds_values == 0).any(axis=0).all()):
                    import warnings

                    warnings.warn(
                        "Every gene contains at least one zero, cannot compute log geometric means. "
                        "Switching to iterative mode.",
                        UserWarning,
                        stacklevel=2,
                    )
                if malformed == "backend_gene_mutation":
                    self.adata.var_names = [f"MUT{index}" for index in range(self.adata.n_vars)]
                elif malformed == "backend_count_mutation":
                    self.adata.X[0, 0] += 1
                elif malformed == "backend_fractional_count_mutation":
                    values = science.np.asarray(self.adata.X, dtype=float).copy()
                    values[0, 0] += 0.5
                    self.adata.X = values
                elif malformed == "backend_design_mutation":
                    self.design.iloc[0, 0] += 1
                elif malformed == "dds_count_mutation":
                    self.dds.X[0, 0] += 1
                elif malformed == "dds_large_integer_count_mutation":
                    self.dds.X[0, 0] += 1
                elif malformed == "dds_sample_axis_mutation":
                    self.dds.obs_names = science.pd.Index(
                        ["DDS_MUTATED", *self.dds.obs_names[1:]],
                    )
                elif malformed == "dds_gene_axis_mutation":
                    self.dds.var_names = science.pd.Index(
                        ["DDS_MUTATED", *self.dds.var_names[1:]],
                    )
                elif malformed == "dds_obs_metadata_mutation":
                    self.dds.obs.iloc[0, 0] = "DDS_MUTATED"
                elif malformed == "dds_var_metadata_mutation":
                    self.dds.var.loc[self.dds.var_names[0], "__openbio_feature_identity__"] = "DDS_MUTATED"
                elif malformed == "dds_design_mutation":
                    self.dds.obsm["design_matrix"].iloc[0, 0] += 1
                elif malformed == "dds_fixed_control_mutation":
                    self.dds.min_replicates = 8
                elif malformed == "dispersion_fallback":
                    import warnings

                    self.dds.uns["disp_function_type"] = "mean"
                    warnings.warn(
                        "The dispersion trend curve fitting did not converge. Switching to a mean-based dispersion trend.",
                        UserWarning,
                        stacklevel=2,
                    )

            def test_contrasts(self, contrasts, **kwargs):
                module.calls.append(("test", science.np.asarray(contrasts).copy(), dict(kwargs)))
                if malformed == "backend_contrast_mutation":
                    contrasts[0] += 1
                genes = self.adata.var_names.tolist()
                frame = science.pd.DataFrame(
                    {
                        "variable": genes,
                        "baseMean": [100.0, 80.0, 20.0, 5.0][: len(genes)],
                        "log_fc": [1.4, -1.1, 0.2, science.np.nan][: len(genes)],
                        "lfcSE": [0.2, 0.3, 0.4, science.np.nan][: len(genes)],
                        "stat": [7.0, -3.7, 0.5, science.np.nan][: len(genes)],
                        "p_value": [0.001, 0.02, 0.5, science.np.nan][: len(genes)],
                        "adj_p_value": [0.002, 0.02, science.np.nan, science.np.nan][: len(genes)],
                        "contrast": None,
                    }
                )
                if malformed == "missing_gene":
                    frame = frame.iloc[:-1].copy()
                elif malformed == "negative_base":
                    frame.loc[0, "baseMean"] = -1
                elif malformed == "infinite":
                    frame.loc[0, "stat"] = science.np.inf
                elif malformed == "wrong_bh":
                    frame.loc[0, "adj_p_value"] = 0.9
                elif malformed == "schema":
                    frame["unexpected"] = 1
                return frame

        module.tl = types.SimpleNamespace(PyDESeq2=PyDESeq2)
    return module


@pytest.mark.parametrize(
    ("operation", "engine", "expected_columns"),
    [
        (pseudobulk_edger_operation, "edger", EDGER_COLUMNS),
        (pseudobulk_deseq2_operation, "pydeseq2", PYDESEQ2_COLUMNS),
    ],
)
def test_pseudobulk_engine_worker_operations_publish_jsonl_tables(
    science,
    monkeypatch,
    tmp_path,
    operation,
    engine,
    expected_columns,
):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    input_root = tmp_path / "pseudobulk"
    input_root.mkdir()
    write_pseudobulk(input_root, artifact)
    monkeypatch.setitem(sys.modules, "pertpy", _fake_pertpy(science, engine=engine))
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))
    parameters = {
        "population": "T",
        "reference_condition": "control",
        "comparison_condition": "treated",
        "categorical_covariate_keys": "",
        "continuous_covariate_keys": "",
        "fdr_threshold": 0.05,
        "min_abs_log2_fold_change": 0.0,
        "min_count": 0,
        "min_total_count": 1,
        "large_n": 10,
        "min_prop": 0.7,
    }
    if engine == "pydeseq2":
        parameters["n_cpus"] = 1

    records = operation(
        context,
        {
            "pseudobulk": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_SINGLE_CELL_PSEUDOBULK",
                "codec": "pseudobulk-h5ad-v1",
            }
        },
        parameters,
    )

    assert [record["name"] for record in records] == ["table", "summary", "code"]
    assert records[0]["kind"] == "OPENBIO_SINGLE_CELL_TABLE"
    assert records[0]["codec"] == "table-jsonl-v1"
    table, metadata = read_table(staging / records[0]["payload"])
    assert table.columns.tolist() == expected_columns
    assert metadata["source"]["operation"] in {"pseudobulk_edger", "pseudobulk_deseq2"}


def test_pseudobulk_schemas_are_typed():
    aggregation = _PseudobulkSchema.define_schema()
    edger = _PseudobulkEdgeRSchema.define_schema()
    pydeseq2 = _PseudobulkDESeq2Schema.define_schema()

    assert [item.id for item in aggregation.inputs] == [
        "adata",
        "sample_key",
        "population_key",
        "condition_key",
        "technical_batch_key",
        "categorical_covariate_keys",
        "continuous_covariate_keys",
        "source",
        "inference_mode",
        "min_cells",
        "min_counts",
    ]
    assert [(item.display_name, item.io_type) for item in aggregation.outputs] == [
        ("pseudobulk", PseudobulkType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert edger.inputs[0].io_type == PseudobulkType.io_type
    assert pydeseq2.inputs[0].io_type == PseudobulkType.io_type
    assert [item.display_name for item in edger.outputs] == ["table", "summary", "code"]
    assert [item.display_name for item in pydeseq2.outputs] == ["table", "summary", "code"]
    assert pydeseq2.display_name == "Pseudobulk PyDESeq2 Contrast"


@pytest.mark.parametrize(
    ("sparse_source", "qc_name"),
    [(False, "psbulk_cells"), (True, "psbulk_n_cells"), ("csc", "psbulk_cells")],
)
def test_pseudobulk_exact_sums_missing_profiles_filter_summary_code_and_immutability(
    science, monkeypatch, sparse_source, qc_name
):
    adata = _fixture_adata(science, sparse_source=sparse_source)
    original = adata.copy()
    adata, backend, artifact, report, code = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        qc_name=qc_name,
    )

    output, metadata = validate_pseudobulk_artifact(artifact)
    assert output.obs[["sample", "cell_type"]].to_dict("records") == [
        {"sample": "s1", "cell_type": "T"},
        {"sample": "s2", "cell_type": "T"},
        {"sample": "s3", "cell_type": "T"},
        {"sample": "s4", "cell_type": "T"},
        {"sample": "s5", "cell_type": "T"},
        {"sample": "s6", "cell_type": "T"},
        {"sample": "s2", "cell_type": "B"},
        {"sample": "s5", "cell_type": "B"},
    ]
    expected_t_s1 = _matrix_values(adata.raw.X, science)[:2].sum(axis=0)
    science.np.testing.assert_array_equal(_matrix_values(output.X, science)[0], expected_t_s1)
    assert output.obs["openbio_n_cells"].tolist() == [2, 2, 2, 2, 2, 2, 2, 2]
    assert metadata["count_state"] == "raw_counts"
    assert metadata["aggregation_mode"] == "sum"
    assert metadata["backend"]["cell_qc_column"] == qc_name
    assert report.summary["key_results"]["profiles_before_qc"] == 10
    assert report.summary["key_results"]["profiles_retained"] == 8
    assert report.summary["key_results"]["profiles_removed"] == 2
    assert report.summary["key_results"]["missing_combination_count"] == 2
    assert report.summary["key_results"]["backend_synthetic_empty_profiles_removed"] == 2
    assert len(report.summary["key_results"]["retained_profile_records_preview"]) == 8
    assert len(report.summary["key_results"]["removed_profile_records_preview"]) == 2
    assert report.summary["key_results"]["retained_profile_records_preview_truncated"] is False
    assert report.summary["key_results"]["removed_profile_records_preview_truncated"] is False
    assert report.summary["key_results"]["strata_preview_truncated"] is False
    assert report.summary["parameters"]["mode"] == "sum"
    assert report.summary["parameters"]["synthetic_empty_profiles"] is False
    assert "Sample, not cell" in report.summary["methods"]
    assert "did not test a Condition effect" in report.summary["results"]
    assert set(report.summary["software_versions"]) == {
        "python",
        "openbio-singlecell",
        "decoupler",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
    }
    json.dumps(report.summary, allow_nan=False)
    assert backend.calls[0]["raw"] is False
    assert backend.calls[0]["empty"] is False
    assert backend.calls[0]["mode"] == "sum"
    assert backend.calls[0]["skip_checks"] is False
    assert backend.calls[1] == {
        "filter_samples": True,
        "min_cells": 2,
        "min_counts": 1,
        "inplace": False,
    }
    compile(code, "<pseudobulk-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk"](adata)
    reproduced_adata, reproduced_metadata = validate_pseudobulk_artifact(reproduced)
    science.np.testing.assert_array_equal(_matrix_values(reproduced_adata.X, science), _matrix_values(output.X, science))
    science.pd.testing.assert_frame_equal(reproduced_adata.obs, output.obs)
    assert reproduced_metadata == metadata
    assert reproduced_summary == report.summary
    science.np.testing.assert_array_equal(_matrix_values(adata.X, science), _matrix_values(original.X, science))
    science.np.testing.assert_array_equal(
        _matrix_values(adata.layers["counts"], science),
        _matrix_values(original.layers["counts"], science),
    )
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert adata.uns == original.uns


def test_pseudobulk_explicit_x_source_is_independent_of_raw_and_formal_status_is_disclosed(
    science,
    monkeypatch,
):
    adata = _fixture_adata(science)
    selected_x = _matrix_values(adata.X, science).astype(science.np.int64) + 7
    adata.X = selected_x.copy()
    _, _, artifact, x_report, x_code = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        source={"source": "X"},
        min_cells=1,
    )
    output, metadata = validate_pseudobulk_artifact(artifact)
    science.np.testing.assert_array_equal(_matrix_values(output.X, science)[0], selected_x[:2].sum(axis=0))
    assert metadata["count_source"]["raw_equivalence_required"] is False
    assert metadata["count_source"]["history_provenance_required"] is False
    assert metadata["count_source"]["feature_axis"] == "adata.var_names"
    assert x_report.summary["key_results"]["declarations"]["count_source"] == (
        "caller_selected_and_declared_raw_counts"
    )
    namespace = {}
    exec(x_code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk"](adata)
    reproduced_adata, reproduced_metadata = validate_pseudobulk_artifact(reproduced)
    science.pd.testing.assert_frame_equal(reproduced_adata.obs, output.obs)
    science.np.testing.assert_array_equal(_matrix_values(reproduced_adata.X, science), _matrix_values(output.X, science))
    assert reproduced_metadata == metadata
    assert reproduced_summary == x_report.summary

    provisional = _fixture_adata(science, annotation_status="provisional")
    _, _, _, report, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=provisional,
        min_cells=1,
        inference_mode="formal",
    )
    assert report.summary["key_results"]["annotation_status"] == "provisional"
    assert report.summary["key_results"]["formal_interpretation_invalid"] is True
    assert "population_annotation_not_curated" in report.summary["key_results"][
        "formal_interpretation_invalid_reasons"
    ]
    assert any("formal interpretation" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
def test_pseudobulk_exact_profile_total_avoids_false_int64_overflow_and_rejects_true_overflow(
    science,
    monkeypatch,
    storage,
):
    adata = _fixture_adata(science)
    selected = _matrix_values(adata.X, science).astype(science.np.int64, copy=True)
    selected[0] = [2**62, 0, 0, 0]
    selected[1] = [0, 1, 0, 0]
    single_cell_profile = science.np.flatnonzero(
        (adata.obs["sample"].to_numpy() == "s1") & (adata.obs["cell_type"].to_numpy() == "B")
    )[0]
    selected[single_cell_profile] = [science.np.iinfo(science.np.int64).max, 0, 0, 0]
    adata.X = (
        science.sparse.csr_matrix(selected)
        if storage == "csr"
        else science.sparse.csc_matrix(selected)
        if storage == "csc"
        else selected
    )
    _, backend, artifact, report, code = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        source={"source": "X"},
        min_cells=1,
    )
    output, _ = validate_pseudobulk_artifact(artifact)
    row_index = science.np.flatnonzero(
        (output.obs["sample"].to_numpy() == "s1") & (output.obs["cell_type"].to_numpy() == "T")
    )[0]
    exact_total = sum(int(value) for value in _matrix_values(output.X, science)[row_index])
    assert exact_total == 2**62 + 1
    max_profile_index = science.np.flatnonzero(
        (output.obs["sample"].to_numpy() == "s1") & (output.obs["cell_type"].to_numpy() == "B")
    )[0]
    max_profile_total = sum(int(value) for value in _matrix_values(output.X, science)[max_profile_index])
    assert max_profile_total == science.np.iinfo(science.np.int64).max
    assert report.summary["parameters"]["summation_policy"] == {
        "output_dtype": "int64",
        "overflow_check": "exact_nonnegative_profile_total_before_backend",
        "backend_exact_integer_limit": 2**53,
        "backend_above_limit_check": "float64_rounding_equivalence_with_disclosure",
    }
    backend_verification = report.summary["key_results"]["backend_count_verification"]
    assert backend_verification["precision_limited_profile_count"] == 2
    assert backend_verification["backend_matrix_dtype"] == "float64"
    assert any("rounding-equivalent rather than exact" in warning for warning in report.summary["warnings"])

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk"](adata)
    reproduced_output, _ = validate_pseudobulk_artifact(reproduced)
    science.np.testing.assert_array_equal(
        _matrix_values(reproduced_output.X, science), _matrix_values(output.X, science)
    )
    assert reproduced_summary == report.summary

    overflowing = adata.copy()
    overflowing_values = selected.copy()
    overflowing_values[1] = [2**62, 0, 0, 0]
    overflowing.X = (
        science.sparse.csr_matrix(overflowing_values)
        if storage == "csr"
        else science.sparse.csc_matrix(overflowing_values)
        if storage == "csc"
        else overflowing_values
    )
    calls_before = len(backend.calls)
    with pytest.raises(OverflowError) as runtime_error:
        run_pseudobulk_owned(
            overflowing,
            sample_key="sample",
            population_key="cell_type",
            condition_key="condition",
            technical_batch_key="batch",
            categorical_covariate_keys="",
            continuous_covariate_keys="age",
            source={"source": "X"},
            inference_mode="formal",
            min_cells=1,
            min_counts=1,
        )
    assert len(backend.calls) == calls_before
    with pytest.raises(OverflowError) as generated_error:
        namespace["run_pseudobulk"](overflowing)
    assert str(generated_error.value) == str(runtime_error.value)
    assert len(backend.calls) == calls_before


def test_pseudobulk_named_counts_layer_is_explicit_and_exact(science, monkeypatch):
    adata, _, artifact, report, _ = _run_aggregation(
        science,
        monkeypatch,
        source={"source": "layer", "layer_name": "counts"},
    )
    output, metadata = validate_pseudobulk_artifact(artifact)
    assert metadata["count_source"] == {
        "kind": "layer",
        "layer_name": "counts",
        "label": "layers['counts']",
        "state": "raw_counts",
        "state_basis": "caller_declared_with_numeric_count_validation",
        "feature_axis": "adata.var_names",
        "history_provenance_required": False,
        "raw_equivalence_required": False,
    }
    assert report.summary["parameters"]["count_source"] == metadata["count_source"]
    assert output.var_names.equals(adata.raw.var_names)


def test_pseudobulk_uses_the_selected_raw_axis_and_rejects_duplicate_selected_identifiers(
    science,
    monkeypatch,
):
    backend = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", backend)

    duplicate_cells = _fixture_adata(science)
    duplicate_names = duplicate_cells.obs_names.tolist()
    duplicate_names[1] = duplicate_names[0]
    duplicate_cells.obs_names = duplicate_names
    with pytest.raises(ValueError, match="unique cell and gene identifiers"):
        run_pseudobulk_owned(
            duplicate_cells,
            population_key="cell_type",
            condition_key="condition",
            source={"source": "raw"},
            min_cells=1,
            min_counts=1,
        )

    duplicate_genes = _fixture_adata(science)
    duplicate_raw = duplicate_genes.raw.to_adata()
    duplicate_gene_names = duplicate_raw.var_names.tolist()
    duplicate_gene_names[1] = duplicate_gene_names[0]
    duplicate_raw.var_names = duplicate_gene_names
    duplicate_genes.raw = duplicate_raw
    with pytest.raises(ValueError, match="unique cell and gene identifiers"):
        run_pseudobulk_owned(
            duplicate_genes,
            population_key="cell_type",
            condition_key="condition",
            source={"source": "raw"},
            min_cells=1,
            min_counts=1,
        )

    selected_raw = _fixture_adata(science)
    raw_values = _matrix_values(selected_raw.raw.X, science)[:, [1, 0, 2]]
    selected_raw.raw = science.ad.AnnData(
        raw_values,
        obs=selected_raw.obs.copy(),
        var=science.pd.DataFrame(index=["RG2", "RG1", "RG3"]),
    )
    _, _, artifact, report, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=selected_raw,
        source={"source": "raw"},
        min_cells=1,
    )
    output, metadata = validate_pseudobulk_artifact(artifact)
    assert output.var_names.tolist() == ["RG2", "RG1", "RG3"]
    assert metadata["input_dimensions"]["genes"] == 3
    assert report.summary["parameters"]["count_source"]["feature_axis"] == "raw.var_names"


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda values: values.__setitem__((0, 0), 0.25), "integer-like"),
        (lambda values: values.__setitem__((0, 0), -1), "negative"),
        (lambda values: values.__setitem__((0, 0), float("nan")), "non-finite"),
    ],
)
def test_pseudobulk_rejects_invalid_raw_counts_before_backend(science, monkeypatch, mutation, message):
    adata = _fixture_adata(science)
    values = _matrix_values(adata.raw.X, science).astype(float)
    mutation(values)
    raw_source = adata.copy()
    raw_source.X = values
    adata.raw = raw_source
    backend = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", backend)
    with pytest.raises(ValueError, match=message):
        run_pseudobulk_owned(
            adata,
            source={"source": "raw"},
            population_key="cell_type",
            condition_key="condition",
            min_cells=1,
            min_counts=1,
        )
    assert backend.calls == []


def test_pseudobulk_rejects_sample_mapping_conflicts_and_string_collisions(science, monkeypatch):
    backend = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", backend)
    conflicting = _fixture_adata(science)
    conflicting.obs.iloc[1, conflicting.obs.columns.get_loc("condition")] = "treated"
    with pytest.raises(ValueError, match="one Condition/batch/covariate mapping.*s1"):
        run_pseudobulk_owned(
            conflicting,
            source={"source": "raw"},
            population_key="cell_type",
            condition_key="condition",
            min_cells=1,
            min_counts=1,
        )

    collision = _fixture_adata(science)
    collision.obs["sample"] = collision.obs["sample"].astype(object)
    collision.obs.iloc[0, collision.obs.columns.get_loc("sample")] = 1
    collision.obs.iloc[1, collision.obs.columns.get_loc("sample")] = "1"
    with pytest.raises(ValueError, match="collapse after string conversion"):
        run_pseudobulk_owned(
            collision,
            source={"source": "raw"},
            population_key="cell_type",
            condition_key="condition",
            min_cells=1,
            min_counts=1,
        )
    assert backend.calls == []


def test_pseudobulk_reports_qc_that_erases_a_population_condition_stratum_without_mutation(
    science,
    monkeypatch,
):
    adata = _fixture_adata(science)
    adata.obs.loc[(adata.obs["sample"] == "s5") & (adata.obs["cell_type"] == "B"), "cell_type"] = "T"
    original = adata.copy()
    backend = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", backend)
    artifact, report, _ = run_pseudobulk_owned(
        adata,
        sample_key="sample",
        population_key="cell_type",
        condition_key="condition",
        technical_batch_key="batch",
        source={"source": "raw"},
        min_cells=2,
        min_counts=1,
    )
    validate_pseudobulk_artifact(artifact)
    assert report.summary["key_results"]["lost_strata"] == [["B", "treated"]]
    assert any("removed an entire" in warning for warning in report.summary["warnings"])
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    science.np.testing.assert_array_equal(_matrix_values(adata.raw.X, science), _matrix_values(original.raw.X, science))


@pytest.mark.parametrize("malformed", ["counts", "identity", "qc", "genes", "filter"])
def test_pseudobulk_rejects_malformed_backend_without_mutating_input(science, monkeypatch, malformed):
    adata = _fixture_adata(science)
    original = adata.copy()
    backend = _fake_decoupler(science, malformed=malformed)
    monkeypatch.setitem(sys.modules, "decoupler", backend)
    with pytest.raises(RuntimeError, match="backend|profile-filter"):
        run_pseudobulk_owned(
            adata,
            source={"source": "raw"},
            population_key="cell_type",
            condition_key="condition",
            technical_batch_key="batch",
            min_cells=1,
            min_counts=1,
        )
    science.np.testing.assert_array_equal(_matrix_values(adata.raw.X, science), _matrix_values(original.raw.X, science))
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)


@pytest.mark.parametrize("malformed", ["counts", "filter"])
def test_pseudobulk_generated_code_matches_runtime_backend_failure_and_preserves_input(
    science,
    monkeypatch,
    malformed,
):
    adata = _fixture_adata(science)
    original = adata.copy()
    good = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", good)
    _, _, code = run_pseudobulk_owned(
        adata,
        sample_key="sample",
        population_key="cell_type",
        condition_key="condition",
        technical_batch_key="batch",
        continuous_covariate_keys="age",
        source={"source": "raw"},
        min_cells=2,
        min_counts=1,
    )
    bad = _fake_decoupler(science, malformed=malformed)
    monkeypatch.setitem(sys.modules, "decoupler", bad)
    with pytest.raises(RuntimeError) as runtime_error:
        run_pseudobulk_owned(
            adata,
            sample_key="sample",
            population_key="cell_type",
            condition_key="condition",
            technical_batch_key="batch",
            continuous_covariate_keys="age",
            source={"source": "raw"},
            min_cells=2,
            min_counts=1,
        )
    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError) as generated_error:
        namespace["run_pseudobulk"](adata)
    assert str(generated_error.value) == str(runtime_error.value)
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    science.np.testing.assert_array_equal(_matrix_values(adata.raw.X, science), _matrix_values(original.raw.X, science))


def test_pseudobulk_backend_version_and_public_signature_drift_fail_closed(science, monkeypatch):
    adata = _fixture_adata(science)
    wrong_version = _fake_decoupler(science)
    wrong_version.__version__ = "1.9.2"
    monkeypatch.setitem(sys.modules, "decoupler", wrong_version)
    with pytest.raises(RuntimeError, match="requires Decoupler 2.x.*1.9.2"):
        run_pseudobulk_owned(
            adata,
            population_key="cell_type",
            condition_key="condition",
            source={"source": "raw"},
            min_cells=1,
            min_counts=1,
        )

    wrong_signature = _fake_decoupler(science)
    wrong_signature.pp.filter_samples = lambda adata: adata.obs_names.to_numpy(dtype=str)
    monkeypatch.setitem(sys.modules, "decoupler", wrong_signature)
    with pytest.raises(RuntimeError, match="filter_samples interface is incompatible"):
        run_pseudobulk_owned(
            adata,
            population_key="cell_type",
            condition_key="condition",
            source={"source": "raw"},
            min_cells=1,
            min_counts=1,
        )


def test_pseudobulk_artifact_is_defensive_and_tamper_evident(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    first = artifact.to_adata()
    first.X[0, 0] += 100
    first.obs.iloc[0, first.obs.columns.get_loc("openbio_total_counts")] += 100
    second = artifact.to_adata()
    assert _matrix_values(second.X, science)[0, 0] != _matrix_values(first.X, science)[0, 0]
    forged = PseudobulkArtifact(first, artifact.metadata)
    with pytest.raises(ValueError, match="retained profile rows differ|content fingerprint changed"):
        validate_pseudobulk_artifact(forged)
    with pytest.raises(AttributeError, match="immutable"):
        artifact.new_field = 1


def test_shared_design_preflight_is_one_population_sample_level_ranked_and_directional(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    prepared = prepare_pseudobulk_design(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        continuous_covariate_keys="age",
        min_count=0,
        min_total_count=1,
    )
    diagnostics = prepared["diagnostics"]
    assert prepared["adata"].n_obs == 6
    assert diagnostics["condition_counts"] == {"control": 3, "treated": 3}
    assert diagnostics["design"]["columns"] == [
        "intercept",
        "batch[b2]",
        "age",
        "condition[treated vs control]",
    ]
    assert diagnostics["design"]["rank"] == 4
    assert diagnostics["design"]["residual_df"] == 2
    assert diagnostics["design"]["contrast_vector"] == [0.0, 0.0, 0.0, 1.0]
    assert "higher in comparison versus reference" in diagnostics["effect_direction"]
    assert diagnostics["expression_filter"]["genes_before"] == 4
    assert diagnostics["expression_filter"]["genes_retained"] == 4
    assert diagnostics["formal_interpretation_invalid"] is False
    assert diagnostics["formal_interpretation_invalid_reasons"] == []


def test_shared_design_allows_low_replication_but_rejects_exact_confounding(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch, min_cells=1)
    prepared = prepare_pseudobulk_design(
        artifact,
        population="B",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
    )
    assert prepared["diagnostics"]["condition_counts"] == {"control": 2, "treated": 2}
    assert prepared["diagnostics"]["formal_interpretation_invalid"] is True
    assert "condition_replication_below_three" in prepared["diagnostics"][
        "formal_interpretation_invalid_reasons"
    ]

    adata = _fixture_adata(science)
    adata.obs["batch"] = adata.obs["condition"].map({"control": "b1", "treated": "b2"})
    _, _, confounded, _, _ = _run_aggregation(science, monkeypatch, adata=adata)
    with pytest.raises(ValueError, match="rank deficient.*confounded"):
        prepare_pseudobulk_design(
            confounded,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
        )


def test_role_aliases_and_constant_nuisance_are_disclosed_without_preemptive_rejection(science, monkeypatch):
    adata = _fixture_adata(science)
    adata.obs["constant_batch"] = "only_batch"
    adata.obs["age"] = 42.0
    _, _, artifact, report, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        technical_batch_key="constant_batch",
        categorical_covariate_keys="constant_batch",
    )
    assert report.summary["key_results"]["formal_interpretation_invalid"] is True
    assert report.summary["key_results"]["role_aliases"] == {
        "constant_batch": ["technical_batch", "categorical_covariate"]
    }
    prepared = prepare_pseudobulk_design(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        categorical_covariate_keys="constant_batch",
        continuous_covariate_keys="age",
        min_count=0,
        min_total_count=1,
    )
    assert prepared["diagnostics"]["design"]["omitted_constant_terms"] == ["constant_batch", "age"]
    assert prepared["diagnostics"]["design"]["columns"] == [
        "intercept",
        "condition[treated vs control]",
    ]


def test_incomplete_but_full_rank_batch_overlap_warns_instead_of_blocking(science, monkeypatch):
    adata = _fixture_adata(science)
    batch_by_sample = {"s1": "b1", "s2": "b1", "s3": "b2", "s4": "b2", "s5": "b2", "s6": "b3"}
    adata.obs["batch"] = adata.obs["sample"].map(batch_by_sample)
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
    )
    prepared = prepare_pseudobulk_design(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
    )
    assert prepared["diagnostics"]["design"]["rank"] == 4
    assert prepared["diagnostics"]["formal_interpretation_invalid"] is True
    assert "incomplete_categorical_condition_overlap" in prepared["diagnostics"][
        "formal_interpretation_invalid_reasons"
    ]
    assert prepared["diagnostics"]["design"]["categorical_condition_overlap"][0]["key"] == "batch"
    backend = _fake_pertpy(science, engine="edger")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    _, report, _ = run_pseudobulk_edger_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
    )
    assert any("incomplete Condition overlap" in warning for warning in report.summary["warnings"])
    assert report.summary["key_results"]["formal_interpretation_invalid"] is True


def test_one_versus_one_samples_remains_a_zero_residual_hard_failure(science, monkeypatch):
    adata = _fixture_adata(science)
    keep_b_samples = {"s2", "s4"}
    remove = (adata.obs["cell_type"] == "B") & ~adata.obs["sample"].isin(keep_b_samples)
    adata = adata[~remove].copy()
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        min_cells=1,
        technical_batch_key="",
        continuous_covariate_keys="",
    )
    with pytest.raises(ValueError, match="no positive residual degrees of freedom.*Samples=2, rank=2"):
        prepare_pseudobulk_design(
            artifact,
            population="B",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
        )


def test_shared_design_rejects_full_rank_design_with_zero_residual_degrees_of_freedom(science, monkeypatch):
    adata = _fixture_adata(science)
    basis_samples = ["s1", "s2", "s4", "s5"]
    for index, sample in enumerate(basis_samples, start=1):
        adata.obs[f"c{index}"] = (adata.obs["sample"] == sample).astype(float)
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        technical_batch_key="",
        continuous_covariate_keys="c1,c2,c3,c4",
    )
    with pytest.raises(ValueError, match="no positive residual degrees of freedom.*Samples=6, rank=6"):
        prepare_pseudobulk_design(
            artifact,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            continuous_covariate_keys="c1,c2,c3,c4",
            min_count=0,
            min_total_count=1,
        )


def test_edger_full_table_summary_backend_controls_generated_parity_and_input_immutability(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    original, original_metadata = validate_pseudobulk_artifact(artifact)
    backend = _fake_pertpy(science, engine="edger")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    result, report, code = run_pseudobulk_edger_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        fdr_threshold=0.05,
        min_abs_log2_fold_change=1.0,
        min_count=0,
        min_total_count=1,
    )
    assert result.table.columns.tolist() == EDGER_COLUMNS
    assert len(result.table) == 4
    assert backend.calls[1] == ("fit", {"robust": True})
    science.np.testing.assert_array_equal(backend.calls[2][1], [0.0, 0.0, 1.0])
    assert report.summary["key_results"]["condition_counts"] == {"control": 3, "treated": 3}
    assert report.summary["key_results"]["significant_total"] == 2
    assert report.summary["key_results"]["significant_up"] == 1
    assert report.summary["key_results"]["significant_down"] == 1
    assert report.summary["parameters"]["fixed_policy"]["ql_fit"] == "glmQLFit robust=True"
    assert report.summary["key_results"]["formal_interpretation_invalid"] is False
    assert any("exactly three" in warning for warning in report.summary["warnings"])
    assert "comparison-minus-reference" in report.summary["methods"]
    assert set(report.summary["software_versions"]) >= {
        "python",
        "openbio-singlecell",
        "pertpy",
        "decoupler",
        "rpy2",
        "R",
        "edgeR",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
    }
    json.dumps(report.summary, allow_nan=False)
    assert "_standalone_validate_pydeseq2_dds_snapshot" not in code
    compile(code, "<edger-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk_edger"](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary
    after, after_metadata = validate_pseudobulk_artifact(artifact)
    science.np.testing.assert_array_equal(_matrix_values(after.X, science), _matrix_values(original.X, science))
    assert after_metadata == original_metadata


@pytest.mark.parametrize(
    "malformed",
    [
        "missing_gene",
        "nonfinite",
        "p_range",
        "wrong_bh",
        "schema",
        "backend_gene_mutation",
        "backend_count_mutation",
        "backend_fractional_count_mutation",
        "backend_design_mutation",
        "backend_contrast_mutation",
    ],
)
def test_edger_malformed_backend_fails_runtime_and_generated(science, monkeypatch, malformed):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    good = _fake_pertpy(science, engine="edger")
    monkeypatch.setitem(sys.modules, "pertpy", good)
    _, _, code = run_pseudobulk_edger_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
    )
    bad = _fake_pertpy(science, engine="edger", malformed=malformed)
    monkeypatch.setitem(sys.modules, "pertpy", bad)
    with pytest.raises(RuntimeError) as runtime_error:
        run_pseudobulk_edger_owned(
            artifact,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
        )
    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError) as generated_error:
        namespace["run_pseudobulk_edger"](artifact)
    assert str(generated_error.value) == str(runtime_error.value)


def test_pydeseq2_null_preservation_summary_controls_generated_parity_and_no_r_claim(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    backend = _fake_pertpy(science, engine="pydeseq2")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    result, report, code = run_pseudobulk_deseq2_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        fdr_threshold=0.05,
        min_abs_log2_fold_change=1.0,
        min_count=0,
        min_total_count=1,
        n_cpus=1,
    )
    assert result.table.columns.tolist() == PYDESEQ2_COLUMNS
    assert len(result.table) == 4
    assert result.table["p_value"].isna().sum() == 1
    assert result.table["p_adjusted"].isna().sum() == 2
    assert backend.calls[1] == (
        "fit",
        {
            "fit_type": "parametric",
            "size_factors_fit_type": "ratio",
            "control_genes": None,
            "min_mu": 0.5,
            "min_disp": 1e-8,
            "max_disp": 10.0,
            "min_replicates": 7,
            "beta_tol": 1e-8,
            "n_cpus": 1,
            "quiet": False,
            "low_memory": False,
        },
    )
    assert backend.calls[2][2] == {
        "alpha": 0.05,
        "lfc_shrink": None,
        "cooks_filter": True,
        "independent_filter": True,
        "prior_LFC_var": None,
        "lfc_null": 0.0,
        "alt_hypothesis": None,
        "quiet": False,
        "n_cpus": 1,
    }
    assert report.summary["key_results"]["null_p_value_count"] == 1
    assert report.summary["key_results"]["null_p_adjusted_count"] == 2
    assert report.summary["key_results"]["independent_filter_null_count"] == 1
    assert report.summary["key_results"]["significant_total"] == 2
    assert report.summary["parameters"]["fixed_policy"]["lfc_shrinkage"] is None
    assert report.summary["parameters"]["fixed_policy"]["control_genes"] is None
    assert report.summary["parameters"]["fixed_policy"]["min_replicates"] == 7
    assert report.summary["parameters"]["fixed_policy"]["lfc_null"] == 0.0
    assert report.summary["parameters"]["fixed_policy"]["alt_hypothesis"] is None
    assert report.summary["parameters"]["realized_policy"]["max_disp"] == 10.0
    assert report.summary["key_results"]["backend_input_postconditions"]["fitted_dds_count_values_preserved"] is True
    assert report.summary["key_results"]["backend_input_postconditions"]["fitted_dds_design_preserved"] is True
    assert report.summary["key_results"]["backend_input_postconditions"]["fitted_dds_fixed_controls_preserved"] is True
    assert report.summary["key_results"]["formal_interpretation_invalid"] is False
    assert any("exactly three" in warning for warning in report.summary["warnings"])
    assert "PyDESeq2" in report.summary["methods"]
    assert "R" not in report.summary["software_versions"]
    assert "DESeq2" not in report.summary["software_versions"]
    json.dumps(report.summary, allow_nan=False)
    assert "_standalone_validate_pydeseq2_dds_snapshot" in code
    compile(code, "<pydeseq2-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk_pydeseq2"](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary


def test_pydeseq2_minimal_backend_state_avoids_reserved_obs_var_collisions(science, monkeypatch):
    adata = _fixture_adata(science)
    adata.var["dispersions"] = science.np.linspace(1.0, 2.0, adata.n_vars)
    size_factors_by_sample = {sample: 1.0 + index / 10 for index, sample in enumerate(sorted(adata.obs["sample"]))}
    adata.obs["size_factors"] = adata.obs["sample"].map(size_factors_by_sample)
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        source={"source": "X"},
        min_cells=1,
        continuous_covariate_keys="age,size_factors",
    )
    artifact_data, _ = validate_pseudobulk_artifact(artifact)
    assert "dispersions" in artifact_data.var
    assert "size_factors" in artifact_data.obs

    backend = _fake_pertpy(science, engine="pydeseq2")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    result, report, code = run_pseudobulk_deseq2_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
        n_cpus=1,
    )
    backend_input = backend.calls[0][1]
    assert backend_input.obs.columns.tolist() == ["__openbio_sample_identity__"]
    assert backend_input.var.columns.tolist() == ["__openbio_feature_identity__"]

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk_pydeseq2"](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary


def test_pydeseq2_reports_realized_iterative_size_factor_fallback_with_generated_parity(science, monkeypatch):
    adata = _fixture_adata(science)
    selected = _matrix_values(adata.X, science).astype(science.np.int64, copy=True)
    for gene_index, sample in enumerate(("s1", "s2", "s3", "s4")):
        sample_cells = (adata.obs["sample"] == sample) & (adata.obs["cell_type"] == "T")
        selected[sample_cells.to_numpy(), gene_index] = 0
    adata.X = selected
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        source={"source": "X"},
        min_cells=1,
    )
    backend = _fake_pertpy(science, engine="pydeseq2")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    result, report, code = run_pseudobulk_deseq2_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
        n_cpus=1,
    )

    assert report.summary["parameters"]["requested_policy"]["size_factors_fit_type"] == "ratio"
    assert report.summary["parameters"]["realized_policy"]["size_factors_fit_type"] == "iterative"
    assert report.summary["key_results"]["every_retained_gene_contains_zero"] is True
    assert report.summary["key_results"]["realized_size_factors_fit_type"] == "iterative"
    assert "realized iterative size-factor fitting" in report.summary["methods"]
    assert any(
        "Switching to iterative mode" in record["message"]
        for record in report.summary["key_results"]["backend_warning_records"]
    )

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk_pydeseq2"](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary


@pytest.mark.parametrize(
    ("owned_operation", "engine", "function_name"),
    [
        (run_pseudobulk_edger_owned, "edger", "run_pseudobulk_edger"),
        (run_pseudobulk_deseq2_owned, "pydeseq2", "run_pseudobulk_pydeseq2"),
    ],
)
def test_engines_run_one_versus_two_samples_with_warning_and_generated_parity(
    science,
    monkeypatch,
    owned_operation,
    engine,
    function_name,
):
    adata = _fixture_adata(science)
    remove = (adata.obs["sample"] == "s1") & (adata.obs["cell_type"] == "B")
    adata = adata[~remove].copy()
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        min_cells=1,
        technical_batch_key="",
        continuous_covariate_keys="",
    )
    backend = _fake_pertpy(science, engine=engine)
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    kwargs = {
        "population": "B",
        "reference_condition": "control",
        "comparison_condition": "treated",
        "min_count": 0,
        "min_total_count": 1,
    }
    if engine == "pydeseq2":
        kwargs["n_cpus"] = 1
    result, report, code = owned_operation(artifact, **kwargs)
    assert report.summary["key_results"]["condition_counts"] == {"control": 1, "treated": 2}
    assert report.summary["key_results"]["design"]["residual_df"] == 1
    assert report.summary["key_results"]["formal_interpretation_invalid"] is True
    assert "condition_replication_below_three" in report.summary["key_results"][
        "formal_interpretation_invalid_reasons"
    ]
    assert any("one independent Sample" in warning for warning in report.summary["warnings"])
    if engine == "pydeseq2":
        assert any("Cook" in warning and "seven" in warning for warning in report.summary["warnings"])
        assert "cooks_refit_below_seven" not in report.summary["key_results"][
            "formal_interpretation_invalid_reasons"
        ]
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace[function_name](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary


@pytest.mark.parametrize(
    "malformed",
    [
        "missing_gene",
        "negative_base",
        "infinite",
        "wrong_bh",
        "schema",
        "backend_gene_mutation",
        "backend_count_mutation",
        "backend_fractional_count_mutation",
        "backend_design_mutation",
        "backend_contrast_mutation",
        "dds_count_mutation",
        "dds_large_integer_count_mutation",
        "dds_sample_axis_mutation",
        "dds_gene_axis_mutation",
        "dds_obs_metadata_mutation",
        "dds_var_metadata_mutation",
        "dds_design_mutation",
        "dds_fixed_control_mutation",
    ],
)
def test_pydeseq2_malformed_backend_fails_runtime_and_generated(science, monkeypatch, malformed):
    adata = _fixture_adata(science)
    adata.var["feature_class"] = ["coding", "coding", "control", "control"]
    if malformed == "dds_large_integer_count_mutation":
        selected = _matrix_values(adata.X, science).astype(science.np.int64, copy=True)
        selected[0, 0] = 2**62
        adata.X = selected
    _, _, artifact, _, _ = _run_aggregation(
        science,
        monkeypatch,
        adata=adata,
        source={"source": "X"},
    )
    good = _fake_pertpy(science, engine="pydeseq2")
    monkeypatch.setitem(sys.modules, "pertpy", good)
    _, _, code = run_pseudobulk_deseq2_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
        n_cpus=1,
    )
    bad = _fake_pertpy(science, engine="pydeseq2", malformed=malformed)
    monkeypatch.setitem(sys.modules, "pertpy", bad)
    with pytest.raises(RuntimeError) as runtime_error:
        run_pseudobulk_deseq2_owned(
            artifact,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
            n_cpus=1,
        )
    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError) as generated_error:
        namespace["run_pseudobulk_pydeseq2"](artifact)
    assert str(generated_error.value) == str(runtime_error.value)


def test_active_engines_reject_generic_anndata_before_backend(science, monkeypatch):
    generic = science.ad.AnnData(science.np.ones((6, 3), dtype=int))
    backend = _fake_pertpy(science, engine="edger")
    monkeypatch.setitem(sys.modules, "pertpy", backend)
    with pytest.raises(TypeError, match="OPENBIO_SINGLE_CELL_PSEUDOBULK"):
        run_pseudobulk_edger_owned(
            generic,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
        )
    assert backend.calls == []


def test_active_engines_reject_unaudited_pertpy_and_pydeseq2_versions(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)

    wrong_pertpy = _fake_pertpy(science, engine="edger")
    wrong_pertpy.__version__ = "1.4.0"
    monkeypatch.setitem(sys.modules, "pertpy", wrong_pertpy)
    with pytest.raises(RuntimeError, match="requires audited Pertpy >=1.3,<1.4"):
        run_pseudobulk_edger_owned(
            artifact,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
        )

    wrong_pydeseq2 = _fake_pertpy(science, engine="pydeseq2")
    wrong_pydeseq2._openbio_backend_versions["pydeseq2"] = "0.6.0"
    monkeypatch.setitem(sys.modules, "pertpy", wrong_pydeseq2)
    with pytest.raises(RuntimeError, match="requires audited pydeseq2 >=0.5,<0.6"):
        run_pseudobulk_deseq2_owned(
            artifact,
            population="T",
            reference_condition="control",
            comparison_condition="treated",
            min_count=0,
            min_total_count=1,
            n_cpus=1,
        )


def test_pseudobulk_rejects_non_scalar_role_metadata_before_backend(science, monkeypatch):
    adata = _fixture_adata(science)
    adata.obs["sample"] = adata.obs["sample"].astype(object)
    adata.obs.at[adata.obs_names[0], "sample"] = ["not", "a", "scalar"]
    backend = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", backend)

    with pytest.raises(TypeError, match="Sample labels value at position 0 must be scalar"):
        run_pseudobulk_owned(
            adata,
            sample_key="sample",
            population_key="cell_type",
            condition_key="condition",
            technical_batch_key="batch",
            source={"source": "layer", "layer_name": "counts"},
            min_cells=1,
            min_counts=1,
        )

    assert backend.calls == []


def test_pydeseq2_reports_requested_and_realized_dispersion_and_captured_warnings(science, monkeypatch):
    _, _, artifact, _, _ = _run_aggregation(science, monkeypatch)
    backend = _fake_pertpy(science, engine="pydeseq2", malformed="dispersion_fallback")
    monkeypatch.setitem(sys.modules, "pertpy", backend)

    result, report, code = run_pseudobulk_deseq2_owned(
        artifact,
        population="T",
        reference_condition="control",
        comparison_condition="treated",
        min_count=0,
        min_total_count=1,
        n_cpus=1,
    )

    assert report.summary["parameters"]["requested_policy"]["dispersion_fit_type"] == "parametric"
    assert report.summary["parameters"]["realized_policy"]["dispersion_fit_type"] == "mean"
    assert report.summary["key_results"]["realized_dispersion_fit_type"] == "mean"
    assert report.summary["key_results"]["backend_warning_records"] == [
        {
            "category": "UserWarning",
            "message": (
                "The dispersion trend curve fitting did not converge. "
                "Switching to a mean-based dispersion trend."
            ),
            "count": 1,
        }
    ]
    assert "realized mean trend" in report.summary["methods"]
    assert any("mean-based dispersion trend" in warning for warning in report.summary["warnings"])

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_pseudobulk_pydeseq2"](artifact)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary
