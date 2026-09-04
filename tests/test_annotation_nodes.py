from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
import time
import types
from dataclasses import replace

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.annotation_core import (
    _celltypist_axis_fingerprint,
    _celltypist_probability_fingerprint,
    _celltypist_selected_fingerprint,
)
from openbio_singlecell.files import resolve_input_path
from openbio_singlecell.marker_evidence import MARKER_COLUMNS, MARKER_UNIVERSE_COLUMNS
from openbio_singlecell.node_types import AnnDataType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_annotation import (
    ANNOTATION_NODE_CLASSES,
    OpenBioSingleCellCellTypistAnnotation,
    OpenBioSingleCellMapClusterAnnotations,
    OpenBioSingleCellMarkerORAEvidence,
)
from openbio_singlecell.operations_annotation import (
    celltypist_owned,
    map_cluster_annotations_owned,
    marker_ora_owned,
)
from openbio_singlecell.ora_evidence import ORA_EVIDENCE_COLUMNS


def _marker_ora_owned(table, universe, *, resource_csv: str = "", **kwargs):
    return marker_ora_owned(
        table,
        universe,
        resource_path=resolve_input_path(resource_csv, extensions=(".csv",)),
        requested_resource_path=resource_csv.strip(),
        **kwargs,
    )


def _matrix_values(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _annotation_adata(science, *, cells=60, logged=False):
    counts = science.np.tile(science.np.array([1.0, 2.0, 3.0, 4.0]), (cells, 1))
    matrix = science.sparse.csr_matrix(counts)
    obs = science.pd.DataFrame(
        {
            "over": science.pd.Categorical(
                ["g1"] * (cells // 2) + ["g2"] * (cells - cells // 2),
                categories=["g2", "g1", "unused"],
            )
        },
        index=[f"cell_{index}" for index in range(cells)],
    )
    var = science.pd.DataFrame(index=["G1", "G2", "G3", "G4"])
    adata = science.ad.AnnData(matrix.copy(), obs=obs, var=var)
    if logged:
        science.sc.pp.normalize_total(adata, target_sum=10_000.0, inplace=True)
        science.sc.pp.log1p(adata)
    return adata


def _fake_celltypist(science, model_path, *, malformed=None, model_features=None, model_classes=None):
    module = types.ModuleType("celltypist")
    module.calls = []
    module.loads = []
    features = model_features or ["G1", "G2", "G3", "G4"]
    classes = model_classes or ["A", "B"]

    class Model:
        @staticmethod
        def load(path):
            module.loads.append(path)
            return types.SimpleNamespace(
                cell_types=science.np.asarray(classes),
                features=science.np.asarray(features),
                description={
                    "date": "2026-01-02",
                    "details": "deterministic fake",
                    "source": "test",
                    "version": "v-test",
                },
            )

    module.models = types.SimpleNamespace(
        Model=Model,
        get_all_models=lambda: ["fake.pkl"],
        get_model_path=lambda name: str(model_path),
    )

    def annotate(filename, **kwargs):
        matrix = _matrix_values(filename.X, science)
        module.calls.append(
            {
                **kwargs,
                "inverse_totals": science.np.expm1(matrix).sum(axis=1),
                "obs_names": filename.obs_names.copy(),
                "var_names": filename.var_names.copy(),
            }
        )
        n_obs = filename.n_obs
        probabilities = science.np.full((n_obs, len(classes)), 0.05, dtype=float)
        probabilities[:, :2] = [0.8, 0.2]
        first_half = n_obs // 2
        probabilities[max(first_half - 10, 0) : first_half, :2] = [0.1, 0.9]
        probabilities[first_half : min(first_half + 10, n_obs), :2] = [0.8, 0.2]
        probabilities[min(first_half + 10, n_obs) :, :2] = [0.1, 0.9]
        individual = science.np.where(probabilities[:, 0] >= probabilities[:, 1], "A", "B")
        predicted = science.pd.DataFrame(
            {"predicted_labels": science.pd.Categorical(individual)},
            index=filename.obs_names,
        )
        over_clustering = kwargs["over_clustering"]
        if kwargs["majority_voting"]:
            over = science.np.asarray(over_clustering, dtype=object)
            majority = science.np.where(over == "g1", "A", "B")
            predicted["over_clustering"] = over
            predicted["majority_voting"] = science.pd.Categorical(majority)
        probability = science.pd.DataFrame(
            probabilities,
            index=filename.obs_names,
            columns=classes,
        )
        decision = science.pd.DataFrame(
            science.np.log(probabilities / (1.0 - probabilities)),
            index=filename.obs_names,
            columns=classes,
        )
        if malformed == "index":
            probability.index = probability.index[::-1]
        elif malformed == "class_order":
            probability = probability[["B", "A"]]
        elif malformed == "range":
            probability.iloc[0, 0] = 1.1
        elif malformed == "labels":
            predicted["predicted_labels"] = predicted["predicted_labels"].astype(object)
            predicted.iloc[0, predicted.columns.get_loc("predicted_labels")] = "unknown"
        elif malformed == "majority":
            predicted.iloc[0, predicted.columns.get_loc("majority_voting")] = "B"
        return types.SimpleNamespace(
            predicted_labels=predicted,
            probability_matrix=probability,
            decision_matrix=decision,
        )

    module.annotate = annotate
    return module


def _mapping_adata(science):
    obs = science.pd.DataFrame(
        {
            "leiden": science.pd.Categorical(
                ["0", "1", "2", "0", None],
                categories=["2", "0", "1", "unused"],
                ordered=True,
            ),
            "unrelated": [1, 2, 3, 4, 5],
        },
        index=[f"cell_{index}" for index in range(5)],
    )
    adata = science.ad.AnnData(
        science.np.arange(15, dtype=float).reshape(5, 3),
        obs=obs,
        var=science.pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["X_umap"] = science.np.arange(10, dtype=float).reshape(5, 2)
    adata.obsp["connectivities"] = science.sparse.eye(5, format="csr")
    return adata


def _frame_content_fingerprint(frame, *, artifact, numeric_columns):
    rows = []
    for row in frame.itertuples(index=False, name=None):
        rows.append(
            [
                float(value).hex() if column in numeric_columns else value
                for column, value in zip(frame.columns, row, strict=True)
            ]
        )
    payload = {
        "schema": f"openbio-singlecell/{artifact}/v2",
        "columns": list(frame.columns),
        "rows": rows,
    }
    return hashlib.sha256(
        json.dumps(
            payload,
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()


def _marker_artifacts(science, *, ranking_truncated=False):
    genes = [f"G{index}" for index in range(1, 9)]
    universe_fingerprint = hashlib.sha256(
        json.dumps(
            {
                "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
                "ordered_genes": genes,
            },
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    analysis_fingerprint = "a" * 64
    ranking_fingerprint = "b" * 64
    shared = {
        "marker_evidence_schema_version": 2,
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": ranking_fingerprint,
        "universe_fingerprint": universe_fingerprint,
        "marker_groupby": "cluster",
        "marker_method": "wilcoxon",
        "marker_source": "X",
        "group_labels": ["A", "B"],
        "group_sizes": {"A": 10, "B": 10},
        "source_gene_count": len(genes),
        "requested_n_genes": 0,
        "actual_n_genes_per_group": len(genes),
        "ranking_truncated": ranking_truncated,
        "tie_correct": True,
        "max_output_rows": 10_000,
        "max_working_memory_gib": 1.0,
        "constant_gene_count": 0,
        "variable_gene_count": len(genes),
        "bh_hypotheses_per_group": len(genes),
        "marker_output_rows": 2 * len(genes),
        "universe_output_rows": len(genes),
        "total_output_rows": 3 * len(genes),
        "groups": "all",
        "reference": "rest",
        "rankby_abs": False,
        "pts": True,
        "corr_method": "benjamini-hochberg",
        "expression_state": "logged",
        "expression_interpretation": "verified logarithmized abundance",
    }
    rows = []
    for group, selected in (("A", ["G1", "G2", "G3"]), ("B", ["G4", "G5", "G6"])):
        for rank, gene in enumerate(selected, start=1):
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "rank": rank,
                    "score": 5.0 - rank,
                    "log2_fold_change_approx": 2.0,
                    "p_value": 0.001 * rank,
                    "p_adjusted": 0.01 * rank,
                    "fraction_in_group": 0.8,
                    "fraction_reference": 0.1,
                }
            )
    marker_frame = science.pd.DataFrame(rows, columns=MARKER_COLUMNS)
    universe_frame = science.pd.DataFrame(
        {"gene": genes, "universe_rank": list(range(1, len(genes) + 1))},
        columns=MARKER_UNIVERSE_COLUMNS,
    )
    marker_parameters = {
        **shared,
        "producer_node_id": "OpenBioSingleCellFilterMarkerGenes",
        "upstream_producer_node_id": "OpenBioSingleCellMarkerGenes",
        "artifact_role": "marker_table",
        "content_fingerprint": _frame_content_fingerprint(
            marker_frame,
            artifact="marker-table",
            numeric_columns=set(MARKER_COLUMNS[2:]),
        ),
        "upstream_content_fingerprint": "c" * 64,
        "min_log2_fold_change": 0.5,
        "min_fraction_in_group": 0.2,
        "max_fraction_reference": 0.5,
        "max_p_adjusted": 0.05,
        "comparison_semantics": "inclusive",
    }
    universe_parameters = {
        **shared,
        "producer_node_id": "OpenBioSingleCellMarkerGenes",
        "artifact_role": "tested_gene_universe",
        "content_fingerprint": _frame_content_fingerprint(
            universe_frame,
            artifact="tested-gene-universe",
            numeric_columns={"universe_rank"},
        ),
    }
    started_at = time.perf_counter()
    marker = make_table_result(
        table=marker_frame,
        title="filtered markers",
        operation="filter_marker_genes",
        parameters=marker_parameters,
        description="fixture",
        warnings=[],
        input_cells=20,
        input_genes=len(genes),
        started_at=started_at,
    )
    universe = make_table_result(
        table=universe_frame,
        title="universe",
        operation="marker_genes",
        parameters=universe_parameters,
        description="fixture",
        warnings=[],
        input_cells=20,
        input_genes=len(genes),
        started_at=started_at,
    )
    return marker, universe


def _replace_filtered_marker_table(marker, frame):
    parameters = {
        **marker.parameters,
        "content_fingerprint": _frame_content_fingerprint(
            frame,
            artifact="marker-table",
            numeric_columns=set(MARKER_COLUMNS[2:]),
        ),
    }
    return make_table_result(
        table=frame,
        title=marker.title,
        operation="filter_marker_genes",
        parameters=parameters,
        description=marker.description,
        warnings=list(marker.warnings),
        input_cells=marker.input_cells,
        input_genes=marker.input_genes,
        started_at=time.perf_counter(),
    )


def _resource_metadata():
    return json.dumps(
        {
            "name": "Test marker sets",
            "version": "1.0",
            "date": "2026-08-28",
            "organism": "human",
            "identifier_namespace": "gene_symbol",
            "scope": "unit-test cell identities",
            "license": "CC0-1.0",
            "citation": "Synthetic resource for deterministic tests",
        },
        separators=(",", ":"),
    )


def _write_ora_resource(input_dir):
    path = input_dir / "ora_resource.csv"
    path.write_text(
        "source,target\n"
        "TypeA,G1\n"
        "TypeA,G2\n"
        "TypeA,G7\n"
        "TypeA,G8\n"
        "TypeA,G1\n"
        "TypeA_tie,G1\n"
        "TypeA_tie,G2\n"
        "TypeA_tie,G7\n"
        "TypeA_tie,G8\n"
        "TypeB,G4\n"
        "TypeB,G5\n"
        "TypeB,G7\n"
        "TypeB,G8\n"
        "TooSmall,G3\n"
        "Outside,NOT_IN_UNIVERSE\n",
        encoding="utf-8",
    )
    return path


def _fake_decoupler(science, *, malformed=None, version="2.2.0"):
    from scipy import stats

    module = types.ModuleType("decoupler")
    module.__version__ = version
    module.calls = []

    def query_set(features, net, alternative, n_bg, ha_corr, tmin, verbose):
        module.calls.append(
            {
                "features": list(features),
                "net": net.copy(),
                "alternative": alternative,
                "n_bg": n_bg,
                "ha_corr": ha_corr,
                "tmin": tmin,
                "verbose": verbose,
            }
        )
        selected = set(features)
        rows = []
        for source in net["source"].drop_duplicates().tolist():
            targets = set(net.loc[net["source"] == source, "target"].tolist())
            a = len(selected & targets)
            b = len(targets - selected)
            c = len(selected - targets)
            d = n_bg - a - b - c
            statistic = float(science.np.log(((a + ha_corr) * (d + ha_corr)) / ((b + ha_corr) * (c + ha_corr))))
            p_value = float(stats.fisher_exact([[a, b], [c, d]], alternative=alternative).pvalue)
            rows.append({"source": source, "stat": statistic, "pval": p_value})
        result = science.pd.DataFrame(rows)
        result["padj"] = stats.false_discovery_control(result["pval"].to_numpy(), method="bh")
        if malformed == "missing_source":
            result = result.iloc[:-1].copy()
        elif malformed == "stat":
            result.loc[0, "stat"] += 1.0
        elif malformed == "pval":
            result.loc[0, "pval"] = 0.9
        elif malformed == "padj":
            result.loc[0, "padj"] = 0.9
        elif malformed == "nonfinite":
            result.loc[0, "stat"] = science.np.nan
        elif malformed == "dtype":
            result["stat"] = result["stat"].astype(str)
        return result.sort_values(["padj", "pval"], kind="mergesort").reset_index(drop=True)

    module.mt = types.SimpleNamespace(query_set=query_set)
    return module


def test_annotation_schemas_expose_explicit_semantics_summary_and_code():
    celltypist = OpenBioSingleCellCellTypistAnnotation.define_schema()
    mapping = OpenBioSingleCellMapClusterAnnotations.define_schema()

    assert [(item.display_name, item.io_type) for item in celltypist.outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    cell_inputs = {item.id: item for item in celltypist.inputs}
    assert cell_inputs["majority_voting"].default is False
    assert cell_inputs["expression_state"].default == "verified_cp10k_log1p"
    assert [option.key for option in cell_inputs["source"].options] == ["X", "raw", "layer"]
    assert all(
        "execute" not in node_class.__dict__
        for node_class in (
            OpenBioSingleCellCellTypistAnnotation,
            OpenBioSingleCellMapClusterAnnotations,
            OpenBioSingleCellMarkerORAEvidence,
        )
    )

    assert [item.display_name for item in mapping.outputs] == ["adata", "summary", "code"]
    map_inputs = {item.id: item for item in mapping.inputs}
    assert map_inputs["unmapped_policy"].default == "error"
    assert map_inputs["annotation_status"].default == "provisional"
    assert map_inputs["overwrite_existing"].default is False


@pytest.mark.parametrize("logged", [False, True])
def test_celltypist_transforms_counts_once_and_preserves_verified_logged_input(science, tmp_path, monkeypatch, logged):
    adata = _annotation_adata(science, logged=logged)
    original = _matrix_values(adata.X, science).copy()
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    output, report, code = celltypist_owned(
        adata.copy(),
        source={"source": "X"},
        expression_state="verified_cp10k_log1p" if logged else "verified_counts",
        model=str(model_path),
        store_decision_matrix=True,
    )

    science.np.testing.assert_allclose(fake.calls[0]["inverse_totals"], 10_000.0, rtol=5e-3)
    science.np.testing.assert_array_equal(_matrix_values(adata.X, science), original)
    assert "celltypist_cell_type" not in adata.obs
    assert output.obsm["celltypist_probabilities"].shape == (60, 2)
    assert output.obsm["celltypist_decision_scores"].shape == (60, 2)
    assert output.uns["celltypist"]["annotation_status"] == "provisional"
    assert output.uns["celltypist"]["classes"] == ["A", "B"]
    provenance = output.uns["celltypist"]
    assert provenance["observation_axis_fingerprint_sha256"] == _celltypist_axis_fingerprint(output.obs_names)
    assert provenance["probability_content_fingerprint_sha256"] == _celltypist_probability_fingerprint(
        output.obsm["celltypist_probabilities"],
        observation_ids=output.obs_names,
        classes=["A", "B"],
    )
    assert provenance["selected_annotation_fingerprint_sha256"] == _celltypist_selected_fingerprint(
        output.obs["celltypist_cell_type"].astype(str),
        output.obs["celltypist_confidence"],
        observation_ids=output.obs_names,
    )
    assert report.summary["key_results"]["expression_transform"] == (
        "none" if logged else "normalize_total_10000_then_log1p"
    )
    assert "Provisional annotation" in report.summary["results"]
    assert "ground truth" in " ".join(report.summary["limitations"])
    assert set(report.summary["software_versions"]) == {
        "python",
        "openbio-singlecell",
        "celltypist",
        "scikit-learn",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
    }
    assert any("pickle artifacts" in warning for warning in report.summary["warnings"])
    assert any("does not prove safety" in warning for warning in report.summary["warnings"])
    assert report.summary["key_results"]["model"]["resolved_path"] == str(model_path.resolve())
    assert len(report.summary["key_results"]["model"]["sha256"]) == 64
    assert output.uns["celltypist"]["software_versions"] == report.summary["software_versions"]
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<celltypist-code>", "exec")

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_celltypist_annotation"](adata)
    science.pd.testing.assert_series_equal(reproduced.obs["celltypist_cell_type"], output.obs["celltypist_cell_type"])
    science.np.testing.assert_allclose(
        reproduced.obsm["celltypist_probabilities"],
        output.obsm["celltypist_probabilities"],
    )
    assert reproduced.uns["celltypist"] == output.uns["celltypist"]
    assert reproduced_summary == report.summary


def test_celltypist_named_model_resolution_never_enumerates_or_downloads(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)

    def forbidden_network_side_effect(*_args, **_kwargs):
        raise AssertionError("CellTypist model resolution must not enumerate or download models")

    fake.models.get_all_models = forbidden_network_side_effect
    fake.models.download_if_required = forbidden_network_side_effect
    fake.models.download_models = forbidden_network_side_effect
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    output, report, code = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model="fake.pkl",
    )

    assert fake.loads == [str(model_path.resolve())]
    assert output.uns["celltypist"]["model"]["resolved_path"] == str(model_path.resolve())
    assert report.summary["parameters"]["model"] == "fake.pkl"


def test_celltypist_model_fingerprint_tracks_explicit_and_named_local_artifact(science, tmp_path, monkeypatch):
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model v1")
    fake = _fake_celltypist(science, model_path)

    def forbidden_network_side_effect(*_args, **_kwargs):
        raise AssertionError("CellTypist fingerprinting must not enumerate or download models")

    fake.models.get_all_models = forbidden_network_side_effect
    fake.models.download_if_required = forbidden_network_side_effect
    fake.models.download_models = forbidden_network_side_effect
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    explicit = OpenBioSingleCellCellTypistAnnotation.fingerprint_inputs(str(model_path))
    named = OpenBioSingleCellCellTypistAnnotation.fingerprint_inputs("fake.pkl")
    assert explicit == named
    assert explicit[0] == "openbio-celltypist-model-v1"
    assert explicit[1] == str(model_path.resolve()).lower()
    assert explicit[2] == model_path.stat().st_size
    assert explicit[3] == model_path.stat().st_mtime_ns
    assert len(explicit[4]) == 64

    model_path.write_bytes(b"trusted fake model v2")
    changed = OpenBioSingleCellCellTypistAnnotation.fingerprint_inputs("fake.pkl")
    assert changed != named
    assert changed[4] != named[4]


def test_celltypist_missing_model_fingerprint_is_stable_and_execution_is_actionable(science, tmp_path, monkeypatch):
    missing_path = tmp_path / "missing.pkl"
    fake = _fake_celltypist(science, missing_path)

    def forbidden_network_side_effect(*_args, **_kwargs):
        raise AssertionError("CellTypist node must not enumerate or download models")

    fake.models.get_all_models = forbidden_network_side_effect
    fake.models.download_if_required = forbidden_network_side_effect
    fake.models.download_models = forbidden_network_side_effect
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    first = OpenBioSingleCellCellTypistAnnotation.fingerprint_inputs("missing.pkl")
    second = OpenBioSingleCellCellTypistAnnotation.fingerprint_inputs("missing.pkl")
    assert first == second == ("openbio-celltypist-model-v1", "missing", "missing.pkl")
    with pytest.raises(FileNotFoundError) as error:
        celltypist_owned(
            _annotation_adata(science),
            expression_state="verified_counts",
            model="missing.pkl",
        )
    message = str(error.value)
    assert "never downloads models" in message
    assert "network-enabled" in message
    assert "celltypist.models.download_models(model='missing.pkl')" in message
    assert "trusted local .pkl path" in message
    assert fake.loads == []
    assert fake.calls == []


def test_celltypist_model_load_failure_discloses_pickle_trust_boundary(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    original_obs = adata.obs.copy(deep=True)
    model_path = tmp_path / "unloadable.pkl"
    model_path.write_bytes(b"not a model")
    fake = _fake_celltypist(science, model_path)

    def fail_load(_path):
        raise ValueError("malformed pickle")

    fake.models.Model.load = staticmethod(fail_load)
    monkeypatch.setitem(sys.modules, "celltypist", fake)
    with pytest.raises(RuntimeError) as error:
        celltypist_owned(
            adata.copy(),
            expression_state="verified_counts",
            model=str(model_path),
        )
    message = str(error.value)
    assert "failed to load model" in message
    assert "pickle artifacts that can execute code" in message
    assert "official or explicitly user-trusted sources" in message
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert "celltypist" not in adata.uns
    assert fake.calls == []


def test_celltypist_rejects_contradictory_expression_state_before_backend(science, tmp_path, monkeypatch):
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    logged = _annotation_adata(science, logged=True)
    with pytest.raises(ValueError, match="verified_counts.*non-integer"):
        celltypist_owned(
            logged,
            expression_state="verified_counts",
            model=str(model_path),
        )
    counts = _annotation_adata(science, logged=False)
    with pytest.raises(ValueError, match="verified_cp10k_log1p"):
        celltypist_owned(
            counts,
            expression_state="verified_cp10k_log1p",
            model=str(model_path),
        )
    assert fake.calls == []
    assert fake.loads == []


def test_celltypist_cp10k_validation_uses_official_absolute_one_count_boundary(science, tmp_path, monkeypatch):
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    invalid = _annotation_adata(science)
    invalid.X = invalid.X.astype(science.np.float32)
    science.sc.pp.normalize_total(invalid, target_sum=9_950.0, inplace=True)
    science.sc.pp.log1p(invalid)
    invalid_values = _matrix_values(invalid.X, science).copy()
    with pytest.raises(ValueError, match="every cell must be within 1 count of 10,000"):
        celltypist_owned(
            invalid,
            expression_state="verified_cp10k_log1p",
            model=str(model_path),
        )
    science.np.testing.assert_array_equal(_matrix_values(invalid.X, science), invalid_values)
    assert "celltypist" not in invalid.uns
    assert fake.loads == []
    assert fake.calls == []

    valid = _annotation_adata(science)
    valid.X = valid.X.astype(science.np.float32)
    science.sc.pp.normalize_total(valid, target_sum=10_000.0, inplace=True)
    science.sc.pp.log1p(valid)
    output, report, _ = celltypist_owned(
        valid,
        expression_state="verified_cp10k_log1p",
        model=str(model_path),
    )
    assert report.summary["parameters"]["target_sum"] is None
    assert report.summary["key_results"]["expression_transform"] == "none"
    assert output.uns["celltypist"]["expression"]["target_sum"] == 10_000.0


def test_celltypist_feature_overlap_model_provenance_and_collision_preflight(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path, model_features=["G2", "absent"])
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    output, report, _ = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model=str(model_path),
    )
    overlap = report.summary["key_results"]["feature_overlap"]
    assert overlap == {
        "matched_features": 1,
        "query_features": 4,
        "model_features": 2,
        "query_fraction": 0.25,
        "model_fraction": 0.5,
    }
    assert len(output.uns["celltypist"]["model"]["sha256"]) == 64
    assert output.uns["celltypist"]["model"]["description"]["version"] == "v-test"

    no_overlap = _fake_celltypist(science, model_path, model_features=["absent"])
    monkeypatch.setitem(sys.modules, "celltypist", no_overlap)
    with pytest.raises(ValueError, match="zero feature overlap"):
        celltypist_owned(
            adata.copy(),
            expression_state="verified_counts",
            model=str(model_path),
        )
    assert no_overlap.calls == []

    collision = adata.copy()
    collision.obs["celltypist_cell_type"] = "prior"
    no_overlap.loads.clear()
    with pytest.raises(ValueError, match="output keys already exist"):
        celltypist_owned(
            collision,
            expression_state="verified_counts",
            model=str(model_path),
        )
    assert no_overlap.loads == []


def test_celltypist_complete_label_count_maps_include_unobserved_model_classes(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    model_path = tmp_path / "three-classes.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path, model_classes=["A", "B", "C"])
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    output, report, code = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model=str(model_path),
    )
    key_results = report.summary["key_results"]
    assert key_results["class_order"] == ["A", "B", "C"]
    assert key_results["individual_label_counts"] == {"A": 30, "B": 30, "C": 0}
    assert key_results["selected_label_counts"] == {"A": 30, "B": 30, "C": 0}
    assert output.obs["celltypist_cell_type"].cat.categories.tolist() == ["A", "B", "C"]

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_celltypist_annotation"](adata)
    assert reproduced.obs["celltypist_cell_type"].cat.categories.tolist() == ["A", "B", "C"]
    assert reproduced_summary == report.summary


def test_celltypist_majority_labels_support_and_selected_probability_are_distinct(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    output, report, code = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model=str(model_path),
        majority_voting=True,
        over_clustering_key="over",
        min_prop=0.6,
    )

    # cell_20 is individually B in a majority-A group: selected probability must look up A (0.1), not row max (0.9).
    assert output.obs.loc["cell_20", "celltypist_individual_label"] == "B"
    assert output.obs.loc["cell_20", "celltypist_majority_label"] == "A"
    assert output.obs.loc["cell_20", "celltypist_individual_probability"] == pytest.approx(0.9)
    assert output.obs.loc["cell_20", "celltypist_confidence"] == pytest.approx(0.1)
    assert output.obs.loc["cell_20", "celltypist_majority_support"] == pytest.approx(2 / 3)
    key_results = report.summary["key_results"]
    assert key_results["majority_voting_completed"] is True
    assert key_results["group_order_preview"] == ["g2", "g1"]
    assert key_results["individual_majority_agreement"] == pytest.approx(2 / 3)

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_celltypist_annotation"](adata)
    science.pd.testing.assert_series_equal(
        reproduced.obs["celltypist_majority_label"],
        output.obs["celltypist_majority_label"],
    )
    science.np.testing.assert_allclose(
        reproduced.obs["celltypist_majority_support"],
        output.obs["celltypist_majority_support"],
    )
    assert reproduced_summary == report.summary


def test_celltypist_overwrite_removes_only_verified_stale_optional_artifacts(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    first, _, _ = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model=str(model_path),
        majority_voting=True,
        over_clustering_key="over",
        min_prop=0.6,
        store_decision_matrix=True,
    )
    first_obs = first.obs.copy(deep=True)
    first_uns = copy.deepcopy(first.uns)
    first_obsm = {key: science.np.asarray(value).copy() for key, value in first.obsm.items() if key is not None}
    first_x = _matrix_values(first.X, science).copy()

    second, report, code = celltypist_owned(
        first.copy(),
        expression_state="verified_counts",
        model=str(model_path),
        majority_voting=False,
        label_column="celltypist_cell_type_v2",
        confidence_column="celltypist_confidence_v2",
        probability_key="celltypist_probabilities_v2",
        store_decision_matrix=False,
        decision_key="celltypist_decision_scores_v2",
        overwrite_existing=True,
    )
    assert "celltypist_cell_type" not in second.obs
    assert "celltypist_confidence" not in second.obs
    assert "celltypist_majority_label" not in second.obs
    assert "celltypist_majority_support" not in second.obs
    assert "celltypist_cell_type_v2" in second.obs
    assert "celltypist_confidence_v2" in second.obs
    assert "celltypist_probabilities" not in second.obsm
    assert "celltypist_probabilities_v2" in second.obsm
    assert "celltypist_decision_scores" not in second.obsm
    assert second.uns["celltypist"]["majority_voting"] is False
    assert second.uns["celltypist"]["matrices"]["decision_key"] is None
    annotations = second.uns["openbio_singlecell"]["annotations"]
    assert "celltypist_cell_type" not in annotations
    assert annotations["celltypist_cell_type_v2"] == second.uns["celltypist"]
    replacement = report.summary["key_results"]["artifact_replacement"]
    assert replacement["prior_artifact_verified"] is True
    assert replacement["removed_stale_obs_keys"] == [
        "celltypist_cell_type",
        "celltypist_confidence",
        "celltypist_majority_label",
        "celltypist_majority_support",
    ]
    assert replacement["removed_stale_obsm_keys"] == [
        "celltypist_probabilities",
        "celltypist_decision_scores",
    ]
    assert replacement["removed_stale_annotation_keys"] == ["celltypist_cell_type"]

    science.pd.testing.assert_frame_equal(first.obs, first_obs)
    assert first.uns == first_uns
    science.np.testing.assert_array_equal(_matrix_values(first.X, science), first_x)
    assert set(key for key in first.obsm if key is not None) == set(first_obsm)
    for key, expected in first_obsm.items():
        science.np.testing.assert_array_equal(science.np.asarray(first.obsm[key]), expected)

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_celltypist_annotation"](first.copy())
    assert "celltypist_cell_type" not in reproduced.obs
    assert "celltypist_confidence" not in reproduced.obs
    assert "celltypist_majority_label" not in reproduced.obs
    assert "celltypist_majority_support" not in reproduced.obs
    assert "celltypist_probabilities" not in reproduced.obsm
    assert "celltypist_decision_scores" not in reproduced.obsm
    assert "celltypist_cell_type" not in reproduced.uns["openbio_singlecell"]["annotations"]
    assert "celltypist_cell_type_v2" in reproduced.uns["openbio_singlecell"]["annotations"]
    assert reproduced.uns["celltypist"] == second.uns["celltypist"]
    assert reproduced_summary == report.summary
    science.pd.testing.assert_frame_equal(first.obs, first_obs)
    assert first.uns == first_uns
    science.np.testing.assert_array_equal(_matrix_values(first.X, science), first_x)

    tampered = first.copy()
    tampered.uns["openbio_singlecell"]["annotations"]["celltypist_cell_type"] = {
        "schema_version": 1,
        "operation": "user_annotation",
    }
    tampered_obs = tampered.obs.copy(deep=True)
    loads_before = len(fake.loads)
    calls_before = len(fake.calls)
    with pytest.raises(ValueError, match="value does not match the verified prior CellTypist provenance"):
        celltypist_owned(
            tampered,
            expression_state="verified_counts",
            model=str(model_path),
            label_column="celltypist_cell_type_v2",
            probability_key="celltypist_probabilities_v2",
            overwrite_existing=True,
        )
    science.pd.testing.assert_frame_equal(tampered.obs, tampered_obs)
    assert len(fake.loads) == loads_before
    assert len(fake.calls) == calls_before


def test_celltypist_overwrite_refuses_unowned_collisions_before_backend(science, tmp_path, monkeypatch):
    adata = _annotation_adata(science)
    adata.obs["celltypist_cell_type"] = "user-owned"
    original_obs = adata.obs.copy(deep=True)
    original_uns = copy.deepcopy(adata.uns)
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)

    with pytest.raises(ValueError, match="ownership of existing outputs could not be verified"):
        celltypist_owned(
            adata.copy(),
            expression_state="verified_counts",
            model=str(model_path),
            overwrite_existing=True,
        )
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns
    assert fake.loads == []
    assert fake.calls == []


def test_celltypist_majority_preflight_and_malformed_backend_are_atomic(science, tmp_path, monkeypatch):
    model_path = tmp_path / "fake.pkl"
    model_path.write_bytes(b"trusted fake model")
    small = _annotation_adata(science, cells=50)
    fake = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", fake)
    with pytest.raises(ValueError, match="more than 50 cells"):
        celltypist_owned(
            small,
            expression_state="verified_counts",
            model=str(model_path),
            majority_voting=True,
            over_clustering_key="over",
        )
    with pytest.raises(ValueError, match="explicit categorical over_clustering_key"):
        celltypist_owned(
            _annotation_adata(science),
            expression_state="verified_counts",
            model=str(model_path),
            majority_voting=True,
        )
    assert fake.calls == []


@pytest.mark.parametrize(
    ("malformed", "message"),
    [
        ("index", "misaligned probability_matrix"),
        ("class_order", "class columns"),
        ("range", "within \\[0, 1\\]"),
        ("labels", "invalid individual label"),
        ("majority", "inconsistent majority labels"),
    ],
)
def test_celltypist_rejects_malformed_backend_without_partial_writes(
    science, tmp_path, monkeypatch, malformed, message
):
    adata = _annotation_adata(science)
    original_obs = adata.obs.copy(deep=True)
    original_uns = copy.deepcopy(adata.uns)
    model_path = tmp_path / f"{malformed}.pkl"
    model_path.write_bytes(b"trusted fake model")

    kwargs = {}
    if malformed == "majority":
        kwargs = {"majority_voting": True, "over_clustering_key": "over", "min_prop": 0.6}
    good = _fake_celltypist(science, model_path)
    monkeypatch.setitem(sys.modules, "celltypist", good)
    _, _, code = celltypist_owned(
        adata.copy(),
        expression_state="verified_counts",
        model=str(model_path),
        **kwargs,
    )
    fake = _fake_celltypist(science, model_path, malformed=malformed)
    monkeypatch.setitem(sys.modules, "celltypist", fake)
    with pytest.raises(RuntimeError, match=message):
        celltypist_owned(
            adata,
            expression_state="verified_counts",
            model=str(model_path),
            **kwargs,
        )
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns
    assert "celltypist_probabilities" not in adata.obsm

    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError, match=message):
        namespace["run_celltypist_annotation"](adata)
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns
    assert "celltypist_probabilities" not in adata.obsm


@pytest.mark.parametrize(
    ("mapping_json", "message"),
    [
        ('{"0":"T","0":"B"}', "duplicate key"),
        ('{"0":"T"," 0 ":"B"}', "collide after trimming"),
        ('{"0":null}', "must be a string"),
        ('{"0":1}', "must be a string"),
        ('{"0":true}', "must be a string"),
        ('{"0":["T"]}', "must be a string"),
        ('{"0":{"label":"T"}}', "must be a string"),
        ('{"0":NaN}', "non-standard JSON constant"),
        ("[]", "non-empty JSON object"),
        ("{}", "non-empty JSON object"),
    ],
)
def test_map_cluster_annotations_strict_json(science, mapping_json, message):
    with pytest.raises((TypeError, ValueError), match=message):
        map_cluster_annotations_owned(
            _mapping_adata(science),
            mapping_json=mapping_json,
            unmapped_policy="set_missing",
        )


def test_map_cluster_annotation_policies_order_categories_provenance_and_code(science):
    adata = _mapping_adata(science)
    original_x = adata.X.copy()
    original_obs = adata.obs.copy(deep=True)
    original_uns = copy.deepcopy(adata.uns)
    mapping_json = '{"0":"T cell","1":"T cell","unused":"Unused target"}'

    with pytest.raises(ValueError, match="observed unmapped levels"):
        map_cluster_annotations_owned(adata, mapping_json=mapping_json)

    output, report, code = map_cluster_annotations_owned(
        adata.copy(),
        mapping_json=mapping_json,
        unmapped_policy="preserve_cluster_label",
        annotation_status="curated",
    )
    assert output.obs["cell_type"].cat.categories.tolist() == ["2", "T cell"]
    assert output.obs["cell_type"].cat.ordered is False
    assert output.obs["cell_type"].astype("string").tolist() == [
        "T cell",
        "T cell",
        "2",
        "T cell",
        science.pd.NA,
    ]
    provenance = output.uns["openbio_singlecell"]["annotations"]["cell_type"]
    assert provenance["annotation_status"] == "curated"
    assert provenance["curation_assertion"] == "caller_declared"
    assert provenance["unused_declared_level_mapping"] == {"unused": "Unused target"}
    assert provenance["many_to_one_merges"] == {"T cell": ["0", "1"]}
    assert provenance["output_categories"] == ["2", "T cell"]
    key_results = report.summary["key_results"]
    assert key_results["missing_source_cells"] == 1
    assert key_results["preserved_cells"] == 1
    assert key_results["unused_mapping_entries"] == 1
    assert "caller-declared Curated annotation" in report.summary["results"]
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<map-code>", "exec")

    namespace = {}
    exec(code, namespace)
    reproduced = namespace["map_cluster_annotations"](adata.copy())
    science.pd.testing.assert_series_equal(reproduced.obs["cell_type"], output.obs["cell_type"])
    assert (
        reproduced.uns["openbio_singlecell"]["annotations"]["cell_type"]
        == output.uns["openbio_singlecell"]["annotations"]["cell_type"]
    )
    science.np.testing.assert_array_equal(adata.X, original_x)
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns

    missing_output, _, _ = map_cluster_annotations_owned(
        adata.copy(),
        mapping_json=mapping_json,
        unmapped_policy="set_missing",
    )
    assert missing_output.obs["cell_type"].cat.categories.tolist() == ["T cell"]
    assert int(missing_output.obs["cell_type"].isna().sum()) == 2


@pytest.mark.parametrize(
    ("unmapped_policy", "mapping_json", "expected_fragment"),
    [
        ("error", '{"0":"A","1":"B","2":"C"}', "preserved 0 cells and set 0 cells missing"),
        (
            "preserve_cluster_label",
            '{"0":"A","1":"B"}',
            "preserved 1 cells and set 0 cells missing",
        ),
        ("set_missing", '{"0":"A","1":"B"}', "preserved 0 cells and set 1 cells missing"),
    ],
)
def test_map_cluster_annotation_results_disclose_each_unmapped_policy(
    science, unmapped_policy, mapping_json, expected_fragment
):
    _, report, _ = map_cluster_annotations_owned(
        _mapping_adata(science),
        mapping_json=mapping_json,
        unmapped_policy=unmapped_policy,
    )

    results = report.summary["results"]
    assert f"Policy {unmapped_policy!r}" in results
    assert expected_fragment in results
    assert "1 cells were already missing in the source" in results
    assert "with counts" in results
    assert "many-to-one target merges" in results
    assert "overwrite_existing was false" in results


def test_map_cluster_annotation_rejects_unknown_colliding_and_ambiguous_labels(science):
    adata = _mapping_adata(science)
    with pytest.raises(ValueError, match="unknown to observed or declared"):
        map_cluster_annotations_owned(
            adata.copy(),
            mapping_json='{"typo":"T"}',
            unmapped_policy="set_missing",
        )
    with pytest.raises(ValueError, match="collide with preserved"):
        map_cluster_annotations_owned(
            adata,
            mapping_json='{"0":"2"}',
            unmapped_policy="preserve_cluster_label",
        )
    with pytest.raises(ValueError, match="cannot overwrite the groupby"):
        map_cluster_annotations_owned(
            adata,
            mapping_json='{"0":"T"}',
            output_column="leiden",
            unmapped_policy="set_missing",
        )

    collision = adata.copy()
    collision.obs["cell_type"] = "prior"
    with pytest.raises(ValueError, match="already exists"):
        map_cluster_annotations_owned(
            collision,
            mapping_json='{"0":"T"}',
            unmapped_policy="set_missing",
        )

    ambiguous = science.ad.AnnData(
        science.np.ones((2, 1)),
        obs=science.pd.DataFrame({"leiden": science.pd.Series([1, "1"], dtype=object)}, index=["a", "b"]),
        var=science.pd.DataFrame(index=["G"]),
    )
    # Build with aligned Series data; the constructor above aligns a default integer index to a/b as missing.
    ambiguous.obs["leiden"] = science.pd.Series([1, "1"], index=ambiguous.obs_names, dtype=object)
    with pytest.raises(ValueError, match="distinct source labels.*collapse"):
        map_cluster_annotations_owned(
            ambiguous,
            mapping_json='{"1":"T"}',
            unmapped_policy="set_missing",
        )


def test_map_cluster_annotation_overwrite_is_explicit_and_replaces_provenance(science):
    adata = _mapping_adata(science)
    adata.obs["cell_type"] = "prior"
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "version": "0.2.0",
        "display_name": "test",
        "source": {},
        "random_seed": 0,
        "warnings": [],
        "analysis_history": {},
        "annotations": {"cell_type": {"operation": "prior"}},
    }
    output, report, _ = map_cluster_annotations_owned(
        adata.copy(),
        mapping_json='{"0":"A","1":"B","2":"C"}',
        overwrite_existing=True,
    )
    provenance = output.uns["openbio_singlecell"]["annotations"]["cell_type"]
    assert provenance["operation"] == "map_cluster_annotations"
    assert provenance["output_existed"] is True
    assert "output column existed before execution: yes" in report.summary["results"]
    assert "overwrite_existing was true" in report.summary["results"]
    assert adata.uns["openbio_singlecell"]["annotations"]["cell_type"] == {"operation": "prior"}


def test_marker_ora_evidence_schema_is_explicit():
    evidence = OpenBioSingleCellMarkerORAEvidence.define_schema()

    assert [item.id for item in evidence.inputs] == [
        "table",
        "universe",
        "resource_csv",
        "resource_metadata_json",
        "source_column",
        "target_column",
        "min_targets",
        "min_overlap",
        "max_p_adjusted",
    ]
    assert [(item.display_name, item.io_type) for item in evidence.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert OpenBioSingleCellMarkerORAEvidence in ANNOTATION_NODE_CLASSES


def test_marker_ora_evidence_gold_statistics_resource_accounting_and_code(science, comfy_directories, monkeypatch):
    input_dir, _, _ = comfy_directories
    resource_path = _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    original_marker = marker.table.copy(deep=True)
    original_universe = universe.table.copy(deep=True)
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    result, report, code = _marker_ora_owned(
        marker,
        universe,
        resource_csv="ora_resource.csv",
        resource_metadata_json=_resource_metadata(),
        min_targets=2,
        min_overlap=2,
        max_p_adjusted=1.0,
    )

    evidence = result.table
    assert list(evidence.columns) == ORA_EVIDENCE_COLUMNS
    assert len(evidence) == 2 * 3
    assert len(fake.calls) == 2
    assert [call["features"] for call in fake.calls] == [["G1", "G2", "G3"], ["G4", "G5", "G6"]]
    assert all(call["n_bg"] == 8 for call in fake.calls)
    assert all(call["alternative"] == "greater" for call in fake.calls)
    assert all(call["ha_corr"] == 0.5 and call["tmin"] == 2 for call in fake.calls)
    assert all(call["verbose"] is False for call in fake.calls)
    assert all(set(call["net"]["source"]) == {"TypeA", "TypeA_tie", "TypeB"} for call in fake.calls)

    row = evidence.loc[(evidence["group"] == "A") & (evidence["source"] == "TypeA")].iloc[0]
    assert (row["a"], row["b"], row["c"], row["d"]) == (2, 2, 1, 3)
    assert row["overlap_genes"] == '["G1","G2"]'
    expected_log_odds = science.np.log((2.5 * 3.5) / (2.5 * 1.5))
    assert row["log_odds_ratio"] == pytest.approx(expected_log_odds)
    from scipy import stats

    expected_p = stats.fisher_exact([[2, 2], [1, 3]], alternative="greater").pvalue
    assert row["p_value"] == pytest.approx(expected_p)
    assert row["candidate_status"] == "tied_best"
    tied = evidence.loc[(evidence["group"] == "A") & evidence["source"].isin(["TypeA", "TypeA_tie"])]
    assert tied["rank_within_group"].tolist() == [1, 1]
    assert tied["candidate_status"].tolist() == ["tied_best", "tied_best"]
    assert bool((evidence["p_adj_global"] >= evidence["p_value"]).all())
    assert result.parameters["n_bg"] == 8
    assert result.parameters["alternative"] == "greater"

    key_results = report.summary["key_results"]
    assert key_results["resource"]["duplicate_rows_removed"] == 1
    assert key_results["resource"]["source_count_before"] == 5
    assert key_results["resource"]["source_count_after"] == 3
    assert key_results["resource"]["sources_removed_by_min_targets_count"] == 2
    assert key_results["resource"]["targets_in_universe_count"] == 7
    assert key_results["resource"]["tested_targets_in_universe_count"] == 6
    assert key_results["unmatched_selected_counts_by_group_preview"] == {"A": 0, "B": 1}
    assert key_results["backend"]["decoupler_version"] == "2.2.0"
    assert key_results["fixed_policy"]["n_bg"] == 8
    assert set(report.summary) == {
        "schema_version",
        "node_id",
        "methods",
        "results",
        "key_results",
        "parameters",
        "warnings",
        "limitations",
        "references",
        "software_versions",
    }
    assert "decoupler 2.x mt.query_set" in report.summary["methods"]
    assert "No annotation label was committed" in report.summary["results"]
    assert report.summary["parameters"] == result.parameters
    assert set(report.summary["software_versions"]) == {
        "python",
        "openbio-singlecell",
        "decoupler",
        "scipy",
        "pandas",
        "numpy",
        "scanpy",
        "anndata",
    }
    assert report.summary["software_versions"]["decoupler"] == "2.2.0"
    assert report.summary["software_versions"]["scipy"] == key_results["backend"]["scipy_version"]
    assert report.summary["warnings"]
    assert report.summary["limitations"]
    resource_reference = report.summary["references"][-1]
    assert resource_reference["citation"] == (
        "Test marker sets 1.0 (2026-08-28): Synthetic resource for deterministic tests"
    )
    assert resource_reference["url"] == f"urn:sha256:{key_results['resource']['sha256']}"
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<marker-ora-code>", "exec")

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_marker_ora_evidence"](marker, universe)
    science.pd.testing.assert_frame_equal(reproduced, evidence)
    assert reproduced_summary == report.summary
    tampered_marker_frame = marker.table.copy(deep=True)
    tampered_marker_frame.loc[0, "score"] += 0.25
    with pytest.raises(ValueError, match="current-content fingerprint"):
        _marker_ora_owned(
            replace(marker, table=tampered_marker_frame),
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
            min_overlap=2,
            max_p_adjusted=1.0,
        )
    with pytest.raises(ValueError, match="marker table content does not match its pinned"):
        namespace["run_marker_ora_evidence"](tampered_marker_frame, universe.table)
    resource_path.write_text(resource_path.read_text(encoding="utf-8") + "TypeB,G6\n", encoding="utf-8")
    with pytest.raises(ValueError, match="resource fingerprint changed"):
        namespace["run_marker_ora_evidence"](marker, universe)
    science.pd.testing.assert_frame_equal(marker.table, original_marker)
    science.pd.testing.assert_frame_equal(universe.table, original_universe)


def test_marker_ora_retains_groups_without_selected_markers_as_unresolved(science, comfy_directories, monkeypatch):
    input_dir, _, _ = comfy_directories
    _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    marker = _replace_filtered_marker_table(
        marker,
        marker.table.loc[marker.table["group"] == "A"].reset_index(drop=True),
    )
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    result, report, code = _marker_ora_owned(
        marker,
        universe,
        resource_csv="ora_resource.csv",
        resource_metadata_json=_resource_metadata(),
        min_targets=2,
        min_overlap=2,
        max_p_adjusted=1.0,
    )

    unresolved = result.table.loc[result.table["group"] == "B"]
    assert len(unresolved) == 3
    assert unresolved["selected_count"].tolist() == [0, 0, 0]
    assert unresolved["a"].tolist() == [0, 0, 0]
    assert unresolved["c"].tolist() == [0, 0, 0]
    assert unresolved["p_value"].tolist() == [1.0, 1.0, 1.0]
    assert not bool(unresolved["positive_enrichment"].any())
    assert set(unresolved["candidate_status"]) == {"unresolved_no_selected_markers"}
    assert [call["features"] for call in fake.calls] == [["G1", "G2", "G3"], []]
    key_results = report.summary["key_results"]
    assert key_results["groups_without_selected_markers_count"] == 1
    assert key_results["groups_without_selected_markers_preview"] == ["B"]
    assert key_results["positive_rows"] == int(
        result.table.loc[result.table["selected_count"] > 0, "positive_enrichment"].sum()
    )
    assert "including 1 with no selected marker" in report.summary["results"]
    assert any("complete ORA rows are retained as unresolved" in warning for warning in report.summary["warnings"])

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_marker_ora_evidence"](marker, universe)
    science.pd.testing.assert_frame_equal(reproduced, result.table)
    assert reproduced_summary == report.summary


@pytest.mark.parametrize(
    ("metadata", "message"),
    [
        ("{}", "exactly the required fields"),
        (_resource_metadata().replace('"name"', '"name": "duplicate", "name"', 1), "duplicate key"),
        (_resource_metadata().replace('"human"', "null"), "must be a nonblank string"),
        (_resource_metadata()[:-1] + ',"unknown":"x"}', "unknown=.*unknown"),
    ],
)
def test_marker_ora_resource_metadata_is_strict_before_backend(
    science, comfy_directories, monkeypatch, metadata, message
):
    input_dir, _, _ = comfy_directories
    _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    with pytest.raises((TypeError, ValueError), match=message):
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=metadata,
            min_targets=2,
        )
    assert fake.calls == []


def test_marker_ora_resource_fingerprint_includes_exact_bytes(comfy_directories):
    input_dir, _, _ = comfy_directories
    resource_path = _write_ora_resource(input_dir)
    original_stat = resource_path.stat()
    first = OpenBioSingleCellMarkerORAEvidence.fingerprint_inputs("ora_resource.csv")

    changed = resource_path.read_bytes().replace(b"TypeA,G1", b"TypeZ,G1", 1)
    assert len(changed) == original_stat.st_size
    resource_path.write_bytes(changed)
    os.utime(resource_path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    second = OpenBioSingleCellMarkerORAEvidence.fingerprint_inputs("ora_resource.csv")

    assert first[1:-1] == second[1:-1]
    assert first[-1] != second[-1]


@pytest.mark.parametrize(
    ("malformed_csv", "message"),
    [
        ("source,target\njunk,A,G1\n", r"data row 1 has 3 fields; expected exactly 2"),
        ("source,target\nA,G1,junk\n", r"data row 1 has 3 fields; expected exactly 2"),
        ("source,target\nA\n", r"data row 1 has 1 fields; expected exactly 2"),
    ],
)
def test_marker_ora_rejects_ragged_resource_rows_in_runtime_and_generated_code(
    science, comfy_directories, monkeypatch, malformed_csv, message
):
    input_dir, _, _ = comfy_directories
    resource_path = _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)
    _, report, code = _marker_ora_owned(
        marker,
        universe,
        resource_csv="ora_resource.csv",
        resource_metadata_json=_resource_metadata(),
        min_targets=2,
    )
    old_sha256 = report.summary["key_results"]["resource"]["sha256"]

    resource_path.write_text(malformed_csv, encoding="utf-8")
    new_sha256 = hashlib.sha256(resource_path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match=message) as runtime_error:
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )

    namespace = {}
    exec(code.replace(old_sha256, new_sha256), namespace)
    with pytest.raises(ValueError, match=message) as generated_error:
        namespace["run_marker_ora_evidence"](marker, universe)
    assert str(generated_error.value) == str(runtime_error.value)
    assert fake.calls


def test_marker_ora_accepts_quoted_commas_without_treating_the_row_as_ragged(science, comfy_directories, monkeypatch):
    input_dir, _, _ = comfy_directories
    resource_path = input_dir / "quoted_ora_resource.csv"
    resource_path.write_text(
        'source,target\n"Type,A",G1\n"Type,A",G2\n"Type,A",G3\n',
        encoding="utf-8",
    )
    marker, universe = _marker_artifacts(science)
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    result, _, _ = _marker_ora_owned(
        marker,
        universe,
        resource_csv="quoted_ora_resource.csv",
        resource_metadata_json=_resource_metadata(),
        min_targets=1,
        min_overlap=1,
    )

    assert result.table["source"].unique().tolist() == ["Type,A"]


def test_marker_ora_rejects_upstream_and_resource_contract_failures_before_backend(
    science, comfy_directories, monkeypatch
):
    input_dir, _, _ = comfy_directories
    resource_path = _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    changed_marker = marker.table.copy()
    changed_marker.loc[0, "score"] += 0.25
    with pytest.raises(ValueError, match="current-content fingerprint"):
        _marker_ora_owned(
            replace(marker, table=changed_marker),
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    changed_universe = universe.table.copy()
    changed_universe.loc[0, "gene"] = "DIFFERENT"
    with pytest.raises(ValueError, match="universe content does not match"):
        _marker_ora_owned(
            marker,
            replace(universe, table=changed_universe),
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    foreign_universe_parameters = {**universe.parameters, "ranking_fingerprint": "d" * 64}
    foreign_universe = make_table_result(
        table=universe.table.copy(deep=True),
        title="foreign universe",
        operation="marker_genes",
        parameters=foreign_universe_parameters,
        description="cross-pair fixture",
        warnings=[],
        input_cells=universe.input_cells,
        input_genes=universe.input_genes,
        started_at=time.perf_counter(),
    )
    with pytest.raises(ValueError, match="ranking fingerprints do not match"):
        _marker_ora_owned(
            marker,
            foreign_universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    truncated_marker, truncated_universe = _marker_artifacts(science, ranking_truncated=True)
    with pytest.raises(ValueError, match="non-truncated upstream ranking"):
        _marker_ora_owned(
            truncated_marker,
            truncated_universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    unfiltered_source = dict(marker.source)
    unfiltered_source["operation"] = "marker_genes"
    with pytest.raises(ValueError, match="requires table operation.*filter_marker_genes"):
        _marker_ora_owned(
            replace(marker, source=unfiltered_source),
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    resource_path.write_text("source,target\nWrong,NOT_IN_UNIVERSE\n", encoding="utf-8")
    with pytest.raises(ValueError, match="zero target overlap"):
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=1,
        )
    assert fake.calls == []


@pytest.mark.parametrize(
    ("malformed", "message"),
    [
        ("missing_source", "source set.*differs"),
        ("stat", "log odds disagree"),
        ("pval", "p-values disagree"),
        ("padj", "adjusted p-values disagree"),
        ("nonfinite", "non-finite"),
        ("dtype", "must be numeric and non-boolean"),
    ],
)
def test_marker_ora_rejects_malformed_decoupler_results(science, comfy_directories, monkeypatch, malformed, message):
    input_dir, _, _ = comfy_directories
    _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)
    original_marker = marker.table.copy(deep=True)
    original_universe = universe.table.copy(deep=True)
    good_fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", good_fake)
    _, _, code = _marker_ora_owned(
        marker,
        universe,
        resource_csv="ora_resource.csv",
        resource_metadata_json=_resource_metadata(),
        min_targets=2,
    )
    fake = _fake_decoupler(science, malformed=malformed)
    monkeypatch.setitem(sys.modules, "decoupler", fake)

    with pytest.raises(RuntimeError, match=message) as runtime_error:
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError, match=message) as generated_error:
        namespace["run_marker_ora_evidence"](marker, universe)
    assert str(generated_error.value) == str(runtime_error.value)
    science.pd.testing.assert_frame_equal(marker.table, original_marker)
    science.pd.testing.assert_frame_equal(universe.table, original_universe)


def test_marker_ora_requires_decoupler_2_query_set_after_preflight(science, comfy_directories, monkeypatch):
    input_dir, _, _ = comfy_directories
    _write_ora_resource(input_dir)
    marker, universe = _marker_artifacts(science)

    legacy = _fake_decoupler(science, version="1.9.2")
    monkeypatch.setitem(sys.modules, "decoupler", legacy)
    with pytest.raises(RuntimeError, match="requires decoupler 2.x"):
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )

    missing_api = types.ModuleType("decoupler")
    missing_api.__version__ = "2.2.0"
    missing_api.mt = types.SimpleNamespace()
    monkeypatch.setitem(sys.modules, "decoupler", missing_api)
    with pytest.raises(RuntimeError, match="mt.query_set API"):
        _marker_ora_owned(
            marker,
            universe,
            resource_csv="ora_resource.csv",
            resource_metadata_json=_resource_metadata(),
            min_targets=2,
        )
