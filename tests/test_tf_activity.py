from __future__ import annotations

import copy
import importlib
import json
from types import SimpleNamespace

import pytest
from scipy import stats

from openbio_singlecell.collectri_ulm import (
    COLLECTRI_PADJ_KEY,
    COLLECTRI_SCORE_KEY,
    COLLECTRI_UNS_KEY,
    collectri_ulm_code,
    run_collectri_ulm,
)
from openbio_singlecell.tf_activity_artifact import (
    TF_ACTIVITY_ARTIFACT_TYPE,
    TFActivityArtifact,
    build_tf_activity_artifact,
    validate_tf_activity_artifact,
)
from openbio_singlecell.tf_activity_ranking import (
    TF_RANKING_COLUMNS,
    rank_tf_activities,
    rank_tf_activities_code,
)


def _artifact_provenance() -> dict[str, object]:
    return {
        "method": "CollecTRI ULM",
        "expression": {
            "source": "X",
            "feature_axis_sha256": "f" * 64,
        },
        "resource": {
            "name": "toy-collectri",
            "organism": "human",
            "identifier_namespace": "HGNC symbol",
            "canonical_sha256": "a" * 64,
        },
        "software_versions": {"decoupler": "2.2.0"},
    }


def _activity_frames(science):
    observations = [f"cell_{index}" for index in range(6)]
    regulators = ["TF_A", "TF_B", "TF_C"]
    scores = science.pd.DataFrame(
        [
            [3.0, -1.0, 0.0],
            [2.0, -2.0, 0.5],
            [1.0, -3.0, 1.0],
            [-1.0, 3.0, 0.0],
            [-2.0, 2.0, -0.5],
            [-3.0, 1.0, -1.0],
        ],
        index=observations,
        columns=regulators,
    )
    adjusted = science.pd.DataFrame(
        [
            [0.01, 0.30, 1.0],
            [0.02, 0.20, 0.8],
            [0.05, 0.10, 0.7],
            [0.05, 0.10, 1.0],
            [0.02, 0.20, 0.8],
            [0.01, 0.30, 0.7],
        ],
        index=observations,
        columns=regulators,
    )
    return scores, adjusted


def _activity_artifact(science):
    scores, adjusted = _activity_frames(science)
    return build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance=_artifact_provenance(),
        numpy=science.np,
        pandas=science.pd,
    )


def test_tf_activity_file_codec_round_trip_uses_h5ad_and_json(science, tmp_path):
    from openbio_singlecell.regulatory_artifact_codecs import (
        read_tf_activity,
        write_tf_activity,
    )

    artifact = _activity_artifact(science)
    root = tmp_path / "tf-activity"
    root.mkdir()

    write_tf_activity(root, artifact)
    portable = read_tf_activity(root)
    scores, adjusted, provenance, metadata = validate_tf_activity_artifact(
        portable,
        exact_type=False,
        numpy=science.np,
        pandas=science.pd,
        copy_result=False,
    )

    expected_scores, expected_adjusted = _activity_frames(science)
    science.pd.testing.assert_frame_equal(scores, expected_scores)
    science.pd.testing.assert_frame_equal(adjusted, expected_adjusted)
    assert provenance == artifact.provenance
    assert metadata == artifact.metadata
    assert (root / "data.h5ad").is_file()
    assert (root / "result.json").is_file()


def _resource_metadata_json() -> str:
    return json.dumps(
        {
            "name": "toy-signed-collectri",
            "version": "2026.08",
            "date": "2026-08-28",
            "organism": "human",
            "identifier_namespace": "HGNC symbol",
            "scope": "deterministic signed unit-test network",
            "license": "CC-BY-4.0",
            "citation": "OpenBio test fixture, 2026.",
        }
    )


def _write_signed_network(path, *, conflicting_duplicate: bool = False) -> None:
    rows = [
        "source,target,weight",
        "TF_A,G0,1",
        "TF_A,G1,1",
        "TF_A,G2,-1",
        "TF_B,G2,1",
        "TF_B,G3,-1",
        "TF_B,G4,1",
        "AP1,G0,1",
        "AP1,G4,-1",
        "AP1,G5,1",
        "NFKB,G1,-1",
        "NFKB,G5,1",
        "NFKB,G6,1",
        "TF_A,G0,1",
    ]
    if conflicting_duplicate:
        rows.append("TF_A,G0,-1")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _logged_adata(science, *, n_obs: int = 12):
    rng = science.np.random.default_rng(20260828)
    counts = rng.poisson(5.0, size=(n_obs, 8)).astype(float)
    counts[:, 0] += science.np.arange(n_obs) % 4
    counts[:, 2] += science.np.arange(n_obs)[::-1] % 5
    normalized = counts / counts.sum(axis=1, keepdims=True) * 10_000.0
    logged = science.np.log1p(normalized)
    obs = science.pd.DataFrame(index=[f"cell_{index:02d}" for index in range(n_obs)])
    var = science.pd.DataFrame(index=[f"G{index}" for index in range(8)])
    adata = science.ad.AnnData(science.sparse.csr_matrix(logged), obs=obs, var=var)
    adata.uns["log1p"] = {"base": None}
    return adata


def _ranking_input(science):
    groups = ["A"] * 4 + ["B"] * 4 + ["C"] * 4 + ["REF"] * 4
    observations = [f"rank_cell_{index:02d}" for index in range(len(groups))]
    obs = science.pd.DataFrame(
        {
            "cell_type": science.pd.Categorical(groups, categories=["unused", "C", "A", "B", "REF"]),
        },
        index=observations,
    )
    adata = science.ad.AnnData(
        science.np.zeros((len(groups), 2), dtype=float),
        obs=obs,
        var=science.pd.DataFrame(index=["G0", "G1"]),
    )
    scores = science.pd.DataFrame(
        {
            "TF_UP": [8.0, 9.0, 7.0, 10.0, 1.0, 2.0, 0.0, 3.0, 4.0, 5.0, 3.0, 6.0, -2.0, -1.0, -3.0, 0.0],
            "TF_DOWN": [-8.0, -9.0, -7.0, -10.0, -1.0, -2.0, 0.0, -3.0, -4.0, -5.0, -3.0, -6.0, 2.0, 1.0, 3.0, 0.0],
            "TF_TIE_A": [1.0, 3.0, 2.0, 4.0, 0.0, 1.0, 2.0, 3.0, 2.0, 4.0, 3.0, 5.0, -1.0, 0.0, 1.0, 2.0],
            "TF_TIE_B": [1.0, 3.0, 2.0, 4.0, 0.0, 1.0, 2.0, 3.0, 2.0, 4.0, 3.0, 5.0, -1.0, 0.0, 1.0, 2.0],
        },
        index=observations,
    )
    adjusted = science.pd.DataFrame(0.5, index=observations, columns=scores.columns)
    artifact = build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance=_artifact_provenance(),
        numpy=science.np,
        pandas=science.pd,
    )
    return adata, artifact


def _adata_snapshot(adata):
    return {
        "X": adata.X.copy(),
        "obs": adata.obs.copy(deep=True),
        "var": adata.var.copy(deep=True),
        "uns": copy.deepcopy(adata.uns),
        "obsm": {key: value.copy() for key, value in adata.obsm.items()},
        "layers": {key: value.copy() for key, value in adata.layers.items()},
    }


def _assert_adata_unchanged(adata, snapshot, science) -> None:
    if science.sparse.issparse(adata.X):
        assert (adata.X != snapshot["X"]).nnz == 0
    else:
        science.np.testing.assert_array_equal(adata.X, snapshot["X"])
    science.pd.testing.assert_frame_equal(adata.obs, snapshot["obs"])
    science.pd.testing.assert_frame_equal(adata.var, snapshot["var"])
    assert adata.uns == snapshot["uns"]
    assert set(adata.obsm) == set(snapshot["obsm"])
    assert set(adata.layers) == set(snapshot["layers"])


def _bh_adjust(pvalues, np):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _fake_decoupler(science, *, network=None, mutation: str | None = None):
    calls: list[dict[str, object]] = []

    class FakeULM:
        def __call__(
            self,
            data,
            net,
            tmin=5,
            raw=False,
            empty=True,
            bsize=250_000,
            verbose=False,
            **kwargs,
        ):
            calls.append(
                {
                    "operation": "ulm",
                    "tmin": tmin,
                    "raw": raw,
                    "empty": empty,
                    "bsize": bsize,
                    "verbose": verbose,
                    "kwargs": dict(kwargs),
                    "net": net.copy(deep=True),
                }
            )
            features = data.var_names.tolist()
            sources = list(dict.fromkeys(net["source"].tolist()))
            matrix = data.X.toarray() if science.sparse.issparse(data.X) else science.np.asarray(data.X)
            score_values = science.np.empty((data.n_obs, len(sources)), dtype=float)
            raw_pvalues = science.np.empty_like(score_values)
            degrees_of_freedom = len(features) - 2
            for source_index, source in enumerate(sources):
                weights_by_target = dict(
                    net.loc[net["source"] == source, ["target", "weight"]].itertuples(index=False, name=None)
                )
                weights = science.np.asarray([weights_by_target.get(feature, 0.0) for feature in features], dtype=float)
                for observation_index, values in enumerate(matrix):
                    correlation = float(science.np.corrcoef(weights, values)[0, 1])
                    statistic = correlation * science.np.sqrt(
                        degrees_of_freedom / max(1.0 - correlation * correlation, science.np.finfo(float).tiny)
                    )
                    score_values[observation_index, source_index] = statistic
                    raw_pvalues[observation_index, source_index] = 2.0 * stats.t.sf(
                        abs(statistic), degrees_of_freedom
                    )
            adjusted_values = science.np.vstack(
                [_bh_adjust(row, science.np) for row in raw_pvalues]
            )
            data.obsm["score_ulm"] = science.pd.DataFrame(
                score_values,
                index=data.obs_names.copy(),
                columns=sources,
            )
            data.obsm["padj_ulm"] = science.pd.DataFrame(
                adjusted_values,
                index=data.obs_names.copy(),
                columns=sources,
            )
            if mutation == "ulm_obs":
                data.obs["unexpected_backend_mutation"] = "changed"
            elif mutation == "ulm_extra_output":
                data.obsm["unexpected_backend_output"] = data.obsm["score_ulm"].copy(deep=True)
            elif mutation == "ulm_nonfinite":
                data.obsm["score_ulm"].iloc[0, 0] = float("nan")
            elif mutation == "ulm_raw":
                data.raw = data.copy()
            elif mutation == "ulm_obsp":
                data.obsp["unexpected_backend_mutation"] = science.sparse.eye(data.n_obs, format="csr")
            elif mutation == "ulm_network":
                net.iloc[0, net.columns.get_loc("weight")] *= -1.0
            return None

        def func(self, mat, adj, tval=True, verbose=False):
            raise AssertionError("The adapter must call the reviewed public Method, not its implementation function.")

    def collectri(organism="human", remove_complexes=False, license="academic", verbose=False):
        calls.append(
            {
                "operation": "collectri",
                "organism": organism,
                "remove_complexes": remove_complexes,
                "license": license,
                "verbose": verbose,
            }
        )
        if network is None:
            raise AssertionError("This fake official loader requires a supplied network fixture.")
        result = network.copy(deep=True)
        if remove_complexes:
            result = result.loc[~result["source"].isin(["AP1", "NFKB"])].reset_index(drop=True)
        return result

    def get_obsm(adata, key):
        calls.append({"operation": "get_obsm", "key": key})
        if mutation == "get_obsm_obs":
            adata.obs["unexpected_backend_mutation"] = "changed"
        elif mutation == "get_obsm_var":
            adata.var["unexpected_backend_mutation"] = "changed"
        frame = adata.obsm[key]
        return science.ad.AnnData(
            frame.to_numpy(dtype=float, copy=True),
            obs=adata.obs.copy(deep=True),
            var=science.pd.DataFrame(index=frame.columns.copy()),
        )

    def rankby_group(adata, groupby, reference="rest", method="t-test_overestim_var"):
        calls.append(
            {
                "operation": "rankby_group",
                "groupby": groupby,
                "reference": copy.deepcopy(reference),
                "method": method,
            }
        )
        labels = adata.obs[groupby]
        groups = labels.unique().tolist()
        rows: list[list[object]] = []
        for group in groups:
            group_mask = (labels == group).to_numpy()
            if reference == "rest":
                reference_mask = ~group_mask
                reference_display = "rest"
            elif isinstance(reference, str):
                if group == reference:
                    continue
                reference_mask = (labels == reference).to_numpy()
                reference_display = reference
            else:
                reference_mask = labels.isin(reference).to_numpy()
                reference_display = ", ".join(reference)
            group_values = science.np.asarray(adata.X[group_mask], dtype=float)
            reference_values = science.np.asarray(adata.X[reference_mask], dtype=float)
            group_rows: list[list[object]] = []
            for feature_index, feature in enumerate(adata.var_names):
                selected = group_values[:, feature_index]
                baseline = reference_values[:, feature_index]
                if method == "wilcoxon":
                    statistic, pvalue = stats.ranksums(selected, baseline)
                else:
                    reference_nobs = selected.size if method == "t-test_overestim_var" else baseline.size
                    statistic, pvalue = stats.ttest_ind_from_stats(
                        mean1=selected.mean(),
                        std1=selected.std(ddof=1),
                        nobs1=selected.size,
                        mean2=baseline.mean(),
                        std2=baseline.std(ddof=1),
                        nobs2=reference_nobs,
                        equal_var=False,
                    )
                group_rows.append(
                    [
                        group,
                        reference_display,
                        feature,
                        float(statistic),
                        float(selected.mean() - baseline.mean()),
                        float(pvalue),
                    ]
                )
            adjusted = _bh_adjust([row[-1] for row in group_rows], science.np)
            for row, adjusted_pvalue in zip(group_rows, adjusted, strict=True):
                rows.append([*row, float(adjusted_pvalue)])
        result = science.pd.DataFrame(
            rows,
            columns=["group", "reference", "name", "stat", "meanchange", "pval", "padj"],
        )
        if mutation == "rank_drop_row":
            result = result.iloc[1:].reset_index(drop=True)
        elif mutation == "rank_bad_padj":
            result.loc[0, "padj"] = 0.0
        elif mutation == "rank_nonfinite":
            result.loc[0, "stat"] = float("nan")
        elif mutation == "rank_carrier_obs":
            adata.obs["unexpected_backend_mutation"] = "changed"
        elif mutation == "rank_var":
            adata.var["unexpected_backend_mutation"] = "changed"
        elif mutation == "rank_uns":
            adata.uns["unexpected_backend_mutation"] = "changed"
        return result

    module = SimpleNamespace(
        __version__="2.2.0",
        mt=SimpleNamespace(ulm=FakeULM()),
        op=SimpleNamespace(collectri=collectri),
        pp=SimpleNamespace(get_obsm=get_obsm),
        tl=SimpleNamespace(rankby_group=rankby_group),
        calls=calls,
    )
    return module


def _run_local_ulm(adata, network_path, science, *, decoupler_module=None, **overrides):
    parameters = {
        "resource_mode": "local_network",
        "network_path": str(network_path),
        "resource_metadata_json": _resource_metadata_json(),
        "source_kind": "X",
        "layer_name": None,
        "affiliation_license": "academic",
        "allow_network_access": False,
        "complex_policy": "retain",
        "min_targets": 3,
        "batch_size": 17,
        "max_output_rows": 10_000,
        "max_working_memory_gib": 1.0,
        "overwrite_existing": False,
        "openbio_version": "test-version",
    }
    parameters.update(overrides)
    return run_collectri_ulm(
        adata,
        decoupler_module=_fake_decoupler(science) if decoupler_module is None else decoupler_module,
        **parameters,
    )


def test_worker_owned_collectri_updates_one_private_anndata_without_full_copy(science, tmp_path):
    from openbio_singlecell.operations_regulatory import collectri_ulm_owned

    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)

    output, artifact, report, code = collectri_ulm_owned(
        adata,
        resource_mode="local_network",
        network_path=network_path,
        resource_metadata_json=_resource_metadata_json(),
        source={"source": "X"},
        affiliation_license="academic",
        allow_network_access=False,
        complex_policy="retain",
        min_targets=3,
        batch_size=17,
        max_output_rows=10_000,
        max_working_memory_gib=1.0,
        overwrite_existing=False,
        decoupler_module=_fake_decoupler(science),
    )

    assert output is adata
    validate_tf_activity_artifact(artifact, copy_result=False)
    assert report.summary["node_id"] == "OpenBioSingleCellCollecTRIULM"
    assert "decoupler.mt.ulm" in code


def test_collectri_worker_operation_writes_new_files_without_rewriting_input(
    science,
    tmp_path,
    monkeypatch,
):
    import uuid

    import openbio_singlecell.operations_regulatory as operations
    from openbio_singlecell.artifact_codecs import read_anndata, write_anndata
    from openbio_singlecell.regulatory_artifact_codecs import read_tf_activity
    from openbio_singlecell.worker_protocol import OperationContext

    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    owned = _logged_adata(science)
    result = operations.collectri_ulm_owned(
        owned,
        resource_mode="local_network",
        network_path=network_path,
        resource_metadata_json=_resource_metadata_json(),
        source={"source": "X"},
        affiliation_license="academic",
        allow_network_access=False,
        complex_policy="retain",
        min_targets=3,
        batch_size=17,
        max_output_rows=10_000,
        max_working_memory_gib=1.0,
        overwrite_existing=False,
        decoupler_module=_fake_decoupler(science),
    )
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _logged_adata(science))
    input_bytes = (input_root / "data.h5ad").read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))
    monkeypatch.setattr(operations, "collectri_ulm_owned", lambda *_args, **_kwargs: result)

    records = operations.collectri_ulm(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
            "network_csv": {
                "type": "file",
                "path": str(network_path.resolve()),
                "provenance": {},
            },
        },
        {
            "resource_mode": "local_network",
            "resource_metadata_json": _resource_metadata_json(),
            "source": {"source": "X"},
            "affiliation_license": "academic",
            "allow_network_access": False,
            "complex_policy": "retain",
            "min_targets": 3,
            "batch_size": 17,
            "max_output_rows": 10_000,
            "max_working_memory_gib": 1.0,
            "overwrite_existing": False,
        },
    )

    assert [record["name"] for record in records] == ["adata", "activities", "summary", "code"]
    assert (input_root / "data.h5ad").read_bytes() == input_bytes
    written = read_anndata(staging / records[0]["payload"])
    portable = read_tf_activity(staging / records[1]["payload"])
    assert COLLECTRI_SCORE_KEY in written.obsm
    validate_tf_activity_artifact(portable, exact_type=False, copy_result=False)


def test_collectri_ulm_local_contract_exact_call_full_bh_and_generated_parity(science, tmp_path):
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)
    before = _adata_snapshot(adata)
    generated_input = _logged_adata(science)
    fake = _fake_decoupler(science)

    output, artifact, summary = _run_local_ulm(
        adata,
        network_path,
        science,
        decoupler_module=fake,
    )

    assert [call["operation"] for call in fake.calls] == ["ulm"]
    call = fake.calls[0]
    assert call["tmin"] == 3
    assert call["raw"] is False
    assert call["empty"] is False
    assert call["bsize"] == 17
    assert call["verbose"] is False
    assert call["kwargs"] == {"tval": True}
    assert call["net"].columns.tolist() == ["source", "target", "weight"]
    assert not call["net"].duplicated(["source", "target"]).any()
    assert call["net"]["source"].tolist() == sorted(call["net"]["source"].tolist())

    assert type(artifact) is TFActivityArtifact
    science.pd.testing.assert_frame_equal(output.obsm[COLLECTRI_SCORE_KEY], artifact.scores)
    science.pd.testing.assert_frame_equal(output.obsm[COLLECTRI_PADJ_KEY], artifact.adjusted_pvalues)
    assert output.obsm[COLLECTRI_SCORE_KEY].index.equals(adata.obs_names)
    assert output.obsm[COLLECTRI_SCORE_KEY].columns.tolist() == ["AP1", "NFKB", "TF_A", "TF_B"]
    assert science.np.isfinite(artifact.scores.to_numpy(dtype=float)).all()
    assert science.np.isfinite(artifact.adjusted_pvalues.to_numpy(dtype=float)).all()
    assert ((artifact.adjusted_pvalues >= 0.0) & (artifact.adjusted_pvalues <= 1.0)).all().all()
    degrees_freedom = adata.n_vars - 2
    independent_raw = 2.0 * stats.t.sf(science.np.abs(artifact.scores.to_numpy()), degrees_freedom)
    independent_adjusted = science.np.vstack(
        [_bh_adjust(row, science.np) for row in independent_raw]
    )
    science.np.testing.assert_allclose(
        artifact.adjusted_pvalues.to_numpy(),
        independent_adjusted,
        rtol=5e-5,
        atol=5e-6,
    )

    resource = summary["key_results"]["resource"]
    assert resource["network_access"] is False
    assert len(resource["file_sha256"]) == 64
    assert len(resource["canonical_network_sha256"]) == 64
    assert resource["accounting"]["duplicate_rows_removed"] == 1
    expression = summary["key_results"]["expression"]
    assert expression["source"] == "X"
    assert expression["feature_axis_sha256"]
    assert "state" not in expression
    assert "state_evidence" not in expression
    assert summary["key_results"]["degrees_freedom"] == 6
    assert summary["parameters"]["raw"] is False
    assert summary["parameters"]["empty"] is False
    assert summary["parameters"]["tval"] is True
    assert summary["status"] == "exploratory_cell_level_tf_activity"
    assert len(summary["references"]) >= 3
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    assert output.uns[COLLECTRI_UNS_KEY]["artifact_fingerprint_sha256"] == artifact.fingerprint
    assert output is adata
    if science.sparse.issparse(adata.X):
        assert (adata.X != before["X"]).nnz == 0
    else:
        science.np.testing.assert_array_equal(adata.X, before["X"])
    science.pd.testing.assert_frame_equal(adata.obs, before["obs"])
    science.pd.testing.assert_frame_equal(adata.var, before["var"])
    assert adata.uns["log1p"] == before["uns"]["log1p"]

    parameters = {
        **summary["parameters"],
        "network_path": str(network_path),
        "resource_metadata_json": _resource_metadata_json(),
        "expected_resource_sha256": resource["file_sha256"],
        "expected_canonical_network_sha256": resource["canonical_network_sha256"],
        "openbio_version": "test-version",
    }
    parameters["source_kind"] = parameters.pop("source")
    code = collectri_ulm_code(parameters=parameters)
    compile(code, "<collectri-ulm-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_output, generated_artifact, generated_summary = namespace["infer_collectri_ulm"](
        generated_input,
        decoupler_module=_fake_decoupler(science),
    )
    science.pd.testing.assert_frame_equal(
        generated_output.obsm[COLLECTRI_SCORE_KEY], output.obsm[COLLECTRI_SCORE_KEY]
    )
    science.pd.testing.assert_frame_equal(generated_artifact.scores, artifact.scores)
    science.pd.testing.assert_frame_equal(
        generated_artifact.adjusted_pvalues, artifact.adjusted_pvalues
    )
    assert generated_artifact.fingerprint == artifact.fingerprint
    assert generated_summary == summary


def test_collectri_ulm_complex_policy_official_access_and_license_are_explicit(science, tmp_path):
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)
    local_fake = _fake_decoupler(science)
    output, artifact, summary = _run_local_ulm(
        adata,
        network_path,
        science,
        decoupler_module=local_fake,
        complex_policy="remove",
    )
    assert output.obsm[COLLECTRI_SCORE_KEY].columns.tolist() == ["TF_A", "TF_B"]
    assert artifact.provenance["retained_regulator_count"] == 2
    assert summary["key_results"]["resource"]["accounting"]["complex_edges_removed"] == 6

    official = science.pd.read_csv(network_path).drop_duplicates(["source", "target"])
    official["resources"] = "fixture"
    official["references"] = "fixture:1"
    official["sign_decision"] = "consensus"
    offline_fake = _fake_decoupler(science, network=official)
    with pytest.raises(RuntimeError, match="direct network access"):
        run_collectri_ulm(
            _logged_adata(science),
            resource_mode="official_collectri",
            affiliation_license="commercial",
            allow_network_access=False,
            min_targets=3,
            decoupler_module=offline_fake,
        )
    assert offline_fake.calls == []

    online_fake = _fake_decoupler(science, network=official)
    _, _, official_summary = run_collectri_ulm(
        _logged_adata(science),
        resource_mode="official_collectri",
        affiliation_license="commercial",
        allow_network_access=True,
        complex_policy="remove",
        min_targets=3,
        decoupler_module=online_fake,
        openbio_version="test-version",
    )
    assert online_fake.calls[0] == {
        "operation": "collectri",
        "organism": "human",
        "remove_complexes": True,
        "license": "commercial",
        "verbose": False,
    }
    assert online_fake.calls[1]["operation"] == "ulm"
    official_resource = official_summary["key_results"]["resource"]
    assert official_resource["network_access"] is True
    assert official_resource["cache_state"] == "not_supported_by_decoupler_2_2"
    assert official_resource["affiliation_license_declaration"] == "commercial"
    assert official_resource["affiliation_filter_enforced_by_decoupler_2_2"] is False
    assert official_summary["software_versions"]["requests"] != "not-installed"


def test_official_collectri_generated_code_uses_materialized_network_without_download_and_has_exact_parity(
    science,
    tmp_path,
):
    source_path = tmp_path / "official_source.csv"
    _write_signed_network(source_path)
    official = science.pd.read_csv(source_path).drop_duplicates(["source", "target"])
    official["resources"] = "fixture"
    official["references"] = "fixture:1"
    official["sign_decision"] = "consensus"
    adata = _logged_adata(science)
    online_fake = _fake_decoupler(science, network=official)
    runtime_output, runtime_artifact, runtime_summary = run_collectri_ulm(
        adata,
        resource_mode="official_collectri",
        affiliation_license="nonprofit",
        allow_network_access=True,
        complex_policy="remove",
        min_targets=3,
        batch_size=19,
        max_output_rows=10_000,
        max_working_memory_gib=1.0,
        decoupler_module=online_fake,
        openbio_version="test-version",
    )
    assert [call["operation"] for call in online_fake.calls] == ["collectri", "ulm"]

    materialized_path = tmp_path / "materialized_post_policy_collectri.csv"
    materialized = (
        official.loc[~official["source"].isin(["AP1", "NFKB"]), ["source", "target", "weight"]]
        .sort_values(["source", "target"], kind="mergesort", ignore_index=True)
    )
    materialized.to_csv(materialized_path, index=False)
    resource = runtime_summary["key_results"]["resource"]
    parameters = {
        **runtime_summary["parameters"],
        "network_path": None,
        "resource_metadata_json": "{}",
        "expected_resource_sha256": None,
        "expected_canonical_network_sha256": resource["canonical_network_sha256"],
        "openbio_version": "test-version",
    }
    parameters["source_kind"] = parameters.pop("source")
    code = collectri_ulm_code(
        parameters=parameters,
        resolved_resource_provenance=resource,
    )
    compile(code, "<official-collectri-ulm-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_fake = _fake_decoupler(science)
    generated_output, generated_artifact, generated_summary = namespace["infer_collectri_ulm"](
        _logged_adata(science),
        str(materialized_path),
        decoupler_module=generated_fake,
    )

    assert [call["operation"] for call in generated_fake.calls] == ["ulm"]
    science.pd.testing.assert_frame_equal(
        generated_output.obsm[COLLECTRI_SCORE_KEY], runtime_output.obsm[COLLECTRI_SCORE_KEY]
    )
    science.pd.testing.assert_frame_equal(
        generated_output.obsm[COLLECTRI_PADJ_KEY], runtime_output.obsm[COLLECTRI_PADJ_KEY]
    )
    science.pd.testing.assert_frame_equal(generated_artifact.scores, runtime_artifact.scores)
    science.pd.testing.assert_frame_equal(
        generated_artifact.adjusted_pvalues,
        runtime_artifact.adjusted_pvalues,
    )
    assert generated_artifact.fingerprint == runtime_artifact.fingerprint
    assert generated_summary == runtime_summary

    changed_path = tmp_path / "changed_materialized_collectri.csv"
    changed = materialized.copy(deep=True)
    changed.loc[0, "weight"] *= -1.0
    changed.to_csv(changed_path, index=False)
    changed_fake = _fake_decoupler(science)
    with pytest.raises(ValueError, match="does not match the runtime-resolved canonical SHA-256"):
        namespace["infer_collectri_ulm"](
            _logged_adata(science),
            str(changed_path),
            decoupler_module=changed_fake,
        )
    assert changed_fake.calls == []


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("ulm_obs", "unsupported private AnnData state"),
        ("ulm_extra_output", "output keys changed"),
        ("ulm_nonfinite", "non-finite"),
        ("ulm_raw", "unsupported private AnnData state"),
        ("ulm_obsp", "unsupported private AnnData state"),
        ("ulm_network", "caller-owned network table"),
    ],
)
def test_collectri_ulm_rejects_backend_tampering(science, tmp_path, mutation, message):
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)
    with pytest.raises(RuntimeError, match=message):
        _run_local_ulm(
            adata,
            network_path,
            science,
            decoupler_module=_fake_decoupler(science, mutation=mutation),
        )


def test_collectri_ulm_rejects_version_signature_resource_and_size_guards(science, tmp_path):
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)
    wrong_version = _fake_decoupler(science)
    wrong_version.__version__ = "2.1.4"
    with pytest.raises(RuntimeError, match="exactly 2.2.0"):
        _run_local_ulm(adata, network_path, science, decoupler_module=wrong_version)

    wrong_signature = _fake_decoupler(science)

    def obsolete_ulm(data, net):
        return None

    wrong_signature.mt.ulm = obsolete_ulm
    with pytest.raises(RuntimeError, match="public mt.ulm signature changed"):
        _run_local_ulm(adata, network_path, science, decoupler_module=wrong_signature)

    conflicting_path = tmp_path / "conflicting.csv"
    _write_signed_network(conflicting_path, conflicting_duplicate=True)
    with pytest.raises(ValueError, match="conflicting weights"):
        _run_local_ulm(adata, conflicting_path, science)

    with pytest.raises(ValueError, match="max_output_rows"):
        _run_local_ulm(adata, network_path, science, max_output_rows=10)
    with pytest.raises(MemoryError, match="working memory"):
        _run_local_ulm(adata, network_path, science, max_working_memory_gib=1e-12)


def test_collectri_ulm_rejects_bad_axes_nonfinite_and_output_collision(science, tmp_path):
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    duplicate_axis = _logged_adata(science)
    duplicate_axis.obs_names = ["duplicate"] * duplicate_axis.n_obs
    with pytest.raises(ValueError, match="unique"):
        _run_local_ulm(duplicate_axis, network_path, science)

    nonfinite = _logged_adata(science)
    nonfinite.X.data[0] = float("nan")
    with pytest.raises(ValueError, match="finite numeric"):
        _run_local_ulm(nonfinite, network_path, science)

    counts = _logged_adata(science)
    counts.X = science.sparse.csr_matrix(science.np.arange(counts.n_obs * counts.n_vars).reshape(counts.shape))
    counts.uns.pop("log1p", None)
    _output, _artifact, count_summary = _run_local_ulm(counts, network_path, science)
    assert any("count-like" in warning for warning in count_summary["warnings"])

    output, _, _ = _run_local_ulm(_logged_adata(science), network_path, science)
    with pytest.raises(ValueError, match="overwrite_existing"):
        _run_local_ulm(output, network_path, science)


def test_real_decoupler_2_2_collectri_ulm_toy_smoke_and_generated_parity(science, tmp_path):
    try:
        decoupler = importlib.import_module("decoupler")
    except (ImportError, OSError) as exc:
        pytest.skip(f"Real decoupler smoke unavailable in this interpreter: {exc}")
    if getattr(decoupler, "__version__", None) != "2.2.0":
        pytest.skip("Real CollecTRI ULM smoke requires exact decoupler 2.2.0.")
    network_path = tmp_path / "toy_collectri.csv"
    _write_signed_network(network_path)
    adata = _logged_adata(science)
    output, artifact, summary = _run_local_ulm(
        adata,
        network_path,
        science,
        decoupler_module=decoupler,
    )
    assert artifact.scores.shape == (adata.n_obs, 4)
    assert science.np.isfinite(artifact.scores.to_numpy()).all()
    assert summary["software_versions"]["decoupler"] == "2.2.0"
    parameters = {
        **summary["parameters"],
        "network_path": str(network_path),
        "resource_metadata_json": _resource_metadata_json(),
        "expected_resource_sha256": summary["key_results"]["resource"]["file_sha256"],
        "expected_canonical_network_sha256": summary["key_results"]["resource"]["canonical_network_sha256"],
        "openbio_version": "test-version",
    }
    parameters["source_kind"] = parameters.pop("source")
    namespace: dict[str, object] = {}
    exec(collectri_ulm_code(parameters=parameters), namespace)
    generated_output, generated_artifact, generated_summary = namespace["infer_collectri_ulm"](
        _logged_adata(science),
        decoupler_module=decoupler,
    )
    science.pd.testing.assert_frame_equal(
        generated_output.obsm[COLLECTRI_SCORE_KEY], output.obsm[COLLECTRI_SCORE_KEY]
    )
    science.pd.testing.assert_frame_equal(generated_artifact.scores, artifact.scores)
    assert generated_summary == summary


def _run_ranking(adata, artifact, science, *, decoupler_module=None, **overrides):
    parameters = {
        "annotation_key": "cell_type",
        "annotation_status": "curated",
        "reference": "rest",
        "method": "t-test_overestim_var",
        "report_p_adjusted": 0.05,
        "max_output_rows": 10_000,
        "openbio_version": "test-version",
    }
    parameters.update(overrides)
    return rank_tf_activities(
        adata,
        artifact,
        decoupler_module=_fake_decoupler(science) if decoupler_module is None else decoupler_module,
        **parameters,
    )


@pytest.mark.parametrize("method", ["t-test_overestim_var", "t-test", "wilcoxon"])
def test_tf_activity_ranking_all_methods_complete_bh_ties_and_generated_parity(
    science,
    method,
):
    adata, artifact = _ranking_input(science)
    before = _adata_snapshot(adata)
    artifact_scores_before = artifact.scores
    fake = _fake_decoupler(science)

    table, summary = _run_ranking(
        adata,
        artifact,
        science,
        decoupler_module=fake,
        method=method,
    )

    assert fake.calls == [
        {"operation": "get_obsm", "key": "openbio_tf_activity"},
        {
            "operation": "rankby_group",
            "groupby": "openbio_annotation",
            "reference": "rest",
            "method": method,
        },
    ]
    assert table.columns.tolist() == list(TF_RANKING_COLUMNS)
    assert len(table) == 4 * 4
    assert table["group"].drop_duplicates().tolist() == ["A", "B", "C", "REF"]
    assert set(table["regulator"]) == {"TF_UP", "TF_DOWN", "TF_TIE_A", "TF_TIE_B"}
    assert science.np.isfinite(
        table[["statistic", "mean_change", "p_value", "p_adjusted"]].to_numpy(dtype=float)
    ).all()
    assert table["significant_adjusted"].equals(table["p_adjusted"] <= 0.05)
    for group, group_frame in table.groupby("group", sort=False, observed=True):
        independent = _bh_adjust(group_frame["p_value"].to_numpy(), science.np)
        science.np.testing.assert_allclose(
            group_frame["p_adjusted"].to_numpy(),
            independent,
            rtol=1e-12,
            atol=1e-14,
            err_msg=f"complete BH family failed for {group}",
        )
        ties = group_frame.loc[group_frame["regulator"].isin(["TF_TIE_A", "TF_TIE_B"])]
        assert ties["rank_within_group"].nunique() == 1
    assert set(table["direction"]) <= {
        "higher_in_group",
        "lower_in_group",
        "no_signed_difference",
    }
    assert summary["status"] == "cluster_annotation_associated_tf_activity_evidence"
    assert summary["key_results"]["complete_family_rows"] == len(table)
    assert summary["key_results"]["unused_categorical_levels"] == ["unused"]
    assert summary["key_results"]["scientific_label"] == (
        "Cluster/annotation-associated TF activity evidence"
    )
    assert "Condition inference" not in summary["results"]
    assert len(summary["references"]) >= 4
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    _assert_adata_unchanged(adata, before, science)
    science.pd.testing.assert_frame_equal(artifact.scores, artifact_scores_before)

    parameters = {
        "annotation_key": "cell_type",
        "annotation_status": "curated",
        "reference": "rest",
        "method": method,
        "report_p_adjusted": 0.05,
        "max_output_rows": 10_000,
        "openbio_version": "test-version",
        "_portable_artifact": True,
    }
    code = rank_tf_activities_code(parameters=parameters)
    compile(code, "<rank-tf-activities-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_table, generated_summary = namespace["rank_tf_activity_artifact"](
        adata,
        artifact.portable(),
        decoupler_module=_fake_decoupler(science),
    )
    science.pd.testing.assert_frame_equal(generated_table, table)
    assert generated_summary == summary


def test_worker_owned_tf_ranking_reuses_decoded_activity_frame(science, monkeypatch):
    adata, artifact = _ranking_input(science)
    portable = artifact.portable()
    scores = portable["scores"]
    copied_shapes = []
    original_copy = science.pd.DataFrame.copy

    def tracked_copy(frame, *args, **kwargs):
        if frame is scores:
            copied_shapes.append(frame.shape)
        return original_copy(frame, *args, **kwargs)

    monkeypatch.setattr(science.pd.DataFrame, "copy", tracked_copy)
    rank_tf_activities(
        adata,
        portable,
        annotation_key="cell_type",
        annotation_status="curated",
        reference="rest",
        method="t-test_overestim_var",
        report_p_adjusted=0.05,
        max_output_rows=10_000,
        decoupler_module=_fake_decoupler(science),
        openbio_version="test-version",
        _portable_artifact=True,
        _worker_owned=True,
    )

    assert copied_shapes == []


@pytest.mark.parametrize(
    ("reference", "backend_reference", "expected_groups"),
    [
        ("REF", "REF", ["A", "B", "C"]),
        ("B, REF", ["B", "REF"], ["A", "C"]),
    ],
)
def test_tf_activity_ranking_single_and_list_references_exclude_reference_groups(
    science,
    reference,
    backend_reference,
    expected_groups,
):
    adata, artifact = _ranking_input(science)
    fake = _fake_decoupler(science)
    table, summary = _run_ranking(
        adata,
        artifact,
        science,
        decoupler_module=fake,
        reference=reference,
        method="t-test",
    )
    assert fake.calls[1]["reference"] == backend_reference
    assert table["group"].drop_duplicates().tolist() == expected_groups
    assert len(table) == len(expected_groups) * artifact.scores.shape[1]
    expected_reference_display = "REF" if reference == "REF" else "B, REF"
    assert set(table["reference"]) == {expected_reference_display}
    assert summary["key_results"]["tested_groups"] == expected_groups
    assert summary["key_results"]["reference_groups"] == (
        ["REF"] if reference == "REF" else ["B", "REF"]
    )
    if isinstance(backend_reference, list):
        assert set(table["group"]).isdisjoint(backend_reference)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("get_obsm_obs", "modified private carrier"),
        ("get_obsm_var", "modified private carrier"),
        ("rank_carrier_obs", "ranking modified private input"),
        ("rank_var", "ranking modified private input"),
        ("rank_uns", "ranking modified private input"),
        ("rank_drop_row", "complete expected group-regulator family"),
        ("rank_bad_padj", "differ from independent"),
        ("rank_nonfinite", "non-finite canonical statistics"),
    ],
)
def test_tf_activity_ranking_rejects_backend_tampering(science, mutation, message):
    adata, artifact = _ranking_input(science)
    with pytest.raises(RuntimeError, match=message):
        _run_ranking(
            adata,
            artifact,
            science,
            decoupler_module=_fake_decoupler(science, mutation=mutation),
        )


def test_tf_activity_ranking_exact_version_signatures_and_row_guard(science):
    adata, artifact = _ranking_input(science)
    wrong_version = _fake_decoupler(science)
    wrong_version.__version__ = "2.3.0"
    with pytest.raises(RuntimeError, match="exactly 2.2.0"):
        _run_ranking(adata, artifact, science, decoupler_module=wrong_version)

    wrong_get_obsm = _fake_decoupler(science)
    wrong_get_obsm.pp.get_obsm = lambda adata: adata
    with pytest.raises(RuntimeError, match="pp.get_obsm signature changed"):
        _run_ranking(adata, artifact, science, decoupler_module=wrong_get_obsm)

    wrong_rank = _fake_decoupler(science)
    wrong_rank.tl.rankby_group = lambda adata, groupby: science.pd.DataFrame()
    with pytest.raises(RuntimeError, match="tl.rankby_group signature changed"):
        _run_ranking(adata, artifact, science, decoupler_module=wrong_rank)

    guard_fake = _fake_decoupler(science)
    with pytest.raises(ValueError, match="max_output_rows"):
        _run_ranking(
            adata,
            artifact,
            science,
            decoupler_module=guard_fake,
            max_output_rows=10,
        )
    assert guard_fake.calls == []


@pytest.mark.parametrize(
    ("reference", "message"),
    [
        ("missing", "not observed"),
        ("B, B", "distinct"),
        ("A, B, C, REF", "every observed group"),
        ("B, ", "blank entry"),
    ],
)
def test_tf_activity_ranking_rejects_invalid_reference_design(science, reference, message):
    adata, artifact = _ranking_input(science)
    with pytest.raises(ValueError, match=message):
        _run_ranking(adata, artifact, science, reference=reference)


def test_tf_activity_ranking_rejects_misaligned_colliding_missing_and_small_groups(science):
    adata, artifact = _ranking_input(science)
    misaligned = adata[adata.obs_names[::-1]].copy()
    with pytest.raises(ValueError, match="do not align exactly"):
        _run_ranking(misaligned, artifact, science)

    colliding = adata.copy()
    colliding.obs["cell_type"] = [1] * 4 + ["1"] * 4 + ["C"] * 4 + ["REF"] * 4
    with pytest.raises(ValueError, match="collide after display normalization"):
        _run_ranking(colliding, artifact, science)

    missing = adata.copy()
    missing.obs["cell_type"] = missing.obs["cell_type"].astype(object)
    missing.obs.iloc[0, missing.obs.columns.get_loc("cell_type")] = None
    with pytest.raises(ValueError, match="missing labels"):
        _run_ranking(missing, artifact, science)

    small = adata.copy()
    small.obs["cell_type"] = ["tiny"] + ["reference"] * (small.n_obs - 1)
    with pytest.raises(ValueError, match="at least two cells"):
        _run_ranking(small, artifact, science, reference="rest")


def test_tf_activity_ranking_rejects_degenerate_statistics_and_tampered_artifact(science):
    adata, artifact = _ranking_input(science)
    scores = artifact.scores
    scores["TF_UP"] = 1.0
    degenerate = build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=artifact.adjusted_pvalues,
        provenance=artifact.provenance,
        numpy=science.np,
        pandas=science.pd,
    )
    with pytest.raises(ValueError, match="non-finite"):
        _run_ranking(adata, degenerate, science, method="t-test")

    private_scores = object.__getattribute__(artifact, "_scores")
    private_scores.iloc[0, 0] = 999.0
    fake = _fake_decoupler(science)
    with pytest.raises(ValueError, match="fingerprint"):
        _run_ranking(adata, artifact, science, decoupler_module=fake)
    assert fake.calls == []


def test_tf_activity_artifact_is_exact_typed_defensive_and_strict_json(science):
    scores, adjusted = _activity_frames(science)
    provenance = _artifact_provenance()
    artifact = build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance=provenance,
        numpy=science.np,
        pandas=science.pd,
    )

    assert type(artifact) is TFActivityArtifact
    assert artifact.artifact_type == TF_ACTIVITY_ARTIFACT_TYPE
    assert len(artifact.fingerprint) == 64
    json.dumps(artifact.provenance, ensure_ascii=False, allow_nan=False)
    json.dumps(artifact.metadata, ensure_ascii=False, allow_nan=False)

    scores.iloc[0, 0] = 999.0
    adjusted.iloc[0, 0] = 0.99
    provenance["method"] = "changed"
    assert artifact.scores.iloc[0, 0] == 3.0
    assert artifact.adjusted_pvalues.iloc[0, 0] == 0.01
    assert artifact.provenance["method"] == "CollecTRI ULM"

    returned_scores = artifact.scores
    returned_adjusted = artifact.adjusted_pvalues
    returned_provenance = artifact.provenance
    returned_metadata = artifact.metadata
    returned_scores.iloc[0, 0] = -999.0
    returned_adjusted.iloc[0, 0] = 0.9
    returned_provenance["method"] = "changed"
    returned_metadata["schema_version"] = 999
    assert artifact.scores.iloc[0, 0] == 3.0
    assert artifact.adjusted_pvalues.iloc[0, 0] == 0.01
    assert artifact.provenance["method"] == "CollecTRI ULM"
    assert artifact.metadata["schema_version"] == 1

    with pytest.raises(AttributeError, match="immutable"):
        artifact.extra = "forbidden"

    validated_scores, validated_adjusted, validated_provenance, validated_metadata = (
        validate_tf_activity_artifact(artifact, numpy=science.np, pandas=science.pd)
    )
    science.pd.testing.assert_frame_equal(validated_scores, artifact.scores)
    science.pd.testing.assert_frame_equal(validated_adjusted, artifact.adjusted_pvalues)
    assert validated_provenance == artifact.provenance
    assert validated_metadata == artifact.metadata


def test_owned_tf_activity_fingerprint_does_not_request_full_frame_copies(
    science,
    monkeypatch,
):
    scores, adjusted = _activity_frames(science)
    copy_requests = []
    original_to_numpy = science.pd.DataFrame.to_numpy

    def tracked_to_numpy(frame, *args, **kwargs):
        if kwargs.get("copy") is True:
            copy_requests.append(frame.shape)
        return original_to_numpy(frame, *args, **kwargs)

    monkeypatch.setattr(science.pd.DataFrame, "to_numpy", tracked_to_numpy)
    artifact = build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance=_artifact_provenance(),
        numpy=science.np,
        pandas=science.pd,
        copy_frames=False,
    )
    validate_tf_activity_artifact(
        artifact,
        numpy=science.np,
        pandas=science.pd,
        copy_result=False,
    )

    assert copy_requests == []


def test_tf_activity_artifact_detects_private_content_and_provenance_tampering(science):
    artifact = _activity_artifact(science)
    private_scores = object.__getattribute__(artifact, "_scores")
    private_scores.iloc[0, 0] = 99.0
    with pytest.raises(ValueError, match="fingerprint"):
        validate_tf_activity_artifact(artifact, numpy=science.np, pandas=science.pd)

    artifact = _activity_artifact(science)
    private_provenance = object.__getattribute__(artifact, "_provenance")
    private_provenance["method"] = "tampered"
    with pytest.raises(ValueError, match="provenance fingerprint"):
        validate_tf_activity_artifact(artifact, numpy=science.np, pandas=science.pd)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda scores, adjusted: scores.rename(index={"cell_0": " cell_0"}, inplace=True), "canonical"),
        (lambda scores, adjusted: scores.rename(columns={"TF_A": "TF_B"}, inplace=True), "unique"),
        (lambda scores, adjusted: scores.__setitem__("TF_A", float("nan")), "non-finite"),
        (lambda scores, adjusted: adjusted.__setitem__("TF_A", 1.1), r"\[0, 1\]"),
        (lambda scores, adjusted: adjusted.rename(columns={"TF_A": "TF_X"}, inplace=True), "regulator axes differ"),
    ],
)
def test_tf_activity_artifact_rejects_invalid_axes_and_values(science, mutate, message):
    scores, adjusted = _activity_frames(science)
    mutate(scores, adjusted)
    with pytest.raises((TypeError, ValueError), match=message):
        build_tf_activity_artifact(
            scores=scores,
            adjusted_pvalues=adjusted,
            provenance=_artifact_provenance(),
            numpy=science.np,
            pandas=science.pd,
        )


def test_tf_activity_runtime_validation_rejects_portable_substitutes(science):
    artifact = _activity_artifact(science)
    portable = artifact.portable()
    with pytest.raises(TypeError, match="exact OPENBIO_TF_ACTIVITY"):
        validate_tf_activity_artifact(portable, numpy=science.np, pandas=science.pd)

    scores, adjusted, provenance, metadata = validate_tf_activity_artifact(
        portable,
        exact_type=False,
        numpy=science.np,
        pandas=science.pd,
    )
    science.pd.testing.assert_frame_equal(scores, artifact.scores)
    science.pd.testing.assert_frame_equal(adjusted, artifact.adjusted_pvalues)
    assert provenance == artifact.provenance
    assert metadata == artifact.metadata
