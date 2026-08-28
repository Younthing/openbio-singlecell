from __future__ import annotations

import json
import math
from types import SimpleNamespace

import pytest

import openbio_singlecell.cnv_analysis as cnv_analysis_module
from openbio_singlecell.cnv_analysis import (
    CNV_STATE_TYPE,
    CNVState,
    analyze_cnv_pca,
    analyze_cnv_score,
    analyze_infer_cnv,
    cnv_pca_code,
    cnv_score_code,
    infer_cnv_code,
)
from openbio_singlecell.node_types import AnnDataType, CNVStateType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_cnv import (
    OpenBioSingleCellCNVPCA,
    OpenBioSingleCellCNVScore,
    OpenBioSingleCellInferCNV,
)


def _cnv_adata(science, *, sparse=False):
    rng = science.np.random.default_rng(42)
    n_cells, n_genes = 48, 40
    counts = rng.poisson(3.0, size=(n_cells, n_genes)).astype(float)
    counts[24:, 20:] += rng.poisson(4.0, size=(24, 20))
    cells = [f"cell_{index:03d}" for index in range(n_cells)]
    genes = [f"gene_{index:03d}" for index in range(n_genes)]
    obs = science.pd.DataFrame(
        {
            "cell_type": science.pd.Categorical(["normal"] * 24 + ["target"] * 24),
            "sample": ["s1"] * 12 + ["s2"] * 12 + ["s3"] * 12 + ["s4"] * 12,
            "partition": science.pd.Categorical(["A"] * 24 + ["B"] * 24, categories=["A", "B", "unused"]),
        },
        index=cells,
    )
    positions = science.np.tile(science.np.arange(20, dtype="int64") * 100, 2)
    var = science.pd.DataFrame(
        {
            "chromosome": ["chr1"] * 20 + ["chr2"] * 20,
            "start": positions,
            "end": positions + 50,
        },
        index=genes,
    )
    count_matrix = science.sparse.csr_matrix(counts) if sparse else counts
    logged = science.np.log1p(counts / counts.sum(axis=1, keepdims=True) * 10_000.0)
    logged_matrix = science.sparse.csr_matrix(logged) if sparse else logged
    adata = science.ad.AnnData(count_matrix.copy(), obs=obs, var=var)
    adata.layers["counts"] = count_matrix.copy()
    adata.raw = adata.copy()
    adata.X = logged_matrix.copy()
    adata.layers["log1p_norm"] = logged_matrix.copy()
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {"destination": "layer_and_raw", "layer_name": "counts"},
            },
            "000001": {
                "operation": "normalize_to_layer",
                "parameters": {"output_layer": "log1p_norm", "transform": "log1p"},
            },
        },
    }
    adata.uns["sentinel"] = {"unchanged": [1, 2, 3]}
    return adata


def _install_fake_backend(
    monkeypatch,
    science,
    *,
    malicious_pca=False,
    malicious_score=False,
    patch_installed_module=False,
):
    np = science.np

    def infercnv(
        adata,
        reference_key=None,
        reference_cat=None,
        reference=None,
        lfc_clip=3,
        window_size=100,
        step=10,
        dynamic_threshold=1.5,
        exclude_chromosomes=("chrX", "chrY"),
        chunksize=5000,
        n_jobs=None,
        inplace=True,
        layer=None,
        key_added="cnv",
        calculate_gene_values=False,
    ):
        del reference, lfc_clip, dynamic_threshold, chunksize, n_jobs, layer, calculate_gene_values
        dense = adata.X.toarray() if science.sparse.issparse(adata.X) else np.asarray(adata.X)
        reference_mask = np.asarray(adata.obs[reference_key].isin(reference_cat))
        columns = []
        chr_pos = {}
        for chromosome in sorted(set(adata.var["chromosome"]), key=lambda value: int(value.removeprefix("chr"))):
            if chromosome in exclude_chromosomes:
                continue
            indices = np.flatnonzero(np.asarray(adata.var["chromosome"] == chromosome))
            chr_pos[chromosome] = len(columns)
            window_count = math.ceil((len(indices) - window_size + 1) / step)
            for window_index in range(window_count):
                start = window_index * step
                selected = indices[start : start + window_size]
                signal = dense[:, selected].mean(axis=1)
                signal -= float(dense[reference_mask][:, selected].mean())
                columns.append(signal)
        matrix = np.column_stack(columns).astype(float)
        if inplace:
            adata.obsm[f"X_{key_added}"] = matrix
            adata.uns[key_added] = {"chr_pos": chr_pos}
            return None
        return matrix

    def pca(
        adata,
        svd_solver="arpack",
        zero_center=False,
        inplace=True,
        use_rep="cnv",
        key_added="cnv_pca",
        **kwargs,
    ):
        del svd_solver, zero_center
        matrix = adata.obsm[f"X_{use_rep}"]
        matrix = matrix.toarray() if science.sparse.issparse(matrix) else np.asarray(matrix)
        n_comps = int(kwargs["n_comps"])
        if malicious_pca:
            scores = np.arange(matrix.shape[0] * n_comps, dtype=float).reshape(matrix.shape[0], n_comps)
        else:
            u, singular_values, _vh = np.linalg.svd(matrix, full_matrices=False)
            scores = (u[:, :n_comps] * singular_values[:n_comps]).astype("float32")
        if inplace:
            adata.obsm[f"X_{key_added}"] = scores
            return None
        return scores

    def cnv_score(adata, groupby="cnv_leiden", use_rep="cnv", key_added="cnv_score", inplace=True, obs_key=None):
        del obs_key
        matrix = adata.obsm[f"X_{use_rep}"]
        scores = {}
        for category in adata.obs[groupby].unique():
            selected = matrix[np.asarray(adata.obs[groupby] == category), :]
            value = float(np.abs(selected).mean())
            scores[category] = value + (1.0 if malicious_score else 0.0)
        if inplace:
            adata.obs[key_added] = [scores[value] for value in adata.obs[groupby]]
            return None
        return scores

    if patch_installed_module:
        infercnvpy = pytest.importorskip("infercnvpy")
        monkeypatch.setattr(infercnvpy.tl, "infercnv", infercnv)
        monkeypatch.setattr(infercnvpy.tl, "pca", pca)
        monkeypatch.setattr(infercnvpy.tl, "cnv_score", cnv_score)
    else:
        infercnvpy = SimpleNamespace(tl=SimpleNamespace(infercnv=infercnv, pca=pca, cnv_score=cnv_score))
        monkeypatch.setattr(
            cnv_analysis_module,
            "_cnv_imports",
            lambda _operation: (
                infercnvpy,
                science.ad,
                science.np,
                science.pd,
                science.sparse,
                SimpleNamespace(natsorted=lambda values: sorted(values, key=lambda value: int(value.removeprefix("chr")))),
                "0.6.1",
            ),
        )
    return SimpleNamespace(infercnv=infercnv, pca=pca, cnv_score=cnv_score)


def _infer_parameters():
    return {
        "source_kind": "layer",
        "layer_name": "log1p_norm",
        "reference_key": "cell_type",
        "reference_categories": "normal",
        "sample_key": "sample",
        "genome_assembly": "GRCh38",
        "window_size": 10,
        "step": 5,
        "lfc_clip": 3.0,
        "dynamic_threshold": 1.5,
        "exclude_chromosomes": "chrX,chrY,chrM",
        "output_key": "cnv",
        "minimum_reference_cells": 20,
        "chunksize": 100,
        "n_jobs": 1,
        "max_output_gib": 1.0,
        "overwrite_existing": False,
    }


def test_cnv_schemas_are_atomic_and_typed():
    infer_schema = OpenBioSingleCellInferCNV.GET_SCHEMA()
    assert [item.id for item in infer_schema.inputs][:7] == [
        "adata",
        "source",
        "reference_key",
        "reference_categories",
        "sample_key",
        "genome_assembly",
        "window_size",
    ]
    assert [(item.display_name, item.io_type) for item in infer_schema.outputs] == [
        ("cnv_state", CNVStateType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert [(item.display_name, item.io_type) for item in OpenBioSingleCellCNVPCA.GET_SCHEMA().outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    score_schema = OpenBioSingleCellCNVScore.GET_SCHEMA()
    assert [item.id for item in score_schema.inputs] == [
        "cnv_state",
        "adata",
        "groupby",
        "output_key",
        "overwrite_existing",
    ]
    assert [(item.display_name, item.io_type) for item in score_schema.outputs] == [
        ("adata", AnnDataType.io_type),
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


@pytest.mark.parametrize("sparse", [False, True])
def test_cnv_fake_backend_chain_is_typed_reportable_and_immutable(science, monkeypatch, sparse):
    _install_fake_backend(monkeypatch, science)
    adata = _cnv_adata(science, sparse=sparse)
    before = adata.copy()
    state, infer_summary = analyze_infer_cnv(adata, **_infer_parameters())
    assert state.artifact_type == CNV_STATE_TYPE
    assert state.to_adata().obsm["X_cnv"].shape == (48, 6)
    assert infer_summary["key_results"]["reference"]["reference_samples"] == 2
    assert infer_summary["key_results"]["cnv_matrix"]["window_metadata"] == {"chr1": 0, "chr2": 3}
    assert infer_summary["references"] and "infercnvpy" in infer_summary["software_versions"]
    json.dumps(infer_summary, allow_nan=False)
    science.np.testing.assert_allclose(adata.layers["log1p_norm"].toarray() if sparse else adata.layers["log1p_norm"], before.layers["log1p_norm"].toarray() if sparse else before.layers["log1p_norm"])
    state_copy = state.to_adata()
    state_copy.obsm["X_cnv"][0, 0] += 100
    assert state.to_adata().obsm["X_cnv"][0, 0] != state_copy.obsm["X_cnv"][0, 0]

    pca, pca_summary = analyze_cnv_pca(state, n_comps=2, random_seed=7)
    assert pca.obsm["X_cnv_pca"].shape == (48, 2)
    assert pca_summary["key_results"]["output"]["independent_scientific_validation"]
    pca.obs["partition"] = pca.obs["partition"].cat.set_categories(["A", "B", "unused"])
    scored, table, score_summary = analyze_cnv_score(state, pca, groupby="partition")
    assert table["group_display"].tolist() == ["B", "A"]
    assert table["cell_count"].tolist() == [24, 24]
    assert score_summary["key_results"]["input_cnv_state"]["downstream_adata_fingerprints_verified"] is True
    assert "unused" not in scored.obs["partition"].cat.categories
    assert "unused" in [item["display"] for item in score_summary["key_results"]["partition"]["unused_categories_removed"]]
    json.dumps(pca_summary, allow_nan=False)
    json.dumps(score_summary, allow_nan=False)


def test_cnv_rejects_invalid_sources_coordinates_and_reference_design(science, monkeypatch):
    _install_fake_backend(monkeypatch, science)
    parameters = _infer_parameters()
    adata = _cnv_adata(science)
    adata.var["start"] = adata.var["start"].astype(float)
    with pytest.raises(TypeError, match="integer coordinates"):
        analyze_infer_cnv(adata, **parameters)

    adata = _cnv_adata(science)
    adata.var.iloc[1, adata.var.columns.get_loc("start")] = adata.var.iloc[0]["start"]
    adata.var.iloc[1, adata.var.columns.get_loc("end")] = adata.var.iloc[0]["end"]
    with pytest.raises(ValueError, match="duplicate chromosome/start/end"):
        analyze_infer_cnv(adata, **parameters)

    adata = _cnv_adata(science)
    with pytest.raises(ValueError, match="at least 30 reference cells"):
        analyze_infer_cnv(adata, **(parameters | {"minimum_reference_cells": 30}))


def test_cnv_artifact_and_downstream_fingerprint_tampering_fail_closed(science, monkeypatch):
    _install_fake_backend(monkeypatch, science)
    state, _summary = analyze_infer_cnv(_cnv_adata(science), **_infer_parameters())
    forged_adata = state.to_adata()
    forged_adata.obsm["X_cnv"][0, 0] += 1.0
    forged = CNVState(forged_adata, state.metadata)
    with pytest.raises(ValueError, match="content fingerprint"):
        analyze_cnv_pca(forged, n_comps=2)

    downstream = state.to_adata()
    downstream.obs["partition"] = downstream.obs["partition"].cat.remove_unused_categories()
    downstream.layers["log1p_norm"][0, 0] += 0.5
    with pytest.raises(ValueError, match="does not match"):
        analyze_cnv_score(state, downstream, groupby="partition")


def test_cnv_independent_postconditions_reject_scientifically_unrelated_results(science, monkeypatch):
    _install_fake_backend(monkeypatch, science, malicious_pca=True)
    state, _summary = analyze_infer_cnv(_cnv_adata(science), **_infer_parameters())
    with pytest.raises(RuntimeError, match="orthogonal|eigenpair|spectrum"):
        analyze_cnv_pca(state, n_comps=2)

    _install_fake_backend(monkeypatch, science, malicious_score=True)
    downstream = state.to_adata()
    downstream.obs["partition"] = downstream.obs["partition"].cat.remove_unused_categories()
    with pytest.raises(RuntimeError, match="independently recomputed"):
        analyze_cnv_score(state, downstream, groupby="partition")


def test_cnv_generated_code_compiles_has_no_openbio_import_and_matches(science, monkeypatch):
    _install_fake_backend(monkeypatch, science, patch_installed_module=True)
    adata = _cnv_adata(science)
    infer_parameters = _infer_parameters()
    expected_state, expected_summary = analyze_infer_cnv(adata, **infer_parameters)
    infer_source = infer_cnv_code(**infer_parameters)
    assert "from openbio_singlecell" not in infer_source
    infer_namespace = {}
    exec(compile(infer_source, "<infer-cnv-code>", "exec"), infer_namespace)
    observed_state, observed_summary = infer_namespace["run_infer_cnv"](adata)
    assert observed_state.fingerprint == expected_state.fingerprint
    assert observed_summary == expected_summary

    pca_parameters = {
        "n_comps": 2,
        "output_key": "X_cnv_pca",
        "overwrite_existing": False,
        "max_output_gib": 1.0,
        "random_seed": 7,
    }
    expected_pca, expected_pca_summary = analyze_cnv_pca(expected_state, **pca_parameters)
    pca_source = cnv_pca_code(**pca_parameters)
    assert "from openbio_singlecell" not in pca_source
    pca_namespace = {}
    exec(compile(pca_source, "<cnv-pca-code>", "exec"), pca_namespace)
    observed_pca, observed_pca_summary = pca_namespace["run_cnv_pca"](expected_state)
    science.np.testing.assert_allclose(observed_pca.obsm["X_cnv_pca"], expected_pca.obsm["X_cnv_pca"])
    assert observed_pca_summary == expected_pca_summary

    expected_pca.obs["partition"] = expected_pca.obs["partition"].cat.remove_unused_categories()
    score_parameters = {"groupby": "partition", "output_key": "cnv_score", "overwrite_existing": False}
    expected_scored, expected_table, expected_score_summary = analyze_cnv_score(
        expected_state, expected_pca, **score_parameters
    )
    score_source = cnv_score_code(**score_parameters)
    assert "from openbio_singlecell" not in score_source
    score_namespace = {}
    exec(compile(score_source, "<cnv-score-code>", "exec"), score_namespace)
    observed_scored, observed_table, observed_score_summary = score_namespace["run_cnv_score"](
        expected_state, expected_pca
    )
    science.pd.testing.assert_frame_equal(observed_table, expected_table)
    science.np.testing.assert_allclose(observed_scored.obs["cnv_score"], expected_scored.obs["cnv_score"])
    assert observed_score_summary == expected_score_summary


def test_real_infercnvpy_061_smoke(science, monkeypatch):
    infercnvpy = pytest.importorskip("infercnvpy")
    if getattr(infercnvpy, "__version__", "0.6.1") not in {"0.6.1", None}:
        pytest.skip("The audited infercnvpy 0.6.1 backend is not installed.")
    import infercnvpy.tl._infercnv as infer_module

    monkeypatch.setattr(
        infer_module,
        "process_map",
        lambda function, *iterables, **_kwargs: list(map(function, *iterables)),
    )
    state, summary = analyze_infer_cnv(_cnv_adata(science), **_infer_parameters())
    assert state.to_adata().obsm["X_cnv"].shape == (48, 6)
    assert summary["software_versions"]["infercnvpy"] == "0.6.1"
    pca, pca_summary = analyze_cnv_pca(state, n_comps=2)
    assert pca.obsm["X_cnv_pca"].shape == (48, 2)
    assert pca_summary["key_results"]["output"]["independent_scientific_validation"]
