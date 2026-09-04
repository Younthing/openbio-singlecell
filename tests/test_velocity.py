from __future__ import annotations

import copy
import json
import os
import random
import sys
import uuid
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from openbio_singlecell.artifact_codecs import read_plot, read_table, write_anndata
from openbio_singlecell.nodes_velocity import VELOCITY_NODE_CLASSES
from openbio_singlecell.operations_velocity import (
    estimate_velocity,
    recover_dynamics,
    velocity_filter_and_normalize,
    velocity_gene_ranking,
    velocity_graph,
    velocity_moments,
    velocity_stream_plot,
)
from openbio_singlecell.staged_state_codec import (
    VELOCITY_STATE_CODEC,
    read_velocity_state,
    write_velocity_state,
)
from openbio_singlecell.velocity_analysis import (
    VelocityState,
    run_velocity_estimate,
    run_velocity_graph,
    run_velocity_moments,
    run_velocity_prepare,
    run_velocity_ranking,
    run_velocity_recover,
    run_velocity_stream,
    validate_velocity_state,
    velocity_code,
)
from openbio_singlecell.velocity_portable import VELOCITY_STATE_KEY
from openbio_singlecell.worker_protocol import OperationContext


def _row_normalize(matrix):
    totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
    target = float(np.median(totals))
    factors = target / totals
    if sparse.issparse(matrix):
        return (sparse.diags(factors) @ sparse.csr_matrix(matrix)).astype(np.float32)
    return (np.asarray(matrix, dtype=np.float32) * factors[:, None]).astype(np.float32)


def _normalize_per_cell(
    data,
    counts_per_cell_after=None,
    counts_per_cell=None,
    key_n_counts=None,
    max_proportion_per_cell=None,
    use_initial_size=True,
    layers=None,
    enforce=None,
    copy=False,
):
    assert counts_per_cell_after is counts_per_cell is key_n_counts is None
    assert max_proportion_per_cell is None
    assert use_initial_size is enforce is True
    assert layers == ["spliced", "unspliced"]
    assert copy is False
    data.X = _row_normalize(data.X)
    for key in layers:
        data.layers[key] = _row_normalize(data.layers[key])
    data.obs["n_counts"] = np.asarray(data.X.sum(axis=1), dtype=float).ravel()
    np.random.seed(987)
    random.seed(987)
    np.set_printoptions(precision=1)


def _moments(
    data,
    n_neighbors=None,
    n_pcs=None,
    mode="connectivities",
    method="umap",
    use_rep=None,
    use_highly_variable=True,
    copy=False,
):
    assert n_neighbors is n_pcs is use_rep is None
    assert mode in {"connectivities", "distances"}
    assert method == "umap" and use_highly_variable is True and copy is False
    graph = sparse.csr_matrix(data.obsp[mode]) + sparse.eye(data.n_obs, format="csr")
    graph = sparse.diags(1.0 / np.asarray(graph.sum(axis=1)).ravel()) @ graph
    for source, target in (("spliced", "Ms"), ("unspliced", "Mu")):
        values = graph @ data.layers[source]
        data.layers[target] = np.asarray(
            values.toarray() if sparse.issparse(values) else values, dtype=np.float32
        )


def _velocity(
    data,
    vkey="velocity",
    mode="stochastic",
    fit_offset=False,
    fit_offset2=False,
    filter_genes=False,
    groups=None,
    groupby=None,
    groups_for_fit=None,
    constrain_ratio=None,
    use_raw=False,
    use_latent_time=None,
    perc=None,
    min_r2=0.01,
    min_likelihood=0.001,
    r2_adjusted=None,
    use_highly_variable=True,
    diff_kinetics=None,
    copy=False,
    **kwargs,
):
    assert fit_offset is fit_offset2 is filter_genes is use_raw is copy is False
    assert groups is groupby is groups_for_fit is constrain_ratio is None
    assert use_latent_time is perc is r2_adjusted is diff_kinetics is None
    assert use_highly_variable is True and not kwargs
    data.var["annotation"] = data.var["annotation"].astype("category")
    mask = np.ones(data.n_vars, dtype=bool)
    if mode == "dynamical":
        mask = data.var["openbio_dynamics_fit_success"].to_numpy(dtype=bool) & (
            data.var["fit_likelihood"].to_numpy(dtype=float) > min_likelihood
        )
    data.var[f"{vkey}_genes"] = mask
    if mode != "dynamical":
        data.var[f"{vkey}_offset"] = np.where(mask, 0.0, np.nan)
        data.var[f"{vkey}_offset2"] = np.where(mask, 0.0, np.nan)
        data.var[f"{vkey}_beta"] = np.where(mask, 0.4, np.nan)
        data.var[f"{vkey}_gamma"] = np.where(mask, 0.2, np.nan)
        data.var[f"{vkey}_qreg_ratio"] = np.where(mask, 0.5, np.nan)
        data.var[f"{vkey}_r2"] = np.where(mask, np.linspace(0.1, 0.9, data.n_vars), np.nan)
    velocity = data.layers["Mu"].astype(float) - 0.5 * data.layers["Ms"].astype(float)
    velocity[:, ~mask] = np.nan
    data.layers[vkey] = velocity
    data.uns[f"{vkey}_params"] = {"mode": mode, "fit_offset": False}
    data.obs["backend_temporary"] = 1


def _recover_dynamics(
    data,
    var_names="velocity_genes",
    n_top_genes=None,
    max_iter=10,
    assignment_mode="projection",
    t_max=None,
    fit_time=True,
    fit_scaling=True,
    fit_steady_states=True,
    fit_connected_states=True,
    fit_basal_transcription=False,
    use_raw=False,
    load_pars=False,
    return_model=None,
    plot_results=False,
    steady_state_prior=None,
    add_key="fit",
    copy=False,
    n_jobs=None,
    backend="loky",
    show_progress_bar=True,
    **kwargs,
):
    assert isinstance(var_names, list) and n_top_genes is None
    assert assignment_mode == "projection" and max_iter >= 1
    assert t_max is steady_state_prior is None
    assert fit_time is fit_scaling is fit_steady_states is fit_connected_states is True
    assert fit_basal_transcription is use_raw is load_pars is return_model is plot_results is copy is False
    assert add_key == "fit" and n_jobs >= 1 and backend == "loky"
    assert show_progress_bar is False and not kwargs
    selected = data.var_names.isin(var_names)
    fitted = selected.copy()
    selected_positions = np.flatnonzero(selected)
    if len(selected_positions) > 5:
        fitted[selected_positions[-1]] = False
    for index, key in enumerate(
        (
            "fit_alpha",
            "fit_beta",
            "fit_gamma",
            "fit_t_",
            "fit_scaling",
            "fit_std_u",
            "fit_std_s",
            "fit_likelihood",
            "fit_u0",
            "fit_s0",
            "fit_pval_steady",
            "fit_steady_u",
            "fit_steady_s",
            "fit_variance",
            "fit_alignment_scaling",
        ),
        start=1,
    ):
        values = np.full(data.n_vars, np.nan)
        values[fitted] = np.linspace(0.5, 2.5, fitted.sum()) if key == "fit_likelihood" else index / 10
        data.var[key] = values
    for key in ("fit_t", "fit_tau", "fit_tau_"):
        values = np.full(data.shape, np.nan)
        values[:, fitted] = np.arange(data.n_obs, dtype=float)[:, None] / max(data.n_obs - 1, 1)
        data.layers[key] = values
    data.varm["loss"] = np.full((data.n_vars, 2), np.nan)
    data.varm["loss"][fitted] = 0.1
    data.uns["recover_dynamics"] = {
        "fit_connected_states": True,
        "fit_basal_transcription": False,
        "use_raw": False,
    }


def _velocity_graph(
    data,
    vkey="velocity",
    xkey="Ms",
    tkey=None,
    basis=None,
    n_neighbors=None,
    n_recurse_neighbors=None,
    random_neighbors_at_max=None,
    sqrt_transform=False,
    variance_stabilization=None,
    gene_subset=None,
    compute_uncertainties=None,
    approx=None,
    mode_neighbors="distances",
    copy=False,
    n_jobs=None,
    backend="loky",
    show_progress_bar=True,
):
    assert xkey == "Ms" and tkey is basis is n_neighbors is None
    assert n_recurse_neighbors == 1
    assert random_neighbors_at_max is variance_stabilization is gene_subset is None
    assert sqrt_transform is compute_uncertainties is approx is copy is False
    assert mode_neighbors in {"distances", "connectivities"}
    assert n_jobs >= 1 and backend == "loky" and show_progress_bar is False
    rows = np.arange(data.n_obs - 1)
    data.uns[f"{vkey}_graph"] = sparse.csr_matrix(
        (np.full(data.n_obs - 1, 0.7), (rows, rows + 1)), shape=(data.n_obs, data.n_obs)
    )
    data.uns[f"{vkey}_graph_neg"] = sparse.csr_matrix(
        (np.full(data.n_obs - 1, -0.2), (rows + 1, rows)), shape=(data.n_obs, data.n_obs)
    )
    data.obs[f"{vkey}_self_transition"] = np.linspace(0.0, 1.0, data.n_obs)
    data.uns[f"{vkey}_params"]["mode_neighbors"] = mode_neighbors
    data.uns[f"{vkey}_params"]["n_recurse_neighbors"] = n_recurse_neighbors


def _velocity_embedding_stream(
    adata,
    basis=None,
    vkey="velocity",
    density=2,
    smooth=None,
    min_mass=None,
    cutoff_perc=None,
    arrow_color=None,
    arrow_size=1,
    arrow_style="-|>",
    max_length=4,
    integration_direction="both",
    linewidth=None,
    n_neighbors=None,
    recompute=None,
    color=None,
    use_raw=None,
    layer=None,
    color_map=None,
    colorbar=True,
    palette=None,
    size=None,
    alpha=0.3,
    perc=None,
    X=None,
    V=None,
    X_grid=None,
    V_grid=None,
    sort_order=True,
    groups=None,
    components=None,
    legend_loc="on data",
    legend_fontsize=None,
    legend_fontweight=None,
    xlabel=None,
    ylabel=None,
    title=None,
    fontsize=None,
    figsize=None,
    dpi=None,
    frameon=None,
    show=None,
    save=None,
    ax=None,
    ncols=None,
    **kwargs,
):
    assert basis == "umap" and vkey == "velocity"
    assert density > 0 and smooth > 0 and min_mass >= 0
    assert show is False and save is None and ax is not None
    assert not kwargs
    embedding = np.asarray(adata.obsm[f"X_{basis}"])
    ax.scatter(embedding[:, 0], embedding[:, 1], s=6)
    ax.plot(embedding[:, 0], embedding[:, 1], linewidth=0.5)
    adata.obsm[f"{vkey}_{basis}"] = np.zeros_like(embedding)


class _Settings:
    verbosity = 1
    presenter_view = False
    autoshow = True
    autosave = False
    figdir = "figures"
    plot_prefix = ""
    dpi = 80
    dpi_save = 150
    frameon = None
    vector_friendly = False
    file_format_figs = "png"
    transparent = False
    color_map = "viridis"


def _fake_scvelo():
    return SimpleNamespace(
        __version__="0.3.4",
        settings=_Settings(),
        pp=SimpleNamespace(normalize_per_cell=_normalize_per_cell, moments=_moments),
        tl=SimpleNamespace(
            velocity=_velocity,
            velocity_graph=_velocity_graph,
            recover_dynamics=_recover_dynamics,
        ),
        pl=SimpleNamespace(velocity_embedding_stream=_velocity_embedding_stream),
    )


def _adata(matrix_kind: str = "dense") -> AnnData:
    rng = np.random.default_rng(41)
    n_obs, n_vars = 24, 16
    spliced = rng.poisson(4, (n_obs, n_vars)).astype(np.int64) + 1
    unspliced = rng.poisson(2, (n_obs, n_vars)).astype(np.int64) + 1
    if matrix_kind == "csr":
        spliced, unspliced = sparse.csr_matrix(spliced), sparse.csr_matrix(unspliced)
    elif matrix_kind == "csc":
        spliced, unspliced = sparse.csc_matrix(spliced), sparse.csc_matrix(unspliced)
    adata = AnnData(
        X=np.log1p(np.asarray(spliced.toarray() if sparse.issparse(spliced) else spliced)),
        obs=pd.DataFrame(index=[f"cell-{index:02d}" for index in range(n_obs)]),
        var=pd.DataFrame(index=[f"gene-{index:02d}" for index in range(n_vars)]),
    )
    adata.var["annotation"] = [f"class-{index % 2}" for index in range(n_vars)]
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.layers["protected"] = np.full(adata.shape, 7.0)
    rows = np.arange(n_obs)
    cols = (rows + 1) % n_obs
    ring = sparse.csr_matrix((np.ones(n_obs), (rows, cols)), shape=(n_obs, n_obs))
    connectivities = (ring + ring.T).tocsr()
    distances = connectivities.copy()
    adata.obsp["vel_connectivities"] = connectivities
    adata.obsp["vel_distances"] = distances
    adata.uns["vel_neighbors"] = {
        "connectivities_key": "vel_connectivities",
        "distances_key": "vel_distances",
        "params": {"method": "umap", "n_neighbors": 2},
    }
    adata.obsm["X_umap"] = rng.normal(size=(n_obs, 2))
    adata.obs["leiden"] = pd.Categorical([str(index % 3) for index in range(n_obs)])
    return adata


def _matrix_equal(left, right) -> bool:
    difference = sparse.csr_matrix(left) - sparse.csr_matrix(right)
    difference.eliminate_zeros()
    return difference.nnz == 0


def _assert_adata_unchanged(observed: AnnData, expected: AnnData) -> None:
    assert observed.obs_names.equals(expected.obs_names)
    assert observed.var_names.equals(expected.var_names)
    pd.testing.assert_frame_equal(observed.obs, expected.obs)
    pd.testing.assert_frame_equal(observed.var, expected.var)
    assert _matrix_equal(observed.X, expected.X)
    assert set(observed.layers) == set(expected.layers)
    for key in observed.layers:
        assert _matrix_equal(observed.layers[key], expected.layers[key])
    assert set(observed.obsp) == set(expected.obsp)
    for key in observed.obsp:
        assert _matrix_equal(observed.obsp[key], expected.obsp[key])
    assert set(observed.uns) == set(expected.uns)


def _chain():
    fake = _fake_scvelo()
    source = _adata()
    prepared = run_velocity_prepare(source, min_shared_counts=0, scvelo_module=fake)
    moments = run_velocity_moments(
        prepared, neighbors_key="vel_neighbors", scvelo_module=fake
    )
    steady = run_velocity_estimate(moments, mode="deterministic", scvelo_module=fake)
    recovered = run_velocity_recover(
        steady, gene_selection="velocity_genes", scvelo_module=fake
    )
    dynamical = run_velocity_estimate(recovered, mode="dynamical", scvelo_module=fake)
    graph = run_velocity_graph(dynamical, scvelo_module=fake)
    return fake, source, prepared, moments, steady, recovered, dynamical, graph


def _operation_context(root: Path) -> OperationContext:
    root.mkdir()
    return OperationContext.from_request_path(root / "request.json", str(uuid.uuid4()))


def _artifact_descriptor(root: Path, *, kind: str, codec: str) -> dict[str, object]:
    return {"type": "artifact", "path": str(root.resolve()), "kind": kind, "codec": codec}


def _file_snapshot(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_velocity_worker_operations_chain_uses_new_artifacts_and_never_rewrites_inputs(
    monkeypatch, tmp_path: Path
):
    monkeypatch.setitem(sys.modules, "scvelo", _fake_scvelo())
    source_root = tmp_path / "source"
    source_root.mkdir()
    write_anndata(source_root, _adata())
    source_before = _file_snapshot(source_root)

    def run_state(operation, input_root: Path, parameters: dict[str, object], name: str):
        before = _file_snapshot(input_root)
        context = _operation_context(tmp_path / name)
        records = operation(
            context,
            {
                "adata" if operation is velocity_filter_and_normalize else "velocity_state":
                    _artifact_descriptor(
                        input_root,
                        kind=(
                            "OPENBIO_ANNDATA"
                            if operation is velocity_filter_and_normalize
                            else "OPENBIO_VELOCITY_STATE"
                        ),
                        codec=(
                            "anndata-h5ad-v1"
                            if operation is velocity_filter_and_normalize
                            else VELOCITY_STATE_CODEC
                        ),
                    )
            },
            parameters,
        )
        assert [record["name"] for record in records] == ["velocity_state", "summary", "code"]
        assert records[0]["codec"] == VELOCITY_STATE_CODEC
        assert _file_snapshot(input_root) == before
        output_root = context.output_root / records[0]["payload"]
        read_velocity_state(output_root)
        return output_root

    prepared = run_state(
        velocity_filter_and_normalize,
        source_root,
        {
            "spliced_layer": "spliced",
            "unspliced_layer": "unspliced",
            "min_shared_counts": 0,
            "min_shared_cells": 0,
            "normalization_target": "median_library",
            "overwrite_existing": False,
        },
        "prepare",
    )
    moments = run_state(
        velocity_moments,
        prepared,
        {
            "neighbors_key": "vel_neighbors",
            "mode": "connectivities",
            "max_dense_gib": 2.0,
            "overwrite_existing": False,
        },
        "moments",
    )
    steady = run_state(
        estimate_velocity,
        moments,
        {
            "mode": "deterministic",
            "vkey": "velocity",
            "min_r2": 0.01,
            "min_likelihood": 0.001,
            "overwrite_existing": False,
        },
        "estimate",
    )
    recovered = run_state(
        recover_dynamics,
        steady,
        {
            "gene_selection": "velocity_genes",
            "n_top_genes": 0,
            "max_iter": 10,
            "n_jobs": 1,
            "max_dense_gib": 2.0,
            "overwrite_existing": False,
        },
        "recover",
    )
    dynamical = run_state(
        estimate_velocity,
        recovered,
        {
            "mode": "dynamical",
            "vkey": "velocity",
            "min_r2": 0.01,
            "min_likelihood": 0.001,
            "overwrite_existing": True,
        },
        "dynamical",
    )
    graph = run_state(
        velocity_graph,
        dynamical,
        {
            "vkey": "velocity",
            "xkey": "Ms",
            "mode_neighbors": "distances",
            "n_jobs": 1,
            "overwrite_existing": False,
        },
        "graph",
    )

    ranking_before = _file_snapshot(recovered)
    ranking_context = _operation_context(tmp_path / "ranking")
    ranking_records = velocity_gene_ranking(
        ranking_context,
        {
            "velocity_state": _artifact_descriptor(
                recovered, kind="OPENBIO_VELOCITY_STATE", codec=VELOCITY_STATE_CODEC
            )
        },
        {"top_n": 8, "include_failed": True, "max_output_rows": 100_000},
    )
    assert [record["name"] for record in ranking_records] == ["table", "summary", "code"]
    table, _metadata = read_table(ranking_context.output_root / ranking_records[0]["payload"])
    assert not table.empty
    assert _file_snapshot(recovered) == ranking_before

    graph_before = _file_snapshot(graph)
    plot_context = _operation_context(tmp_path / "stream")
    plot_records = velocity_stream_plot(
        plot_context,
        {
            "velocity_state": _artifact_descriptor(
                graph, kind="OPENBIO_VELOCITY_STATE", codec=VELOCITY_STATE_CODEC
            )
        },
        {"basis": "umap", "color_key": "leiden", "density": 2.0, "smooth": 0.5, "min_mass": 1.0},
    )
    assert [record["name"] for record in plot_records] == ["plot", "summary", "code"]
    png, _metadata = read_plot(plot_context.output_root / plot_records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert _file_snapshot(graph) == graph_before
    assert _file_snapshot(source_root) == source_before


def test_velocity_state_codec_roundtrip_is_h5ad_plus_strict_json(tmp_path: Path):
    state = run_velocity_prepare(
        _adata(), min_shared_counts=0, scvelo_module=_fake_scvelo()
    )
    root = tmp_path / "velocity-state"
    root.mkdir()

    descriptors = write_velocity_state(root, state)

    assert VELOCITY_STATE_CODEC == "velocity-state-h5ad-json-v1"
    assert {path.name for path in root.iterdir()} == {"data.h5ad", "state.json"}
    assert {item["path"] for item in descriptors} == {"data.h5ad", "state.json"}
    json.loads((root / "state.json").read_text(encoding="utf-8"))
    restored = read_velocity_state(root)
    payload, summary, metadata, portable = validate_velocity_state(restored)
    assert restored.fingerprint == state.fingerprint
    assert summary == state.summary
    assert metadata == state.metadata
    assert portable["state_fingerprint_sha256"] == metadata["state_fingerprint_sha256"]
    assert payload.shape == state.portable_adata().shape


def test_velocity_nodes_are_schema_only_and_operation_module_is_worker_pure():
    assert all("execute" not in node.__dict__ for node in VELOCITY_NODE_CLASSES)
    assert [node.define_schema().node_id for node in VELOCITY_NODE_CLASSES] == [
        "OpenBioSingleCellVelocityFilterAndNormalize",
        "OpenBioSingleCellVelocityMoments",
        "OpenBioSingleCellEstimateVelocity",
        "OpenBioSingleCellVelocityGraph",
        "OpenBioSingleCellRecoverDynamics",
        "OpenBioSingleCellVelocityGeneRanking",
        "OpenBioSingleCellVelocityDynamicsPlot",
        "OpenBioSingleCellVelocityGeneRankingPlot",
        "OpenBioSingleCellVelocityStreamPlot",
    ]
    source = Path(__file__).parents[1] / "openbio_singlecell" / "operations_velocity.py"
    text = source.read_text(encoding="utf-8")
    assert "comfy_api" not in text
    assert "folder_paths" not in text
    assert "nodes_velocity" not in text


@pytest.mark.parametrize("matrix_kind", ["dense", "csr", "csc"])
def test_prepare_is_sparse_safe_strict_and_does_not_mutate_input(matrix_kind):
    adata = _adata(matrix_kind)
    before = adata.copy()
    state = run_velocity_prepare(
        adata,
        min_shared_counts=0,
        min_shared_cells=1,
        scvelo_module=_fake_scvelo(),
    )
    payload, summary, metadata, portable = validate_velocity_state(
        state, allowed_stages=("prepared",)
    )
    _assert_adata_unchanged(adata, before)
    assert isinstance(state, VelocityState)
    assert metadata["state_fingerprint_sha256"] == portable["state_fingerprint_sha256"]
    assert summary["software_versions"]["scvelo"] == "0.3.4"
    assert summary["key_results"]["expression_X_restored_exactly"] is True
    assert _matrix_equal(payload.X, before.X)
    assert _matrix_equal(payload.layers["protected"], before.layers["protected"])
    json.dumps(summary, allow_nan=False)


def test_complete_fake_chain_has_atomic_stages_views_and_private_plotting():
    fake, source, prepared, moments, steady, recovered, dynamical, graph = _chain()
    assert [item.stage for item in (prepared, moments, steady, recovered, dynamical, graph)] == [
        "prepared",
        "moments",
        "velocity_estimated",
        "dynamics_recovered",
        "velocity_estimated",
        "velocity_graph",
    ]
    steady_payload = steady.portable_adata()
    pd.testing.assert_series_equal(
        steady_payload.var["annotation"], source.var["annotation"]
    )
    assert {
        "velocity_offset",
        "velocity_offset2",
        "velocity_beta",
        "velocity_gamma",
        "velocity_qreg_ratio",
        "velocity_r2",
        "velocity_genes",
    }.issubset(steady_payload.var.columns)
    table, ranking_summary = run_velocity_ranking(
        recovered, top_n=8, include_failed=True
    )
    assert list(table.columns) == [
        "rank",
        "gene",
        "fit_status",
        "fit_likelihood",
        "fit_r2",
        "velocity_gene",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_switch_time",
    ]
    assert table.iloc[:8]["fit_likelihood"].is_monotonic_decreasing
    assert table.iloc[-1]["fit_status"] == "failed"
    json.dumps(ranking_summary, allow_nan=False)
    before_graph = graph.portable_adata()
    png, plot_summary = run_velocity_stream(graph, scvelo_module=fake)
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert plot_summary["key_results"]["private_embedding_cache_created"] is True
    assert "velocity_umap" not in graph.portable_adata().obsm
    assert "velocity_umap" not in before_graph.obsm
    _assert_adata_unchanged(source, _adata())


def test_artifact_is_defensive_immutable_and_detects_internal_tampering():
    state = run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=_fake_scvelo())
    with pytest.raises(AttributeError, match="immutable"):
        state.stage = "moments"
    payload = state.portable_adata()
    payload.layers["spliced"][0, 0] = 999
    validate_velocity_state(state)
    state._adata.layers["spliced"][0, 0] = 999
    with pytest.raises(ValueError, match="fingerprint"):
        validate_velocity_state(state)

    fake = _fake_scvelo()
    prepared = run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=fake)
    moments = run_velocity_moments(
        prepared, neighbors_key="vel_neighbors", scvelo_module=fake
    )
    estimate = run_velocity_estimate(moments, mode="deterministic", scvelo_module=fake)
    estimate._adata.var.loc[estimate._adata.var_names[0], "velocity_offset"] = 99.0
    with pytest.raises(ValueError, match="fingerprint"):
        validate_velocity_state(estimate)


def test_stage_memory_and_prerequisite_guards_fail_before_backend():
    fake = _fake_scvelo()
    source = _adata()
    prepared = run_velocity_prepare(source, min_shared_counts=0, scvelo_module=fake)
    with pytest.raises(MemoryError, match="max_dense_gib"):
        run_velocity_moments(
            prepared,
            neighbors_key="vel_neighbors",
            max_dense_gib=1e-12,
            scvelo_module=fake,
        )
    with pytest.raises(ValueError, match="stage"):
        run_velocity_estimate(prepared, scvelo_module=fake)
    with pytest.raises(ValueError, match="stage"):
        run_velocity_graph(prepared, scvelo_module=fake)
    moments = run_velocity_moments(
        prepared, neighbors_key="vel_neighbors", scvelo_module=fake
    )
    with pytest.raises(ValueError, match="velocity-gene mask"):
        run_velocity_recover(moments, gene_selection="velocity_genes", scvelo_module=fake)
    with pytest.raises(ValueError, match="stage"):
        run_velocity_ranking(moments)

    collision_source = _adata()
    collision_source.var["custom_offset"] = 0.0
    collision_prepared = run_velocity_prepare(
        collision_source, min_shared_counts=0, scvelo_module=fake
    )
    collision_moments = run_velocity_moments(
        collision_prepared, neighbors_key="vel_neighbors", scvelo_module=fake
    )
    with pytest.raises(ValueError, match="custom_offset"):
        run_velocity_estimate(
            collision_moments,
            mode="deterministic",
            vkey="custom",
            scvelo_module=fake,
        )


def test_backend_version_signature_fallback_and_malformed_output_guards():
    bad_version = _fake_scvelo()
    bad_version.__version__ = "0.3.3"
    with pytest.raises(RuntimeError, match="exact scVelo 0.3.4"):
        run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=bad_version)

    bad_signature = _fake_scvelo()
    bad_signature.pp.normalize_per_cell = lambda data: None
    with pytest.raises(RuntimeError, match="interface is incompatible"):
        run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=bad_signature)

    fake = _fake_scvelo()
    prepared = run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=fake)
    moments = run_velocity_moments(
        prepared, neighbors_key="vel_neighbors", scvelo_module=fake
    )
    fallback = _fake_scvelo()

    def fallback_velocity(*args, **kwargs):
        _velocity(*args, **kwargs)
        data = args[0] if args else kwargs["data"]
        data.uns["velocity_params"]["mode"] = "deterministic"

    fallback_velocity.__signature__ = __import__("inspect").signature(_velocity)
    fallback.tl.velocity = fallback_velocity
    with pytest.raises(RuntimeError, match="fell back"):
        run_velocity_estimate(moments, mode="stochastic", scvelo_module=fallback)

    steady = run_velocity_estimate(
        moments, mode="deterministic", scvelo_module=fake
    )
    recovered = run_velocity_recover(
        steady, gene_selection="velocity_genes", scvelo_module=fake
    )
    high_density_threshold = run_velocity_estimate(
        recovered,
        mode="dynamical",
        min_likelihood=1.1,
        scvelo_module=fake,
    )
    assert high_density_threshold.summary["key_results"]["selected_velocity_genes"] >= 10


def test_global_random_print_matplotlib_and_backend_settings_are_restored():
    import matplotlib

    fake = _fake_scvelo()
    np.random.seed(111)
    random.seed(222)
    np.set_printoptions(precision=6)
    matplotlib.rcParams["lines.linewidth"] = 3.25
    numpy_state = copy.deepcopy(np.random.get_state())
    python_state = random.getstate()
    print_options = np.get_printoptions()
    rc_value = matplotlib.rcParams["lines.linewidth"]
    setting_names = (
        "verbosity",
        "presenter_view",
        "autoshow",
        "autosave",
        "figdir",
        "plot_prefix",
        "dpi",
        "dpi_save",
        "frameon",
        "vector_friendly",
        "file_format_figs",
        "transparent",
        "color_map",
    )
    settings = {name: copy.deepcopy(getattr(fake.settings, name)) for name in setting_names}
    run_velocity_prepare(_adata(), min_shared_counts=0, scvelo_module=fake)
    observed_numpy = np.random.get_state()
    assert observed_numpy[0] == numpy_state[0]
    assert np.array_equal(observed_numpy[1], numpy_state[1])
    assert observed_numpy[2:] == numpy_state[2:]
    assert random.getstate() == python_state
    assert np.get_printoptions() == print_options
    assert matplotlib.rcParams["lines.linewidth"] == rc_value
    assert {name: getattr(fake.settings, name) for name in setting_names} == settings


def test_generated_sources_compile_are_standalone_and_match_portable_results():
    fake, _, prepared, moments, steady, recovered, dynamical, graph = _chain()
    cases = [
        (
            "prepare",
            _adata(),
            {"min_shared_counts": 0},
            "run_velocity_filter_and_normalize",
            prepared,
        ),
        (
            "moments",
            prepared.portable_adata(),
            {"neighbors_key": "vel_neighbors"},
            "run_velocity_moments",
            moments,
        ),
        (
            "estimate",
            moments.portable_adata(),
            {"mode": "deterministic"},
            "run_velocity_estimation",
            steady,
        ),
        (
            "recover",
            steady.portable_adata(),
            {"gene_selection": "velocity_genes"},
            "run_velocity_dynamics_recovery",
            recovered,
        ),
        (
            "estimate",
            recovered.portable_adata(),
            {"mode": "dynamical"},
            "run_velocity_estimation",
            dynamical,
        ),
        (
            "graph",
            dynamical.portable_adata(),
            {},
            "run_velocity_graph",
            graph,
        ),
    ]
    for operation, adata, parameters, function_name, expected in cases:
        source = velocity_code(operation, **parameters)
        assert "openbio_singlecell" not in source
        namespace = {}
        exec(compile(source, f"<{operation}>", "exec"), namespace)
        portable, summary = namespace[function_name](adata, scvelo_module=fake)
        expected_adata, expected_summary, _, _ = validate_velocity_state(expected)
        assert portable.uns[VELOCITY_STATE_KEY]["state_fingerprint_sha256"] == expected_adata.uns[
            VELOCITY_STATE_KEY
        ]["state_fingerprint_sha256"]
        assert summary == expected_summary

    ranking_source = velocity_code(
        "ranking", top_n=8, include_failed=True, max_output_rows=100_000
    )
    ranking_namespace = {}
    exec(compile(ranking_source, "<ranking>", "exec"), ranking_namespace)
    observed_table, observed_summary = ranking_namespace["run_velocity_fit_ranking"](
        recovered.portable_adata()
    )
    expected_table, expected_summary = run_velocity_ranking(
        recovered, top_n=8, include_failed=True
    )
    pd.testing.assert_frame_equal(observed_table, expected_table)
    assert observed_summary == expected_summary

    stream_source = velocity_code(
        "stream", basis="umap", color_key="leiden", density=2.0, smooth=0.5, min_mass=1.0
    )
    stream_namespace = {}
    exec(compile(stream_source, "<stream>", "exec"), stream_namespace)
    observed_png, observed_summary = stream_namespace["run_velocity_stream_plot"](
        graph.portable_adata(), scvelo_module=fake
    )
    expected_png, expected_summary = run_velocity_stream(graph, scvelo_module=fake)
    assert observed_png == expected_png
    assert observed_summary == expected_summary


def test_generated_sources_exclude_unrelated_velocity_operations():
    implementations = {
        "prepare": "prepare_velocity_abundances",
        "moments": "compute_velocity_moments",
        "estimate": "estimate_rna_velocity",
        "recover": "recover_velocity_dynamics",
        "graph": "build_velocity_graph",
        "ranking": "rank_recovered_dynamics",
        "stream": "render_velocity_stream",
    }
    for operation, implementation in implementations.items():
        source = velocity_code(operation)
        assert f"def {implementation}(" in source
        for unrelated in set(implementations.values()) - {implementation}:
            assert f"def {unrelated}(" not in source


def test_real_scvelo_034_full_cpu_smoke_when_wheel_is_provided(monkeypatch):
    wheel = os.environ.get("OPENBIO_SCVELO_WHEEL")
    if not wheel:
        pytest.skip("Set OPENBIO_SCVELO_WHEEL to the audited scVelo 0.3.4 wheel for the real smoke.")
    monkeypatch.syspath_prepend(wheel)
    sys.modules.pop("scvelo", None)
    import scanpy as sc
    import scvelo

    assert scvelo.__version__ == "0.3.4"
    adata = scvelo.datasets.simulation(n_obs=60, n_vars=24, random_seed=4)
    adata.layers["spliced"] = np.rint(adata.layers["spliced"]).astype(int) + 1
    adata.layers["unspliced"] = np.rint(adata.layers["unspliced"]).astype(int) + 1
    adata.X = np.asarray(adata.layers["spliced"], dtype=float)
    sc.pp.pca(adata, n_comps=12)
    sc.pp.neighbors(adata, n_neighbors=8, n_pcs=12, key_added="vel_neighbors")
    adata.obsm["X_umap"] = np.column_stack(
        [np.linspace(-1.0, 1.0, adata.n_obs), np.sin(np.linspace(0.0, 4.0, adata.n_obs))]
    )
    adata.obs["leiden"] = pd.Categorical(
        np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "early", "late")
    )
    recover_globals = scvelo.tl.recover_dynamics.__globals__
    original_unique = recover_globals["make_unique_list"]
    original_read = recover_globals["_read_pars"]
    divergence_module = __import__(
        "scvelo.tools._em_model_utils", fromlist=["compute_divergence"]
    )
    original_divergence = divergence_module.compute_divergence
    prepared = run_velocity_prepare(
        adata, min_shared_counts=0, scvelo_module=scvelo
    )
    moments = run_velocity_moments(
        prepared, neighbors_key="vel_neighbors", scvelo_module=scvelo
    )
    steady = run_velocity_estimate(
        moments, mode="deterministic", min_r2=0.0, scvelo_module=scvelo
    )
    recovered = run_velocity_recover(
        steady,
        n_top_genes=15,
        max_iter=3,
        n_jobs=1,
        scvelo_module=scvelo,
    )
    dynamical = run_velocity_estimate(
        recovered,
        mode="dynamical",
        min_likelihood=0.0,
        scvelo_module=scvelo,
    )
    graph = run_velocity_graph(dynamical, scvelo_module=scvelo)
    table, ranking_summary = run_velocity_ranking(recovered, top_n=10)
    png, plot_summary = run_velocity_stream(graph, scvelo_module=scvelo)
    payload, summary, _, _ = validate_velocity_state(graph, allowed_stages=("velocity_graph",))
    assert payload.layers["Ms"].dtype == np.float32
    assert payload.layers["Mu"].dtype == np.float32
    assert summary["software_versions"]["scvelo"] == "0.3.4"
    assert recovered.summary["key_results"]["successful_fits"] >= 10
    assert dynamical.summary["key_results"]["selected_velocity_genes"] >= 10
    assert not table.empty and ranking_summary["key_results"]["successful_fits"] >= 10
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert plot_summary["key_results"]["png_bytes"] == len(png)
    assert recover_globals["make_unique_list"] is original_unique
    assert recover_globals["_read_pars"] is original_read
    assert divergence_module.compute_divergence is original_divergence
