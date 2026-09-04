from __future__ import annotations

import json
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.artifact_codecs import read_plot, write_anndata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.node_types import AnnDataType, PlotResultType
from openbio_singlecell.nodes_trajectory import (
    OpenBioSingleCellDiffusionSpectrumPlot,
    OpenBioSingleCellDPTGeneTrendPlot,
    OpenBioSingleCellPAGAPlot,
)
from openbio_singlecell.operations_trajectory import (
    diffusion_map_owned,
    diffusion_spectrum_plot,
    diffusion_spectrum_plot_owned,
    dpt_gene_trend_plot,
    dpt_gene_trend_plot_owned,
    dpt_owned,
    paga_owned,
    paga_plot,
    paga_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _trajectory_input(science):
    rng = science.np.random.default_rng(73)
    adata = science.ad.AnnData(rng.normal(size=(36, 7)))
    adata.obs_names = [f"cell_{index:03d}" for index in range(36)]
    adata.var_names = [f"gene_{index}" for index in range(7)]
    adata.obs["state"] = science.pd.Categorical(
        ["early"] * 12 + ["middle"] * 12 + ["late"] * 12,
        categories=["early", "middle", "late"],
        ordered=True,
    )
    science.sc.pp.neighbors(
        adata,
        n_neighbors=6,
        use_rep="X",
        metric="euclidean",
        method="umap",
        random_state=11,
    )
    return adata


def _context(root: Path) -> OperationContext:
    root.mkdir()
    return OperationContext.from_request_path(root / "request.json", str(uuid.uuid4()))


def _descriptor(root: Path) -> dict[str, object]:
    return {
        "type": "artifact",
        "path": str(root.resolve()),
        "kind": "OPENBIO_ANNDATA",
        "codec": "anndata-h5ad-v1",
    }


def _snapshot(root: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in root.iterdir() if path.is_file()}


def test_trajectory_plot_schemas_are_domain_specific_read_only_companions():
    paga = OpenBioSingleCellPAGAPlot.define_schema()
    diffusion = OpenBioSingleCellDiffusionSpectrumPlot.define_schema()
    dpt = OpenBioSingleCellDPTGeneTrendPlot.define_schema()

    assert paga.node_id == "OpenBioSingleCellPAGAPlot"
    assert paga.display_name == "PAGA Plot"
    assert paga.category == "openbio/single-cell/trajectory"
    assert [item.id for item in paga.inputs] == ["adata", "min_connectivity"]
    assert paga.inputs[0].io_type == AnnDataType.io_type
    assert [item.io_type for item in paga.outputs] == [
        PlotResultType.io_type,
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]

    assert diffusion.node_id == "OpenBioSingleCellDiffusionSpectrumPlot"
    assert diffusion.display_name == "Diffusion Spectrum Plot"
    assert diffusion.category == "openbio/single-cell/trajectory"
    assert [item.id for item in diffusion.inputs] == ["adata"]
    assert [item.io_type for item in diffusion.outputs] == [
        PlotResultType.io_type,
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]

    assert dpt.node_id == "OpenBioSingleCellDPTGeneTrendPlot"
    assert dpt.display_name == "DPT Gene Trend Plot"
    assert dpt.category == "openbio/single-cell/trajectory"
    assert [item.id for item in dpt.inputs] == ["adata", "genes", "source", "n_bins"]
    assert [option.key for option in dpt.inputs[2].options] == ["layer", "X", "raw"]
    assert [item.io_type for item in dpt.outputs] == [
        PlotResultType.io_type,
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]


def test_paga_plot_renders_verified_stored_abstraction_without_direction_claim(science):
    output, _report, _code = paga_owned(_trajectory_input(science), groupby="state")
    before_connectivities = output.uns["paga"]["connectivities"].copy()
    before_tree = output.uns["paga"]["connectivities_tree"].copy()

    plotted, report, code = paga_plot_owned(output, min_connectivity=0.0)

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["groupby"] == "state"
    assert details["groups"] == ["early", "middle", "late"]
    assert details["group_sizes"] == [12, 12, 12]
    assert details["available_edges"] >= details["plotted_edges"]
    assert details["graph_fingerprint_sha256"] == output.uns["openbio_paga"][
        "graph_fingerprint_sha256"
    ]
    assert details["edge_style_semantics"] == {
        "solid_blue": "stored spanning-forest edge",
        "dashed_gray": "other retained stored connectivity edge",
    }
    assert "direction" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<paga-plot-code>", "exec"), namespace)
    assert namespace["plot_paga"](output) == plotted.png
    difference = output.uns["paga"]["connectivities"] - before_connectivities
    difference.eliminate_zeros()
    assert difference.nnz == 0
    difference = output.uns["paga"]["connectivities_tree"] - before_tree
    difference.eliminate_zeros()
    assert difference.nnz == 0


def test_diffusion_spectrum_plot_uses_graph_bound_stored_eigenvalues(science):
    output, _report, _code = diffusion_map_owned(_trajectory_input(science), n_comps=6, random_seed=17)
    coordinates_before = output.obsm["X_diffmap"].copy()
    eigenvalues_before = output.uns["diffmap_evals"].copy()

    plotted, report, code = diffusion_spectrum_plot_owned(output)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["component_indices"] == [0, 1, 2, 3, 4, 5]
    science.np.testing.assert_allclose(details["eigenvalues"], eigenvalues_before)
    assert details["stationary_component_index"] == 0
    assert details["informative_components"] == 5
    assert details["color_semantics"] == {
        "orange": "stationary component index 0",
        "blue": "informative diffusion components",
    }
    assert details["graph_fingerprint_sha256"] == output.uns["openbio_diffusion_map"][
        "graph_fingerprint_sha256"
    ]
    assert "sign" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<diffusion-spectrum-plot-code>", "exec"), namespace)
    assert namespace["plot_diffusion_spectrum"](output) == plotted.png
    science.np.testing.assert_array_equal(output.obsm["X_diffmap"], coordinates_before)
    science.np.testing.assert_array_equal(output.uns["diffmap_evals"], eigenvalues_before)


def test_dpt_gene_trend_plot_bins_explicit_expression_along_verified_pseudotime(science):
    diffmapped, _report, _code = diffusion_map_owned(_trajectory_input(science), n_comps=6, random_seed=17)
    output, _report, _code = dpt_owned(
        diffmapped,
        root_mode={"root_mode": "cell_id", "root_cell_id": "cell_000"},
        n_dcs=5,
    )
    expression_before = output.X.copy()
    pseudotime_before = output.obs["dpt_pseudotime"].copy()

    plotted, report, code = dpt_gene_trend_plot_owned(
        output,
        genes="gene_0,gene_1",
        source={"source": "X"},
        n_bins=10,
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["genes"] == ["gene_0", "gene_1"]
    assert details["source"] == {"source": "X"}
    assert details["requested_bins"] == 10
    assert 3 <= details["nonempty_bins"] <= 10
    assert sum(details["bin_cell_counts"]) == output.n_obs
    assert details["pseudotime_fingerprint_sha256"] == output.uns["openbio_dpt"][
        "output_fingerprint_sha256"
    ]
    assert "causal" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<dpt-gene-trend-plot-code>", "exec"), namespace)
    assert namespace["plot_dpt_gene_trends"](output) == plotted.png
    science.np.testing.assert_array_equal(output.X, expression_before)
    science.pd.testing.assert_series_equal(output.obs["dpt_pseudotime"], pseudotime_before)


def test_trajectory_plot_workers_read_ann_data_without_rewriting_inputs(science, tmp_path):
    paga_adata, _report, _code = paga_owned(_trajectory_input(science), groupby="state")
    paga_root = tmp_path / "paga-input"
    paga_root.mkdir()
    write_anndata(paga_root, paga_adata)
    paga_before = _snapshot(paga_root)
    paga_context = _context(tmp_path / "paga-plot")
    paga_records = paga_plot(
        paga_context,
        {"adata": _descriptor(paga_root)},
        {"min_connectivity": 0.0},
    )
    assert [record["name"] for record in paga_records] == ["plot", "summary", "code"]
    assert read_plot(paga_context.output_root / paga_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(paga_root) == paga_before

    diffmapped, _report, _code = diffusion_map_owned(_trajectory_input(science), n_comps=6, random_seed=17)
    dpt_adata, _report, _code = dpt_owned(
        diffmapped,
        root_mode={"root_mode": "cell_id", "root_cell_id": "cell_000"},
        n_dcs=5,
    )
    dpt_root = tmp_path / "dpt-input"
    dpt_root.mkdir()
    write_anndata(dpt_root, dpt_adata)
    dpt_before = _snapshot(dpt_root)
    diffusion_context = _context(tmp_path / "diffusion-plot")
    diffusion_records = diffusion_spectrum_plot(
        diffusion_context,
        {"adata": _descriptor(dpt_root)},
        {},
    )
    trend_context = _context(tmp_path / "trend-plot")
    trend_records = dpt_gene_trend_plot(
        trend_context,
        {"adata": _descriptor(dpt_root)},
        {"genes": "gene_0", "source": {"source": "X"}, "n_bins": 10},
    )
    assert [record["name"] for record in diffusion_records] == ["plot", "summary", "code"]
    assert [record["name"] for record in trend_records] == ["plot", "summary", "code"]
    assert read_plot(diffusion_context.output_root / diffusion_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert read_plot(trend_context.output_root / trend_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(dpt_root) == dpt_before


def test_trajectory_plots_reject_tampered_graph_diffusion_and_pseudotime_evidence(science):
    paga_adata, _report, _code = paga_owned(_trajectory_input(science), groupby="state")
    paga_adata.uns["paga"]["connectivities"].data[0] += 0.1
    with pytest.raises((ValueError, RuntimeError), match="symmetric|fingerprint|weights"):
        paga_plot_owned(paga_adata)

    diffmapped, _report, _code = diffusion_map_owned(_trajectory_input(science), n_comps=6, random_seed=17)
    diffmapped.uns["diffmap_evals"][1] += 0.05
    with pytest.raises((ValueError, RuntimeError), match="fingerprint|eigenpair|non-increasing"):
        diffusion_spectrum_plot_owned(diffmapped)

    diffmapped, _report, _code = diffusion_map_owned(_trajectory_input(science), n_comps=6, random_seed=17)
    dpt_adata, _report, _code = dpt_owned(
        diffmapped,
        root_mode={"root_mode": "cell_id", "root_cell_id": "cell_000"},
        n_dcs=5,
    )
    dpt_adata.obs.loc[dpt_adata.obs_names[-1], "dpt_pseudotime"] *= 0.5
    with pytest.raises(ValueError, match="fingerprint"):
        dpt_gene_trend_plot_owned(dpt_adata, genes="gene_0", source={"source": "X"}, n_bins=10)
