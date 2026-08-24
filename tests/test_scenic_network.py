from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from openbio_singlecell import nodes_regulatory
from openbio_singlecell.node_types import AnnDataType, ScenicNetworkType
from openbio_singlecell.nodes_regulatory import OpenBioSingleCellRunPySCENIC, OpenBioSingleCellSCENICTFModules
from openbio_singlecell.scenic_network import ScenicNetwork


def output_values(node_output):
    return node_output.result


@pytest.fixture
def adata(science):
    obs = science.pd.DataFrame(index=["cell_1", "cell_2"])
    var = science.pd.DataFrame(index=["TF1", "G1", "G2"])
    return science.ad.AnnData(science.np.asarray([[2.0, 1.0, 5.0], [3.0, 4.0, 6.0]]), obs=obs, var=var)


def adjacency_frame(science):
    return science.pd.DataFrame(
        {
            "TF": ["TF1", "TF1"],
            "target": ["G1", "TF1"],
            "importance": [0.8, 0.2],
            "rho": [0.5, -0.1],
        }
    )


def test_scenic_node_schemas_use_the_network_contract():
    run_schema = OpenBioSingleCellRunPySCENIC.GET_SCHEMA()
    module_schema = OpenBioSingleCellSCENICTFModules.GET_SCHEMA()
    run_inputs = {input_.id: input_ for input_ in run_schema.inputs}

    assert [output.display_name for output in run_schema.outputs] == ["adata", "network"]
    assert [output.io_type for output in run_schema.outputs] == [AnnDataType.io_type, ScenicNetworkType.io_type]
    assert run_inputs["use_highly_variable"].default is False
    assert not run_inputs["use_highly_variable"].advanced
    assert run_inputs["highly_variable_key"].default == "highly_variable"
    assert run_inputs["highly_variable_key"].advanced is True
    assert run_inputs["grn_method"].default == "grnboost2"
    assert run_inputs["grn_method"].advanced is True
    assert run_inputs["mask_dropouts"].default is True
    assert run_inputs["mask_dropouts"].advanced is True
    assert run_inputs["auc_threshold"].default == 0.05
    assert not run_inputs["auc_threshold"].advanced
    assert run_inputs["source"].get_io_type() == "COMFY_DYNAMICCOMBO_V3"
    assert [(option.key, [item.id for item in option.inputs]) for option in run_inputs["source"].options] == [
        ("X", []),
        ("raw", []),
        ("layer", ["layer_name"]),
    ]
    assert [input_.id for input_ in module_schema.inputs] == [
        "adata",
        "network",
        "transcription_factor",
        "source",
    ]
    assert module_schema.inputs[1].io_type == ScenicNetworkType.io_type
    assert "validate_inputs" in OpenBioSingleCellRunPySCENIC.__dict__
    assert "fingerprint_inputs" in OpenBioSingleCellRunPySCENIC.__dict__
    assert "validate_inputs" not in OpenBioSingleCellSCENICTFModules.__dict__
    assert "fingerprint_inputs" not in OpenBioSingleCellSCENICTFModules.__dict__


def test_scenic_network_validates_and_owns_its_adjacency(science):
    adjacency = adjacency_frame(science)
    provenance = {"operation": "test", "resource": {"path": "ranking.feather"}}
    network = ScenicNetwork(
        adjacency=adjacency,
        gene_names=("TF1", "G1"),
        grn_method=" grnboost2 ",
        provenance=provenance,
    )

    adjacency.loc[0, "TF"] = "changed"
    provenance["resource"]["path"] = "changed"
    assert network.adjacency.loc[0, "TF"] == "TF1"
    assert network.gene_names == ("TF1", "G1")
    assert network.grn_method == "grnboost2"
    assert network.provenance["resource"]["path"] == "ranking.feather"
    assert list(network.adjacency.columns) == ["TF", "target", "importance", "rho"]

    with pytest.raises(TypeError, match="pandas DataFrame"):
        ScenicNetwork(
            adjacency=object(),
            gene_names=("TF1", "G1"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="cannot be empty"):
        ScenicNetwork(
            adjacency=adjacency_frame(science).iloc[0:0],
            gene_names=("TF1", "G1"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="missing required columns"):
        ScenicNetwork(
            adjacency=adjacency_frame(science).drop(columns="importance"),
            gene_names=("TF1", "G1"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="finite numbers"):
        invalid_importance = adjacency_frame(science)
        invalid_importance.loc[0, "importance"] = float("nan")
        ScenicNetwork(
            adjacency=invalid_importance,
            gene_names=("TF1", "G1"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="must be unique"):
        ScenicNetwork(
            adjacency=adjacency_frame(science),
            gene_names=("TF1", "G1", "G1"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="outside the inference gene set"):
        ScenicNetwork(
            adjacency=adjacency_frame(science),
            gene_names=("TF1", "G2"),
            grn_method="grnboost2",
            provenance={"operation": "test"},
        )
    with pytest.raises(ValueError, match="provenance cannot be empty"):
        ScenicNetwork(
            adjacency=adjacency_frame(science),
            gene_names=("TF1", "G1"),
            grn_method="grnboost2",
            provenance={},
        )


def test_run_pyscenic_materializes_network_before_temporary_directory_is_removed(
    adata,
    comfy_directories,
    science,
    monkeypatch,
):
    adata.var["highly_variable"] = [False, True, False]
    input_dir, _, _ = comfy_directories
    (input_dir / "tfs.txt").write_text("TF1\n", encoding="utf-8")
    (input_dir / "ranking.feather").write_bytes(b"ranking")
    (input_dir / "motifs.tbl").write_text("motifs\n", encoding="utf-8")
    run_parameters = {
        "tf_list_file": "tfs.txt",
        "ranking_database_files": "ranking.feather",
        "motif_annotations_file": "motifs.tbl",
        "grn_method": "grnboost2",
        "num_workers": 1,
        "random_seed": 7,
    }

    missing_key = adata.copy()
    del missing_key.var["highly_variable"]
    with pytest.raises(ValueError, match="Highly-variable gene key not found"):
        OpenBioSingleCellRunPySCENIC.execute(missing_key, use_highly_variable=True, **run_parameters)

    no_available_genes = adata[:, ["G1", "G2"]].copy()
    no_available_genes.var["highly_variable"] = False
    with pytest.raises(ValueError, match="No highly-variable or TF-list genes"):
        OpenBioSingleCellRunPySCENIC.execute(no_available_genes, use_highly_variable=True, **run_parameters)

    loom_input = {}

    class FakeLoompy:
        @staticmethod
        def create(path, matrix, row_attributes, column_attributes):
            loom_input["genes"] = tuple(row_attributes["Gene"])

    def require_optional_dependency(name):
        if name == "loompy":
            return FakeLoompy
        if name == "pyscenic.cli.pyscenic":
            return object()
        raise AssertionError(f"unexpected optional dependency: {name}")

    cli_calls = []

    def run_cli(arguments, working_directory):
        cli_calls.append((arguments, Path(working_directory)))
        if arguments[0] == "grn":
            adjacency_path = Path(arguments[arguments.index("-o") + 1])
            adjacency_frame(science).to_csv(adjacency_path, index=False)

    monkeypatch.setattr(nodes_regulatory, "_require_optional_dependency", require_optional_dependency)
    monkeypatch.setattr(nodes_regulatory, "_run_pyscenic_cli", run_cli)
    monkeypatch.setattr(nodes_regulatory, "_attach_pyscenic_results", lambda *args: None)

    output, network = output_values(
        OpenBioSingleCellRunPySCENIC.execute(
            adata,
            use_highly_variable=True,
            **run_parameters,
        )
    )

    assert output is not adata
    science.pd.testing.assert_frame_equal(network.adjacency, adjacency_frame(science))
    assert loom_input["genes"] == ("TF1", "G1")
    assert network.gene_names == ("TF1", "G1")
    assert network.grn_method == "grnboost2"
    assert network.provenance["expression"]["use_highly_variable"] is True
    assert network.provenance["tf_list"]["path"] == "tfs.txt"
    assert network.provenance["ranking_databases"][0]["path"] == "ranking.feather"
    assert [arguments[0] for arguments, _ in cli_calls] == ["grn", "ctx", "aucell"]
    assert all(not working_directory.exists() for _, working_directory in cli_calls)


def test_scenic_tf_modules_consumes_network_without_reading_a_file(adata, science, monkeypatch):
    network = ScenicNetwork(
        adjacency=adjacency_frame(science),
        gene_names=("G1", "TF1"),
        grn_method="grnboost2",
        provenance={"operation": "test"},
    )
    received = {}

    def modules_from_adjacencies(adjacency, expression):
        received["adjacency"] = adjacency
        received["expression"] = expression
        return [SimpleNamespace(transcription_factor="TF1", genes=["G1"])]

    def fail_file_access(*args, **kwargs):
        raise AssertionError("SCENIC TF Modules must not access adjacency files")

    monkeypatch.setattr(
        nodes_regulatory,
        "_require_optional_dependency",
        lambda name: SimpleNamespace(modules_from_adjacencies=modules_from_adjacencies),
    )
    monkeypatch.setattr(nodes_regulatory, "resolve_input_path", fail_file_access)
    monkeypatch.setattr(nodes_regulatory, "input_file_fingerprint", fail_file_access)
    monkeypatch.setattr(science.pd, "read_csv", fail_file_access)

    with pytest.raises(ValueError, match="missing SCENIC network genes.*G1"):
        OpenBioSingleCellSCENICTFModules.execute(
            adata[:, ["TF1", "G2"]].copy(), network, "TF1", {"source": "X"}
        )

    (result,) = output_values(
        OpenBioSingleCellSCENICTFModules.execute(adata, network, "TF1", {"source": "X"})
    )

    science.pd.testing.assert_frame_equal(received["adjacency"], network.adjacency)
    assert list(received["expression"].columns) == ["G1", "TF1"]
    assert result.parameters == {
        "grn_method": "grnboost2",
        "transcription_factor": "TF1",
        "source": "X",
    }
    assert result.table.to_dict("records") == [{"transcription_factor": "TF1", "module": 0, "gene": "G1"}]
