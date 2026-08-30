from __future__ import annotations

import asyncio
import gc
import json
import sys
import tempfile
import uuid
from pathlib import Path
from types import SimpleNamespace

import pytest

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION
from openbio_singlecell.artifact_codecs import read_anndata, write_anndata
from openbio_singlecell.artifact_envelope import summary_from_metadata
from openbio_singlecell.artifact_runtime import ArtifactTicket
from openbio_singlecell.artifact_service import (
    current_artifact_runtime,
    execute_artifact_node,
    initialize_artifact_service,
)
from openbio_singlecell.contracts import ensure_metadata, record_history
from openbio_singlecell.files import (
    input_file_provenance,
    resolve_input_path,
    tenx_mtx_provenance,
)
from openbio_singlecell.nodes_input import (
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoad10xStudy,
    OpenBioSingleCellLoadH5AD,
)
from openbio_singlecell.operations_input import (
    anndata_summary,
    load_10x_h5,
    load_10x_mtx,
    load_10x_study,
    load_h5ad,
)
from openbio_singlecell.payload import result_to_payload
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError


def output_value(node_output):
    return node_output.result[0]


def _json_default(value):
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"Unsupported JSON test value: {type(value).__name__}")


def _run_input_operation(operation, inputs, parameters):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        context = OperationContext(root, str(uuid.uuid4()))
        records = operation(context, inputs, parameters)
        values = []
        for record in records:
            if record["type"] == "artifact":
                values.append(read_anndata(root / record["payload"]))
            elif record["type"] == "summary":
                values.append(summary_from_metadata(record["value"]))
            else:
                values.append(record["value"])
        return SimpleNamespace(result=tuple(values))


def _file_descriptor(path, extensions):
    resolved = Path(resolve_input_path(path, extensions=extensions))
    return {
        "type": "file",
        "path": str(resolved),
        "provenance": input_file_provenance(path, extensions),
    }


def _execute_load_h5ad(path, make_var_names_unique=True):
    return _run_input_operation(
        load_h5ad,
        {"path": _file_descriptor(path, (".h5ad",))},
        {"make_var_names_unique": make_var_names_unique},
    )


def _execute_load_10x_mtx(
    directory,
    var_names="gene_symbols",
    make_unique=True,
    gex_only=True,
):
    resolved = Path(resolve_input_path(directory, kind="directory"))
    return _run_input_operation(
        load_10x_mtx,
        {
            "directory": {
                "type": "directory",
                "path": str(resolved),
                "provenance": tenx_mtx_provenance(directory),
            }
        },
        {"var_names": var_names, "make_unique": make_unique, "gex_only": gex_only},
    )


def _execute_load_10x_study(
    directory,
    sample_key="sample",
    var_names="gene_symbols",
    join="inner",
    make_unique=True,
    gex_only=True,
):
    resolved = Path(resolve_input_path(directory, kind="directory"))
    portable_root = Path(directory).as_posix().rstrip("/")
    samples = {
        child.name: tenx_mtx_provenance(f"{portable_root}/{child.name}")
        for child in sorted(resolved.iterdir(), key=lambda item: item.name.casefold())
        if child.is_dir() or child.is_symlink()
    }
    return _run_input_operation(
        load_10x_study,
        {
            "directory": {
                "type": "directory",
                "path": str(resolved),
                "provenance": {"path": portable_root, "samples": samples},
            }
        },
        {
            "sample_key": sample_key,
            "var_names": var_names,
            "join": join,
            "make_unique": make_unique,
            "gex_only": gex_only,
        },
    )


def _execute_load_10x_h5(path, genome="", gex_only=True, make_unique=True):
    return _run_input_operation(
        load_10x_h5,
        {"path": _file_descriptor(path, (".h5", ".hdf5"))},
        {"genome": genome, "gex_only": gex_only, "make_unique": make_unique},
    )


def _execute_anndata_summary(adata):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        artifact_root = root / "adata"
        artifact_root.mkdir()
        write_anndata(artifact_root, adata)
        return _run_input_operation(
            anndata_summary,
            {
                "adata": {
                    "type": "artifact",
                    "path": str(artifact_root.resolve()),
                    "kind": "OPENBIO_ANNDATA",
                    "codec": "anndata-h5ad-v1",
                }
            },
            {},
        )


def write_10x_mtx(directory, science, matrix, barcodes, features):
    directory.mkdir(parents=True)
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.mmwrite(directory / "matrix.mtx", science.sparse.csr_matrix(matrix))
    (directory / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")
    (directory / "features.tsv").write_text(
        "".join("\t".join(feature) + "\n" for feature in features),
        encoding="utf-8",
    )


def write_10x_h5(path, science, matrix, barcodes, gene_ids, gene_names, feature_types=None):
    h5py = pytest.importorskip("h5py")
    feature_types = feature_types or ["Gene Expression"] * len(gene_ids)
    stored = science.sparse.csc_matrix(matrix)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=stored.data)
        group.create_dataset("indices", data=stored.indices)
        group.create_dataset("indptr", data=stored.indptr)
        group.create_dataset("shape", data=stored.shape)
        group.create_dataset("barcodes", data=science.np.asarray(barcodes, dtype="S"))
        features = group.create_group("features")
        features.create_dataset("id", data=science.np.asarray(gene_ids, dtype="S"))
        features.create_dataset("name", data=science.np.asarray(gene_names, dtype="S"))
        features.create_dataset("feature_type", data=science.np.asarray(feature_types, dtype="S"))
        features.create_dataset("genome", data=science.np.asarray(["test"] * len(gene_ids), dtype="S"))
        features.create_dataset("_all_tag_keys", data=science.np.asarray(["genome"], dtype="S"))


def test_load_10x_mtx_operation_publishes_anndata_artifact(tmp_path, science):
    source = tmp_path / "tenx"
    write_10x_mtx(
        source,
        science,
        science.np.asarray([[1, 0], [0, 2]]),
        ["cell-a", "cell-b"],
        [("id-a", "A", "Gene Expression"), ("id-b", "B", "Gene Expression")],
    )
    output_root = tmp_path / "run"
    output_root.mkdir()
    records = load_10x_mtx(
        OperationContext(output_root, "00000000-0000-0000-0000-000000000001"),
        {
            "directory": {
                "type": "directory",
                "path": str(source.resolve()),
                "provenance": {"path": "tenx", "files": {}},
            }
        },
        {"var_names": "gene_symbols", "make_unique": True, "gex_only": True},
    )

    assert records == [
        {
            "type": "artifact",
            "name": "adata",
            "kind": "OPENBIO_ANNDATA",
            "codec": "anndata-h5ad-v1",
            "payload": "outputs/adata",
        }
    ]
    loaded = read_anndata(output_root / "outputs" / "adata")
    assert loaded.shape == (2, 2)
    assert loaded.uns["openbio_singlecell"]["source"]["path"] == "tenx"


def test_load_10x_study_operation_discovers_sorted_samples_and_preserves_provenance(tmp_path, science):
    study = tmp_path / "study"
    for sample_name, value in (("sample-b", 2), ("sample-a", 1)):
        write_10x_mtx(
            study / sample_name,
            science,
            science.np.asarray([[value]]),
            ["cell"],
            [("id-a", "A", "Gene Expression")],
        )
    output_root = tmp_path / "run"
    output_root.mkdir()
    sample_provenance = {name: {"path": f"study/{name}", "files": {}} for name in ("sample-a", "sample-b")}
    records = load_10x_study(
        OperationContext(output_root, "00000000-0000-0000-0000-000000000002"),
        {
            "directory": {
                "type": "directory",
                "path": str(study.resolve()),
                "provenance": {"path": "study", "samples": sample_provenance},
            }
        },
        {
            "sample_key": "sample",
            "var_names": "gene_symbols",
            "join": "inner",
            "make_unique": True,
            "gex_only": True,
        },
    )

    loaded = read_anndata(output_root / records[0]["payload"])
    assert loaded.obs["sample"].astype(str).tolist() == ["sample-a", "sample-b"]
    source = loaded.uns["openbio_singlecell"]["source"]
    assert list(source["samples"]) == ["sample-a", "sample-b"]
    assert source["samples"]["sample-a"]["path"] == "study/sample-a"


def test_load_10x_h5_operation_filters_to_gene_expression_and_writes_new_artifact(tmp_path, science):
    source = tmp_path / "matrix.h5"
    write_10x_h5(
        source,
        science,
        science.np.asarray([[1, 0], [0, 2]]),
        ["cell-a", "cell-b"],
        ["id-a", "id-b"],
        ["A", "B"],
        ["Gene Expression", "Antibody Capture"],
    )
    output_root = tmp_path / "run"
    output_root.mkdir()
    records = load_10x_h5(
        OperationContext(output_root, "00000000-0000-0000-0000-000000000003"),
        {
            "path": {
                "type": "file",
                "path": str(source.resolve()),
                "provenance": {"path": "matrix.h5", "size": source.stat().st_size, "mtime_ns": 1},
            }
        },
        {"genome": "", "gex_only": True, "make_unique": True},
    )

    loaded = read_anndata(output_root / records[0]["payload"])
    assert loaded.shape == (2, 1)
    assert loaded.var["feature_types"].tolist() == ["Gene Expression"]
    assert loaded.uns["openbio_singlecell"]["source"]["path"] == "matrix.h5"


def test_anndata_summary_operation_returns_only_small_json_and_code(tmp_path, adata):
    source_root = tmp_path / "input-artifact"
    source_root.mkdir()
    write_anndata(source_root, adata)
    before = (source_root / "data.h5ad").read_bytes()
    output_root = tmp_path / "run"
    output_root.mkdir()

    records = anndata_summary(
        OperationContext(output_root, "00000000-0000-0000-0000-000000000004"),
        {
            "adata": {
                "type": "artifact",
                "path": str(source_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            }
        },
        {},
    )

    assert [record["type"] for record in records] == ["summary", "string"]
    assert [record["name"] for record in records] == ["summary", "code"]
    assert records[0]["value"]["summary"]["node_id"] == "OpenBioSingleCellAnnDataSummary"
    compile(records[1]["value"], "<anndata-summary>", "exec")
    assert not (output_root / "outputs").exists()
    assert (source_root / "data.h5ad").read_bytes() == before


def test_every_input_operation_rejects_extended_parameter_objects(tmp_path, adata, science):
    h5ad_path = tmp_path / "input.h5ad"
    adata.write_h5ad(h5ad_path)
    mtx = tmp_path / "tenx"
    write_10x_mtx(
        mtx,
        science,
        [[1]],
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    study = tmp_path / "study"
    write_10x_mtx(
        study / "sample",
        science,
        [[1]],
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    h5_path = tmp_path / "matrix.h5"
    write_10x_h5(h5_path, science, [[1]], ["cell"], ["id"], ["gene"])
    artifact = tmp_path / "artifact"
    artifact.mkdir()
    write_anndata(artifact, adata)
    cases = [
        (
            load_h5ad,
            {
                "path": {
                    "type": "file",
                    "path": str(h5ad_path.resolve()),
                    "provenance": {"path": "input.h5ad"},
                }
            },
            {"make_var_names_unique": True},
        ),
        (
            load_10x_mtx,
            {
                "directory": {
                    "type": "directory",
                    "path": str(mtx.resolve()),
                    "provenance": {"path": "tenx", "files": {}},
                }
            },
            {"var_names": "gene_symbols", "make_unique": True, "gex_only": True},
        ),
        (
            load_10x_study,
            {
                "directory": {
                    "type": "directory",
                    "path": str(study.resolve()),
                    "provenance": {
                        "path": "study",
                        "samples": {"sample": {"path": "study/sample", "files": {}}},
                    },
                }
            },
            {
                "sample_key": "sample",
                "var_names": "gene_symbols",
                "join": "inner",
                "make_unique": True,
                "gex_only": True,
            },
        ),
        (
            load_10x_h5,
            {
                "path": {
                    "type": "file",
                    "path": str(h5_path.resolve()),
                    "provenance": {"path": "matrix.h5"},
                }
            },
            {"genome": "", "gex_only": True, "make_unique": True},
        ),
        (
            anndata_summary,
            {
                "adata": {
                    "type": "artifact",
                    "path": str(artifact.resolve()),
                    "kind": "OPENBIO_ANNDATA",
                    "codec": "anndata-h5ad-v1",
                }
            },
            {},
        ),
    ]

    for index, (operation, inputs, parameters) in enumerate(cases):
        output_root = tmp_path / f"strict-{index}"
        output_root.mkdir()
        with pytest.raises(ProtocolError, match="parameters must be exactly"):
            operation(
                OperationContext(output_root, str(uuid.uuid4())),
                inputs,
                {**parameters, "unexpected": True},
            )


def test_input_operation_module_is_worker_only():
    source = (Path(__file__).parents[1] / "openbio_singlecell" / "operations_input.py").read_text(encoding="utf-8")
    assert "comfy_api" not in source
    assert "folder_paths" not in source
    assert "nodes_" not in source


def test_10x_directory_loaders_run_through_artifact_service_and_real_workers(
    comfy_directories,
    science,
):
    input_dir, _, temp_dir = comfy_directories
    mtx_files = input_dir / "service-tenx"
    write_10x_mtx(
        mtx_files,
        science,
        [[1, 0], [0, 2]],
        ["cell-a", "cell-b"],
        [("id-a", "A", "Gene Expression"), ("id-b", "B", "Gene Expression")],
    )
    study_files = input_dir / "service-study"
    for sample_name, value in (("sample-b", 2), ("sample-a", 1)):
        write_10x_mtx(
            study_files / sample_name,
            science,
            [[value]],
            ["cell"],
            [("id-a", "A", "Gene Expression")],
        )
    source_bytes = {
        path.relative_to(input_dir).as_posix(): path.read_bytes()
        for path in [*mtx_files.iterdir(), *study_files.rglob("*")]
        if path.is_file()
    }

    async def execute_loaders():
        gc.collect()
        await initialize_artifact_service(temp_dir, sys.executable)
        mtx_output = await execute_artifact_node(
            OpenBioSingleCellLoad10xMTX,
            None,
            {
                "directory": "service-tenx",
                "var_names": "gene_symbols",
                "make_unique": True,
                "gex_only": True,
            },
        )
        study_output = await execute_artifact_node(
            OpenBioSingleCellLoad10xStudy,
            None,
            {
                "directory": "service-study",
                "sample_key": "sample",
                "var_names": "gene_symbols",
                "join": "inner",
                "make_unique": True,
                "gex_only": True,
            },
        )
        return mtx_output[0], study_output[0]

    mtx_ticket, study_ticket = asyncio.run(execute_loaders())
    assert isinstance(mtx_ticket, ArtifactTicket)
    assert isinstance(study_ticket, ArtifactTicket)
    runtime = current_artifact_runtime()
    mtx = read_anndata(runtime.resolve(mtx_ticket))
    study = read_anndata(runtime.resolve(study_ticket))
    assert mtx.shape == (2, 2)
    assert study.obs["sample"].astype(str).tolist() == ["sample-a", "sample-b"]
    assert mtx.uns["openbio_singlecell"]["source"]["path"] == "service-tenx"
    assert study.uns["openbio_singlecell"]["source"]["path"] == "service-study"
    assert {
        path.relative_to(input_dir).as_posix(): path.read_bytes()
        for path in [*mtx_files.iterdir(), *study_files.rglob("*")]
        if path.is_file()
    } == source_bytes


def test_binary_file_loaders_use_the_web_upload_widget():
    expected = {
        OpenBioSingleCellLoadH5AD: [".h5ad"],
        OpenBioSingleCellLoad10xH5: [".h5", ".hdf5"],
    }

    for node_class, extensions in expected.items():
        path_input = next(input_ for input_ in node_class.define_schema().inputs if input_.id == "path")
        assert path_input.extra_dict["widgetType"] == "OPENBIO_INPUT_FILE_UPLOAD_WIDGET"
        assert path_input.extra_dict["allowed_extensions"] == extensions
        assert path_input.extra_dict["upload_subfolder"] == "openbio-singlecell"


@pytest.fixture
def adata(science):
    counts = science.sparse.csr_matrix(science.np.arange(24).reshape(6, 4))
    obs = science.pd.DataFrame(index=[f"cell_{index}" for index in range(6)])
    var = science.pd.DataFrame(index=["MT-A", "B", "C", "D"])
    value = science.ad.AnnData(counts, obs=obs, var=var)
    ensure_metadata(value, display_name="synthetic", source={"kind": "test"})
    return value


def test_h5ad_and_10x_mtx_loaders(comfy_directories, adata, science):
    input_dir, _, _ = comfy_directories
    h5ad_path = input_dir / "sample.h5ad"
    adata.write_h5ad(h5ad_path)

    loaded = output_value(_execute_load_h5ad("sample.h5ad", True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "sample"
    h5ad_source = loaded.uns["openbio_singlecell"]["source"]
    assert h5ad_source["kind"] == "h5ad"
    assert h5ad_source["path"] == "sample.h5ad"
    assert "fingerprint" not in h5ad_source
    assert str(input_dir.resolve()).lower() not in json.dumps(h5ad_source, default=_json_default).lower()

    tenx = input_dir / "tenx"
    write_10x_mtx(
        tenx,
        science,
        adata.X.transpose(),
        list(adata.obs_names),
        [(f"id_{name}", name, "Gene Expression") for name in adata.var_names],
    )

    loaded_10x = output_value(_execute_load_10x_mtx("tenx", "gene_symbols", True, True))

    assert loaded_10x.shape == adata.shape
    assert science.sparse.issparse(loaded_10x.X)
    tenx_source = loaded_10x.uns["openbio_singlecell"]["source"]
    assert loaded_10x.uns["openbio_singlecell"]["display_name"] == "tenx"
    assert tenx_source["kind"] == "10x_mtx"
    assert tenx_source["path"] == "tenx"
    assert tenx_source["files"]["matrix"]["path"] == "tenx/matrix.mtx"
    assert "fingerprint" not in tenx_source
    assert str(input_dir.resolve()).lower() not in json.dumps(tenx_source, default=_json_default).lower()


def test_h5ad_loader_boundaries(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    assert "not found" in OpenBioSingleCellLoadH5AD.validate_inputs("missing.h5ad", True).lower()

    (input_dir / "corrupt.h5ad").write_bytes(b"not an hdf5 file")
    with pytest.raises(OSError):
        _execute_load_h5ad("corrupt.h5ad", True)

    empty = science.ad.AnnData(science.np.empty((0, 2)))
    empty.write_h5ad(input_dir / "empty.h5ad")
    empty_loaded = output_value(_execute_load_h5ad("empty.h5ad", True))
    assert empty_loaded.shape == (0, 2)
    assert any(
        "empty observation axis" in warning for warning in empty_loaded.uns["openbio_singlecell"]["source"]["warnings"]
    )

    duplicate_names = science.ad.AnnData(science.np.eye(2))
    duplicate_names.var_names = ["gene", "gene"]
    duplicate_names.write_h5ad(input_dir / "duplicate.h5ad")
    loaded = output_value(_execute_load_h5ad("duplicate.h5ad", True))
    assert loaded.var_names.is_unique
    assert loaded.uns["openbio_singlecell"]["source"]["axis_names"]["var_names_repaired"] == 1

    preserved = output_value(_execute_load_h5ad("duplicate.h5ad", False))
    assert not preserved.var_names.is_unique
    assert preserved.uns["openbio_singlecell"]["source"]["warnings"]

    duplicate_cells = science.ad.AnnData(science.np.eye(2))
    duplicate_cells.obs_names = ["cell", "cell"]
    duplicate_cells.write_h5ad(input_dir / "duplicate_cells.h5ad")
    duplicate_cells_loaded = output_value(_execute_load_h5ad("duplicate_cells.h5ad", True))
    duplicate_source = duplicate_cells_loaded.uns["openbio_singlecell"]["source"]
    assert not duplicate_cells_loaded.obs_names.is_unique
    assert any("duplicate observation-name" in warning for warning in duplicate_source["warnings"])


def test_summary_normalizes_metadata(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    value = science.ad.AnnData(science.np.eye(2))
    value.uns["openbio_singlecell"] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "display_name": "sample",
        "warnings": 1,
        "analysis_history": {},
        "random_seed": "bad",
    }
    value.write_h5ad(input_dir / "metadata.h5ad")

    loaded = output_value(_execute_load_h5ad("metadata.h5ad", True))
    summary = output_value(_execute_anndata_summary(loaded))
    record_history(loaded, "test", {}, loaded.n_obs, loaded.n_vars)
    metadata = loaded.uns["openbio_singlecell"]

    assert summary.title == "metadata summary"
    assert summary.summary["key_results"]["shape"] == [2, 2]
    assert summary.summary["node_id"] == "OpenBioSingleCellAnnDataSummary"
    assert summary.warnings == []
    assert summary.random_seed == 0
    assert metadata["warnings"] == []
    assert metadata["analysis_history"]["000000"]["operation"] == "test"
    assert metadata["source"]["kind"] == "h5ad"
    assert metadata["source"]["path"] == "metadata.h5ad"


def test_h5ad_loader_redacts_absolute_paths_in_embedded_provenance(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    value = science.ad.AnnData(science.np.eye(2))
    ensure_metadata(
        value,
        source={
            "kind": "external",
            "windows_path": r"C:\private\donor.h5ad",
            "posix_path": "/private/donor.h5ad",
            "relative_path": "portable/donor.h5ad",
            "nested": {"file": "file:///private/donor.h5ad"},
        },
    )
    value.write_h5ad(input_dir / "embedded.h5ad")

    loaded = output_value(_execute_load_h5ad("embedded.h5ad"))
    embedded = loaded.uns["openbio_singlecell"]["source"]["embedded_source"]

    assert embedded["windows_path"] == "<redacted-absolute-path>"
    assert embedded["posix_path"] == "<redacted-absolute-path>"
    assert embedded["nested"]["file"] == "<redacted-absolute-path>"
    assert embedded["relative_path"] == "portable/donor.h5ad"

    array_value = science.ad.AnnData(science.np.eye(1))
    array_value.uns["openbio_singlecell"] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "display_name": "array source",
        "warnings": [],
        "analysis_history": {},
        "random_seed": 0,
        "source": {
            "paths": science.np.asarray(["/private/a.h5ad", r"C:\private\b.h5ad"]),
            "byte_paths": science.np.asarray([b"/private/c.h5ad", b"C:\\private\\d.h5ad"]),
        },
    }
    array_value.write_h5ad(input_dir / "embedded_array.h5ad")
    array_loaded = output_value(_execute_load_h5ad("embedded_array.h5ad"))
    assert array_loaded.uns["openbio_singlecell"]["source"]["embedded_source"]["paths"].tolist() == [
        "<redacted-absolute-path>",
        "<redacted-absolute-path>",
    ]
    assert array_loaded.uns["openbio_singlecell"]["source"]["embedded_source"]["byte_paths"].tolist() == [
        "<redacted-absolute-path>",
        "<redacted-absolute-path>",
    ]


def test_10x_mtx_loader_boundaries(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    missing = input_dir / "missing"
    missing.mkdir()
    (missing / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n0 0 0\n",
        encoding="utf-8",
    )
    validation = OpenBioSingleCellLoad10xMTX.validate_inputs("missing")
    assert "missing barcodes" in validation

    mismatched = input_dir / "mismatched"
    mismatched.mkdir()
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.mmwrite(mismatched / "matrix.mtx", science.sparse.csr_matrix([[1], [2]]))
    (mismatched / "barcodes.tsv").write_text("cell\n", encoding="utf-8")
    (mismatched / "features.tsv").write_text("id\tgene\tGene Expression\n", encoding="utf-8")
    with pytest.raises(ValueError, match="dimensions do not match"):
        _execute_load_10x_mtx("mismatched", "gene_symbols", True, True)


@pytest.mark.parametrize(
    ("name", "matrix", "barcodes", "features", "message"),
    [
        (
            "duplicate_barcodes",
            [[1, 2]],
            ["cell", "cell"],
            [("id", "gene", "Gene Expression")],
            "duplicate barcode",
        ),
        (
            "blank_barcode",
            [[1, 2]],
            ["cell", ""],
            [("id", "gene", "Gene Expression")],
            "blank barcode",
        ),
        (
            "blank_gene_id",
            [[1]],
            ["cell"],
            [("", "gene", "Gene Expression")],
            "blank stable feature ID",
        ),
        (
            "duplicate_gene_id",
            [[1], [2]],
            ["cell"],
            [("id", "A", "Gene Expression"), ("id", "B", "Gene Expression")],
            "duplicate stable feature ID",
        ),
        (
            "all_zero",
            [[0]],
            ["cell"],
            [("id", "gene", "Gene Expression")],
            "no positive expression",
        ),
    ],
)
def test_10x_mtx_discloses_identity_and_empty_count_advisories(
    comfy_directories, science, name, matrix, barcodes, features, message
):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(input_dir / name, science, matrix, barcodes, features)
    loaded = output_value(_execute_load_10x_mtx(name, "gene_symbols", True, True))
    source = loaded.uns["openbio_singlecell"]["source"]
    assert any(message in warning for warning in source["warnings"])


@pytest.mark.parametrize(
    ("name", "matrix", "message"),
    [
        ("negative", [[-1]], "negative expression"),
        ("fractional", [[1.5]], "non-integer"),
    ],
)
def test_10x_mtx_rejects_values_incompatible_with_count_format(comfy_directories, science, name, matrix, message):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / name,
        science,
        matrix,
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    with pytest.raises(ValueError, match=message):
        _execute_load_10x_mtx(name, "gene_symbols", True, True)


def test_10x_mtx_repairs_only_duplicate_feature_names(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / "duplicate_symbols",
        science,
        [[1], [2]],
        ["cell"],
        [("id_1", "gene", "Gene Expression"), ("id_2", "gene", "Gene Expression")],
    )

    repaired = output_value(_execute_load_10x_mtx("duplicate_symbols", "gene_symbols", True, True))
    preserved = output_value(_execute_load_10x_mtx("duplicate_symbols", "gene_symbols", False, True))

    assert repaired.obs_names.tolist() == ["cell"]
    assert repaired.var_names.tolist() == ["gene", "gene-1"]
    assert not preserved.var_names.is_unique
    assert repaired.uns["openbio_singlecell"]["source"]["axis_names"]["var_names_repaired"] == 1


def test_10x_mtx_rejects_ambiguous_file_roles(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / "ambiguous",
        science,
        [[1]],
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    (input_dir / "ambiguous" / "barcodes.tsv.gz").write_bytes(b"placeholder")
    message = OpenBioSingleCellLoad10xMTX.validate_inputs("ambiguous")
    assert "ambiguous barcodes" in message


def test_10x_study_never_silently_omits_invalid_samples(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1]],
        ["cell_a"],
        [("id", "gene", "Gene Expression")],
    )
    broken = root / "Sample_B"
    broken.mkdir()
    (broken / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n0 0 0\n")

    validation = OpenBioSingleCellLoad10xStudy.validate_inputs("study")
    assert "Sample_B" in validation
    assert "missing barcodes" in validation
    with pytest.raises(ValueError, match="Sample_B"):
        OpenBioSingleCellLoad10xStudy.fingerprint_inputs("study")
    with pytest.raises((ValueError, FileNotFoundError), match="Sample_B|missing barcodes"):
        _execute_load_10x_study("study")


def test_10x_study_reports_inner_and_outer_alignment(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1], [2]],
        ["barcode"],
        [("id_common", "common", "Gene Expression"), ("id_a", "only_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[3], [4]],
        ["barcode"],
        [("id_common", "common", "Gene Expression"), ("id_b", "only_b", "Gene Expression")],
    )

    inner = output_value(_execute_load_10x_study("study", join="inner"))
    outer = output_value(_execute_load_10x_study("study", join="outer"))

    assert inner.shape == (2, 1)
    assert outer.shape == (2, 3)
    assert science.np.asarray(outer.X.toarray()).tolist() == [[1, 2, 0], [3, 0, 4]]
    assert inner.obs["sample"].astype(str).tolist() == ["Sample_A", "Sample_B"]
    assert inner.obs_names.tolist() == ["barcode-Sample_A", "barcode-Sample_B"]
    inner_source = inner.uns["openbio_singlecell"]["source"]
    outer_source = outer.uns["openbio_singlecell"]["source"]
    assert inner_source["sample_count"] == 2
    assert set(inner_source["samples"]) == {"Sample_A", "Sample_B"}
    assert inner_source["feature_alignment"]["Sample_A"]["dropped_features"] == 1
    assert outer_source["feature_alignment"]["Sample_A"]["zero_filled_features"] == 1
    assert outer_source["warnings"]
    assert outer.var["gene_ids"].tolist() == ["id_common", "id_a", "id_b"]
    assert outer.var["gene_symbols"].tolist() == ["common", "only_a", "only_b"]
    assert outer.var["feature_types"].tolist() == ["Gene Expression"] * 3
    assert str(input_dir.resolve()).lower() not in json.dumps(outer_source, default=_json_default).lower()


def test_10x_study_returns_empty_inner_feature_intersection_with_advisory(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    root = input_dir / "disjoint"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1]],
        ["a"],
        [("id_a", "gene_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[2]],
        ["b"],
        [("id_b", "gene_b", "Gene Expression")],
    )

    output = output_value(_execute_load_10x_study("disjoint", join="inner"))
    source = output.uns["openbio_singlecell"]["source"]
    assert output.shape == (2, 0)
    assert source["retained_features"] == 0
    assert any("empty feature intersection" in warning for warning in source["warnings"])


def test_10x_study_rejects_conflicting_symbol_identity(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "A",
        science,
        [[1]],
        ["a"],
        [("id_1", "shared_name", "Gene Expression")],
    )
    write_10x_mtx(
        root / "B",
        science,
        [[1]],
        ["b"],
        [("id_2", "shared_name", "Gene Expression")],
    )
    with pytest.raises(ValueError, match="conflicting identity metadata"):
        _execute_load_10x_study("study")


def test_10x_study_aggregates_sample_payload_errors(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[-1]],
        ["a"],
        [("id_a", "gene_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[1, 2]],
        ["b", "b"],
        [("id_b", "gene_b", "Gene Expression")],
    )

    with pytest.raises(ValueError, match="Invalid 10x Study Sample payloads") as error:
        _execute_load_10x_study("study")

    message = str(error.value)
    assert "Sample_A" in message and "negative expression" in message
    assert "Sample_B" not in message


def test_10x_h5_loader(comfy_directories, adata, science):
    input_dir, _, _ = comfy_directories
    path = input_dir / "matrix.h5"
    write_10x_h5(
        path,
        science,
        adata.X.transpose(),
        list(adata.obs_names),
        [f"id_{value}" for value in adata.var_names],
        list(adata.var_names),
    )

    loaded = output_value(_execute_load_10x_h5("matrix.h5", "", True, True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "matrix"
    source = loaded.uns["openbio_singlecell"]["source"]
    assert source["kind"] == "10x_h5"
    assert source["path"] == "matrix.h5"
    assert "fingerprint" not in source
    assert str(input_dir.resolve()).lower() not in json.dumps(source, default=_json_default).lower()


def test_10x_h5_rejects_duplicate_barcodes_and_invalid_counts(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_h5(
        input_dir / "duplicate.h5",
        science,
        [[1, 2]],
        ["cell", "cell"],
        ["id"],
        ["gene"],
    )
    duplicate = output_value(_execute_load_10x_h5("duplicate.h5"))
    assert not duplicate.obs_names.is_unique
    assert any(
        "duplicate observation-name" in warning for warning in duplicate.uns["openbio_singlecell"]["source"]["warnings"]
    )

    write_10x_h5(
        input_dir / "fractional.h5",
        science,
        [[1.5]],
        ["cell"],
        ["id"],
        ["gene"],
    )
    with pytest.raises(ValueError, match="non-integer"):
        _execute_load_10x_h5("fractional.h5")


def test_10x_h5_discloses_feature_filtering(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_h5(
        input_dir / "multimodal.h5",
        science,
        [[1], [2]],
        ["cell"],
        ["gene_id", "antibody_id"],
        ["gene", "antibody"],
        ["Gene Expression", "Antibody Capture"],
    )

    loaded = output_value(_execute_load_10x_h5("multimodal.h5", gex_only=True))
    source = loaded.uns["openbio_singlecell"]["source"]

    assert loaded.var_names.tolist() == ["gene"]
    assert source["input_features"] == 2
    assert source["retained_features"] == 1
    assert source["filtered_features"] == 1
    assert source["feature_types_before_filter"] == {
        "Gene Expression": 1,
        "Antibody Capture": 1,
    }


def test_anndata_summary_is_strict_bounded_read_only_and_code_equivalent(science):
    value = science.ad.AnnData(
        X=science.sparse.csr_matrix(science.np.eye(3)),
        obs=science.pd.DataFrame(index=["first", "second", "third"]),
        var=science.pd.DataFrame(index=["a", "b", "c"]),
    )
    value.layers["counts"] = value.X.copy()
    value.obsm["X_test"] = science.np.ones((3, 2))
    value.varm["loadings"] = science.np.ones((3, 2))
    value.obsp["connectivities"] = science.sparse.eye(3, format="csr")
    value.varp["similarity"] = science.sparse.eye(3, format="csr")
    value.uns["method"] = {"value": 1}
    value.raw = value.copy()
    value.obs_names = ["cell", "cell", "third"]
    before_x = value.X.copy()
    before_obs = value.obs.copy(deep=True)

    report, code = _execute_anndata_summary(value).result
    key_results = report.summary["key_results"]

    json.dumps(report.summary, allow_nan=False)
    payload = result_to_payload(report)
    assert key_results["shape"] == [3, 3]
    assert key_results["X"]["storage"] == "sparse_csr"
    assert key_results["obs_index"]["is_unique"] is False
    assert key_results["raw"]["present"] is True
    assert key_results["layers"]["names"] == ["counts"]
    assert payload["summary"]["key_results"]["layers"]["names"] == ["counts"]
    assert key_results["obsp"]["names"] == ["connectivities"]
    assert key_results["varp"]["names"] == ["similarity"]
    assert any("not unique" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(compile(code, "<anndata-summary>", "exec"), namespace)
    assert namespace["summarize_anndata"](value) == key_results
    assert (value.X != before_x).nnz == 0
    assert value.obs.equals(before_obs)


def test_anndata_summary_bounds_names_without_mutating_invalid_metadata(science):
    value = science.ad.AnnData(science.np.eye(2))
    for index in range(70):
        value.uns[f"key_{index:02d}_" + "x" * 300] = index
    value.uns["openbio_singlecell"] = ["invalid"]
    before = list(value.uns["openbio_singlecell"])

    report, _ = _execute_anndata_summary(value).result
    inventory = report.summary["key_results"]["uns"]

    assert inventory["total"] == 71
    assert len(inventory["names"]) == 64
    assert inventory["truncated"] is True
    assert max(map(len, inventory["names"])) <= 256
    assert report.summary["key_results"]["openbio_metadata"]["valid_mapping"] is False
    assert value.uns["openbio_singlecell"] == before


def test_anndata_summary_bounds_embedded_display_name_and_warnings(science):
    value = science.ad.AnnData(science.np.eye(2))
    value.uns["openbio_singlecell"] = {
        "display_name": "d" * 500,
        "warnings": [f"warning-{index}-" + "x" * 500 for index in range(40)],
        "source": {"kind": "source-" + "s" * 500},
        "analysis_history": {},
    }

    report, _ = _execute_anndata_summary(value).result
    metadata = report.summary["key_results"]["openbio_metadata"]

    assert len(metadata["display_name"]) == 256
    assert len(metadata["source_kind"]) == 256
    assert len(report.title) == 264
    assert len(report.summary["warnings"]) == 33
    assert all(len(warning) <= 256 for warning in report.summary["warnings"][:-1])
    assert "limited to the first 32" in report.summary["warnings"][-1]
