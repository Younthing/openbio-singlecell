from __future__ import annotations

import os

import pytest

from openbio_singlecell.files import (
    file_fingerprint,
    input_file_provenance,
    prepare_output_target,
    resolve_10x_mtx_files,
    resolve_input_path,
    tenx_mtx_fingerprint,
    tenx_mtx_provenance,
)


def test_input_path_is_contained_and_rejects_annotations(comfy_directories):
    input_dir, output_dir, _ = comfy_directories
    nested = input_dir / "nested"
    nested.mkdir()
    source = nested / "sample.h5ad"
    source.write_bytes(b"h5ad")
    (output_dir / "outside.h5ad").write_bytes(b"outside")

    assert resolve_input_path("nested/sample.h5ad", extensions=(".h5ad",)) == os.path.realpath(source)
    with pytest.raises(ValueError):
        resolve_input_path("../output/outside.h5ad", extensions=(".h5ad",))
    with pytest.raises(ValueError):
        resolve_input_path("outside.h5ad [output]", extensions=(".h5ad",))


def test_input_symlink_escape_is_rejected(comfy_directories, tmp_path):
    input_dir, _, _ = comfy_directories
    outside = tmp_path / "outside.h5ad"
    outside.write_bytes(b"outside")
    link = input_dir / "linked.h5ad"
    try:
        link.symlink_to(outside)
    except OSError:
        pytest.skip("Symbolic links are unavailable on this platform.")

    with pytest.raises(ValueError):
        resolve_input_path("linked.h5ad", extensions=(".h5ad",))


def test_file_fingerprint_changes_with_stat(comfy_directories):
    input_dir, _, _ = comfy_directories
    source = input_dir / "sample.h5"
    source.write_bytes(b"one")
    first = file_fingerprint(str(source))
    provenance = input_file_provenance("sample.h5", (".h5",))
    assert os.path.isabs(first[0])
    assert provenance["path"] == "sample.h5"
    assert str(input_dir) not in str(provenance)
    source.write_bytes(b"two-two")
    stat = source.stat()
    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
    second = file_fingerprint(str(source))
    assert first != second


def test_10x_selection_and_fingerprint_use_the_same_actual_files(comfy_directories):
    input_dir, _, _ = comfy_directories
    directory = input_dir / "tenx"
    directory.mkdir()
    (directory / "matrix.mtx").write_bytes(b"plain-matrix")
    (directory / "matrix.mtx.gz").write_bytes(b"compressed-matrix")
    (directory / "barcodes.tsv").write_text("cell\n", encoding="utf-8")
    (directory / "features.tsv").write_text("id\tgene\n", encoding="utf-8")

    selected = resolve_10x_mtx_files("tenx")
    fingerprint = tenx_mtx_fingerprint("tenx")
    provenance = tenx_mtx_provenance("tenx")
    assert selected["matrix"].endswith("matrix.mtx")
    assert fingerprint[0][1] == file_fingerprint(selected["matrix"])[0]
    assert provenance["files"]["matrix"]["path"] == "tenx/matrix.mtx"
    assert str(input_dir) not in str(provenance)

    (directory / "features.tsv").write_text("id\tgene\tGene Expression\n", encoding="utf-8")
    stat = (directory / "features.tsv").stat()
    os.utime(directory / "features.tsv", ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
    assert tenx_mtx_fingerprint("tenx") != fingerprint


def test_output_increment_overwrite_and_containment(comfy_directories):
    _, _, _ = comfy_directories
    first = prepare_output_target("nested/table", "csv")
    os.makedirs(os.path.dirname(first.path), exist_ok=True)
    with open(first.path, "wb") as handle:
        handle.write(b"first")
    second = prepare_output_target("nested/table", "csv")
    overwrite = prepare_output_target("nested/table", "csv", overwrite=True)

    assert first.filename == "table_00001.csv"
    assert second.filename == "table_00002.csv"
    assert overwrite.filename == "table.csv"
    assert first.subfolder == "openbio-singlecell/nested"
    with pytest.raises(ValueError):
        prepare_output_target("../../escape", "csv")
