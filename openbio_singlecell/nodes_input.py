from __future__ import annotations

import gzip
import os
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .contracts import SingleCellResult, ensure_metadata, make_analysis_source
from .files import input_file_fingerprint, resolve_10x_mtx_files, resolve_input_path, tenx_mtx_fingerprint
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/input"


def _validate_nonempty(adata: AnnData) -> None:
    if int(adata.n_obs) == 0 or int(adata.n_vars) == 0:
        raise ValueError("The loaded single-cell matrix must contain at least one cell and one gene.")


def _fingerprint_metadata(fingerprint: Any) -> dict[str, Any]:
    if len(fingerprint) == 3 and isinstance(fingerprint[0], str):
        return {"path": fingerprint[0], "size": fingerprint[1], "mtime_ns": fingerprint[2]}
    return {
        role: {"path": canonical, "size": size, "mtime_ns": mtime_ns} for role, canonical, size, mtime_ns in fingerprint
    }


def _source_metadata(kind: str, path: str, fingerprint: Any) -> dict[str, Any]:
    return {"kind": kind, "path": path, "fingerprint": _fingerprint_metadata(fingerprint)}


def read_10x_mtx(
    relative_directory: str,
    var_names: str = "gene_symbols",
    make_unique: bool = True,
    gex_only: bool = True,
) -> AnnData:
    dependencies.require_scientific_dependencies()
    files = resolve_10x_mtx_files(relative_directory)

    opener = gzip.open if files["matrix"].lower().endswith(".gz") else open
    with opener(files["matrix"], "rb") as matrix_file:
        matrix = dependencies.mmread(matrix_file).tocsr()

    barcodes = dependencies.pd.read_csv(files["barcodes"], sep="\t", header=None, dtype=str, compression="infer")
    features = dependencies.pd.read_csv(files["features"], sep="\t", header=None, dtype=str, compression="infer")
    if features.shape[1] == 0:
        raise ValueError("The 10x features file is empty.")

    gene_ids = features.iloc[:, 0].astype(str)
    gene_symbols = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else gene_ids.copy()
    feature_types = features.iloc[:, 2].astype(str) if features.shape[1] > 2 else None

    if matrix.shape[0] != len(features) or matrix.shape[1] != len(barcodes):
        raise ValueError(
            "10x matrix dimensions do not match the features and barcodes files: "
            f"matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}."
        )

    if gex_only and feature_types is not None:
        mask = feature_types.to_numpy() == "Gene Expression"
        matrix = matrix[mask, :]
        gene_ids = gene_ids[mask].reset_index(drop=True)
        gene_symbols = gene_symbols[mask].reset_index(drop=True)
        feature_types = feature_types[mask].reset_index(drop=True)

    selected_names = gene_ids if var_names == "gene_ids" else gene_symbols
    obs = dependencies.pd.DataFrame(index=barcodes.iloc[:, 0].astype(str).to_numpy())
    var = dependencies.pd.DataFrame(index=selected_names.astype(str).to_numpy())
    var["gene_ids"] = gene_ids.astype(str).to_numpy()
    var["gene_symbols"] = gene_symbols.astype(str).to_numpy()
    if feature_types is not None:
        var["feature_types"] = feature_types.astype(str).to_numpy()

    adata = dependencies.ad.AnnData(X=matrix.transpose().tocsr(), obs=obs, var=var)
    if make_unique:
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
    _validate_nonempty(adata)
    return adata


class OpenBioSingleCellLoadH5AD(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoadH5AD",
            display_name="Load H5AD",
            category=CATEGORY,
            description="Load an AnnData H5AD file from the ComfyUI input directory.",
            inputs=[
                io.String.Input("path", default="openbio-singlecell/openbio_singlecell_demo.h5ad"),
                io.Boolean.Input("make_var_names_unique", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, path: str, make_var_names_unique: bool = True) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5ad",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, path: str, make_var_names_unique: bool = True) -> Any:
        return input_file_fingerprint(path, (".h5ad",))

    @classmethod
    def execute(cls, path: str, make_var_names_unique: bool = True) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        resolved = resolve_input_path(path, extensions=(".h5ad",))
        adata = dependencies.ad.read_h5ad(resolved)
        if make_var_names_unique:
            adata.var_names_make_unique()
        _validate_nonempty(adata)
        ensure_metadata(
            adata,
            display_name=os.path.splitext(os.path.basename(resolved))[0],
            source=_source_metadata("h5ad", path, input_file_fingerprint(path, (".h5ad",))),
        )
        return io.NodeOutput(adata)


class OpenBioSingleCellLoad10xMTX(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xMTX",
            display_name="Load 10x MTX",
            category=CATEGORY,
            description="Load a 10x matrix, barcodes, and features directory under ComfyUI input.",
            inputs=[
                io.String.Input("directory", default="openbio-singlecell/10x"),
                io.Combo.Input("var_names", options=["gene_symbols", "gene_ids"], default="gene_symbols"),
                io.Boolean.Input("make_unique", default=True, advanced=True),
                io.Boolean.Input("gex_only", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, directory: str, **kwargs) -> bool | str:
        try:
            resolve_10x_mtx_files(directory)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, directory: str, **kwargs) -> Any:
        return tenx_mtx_fingerprint(directory)

    @classmethod
    def execute(
        cls,
        directory: str,
        var_names: str = "gene_symbols",
        make_unique: bool = True,
        gex_only: bool = True,
    ) -> io.NodeOutput:
        adata = read_10x_mtx(directory, var_names, make_unique, gex_only)
        ensure_metadata(
            adata,
            display_name=os.path.basename(os.path.normpath(directory)) or "10x",
            source=_source_metadata("10x_mtx", directory, tenx_mtx_fingerprint(directory)),
        )
        return io.NodeOutput(adata)


class OpenBioSingleCellLoad10xH5(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xH5",
            display_name="Load 10x H5",
            category=CATEGORY,
            description="Load a 10x Genomics HDF5 matrix from the ComfyUI input directory.",
            inputs=[
                io.String.Input("path", default="openbio-singlecell/filtered_feature_bc_matrix.h5"),
                io.String.Input("genome", default="", advanced=True),
                io.Boolean.Input("gex_only", default=True, advanced=True),
                io.Boolean.Input("make_unique", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, path: str, **kwargs) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5", ".hdf5"))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, path: str, **kwargs) -> Any:
        return input_file_fingerprint(path, (".h5", ".hdf5"))

    @classmethod
    def execute(
        cls,
        path: str,
        genome: str = "",
        gex_only: bool = True,
        make_unique: bool = True,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        resolved = resolve_input_path(path, extensions=(".h5", ".hdf5"))
        adata = dependencies.sc.read_10x_h5(resolved, genome=genome or None, gex_only=gex_only)
        if make_unique:
            adata.obs_names_make_unique()
            adata.var_names_make_unique()
        _validate_nonempty(adata)
        ensure_metadata(
            adata,
            display_name=os.path.splitext(os.path.basename(resolved))[0],
            source=_source_metadata("10x_h5", path, input_file_fingerprint(path, (".h5", ".hdf5"))),
        )
        return io.NodeOutput(adata)


def _limited_names(values: Any) -> dict[str, Any]:
    names = [str(value) for value in list(values)]
    return {"names": names[:64], "total": len(names), "truncated": len(names) > 64}


class OpenBioSingleCellAnnDataSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAnnDataSummary",
            display_name="AnnData Summary",
            category=CATEGORY,
            description="Create a bounded structural summary of AnnData.",
            inputs=[AnnDataType.Input("adata")],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        start = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        summary = {
            "representation": repr(adata),
            "shape": [cells, genes],
            "obs": _limited_names(adata.obs.columns),
            "var": _limited_names(adata.var.columns),
            "layers": _limited_names(adata.layers.keys()),
            "obsm": _limited_names(adata.obsm.keys()),
            "varm": _limited_names(adata.varm.keys()),
            "uns": _limited_names(adata.uns.keys()),
        }
        elapsed = time.perf_counter() - start
        source = make_analysis_source("anndata_summary", {}, cells, genes, 0, elapsed)
        metadata = adata.uns.get("openbio_singlecell", {})
        result = SingleCellResult(
            kind="summary",
            title=f"{metadata.get('display_name', 'AnnData')} summary",
            parameters={},
            description="AnnData structure and annotations.",
            warnings=[str(value) for value in metadata.get("warnings", [])],
            input_cells=cells,
            input_genes=genes,
            random_seed=int(metadata.get("random_seed", 0)),
            elapsed_seconds=elapsed,
            source=source,
            summary=summary,
        )
        return io.NodeOutput(result)


INPUT_NODE_CLASSES = [
    OpenBioSingleCellLoadH5AD,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellAnnDataSummary,
]
