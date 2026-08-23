from comfy_api.latest import io

AnnDataType = io.Custom("OPENBIO_ANNDATA")
SingleCellResultType = io.Custom("OPENBIO_SINGLE_CELL_RESULT")

__all__ = ["AnnDataType", "SingleCellResultType"]
