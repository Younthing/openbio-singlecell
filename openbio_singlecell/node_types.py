from comfy_api.latest import io

AnnDataType = io.Custom("OPENBIO_ANNDATA")
SummaryResultType = io.Custom("OPENBIO_SINGLE_CELL_SUMMARY")
TableResultType = io.Custom("OPENBIO_SINGLE_CELL_TABLE")
PlotResultType = io.Custom("OPENBIO_SINGLE_CELL_PLOT")

__all__ = ["AnnDataType", "PlotResultType", "SummaryResultType", "TableResultType"]
