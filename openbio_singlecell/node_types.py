from comfy_api.latest import io

AnnDataType = io.Custom("OPENBIO_ANNDATA")
SummaryResultType = io.Custom("OPENBIO_SINGLE_CELL_SUMMARY")
TableResultType = io.Custom("OPENBIO_SINGLE_CELL_TABLE")
PlotResultType = io.Custom("OPENBIO_SINGLE_CELL_PLOT")
SCVIModelType = io.Custom("OPENBIO_SCVI_MODEL")
ScenicNetworkType = io.Custom("OPENBIO_SCENIC_NETWORK")
CassiopeiaTreeType = io.Custom("OPENBIO_CASSIOPEIA_TREE")

__all__ = [
    "AnnDataType",
    "CassiopeiaTreeType",
    "PlotResultType",
    "SCVIModelType",
    "ScenicNetworkType",
    "SummaryResultType",
    "TableResultType",
]
