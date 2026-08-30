from comfy_api.latest import io

AnnDataType = io.Custom("OPENBIO_ANNDATA")
AugurResultType = io.Custom("OPENBIO_AUGUR_RESULT")
DGIdbResourceType = io.Custom("OPENBIO_DGIDB_RESOURCE")
LianaResultType = io.Custom("OPENBIO_LIANA_RESULT")
SummaryResultType = io.Custom("OPENBIO_SINGLE_CELL_SUMMARY")
TableResultType = io.Custom("OPENBIO_SINGLE_CELL_TABLE")
TFActivityArtifactType = io.Custom("OPENBIO_TF_ACTIVITY")
PlotResultType = io.Custom("OPENBIO_SINGLE_CELL_PLOT")
PseudobulkType = io.Custom("OPENBIO_SINGLE_CELL_PSEUDOBULK")
CNMFRunType = io.Custom("OPENBIO_CNMF_RUN")
SCVIModelType = io.Custom("OPENBIO_SCVI_MODEL")
SCENICResultArtifactType = io.Custom("OPENBIO_SCENIC_RESULT")
CassiopeiaCharactersType = io.Custom("OPENBIO_CASSIOPEIA_CHARACTERS")
CassiopeiaTreeType = io.Custom("OPENBIO_CASSIOPEIA_TREE")
CNVStateType = io.Custom("OPENBIO_CNV_STATE")
VelocityStateType = io.Custom("OPENBIO_VELOCITY_STATE")
WorkerType = io.Custom("OPENBIO_WORKER")


def analysis_outputs(*primary: io.Output) -> list[io.Output]:
    return [*primary, SummaryResultType.Output(display_name="summary"), io.String.Output("code")]


__all__ = [
    "analysis_outputs",
    "AnnDataType",
    "AugurResultType",
    "CassiopeiaCharactersType",
    "CassiopeiaTreeType",
    "CNMFRunType",
    "CNVStateType",
    "DGIdbResourceType",
    "LianaResultType",
    "PlotResultType",
    "PseudobulkType",
    "SCVIModelType",
    "SCENICResultArtifactType",
    "SummaryResultType",
    "TableResultType",
    "TFActivityArtifactType",
    "VelocityStateType",
    "WorkerType",
]
