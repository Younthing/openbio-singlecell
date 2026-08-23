from __future__ import annotations

import time
from typing import TYPE_CHECKING, Literal

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/correction"
MAX_RANDOM_SEED = 2**31 - 1


def _comma_separated_metrics(value: str) -> list[str]:
    metrics = list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))
    if not metrics:
        raise ValueError("At least one QC metric is required.")
    return metrics


def _mad_outliers(
    values: object,
    *,
    nmads: float,
    direction: Literal["both", "upper", "lower"],
    scale_mad: bool,
) -> object:
    science = dependencies.require_scientific_dependencies()
    from scipy.stats import median_abs_deviation

    array = science.np.asarray(values, dtype=float)
    median = science.np.median(array)
    scale = 1.4826 if scale_mad else 1.0
    mad = median_abs_deviation(array, scale=scale)

    lower = -science.np.inf if direction == "upper" else median - nmads * mad
    upper = science.np.inf if direction == "lower" else median + nmads * mad
    return (array < lower) | (array > upper)


class OpenBioSingleCellMarkMADOutliers(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkMADOutliers",
            display_name="Mark MAD Outliers",
            category=CATEGORY,
            description="Mark metric outliers, optionally within each batch.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input(
                    "metrics",
                    default="total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
                ),
                io.String.Input("batch_key", default="batch"),
                io.Float.Input("nmads", default=5.0, min=0.1, max=100.0, step=0.1),
                io.Combo.Input("direction", options=["both", "upper", "lower"], default="both"),
                io.String.Input("output_column", default="outlier", advanced=True),
                io.Boolean.Input("scale_mad", default=False, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        metrics: str = "total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
        batch_key: str = "batch",
        nmads: float = 5.0,
        direction: Literal["both", "upper", "lower"] = "both",
        output_column: str = "outlier",
        scale_mad: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        selected_metrics = _comma_separated_metrics(metrics)
        missing_metrics = [metric for metric in selected_metrics if metric not in output.obs]
        if missing_metrics:
            raise ValueError(f"QC metrics not found in obs: {missing_metrics}")
        if batch_key and batch_key not in output.obs:
            raise ValueError(f"QC batch column not found in obs: {batch_key!r}")

        flags = science.pd.Series(False, index=output.obs_names)
        for metric in selected_metrics:
            metric_flags = science.pd.Series(False, index=output.obs_names)
            if batch_key:
                groups = output.obs.groupby(batch_key, observed=False).groups.values()
                for indices in groups:
                    metric_flags.loc[indices] = _mad_outliers(
                        output.obs.loc[indices, metric],
                        nmads=nmads,
                        direction=direction,
                        scale_mad=scale_mad,
                    )
            else:
                metric_flags[:] = _mad_outliers(
                    output.obs[metric],
                    nmads=nmads,
                    direction=direction,
                    scale_mad=scale_mad,
                )
            flags |= metric_flags
        output.obs[output_column] = flags

        parameters = {
            "metrics": selected_metrics,
            "batch_key": batch_key,
            "nmads": nmads,
            "direction": direction,
            "output_column": output_column,
            "scale_mad": scale_mad,
        }
        finish_adata(output, "mark_mad_outliers", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellScrublet(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScrublet",
            display_name="Run Scrublet",
            category=CATEGORY,
            description="Add doublet_score and predicted_doublet to obs.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("batch_key", default="batch"),
                io.Int.Input("random_seed", default=123, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        batch_key: str = "batch",
        random_seed: int = 123,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        output.obs_names_make_unique()
        science.sc.pp.scrublet(
            output,
            batch_key=batch_key or None,
            random_state=random_seed,
        )
        parameters = {"batch_key": batch_key, "random_seed": random_seed}
        finish_adata(
            output,
            "scrublet",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellFilterDoublets(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterDoublets",
            display_name="Filter Predicted Doublets",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("prediction_column", default="predicted_doublet", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        prediction_column: str = "predicted_doublet",
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata[~adata.obs[prediction_column].astype(bool)].copy()
        parameters = {"prediction_column": prediction_column}
        finish_adata(output, "filter_doublets", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


CORRECTION_NODE_CLASSES = [
    OpenBioSingleCellMarkMADOutliers,
    OpenBioSingleCellScrublet,
    OpenBioSingleCellFilterDoublets,
]


__all__ = ["CORRECTION_NODE_CLASSES"]
