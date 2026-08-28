from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import PLUGIN_VERSION
from .analysis_utils import finish_adata, make_summary_result
from .node_types import AnnDataType, SummaryResultType
from .schist_analysis import run_schist_nested_model, schist_nested_model_code

if TYPE_CHECKING:
    from anndata import AnnData


SCHIST_CATEGORY = "openbio/single-cell/clustering"


class OpenBioSingleCellSchistNestedModel(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSchistNestedModel",
            display_name="Schist Nested-SBM Hierarchy",
            category=SCHIST_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("random_seed", default=123, min=1, max=2**31 - 1, advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.String.Input("key_added", default="nsbm", advanced=True),
                io.Int.Input("posterior_samples", default=100, min=100, max=2**31 - 1, advanced=True),
                io.Boolean.Input("degree_correction", default=True, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_working_memory_gib", default=8.0, min=0.001, step=0.5, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        random_seed: int = 123,
        neighbors_key: str = "neighbors",
        key_added: str = "nsbm",
        posterior_samples: int = 100,
        degree_correction: bool = True,
        overwrite_existing: bool = False,
        max_working_memory_gib: float = 8.0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        output, summary = run_schist_nested_model(
            adata,
            random_seed=random_seed,
            neighbors_key=neighbors_key,
            key_added=key_added,
            posterior_samples=posterior_samples,
            degree_correction=degree_correction,
            overwrite_existing=overwrite_existing,
            max_working_memory_gib=max_working_memory_gib,
            openbio_version=PLUGIN_VERSION,
        )
        code = schist_nested_model_code(
            random_seed=random_seed,
            neighbors_key=neighbors_key,
            key_added=key_added,
            posterior_samples=posterior_samples,
            degree_correction=degree_correction,
            overwrite_existing=overwrite_existing,
            max_working_memory_gib=max_working_memory_gib,
            openbio_version=PLUGIN_VERSION,
        )
        finish_adata(
            output,
            "schist_nested_sbm_hierarchy",
            dict(summary["parameters"]),
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
            random_seed=random_seed,
            warnings=list(summary["warnings"]),
        )
        report = make_summary_result(
            summary=summary,
            title="Schist Nested-SBM hierarchy summary",
            operation="schist_nested_sbm_hierarchy",
            parameters=dict(summary["parameters"]),
            description=str(summary["results"]),
            warnings=list(summary["warnings"]),
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output, report, code)


SCHIST_NODE_CLASSES = [OpenBioSingleCellSchistNestedModel]


__all__ = [
    "OpenBioSingleCellSchistNestedModel",
    "SCHIST_NODE_CLASSES",
]
