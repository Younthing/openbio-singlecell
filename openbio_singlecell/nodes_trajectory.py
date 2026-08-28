from __future__ import annotations

import math
import time
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from .analysis_utils import finish_adata, make_summary_result
from .cell_cycle import analyze_cell_cycle_score, cell_cycle_code
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, SummaryResultType
from .trajectory_analysis import (
    analyze_diffusion_map,
    analyze_dpt,
    analyze_paga,
    diffusion_map_code,
    dpt_code,
    paga_code,
)

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/trajectory"
ANNOTATION_CATEGORY = "openbio/single-cell/annotation"


def _parse_dpt_root_value(value_type: object, value: object) -> str | int | float:
    if value_type not in {"string", "integer", "number"}:
        raise ValueError("DPT root_value_type must be 'string', 'integer', or 'number'.")
    if not isinstance(value, str):
        raise TypeError("DPT root_value must be entered as text and interpreted by root_value_type.")
    if not value or value != value.strip():
        raise ValueError("DPT root_value must be nonempty and free of surrounding whitespace.")
    if value_type == "string":
        return value
    if value_type == "integer":
        try:
            parsed = int(value, 10)
        except ValueError as error:
            raise ValueError("DPT integer root_value must use canonical base-10 syntax.") from error
        if str(parsed) != value:
            raise ValueError("DPT integer root_value must use canonical base-10 syntax.")
        return parsed
    try:
        parsed_number = float(value)
    except ValueError as error:
        raise ValueError("DPT number root_value must be a finite numeric literal.") from error
    if not math.isfinite(parsed_number):
        raise ValueError("DPT number root_value must be a finite numeric literal.")
    return parsed_number


class OpenBioSingleCellCellCycleScore(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Cell-cycle expression source; full-gene log-normalized values are recommended",
        default="layer",
        include_raw=True,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellCycleScore",
            display_name="Cell Cycle Score",
            category=ANNOTATION_CATEGORY,
            description=(
                "Score one versioned human or custom cell-cycle program on an explicit Raw, X, or layer source; "
                "the report discloses departures from recommended full-gene log-normalized input."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.DynamicCombo.Input(
                    "gene_set_source",
                    options=[
                        io.DynamicCombo.Option("regev_human_97", []),
                        io.DynamicCombo.Option(
                            "custom",
                            [
                                io.String.Input("s_genes", default=""),
                                io.String.Input("g2m_genes", default=""),
                            ],
                        ),
                    ],
                ),
                io.Combo.Input("organism", options=["human", "mouse", "other"], default="human"),
                io.String.Input("output_prefix", default="cell_cycle", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
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
        source: DynamicExpressionSource | None = None,
        gene_set_source: Mapping[str, object] | None = None,
        organism: str = "human",
        output_prefix: str = "cell_cycle",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        if gene_set_source is None:
            gene_set_source = {"gene_set_source": "regev_human_97"}
        if not isinstance(gene_set_source, Mapping):
            raise TypeError("Cell-cycle gene_set_source must be a DynamicCombo value.")
        mode = gene_set_source.get("gene_set_source")
        if mode not in {"regev_human_97", "custom"}:
            raise ValueError(f"Unsupported cell-cycle gene-set source: {mode!r}.")
        custom_s: Any = gene_set_source.get("s_genes") if mode == "custom" else None
        custom_g2m: Any = gene_set_source.get("g2m_genes") if mode == "custom" else None
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, summary = analyze_cell_cycle_score(
            adata,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            gene_set_source=str(mode),
            organism=organism,
            s_genes=custom_s,
            g2m_genes=custom_g2m,
            output_prefix=output_prefix,
            overwrite_existing=overwrite_existing,
            random_seed=random_seed,
        )
        parameters = dict(summary["parameters"])
        warnings = [str(warning) for warning in summary["warnings"]]
        finish_adata(
            output,
            "cell_cycle_score",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=warnings,
        )
        report = make_summary_result(
            summary=summary,
            title="Cell-cycle program score summary",
            operation="cell_cycle_score",
            parameters=parameters,
            description=str(summary["results"]),
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = cell_cycle_code(
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            gene_set_source=str(mode),
            organism=organism,
            s_genes=custom_s,
            g2m_genes=custom_g2m,
            output_prefix=output_prefix,
            overwrite_existing=overwrite_existing,
            random_seed=random_seed,
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellDiffusionMap(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDiffusionMap",
            display_name="Diffusion Map",
            category=CATEGORY,
            description=(
                "Compute and fingerprint one diffusion basis from an exact named neighbor graph; no root or "
                "pseudotime is selected."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Int.Input("n_comps", default=15, min=3, max=4096),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
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
        neighbors_key: str = "neighbors",
        n_comps: int = 15,
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, summary = analyze_diffusion_map(
            adata,
            neighbors_key=neighbors_key,
            n_comps=n_comps,
            overwrite_existing=overwrite_existing,
            random_seed=random_seed,
        )
        parameters = dict(summary["parameters"])
        report_warnings = [str(value) for value in summary["warnings"]]
        finish_adata(
            output,
            "diffusion_map",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=report_warnings,
        )
        report = make_summary_result(
            summary=summary,
            title="Diffusion Map summary",
            operation="diffusion_map",
            parameters=parameters,
            description=str(summary["results"]),
            warnings=report_warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = diffusion_map_code(
            neighbors_key=str(parameters["neighbors_key"]),
            n_comps=int(parameters["n_comps"]),
            overwrite_existing=bool(parameters["overwrite_existing"]),
            random_seed=int(parameters["random_seed"]),
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellPAGA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPAGA",
            display_name="PAGA",
            category=CATEGORY,
            description=("Coarse-grain one named undirected neighbor graph over one explicit categorical partition."),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
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
        groupby: str = "leiden",
        neighbors_key: str = "neighbors",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, summary = analyze_paga(
            adata,
            groupby=groupby,
            neighbors_key=neighbors_key,
            overwrite_existing=overwrite_existing,
        )
        parameters = dict(summary["parameters"])
        report_warnings = [str(value) for value in summary["warnings"]]
        finish_adata(output, "paga", parameters, cells, genes, started_at, warnings=report_warnings)
        report = make_summary_result(
            summary=summary,
            title="PAGA connectivity summary",
            operation="paga",
            parameters=parameters,
            description=str(summary["results"]),
            warnings=report_warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
        )
        code = paga_code(
            groupby=str(parameters["groupby"]),
            neighbors_key=str(parameters["neighbors_key"]),
            overwrite_existing=bool(parameters["overwrite_existing"]),
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellDPT(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDPT",
            display_name="Diffusion Pseudotime",
            category=CATEGORY,
            description=(
                "Orient an already graph-bound Diffusion Map from one exact cell or deterministic nominated-"
                "population medoid."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.DynamicCombo.Input(
                    "root_mode",
                    options=[
                        io.DynamicCombo.Option(
                            "cell_id",
                            [io.String.Input("root_cell_id", default="")],
                        ),
                        io.DynamicCombo.Option(
                            "group_medoid",
                            [
                                io.String.Input("root_column", default="leiden"),
                                io.Combo.Input(
                                    "root_value_type",
                                    options=["string", "integer", "number"],
                                    default="string",
                                ),
                                io.String.Input("root_value", default=""),
                            ],
                        ),
                    ],
                ),
                io.Int.Input("n_dcs", default=10, min=2, max=4096),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
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
        neighbors_key: str = "neighbors",
        root_mode: Mapping[str, object] | None = None,
        n_dcs: int = 10,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        if root_mode is None:
            root_mode = {"root_mode": "cell_id", "root_cell_id": ""}
        if not isinstance(root_mode, Mapping):
            raise TypeError("DPT root_mode must be a DynamicCombo value.")
        mode = root_mode.get("root_mode")
        if mode not in {"cell_id", "group_medoid"}:
            raise ValueError(f"Unsupported DPT root mode: {mode!r}.")
        root_cell_id = root_mode.get("root_cell_id", "") if mode == "cell_id" else ""
        root_column = root_mode.get("root_column", "") if mode == "group_medoid" else ""
        root_value: Any = ""
        if mode == "group_medoid":
            root_value = _parse_dpt_root_value(
                root_mode.get("root_value_type"),
                root_mode.get("root_value"),
            )
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, summary = analyze_dpt(
            adata,
            neighbors_key=neighbors_key,
            root_mode=str(mode),
            root_cell_id=root_cell_id,
            root_column=root_column,
            root_value=root_value,
            n_dcs=n_dcs,
            overwrite_existing=overwrite_existing,
        )
        parameters = dict(summary["parameters"])
        report_warnings = [str(value) for value in summary["warnings"]]
        finish_adata(output, "dpt", parameters, cells, genes, started_at, warnings=report_warnings)
        report = make_summary_result(
            summary=summary,
            title="Diffusion Pseudotime summary",
            operation="dpt",
            parameters=parameters,
            description=str(summary["results"]),
            warnings=report_warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
        )
        code = dpt_code(
            neighbors_key=str(parameters["neighbors_key"]),
            root_mode=str(parameters["root_mode"]),
            root_cell_id=str(parameters["root_cell_id"]),
            root_column=str(parameters["root_column"] or ""),
            root_value=root_value,
            n_dcs=int(parameters["n_dcs"]),
            overwrite_existing=bool(parameters["overwrite_existing"]),
        )
        return io.NodeOutput(output, report, code)


TRAJECTORY_NODE_CLASSES = [
    OpenBioSingleCellCellCycleScore,
    OpenBioSingleCellDiffusionMap,
    OpenBioSingleCellPAGA,
    OpenBioSingleCellDPT,
]
