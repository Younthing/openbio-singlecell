from __future__ import annotations

import hashlib
import json

from comfy_api.latest import io

from .monocle2 import MONOCLE2_COMPAT_SCRIPT, MONOCLE2_PACKAGES, MONOCLE2_SCRIPT, MONOCLE2_SOURCE
from .node_types import AnnDataType, Monocle2CDSType, PlotResultType, TableResultType, analysis_outputs
from .r_runtime import PROBE_SCRIPT, runtime_fingerprint, validate_r_runtime

CATEGORY = "openbio/single-cell/trajectory/monocle2"


def _script_fingerprint() -> tuple[str, ...]:
    return (
        hashlib.sha256(MONOCLE2_SCRIPT.read_bytes()).hexdigest(),
        hashlib.sha256(MONOCLE2_COMPAT_SCRIPT.read_bytes()).hexdigest(),
        hashlib.sha256(PROBE_SCRIPT.read_bytes()).hexdigest(),
    )


def _analysis_fingerprint(cls, r_runtime: str = "", **kwargs):
    runtime = validate_r_runtime(json.loads(r_runtime)) if r_runtime else r_runtime
    return json.dumps(runtime, sort_keys=True, separators=(",", ":")), _script_fingerprint()


class OpenBioSingleCellMonocle2Runtime(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2Runtime",
            display_name="Monocle 2 Runtime",
            category=CATEGORY,
            description=(
                "Probe Rscript and the installed Monocle 2 packages for downstream R nodes. Optional library_paths "
                "and path_prefix accept one directory per line and affect only child R processes."
            ),
            inputs=[
                io.String.Input("rscript", default="Rscript"),
                io.String.Input("r_home", default="", advanced=True),
                io.String.Input("library_paths", default="", advanced=True),
                io.String.Input("path_prefix", default="", advanced=True),
            ],
            outputs=analysis_outputs(io.String.Output("r_runtime")),
        )

    @classmethod
    def fingerprint_inputs(cls, rscript="Rscript", r_home="", library_paths="", path_prefix="", **kwargs):
        return runtime_fingerprint(
            rscript, r_home=r_home, library_paths=library_paths, path_prefix=path_prefix, packages=MONOCLE2_PACKAGES,
        ), _script_fingerprint()


class OpenBioSingleCellMonocle2Prepare(io.ComfyNode):
    EXPRESSION_SOURCE = MONOCLE2_SOURCE
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2Prepare",
            display_name="Monocle 2 Prepare",
            category=CATEGORY,
            description=(
                "Create a native R CellDataSet from the chosen expression source. Choose the expression family "
                "and preprocessing explicitly; Raw does not imply counts. Empty gene_short_name_column uses feature IDs. "
                "Advanced JSON arguments override the visible native defaults."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("r_runtime", force_input=True),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("gene_short_name_column", default=""),
                io.Combo.Input(
                    "expression_family",
                    options=["negbinomial.size", "negbinomial", "tobit", "gaussianff"],
                    default="negbinomial.size",
                ),
                io.Float.Input("lower_detection_limit", default=0.1),
                io.Boolean.Input("estimate_size_factors", default=True),
                io.Boolean.Input("estimate_dispersions", default=True),
                io.Boolean.Input("detect_genes", default=True),
                io.String.Input("size_factor_parameters_json", default="{}", advanced=True),
                io.String.Input("dispersion_parameters_json", default="{}", advanced=True),
                io.String.Input("detect_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(Monocle2CDSType.Output(display_name="cds")),
        )


class OpenBioSingleCellMonocle2OrderingGenes(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2OrderingGenes",
            display_name="Monocle 2 Ordering Genes",
            category=CATEGORY,
            description=(
                "Mark ordering genes without dropping other features. Select native dispersion evidence "
                "(mean_expression >= threshold and dispersion_empirical >= dispersion_fit × fold), "
                "an explicit JSON array of feature IDs, or a Boolean feature metadata column."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.Combo.Input("method", options=["dispersion", "explicit", "var_column"], default="dispersion"),
                io.String.Input("genes_json", default="[]"),
                io.String.Input("var_column", default="highly_variable"),
                io.Float.Input("min_mean_expression", default=0.5),
                io.Float.Input("dispersion_fold", default=1.0),
            ],
            outputs=analysis_outputs(
                Monocle2CDSType.Output(display_name="cds"), TableResultType.Output(display_name="table"),
            ),
        )


class OpenBioSingleCellMonocle2DDRTree(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2DDRTree",
            display_name="Monocle 2 DDRTree",
            category=CATEGORY,
            description=(
                "Recompute a DDRTree trajectory from the selected ordering genes using native reduceDimension. "
                "Advanced JSON arguments override visible native defaults; the algorithm remains DDRTree. "
                "Monocle 2 sets its own fixed seed (2016)."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.Int.Input("num_components", default=2, min=1),
                io.Combo.Input("norm_method", options=["log", "vstExprs", "none"], default="log"),
                io.String.Input("residual_formula", default=""),
                io.Float.Input("pseudo_expr", default=1.0),
                io.Boolean.Input("auto_param_selection", default=True),
                io.Boolean.Input("scaling", default=True),
                io.String.Input("extra_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(Monocle2CDSType.Output(display_name="cds")),
        )


class OpenBioSingleCellMonocle2OrderCells(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2OrderCells",
            display_name="Monocle 2 Order Cells",
            category=CATEGORY,
            description=(
                "Run native orderCells. Automatic orientation is not a biological root: inspect State and known "
                "annotations first, then select a root State explicitly. State labels are trajectory segments, "
                "not cell types."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.DynamicCombo.Input(
                    "root",
                    options=[
                        io.DynamicCombo.Option("automatic", []),
                        io.DynamicCombo.Option("state", [io.String.Input("state", default="")]),
                    ],
                ),
                io.Boolean.Input("reverse", default=False),
            ],
            outputs=analysis_outputs(
                Monocle2CDSType.Output(display_name="cds"), TableResultType.Output(display_name="table"),
            ),
        )


class OpenBioSingleCellMonocle2TrajectoryPlot(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2TrajectoryPlot",
            display_name="Monocle 2 Trajectory Plot",
            category=CATEGORY,
            description=(
                "Render the stored trajectory with native plot_cell_trajectory. Color by State, Pseudotime "
                "or a cell metadata column; advanced native arguments also expose gene coloring and override "
                "visible defaults. State labels identify trajectory segments."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.String.Input("color_by", default="State"),
                io.Int.Input("x_dimension", default=1, min=1),
                io.Int.Input("y_dimension", default=2, min=1),
                io.Boolean.Input("show_tree", default=True),
                io.Boolean.Input("show_branch_points", default=True),
                io.Boolean.Input("show_state_number", default=False),
                io.Float.Input("cell_size", default=1.5),
                io.Float.Input("width", default=8.0),
                io.Float.Input("height", default=6.0),
                io.Int.Input("dpi", default=150, min=1),
                io.String.Input("extra_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellMonocle2DifferentialTest(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2DifferentialTest",
            display_name="Monocle 2 Differential Test",
            category=CATEGORY,
            description=(
                "Run native differentialGeneTest with explicit full and reduced formulas. An empty gene array "
                "tests all features. Pseudotime association is exploratory evidence, not replicate-aware Condition "
                "inference. Advanced JSON arguments override visible native defaults."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.String.Input("genes_json", default="[]"),
                io.String.Input("full_formula", default="~sm.ns(Pseudotime, df=3)"),
                io.String.Input("reduced_formula", default="~1"),
                io.Boolean.Input("relative_expr", default=True),
                io.Int.Input("cores", default=1, min=1),
                io.String.Input("extra_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellMonocle2BEAM(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2BEAM",
            display_name="Monocle 2 BEAM",
            category=CATEGORY,
            description=(
                "Test branch-dependent expression with native BEAM at a user-selected branch point. "
                "An empty gene array tests all features. Branch evidence is not replicate-aware Condition inference. "
                "Advanced JSON arguments override visible native defaults."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.String.Input("genes_json", default="[]"),
                io.Int.Input("branch_point", default=1, min=1),
                io.Combo.Input(
                    "progenitor_method", options=["duplicate", "sequential_split"], default="duplicate",
                ),
                io.Int.Input("cores", default=1, min=1),
                io.String.Input("extra_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellMonocle2GeneTrends(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2GeneTrends",
            display_name="Monocle 2 Gene Trends",
            category=CATEGORY,
            description=(
                "Plot genes in pseudotime with the native trend formula. Specify feature IDs as a JSON array, "
                "or connect a differential result table and select its top_n genes by q-value. The report lists "
                "the plotted genes. Advanced JSON arguments override visible native defaults."
            ),
            inputs=[
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                TableResultType.Input("table", optional=True),
                io.String.Input("genes_json", default="[]"),
                io.Int.Input("top_n", default=6, min=1),
                io.String.Input("color_by", default="State"),
                io.String.Input("trend_formula", default="~ sm.ns(Pseudotime, df=3)"),
                io.Boolean.Input("relative_expr", default=True),
                io.Int.Input("ncol", default=2, min=1),
                io.Float.Input("width", default=10.0),
                io.Float.Input("height", default=7.0),
                io.Int.Input("dpi", default=150, min=1),
                io.String.Input("extra_parameters_json", default="{}", advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellMonocle2Export(io.ComfyNode):
    fingerprint_inputs = classmethod(_analysis_fingerprint)

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMonocle2Export",
            display_name="Monocle 2 Export to AnnData",
            category=CATEGORY,
            description=(
                "Attach DDRTree coordinates, Pseudotime and State to the selected AnnData by exact cell identity. "
                "Keep its expression, feature axis, Raw snapshot, layers and unrelated metadata intact."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                Monocle2CDSType.Input("cds"),
                io.String.Input("r_runtime", force_input=True),
                io.String.Input("embedding_key", default="X_monocle2"),
                io.String.Input("pseudotime_key", default="monocle2_pseudotime"),
                io.String.Input("state_key", default="monocle2_state"),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"), TableResultType.Output(display_name="table"),
            ),
        )


MONOCLE2_NODE_CLASSES = [
    OpenBioSingleCellMonocle2Runtime,
    OpenBioSingleCellMonocle2Prepare,
    OpenBioSingleCellMonocle2OrderingGenes,
    OpenBioSingleCellMonocle2DDRTree,
    OpenBioSingleCellMonocle2OrderCells,
    OpenBioSingleCellMonocle2TrajectoryPlot,
    OpenBioSingleCellMonocle2DifferentialTest,
    OpenBioSingleCellMonocle2BEAM,
    OpenBioSingleCellMonocle2GeneTrends,
    OpenBioSingleCellMonocle2Export,
]
