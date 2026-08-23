from __future__ import annotations

import importlib
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_result, matrix_totals_and_nonzero
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/regulatory"
EXPRESSION_SOURCES = ["X", "raw", "layer"]
RANKING_METHODS = ["wilcoxon", "t-test_overestim_var"]
PYSCENIC_GRN_METHODS = ["grnboost2", "genie3"]


def _require_optional_dependency(name: str) -> Any:
    try:
        return importlib.import_module(name)
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"The {name!r} package is required for this regulatory node ({exc}).") from exc


def _required_name(value: str, description: str) -> str:
    name = value.strip()
    if not name:
        raise ValueError(f"{description} cannot be empty.")
    return name


def _expression_adata(adata: AnnData, source: str, layer_name: str) -> AnnData:
    if source == "raw":
        if adata.raw is None:
            raise ValueError("Expression source 'raw' was selected, but adata.raw is unavailable.")
        work = adata.raw.to_adata()
        work.obs = adata.obs.copy()
        return work

    work = adata.copy()
    if source == "layer":
        layer_name = _required_name(layer_name, "Expression layer")
        if layer_name not in adata.layers:
            raise ValueError(f"Expression layer not found: {layer_name!r}")
        work.X = adata.layers[layer_name].copy()
    return work


def _resource_names(value: str, description: str) -> tuple[str, ...]:
    names = tuple(name.strip() for line in value.splitlines() for name in line.split(",") if name.strip())
    if not names:
        raise ValueError(f"{description} cannot be empty.")
    return names


def _resolve_pyscenic_resources(
    tf_list_file: str,
    ranking_database_files: str,
    motif_annotations_file: str,
) -> tuple[str, tuple[str, ...], str]:
    tf_path = resolve_input_path(tf_list_file, extensions=(".txt",))
    ranking_paths = tuple(
        resolve_input_path(name, extensions=(".feather",))
        for name in _resource_names(ranking_database_files, "Ranking database files")
    )
    motif_path = resolve_input_path(motif_annotations_file, extensions=(".tbl",))
    return tf_path, ranking_paths, motif_path


def _run_pyscenic_cli(arguments: list[str], working_directory: str) -> None:
    command = [sys.executable, "-m", "pyscenic.cli.pyscenic", *arguments]
    completed = subprocess.run(
        command,
        cwd=working_directory,
        check=False,
        capture_output=True,
        text=True,
        errors="replace",
    )
    if completed.returncode == 0:
        return

    details = completed.stderr.strip() or completed.stdout.strip() or "pySCENIC produced no diagnostic output."
    raise RuntimeError(f"pySCENIC {arguments[0]} failed (exit {completed.returncode}):\n{details[-4000:]}")


def _attach_pyscenic_results(
    adata: AnnData,
    loom_path: str,
    regulons_path: str,
    activity_key: str,
    science: dependencies.ScientificDependencies,
) -> None:
    loompy = _require_optional_dependency("loompy")
    scenic_utils = _require_optional_dependency("pyscenic.utils")
    scenic_transform = _require_optional_dependency("pyscenic.transform")
    scenic_export = _require_optional_dependency("pyscenic.export")

    connection = loompy.connect(loom_path, mode="r", validate=False)
    try:
        auc_matrix = science.pd.DataFrame(
            connection.ca.RegulonsAUC,
            index=science.pd.Index(connection.ca.CellID.astype(str)),
        )
    finally:
        connection.close()

    observation_ids = science.pd.Index(adata.obs_names.astype(str))
    if not observation_ids.isin(auc_matrix.index).all():
        missing_count = int((~observation_ids.isin(auc_matrix.index)).sum())
        raise ValueError(f"pySCENIC output is missing {missing_count} AnnData observations.")
    auc_matrix = auc_matrix.reindex(observation_ids)
    auc_matrix.index = adata.obs_names

    motifs = scenic_utils.load_motifs(regulons_path)
    regulons = scenic_transform.df2regulons(motifs)
    scenic_export.add_scenic_metadata(adata, auc_matrix, regulons)
    adata.obsm[activity_key] = auc_matrix.copy()


class OpenBioSingleCellCollecTRIULM(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCollecTRIULM",
            display_name="CollecTRI ULM Activities",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("organism", options=["human", "mouse"], default="human"),
                io.Boolean.Input("split_complexes", default=False),
                io.Combo.Input("source", options=EXPRESSION_SOURCES, default="X"),
                io.String.Input("layer_name", default="log1p_norm"),
                io.String.Input("activity_key", default="collectri_ulm_estimate", advanced=True),
                io.String.Input("pvalue_key", default="collectri_ulm_pvals", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        organism: str = "human",
        split_complexes: bool = False,
        source: str = "X",
        layer_name: str = "log1p_norm",
        activity_key: str = "collectri_ulm_estimate",
        pvalue_key: str = "collectri_ulm_pvals",
    ) -> io.NodeOutput:
        activity_key = _required_name(activity_key, "TF activity output key")
        pvalue_key = _required_name(pvalue_key, "TF activity p-value output key")
        if activity_key == pvalue_key:
            raise ValueError("TF activity and p-value output keys must be different.")

        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = _expression_adata(output, source, layer_name)
        network = decoupler.get_collectri(organism=organism, split_complexes=split_complexes)
        decoupler.run_ulm(
            mat=work,
            net=network,
            source="source",
            target="target",
            weight="weight",
            verbose=False,
        )
        output.obsm[activity_key] = work.obsm["ulm_estimate"].copy()
        output.obsm[pvalue_key] = work.obsm["ulm_pvals"].copy()

        parameters = {
            "organism": organism,
            "split_complexes": split_complexes,
            "source": source,
            "layer_name": layer_name,
            "activity_key": activity_key,
            "pvalue_key": pvalue_key,
        }
        finish_adata(output, "collectri_ulm", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellRankTFActivities(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRankTFActivities",
            display_name="Rank TF Activities",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.String.Input("reference", default="rest"),
                io.Combo.Input("method", options=RANKING_METHODS, default="wilcoxon"),
                io.Int.Input("top_n", default=3, min=1, max=2**31 - 1),
                io.Float.Input("max_pvalue", default=0.05, min=0.0, max=1.0, step=0.01),
                io.String.Input("activity_key", default="collectri_ulm_estimate", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        reference: str = "rest",
        method: str = "wilcoxon",
        top_n: int = 3,
        max_pvalue: float = 0.05,
        activity_key: str = "collectri_ulm_estimate",
    ) -> io.NodeOutput:
        groupby = _required_name(groupby, "TF activity groupby column")
        reference = _required_name(reference, "TF activity reference group")
        activity_key = _required_name(activity_key, "TF activity key")
        if groupby not in adata.obs:
            raise ValueError(f"TF activity groupby column not found in obs: {groupby!r}")
        if activity_key not in adata.obsm:
            raise ValueError(f"TF activities not found in obsm: {activity_key!r}")

        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        activities = decoupler.get_acts(adata, obsm_key=activity_key)
        ranked = decoupler.rank_sources_groups(
            activities,
            groupby=groupby,
            reference=reference,
            method=method,
        )
        ranked = ranked.dropna(subset=["group", "names", "pvals"]).copy()
        ranked = ranked[ranked["pvals"] < max_pvalue]
        ranked = (
            ranked.sort_values(["group", "pvals"])
            .groupby("group", sort=False, observed=True)
            .head(top_n)
            .reset_index(drop=True)
        )

        parameters = {
            "groupby": groupby,
            "reference": reference,
            "method": method,
            "top_n": top_n,
            "max_pvalue": max_pvalue,
            "activity_key": activity_key,
        }
        result = make_result(
            kind="table",
            title=f"TF activities by {groupby}",
            operation="rank_tf_activities",
            parameters=parameters,
            description="Top transcription-factor activities for each observation group.",
            warnings=[] if not ranked.empty else ["No TF activities passed the requested p-value threshold."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=ranked,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellRunPySCENIC(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRunPySCENIC",
            display_name="Run pySCENIC",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("source", options=EXPRESSION_SOURCES, default="X"),
                io.String.Input("layer_name", default="log1p_norm"),
                io.String.Input(
                    "tf_list_file",
                    display_name="TF list",
                    default="",
                    placeholder="openbio-singlecell/allTFs_hg38.txt",
                    tooltip="TF list under the ComfyUI input directory.",
                ),
                io.String.Input(
                    "ranking_database_files",
                    display_name="Ranking databases",
                    default="",
                    multiline=True,
                    placeholder="One .feather file per line",
                    tooltip="One or more cisTarget ranking databases under the ComfyUI input directory.",
                ),
                io.String.Input(
                    "motif_annotations_file",
                    display_name="Motif annotations",
                    default="",
                    placeholder="openbio-singlecell/motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
                    tooltip="Motif-to-TF annotations under the ComfyUI input directory.",
                ),
                io.Combo.Input("grn_method", options=PYSCENIC_GRN_METHODS, default="grnboost2"),
                io.Boolean.Input("mask_dropouts", default=True),
                io.Float.Input("auc_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Int.Input("num_workers", default=4, min=1, max=1024, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
                io.String.Input("activity_key", default="scenic_auc", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(
        cls,
        tf_list_file: str,
        ranking_database_files: str,
        motif_annotations_file: str,
        **kwargs: Any,
    ) -> bool | str:
        try:
            _resolve_pyscenic_resources(tf_list_file, ranking_database_files, motif_annotations_file)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(
        cls,
        tf_list_file: str,
        ranking_database_files: str,
        motif_annotations_file: str,
        **kwargs: Any,
    ) -> Any:
        ranking_files = _resource_names(ranking_database_files, "Ranking database files")
        return (
            input_file_fingerprint(tf_list_file, (".txt",)),
            tuple(input_file_fingerprint(name, (".feather",)) for name in ranking_files),
            input_file_fingerprint(motif_annotations_file, (".tbl",)),
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: str = "X",
        layer_name: str = "log1p_norm",
        tf_list_file: str = "",
        ranking_database_files: str = "",
        motif_annotations_file: str = "",
        grn_method: str = "grnboost2",
        mask_dropouts: bool = True,
        auc_threshold: float = 0.05,
        num_workers: int = 4,
        random_seed: int = 0,
        activity_key: str = "scenic_auc",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        activity_key = _required_name(activity_key, "pySCENIC activity key")
        tf_path, ranking_paths, motif_path = _resolve_pyscenic_resources(
            tf_list_file,
            ranking_database_files,
            motif_annotations_file,
        )
        _require_optional_dependency("pyscenic.cli.pyscenic")
        loompy = _require_optional_dependency("loompy")

        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError("pySCENIC requires a non-empty AnnData object.")
        if not adata.obs_names.is_unique or not adata.var_names.is_unique:
            raise ValueError("pySCENIC requires unique observation and variable names.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = _expression_adata(output, source, layer_name)
        if work.X is None:
            raise ValueError("The selected pySCENIC expression source has no matrix.")

        tf_names = {
            line.strip()
            for line in Path(tf_path).read_text(encoding="utf-8-sig").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        if not tf_names:
            raise ValueError("The pySCENIC TF list is empty.")
        if not work.var_names.isin(tf_names).any():
            raise ValueError("No genes in the selected expression source occur in the pySCENIC TF list.")

        total_counts, detected_genes = matrix_totals_and_nonzero(work.X, axis=1)
        row_attributes = {"Gene": science.np.asarray(work.var_names.astype(str))}
        column_attributes = {
            "CellID": science.np.asarray(work.obs_names.astype(str)),
            "nGene": science.np.asarray(detected_genes).ravel(),
            "nUMI": science.np.asarray(total_counts).ravel(),
        }

        with tempfile.TemporaryDirectory(prefix="openbio_pyscenic_") as temporary_directory:
            temporary_path = Path(temporary_directory)
            input_loom = temporary_path / "expression.loom"
            adjacency_csv = temporary_path / "adjacencies.csv"
            regulons_csv = temporary_path / "regulons.csv"
            final_loom = temporary_path / "aucell.loom"
            loompy.create(str(input_loom), work.X.transpose(), row_attributes, column_attributes)

            _run_pyscenic_cli(
                [
                    "grn",
                    str(input_loom),
                    tf_path,
                    "-o",
                    str(adjacency_csv),
                    "--method",
                    grn_method,
                    "--num_workers",
                    str(num_workers),
                    "--seed",
                    str(random_seed),
                ],
                temporary_directory,
            )

            ctx_arguments = [
                "ctx",
                str(adjacency_csv),
                *ranking_paths,
                "--annotations_fname",
                motif_path,
                "--expression_mtx_fname",
                str(input_loom),
                "--output",
                str(regulons_csv),
                "--num_workers",
                str(num_workers),
            ]
            if mask_dropouts:
                ctx_arguments.append("--mask_dropouts")
            _run_pyscenic_cli(ctx_arguments, temporary_directory)

            _run_pyscenic_cli(
                [
                    "aucell",
                    str(input_loom),
                    str(regulons_csv),
                    "--auc_threshold",
                    str(auc_threshold),
                    "--output",
                    str(final_loom),
                    "--num_workers",
                    str(num_workers),
                    "--seed",
                    str(random_seed),
                ],
                temporary_directory,
            )

            _attach_pyscenic_results(output, str(final_loom), str(regulons_csv), activity_key, science)

        parameters = {
            "source": source,
            "layer_name": layer_name,
            "tf_list_file": tf_list_file,
            "ranking_database_files": list(_resource_names(ranking_database_files, "Ranking database files")),
            "motif_annotations_file": motif_annotations_file,
            "grn_method": grn_method,
            "mask_dropouts": mask_dropouts,
            "auc_threshold": auc_threshold,
            "num_workers": num_workers,
            "activity_key": activity_key,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "run_pyscenic",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellImportPySCENICResults(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellImportPySCENICResults",
            display_name="Import pySCENIC Results",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("final_loom_file", default="openbio-singlecell/pyscenic_output.loom"),
                io.String.Input("regulons_csv", default="openbio-singlecell/regulons.csv"),
                io.String.Input("activity_key", default="scenic_auc", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, final_loom_file: str, regulons_csv: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(final_loom_file, extensions=(".loom",))
            resolve_input_path(regulons_csv, extensions=(".csv",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, final_loom_file: str, regulons_csv: str, **kwargs: Any) -> Any:
        return (
            input_file_fingerprint(final_loom_file, (".loom",)),
            input_file_fingerprint(regulons_csv, (".csv",)),
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        final_loom_file: str = "openbio-singlecell/pyscenic_output.loom",
        regulons_csv: str = "openbio-singlecell/regulons.csv",
        activity_key: str = "scenic_auc",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        activity_key = _required_name(activity_key, "pySCENIC activity key")
        loom_path = resolve_input_path(final_loom_file, extensions=(".loom",))
        regulons_path = resolve_input_path(regulons_csv, extensions=(".csv",))
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        _attach_pyscenic_results(output, loom_path, regulons_path, activity_key, science)

        parameters = {
            "final_loom_file": final_loom_file,
            "regulons_csv": regulons_csv,
            "activity_key": activity_key,
        }
        finish_adata(output, "import_pyscenic_results", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellSCENICRegulonSpecificity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICRegulonSpecificity",
            display_name="SCENIC Regulon Specificity",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.String.Input("activity_key", default="scenic_auc", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        activity_key: str = "scenic_auc",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        groupby = _required_name(groupby, "SCENIC groupby column")
        activity_key = _required_name(activity_key, "SCENIC activity key")
        if groupby not in adata.obs:
            raise ValueError(f"SCENIC groupby column not found in obs: {groupby!r}")
        if activity_key not in adata.obsm:
            raise ValueError(f"SCENIC activities not found in obsm: {activity_key!r}")
        auc_matrix = adata.obsm[activity_key]
        if not isinstance(auc_matrix, science.pd.DataFrame):
            raise ValueError(f"obsm[{activity_key!r}] must be a named regulon activity DataFrame.")

        scenic_rss = _require_optional_dependency("pyscenic.rss")
        started_at = time.perf_counter()
        rss = scenic_rss.regulon_specificity_scores(auc_matrix, adata.obs[groupby])
        rss.index.name = "regulon"
        table = rss.reset_index().melt(id_vars="regulon", var_name="group", value_name="rss")
        table = table.sort_values(["group", "rss"], ascending=[True, False]).reset_index(drop=True)
        parameters = {"groupby": groupby, "activity_key": activity_key}
        result = make_result(
            kind="table",
            title=f"SCENIC regulon specificity by {groupby}",
            operation="scenic_regulon_specificity",
            parameters=parameters,
            description="Regulon specificity scores calculated from the pySCENIC AUCell matrix.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellSCENICActivityBinarization(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICActivityBinarization",
            display_name="SCENIC Activity Binarization",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("activity_key", default="scenic_auc", advanced=True),
                io.String.Input("binary_key", default="scenic_binary", advanced=True),
                io.String.Input("threshold_key", default="scenic_thresholds", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        activity_key: str = "scenic_auc",
        binary_key: str = "scenic_binary",
        threshold_key: str = "scenic_thresholds",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        activity_key = _required_name(activity_key, "SCENIC activity key")
        binary_key = _required_name(binary_key, "SCENIC binary output key")
        threshold_key = _required_name(threshold_key, "SCENIC threshold output key")
        if activity_key not in adata.obsm:
            raise ValueError(f"SCENIC activities not found in obsm: {activity_key!r}")
        activities = adata.obsm[activity_key]
        if not isinstance(activities, science.pd.DataFrame):
            raise ValueError(f"obsm[{activity_key!r}] must be a named regulon activity DataFrame.")

        scenic_binarization = _require_optional_dependency("pyscenic.binarization")
        started_at = time.perf_counter()
        binary, thresholds = scenic_binarization.binarize(activities)
        output = adata.copy()
        output.obsm[binary_key] = science.pd.DataFrame(binary).reindex(output.obs_names)
        threshold_series = science.pd.Series(thresholds, dtype=float)
        output.uns[threshold_key] = {str(name): float(value) for name, value in threshold_series.items()}
        parameters = {
            "activity_key": activity_key,
            "binary_key": binary_key,
            "threshold_key": threshold_key,
        }
        finish_adata(
            output,
            "scenic_activity_binarization",
            parameters,
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellSCENICTFModules(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICTFModules",
            display_name="SCENIC TF Modules",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("adjacency_csv", default="openbio-singlecell/adjacencies.csv"),
                io.String.Input("transcription_factor", default=""),
                io.Combo.Input("source", options=EXPRESSION_SOURCES, default="X"),
                io.String.Input("layer_name", default="log1p_norm"),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def validate_inputs(cls, adjacency_csv: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(adjacency_csv, extensions=(".csv",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, adjacency_csv: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(adjacency_csv, (".csv",))

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        adjacency_csv: str = "openbio-singlecell/adjacencies.csv",
        transcription_factor: str = "",
        source: str = "X",
        layer_name: str = "log1p_norm",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        path = resolve_input_path(adjacency_csv, extensions=(".csv",))
        adjacency = science.pd.read_csv(path)
        if adjacency.empty:
            raise ValueError("SCENIC adjacency file is empty.")
        work = _expression_adata(adata, source, layer_name)
        matrix = work.X.toarray() if science.sparse.issparse(work.X) else science.np.asarray(work.X)
        expression = science.pd.DataFrame(matrix, index=work.obs_names, columns=work.var_names)

        scenic_utils = _require_optional_dependency("pyscenic.utils")
        started_at = time.perf_counter()
        modules = list(scenic_utils.modules_from_adjacencies(adjacency, expression))
        transcription_factor = transcription_factor.strip()
        if transcription_factor:
            modules = [module for module in modules if str(module.transcription_factor) == transcription_factor]
        rows = []
        for index, module in enumerate(modules):
            rows.extend(
                {
                    "transcription_factor": str(module.transcription_factor),
                    "module": index,
                    "gene": str(gene),
                }
                for gene in module.genes
            )
        table = science.pd.DataFrame.from_records(
            rows,
            columns=["transcription_factor", "module", "gene"],
        )
        warnings = [] if not table.empty else ["No SCENIC modules matched the requested transcription factor."]
        parameters = {
            "adjacency_csv": adjacency_csv,
            "transcription_factor": transcription_factor,
            "source": source,
            "layer_name": layer_name,
        }
        result = make_result(
            kind="table",
            title=f"SCENIC modules: {transcription_factor or 'all TFs'}",
            operation="scenic_tf_modules",
            parameters=parameters,
            description="TF-module target genes generated from SCENIC adjacencies and the selected expression matrix.",
            warnings=warnings,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


REGULATORY_NODE_CLASSES = [
    OpenBioSingleCellCollecTRIULM,
    OpenBioSingleCellRankTFActivities,
    OpenBioSingleCellRunPySCENIC,
    OpenBioSingleCellImportPySCENICResults,
    OpenBioSingleCellSCENICRegulonSpecificity,
    OpenBioSingleCellSCENICActivityBinarization,
    OpenBioSingleCellSCENICTFModules,
]


__all__ = [
    "REGULATORY_NODE_CLASSES",
    "OpenBioSingleCellCollecTRIULM",
    "OpenBioSingleCellImportPySCENICResults",
    "OpenBioSingleCellRankTFActivities",
    "OpenBioSingleCellRunPySCENIC",
    "OpenBioSingleCellSCENICActivityBinarization",
    "OpenBioSingleCellSCENICRegulonSpecificity",
    "OpenBioSingleCellSCENICTFModules",
]
