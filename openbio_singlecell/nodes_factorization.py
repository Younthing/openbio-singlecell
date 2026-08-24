from __future__ import annotations

import copy
import importlib
import os
import time
from typing import TYPE_CHECKING, Any

import folder_paths
from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .expression_source import (
    DynamicExpressionSource,
    ExpressionSource,
    ExpressionSourceSpec,
)
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/factorization"


def _require_omicverse() -> Any:
    try:
        return importlib.import_module("omicverse")
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"The 'omicverse' package is required for cNMF ({exc}).") from exc


def _expression_adata(adata: AnnData, source: ExpressionSource) -> AnnData:
    if source.kind == "raw":
        work = adata.raw.to_adata()
        work.obs = adata.obs.copy()
        return work

    work = adata.copy()
    if source.kind == "layer":
        work.X = adata.layers[source.layer_name].copy()
    return work


def _working_directory(relative_directory: str) -> str:
    relative_directory = relative_directory.strip()
    if not relative_directory:
        raise ValueError("cNMF working directory cannot be empty.")
    if os.path.isabs(relative_directory) or os.path.splitdrive(relative_directory)[0]:
        raise ValueError("cNMF working directory must be relative to the ComfyUI temp directory.")

    root = folder_paths.get_temp_directory()
    target = os.path.abspath(os.path.join(root, os.path.normpath(relative_directory)))
    if not folder_paths.is_within_directory(root, target):
        raise ValueError("cNMF working directory escapes the ComfyUI temp directory.")
    os.makedirs(target, exist_ok=True)
    return target


def _analysis_name(value: str) -> str:
    name = value.strip()
    if not name or name in {".", ".."}:
        raise ValueError("cNMF analysis name cannot be empty.")
    if os.path.basename(name) != name or os.path.splitdrive(name)[0]:
        raise ValueError("cNMF analysis name must be a single directory name.")
    return name


def _align_frame(frame: Any, index: Any, science: dependencies.ScientificDependencies, description: str) -> Any:
    table = science.pd.DataFrame(frame).copy()
    table.index = table.index.astype(str)
    target_index = science.pd.Index(index.astype(str))
    missing = target_index[~target_index.isin(table.index)]
    if len(missing):
        raise RuntimeError(f"cNMF {description} is missing {len(missing)} AnnData entries.")
    table = table.reindex(target_index)
    table.index = index
    return table


class OpenBioSingleCellCNMF(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="cNMF expression source",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNMF",
            display_name="cNMF Gene Programs",
            category=CATEGORY,
            description="Run the complete OmicVerse cNMF workflow and attach the selected consensus programs.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("selected_k", default=7, min=2, max=1024),
                io.Combo.Input("cluster_assignment", options=["rfc", "usage_argmax"], default="rfc"),
                io.String.Input("rfc_use_rep", default="scaled|original|X_pca"),
                io.Float.Input("rfc_threshold", default=0.5, min=0.0, max=1.0, step=0.05),
                io.Int.Input("components_min", default=3, min=2, max=1024, advanced=True),
                io.Int.Input("components_max", default=19, min=2, max=1024, advanced=True),
                io.Int.Input("n_iter", default=200, min=1, max=100000, advanced=True),
                io.Int.Input("num_highvar_genes", default=2000, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input("density_threshold", default=0.1, min=0.0, max=2.0, step=0.05),
                io.Int.Input("n_top_genes", default=100, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("workers", default=1, min=1, max=1024, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
                io.String.Input(
                    "working_directory",
                    default="openbio-singlecell/cnmf",
                    advanced=True,
                ),
                io.String.Input("name", default="cnmf", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        selected_k: int = 7,
        cluster_assignment: str = "rfc",
        rfc_use_rep: str = "scaled|original|X_pca",
        rfc_threshold: float = 0.5,
        components_min: int = 3,
        components_max: int = 19,
        n_iter: int = 200,
        num_highvar_genes: int = 2000,
        density_threshold: float = 0.1,
        n_top_genes: int = 100,
        workers: int = 1,
        random_seed: int = 123,
        working_directory: str = "openbio-singlecell/cnmf",
        name: str = "cnmf",
    ) -> io.NodeOutput:
        if components_min > components_max:
            raise ValueError("cNMF components_min cannot exceed components_max.")
        if not components_min <= selected_k <= components_max:
            raise ValueError("cNMF selected_k must be within the requested component range.")
        if cluster_assignment not in {"rfc", "usage_argmax"}:
            raise ValueError(f"Unsupported cNMF cluster assignment method: {cluster_assignment!r}")
        rfc_use_rep = rfc_use_rep.strip()
        if cluster_assignment == "rfc" and not rfc_use_rep:
            raise ValueError("cNMF RFC representation cannot be empty.")

        science = dependencies.require_scientific_dependencies()
        omicverse = _require_omicverse()
        resolved_working_directory = _working_directory(working_directory)
        name = _analysis_name(name)
        components = science.np.arange(components_min, components_max + 1, dtype=int)

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        work = _expression_adata(adata, expression)
        cnmf = omicverse.single.cNMF(
            work,
            components=components,
            n_iter=n_iter,
            seed=random_seed,
            num_highvar_genes=num_highvar_genes,
            output_dir=resolved_working_directory,
            name=name,
        )
        for worker_i in range(workers):
            cnmf.factorize(worker_i=worker_i, total_workers=workers)
        cnmf.combine(skip_missing_files=True)
        cnmf.consensus(
            k=selected_k,
            density_threshold=density_threshold,
            show_clustering=False,
            close_clustergram_fig=True,
        )
        result_dict = cnmf.load_results(
            K=selected_k,
            density_threshold=density_threshold,
            n_top_genes=n_top_genes,
        )

        obs_columns_before = set(work.obs.columns)
        var_columns_before = set(work.var.columns)
        obsm_keys_before = set(work.obsm)
        uns_keys_before = set(work.uns)
        cnmf.get_results(work, result_dict)
        if cluster_assignment == "rfc":
            cnmf.get_results_rfc(
                work,
                result_dict,
                use_rep=rfc_use_rep,
                cNMF_threshold=rfc_threshold,
            )

        for column in work.obs.columns:
            if column not in obs_columns_before or str(column).lower().startswith(("cnmf", "nmf")):
                output.obs[column] = work.obs[column].reindex(output.obs_names)
        for column in work.var.columns:
            if column not in var_columns_before or str(column).lower().startswith(("cnmf", "nmf")):
                output.var[column] = work.var[column].reindex(output.var_names)
        for key in work.obsm:
            if key not in obsm_keys_before or str(key).lower().startswith(("cnmf", "nmf")):
                value = work.obsm[key]
                output.obsm[key] = value.copy() if hasattr(value, "copy") else copy.deepcopy(value)
        for key in work.uns:
            if key not in uns_keys_before or str(key).lower().startswith(("cnmf", "nmf")):
                output.uns[key] = copy.deepcopy(work.uns[key])

        if "usage_norm" not in result_dict or "gep_scores" not in result_dict:
            raise RuntimeError("cNMF did not return usage_norm and gep_scores results.")
        usage = _align_frame(result_dict["usage_norm"], output.obs_names, science, "cell usage matrix")
        gene_scores = _align_frame(result_dict["gep_scores"], output.var_names, science, "gene score matrix")
        for column in usage.columns:
            output.obs[str(column)] = usage[column].to_numpy()
        if cluster_assignment == "rfc":
            if "cNMF_cluster_rfc" not in work.obs:
                raise RuntimeError("cNMF RFC assignment did not produce obs['cNMF_cluster_rfc'].")
            assigned = work.obs["cNMF_cluster_rfc"].reindex(output.obs_names).astype(str)
            assigned = assigned.map(lambda value: value if value.startswith("cNMF_") else f"cNMF_{value}")
        else:
            assigned = usage.idxmax(axis=1).astype(str)
        output.obs["cNMF_cluster"] = science.pd.Categorical(assigned.to_numpy())
        output.obsm["cNMF_usage"] = usage
        for column in gene_scores.columns:
            output.var[str(column)] = gene_scores[column].to_numpy()

        top_genes = result_dict.get("top_genes")
        if top_genes is None:
            top_gene_records: dict[str, list[Any]] = {}
        else:
            top_gene_records = science.pd.DataFrame(top_genes).to_dict(orient="list")
        output.uns["cnmf"] = {
            "selected_k": int(selected_k),
            "component_range": [int(components_min), int(components_max)],
            "usage_key": "cNMF_usage",
            "cluster_key": "cNMF_cluster",
            "cluster_assignment": cluster_assignment,
            "top_genes": top_gene_records,
        }

        parameters = {
            **expression.parameters(),
            "selected_k": selected_k,
            "cluster_assignment": cluster_assignment,
            "rfc_use_rep": rfc_use_rep,
            "rfc_threshold": rfc_threshold,
            "components_min": components_min,
            "components_max": components_max,
            "n_iter": n_iter,
            "num_highvar_genes": num_highvar_genes,
            "density_threshold": density_threshold,
            "n_top_genes": n_top_genes,
            "workers": workers,
            "random_seed": random_seed,
            "working_directory": working_directory,
            "name": name,
        }
        finish_adata(
            output,
            "cnmf_gene_programs",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


FACTORIZATION_NODE_CLASSES = [OpenBioSingleCellCNMF]


__all__ = ["FACTORIZATION_NODE_CLASSES", "OpenBioSingleCellCNMF"]
