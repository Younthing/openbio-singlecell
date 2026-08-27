from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, SCVIModelType
from .scvi_model import SCVIModel

if TYPE_CHECKING:
    from anndata import AnnData


INTEGRATION_CATEGORY = "openbio/single-cell/batch-integration"
CLUSTERING_CATEGORY = "openbio/single-cell/clustering"


def _parse_resolutions(value: str) -> list[float]:
    resolutions = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        resolution = float(item)
        if resolution <= 0:
            raise ValueError("Leiden resolutions must be greater than zero.")
        if resolution not in resolutions:
            resolutions.append(resolution)
    if not resolutions:
        raise ValueError("At least one Leiden resolution is required.")
    return resolutions


def _resolution_key(prefix: str, resolution: float) -> str:
    suffix = f"{resolution:g}".replace("-", "m").replace(".", "_")
    return f"{prefix}_{suffix}"


def _comma_separated_keys(value: str) -> list[str]:
    return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))


class OpenBioSingleCellHarmonyIntegration(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHarmonyIntegration",
            display_name="Harmony Integration",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("batch_key", default="batch"),
                io.String.Input("basis", default="X_pca"),
                io.String.Input("adjusted_basis", default="X_pca_harmony", advanced=True),
                io.Int.Input("max_iter_harmony", default=10, min=1, max=1000, advanced=True),
                io.Int.Input("max_iter_kmeans", default=100, min=1, max=1000, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        batch_key: str = "batch",
        basis: str = "X_pca",
        adjusted_basis: str = "X_pca_harmony",
        max_iter_harmony: int = 10,
        max_iter_kmeans: int = 100,
    ) -> io.NodeOutput:
        if batch_key not in adata.obs:
            raise ValueError(f"Harmony batch column not found in obs: {batch_key!r}")
        if basis not in adata.obsm:
            raise ValueError(f"Harmony basis not found in obsm: {basis!r}")
        try:
            import scanpy.external as sce
        except ImportError as error:
            raise RuntimeError("Harmony Integration requires scanpy external dependencies and harmonypy.") from error

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        sce.pp.harmony_integrate(
            output,
            key=batch_key,
            basis=basis,
            adjusted_basis=adjusted_basis,
            max_iter_harmony=max_iter_harmony,
            max_iter_kmeans=max_iter_kmeans,
        )
        parameters = {
            "batch_key": batch_key,
            "basis": basis,
            "adjusted_basis": adjusted_basis,
            "max_iter_harmony": max_iter_harmony,
            "max_iter_kmeans": max_iter_kmeans,
        }
        finish_adata(output, "harmony_integration", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellSCVIIntegration(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="scVI counts source",
        default="layer",
        layer_input_id="counts_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVIIntegration",
            display_name="scVI Integration",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("batch_key", default=""),
                io.String.Input("size_factor_key", default="", advanced=True),
                io.String.Input("categorical_covariates", default=""),
                io.String.Input("continuous_covariates", default=""),
                io.Int.Input("n_latent", default=10, min=1, max=4096),
                io.Combo.Input("gene_likelihood", options=["zinb", "nb", "poisson"], default="zinb"),
                io.Int.Input("n_layers", default=2, min=1, max=20, advanced=True),
                io.Combo.Input(
                    "dispersion",
                    options=["gene-batch", "gene", "gene-label", "gene-cell"],
                    default="gene-batch",
                    advanced=True,
                ),
                io.Float.Input("dropout_rate", default=0.3, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("max_epochs", default=0, min=0, max=100000, advanced=True),
                io.Boolean.Input("early_stopping", default=True, advanced=True),
                io.Boolean.Input("compute_mde", default=False, advanced=True),
                io.Boolean.Input("store_latent_distribution", default=False, advanced=True),
                io.String.Input("output_key", default="X_scVI", advanced=True),
                io.String.Input("qzm_key", default="X_latent_qzm", advanced=True),
                io.String.Input("qzv_key", default="X_latent_qzv", advanced=True),
                io.String.Input("mde_key", default="X_mde", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SCVIModelType.Output(display_name="model"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        batch_key: str = "",
        size_factor_key: str = "",
        categorical_covariates: str = "",
        continuous_covariates: str = "",
        n_latent: int = 10,
        gene_likelihood: str = "zinb",
        n_layers: int = 2,
        dispersion: str = "gene-batch",
        dropout_rate: float = 0.3,
        max_epochs: int = 0,
        early_stopping: bool = True,
        compute_mde: bool = False,
        store_latent_distribution: bool = False,
        output_key: str = "X_scVI",
        qzm_key: str = "X_latent_qzm",
        qzv_key: str = "X_latent_qzv",
        mde_key: str = "X_mde",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        batch_key = batch_key.strip()
        size_factor_key = size_factor_key.strip()
        output_key = output_key.strip()
        qzm_key = qzm_key.strip()
        qzv_key = qzv_key.strip()
        mde_key = mde_key.strip()
        categorical_keys = _comma_separated_keys(categorical_covariates)
        continuous_keys = _comma_separated_keys(continuous_covariates)
        obs_keys = [batch_key, size_factor_key, *categorical_keys, *continuous_keys]
        missing_keys = [key for key in obs_keys if key and key not in adata.obs]
        if missing_keys:
            raise ValueError(f"scVI observation columns not found: {missing_keys}")
        if not output_key:
            raise ValueError("scVI output_key cannot be empty.")
        if store_latent_distribution and (not qzm_key or not qzv_key):
            raise ValueError("scVI qzm_key and qzv_key cannot be empty when storing the latent distribution.")
        if compute_mde and not mde_key:
            raise ValueError("scVI mde_key cannot be empty when computing MDE.")

        try:
            import scvi
            import torch
        except (ImportError, OSError) as error:
            raise RuntimeError("scVI Integration requires the scvi-tools package.") from error

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        # ComfyUI executes nodes under torch.inference_mode(), but scVI creates and
        # trains a gradient-based model. Constructing the model inside inference
        # mode also turns its parameters into inference tensors, so the entire
        # scVI lifecycle must opt out locally.
        with torch.inference_mode(False):
            scvi.settings.seed = random_seed
            scvi.model.SCVI.setup_anndata(
                output,
                layer=expression.scanpy_layer,
                batch_key=batch_key or None,
                size_factor_key=size_factor_key or None,
                categorical_covariate_keys=categorical_keys or None,
                continuous_covariate_keys=continuous_keys or None,
            )
            model = scvi.model.SCVI(
                output,
                n_layers=n_layers,
                n_latent=n_latent,
                gene_likelihood=gene_likelihood,
                dispersion=dispersion,
                dropout_rate=dropout_rate,
            )
            train_kwargs = {"early_stopping": early_stopping}
            if max_epochs:
                train_kwargs["max_epochs"] = max_epochs
            model.train(**train_kwargs)
            output.obsm[output_key] = model.get_latent_representation()
            if store_latent_distribution or compute_mde:
                qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
                if store_latent_distribution:
                    output.obsm[qzm_key] = qzm
                    output.obsm[qzv_key] = qzv
                if compute_mde:
                    output.obsm[mde_key] = scvi.model.utils.mde(qzm)

        parameters = {
            **expression.parameters(),
            "batch_key": batch_key,
            "size_factor_key": size_factor_key,
            "categorical_covariates": categorical_keys,
            "continuous_covariates": continuous_keys,
            "n_latent": n_latent,
            "gene_likelihood": gene_likelihood,
            "n_layers": n_layers,
            "dispersion": dispersion,
            "dropout_rate": dropout_rate,
            "max_epochs": max_epochs,
            "early_stopping": early_stopping,
            "compute_mde": compute_mde,
            "store_latent_distribution": store_latent_distribution,
            "output_key": output_key,
            "qzm_key": qzm_key if store_latent_distribution else None,
            "qzv_key": qzv_key if store_latent_distribution else None,
            "mde_key": mde_key if compute_mde else None,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "scvi_integration",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        trained_model = SCVIModel(model, output, parameters)
        return io.NodeOutput(output, trained_model)


class OpenBioSingleCellLeidenResolutionSweep(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLeidenResolutionSweep",
            display_name="Leiden Resolution Sweep",
            category=CLUSTERING_CATEGORY,
            description="Run Leiden clustering for a comma-separated list of resolutions.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("resolutions", default="0.25,0.5,1.0,2.0"),
                io.String.Input("key_prefix", default="leiden", advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resolutions: str = "0.25,0.5,1.0,2.0",
        key_prefix: str = "leiden",
        neighbors_key: str = "neighbors",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not key_prefix.strip():
            raise ValueError("Leiden key_prefix cannot be empty.")
        values = _parse_resolutions(resolutions)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        keys = []
        for resolution in values:
            key = _resolution_key(key_prefix.strip(), resolution)
            science.sc.tl.leiden(
                output,
                resolution=resolution,
                key_added=key,
                neighbors_key=neighbors_key,
                random_state=random_seed,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            keys.append(key)

        parameters = {
            "resolutions": values,
            "keys": keys,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "leiden_resolution_sweep",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


INTEGRATION_NODE_CLASSES = [
    OpenBioSingleCellHarmonyIntegration,
    OpenBioSingleCellSCVIIntegration,
    OpenBioSingleCellLeidenResolutionSweep,
]
