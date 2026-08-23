from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_result
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


ABUNDANCE_CATEGORY = "openbio/single-cell/differential-abundance"


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError(
            "Milo and tascCODA nodes require the pertpy package and the relevant optional extras."
        ) from error
    return pertpy


def _require_schist() -> Any:
    try:
        import schist
    except (ImportError, OSError) as error:
        raise RuntimeError("Schist Nested Model requires the schist package.") from error
    return schist


def _comma_separated_columns(value: str) -> list[str]:
    return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))


def _table_with_index(frame: Any, index_name: str, science: dependencies.ScientificDependencies) -> Any:
    table = science.pd.DataFrame(frame).copy()
    if isinstance(table.index, science.pd.MultiIndex):
        table.index = table.index.set_names(
            [name or f"{index_name}_{position + 1}" for position, name in enumerate(table.index.names)]
        )
    elif table.index.name is None:
        table.index.name = index_name
    return table.reset_index()


class OpenBioSingleCellSchistNestedModel(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSchistNestedModel",
            display_name="Schist Nested Model",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, random_seed: int = 123) -> io.NodeOutput:
        schist = _require_schist()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        schist.inference.nested_model(output, random_seed=random_seed)

        parameters = {"random_seed": random_seed}
        finish_adata(
            output,
            "schist_nested_model",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellMiloDifferentialAbundance(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMiloDifferentialAbundance",
            display_name="Milo Differential Abundance",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("design", default="~ group"),
                io.String.Input("contrast", default=""),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("use_rep", default="X_pca"),
                io.Int.Input("n_neighbors", default=30, min=2, max=2**31 - 1),
                io.Float.Input(
                    "neighborhood_proportion",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "mixed_annotation_threshold",
                    default=0.6,
                    min=0.0,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.String.Input("neighbors_key", default="openbio_milo", advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        design: str = "~ group",
        contrast: str = "",
        annotation_key: str = "cell_type",
        use_rep: str = "X_pca",
        n_neighbors: int = 30,
        neighborhood_proportion: float = 0.1,
        mixed_annotation_threshold: float = 0.6,
        neighbors_key: str = "openbio_milo",
        random_seed: int = 123,
    ) -> io.NodeOutput:
        if sample_key not in adata.obs:
            raise ValueError(f"Milo sample column not found in obs: {sample_key!r}")
        if annotation_key not in adata.obs:
            raise ValueError(f"Milo annotation column not found in obs: {annotation_key!r}")
        if use_rep not in adata.obsm:
            raise ValueError(f"Milo representation not found in obsm: {use_rep!r}")
        if not design.strip():
            raise ValueError("Milo design cannot be empty.")
        if not neighbors_key.strip():
            raise ValueError("Milo neighbors_key cannot be empty.")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        milo = pertpy.tl.Milo()
        mdata = milo.load(adata.copy())
        rna = mdata["rna"]
        science.sc.pp.neighbors(
            rna,
            n_neighbors=n_neighbors,
            use_rep=use_rep,
            key_added=neighbors_key,
            random_state=random_seed,
        )
        milo.make_nhoods(mdata, neighbors_key=neighbors_key, prop=neighborhood_proportion)
        counted = milo.count_nhoods(mdata, sample_col=sample_key)
        if counted is not None:
            mdata = counted

        da_options = {"design": design.strip()}
        if contrast.strip():
            da_options["model_contrasts"] = contrast.strip()
        milo.da_nhoods(mdata, **da_options)
        milo.annotate_nhoods(mdata, anno_col=annotation_key)

        abundance = mdata["milo"].var.copy()
        if "nhood_annotation" in abundance and "nhood_annotation_frac" in abundance:
            mixed = abundance["nhood_annotation_frac"] < mixed_annotation_threshold
            abundance["nhood_annotation"] = abundance["nhood_annotation"].astype(object)
            abundance.loc[mixed, "nhood_annotation"] = "Mixed"
        table = _table_with_index(abundance, "neighborhood", science)

        parameters = {
            "sample_key": sample_key,
            "design": design.strip(),
            "contrast": contrast.strip(),
            "annotation_key": annotation_key,
            "use_rep": use_rep,
            "n_neighbors": n_neighbors,
            "neighborhood_proportion": neighborhood_proportion,
            "mixed_annotation_threshold": mixed_annotation_threshold,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title="Milo differential abundance",
            operation="milo_differential_abundance",
            parameters=parameters,
            description="Differential abundance statistics and majority annotation for each cell neighborhood.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellTasccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellTasccodaDifferentialComposition",
            display_name="tascCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="batch"),
                io.String.Input("annotation_key", default="nsbm_level_1"),
                io.String.Input("covariate_keys", default="group"),
                io.String.Input("formula", default="group"),
                io.String.Input(
                    "hierarchy_levels",
                    default="nsbm_level_4,nsbm_level_3,nsbm_level_2,nsbm_level_1",
                ),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("phi", default=0.0, min=0.0, step=0.1, advanced=True),
                io.Float.Input("lambda_1", default=3.5, min=0.0, step=0.1, advanced=True),
                io.Float.Input("estimated_fdr", default=0.05, min=0.0, max=1.0, step=0.01, advanced=True),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=1234, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "batch",
        annotation_key: str = "nsbm_level_1",
        covariate_keys: str = "group",
        formula: str = "group",
        hierarchy_levels: str = "nsbm_level_4,nsbm_level_3,nsbm_level_2,nsbm_level_1",
        reference_cell_type: str = "automatic",
        phi: float = 0.0,
        lambda_1: float = 3.5,
        estimated_fdr: float = 0.05,
        num_samples: int = 10000,
        num_warmup: int = 1000,
        random_seed: int = 1234,
    ) -> io.NodeOutput:
        covariates = _comma_separated_columns(covariate_keys)
        levels = _comma_separated_columns(hierarchy_levels)
        required_obs = [sample_key, annotation_key, *covariates, *levels]
        missing_obs = [key for key in required_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"tascCODA observation columns not found: {missing_obs}")
        if not covariates:
            raise ValueError("tascCODA requires at least one covariate column.")
        if not levels:
            raise ValueError("tascCODA requires at least one hierarchy level.")
        if not formula.strip():
            raise ValueError("tascCODA formula cannot be empty.")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        model = pertpy.tl.Tasccoda()
        mdata = model.load(
            adata.copy(),
            type="cell_level",
            cell_type_identifier=annotation_key,
            sample_identifier=sample_key,
            covariate_obs=covariates,
            levels_orig=levels,
            add_level_name=True,
        )
        prepared = model.prepare(
            mdata,
            modality_key="coda",
            reference_cell_type=reference_cell_type or "automatic",
            formula=formula.strip(),
            pen_args={"phi": phi, "lambda_1": lambda_1},
            tree_key="tree",
        )
        if prepared is not None:
            mdata = prepared
        model.run_nuts(
            mdata,
            modality_key="coda",
            rng_key=random_seed,
            num_samples=num_samples,
            num_warmup=num_warmup,
        )

        summaries = model.summary_prepare(mdata["coda"], est_fdr=estimated_fdr)
        sections = ("intercept", "effect", "node_effect")
        tables = []
        for section, frame in zip(sections, summaries, strict=False):
            table = _table_with_index(frame, "parameter", science)
            table.insert(0, "summary_section", section)
            tables.append(table)
        table = science.pd.concat(tables, ignore_index=True, sort=False)

        parameters = {
            "sample_key": sample_key,
            "annotation_key": annotation_key,
            "covariate_keys": covariates,
            "formula": formula.strip(),
            "hierarchy_levels": levels,
            "reference_cell_type": reference_cell_type or "automatic",
            "phi": phi,
            "lambda_1": lambda_1,
            "estimated_fdr": estimated_fdr,
            "num_samples": num_samples,
            "num_warmup": num_warmup,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title="tascCODA differential composition",
            operation="tasccoda_differential_composition",
            parameters=parameters,
            description="Intercept, cell-type effect, and hierarchy-node effect summaries from tascCODA.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


ABUNDANCE_NODE_CLASSES = [
    OpenBioSingleCellSchistNestedModel,
    OpenBioSingleCellMiloDifferentialAbundance,
    OpenBioSingleCellTasccodaDifferentialComposition,
]


__all__ = [
    "ABUNDANCE_NODE_CLASSES",
    "OpenBioSingleCellMiloDifferentialAbundance",
    "OpenBioSingleCellSchistNestedModel",
    "OpenBioSingleCellTasccodaDifferentialComposition",
]
