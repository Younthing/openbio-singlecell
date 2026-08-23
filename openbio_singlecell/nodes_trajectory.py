from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/trajectory"


def _parse_gene_list(value: str, name: str) -> list[str]:
    genes = list(dict.fromkeys(gene.strip() for gene in value.split(",") if gene.strip()))
    if not genes:
        raise ValueError(f"{name} must contain at least one comma-separated gene.")
    return genes


class OpenBioSingleCellCellCycleScore(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellCycleScore",
            display_name="Cell Cycle Score",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("s_genes", default="MCM5,PCNA,TYMS,FEN1,MCM2,MCM4,RRM1"),
                io.String.Input("g2m_genes", default="HMGB2,CDK1,NUSAP1,UBE2C,BIRC5,TPX2,TOP2A"),
                io.Boolean.Input("use_raw", default=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        s_genes: str,
        g2m_genes: str,
        use_raw: bool = True,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        parsed_s_genes = _parse_gene_list(s_genes, "s_genes")
        parsed_g2m_genes = _parse_gene_list(g2m_genes, "g2m_genes")
        if use_raw and adata.raw is None:
            raise ValueError("Cell Cycle Score with use_raw enabled requires adata.raw.")

        available_genes = set(adata.raw.var_names if use_raw else adata.var_names)
        if not available_genes.intersection(parsed_s_genes):
            raise ValueError("None of the requested S-phase genes are present in the selected expression source.")
        if not available_genes.intersection(parsed_g2m_genes):
            raise ValueError("None of the requested G2/M-phase genes are present in the selected expression source.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.score_genes_cell_cycle(
            output,
            s_genes=parsed_s_genes,
            g2m_genes=parsed_g2m_genes,
            use_raw=use_raw,
            random_state=random_seed,
        )
        parameters = {
            "s_genes": parsed_s_genes,
            "g2m_genes": parsed_g2m_genes,
            "use_raw": use_raw,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "cell_cycle_score",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellDiffusionMap(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDiffusionMap",
            display_name="Diffusion Map",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_comps", default=15, min=2, max=4096),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_comps: int = 15,
        neighbors_key: str = "neighbors",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if neighbors_key not in adata.uns:
            raise ValueError(f"Diffusion Map requires uns[{neighbors_key!r}]; run Neighbors first.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.diffmap(
            output,
            n_comps=n_comps,
            neighbors_key=neighbors_key,
            random_state=random_seed,
        )
        parameters = {
            "n_comps": n_comps,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "diffusion_map",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellPAGA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPAGA",
            display_name="PAGA",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "leiden",
        neighbors_key: str = "neighbors",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not groupby:
            raise ValueError("PAGA groupby cannot be empty.")
        if groupby not in adata.obs:
            raise ValueError(f"PAGA groupby column not found in obs: {groupby!r}")
        if neighbors_key not in adata.uns:
            raise ValueError(f"PAGA requires uns[{neighbors_key!r}]; run Neighbors first.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.paga(output, groups=groupby, neighbors_key=neighbors_key)
        parameters = {"groupby": groupby, "neighbors_key": neighbors_key}
        finish_adata(output, "paga", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellDPT(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDPT",
            display_name="Diffusion Pseudotime",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("root_column", default="leiden"),
                io.String.Input("root_value", default=""),
                io.Int.Input("n_dcs", default=10, min=1, max=4096, advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        root_column: str = "leiden",
        root_value: str = "",
        n_dcs: int = 10,
        neighbors_key: str = "neighbors",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not root_column:
            raise ValueError("DPT root_column cannot be empty.")
        if root_column not in adata.obs:
            raise ValueError(f"DPT root column not found in obs: {root_column!r}")
        root_value = root_value.strip()
        if not root_value:
            raise ValueError("DPT root_value cannot be empty.")
        if neighbors_key not in adata.uns:
            raise ValueError(f"DPT requires uns[{neighbors_key!r}]; run Neighbors first.")

        matches = science.np.flatnonzero(adata.obs[root_column].astype(str).to_numpy() == root_value)
        if matches.size == 0:
            raise ValueError(f"DPT root value {root_value!r} was not found in obs[{root_column!r}].")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        output.uns["iroot"] = int(matches[0])
        science.sc.tl.dpt(output, n_dcs=n_dcs, neighbors_key=neighbors_key)
        parameters = {
            "root_column": root_column,
            "root_value": root_value,
            "n_dcs": n_dcs,
            "neighbors_key": neighbors_key,
        }
        finish_adata(output, "dpt", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


TRAJECTORY_NODE_CLASSES = [
    OpenBioSingleCellCellCycleScore,
    OpenBioSingleCellDiffusionMap,
    OpenBioSingleCellPAGA,
    OpenBioSingleCellDPT,
]
