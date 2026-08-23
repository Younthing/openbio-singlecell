from __future__ import annotations

import gzip
import re
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/data"
GTF_EXTENSIONS = (".gtf", ".gtf.gz")
_GTF_ATTRIBUTE_PATTERN = re.compile(r"([^\s;]+)\s+\"([^\"]*)\"")


def _comma_separated_values(values: str) -> tuple[str, ...]:
    parsed = tuple(dict.fromkeys(value.strip() for value in values.split(",") if value.strip()))
    if not parsed:
        raise ValueError("Observation values must contain at least one non-empty value.")
    return parsed


def _stable_gene_id(value: str) -> str:
    return value.split(".", 1)[0]


def _read_gtf_gene_mapping(path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(_GTF_ATTRIBUTE_PATTERN.findall(fields[8]))
            gene_id = attributes.get("gene_id")
            gene_name = attributes.get("gene_name")
            if gene_id and gene_name:
                mapping.setdefault(_stable_gene_id(gene_id), gene_name)
    if not mapping:
        raise ValueError("The GTF file does not contain gene entries with gene_id and gene_name attributes.")
    return mapping


class OpenBioSingleCellUseExpressionLayer(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellUseExpressionLayer",
            display_name="Use Expression Layer",
            category=CATEGORY,
            description="Copy a named AnnData expression layer into X.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("layer_name", default="counts"),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, layer_name: str = "counts") -> io.NodeOutput:
        if not layer_name:
            raise ValueError("Expression layer name cannot be empty.")
        if layer_name not in adata.layers:
            raise ValueError(f"Expression layer not found: {layer_name!r}")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        output.X = output.layers[layer_name].copy()
        finish_adata(output, "use_expression_layer", {"layer_name": layer_name}, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellSnapshotExpression(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSnapshotExpression",
            display_name="Snapshot Expression",
            category=CATEGORY,
            description="Copy the current X matrix into an AnnData layer, raw, or both.",
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("destination", options=["layer", "raw", "layer_and_raw"], default="layer"),
                io.String.Input("layer_name", default="counts", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        destination: str = "layer",
        layer_name: str = "counts",
    ) -> io.NodeOutput:
        if destination not in {"layer", "raw", "layer_and_raw"}:
            raise ValueError(f"Unsupported expression snapshot destination: {destination!r}")
        layer_name = layer_name.strip()
        if destination in {"layer", "layer_and_raw"} and not layer_name:
            raise ValueError("Expression snapshot layer name cannot be empty.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        if destination in {"layer", "layer_and_raw"}:
            output.layers[layer_name] = output.X.copy()
        if destination in {"raw", "layer_and_raw"}:
            output.raw = output.copy()
        parameters = {"destination": destination, "layer_name": layer_name}
        finish_adata(output, "snapshot_expression", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellSubsetObservations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSubsetObservations",
            display_name="Subset Observations",
            category=CATEGORY,
            description="Subset observations by string values in an obs column.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("column", default="sample"),
                io.String.Input("values", default=""),
                io.Boolean.Input("invert", default=False),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        column: str = "sample",
        values: str = "",
        invert: bool = False,
    ) -> io.NodeOutput:
        if column not in adata.obs:
            raise ValueError(f"Observation annotation not found in obs: {column!r}")
        selected_values = _comma_separated_values(values)

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        mask = adata.obs[column].astype(str).isin(selected_values)
        if invert:
            mask = ~mask
        output = adata[mask.to_numpy()].copy()
        parameters = {"column": column, "values": list(selected_values), "invert": invert}
        finish_adata(output, "subset_observations", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellMergeObservationAnnotations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMergeObservationAnnotations",
            display_name="Merge Observation Annotations",
            category=CATEGORY,
            description="Merge an observation annotation from a subset AnnData by obs_names.",
            inputs=[
                AnnDataType.Input("adata"),
                AnnDataType.Input("subset_adata"),
                io.String.Input("source_column", default="cell_type"),
                io.String.Input("target_column", default="cell_type", advanced=True),
                io.String.Input("prefix", default="", advanced=True),
                io.Boolean.Input("overwrite", default=False, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        subset_adata: AnnData,
        source_column: str = "cell_type",
        target_column: str = "cell_type",
        prefix: str = "",
        overwrite: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if source_column not in subset_adata.obs:
            raise ValueError(f"Source observation annotation not found in subset_adata.obs: {source_column!r}")
        if not target_column:
            raise ValueError("Target observation annotation name cannot be empty.")
        if not adata.obs_names.is_unique or not subset_adata.obs_names.is_unique:
            raise ValueError("Merging observation annotations requires unique obs_names in both AnnData inputs.")

        shared_names = adata.obs_names[adata.obs_names.isin(subset_adata.obs_names)]
        if len(shared_names) == 0:
            raise ValueError("The AnnData inputs do not share any obs_names.")

        source_values = subset_adata.obs.loc[shared_names, source_column].astype("object")
        if prefix:
            source_values = source_values.map(lambda value: f"{prefix}{value}" if science.pd.notna(value) else value)
        source_values = source_values[source_values.notna()]
        if source_values.empty:
            raise ValueError("The source annotation has no non-missing values for shared observations.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        if target_column in output.obs:
            target_values = output.obs[target_column].astype("object").copy()
        else:
            target_values = science.pd.Series(science.pd.NA, index=output.obs_names, dtype="object")

        values_to_merge = source_values
        if not overwrite:
            values_to_merge = values_to_merge[target_values.loc[values_to_merge.index].isna()]
        target_values.loc[values_to_merge.index] = values_to_merge
        output.obs[target_column] = target_values

        warnings = []
        if values_to_merge.empty:
            warnings.append("No observation annotations changed because overwrite is disabled.")
        parameters = {
            "source_column": source_column,
            "target_column": target_column,
            "prefix": prefix,
            "overwrite": overwrite,
        }
        finish_adata(
            output,
            "merge_observation_annotations",
            parameters,
            cells,
            genes,
            started_at,
            warnings=warnings,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellMapGeneIdsFromGTF(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMapGeneIdsFromGTF",
            display_name="Map Gene IDs from GTF",
            category=CATEGORY,
            description="Map AnnData var_names from gene IDs to gene names using a GTF file.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gtf_path", default="annotations/genes.gtf.gz"),
                io.String.Input("original_id_column", default="ensembl_id", advanced=True),
                io.Boolean.Input("drop_unmapped", default=False),
                io.Boolean.Input("make_unique", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, gtf_path: str, **kwargs) -> bool | str:
        try:
            resolve_input_path(gtf_path, extensions=GTF_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gtf_path: str, **kwargs) -> Any:
        return input_file_fingerprint(gtf_path, GTF_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gtf_path: str = "annotations/genes.gtf.gz",
        original_id_column: str = "ensembl_id",
        drop_unmapped: bool = False,
        make_unique: bool = True,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not original_id_column:
            raise ValueError("Original gene ID annotation name cannot be empty.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        resolved = resolve_input_path(gtf_path, extensions=GTF_EXTENSIONS)
        mapping = _read_gtf_gene_mapping(resolved)
        original_ids = [str(value) for value in adata.var_names]
        mapped_names = [mapping.get(_stable_gene_id(value)) for value in original_ids]
        mapped_mask = science.np.asarray([value is not None for value in mapped_names], dtype=bool)
        mapped_count = int(mapped_mask.sum())
        if mapped_count == 0:
            raise ValueError("None of the AnnData var_names matched a gene_id in the GTF file.")

        output = adata.copy()
        output.var[original_id_column] = original_ids
        if drop_unmapped:
            output = output[:, mapped_mask].copy()
            output.var_names = science.pd.Index([mapped_names[index] for index, keep in enumerate(mapped_mask) if keep])
        else:
            output.var_names = science.pd.Index(
                [
                    mapped_name if mapped_name is not None else original_id
                    for original_id, mapped_name in zip(original_ids, mapped_names, strict=True)
                ]
            )
        if make_unique:
            output.var_names_make_unique()

        unmapped_count = genes - mapped_count
        warnings = []
        if unmapped_count:
            action = "dropped" if drop_unmapped else "left unchanged"
            warnings.append(f"{unmapped_count} gene IDs were not mapped and were {action}.")
        parameters = {
            "gtf_path": gtf_path,
            "original_id_column": original_id_column,
            "drop_unmapped": drop_unmapped,
            "make_unique": make_unique,
        }
        finish_adata(
            output,
            "map_gene_ids_from_gtf",
            parameters,
            cells,
            genes,
            started_at,
            warnings=warnings,
        )
        return io.NodeOutput(output)


DATA_NODE_CLASSES = [
    OpenBioSingleCellUseExpressionLayer,
    OpenBioSingleCellSnapshotExpression,
    OpenBioSingleCellSubsetObservations,
    OpenBioSingleCellMergeObservationAnnotations,
    OpenBioSingleCellMapGeneIdsFromGTF,
]
