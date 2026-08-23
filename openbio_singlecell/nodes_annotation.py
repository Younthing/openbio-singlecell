from __future__ import annotations

import json
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/annotation"


def _required_name(value: str, description: str) -> str:
    name = value.strip()
    if not name:
        raise ValueError(f"{description} cannot be empty.")
    return name


def _require_optional_dependency(name: str) -> Any:
    try:
        return __import__(name)
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            f"The {name!r} package is required for this annotation node but is unavailable ({exc})."
        ) from exc


class OpenBioSingleCellCellTypistAnnotation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypistAnnotation",
            display_name="CellTypist Annotation",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("model", default=""),
                io.Boolean.Input("use_raw", default=True),
                io.Boolean.Input("majority_voting", default=True),
                io.String.Input("label_column", default="celltypist_cell_type", advanced=True),
                io.String.Input("confidence_column", default="celltypist_confidence", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        model: str = "",
        use_raw: bool = True,
        majority_voting: bool = True,
        label_column: str = "celltypist_cell_type",
        confidence_column: str = "celltypist_confidence",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        model = _required_name(model, "CellTypist model")
        label_column = _required_name(label_column, "CellTypist label column")
        confidence_column = _required_name(confidence_column, "CellTypist confidence column")
        if label_column == confidence_column:
            raise ValueError("CellTypist label and confidence columns must be different.")
        if use_raw and adata.raw is None:
            raise ValueError("CellTypist use_raw was enabled, but adata.raw is unavailable.")

        celltypist = _require_optional_dependency("celltypist")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)

        prediction_input = adata.raw.to_adata() if use_raw else adata.copy()
        science.sc.pp.normalize_total(prediction_input, target_sum=10_000.0, inplace=True)
        science.sc.pp.log1p(prediction_input)

        predictions = celltypist.annotate(
            prediction_input,
            model=model,
            majority_voting=majority_voting,
        )
        prediction_obs = predictions.to_adata().obs.reindex(adata.obs_names)
        source_label = "majority_voting" if majority_voting else "predicted_labels"

        output = adata.copy()
        output.obs[label_column] = science.pd.Categorical(prediction_obs[source_label].astype(str).to_numpy())
        output.obs[confidence_column] = science.pd.to_numeric(prediction_obs["conf_score"]).to_numpy()
        parameters = {
            "model": model,
            "use_raw": use_raw,
            "majority_voting": majority_voting,
            "label_column": label_column,
            "confidence_column": confidence_column,
            "target_sum": 10_000.0,
            "log1p": True,
        }
        finish_adata(output, "celltypist_annotation", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellMarkerORAAnnotation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerORAAnnotation",
            display_name="Marker ORA Annotation",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("resource_csv", default="openbio-singlecell/markers.csv"),
                io.String.Input("groupby", default="leiden"),
                io.Boolean.Input("use_raw", default=True),
                io.String.Input("source_column", default="cell_type", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.String.Input("output_column", default="ora_annotation", advanced=True),
                io.Int.Input("min_n", default=3, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, resource_csv: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(resource_csv, extensions=(".csv",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource_csv: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(resource_csv, (".csv",))

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resource_csv: str = "openbio-singlecell/markers.csv",
        source_column: str = "cell_type",
        target_column: str = "genesymbol",
        groupby: str = "leiden",
        use_raw: bool = True,
        output_column: str = "ora_annotation",
        min_n: int = 3,
        seed: int = 123,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        source_column = _required_name(source_column, "ORA source column")
        target_column = _required_name(target_column, "ORA target column")
        groupby = _required_name(groupby, "ORA groupby column")
        output_column = _required_name(output_column, "ORA output column")
        if source_column == target_column:
            raise ValueError("ORA source and target columns must be different.")
        if groupby not in adata.obs:
            raise ValueError(f"ORA groupby column not found in obs: {groupby!r}")
        if use_raw and adata.raw is None:
            raise ValueError("Marker ORA with use_raw enabled requires adata.raw.")

        resolved = resolve_input_path(resource_csv, extensions=(".csv",))
        resource = science.pd.read_csv(resolved)
        missing_columns = [column for column in (source_column, target_column) if column not in resource]
        if missing_columns:
            raise ValueError(f"ORA resource CSV is missing columns: {missing_columns}")
        resource = resource.dropna(subset=[source_column, target_column]).copy()
        resource[source_column] = resource[source_column].astype(str)
        resource[target_column] = resource[target_column].astype(str)

        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        if use_raw:
            raw = output.raw.to_adata()
            raw.X = raw.X.astype(science.np.float64)
            if science.sparse.issparse(raw.X):
                raw.X = raw.X.tocsr()
            output.raw = raw
        decoupler.run_ora(
            mat=output,
            net=resource,
            source=source_column,
            target=target_column,
            min_n=min_n,
            seed=seed,
            use_raw=use_raw,
            verbose=False,
        )
        acts = decoupler.get_acts(output, obsm_key="ora_estimate")
        ranked = decoupler.rank_sources_groups(
            acts,
            groupby=groupby,
            reference="rest",
            method="t-test_overestim_var",
        )
        ranked = ranked.dropna(subset=["group", "names"]).copy()
        ranked["_openbio_group"] = ranked["group"].astype(str)
        top = ranked.groupby("_openbio_group", sort=False, observed=True).head(1)
        annotation_map = dict(zip(top["_openbio_group"], top["names"].astype(str), strict=True))

        group_values = output.obs[groupby].astype("string")
        mapped = group_values.map(annotation_map)
        output.obs[output_column] = science.pd.Categorical(mapped.to_numpy())

        parameters = {
            "resource_csv": resource_csv,
            "source_column": source_column,
            "target_column": target_column,
            "groupby": groupby,
            "output_column": output_column,
            "min_n": min_n,
            "seed": seed,
            "use_raw": use_raw,
        }
        finish_adata(
            output,
            "marker_ora_annotation",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellMapClusterAnnotations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMapClusterAnnotations",
            display_name="Map Cluster Annotations",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("mapping_json", default="{}"),
                io.String.Input("output_column", default="cell_type", advanced=True),
                io.Boolean.Input("keep_unmapped", default=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "leiden",
        mapping_json: str = "{}",
        output_column: str = "cell_type",
        keep_unmapped: bool = True,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        groupby = _required_name(groupby, "Cluster groupby column")
        output_column = _required_name(output_column, "Cluster annotation output column")
        if groupby not in adata.obs:
            raise ValueError(f"Cluster groupby column not found in obs: {groupby!r}")

        try:
            raw_mapping = json.loads(mapping_json)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Cluster annotation mapping is not valid JSON ({exc.msg}).") from exc
        if not isinstance(raw_mapping, dict) or not raw_mapping:
            raise ValueError("Cluster annotation mapping must be a non-empty JSON object.")
        mapping = {str(key): str(value) for key, value in raw_mapping.items()}

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        group_values = output.obs[groupby].astype("string")
        mapped = group_values.map(mapping)
        if keep_unmapped:
            mapped = mapped.fillna(group_values)
        output.obs[output_column] = science.pd.Categorical(mapped.to_numpy())

        parameters = {
            "groupby": groupby,
            "mapping": mapping,
            "output_column": output_column,
            "keep_unmapped": keep_unmapped,
        }
        finish_adata(output, "map_cluster_annotations", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


ANNOTATION_NODE_CLASSES = [
    OpenBioSingleCellCellTypistAnnotation,
    OpenBioSingleCellMarkerORAAnnotation,
    OpenBioSingleCellMapClusterAnnotations,
]


__all__ = [
    "ANNOTATION_NODE_CLASSES",
    "OpenBioSingleCellCellTypistAnnotation",
    "OpenBioSingleCellMapClusterAnnotations",
    "OpenBioSingleCellMarkerORAAnnotation",
]
