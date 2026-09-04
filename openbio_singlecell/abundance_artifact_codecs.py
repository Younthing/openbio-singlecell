from __future__ import annotations

from pathlib import Path

from .artifact_codecs import read_anndata, read_table, write_anndata, write_table
from .composition_result import (
    COMPOSITION_MODEL_ARTIFACT_TYPE,
    CompositionModelResult,
    build_composition_model_result,
    validate_composition_model_result,
)
from .milo_result import (
    MILO_RESULT_ARTIFACT_TYPE,
    MiloResult,
    build_milo_result,
    validate_milo_result,
)
from .worker_protocol import read_json, write_json

MILO_RESULT_CODEC = "milo-result-h5ad-v1"
MILO_RESULT_KIND = MILO_RESULT_ARTIFACT_TYPE
MILO_RESULT_METADATA = "metadata.json"
COMPOSITION_MODEL_CODEC = "composition-inferencedata-zarr-v1"
COMPOSITION_MODEL_KIND = COMPOSITION_MODEL_ARTIFACT_TYPE
COMPOSITION_MODEL_METADATA = "metadata.json"
CODEC_VERSION = 1


def _root(directory: str | Path) -> Path:
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Abundance artifact root is not a directory: {root}")
    return root


def _descriptors(root: Path) -> list[dict[str, str | int]]:
    return [
        {"path": path.relative_to(root).as_posix(), "size": path.stat().st_size}
        for path in sorted(root.rglob("*"))
        if path.is_file()
    ]


def write_milo_result(directory: str | Path, result: MiloResult) -> list[dict[str, str | int]]:
    import anndata as ad
    import pandas as pd

    root = _root(directory)
    table, membership, graph, coordinates, observations, neighborhoods, provenance, metadata = (
        validate_milo_result(result, copy_result=False)
    )
    payload = ad.AnnData(
        X=membership,
        obs=pd.DataFrame(index=observations),
        var=pd.DataFrame(index=neighborhoods),
    )
    payload.varp["nhood_connectivities"] = graph
    payload.varm["X_milo_graph"] = coordinates
    write_anndata(root, payload)
    table_root = root / "table"
    table_root.mkdir(exist_ok=False)
    write_table(table_root, table, {"codec": MILO_RESULT_CODEC, "version": CODEC_VERSION})
    write_json(
        root / MILO_RESULT_METADATA,
        {
            "codec": MILO_RESULT_CODEC,
            "version": CODEC_VERSION,
            "artifact_type": MILO_RESULT_ARTIFACT_TYPE,
            "provenance": provenance,
            "metadata": metadata,
        },
    )
    return _descriptors(root)


def read_milo_result(directory: str | Path) -> MiloResult:
    root = _root(directory)
    envelope = read_json(root / MILO_RESULT_METADATA)
    if (
        not isinstance(envelope, dict)
        or set(envelope) != {"codec", "version", "artifact_type", "provenance", "metadata"}
        or envelope["codec"] != MILO_RESULT_CODEC
        or envelope["version"] != CODEC_VERSION
        or envelope["artifact_type"] != MILO_RESULT_ARTIFACT_TYPE
    ):
        raise ValueError("Milo result artifact envelope is invalid.")
    table, table_metadata = read_table(root / "table")
    if table_metadata != {"codec": MILO_RESULT_CODEC, "version": CODEC_VERSION}:
        raise ValueError("Milo result table envelope is invalid.")
    payload = read_anndata(root)
    if (
        len(payload.obs.columns)
        or len(payload.var.columns)
        or payload.layers
        or payload.obsm
        or payload.obsp
        or payload.uns
        or set(payload.varm) != {"X_milo_graph"}
        or set(payload.varp) != {"nhood_connectivities"}
    ):
        raise ValueError("Milo result H5AD has an invalid exact slot schema.")
    provenance = envelope["provenance"]
    if not isinstance(provenance, dict):
        raise ValueError("Milo result provenance must be a JSON object.")
    result = build_milo_result(
        table=table,
        membership=payload.X,
        graph=payload.varp["nhood_connectivities"],
        representative_coordinates=payload.varm["X_milo_graph"],
        observation_names=payload.obs_names.tolist(),
        neighborhood_names=payload.var_names.tolist(),
        representation_key=provenance.get("representation_key"),
        representation_sha256=provenance.get("representation_sha256"),
        condition_key=provenance.get("condition_key"),
        annotation_key=provenance.get("annotation_key"),
        annotation_status=provenance.get("annotation_status"),
        n_neighbors=provenance.get("n_neighbors"),
        neighborhood_proportion=provenance.get("neighborhood_proportion"),
        mixed_annotation_threshold=provenance.get("mixed_annotation_threshold"),
        spatial_fdr_threshold=provenance.get("spatial_fdr_threshold"),
        min_abs_log2_fold_change=provenance.get("min_abs_log2_fold_change"),
        random_seed=provenance.get("random_seed"),
    )
    if result.provenance != provenance or result.metadata != envelope["metadata"]:
        raise ValueError("Milo result artifact failed its current-content fingerprint check.")
    return result


def write_composition_model_result(
    directory: str | Path,
    result: CompositionModelResult,
) -> list[dict[str, str | int]]:
    root = _root(directory)
    table, _posterior, _sample_stats, cell_types, hierarchy, model_metadata, metadata = (
        validate_composition_model_result(result, copy_result=False)
    )
    table_root = root / "table"
    table_root.mkdir(exist_ok=False)
    write_table(table_root, table, {"codec": COMPOSITION_MODEL_CODEC, "version": CODEC_VERSION})
    inference_root = root / "inference_data"
    result.inference_data().to_zarr(inference_root, mode="w-", consolidated=False)
    write_json(
        root / COMPOSITION_MODEL_METADATA,
        {
            "codec": COMPOSITION_MODEL_CODEC,
            "version": CODEC_VERSION,
            "artifact_type": COMPOSITION_MODEL_ARTIFACT_TYPE,
            "method": result.method,
            "cell_types": cell_types,
            "hierarchy": hierarchy,
            "model_metadata": model_metadata,
            "metadata": metadata,
        },
    )
    return _descriptors(root)


def read_composition_model_result(directory: str | Path) -> CompositionModelResult:
    import arviz as az

    root = _root(directory)
    envelope = read_json(root / COMPOSITION_MODEL_METADATA)
    expected_fields = {
        "codec",
        "version",
        "artifact_type",
        "method",
        "cell_types",
        "hierarchy",
        "model_metadata",
        "metadata",
    }
    if (
        not isinstance(envelope, dict)
        or set(envelope) != expected_fields
        or envelope["codec"] != COMPOSITION_MODEL_CODEC
        or envelope["version"] != CODEC_VERSION
        or envelope["artifact_type"] != COMPOSITION_MODEL_ARTIFACT_TYPE
    ):
        raise ValueError("Composition-model artifact envelope is invalid.")
    table, table_metadata = read_table(root / "table")
    if table_metadata != {"codec": COMPOSITION_MODEL_CODEC, "version": CODEC_VERSION}:
        raise ValueError("Composition-model table envelope is invalid.")
    inference_data = az.from_zarr(root / "inference_data", backend_kwargs={"consolidated": False})
    if set(inference_data.children) != {"posterior", "sample_stats"}:
        raise ValueError("Composition-model InferenceData groups are invalid.")
    posterior_dataset = inference_data["posterior"].to_dataset()
    sample_stats_dataset = inference_data["sample_stats"].to_dataset()
    metadata = envelope["metadata"]
    if not isinstance(metadata, dict):
        raise ValueError("Composition-model artifact metadata must be a JSON object.")
    draw_count = metadata.get("draw_count")
    if posterior_dataset.sizes.get("chain") != 1 or posterior_dataset.sizes.get("draw") != draw_count:
        raise ValueError("Composition-model posterior chain/draw axes are invalid.")
    if sample_stats_dataset.sizes.get("chain") != 1 or sample_stats_dataset.sizes.get("draw") != draw_count:
        raise ValueError("Composition-model sample-stat chain/draw axes are invalid.")
    method = envelope["method"]
    posterior_axes = {
        "sccoda": {"intercept": "cell_type", "condition_effect": "cell_type"},
        "tasccoda": {
            "intercept": "cell_type",
            "hierarchy_node_effect": "hierarchy_node",
            "derived_leaf_effect": "cell_type",
            "theta": None,
        },
    }.get(method)
    if posterior_axes is None or set(posterior_dataset.data_vars) != set(posterior_axes):
        raise ValueError("Composition-model posterior variables are invalid.")
    for name, axis in posterior_axes.items():
        expected_dims = ("chain", "draw", axis) if axis is not None else ("chain", "draw")
        if posterior_dataset[name].dims != expected_dims:
            raise ValueError(f"Composition-model posterior variable {name!r} dimensions are invalid.")
    if set(sample_stats_dataset.data_vars) != {"potential_energy", "num_steps", "step_size"} or any(
        sample_stats_dataset[name].dims != ("chain", "draw") for name in sample_stats_dataset.data_vars
    ):
        raise ValueError("Composition-model sample-stat variables are invalid.")
    hierarchy = envelope["hierarchy"]
    expected_posterior_coords = {"chain", "draw", "cell_type"}
    if method == "tasccoda":
        expected_posterior_coords.add("hierarchy_node")
    if set(posterior_dataset.coords) != expected_posterior_coords or set(sample_stats_dataset.coords) != {
        "chain",
        "draw",
    }:
        raise ValueError("Composition-model InferenceData coordinate schema is invalid.")
    if (
        posterior_dataset.coords["chain"].values.tolist() != [0]
        or posterior_dataset.coords["draw"].values.tolist() != list(range(draw_count))
        or posterior_dataset.coords["cell_type"].values.tolist() != envelope["cell_types"]
        or sample_stats_dataset.coords["chain"].values.tolist() != [0]
        or sample_stats_dataset.coords["draw"].values.tolist() != list(range(draw_count))
        or (
            method == "tasccoda"
            and (
                not isinstance(hierarchy, dict)
                or posterior_dataset.coords["hierarchy_node"].values.tolist() != hierarchy.get("node_names")
            )
        )
    ):
        raise ValueError("Composition-model InferenceData axes do not match the artifact envelope.")
    posterior = {
        name: posterior_dataset[name].values[0].copy()
        for name in posterior_dataset.data_vars
    }
    sample_stats = {
        name: sample_stats_dataset[name].values[0].copy()
        for name in sample_stats_dataset.data_vars
    }
    inference_data.close()
    result = build_composition_model_result(
        table=table,
        method=method,
        posterior=posterior,
        sample_stats=sample_stats,
        cell_types=envelope["cell_types"],
        hierarchy=envelope["hierarchy"],
        model_metadata=envelope["model_metadata"],
    )
    if result.metadata != metadata:
        raise ValueError("Composition-model artifact failed its current-content fingerprint check.")
    return result


__all__ = [
    "COMPOSITION_MODEL_CODEC",
    "COMPOSITION_MODEL_KIND",
    "MILO_RESULT_CODEC",
    "MILO_RESULT_KIND",
    "read_composition_model_result",
    "read_milo_result",
    "write_composition_model_result",
    "write_milo_result",
]
