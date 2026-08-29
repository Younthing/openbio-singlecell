from __future__ import annotations

from pathlib import Path
from typing import Any

from .artifact_codecs import read_anndata, read_table, write_anndata, write_table
from .artifact_envelope import write_strict_json
from .scenic_artifact import (
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    validate_scenic_result_artifact,
)
from .tf_activity_artifact import (
    TF_ACTIVITY_ARTIFACT_TYPE,
    TFActivityArtifact,
    validate_tf_activity_artifact,
)
from .worker_protocol import read_json

TF_ACTIVITY_CODEC = "tf-activity-h5ad-v1"
SCENIC_CODEC = "scenic-activity-membership-v1"
CODEC_VERSION = 1


def _root(directory: str | Path) -> Path:
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Artifact root is not a directory: {root}")
    return root


def _frame(matrix: Any, *, index: Any, columns: Any) -> Any:
    import pandas as pd
    from scipy import sparse

    values = matrix.toarray() if sparse.issparse(matrix) else matrix
    return pd.DataFrame(values, index=index, columns=columns)


def write_tf_activity(directory: str | Path, artifact: TFActivityArtifact) -> None:
    import anndata as ad
    import pandas as pd

    root = _root(directory)
    scores, adjusted, provenance, metadata = validate_tf_activity_artifact(
        artifact,
        exact_type=True,
        copy_result=False,
    )
    activity = ad.AnnData(
        X=scores.to_numpy(dtype=float, copy=False),
        obs=pd.DataFrame(index=scores.index),
        var=pd.DataFrame(index=scores.columns),
    )
    activity.layers["adjusted_pvalues"] = adjusted.to_numpy(dtype=float, copy=False)
    write_anndata(root, activity)
    write_strict_json(
        root / "result.json",
        {
            "codec": TF_ACTIVITY_CODEC,
            "version": CODEC_VERSION,
            "artifact_type": TF_ACTIVITY_ARTIFACT_TYPE,
            "provenance": provenance,
            "metadata": metadata,
        },
    )


def read_tf_activity(directory: str | Path) -> dict[str, Any]:
    root = _root(directory)
    activity = read_anndata(root)
    if set(activity.layers) - {None} != {"adjusted_pvalues"}:
        raise ValueError("TF activity H5AD must contain exactly the adjusted_pvalues layer.")
    envelope = read_json(root / "result.json")
    if (
        not isinstance(envelope, dict)
        or set(envelope) != {"codec", "version", "artifact_type", "provenance", "metadata"}
        or envelope["codec"] != TF_ACTIVITY_CODEC
        or envelope["version"] != CODEC_VERSION
        or envelope["artifact_type"] != TF_ACTIVITY_ARTIFACT_TYPE
    ):
        raise ValueError("TF activity artifact envelope is invalid.")
    result = {
        "artifact_type": TF_ACTIVITY_ARTIFACT_TYPE,
        "scores": _frame(activity.X, index=activity.obs_names, columns=activity.var_names),
        "adjusted_pvalues": _frame(
            activity.layers["adjusted_pvalues"],
            index=activity.obs_names,
            columns=activity.var_names,
        ),
        "provenance": envelope["provenance"],
        "metadata": envelope["metadata"],
    }
    validate_tf_activity_artifact(result, exact_type=False, copy_result=False)
    return result


def write_scenic(directory: str | Path, artifact: Any) -> None:
    import anndata as ad
    import pandas as pd

    root = _root(directory)
    activity, membership, provenance, metadata = validate_scenic_result_artifact(
        artifact,
        exact_type=False,
        numpy=__import__("numpy"),
        pandas=pd,
        copy_result=False,
    )
    activity_root = root / "activity"
    membership_root = root / "membership"
    activity_root.mkdir(exist_ok=False)
    membership_root.mkdir(exist_ok=False)
    activity_adata = ad.AnnData(
        X=activity.to_numpy(dtype=float, copy=False),
        obs=pd.DataFrame(index=activity.index),
        var=pd.DataFrame(index=activity.columns),
    )
    write_anndata(activity_root, activity_adata)
    write_table(
        membership_root,
        membership,
        {"codec": SCENIC_CODEC, "version": CODEC_VERSION},
        json_list_columns=("motif_ids",),
    )
    write_strict_json(
        root / "result.json",
        {
            "codec": SCENIC_CODEC,
            "version": CODEC_VERSION,
            "artifact_type": SCENIC_RESULT_ARTIFACT_TYPE,
            "artifact_schema_version": SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
            "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
            "provenance": provenance,
        },
    )


def read_scenic(directory: str | Path) -> dict[str, Any]:
    import pandas as pd

    root = _root(directory)
    activity_adata = read_anndata(root / "activity")
    membership, membership_metadata = read_table(root / "membership")
    if membership_metadata != {"codec": SCENIC_CODEC, "version": CODEC_VERSION}:
        raise ValueError("SCENIC membership artifact envelope is invalid.")
    envelope = read_json(root / "result.json")
    if (
        not isinstance(envelope, dict)
        or set(envelope)
        != {
            "codec",
            "version",
            "artifact_type",
            "artifact_schema_version",
            "artifact_fingerprint_sha256",
            "provenance",
        }
        or envelope["codec"] != SCENIC_CODEC
        or envelope["version"] != CODEC_VERSION
        or envelope["artifact_type"] != SCENIC_RESULT_ARTIFACT_TYPE
        or envelope["artifact_schema_version"] != SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION
    ):
        raise ValueError("SCENIC artifact envelope is invalid.")
    observations = activity_adata.obs_names.tolist()
    regulons = activity_adata.var_names.tolist()
    result = {
        "artifact_type": SCENIC_RESULT_ARTIFACT_TYPE,
        "artifact_schema_version": SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
        "observations": observations,
        "regulons": regulons,
        "activity": _frame(activity_adata.X, index=observations, columns=regulons),
        "membership": membership,
        "provenance": envelope["provenance"],
        "artifact_fingerprint_sha256": envelope["artifact_fingerprint_sha256"],
    }
    validate_scenic_result_artifact(
        result,
        exact_type=False,
        numpy=__import__("numpy"),
        pandas=pd,
        copy_result=False,
    )
    return result


__all__ = [
    "SCENIC_CODEC",
    "TF_ACTIVITY_CODEC",
    "read_scenic",
    "read_tf_activity",
    "write_scenic",
    "write_tf_activity",
]
