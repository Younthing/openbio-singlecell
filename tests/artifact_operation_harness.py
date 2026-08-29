from __future__ import annotations

import tempfile
import uuid
import warnings
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from openbio_singlecell.artifact_codecs import read_anndata, write_anndata
from openbio_singlecell.artifact_envelope import summary_from_metadata
from openbio_singlecell.worker_protocol import OperationContext

_ANNDATA_DESCRIPTOR = {"type": "artifact", "kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1"}


def run_anndata_operation(
    operation: Any,
    adata: Any,
    parameters: dict[str, Any],
    *,
    artifact_output: str | None = None,
) -> SimpleNamespace:
    with tempfile.TemporaryDirectory(prefix="openbio-owned-operation-") as directory:
        root = Path(directory)
        input_root = root / "input"
        input_root.mkdir()
        write_anndata(input_root, adata)
        staging = root / "run.partial"
        staging.mkdir()
        context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))
        records = operation(
            context,
            {"adata": {**_ANNDATA_DESCRIPTOR, "path": str(input_root.resolve())}},
            parameters,
        )
        by_name = {record["name"]: record for record in records}
        output = read_anndata(staging / by_name["adata"]["payload"])
        report = summary_from_metadata(by_name["summary"]["value"])
        for message in report.warnings:
            warnings.warn(message, UserWarning, stacklevel=2)
        result = [output]
        if artifact_output is not None:
            record = by_name[artifact_output]
            artifact_root = staging / record["payload"]
            result.append(
                {
                    "kind": record["kind"],
                    "codec": record["codec"],
                    "members": tuple(
                        sorted(
                            path.relative_to(artifact_root).as_posix()
                            for path in artifact_root.rglob("*")
                            if path.is_file()
                        )
                    ),
                }
            )
        result.extend((report, by_name["code"]["value"]))
        return SimpleNamespace(result=tuple(result))


__all__ = ["run_anndata_operation"]
