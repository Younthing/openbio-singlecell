from __future__ import annotations

from pathlib import Path
from typing import Any

from .artifact_codecs import read_table, write_table
from .liana_result import LIANA_RESULT_ARTIFACT_TYPE, LianaResult, validate_liana_result

LIANA_CODEC = "liana-table-jsonl-v1"
CODEC_VERSION = 1


def write_liana_result(directory: str | Path, result: LianaResult) -> None:
    table, provenance, metadata = validate_liana_result(result, copy_result=False)
    write_table(
        directory,
        table,
        {
            "codec": LIANA_CODEC,
            "version": CODEC_VERSION,
            "artifact_type": LIANA_RESULT_ARTIFACT_TYPE,
            "provenance": provenance,
            "metadata": metadata,
        },
    )


def read_liana_result(directory: str | Path) -> dict[str, Any]:
    table, envelope = read_table(directory)
    if (
        not isinstance(envelope, dict)
        or set(envelope) != {
            "codec",
            "version",
            "artifact_type",
            "provenance",
            "metadata",
        }
        or envelope["codec"] != LIANA_CODEC
        or envelope["version"] != CODEC_VERSION
        or envelope["artifact_type"] != LIANA_RESULT_ARTIFACT_TYPE
    ):
        raise ValueError("LIANA artifact envelope is invalid.")
    result = {
        "artifact_type": LIANA_RESULT_ARTIFACT_TYPE,
        "table": table,
        "provenance": envelope["provenance"],
        "metadata": envelope["metadata"],
    }
    validate_liana_result(result, exact_type=False, copy_result=False)
    return result


__all__ = ["LIANA_CODEC", "read_liana_result", "write_liana_result"]
