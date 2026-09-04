from __future__ import annotations

from typing import Any


def differential_table_content_fingerprint(table: Any, *, contract: str) -> str:
    """Hash one canonical differential-evidence table, preserving scientific nulls."""
    import hashlib
    import json
    import math

    import pandas as pd

    if not isinstance(contract, str) or not contract or contract != contract.strip():
        raise ValueError("Differential evidence contract must be canonical nonblank text.")
    if not isinstance(table, pd.DataFrame):
        raise TypeError("Differential evidence must be a pandas DataFrame.")

    def scalar(value):
        if value is None or bool(pd.isna(value)):
            return ["null", None]
        if isinstance(value, bool):
            return ["bool", value]
        if isinstance(value, str):
            return ["str", value]
        if hasattr(value, "item"):
            return scalar(value.item())
        if isinstance(value, int):
            return ["int", value]
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError("Differential evidence fingerprint rejects infinite values.")
            return ["float", value.hex()]
        raise TypeError(f"Differential evidence fingerprint does not support {type(value).__name__} values.")

    payload = {
        "schema": "openbio-singlecell/differential-table-evidence/v1",
        "contract": contract,
        "columns": table.columns.tolist(),
        "rows": [[scalar(value) for value in row] for row in table.itertuples(index=False, name=None)],
    }
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


__all__ = ["differential_table_content_fingerprint"]
