from __future__ import annotations

import hashlib
import json
import math
from typing import Any

PLOT_EVIDENCE_SCHEMA = "openbio-singlecell/plot-evidence/v1"


def table_evidence_fingerprint(frame: Any, *, contract: str) -> str:
    """Hash one canonical result table without depending on pandas dtype round-trips."""
    import pandas as pd

    if not isinstance(contract, str) or not contract or contract != contract.strip():
        raise ValueError("Plot evidence contract must be canonical nonblank text.")
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Plot evidence must be a pandas DataFrame.")

    def plain(value: Any) -> Any:
        if value is None:
            return ["none", None]
        if isinstance(value, bool):
            return ["bool", value]
        if isinstance(value, str):
            return ["str", value]
        if hasattr(value, "item"):
            return plain(value.item())
        if isinstance(value, int):
            return ["int", value]
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError("Plot evidence fingerprint rejects non-finite values.")
            return ["float", value.hex()]
        if isinstance(value, (list, tuple)):
            return ["list", [plain(item) for item in value]]
        raise TypeError(f"Plot evidence fingerprint does not support {type(value).__name__} values.")

    payload = {
        "schema": PLOT_EVIDENCE_SCHEMA,
        "contract": contract,
        "columns": frame.columns.tolist(),
        "rows": [[plain(value) for value in row] for row in frame.itertuples(index=False, name=None)],
    }
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


__all__ = ["PLOT_EVIDENCE_SCHEMA", "table_evidence_fingerprint"]
