from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from .expression_source import ExpressionSource

if TYPE_CHECKING:
    from anndata import AnnData


def _history_entries(adata: AnnData) -> list[Mapping[str, Any]]:
    metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(metadata, Mapping):
        return []
    history = metadata.get("analysis_history")
    if not isinstance(history, Mapping):
        return []
    return [entry for entry in history.values() if isinstance(entry, Mapping)]


def resolve_expression_state(adata: AnnData, expression: ExpressionSource) -> tuple[str, str | None]:
    """Resolve the latest known semantic state of one current-axis expression representation."""
    x_state = "unknown"
    x_evidence = None
    layer_states: dict[str, tuple[str, str]] = {}
    for entry in _history_entries(adata):
        operation = entry.get("operation")
        parameters = entry.get("parameters")
        parameters = parameters if isinstance(parameters, Mapping) else {}
        if operation == "snapshot_expression":
            layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
            if parameters.get("source") == "X":
                x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
        elif operation == "normalize_to_layer":
            output_layer = parameters.get("output_layer")
            transform = parameters.get("transform")
            if isinstance(output_layer, str):
                state = "logged" if transform == "log1p" else "transformed" if transform == "sqrt" else "normalized"
                layer_states[output_layer] = (state, "OpenBio Normalize To Layer history")
        elif operation == "pearson_residuals_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("pearson_residuals", "OpenBio Pearson Residuals history")
        elif operation == "scale_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("scaled", "OpenBio Scale history")

        if operation == "normalize_total":
            x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
        elif operation == "log1p":
            x_state, x_evidence = "logged", "OpenBio Log1p history"

    if expression.kind == "layer":
        return layer_states.get(str(expression.layer_name), ("unknown", None))
    if x_state == "unknown" and isinstance(adata.uns.get("log1p"), Mapping):
        return "logged", "AnnData uns['log1p'] marker"
    return x_state, x_evidence


__all__ = ["resolve_expression_state"]
