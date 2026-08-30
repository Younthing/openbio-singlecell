from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal, cast

if TYPE_CHECKING:
    from anndata import AnnData
    from comfy_api.latest import io


ExpressionSourceKind = Literal["X", "raw", "layer"]
DynamicExpressionSource = Mapping[str, object]


@dataclass(frozen=True, slots=True)
class ExpressionSource:
    kind: ExpressionSourceKind
    layer_name: str | None = None
    layer_input_id: str = "layer_name"

    def __post_init__(self) -> None:
        if not isinstance(self.kind, str) or self.kind not in ("X", "raw", "layer"):
            raise ValueError(f"Unsupported expression source: {self.kind!r}")
        if not isinstance(self.layer_input_id, str):
            raise TypeError("Expression source layer input id must be a string.")
        layer_input_id = self.layer_input_id.strip()
        if not layer_input_id:
            raise ValueError("Expression source layer input id cannot be empty.")
        object.__setattr__(self, "layer_input_id", layer_input_id)

        if self.kind == "layer":
            if not isinstance(self.layer_name, str):
                raise TypeError("A layer expression source requires a string layer name.")
            layer_name = self.layer_name.strip()
            if not layer_name:
                raise ValueError("A layer expression source requires a layer name.")
            object.__setattr__(self, "layer_name", layer_name)
        elif self.layer_name is not None:
            raise ValueError(f"Expression source {self.kind!r} cannot carry a layer name.")

    @property
    def use_raw(self) -> bool:
        return self.kind == "raw"

    @property
    def scanpy_layer(self) -> str | None:
        return self.layer_name if self.kind == "layer" else None

    def matrix(self, adata: AnnData) -> Any:
        if self.kind == "raw":
            return adata.raw.X
        if self.kind == "layer":
            return adata.layers[self.layer_name]
        return adata.X

    def parameters(self) -> dict[str, str]:
        parameters = {"source": self.kind}
        if self.kind == "layer":
            parameters[self.layer_input_id] = cast(str, self.layer_name)
        return parameters


@dataclass(frozen=True, slots=True)
class ExpressionSourceSpec:
    """Own the UI shape, default, validation, and provenance for one source input."""

    description: str
    default: ExpressionSourceKind = "X"
    include_raw: bool = False
    layer_input_id: str = "layer_name"
    layer_default: str = "counts"

    def __post_init__(self) -> None:
        if not isinstance(self.description, str):
            raise TypeError("Expression source description must be a string.")
        description = self.description.strip()
        if not description:
            raise ValueError("Expression source description cannot be empty.")
        object.__setattr__(self, "description", description)

        if not isinstance(self.include_raw, bool):
            raise TypeError("Expression source include_raw must be a boolean.")
        if not isinstance(self.layer_input_id, str):
            raise TypeError("Expression source layer input id must be a string.")
        layer_input_id = self.layer_input_id.strip()
        if not layer_input_id:
            raise ValueError("Expression source layer input id cannot be empty.")
        object.__setattr__(self, "layer_input_id", layer_input_id)

        if not isinstance(self.layer_default, str):
            raise TypeError("Expression source layer default must be a string.")
        layer_default = self.layer_default.strip()
        if not layer_default:
            raise ValueError("Expression source layer default cannot be empty.")
        object.__setattr__(self, "layer_default", layer_default)

        if self.default not in self._choices():
            raise ValueError(f"Unsupported default expression source: {self.default!r}")

    def input(self, *, optional: bool = False) -> io.DynamicCombo.Input:
        from comfy_api.latest import io

        choices = self._choices()
        ordered_choices = (self.default, *(choice for choice in choices if choice != self.default))
        return io.DynamicCombo.Input(
            "source",
            optional=optional,
            options=[
                io.DynamicCombo.Option(
                    choice,
                    [io.String.Input(self.layer_input_id, default=self.layer_default)] if choice == "layer" else [],
                )
                for choice in ordered_choices
            ],
        )

    def resolve(
        self,
        adata: AnnData,
        value: DynamicExpressionSource | None,
    ) -> ExpressionSource:
        if value is None:
            value = self._default_value()
        if not isinstance(value, Mapping):
            raise TypeError(f"{self.description} must be a DynamicCombo value.")

        kind_value = value.get("source")
        kind = kind_value.strip() if isinstance(kind_value, str) else kind_value
        if kind not in self._choices():
            raise ValueError(f"Unsupported {self.description.lower()}: {kind!r}")

        if kind == "raw":
            if adata.raw is None:
                raise ValueError(f"{self.description} 'raw' was selected, but adata.raw is unavailable.")
            return ExpressionSource("raw", layer_input_id=self.layer_input_id)

        if kind == "layer":
            layer_value = value.get(self.layer_input_id)
            if layer_value is not None and not isinstance(layer_value, str):
                raise TypeError(f"{self.description} layer must be a string.")
            layer_name = layer_value.strip() if isinstance(layer_value, str) else ""
            if not layer_name:
                raise ValueError(f"{self.description} layer cannot be empty.")
            if layer_name not in adata.layers:
                raise ValueError(f"{self.description} layer not found: {layer_name!r}")
            return ExpressionSource("layer", layer_name, self.layer_input_id)

        return ExpressionSource("X", layer_input_id=self.layer_input_id)

    def _choices(self) -> tuple[ExpressionSourceKind, ...]:
        if self.include_raw:
            return ("X", "raw", "layer")
        return ("X", "layer")

    def _default_value(self) -> dict[str, object]:
        value: dict[str, object] = {"source": self.default}
        if self.default == "layer":
            value[self.layer_input_id] = self.layer_default
        return value
