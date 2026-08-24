from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

DYNAMIC_COMBO_TYPE = "COMFY_DYNAMICCOMBO_V3"


def _selected_option(item: Any, named_values: Mapping[str, object]) -> Any:
    if item.id not in named_values:
        raise ValueError(f"workflow is missing required widget value {item.id!r}")
    selected = named_values[item.id]
    option = next((option for option in item.options if option.key == selected), None)
    if option is None:
        raise ValueError(f"{item.id} has no dynamic option named {selected!r}")
    return option


def _is_widget_input(item: Any) -> bool:
    return item.get_io_type() == DYNAMIC_COMBO_TYPE or hasattr(item, "default")


def selected_widget_names(widget_inputs: Iterable[Any], named_values: Mapping[str, object]) -> list[str]:
    names = []
    for item in widget_inputs:
        names.append(item.id)
        if item.get_io_type() == DYNAMIC_COMBO_TYPE:
            option = _selected_option(item, named_values)
            names.extend(f"{item.id}.{nested.id}" for nested in option.inputs)
    return names


def workflow_execute_kwargs(node_class: Any, workflow_node: Mapping[str, Any]) -> dict[str, Any]:
    named_values = workflow_node["widgets_values_named"]
    widget_inputs = [item for item in node_class.GET_SCHEMA().inputs if _is_widget_input(item)]
    expected_names = selected_widget_names(widget_inputs, named_values)
    expected = set(expected_names)
    actual = set(named_values)
    missing = expected - actual
    if missing:
        raise ValueError(f"workflow is missing active widget values {sorted(missing)}")
    unexpected = actual - expected
    if unexpected:
        raise ValueError(f"workflow has unknown or inactive widget values {sorted(unexpected)}")

    kwargs = {}
    for item in widget_inputs:
        if item.get_io_type() != DYNAMIC_COMBO_TYPE:
            kwargs[item.id] = named_values[item.id]
            continue

        option = _selected_option(item, named_values)
        kwargs[item.id] = {
            item.id: named_values[item.id],
            **{nested.id: named_values[f"{item.id}.{nested.id}"] for nested in option.inputs},
        }
    return kwargs
