from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping, MutableMapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData


ANNDATA_PAYLOAD = "data.h5ad"
ANNDATA_CODEC = "anndata-h5ad-v1"
TABLE_PAYLOAD = "data.jsonl"
TABLE_SCHEMA = "schema.json"
TABLE_CODEC = "table-jsonl-v1"
RESULT_METADATA = "result.json"
PLOT_PAYLOAD = "plot.png"
PLOT_CODEC = "plot-png-v1"
_JSON_UNS_MANIFEST = "__openbio_json_uns_v1__"


def _artifact_root(root: str | os.PathLike[str]) -> Path:
    path = Path(root).resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(f"Artifact root is not a directory: {path}")
    return path


def _payload_path(root: str | os.PathLike[str], name: str) -> Path:
    artifact_root = _artifact_root(root)
    path = artifact_root / name
    if path.parent != artifact_root:
        raise ValueError("Artifact payload path escapes its root.")
    if os.path.lexists(path) and path.resolve().parent != artifact_root:
        raise ValueError("Artifact payload path escapes its root.")
    return path


def _payload_descriptors(root: str | os.PathLike[str], names: tuple[str, ...]) -> list[dict[str, str | int]]:
    return [{"path": name, "size": _payload_path(root, name).stat().st_size} for name in names]


def _nested_json_sequence(value: Any) -> bool:
    return isinstance(value, (list, tuple)) and any(
        item is None or isinstance(item, (Mapping, list, tuple)) for item in value
    )


def _encode_nested_uns(
    mapping: MutableMapping[str, Any],
    *,
    path: tuple[str, ...],
    restored: list[tuple[MutableMapping[str, Any], str, Any]],
    encoded_paths: list[tuple[str, ...]],
    mapping_orders: list[tuple[tuple[str, ...], tuple[str, ...]]],
) -> None:
    keys = tuple(mapping)
    if len(keys) > 1:
        mapping_orders.append((path, keys))
    for key, value in tuple(mapping.items()):
        if not isinstance(key, str):
            raise TypeError("AnnData uns mapping keys must be strings.")
        item_path = (*path, key)
        if isinstance(value, MutableMapping):
            _encode_nested_uns(
                value,
                path=item_path,
                restored=restored,
                encoded_paths=encoded_paths,
                mapping_orders=mapping_orders,
            )
        elif _nested_json_sequence(value):
            payload = json.dumps(
                _strict_json(value, location=f"uns.{'.'.join(item_path)}"),
                allow_nan=False,
                ensure_ascii=False,
                separators=(",", ":"),
            )
            restored.append((mapping, key, value))
            mapping[key] = payload
            encoded_paths.append(item_path)


def _decode_nested_uns(mapping: MutableMapping[str, Any]) -> None:
    manifest = mapping.pop(_JSON_UNS_MANIFEST, None)
    if manifest is None:
        return
    if not isinstance(manifest, Mapping):
        raise TypeError("AnnData nested-uns manifest must be a mapping.")
    encoded_paths: list[list[str]] = []
    mapping_orders: list[tuple[list[str], list[str]]] = []
    for index in range(len(manifest)):
        record_text = manifest.get(str(index))
        if not isinstance(record_text, str):
            raise ValueError("AnnData nested-uns manifest entries must be consecutive strings.")
        record = _strict_json(json.loads(record_text), location="AnnData nested-uns manifest entry")
        if isinstance(record, dict):
            if set(record) != {"mapping", "keys"}:
                raise ValueError("AnnData nested-uns mapping-order record is invalid.")
            order_path = record["mapping"]
            keys = record["keys"]
            if (
                not isinstance(order_path, list)
                or not all(isinstance(key, str) for key in order_path)
                or not isinstance(keys, list)
                or len(keys) < 2
                or not all(isinstance(key, str) for key in keys)
                or len(keys) != len(set(keys))
            ):
                raise ValueError("AnnData nested-uns mapping-order record is invalid.")
            mapping_orders.append((order_path, keys))
            continue
        path = record
        if not isinstance(path, list) or not path or not all(isinstance(key, str) for key in path):
            raise ValueError("AnnData nested-uns paths must be non-empty string lists.")
        encoded_paths.append(path)

    for path in encoded_paths:
        parent: MutableMapping[str, Any] = mapping
        for key in path[:-1]:
            child = parent.get(key)
            if not isinstance(child, MutableMapping):
                raise ValueError("AnnData nested-uns path does not resolve to a mapping.")
            parent = child
        payload = parent.get(path[-1])
        if not isinstance(payload, str):
            raise ValueError("AnnData nested-uns payload must be a JSON string.")
        parent[path[-1]] = _strict_json(
            json.loads(payload, parse_constant=_reject_json_constant),
            location=f"uns.{'.'.join(path)}",
        )

    for path, keys in mapping_orders:
        ordered: MutableMapping[str, Any] = mapping
        for key in path:
            child = ordered.get(key)
            if not isinstance(child, MutableMapping):
                raise ValueError("AnnData nested-uns mapping-order path does not resolve to a mapping.")
            ordered = child
        if set(ordered) != set(keys):
            raise ValueError("AnnData nested-uns mapping-order keys do not match the stored mapping.")
        values = dict(ordered)
        ordered.clear()
        ordered.update((key, values[key]) for key in keys)


def write_anndata(root: str | os.PathLike[str], adata: AnnData) -> list[dict[str, str | int]]:
    from anndata import AnnData, read_h5ad

    if not isinstance(adata, AnnData):
        raise TypeError("Expected an AnnData value.")
    path = _payload_path(root, ANNDATA_PAYLOAD)
    if os.path.lexists(path):
        raise FileExistsError(f"Artifact payload already exists: {ANNDATA_PAYLOAD}")

    if _JSON_UNS_MANIFEST in adata.uns:
        raise ValueError(f"AnnData uns key {_JSON_UNS_MANIFEST!r} is reserved by the artifact codec.")
    restored_uns: list[tuple[MutableMapping[str, Any], str, Any]] = []
    encoded_paths: list[tuple[str, ...]] = []
    mapping_orders: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    try:
        _encode_nested_uns(
            adata.uns,
            path=(),
            restored=restored_uns,
            encoded_paths=encoded_paths,
            mapping_orders=mapping_orders,
        )
        records = [json.dumps(item_path, ensure_ascii=False, separators=(",", ":")) for item_path in encoded_paths]
        records.extend(
            json.dumps(
                {"mapping": item_path, "keys": keys},
                ensure_ascii=False,
                separators=(",", ":"),
            )
            for item_path, keys in mapping_orders
        )
        if records:
            adata.uns[_JSON_UNS_MANIFEST] = {str(index): record for index, record in enumerate(records)}
        adata.write_h5ad(path, compression=None, convert_strings_to_categoricals=False)
    finally:
        adata.uns.pop(_JSON_UNS_MANIFEST, None)
        for mapping, key, value in reversed(restored_uns):
            mapping[key] = value
    stored = read_h5ad(path, backed="r")
    try:
        if stored.shape != adata.shape:
            raise ValueError("Written AnnData shape does not match its owned state.")
        if not stored.obs_names.equals(adata.obs_names) or not stored.var_names.equals(adata.var_names):
            raise ValueError("Written AnnData axes do not match its owned state.")
    finally:
        stored.file.close()
    return _payload_descriptors(root, (ANNDATA_PAYLOAD,))


def read_anndata(root: str | os.PathLike[str]) -> AnnData:
    from anndata import read_h5ad

    adata = read_h5ad(_payload_path(root, ANNDATA_PAYLOAD))
    _decode_nested_uns(adata.uns)
    return adata


def _strict_json(value: Any, *, location: str) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{location} contains a non-finite JSON number.")
        return value
    if isinstance(value, (list, tuple)):
        return [_strict_json(item, location=f"{location}[]") for item in value]
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError(f"{location} JSON object keys must be strings.")
        return {key: _strict_json(item, location=f"{location}.{key}") for key, item in value.items()}
    if type(value).__module__.startswith("numpy") and hasattr(value, "item"):
        return _strict_json(value.item(), location=location)
    raise TypeError(f"{location} contains a non-JSON value of type {type(value).__name__}.")


def _write_json(path: Path, value: Any) -> None:
    if os.path.lexists(path):
        raise FileExistsError(f"Artifact payload already exists: {path.name}")
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(
            _strict_json(value, location=path.name),
            handle,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
        )
        handle.write("\n")


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains a non-finite number: {value}")


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_json_constant)
    return _strict_json(value, location=path.name)


def _scalar(value: Any, *, location: str) -> None | bool | int | float | str:
    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{location} contains a non-finite JSON number.")
        return value
    raise TypeError(f"{location} contains a non-scalar value of type {type(value).__name__}.")


def _column_schema(column: Any, *, json_list: bool) -> dict[str, Any]:
    import pandas as pd
    from pandas.api import types as pandas_types

    dtype = column.dtype
    if json_list:
        if not pandas_types.is_object_dtype(dtype):
            raise TypeError("JSON-list columns must use object dtype.")
        return {"kind": "json-list", "dtype": "object"}
    if isinstance(dtype, pd.CategoricalDtype):
        categories = [
            _scalar(value, location=f"category {value!r}")
            for value in dtype.categories.tolist()
        ]
        if any(value is None for value in categories):
            raise ValueError("Categorical levels cannot be null.")
        return {
            "kind": "categorical",
            "categories": categories,
            "ordered": bool(dtype.ordered),
        }
    if pandas_types.is_extension_array_dtype(dtype):
        if not any(
            predicate(dtype)
            for predicate in (
                pandas_types.is_bool_dtype,
                pandas_types.is_integer_dtype,
                pandas_types.is_float_dtype,
                pandas_types.is_string_dtype,
            )
        ):
            raise TypeError(f"Unsupported table dtype: {dtype}")
        return {"kind": "nullable", "dtype": str(dtype)}
    if any(
        predicate(dtype)
        for predicate in (
            pandas_types.is_bool_dtype,
            pandas_types.is_integer_dtype,
            pandas_types.is_float_dtype,
            pandas_types.is_object_dtype,
        )
    ):
        return {"kind": "scalar", "dtype": str(dtype)}
    raise TypeError(f"Unsupported table dtype: {dtype}")


def _index_schema(index: Any) -> dict[str, Any]:
    import pandas as pd
    from pandas.api import types as pandas_types

    if isinstance(index, pd.MultiIndex):
        raise TypeError("Table MultiIndex values are not portable.")
    name = _scalar(index.name, location="table index name")
    if isinstance(index, pd.RangeIndex):
        return {
            "kind": "range",
            "name": name,
            "start": int(index.start),
            "stop": int(index.stop),
            "step": int(index.step),
        }
    if not any(
        predicate(index.dtype)
        for predicate in (
            pandas_types.is_bool_dtype,
            pandas_types.is_integer_dtype,
            pandas_types.is_float_dtype,
            pandas_types.is_string_dtype,
            pandas_types.is_object_dtype,
        )
    ):
        raise TypeError(f"Unsupported table index dtype: {index.dtype}")
    return {"kind": "scalar", "name": name, "dtype": str(index.dtype)}


def _encode_cell(value: Any, schema: Mapping[str, Any], *, location: str) -> Any:
    from pandas.api.types import is_float_dtype

    if schema["kind"] == "json-list":
        if value is None or type(value).__name__ in {"NAType", "NaTType"}:
            return None
        if not isinstance(value, list):
            raise TypeError(f"{location} must be a JSON list or null.")
        return _strict_json(value, location=location)
    if schema["kind"] == "categorical":
        try:
            if bool(value != value):
                return None
        except (TypeError, ValueError):
            pass
    if schema["kind"] in {"scalar", "nullable"} and is_float_dtype(schema["dtype"]):
        scalar = value.item() if hasattr(value, "item") else value
        if isinstance(scalar, float) and math.isinf(scalar):
            # The existing dtype-based reader restores these quoted values as floats.
            # Keep JSON valid without turning native R infinities into missing data.
            return "Infinity" if scalar > 0 else "-Infinity"
    return _scalar(value, location=location)


def write_table(
    root: str | os.PathLike[str],
    table: Any,
    result_metadata: Mapping[str, Any],
    *,
    json_list_columns: Sequence[str] = (),
) -> list[dict[str, str | int]]:
    import pandas as pd

    if not isinstance(table, pd.DataFrame):
        raise TypeError("Expected a pandas DataFrame.")
    if isinstance(table.columns, pd.MultiIndex):
        raise TypeError("Table MultiIndex columns are not portable.")
    columns = list(table.columns)
    if not all(isinstance(column, str) for column in columns):
        raise TypeError("Table column names must be strings.")
    if len(set(columns)) != len(columns):
        raise ValueError("Table column names must be unique.")
    declared_lists = set(json_list_columns)
    unknown_lists = declared_lists.difference(columns)
    if unknown_lists:
        raise ValueError(f"Unknown JSON-list columns: {sorted(unknown_lists)!r}")

    column_schemas = [
        {"name": name, **_column_schema(table[name], json_list=name in declared_lists)}
        for name in columns
    ]
    schema = {
        "version": 1,
        "format": "openbio-table-jsonl",
        "index": _index_schema(table.index),
        "columns": column_schemas,
    }
    data_path = _payload_path(root, TABLE_PAYLOAD)
    if os.path.lexists(data_path):
        raise FileExistsError(f"Artifact payload already exists: {TABLE_PAYLOAD}")
    with data_path.open("x", encoding="utf-8", newline="\n") as handle:
        for row_number, row in enumerate(table.itertuples(index=True, name=None)):
            encoded = [_scalar(row[0], location=f"table index row {row_number}")]
            encoded.extend(
                _encode_cell(value, column_schema, location=f"table column {column_schema['name']!r} row {row_number}")
                for value, column_schema in zip(row[1:], column_schemas, strict=True)
            )
            json.dump(encoded, handle, allow_nan=False, ensure_ascii=False, separators=(",", ":"))
            handle.write("\n")

    _write_json(_payload_path(root, TABLE_SCHEMA), schema)
    _write_json(_payload_path(root, RESULT_METADATA), dict(result_metadata))
    return _payload_descriptors(root, (TABLE_PAYLOAD, TABLE_SCHEMA, RESULT_METADATA))


def _restore_column(values: list[Any], schema: Mapping[str, Any]) -> Any:
    import pandas as pd

    kind = schema["kind"]
    if kind == "categorical":
        return pd.Categorical(values, categories=schema["categories"], ordered=schema["ordered"])
    if kind == "nullable":
        return pd.array(values, dtype=schema["dtype"])
    if kind == "json-list" or schema["dtype"] == "object":
        return pd.Series(values, dtype=object)
    return pd.Series(values, dtype=schema["dtype"])


def read_table(root: str | os.PathLike[str]) -> tuple[Any, dict[str, Any]]:
    import pandas as pd

    schema = _read_json(_payload_path(root, TABLE_SCHEMA))
    result_metadata = _read_json(_payload_path(root, RESULT_METADATA))
    column_schemas = schema["columns"]
    rows: list[list[Any]] = []
    with _payload_path(root, TABLE_PAYLOAD).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            row = _strict_json(
                json.loads(line, parse_constant=_reject_json_constant),
                location=f"{TABLE_PAYLOAD} line {line_number}",
            )
            if not isinstance(row, list) or len(row) != len(column_schemas) + 1:
                raise ValueError(f"Invalid table JSONL row at line {line_number}.")
            rows.append(row)

    index_values = [row[0] for row in rows]
    index_schema = schema["index"]
    if index_schema["kind"] == "range":
        index = pd.RangeIndex(
            start=index_schema["start"],
            stop=index_schema["stop"],
            step=index_schema["step"],
            name=index_schema["name"],
        )
        if list(index) != index_values:
            raise ValueError("Table JSONL index does not match its schema.")
    else:
        index = pd.Index(index_values, dtype=index_schema["dtype"], name=index_schema["name"])

    data = {
        column_schema["name"]: _restore_column(
            [row[column_number + 1] for row in rows],
            column_schema,
        )
        for column_number, column_schema in enumerate(column_schemas)
    }
    table = pd.DataFrame(data, columns=[column["name"] for column in column_schemas])
    table.index = index
    if not isinstance(result_metadata, dict):
        raise TypeError("Table result metadata must be a JSON object.")
    return table, result_metadata


def _validate_png(path: Path) -> None:
    from PIL import Image, UnidentifiedImageError

    try:
        with Image.open(path) as image:
            if image.format != "PNG":
                raise ValueError("Plot payload must be a PNG image.")
            image.verify()
    except (OSError, UnidentifiedImageError) as error:
        raise ValueError("Plot payload is not a valid PNG image.") from error


def write_plot(
    root: str | os.PathLike[str],
    png: bytes,
    result_metadata: Mapping[str, Any],
) -> list[dict[str, str | int]]:
    if not isinstance(png, bytes):
        raise TypeError("Plot PNG payload must be bytes.")
    path = _payload_path(root, PLOT_PAYLOAD)
    if os.path.lexists(path):
        raise FileExistsError(f"Artifact payload already exists: {PLOT_PAYLOAD}")
    with path.open("xb") as handle:
        handle.write(png)
    _validate_png(path)
    _write_json(_payload_path(root, RESULT_METADATA), dict(result_metadata))
    return _payload_descriptors(root, (PLOT_PAYLOAD, RESULT_METADATA))


def read_plot(root: str | os.PathLike[str]) -> tuple[bytes, dict[str, Any]]:
    path = _payload_path(root, PLOT_PAYLOAD)
    _validate_png(path)
    metadata = _read_json(_payload_path(root, RESULT_METADATA))
    if not isinstance(metadata, dict):
        raise TypeError("Plot result metadata must be a JSON object.")
    return path.read_bytes(), metadata
