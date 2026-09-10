"""File exchange for native Monocle 2; Python retains ownership of AnnData."""

from __future__ import annotations

import hashlib
import inspect
import tempfile
from pathlib import Path

from . import expression_source, r_runtime
from .expression_source import ExpressionSourceSpec
from .r_runtime import run_r_script
from .worker_protocol import read_json, write_json

MONOCLE2_KIND = "OPENBIO_MONOCLE2_CDS"
MONOCLE2_CODEC = "monocle2-cds-rds-v1"
MONOCLE2_PACKAGES = ("monocle", "DDRTree", "Matrix", "jsonlite", "Biobase", "VGAM", "igraph", "dplyr", "ggplot2")
MONOCLE2_SCRIPT = Path(__file__).with_name("r") / "monocle2.R"
MONOCLE2_COMPAT_SCRIPT = MONOCLE2_SCRIPT.with_name("monocle2_compat.R")
MONOCLE2_SOURCE = ExpressionSourceSpec(description="Monocle 2 expression source", include_raw=True)


def _hash_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _frame_payload(frame):
    import numpy as np
    import pandas as pd

    def scalar(value):
        if pd.isna(value):
            return None
        if isinstance(value, (float, np.floating)) and np.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return value.item() if isinstance(value, np.generic) else value

    columns = {}
    for name in frame.columns:
        column = frame[name]
        if isinstance(column.dtype, pd.CategoricalDtype):
            record = {"type": "category", "levels": [scalar(value) for value in column.cat.categories],
                      "ordered": bool(column.cat.ordered)}
        elif pd.api.types.is_bool_dtype(column.dtype):
            record = {"type": "boolean"}
        elif pd.api.types.is_integer_dtype(column.dtype):
            record = {"type": "integer"}
        elif pd.api.types.is_numeric_dtype(column.dtype):
            record = {"type": "number"}
        else:
            record = {"type": "string"}
        record["values"] = [scalar(value) for value in column]
        columns[str(name)] = record
    return {"index": frame.index.tolist(), "columns": columns}


def _read_frame(path):
    import pandas as pd

    payload = read_json(path)
    data, aliases = {}, {}
    dtypes = {"string": "string", "integer": "Int64", "number": "Float64", "boolean": "boolean"}
    for name, column in payload["columns"].items():
        if column["type"] == "category":
            data[name] = pd.Categorical(column["values"], categories=column["levels"], ordered=column["ordered"])
        else:
            data[name] = pd.array(column["values"], dtype=dtypes[column["type"]])
        if "original_name" in column:
            aliases[name] = column["original_name"]
    frame = pd.DataFrame(data, index=pd.Index(payload["index"], dtype="object"))
    frame.attrs["r_column_aliases"] = aliases
    return frame


def read_cds(root):
    root = Path(root).resolve(strict=True)
    metadata = read_json(root / "metadata.json")
    if not isinstance(metadata, dict) or metadata.get("kind") != MONOCLE2_KIND or metadata.get("codec") != MONOCLE2_CODEC:
        raise ValueError("Expected a Monocle 2 CellDataSet artifact.")
    if metadata.get("schema_version") != 1 or not isinstance(metadata.get("files"), dict):
        raise ValueError("Monocle 2 artifact metadata is invalid.")
    files = metadata["files"]
    if not {"state.rds", "cell_metadata.json", "ordering_genes.json", "summary.json"}.issubset(files):
        raise ValueError("Monocle 2 artifact is missing native state or result sidecars.")
    for name, digest in files.items():
        path = root / name
        if Path(name).name != name or path.resolve().parent != root or _hash_file(path) != digest:
            raise ValueError(f"Monocle 2 artifact file does not match its recorded content: {name!r}")
    return metadata


def _seal_cds(root, summary, runtime, source):
    names = ["state.rds", "cell_metadata.csv", "cell_metadata.json", "ordering_genes.csv", "ordering_genes.json", "summary.json"]
    if (root / "embedding.csv").exists():
        names.extend(["embedding.csv", "embedding.json"])
    metadata = {
        "schema_version": 1, "kind": MONOCLE2_KIND, "codec": MONOCLE2_CODEC,
        "files": {name: _hash_file(root / name) for name in names},
        "n_obs": summary["n_obs"], "n_vars": summary["n_vars"], "operation": summary["operation"],
        "versions": summary["versions"], "source": source,
        "runtime_version": runtime["r_version"],
    }
    write_json(root / "metadata.json", metadata)
    return metadata


def _attach_results(adata, cells, embedding, summary, parameters):
    import pandas as pd

    if not {"Pseudotime", "State"}.issubset(cells.columns):
        raise ValueError("Monocle 2 export requires orderCells results (State and Pseudotime).")
    cells = cells.set_index("cell_id")
    embedding = embedding.set_index("cell_id")
    for frame in (cells, embedding):
        if not adata.obs_names.is_unique or not frame.index.is_unique:
            raise ValueError("Monocle 2 attachment requires unique cell identifiers.")
        if len(frame) != adata.n_obs or len(frame.index.difference(adata.obs_names)):
            raise ValueError("Monocle 2 results must identify the same cells as the selected AnnData.")
    embedding_key = parameters["embedding_key"]
    pseudotime_key = parameters["pseudotime_key"]
    state_key = parameters["state_key"]
    if not all(isinstance(key, str) and key for key in (embedding_key, pseudotime_key, state_key)):
        raise ValueError("Monocle 2 output names cannot be empty.")
    if pseudotime_key == state_key:
        raise ValueError("Pseudotime and State require distinct observation columns.")
    if not parameters["overwrite_existing"]:
        existing = [key for key in (pseudotime_key, state_key) if key in adata.obs]
        if embedding_key in adata.obsm:
            existing.append(embedding_key)
        if existing:
            raise ValueError(f"Monocle 2 output fields already exist: {existing}; enable overwrite_existing to replace.")
    aligned = cells.reindex(adata.obs_names)
    adata.obs[pseudotime_key] = aligned["Pseudotime"].to_numpy(dtype=float)
    adata.obs[state_key] = pd.Categorical(aligned["State"], categories=list(summary["state_counts"]))
    adata.obsm[embedding_key] = embedding.reindex(adata.obs_names).to_numpy(dtype=float)
    return adata


def run_analysis(operation, runtime, output_dir, parameters, *, adata=None, cds_path=None, script_path=None):
    from scipy import sparse
    from scipy.io import mmwrite

    output_dir = Path(output_dir).resolve()
    parameters = dict(parameters)
    attachment = parameters if operation == "export" else None
    if attachment is not None:
        parameters = {}
    request = {"operation": operation, "output_dir": str(output_dir)}
    source = None
    with tempfile.TemporaryDirectory(prefix=".monocle-exchange-", dir=output_dir.parent) as exchange:
        exchange = Path(exchange)
        if operation == "create":
            selected = MONOCLE2_SOURCE.resolve(adata, parameters.pop("source", None))
            matrix = selected.matrix(adata)
            var = (adata.raw.var if selected.use_raw else adata.var).copy()
            symbol_column = parameters.pop("gene_short_name_column", "")
            var["gene_short_name"] = var[symbol_column] if symbol_column else var.index
            mmwrite(exchange / "expression.mtx", sparse.coo_matrix(matrix.T), precision=17, symmetry="general")
            write_json(exchange / "obs.json", _frame_payload(adata.obs))
            write_json(exchange / "var.json", _frame_payload(var))
            request.update(matrix_path=str(exchange / "expression.mtx"), obs_path=str(exchange / "obs.json"),
                           var_path=str(exchange / "var.json"))
            source = {"expression": selected.parameters(), "gene_short_name_column": symbol_column}
        else:
            metadata = read_cds(cds_path)
            source = metadata.get("source")
            request["input_rds"] = str(Path(cds_path).resolve() / "state.rds")
        request["parameters"] = parameters
        summary = run_r_script(runtime, Path(script_path) if script_path else MONOCLE2_SCRIPT, request, output_dir)
    result = {"summary": summary, "output_dir": output_dir}
    if (output_dir / "state.rds").exists():
        _seal_cds(output_dir, summary, runtime, source)
        result["cds_path"] = output_dir
    table_file = {
        "ordering_genes": "ordering_genes.json", "order_cells": "cell_metadata.json",
        "differential_test": "table.json", "beam": "table.json", "export": "cell_metadata.json",
    }.get(operation)
    if table_file:
        result["table"] = _read_frame(output_dir / table_file)
        if result["table"].attrs["r_column_aliases"]:
            summary["table_column_aliases"] = result["table"].attrs["r_column_aliases"]
    if (output_dir / "plot.png").exists():
        result["plot_path"] = output_dir / "plot.png"
    if attachment is not None:
        result["adata"] = _attach_results(
            adata, result["table"], _read_frame(output_dir / "embedding.json"), summary, attachment,
        )
    return result


_PORTABLE_HEADER = '''from __future__ import annotations
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Literal, cast

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, allow_nan=False), encoding="utf-8")

OUTPUT_TAIL_BYTES = 64 * 1024
ExpressionSourceKind = Literal["X", "raw", "layer"]
DynamicExpressionSource = Mapping[str, object]
'''


def _runtime_source():
    return "\n\n".join(inspect.getsource(function) for function in (
        r_runtime.RRuntimeError, r_runtime._child_environment, r_runtime._description_identity,
        r_runtime._invoke, r_runtime.validate_r_runtime, r_runtime.run_r_script,
    ))


def reproduction_code(operation, parameters, runtime):
    parts = [_PORTABLE_HEADER, _runtime_source(), inspect.getsource(expression_source.ExpressionSource),
             inspect.getsource(expression_source.ExpressionSourceSpec),
             f"MONOCLE2_KIND = {MONOCLE2_KIND!r}\nMONOCLE2_CODEC = {MONOCLE2_CODEC!r}\n"
             "MONOCLE2_SOURCE = ExpressionSourceSpec(description='Monocle 2 expression source', include_raw=True)\n"]
    parts.extend(inspect.getsource(function) for function in (
        _hash_file, _frame_payload, _read_frame, read_cds, _seal_cds, _attach_results, run_analysis,
    ))
    parts.append(f'''
def run_monocle2(adata=None, cds_path=None, *, output_dir, r_runtime=None):
    runtime = {runtime!r} if r_runtime is None else r_runtime
    if isinstance(runtime, str):
        runtime = json.loads(runtime)
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    if {operation!r} == "export":
        adata = adata.copy()
    with tempfile.TemporaryDirectory(prefix="monocle2-source-") as bundle:
        script = Path(bundle) / "monocle2.R"
        script.write_text({MONOCLE2_SCRIPT.read_text(encoding="utf-8")!r}, encoding="utf-8")
        script.with_name("monocle2_compat.R").write_text(
            {MONOCLE2_COMPAT_SCRIPT.read_text(encoding="utf-8")!r}, encoding="utf-8")
        return run_analysis({operation!r}, runtime, destination, {parameters!r},
                            adata=adata, cds_path=cds_path, script_path=script)
''')
    return "\n\n".join(parts) + "\n"


def runtime_reproduction_code(parameters):
    return _PORTABLE_HEADER + "\n" + _runtime_source() + f'''
def probe_monocle2_runtime():
    options = {parameters!r}
    runtime = {{key: options[key] for key in ("r_home", "library_paths", "path_prefix")}}
    executable = shutil.which(options["rscript"], path=_child_environment(runtime).get("PATH"))
    if executable is None:
        raise ValueError("Rscript was not found in the selected environment.")
    executable = os.path.normcase(os.path.realpath(executable))
    stat = Path(executable).stat()
    runtime.update(executable=executable, size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    with tempfile.TemporaryDirectory(prefix="monocle2-probe-") as bundle:
        root = Path(bundle)
        script = root / "probe.R"
        script.write_text({r_runtime.PROBE_SCRIPT.read_text(encoding="utf-8")!r}, encoding="utf-8")
        result = _invoke(runtime, script, {{"packages": {list(MONOCLE2_PACKAGES)!r}, "output_dir": bundle}}, root, timeout=60)
    for package in result["packages"].values():
        package["description"] = _description_identity(Path(package["path"]) / "DESCRIPTION")
    return {{**runtime, **result}}
'''
