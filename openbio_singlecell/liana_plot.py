from __future__ import annotations

import copy
import importlib
import inspect
import io as bytes_io
import json
import math
import platform
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from typing import Any

from .liana_communication import LIANA_AUDITED_VERSION, LIANA_METHOD_PROFILES
from .liana_result import (
    _liana_result_table_identity,
    liana_result_portable_code,
    validate_liana_result,
)

LIANA_DOT_PLOT_NODE_ID = "OpenBioSingleCellLianaDotPlot"
LIANA_DOT_PLOT_SIGNATURE = (
    "adata",
    "uns_key",
    "liana_res",
    "sample_key",
    "colour",
    "size",
    "inverse_colour",
    "inverse_size",
    "source_labels",
    "target_labels",
    "ligand_complex",
    "receptor_complex",
    "size_range",
    "cmap",
    "figure_size",
    "return_fig",
)
LIANA_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
LIANA_PLOT_REPORT_LIMIT = 20


def _liana_plot_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    result = value.strip()
    if not result or result != value:
        raise ValueError(f"{label} must be a canonical nonblank string.")
    if any(ord(character) < 32 or ord(character) == 127 for character in result):
        raise ValueError(f"{label} cannot contain ASCII control characters.")
    return result


def _liana_plot_labels(values: Sequence[str] | None, *, label: str) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{label} must be a sequence of exact labels, not a comma-delimited string.")
    labels = [_liana_plot_text(value, label=label) for value in values]
    if len(labels) != len(set(labels)):
        raise ValueError(f"{label} cannot contain duplicate labels.")
    return labels


def _liana_plot_number(
    value: Any,
    *,
    label: str,
    minimum: float,
    maximum: float,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a real number.")
    result = float(value)
    if not math.isfinite(result) or result < minimum or result > maximum:
        raise ValueError(f"{label} must be finite and lie in [{minimum}, {maximum}].")
    return result


def _liana_plot_integer(value: Any, *, label: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer.")
    if value < minimum or value > maximum:
        raise ValueError(f"{label} must lie in [{minimum}, {maximum}].")
    return value


def _liana_plot_signature(callable_object: Any) -> None:
    if not callable(callable_object):
        raise RuntimeError("LIANA 1.9.0 public pl.dotplot_by_sample is not callable.")
    try:
        observed = tuple(inspect.signature(callable_object).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Cannot inspect LIANA 1.9.0 public pl.dotplot_by_sample signature.") from exc
    if observed != LIANA_DOT_PLOT_SIGNATURE:
        raise RuntimeError(
            "LIANA 1.9.0 public pl.dotplot_by_sample signature drifted: expected "
            f"{LIANA_DOT_PLOT_SIGNATURE!r}, observed {observed!r}."
        )


def _liana_plot_runtime(liana_module: Any | None) -> tuple[Any, Any]:
    injected = liana_module is not None
    if liana_module is None:
        try:
            liana_module = importlib.import_module("liana")
        except (ImportError, OSError) as exc:
            raise RuntimeError(
                "LIANA Dot Plot requires optional liana==1.9.0 with pandas<3 and plotnine ({exc})."
            ) from exc
    version = getattr(liana_module, "__version__", None)
    if version != LIANA_AUDITED_VERSION:
        raise RuntimeError(f"LIANA Dot Plot requires exact liana=={LIANA_AUDITED_VERSION}; observed {version!r}.")
    pandas = importlib.import_module("pandas")
    if not injected:
        pandas_version = getattr(pandas, "__version__", "unknown")
        try:
            pandas_major = int(str(pandas_version).split(".", 1)[0])
        except ValueError as exc:
            raise RuntimeError(f"Cannot interpret pandas version {pandas_version!r} for LIANA 1.9.0.") from exc
        if pandas_major >= 3:
            raise RuntimeError(
                f"LIANA 1.9.0 declares pandas<3 and this audited plot adapter rejects pandas {pandas_version}."
            )
    renderer = getattr(getattr(liana_module, "pl", None), "dotplot_by_sample", None)
    _liana_plot_signature(renderer)
    return liana_module, renderer


def _liana_select_plot_rows(
    table: Any,
    *,
    method: str,
    source_labels: Sequence[str] | None,
    target_labels: Sequence[str] | None,
    selection_threshold: float,
    top_n: int,
    max_plot_rows: int,
    numpy: Any,
    pandas: Any,
) -> tuple[Any, dict[str, Any]]:
    source_labels = _liana_plot_labels(source_labels, label="LIANA source label")
    target_labels = _liana_plot_labels(target_labels, label="LIANA target label")
    selection_threshold = _liana_plot_number(
        selection_threshold,
        label="LIANA method-specific selection threshold",
        minimum=0.0,
        maximum=1.0,
    )
    top_n = _liana_plot_integer(top_n, label="LIANA top_n", minimum=1, maximum=2**31 - 1)
    max_plot_rows = _liana_plot_integer(
        max_plot_rows,
        label="LIANA max_plot_rows",
        minimum=1,
        maximum=2**31 - 1,
    )
    if method not in LIANA_METHOD_PROFILES:
        raise ValueError(f"Unsupported LIANA plot method: {method!r}.")
    _liana_result_table_identity(table, method=method, numpy=numpy, pandas=pandas)
    available_sources = list(dict.fromkeys(table["source"].tolist()))
    available_targets = list(dict.fromkeys(table["target"].tolist()))
    missing_sources = [label for label in source_labels if label not in set(available_sources)]
    missing_targets = [label for label in target_labels if label not in set(available_targets)]
    if missing_sources:
        raise ValueError(f"LIANA source labels are absent from the typed result: {missing_sources!r}.")
    if missing_targets:
        raise ValueError(f"LIANA target labels are absent from the typed result: {missing_targets!r}.")
    work = table.copy(deep=True)
    input_rows = len(work)
    if source_labels:
        work = work.loc[work["source"].isin(source_labels)].copy()
    if target_labels:
        work = work.loc[work["target"].isin(target_labels)].copy()
    rows_after_labels = len(work)
    specificity_field = str(LIANA_METHOD_PROFILES[method]["specificity_field"])
    work = work.loc[work[specificity_field] <= selection_threshold].copy()
    rows_after_threshold = len(work)
    if work.empty:
        raise ValueError("LIANA visual selection is empty after exact label and method-specific threshold filters.")
    sort_columns = ["sample", "source", "target"]
    ascending = [True, True, True]
    for column, direction in LIANA_METHOD_PROFILES[method]["sort"]:
        sort_columns.append(column)
        ascending.append(direction)
    sort_columns.extend(["ligand_complex", "receptor_complex"])
    ascending.extend([True, True])
    work = work.sort_values(sort_columns, ascending=ascending, kind="mergesort")
    work = (
        work.groupby(["sample", "source", "target"], sort=False, observed=True, group_keys=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    if work.empty:
        raise ValueError("LIANA visual top-N selection unexpectedly produced no rows.")
    if len(work) > max_plot_rows:
        raise ValueError(f"LIANA visual selection has {len(work):,} rows, exceeding max_plot_rows={max_plot_rows:,}.")
    selected_samples = list(dict.fromkeys(work["sample"].tolist()))
    selection = {
        "input_rows": input_rows,
        "rows_after_label_filters": rows_after_labels,
        "rows_after_method_threshold": rows_after_threshold,
        "plotted_rows": len(work),
        "source_labels": source_labels,
        "target_labels": target_labels,
        "threshold_field": specificity_field,
        "threshold_operator": "<=",
        "threshold_value": selection_threshold,
        "top_n": top_n,
        "top_n_scope": "within each Sample/source/target tuple after threshold filtering",
        "selected_samples": selected_samples,
        "selected_sample_count": len(selected_samples),
        "interaction_count": int(work[["ligand_complex", "receptor_complex"]].drop_duplicates().shape[0]),
    }
    return work, selection


def _liana_plot_package_version(package: str) -> str:
    try:
        return importlib_metadata.version(package)
    except importlib_metadata.PackageNotFoundError:
        return "not-installed"


def _liana_png(figure: Any, *, max_image_pixels: int) -> bytes:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from PIL import Image

    if not isinstance(figure, Figure):
        raise RuntimeError("LIANA plotnine draw() did not return a Matplotlib Figure.")
    original_canvas = figure.canvas
    buffer = bytes_io.BytesIO()
    try:
        FigureCanvasAgg(figure)
        figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    finally:
        if original_canvas is not None:
            figure.set_canvas(original_canvas)
    png = buffer.getvalue()
    if not png.startswith(LIANA_PNG_SIGNATURE):
        raise RuntimeError("LIANA renderer did not produce a valid PNG signature.")
    with Image.open(bytes_io.BytesIO(png)) as image:
        width, height = image.size
        if width < 1 or height < 1 or width * height > max_image_pixels:
            raise MemoryError(
                f"LIANA rendered image has {width * height:,} pixels, exceeding max_image_pixels={max_image_pixels:,}."
            )
        if image.format != "PNG":
            raise RuntimeError("LIANA renderer returned bytes that Pillow does not identify as PNG.")
        image.verify()
    return png


def _liana_plot_references(method: str, provenance: Mapping[str, Any]) -> list[dict[str, Any]]:
    references = [
        {
            "citation": (
                "Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference "
                "from single-cell RNA-Seq data. Nature Communications. 2022;13:3224."
            ),
            "doi": "10.1038/s41467-022-30755-0",
            "url": "https://doi.org/10.1038/s41467-022-30755-0",
            "kind": "method",
        },
        {
            "citation": (
                "Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication "
                "inference. Nature Cell Biology. 2024;26:1613-1622."
            ),
            "doi": "10.1038/s41556-024-01469-w",
            "url": "https://doi.org/10.1038/s41556-024-01469-w",
            "kind": "software",
        },
    ]
    if method == "rank_aggregate":
        references.append(
            {
                "citation": (
                    "Kolde R, et al. Robust rank aggregation for gene list integration and meta-analysis. "
                    "Bioinformatics. 2012;28:573-580."
                ),
                "doi": "10.1093/bioinformatics/btr709",
                "url": "https://doi.org/10.1093/bioinformatics/btr709",
                "kind": "method",
            }
        )
    else:
        references.append(
            {
                "citation": (
                    "Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression "
                    "of multi-subunit ligand-receptor complexes. Nature Protocols. 2020;15:1484-1506."
                ),
                "doi": "10.1038/s41596-020-0292-x",
                "url": "https://doi.org/10.1038/s41596-020-0292-x",
                "kind": "method",
            }
        )
    metadata = provenance["resource"]["metadata"]
    references.append(
        {
            "citation": metadata["citation"],
            "doi": None,
            "url": metadata["download_url"],
            "kind": "resource",
        }
    )
    return references


def _liana_plot_impl(
    table: Any,
    provenance: Mapping[str, Any],
    *,
    source_labels: Sequence[str] | None,
    target_labels: Sequence[str] | None,
    selection_threshold: float,
    top_n: int,
    figure_width: float,
    figure_height: float,
    max_plot_rows: int,
    max_image_pixels: int,
    openbio_version: str,
    liana_module: Any | None = None,
) -> tuple[bytes, dict[str, Any]]:
    expected_provenance = {
        "schema_version",
        "method",
        "method_profile",
        "roles",
        "organism",
        "expression",
        "design",
        "resource",
        "backend",
        "parameters",
        "result",
        "summary_sha256",
    }
    if not isinstance(provenance, Mapping) or set(provenance) != expected_provenance:
        raise ValueError("LIANA typed result provenance schema is invalid for plotting.")
    if provenance["schema_version"] != 1:
        raise ValueError("LIANA typed result provenance schema version is unsupported.")
    method = provenance["method"]
    if method not in LIANA_METHOD_PROFILES:
        raise ValueError("LIANA typed result carries an unsupported method.")
    if provenance["backend"].get("version") != LIANA_AUDITED_VERSION:
        raise ValueError("LIANA typed result was not produced by the audited LIANA 1.9.0 backend.")
    if provenance["method_profile"] != {
        key: LIANA_METHOD_PROFILES[method][key]
        for key in (
            "magnitude_field",
            "magnitude_direction",
            "specificity_field",
            "specificity_direction",
        )
    }:
        raise ValueError("LIANA typed result method score semantics do not match the audited profile.")
    openbio_version = _liana_plot_text(openbio_version, label="OpenBio version")
    figure_width = _liana_plot_number(figure_width, label="LIANA figure_width", minimum=1.0, maximum=30.0)
    figure_height = _liana_plot_number(figure_height, label="LIANA figure_height", minimum=1.0, maximum=30.0)
    max_image_pixels = _liana_plot_integer(
        max_image_pixels,
        label="LIANA max_image_pixels",
        minimum=1,
        maximum=2**31 - 1,
    )
    requested_pixels = math.ceil(figure_width * 120) * math.ceil(figure_height * 120)
    if requested_pixels > max_image_pixels:
        raise MemoryError(
            f"LIANA requested image needs at least {requested_pixels:,} nominal pixels, exceeding "
            f"max_image_pixels={max_image_pixels:,}."
        )
    numpy = importlib.import_module("numpy")
    pandas = importlib.import_module("pandas")
    matplotlib = importlib.import_module("matplotlib")
    pyplot = importlib.import_module("matplotlib.pyplot")
    liana_module, renderer = _liana_plot_runtime(liana_module)
    selected, selection = _liana_select_plot_rows(
        table,
        method=method,
        source_labels=source_labels,
        target_labels=target_labels,
        selection_threshold=selection_threshold,
        top_n=top_n,
        max_plot_rows=max_plot_rows,
        numpy=numpy,
        pandas=pandas,
    )
    selected_identity = _liana_result_table_identity(selected, method=method, numpy=numpy, pandas=pandas)
    render_frame = selected.copy(deep=True)
    render_identity = _liana_result_table_identity(render_frame, method=method, numpy=numpy, pandas=pandas)
    profile = LIANA_METHOD_PROFILES[method]
    inverse_colour = method == "rank_aggregate"
    inverse_size = True
    existing_figures = set(pyplot.get_fignums())
    rc_before = dict(matplotlib.rcParams)
    rng_before = copy.deepcopy(numpy.random.get_state())
    figure = None
    try:
        plot = renderer(
            adata=None,
            uns_key="liana_res",
            liana_res=render_frame,
            sample_key="sample",
            colour=profile["magnitude_field"],
            size=profile["specificity_field"],
            inverse_colour=inverse_colour,
            inverse_size=inverse_size,
            source_labels=list(source_labels) if source_labels else None,
            target_labels=list(target_labels) if target_labels else None,
            ligand_complex=None,
            receptor_complex=None,
            size_range=(2, 9),
            cmap="viridis",
            figure_size=(figure_width, figure_height),
            return_fig=True,
        )
        if not callable(getattr(plot, "draw", None)):
            raise RuntimeError("LIANA 1.9.0 dotplot_by_sample did not return the audited plotnine draw object.")
        figure = plot.draw()
        png = _liana_png(figure, max_image_pixels=max_image_pixels)
    finally:
        numpy.random.set_state(rng_before)
        matplotlib.rcParams.update(rc_before)
        for figure_number in set(pyplot.get_fignums()) - existing_figures:
            pyplot.close(figure_number)
        if figure is not None:
            pyplot.close(figure)
    if _liana_result_table_identity(render_frame, method=method, numpy=numpy, pandas=pandas) != render_identity:
        raise RuntimeError("LIANA plotting backend modified its private result-table copy.")
    if _liana_result_table_identity(selected, method=method, numpy=numpy, pandas=pandas) != selected_identity:
        raise RuntimeError("LIANA visual selection changed during rendering.")
    warnings_list: list[str] = []
    all_samples = list(provenance["design"]["sample_order"])
    omitted_samples = [sample for sample in all_samples if sample not in selection["selected_samples"]]
    if omitted_samples:
        warnings_list.append(
            f"Visual filters omitted {len(omitted_samples):,} Samples with no selected rows: {omitted_samples!r}."
        )
    annotation_status = provenance["roles"]["annotation_status"]
    if annotation_status != "curated":
        warnings_list.append(
            f"Identity labels are declared {annotation_status}; the figure requires annotation review."
        )
    parameters = {
        "method": method,
        "source_labels": list(source_labels or []),
        "target_labels": list(target_labels or []),
        "selection_threshold": float(selection_threshold),
        "selection_field": profile["specificity_field"],
        "top_n": top_n,
        "top_n_scope": selection["top_n_scope"],
        "figure_width": figure_width,
        "figure_height": figure_height,
        "max_plot_rows": max_plot_rows,
        "max_image_pixels": max_image_pixels,
        "colour_field": profile["magnitude_field"],
        "size_field": profile["specificity_field"],
        "inverse_colour": inverse_colour,
        "inverse_size": inverse_size,
        "sample_key": "sample",
        "return_fig": True,
        "renderer": "plotnine ggplot.draw to Matplotlib PNG",
    }
    key_columns = [
        "sample",
        "condition",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        profile["magnitude_field"],
        profile["specificity_field"],
    ]
    key_results = {
        "selection": selection,
        "visual_encoding": {
            "colour_field": profile["magnitude_field"],
            "colour_direction": profile["magnitude_direction"],
            "colour_inverted_for_display": inverse_colour,
            "size_field": profile["specificity_field"],
            "size_direction": profile["specificity_direction"],
            "size_inverted_for_display": inverse_size,
        },
        "plotted_interaction_limit": LIANA_PLOT_REPORT_LIMIT,
        "plotted_interactions": selected[key_columns].head(LIANA_PLOT_REPORT_LIMIT).to_dict(orient="records"),
        "resource": copy.deepcopy(provenance["resource"]),
        "producer_summary_sha256": provenance["summary_sha256"],
        "png_bytes": len(png),
        "nominal_pixels": requested_pixels,
    }
    limitations = [
        "This is a deterministic view of an existing typed LIANA result and performs no new communication inference.",
        "Within-Sample cell-label permutation p-values and aggregate ranks are not replicate-aware Condition tests.",
        "Point size, colour, and Sample facets do not establish secretion, binding, spatial contact, activation, "
        "causality, or differential communication.",
        "Visual top-N selection is descriptive and must not be treated as a multiple-testing procedure.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": LIANA_DOT_PLOT_NODE_ID,
        "status": "ok",
        "methods": (
            f"Validated the tamper-evident {method} result, applied the disclosed method-specific threshold and "
            "stable per-Sample/source/target top-N selection, then rendered LIANA 1.9.0 "
            "pl.dotplot_by_sample without rerunning inference."
        ),
        "results": (
            f"Rendered {selection['plotted_rows']:,} selected candidate-interaction rows across "
            f"{selection['selected_sample_count']:,} Sample facets."
        ),
        "key_results": key_results,
        "parameters": parameters,
        "references": _liana_plot_references(method, provenance),
        "software_versions": {
            "python": platform.python_version(),
            "openbio-singlecell": openbio_version,
            "liana": str(liana_module.__version__),
            "plotnine": _liana_plot_package_version("plotnine"),
            "matplotlib": str(matplotlib.__version__),
            "pandas": str(pandas.__version__),
            "numpy": str(numpy.__version__),
            "Pillow": _liana_plot_package_version("Pillow"),
        },
        "warnings": warnings_list,
        "limitations": limitations,
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return png, summary


def render_liana_dot_plot(
    result: Any,
    *,
    exact_type: bool = True,
    **parameters: Any,
) -> tuple[bytes, dict[str, Any]]:
    table, provenance, metadata = validate_liana_result(
        result,
        exact_type=exact_type,
        copy_result=False,
    )
    fingerprint = metadata["artifact_fingerprint_sha256"]
    png, summary = _liana_plot_impl(table, provenance, **parameters)
    _table_after, _provenance_after, metadata_after = validate_liana_result(
        result,
        exact_type=exact_type,
        copy_result=False,
    )
    if metadata_after["artifact_fingerprint_sha256"] != fingerprint:
        raise RuntimeError("LIANA typed result changed during plotting.")
    return png, summary


def liana_dot_plot_code(*, parameters: Mapping[str, Any]) -> str:
    if not isinstance(parameters, Mapping):
        raise TypeError("LIANA dot-plot generated-code parameters must be a mapping.")
    payload = copy.deepcopy(dict(parameters))
    json.dumps(payload, ensure_ascii=False, allow_nan=False)
    helpers = (
        _liana_plot_text,
        _liana_plot_labels,
        _liana_plot_number,
        _liana_plot_integer,
        _liana_plot_signature,
        _liana_plot_runtime,
        _liana_select_plot_rows,
        _liana_plot_package_version,
        _liana_png,
        _liana_plot_references,
        _liana_plot_impl,
    )
    helper_source = "\n\n".join(inspect.getsource(helper) for helper in helpers)
    portable_result_source = liana_result_portable_code().replace("from __future__ import annotations\n\n", "", 1)
    code = f'''from __future__ import annotations

import copy
import importlib
import inspect
import io as bytes_io
import json
import math
import platform
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from typing import Any

{portable_result_source}

LIANA_AUDITED_VERSION = {LIANA_AUDITED_VERSION!r}
LIANA_METHOD_PROFILES = {LIANA_METHOD_PROFILES!r}
LIANA_DOT_PLOT_NODE_ID = {LIANA_DOT_PLOT_NODE_ID!r}
LIANA_DOT_PLOT_SIGNATURE = {LIANA_DOT_PLOT_SIGNATURE!r}
LIANA_PNG_SIGNATURE = {LIANA_PNG_SIGNATURE!r}
LIANA_PLOT_REPORT_LIMIT = {LIANA_PLOT_REPORT_LIMIT!r}

{helper_source}

_OPENBIO_LIANA_PLOT_PARAMETERS = {payload!r}

def run_liana_dot_plot(portable_result, liana_module=None):
    """Render a portable typed LIANA result without rerunning communication inference."""
    table, provenance, _metadata = validate_liana_result(portable_result, exact_type=False)
    parameters = copy.deepcopy(_OPENBIO_LIANA_PLOT_PARAMETERS)
    parameters["liana_module"] = liana_module
    return _liana_plot_impl(table, provenance, **parameters)
'''
    compile(code, "<openbio_liana_dot_plot_code>", "exec")
    return code


__all__ = [
    "LIANA_DOT_PLOT_NODE_ID",
    "liana_dot_plot_code",
    "render_liana_dot_plot",
]
