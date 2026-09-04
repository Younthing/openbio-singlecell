from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .analysis_reporting import _plain_json
from .differential_evidence import differential_table_content_fingerprint
from .pseudobulk import _standalone_artifact_fingerprint, _standalone_validate_pseudobulk_artifact


def _standalone_pseudobulk_qc_plot(artifact, *, view=None, _return_diagnostics=False):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Pseudobulk QC Plot"
    if view is None:
        mode = {"view": "profile_qc"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") != "profile_qc":
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")

    adata, metadata = _standalone_validate_pseudobulk_artifact(artifact, copy_result=True)
    if adata.n_obs > 100_000:
        raise ValueError(f"{operation} supports at most 100,000 retained profiles in one static PNG.")
    roles = metadata["role_keys"]
    population_key = roles["population"]
    profile_records = metadata["profile_records"]
    samples = list(dict.fromkeys(record["sample"] for record in profile_records))
    populations = list(dict.fromkeys(record["population"] for record in profile_records))
    conditions = list(dict.fromkeys(record["condition"] for record in profile_records))
    sample_positions = {value: index for index, value in enumerate(samples)}
    population_positions = {value: index for index, value in enumerate(populations)}

    retained_records = [record for record in profile_records if record["retained"]]
    removed_records = [record for record in profile_records if not record["retained"]]
    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "samples": samples,
        "populations": populations,
        "conditions": conditions,
        "plotted_profiles": int(adata.n_obs),
        "removed_profiles": len(removed_records),
        "plotted_samples": len(samples),
        "plotted_genes": int(adata.n_vars),
        "count_source": metadata["count_source"],
        "inference_status": metadata["inference_status"],
    }

    if len(samples) > 200 or len(populations) > 100:
        raise ValueError(
            f"{operation} profile_qc supports at most 200 Samples and 100 populations in one static PNG."
        )
    coverage = np.zeros((len(samples), len(populations)), dtype=float)
    for record in retained_records:
        coverage[sample_positions[record["sample"]], population_positions[record["population"]]] = int(
            record["n_cells"]
        )
    figure = Figure(
        figsize=(max(9.0, min(20.0, 0.45 * len(populations) + 7.0)), 5.5),
        constrained_layout=True,
    )
    FigureCanvasAgg(figure)
    scatter_axis, coverage_axis = figure.subplots(1, 2)
    colors = [f"C{index % 10}" for index in range(len(populations))]
    for population, color in zip(populations, colors, strict=True):
        selected = adata.obs[population_key].astype(str).to_numpy() == population
        scatter_axis.scatter(
            adata.obs.loc[selected, "openbio_n_cells"].to_numpy(dtype=float),
            adata.obs.loc[selected, "openbio_total_counts"].to_numpy(dtype=float),
            label=population,
            color=color,
            alpha=0.8,
        )
    scatter_axis.set_xlabel("Contributing cells per Sample–population profile")
    scatter_axis.set_ylabel("Raw-count library size")
    scatter_axis.set_yscale("log")
    scatter_axis.set_title("Retained pseudobulk profiles")
    scatter_axis.grid(alpha=0.2)
    scatter_axis.legend(title="Population")
    image = coverage_axis.imshow(coverage, aspect="auto", interpolation="nearest", cmap="Blues")
    coverage_axis.set_xticks(np.arange(len(populations)), populations, rotation=45, ha="right")
    coverage_axis.set_yticks(np.arange(len(samples)), samples)
    coverage_axis.set_xlabel("Population")
    coverage_axis.set_ylabel("Sample")
    coverage_axis.set_title("Contributing-cell coverage")
    figure.colorbar(image, ax=coverage_axis, label="Cells")
    title = "Pseudobulk profile quality and coverage"
    details["profile_cell_counts"] = [int(value) for value in adata.obs["openbio_n_cells"]]
    details["profile_library_sizes"] = [int(value) for value in adata.obs["openbio_total_counts"]]
    details["coverage_cells"] = coverage.astype(int).tolist()

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details["title"] = title
    if _return_diagnostics:
        return png, details
    return png


def run_pseudobulk_qc_plot(artifact: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_pseudobulk_qc_plot(artifact, view=view, _return_diagnostics=True)


def pseudobulk_qc_plot_code(*, view: dict[str, Any]) -> str:
    functions = (
        _plain_json,
        _standalone_artifact_fingerprint,
        _standalone_validate_pseudobulk_artifact,
        _standalone_pseudobulk_qc_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

{implementation}


def plot_pseudobulk_qc(pseudobulk):
    return _standalone_pseudobulk_qc_plot(pseudobulk, view={view!r})
"""


def _standalone_pseudobulk_condition_contrast_plot(
    table,
    *,
    producer,
    view=None,
    _return_diagnostics=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Pseudobulk Condition Contrast Plot"
    if view is None:
        mode = {"view": "volcano"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant in {"volcano", "ma"}:
        expected_view = {"view"}
    elif variant == "top_genes":
        expected_view = {"view", "n_genes"}
    else:
        expected_view = set()
    if set(mode) != expected_view:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if variant == "top_genes" and (
        isinstance(mode["n_genes"], (bool, np.bool_))
        or not isinstance(mode["n_genes"], (int, np.integer))
        or not 1 <= int(mode["n_genes"]) <= 100
    ):
        raise ValueError(f"{operation} n_genes must be an integer from 1 through 100.")

    expected_producer = {"operation", "parameters", "input_cells", "input_genes", "random_seed"}
    if not isinstance(producer, Mapping) or set(producer) != expected_producer:
        raise ValueError(f"{operation} requires exact pseudobulk contrast producer provenance.")
    producer_operation = producer["operation"]
    if producer_operation not in {"pseudobulk_edger", "pseudobulk_deseq2"}:
        raise ValueError(f"{operation} requires an edgeR or PyDESeq2 pseudobulk contrast table.")
    parameters = producer["parameters"]
    required_parameters = {
        "engine",
        "population",
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "fdr_threshold",
        "min_abs_log2_fold_change",
        "table_content_fingerprint_sha256",
    }
    if not isinstance(parameters, Mapping) or not required_parameters.issubset(parameters):
        raise ValueError(f"{operation} producer provenance is missing the Sample-level contrast definition.")
    expected_engine = "edgeR" if producer_operation == "pseudobulk_edger" else "PyDESeq2"
    if parameters["engine"] != expected_engine:
        raise ValueError(f"{operation} producer operation and engine identity disagree.")
    for name in ("population", "condition_key", "reference_condition", "comparison_condition"):
        value = parameters[name]
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{operation} producer {name} must be canonical nonblank text.")
    if parameters["reference_condition"] == parameters["comparison_condition"]:
        raise ValueError(f"{operation} producer Conditions must differ.")
    for name in ("input_cells", "input_genes", "random_seed"):
        value = producer[name]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) < 0:
            raise ValueError(f"{operation} producer {name} must be a nonnegative integer.")
    fdr = parameters["fdr_threshold"]
    minimum_effect = parameters["min_abs_log2_fold_change"]
    if (
        isinstance(fdr, (bool, np.bool_))
        or not isinstance(fdr, (int, float, np.integer, np.floating))
        or not np.isfinite(fdr)
        or not 0 <= float(fdr) <= 1
    ):
        raise ValueError(f"{operation} producer FDR threshold is invalid.")
    if (
        isinstance(minimum_effect, (bool, np.bool_))
        or not isinstance(minimum_effect, (int, float, np.integer, np.floating))
        or not np.isfinite(minimum_effect)
        or float(minimum_effect) < 0
    ):
        raise ValueError(f"{operation} producer minimum absolute log2 fold change is invalid.")

    edger_columns = [
        "gene",
        "log2_fold_change",
        "log_counts_per_million",
        "quasi_likelihood_f",
        "p_value",
        "p_adjusted",
    ]
    deseq2_columns = [
        "gene",
        "base_mean",
        "log2_fold_change",
        "log2_fold_change_standard_error",
        "wald_statistic",
        "p_value",
        "p_adjusted",
    ]
    columns = edger_columns if producer_operation == "pseudobulk_edger" else deseq2_columns
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} table does not match the canonical {expected_engine} result schema.")
    if table.empty or len(table) > 500_000:
        raise ValueError(f"{operation} requires 1 through 500,000 tested-gene rows.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} table must use the canonical zero-based row index.")
    genes = table["gene"].tolist()
    if (
        any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes)
        or len(genes) != len(set(genes))
    ):
        raise ValueError(f"{operation} gene axis must be unique canonical nonblank text.")
    numeric = {}
    for column in columns[1:]:
        if pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_numeric_dtype(table[column].dtype):
            raise TypeError(f"{operation} column {column!r} must be numeric and non-Boolean.")
        numeric[column] = table[column].to_numpy(dtype=float)
    if producer_operation == "pseudobulk_edger":
        if any(not bool(np.isfinite(values).all()) for values in numeric.values()):
            raise ValueError(f"{operation} edgeR columns must contain finite values.")
        if bool((numeric["quasi_likelihood_f"] < 0).any()):
            raise ValueError(f"{operation} edgeR quasi-likelihood F values cannot be negative.")
    else:
        if not bool(np.isfinite(numeric["base_mean"]).all()) or bool((numeric["base_mean"] < 0).any()):
            raise ValueError(f"{operation} PyDESeq2 base means must be finite and nonnegative.")
        if any(bool(np.isinf(numeric[column]).any()) for column in columns[2:]):
            raise ValueError(f"{operation} PyDESeq2 statistical estimates cannot contain infinity.")
        finite_standard_errors = numeric["log2_fold_change_standard_error"][
            np.isfinite(numeric["log2_fold_change_standard_error"])
        ]
        if bool((finite_standard_errors < 0).any()):
            raise ValueError(f"{operation} PyDESeq2 finite standard errors cannot be negative.")
    for column in ("p_value", "p_adjusted"):
        finite = numeric[column][np.isfinite(numeric[column])]
        if bool(((finite < 0) | (finite > 1)).any()):
            raise ValueError(f"{operation} finite {column} values must lie in [0, 1].")
    adjusted_mask = np.isfinite(numeric["p_adjusted"])
    if bool((adjusted_mask & ~np.isfinite(numeric["p_value"])).any()):
        raise ValueError(f"{operation} cannot retain adjusted p-values where raw p-values are null.")

    def bh_adjust(values):
        order = np.argsort(values, kind="mergesort")
        ranked = values[order] * float(values.size) / np.arange(1, values.size + 1, dtype=float)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adjusted = np.empty_like(ranked)
        adjusted[order] = np.clip(ranked, 0.0, 1.0)
        return adjusted

    observed_adjusted = numeric["p_adjusted"][adjusted_mask]
    expected_adjusted = bh_adjust(numeric["p_value"][adjusted_mask])
    if not bool(np.allclose(observed_adjusted, expected_adjusted, rtol=1e-8, atol=1e-12)):
        raise ValueError(f"{operation} adjusted p-values disagree with the retained Benjamini–Hochberg universe.")
    content_fingerprint = differential_table_content_fingerprint(table, contract=producer_operation)
    if parameters["table_content_fingerprint_sha256"] != content_fingerprint:
        raise ValueError(f"{operation} table failed its producer-generated current-content fingerprint check.")

    effects = numeric["log2_fold_change"]
    valid = adjusted_mask & np.isfinite(effects)
    calls = valid & (numeric["p_adjusted"] <= float(fdr)) & (np.abs(effects) >= float(minimum_effect))
    up = calls & (effects > 0)
    down = calls & (effects < 0)
    colors = np.full(len(table), "#B8B8B8", dtype=object)
    colors[up] = "#D64A4A"
    colors[down] = "#3D6FB6"
    if variant == "volcano":
        selected = valid
        x = effects[selected]
        y = -np.log10(np.maximum(numeric["p_adjusted"][selected], 1e-300))
        figure = Figure(figsize=(7.2, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.scatter(x, y, c=colors[selected], alpha=0.75, s=18, linewidths=0)
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.axhline(-np.log10(max(float(fdr), 1e-300)), color="#555555", linestyle="--", linewidth=0.9)
        if float(minimum_effect) > 0:
            axis.axvline(float(minimum_effect), color="#777777", linestyle="--", linewidth=0.8)
            axis.axvline(-float(minimum_effect), color="#777777", linestyle="--", linewidth=0.8)
        axis.set_xlabel(f"log2 fold change ({parameters['comparison_condition']} − {parameters['reference_condition']})")
        axis.set_ylabel("−log10(Benjamini–Hochberg adjusted p-value)")
        axis.set_title(f"{expected_engine} Sample-level Condition contrast")
        axis.grid(alpha=0.2)
        title = f"{expected_engine} pseudobulk Condition contrast"
    elif variant == "ma":
        selected = valid
        abundance = (
            numeric["log_counts_per_million"]
            if producer_operation == "pseudobulk_edger"
            else np.log2(numeric["base_mean"] + 1.0)
        )
        selected &= np.isfinite(abundance)
        figure = Figure(figsize=(7.2, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.scatter(abundance[selected], effects[selected], c=colors[selected], alpha=0.75, s=18, linewidths=0)
        axis.axhline(0.0, color="#555555", linewidth=0.8)
        if float(minimum_effect) > 0:
            axis.axhline(float(minimum_effect), color="#777777", linestyle="--", linewidth=0.8)
            axis.axhline(-float(minimum_effect), color="#777777", linestyle="--", linewidth=0.8)
        axis.set_xlabel("log counts per million" if producer_operation == "pseudobulk_edger" else "log2(base mean + 1)")
        axis.set_ylabel(f"log2 fold change ({parameters['comparison_condition']} − {parameters['reference_condition']})")
        axis.set_title(f"{expected_engine} mean–difference evidence")
        axis.grid(alpha=0.2)
        title = f"{expected_engine} pseudobulk mean–difference plot"
    else:
        finite_rows = np.flatnonzero(valid)
        order = sorted(
            finite_rows.tolist(),
            key=lambda index: (numeric["p_adjusted"][index], -abs(effects[index]), index),
        )[: int(mode["n_genes"])]
        if not order:
            raise ValueError(f"{operation} top_genes has no finite adjusted-p-value and fold-change rows.")
        selected_genes = [genes[index] for index in order]
        selected_effects = effects[order]
        figure = Figure(figsize=(8.0, max(4.5, min(20.0, 0.35 * len(order) + 2.5))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(order), dtype=float)
        axis.barh(positions, selected_effects, color=colors[order])
        axis.set_yticks(positions, selected_genes)
        axis.invert_yaxis()
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel(f"log2 fold change ({parameters['comparison_condition']} − {parameters['reference_condition']})")
        axis.set_title(f"Top {expected_engine} genes by adjusted p-value")
        axis.grid(axis="x", alpha=0.2)
        title = f"Top {expected_engine} pseudobulk contrast genes"

    color_legend = {
        "Significant: comparison higher": "#D64A4A",
        "Significant: comparison lower": "#3D6FB6",
        "Not called significant": "#B8B8B8",
    }
    for label, color in color_legend.items():
        axis.scatter([], [], color=color, s=24, linewidths=0, label=label)
    axis.legend(frameon=False, bbox_to_anchor=(1.02, 1.0), loc="upper left")
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "engine": expected_engine,
        "population": parameters["population"],
        "condition_key": parameters["condition_key"],
        "reference_condition": parameters["reference_condition"],
        "comparison_condition": parameters["comparison_condition"],
        "tested_genes": len(table),
        "plotted_genes": int(valid.sum()) if variant != "top_genes" else len(order),
        "null_adjusted_p_values": int((~adjusted_mask).sum()),
        "fdr_threshold": float(fdr),
        "minimum_absolute_log2_fold_change": float(minimum_effect),
        "significant_up": int(up.sum()),
        "significant_down": int(down.sum()),
        "color_legend": color_legend,
        "table_content_fingerprint_sha256": content_fingerprint,
        "title": title,
    }
    if variant == "top_genes":
        details["top_genes"] = [
            {
                "gene": genes[index],
                "log2_fold_change": float(effects[index]),
                "p_adjusted": float(numeric["p_adjusted"][index]),
            }
            for index in order
        ]
    if _return_diagnostics:
        return png, details
    return png


def run_pseudobulk_condition_contrast_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_pseudobulk_condition_contrast_plot(
        table,
        producer=producer,
        view=view,
        _return_diagnostics=True,
    )


def pseudobulk_condition_contrast_plot_code(
    *,
    producer: dict[str, Any],
    view: dict[str, Any],
) -> str:
    implementation = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (differential_table_content_fingerprint, _standalone_pseudobulk_condition_contrast_plot)
    )
    return f"""from __future__ import annotations

{implementation}


def plot_pseudobulk_condition_contrast(table):
    return _standalone_pseudobulk_condition_contrast_plot(table, producer={producer!r}, view={view!r})
"""


def _standalone_scvi_population_de_evidence_plot(
    table,
    *,
    producer,
    view=None,
    _return_diagnostics=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "scVI Population DE Evidence Plot"
    if view is None:
        mode = {"view": "probability_lfc"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant in {"probability_lfc", "bayes_factor"}:
        expected_view = {"view"}
    elif variant == "top_genes":
        expected_view = {"view", "n_genes"}
    else:
        expected_view = set()
    if set(mode) != expected_view:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if variant == "top_genes" and (
        isinstance(mode["n_genes"], (bool, np.bool_))
        or not isinstance(mode["n_genes"], (int, np.integer))
        or not 1 <= int(mode["n_genes"]) <= 100
    ):
        raise ValueError(f"{operation} n_genes must be an integer from 1 through 100.")

    expected_producer = {"operation", "parameters", "input_cells", "input_genes", "random_seed"}
    if not isinstance(producer, Mapping) or set(producer) != expected_producer:
        raise ValueError(f"{operation} requires exact scVI producer provenance.")
    if producer["operation"] != "scvi_model_de_evidence":
        raise ValueError(f"{operation} requires a scVI Model DE Evidence table.")
    parameters = producer["parameters"]
    required_parameters = {
        "groupby",
        "group1",
        "group2",
        "population_scope",
        "mode",
        "delta",
        "fdr_target",
        "batch_handling",
        "n_samples_overall",
        "random_seed",
        "table_content_fingerprint_sha256",
    }
    if not isinstance(parameters, Mapping) or not required_parameters.issubset(parameters):
        raise ValueError(f"{operation} producer provenance is missing the scVI comparison definition.")
    for name in ("groupby", "group1", "group2"):
        value = parameters[name]
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{operation} producer {name} must be canonical nonblank text.")
    if parameters["group1"] == parameters["group2"]:
        raise ValueError(f"{operation} producer populations must differ.")
    analysis_mode = parameters["mode"]
    if analysis_mode not in {"change", "vanilla"}:
        raise ValueError(f"{operation} producer mode is unsupported.")
    for name in ("input_cells", "input_genes", "random_seed"):
        value = producer[name]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) < 0:
            raise ValueError(f"{operation} producer {name} must be a nonnegative integer.")

    common_numeric = [
        "scale1",
        "scale2",
        "raw_mean1",
        "raw_mean2",
        "non_zeros_proportion1",
        "non_zeros_proportion2",
        "raw_normalized_mean1",
        "raw_normalized_mean2",
    ]
    change_numeric = [
        "proba_de",
        "proba_not_de",
        "bayes_factor",
        "scale1",
        "scale2",
        "pseudocounts",
        "delta",
        "lfc_mean",
        "lfc_median",
        "lfc_std",
        "lfc_min",
        "lfc_max",
        "raw_mean1",
        "raw_mean2",
        "non_zeros_proportion1",
        "non_zeros_proportion2",
        "raw_normalized_mean1",
        "raw_normalized_mean2",
    ]
    vanilla_numeric = ["proba_m1", "proba_m2", "bayes_factor", *common_numeric]
    scientific_columns = change_numeric + ["is_de_fdr"] if analysis_mode == "change" else vanilla_numeric
    columns = ["gene", "comparison", "group1", "group2", *scientific_columns]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} table does not match the canonical scVI {analysis_mode} schema.")
    if table.empty or len(table) > 500_000:
        raise ValueError(f"{operation} requires 1 through 500,000 fitted-feature rows.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} table must use the canonical zero-based row index.")
    genes = table["gene"].tolist()
    if (
        any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes)
        or len(genes) != len(set(genes))
    ):
        raise ValueError(f"{operation} fitted-feature axis must be unique canonical nonblank text.")
    comparison = f"{parameters['group1']} vs {parameters['group2']}"
    for column, expected in (
        ("comparison", comparison),
        ("group1", parameters["group1"]),
        ("group2", parameters["group2"]),
    ):
        if table[column].tolist() != [expected] * len(table):
            raise ValueError(f"{operation} table {column} values disagree with producer provenance.")
    numeric = {}
    for column in (change_numeric if analysis_mode == "change" else vanilla_numeric):
        if pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_numeric_dtype(table[column].dtype):
            raise TypeError(f"{operation} column {column!r} must be numeric and non-Boolean.")
        values = table[column].to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} column {column!r} must contain finite values.")
        numeric[column] = values
    for column in ("scale1", "scale2", "raw_mean1", "raw_mean2", "raw_normalized_mean1", "raw_normalized_mean2"):
        if bool((numeric[column] < 0).any()):
            raise ValueError(f"{operation} column {column!r} cannot contain negative values.")
    for column in ("non_zeros_proportion1", "non_zeros_proportion2"):
        if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
            raise ValueError(f"{operation} column {column!r} must lie in [0, 1].")
    if analysis_mode == "change":
        for column in ("proba_de", "proba_not_de"):
            if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
                raise ValueError(f"{operation} posterior probability {column!r} must lie in [0, 1].")
        if not bool(np.allclose(numeric["proba_de"] + numeric["proba_not_de"], 1.0, rtol=0, atol=1e-6)):
            raise ValueError(f"{operation} posterior probabilities must be complementary.")
        delta = parameters["delta"]
        fdr_target = parameters["fdr_target"]
        if (
            isinstance(delta, (bool, np.bool_))
            or not isinstance(delta, (int, float, np.integer, np.floating))
            or not np.isfinite(delta)
            or float(delta) <= 0
            or not bool(np.allclose(numeric["delta"], float(delta), rtol=0, atol=1e-12))
        ):
            raise ValueError(f"{operation} practical-change threshold disagrees with producer provenance.")
        if (
            isinstance(fdr_target, (bool, np.bool_))
            or not isinstance(fdr_target, (int, float, np.integer, np.floating))
            or not 0 < float(fdr_target) < 1
        ):
            raise ValueError(f"{operation} posterior expected-FDR target is invalid.")
        if bool((numeric["pseudocounts"] <= 0).any()) or not bool(
            np.allclose(numeric["pseudocounts"], numeric["pseudocounts"][0], rtol=0, atol=0)
        ):
            raise ValueError(f"{operation} realized pseudocount evidence is invalid.")
        if bool((numeric["lfc_std"] < 0).any()) or bool((numeric["lfc_min"] > numeric["lfc_max"]).any()):
            raise ValueError(f"{operation} posterior log2-fold-change ranges are invalid.")
        if any(
            bool(((numeric[column] < numeric["lfc_min"]) | (numeric[column] > numeric["lfc_max"])).any())
            for column in ("lfc_mean", "lfc_median")
        ):
            raise ValueError(f"{operation} posterior log2-fold-change centers lie outside stored ranges.")
        if not pd.api.types.is_bool_dtype(table["is_de_fdr"].dtype) or bool(table["is_de_fdr"].isna().any()):
            raise TypeError(f"{operation} posterior-FDR tags must be complete Boolean values.")
        if bool(np.any(numeric["proba_de"][:-1] < numeric["proba_de"][1:])):
            raise ValueError(f"{operation} change-mode rows must retain upstream posterior-probability order.")
        effect = numeric["lfc_mean"]
        probability = numeric["proba_de"]
        tagged = table["is_de_fdr"].to_numpy(dtype=bool)
        colors = np.where(tagged & (effect > 0), "#D64A4A", np.where(tagged & (effect < 0), "#3D6FB6", "#B8B8B8"))
        effect_label = f"Posterior mean decoded log2 fold change ({parameters['group1']} − {parameters['group2']})"
        probability_label = "Posterior probability of practical differential expression"
    else:
        for column in ("proba_m1", "proba_m2"):
            if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
                raise ValueError(f"{operation} directional probability {column!r} must lie in [0, 1].")
        if not bool(np.allclose(numeric["proba_m1"] + numeric["proba_m2"], 1.0, rtol=0, atol=1e-6)):
            raise ValueError(f"{operation} directional probabilities must be complementary.")
        if bool(np.any(numeric["bayes_factor"][:-1] < numeric["bayes_factor"][1:])):
            raise ValueError(f"{operation} vanilla rows must retain upstream Bayes-factor order.")
        effect = np.log2((numeric["raw_normalized_mean1"] + 1.0) / (numeric["raw_normalized_mean2"] + 1.0))
        probability = numeric["proba_m1"]
        tagged = np.zeros(len(table), dtype=bool)
        colors = np.where(probability > 0.5, "#D64A4A", np.where(probability < 0.5, "#3D6FB6", "#B8B8B8"))
        effect_label = f"Decoded normalized-mean log2 ratio ({parameters['group1']} − {parameters['group2']}; +1 display offset)"
        probability_label = f"Posterior probability {parameters['group1']} is higher"

    content_fingerprint = differential_table_content_fingerprint(
        table,
        contract=f"scvi_model_de_evidence/{analysis_mode}",
    )
    if parameters["table_content_fingerprint_sha256"] != content_fingerprint:
        raise ValueError(f"{operation} table failed its producer-generated current-content fingerprint check.")

    if variant == "probability_lfc":
        figure = Figure(figsize=(7.2, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.scatter(effect, probability, c=colors, alpha=0.75, s=18, linewidths=0)
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_ylim(0.0, 1.0)
        axis.set_xlabel(effect_label)
        axis.set_ylabel(probability_label)
        axis.set_title("scVI cell/model-conditional population evidence")
        axis.grid(alpha=0.2)
        title = "scVI population probability and effect evidence"
        plotted_rows = len(table)
    elif variant == "bayes_factor":
        figure = Figure(figsize=(7.2, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.scatter(effect, numeric["bayes_factor"], c=colors, alpha=0.75, s=18, linewidths=0)
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel(effect_label)
        axis.set_ylabel("scVI Bayes factor")
        axis.set_title("scVI Bayes-factor evidence")
        axis.grid(alpha=0.2)
        title = "scVI population Bayes-factor evidence"
        plotted_rows = len(table)
    else:
        count = min(int(mode["n_genes"]), len(table))
        selected_genes = genes[:count]
        selected_effect = effect[:count]
        figure = Figure(figsize=(8.0, max(4.5, min(20.0, 0.35 * count + 2.5))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(count, dtype=float)
        axis.barh(positions, selected_effect, color=colors[:count])
        axis.set_yticks(positions, selected_genes)
        axis.invert_yaxis()
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel(effect_label)
        axis.set_title("Top upstream-ranked scVI population evidence")
        axis.grid(axis="x", alpha=0.2)
        title = "Top scVI population evidence"
        plotted_rows = count

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "mode": analysis_mode,
        "groupby": parameters["groupby"],
        "group1": parameters["group1"],
        "group2": parameters["group2"],
        "population_scope": parameters["population_scope"],
        "batch_handling": parameters["batch_handling"],
        "tested_fitted_features": len(table),
        "plotted_features": plotted_rows,
        "posterior_fdr_tagged": int(tagged.sum()) if analysis_mode == "change" else None,
        "statistical_unit": "cells and fitted-model posterior draws",
        "table_content_fingerprint_sha256": content_fingerprint,
        "title": title,
    }
    if variant == "top_genes":
        details["top_genes"] = [
            {
                "gene": genes[index],
                "display_effect": float(effect[index]),
                "bayes_factor": float(numeric["bayes_factor"][index]),
                "posterior_probability": float(probability[index]),
                "posterior_fdr_tagged": bool(tagged[index]) if analysis_mode == "change" else None,
            }
            for index in range(plotted_rows)
        ]
    if _return_diagnostics:
        return png, details
    return png


def run_scvi_population_de_evidence_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_scvi_population_de_evidence_plot(
        table,
        producer=producer,
        view=view,
        _return_diagnostics=True,
    )


def scvi_population_de_evidence_plot_code(
    *,
    producer: dict[str, Any],
    view: dict[str, Any],
) -> str:
    implementation = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (differential_table_content_fingerprint, _standalone_scvi_population_de_evidence_plot)
    )
    return f"""from __future__ import annotations

{implementation}


def plot_scvi_population_de_evidence(table):
    return _standalone_scvi_population_de_evidence_plot(table, producer={producer!r}, view={view!r})
"""


__all__ = [
    "pseudobulk_condition_contrast_plot_code",
    "pseudobulk_qc_plot_code",
    "run_pseudobulk_condition_contrast_plot",
    "run_pseudobulk_qc_plot",
    "run_scvi_population_de_evidence_plot",
    "scvi_population_de_evidence_plot_code",
]
