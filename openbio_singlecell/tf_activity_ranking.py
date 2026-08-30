from __future__ import annotations

import copy
import importlib
import inspect
import json
import math
from collections.abc import Mapping, Sequence
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import _package_version, collect_software_versions
from .collectri_ulm import DECOUPLER_VERSION, _bh_adjust, _signature_shape
from .tf_activity_artifact import (
    TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION,
    TF_ACTIVITY_ARTIFACT_TYPE,
    TF_ACTIVITY_PRODUCER_NODE_ID,
    TF_ACTIVITY_PRODUCER_SCHEMA,
    TFActivityArtifact,
    _activity_frame_fingerprint,
    _canonical_axis,
    _canonical_json_sha256,
    _portable_tf_activity_payload,
    validate_tf_activity_artifact,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


TF_RANKING_METHODS = ("t-test_overestim_var", "t-test", "wilcoxon")
TF_RANKING_SUMMARY_SCHEMA = "openbio-singlecell/tf-activity-ranking-summary/v1"
TF_RANKING_COLUMNS = (
    "group",
    "reference",
    "regulator",
    "statistic",
    "mean_change",
    "p_value",
    "p_adjusted",
    "direction",
    "significant_adjusted",
    "rank_within_group",
)
TF_RANKING_REFERENCES = [
    {
        "citation": (
            "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological "
            "activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
        ),
        "doi": "10.1093/bioadv/vbac016",
        "url": "https://doi.org/10.1093/bioadv/vbac016",
        "kind": "software",
    },
    {
        "citation": (
            "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach "
            "to multiple testing. Journal of the Royal Statistical Society B. 1995;57:289-300."
        ),
        "doi": "10.1111/j.2517-6161.1995.tb02031.x",
        "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
        "kind": "method",
    },
    {
        "citation": "decoupler 2.2.0 get_obsm and rankby_group public interfaces and tagged implementation.",
        "doi": None,
        "url": "https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/tl/_rankby_group.py",
        "kind": "software_documentation",
    },
]


def _validate_decoupler_ranking(decoupler: Any) -> None:
    if getattr(decoupler, "__version__", None) != DECOUPLER_VERSION:
        raise RuntimeError(
            f"Rank TF Activities requires decoupler exactly {DECOUPLER_VERSION}; "
            f"observed {getattr(decoupler, '__version__', None)!r}."
        )
    try:
        get_obsm = decoupler.pp.get_obsm
        rankby_group = decoupler.tl.rankby_group
    except AttributeError as exc:
        raise RuntimeError("decoupler 2.2.0 public get_obsm/rankby_group interfaces are unavailable.") from exc
    if _signature_shape(get_obsm) != [
        ("adata", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("key", "POSITIONAL_OR_KEYWORD", "<required>"),
    ]:
        raise RuntimeError("decoupler 2.2.0 public pp.get_obsm signature changed.")
    if _signature_shape(rankby_group) != [
        ("adata", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("groupby", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("reference", "POSITIONAL_OR_KEYWORD", "rest"),
        ("method", "POSITIONAL_OR_KEYWORD", "t-test_overestim_var"),
    ]:
        raise RuntimeError("decoupler 2.2.0 public tl.rankby_group signature changed.")


def _canonical_annotation(series: Any, *, pandas: Any) -> tuple[list[str], list[str], list[str]]:
    labels: list[str] = []
    observed: list[str] = []
    identities_by_display: dict[str, tuple[str, str]] = {}
    for position, value in enumerate(series.tolist()):
        missing = pandas.isna(value)
        if not isinstance(missing, bool):
            try:
                missing = bool(missing)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"TF annotation at row {position} must be scalar.") from exc
        if missing:
            raise ValueError("TF annotation contains missing labels.")
        if isinstance(value, str):
            if not value or value != value.strip():
                raise ValueError("TF annotation string labels must be canonical and nonblank.")
            display = value
            identity = ("str", value)
        elif isinstance(value, (bool, int, float)):
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError("TF annotation numeric labels must be finite.")
            display = str(value)
            identity = (type(value).__name__, repr(value))
        else:
            raise TypeError("TF annotation labels must be scalar strings, booleans, or finite numbers.")
        previous = identities_by_display.get(display)
        if previous is not None and previous != identity:
            raise ValueError(
                f"TF annotation labels collide after display normalization at {display!r}: {previous!r} and {identity!r}."
            )
        identities_by_display[display] = identity
        labels.append(display)
        if display not in observed:
            observed.append(display)
    unused: list[str] = []
    if isinstance(series.dtype, pandas.CategoricalDtype):
        for category in series.cat.categories.tolist():
            if pandas.isna(category):
                continue
            display = str(category)
            if display not in observed:
                unused.append(display)
    if len(observed) < 2:
        raise ValueError("Rank TF Activities requires at least two observed annotation groups.")
    return labels, observed, unused


def _parse_reference(reference: str, *, observed_groups: Sequence[str]) -> tuple[str | list[str], list[str]]:
    if not isinstance(reference, str) or not reference:
        raise ValueError("TF activity reference must be a canonical nonblank string.")
    observed_set = set(observed_groups)
    if "," in reference:
        values = [value.strip() for value in reference.split(",")]
        if not values or any(not value for value in values):
            raise ValueError("TF activity reference list contains a blank entry.")
        if len(values) != len(set(values)):
            raise ValueError("TF activity reference list entries must be distinct.")
        unknown = [value for value in values if value not in observed_set]
        if unknown:
            raise ValueError(f"TF activity reference list contains unknown annotation labels: {unknown!r}.")
        if len(values) == len(observed_groups):
            raise ValueError("TF activity reference list cannot contain every observed group.")
        return values, values
    if reference != reference.strip():
        raise ValueError("TF activity reference must be a canonical nonblank string.")
    if reference == "rest":
        return "rest", []
    if reference not in observed_set:
        raise ValueError(f"TF activity reference label is not observed: {reference!r}.")
    return reference, [reference]


def _independent_comparison(
    group_values: Any,
    reference_values: Any,
    *,
    method: str,
    scipy_stats: Any,
    numpy: Any,
) -> tuple[float, float, float]:
    mean_change = float(numpy.mean(group_values) - numpy.mean(reference_values))
    if method == "wilcoxon":
        statistic, p_value = scipy_stats.ranksums(group_values, reference_values)
    else:
        statistic, p_value = scipy_stats.ttest_ind_from_stats(
            mean1=numpy.mean(group_values),
            std1=numpy.std(group_values, ddof=1),
            nobs1=group_values.size,
            mean2=numpy.mean(reference_values),
            std2=numpy.std(reference_values, ddof=1),
            nobs2=group_values.size if method == "t-test_overestim_var" else reference_values.size,
            equal_var=False,
        )
    return float(statistic), mean_change, float(p_value)


def rank_tf_activities(
    adata: AnnData,
    activities: TFActivityArtifact,
    *,
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    reference: str = "rest",
    method: str = "t-test_overestim_var",
    report_p_adjusted: float = 0.05,
    max_output_rows: int = 100_000,
    decoupler_module: Any | None = None,
    openbio_version: str = PLUGIN_VERSION,
    _portable_artifact: bool = False,
    _worker_owned: bool = False,
) -> tuple[DataFrame, dict[str, Any]]:
    """Characterize one validated TF-activity artifact across annotation groups."""

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy

    if not isinstance(annotation_key, str) or not annotation_key or annotation_key != annotation_key.strip():
        raise ValueError("TF activity annotation_key must be a canonical nonblank string.")
    if annotation_key not in adata.obs:
        raise ValueError(f"TF activity annotation column not found in obs: {annotation_key!r}.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError("TF activity annotation_status must be unknown, provisional, or curated.")
    if method not in TF_RANKING_METHODS:
        raise ValueError(f"TF activity method must be one of {list(TF_RANKING_METHODS)!r}.")
    if isinstance(report_p_adjusted, bool) or not isinstance(report_p_adjusted, (int, float)):
        raise TypeError("TF activity report_p_adjusted must be numeric.")
    report_p_adjusted = float(report_p_adjusted)
    if not math.isfinite(report_p_adjusted) or not 0.0 <= report_p_adjusted <= 1.0:
        raise ValueError("TF activity report_p_adjusted must be finite and lie in [0, 1].")
    if isinstance(max_output_rows, bool) or not isinstance(max_output_rows, int) or max_output_rows < 1:
        raise TypeError("TF activity max_output_rows must be a positive integer.")
    scores, _activity_padj, provenance, artifact_metadata = validate_tf_activity_artifact(
        activities,
        exact_type=not _portable_artifact,
        numpy=np,
        pandas=pd,
        copy_result=not _portable_artifact,
    )
    observation_names = []
    for position, value in enumerate(adata.obs_names.tolist()):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(
                f"Rank TF Activities observation identifier at position {position} must be canonical and nonblank."
            )
        observation_names.append(value)
    if len(observation_names) != len(set(observation_names)):
        raise ValueError("Rank TF Activities observation identifiers must be unique.")
    if scores.index.tolist() != observation_names:
        raise ValueError("TF activity artifact observation IDs do not align exactly to AnnData observations.")
    regulators = scores.columns.tolist()
    labels, observed_groups, unused_categories = _canonical_annotation(
        adata.obs[annotation_key], pandas=pd
    )
    backend_reference, explicit_reference_groups = _parse_reference(
        reference, observed_groups=observed_groups
    )
    explicit_set = set(explicit_reference_groups)
    tested_groups = [group for group in observed_groups if group not in explicit_set]
    if not tested_groups:
        raise ValueError("TF activity reference leaves no annotation group to test.")
    label_array = np.asarray(labels, dtype=object)
    comparison_masks: dict[str, tuple[Any, Any]] = {}
    group_support: list[dict[str, Any]] = []
    for group in tested_groups:
        group_mask = label_array == group
        if backend_reference == "rest":
            reference_mask = ~group_mask
            reference_display = "rest"
        elif isinstance(backend_reference, list):
            reference_mask = np.isin(label_array, backend_reference)
            reference_display = ", ".join(backend_reference)
        else:
            reference_mask = label_array == backend_reference
            reference_display = backend_reference
        group_size = int(group_mask.sum())
        reference_size = int(reference_mask.sum())
        if group_size < 2 or reference_size < 2:
            raise ValueError(
                f"TF activity comparison {group!r} versus {reference_display!r} requires at least two cells "
                f"on each side; observed {group_size} and {reference_size}."
            )
        comparison_masks[group] = (group_mask, reference_mask)
        group_support.append(
            {
                "group": group,
                "group_cells": group_size,
                "reference": reference_display,
                "reference_cells": reference_size,
            }
        )
    backend_groups = (
        observed_groups
        if backend_reference == "rest" or isinstance(backend_reference, list)
        else [group for group in observed_groups if group != backend_reference]
    )
    backend_rows = len(backend_groups) * len(regulators)
    canonical_rows = len(tested_groups) * len(regulators)
    if backend_rows > max_output_rows or canonical_rows > max_output_rows:
        raise ValueError(
            f"TF activity ranking requires up to {backend_rows:,} backend rows and {canonical_rows:,} canonical "
            f"rows, exceeding max_output_rows={max_output_rows:,}."
        )
    decoupler = importlib.import_module("decoupler") if decoupler_module is None else decoupler_module
    _validate_decoupler_ranking(decoupler)
    carrier = ad.AnnData(
        X=np.zeros((len(observation_names), 1), dtype=float),
        obs=pd.DataFrame({"openbio_annotation": labels}, index=observation_names),
        var=pd.DataFrame(index=["carrier_placeholder"]),
    )
    carrier.obsm["openbio_tf_activity"] = scores if _worker_owned else scores.copy(deep=True)
    carrier_x_before = carrier.X.copy()
    carrier_obs_before = carrier.obs.copy(deep=True)
    carrier_var_before = carrier.var.copy(deep=True)
    carrier_obsm_before = (
        None
        if _worker_owned
        else carrier.obsm["openbio_tf_activity"].copy(deep=True)
    )
    carrier_uns_before = copy.deepcopy(dict(carrier.uns))
    carrier_layer_keys_before = list(carrier.layers.keys())
    carrier_obsp_keys_before = list(carrier.obsp.keys())
    carrier_varm_keys_before = list(carrier.varm.keys())
    carrier_varp_keys_before = list(carrier.varp.keys())
    carrier_raw_before = carrier.raw
    extracted = decoupler.pp.get_obsm(carrier, key="openbio_tf_activity")
    if not isinstance(extracted, ad.AnnData):
        raise RuntimeError("decoupler pp.get_obsm returned a non-AnnData object.")
    if not np.array_equal(carrier.X, carrier_x_before):
        raise RuntimeError("decoupler pp.get_obsm modified private carrier X.")
    try:
        pd.testing.assert_frame_equal(carrier.obs, carrier_obs_before, check_exact=True)
        pd.testing.assert_frame_equal(carrier.var, carrier_var_before, check_exact=True)
        if carrier_obsm_before is not None:
            pd.testing.assert_frame_equal(
                carrier.obsm["openbio_tf_activity"], carrier_obsm_before, check_exact=True
            )
    except AssertionError as exc:
        raise RuntimeError("decoupler pp.get_obsm modified private carrier annotations or activities.") from exc
    if (
        dict(carrier.uns) != carrier_uns_before
        or list(carrier.layers.keys()) != carrier_layer_keys_before
        or list(carrier.obsp.keys()) != carrier_obsp_keys_before
        or list(carrier.varm.keys()) != carrier_varm_keys_before
        or list(carrier.varp.keys()) != carrier_varp_keys_before
        or carrier.raw is not carrier_raw_before
    ):
        raise RuntimeError("decoupler pp.get_obsm modified private carrier state.")
    if extracted.obs_names.tolist() != observation_names or extracted.var_names.tolist() != regulators:
        raise RuntimeError("decoupler pp.get_obsm returned misaligned activity axes.")
    if not np.array_equal(np.asarray(extracted.X, dtype=float), scores.to_numpy(dtype=float)):
        raise RuntimeError("decoupler pp.get_obsm returned activity values that differ from the typed artifact.")
    score_work = extracted if _worker_owned else extracted.copy()
    score_work.obs["openbio_annotation"] = labels
    work_x_before = None if _worker_owned else np.asarray(score_work.X, dtype=float).copy()
    work_obs_before = score_work.obs.copy(deep=True)
    work_var_before = score_work.var.copy(deep=True)
    work_uns_before = copy.deepcopy(dict(score_work.uns))
    work_layer_keys_before = list(score_work.layers.keys())
    work_obsp_keys_before = list(score_work.obsp.keys())
    work_varm_keys_before = list(score_work.varm.keys())
    work_varp_keys_before = list(score_work.varp.keys())
    work_raw_before = score_work.raw
    work_obsm_keys_before = list(score_work.obsm.keys())
    work_activity_before = (
        score_work.obsm["openbio_tf_activity"].copy()
        if not _worker_owned and "openbio_tf_activity" in score_work.obsm
        else None
    )
    ranked = decoupler.tl.rankby_group(
        adata=score_work,
        groupby="openbio_annotation",
        reference=backend_reference,
        method=method,
    )
    if work_x_before is not None and not np.array_equal(
        np.asarray(score_work.X, dtype=float), work_x_before
    ):
        raise RuntimeError("decoupler rankby_group modified the private activity matrix.")
    try:
        pd.testing.assert_frame_equal(score_work.obs, work_obs_before, check_exact=True)
        pd.testing.assert_frame_equal(score_work.var, work_var_before, check_exact=True)
        pd.testing.assert_frame_equal(carrier.obs, carrier_obs_before, check_exact=True)
        pd.testing.assert_frame_equal(carrier.var, carrier_var_before, check_exact=True)
        if carrier_obsm_before is not None:
            pd.testing.assert_frame_equal(
                carrier.obsm["openbio_tf_activity"], carrier_obsm_before, check_exact=True
            )
    except AssertionError as exc:
        raise RuntimeError("decoupler ranking modified private input annotations or carrier activities.") from exc
    if (
        dict(score_work.uns) != work_uns_before
        or list(score_work.layers.keys()) != work_layer_keys_before
        or list(score_work.obsp.keys()) != work_obsp_keys_before
        or list(score_work.varm.keys()) != work_varm_keys_before
        or list(score_work.varp.keys()) != work_varp_keys_before
        or score_work.raw is not work_raw_before
        or list(score_work.obsm.keys()) != work_obsm_keys_before
        or dict(carrier.uns) != carrier_uns_before
        or list(carrier.layers.keys()) != carrier_layer_keys_before
        or list(carrier.obsp.keys()) != carrier_obsp_keys_before
        or list(carrier.varm.keys()) != carrier_varm_keys_before
        or list(carrier.varp.keys()) != carrier_varp_keys_before
        or carrier.raw is not carrier_raw_before
    ):
        raise RuntimeError("decoupler ranking modified private input or carrier state.")
    if work_activity_before is not None:
        try:
            pd.testing.assert_frame_equal(
                score_work.obsm["openbio_tf_activity"], work_activity_before, check_exact=True
            )
        except AssertionError as exc:
            raise RuntimeError("decoupler ranking modified private activity storage.") from exc
    if not isinstance(ranked, pd.DataFrame):
        raise RuntimeError("decoupler rankby_group returned a non-DataFrame result.")
    backend_columns = ["group", "reference", "name", "stat", "meanchange", "pval", "padj"]
    if ranked.columns.tolist() != backend_columns:
        raise RuntimeError(
            f"decoupler rankby_group output schema changed: observed={ranked.columns.tolist()!r}."
        )
    backend = ranked if _worker_owned else ranked.copy(deep=True)
    for column in ("group", "reference", "name"):
        backend[column] = backend[column].astype(object)
    if bool(backend.duplicated(["group", "name"]).any()):
        raise RuntimeError("decoupler rankby_group returned duplicate group-regulator rows.")
    expected_keys = {(group, regulator) for group in backend_groups for regulator in regulators}
    observed_keys = set(backend[["group", "name"]].itertuples(index=False, name=None))
    if observed_keys != expected_keys:
        raise RuntimeError(
            "decoupler rankby_group rows do not form the complete expected group-regulator family."
        )
    expected_reference_display = (
        "rest" if backend_reference == "rest" else ", ".join(backend_reference)
        if isinstance(backend_reference, list)
        else backend_reference
    )
    if set(backend["reference"].tolist()) != {expected_reference_display}:
        raise RuntimeError("decoupler rankby_group reference labels differ from the declared comparator.")
    backend_indexed = backend.set_index(["group", "name"])
    score_values = scores.to_numpy(dtype=float, copy=not _worker_owned)
    regulator_position = {regulator: position for position, regulator in enumerate(regulators)}
    independent_by_group: dict[str, list[dict[str, Any]]] = {}
    for group in backend_groups:
        group_mask = label_array == group
        if backend_reference == "rest":
            reference_mask = ~group_mask
        elif isinstance(backend_reference, list):
            reference_mask = np.isin(label_array, backend_reference)
        else:
            reference_mask = label_array == backend_reference
        rows: list[dict[str, Any]] = []
        raw_pvalues: list[float] = []
        raw_statistics: list[tuple[float, float, float]] = []
        for regulator in regulators:
            position = regulator_position[regulator]
            values = _independent_comparison(
                score_values[group_mask, position],
                score_values[reference_mask, position],
                method=method,
                scipy_stats=scipy.stats,
                numpy=np,
            )
            if not all(math.isfinite(value) for value in values):
                raise ValueError(
                    f"TF activity comparison {group!r}/{regulator!r} is non-finite; "
                    "constant or degenerate group values are not reportable."
                )
            raw_statistics.append(values)
            raw_pvalues.append(values[2])
        adjusted_values = _bh_adjust(raw_pvalues, numpy=np)
        for regulator, (statistic, mean_change, p_value), p_adjusted in zip(
            regulators, raw_statistics, adjusted_values, strict=True
        ):
            observed = backend_indexed.loc[(group, regulator)]
            observed_values = np.asarray(
                [observed["stat"], observed["meanchange"], observed["pval"], observed["padj"]],
                dtype=float,
            )
            expected_values = np.asarray([statistic, mean_change, p_value, p_adjusted], dtype=float)
            if not bool(np.isfinite(observed_values).all()):
                raise RuntimeError("decoupler rankby_group returned non-finite canonical statistics.")
            if bool(((observed_values[2:] < 0.0) | (observed_values[2:] > 1.0)).any()):
                raise RuntimeError("decoupler rankby_group returned p-values outside [0, 1].")
            if not bool(np.allclose(observed_values, expected_values, rtol=1e-11, atol=1e-13)):
                raise RuntimeError(
                    f"decoupler rankby_group values differ from independent {method} and BH calculation "
                    f"for {group!r}/{regulator!r}."
                )
            rows.append(
                {
                    "group": group,
                    "reference": expected_reference_display,
                    "regulator": regulator,
                    "statistic": statistic,
                    "mean_change": mean_change,
                    "p_value": p_value,
                    "p_adjusted": float(p_adjusted),
                }
            )
        independent_by_group[group] = rows
    canonical_rows_data: list[dict[str, Any]] = []
    sign_discordance_count = 0
    for group in tested_groups:
        rows = independent_by_group[group]
        rows.sort(
            key=lambda row: (
                row["p_adjusted"],
                row["p_value"],
                -row["statistic"],
                row["regulator"],
            )
        )
        previous_key: tuple[float, float, float] | None = None
        previous_rank = 0
        for position, row in enumerate(rows, start=1):
            scientific_key = (row["p_adjusted"], row["p_value"], -row["statistic"])
            if scientific_key != previous_key:
                previous_rank = position
                previous_key = scientific_key
            statistic = float(row["statistic"])
            mean_change = float(row["mean_change"])
            if statistic > 0.0:
                direction = "higher_in_group"
            elif statistic < 0.0:
                direction = "lower_in_group"
            else:
                direction = "no_signed_difference"
            if statistic != 0.0 and mean_change != 0.0 and (statistic > 0.0) != (mean_change > 0.0):
                sign_discordance_count += 1
            canonical_rows_data.append(
                {
                    **row,
                    "direction": direction,
                    "significant_adjusted": bool(row["p_adjusted"] <= report_p_adjusted),
                    "rank_within_group": previous_rank,
                }
            )
    table = pd.DataFrame(canonical_rows_data, columns=TF_RANKING_COLUMNS)
    table["rank_within_group"] = table["rank_within_group"].astype(np.int64)
    if table.columns.tolist() != list(TF_RANKING_COLUMNS) or len(table) != canonical_rows:
        raise RuntimeError("TF activity canonical ranking table construction failed.")
    warnings = [
        "TF ranking treats cells as independent descriptive units; it is not biological-Sample-level Condition inference.",
        "Adjusted p-values are BH corrections across the complete regulator family within each tested annotation group.",
    ]
    if annotation_status == "unknown":
        warnings.append("Annotation status is unknown; group labels must not be presented as curated identities.")
    elif annotation_status == "provisional":
        warnings.append("Annotation labels are provisional and require independent expert review.")
    if method == "wilcoxon":
        warnings.append(
            "decoupler 2.2.0 'wilcoxon' uses SciPy ranksums, an asymptotic two-sided test without tie correction."
        )
    if method == "t-test_overestim_var":
        warnings.append(
            "t-test_overestim_var substitutes the tested-group size for reference nobs and is conservative only "
            "when this decreases the reference size."
        )
    if sign_discordance_count:
        warnings.append(
            f"Statistic and arithmetic mean-change signs disagree for {sign_discordance_count} row(s); "
            "direction follows the test statistic."
        )
    references = copy.deepcopy(TF_RANKING_REFERENCES)
    if method == "wilcoxon":
        references.append(
            {
                "citation": "Wilcoxon F. Individual comparisons by ranking methods. Biometrics Bulletin. 1945.",
                "doi": "10.2307/3001968",
                "url": "https://doi.org/10.2307/3001968",
                "kind": "method",
            }
        )
    else:
        references.append(
            {
                "citation": (
                    "Welch BL. The generalization of Student's problem when several different population "
                    "variances are involved. Biometrika. 1947;34:28-35."
                ),
                "doi": "10.1093/biomet/34.1-2.28",
                "url": "https://doi.org/10.1093/biomet/34.1-2.28",
                "kind": "method",
            }
        )
    leading_by_group = []
    for group in tested_groups:
        leading_by_group.append(
            {
                "group": group,
                "leading_regulators": table.loc[table["group"] == group].head(10).to_dict(orient="records"),
            }
        )
    significant = int(table["significant_adjusted"].sum())
    parameters = {
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "reference": reference,
        "resolved_backend_reference": backend_reference,
        "method": method,
        "report_p_adjusted": report_p_adjusted,
        "max_output_rows": max_output_rows,
    }
    summary = {
        "schema_version": TF_RANKING_SUMMARY_SCHEMA,
        "node_id": "OpenBioSingleCellRankTFActivities",
        "status": "cluster_annotation_associated_tf_activity_evidence",
        "methods": (
            f"decoupler {DECOUPLER_VERSION} rankby_group characterized the validated cell-by-regulator activity "
            f"artifact using {method}. Independent statistics and per-group complete-family BH adjustment were "
            "recomputed before canonical release."
        ),
        "results": (
            f"The complete table contains {len(table):,} regulator comparisons across {len(tested_groups):,} "
            f"tested annotation groups; {significant:,} meet adjusted p <= {report_p_adjusted:g}."
        ),
        "key_results": {
            "scientific_label": "Cluster/annotation-associated TF activity evidence",
            "observations": len(observation_names),
            "regulators": len(regulators),
            "tested_groups": tested_groups,
            "reference_groups": explicit_reference_groups,
            "group_support": group_support,
            "complete_family_rows": len(table),
            "significant_adjusted_rows": significant,
            "positive_statistic_rows": int((table["statistic"] > 0.0).sum()),
            "negative_statistic_rows": int((table["statistic"] < 0.0).sum()),
            "statistic_mean_change_sign_discordance_rows": sign_discordance_count,
            "unused_categorical_levels": unused_categories,
            "leading_by_group": leading_by_group,
            "activity_artifact_fingerprint_sha256": artifact_metadata["artifact_fingerprint_sha256"],
            "activity_provenance": provenance,
        },
        "parameters": parameters,
        "references": references,
        "software_versions": collect_software_versions(
            ["anndata", "decoupler", "numpy", "pandas", "scipy"],
            openbio_version=openbio_version,
        ),
        "warnings": warnings,
        "limitations": [
            "Cell-level tests do not use biological Sample as the replicate and cannot support a Condition claim.",
            "Ranking describes activity-score association with supplied annotations, not TF binding or causality.",
            "Results inherit expression-source, feature-axis, identifier, and CollecTRI coverage limitations from the activity artifact.",
            "Technical batch, repeated measures, and donor effects are not modeled by rankby_group.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return table, summary


def rank_tf_activities_code(
    *,
    function_name: str = "rank_tf_activity_artifact",
    parameters: Mapping[str, Any],
) -> str:
    """Return importable equivalent ranking source that never reruns activity inference."""

    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated TF ranking function_name must be a Python identifier.")
    portable_parameters = dict(parameters)
    portable_parameters["_portable_artifact"] = True
    helper_sources = "\n\n".join(
        dedent(inspect.getsource(helper)).strip()
        for helper in (
            _canonical_json_sha256,
            _canonical_axis,
            _activity_frame_fingerprint,
            _portable_tf_activity_payload,
            validate_tf_activity_artifact,
            _signature_shape,
            _bh_adjust,
            _package_version,
            collect_software_versions,
            _validate_decoupler_ranking,
            _canonical_annotation,
            _parse_reference,
            _independent_comparison,
            rank_tf_activities,
        )
    )
    argument_lines = [
        f"        {name}={value!r}," for name, value in portable_parameters.items()
    ]
    argument_lines.append("        decoupler_module=decoupler_module,")
    arguments = "\n".join(argument_lines)
    return f'''from __future__ import annotations

import copy
import hashlib
import importlib
import inspect
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
DECOUPLER_VERSION = {DECOUPLER_VERSION!r}
TF_ACTIVITY_ARTIFACT_TYPE = {TF_ACTIVITY_ARTIFACT_TYPE!r}
TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION = {TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION!r}
TF_ACTIVITY_PRODUCER_NODE_ID = {TF_ACTIVITY_PRODUCER_NODE_ID!r}
TF_ACTIVITY_PRODUCER_SCHEMA = {TF_ACTIVITY_PRODUCER_SCHEMA!r}
TF_RANKING_METHODS = {TF_RANKING_METHODS!r}
TF_RANKING_SUMMARY_SCHEMA = {TF_RANKING_SUMMARY_SCHEMA!r}
TF_RANKING_COLUMNS = {TF_RANKING_COLUMNS!r}
TF_RANKING_REFERENCES = {TF_RANKING_REFERENCES!r}

{helper_sources}


def {function_name}(adata, activities, decoupler_module=None):
    """Return the complete canonical ranking and strict summary without rerunning ULM."""
    return rank_tf_activities(
        adata,
        activities,
{arguments}
    )
'''


__all__ = [
    "TF_RANKING_COLUMNS",
    "TF_RANKING_METHODS",
    "TF_RANKING_SUMMARY_SCHEMA",
    "rank_tf_activities",
    "rank_tf_activities_code",
]
