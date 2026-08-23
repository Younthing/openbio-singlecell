from __future__ import annotations

import math

import pytest

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION
from openbio_singlecell.contracts import METADATA_KEY, ensure_metadata, record_history


def test_metadata_uses_anndata_display_name_and_normalizes_values(science):
    adata = science.ad.AnnData(science.np.eye(2))

    metadata = ensure_metadata(
        adata,
        display_name="sample",
        source={"path": "sample.h5ad", "ignored": None},
    )

    assert metadata["display_name"] == "sample"
    assert metadata["source"] == {"path": "sample.h5ad"}
    assert metadata["analysis_history"] == {}
    assert adata.uns[METADATA_KEY] == metadata


def test_history_is_append_only_and_json_safe(science):
    adata = science.ad.AnnData(science.np.eye(1))
    adata.uns[METADATA_KEY] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "analysis_history": {
            "000001": {"operation": "first"},
            "000004": {"operation": "fourth"},
        },
    }

    entry = record_history(
        adata,
        "normalize_total",
        {"target_sum": science.np.float64(10_000), "invalid": math.nan},
        1,
        1,
    )

    metadata = ensure_metadata(adata)
    assert metadata["analysis_history"]["000005"] == entry
    assert entry["parameters"] == {"target_sum": 10_000.0, "invalid": "nan"}
    assert metadata["display_name"] == "AnnData"


@pytest.mark.parametrize(
    "history_type",
    ["list", "tuple", "ndarray", "scalar"],
)
def test_analysis_history_rejects_non_mapping_values(history_type, science):
    values = [{"operation": "old"}]
    if history_type == "list":
        history = values
    elif history_type == "tuple":
        history = tuple(values)
    elif history_type == "ndarray":
        history = science.np.asarray(["old"])
    else:
        history = 42

    adata = science.ad.AnnData(science.np.eye(1))
    adata.uns[METADATA_KEY] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "analysis_history": history,
    }

    with pytest.raises(ValueError, match="analysis_history must be a mapping"):
        ensure_metadata(adata)


@pytest.mark.parametrize("schema_version", [None, 0, SCHEMA_VERSION + 1])
def test_metadata_rejects_missing_or_unsupported_schema_versions(schema_version, science):
    adata = science.ad.AnnData(science.np.eye(1))
    adata.uns[METADATA_KEY] = {
        "schema_version": schema_version,
        "version": PLUGIN_VERSION,
        "analysis_history": {},
    }

    with pytest.raises(ValueError, match="unsupported schema_version"):
        ensure_metadata(adata)


def test_metadata_rejects_non_mapping_payload(science):
    adata = science.ad.AnnData(science.np.eye(1))
    adata.uns[METADATA_KEY] = "invalid"

    with pytest.raises(ValueError, match="must be a mapping"):
        ensure_metadata(adata)
