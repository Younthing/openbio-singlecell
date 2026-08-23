from __future__ import annotations

import math

from openbio_singlecell import dependencies
from openbio_singlecell.contracts import METADATA_KEY, ensure_metadata, record_history


def test_metadata_uses_anndata_display_name_and_normalizes_values():
    adata = dependencies.ad.AnnData(dependencies.np.eye(2))

    metadata = ensure_metadata(
        adata,
        display_name="sample",
        source={"path": "sample.h5ad", "ignored": None},
    )

    assert metadata["display_name"] == "sample"
    assert metadata["source"] == {"path": "sample.h5ad"}
    assert metadata["analysis_history"] == {}
    assert adata.uns[METADATA_KEY] == metadata


def test_history_is_append_only_and_json_safe():
    adata = dependencies.ad.AnnData(dependencies.np.eye(1))
    adata.uns[METADATA_KEY] = {
        "analysis_history": {
            "000001": {"operation": "first"},
            "000004": {"operation": "fourth"},
        }
    }

    entry = record_history(
        adata,
        "normalize_total",
        {"target_sum": dependencies.np.float64(10_000), "invalid": math.nan},
        1,
        1,
    )

    metadata = ensure_metadata(adata)
    assert metadata["analysis_history"]["000005"] == entry
    assert entry["parameters"] == {"target_sum": 10_000.0, "invalid": "nan"}
    assert metadata["display_name"] == "AnnData"


def test_sequence_history_is_normalized():
    adata = dependencies.ad.AnnData(dependencies.np.eye(1))
    adata.uns[METADATA_KEY] = {"analysis_history": [{"operation": "legacy"}]}

    metadata = ensure_metadata(adata)

    assert metadata["analysis_history"] == {"000000": {"operation": "legacy"}}
    assert metadata["display_name"] == "AnnData"
