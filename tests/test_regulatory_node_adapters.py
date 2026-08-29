from __future__ import annotations

from openbio_singlecell.nodes_regulatory import (
    OpenBioSingleCellCollecTRIULM,
    OpenBioSingleCellSCENICActivityBinarization,
)


def test_collectri_node_flattens_only_the_selected_resource_for_worker_transport():
    local = OpenBioSingleCellCollecTRIULM.prepare_worker_arguments(
        {
            "adata": "ticket",
            "resource": {
                "resource": "local_network",
                "network_csv": "collectri.csv",
                "resource_metadata_json": '{"license":"fixture"}',
            },
        }
    )
    official = OpenBioSingleCellCollecTRIULM.prepare_worker_arguments(
        {
            "adata": "ticket",
            "resource": {
                "resource": "official_collectri",
                "affiliation_license": "nonprofit",
                "allow_network_access": True,
            },
        }
    )

    assert local == {
        "adata": "ticket",
        "resource_mode": "local_network",
        "network_csv": "collectri.csv",
        "resource_metadata_json": '{"license":"fixture"}',
        "affiliation_license": "academic",
        "allow_network_access": False,
    }
    assert official == {
        "adata": "ticket",
        "resource_mode": "official_collectri",
        "resource_metadata_json": "{}",
        "affiliation_license": "nonprofit",
        "allow_network_access": True,
    }


def test_scenic_binarization_omits_only_an_unconnected_threshold_override():
    assert OpenBioSingleCellSCENICActivityBinarization.prepare_worker_arguments(
        {"scenic_result": "ticket", "threshold_overrides": None, "random_seed": 7}
    ) == {"scenic_result": "ticket", "random_seed": 7}

    threshold_ticket = object()
    assert (
        OpenBioSingleCellSCENICActivityBinarization.prepare_worker_arguments(
            {
                "scenic_result": "ticket",
                "threshold_overrides": threshold_ticket,
                "random_seed": 7,
            }
        )["threshold_overrides"]
        is threshold_ticket
    )
