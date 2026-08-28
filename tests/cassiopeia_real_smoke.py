"""Standalone real-backend smoke; run in the pinned Linux lineage environment."""

from __future__ import annotations

import json
import tempfile
from importlib import metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from openbio_singlecell.cassiopeia_tree import (
    CASSIOPEIA_DISTRIBUTION,
    CASSIOPEIA_VERSION,
    add_cassiopeia_plasticity,
    compute_cassiopeia_expansions,
    prepare_cassiopeia_characters,
    reconstruct_cassiopeia_tree,
)


def main() -> None:
    assert metadata.version(CASSIOPEIA_DISTRIBUTION) == CASSIOPEIA_VERSION
    rows = pd.DataFrame(
        [
            {"Tumor": "T1", "cellBC": "c1", "intBC": "i1", "r1": "A", "r2": "NONE"},
            {"Tumor": "T1", "cellBC": "c2", "intBC": "i1", "r1": "A", "r2": "B"},
            {"Tumor": "T1", "cellBC": "c3", "intBC": "i1", "r1": "C", "r2": "B"},
            {"Tumor": "T1", "cellBC": "c4", "intBC": "i1", "r1": "C", "r2": "NONE"},
        ]
    )
    with tempfile.TemporaryDirectory(prefix="openbio-cassiopeia-smoke-") as directory:
        path = Path(directory) / "alleles.tsv"
        rows.to_csv(path, sep="\t", index=False)
        characters, qc, qc_summary, qc_code = prepare_cassiopeia_characters(
            str(path),
            first_column_as_index=False,
            cut_site_columns="r1,r2",
            prior_grouping_columns="Tumor,intBC",
            minimum_unique_fraction=0.0,
            minimum_informative_character_fraction=0.0,
            maximum_missing_fraction=1.0,
            maximum_uncut_fraction=1.0,
            openbio_version="real-smoke",
        )
        assert qc.loc[0, "status"] == "pass"
        first_tree, tree_summary, tree_code = reconstruct_cassiopeia_tree(
            characters, "T1", openbio_version="real-smoke"
        )
        second_tree, _, _ = reconstruct_cassiopeia_tree(
            characters, "T1", openbio_version="real-smoke"
        )
        assert first_tree.topology_fingerprint == second_tree.topology_fingerprint
        expansion, expansion_summary, expansion_code = compute_cassiopeia_expansions(
            first_tree,
            minimum_clade_size=2,
            minimum_depth=1,
            openbio_version="real-smoke",
        )
        assert expansion.shape[0] == tree_summary["key_results"]["node_count"]
        adata = ad.AnnData(np.zeros((4, 1), dtype=float))
        adata.obs_names = ["c1", "c2", "c3", "c4"]
        adata.var_names = ["g1"]
        adata.obs["cell_type"] = pd.Categorical(["A", "A", "B", "B"])
        output, cells, plasticity_summary, plasticity_code = add_cassiopeia_plasticity(
            adata,
            first_tree,
            annotation_key="cell_type",
            minimum_state_fraction=0.025,
            openbio_version="real-smoke",
        )
        assert cells.shape[0] == 4
        assert output.obs["sc_effective_plasticity"].notna().all()
        assert plasticity_summary["key_results"]["global_small_parsimony"] == 1
        for summary in (qc_summary, tree_summary, expansion_summary, plasticity_summary):
            json.dumps(summary, allow_nan=False)
            assert summary["software_versions"][CASSIOPEIA_DISTRIBUTION] == CASSIOPEIA_VERSION
        for source in (qc_code, tree_code, expansion_code, plasticity_code):
            compile(source, "<cassiopeia-real-code>", "exec")
    print("cassiopeia-mt 2.1.3 real smoke: PASS")


if __name__ == "__main__":
    main()
