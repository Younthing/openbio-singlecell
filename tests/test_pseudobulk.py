from __future__ import annotations

from types import SimpleNamespace

from openbio_singlecell import nodes_differential
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellPseudobulk,
    OpenBioSingleCellPseudobulkDESeq2,
)


def output_value(node_output):
    return node_output.result[0]


def test_pseudobulk_filters_aggregates_with_node_thresholds(science, monkeypatch):
    obs = science.pd.DataFrame(
        {
            "sample": ["s1", "s1", "s1", "s2", "s2", "s2", "s2", "s2"],
            "cell_type": ["A", "A", "B", "A", "A", "B", "B", "B"],
        },
        index=[f"cell_{index}" for index in range(8)],
    )
    var = science.pd.DataFrame(index=["G1", "G2"])
    counts = science.np.asarray(
        [
            [1, 1],
            [2, 2],
            [10, 0],
            [1, 0],
            [1, 0],
            [2, 1],
            [2, 1],
            [2, 1],
        ],
        dtype=float,
    )
    adata = science.ad.AnnData(counts.copy(), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    compute_calls = []

    class FakePseudobulkSpace:
        def compute(
            self,
            adata,
            *,
            target_col,
            groups_col,
            layer_key=None,
            embedding_key=None,
            mode="sum",
        ):
            compute_calls.append(
                {
                    "target_col": target_col,
                    "groups_col": groups_col,
                    "layer_key": layer_key,
                    "embedding_key": embedding_key,
                    "mode": mode,
                }
            )
            matrix = adata.layers[layer_key] if layer_key is not None else adata.X
            rows = []
            aggregate_obs = []
            aggregate_names = []
            grouped = adata.obs.groupby([target_col, groups_col], sort=False, observed=True)
            for (sample, population), labels in grouped.groups.items():
                positions = adata.obs.index.get_indexer(labels)
                values = matrix[positions]
                if science.sparse.issparse(values):
                    values = values.toarray()
                values = science.np.asarray(values)
                rows.append(values.sum(axis=0) if mode == "sum" else values.mean(axis=0))
                aggregate_obs.append({target_col: sample, groups_col: population})
                aggregate_names.append(f"{sample}|{population}")
            return science.ad.AnnData(
                science.np.asarray(rows),
                obs=science.pd.DataFrame(aggregate_obs, index=aggregate_names),
                var=adata.var.copy(),
            )

    fake_pertpy = SimpleNamespace(tl=SimpleNamespace(PseudobulkSpace=FakePseudobulkSpace))
    monkeypatch.setattr(nodes_differential, "_require_pertpy", lambda: fake_pertpy)

    output = output_value(
        OpenBioSingleCellPseudobulk.execute(
            adata,
            sample_key="sample",
            groupby="cell_type",
            source={"source": "layer", "layer_name": "counts"},
            mode="sum",
            min_cells=2,
            min_counts=5,
        )
    )

    assert compute_calls == [
        {
            "target_col": "sample",
            "groups_col": "cell_type",
            "layer_key": "counts",
            "embedding_key": None,
            "mode": "sum",
        }
    ]
    assert output.obs[["sample", "cell_type"]].to_dict("records") == [
        {"sample": "s1", "cell_type": "A"},
        {"sample": "s2", "cell_type": "B"},
    ]
    assert output.obs["n_cells"].tolist() == [2, 3]
    assert output.obs["total_counts"].tolist() == [6.0, 9.0]
    science.np.testing.assert_array_equal(output.X, science.np.asarray([[3, 3], [6, 3]]))
    science.np.testing.assert_array_equal(output.layers["counts"], output.X)


def test_pseudobulk_deseq2_uses_pertpy_contrast_object(science, monkeypatch):
    obs = science.pd.DataFrame(
        {"group": ["Normal", "Normal", "nonDM_ED", "nonDM_ED"]},
        index=["s1", "s2", "s3", "s4"],
    )
    adata = science.ad.AnnData(
        science.np.asarray([[10, 2], [12, 3], [21, 8], [19, 7]], dtype=int),
        obs=obs,
        var=science.pd.DataFrame(index=["G1", "G2"]),
    )
    contrast = object()
    calls = []

    class FakePyDESeq2:
        def __init__(self, *, adata, design):
            calls.append(("init", adata, design))

        def fit(self):
            calls.append(("fit",))

        def contrast(self, *, column, baseline, group_to_compare):
            calls.append(("contrast", column, baseline, group_to_compare))
            return contrast

        def test_contrasts(self, selected_contrast):
            assert selected_contrast is contrast, "test_contrasts requires the object returned by contrast()"
            calls.append(("test_contrasts", selected_contrast))
            return science.pd.DataFrame(
                {"variable": ["G1"], "log2FoldChange": [1.25], "padj": [0.01]},
            )

    fake_pertpy = SimpleNamespace(tl=SimpleNamespace(PyDESeq2=FakePyDESeq2))
    monkeypatch.setattr(nodes_differential, "_require_pertpy", lambda: fake_pertpy)

    result = output_value(
        OpenBioSingleCellPseudobulkDESeq2.execute(
            adata,
            design="~group",
            contrast_column="group",
            baseline="Normal",
            comparison="nonDM_ED",
        )
    )

    assert calls == [
        ("init", adata, "~group"),
        ("fit",),
        ("contrast", "group", "Normal", "nonDM_ED"),
        ("test_contrasts", contrast),
    ]
    assert result.title == "DESeq2: nonDM_ED vs Normal"
    assert result.table.columns.tolist().count("gene") == 1
    assert result.table.to_dict("records") == [
        {"gene": "G1", "log2FoldChange": 1.25, "padj": 0.01}
    ]
