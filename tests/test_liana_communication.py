from __future__ import annotations

import copy
import importlib.metadata
import json
import types
import uuid

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from matplotlib.figure import Figure

from openbio_singlecell.liana_artifact_codec import LIANA_CODEC, read_liana_result, write_liana_result
from openbio_singlecell.liana_communication import (
    liana_communication_code,
    liana_resource_cache_fingerprint,
    run_liana_communication,
)
from openbio_singlecell.liana_plot import liana_dot_plot_code, render_liana_dot_plot
from openbio_singlecell.liana_result import LianaResult, build_liana_result, validate_liana_result
from openbio_singlecell.node_types import (
    LianaResultType,
    PlotResultType,
    SummaryResultType,
)
from openbio_singlecell.nodes_communication import (
    OpenBioSingleCellLianaCommunication,
    OpenBioSingleCellLianaDotPlot,
)


def _resource_metadata(name: str = "consensus", organism: str = "Homo sapiens") -> str:
    return json.dumps(
        {
            "name": name,
            "version": "test-snapshot-1",
            "release_date": "2026-08-20",
            "download_url": "https://example.org/licensed-liana-resource.csv",
            "organism": organism,
            "gene_identifier_namespace": "HGNC symbol",
            "scope": "two-pair test snapshot",
            "license": "CC0 test fixture",
            "source_license_review": "reviewed for this test",
            "citation": "Test ligand-receptor resource",
        },
        separators=(",", ":"),
    )


def _adata() -> AnnData:
    observations = pd.DataFrame(
        {
            "sample": ["S1"] * 4 + ["S2"] * 4,
            "condition": ["control"] * 4 + ["treated"] * 4,
            "cell_type": ["A", "A", "B", "B"] * 2,
        },
        index=[f"cell_{index}" for index in range(8)],
    )
    values = np.log1p(np.arange(32, dtype=float).reshape(8, 4) / 7.0 + 1.0)
    adata = AnnData(
        X=values,
        obs=observations,
        var=pd.DataFrame(index=["L1", "L2", "R1", "R2"]),
    )
    adata.raw = adata.copy()
    adata.uns["openbio_singlecell"] = {
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {"source": "X", "layer_name": "counts"},
            },
            "000001": {"operation": "log1p", "parameters": {}},
        }
    }
    return adata


class _FakeResourceSelector:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def select_resource(self, resource_name="consensus"):
        self.calls.append(resource_name)
        return pd.DataFrame({"ligand": ["L1", "L2"], "receptor": ["R1", "R2"]})


class _FakeRankAggregate:
    def __init__(self, *, original=None, attack=None, transform=None) -> None:
        self.original = original
        self.attack = attack
        self.transform = transform
        self.calls = []

    def __call__(
        self,
        adata,
        groupby,
        resource_name="consensus",
        expr_prop=0.1,
        min_cells=5,
        groupby_pairs=None,
        base=np.e,
        aggregate_method="rra",
        consensus_opts=None,
        return_all_lrs=False,
        key_added="liana_res",
        use_raw=True,
        layer=None,
        de_method="t-test",
        n_perms=1000,
        seed=1337,
        n_jobs=1,
        resource=None,
        interactions=None,
        mdata_kwargs=None,
        spatial_key=None,
        spatial_kwargs=None,
        inplace=True,
        verbose=False,
    ):
        raise AssertionError("The adapter must call by_sample, not the pooled method.")

    def by_sample(
        self,
        adata,
        sample_key,
        key_added="liana_res",
        inplace=True,
        verbose=False,
        **kwargs,
    ):
        self.calls.append((adata, sample_key, key_added, inplace, verbose, kwargs))
        adata.obs[sample_key] = adata.obs[sample_key].astype("category")
        adata.uns[key_added] = {"expected_private_scratch": True}
        if self.attack == "caller":
            self.original.X[0, 0] += 1.0
        if self.attack == "resource":
            kwargs["resource"].iloc[0, 0] = "MUTATED"
        rows = []
        for sample in ("S1", "S2"):
            rows.extend(
                [
                    [sample, "A", "B", "L1", "R1", 2.0, 0.04, 2.0, 0.3, 0.5, 0.6, 0.7, 0.04, 0.10],
                    [sample, "A", "B", "L2", "R2", 1.0, 0.20, 1.0, 0.1, -0.2, 0.3, 0.4, 0.20, 0.30],
                ]
            )
        table = pd.DataFrame(
            rows,
            columns=[
                "sample",
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "lr_means",
                "cellphone_pvals",
                "expr_prod",
                "scaled_weight",
                "lr_logfc",
                "spec_weight",
                "lrscore",
                "specificity_rank",
                "magnitude_rank",
            ],
        )
        return table if self.transform is None else self.transform(table)


class _FakeCellPhoneDB:
    def __init__(self, *, transform=None) -> None:
        self.transform = transform
        self.calls = []

    def __call__(
        self,
        adata,
        groupby,
        resource_name="consensus",
        expr_prop=0.1,
        min_cells=5,
        groupby_pairs=None,
        base=np.e,
        supp_columns=None,
        return_all_lrs=False,
        key_added="liana_res",
        use_raw=True,
        layer=None,
        de_method="t-test",
        n_perms=1000,
        seed=1337,
        n_jobs=1,
        resource=None,
        interactions=None,
        spatial_key="spatial",
        spatial_kwargs=None,
        mdata_kwargs=None,
        inplace=True,
        verbose=False,
    ):
        raise AssertionError("The adapter must call by_sample, not the pooled method.")

    def by_sample(
        self,
        adata,
        sample_key,
        key_added="liana_res",
        inplace=True,
        verbose=False,
        **kwargs,
    ):
        self.calls.append((adata, sample_key, key_added, inplace, verbose, kwargs))
        adata.obs[sample_key] = adata.obs[sample_key].astype("category")
        adata.uns[key_added] = {}
        rows = []
        for sample in ("S1", "S2"):
            rows.extend(
                [
                    [sample, "L1", "L1", 2.0, 0.8, "R1", "R1", 3.0, 0.9, "A", "B", 2.5, 0.04],
                    [sample, "L2", "L2", 1.0, 0.7, "R2", "R2", 1.5, 0.8, "A", "B", 1.25, 0.20],
                ]
            )
        table = pd.DataFrame(
            rows,
            columns=[
                "sample",
                "ligand",
                "ligand_complex",
                "ligand_means",
                "ligand_props",
                "receptor",
                "receptor_complex",
                "receptor_means",
                "receptor_props",
                "source",
                "target",
                "lr_means",
                "cellphone_pvals",
            ],
        )
        return table if self.transform is None else self.transform(table)


def _fake_liana(method: str, method_object=None, *, version="1.9.0"):
    method_object = method_object or (_FakeRankAggregate() if method == "rank_aggregate" else _FakeCellPhoneDB())
    selector = _FakeResourceSelector()
    return (
        types.SimpleNamespace(
            __version__=version,
            mt=types.SimpleNamespace(**{method: method_object}),
            rs=selector,
        ),
        method_object,
        selector,
    )


def _parameters(method: str, fake, **overrides):
    values = {
        "sample_key": "sample",
        "condition_key": "condition",
        "identity_key": "cell_type",
        "annotation_status": "curated",
        "organism": "Homo sapiens",
        "method": method,
        "resource_mode": "bundled_human",
        "resource_name": "consensus",
        "resource_path": None,
        "resource_metadata_json": _resource_metadata(),
        "source_kind": "X",
        "layer_name": None,
        "expression_proportion": 0.1,
        "min_cells_per_identity_sample": 2,
        "permutations": 1000,
        "random_seed": 7,
        "jobs": 1,
        "max_output_rows": 100,
        "max_working_memory_gib": 1.0,
        "openbio_version": "0.2.0",
        "liana_module": fake,
    }
    values.update(overrides)
    return values


def _run(method="rank_aggregate", *, adata=None, method_object=None, **overrides):
    adata = _adata() if adata is None else adata
    fake, backend, selector = _fake_liana(method, method_object)
    artifact, summary = run_liana_communication(
        adata,
        **_parameters(method, fake, **overrides),
    )
    return adata, fake, backend, selector, artifact, summary


@pytest.mark.parametrize("method", ["rank_aggregate", "cellphonedb"])
def test_communication_exact_by_sample_complete_family_summary_and_immutability(method):
    adata = _adata()
    original = adata.copy()
    adata, _fake, backend, selector, artifact, summary = _run(method, adata=adata)
    table, provenance, metadata = validate_liana_result(artifact)

    assert isinstance(artifact, LianaResult)
    assert list(table.columns) == list(
        {
            "rank_aggregate": (
                "sample",
                "condition",
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "lr_means",
                "cellphone_pvals",
                "expr_prod",
                "scaled_weight",
                "lr_logfc",
                "spec_weight",
                "lrscore",
                "specificity_rank",
                "magnitude_rank",
            ),
            "cellphonedb": (
                "sample",
                "condition",
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "ligand",
                "receptor",
                "ligand_means",
                "ligand_props",
                "receptor_means",
                "receptor_props",
                "lr_means",
                "cellphone_pvals",
            ),
        }[method]
    )
    assert len(table) == 4
    assert table.groupby("sample", sort=False).size().to_dict() == {"S1": 2, "S2": 2}
    assert table[["sample", "condition"]].drop_duplicates().to_dict("records") == [
        {"sample": "S1", "condition": "control"},
        {"sample": "S2", "condition": "treated"},
    ]
    assert summary["key_results"]["result"]["complete_backend_family_returned"] is True
    assert summary["key_results"]["resource"]["accounting"]["retained_rows_on_exact_feature_axis"] == 2
    assert summary["software_versions"]["liana"] == "1.9.0"
    assert "not tested" in summary["results"]
    assert all("adjusted" not in text.lower() for text in [summary["methods"], summary["results"]])
    json.dumps(summary, allow_nan=False)
    assert len(summary["references"]) >= 4
    assert provenance["resource"]["metadata"]["license"] == "CC0 test fixture"
    assert len(metadata["table_identity"]["content_sha256"]) == 64
    assert selector.calls == ["consensus"]
    assert len(backend.calls) == 1
    work, sample_key, key_added, inplace, verbose, kwargs = backend.calls[0]
    assert work is adata
    assert sample_key == "sample"
    assert key_added == "__openbio_liana_private_scratch__"
    assert inplace is False and verbose is False
    assert kwargs["resource_name"] == "consensus"
    assert kwargs["return_all_lrs"] is False
    assert kwargs["use_raw"] is False
    assert kwargs["layer"] is None
    assert kwargs["n_perms"] == 1000
    assert kwargs["seed"] == 7
    assert kwargs["n_jobs"] == 1
    assert kwargs["spatial_key"] is None
    if method == "rank_aggregate":
        assert kwargs["aggregate_method"] == "rra"
        assert kwargs["consensus_opts"] is None
    else:
        assert kwargs["supp_columns"] is None
    assert np.array_equal(adata.X, original.X)
    assert np.array_equal(adata.obs.to_numpy(dtype=str), original.obs.to_numpy(dtype=str))
    assert np.array_equal(adata.raw.X, original.raw.X)


def test_liana_result_is_defensive_exact_typed_and_detects_private_tampering():
    *_prefix, artifact, _summary = _run()
    exported = artifact.table
    exported.iloc[0, exported.columns.get_loc("magnitude_rank")] = 0.99
    assert artifact.table.loc[0, "magnitude_rank"] == 0.10
    portable = artifact.portable()
    portable["provenance"]["organism"] = "tampered"
    assert artifact.provenance["organism"] == "Homo sapiens"
    with pytest.raises(TypeError, match="exact OPENBIO_LIANA_RESULT"):
        validate_liana_result(artifact.portable())
    validate_liana_result(artifact.portable(), exact_type=False)
    private_table = object.__getattribute__(artifact, "_table")
    private_table.loc[0, "magnitude_rank"] = 0.99
    with pytest.raises(ValueError, match="current-content fingerprint"):
        validate_liana_result(artifact)


def test_liana_artifact_round_trip_preserves_validated_result_without_copying_owned_table(tmp_path):
    *_prefix, artifact, _summary = _run()
    owned_table = object.__getattribute__(artifact, "_table")

    assert validate_liana_result(artifact, copy_result=False)[0] is owned_table

    artifact_root = tmp_path / "liana"
    artifact_root.mkdir()
    write_liana_result(artifact_root, artifact)
    restored = read_liana_result(artifact_root)
    table, provenance, metadata = validate_liana_result(
        restored,
        exact_type=False,
        copy_result=False,
    )

    assert table.equals(owned_table)
    assert provenance == artifact.provenance
    assert metadata == artifact.metadata


def test_liana_result_builder_can_adopt_worker_owned_table_without_copying():
    *_prefix, artifact, _summary = _run()
    owned_table = artifact.table

    adopted = build_liana_result(
        table=owned_table,
        method="rank_aggregate",
        provenance=artifact.provenance,
        numpy=np,
        pandas=pd,
        copy_table=False,
    )

    assert object.__getattribute__(adopted, "_table") is owned_table


def test_worker_owned_liana_uses_the_worker_private_anndata_in_place():
    from openbio_singlecell.operations_communication import liana_communication_owned

    adata = _adata()
    fake, backend, _selector = _fake_liana("rank_aggregate")

    result, report, code = liana_communication_owned(
        adata,
        sample_key="sample",
        condition_key="condition",
        identity_key="cell_type",
        annotation_status="curated",
        organism="Homo sapiens",
        method="rank_aggregate",
        resource_mode="bundled_human",
        resource_name="consensus",
        resource_path=None,
        resource_metadata_json=_resource_metadata(),
        source={"source": "X"},
        expression_proportion=0.1,
        min_cells_per_identity_sample=2,
        permutations=1000,
        random_seed=7,
        jobs=1,
        max_output_rows=100,
        max_working_memory_gib=1.0,
        liana_module=fake,
    )

    assert backend.calls[0][0] is adata
    validate_liana_result(result, copy_result=False)
    assert report.summary["node_id"] == "OpenBioSingleCellLianaCommunication"
    assert "method_object.by_sample" in code


def test_worker_owned_liana_copies_resource_only_for_backend_workspace(monkeypatch):
    import openbio_singlecell.liana_communication as communication

    adata = _adata()
    fake, _backend, _selector = _fake_liana("rank_aggregate")
    resolved_resource_id = None
    copy_count = 0
    original_resolve = communication._liana_resolve_resource
    original_copy = pd.DataFrame.copy

    def tracked_resolve(*args, **kwargs):
        nonlocal resolved_resource_id
        result = original_resolve(*args, **kwargs)
        resolved_resource_id = id(result[0])
        return result

    def tracked_copy(frame, *args, **kwargs):
        nonlocal copy_count
        if id(frame) == resolved_resource_id:
            copy_count += 1
        return original_copy(frame, *args, **kwargs)

    monkeypatch.setattr(communication, "_liana_resolve_resource", tracked_resolve)
    monkeypatch.setattr(pd.DataFrame, "copy", tracked_copy)
    communication._liana_run_impl(adata, **_parameters("rank_aggregate", fake))

    assert copy_count == 1


@pytest.mark.parametrize(
    ("mutation", "match"),
    [("caller", "worker-owned selected expression"), ("resource", "private pinned resource")],
)
def test_adversarial_backend_mutation_is_detected(mutation, match):
    adata = _adata()
    backend = _FakeRankAggregate(original=adata, attack=mutation)
    fake, _backend, _selector = _fake_liana("rank_aggregate", backend)
    with pytest.raises(RuntimeError, match=match):
        run_liana_communication(adata, **_parameters("rank_aggregate", fake))


@pytest.mark.parametrize(
    ("transform", "match"),
    [
        (lambda table: table.assign(unexpected=1.0), "complete column family"),
        (
            lambda table: pd.concat([table, table.iloc[[0]]], ignore_index=True),
            "duplicate rows",
        ),
        (lambda table: table.assign(magnitude_rank=np.nan), "non-finite"),
        (lambda table: table.assign(sample="S1"), "Sample family"),
        (lambda table: table.assign(ligand_complex="OUTSIDE"), "outside the pinned"),
    ],
)
def test_malformed_backend_outputs_fail_closed(transform, match):
    backend = _FakeRankAggregate(transform=transform)
    fake, _backend, _selector = _fake_liana("rank_aggregate", backend)
    with pytest.raises(RuntimeError, match=match):
        run_liana_communication(_adata(), **_parameters("rank_aggregate", fake))


def test_exact_version_signature_and_parameter_guards_fail_before_backend():
    bad_version, backend, _selector = _fake_liana("rank_aggregate", version="1.10.0")
    with pytest.raises(RuntimeError, match="exact liana==1.9.0"):
        run_liana_communication(_adata(), **_parameters("rank_aggregate", bad_version))
    assert backend.calls == []

    class Drifted(_FakeRankAggregate):
        def by_sample(self, adata, sample_key, **kwargs):
            return None

    fake, drifted, _selector = _fake_liana("rank_aggregate", Drifted())
    with pytest.raises(RuntimeError, match="signature drifted"):
        run_liana_communication(_adata(), **_parameters("rank_aggregate", fake))
    assert drifted.calls == []

    fake, backend, _selector = _fake_liana("rank_aggregate")
    with pytest.raises(ValueError, match="permutations"):
        run_liana_communication(_adata(), **_parameters("rank_aggregate", fake, permutations=999))
    assert backend.calls == []


def test_sample_condition_and_memory_guards():
    fake, backend, _selector = _fake_liana("rank_aggregate")
    conflict = _adata()
    conflict.obs.loc["cell_1", "condition"] = "treated"
    with pytest.raises(ValueError, match="more than one Condition"):
        run_liana_communication(conflict, **_parameters("rank_aggregate", fake))
    assert backend.calls == []

    large_features = ["L1", "L2", "R1", "R2"] + [f"G{i}" for i in range(50_000)]
    large = AnnData(
        X=np.full((8, len(large_features)), 0.25),
        obs=_adata().obs.copy(),
        var=pd.DataFrame(index=large_features),
    )
    large.raw = large.copy()
    large.uns["openbio_singlecell"] = copy.deepcopy(_adata().uns["openbio_singlecell"])
    with pytest.raises(MemoryError, match="working memory"):
        run_liana_communication(
            large,
            **_parameters("rank_aggregate", fake, max_working_memory_gib=0.001),
        )


def test_generated_communication_code_has_full_runtime_parity_and_pinned_resource():
    adata, fake, _backend, _selector, artifact, summary = _run()
    accounting = summary["key_results"]["resource"]["accounting"]
    parameters = _parameters("rank_aggregate", fake)
    parameters.pop("liana_module")
    parameters.update(
        {
            "expected_file_sha256": accounting["raw_file_sha256"],
            "expected_canonical_resource_sha256": accounting["canonical_resource_sha256"],
            "expected_effective_resource_sha256": accounting["effective_resource_sha256"],
        }
    )
    code = liana_communication_code(parameters=parameters)
    namespace = {}
    exec(code, namespace)
    generated_table, generated_summary = namespace["run_liana_communication"](adata, liana_module=fake)
    pd.testing.assert_frame_equal(generated_table, artifact.table, check_exact=True)
    assert generated_summary == summary
    assert "from openbio_singlecell" not in code
    assert "import openbio_singlecell" not in code
    assert "import liana" not in code


def test_local_resource_raw_and_effective_hashes_are_pinned(tmp_path):
    resource_path = tmp_path / "resource.csv"
    resource_path.write_text("ligand,receptor\nL1,R1\nL1,R1\nL2,R2\n", encoding="utf-8")
    fake, backend, _selector = _fake_liana("rank_aggregate")
    artifact, summary = run_liana_communication(
        _adata(),
        **_parameters(
            "rank_aggregate",
            fake,
            resource_mode="local_resource",
            resource_name="local",
            resource_path=str(resource_path),
            resource_metadata_json=_resource_metadata("local"),
        ),
    )
    accounting = summary["key_results"]["resource"]["accounting"]
    cache_fingerprint = liana_resource_cache_fingerprint(str(resource_path), _resource_metadata("local"))
    assert cache_fingerprint[4] == accounting["raw_file_sha256"]
    with pytest.raises(ValueError, match="regular file"):
        liana_resource_cache_fingerprint(str(tmp_path), _resource_metadata("local"))
    assert accounting["duplicate_rows_collapsed"] == 1
    assert len(accounting["raw_file_sha256"]) == 64
    assert len(accounting["canonical_resource_sha256"]) == 64
    assert len(accounting["effective_resource_sha256"]) == 64
    validate_liana_result(artifact)
    with pytest.raises(ValueError, match="SHA-256 changed"):
        run_liana_communication(
            _adata(),
            **_parameters(
                "rank_aggregate",
                fake,
                resource_mode="local_resource",
                resource_name="local",
                resource_path=str(resource_path),
                resource_metadata_json=_resource_metadata("local"),
                expected_file_sha256="0" * 64,
            ),
        )
    assert len(backend.calls) == 1


class _DrawResult:
    def __init__(self) -> None:
        self.draw_calls = 0

    def draw(self):
        self.draw_calls += 1
        figure = Figure(figsize=(3, 2))
        axis = figure.subplots()
        axis.scatter([0.0, 1.0], [1.0, 0.0])
        return figure


class _FakePlotting:
    def __init__(self, *, mutate=False, wrong_return=False) -> None:
        self.calls = []
        self.mutate = mutate
        self.wrong_return = wrong_return

    def dotplot_by_sample(
        self,
        adata=None,
        uns_key="liana_res",
        liana_res=None,
        sample_key="sample",
        colour=None,
        size=None,
        inverse_colour=False,
        inverse_size=False,
        source_labels=None,
        target_labels=None,
        ligand_complex=None,
        receptor_complex=None,
        size_range=(2, 9),
        cmap="viridis",
        figure_size=(8, 6),
        return_fig=True,
    ):
        self.calls.append(
            {
                "adata": adata,
                "uns_key": uns_key,
                "liana_res": liana_res,
                "sample_key": sample_key,
                "colour": colour,
                "size": size,
                "inverse_colour": inverse_colour,
                "inverse_size": inverse_size,
                "source_labels": source_labels,
                "target_labels": target_labels,
                "ligand_complex": ligand_complex,
                "receptor_complex": receptor_complex,
                "size_range": size_range,
                "cmap": cmap,
                "figure_size": figure_size,
                "return_fig": return_fig,
            }
        )
        if self.mutate:
            liana_res.loc[0, colour] = 0.99
        return object() if self.wrong_return else _DrawResult()


def _plot_fake(*, mutate=False, wrong_return=False, version="1.9.0"):
    plotting = _FakePlotting(mutate=mutate, wrong_return=wrong_return)
    return types.SimpleNamespace(__version__=version, pl=plotting), plotting


def _plot_parameters(fake, **overrides):
    values = {
        "source_labels": ["A"],
        "target_labels": ["B"],
        "selection_threshold": 0.25,
        "top_n": 1,
        "figure_width": 5.0,
        "figure_height": 4.0,
        "max_plot_rows": 100,
        "max_image_pixels": 5_000_000,
        "openbio_version": "0.2.0",
        "liana_module": fake,
    }
    values.update(overrides)
    return values


def test_worker_owned_liana_plot_consumes_portable_artifact(tmp_path):
    from openbio_singlecell.operations_communication import liana_dot_plot_owned

    *_prefix, artifact, _summary = _run()
    root = tmp_path / "result"
    root.mkdir()
    write_liana_result(root, artifact)
    portable = read_liana_result(root)
    fake, plotting = _plot_fake()

    plotted, report, code = liana_dot_plot_owned(
        portable,
        source_labels="A",
        target_labels="B",
        selection_method="rank_aggregate",
        selection_threshold=0.25,
        top_n=1,
        figure_width=5.0,
        figure_height=4.0,
        max_plot_rows=100,
        max_image_pixels=5_000_000,
        liana_module=fake,
    )

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["node_id"] == "OpenBioSingleCellLianaDotPlot"
    assert len(plotting.calls) == 1
    assert "dotplot_by_sample" in code


def test_liana_worker_operations_publish_portable_artifacts_without_rewriting_input(
    tmp_path,
    monkeypatch,
):
    import openbio_singlecell.operations_communication as operations
    from openbio_singlecell.artifact_codecs import read_plot, write_anndata
    from openbio_singlecell.worker_protocol import OperationContext

    adata = _adata()
    fake, _backend, _selector = _fake_liana("rank_aggregate")
    result, report, code = operations.liana_communication_owned(
        adata,
        sample_key="sample",
        condition_key="condition",
        identity_key="cell_type",
        annotation_status="curated",
        organism="Homo sapiens",
        method="rank_aggregate",
        resource_mode="bundled_human",
        resource_name="consensus",
        resource_path=None,
        resource_metadata_json=_resource_metadata(),
        source={"source": "X"},
        expression_proportion=0.1,
        min_cells_per_identity_sample=2,
        permutations=1000,
        random_seed=7,
        jobs=1,
        max_output_rows=100,
        max_working_memory_gib=1.0,
        liana_module=fake,
    )
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _adata())
    input_bytes = (input_root / "data.h5ad").read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))
    monkeypatch.setattr(
        operations,
        "liana_communication_owned",
        lambda *_args, **_kwargs: (result, report, code),
    )

    records = operations.liana_communication(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            }
        },
        {
            "sample_key": "sample",
            "condition_key": "condition",
            "identity_key": "cell_type",
            "annotation_status": "curated",
            "organism": "Homo sapiens",
            "method": "rank_aggregate",
            "resource_mode": "bundled_human",
            "resource_name": "consensus",
            "resource_metadata_json": _resource_metadata(),
            "source": {"source": "X"},
            "expression_proportion": 0.1,
            "min_cells_per_identity_sample": 2,
            "permutations": 1000,
            "random_seed": 7,
            "jobs": 1,
            "max_output_rows": 100,
            "max_working_memory_gib": 1.0,
        },
    )
    portable = read_liana_result(staging / records[0]["payload"])
    fake_plot, _plotting = _plot_fake()
    plotted, plot_report, plot_code = operations.liana_dot_plot_owned(
        portable,
        source_labels="A",
        target_labels="B",
        selection_method="rank_aggregate",
        selection_threshold=0.25,
        top_n=1,
        figure_width=5.0,
        figure_height=4.0,
        max_plot_rows=100,
        max_image_pixels=5_000_000,
        liana_module=fake_plot,
    )
    monkeypatch.setattr(
        operations,
        "liana_dot_plot_owned",
        lambda *_args, **_kwargs: (plotted, plot_report, plot_code),
    )
    plot_records = operations.liana_dot_plot(
        context,
        {
            "result": {
                "type": "artifact",
                "path": str((staging / records[0]["payload"]).resolve()),
                "kind": "OPENBIO_LIANA_RESULT",
                "codec": LIANA_CODEC,
            }
        },
        {
            "source_labels": "A",
            "target_labels": "B",
            "selection_method": "rank_aggregate",
            "selection_threshold": 0.25,
            "top_n": 1,
            "figure_width": 5.0,
            "figure_height": 4.0,
            "max_plot_rows": 100,
            "max_image_pixels": 5_000_000,
        },
    )

    validate_liana_result(portable, exact_type=False, copy_result=False)
    assert [record["name"] for record in records] == ["result", "summary", "code"]
    assert [record["name"] for record in plot_records] == ["plot", "summary", "code"]
    assert read_plot(staging / plot_records[0]["payload"])[0].startswith(b"\x89PNG\r\n\x1a\n")
    assert (input_root / "data.h5ad").read_bytes() == input_bytes


@pytest.mark.parametrize(
    ("method", "colour", "size", "inverse_colour"),
    [
        ("rank_aggregate", "magnitude_rank", "specificity_rank", True),
        ("cellphonedb", "lr_means", "cellphone_pvals", False),
    ],
)
def test_dot_plot_consumes_only_typed_result_selects_stably_and_reports(method, colour, size, inverse_colour):
    *_prefix, artifact, communication_summary = _run(method)
    fake, plotting = _plot_fake()
    rng_before = copy.deepcopy(np.random.get_state())
    png, summary = render_liana_dot_plot(artifact, **_plot_parameters(fake))
    rng_after = np.random.get_state()
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert summary["key_results"]["selection"]["input_rows"] == 4
    assert summary["key_results"]["selection"]["plotted_rows"] == 2
    assert summary["key_results"]["selection"]["top_n_scope"].startswith("within each Sample")
    assert summary["key_results"]["producer_summary_sha256"] == artifact.provenance["summary_sha256"]
    assert communication_summary["key_results"]["result"]["rows"] == 4
    if method == "rank_aggregate":
        assert any(reference["doi"] == "10.1093/bioinformatics/btr709" for reference in summary["references"])
    json.dumps(summary, allow_nan=False)
    assert len(plotting.calls) == 1
    call = plotting.calls[0]
    assert call["adata"] is None
    assert isinstance(call["liana_res"], pd.DataFrame)
    assert len(call["liana_res"]) == 2
    assert call["sample_key"] == "sample"
    assert call["colour"] == colour
    assert call["size"] == size
    assert call["inverse_colour"] is inverse_colour
    assert call["inverse_size"] is True
    assert call["return_fig"] is True
    assert artifact.table.equals(validate_liana_result(artifact)[0])
    assert all(np.array_equal(left, right) for left, right in zip(rng_before, rng_after, strict=True))


def test_dot_plot_generated_code_matches_png_and_summary_exactly():
    *_prefix, artifact, _communication_summary = _run()
    fake, _plotting = _plot_fake()
    parameters = _plot_parameters(fake)
    parameters.pop("liana_module")
    runtime_png, runtime_summary = render_liana_dot_plot(artifact, **{**parameters, "liana_module": fake})
    code = liana_dot_plot_code(parameters=parameters)
    namespace = {}
    exec(code, namespace)
    generated_png, generated_summary = namespace["run_liana_dot_plot"](artifact.portable(), liana_module=fake)
    assert generated_png == runtime_png
    assert generated_summary == runtime_summary
    assert "run_liana_communication" not in code
    assert "from openbio_singlecell" not in code
    assert "import openbio_singlecell" not in code


def test_dot_plot_rejects_tamper_wrong_version_signature_return_and_backend_mutation():
    *_prefix, artifact, _summary = _run()
    fake, _plotting = _plot_fake(version="1.10.0")
    with pytest.raises(RuntimeError, match="exact liana==1.9.0"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake))

    class Drifted:
        def dotplot_by_sample(self, liana_res, sample_key="sample"):
            return _DrawResult()

    fake = types.SimpleNamespace(__version__="1.9.0", pl=Drifted())
    with pytest.raises(RuntimeError, match="signature drifted"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake))

    fake, _plotting = _plot_fake(wrong_return=True)
    with pytest.raises(RuntimeError, match="plotnine draw object"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake))

    fake, _plotting = _plot_fake(mutate=True)
    with pytest.raises(RuntimeError, match="modified its private result-table copy"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake))

    portable = artifact.portable()
    portable["table"].loc[0, "magnitude_rank"] = 0.99
    code = liana_dot_plot_code(parameters={k: v for k, v in _plot_parameters(fake).items() if k != "liana_module"})
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        namespace["run_liana_dot_plot"](portable, liana_module=fake)


def test_dot_plot_selection_contract_and_empty_or_missing_filters_fail_before_render():
    *_prefix, artifact, _summary = _run()
    fake, plotting = _plot_fake()
    with pytest.raises(ValueError, match="absent from the typed result"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake, source_labels=["missing"]))
    with pytest.raises(ValueError, match="selection is empty"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake, selection_threshold=0.001))
    with pytest.raises(MemoryError, match="nominal pixels"):
        render_liana_dot_plot(artifact, **_plot_parameters(fake, max_image_pixels=10))
    assert plotting.calls == []


def test_node_schemas_have_typed_ports_and_dynamic_branches():
    communication = OpenBioSingleCellLianaCommunication.define_schema()
    plot = OpenBioSingleCellLianaDotPlot.define_schema()
    assert [item.id for item in communication.inputs[:7]] == [
        "adata",
        "sample_key",
        "condition_key",
        "identity_key",
        "annotation_status",
        "organism",
        "method",
    ]
    assert [(item.display_name, item.io_type) for item in communication.outputs] == [
        ("result", LianaResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert [(item.id, item.io_type) for item in plot.inputs[:4]] == [
        ("result", LianaResultType.io_type),
        ("source_labels", "STRING"),
        ("target_labels", "STRING"),
        ("selection", "COMFY_DYNAMICCOMBO_V3"),
    ]
    assert [(item.display_name, item.io_type) for item in plot.outputs] == [
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_liana_nodes_flatten_dynamic_inputs_for_the_worker_protocol():
    communication = OpenBioSingleCellLianaCommunication.prepare_worker_arguments(
        {
            "adata": "ticket",
            "resource": {
                "resource": "local_resource",
                "resource_name": "snapshot",
                "resource_csv": "resource.csv",
                "resource_metadata_json": _resource_metadata("snapshot"),
            },
        }
    )
    plot = OpenBioSingleCellLianaDotPlot.prepare_worker_arguments(
        {
            "result": "ticket",
            "selection": {
                "selection": "cellphonedb",
                "max_cellphone_pvalue": 0.025,
            },
        }
    )

    assert communication == {
        "adata": "ticket",
        "resource_mode": "local_resource",
        "resource_name": "snapshot",
        "resource_csv": "resource.csv",
        "resource_metadata_json": _resource_metadata("snapshot"),
    }
    assert plot == {
        "result": "ticket",
        "selection_method": "cellphonedb",
        "selection_threshold": 0.025,
    }



def test_dot_plot_node_rejects_method_branch_mismatch_before_renderer_import():
    from openbio_singlecell.operations_communication import liana_dot_plot_owned

    *_prefix, artifact, _summary = _run("rank_aggregate")
    with pytest.raises(ValueError, match="does not match the typed result method"):
        liana_dot_plot_owned(
            artifact.portable(),
            selection_method="cellphonedb",
            selection_threshold=0.05,
        )


def test_current_liana_schemas_are_exact():
    current_communication = tuple(item.id for item in OpenBioSingleCellLianaCommunication.define_schema().inputs)
    assert current_communication == (
        "adata",
        "sample_key",
        "condition_key",
        "identity_key",
        "annotation_status",
        "organism",
        "method",
        "resource",
        "source",
        "expression_proportion",
        "min_cells_per_identity_sample",
        "permutations",
        "random_seed",
        "jobs",
        "max_output_rows",
        "max_working_memory_gib",
    )
    current_plot = tuple(item.id for item in OpenBioSingleCellLianaDotPlot.define_schema().inputs)
    assert current_plot == (
        "result",
        "source_labels",
        "target_labels",
        "selection",
        "top_n",
        "figure_width",
        "figure_height",
        "max_plot_rows",
        "max_image_pixels",
    )


def test_real_liana_1_9_sample_resolved_smoke_for_both_methods_and_plot():
    try:
        installed = importlib.metadata.version("liana")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("optional exact liana==1.9.0 smoke dependency is not installed")
    if installed != "1.9.0":
        pytest.skip(f"optional smoke requires exact liana==1.9.0, found {installed}")
    if int(pd.__version__.split(".", maxsplit=1)[0]) >= 3:
        pytest.skip("liana==1.9.0 declares pandas<3")

    import liana as li

    samples = np.repeat(["S1", "S2"], 16)
    identities = np.tile(np.repeat(["sender", "receiver"], 8), 2)
    observations = pd.DataFrame(
        {
            "sample": samples,
            "condition": np.where(samples == "S1", "control", "treated"),
            "cell_type": identities,
        },
        index=[f"real_cell_{index}" for index in range(32)],
    )
    rng = np.random.default_rng(9341)
    values = rng.gamma(shape=2.0, scale=0.8, size=(32, 4)) + 0.2
    values[identities == "sender", :2] += 2.5
    values[identities == "receiver", 2:] += 2.5
    adata = AnnData(
        X=np.log1p(values),
        obs=observations,
        var=pd.DataFrame(index=["LGALS9", "MET", "PTPRC", "CD44"]),
    )
    adata.raw = adata.copy()
    adata.uns["openbio_singlecell"] = {
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {"source": "X", "layer_name": "counts"},
            },
            "000001": {"operation": "log1p", "parameters": {}},
        }
    }
    metadata = _resource_metadata("consensus")

    artifacts = {}
    for method in ("rank_aggregate", "cellphonedb"):
        artifact, summary = run_liana_communication(
            adata,
            **_parameters(
                method,
                li,
                resource_metadata_json=metadata,
                min_cells_per_identity_sample=5,
                max_output_rows=10_000,
            ),
        )
        table, _provenance, _artifact_metadata = validate_liana_result(artifact)
        assert not table.empty
        assert set(table["sample"]) == {"S1", "S2"}
        assert summary["software_versions"]["liana"] == "1.9.0"
        artifacts[method] = artifact

    png, plot_summary = render_liana_dot_plot(
        artifacts["rank_aggregate"],
        source_labels=None,
        target_labels=None,
        selection_threshold=1.0,
        top_n=1,
        figure_width=5.0,
        figure_height=4.0,
        max_plot_rows=100,
        max_image_pixels=5_000_000,
        openbio_version="0.2.0",
        liana_module=li,
    )
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert plot_summary["key_results"]["selection"]["plotted_rows"] > 0
