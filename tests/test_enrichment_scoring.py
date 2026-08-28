from __future__ import annotations

import copy
import importlib
import json
import os
from types import SimpleNamespace

import pytest

from openbio_singlecell.node_types import AnnDataType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_enrichment import (
    OpenBioSingleCellAUCellScores,
    OpenBioSingleCellGenePanelScores,
    OpenBioSingleCellGSVAScores,
    OpenBioSingleCellPathwayScoreTTest,
)
from openbio_singlecell.pathway_score_contrast import PATHWAY_SCORE_CONTRAST_COLUMNS
from openbio_singlecell.score_artifact import build_score_artifact, store_score_artifact


def _outputs(node_output):
    return node_output.result


def _metadata_json():
    return json.dumps(
        {
            "name": "reviewed-test-sets",
            "version": "2026.08",
            "date": "2026-08-01",
            "organism": "human",
            "identifier_namespace": "HGNC symbol",
            "scope": "unit-test scientific contract",
            "license": "CC0-1.0",
            "citation": "OpenBio test fixture, 2026.",
        }
    )


def _logged_adata(science, *, cells=16, genes=100):
    rng = science.np.random.default_rng(7)
    counts = rng.poisson(lam=4.0, size=(cells, genes)).astype(float)
    counts[:, 0] += science.np.arange(cells) % 5
    obs = science.pd.DataFrame(index=[f"cell_{index:03d}" for index in range(cells)])
    var = science.pd.DataFrame(index=[f"G{index:03d}" for index in range(genes)])
    adata = science.ad.AnnData(science.sparse.csr_matrix(counts), obs=obs, var=var)
    science.sc.pp.normalize_total(adata, target_sum=10_000.0)
    science.sc.pp.log1p(adata)
    return adata


def _write_resource(path, *, genes=100):
    rows = ["geneset,genesymbol"]
    rows.extend(f"set_a,G{index:03d}" for index in range(0, min(8, genes)))
    rows.extend(f"set_b,G{index:03d}" for index in range(8, min(18, genes)))
    rows.extend(f"panel_x,G{index:03d}" for index in range(20, min(27, genes)))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_scoring_and_pathway_schemas_are_atomic():
    aucell = OpenBioSingleCellAUCellScores.GET_SCHEMA()
    assert [item.id for item in aucell.inputs] == [
        "adata",
        "gene_sets_file",
        "resource_metadata_json",
        "source",
        "source_column",
        "target_column",
        "min_targets",
        "n_top_features",
        "batch_size",
        "max_output_rows",
        "max_working_memory_gib",
        "output_key",
        "overwrite_existing",
    ]
    assert [(item.display_name, item.io_type) for item in aucell.outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    gsva = OpenBioSingleCellGSVAScores.GET_SCHEMA()
    assert [item.id for item in gsva.inputs] == [
        "adata",
        "gene_sets_file",
        "resource_metadata_json",
        "source",
        "source_column",
        "target_column",
        "min_targets",
        "kernel",
        "maxdiff",
        "absrnk",
        "tau",
        "max_working_memory_gib",
        "batch_size",
        "max_output_rows",
        "output_key",
        "overwrite_existing",
    ]
    panel = OpenBioSingleCellGenePanelScores.GET_SCHEMA()
    assert [item.id for item in panel.inputs] == [
        "adata",
        "gene_sets_file",
        "resource_metadata_json",
        "panel",
        "source",
        "source_column",
        "target_column",
        "output_key",
        "ctrl_size",
        "n_bins",
        "random_seed",
        "overwrite_existing",
    ]
    pathway = OpenBioSingleCellPathwayScoreTTest.GET_SCHEMA()
    assert [item.id for item in pathway.inputs] == [
        "adata",
        "sample_key",
        "condition_key",
        "annotation_key",
        "population",
        "condition_a",
        "condition_b",
        "score_key",
        "min_cells_per_sample_population",
        "min_samples_per_condition",
        "technical_batch_key",
        "confidence_level",
        "annotation_status",
    ]
    assert [(item.display_name, item.io_type) for item in pathway.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


@pytest.mark.parametrize(
    "node",
    [
        OpenBioSingleCellAUCellScores,
        OpenBioSingleCellGSVAScores,
        OpenBioSingleCellGenePanelScores,
    ],
)
def test_scoring_cache_fingerprint_tracks_same_size_same_mtime_resource_bytes(
    node, comfy_directories
):
    input_dir, _, _ = comfy_directories
    resource_path = input_dir / "cache-sets.csv"
    _write_resource(resource_path)
    before = node.fingerprint_inputs(resource_path.name)
    payload = resource_path.read_bytes().replace(b"set_a", b"set_z", 1)
    stat = resource_path.stat()
    resource_path.write_bytes(payload)
    os.utime(resource_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    after = node.fingerprint_inputs(resource_path.name)
    assert before != after


def test_real_aucell_22_generated_parity_provenance_and_immutability(
    science, comfy_directories
):
    pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    original_x = adata.X.copy()
    original_uns = copy.deepcopy(adata.uns)

    output, report, code = _outputs(
        OpenBioSingleCellAUCellScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
            n_top_features=20,
            output_key="aucell_test",
        )
    )
    assert output.obsm["aucell_test"].columns.tolist() == ["set_a", "set_b", "panel_x"]
    assert output.obsm["aucell_test"].index.equals(adata.obs_names)
    values = output.obsm["aucell_test"].to_numpy(dtype=float)
    assert science.np.isfinite(values).all()
    assert ((0.0 <= values) & (values <= 1.0)).all()
    assert report.summary["key_results"]["p_values_produced"] is False
    assert report.summary["parameters"]["resolved_n_up"] == 20
    assert report.summary["parameters"]["tie_policy"] == (
        "decoupler 2.2 fixed seed-0 feature permutation followed by scipy ordinal ranking; "
        "ties then follow the permuted feature order"
    )
    assert "seed-0 feature permutation" in json.dumps(report.summary)
    assert report.summary["software_versions"]["decoupler"] == "2.2.0"
    assert len(report.summary["key_results"]["resource"]["sha256"]) == 64
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<aucell-code>", "exec")
    namespace = {}
    exec(code, namespace)
    generated, generated_summary = namespace["score_aucell"](adata)
    science.pd.testing.assert_frame_equal(generated.obsm["aucell_test"], output.obsm["aucell_test"])
    assert generated.uns["openbio_singlecell_score_artifacts"] == output.uns[
        "openbio_singlecell_score_artifacts"
    ]
    assert generated_summary == report.summary
    assert "seed-0 feature permutation" in code
    assert (adata.X != original_x).nnz == 0
    assert adata.uns == original_uns
    assert "aucell_test" not in adata.obsm


def test_aucell_uses_exact_reviewed_decoupler_public_arguments(
    science, comfy_directories, monkeypatch
):
    decoupler = pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    calls = []

    def fake_aucell(data, net, tmin=5, raw=False, empty=True, bsize=250_000, verbose=False, n_up=None):
        calls.append(
            {
                "data": data,
                "net": net.copy(),
                "tmin": tmin,
                "raw": raw,
                "empty": empty,
                "bsize": bsize,
                "verbose": verbose,
                "n_up": n_up,
            }
        )
        sources = list(dict.fromkeys(net["source"].tolist()))
        data.obsm["score_aucell"] = science.pd.DataFrame(
            science.np.full((data.n_obs, len(sources)), 0.25),
            index=data.obs_names.copy(),
            columns=sources,
        )

    monkeypatch.setattr(decoupler.mt, "aucell", fake_aucell)
    output, report, _ = _outputs(
        OpenBioSingleCellAUCellScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            n_top_features=0,
            batch_size=17,
        )
    )
    assert len(calls) == 1
    assert calls[0]["tmin"] == 3
    assert calls[0]["raw"] is False
    assert calls[0]["empty"] is False
    assert calls[0]["bsize"] == 17
    assert calls[0]["verbose"] is False
    assert calls[0]["n_up"] is None
    assert calls[0]["net"].columns.tolist() == ["source", "target"]
    assert output.obsm["aucell_scores"].columns.tolist() == ["set_a", "set_b", "panel_x"]
    assert report.summary["parameters"]["resolved_n_up"] == 5


def test_aucell_public_wrapper_applies_decoupler_seed_zero_feature_tie_break(
    science, comfy_directories, monkeypatch
):
    pytest.importorskip("decoupler")
    decoupler_data = importlib.import_module("decoupler.pp.data")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    original_break_ties = decoupler_data._break_ties
    observed_orders = []

    def observe_break_ties(mat, features):
        shuffled_mat, shuffled_features = original_break_ties(mat, features)
        expected_indices = science.np.random.default_rng(seed=0).choice(
            science.np.arange(features.size), features.size, replace=False
        )
        assert science.np.array_equal(shuffled_features, features[expected_indices])
        observed_orders.append(tuple(map(str, shuffled_features)))
        return shuffled_mat, shuffled_features

    monkeypatch.setattr(decoupler_data, "_break_ties", observe_break_ties)
    _, report, _ = _outputs(
        OpenBioSingleCellAUCellScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            n_top_features=20,
        )
    )

    assert len(observed_orders) == 1
    assert report.summary["parameters"]["tie_policy"].startswith(
        "decoupler 2.2 fixed seed-0 feature permutation"
    )


@pytest.mark.parametrize(
    ("node", "method_name"),
    [
        (OpenBioSingleCellAUCellScores, "aucell"),
        (OpenBioSingleCellGSVAScores, "gsva"),
    ],
)
def test_decoupler_scoring_rejects_unreviewed_version_and_missing_public_method(
    science, comfy_directories, monkeypatch, node, method_name
):
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    original_import_module = importlib.import_module

    def replace_decoupler(module):
        monkeypatch.setattr(
            importlib,
            "import_module",
            lambda name: module if name == "decoupler" else original_import_module(name),
        )

    replace_decoupler(SimpleNamespace(__version__="1.9.9", mt=SimpleNamespace()))
    with pytest.raises(RuntimeError, match="requires the reviewed decoupler 2.2 public mt API"):
        node.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
        )

    replace_decoupler(SimpleNamespace(__version__="2.2.0", mt=SimpleNamespace()))
    with pytest.raises(RuntimeError, match=rf"requires callable decoupler\.mt\.{method_name}"):
        node.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
        )


def test_real_gsva_22_official_defaults_memory_guard_and_generated_parity(
    science, comfy_directories
):
    pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)

    output, report, code = _outputs(
        OpenBioSingleCellGSVAScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
            kernel="gaussian_normalized",
            output_key="gsva_test",
        )
    )
    scores = output.obsm["gsva_test"]
    assert scores.columns.tolist() == ["set_a", "set_b", "panel_x"]
    assert science.np.isfinite(scores.to_numpy(dtype=float)).all()
    assert report.summary["parameters"]["maxdiff"] is True
    assert report.summary["parameters"]["absrnk"] is False
    assert report.summary["parameters"]["kcdf"] == "gaussian"
    assert report.summary["key_results"]["resource"]["duplicate_pairs_removed"] == 0
    namespace = {}
    exec(compile(code, "<gsva-code>", "exec"), namespace)
    generated, generated_summary = namespace["score_gsva"](adata)
    science.pd.testing.assert_frame_equal(generated.obsm["gsva_test"], scores)
    assert generated_summary == report.summary

    with pytest.raises(ValueError, match="working-memory preflight"):
        OpenBioSingleCellGSVAScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
            max_working_memory_gib=1e-9,
        )


def test_gsva_uses_exact_reviewed_decoupler_arguments_and_empirical_kernel(
    science, comfy_directories, monkeypatch
):
    decoupler = pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    calls = []

    def fake_gsva(
        data,
        net,
        tmin=5,
        raw=False,
        empty=True,
        bsize=250_000,
        verbose=False,
        layer=None,
        kcdf="gaussian",
        maxdiff=True,
        absrnk=False,
        tau=1.0,
    ):
        calls.append(locals().copy())
        sources = list(dict.fromkeys(net["source"].tolist()))
        data.obsm["score_gsva"] = science.pd.DataFrame(
            science.np.linspace(-0.5, 0.5, data.n_obs * len(sources)).reshape(data.n_obs, len(sources)),
            index=data.obs_names.copy(),
            columns=sources,
        )

    monkeypatch.setattr(decoupler.mt, "gsva", fake_gsva)
    output, report, _ = _outputs(
        OpenBioSingleCellGSVAScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            kernel="empirical",
            maxdiff=False,
            absrnk=True,
            tau=1.5,
            batch_size=31,
        )
    )
    assert len(calls) == 1
    call = calls[0]
    assert call["layer"] is None
    assert call["kcdf"] is None
    assert call["maxdiff"] is False
    assert call["absrnk"] is True
    assert call["tau"] == 1.5
    assert call["raw"] is False and call["empty"] is False and call["verbose"] is False
    assert call["bsize"] == 31
    assert science.np.isfinite(output.obsm["gsva_scores"].to_numpy()).all()
    assert report.summary["parameters"]["kcdf"] is None


def test_real_poisson_gsva_requires_and_uses_explicit_count_state(science, comfy_directories):
    pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    rng = science.np.random.default_rng(11)
    counts = rng.poisson(4.0, size=(16, 100)).astype(float)
    adata = science.ad.AnnData(
        science.sparse.csr_matrix(counts),
        obs=science.pd.DataFrame(index=[f"cell_{index}" for index in range(16)]),
        var=science.pd.DataFrame(index=[f"G{index:03d}" for index in range(100)]),
    )
    adata.raw = adata.copy()
    science.sc.pp.normalize_total(adata, target_sum=10_000.0)
    science.sc.pp.log1p(adata)
    output, report, code = _outputs(
        OpenBioSingleCellGSVAScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "raw"},
            min_targets=3,
            kernel="poisson_counts",
            output_key="poisson_gsva",
        )
    )
    assert report.summary["key_results"]["expression"]["state"] == "counts"
    assert report.summary["parameters"]["kcdf"] == "poisson"
    assert science.np.isfinite(output.obsm["poisson_gsva"].to_numpy(dtype=float)).all()
    namespace = {}
    exec(compile(code, "<poisson-gsva-code>", "exec"), namespace)
    generated, generated_summary = namespace["score_gsva"](adata)
    science.pd.testing.assert_frame_equal(generated.obsm["poisson_gsva"], output.obsm["poisson_gsva"])
    assert generated_summary == report.summary

    with pytest.raises(ValueError, match="finite nonnegative integer count values"):
        OpenBioSingleCellGSVAScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": "X"},
            min_targets=3,
            kernel="poisson_counts",
        )


def test_real_scanpy_one_panel_exact_arguments_and_generated_parity(science, comfy_directories):
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)

    output, report, code = _outputs(
        OpenBioSingleCellGenePanelScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            panel="panel_x",
            source={"source": "X"},
            output_key="panel_test",
            ctrl_size=0,
            n_bins=5,
            random_seed=19,
        )
    )
    assert "panel_test" in output.obs
    assert science.np.isfinite(output.obs["panel_test"].to_numpy(dtype=float)).all()
    assert report.summary["parameters"]["panel"] == "panel_x"
    assert report.summary["parameters"]["resolved_ctrl_size"] == 7
    assert report.summary["parameters"]["scanpy_fixed_arguments"] == {
        "ctrl_as_ref": False,
        "copy": False,
        "use_raw": False,
        "layer": None,
    }
    namespace = {}
    exec(compile(code, "<panel-code>", "exec"), namespace)
    generated, generated_summary = namespace["score_gene_panel"](adata)
    science.pd.testing.assert_series_equal(generated.obs["panel_test"], output.obs["panel_test"])
    assert generated_summary == report.summary
    assert "panel_test" not in adata.obs


def test_panel_uses_exact_scanpy_arguments_and_isolates_global_rng(
    science, comfy_directories, monkeypatch
):
    import scanpy

    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    calls = []

    def fake_score_genes(
        adata,
        gene_list,
        *,
        ctrl_as_ref=True,
        ctrl_size=50,
        gene_pool=None,
        n_bins=25,
        score_name="score",
        random_state=0,
        copy=False,
        use_raw=None,
        layer=None,
    ):
        calls.append(
            {
                "gene_list": list(gene_list),
                "ctrl_as_ref": ctrl_as_ref,
                "ctrl_size": ctrl_size,
                "gene_pool": list(gene_pool),
                "n_bins": n_bins,
                "score_name": score_name,
                "random_state": random_state,
                "copy": copy,
                "use_raw": use_raw,
                "layer": layer,
            }
        )
        science.np.random.seed(999)
        adata.obs[score_name] = science.np.arange(adata.n_obs, dtype=float)

    monkeypatch.setattr(scanpy.tl, "score_genes", fake_score_genes)
    before = science.np.random.get_state()
    output, report, _ = _outputs(
        OpenBioSingleCellGenePanelScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            panel="panel_x",
            ctrl_size=0,
            n_bins=7,
            random_seed=23,
        )
    )
    after = science.np.random.get_state()
    assert before[0] == after[0] and science.np.array_equal(before[1], after[1]) and before[2:] == after[2:]
    assert calls == [
        {
            "gene_list": [f"G{index:03d}" for index in range(20, 27)],
            "ctrl_as_ref": False,
            "ctrl_size": 7,
            "gene_pool": [f"G{index:03d}" for index in range(100)],
            "n_bins": 7,
            "score_name": "__openbio_gene_panel_score__",
            "random_state": 23,
            "copy": False,
            "use_raw": False,
            "layer": None,
        }
    ]
    assert output.obs["panel_score"].tolist() == list(map(float, range(16)))
    assert report.summary["parameters"]["scanpy_changed_numpy_global_rng_state_inside_isolation"] is True


def test_score_output_collision_requires_opt_in_and_generated_code_pins_resource_bytes(
    science, comfy_directories
):
    pytest.importorskip("decoupler")
    input_dir, _, _ = comfy_directories
    resource_path = input_dir / "sets.csv"
    _write_resource(resource_path)
    adata = _logged_adata(science)
    output, _, code = _outputs(
        OpenBioSingleCellAUCellScores.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            output_key="scores",
        )
    )
    with pytest.raises(ValueError, match="overwrite_existing"):
        OpenBioSingleCellAUCellScores.execute(
            output,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            output_key="scores",
        )
    replaced, replacement_report, _ = _outputs(
        OpenBioSingleCellAUCellScores.execute(
            output,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            min_targets=3,
            output_key="scores",
            overwrite_existing=True,
        )
    )
    assert any("Explicitly replaced" in warning for warning in replacement_report.summary["warnings"])
    assert "scores" in replaced.obsm

    original_text = resource_path.read_text(encoding="utf-8")
    changed_text = original_text.replace("set_a,G000", "set_a,G099")
    assert len(changed_text.encode()) == len(original_text.encode())
    resource_path.write_text(changed_text, encoding="utf-8")
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="resource fingerprint changed"):
        namespace["score_aucell"](adata)


@pytest.mark.parametrize(
    "node,extra,entrypoint,storage",
    [
        (OpenBioSingleCellAUCellScores, {"min_targets": 3}, "score_aucell", "obsm"),
        (
            OpenBioSingleCellGSVAScores,
            {"min_targets": 3, "kernel": "gaussian_normalized"},
            "score_gsva",
            "obsm",
        ),
        (
            OpenBioSingleCellGenePanelScores,
            {"panel": "panel_x"},
            "score_gene_panel",
            "obs",
        ),
    ],
)
@pytest.mark.parametrize("source_kind", ["X", "raw"])
def test_scoring_allows_explicit_count_like_expression_with_cautious_generated_parity(
    science, comfy_directories, node, extra, entrypoint, storage, source_kind
):
    input_dir, _, _ = comfy_directories
    _write_resource(input_dir / "sets.csv")
    adata = _logged_adata(science)
    counts = science.np.rint(science.np.expm1(adata.X.toarray())).astype(int)
    adata.raw = science.ad.AnnData(
        counts.copy(),
        obs=science.pd.DataFrame(index=adata.obs_names.copy()),
        var=science.pd.DataFrame(index=adata.var_names.copy()),
    )
    if source_kind == "X":
        adata.X = counts.copy()
    adata.uns.pop("log1p", None)
    output, report, code = _outputs(
        node.execute(
            adata,
            gene_sets_file="sets.csv",
            resource_metadata_json=_metadata_json(),
            source={"source": source_kind},
            **extra,
        )
    )
    assert any("count-like" in warning for warning in report.summary["warnings"])
    assert report.summary["key_results"]["expression"]["state"] in {"counts", "unknown"}
    namespace = {}
    exec(compile(code, f"<{entrypoint}-{source_kind}-counts>", "exec"), namespace)
    generated, generated_summary = namespace[entrypoint](adata)
    assert generated_summary == report.summary
    if storage == "obs":
        science.pd.testing.assert_series_equal(generated.obs["panel_score"], output.obs["panel_score"])
    else:
        key = "aucell_scores" if entrypoint == "score_aucell" else "gsva_scores"
        science.pd.testing.assert_frame_equal(generated.obsm[key], output.obsm[key])


def _pathway_adata_with_artifact(science):
    samples = [f"s{index}" for index in range(8)]
    conditions = {sample: "control" if index < 4 else "treated" for index, sample in enumerate(samples)}
    batches = {sample: f"b{index % 2}" for index, sample in enumerate(samples)}
    obs_rows = []
    score_rows = []
    obs_names = []
    sample_values = {
        "p1": [0.0, 1.0, 2.5, 3.0, 4.0, 5.5, 6.0, 8.0],
        "p2": [8.0, 6.0, 5.0, 3.0, 4.0, 3.0, 1.0, 0.0],
    }
    for sample_index, sample in enumerate(samples):
        for cell_index in range(10):
            obs_names.append(f"{sample}_c{cell_index}")
            obs_rows.append(
                {
                    "sample": sample,
                    "condition": conditions[sample],
                    "cell_type": "T",
                    "batch": batches[sample],
                }
            )
            score_rows.append(
                {
                    "p1": sample_values["p1"][sample_index] + (cell_index - 4.5) * 0.01,
                    "p2": sample_values["p2"][sample_index] - (cell_index - 4.5) * 0.01,
                }
            )
    obs = science.pd.DataFrame(obs_rows, index=obs_names)
    adata = science.ad.AnnData(
        science.np.ones((len(obs), 3), dtype=float),
        obs=obs,
        var=science.pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    scores = science.pd.DataFrame(score_rows, index=obs_names)
    adata.obsm["scores"] = scores
    references = [
        {
            "citation": "Aibar S, et al. SCENIC. Nature Methods. 2017.",
            "url": "https://doi.org/10.1038/nmeth.4463",
            "kind": "method",
            "doi": "10.1038/nmeth.4463",
        }
    ]
    artifact = build_score_artifact(
        frame=scores,
        feature_names=list(adata.var_names),
        producer_node="OpenBioSingleCellAUCellScores",
        method="AUCell",
        storage="obsm",
        score_key="scores",
        resource={"sha256": "a" * 64, "metadata": {"name": "fixture"}},
        expression={"source": "X", "state": "logged", "content_sha256": "b" * 64},
        parameters={"min_targets": 2},
        references=references,
        np=science.np,
        pd=science.pd,
    )
    store_score_artifact(adata, artifact)
    return adata


def test_pathway_contrast_uses_sample_means_full_bh_family_and_generated_parity(science):
    adata = _pathway_adata_with_artifact(science)
    original = adata.copy()
    result, report, code = _outputs(
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
            min_cells_per_sample_population=10,
            min_samples_per_condition=3,
            technical_batch_key="batch",
            annotation_status="curated",
        )
    )
    table = result.table
    assert table.columns.tolist() == PATHWAY_SCORE_CONTRAST_COLUMNS
    assert table["pathway"].nunique() == 2
    assert table["n_samples_a"].tolist() == [4, 4]
    assert table["n_samples_b"].tolist() == [4, 4]
    assert table["n_cells_a"].tolist() == [40, 40]
    assert table["n_cells_b"].tolist() == [40, 40]
    assert science.np.isfinite(table.select_dtypes("number").to_numpy(dtype=float)).all()
    expected_bh = science.np.asarray(
        __import__("scipy").stats.false_discovery_control(table["p_value"].to_numpy(), method="bh")
    )
    science.np.testing.assert_allclose(table["p_adjusted"], expected_bh)
    p1 = table.set_index("pathway").loc["p1"]
    expected = __import__("scipy").stats.ttest_ind(
        [0.0, 1.0, 2.5, 3.0],
        [4.0, 5.5, 6.0, 8.0],
        equal_var=False,
        nan_policy="raise",
        alternative="two-sided",
    )
    expected_interval = expected.confidence_interval(confidence_level=0.95)
    assert p1["mean_difference"] == pytest.approx(-4.25)
    assert p1["t_statistic"] == pytest.approx(expected.statistic)
    assert p1["degrees_of_freedom"] == pytest.approx(expected.df)
    assert p1["p_value"] == pytest.approx(expected.pvalue)
    assert p1["ci_low"] == pytest.approx(expected_interval.low)
    assert p1["ci_high"] == pytest.approx(expected_interval.high)
    assert report.summary["key_results"]["inference_unit"] == "Sample"
    assert report.summary["key_results"]["technical_batch_audit"]["status"] == (
        "overlap_present_no_adjustment"
    )
    assert report.summary["key_results"]["retained_sample_support_preview"] == [
        {
            "sample": f"s{index}",
            "condition": "control" if index < 4 else "treated",
            "cells": 10,
            "technical_batches": [f"b{index % 2}"],
        }
        for index in range(8)
    ]
    assert "biological Sample" in report.summary["methods"]
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(compile(code, "<pathway-score-code>", "exec"), namespace)
    generated_table, generated_summary = namespace["contrast_pathway_scores"](adata)
    science.pd.testing.assert_frame_equal(generated_table, table)
    assert generated_summary == report.summary
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    science.pd.testing.assert_frame_equal(adata.obsm["scores"], original.obsm["scores"])
    assert adata.uns == original.uns


def test_pathway_contrast_rejects_tampered_scores(science):
    adata = _pathway_adata_with_artifact(science)
    adata.obsm["scores"].iloc[0, 0] += 1.0
    with pytest.raises(ValueError, match="immutable provenance"):
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
            technical_batch_key="batch",
        )



def test_pathway_contrast_reports_batch_confounding_and_sample_batch_multiplicity(science):
    adata = _pathway_adata_with_artifact(science)
    adata.obs["confounded"] = adata.obs["condition"].map({"control": "b1", "treated": "b2"})
    result, report, code = _outputs(
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
            technical_batch_key="confounded",
        )
    )
    audit = report.summary["key_results"]["technical_batch_audit"]
    assert audit["formal_interpretation_invalid"] is True
    assert audit["formal_interpretation_invalid_reasons"] == [
        "perfect_condition_technical_batch_confounding"
    ]
    assert any("formal Condition-effect interpretation is invalid" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(compile(code, "<pathway-score-perfect-confounding>", "exec"), namespace)
    generated_table, generated_summary = namespace["contrast_pathway_scores"](adata)
    science.pd.testing.assert_frame_equal(generated_table, result.table)
    assert generated_summary == report.summary

    adata = _pathway_adata_with_artifact(science)
    adata.obs.loc[adata.obs["sample"] == "s0", "batch"] = ["b0"] * 9 + ["b9"]
    result, report, code = _outputs(
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
            technical_batch_key="batch",
        )
    )
    audit = report.summary["key_results"]["technical_batch_audit"]
    assert audit["formal_interpretation_invalid"] is True
    assert audit["samples_spanning_batches"] == [
        {"sample": "s0", "technical_batches": ["b0", "b9"]}
    ]
    assert "sample_spans_multiple_technical_batches" in audit["formal_interpretation_invalid_reasons"]
    namespace = {}
    exec(compile(code, "<pathway-score-batch-multiplicity>", "exec"), namespace)
    generated_table, generated_summary = namespace["contrast_pathway_scores"](adata)
    science.pd.testing.assert_frame_equal(generated_table, result.table)
    assert generated_summary == report.summary


def test_pathway_contrast_allows_two_samples_per_arm_with_warning_and_generated_parity(science):
    adata = _pathway_adata_with_artifact(science)
    keep_samples = {"s0", "s1", "s4", "s5"}
    adata = adata[adata.obs["sample"].isin(keep_samples)].copy()
    scores = adata.obsm["scores"]
    artifact = build_score_artifact(
        frame=scores,
        feature_names=list(adata.var_names),
        producer_node="OpenBioSingleCellAUCellScores",
        method="AUCell",
        storage="obsm",
        score_key="scores",
        resource={"sha256": "a" * 64, "metadata": {"name": "fixture"}},
        expression={"source": "X", "state": "logged", "content_sha256": "b" * 64},
        parameters={"min_targets": 2},
        references=[
            {
                "citation": "Aibar S, et al. SCENIC. Nature Methods. 2017.",
                "url": "https://doi.org/10.1038/nmeth.4463",
                "kind": "method",
                "doi": "10.1038/nmeth.4463",
            }
        ],
        np=science.np,
        pd=science.pd,
    )
    store_score_artifact(adata, artifact)
    result, report, code = _outputs(
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
            min_samples_per_condition=2,
            technical_batch_key="batch",
        )
    )
    table = result.table
    assert table["n_samples_a"].eq(2).all()
    assert table["n_samples_b"].eq(2).all()
    assert any("only two retained biological Samples" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(compile(code, "<pathway-score-two-samples>", "exec"), namespace)
    generated_table, generated_summary = namespace["contrast_pathway_scores"](adata)
    science.pd.testing.assert_frame_equal(generated_table, table)
    assert generated_summary == report.summary


def test_pathway_contrast_fails_closed_on_zero_variance_family_member(science):
    adata = _pathway_adata_with_artifact(science)
    adata.obsm["scores"]["constant"] = 1.0
    artifact = build_score_artifact(
        frame=adata.obsm["scores"],
        feature_names=list(adata.var_names),
        producer_node="OpenBioSingleCellAUCellScores",
        method="AUCell",
        storage="obsm",
        score_key="scores",
        resource={"sha256": "a" * 64, "metadata": {"name": "fixture"}},
        expression={"source": "X", "state": "logged", "content_sha256": "b" * 64},
        parameters={"min_targets": 2},
        references=[
            {
                "citation": "Aibar S et al. SCENIC.",
                "url": "https://doi.org/10.1038/nmeth.4463",
                "kind": "method",
                "doi": "10.1038/nmeth.4463",
            }
        ],
        np=science.np,
        pd=science.pd,
    )
    store_score_artifact(adata, artifact)
    with pytest.raises(ValueError, match="zero Sample-level variance in both Conditions"):
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
        )


def test_pathway_contrast_rejects_sample_mapping_conflicts_and_insufficient_population_support(science):
    adata = _pathway_adata_with_artifact(science)
    adata.obs.loc["s0_c0", "condition"] = "treated"
    with pytest.raises(ValueError, match="one Condition per Sample"):
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
        )

    adata = _pathway_adata_with_artifact(science)
    adata.obs.loc[adata.obs["sample"].isin(["s0", "s1"]), "cell_type"] = "B"
    with pytest.raises(ValueError, match="at least 3 retained Samples per Condition"):
        OpenBioSingleCellPathwayScoreTTest.execute(
            adata,
            population="T",
            condition_a="control",
            condition_b="treated",
            score_key="scores",
        )
