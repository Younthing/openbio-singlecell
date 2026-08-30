from __future__ import annotations

import hashlib
import json

import pytest

from openbio_singlecell.cell_cycle import (
    REGEV_HUMAN_G2M_GENES,
    REGEV_HUMAN_RESOURCE_SHA256,
    REGEV_HUMAN_S_GENES,
    analyze_cell_cycle_score,
)
from openbio_singlecell.node_types import AnnDataType, SummaryResultType
from openbio_singlecell.nodes_trajectory import OpenBioSingleCellCellCycleScore
from openbio_singlecell.operations_trajectory import cell_cycle_score_owned


def _cell_cycle_adata(science, *, sparse=False, genes=None):
    gene_names = list(genes or [*REGEV_HUMAN_S_GENES, *REGEV_HUMAN_G2M_GENES, *[f"CTRL{i}" for i in range(120)]])
    rng = science.np.random.default_rng(42)
    counts = rng.poisson(2.0, size=(36, len(gene_names))).astype(float)
    counts[:, : min(12, len(gene_names))] += science.np.repeat([0.0, 2.0, 5.0], 12)[:, None]
    counts[:, min(43, len(gene_names)) : min(55, len(gene_names))] += science.np.repeat(
        [5.0, 2.0, 0.0], 12
    )[:, None]
    counts[:, 0] += 1.0
    totals = counts.sum(axis=1, keepdims=True)
    logged = science.np.log1p(counts / totals * 10_000.0)
    obs = science.pd.DataFrame(index=[f"cell_{index:03d}" for index in range(counts.shape[0])])
    var = science.pd.DataFrame(index=gene_names)
    count_matrix = science.sparse.csr_matrix(counts) if sparse else counts
    logged_matrix = science.sparse.csr_matrix(logged) if sparse else logged
    adata = science.ad.AnnData(count_matrix.copy(), obs=obs, var=var)
    adata.layers["counts"] = count_matrix.copy()
    adata.raw = adata.copy()
    adata.X = logged_matrix.copy()
    adata.layers["log1p_norm"] = logged_matrix.copy()
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "analysis_history": {
            "000000": {
                "operation": "snapshot_expression",
                "parameters": {"destination": "layer_and_raw", "layer_name": "counts"},
            },
            "000001": {
                "operation": "normalize_to_layer",
                "parameters": {"output_layer": "log1p_norm", "transform": "log1p"},
            },
        }
    }
    adata.uns["sentinel"] = {"unchanged": [1, 2, 3]}
    return adata


def _run(adata, **overrides):
    parameters = {
        "source_kind": "layer",
        "layer_name": "log1p_norm",
        "gene_set_source": "regev_human_97",
        "organism": "human",
        "s_genes": None,
        "g2m_genes": None,
        "output_prefix": "cell_cycle",
        "overwrite_existing": False,
        "random_seed": 7,
    }
    parameters.update(overrides)
    return analyze_cell_cycle_score(adata, **parameters)


def test_bundled_cell_cycle_resource_has_exact_43_54_identity():
    assert len(REGEV_HUMAN_S_GENES) == 43
    assert len(REGEV_HUMAN_G2M_GENES) == 54
    assert len(set(REGEV_HUMAN_S_GENES) | set(REGEV_HUMAN_G2M_GENES)) == 97
    payload = "\n".join([*REGEV_HUMAN_S_GENES, *REGEV_HUMAN_G2M_GENES]) + "\n"
    assert hashlib.sha256(payload.encode()).hexdigest() == REGEV_HUMAN_RESOURCE_SHA256


def test_cell_cycle_schema_is_atomic_and_uses_transformed_source():
    schema = OpenBioSingleCellCellCycleScore.define_schema()
    assert [input_.id for input_ in schema.inputs] == [
        "adata",
        "source",
        "gene_set_source",
        "organism",
        "output_prefix",
        "overwrite_existing",
        "random_seed",
    ]
    assert [(output.display_name, output.io_type) for output in schema.outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_cell_cycle_matches_explicit_scanpy_call_on_worker_owned_input(science):
    adata = _cell_cycle_adata(science)
    before = adata.copy()
    output, summary = _run(adata)
    assert output is adata

    work = science.ad.AnnData(
        adata.layers["log1p_norm"].copy(),
        obs=science.pd.DataFrame(index=adata.obs_names.copy()),
        var=science.pd.DataFrame(index=adata.var_names.copy()),
    )
    science.sc.tl.score_genes_cell_cycle(
        work,
        s_genes=list(REGEV_HUMAN_S_GENES),
        g2m_genes=list(REGEV_HUMAN_G2M_GENES),
        gene_pool=list(adata.var_names),
        n_bins=25,
        ctrl_as_ref=True,
        use_raw=False,
        random_state=7,
        copy=False,
    )
    science.np.testing.assert_allclose(output.obs["cell_cycle_s_score"], work.obs["S_score"])
    science.np.testing.assert_allclose(output.obs["cell_cycle_g2m_score"], work.obs["G2M_score"])
    assert output.obs["cell_cycle_phase"].astype(str).tolist() == work.obs["phase"].astype(str).tolist()
    assert summary["status"] == "descriptive_cell_cycle_program_score"
    assert summary["key_results"]["resource"]["sha256"] == REGEV_HUMAN_RESOURCE_SHA256
    assert sum(summary["key_results"]["phase_counts"].values()) == adata.n_obs
    assert summary["parameters"]["control_size"] == 43
    json.dumps(summary, allow_nan=False)

    science.np.testing.assert_array_equal(adata.X, before.X)
    science.np.testing.assert_array_equal(adata.layers["counts"], before.layers["counts"])
    science.np.testing.assert_array_equal(adata.layers["log1p_norm"], before.layers["log1p_norm"])
    assert adata.uns["sentinel"] == before.uns["sentinel"]


def test_cell_cycle_node_and_generated_code_are_equivalent(science):
    adata = _cell_cycle_adata(science)
    reproduction_input = adata.copy()
    output, report, code = cell_cycle_score_owned(
        adata,
        source={"source": "layer", "layer_name": "log1p_norm"},
        gene_set_source={"gene_set_source": "regev_human_97"},
        organism="human",
        output_prefix="cc",
        overwrite_existing=False,
        random_seed=11,
    )
    assert "from openbio_singlecell" not in code
    compile(code, "<cell-cycle-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["score_cell_cycle"](reproduction_input)
    science.pd.testing.assert_series_equal(output.obs["cc_s_score"], reproduced.obs["cc_s_score"])
    science.pd.testing.assert_series_equal(output.obs["cc_g2m_score"], reproduced.obs["cc_g2m_score"])
    science.pd.testing.assert_series_equal(output.obs["cc_phase"], reproduced.obs["cc_phase"])
    assert reproduced_summary == report.summary
    assert report.summary["software_versions"]["scanpy"] == "1.12.3"


def test_cell_cycle_sparse_and_dense_are_numerically_equal(science):
    dense_output, dense_summary = _run(_cell_cycle_adata(science))
    sparse_output, sparse_summary = _run(_cell_cycle_adata(science, sparse=True))
    science.np.testing.assert_allclose(dense_output.obs["cell_cycle_s_score"], sparse_output.obs["cell_cycle_s_score"])
    science.np.testing.assert_allclose(
        dense_output.obs["cell_cycle_g2m_score"], sparse_output.obs["cell_cycle_g2m_score"]
    )
    assert dense_summary["key_results"]["phase_counts"] == sparse_summary["key_results"]["phase_counts"]


def test_cell_cycle_custom_program_and_collision_policy(science):
    adata = _cell_cycle_adata(science)
    s_genes = ",".join(REGEV_HUMAN_S_GENES[:10])
    g2m_genes = ",".join(REGEV_HUMAN_G2M_GENES[:10])
    output, summary = _run(
        adata,
        gene_set_source="custom",
        organism="mouse",
        s_genes=s_genes,
        g2m_genes=g2m_genes,
    )
    assert summary["key_results"]["resource"]["organism"] == "mouse"
    assert output.obs["cell_cycle_phase"].notna().all()

    _, override_summary = _run(_cell_cycle_adata(science), organism="mouse")
    assert override_summary["key_results"]["resource"]["organism_compatibility"] == "nonhuman_expert_override"
    assert any("no orthology mapping" in warning for warning in override_summary["warnings"])

    adata = _cell_cycle_adata(science)
    adata.obs["cell_cycle_s_score"] = 999.0
    with pytest.raises(ValueError, match="already exist"):
        _run(adata, gene_set_source="custom", s_genes=s_genes, g2m_genes=g2m_genes)
    overwritten, _ = _run(
        adata,
        gene_set_source="custom",
        s_genes=s_genes,
        g2m_genes=g2m_genes,
        overwrite_existing=True,
    )
    assert not (overwritten.obs["cell_cycle_s_score"] == 999.0).all()


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        (
            {
                "gene_set_source": "custom",
                "s_genes": ",".join([*REGEV_HUMAN_S_GENES[:9], REGEV_HUMAN_S_GENES[0]]),
                "g2m_genes": ",".join(REGEV_HUMAN_G2M_GENES[:10]),
            },
            "duplicate",
        ),
    ],
)
def test_cell_cycle_rejects_invalid_resource_contract(science, overrides, message):
    with pytest.raises(ValueError, match=message):
        _run(_cell_cycle_adata(science), **overrides)


def test_cell_cycle_reports_low_coverage_and_cross_phase_overlap(science):
    genes = [*REGEV_HUMAN_S_GENES[:21], *REGEV_HUMAN_G2M_GENES[:27], *[f"CTRL{i}" for i in range(80)]]
    output, summary = _run(_cell_cycle_adata(science, genes=genes))
    assert output.n_obs == 36
    assert any("coverage is below" in warning for warning in summary["warnings"])

    overlap_gene = REGEV_HUMAN_S_GENES[0]
    _, overlap_summary = _run(
        _cell_cycle_adata(science),
        gene_set_source="custom",
        organism="other",
        s_genes=",".join(REGEV_HUMAN_S_GENES[:10]),
        g2m_genes=",".join([overlap_gene, *REGEV_HUMAN_G2M_GENES[:9]]),
    )
    assert overlap_summary["key_results"]["resource"]["cross_phase_overlap"] == [overlap_gene]
    assert any("overlap" in warning for warning in overlap_summary["warnings"])


def test_cell_cycle_accepts_explicit_raw_and_signed_expert_sources(science):
    adata = _cell_cycle_adata(science)
    raw_output, raw_summary = _run(adata.copy(), source_kind="raw", layer_name=None)
    assert raw_output.n_obs == adata.n_obs
    assert raw_summary["key_results"]["expression_source"] == "raw.X"
    assert "expression_state" not in raw_summary["key_results"]
    assert "expression_state_evidence" not in raw_summary["key_results"]
    assert not any("Raw was selected explicitly" in warning for warning in raw_summary["warnings"])

    reproduction_input = adata.copy()
    node_output, node_report, raw_code = cell_cycle_score_owned(
        adata.copy(),
        source={"source": "raw"},
        gene_set_source={"gene_set_source": "regev_human_97"},
        organism="human",
        output_prefix="raw_cc",
        overwrite_existing=False,
        random_seed=7,
    )
    namespace = {}
    exec(raw_code, namespace)
    reproduced, reproduced_summary = namespace["score_cell_cycle"](reproduction_input)
    science.pd.testing.assert_series_equal(node_output.obs["raw_cc_s_score"], reproduced.obs["raw_cc_s_score"])
    science.pd.testing.assert_series_equal(node_output.obs["raw_cc_phase"], reproduced.obs["raw_cc_phase"])
    assert reproduced_summary == node_report.summary

    signed = _cell_cycle_adata(science)
    signed.layers["signed"] = signed.layers["log1p_norm"] - 4.0
    signed_output, signed_summary = _run(signed, source_kind="layer", layer_name="signed")
    assert signed_output.n_obs == signed.n_obs
    assert "expression_state" not in signed_summary["key_results"]
    assert any("negative values" in warning for warning in signed_summary["warnings"])


def test_cell_cycle_retains_required_value_and_gene_availability_failures(science):
    nonfinite = _cell_cycle_adata(science)
    nonfinite.layers["log1p_norm"][0, 0] = float("nan")
    with pytest.raises(ValueError, match="NaN or infinity"):
        _run(nonfinite)

    genes = [REGEV_HUMAN_S_GENES[0], *[f"CTRL{i}" for i in range(120)]]
    with pytest.raises(ValueError, match="at least one observed gene in each phase"):
        _run(_cell_cycle_adata(science, genes=genes))


def test_cell_cycle_restores_global_numpy_rng_even_if_backend_changes_it(science, monkeypatch):
    original = science.sc.tl.score_genes_cell_cycle

    def mutating_backend(*args, **kwargs):
        science.np.random.seed(999)
        return original(*args, **kwargs)

    monkeypatch.setattr(science.sc.tl, "score_genes_cell_cycle", mutating_backend)
    science.np.random.seed(12345)
    before = science.np.random.get_state()
    _run(_cell_cycle_adata(science))
    after = science.np.random.get_state()
    assert before[0] == after[0]
    science.np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]


def test_cell_cycle_rejects_malformed_backend_output(science, monkeypatch):
    original = science.sc.tl.score_genes_cell_cycle

    def malformed(*args, **kwargs):
        result = original(*args, **kwargs)
        args[0].obs.loc[args[0].obs.index[0], "S_score"] = float("nan")
        return result

    monkeypatch.setattr(science.sc.tl, "score_genes_cell_cycle", malformed)
    with pytest.raises(RuntimeError, match="malformed"):
        _run(_cell_cycle_adata(science))
