from __future__ import annotations

from types import SimpleNamespace

import pytest

from openbio_singlecell.expression_source import ExpressionSource, ExpressionSourceSpec
from openbio_singlecell.extension import NODE_CLASSES

RAW_SENTINEL = object()


def _adata(*, raw=RAW_SENTINEL):
    return SimpleNamespace(X="X", raw=raw, layers={"counts": "counts"})


def test_expression_source_input_only_nests_layer_name_under_layer_option():
    spec = ExpressionSourceSpec(
        description="Test expression source",
        default="X",
        include_raw=True,
        layer_default="counts",
    )
    source = spec.input()

    assert source.get_io_type() == "COMFY_DYNAMICCOMBO_V3"
    assert [(option.key, [item.id for item in option.inputs]) for option in source.options] == [
        ("X", []),
        ("raw", []),
        ("layer", ["layer_name"]),
    ]


@pytest.mark.parametrize(
    ("value", "kind", "layer_name", "matrix"),
    [
        ({"source": "X"}, "X", None, "X"),
        ({"source": "raw"}, "raw", None, "raw"),
        ({"source": "layer", "layer_name": "counts"}, "layer", "counts", "counts"),
    ],
)
def test_resolve_expression_source_exposes_one_coherent_selection(value, kind, layer_name, matrix):
    adata = _adata(raw=SimpleNamespace(X="raw"))
    spec = ExpressionSourceSpec(description="Test expression source", include_raw=True)

    source = spec.resolve(adata, value)

    assert source.kind == kind
    assert source.layer_name == layer_name
    assert source.matrix(adata) == matrix


def test_resolve_expression_source_validates_only_the_selected_branch():
    adata = _adata(raw=None)
    spec = ExpressionSourceSpec(description="Test expression source", include_raw=True)

    assert spec.resolve(adata, {"source": "X"}).kind == "X"
    with pytest.raises(ValueError, match="layer not found"):
        spec.resolve(adata, {"source": "layer", "layer_name": "missing"})
    with pytest.raises(ValueError, match="adata.raw is unavailable"):
        spec.resolve(adata, {"source": "raw"})


def test_expression_source_value_object_enforces_and_canonicalizes_its_invariants():
    source = ExpressionSource("layer", "  counts  ", "  source_layer  ")

    assert source.layer_name == "counts"
    assert source.layer_input_id == "source_layer"
    with pytest.raises(ValueError, match="Unsupported expression source"):
        ExpressionSource("bogus")
    with pytest.raises(ValueError, match="requires a layer name"):
        ExpressionSource("layer", "   ")
    with pytest.raises(TypeError, match="string layer name"):
        ExpressionSource("layer", 42)


def test_resolve_expression_source_rejects_non_string_layer_values():
    spec = ExpressionSourceSpec(description="Test expression source")

    with pytest.raises(TypeError, match="layer must be a string"):
        spec.resolve(_adata(), {"source": "layer", "layer_name": 42})


def test_expression_source_spec_owns_the_runtime_default_and_parameter_shape():
    spec = ExpressionSourceSpec(
        description="Counts source",
        default="layer",
        layer_input_id="counts_layer",
        layer_default="counts",
    )
    adata = _adata()

    source = spec.resolve(adata, None)

    assert source.kind == "layer"
    assert source.matrix(adata) == "counts"
    assert source.parameters() == {"source": "layer", "counts_layer": "counts"}
    assert spec.resolve(adata, {"source": "X"}).parameters() == {"source": "X"}
    raw_spec = ExpressionSourceSpec(description="Raw-capable source", include_raw=True)
    assert raw_spec.resolve(adata, {"source": "raw"}).parameters() == {"source": "raw"}
    with pytest.raises(ValueError, match="Unsupported counts source"):
        spec.resolve(adata, {})


def test_every_expression_source_node_uses_the_native_dynamic_combo_contract():
    source_nodes = [
        (node_class, item)
        for node_class in NODE_CLASSES
        for item in node_class.GET_SCHEMA().inputs
        if item.id == "source"
    ]

    assert len(source_nodes) == 27
    assert {item.get_io_type() for _, item in source_nodes} == {"COMFY_DYNAMICCOMBO_V3"}
    for node_class, source in source_nodes:
        spec = node_class.EXPRESSION_SOURCE
        assert isinstance(spec, ExpressionSourceSpec)
        assert source.options[0].key == spec.default
        assert source.options[0].key != "raw"
        layer_option = next(option for option in source.options if option.key == "layer")
        assert [item.id for item in layer_option.inputs] in [["layer_name"], ["source_layer"], ["counts_layer"]]
        for option in source.options:
            if option.key != "layer":
                assert option.inputs == []


def test_expression_source_defaults_follow_the_project_source_policy():
    count_layer_nodes = {
        "OpenBioSingleCellCalculateQC",
        "OpenBioSingleCellCellTypistAnnotation",
        "OpenBioSingleCellAugur",
        "OpenBioSingleCellCNMFRankSurvey",
        "OpenBioSingleCellFilterCells",
        "OpenBioSingleCellFilterGenes",
        "OpenBioSingleCellNormalizeTotal",
        "OpenBioSingleCellNormalizeToLayer",
        "OpenBioSingleCellPearsonResidualsToLayer",
        "OpenBioSingleCellPseudobulk",
        "OpenBioSingleCellQCPlots",
        "OpenBioSingleCellScrublet",
        "OpenBioSingleCellSCVIIntegration",
        "OpenBioSingleCellSnapshotExpression",
    }
    default_layer_nodes = {
        "OpenBioSingleCellCNMFRankSurvey",
        "OpenBioSingleCellHighlyVariableGenes",
        "OpenBioSingleCellPCA",
        "OpenBioSingleCellCellCycleScore",
        "OpenBioSingleCellInferCNV",
        "OpenBioSingleCellCollecTRIULM",
        "OpenBioSingleCellMarkerGenes",
        "OpenBioSingleCellMarkerExpressionPlot",
        "OpenBioSingleCellPearsonResidualsToLayer",
        "OpenBioSingleCellPseudobulk",
        "OpenBioSingleCellSCVIIntegration",
        "OpenBioSingleCellScale",
    }

    source_nodes = {
        node_class.GET_SCHEMA().node_id: node_class
        for node_class in NODE_CLASSES
        if any(item.id == "source" for item in node_class.GET_SCHEMA().inputs)
    }
    assert len(source_nodes) == 27

    for node_id, node_class in source_nodes.items():
        spec = node_class.EXPRESSION_SOURCE
        assert spec.default == ("layer" if node_id in default_layer_nodes else "X")
        assert spec.layer_default == ("counts" if node_id in count_layer_nodes else "log1p_norm")
