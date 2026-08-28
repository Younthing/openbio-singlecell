from __future__ import annotations

import json
from dataclasses import FrozenInstanceError

import pytest

from openbio_singlecell.nodes_study import OpenBioSingleCellCoreStudyParameters
from openbio_singlecell.study_inputs import (
    DEFAULT_CORE_STUDY_PARAMETERS_JSON,
    CoreStudyParameters,
)


def study_parameters_json(**overrides: object) -> str:
    payload = json.loads(DEFAULT_CORE_STUDY_PARAMETERS_JSON)
    payload.update(overrides)
    return json.dumps(payload)


def test_core_study_parameters_parse_plain_strings():
    parameters = CoreStudyParameters.from_json(
        study_parameters_json(
            batch_column="batch",
            reference="normal",
            comparison="tumor",
        )
    )

    assert parameters == CoreStudyParameters(
        sample_column="sample",
        condition_column="group",
        batch_column="batch",
        annotation_column="cell_type",
        reference="normal",
        comparison="tumor",
    )


def test_core_study_parameters_are_immutable():
    parameters = CoreStudyParameters.from_json(DEFAULT_CORE_STUDY_PARAMETERS_JSON)

    with pytest.raises(FrozenInstanceError):
        parameters.sample_column = "replacement"  # type: ignore[misc]


@pytest.mark.parametrize(
    "overrides",
    [
        {"reference": "", "comparison": ""},
        {"reference": "normal", "comparison": ""},
        {"reference": "", "comparison": "tumor"},
        {"reference": "normal", "comparison": "normal"},
        {"sample_column": "", "annotation_column": ""},
        {"sample_column": "group", "condition_column": "group"},
    ],
)
def test_core_study_parameters_do_not_impose_consumer_semantics(overrides: dict[str, str]):
    parameters = CoreStudyParameters.from_json(study_parameters_json(**overrides))

    for field, value in overrides.items():
        assert getattr(parameters, field) == value


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"reference": 3}, "reference must be a string"),
    ],
)
def test_core_study_parameters_reject_invalid_values(overrides: dict[str, object], message: str):
    with pytest.raises(ValueError, match=message):
        CoreStudyParameters.from_json(study_parameters_json(**overrides))


def test_core_study_parameters_preserve_exact_whitespace_owned_by_user():
    parameters = CoreStudyParameters.from_json(
        study_parameters_json(condition_column=" group ", reference=" control ")
    )

    assert parameters.condition_column == " group "
    assert parameters.reference == " control "


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ([], "object at the top level"),
        ({"schema_version": 2}, "missing fields"),
        (
            {
                **json.loads(DEFAULT_CORE_STUDY_PARAMETERS_JSON),
                "schema_version": 2,
            },
            "schema_version must be 1",
        ),
        (
            {
                **json.loads(DEFAULT_CORE_STUDY_PARAMETERS_JSON),
                "unexpected": "value",
            },
            "unexpected fields",
        ),
    ],
)
def test_core_study_parameters_reject_invalid_json_contract(payload: object, message: str):
    with pytest.raises(ValueError, match=message):
        CoreStudyParameters.from_json(json.dumps(payload))


def test_core_study_parameters_node_has_one_local_widget_and_plain_string_outputs():
    schema = OpenBioSingleCellCoreStudyParameters.define_schema()

    assert [input_.id for input_ in schema.inputs] == ["study_parameters_json"]
    assert schema.inputs[0].default == DEFAULT_CORE_STUDY_PARAMETERS_JSON
    assert schema.inputs[0].socketless is True
    assert schema.inputs[0].extra_dict["widgetType"] == "OPENBIO_CORE_STUDY_PARAMETERS_WIDGET"
    assert [output.id for output in schema.outputs] == [
        "sample_column",
        "condition_column",
        "batch_column",
        "annotation_column",
        "reference",
        "comparison",
    ]
    assert [output.io_type for output in schema.outputs] == ["STRING"] * 6


def test_core_study_parameters_node_validates_and_emits_values():
    payload = study_parameters_json(
        batch_column="batch",
        reference="normal",
        comparison="tumor",
    )

    assert OpenBioSingleCellCoreStudyParameters.validate_inputs(payload) is True
    assert OpenBioSingleCellCoreStudyParameters.execute(payload).result == (
        "sample",
        "group",
        "batch",
        "cell_type",
        "normal",
        "tumor",
    )

    invalid_payload = study_parameters_json(reference=3)
    message = OpenBioSingleCellCoreStudyParameters.validate_inputs(invalid_payload)
    assert isinstance(message, str)
    assert "must be a string" in message
    with pytest.raises(ValueError, match="must be a string"):
        OpenBioSingleCellCoreStudyParameters.execute(invalid_payload)
