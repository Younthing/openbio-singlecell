from __future__ import annotations

from comfy_api.latest import io

from .study_inputs import DEFAULT_CORE_STUDY_PARAMETERS_JSON, CoreStudyParameters

CATEGORY = "openbio/single-cell/study"


class OpenBioSingleCellCoreStudyParameters(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCoreStudyParameters",
            display_name="Core Study Parameters",
            category=CATEGORY,
            description="Validate core study settings and expose each value as a reusable string.",
            inputs=[
                io.String.Input(
                    "study_parameters_json",
                    default=DEFAULT_CORE_STUDY_PARAMETERS_JSON,
                    multiline=True,
                    socketless=True,
                    extra_dict={"widgetType": "OPENBIO_CORE_STUDY_PARAMETERS_WIDGET"},
                )
            ],
            outputs=[
                io.String.Output("sample_column"),
                io.String.Output("condition_column"),
                io.String.Output("batch_column"),
                io.String.Output("annotation_column"),
                io.String.Output("reference"),
                io.String.Output("comparison"),
            ],
        )

    @classmethod
    def validate_inputs(cls, study_parameters_json: str) -> bool | str:
        try:
            CoreStudyParameters.from_json(study_parameters_json)
        except (TypeError, ValueError) as exc:
            return str(exc)
        return True

    @classmethod
    def execute(
        cls,
        study_parameters_json: str = DEFAULT_CORE_STUDY_PARAMETERS_JSON,
    ) -> io.NodeOutput:
        parameters = CoreStudyParameters.from_json(study_parameters_json)
        return io.NodeOutput(
            parameters.sample_column,
            parameters.condition_column,
            parameters.batch_column,
            parameters.annotation_column,
            parameters.reference,
            parameters.comparison,
        )


STUDY_NODE_CLASSES = (OpenBioSingleCellCoreStudyParameters,)


__all__ = [
    "OpenBioSingleCellCoreStudyParameters",
    "STUDY_NODE_CLASSES",
]
