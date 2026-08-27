export const CORE_STUDY_PARAMETERS_WIDGET = "study_parameters_json";

const CORE_STUDY_PARAMETER_FIELDS = Object.freeze([
    "sample_column",
    "condition_column",
    "batch_column",
    "annotation_column",
    "reference",
    "comparison",
]);

const CORE_STUDY_PARAMETERS_DEFAULTS = Object.freeze({
    schema_version: 1,
    sample_column: "sample",
    condition_column: "group",
    batch_column: "",
    annotation_column: "cell_type",
    reference: "",
    comparison: "",
});

function asText(value) {
    return value == null ? "" : String(value);
}

function parseObject(value) {
    if (value && typeof value === "object") return value;
    if (typeof value !== "string" || value.trim() === "") return null;

    try {
        const parsed = JSON.parse(value);
        return parsed && typeof parsed === "object" ? parsed : null;
    } catch {
        return null;
    }
}

export function normalizeCoreStudyParameters(value) {
    const parsed = parseObject(value) ?? {};
    return {
        schema_version: 1,
        sample_column: asText(parsed.sample_column ?? CORE_STUDY_PARAMETERS_DEFAULTS.sample_column).trim(),
        condition_column: asText(parsed.condition_column ?? CORE_STUDY_PARAMETERS_DEFAULTS.condition_column).trim(),
        batch_column: asText(parsed.batch_column).trim(),
        annotation_column: asText(
            parsed.annotation_column ?? CORE_STUDY_PARAMETERS_DEFAULTS.annotation_column,
        ).trim(),
        reference: asText(parsed.reference).trim(),
        comparison: asText(parsed.comparison).trim(),
    };
}

export function canonicalCoreStudyParametersJSON(value) {
    return JSON.stringify(normalizeCoreStudyParameters(value));
}

function controlRow(labelText, control) {
    const row = document.createElement("tr");
    const label = document.createElement("th");
    label.className = "openbio-study-parameters__label";
    label.textContent = labelText;
    const value = document.createElement("td");
    value.className = "openbio-study-parameters__value";
    value.appendChild(control);
    row.appendChild(label);
    row.appendChild(value);
    return row;
}

function textControl(name) {
    const input = document.createElement("input");
    input.type = "text";
    input.className = "openbio-study-parameters__control";
    input.dataset.field = name;
    input.setAttribute("aria-label", name.replaceAll("_", " "));
    return input;
}

export function createCoreStudyParametersEditor({ value = null, onChange = () => {} } = {}) {
    let parameters = normalizeCoreStudyParameters(value);
    const root = document.createElement("section");
    root.className = "openbio-study-parameters";
    const table = document.createElement("table");
    table.className = "openbio-study-parameters__table";
    const body = document.createElement("tbody");
    const controls = {};
    const labels = {
        sample_column: "Analysis unit column",
        condition_column: "Condition column",
        batch_column: "Batch column",
        annotation_column: "Annotation column",
        reference: "Reference",
        comparison: "Comparison",
    };

    for (const field of CORE_STUDY_PARAMETER_FIELDS) {
        const control = textControl(field);
        control.value = parameters[field];
        control.addEventListener("input", () => {
            const previous = parameters;
            parameters = normalizeCoreStudyParameters({ ...parameters, [field]: control.value });
            onChange({ ...parameters }, previous);
        });
        controls[field] = control;
        body.appendChild(controlRow(labels[field], control));
    }

    table.appendChild(body);
    root.appendChild(table);

    return {
        element: root,
        getValue() {
            return { ...parameters };
        },
        setValue(nextValue, { emit = false } = {}) {
            const previous = parameters;
            parameters = normalizeCoreStudyParameters(nextValue);
            for (const field of CORE_STUDY_PARAMETER_FIELDS) controls[field].value = parameters[field];
            if (emit) onChange({ ...parameters }, previous);
        },
    };
}

function notifyWidgetChange(node, widget, inputName, previousValue, nextValue) {
    if (previousValue === nextValue) return;
    widget?.callback?.(nextValue);
    node.onWidgetChanged?.(inputName, nextValue, previousValue, widget);
    node.setDirtyCanvas?.(true, true);
}

function createCoreStudyParametersWidget(node, inputName = CORE_STUDY_PARAMETERS_WIDGET, inputData) {
    let serialized = canonicalCoreStudyParametersJSON(inputData?.[1]?.default);
    let widget;
    const editor = createCoreStudyParametersEditor({
        value: serialized,
        onChange(parameters) {
            const previous = serialized;
            serialized = canonicalCoreStudyParametersJSON(parameters);
            notifyWidgetChange(node, widget, inputName, previous, serialized);
        },
    });
    widget = node.addDOMWidget(inputName, "openbio-core-study-parameters", editor.element, {
        socketless: true,
        hideInPanel: true,
        getMinHeight: () => 240,
        getMaxHeight: () => 420,
        getValue: () => serialized,
        setValue(nextValue) {
            serialized = canonicalCoreStudyParametersJSON(nextValue);
            editor.setValue(serialized);
        },
    });
    widget.serializeValue = () => serialized;
    return { widget, minWidth: 420, minHeight: 280 };
}

export const coreStudyParameterWidgets = Object.freeze({
    OPENBIO_CORE_STUDY_PARAMETERS_WIDGET: createCoreStudyParametersWidget,
});
