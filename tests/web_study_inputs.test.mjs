import assert from "node:assert/strict";
import test from "node:test";

import {
    canonicalCoreStudyParametersJSON,
    createCoreStudyParametersEditor,
    normalizeCoreStudyParameters,
    coreStudyParameterWidgets,
} from "../web/study_input_widgets.mjs";

class FakeElement {
    constructor(tag) {
        this.tag = tag;
        this.type = "";
        this.className = "";
        this.textContent = "";
        this.value = "";
        this.children = [];
        this.dataset = {};
        this.listeners = {};
    }

    appendChild(child) {
        this.children.push(child);
        return child;
    }

    addEventListener(name, listener) {
        this.listeners[name] ??= [];
        this.listeners[name].push(listener);
    }

    setAttribute(name, value) {
        this[name] = value;
    }
}

function withFakeDocument(callback) {
    const originalDocument = globalThis.document;
    globalThis.document = { createElement: (tag) => new FakeElement(tag) };
    try {
        return callback();
    } finally {
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
    }
}

function fakeNode() {
    return {
        widgets: [],
        addDOMWidget(name, type, element, options) {
            const widget = { name, type, element, options };
            Object.defineProperty(widget, "value", {
                get: () => options.getValue(),
                set: (value) => options.setValue(value),
            });
            this.widgets.push(widget);
            return widget;
        },
    };
}

test("core study parameters normalize the six reusable string outputs", () => {
    assert.deepEqual(
        normalizeCoreStudyParameters({
            sample_column: " sample_id ",
            condition_column: " treatment ",
            batch_column: " batch ",
            annotation_column: " cluster ",
            reference: " control ",
            comparison: " treated ",
            ignored: "not serialized",
        }),
        {
            schema_version: 1,
            sample_column: "sample_id",
            condition_column: "treatment",
            batch_column: "batch",
            annotation_column: "cluster",
            reference: "control",
            comparison: "treated",
        },
    );

    assert.deepEqual(normalizeCoreStudyParameters(null), {
        schema_version: 1,
        sample_column: "sample",
        condition_column: "group",
        batch_column: "",
        annotation_column: "cell_type",
        reference: "",
        comparison: "",
    });
});

test("core study parameter editor uses six independent text inputs", () => {
    withFakeDocument(() => {
        const changes = [];
        const editor = createCoreStudyParametersEditor({
            value: { reference: "normal", comparison: "tumor" },
            onChange: (value) => changes.push(value),
        });
        const rows = editor.element.children[0].children[0].children;
        const controls = rows.map((row) => row.children[1].children[0]);

        assert.equal(controls.length, 6);
        assert.ok(controls.every((control) => control.tag === "input" && control.type === "text"));
        assert.deepEqual(
            controls.map((control) => control.dataset.field),
            [
                "sample_column",
                "condition_column",
                "batch_column",
                "annotation_column",
                "reference",
                "comparison",
            ],
        );

        controls[4].value = "";
        controls[4].listeners.input[0]();
        controls[5].value = "";
        controls[5].listeners.input[0]();
        assert.equal(editor.getValue().reference, "");
        assert.equal(editor.getValue().comparison, "");
        assert.equal(changes.length, 2);
    });
});

test("clearing required fields stays visible in the serialized value", () => {
    withFakeDocument(() => {
        const editor = createCoreStudyParametersEditor();
        const rows = editor.element.children[0].children[0].children;
        const controls = rows.map((row) => row.children[1].children[0]);

        for (const index of [0, 1, 3]) {
            controls[index].value = "";
            controls[index].listeners.input[0]();
        }

        assert.deepEqual(editor.getValue(), {
            schema_version: 1,
            sample_column: "",
            condition_column: "",
            batch_column: "",
            annotation_column: "",
            reference: "",
            comparison: "",
        });
        assert.equal(
            canonicalCoreStudyParametersJSON(editor.getValue()),
            '{"schema_version":1,"sample_column":"","condition_column":"","batch_column":"","annotation_column":"","reference":"","comparison":""}',
        );
    });
});

test("the only study custom widget serializes canonical core parameters", () => {
    withFakeDocument(() => {
        assert.deepEqual(Object.keys(coreStudyParameterWidgets), [
            "OPENBIO_CORE_STUDY_PARAMETERS_WIDGET",
        ]);

        const result = coreStudyParameterWidgets.OPENBIO_CORE_STUDY_PARAMETERS_WIDGET(
            fakeNode(),
            "study_parameters_json",
        );
        const { widget } = result;
        widget.value = JSON.stringify({ reference: " normal ", comparison: " tumor " });

        assert.equal(
            widget.serializeValue(),
            canonicalCoreStudyParametersJSON({ reference: "normal", comparison: "tumor" }),
        );
        assert.equal(widget.name, "study_parameters_json");
        assert.equal(widget.type, "openbio-core-study-parameters");
        assert.equal(widget.options.socketless, true);
        assert.equal(widget.options.hideInPanel, true);
        assert.equal(widget.options.hideOnZoom, undefined);
        assert.equal(result.minWidth, 420);
    });
});

test("new study widgets initialize from the backend input default", () => {
    withFakeDocument(() => {
        const backendDefault = canonicalCoreStudyParametersJSON({
            sample_column: "donor",
            condition_column: "treatment",
            batch_column: "run",
            annotation_column: "cluster",
            reference: "vehicle",
            comparison: "drug",
        });
        const { widget } = coreStudyParameterWidgets.OPENBIO_CORE_STUDY_PARAMETERS_WIDGET(
            fakeNode(),
            "study_parameters_json",
            ["STRING", { default: backendDefault }],
        );

        assert.equal(widget.serializeValue(), backendDefault);
    });
});
