import assert from "node:assert/strict";
import { registerHooks } from "node:module";
import test from "node:test";

let registeredExtension;

globalThis.__openbioTestApp = {
    registerExtension(extension) {
        registeredExtension = extension;
    },
};
globalThis.__openbioTestApi = {
    fetchApi() {
        throw new Error("Unexpected upload request in preview lifecycle test.");
    },
};

registerHooks({
    resolve(specifier, context, nextResolve) {
        if (specifier === "/scripts/app.js") {
            return {
                shortCircuit: true,
                url: "data:text/javascript,export const app=globalThis.__openbioTestApp",
            };
        }
        if (specifier === "/scripts/api.js") {
            return {
                shortCircuit: true,
                url: "data:text/javascript,export const api=globalThis.__openbioTestApi",
            };
        }
        return nextResolve(specifier, context);
    },
});

function element(tag) {
    return {
        tag,
        id: "",
        rel: "",
        href: "",
        className: "",
        textContent: "",
        children: [],
        appendChild(child) {
            this.children.push(child);
            return child;
        },
        replaceChildren(...children) {
            this.children = children;
        },
    };
}

test("live execution replaces the preview placeholder with a summary", async () => {
    const originalDocument = globalThis.document;
    const head = element("head");
    globalThis.document = {
        head,
        createElement: element,
        getElementById(id) {
            return head.children.find((child) => child.id === id) ?? null;
        },
    };

    try {
        await import("../web/openbio_singlecell.js?live-execution-test");

        assert.deepEqual(Object.keys(registeredExtension.getCustomWidgets()), [
            "OPENBIO_CORE_STUDY_PARAMETERS_WIDGET",
            "OPENBIO_INPUT_FILE_UPLOAD_WIDGET",
        ]);
        assert.equal(head.children.length, 1);

        const widgets = [];
        const previousMessages = [];
        const node = {
            constructor: { comfyClass: "OpenBioSingleCellPreviewResult" },
            type: "OpenBioSingleCellPreviewResult",
            size: [300, 220],
            addDOMWidget(...args) {
                widgets.push(args);
            },
            onExecuted(message) {
                previousMessages.push(message);
            },
            setSize(size) {
                this.size = size;
            },
            setDirtyCanvas() {},
        };
        const output = {
            openbio_singlecell: [
                {
                    schema_version: 1,
                    kind: "summary",
                    title: "AnnData summary",
                    columns: [],
                    rows: [],
                    total_rows: 0,
                    summary: { shape: [577, 500] },
                    warnings: [],
                    elapsed_seconds: 0.001,
                },
            ],
        };

        registeredExtension.nodeCreated(node);
        registeredExtension.nodeCreated(node);
        node.onExecuted(output);

        assert.deepEqual(previousMessages, [output]);
        assert.equal(widgets.length, 1);
        const root = widgets[0][2];
        assert.equal(root.children[0].textContent, "AnnData summary");
        assert.match(root.children[2].textContent, /"shape": \[/);
        assert.equal(
            root.children.some(
                (child) => child.textContent === "Run this node to preview the result.",
            ),
            false,
        );
    } finally {
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
        delete globalThis.__openbioTestApp;
        delete globalThis.__openbioTestApi;
    }
});
