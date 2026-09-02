import assert from "node:assert/strict";
import { registerHooks } from "node:module";
import test from "node:test";

let registeredExtension;
const apiPaths = [];

const testApp = {
    registerExtension(extension) {
        registeredExtension = extension;
    },
};
const testApi = {
    apiURL(path) {
        apiPaths.push(path);
        return `api:${path}`;
    },
    fetchApi() {
        throw new Error("Unexpected upload request in preview lifecycle test.");
    },
};
globalThis.__openbioTestApp = testApp;
globalThis.__openbioTestApi = testApi;

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
    const listeners = new Map();
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
        setAttribute(name, value) {
            this[name] = String(value);
        },
        addEventListener(name, listener) {
            listeners.set(name, listener);
        },
        click() {
            return listeners.get("click")?.({ preventDefault() {} });
        },
    };
}

function descendants(root) {
    return [root, ...root.children.flatMap(descendants)];
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

        assert.equal(Object.hasOwn(registeredExtension, "beforeConfigureGraph"), false);
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

        testApp.graph = { getNodeById: () => node };
        registeredExtension.onNodeOutputsUpdated({
            7: {
                openbio_singlecell: [{ ...output.openbio_singlecell[0], title: "Restored summary" }],
            },
        });
        assert.equal(root.children[0].textContent, "Restored summary");
    } finally {
        delete testApp.graph;
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
        delete globalThis.__openbioTestApp;
        delete globalThis.__openbioTestApi;
    }
});

test("output receipts update for live execution and restored history", () => {
    const originalDocument = globalThis.document;
    const head = element("head");
    globalThis.document = {
        head,
        createElement: element,
        getElementById(id) {
            return head.children.find((child) => child.id === id) ?? null;
        },
    };

    const scenarios = [
        {
            type: "OpenBioSingleCellSaveH5AD",
            live: { files: [{ filename: "live.h5ad", subfolder: "run one", type: "output" }] },
            restored: {
                files: [{ filename: "restored.h5ad", subfolder: "run two", type: "output" }],
            },
            expected: "restored.h5ad",
        },
        {
            type: "OpenBioSingleCellExportCSV",
            live: { files: [{ filename: "live.csv", subfolder: "run one", type: "output" }] },
            restored: {
                files: [{ filename: "restored.csv", subfolder: "run two", type: "output" }],
            },
            expected: "restored.csv",
        },
        {
            type: "OpenBioSingleCellSavePNG",
            live: { images: [{ filename: "live.png", subfolder: "run one", type: "output" }] },
            restored: {
                images: [{ filename: "restored.png", subfolder: "run two", type: "output" }],
            },
            expected: "restored.png",
        },
        {
            type: "OpenBioSingleCellPersistArtifact",
            live: { text: ["D:/output/live"] },
            restored: { text: ["D:/output/restored"] },
            expected: "D:/output/restored",
        },
    ];
    const nodes = new Map();

    try {
        for (const [index, scenario] of scenarios.entries()) {
            const widgets = [];
            const previousMessages = [];
            const node = {
                widgets,
                constructor: { comfyClass: scenario.type },
                type: scenario.type,
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

            registeredExtension.nodeCreated(node);
            registeredExtension.nodeCreated(node);
            node.onExecuted(scenario.live);

            assert.deepEqual(previousMessages, [scenario.live]);
            assert.equal(widgets.length, 1);
            assert.match(
                descendants(widgets[0][2])
                    .map((item) => item.textContent)
                    .join("\n"),
                /Last saved|Last persisted/,
            );
            nodes.set(String(index + 1), node);
        }

        testApp.graph = { getNodeById: (id) => nodes.get(String(id)) };
        registeredExtension.onNodeOutputsUpdated(
            Object.fromEntries(
                scenarios.map((scenario, index) => [String(index + 1), scenario.restored]),
            ),
        );

        for (const [index, scenario] of scenarios.entries()) {
            const root = nodes.get(String(index + 1)).widgets[0][2];
            assert.equal(
                descendants(root)
                    .map((item) => item.textContent)
                    .join("\n")
                    .includes(scenario.expected),
                true,
            );
        }
        assert.equal(apiPaths.some((path) => path.startsWith("/view?")), true);
    } finally {
        delete testApp.graph;
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
    }
});
