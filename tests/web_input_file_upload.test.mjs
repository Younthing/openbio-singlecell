import assert from "node:assert/strict";
import test from "node:test";

import {
    acceptsInputFile,
    createInputFileUploadWidgetFactory,
    normalizeAllowedExtensions,
    uploadInputFile,
} from "../web/input_file_upload_widget.mjs";

class FakeClassList {
    values = new Set();

    add(value) {
        this.values.add(value);
    }

    remove(value) {
        this.values.delete(value);
    }

    contains(value) {
        return this.values.has(value);
    }
}

class FakeElement {
    constructor(tag) {
        this.tag = tag;
        this.type = "";
        this.className = "";
        this.textContent = "";
        this.value = "";
        this.children = [];
        this.listeners = {};
        this.attributes = {};
        this.classList = new FakeClassList();
        this.disabled = false;
        this.files = [];
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
        this.attributes[name] = value;
    }

    click() {
        this.clicked = true;
    }
}

function withFakeDocument(callback) {
    const originalDocument = globalThis.document;
    globalThis.document = { createElement: (tag) => new FakeElement(tag) };
    return Promise.resolve()
        .then(callback)
        .finally(() => {
            if (originalDocument === undefined) delete globalThis.document;
            else globalThis.document = originalDocument;
        });
}

function fakeNode() {
    return {
        changes: [],
        addDOMWidget(name, type, element, options) {
            const widget = { name, type, element, options };
            this.widget = widget;
            return widget;
        },
        onWidgetChanged(...args) {
            this.changes.push(args);
        },
        setDirtyCanvas() {},
        graph: { setDirtyCanvas() {} },
    };
}

test("file acceptance is extension based for binary single-cell files", () => {
    assert.deepEqual(normalizeAllowedExtensions(["H5AD", ".H5"]), [".h5ad", ".h5"]);
    assert.equal(acceptsInputFile({ name: "study.H5AD" }, [".h5ad"]), true);
    assert.equal(acceptsInputFile({ name: "study.csv" }, [".h5ad"]), false);
});

test("upload uses the ComfyUI input endpoint and returns its relative path", async () => {
    const file = new File([new Uint8Array([1, 2, 3])], "study.h5ad");
    let request;
    const path = await uploadInputFile(file, {
        subfolder: "openbio-singlecell",
        fetchApi: async (...args) => {
            request = args;
            return {
                ok: true,
                async json() {
                    return { name: "study (1).h5ad", subfolder: "openbio-singlecell", type: "input" };
                },
            };
        },
    });

    assert.equal(request[0], "/upload/image");
    assert.equal(request[1].method, "POST");
    assert.equal(request[1].body.get("type"), "input");
    assert.equal(request[1].body.get("subfolder"), "openbio-singlecell");
    assert.equal(request[1].body.get("image").name, "study.h5ad");
    assert.equal(path, "openbio-singlecell/study (1).h5ad");
});

test("widget supports manual paths and uploads an H5AD dropped on the node", async () => {
    await withFakeDocument(async () => {
        const uploads = [];
        const node = fakeNode();
        const createWidget = createInputFileUploadWidgetFactory({
            fetchApi: async (_route, request) => {
                uploads.push(request.body.get("image").name);
                return {
                    ok: true,
                    async json() {
                        return {
                            name: "dropped.h5ad",
                            subfolder: "openbio-singlecell",
                            type: "input",
                        };
                    },
                };
            },
        });
        const { widget } = createWidget(node, "path", [
            "STRING",
            {
                default: "openbio-singlecell/demo.h5ad",
                allowed_extensions: [".h5ad"],
                accept: ".h5ad",
                upload_subfolder: "openbio-singlecell",
            },
        ]);
        const root = widget.element;
        const pathInput = root.children[0].children[0];
        const status = root.children[2];

        assert.equal(widget.options.hideInPanel, true);
        pathInput.value = "openbio-singlecell/existing.h5ad";
        pathInput.listeners.input[0]();
        assert.equal(widget.serializeValue(), "openbio-singlecell/existing.h5ad");

        const file = new File([new Uint8Array([4, 5, 6])], "dropped.h5ad");
        assert.equal(node.onDragOver({ dataTransfer: { files: [file], items: [] } }), true);
        assert.equal(await node.onDragDrop({ dataTransfer: { files: [file], items: [] } }), true);
        assert.deepEqual(uploads, ["dropped.h5ad"]);
        assert.equal(widget.serializeValue(), "openbio-singlecell/dropped.h5ad");
        assert.equal(status.textContent, "Uploaded: openbio-singlecell/dropped.h5ad");
        assert.equal(node.changes.at(-1)[1], "openbio-singlecell/dropped.h5ad");
    });
});
