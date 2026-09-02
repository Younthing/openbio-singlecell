import assert from "node:assert/strict";
import test from "node:test";

import {
    createOutputReceipt,
    createSingleCellPreview,
    normalizeSingleCellPayload,
    updateOutputReceipt,
    updateSingleCellPreview,
} from "../web/single_cell_result_renderer.mjs";

function element(tag) {
    const listeners = new Map();
    return {
        tag,
        className: "",
        textContent: "",
        children: [],
        attributes: {},
        appendChild(child) {
            this.children.push(child);
            return child;
        },
        replaceChildren(...children) {
            this.children = children;
        },
        setAttribute(name, value) {
            this.attributes[name] = String(value);
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

function findClass(root, className) {
    return descendants(root).find((item) => item.className === className);
}

function previewNode() {
    const widgets = [];
    return {
        widgets,
        size: [300, 220],
        addDOMWidget(...args) {
            widgets.push(args);
        },
        setSize(size) {
            this.size = size;
        },
        setDirtyCanvas() {},
    };
}

function withDocument(run) {
    const originalDocument = globalThis.document;
    globalThis.document = { createElement: element };
    try {
        return run();
    } finally {
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
    }
}

test("accepts each discriminated result kind", () => {
    for (const kind of ["summary", "table", "plot"]) {
        const payload = normalizeSingleCellPayload({
            schema_version: 1,
            kind,
            title: "result",
            columns: ["gene"],
            rows: [["CD3D"]],
            total_rows: 1,
            summary: { cells: 10 },
            warnings: [],
            elapsed_seconds: 0.25,
        });

        assert.equal(payload.kind, kind);
        assert.equal(payload.schema_version, 1);
    }
});

test("rejects unknown schemas and kinds", () => {
    assert.equal(normalizeSingleCellPayload(null), null);
    assert.equal(normalizeSingleCellPayload({ schema_version: 2, kind: "summary" }), null);
    assert.equal(normalizeSingleCellPayload({ schema_version: 1, kind: "unknown" }), null);
});

test("normalizes untrusted collection fields", () => {
    const payload = normalizeSingleCellPayload({
        schema_version: 1,
        kind: "table",
        title: null,
        columns: "gene",
        rows: ["not-a-row", ["CD3D", null]],
        total_rows: -1,
        warnings: "warning",
        elapsed_seconds: Number.NaN,
    });

    assert.deepEqual(payload.columns, []);
    assert.deepEqual(payload.rows, [["CD3D", null]]);
    assert.equal(payload.total_rows, 0);
    assert.deepEqual(payload.warnings, []);
    assert.equal(payload.elapsed_seconds, null);
});

test("summary preview leads with conclusions and collapses supporting report sections", () => {
    withDocument(() => {
        const node = previewNode();
        const updated = updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "summary",
            title: "Differential expression",
            columns: [],
            rows: [],
            total_rows: 0,
            warnings: [],
            elapsed_seconds: 1.25,
            summary: {
                schema_version: 1,
                node_id: "OpenBioSingleCellTest",
                methods: "Compared biological replicates.",
                results: "Detected one supported contrast.",
                key_results: { supported_contrasts: 1 },
                parameters: { alpha: 0.05 },
                warnings: [],
                limitations: ["Small cohort."],
                references: [
                    { citation: "Method paper.", url: "https://example.org", kind: "method" },
                ],
                software_versions: { python: "3.12", "openbio-singlecell": "0.2.0" },
            },
        });

        const root = node.widgets[0][2];
        assert.equal(updated, true);
        assert.equal(
            findClass(root, "openbio-sc-preview__results").textContent,
            "Detected one supported contrast.",
        );
        assert.equal(findClass(root, "openbio-sc-preview__key-results")?.tag, "section");
        assert.deepEqual(
            descendants(root)
                .filter((item) => item.tag === "summary")
                .map((item) => item.textContent),
            ["Methods", "Parameters", "Limitations", "References", "Software versions"],
        );
    });
});

test("summary preview deduplicates report warnings", () => {
    withDocument(() => {
        const node = previewNode();
        updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "summary",
            title: "Result",
            columns: [],
            rows: [],
            total_rows: 0,
            warnings: ["Shared warning.", "Preview warning."],
            elapsed_seconds: 0.1,
            summary: {
                methods: "Method.",
                results: "Result.",
                key_results: {},
                parameters: {},
                warnings: ["Shared warning.", "Report warning."],
                limitations: [],
                references: [],
                software_versions: {},
            },
        });

        const warningList = findClass(node.widgets[0][2], "openbio-sc-preview__warning-list");
        assert.deepEqual(
            warningList.children.map((item) => item.textContent),
            ["Shared warning.", "Preview warning.", "Report warning."],
        );
    });
});

test("summary preview only makes HTTP references clickable", () => {
    withDocument(() => {
        const node = previewNode();
        updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "summary",
            title: "Result",
            columns: [],
            rows: [],
            total_rows: 0,
            warnings: [],
            elapsed_seconds: 0.1,
            summary: {
                methods: "Method.",
                results: "Result.",
                key_results: {},
                parameters: {},
                warnings: [],
                limitations: [],
                references: [
                    { citation: "Safe.", url: "https://example.org/paper", kind: "method" },
                    { citation: "Unsafe.", url: "javascript:alert(1)", kind: "method" },
                ],
                software_versions: {},
            },
        });

        const links = descendants(node.widgets[0][2]).filter((item) => item.tag === "a");
        assert.deepEqual(links.map((item) => item.href), ["https://example.org/paper"]);
        assert.equal(links[0].rel, "noopener noreferrer");
        assert.ok(descendants(node.widgets[0][2]).some((item) => item.textContent === "javascript:alert(1)"));
    });
});

test("summary preview copies the complete report JSON", async () => {
    const originalNavigator = globalThis.navigator;
    const copied = [];
    Object.defineProperty(globalThis, "navigator", {
        configurable: true,
        value: { clipboard: { async writeText(value) { copied.push(value); } } },
    });
    try {
        await withDocument(() => {
            const node = previewNode();
            const summary = {
                methods: "Method.",
                results: "Result.",
                key_results: { cells: 42 },
                parameters: {},
                warnings: [],
                limitations: [],
                references: [],
                software_versions: {},
            };
            updateSingleCellPreview(node, {
                schema_version: 1,
                kind: "summary",
                title: "Result",
                columns: [],
                rows: [],
                total_rows: 0,
                warnings: [],
                elapsed_seconds: 0.1,
                summary,
            });

            const button = findClass(node.widgets[0][2], "openbio-sc-preview__copy");
            return Promise.resolve(button.click()).then(() => {
                assert.deepEqual(JSON.parse(copied[0]), summary);
                assert.equal(button.textContent, "Copied");
            });
        });
    } finally {
        if (originalNavigator === undefined) delete globalThis.navigator;
        else Object.defineProperty(globalThis, "navigator", { configurable: true, value: originalNavigator });
    }
});

test("table preview shows row identity and an honest truncation message", () => {
    withDocument(() => {
        const node = previewNode();
        updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "table",
            title: "Markers",
            columns: ["gene", "effect"],
            rows: [["CD3D", 1.25], [null, 0]],
            index_name: "cell_id",
            index_values: ["cell-a", "cell-b"],
            total_rows: 101,
            summary: { description: "Ranked markers." },
            warnings: [],
            elapsed_seconds: 0.5,
        });

        const table = findClass(node.widgets[0][2], "openbio-sc-preview__table");
        const header = table.children[0].children[0];
        const firstRow = table.children[1].children[0];
        assert.deepEqual(header.children.map((cell) => cell.textContent), ["cell_id", "gene", "effect"]);
        assert.equal(firstRow.children[0].tag, "th");
        assert.equal(firstRow.children[0].textContent, "cell-a");
        assert.equal(firstRow.children[0].attributes.scope, "row");
        assert.equal(findClass(node.widgets[0][2], "openbio-sc-preview__hint").textContent,
            "Showing first 2 of 101 rows. Use Export CSV for the complete table.");
    });
});

test("table preview copies visible rows as TSV and preserves missing values", async () => {
    const originalNavigator = globalThis.navigator;
    const copied = [];
    Object.defineProperty(globalThis, "navigator", {
        configurable: true,
        value: { clipboard: { async writeText(value) { copied.push(value); } } },
    });
    try {
        await withDocument(() => {
            const node = previewNode();
            updateSingleCellPreview(node, {
                schema_version: 1,
                kind: "table",
                title: "Markers",
                columns: ["gene", "effect"],
                rows: [["CD3D", 1.25], [null, 0]],
                index_name: "cell_id",
                index_values: ["cell-a", "cell-b"],
                total_rows: 2,
                summary: { description: "Ranked markers." },
                warnings: [],
                elapsed_seconds: 0.5,
            });

            const root = node.widgets[0][2];
            const table = findClass(root, "openbio-sc-preview__table");
            assert.equal(table.children[1].children[1].children[1].textContent, "—");
            const button = descendants(root).find((item) => item.textContent === "Copy preview TSV");
            return Promise.resolve(button.click()).then(() => {
                assert.equal(copied[0], "cell_id\tgene\teffect\ncell-a\tCD3D\t1.25\ncell-b\tNA\t0");
            });
        });
    } finally {
        if (originalNavigator === undefined) delete globalThis.navigator;
        else Object.defineProperty(globalThis, "navigator", { configurable: true, value: originalNavigator });
    }
});

test("table preview has an explicit empty state", () => {
    withDocument(() => {
        const node = previewNode();
        updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "table",
            title: "No markers",
            columns: ["gene"],
            rows: [],
            index_name: "cell_id",
            index_values: [],
            total_rows: 0,
            summary: {
                description: "No rows met the threshold.",
                parameters: {},
                input_cells: 12,
                input_genes: 100,
                random_seed: 0,
                source: {},
            },
            warnings: [],
            elapsed_seconds: 0.2,
        });

        assert.equal(
            findClass(node.widgets[0][2], "openbio-sc-preview__empty-table").textContent,
            "No rows to preview.",
        );
    });
});

test("plot preview explains the native image with its run metadata", () => {
    withDocument(() => {
        const node = previewNode();
        updateSingleCellPreview(node, {
            schema_version: 1,
            kind: "plot",
            title: "UMAP",
            columns: [],
            rows: [],
            total_rows: 0,
            index_name: "",
            index_values: [],
            summary: {
                description: "Colored by curated cell type.",
                parameters: { color: "cell_type" },
                input_cells: 120,
                input_genes: 500,
                random_seed: 0,
                source: { embedding: "X_umap" },
            },
            warnings: [],
            elapsed_seconds: 0.75,
        });

        const root = node.widgets[0][2];
        assert.equal(
            findClass(root, "openbio-sc-preview__description").textContent,
            "Colored by curated cell type.",
        );
        assert.equal(findClass(root, "openbio-sc-preview__input-size").textContent, "120 cells · 500 genes");
        assert.deepEqual(
            descendants(root)
                .filter((item) => item.tag === "summary")
                .map((item) => item.textContent),
            ["Parameters", "Source"],
        );
        assert.equal(descendants(root).some((item) => item.tag === "img"), false);
    });
});

test("saved file receipt provides a download and copies its relative output path", async () => {
    const originalNavigator = globalThis.navigator;
    const copied = [];
    Object.defineProperty(globalThis, "navigator", {
        configurable: true,
        value: { clipboard: { async writeText(value) { copied.push(value); } } },
    });
    try {
        await withDocument(() => {
            const node = previewNode();
            createOutputReceipt(node);
            const updated = updateOutputReceipt(
                node,
                "OpenBioSingleCellExportCSV",
                {
                    files: [{
                        filename: "markers 00001.csv",
                        subfolder: "openbio-singlecell/run 1",
                        type: "output",
                    }],
                },
                (route) => `/api${route}`,
            );

            const root = node.widgets[0][2];
            const link = descendants(root).find((item) => item.tag === "a");
            assert.equal(updated, true);
            assert.equal(findClass(root, "openbio-sc-output__status").textContent, "Last saved");
            assert.equal(
                findClass(root, "openbio-sc-output__path").textContent,
                "openbio-singlecell/run 1/markers 00001.csv",
            );
            assert.equal(
                link.href,
                "/api/view?filename=markers+00001.csv&subfolder=openbio-singlecell%2Frun+1&type=output",
            );
            assert.equal(link.download, "markers 00001.csv");
            const copy = descendants(root).find((item) => item.textContent === "Copy path");
            return Promise.resolve(copy.click()).then(() => {
                assert.deepEqual(copied, ["openbio-singlecell/run 1/markers 00001.csv"]);
            });
        });
    } finally {
        if (originalNavigator === undefined) delete globalThis.navigator;
        else Object.defineProperty(globalThis, "navigator", { configurable: true, value: originalNavigator });
    }
});

test("preview updates are announced without stealing focus", () => {
    withDocument(() => {
        const root = createSingleCellPreview(previewNode());

        assert.equal(root.attributes["aria-live"], "polite");
        assert.equal(root.attributes.role, "region");
    });
});

test("preview state is idempotent without modifying the node instance", () => {
    const originalDocument = globalThis.document;
    const widgets = [];
    const node = {
        addDOMWidget(...args) {
            widgets.push(args);
        },
    };
    const initialKeys = Reflect.ownKeys(node);
    globalThis.document = {
        createElement(tag) {
            return {
                tag,
                className: "",
                textContent: "",
                children: [],
                appendChild(child) {
                    this.children.push(child);
                },
                setAttribute() {},
            };
        },
    };

    try {
        const first = createSingleCellPreview(node);
        const second = createSingleCellPreview(node);

        assert.equal(second, first);
        assert.equal(widgets.length, 1);
        assert.deepEqual(Reflect.ownKeys(node), initialKeys);
    } finally {
        if (originalDocument === undefined) delete globalThis.document;
        else globalThis.document = originalDocument;
    }
});
