const RESULT_KINDS = new Set(["summary", "table", "plot"]);
const WIDGET_NAME = "openbio_singlecell_preview";
const OUTPUT_WIDGET_NAME = "openbio_singlecell_output";
const previewRoots = new WeakMap();
const outputRoots = new WeakMap();
const FILE_OUTPUTS = new Map([
    ["OpenBioSingleCellSaveH5AD", "files"],
    ["OpenBioSingleCellExportCSV", "files"],
    ["OpenBioSingleCellSavePNG", "images"],
]);
const PERSIST_OUTPUT = "OpenBioSingleCellPersistArtifact";

function asText(value) {
    return value == null ? "" : String(value);
}

function asCell(value) {
    return value == null ? null : String(value);
}

function textElement(tag, text, className) {
    const element = document.createElement(tag);
    element.className = className;
    element.textContent = asText(text);
    return element;
}

function copyButton(label, value, className = "openbio-sc-preview__copy") {
    const button = textElement("button", label, className);
    button.type = "button";
    button.addEventListener("click", async () => {
        try {
            await globalThis.navigator.clipboard.writeText(value);
            button.textContent = "Copied";
        } catch {
            button.textContent = "Copy failed";
        }
    });
    return button;
}

export function normalizeSingleCellPayload(value) {
    if (!value || typeof value !== "object") return null;
    if (value.schema_version !== 1 || !RESULT_KINDS.has(value.kind)) return null;

    return {
        schema_version: value.schema_version,
        kind: value.kind,
        title: asText(value.title) || "Single-cell result",
        columns: Array.isArray(value.columns) ? value.columns.map(asText) : [],
        rows: Array.isArray(value.rows)
            ? value.rows.filter(Array.isArray).map((row) => row.map(asCell))
            : [],
        index_name: asText(value.index_name) || "Row",
        index_values: Array.isArray(value.index_values) ? value.index_values.map(asText) : [],
        total_rows: Number.isInteger(value.total_rows) && value.total_rows >= 0 ? value.total_rows : 0,
        summary: value.summary ?? {},
        warnings: Array.isArray(value.warnings) ? value.warnings.map(asText) : [],
        elapsed_seconds: Number.isFinite(value.elapsed_seconds) ? value.elapsed_seconds : null,
    };
}

function renderWarnings(root, warnings) {
    if (warnings.length === 0) return;

    const section = document.createElement("section");
    section.className = "openbio-sc-preview__warnings";
    section.appendChild(textElement("strong", "Warnings", "openbio-sc-preview__section-title"));
    const list = document.createElement("ul");
    list.className = "openbio-sc-preview__warning-list";
    for (const warning of warnings) list.appendChild(textElement("li", warning, ""));
    section.appendChild(list);
    root.appendChild(section);
}

function renderSummary(root, summary) {
    const pre = document.createElement("pre");
    pre.className = "openbio-sc-preview__summary";
    pre.textContent = typeof summary === "string" ? summary : JSON.stringify(summary, null, 2);
    root.appendChild(pre);
}

function renderStructuredValue(value, className = "") {
    if (Array.isArray(value)) {
        const list = document.createElement("ul");
        list.className = className;
        for (const item of value) {
            const entry = document.createElement("li");
            entry.appendChild(renderStructuredValue(item));
            list.appendChild(entry);
        }
        return list;
    }
    if (value && typeof value === "object") {
        const list = document.createElement("dl");
        list.className = className;
        for (const [key, item] of Object.entries(value)) {
            list.appendChild(textElement("dt", key, "openbio-sc-preview__term"));
            const description = document.createElement("dd");
            description.className = "openbio-sc-preview__definition";
            description.appendChild(renderStructuredValue(item));
            list.appendChild(description);
        }
        return list;
    }
    return textElement("span", value, className);
}

function renderDetails(root, title, value, renderer = renderStructuredValue) {
    const details = document.createElement("details");
    details.className = "openbio-sc-preview__details";
    details.appendChild(textElement("summary", title, "openbio-sc-preview__details-title"));
    details.appendChild(renderer(value, "openbio-sc-preview__details-body"));
    root.appendChild(details);
}

function safeHttpUrl(value) {
    try {
        const url = new URL(asText(value));
        return url.protocol === "http:" || url.protocol === "https:" ? url.href : null;
    } catch {
        return null;
    }
}

function renderReferences(references, className) {
    const list = document.createElement("ul");
    list.className = className;
    for (const reference of Array.isArray(references) ? references : []) {
        const entry = document.createElement("li");
        entry.appendChild(textElement("span", reference?.citation, ""));
        const url = safeHttpUrl(reference?.url);
        if (url) {
            const link = textElement("a", "Open reference", "openbio-sc-preview__link");
            link.href = url;
            link.target = "_blank";
            link.rel = "noopener noreferrer";
            entry.appendChild(link);
        } else if (reference?.url) {
            entry.appendChild(textElement("span", reference.url, "openbio-sc-preview__reference-url"));
        }
        list.appendChild(entry);
    }
    return list;
}

function renderReport(root, summary) {
    const results = textElement("p", summary.results, "openbio-sc-preview__results");
    root.appendChild(results);

    const keyResults = document.createElement("section");
    keyResults.className = "openbio-sc-preview__key-results";
    keyResults.appendChild(textElement("strong", "Key results", "openbio-sc-preview__section-title"));
    keyResults.appendChild(renderStructuredValue(summary.key_results));
    root.appendChild(keyResults);
    root.appendChild(copyButton("Copy report JSON", JSON.stringify(summary, null, 2)));

    renderDetails(root, "Methods", summary.methods);
    renderDetails(root, "Parameters", summary.parameters);
    renderDetails(root, "Limitations", summary.limitations);
    renderDetails(root, "References", summary.references, renderReferences);
    renderDetails(root, "Software versions", summary.software_versions);
}

function renderArtifactDetails(root, summary) {
    root.appendChild(textElement("p", summary.description, "openbio-sc-preview__description"));
    root.appendChild(
        textElement(
            "p",
            `${asText(summary.input_cells)} cells · ${asText(summary.input_genes)} genes`,
            "openbio-sc-preview__input-size",
        ),
    );
    root.appendChild(
        textElement("p", `Random seed: ${asText(summary.random_seed)}`, "openbio-sc-preview__seed"),
    );
    renderDetails(root, "Parameters", summary.parameters);
    renderDetails(root, "Source", summary.source);
}

function renderTable(root, payload) {
    if (payload.rows.length === 0) {
        root.appendChild(
            textElement("p", "No rows to preview.", "openbio-sc-preview__empty-table"),
        );
        return;
    }

    const wrapper = document.createElement("div");
    wrapper.className = "openbio-sc-preview__table-wrap";
    const table = document.createElement("table");
    table.className = "openbio-sc-preview__table";
    const head = document.createElement("thead");
    const headerRow = document.createElement("tr");

    const indexHeader = textElement("th", payload.index_name, "openbio-sc-preview__index-header");
    indexHeader.setAttribute("scope", "col");
    headerRow.appendChild(indexHeader);

    for (const column of payload.columns) {
        const header = textElement("th", column, "openbio-sc-preview__header-cell");
        header.setAttribute("scope", "col");
        headerRow.appendChild(header);
    }
    head.appendChild(headerRow);
    table.appendChild(head);

    const body = document.createElement("tbody");
    for (const [rowIndex, row] of payload.rows.entries()) {
        const tableRow = document.createElement("tr");
        const index = textElement(
            "th",
            payload.index_values[rowIndex],
            "openbio-sc-preview__index-cell",
        );
        index.setAttribute("scope", "row");
        tableRow.appendChild(index);
        for (const value of row) {
            tableRow.appendChild(
                textElement(
                    "td",
                    value == null ? "—" : value,
                    value == null
                        ? "openbio-sc-preview__body-cell openbio-sc-preview__body-cell--missing"
                        : "openbio-sc-preview__body-cell",
                ),
            );
        }
        body.appendChild(tableRow);
    }
    table.appendChild(body);
    wrapper.appendChild(table);
    root.appendChild(wrapper);
    const tsvValue = (value) => {
        const text = value == null ? "NA" : asText(value);
        return /[\t\r\n"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
    };
    const lines = [
        [payload.index_name, ...payload.columns],
        ...payload.rows.map((row, index) => [payload.index_values[index], ...row]),
    ];
    root.appendChild(
        copyButton(
            "Copy preview TSV",
            lines.map((row) => row.map(tsvValue).join("\t")).join("\n"),
        ),
    );
    const message = payload.total_rows > payload.rows.length
        ? `Showing first ${payload.rows.length} of ${payload.total_rows} rows. Use Export CSV for the complete table.`
        : `${payload.total_rows} rows.`;
    root.appendChild(textElement("p", message, "openbio-sc-preview__hint"));
}

export function createSingleCellPreview(node) {
    const existing = previewRoots.get(node);
    if (existing) return existing;

    const root = document.createElement("section");
    root.className = "openbio-sc-preview";
    root.setAttribute("role", "region");
    root.setAttribute("aria-live", "polite");
    root.appendChild(
        textElement("p", "Run this node to preview the result.", "openbio-sc-preview__empty"),
    );
    node.addDOMWidget(WIDGET_NAME, "openbio-single-cell", root, {
        serialize: false,
        hideOnZoom: false,
    });
    previewRoots.set(node, root);
    return root;
}

export function createOutputReceipt(node) {
    const existing = outputRoots.get(node);
    if (existing) return existing;

    const root = document.createElement("section");
    root.className = "openbio-sc-output";
    root.setAttribute("aria-live", "polite");
    root.appendChild(textElement("p", "Run this node to save the result.", "openbio-sc-output__empty"));
    node.addDOMWidget(OUTPUT_WIDGET_NAME, "openbio-single-cell-output", root, {
        serialize: false,
        hideOnZoom: false,
    });
    outputRoots.set(node, root);
    return root;
}

function outputReceipt(nodeType, output, apiUrlBuilder) {
    if (nodeType === PERSIST_OUTPUT) {
        const path = Array.isArray(output?.text) ? output.text[0] : null;
        return typeof path === "string" && path.trim()
            ? { status: "Last persisted", path, format: "Persisted artifact", download: null }
            : null;
    }

    const outputKey = FILE_OUTPUTS.get(nodeType);
    const file = outputKey && Array.isArray(output?.[outputKey]) ? output[outputKey][0] : null;
    if (
        !file
        || typeof file.filename !== "string"
        || typeof file.subfolder !== "string"
        || typeof file.type !== "string"
    ) return null;

    const path = [file.subfolder.replace(/^\/+|\/+$/g, ""), file.filename]
        .filter(Boolean)
        .join("/");
    const parameters = new URLSearchParams({
        filename: file.filename,
        subfolder: file.subfolder,
        type: file.type,
    });
    const extension = file.filename.includes(".") ? file.filename.split(".").pop().toUpperCase() : "File";
    return {
        status: "Last saved",
        path,
        format: `${extension} · ${file.type}`,
        download: {
            filename: file.filename,
            url: apiUrlBuilder(`/view?${parameters}`),
        },
    };
}

export function updateOutputReceipt(node, nodeType, output, apiUrlBuilder = (route) => route) {
    const receipt = outputReceipt(nodeType, output, apiUrlBuilder);
    if (!receipt) return false;

    const root = createOutputReceipt(node);
    root.replaceChildren();
    root.appendChild(textElement("strong", receipt.status, "openbio-sc-output__status"));
    root.appendChild(textElement("p", receipt.format, "openbio-sc-output__meta"));
    root.appendChild(textElement("code", receipt.path, "openbio-sc-output__path"));

    const actions = document.createElement("div");
    actions.className = "openbio-sc-output__actions";
    if (receipt.download) {
        const link = textElement("a", "Download", "openbio-sc-output__action");
        link.href = receipt.download.url;
        link.download = receipt.download.filename;
        actions.appendChild(link);
    }
    actions.appendChild(copyButton("Copy path", receipt.path, "openbio-sc-output__action"));
    root.appendChild(actions);

    node.setSize([
        Math.max(node.size?.[0] ?? 0, 360),
        Math.max(node.size?.[1] ?? 0, 170),
    ]);
    node.setDirtyCanvas?.(true, true);
    return true;
}

export function updateSingleCellPreview(node, value) {
    const payload = normalizeSingleCellPayload(value);
    if (!payload) return false;

    const root = createSingleCellPreview(node);
    root.replaceChildren();
    root.appendChild(textElement("strong", payload.title, "openbio-sc-preview__title"));
    const elapsed = payload.elapsed_seconds == null ? "—" : `${payload.elapsed_seconds.toFixed(3)}s`;
    root.appendChild(
        textElement("p", `${payload.kind} · ${elapsed}`, "openbio-sc-preview__meta"),
    );

    if (payload.kind === "table") {
        renderArtifactDetails(root, payload.summary);
        renderTable(root, payload);
    }
    else if (payload.kind === "summary" && typeof payload.summary?.results === "string") {
        renderReport(root, payload.summary);
    } else if (payload.kind === "plot" && typeof payload.summary?.description === "string") {
        renderArtifactDetails(root, payload.summary);
    } else renderSummary(root, payload.summary);
    const reportWarnings = Array.isArray(payload.summary?.warnings)
        ? payload.summary.warnings.map(asText)
        : [];
    renderWarnings(root, [...new Set([...payload.warnings, ...reportWarnings])]);

    node.setSize([
        Math.max(node.size?.[0] ?? 0, 360),
        Math.max(node.size?.[1] ?? 0, 240),
    ]);
    node.setDirtyCanvas?.(true, true);
    return true;
}
