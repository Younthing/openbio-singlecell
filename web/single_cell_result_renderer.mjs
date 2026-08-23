const RESULT_KINDS = new Set(["summary", "table", "plot"]);
const PREVIEW_PROPERTY = "__openbioSingleCellPreview";
const WIDGET_NAME = "openbio_singlecell_preview";

function asText(value) {
    return value == null ? "" : String(value);
}

function textElement(tag, text, className) {
    const element = document.createElement(tag);
    element.className = className;
    element.textContent = asText(text);
    return element;
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
            ? value.rows.filter(Array.isArray).map((row) => row.map(asText))
            : [],
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

function renderTable(root, payload) {
    const wrapper = document.createElement("div");
    wrapper.className = "openbio-sc-preview__table-wrap";
    const table = document.createElement("table");
    table.className = "openbio-sc-preview__table";
    const head = document.createElement("thead");
    const headerRow = document.createElement("tr");

    for (const column of payload.columns) {
        headerRow.appendChild(textElement("th", column, "openbio-sc-preview__header-cell"));
    }
    head.appendChild(headerRow);
    table.appendChild(head);

    const body = document.createElement("tbody");
    for (const row of payload.rows) {
        const tableRow = document.createElement("tr");
        for (const value of row) {
            tableRow.appendChild(textElement("td", value, "openbio-sc-preview__body-cell"));
        }
        body.appendChild(tableRow);
    }
    table.appendChild(body);
    wrapper.appendChild(table);
    root.appendChild(wrapper);
    root.appendChild(
        textElement(
            "p",
            `${payload.total_rows} total rows. Use Export CSV for the complete table.`,
            "openbio-sc-preview__hint",
        ),
    );
}

export function createSingleCellPreview(node) {
    if (node[PREVIEW_PROPERTY]) return node[PREVIEW_PROPERTY];

    const root = document.createElement("section");
    root.className = "openbio-sc-preview";
    root.appendChild(
        textElement("p", "Run this node to preview the result.", "openbio-sc-preview__empty"),
    );
    node.addDOMWidget(WIDGET_NAME, "openbio-single-cell", root, {
        serialize: false,
        hideOnZoom: false,
    });
    node[PREVIEW_PROPERTY] = root;
    return root;
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

    if (payload.kind === "table") renderTable(root, payload);
    else renderSummary(root, payload.summary);
    renderWarnings(root, payload.warnings);

    node.setSize([
        Math.max(node.size?.[0] ?? 0, 360),
        Math.max(node.size?.[1] ?? 0, 240),
    ]);
    node.setDirtyCanvas?.(true, true);
    return true;
}
