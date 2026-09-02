import { app } from "/scripts/app.js";
import { api } from "/scripts/api.js";

import {
    createOutputReceipt,
    createSingleCellPreview,
    updateOutputReceipt,
    updateSingleCellPreview,
} from "./single_cell_result_renderer.mjs";
import { createInputFileUploadWidgets } from "./input_file_upload_widget.mjs";
import { coreStudyParameterWidgets } from "./study_input_widgets.mjs";

const PREVIEW_NODE_TYPE = "OpenBioSingleCellPreviewResult";
const OUTPUT_RECEIPT_NODE_TYPES = new Set([
    "OpenBioSingleCellSaveH5AD",
    "OpenBioSingleCellExportCSV",
    "OpenBioSingleCellSavePNG",
    "OpenBioSingleCellPersistArtifact",
]);
const STYLESHEET_ID = "openbio-single-cell-preview-styles";
const liveResultNodes = new WeakSet();
const inputFileUploadWidgets = createInputFileUploadWidgets({ fetchApi: api.fetchApi.bind(api) });
const apiUrlBuilder = api.apiURL.bind(api);

function ensureStylesheet() {
    if (document.getElementById(STYLESHEET_ID)) return;

    const stylesheet = document.createElement("link");
    stylesheet.id = STYLESHEET_ID;
    stylesheet.rel = "stylesheet";
    stylesheet.href = new URL("./openbio_singlecell.css", import.meta.url).href;
    document.head.appendChild(stylesheet);
}

function resultNodeType(node) {
    for (const type of [node?.constructor?.comfyClass, node?.type]) {
        if (type === PREVIEW_NODE_TYPE || OUTPUT_RECEIPT_NODE_TYPES.has(type)) return type;
    }
    return null;
}

function getNodeByLocatorId(locatorId) {
    const identifiers = [locatorId];
    if (/^\d+$/.test(locatorId)) identifiers.push(Number(locatorId));

    for (const graph of [app.rootGraph, app.graph]) {
        for (const identifier of identifiers) {
            const node = graph?.getNodeById?.(identifier);
            if (node) return node;
        }
    }
    return null;
}

function updateResultNode(node, nodeType, output) {
    if (nodeType === PREVIEW_NODE_TYPE) {
        if (Array.isArray(output?.openbio_singlecell)) {
            updateSingleCellPreview(node, output.openbio_singlecell[0]);
        }
        return;
    }
    updateOutputReceipt(node, nodeType, output, apiUrlBuilder);
}

function attachLiveResult(node, nodeType) {
    if (liveResultNodes.has(node)) return;

    const onExecuted = node.onExecuted;
    node.onExecuted = function (output) {
        onExecuted?.apply(this, [output]);
        updateResultNode(this, nodeType, output);
    };
    liveResultNodes.add(node);
}

app.registerExtension({
    name: "openbio-singlecell",

    getCustomWidgets() {
        ensureStylesheet();
        return { ...coreStudyParameterWidgets, ...inputFileUploadWidgets };
    },

    nodeCreated(node) {
        const nodeType = resultNodeType(node);
        if (!nodeType) return;

        ensureStylesheet();
        if (nodeType === PREVIEW_NODE_TYPE) createSingleCellPreview(node);
        else createOutputReceipt(node);
        attachLiveResult(node, nodeType);
    },

    onNodeOutputsUpdated(nodeOutputs) {
        for (const [locatorId, output] of Object.entries(nodeOutputs)) {
            const node = getNodeByLocatorId(locatorId);
            const nodeType = resultNodeType(node);
            if (nodeType) updateResultNode(node, nodeType, output);
        }
    },
});
