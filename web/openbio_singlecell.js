import { app } from "/scripts/app.js";

import {
    createSingleCellPreview,
    updateSingleCellPreview,
} from "./single_cell_result_renderer.mjs";

const PREVIEW_NODE_TYPE = "OpenBioSingleCellPreviewResult";
const STYLESHEET_ID = "openbio-single-cell-preview-styles";

function ensureStylesheet() {
    if (document.getElementById(STYLESHEET_ID)) return;

    const stylesheet = document.createElement("link");
    stylesheet.id = STYLESHEET_ID;
    stylesheet.rel = "stylesheet";
    stylesheet.href = new URL("./openbio_singlecell.css", import.meta.url).href;
    document.head.appendChild(stylesheet);
}

function isPreviewNode(node) {
    return node?.constructor?.comfyClass === PREVIEW_NODE_TYPE || node?.type === PREVIEW_NODE_TYPE;
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

app.registerExtension({
    name: "openbio-singlecell.preview",

    nodeCreated(node) {
        if (!isPreviewNode(node)) return;

        ensureStylesheet();
        createSingleCellPreview(node);
    },

    onNodeOutputsUpdated(nodeOutputs) {
        for (const [locatorId, output] of Object.entries(nodeOutputs)) {
            if (!Array.isArray(output?.openbio_singlecell)) continue;

            const node = getNodeByLocatorId(locatorId);
            if (!isPreviewNode(node)) continue;

            updateSingleCellPreview(node, output.openbio_singlecell[0]);
        }
    },
});
