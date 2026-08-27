export const INPUT_FILE_UPLOAD_WIDGET = "OPENBIO_INPUT_FILE_UPLOAD_WIDGET";

function asOptions(inputData) {
    return inputData?.[1] && typeof inputData[1] === "object" ? inputData[1] : {};
}

export function normalizeAllowedExtensions(values) {
    const source = Array.isArray(values) ? values : String(values ?? "").split(",");
    return source
        .map((value) => String(value).trim().toLowerCase())
        .filter(Boolean)
        .map((value) => (value.startsWith(".") ? value : `.${value}`));
}

export function acceptsInputFile(file, allowedExtensions) {
    if (!file?.name) return false;
    const extensions = normalizeAllowedExtensions(allowedExtensions);
    if (extensions.length === 0) return true;
    const name = String(file.name).toLowerCase();
    return extensions.some((extension) => name.endsWith(extension));
}

function uploadedRelativePath(payload) {
    if (!payload || typeof payload.name !== "string" || payload.name.trim() === "") {
        throw new Error("The upload response did not contain a file name.");
    }
    const name = payload.name.replaceAll("\\", "/").replace(/^\/+/, "");
    const subfolder = String(payload.subfolder ?? "")
        .replaceAll("\\", "/")
        .replace(/^\/+|\/+$/g, "");
    return subfolder ? `${subfolder}/${name}` : name;
}

export async function uploadInputFile(file, { fetchApi, subfolder = "openbio-singlecell" } = {}) {
    if (typeof fetchApi !== "function") throw new TypeError("A ComfyUI upload function is required.");

    const body = new FormData();
    body.append("image", file, file.name);
    body.append("type", "input");
    if (subfolder) body.append("subfolder", subfolder);

    const response = await fetchApi("/upload/image", { method: "POST", body });
    if (!response?.ok) {
        const status = response?.status ? ` (${response.status} ${response.statusText ?? ""})`.trimEnd() : "";
        throw new Error(`File upload failed${status}.`);
    }
    return uploadedRelativePath(await response.json());
}

function droppedFiles(event) {
    const files = Array.from(event?.dataTransfer?.files ?? []);
    if (files.length > 0) return files;
    return Array.from(event?.dataTransfer?.items ?? [])
        .filter((item) => item.kind === "file")
        .map((item) => item.getAsFile?.())
        .filter(Boolean);
}

function transferContainsFile(event, allowedExtensions) {
    const files = droppedFiles(event);
    if (files.length > 0) return files.some((file) => acceptsInputFile(file, allowedExtensions));
    return Array.from(event?.dataTransfer?.items ?? []).some((item) => item.kind === "file");
}

function notifyWidgetChange(node, widget, inputName, previousValue, nextValue) {
    if (previousValue === nextValue) return;
    widget?.callback?.(nextValue);
    node.onWidgetChanged?.(inputName, nextValue, previousValue, widget);
    node.setDirtyCanvas?.(true, true);
    node.graph?.setDirtyCanvas?.(true, true);
}

export function createInputFileUploadWidgetFactory({ fetchApi }) {
    return function createInputFileUploadWidget(node, inputName, inputData) {
        const options = asOptions(inputData);
        const allowedExtensions = normalizeAllowedExtensions(options.allowed_extensions);
        const accept = String(options.accept ?? allowedExtensions.join(","));
        const subfolder = String(options.upload_subfolder ?? "openbio-singlecell");
        const chooseLabel = String(options.upload_label ?? "Choose file");
        const dropLabel = String(
            options.drop_label ??
                `Drop ${allowedExtensions.length ? allowedExtensions.join(" / ") : "a file"} here`,
        );
        let value = String(options.default ?? "");
        let widget;
        let uploading = false;

        const root = document.createElement("section");
        root.className = "openbio-input-file";
        const pathInput = document.createElement("input");
        pathInput.type = "text";
        pathInput.className = "openbio-input-file__path";
        pathInput.setAttribute("aria-label", `${inputName} relative input path`);
        pathInput.value = value;

        const chooseButton = document.createElement("button");
        chooseButton.type = "button";
        chooseButton.className = "openbio-input-file__button";
        chooseButton.textContent = chooseLabel;

        const fileInput = document.createElement("input");
        fileInput.type = "file";
        fileInput.className = "openbio-input-file__native";
        fileInput.accept = accept;
        fileInput.tabIndex = -1;
        fileInput.setAttribute("aria-hidden", "true");

        const dropZone = document.createElement("div");
        dropZone.className = "openbio-input-file__drop-zone";
        dropZone.textContent = dropLabel;

        const status = document.createElement("p");
        status.className = "openbio-input-file__status";
        status.setAttribute("aria-live", "polite");
        status.textContent = value ? `Selected: ${value}` : "No file selected";

        const controls = document.createElement("div");
        controls.className = "openbio-input-file__controls";
        controls.appendChild(pathInput);
        controls.appendChild(chooseButton);
        root.appendChild(controls);
        root.appendChild(dropZone);
        root.appendChild(status);
        root.appendChild(fileInput);

        const updateValue = (nextValue, { notify = true } = {}) => {
            const previous = value;
            value = String(nextValue ?? "");
            pathInput.value = value;
            if (notify) notifyWidgetChange(node, widget, inputName, previous, value);
        };

        const setUploading = (nextUploading, fileName = "") => {
            uploading = nextUploading;
            chooseButton.disabled = nextUploading;
            pathInput.disabled = nextUploading;
            root.setAttribute("aria-busy", String(nextUploading));
            if (nextUploading) status.textContent = `Uploading ${fileName}…`;
        };

        const handleFiles = async (files) => {
            if (uploading) return true;
            const file = Array.from(files ?? []).find((candidate) =>
                acceptsInputFile(candidate, allowedExtensions),
            );
            if (!file) {
                status.textContent = `Choose a supported file: ${allowedExtensions.join(", ")}`;
                return false;
            }

            setUploading(true, file.name);
            try {
                const uploadedPath = await uploadInputFile(file, { fetchApi, subfolder });
                updateValue(uploadedPath);
                status.textContent = `Uploaded: ${uploadedPath}`;
                return true;
            } catch (error) {
                status.textContent = error instanceof Error ? error.message : String(error);
                return true;
            } finally {
                setUploading(false);
            }
        };

        pathInput.addEventListener("input", () => {
            updateValue(pathInput.value);
            status.textContent = value ? `Selected: ${value}` : "No file selected";
        });
        chooseButton.addEventListener("click", () => {
            fileInput.value = "";
            fileInput.click();
        });
        fileInput.addEventListener("change", () => void handleFiles(fileInput.files));

        for (const element of [root, dropZone]) {
            element.addEventListener("dragover", (event) => {
                if (!transferContainsFile(event, allowedExtensions)) return;
                event.preventDefault();
                event.stopPropagation();
                if (event.dataTransfer) event.dataTransfer.dropEffect = "copy";
                root.classList.add("openbio-input-file--dragging");
            });
            element.addEventListener("dragleave", () => {
                root.classList.remove("openbio-input-file--dragging");
            });
            element.addEventListener("drop", (event) => {
                if (!transferContainsFile(event, allowedExtensions)) return;
                event.preventDefault();
                event.stopPropagation();
                root.classList.remove("openbio-input-file--dragging");
                void handleFiles(droppedFiles(event));
            });
        }

        const previousOnDragOver = node.onDragOver;
        const previousOnDragDrop = node.onDragDrop;
        node.onDragOver = function (event) {
            return (
                transferContainsFile(event, allowedExtensions) ||
                previousOnDragOver?.call(this, event) === true
            );
        };
        node.onDragDrop = async function (event) {
            const files = droppedFiles(event);
            if (files.some((file) => acceptsInputFile(file, allowedExtensions))) {
                await handleFiles(files);
                return true;
            }
            return (await previousOnDragDrop?.call(this, event)) ?? false;
        };

        widget = node.addDOMWidget(inputName, "openbio-input-file", root, {
            hideInPanel: true,
            getMinHeight: () => 126,
            getMaxHeight: () => 176,
            getValue: () => value,
            setValue(nextValue) {
                updateValue(nextValue, { notify: false });
                status.textContent = value ? `Selected: ${value}` : "No file selected";
            },
        });
        widget.serializeValue = () => value;
        return { widget, minWidth: 360, minHeight: 154 };
    };
}

export function createInputFileUploadWidgets(options) {
    return Object.freeze({
        [INPUT_FILE_UPLOAD_WIDGET]: createInputFileUploadWidgetFactory(options),
    });
}
