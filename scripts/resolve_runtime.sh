resolve_openbio_comfy_root() {
    requested=$1
    plugin_root=$2

    if [ -n "$requested" ]; then
        comfy_root=$requested
        comfy_root_source=--comfy-root
    elif [ -n "${OPENBIO_COMFYUI_ROOT:-}" ]; then
        comfy_root=$OPENBIO_COMFYUI_ROOT
        comfy_root_source=OPENBIO_COMFYUI_ROOT
    else
        manager_comfy_root=$(dirname -- "$(dirname -- "$plugin_root")")
        sibling_comfy_root=$(dirname -- "$plugin_root")/ComfyUI
        if [ "$(basename -- "$(dirname -- "$plugin_root")")" = custom_nodes ] && [ -f "$manager_comfy_root/main.py" ]; then
            comfy_root=$manager_comfy_root
        elif [ -f "$sibling_comfy_root/main.py" ]; then
            comfy_root=$sibling_comfy_root
        else
            printf '%s\n' "ComfyUI root was not found. Pass --comfy-root or set OPENBIO_COMFYUI_ROOT." >&2
            return 1
        fi
        comfy_root_source=
    fi

    if [ ! -d "$comfy_root" ] || [ ! -f "$comfy_root/main.py" ]; then
        printf '%s\n' "${comfy_root_source:-Resolved path} does not point to a ComfyUI root containing main.py: $comfy_root" >&2
        return 1
    fi
    CDPATH= cd -- "$comfy_root" && pwd
}

resolve_openbio_python() {
    comfy_root=$1

    if [ -n "${OPENBIO_PYTHON:-}" ]; then
        printf '%s\n' "$OPENBIO_PYTHON"
    elif [ -x "$comfy_root/.venv/bin/python" ]; then
        printf '%s\n' "$comfy_root/.venv/bin/python"
    elif [ -x "$comfy_root/venv/bin/python" ]; then
        printf '%s\n' "$comfy_root/venv/bin/python"
    elif command -v python3 >/dev/null 2>&1; then
        command -v python3
    else
        printf '%s\n' "Python was not found. Install Python 3.12+ or set OPENBIO_PYTHON." >&2
        return 1
    fi
}
