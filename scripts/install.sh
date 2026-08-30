#!/usr/bin/env sh
set -eu

usage() {
    cat <<'EOF'
Usage: OPENBIO_PYTHON=/path/to/python sh scripts/install.sh [--comfy-root PATH] [--force]

Installs only openbio-singlecell Python requirements, checks imports, and
generates the local demonstration AnnData H5AD file. Set OPENBIO_PYTHON to the same
interpreter that runs ComfyUI when automatic detection is not appropriate.
ComfyUI is resolved from --comfy-root, OPENBIO_COMFYUI_ROOT, the Manager install
location, or a sibling ComfyUI repository (in that order).
EOF
}

DEMO_ARG=
COMFY_ROOT_ARG=
while [ "$#" -gt 0 ]; do
    case "$1" in
        --force) DEMO_ARG=--force; shift ;;
        --comfy-root)
            if [ "$#" -lt 2 ]; then
                printf '%s\n' "--comfy-root requires a path." >&2
                exit 2
            fi
            COMFY_ROOT_ARG=$2
            shift 2
            ;;
        --comfy-root=*) COMFY_ROOT_ARG=${1#*=}; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PLUGIN_ROOT=$(dirname -- "$SCRIPT_DIR")
. "$SCRIPT_DIR/resolve_runtime.sh"

COMFY_ROOT=$(resolve_openbio_comfy_root "$COMFY_ROOT_ARG" "$PLUGIN_ROOT")
PYTHON_EXE=$(resolve_openbio_python "$COMFY_ROOT")

printf 'Using Python: %s\n' "$PYTHON_EXE"
"$PYTHON_EXE" -c 'import sys; raise SystemExit(sys.version_info < (3, 12))'
if command -v uv >/dev/null 2>&1; then
    uv pip install --python "$PYTHON_EXE" -r "$PLUGIN_ROOT/requirements.txt"
elif "$PYTHON_EXE" -m pip --version >/dev/null 2>&1; then
    "$PYTHON_EXE" -m pip install --disable-pip-version-check -r "$PLUGIN_ROOT/requirements.txt"
else
    printf '%s\n' 'Neither uv nor pip is available. Install uv or add pip to the ComfyUI Python environment.' >&2
    exit 1
fi
"$PYTHON_EXE" "$SCRIPT_DIR/check_imports.py" "$PLUGIN_ROOT" "$COMFY_ROOT"

set -- "$SCRIPT_DIR/generate_demo.py" --comfy-root "$COMFY_ROOT"
if [ -n "$DEMO_ARG" ]; then
    set -- "$@" "$DEMO_ARG"
fi
"$PYTHON_EXE" "$@"
printf '%s\n' 'openbio-singlecell installation complete.'
