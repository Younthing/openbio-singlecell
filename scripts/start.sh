#!/usr/bin/env sh
set -eu

usage() {
    cat <<'EOF'
Usage: OPENBIO_PYTHON=/path/to/python sh scripts/start.sh [--comfy-root PATH] [ComfyUI arguments...]

Starts ComfyUI with the sibling OpenBio frontend dist, local assets,
and API nodes disabled. Additional arguments are forwarded to ComfyUI.
ComfyUI can also be selected with OPENBIO_COMFYUI_ROOT. The frontend repository
can be selected with OPENBIO_FRONTEND_ROOT.
EOF
}

COMFY_ROOT_ARG=
while [ "$#" -gt 0 ]; do
    case "$1" in
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
        --) shift; break ;;
        *) break ;;
    esac
done

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PLUGIN_ROOT=$(dirname -- "$SCRIPT_DIR")

if [ -n "$COMFY_ROOT_ARG" ]; then
    COMFY_ROOT=$COMFY_ROOT_ARG
    COMFY_ROOT_SOURCE=--comfy-root
elif [ -n "${OPENBIO_COMFYUI_ROOT:-}" ]; then
    COMFY_ROOT=$OPENBIO_COMFYUI_ROOT
    COMFY_ROOT_SOURCE=OPENBIO_COMFYUI_ROOT
else
    MANAGER_COMFY_ROOT=$(dirname -- "$(dirname -- "$PLUGIN_ROOT")")
    SIBLING_COMFY_ROOT=$(dirname -- "$PLUGIN_ROOT")/ComfyUI
    if [ "$(basename -- "$(dirname -- "$PLUGIN_ROOT")")" = custom_nodes ] && [ -f "$MANAGER_COMFY_ROOT/main.py" ]; then
        COMFY_ROOT=$MANAGER_COMFY_ROOT
    elif [ -f "$SIBLING_COMFY_ROOT/main.py" ]; then
        COMFY_ROOT=$SIBLING_COMFY_ROOT
    else
        printf '%s\n' "ComfyUI root was not found. Pass --comfy-root or set OPENBIO_COMFYUI_ROOT." >&2
        exit 1
    fi
    COMFY_ROOT_SOURCE=
fi

if [ ! -d "$COMFY_ROOT" ] || [ ! -f "$COMFY_ROOT/main.py" ]; then
    printf '%s\n' "${COMFY_ROOT_SOURCE:-Resolved path} does not point to a ComfyUI root containing main.py: $COMFY_ROOT" >&2
    exit 1
fi
COMFY_ROOT=$(CDPATH= cd -- "$COMFY_ROOT" && pwd)

if [ -n "${OPENBIO_FRONTEND_ROOT:-}" ]; then
    FRONTEND_ROOT=$OPENBIO_FRONTEND_ROOT
    FRONTEND_ROOT_SOURCE=OPENBIO_FRONTEND_ROOT
else
    SIBLING_FRONTEND_ROOT=$(dirname -- "$PLUGIN_ROOT")/ComfyUI_frontend
    PAIRED_FRONTEND_ROOT=$(dirname -- "$COMFY_ROOT")/ComfyUI_frontend
    if [ -d "$SIBLING_FRONTEND_ROOT" ]; then
        FRONTEND_ROOT=$SIBLING_FRONTEND_ROOT
    elif [ -d "$PAIRED_FRONTEND_ROOT" ]; then
        FRONTEND_ROOT=$PAIRED_FRONTEND_ROOT
    else
        printf '%s\n' "ComfyUI frontend repository was not found. Set OPENBIO_FRONTEND_ROOT." >&2
        exit 1
    fi
    FRONTEND_ROOT_SOURCE=
fi

if [ ! -d "$FRONTEND_ROOT" ]; then
    printf '%s\n' "${FRONTEND_ROOT_SOURCE:-Resolved path} does not point to a frontend repository: $FRONTEND_ROOT" >&2
    exit 1
fi
FRONTEND_ROOT=$(CDPATH= cd -- "$FRONTEND_ROOT" && pwd)
FRONTEND_DIST=$FRONTEND_ROOT/dist

if [ ! -s "$FRONTEND_DIST/index.html" ]; then
    cat >&2 <<EOF
OpenBio frontend dist is incomplete: $FRONTEND_DIST
Missing or empty: index.html
Build it first:
  cd "$FRONTEND_ROOT"
  corepack pnpm install --frozen-lockfile
  corepack pnpm build:openbio
EOF
    exit 1
fi

if [ -n "${OPENBIO_PYTHON:-}" ]; then
    PYTHON_EXE=$OPENBIO_PYTHON
elif [ -x "$COMFY_ROOT/.venv/bin/python" ]; then
    PYTHON_EXE=$COMFY_ROOT/.venv/bin/python
elif [ -x "$COMFY_ROOT/venv/bin/python" ]; then
    PYTHON_EXE=$COMFY_ROOT/venv/bin/python
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_EXE=$(command -v python3)
else
    printf '%s\n' "Python was not found. Install Python 3.12+ or set OPENBIO_PYTHON." >&2
    exit 1
fi

"$PYTHON_EXE" -c 'import sys; raise SystemExit(sys.version_info < (3, 12))'
printf 'Using Python: %s\n' "$PYTHON_EXE"
cd "$COMFY_ROOT"
exec "$PYTHON_EXE" main.py --disable-api-nodes --enable-assets --front-end-root "$FRONTEND_DIST" "$@"
