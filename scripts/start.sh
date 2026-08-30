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
. "$SCRIPT_DIR/resolve_runtime.sh"

COMFY_ROOT=$(resolve_openbio_comfy_root "$COMFY_ROOT_ARG" "$PLUGIN_ROOT")

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

PYTHON_EXE=$(resolve_openbio_python "$COMFY_ROOT")

"$PYTHON_EXE" -c 'import sys; raise SystemExit(sys.version_info < (3, 12))'
printf 'Using Python: %s\n' "$PYTHON_EXE"
cd "$COMFY_ROOT"

CACHE_MODE_SET=false
for ARG in "$@"; do
    case "$ARG" in
        --cache-classic|--cache-none|--cache-lru|--cache-lru=*|--cache-ram|--cache-ram=*|--high-ram) CACHE_MODE_SET=true ;;
    esac
done
if [ "$CACHE_MODE_SET" = false ]; then
    set -- "$@" --cache-classic
fi

exec "$PYTHON_EXE" main.py --disable-api-nodes --enable-assets --front-end-root "$FRONTEND_DIST" "$@"
