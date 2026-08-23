from __future__ import annotations

import importlib
import importlib.metadata
import sys


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: check_imports.py <plugin-root> <comfy-root>")

    sys.path[:0] = sys.argv[1:3]
    importlib.import_module("anndata")
    importlib.import_module("scanpy")
    importlib.import_module("openbio_singlecell.extension")

    print(
        "openbio-singlecell imports OK "
        f"(scanpy {importlib.metadata.version('scanpy')}, "
        f"anndata {importlib.metadata.version('anndata')})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
