from __future__ import annotations

from pathlib import Path
from typing import Any

AVAILABLE = False
IMPORT_ERROR: Exception | None = None

ad: Any = None
np: Any = None
pd: Any = None
sc: Any = None
sparse: Any = None
Figure: Any = None
FigureCanvasAgg: Any = None
mmread: Any = None

REQUIREMENTS_PATH = Path(__file__).resolve().parents[1] / "requirements.txt"
INSTALL_COMMAND = f'python -m pip install -r "{REQUIREMENTS_PATH}"'

try:
    import anndata as _ad
    import numpy as _np
    import pandas as _pd
    import scanpy as _sc
    from matplotlib.backends.backend_agg import FigureCanvasAgg as _FigureCanvasAgg
    from matplotlib.figure import Figure as _Figure
    from scipy import sparse as _sparse
    from scipy.io import mmread as _mmread

    ad = _ad
    np = _np
    pd = _pd
    sc = _sc
    sparse = _sparse
    Figure = _Figure
    FigureCanvasAgg = _FigureCanvasAgg
    mmread = _mmread

    AVAILABLE = True
except (ImportError, OSError) as exc:
    IMPORT_ERROR = exc


def require_scientific_dependencies() -> None:
    if AVAILABLE:
        return

    detail = f" ({IMPORT_ERROR})" if IMPORT_ERROR is not None else ""
    raise RuntimeError(
        f"openbio-singlecell scientific dependencies are unavailable{detail}. Install them with: {INSTALL_COMMAND}"
    )
