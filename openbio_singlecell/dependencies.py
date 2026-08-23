from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

_REQUIREMENTS_PATH = Path(__file__).resolve().parents[1] / "requirements.txt"
_INSTALL_COMMAND = f'python -m pip install -r "{_REQUIREMENTS_PATH}"'


@dataclass(frozen=True, slots=True)
class ScientificDependencies:
    ad: Any
    np: Any
    pd: Any
    sc: Any
    sparse: Any
    Figure: Any
    FigureCanvasAgg: Any
    mmread: Any


@lru_cache(maxsize=1)
def _load_scientific_dependencies() -> ScientificDependencies:
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from scipy import sparse
        from scipy.io import mmread
    except (ImportError, OSError) as error:
        raise RuntimeError(
            "openbio-singlecell scientific dependencies are unavailable "
            f"({error}). Install them with: {_INSTALL_COMMAND}"
        ) from error

    return ScientificDependencies(
        ad=ad,
        np=np,
        pd=pd,
        sc=sc,
        sparse=sparse,
        Figure=Figure,
        FigureCanvasAgg=FigureCanvasAgg,
        mmread=mmread,
    )


def require_scientific_dependencies() -> ScientificDependencies:
    return _load_scientific_dependencies()


__all__ = [
    "ScientificDependencies",
    "require_scientific_dependencies",
]
