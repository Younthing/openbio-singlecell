from __future__ import annotations

import random
import threading
import time
import warnings
from collections.abc import Mapping
from importlib.util import find_spec
from textwrap import dedent, indent
from typing import TYPE_CHECKING, Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_plot_result
from .embedding_diagnostic_plotting import (
    _standalone_neighbor_graph_diagnostics_plot,
    _standalone_pca_loadings_plot,
    neighbor_graph_diagnostics_plot_code,
    pca_loadings_plot_code,
)
from .expression_source import _PCA_SPEC, DynamicExpressionSource, ExpressionSource
from .graph_analysis import (
    assess_leiden_stability,
    graph_diagnostics,
    resolve_named_graph,
    run_leiden_partition,
    validate_graph_result_key,
    validate_leiden_settings,
    validate_random_seed,
)
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .pca_plotting import _pca_variance_plot_impl, pca_variance_plot_code
from .worker_protocol import JSONValue, OperationContext, register_operation

if TYPE_CHECKING:
    from anndata import AnnData


SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "scipy", "scikit-learn")
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
_DRAW_GRAPH_RNG_LOCK = threading.RLock()

SCANPY_REFERENCE = AnalysisReference(
    citation="Wolf FA, Angerer P, Theis FJ. SCANPY. Genome Biology. 2018;19:15.",
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
PCA_REFERENCE = AnalysisReference(
    citation="Pearson K. On lines and planes of closest fit. Philosophical Magazine. 1901;2:559-572.",
    doi="10.1080/14786440109462720",
    url="https://doi.org/10.1080/14786440109462720",
    kind="method",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)
UMAP_REFERENCE = AnalysisReference(
    citation="McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection. 2018.",
    url="https://arxiv.org/abs/1802.03426",
    kind="method",
)
GAUSSIAN_KERNEL_REFERENCES = (
    AnalysisReference(
        citation="Coifman RR et al. Geometric diffusions as a tool for harmonic analysis and structure definition of data. PNAS. 2005;102:7426-7431.",
        doi="10.1073/pnas.0500334102",
        url="https://doi.org/10.1073/pnas.0500334102",
        kind="method",
    ),
    AnalysisReference(
        citation="Haghverdi L et al. Diffusion pseudotime robustly reconstructs lineage branching. Nature Methods. 2016;13:845-848.",
        doi="10.1038/nmeth.3971",
        url="https://doi.org/10.1038/nmeth.3971",
        kind="method",
    ),
)
JACCARD_KERNEL_REFERENCE = AnalysisReference(
    citation="Levine JH et al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell. 2015;162:184-197.",
    doi="10.1016/j.cell.2015.05.047",
    url="https://doi.org/10.1016/j.cell.2015.05.047",
    kind="method",
)
TSNE_REFERENCE = AnalysisReference(
    citation="van der Maaten L, Hinton G. Visualizing Data using t-SNE. JMLR. 2008;9:2579-2605.",
    url="https://www.jmlr.org/papers/v9/vandermaaten08a.html",
    kind="method",
)
FR_REFERENCE = AnalysisReference(
    citation="Fruchterman TMJ, Reingold EM. Graph drawing by force-directed placement. SPE. 1991;21:1129-1164.",
    doi="10.1002/spe.4380211102",
    url="https://doi.org/10.1002/spe.4380211102",
    kind="method",
)
KAMADA_KAWAI_REFERENCE = AnalysisReference(
    citation="Kamada T, Kawai S. An algorithm for drawing general undirected graphs. Information Processing Letters. 1989;31:7-15.",
    doi="10.1016/0020-0190(89)90102-6",
    url="https://doi.org/10.1016/0020-0190(89)90102-6",
    kind="method",
)
FORCEATLAS_REFERENCE = AnalysisReference(
    citation="Jacomy M, Venturini T, Heymann S, Bastian M. ForceAtlas2. PLoS ONE. 2014;9:e98679.",
    doi="10.1371/journal.pone.0098679",
    url="https://doi.org/10.1371/journal.pone.0098679",
    kind="method",
)
IGRAPH_REFERENCE = AnalysisReference(
    citation="igraph Python documentation: graph layouts.",
    url="https://python.igraph.org/en/stable/visualisation.html",
    kind="software_documentation",
)
LEIDEN_REFERENCE = AnalysisReference(
    citation="Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden. Scientific Reports. 2019;9:5233.",
    doi="10.1038/s41598-019-41695-z",
    url="https://doi.org/10.1038/s41598-019-41695-z",
    kind="method",
)
ARI_REFERENCE = AnalysisReference(
    citation="Hubert L, Arabie P. Comparing partitions. Journal of Classification. 1985;2:193-218.",
    doi="10.1007/BF01908075",
    url="https://doi.org/10.1007/BF01908075",
    kind="method",
)
SUPPORTED_GRAPH_LAYOUTS = ("fr", "kk", "fa")


def _require_adata(
    adata: AnnData,
    operation: str,
    *,
    min_obs: int = 2,
    require_variables: bool = True,
) -> None:
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
    if int(adata.n_obs) < min_obs:
        raise ValueError(f"{operation} requires at least {min_obs} observations.")
    if require_variables and int(adata.n_vars) < 1:
        raise ValueError(f"{operation} requires at least one variable.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")


def _matrix_values(matrix: Any) -> Any:
    science = dependencies.require_scientific_dependencies()
    return science.np.asarray(matrix.data if science.sparse.issparse(matrix) else matrix)


def _validate_matrix(matrix: Any, *, operation: str, expected_rows: int) -> None:
    science = dependencies.require_scientific_dependencies()
    if len(matrix.shape) != 2 or int(matrix.shape[0]) != expected_rows or int(matrix.shape[1]) < 1:
        raise ValueError(f"{operation} representation must be a two-dimensional observation-aligned matrix.")
    values = _matrix_values(matrix)
    if values.size and (
        not bool(science.np.issubdtype(values.dtype, science.np.number)) or bool(science.np.iscomplexobj(values))
    ):
        raise TypeError(f"{operation} representation must be real numeric.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{operation} representation contains non-finite values.")


def _positive_float(value: Any, name: str, *, allow_zero: bool = False) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a number.")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a number.") from error
    if not bool(science.np.isfinite(value)) or (value < 0 if allow_zero else value <= 0):
        condition = "non-negative" if allow_zero else "greater than zero"
        raise ValueError(f"{name} must be finite and {condition}.")
    return value


def _positive_int(value: Any, name: str, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")
    if value < 0 if allow_zero else value < 1:
        condition = "non-negative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {condition}.")
    return value


def _resolve_representation(adata: AnnData, use_rep: str, n_dimensions: int, operation: str) -> tuple[Any, int]:
    use_rep = use_rep.strip() if isinstance(use_rep, str) else ""
    if not use_rep:
        raise ValueError(f"{operation} use_rep cannot be empty; select X or an explicit obsm key.")
    if use_rep == "X":
        matrix = adata.X
    elif use_rep in adata.obsm:
        matrix = adata.obsm[use_rep]
    else:
        raise ValueError(f"{operation} representation not found in obsm: {use_rep!r}")
    _validate_matrix(matrix, operation=operation, expected_rows=int(adata.n_obs))
    n_dimensions = _positive_int(n_dimensions, f"{operation} n_dimensions", allow_zero=True)
    available = int(matrix.shape[1])
    if n_dimensions > available:
        raise ValueError(f"{operation} requested {n_dimensions} dimensions, but {use_rep!r} has {available}.")
    return matrix, available if n_dimensions == 0 else n_dimensions


def _validate_metric_rows(matrix: Any, *, n_dimensions: int, metric: str, operation: str) -> list[str]:
    science = dependencies.require_scientific_dependencies()
    selected = (
        matrix[:, :n_dimensions] if science.sparse.issparse(matrix) else science.np.asarray(matrix)[:, :n_dimensions]
    )
    if metric == "cosine":
        if science.sparse.issparse(selected):
            squared_norms = science.np.asarray(selected.multiply(selected).sum(axis=1), dtype=float).ravel()
        else:
            squared_norms = science.np.einsum("ij,ij->i", selected, selected, dtype=float)
        invalid = squared_norms == 0
        if bool(invalid.any()):
            return [
                f"{operation} cosine distance uses the backend convention for {int(invalid.sum())} "
                "zero-norm observation row(s); these rows have no expression direction."
            ]
    elif metric == "correlation":
        if science.sparse.issparse(selected):
            row_maximum = selected.max(axis=1)
            row_minimum = selected.min(axis=1)
            maxima = science.np.asarray(
                row_maximum.toarray() if science.sparse.issparse(row_maximum) else row_maximum,
                dtype=float,
            ).ravel()
            minima = science.np.asarray(
                row_minimum.toarray() if science.sparse.issparse(row_minimum) else row_minimum,
                dtype=float,
            ).ravel()
            constant = maxima == minima
        else:
            constant = science.np.ptp(selected, axis=1) == 0
        if bool(constant.any()):
            raise ValueError(
                f"{operation} correlation metric is undefined for {int(constant.sum())} constant observation row(s)."
            )
    return []


def _coordinate_diagnostics(coordinates: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    values = science.np.asarray(coordinates, dtype=float)
    return {
        "shape": list(values.shape),
        "axis_summaries": [summarize_numeric(values[:, axis]) for axis in range(values.shape[1])],
        "duplicate_coordinate_rows": int(values.shape[0] - science.np.unique(values, axis=0).shape[0]),
    }


def _validate_embedding_coordinates(coordinates: Any, *, n_obs: int, operation: str) -> Any:
    science = dependencies.require_scientific_dependencies()
    values = science.np.asarray(coordinates)
    if values.shape != (n_obs, 2):
        raise RuntimeError(f"{operation} backend must return coordinates with shape (n_obs, 2).")
    if not bool(science.np.issubdtype(values.dtype, science.np.number)) or bool(science.np.iscomplexobj(values)):
        raise TypeError(f"{operation} backend coordinates must be real numeric values.")
    if not bool(science.np.isfinite(values).all()):
        raise RuntimeError(f"{operation} backend returned non-finite coordinates.")
    return values


def _validate_pca_bundle(adata: AnnData, *, n_comps: int, mask: Any) -> tuple[Any, Any, Any, Any]:
    science = dependencies.require_scientific_dependencies()

    def floating_array(value: Any, *, shape: tuple[int, ...], name: str) -> Any:
        array = science.np.asarray(value)
        if array.shape != shape:
            raise RuntimeError(f"PCA backend returned invalid {name} shape {array.shape!r}; expected {shape!r}.")
        if not bool(science.np.issubdtype(array.dtype, science.np.floating)):
            raise RuntimeError(f"PCA backend returned non-floating {name} values.")
        if not bool(science.np.isfinite(array).all()):
            raise RuntimeError(f"PCA backend returned non-finite {name} values.")
        return array

    scores = floating_array(
        adata.obsm.get("X_pca"),
        shape=(int(adata.n_obs), n_comps),
        name="scores",
    )
    if scores.dtype != science.np.dtype("float32"):
        raise RuntimeError("PCA backend did not honor the requested float32 score dtype.")
    loadings = floating_array(
        adata.varm.get("PCs"),
        shape=(int(adata.n_vars), n_comps),
        name="loadings",
    )
    metadata = adata.uns.get("pca")
    if not isinstance(metadata, Mapping):
        raise RuntimeError("PCA backend did not return the canonical uns['pca'] mapping.")
    variance = floating_array(metadata.get("variance"), shape=(n_comps,), name="explained variance")
    variance_ratio = floating_array(
        metadata.get("variance_ratio"),
        shape=(n_comps,),
        name="explained variance ratio",
    )
    if bool((variance < 0).any()) or bool((variance_ratio < 0).any()):
        raise RuntimeError("PCA backend returned negative variance diagnostics.")
    if mask is not None and not bool(
        science.np.allclose(loadings[~science.np.asarray(mask, dtype=bool)], 0.0, rtol=0.0, atol=1e-12)
    ):
        raise RuntimeError("PCA backend returned non-zero loadings for excluded variables.")
    return scores, loadings, variance, variance_ratio


def _validate_paga_initialization(adata: AnnData) -> None:
    science = dependencies.require_scientific_dependencies()
    paga = adata.uns.get("paga")
    if not isinstance(paga, Mapping):
        raise ValueError("PAGA initialization requires an uns['paga'] mapping; run PAGA first.")
    groups = paga.get("groups")
    if not isinstance(groups, str) or not groups or groups not in adata.obs:
        raise ValueError("PAGA initialization requires its groups column in obs.")
    labels = adata.obs[groups]
    if not isinstance(labels.dtype, science.pd.CategoricalDtype):
        raise ValueError("PAGA initialization groups must be categorical.")
    if bool(labels.isna().any()) or len(labels.cat.categories) < 1:
        raise ValueError("PAGA initialization groups must be complete and categorical.")
    n_categories = len(labels.cat.categories)
    positions = science.np.asarray(paga.get("pos"))
    if positions.shape != (n_categories, 2):
        raise ValueError("PAGA initialization positions must align to group categories and have two columns.")
    if not bool(science.np.issubdtype(positions.dtype, science.np.number)) or bool(science.np.iscomplexobj(positions)):
        raise TypeError("PAGA initialization positions must be real numeric values.")
    if not bool(science.np.isfinite(positions).all()):
        raise ValueError("PAGA initialization positions must be finite.")
    connectivities = paga.get("connectivities")
    if not science.sparse.issparse(connectivities):
        raise TypeError("PAGA initialization connectivities must be a SciPy sparse matrix.")
    if tuple(connectivities.shape) != (n_categories, n_categories):
        raise ValueError("PAGA initialization connectivities must align to group categories.")
    connectivities = connectivities.tocsr(copy=True)
    connectivities.sum_duplicates()
    connectivities.eliminate_zeros()
    weights = science.np.asarray(connectivities.data)
    if weights.size and (
        not bool(science.np.issubdtype(weights.dtype, science.np.number)) or bool(science.np.iscomplexobj(weights))
    ):
        raise TypeError("PAGA initialization connectivities must contain real numeric weights.")
    if weights.size and (not bool(science.np.isfinite(weights).all()) or bool((weights < 0).any())):
        raise ValueError("PAGA initialization connectivities must contain finite non-negative weights.")
    if not bool(science.np.allclose(connectivities.diagonal(), 0.0, rtol=0.0, atol=1e-12)):
        raise ValueError("PAGA initialization connectivities must have a zero diagonal.")
    delta = connectivities - connectivities.T
    delta.eliminate_zeros()
    if delta.nnz and float(science.np.max(science.np.abs(delta.data))) > 1e-8:
        raise ValueError("PAGA initialization connectivities must be symmetric.")


def _neighbor_method_references(method: str) -> list[AnalysisReference]:
    if method == "umap":
        return [UMAP_REFERENCE, SCANPY_REFERENCE]
    if method == "gauss":
        return [*GAUSSIAN_KERNEL_REFERENCES, SCANPY_REFERENCE]
    return [JACCARD_KERNEL_REFERENCE, SCANPY_REFERENCE]


def _layout_references(layout: str) -> list[AnalysisReference]:
    if layout == "fa":
        return [FORCEATLAS_REFERENCE, SCANPY_REFERENCE]
    if layout == "kk":
        return [KAMADA_KAWAI_REFERENCE, IGRAPH_REFERENCE, SCANPY_REFERENCE]
    return [FR_REFERENCE, IGRAPH_REFERENCE, SCANPY_REFERENCE]


def _temporary_obsm_key(adata: AnnData, prefix: str) -> str:
    key = prefix
    suffix = 0
    while key in adata.obsm:
        suffix += 1
        key = f"{prefix}_{suffix}"
    return key


def _graph_warnings(diagnostics: dict[str, Any]) -> list[str]:
    warnings = []
    if diagnostics["connected_components"] > 1:
        warnings.append("The selected neighbor graph is disconnected; global layout separation is not identifiable.")
    if diagnostics["isolated_observation_count"]:
        warnings.append("The selected neighbor graph contains isolated observations.")
    connectivities = diagnostics["connectivities"]
    if not connectivities["symmetric"]:
        warnings.append("The selected neighbor graph is asymmetric; UMAP consumed the directed weights as supplied.")
    if connectivities["self_loop_count"]:
        warnings.append(
            f"The selected neighbor graph contains {connectivities['self_loop_count']} nonzero self-loop(s)."
        )
    return warnings


def _generated_graph_helper() -> str:
    return dedent(
        """
        from collections.abc import Mapping

        import numpy as np
        from scipy import sparse


        def _openbio_sparse_graph_matrix(
            matrix, n_obs, operation, matrix_name, require_zero_diagonal=True
        ):
            if not sparse.issparse(matrix):
                raise TypeError(f"{operation} {matrix_name} matrix must be a SciPy sparse matrix.")
            if matrix.shape != (n_obs, n_obs):
                raise ValueError(f"{operation} {matrix_name} matrix must have shape (n_obs, n_obs).")
            matrix = matrix.tocsr(copy=True)
            matrix.sum_duplicates()
            matrix.eliminate_zeros()
            values = np.asarray(matrix.data)
            if values.size and (not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values)):
                raise TypeError(f"{operation} {matrix_name} matrix must contain real numeric weights.")
            if values.size and not np.isfinite(values).all():
                raise ValueError(f"{operation} {matrix_name} matrix contains non-finite weights.")
            if values.size and (values < 0).any():
                raise ValueError(f"{operation} {matrix_name} matrix contains negative weights.")
            if require_zero_diagonal and not np.allclose(matrix.diagonal(), 0.0, rtol=0.0, atol=1e-12):
                raise ValueError(f"{operation} {matrix_name} matrix must have a zero diagonal.")
            return matrix


        def _openbio_validate_metric_rows(matrix, n_dimensions, metric, operation):
            selected = matrix[:, :n_dimensions] if sparse.issparse(matrix) else np.asarray(matrix)[:, :n_dimensions]
            if metric == "cosine":
                if sparse.issparse(selected):
                    squared_norms = np.asarray(selected.multiply(selected).sum(axis=1), dtype=float).ravel()
                else:
                    squared_norms = np.einsum("ij,ij->i", selected, selected, dtype=float)
                invalid = squared_norms == 0
                if invalid.any():
                    warnings.warn(
                        f"{operation} cosine distance uses the backend convention for {int(invalid.sum())} "
                        "zero-norm observation row(s); these rows have no expression direction.",
                        UserWarning,
                        stacklevel=2,
                    )
            elif metric == "correlation":
                if sparse.issparse(selected):
                    row_maximum = selected.max(axis=1)
                    row_minimum = selected.min(axis=1)
                    maxima = np.asarray(
                        row_maximum.toarray() if sparse.issparse(row_maximum) else row_maximum,
                        dtype=float,
                    ).ravel()
                    minima = np.asarray(
                        row_minimum.toarray() if sparse.issparse(row_minimum) else row_minimum,
                        dtype=float,
                    ).ravel()
                    constant = maxima == minima
                else:
                    constant = np.ptp(selected, axis=1) == 0
                if constant.any():
                    raise ValueError(
                        f"{operation} correlation metric is undefined for {int(constant.sum())} constant observation row(s)."
                    )


        def _openbio_canonical_graph(
            adata,
            neighbors_key,
            operation,
            require_distances=False,
            require_undirected=True,
            require_zero_diagonal=True,
        ):
            if adata.isbacked or adata.n_obs < 2:
                raise ValueError(f"{operation} requires an in-memory AnnData with at least two observations.")
            if not adata.obs_names.is_unique:
                raise ValueError(f"{operation} requires unique observation identifiers.")
            metadata = adata.uns.get(neighbors_key)
            if not isinstance(metadata, Mapping):
                raise ValueError(f"{operation} named neighbor graph is missing.")
            connectivity_key = metadata.get("connectivities_key")
            if not isinstance(connectivity_key, str) or not connectivity_key.strip():
                raise ValueError(f"{operation} named connectivity key is invalid.")
            connectivity_key = connectivity_key.strip()
            if connectivity_key not in adata.obsp:
                raise ValueError(f"{operation} named connectivity matrix is missing.")
            matrix = _openbio_sparse_graph_matrix(
                adata.obsp[connectivity_key],
                adata.n_obs,
                operation,
                "connectivity",
                require_zero_diagonal=require_zero_diagonal,
            )
            if require_undirected:
                delta = matrix - matrix.T
                delta.eliminate_zeros()
                if delta.nnz and float(np.max(np.abs(delta.data))) > 1e-8:
                    raise ValueError(f"{operation} requires an undirected graph.")
                matrix = ((matrix + matrix.T) * 0.5).tocsr()
            if require_zero_diagonal:
                matrix.setdiag(0)
                matrix.eliminate_zeros()
            off_diagonal = matrix.copy()
            off_diagonal.setdiag(0)
            off_diagonal.eliminate_zeros()
            if off_diagonal.nnz == 0:
                raise ValueError(f"{operation} graph has no positive edges.")
            adata.obsp[connectivity_key] = matrix
            distances_key = None
            if require_distances:
                distances_key = metadata.get("distances_key")
                if not isinstance(distances_key, str) or not distances_key.strip():
                    raise ValueError(f"{operation} named distance key is invalid.")
                distances_key = distances_key.strip()
                if distances_key not in adata.obsp:
                    raise ValueError(f"{operation} named distance matrix is missing.")
                _openbio_sparse_graph_matrix(
                    adata.obsp[distances_key], adata.n_obs, operation, "distance"
                )
            return distances_key, connectivity_key


        """
    )


def _simple_code(function_name: str, body: str) -> str:
    statements = indent(dedent(body).strip(), "    ")
    return (
        "from collections.abc import Mapping\n\n"
        "import numpy as np\n"
        "import scanpy as sc\n\n\n"
        "from scipy import sparse\n\n\n"
        f"def {function_name}(adata):\n"
        "    output = adata\n"
        f"{statements}\n"
        "    return output\n"
    )


def _pca_code(expression: ExpressionSource, parameters: dict[str, Any]) -> str:
    source_kind = expression.kind
    layer_name = expression.layer_name
    return dedent(
        f"""from collections.abc import Mapping

import numpy as np
import scanpy as sc
from scipy import sparse


def run_pca(adata):
    if adata.isbacked or adata.n_obs < 2 or adata.n_vars < 2:
        raise ValueError("PCA requires an in-memory AnnData with at least two observations and variables.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("PCA requires unique observation and variable identifiers.")
    output = adata
    if {source_kind!r} == "layer":
        if {layer_name!r} not in output.layers:
            raise ValueError("PCA source layer is missing.")
        matrix = output.layers[{layer_name!r}]
        layer = {layer_name!r}
    else:
        matrix, layer = output.X, None
    values = np.asarray(matrix.data if sparse.issparse(matrix) else matrix)
    if values.size and (not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values) or not np.isfinite(values).all()):
        raise ValueError("PCA input must contain finite numeric values.")
    if {parameters["use_hvg"]!r}:
        if "highly_variable" not in output.var or str(output.var["highly_variable"].dtype) not in {{"bool", "boolean"}} or output.var["highly_variable"].isna().any():
            raise TypeError("PCA requires a complete boolean highly_variable mask.")
        mask = np.asarray(output.var["highly_variable"], dtype=bool)
        selected_vars = int(mask.sum())
    else:
        mask, selected_vars = None, output.n_vars
    if selected_vars < 2 or not 1 <= {parameters["n_comps"]!r} < min(output.n_obs, selected_vars):
        raise ValueError("PCA with ARPACK requires n_comps < min(n_obs, selected variables).")
    selected = matrix[:, mask] if mask is not None else matrix
    if sparse.issparse(selected):
        means = np.asarray(selected.mean(axis=0), dtype=float).ravel()
        mean_squares = np.asarray(selected.power(2).mean(axis=0), dtype=float).ravel()
        variances = np.maximum(mean_squares - means**2, 0.0)
    else:
        variances = np.var(np.asarray(selected), axis=0)
    if np.all(variances == 0):
        raise ValueError("PCA cannot run because every selected variable is constant.")
    estimate = output.n_obs * {parameters["n_comps"]!r} * 4 + output.n_vars * {parameters["n_comps"]!r} * 8 + {parameters["n_comps"]!r} * 16
    if estimate > {parameters["max_output_gib"]!r} * 1024**3:
        raise ValueError("PCA dense output estimate exceeds max_output_gib.")
    collisions = [name for name, present in (("X_pca", "X_pca" in output.obsm), ("PCs", "PCs" in output.varm), ("pca", "pca" in output.uns)) if present]
    if collisions and not {parameters["overwrite_existing"]!r}:
        raise ValueError(f"PCA output bundle already exists: {{collisions}}")
    for container, key in ((output.obsm, "X_pca"), (output.varm, "PCs"), (output.uns, "pca")):
        if {parameters["overwrite_existing"]!r} and key in container:
            del container[key]
    sc.pp.pca(output, n_comps={parameters["n_comps"]!r}, layer=layer, mask_var={"'highly_variable'" if parameters["use_hvg"] else "None"}, zero_center=True, svd_solver="arpack", random_state={parameters["random_seed"]!r}, dtype="float32", chunked=False)
    scores = np.asarray(output.obsm.get("X_pca"))
    loadings = np.asarray(output.varm.get("PCs"))
    metadata = output.uns.get("pca")
    if scores.shape != (output.n_obs, {parameters["n_comps"]!r}):
        raise RuntimeError("PCA backend returned invalid score shape.")
    if scores.dtype != np.dtype("float32") or not np.isfinite(scores).all():
        raise RuntimeError("PCA backend returned invalid float32 scores.")
    if loadings.shape != (output.n_vars, {parameters["n_comps"]!r}) or not np.issubdtype(loadings.dtype, np.floating) or not np.isfinite(loadings).all():
        raise RuntimeError("PCA backend returned invalid loadings.")
    if not isinstance(metadata, Mapping):
        raise RuntimeError("PCA backend did not return the canonical metadata mapping.")
    variance = np.asarray(metadata.get("variance"))
    variance_ratio = np.asarray(metadata.get("variance_ratio"))
    for name, values in (("variance", variance), ("variance ratio", variance_ratio)):
        if values.shape != ({parameters["n_comps"]!r},) or not np.issubdtype(values.dtype, np.floating) or not np.isfinite(values).all() or (values < 0).any():
            raise RuntimeError(f"PCA backend returned invalid {{name}} values.")
    if mask is not None and not np.allclose(loadings[~mask], 0.0, rtol=0.0, atol=1e-12):
        raise RuntimeError("PCA backend returned non-zero loadings for excluded variables.")
    return output
"""
    )


class OpenBioSingleCellPCA:
    EXPRESSION_SOURCE = _PCA_SPEC

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_comps: int = 50,
        use_hvg: bool = True,
        source: DynamicExpressionSource | None = None,
        overwrite_existing: bool = False,
        max_output_gib: float = 2.0,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_adata(adata, "PCA")
        if int(adata.n_vars) < 2:
            raise ValueError("PCA requires at least two variables on the current axis.")
        if not bool(adata.var_names.is_unique):
            raise ValueError("PCA requires unique variable identifiers.")
        random_seed = validate_random_seed(random_seed)
        if use_hvg and "highly_variable" not in adata.var:
            raise ValueError(
                "PCA with use_hvg enabled requires var['highly_variable']; run Highly Variable Genes first."
            )
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        matrix = expression.matrix(adata)
        _validate_matrix(matrix, operation="PCA", expected_rows=int(adata.n_obs))
        warnings = []
        mask: Any = None
        selected_vars = int(adata.n_vars)
        if use_hvg:
            if (
                str(adata.var["highly_variable"].dtype) not in {"bool", "boolean"}
                or adata.var["highly_variable"].isna().any()
            ):
                raise TypeError("PCA var['highly_variable'] must be a complete boolean mask.")
            mask = science.np.asarray(adata.var["highly_variable"], dtype=bool)
            selected_vars = int(mask.sum())
            if selected_vars < 2:
                raise ValueError("PCA requires at least two selected highly-variable genes.")
        n_comps = _positive_int(n_comps, "PCA n_comps")
        if n_comps >= min(int(adata.n_obs), selected_vars):
            raise ValueError("PCA with the ARPACK solver requires n_comps < min(n_obs, selected variables).")
        selected = matrix[:, mask] if mask is not None else matrix
        if science.sparse.issparse(selected):
            means = science.np.asarray(selected.mean(axis=0), dtype=float).ravel()
            mean_squares = science.np.asarray(selected.power(2).mean(axis=0), dtype=float).ravel()
            variances = science.np.maximum(mean_squares - means**2, 0.0)
        else:
            variances = science.np.var(science.np.asarray(selected), axis=0)
        constant_mask = science.np.asarray(variances == 0, dtype=bool)
        constant_count = int(constant_mask.sum())
        if constant_count == selected_vars:
            raise ValueError("PCA cannot run because every selected variable is constant.")
        selected_indices = science.np.flatnonzero(mask) if mask is not None else science.np.arange(int(adata.n_vars))
        constant_examples = [str(adata.var_names[index]) for index in selected_indices[constant_mask][:20]]
        if constant_count:
            warnings.append(
                f"{constant_count} selected variable(s) had zero variance and therefore cannot contribute to PCA."
            )
        budget = _positive_float(max_output_gib, "PCA max_output_gib") * 1024**3
        estimated = int(adata.n_obs) * n_comps * 4 + int(adata.n_vars) * n_comps * 8 + n_comps * 16
        if estimated > budget:
            raise ValueError(f"PCA dense output estimate ({estimated / 1024**3:.3f} GiB) exceeds max_output_gib.")
        collisions = [
            key
            for key, present in (
                ("obsm['X_pca']", "X_pca" in adata.obsm),
                ("varm['PCs']", "PCs" in adata.varm),
                ("uns['pca']", "pca" in adata.uns),
            )
            if present
        ]
        if collisions and not overwrite_existing:
            raise ValueError(f"PCA output bundle already exists: {', '.join(collisions)}")
        started_at = time.perf_counter()
        output = adata
        for container, key in ((output.obsm, "X_pca"), (output.varm, "PCs"), (output.uns, "pca")):
            if overwrite_existing and key in container:
                del container[key]
        science.sc.pp.pca(
            output,
            n_comps=n_comps,
            layer=expression.scanpy_layer,
            mask_var="highly_variable" if use_hvg else None,
            zero_center=True,
            svd_solver="arpack",
            random_state=random_seed,
            dtype="float32",
            chunked=False,
        )
        scores, loadings, variance, variance_ratio = _validate_pca_bundle(output, n_comps=n_comps, mask=mask)
        parameters = {
            "n_comps": n_comps,
            "use_hvg": use_hvg,
            **expression.parameters(),
            "zero_center": True,
            "svd_solver": "arpack",
            "dtype": "float32",
            "overwrite_existing": overwrite_existing,
            "max_output_gib": max_output_gib,
            "random_seed": random_seed,
        }
        finish_adata(
            output, "pca", parameters, int(adata.n_obs), int(adata.n_vars), started_at, random_seed=random_seed
        )
        top = []
        for component in range(min(5, n_comps)):
            selected_loadings = science.np.abs(loadings[selected_indices, component])
            indices = selected_indices[science.np.argsort(selected_loadings)[-5:][::-1]]
            top.append({"component": component + 1, "genes": [str(output.var_names[index]) for index in indices]})
        code = _pca_code(expression, parameters)
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellPCA",
            title="PCA",
            operation="pca",
            methods=(
                f"Centered PCA used Scanpy's ARPACK path on {selected_vars} variables from the explicitly selected "
                f"{expression.kind} expression source."
            ),
            results=f"Computed {n_comps} components explaining {float(variance_ratio.sum()):.3%} of selected-expression variance.",
            key_results={
                "components": n_comps,
                "selected_variables": selected_vars,
                "constant_selected_variables": constant_count,
                "constant_selected_variable_examples": constant_examples,
                "variance": variance.tolist(),
                "cumulative_variance_ratio": float(variance_ratio.sum()),
                "variance_ratio": variance_ratio.tolist(),
                "top_absolute_loading_genes": top,
                "scores": {
                    **_coordinate_diagnostics(scores),
                    "dtype": str(scores.dtype),
                },
                "loadings": {"shape": list(loadings.shape), "dtype": str(loadings.dtype)},
                "variance_dtype": str(variance.dtype),
                "variance_ratio_dtype": str(variance_ratio.dtype),
                "estimated_dense_output_bytes": estimated,
            },
            parameters=parameters,
            references=[PCA_REFERENCE, SCANPY_REFERENCE],
            software_packages=SOFTWARE_PACKAGES,
            warnings=warnings,
            limitations=["PCA is a linear projection; component sign is arbitrary and loading rank is descriptive."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


def pca_variance_plot_owned(adata: Any, n_pcs: int = 0) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _pca_variance_plot_impl(adata, n_pcs=n_pcs)
    parameters = {"n_pcs": details["requested_n_pcs"]}
    warnings_list = details["warnings"]
    code = pca_variance_plot_code(n_pcs=parameters["n_pcs"])
    plotted = make_plot_result(
        png=png,
        title="PCA variance",
        operation="pca_variance_plot",
        parameters=parameters,
        description=(
            f"Read-only rendering of explained and cumulative variance for "
            f"{details['plotted_components']:,} stored principal components."
        ),
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPCAVariancePlot",
        title="PCA variance",
        operation="pca_variance_plot",
        methods=(
            "Read the canonical PCA variance ratios stored in adata.uns['pca']['variance_ratio'] without "
            "recomputing PCA, then rendered per-component bars and a cumulative curve with Matplotlib Agg."
        ),
        results=(
            f"Displayed {details['plotted_components']:,} of {details['available_components']:,} stored "
            f"components, reaching {details['cumulative_variance_ratio'][-1]:.3%} cumulative explained variance."
        ),
        key_results=details,
        parameters=parameters,
        references=[PCA_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy"],
        warnings=warnings_list,
        limitations=[
            "Explained variance is descriptive of the stored PCA fit and does not determine biological relevance."
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def pca_loadings_plot_owned(adata: Any, component: int = 1, n_genes: int = 20) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _standalone_pca_loadings_plot(
        adata,
        component=component,
        n_genes=n_genes,
        _return_details=True,
    )
    parameters = {"component": details["component"], "n_genes": details["requested_genes"]}
    warnings_list = []
    if details["plotted_genes"] < details["requested_genes"]:
        warnings_list.append(
            f"Requested {details['requested_genes']:,} genes but only {details['plotted_genes']:,} features exist."
        )
    code = pca_loadings_plot_code(**parameters)
    plotted = make_plot_result(
        png=png,
        title=f"PC{details['component']} loadings",
        operation="pca_loadings_plot",
        parameters=parameters,
        description="Read-only signed loading plot from the stored PCA bundle.",
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPCALoadingsPlot",
        title=f"PC{details['component']} loadings",
        operation="pca_loadings_plot",
        methods=(
            "Read the selected column of adata.varm['PCs'], ranked features by descending absolute loading "
            "with stable feature-order ties, and rendered signed bars without recomputing PCA."
        ),
        results=(
            f"Displayed {details['plotted_genes']:,} features for PC{details['component']}; "
            f"{details['positive_loading_count']:,} shown loadings were positive and "
            f"{details['negative_loading_count']:,} were negative."
        ),
        key_results=details,
        parameters=parameters,
        references=[PCA_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy"],
        warnings=warnings_list,
        limitations=[
            "PCA component signs are arbitrary; loading magnitude is descriptive and is not marker evidence."
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


class OpenBioSingleCellNeighbors:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        use_rep: str = "X_pca",
        n_dimensions: int = 0,
        n_neighbors: int = 15,
        metric: str = "euclidean",
        method: str = "umap",
        key_added: str = "neighbors",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_adata(adata, "Neighbors", require_variables=False)
        random_seed = validate_random_seed(random_seed)
        matrix, used = _resolve_representation(adata, use_rep, n_dimensions, "Neighbors")
        n_neighbors = _positive_int(n_neighbors, "Neighbors n_neighbors")
        if n_neighbors < 2:
            raise ValueError("Neighbors n_neighbors must be at least 2.")
        if metric not in {"euclidean", "cosine", "correlation", "manhattan"}:
            raise ValueError(f"Unsupported Neighbors metric: {metric!r}")
        metric_warnings = _validate_metric_rows(matrix, n_dimensions=used, metric=metric, operation="Neighbors")
        if method not in {"umap", "gauss", "jaccard"}:
            raise ValueError(f"Unsupported Neighbors method: {method!r}")
        key_added = key_added.strip()
        if not key_added:
            raise ValueError("Neighbors key_added cannot be empty.")
        distance_key = "distances" if key_added == "neighbors" else f"{key_added}_distances"
        connectivity_key = "connectivities" if key_added == "neighbors" else f"{key_added}_connectivities"
        collisions = [
            name
            for name, present in (
                (f"uns[{key_added!r}]", key_added in adata.uns),
                (f"obsp[{distance_key!r}]", distance_key in adata.obsp),
                (f"obsp[{connectivity_key!r}]", connectivity_key in adata.obsp),
            )
            if present
        ]
        if collisions and not overwrite_existing:
            raise ValueError(f"Neighbors output bundle already exists: {', '.join(collisions)}")
        started_at = time.perf_counter()
        output = adata
        for container, key in ((output.uns, key_added), (output.obsp, distance_key), (output.obsp, connectivity_key)):
            if overwrite_existing and key in container:
                del container[key]
        backend_rep, backend_dims = use_rep, used
        temp_key = _temporary_obsm_key(output, "__openbio_neighbors_representation")
        if use_rep == "X" and used < int(matrix.shape[1]):
            output.obsm[temp_key] = matrix[:, :used].copy()
            backend_rep, backend_dims = temp_key, None
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            science.sc.pp.neighbors(
                output,
                n_neighbors=n_neighbors,
                n_pcs=backend_dims,
                use_rep=backend_rep,
                metric=metric,
                method=method,
                knn=True,
                transformer=None,
                random_state=random_seed,
                key_added=None if key_added == "neighbors" else key_added,
            )
        if temp_key in output.obsm:
            del output.obsm[temp_key]
        graph = resolve_named_graph(output, key_added, operation="Neighbors", require_distances=True)
        if graph.distances_key != distance_key or graph.connectivities_key != connectivity_key:
            raise RuntimeError("Neighbors backend returned output pointers that do not match the requested key bundle.")
        metadata = output.uns.get(key_added)
        params = metadata.get("params") if isinstance(metadata, Mapping) else None
        if not isinstance(params, Mapping):
            raise RuntimeError("Neighbors backend did not return a valid parameter mapping.")
        effective_value = params.get("n_neighbors")
        if isinstance(effective_value, bool) or not isinstance(effective_value, (int, science.np.integer)):
            raise RuntimeError("Neighbors backend did not disclose an integer effective n_neighbors.")
        effective_n_neighbors = int(effective_value)
        if effective_n_neighbors < 2:
            raise RuntimeError("Neighbors backend returned an invalid effective n_neighbors.")
        metadata["params"] = dict(params)
        metadata["params"]["use_rep"] = use_rep
        metadata["params"]["n_pcs"] = used
        parameters = {
            "use_rep": use_rep,
            "n_dimensions": used,
            "n_neighbors": n_neighbors,
            "effective_n_neighbors": effective_n_neighbors,
            "metric": metric,
            "method": method,
            "knn": True,
            "transformer": "auto (Scanpy default selection)",
            "key_added": key_added,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        finish_adata(
            output, "neighbors", parameters, int(adata.n_obs), int(adata.n_vars), started_at, random_seed=random_seed
        )
        diagnostics = graph_diagnostics(graph)
        code = _simple_code(
            "compute_neighbors",
            f"""            if output.isbacked or output.n_obs < 2 or not output.obs_names.is_unique:\n                raise ValueError("Neighbors requires an in-memory AnnData with unique observations.")\n            use_rep = {use_rep!r}\n            if use_rep != "X" and use_rep not in output.obsm:\n                raise ValueError(f"Neighbors representation not found: {{use_rep!r}}")\n            matrix = output.X if use_rep == "X" else output.obsm[use_rep]\n            values = np.asarray(matrix.data if sparse.issparse(matrix) else matrix)\n            if len(matrix.shape) != 2 or matrix.shape[0] != output.n_obs or matrix.shape[1] < 1:\n                raise ValueError("Neighbors representation must be observation aligned.")\n            if values.size and (not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values)):\n                raise TypeError("Neighbors representation must be real numeric.")\n            if values.size and not np.isfinite(values).all():\n                raise ValueError("Neighbors representation contains non-finite values.")\n            if {used!r} > matrix.shape[1] or {n_neighbors!r} < 2:\n                raise ValueError("Neighbors dimensions or neighbor count exceed input bounds.")\n            _openbio_validate_metric_rows(matrix, {used!r}, {metric!r}, "Neighbors")\n            distance_key = {distance_key!r}\n            connectivity_key = {connectivity_key!r}\n            collisions = [{key_added!r} in output.uns, distance_key in output.obsp, connectivity_key in output.obsp]\n            if any(collisions) and not {overwrite_existing!r}:\n                raise ValueError("Neighbors output bundle already exists.")\n            for container, key in ((output.uns, {key_added!r}), (output.obsp, distance_key), (output.obsp, connectivity_key)):\n                if {overwrite_existing!r} and key in container:\n                    del container[key]\n            backend_rep, backend_dims = use_rep, {used!r}\n            temp_key = "__openbio_neighbors_representation"\n            if use_rep == "X" and {used!r} < matrix.shape[1]:\n                output.obsm[temp_key] = matrix[:, :{used!r}].copy()\n                backend_rep, backend_dims = temp_key, None\n            with warnings.catch_warnings(record=True):\n                warnings.simplefilter("always")\n                sc.pp.neighbors(output, n_neighbors={n_neighbors!r}, n_pcs=backend_dims, use_rep=backend_rep, metric={metric!r}, method={method!r}, knn=True, transformer="auto", random_state={random_seed!r}, key_added={None if key_added == "neighbors" else key_added!r})\n            if temp_key in output.obsm:\n                del output.obsm[temp_key]\n            actual_distance_key, actual_connectivity_key = _openbio_canonical_graph(output, {key_added!r}, "Neighbors", require_distances=True)\n            if (actual_distance_key, actual_connectivity_key) != (distance_key, connectivity_key):\n                raise RuntimeError("Neighbors backend returned inconsistent output pointers.")\n            metadata = output.uns.get({key_added!r})\n            params = metadata.get("params") if isinstance(metadata, Mapping) else None\n            if not isinstance(params, Mapping):\n                raise RuntimeError("Neighbors backend did not return a valid parameter mapping.")\n            effective_value = params.get("n_neighbors")\n            if isinstance(effective_value, bool) or not isinstance(effective_value, (int, np.integer)) or int(effective_value) < 2:\n                raise RuntimeError("Neighbors backend did not disclose a valid integer effective n_neighbors.")\n            effective_n_neighbors = int(effective_value)\n            metadata["params"] = dict(params)\n            metadata["params"].update(use_rep=use_rep, n_pcs={used!r})\n            if effective_n_neighbors != {n_neighbors!r}:\n                warnings.warn(f"Requested n_neighbors={{{n_neighbors!r}}} was resolved to effective n_neighbors={{effective_n_neighbors}} by Scanpy.", UserWarning, stacklevel=2)""",
        )
        code = code.replace('transformer="auto"', "transformer=None")
        code = code.replace(
            "from collections.abc import Mapping", "import warnings\nfrom collections.abc import Mapping", 1
        )
        code = code.replace(
            'temp_key = "__openbio_neighbors_representation"',
            'temp_prefix = "__openbio_neighbors_representation"\n'
            "    temp_key, suffix = temp_prefix, 0\n"
            "    while temp_key in output.obsm:\n"
            "        suffix += 1\n"
            '        temp_key = f"{temp_prefix}_{suffix}"',
        )
        code = _generated_graph_helper() + code
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellNeighbors",
            title="Neighbor graph",
            operation="neighbors",
            methods=f"Constructed a neighbor graph from {used} dimensions of {use_rep!r} using requested k={n_neighbors}, effective k={effective_n_neighbors}, {metric} distance, and Scanpy's {method} connectivity method.",
            results=f"The graph contains {graph.positive_edges} undirected edges in {len(graph.component_sizes)} connected component(s).",
            key_results={
                "actual_dimensions": used,
                "metric_row_warnings": metric_warnings,
                "requested_n_neighbors": n_neighbors,
                "effective_n_neighbors": effective_n_neighbors,
                "graph": diagnostics,
            },
            parameters=parameters,
            references=_neighbor_method_references(method),
            software_packages=SOFTWARE_PACKAGES + ("pynndescent", "umap-learn"),
            warnings=[
                *metric_warnings,
                *(
                    [
                        f"Requested n_neighbors={n_neighbors} was resolved to effective n_neighbors={effective_n_neighbors} by Scanpy."
                    ]
                    if effective_n_neighbors != n_neighbors
                    else []
                ),
                *_graph_warnings(diagnostics),
            ],
            limitations=["Neighborhood structure depends on the chosen representation, metric, and k."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


def neighbor_graph_diagnostics_plot_owned(
    adata: Any,
    neighbors_key: str = "neighbors",
    max_working_memory_gib: float = 4.0,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _standalone_neighbor_graph_diagnostics_plot(
        adata,
        neighbors_key=neighbors_key,
        max_working_memory_gib=max_working_memory_gib,
        _return_details=True,
    )
    parameters = {
        "neighbors_key": details["neighbors_key"],
        "max_working_memory_gib": details["memory_preflight"]["max_working_memory_gib"],
    }
    warnings_list = []
    if details["connected_components"] > 1:
        warnings_list.append(
            f"The stored graph has {details['connected_components']:,} connected components."
        )
    code = neighbor_graph_diagnostics_plot_code(**parameters)
    plotted = make_plot_result(
        png=png,
        title=f"Neighbor graph diagnostics: {details['neighbors_key']}",
        operation="neighbor_graph_diagnostics_plot",
        parameters=parameters,
        description="Read-only graph degree, distance, and connected-component diagnostics.",
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellNeighborGraphDiagnosticsPlot",
        title=f"Neighbor graph diagnostics: {details['neighbors_key']}",
        operation="neighbor_graph_diagnostics_plot",
        methods=(
            "Validated the named stored connectivity and distance matrices, then rendered cell degree, stored "
            "positive distances, and connected-component sizes without rebuilding the neighbor graph."
        ),
        results=(
            f"The graph contains {details['observations']:,} observations across "
            f"{details['connected_components']:,} connected component(s)."
        ),
        key_results=details,
        parameters=parameters,
        references=[SCANPY_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy", "scipy"],
        warnings=warnings_list,
        limitations=[
            "Graph diagnostics describe stored geometry; they do not establish biological populations or an optimal neighbor setting."
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


class OpenBioSingleCellUMAP:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        neighbors_key: str = "neighbors",
        min_dist: float = 0.5,
        spread: float = 1.0,
        key_added: str = "X_umap",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_adata(adata, "UMAP", require_variables=False)
        random_seed = validate_random_seed(random_seed)
        graph = resolve_named_graph(
            adata,
            neighbors_key,
            operation="UMAP",
            require_undirected=False,
            require_zero_diagonal=False,
        )
        min_dist = _positive_float(min_dist, "UMAP min_dist", allow_zero=True)
        spread = _positive_float(spread, "UMAP spread")
        key_added = key_added.strip()
        if not key_added:
            raise ValueError("UMAP key_added cannot be empty.")
        uns_key = "umap" if key_added == "X_umap" else key_added
        if (key_added in adata.obsm or uns_key in adata.uns) and not overwrite_existing:
            raise ValueError("UMAP coordinate/parameter output already exists.")
        started_at = time.perf_counter()
        output = adata
        if overwrite_existing:
            output.obsm.pop(key_added, None)
            output.uns.pop(uns_key, None)
        output.obsp[graph.connectivities_key] = graph.connectivities
        science.sc.tl.umap(
            output,
            min_dist=min_dist,
            spread=spread,
            n_components=2,
            maxiter=None,
            alpha=1.0,
            gamma=1.0,
            negative_sample_rate=5,
            init_pos="spectral",
            a=None,
            b=None,
            method="umap",
            neighbors_key=neighbors_key,
            random_state=random_seed,
            key_added=None if key_added == "X_umap" else key_added,
        )
        coordinate_values = _validate_embedding_coordinates(
            output.obsm.get(key_added), n_obs=int(output.n_obs), operation="UMAP"
        )
        umap_metadata = output.uns.get(uns_key)
        if not isinstance(umap_metadata, Mapping) or not isinstance(umap_metadata.get("params"), Mapping):
            raise RuntimeError("UMAP backend did not return the expected parameter record.")
        parameters = {
            "neighbors_key": neighbors_key,
            "min_dist": min_dist,
            "spread": spread,
            "n_components": 2,
            "maxiter": None,
            "alpha": 1.0,
            "gamma": 1.0,
            "negative_sample_rate": 5,
            "init_pos": "spectral",
            "a": None,
            "b": None,
            "method": "umap",
            "key_added": key_added,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        finish_adata(
            output, "umap", parameters, int(adata.n_obs), int(adata.n_vars), started_at, random_seed=random_seed
        )
        coordinates = _coordinate_diagnostics(coordinate_values)
        graph_report = graph_diagnostics(graph)
        code = _simple_code(
            "run_umap",
            f"""            _openbio_canonical_graph(output, {neighbors_key!r}, "UMAP", require_undirected=False, require_zero_diagonal=False)
            if {min_dist!r} > {spread!r}:
                warnings.warn("UMAP min_dist exceeds spread; this expert setting changes the usual compactness interpretation.", UserWarning, stacklevel=2)
            if ({key_added!r} in output.obsm or {uns_key!r} in output.uns) and not {overwrite_existing!r}:
                raise ValueError("UMAP output already exists.")
            if {overwrite_existing!r}:
                output.obsm.pop({key_added!r}, None)
                output.uns.pop({uns_key!r}, None)
            sc.tl.umap(output, min_dist={min_dist!r}, spread={spread!r}, n_components=2, maxiter=None, alpha=1.0, gamma=1.0, negative_sample_rate=5, init_pos="spectral", a=None, b=None, method="umap", neighbors_key={neighbors_key!r}, random_state={random_seed!r}, key_added={None if key_added == "X_umap" else key_added!r})
            coordinates = np.asarray(output.obsm.get({key_added!r}))
            if coordinates.shape != (output.n_obs, 2):
                raise RuntimeError("UMAP backend returned invalid coordinate shape.")
            if not np.issubdtype(coordinates.dtype, np.number) or np.iscomplexobj(coordinates):
                raise TypeError("UMAP backend coordinates must be real numeric values.")
            if not np.isfinite(coordinates).all():
                raise RuntimeError("UMAP backend returned non-finite coordinates.")
            metadata = output.uns.get({uns_key!r})
            if not isinstance(metadata, Mapping) or not isinstance(metadata.get("params"), Mapping):
                raise RuntimeError("UMAP backend did not return the expected parameter record.")""",
        )
        code = code.replace(
            "from collections.abc import Mapping", "import warnings\nfrom collections.abc import Mapping", 1
        )
        code = _generated_graph_helper() + code
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellUMAP",
            title="UMAP",
            operation="umap",
            methods="Computed a two-dimensional CPU UMAP embedding with spectral initialization from the explicitly named neighbor graph.",
            results=f"Generated {key_added!r} coordinates for {adata.n_obs} observations.",
            key_results={"coordinates": coordinates, "graph": graph_report},
            parameters=parameters,
            references=[UMAP_REFERENCE, SCANPY_REFERENCE],
            software_packages=SOFTWARE_PACKAGES + ("umap-learn",),
            warnings=[
                *(
                    ["UMAP min_dist exceeds spread; this expert setting changes the usual compactness interpretation."]
                    if min_dist > spread
                    else []
                ),
                *_graph_warnings(graph_report),
            ],
            limitations=["UMAP coordinates and global distances are visualization aids, not inferential measurements."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


class OpenBioSingleCellTSNE:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        use_rep: str = "X_pca",
        n_dimensions: int = 0,
        perplexity: float = 30.0,
        metric: str = "euclidean",
        early_exaggeration: float = 12.0,
        learning_rate: float = 1000.0,
        key_added: str = "X_tsne",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_adata(adata, "t-SNE", require_variables=False)
        random_seed = validate_random_seed(random_seed)
        matrix, used = _resolve_representation(adata, use_rep, n_dimensions, "t-SNE")
        perplexity = _positive_float(perplexity, "t-SNE perplexity")
        if perplexity >= int(adata.n_obs):
            raise ValueError("t-SNE requires perplexity < n_obs.")
        early_exaggeration = _positive_float(early_exaggeration, "t-SNE early_exaggeration")
        if metric not in {"euclidean", "cosine", "correlation", "manhattan"}:
            raise ValueError(f"Unsupported t-SNE metric: {metric!r}")
        metric_warnings = _validate_metric_rows(matrix, n_dimensions=used, metric=metric, operation="t-SNE")
        learning_rate = _positive_float(learning_rate, "t-SNE learning_rate")
        key_added = key_added.strip()
        if not key_added:
            raise ValueError("t-SNE key_added cannot be empty.")
        if (key_added in adata.obsm or key_added in adata.uns) and not overwrite_existing:
            raise ValueError(f"t-SNE coordinate/parameter output already exists for key: {key_added!r}")
        started_at = time.perf_counter()
        output = adata
        if overwrite_existing:
            output.obsm.pop(key_added, None)
            output.uns.pop(key_added, None)
        backend_rep, backend_dims = use_rep, used
        temp_key = _temporary_obsm_key(output, "__openbio_tsne_representation")
        if use_rep == "X" and used < int(matrix.shape[1]):
            output.obsm[temp_key] = matrix[:, :used].copy()
            backend_rep, backend_dims = temp_key, None
        science.sc.tl.tsne(
            output,
            n_pcs=backend_dims,
            n_components=2,
            use_rep=backend_rep,
            perplexity=perplexity,
            metric=metric,
            early_exaggeration=early_exaggeration,
            learning_rate=learning_rate,
            random_state=random_seed,
            use_fast_tsne=False,
            n_jobs=None,
            key_added=key_added,
            copy=False,
        )
        output.obsm.pop(temp_key, None)
        coordinate_values = _validate_embedding_coordinates(
            output.obsm.get(key_added), n_obs=int(output.n_obs), operation="t-SNE"
        )
        tsne_metadata = output.uns.get(key_added)
        if not isinstance(tsne_metadata, Mapping) or not isinstance(tsne_metadata.get("params"), Mapping):
            raise RuntimeError("t-SNE backend did not return the expected parameter record.")
        parameters = {
            "use_rep": use_rep,
            "n_dimensions": used,
            "perplexity": perplexity,
            "metric": metric,
            "early_exaggeration": early_exaggeration,
            "learning_rate": learning_rate,
            "n_components": 2,
            "use_fast_tsne": False,
            "n_jobs": None,
            "key_added": key_added,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        finish_adata(
            output, "tsne", parameters, int(adata.n_obs), int(adata.n_vars), started_at, random_seed=random_seed
        )
        diagnostics = _coordinate_diagnostics(coordinate_values)
        code = _simple_code(
            "run_tsne",
            f"""            if output.isbacked or output.n_obs < 2 or not output.obs_names.is_unique:
                raise ValueError("t-SNE requires an in-memory AnnData with unique observations.")
            use_rep = {use_rep!r}
            if use_rep != "X" and use_rep not in output.obsm:
                raise ValueError("t-SNE representation is missing.")
            matrix = output.X if use_rep == "X" else output.obsm[use_rep]
            values = np.asarray(matrix.data if sparse.issparse(matrix) else matrix)
            if len(matrix.shape) != 2 or matrix.shape[0] != output.n_obs or matrix.shape[1] < 1:
                raise ValueError("t-SNE representation must be observation aligned.")
            if values.size and (not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values)):
                raise TypeError("t-SNE representation must be real numeric.")
            if values.size and not np.isfinite(values).all():
                raise ValueError("t-SNE representation contains non-finite values.")
            if {used!r} > matrix.shape[1] or not 0 < {perplexity!r} < output.n_obs:
                raise ValueError("t-SNE dimensions or perplexity exceed input bounds.")
            _openbio_validate_metric_rows(matrix, {used!r}, {metric!r}, "t-SNE")
            if {perplexity!r} < 1:
                warnings.warn("t-SNE perplexity is below 1; this is executable but unusually local.", UserWarning, stacklevel=2)
            if {perplexity!r} >= 0.5 * output.n_obs:
                warnings.warn("t-SNE perplexity is at least half the observation count; neighborhood interpretation may be unstable.", UserWarning, stacklevel=2)
            if ({key_added!r} in output.obsm or {key_added!r} in output.uns) and not {overwrite_existing!r}:
                raise ValueError("t-SNE output already exists.")
            if {overwrite_existing!r}:
                output.obsm.pop({key_added!r}, None)
                output.uns.pop({key_added!r}, None)
            backend_rep, backend_dims = use_rep, {used!r}
            temp_key = "__openbio_tsne_representation"
            if use_rep == "X" and {used!r} < matrix.shape[1]:
                output.obsm[temp_key] = matrix[:, :{used!r}].copy()
                backend_rep, backend_dims = temp_key, None
            sc.tl.tsne(output, n_pcs=backend_dims, n_components=2, use_rep=backend_rep, perplexity={perplexity!r}, metric={metric!r}, early_exaggeration={early_exaggeration!r}, learning_rate={learning_rate!r}, random_state={random_seed!r}, use_fast_tsne=False, n_jobs=None, key_added={key_added!r}, copy=False)
            output.obsm.pop(temp_key, None)
            coordinates = np.asarray(output.obsm.get({key_added!r}))
            if coordinates.shape != (output.n_obs, 2):
                raise RuntimeError("t-SNE backend returned invalid coordinate shape.")
            if not np.issubdtype(coordinates.dtype, np.number) or np.iscomplexobj(coordinates):
                raise TypeError("t-SNE backend coordinates must be real numeric values.")
            if not np.isfinite(coordinates).all():
                raise RuntimeError("t-SNE backend returned non-finite coordinates.")
            metadata = output.uns.get({key_added!r})
            if not isinstance(metadata, Mapping) or not isinstance(metadata.get("params"), Mapping):
                raise RuntimeError("t-SNE backend did not return the expected parameter record.")""",
        )
        code = code.replace(
            'temp_key = "__openbio_tsne_representation"',
            'temp_prefix = "__openbio_tsne_representation"\n'
            "    temp_key, suffix = temp_prefix, 0\n"
            "    while temp_key in output.obsm:\n"
            "        suffix += 1\n"
            '        temp_key = f"{temp_prefix}_{suffix}"',
        )
        code = code.replace(
            "from collections.abc import Mapping", "import warnings\nfrom collections.abc import Mapping", 1
        )
        code = _generated_graph_helper() + code
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellTSNE",
            title="t-SNE",
            operation="tsne",
            methods=f"Computed a two-dimensional t-SNE embedding from {used} dimensions of {use_rep!r} with {metric} distance.",
            results=f"Generated {key_added!r} coordinates for {adata.n_obs} observations.",
            key_results={"actual_dimensions": used, "coordinates": diagnostics, "metric_row_warnings": metric_warnings},
            parameters=parameters,
            references=[TSNE_REFERENCE, SCANPY_REFERENCE],
            software_packages=SOFTWARE_PACKAGES,
            warnings=[
                *metric_warnings,
                *(["t-SNE perplexity is below 1; this is executable but unusually local."] if perplexity < 1 else []),
                *(
                    [
                        "t-SNE perplexity is at least half the observation count; neighborhood interpretation may be unstable."
                    ]
                    if perplexity >= 0.5 * int(adata.n_obs)
                    else []
                ),
            ],
            limitations=[
                "t-SNE emphasizes local neighborhoods; axes and global distances are not directly interpretable."
            ],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


class OpenBioSingleCellForceDirectedGraph:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        layout: str = "fr",
        init_mode: str = "random",
        init_key: str = "X_draw_graph_fr",
        neighbors_key: str = "neighbors",
        key_suffix: str = "",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_adata(adata, "Force-directed graph", require_variables=False)
        random_seed = validate_random_seed(random_seed)
        graph = resolve_named_graph(adata, neighbors_key, operation="Force-directed graph")
        if layout not in SUPPORTED_GRAPH_LAYOUTS:
            raise ValueError(f"Unsupported graph layout: {layout!r}")
        if init_mode not in {"random", "existing", "paga"}:
            raise ValueError(f"Unsupported graph initialization: {init_mode!r}")
        if layout == "fa" and find_spec("fa2_modified") is None:
            raise RuntimeError(
                "ForceAtlas2 layout requires the optional 'fa2-modified' package; no fallback was applied."
            )
        backend = "fa2-modified" if layout == "fa" else "python-igraph"
        layout_kwargs = {"maxiter": 500} if layout in {"fa", "kk"} else {"niter": 500}
        backend_init: str | None = None
        if init_mode == "existing":
            if init_key not in adata.obsm:
                raise ValueError(f"Existing graph initialization not found in obsm: {init_key!r}")
            coords = adata.obsm[init_key]
            _validate_matrix(coords, operation="Force-directed graph initialization", expected_rows=int(adata.n_obs))
            if int(coords.shape[1]) != 2:
                raise ValueError("Existing graph initialization must have exactly two columns.")
            backend_init = init_key
        elif init_mode == "paga":
            _validate_paga_initialization(adata)
            backend_init = "paga"
        if not isinstance(key_suffix, str):
            raise TypeError("Force-directed graph key_suffix must be a string.")
        suffix = key_suffix.strip() or layout
        coord_key = f"X_draw_graph_{suffix}"
        if (coord_key in adata.obsm or "draw_graph" in adata.uns) and not overwrite_existing:
            raise ValueError("Force-directed graph coordinate/parameter output already exists.")
        started_at = time.perf_counter()
        output = adata
        if overwrite_existing:
            output.obsm.pop(coord_key, None)
            output.uns.pop("draw_graph", None)
        output.obsp[graph.connectivities_key] = graph.connectivities
        with _DRAW_GRAPH_RNG_LOCK:
            numpy_state = science.np.random.get_state()
            python_state = random.getstate()
            try:
                science.sc.tl.draw_graph(
                    output,
                    layout=layout,
                    init_pos=backend_init,
                    root=None,
                    neighbors_key=neighbors_key,
                    key_added_ext=suffix,
                    random_state=random_seed,
                    n_jobs=None,
                    adjacency=None,
                    obsp=None,
                    copy=False,
                    **layout_kwargs,
                )
            finally:
                science.np.random.set_state(numpy_state)
                random.setstate(python_state)
        coordinate_values = _validate_embedding_coordinates(
            output.obsm.get(coord_key), n_obs=int(output.n_obs), operation="Force-directed graph"
        )
        draw_metadata = output.uns.get("draw_graph")
        draw_params = draw_metadata.get("params") if isinstance(draw_metadata, Mapping) else None
        if not isinstance(draw_params, Mapping):
            raise RuntimeError("Force-directed graph backend did not return the expected parameter record.")
        actual_layout = draw_params.get("layout")
        if actual_layout != layout:
            raise RuntimeError(f"Requested graph layout {layout!r}, but backend reported {actual_layout!r}.")
        parameters = {
            "layout": layout,
            "backend": backend,
            "init_mode": init_mode,
            "init_key": init_key if init_mode == "existing" else None,
            "neighbors_key": neighbors_key,
            "key_suffix": suffix,
            "root": None,
            "n_jobs": None,
            "layout_kwargs": layout_kwargs,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "force_directed_graph",
            parameters,
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
            random_seed=random_seed,
        )
        diagnostics = _coordinate_diagnostics(coordinate_values)
        graph_report = graph_diagnostics(graph)
        code = dedent(
            f"""import random
import threading
from importlib.util import find_spec

import numpy as np
import pandas as pd
import scanpy as sc


_OPENBIO_DRAW_GRAPH_RNG_LOCK = threading.RLock()


def _openbio_validate_paga(adata):
    paga = adata.uns.get("paga")
    if not isinstance(paga, Mapping):
        raise ValueError("PAGA initialization requires an uns['paga'] mapping; run PAGA first.")
    groups = paga.get("groups")
    if not isinstance(groups, str) or not groups or groups not in adata.obs:
        raise ValueError("PAGA initialization requires its groups column in obs.")
    labels = adata.obs[groups]
    if not isinstance(labels.dtype, pd.CategoricalDtype):
        raise ValueError("PAGA initialization groups must be categorical.")
    if labels.isna().any() or len(labels.cat.categories) < 1:
        raise ValueError("PAGA initialization groups must be complete and categorical.")
    n_categories = len(labels.cat.categories)
    positions = np.asarray(paga.get("pos"))
    if positions.shape != (n_categories, 2):
        raise ValueError("PAGA initialization positions must align to group categories and have two columns.")
    if not np.issubdtype(positions.dtype, np.number) or np.iscomplexobj(positions):
        raise TypeError("PAGA initialization positions must be real numeric values.")
    if not np.isfinite(positions).all():
        raise ValueError("PAGA initialization positions must be finite.")
    connectivities = paga.get("connectivities")
    if not sparse.issparse(connectivities):
        raise TypeError("PAGA initialization connectivities must be a SciPy sparse matrix.")
    if connectivities.shape != (n_categories, n_categories):
        raise ValueError("PAGA initialization connectivities must align to group categories.")
    connectivities = connectivities.tocsr(copy=True)
    connectivities.sum_duplicates()
    connectivities.eliminate_zeros()
    weights = np.asarray(connectivities.data)
    if weights.size and (not np.issubdtype(weights.dtype, np.number) or np.iscomplexobj(weights)):
        raise TypeError("PAGA initialization connectivities must contain real numeric weights.")
    if weights.size and (not np.isfinite(weights).all() or (weights < 0).any()):
        raise ValueError("PAGA initialization connectivities must contain finite non-negative weights.")
    if not np.allclose(connectivities.diagonal(), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("PAGA initialization connectivities must have a zero diagonal.")
    delta = connectivities - connectivities.T
    delta.eliminate_zeros()
    if delta.nnz and float(np.max(np.abs(delta.data))) > 1e-8:
        raise ValueError("PAGA initialization connectivities must be symmetric.")


def run_force_directed_graph(adata):
    output = adata
    _openbio_canonical_graph(output, {neighbors_key!r}, "Force-directed graph")
    if {layout!r} == "fa" and find_spec("fa2_modified") is None:
        raise RuntimeError("ForceAtlas2 requires fa2-modified; no fallback was applied.")
    init_pos = {backend_init!r}
    if init_pos == "paga":
        _openbio_validate_paga(output)
    elif init_pos is not None:
        if init_pos not in output.obsm:
            raise ValueError("Existing initialization is missing.")
        initial = np.asarray(output.obsm[init_pos])
        if initial.shape != (output.n_obs, 2):
            raise ValueError("Existing initialization must have shape (n_obs, 2).")
        if not np.issubdtype(initial.dtype, np.number) or np.iscomplexobj(initial):
            raise TypeError("Existing initialization must contain real numeric coordinates.")
        if not np.isfinite(initial).all():
            raise ValueError("Existing initialization must contain finite coordinates.")
    if ({coord_key!r} in output.obsm or "draw_graph" in output.uns) and not {overwrite_existing!r}:
        raise ValueError("Force-directed graph output already exists.")
    if {overwrite_existing!r}:
        output.obsm.pop({coord_key!r}, None)
        output.uns.pop("draw_graph", None)
    with _OPENBIO_DRAW_GRAPH_RNG_LOCK:
        np_state, py_state = np.random.get_state(), random.getstate()
        try:
            sc.tl.draw_graph(output, layout={layout!r}, init_pos=init_pos, root=None, neighbors_key={neighbors_key!r}, key_added_ext={suffix!r}, random_state={random_seed!r}, n_jobs=None, adjacency=None, obsp=None, copy=False, **{layout_kwargs!r})
        finally:
            np.random.set_state(np_state)
            random.setstate(py_state)
    coordinates = np.asarray(output.obsm.get({coord_key!r}))
    if coordinates.shape != (output.n_obs, 2):
        raise RuntimeError("Force-directed graph backend returned invalid coordinate shape.")
    if not np.issubdtype(coordinates.dtype, np.number) or np.iscomplexobj(coordinates):
        raise TypeError("Force-directed graph backend coordinates must be real numeric values.")
    if not np.isfinite(coordinates).all():
        raise RuntimeError("Force-directed graph backend returned non-finite coordinates.")
    metadata = output.uns.get("draw_graph")
    params = metadata.get("params") if isinstance(metadata, Mapping) else None
    if not isinstance(params, Mapping):
        raise RuntimeError("Force-directed graph backend did not return the expected parameter record.")
    if params.get("layout") != {layout!r}:
        raise RuntimeError("Graph backend did not use the requested layout.")
    return output
"""
        )
        code = _generated_graph_helper() + code
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellForceDirectedGraph",
            title="Force-directed graph",
            operation="force_directed_graph",
            methods=f"Computed the requested {layout!r} layout with {backend} from the explicitly named undirected graph; initialization was {init_mode}, with a fixed 500-iteration policy.",
            results=f"Generated {coord_key!r} coordinates for {adata.n_obs} observations.",
            key_results={
                "actual_layout": actual_layout,
                "backend": backend,
                "coordinates": diagnostics,
                "graph": graph_report,
            },
            parameters=parameters,
            references=_layout_references(layout),
            software_packages=SOFTWARE_PACKAGES + (("fa2-modified",) if layout == "fa" else ("igraph",)),
            warnings=_graph_warnings(graph_report),
            limitations=[
                "Force-directed layouts are stochastic visual summaries; geometric distances are not inferential effect sizes."
            ],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


def _leiden_size_diagnostics(cluster_sizes: tuple[tuple[str, int], ...], *, n_obs: int) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    sizes = dict(cluster_sizes)
    counts = science.np.asarray(list(sizes.values()), dtype=int)
    if counts.size < 1 or int(counts.sum()) != n_obs or bool((counts <= 0).any()):
        raise RuntimeError("Leiden cluster-size diagnostics do not cover every observation exactly once.")
    return {
        "cluster_count": len(sizes),
        "cluster_sizes": sizes,
        "cluster_fractions": {label: float(size / n_obs) for label, size in sizes.items()},
        "min_cluster_size": int(counts.min()),
        "median_cluster_size": float(science.np.median(counts)),
        "max_cluster_size": int(counts.max()),
        "singleton_count": int((counts == 1).sum()),
        "largest_cluster_fraction": float(counts.max() / n_obs),
    }


def _leiden_code(parameters: dict[str, Any]) -> str:
    neighbors_key = parameters["neighbors_key"]
    key_added = parameters["key_added"]
    overwrite = parameters["overwrite_existing"]
    repeats = parameters["stability_repeats"]
    seed = parameters["random_seed"]
    resolution = parameters["resolution"]
    iterations = parameters["n_iterations"]
    return dedent(
        f"""from collections.abc import Mapping

import igraph
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
from scipy import sparse
from sklearn.metrics import adjusted_rand_score


def run_leiden(adata):
    if adata.isbacked or adata.n_obs < 2:
        raise ValueError("Leiden requires an in-memory AnnData with at least two observations.")
    if not adata.obs_names.is_unique:
        raise ValueError("Leiden requires unique observation identifiers.")
    metadata = adata.uns.get({neighbors_key!r})
    if not isinstance(metadata, Mapping):
        raise ValueError("Named neighbor graph is missing.")
    connectivity_key = metadata.get("connectivities_key")
    if not isinstance(connectivity_key, str) or not connectivity_key.strip():
        raise ValueError("Named connectivity pointer is invalid.")
    connectivity_key = connectivity_key.strip()
    if connectivity_key not in adata.obsp:
        raise ValueError("Named connectivity matrix is missing.")
    adjacency = adata.obsp[connectivity_key]
    if not sparse.issparse(adjacency) or adjacency.shape != (adata.n_obs, adata.n_obs):
        raise ValueError("Leiden requires an observation-aligned sparse graph.")
    adjacency = adjacency.tocsr(copy=True)
    adjacency.sum_duplicates()
    adjacency.eliminate_zeros()
    values = np.asarray(adjacency.data)
    if values.size and not np.issubdtype(values.dtype, np.number):
        raise TypeError("Leiden graph weights must be numeric.")
    if values.size and (not np.isfinite(values).all() or (values < 0).any()):
        raise ValueError("Leiden graph has invalid weights.")
    if not np.allclose(adjacency.diagonal(), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Leiden graph diagonal must be zero.")
    delta = adjacency - adjacency.T
    delta.eliminate_zeros()
    if delta.nnz and float(np.max(np.abs(delta.data))) > 1e-8:
        raise ValueError("Leiden requires an undirected graph.")
    adjacency = ((adjacency + adjacency.T) * 0.5).tocsr()
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    if adjacency.nnz == 0:
        raise ValueError("Leiden graph has no positive edges.")
    if {key_added!r} == {neighbors_key!r}:
        raise ValueError("Leiden key_added cannot equal the resolved neighbors_key because that uns key stores graph metadata.")
    if ({key_added!r} in adata.obs or {key_added!r} in adata.uns) and not {overwrite!r}:
        raise ValueError("Leiden output already exists.")
    output = adata
    if {overwrite!r}:
        output.obs.drop(columns=[{key_added!r}], inplace=True, errors="ignore")
        output.uns.pop({key_added!r}, None)
    upper = sparse.triu(adjacency, k=1).tocoo()
    graph = igraph.Graph(n=output.n_obs, edges=list(zip(upper.row.tolist(), upper.col.tolist(), strict=True)), directed=False)
    graph.es["weight"] = [float(value) for value in upper.data]
    def _partition(key, offset):
        try:
            sc.tl.leiden(output, adjacency=adjacency, resolution={resolution!r}, key_added=key, random_state={seed!r} + offset, flavor="igraph", n_iterations={iterations!r}, directed=False, use_weights=True, objective_function="modularity")
            labels = output.obs.get(key)
            if labels is None or not isinstance(labels.dtype, pd.CategoricalDtype):
                raise RuntimeError("Leiden backend did not return categorical memberships.")
            if len(labels) != output.n_obs or not labels.index.equals(output.obs_names):
                raise RuntimeError("Leiden memberships do not align to the graph observations.")
            if labels.isna().any():
                raise RuntimeError("Leiden backend returned missing memberships.")
            observed_values = tuple(str(value) for value in labels.astype(object))
            observed = tuple(dict.fromkeys(observed_values))
            mapping = {{label: str(index) for index, label in enumerate(observed)}}
            membership = tuple(mapping[label] for label in observed_values)
            categories = tuple(str(index) for index in range(len(observed)))
            output.obs[key] = pd.Series(
                pd.Categorical(membership, categories=categories),
                index=output.obs_names,
                name=key,
            )
            labels = output.obs[key]
            counts = labels.value_counts(sort=False)
            if int(counts.sum()) != output.n_obs or (counts <= 0).any():
                raise RuntimeError("Leiden cluster sizes do not cover every observation exactly once.")
            result_metadata = output.uns.get(key)
            backend_value = result_metadata.get("modularity") if isinstance(result_metadata, Mapping) else None
            if isinstance(backend_value, bool):
                raise RuntimeError("Leiden backend modularity is missing or invalid.")
            try:
                backend_modularity = float(backend_value)
            except (TypeError, ValueError) as error:
                raise RuntimeError("Leiden backend modularity is missing or invalid.") from error
            if not np.isfinite(backend_modularity):
                raise RuntimeError("Leiden backend modularity is non-finite.")
            label_ids = {{label: index for index, label in enumerate(dict.fromkeys(membership))}}
            calculated_modularity = float(graph.modularity([label_ids[label] for label in membership], weights=graph.es["weight"], resolution={resolution!r}, directed=False))
            if not np.isfinite(calculated_modularity) or not np.isclose(backend_modularity, calculated_modularity, rtol=1e-9, atol=1e-12):
                raise RuntimeError("Leiden backend modularity is inconsistent with the validated graph partition.")
            return membership, calculated_modularity
        finally:
            if offset:
                output.obs.drop(columns=[key], inplace=True, errors="ignore")
                output.uns.pop(key, None)

    memberships = []
    modularities = []
    for offset in range({repeats!r}):
        key = {key_added!r} if offset == 0 else f"__openbio_leiden_stability_{{offset}}"
        suffix = 0
        while offset and (key in output.obs or key in output.uns):
            suffix += 1
            key = f"__openbio_leiden_stability_{{offset}}_{{suffix}}"
        membership, calculated_modularity = _partition(key, offset)
        memberships.append(membership)
        modularities.append(calculated_modularity)
    pairwise = [float(adjusted_rand_score(memberships[left], memberships[right])) for left in range(len(memberships)) for right in range(left + 1, len(memberships))]
    if pairwise and (not np.isfinite(pairwise).all() or any(not -1.0 <= value <= 1.0 for value in pairwise)):
        raise RuntimeError("Leiden stability ARI diagnostics are invalid.")
    counts = output.obs[{key_added!r}].value_counts(sort=False)
    cluster_sizes = {{str(label): int(count) for label, count in counts.items()}}
    cluster_counts = np.asarray(list(cluster_sizes.values()), dtype=int)
    if cluster_counts.size < 1 or int(cluster_counts.sum()) != output.n_obs or (cluster_counts <= 0).any():
        raise RuntimeError("Leiden cluster-size diagnostics do not cover every observation exactly once.")
    size_diagnostics = {{
        "cluster_count": len(cluster_sizes),
        "cluster_sizes": cluster_sizes,
        "cluster_fractions": {{label: float(size / output.n_obs) for label, size in cluster_sizes.items()}},
        "min_cluster_size": int(cluster_counts.min()),
        "median_cluster_size": float(np.median(cluster_counts)),
        "max_cluster_size": int(cluster_counts.max()),
        "singleton_count": int((cluster_counts == 1).sum()),
        "largest_cluster_fraction": float(cluster_counts.max() / output.n_obs),
    }}
    modularity = modularities[0]
    output.uns[{key_added!r}]["openbio_diagnostics"] = {{
        "modularity": modularity,
        **size_diagnostics,
        "stability": {{"starts": {repeats!r}, "seeds": [{seed!r} + offset for offset in range({repeats!r})], "pairwise_ari": pairwise, "mean_ari": float(np.mean(pairwise)) if pairwise else None, "min_ari": min(pairwise) if pairwise else None, "max_ari": max(pairwise) if pairwise else None}},
    }}
    if {resolution!r} == 0:
        warnings.warn("Leiden resolution=0 is executable but commonly collapses the graph into very few communities.", UserWarning, stacklevel=2)
    if {iterations!r} == 0:
        warnings.warn("Leiden n_iterations=0 requests no optimization iterations; interpret the returned partition cautiously.", UserWarning, stacklevel=2)
    if {repeats!r} == 1:
        warnings.warn("Leiden stability was not assessed because stability_repeats=1.", UserWarning, stacklevel=2)
    return output
"""
    )


class OpenBioSingleCellLeiden:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resolution: float = 1.0,
        key_added: str = "leiden",
        neighbors_key: str = "neighbors",
        n_iterations: int = 2,
        stability_repeats: int = 5,
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        _require_adata(adata, "Leiden", require_variables=False)
        if not isinstance(key_added, str):
            raise TypeError("Leiden key_added must be a string.")
        if not key_added.strip():
            raise ValueError("Leiden key_added cannot be empty.")
        resolution, n_iterations, stability_repeats = validate_leiden_settings(
            resolution, n_iterations, stability_repeats
        )
        random_seed = validate_random_seed(random_seed, stability_repeats=stability_repeats)
        graph = resolve_named_graph(adata, neighbors_key, operation="Leiden")
        if isinstance(key_added, str) and key_added.strip() == graph.neighbors_key:
            raise ValueError(
                "Leiden key_added cannot equal the resolved neighbors_key because that uns key stores graph metadata."
            )
        key_added = validate_graph_result_key(
            adata, key_added, operation="Leiden", obs=True, uns=True, overwrite_existing=overwrite_existing
        )
        started_at = time.perf_counter()
        output = adata
        if overwrite_existing:
            output.obs.drop(columns=[key_added], inplace=True, errors="ignore")
            output.uns.pop(key_added, None)
        partition = run_leiden_partition(
            output,
            graph,
            resolution=resolution,
            key_added=key_added,
            random_seed=random_seed,
            n_iterations=n_iterations,
        )
        stability = assess_leiden_stability(
            output,
            graph,
            base_membership=partition.membership,
            resolution=resolution,
            random_seed=random_seed,
            n_iterations=n_iterations,
            repeats=stability_repeats,
        )
        diagnostics = {
            "starts": stability.starts,
            "seeds": list(stability.seeds),
            "pairwise_ari": list(stability.pairwise_ari),
            "mean_ari": stability.mean_ari,
            "min_ari": stability.min_ari,
            "max_ari": stability.max_ari,
        }
        size_diagnostics = _leiden_size_diagnostics(partition.cluster_sizes, n_obs=int(output.n_obs))
        output.uns[key_added]["openbio_diagnostics"] = {
            "modularity": partition.modularity,
            **size_diagnostics,
            "stability": diagnostics,
        }
        parameters = {
            "resolution": resolution,
            "key_added": key_added,
            "neighbors_key": graph.neighbors_key,
            "flavor": "igraph",
            "directed": False,
            "use_weights": True,
            "objective_function": "modularity",
            "n_iterations": n_iterations,
            "stability_repeats": stability_repeats,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        finish_adata(
            output, "leiden", parameters, int(adata.n_obs), int(adata.n_vars), started_at, random_seed=random_seed
        )
        graph_report = graph_diagnostics(graph)
        code = _leiden_code(parameters)
        summary, code = make_analysis_report(
            node_id="OpenBioSingleCellLeiden",
            title="Leiden clustering",
            operation="leiden",
            methods="Partitioned the weighted, undirected named graph with Scanpy's igraph Leiden backend and the modularity objective.",
            results=f"Identified {size_diagnostics['cluster_count']} graph communities; {size_diagnostics['singleton_count']} were singletons, and the largest contained {size_diagnostics['largest_cluster_fraction']:.1%} of observations.",
            key_results={
                **size_diagnostics,
                "modularity": partition.modularity,
                "stability": diagnostics,
                "graph": graph_report,
            },
            parameters=parameters,
            references=[
                LEIDEN_REFERENCE,
                SCANPY_REFERENCE,
                *([ARI_REFERENCE] if stability_repeats > 1 else []),
            ],
            software_packages=SOFTWARE_PACKAGES + ("pandas", "igraph"),
            warnings=[
                *(
                    ["Leiden resolution=0 is executable but commonly collapses the graph into very few communities."]
                    if resolution == 0
                    else []
                ),
                *(
                    [
                        "Leiden n_iterations=0 requests no optimization iterations; interpret the returned partition cautiously."
                    ]
                    if n_iterations == 0
                    else []
                ),
                *(["Leiden stability was not assessed because stability_repeats=1."] if stability_repeats == 1 else []),
                *_graph_warnings(graph_report),
            ],
            limitations=[
                "Communities are algorithmic graph partitions, not validated biological cell types or ground truth."
            ],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, summary, code


def _run_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
    *,
    name: str,
    expected_parameters: set[str],
    runner: Any,
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation=name)
    require_parameters(parameters, expected_parameters, operation=name)
    adata = read_anndata_input(inputs)
    output, summary, code = runner(adata, **parameters)
    return analysis_outputs(summary, code, write_anndata_output(context, output))


@register_operation("openbio.node.pca")
def pca(context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="PCA",
        expected_parameters={"n_comps", "use_hvg", "source", "overwrite_existing", "max_output_gib", "random_seed"},
        runner=OpenBioSingleCellPCA.execute,
    )


@register_operation("openbio.node.pcavarianceplot")
def pca_variance_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="PCA Variance Plot")
    require_parameters(parameters, {"n_pcs"}, operation="PCA Variance Plot")
    plotted, report, code = pca_variance_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.pcaloadingsplot")
def pca_loadings_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="PCA Loadings Plot")
    require_parameters(parameters, {"component", "n_genes"}, operation="PCA Loadings Plot")
    plotted, report, code = pca_loadings_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.neighbors")
def neighbors(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="Neighbors",
        expected_parameters={
            "use_rep",
            "n_dimensions",
            "n_neighbors",
            "metric",
            "method",
            "key_added",
            "overwrite_existing",
            "random_seed",
        },
        runner=OpenBioSingleCellNeighbors.execute,
    )


@register_operation("openbio.node.neighborgraphdiagnosticsplot")
def neighbor_graph_diagnostics_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Neighbor Graph Diagnostics Plot")
    require_parameters(
        parameters,
        {"neighbors_key", "max_working_memory_gib"},
        operation="Neighbor Graph Diagnostics Plot",
    )
    plotted, report, code = neighbor_graph_diagnostics_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.umap")
def umap(context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="UMAP",
        expected_parameters={"neighbors_key", "min_dist", "spread", "key_added", "overwrite_existing", "random_seed"},
        runner=OpenBioSingleCellUMAP.execute,
    )


@register_operation("openbio.node.tsne")
def tsne(context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="t-SNE",
        expected_parameters={
            "use_rep",
            "n_dimensions",
            "perplexity",
            "metric",
            "early_exaggeration",
            "learning_rate",
            "key_added",
            "overwrite_existing",
            "random_seed",
        },
        runner=OpenBioSingleCellTSNE.execute,
    )


@register_operation("openbio.node.forcedirectedgraph")
def force_directed_graph(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="Force-Directed Graph",
        expected_parameters={
            "layout",
            "init_mode",
            "init_key",
            "neighbors_key",
            "key_suffix",
            "overwrite_existing",
            "random_seed",
        },
        runner=OpenBioSingleCellForceDirectedGraph.execute,
    )


@register_operation("openbio.node.leiden")
def leiden(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    return _run_operation(
        context,
        inputs,
        parameters,
        name="Leiden",
        expected_parameters={
            "resolution",
            "key_added",
            "neighbors_key",
            "n_iterations",
            "stability_repeats",
            "overwrite_existing",
            "random_seed",
        },
        runner=OpenBioSingleCellLeiden.execute,
    )


__all__ = [
    "force_directed_graph",
    "leiden",
    "neighbor_graph_diagnostics_plot",
    "neighbor_graph_diagnostics_plot_owned",
    "neighbors",
    "pca",
    "pca_loadings_plot",
    "pca_loadings_plot_owned",
    "pca_variance_plot",
    "pca_variance_plot_owned",
    "tsne",
    "umap",
]
