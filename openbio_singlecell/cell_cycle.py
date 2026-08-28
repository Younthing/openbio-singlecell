from __future__ import annotations

import hashlib
import inspect
import json
import threading
from collections.abc import Sequence
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import collect_software_versions, summarize_numeric

if TYPE_CHECKING:
    from anndata import AnnData


CELL_CYCLE_NODE_ID = "OpenBioSingleCellCellCycleScore"
CELL_CYCLE_SCHEMA = "openbio-singlecell/cell-cycle-score/v1"
_CELL_CYCLE_RNG_LOCK = threading.RLock()

REGEV_HUMAN_S_GENES = (
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
)

REGEV_HUMAN_G2M_GENES = (
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
)

REGEV_HUMAN_RESOURCE_SHA256 = "eb666ba4f7b7a8b2845e1586d42687feaf813e9fa1fdadba89dacd4093411ac0"

CELL_CYCLE_REFERENCES = [
    {
        "citation": "Satija R, Farrell JA, Gennert D, et al. Spatial reconstruction of single-cell gene expression data. Nature Biotechnology. 2015;33:495-502.",
        "doi": "10.1038/nbt.3192",
        "url": "https://doi.org/10.1038/nbt.3192",
        "kind": "method",
    },
    {
        "citation": "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016;352:189-196.",
        "doi": "10.1126/science.aad0501",
        "url": "https://doi.org/10.1126/science.aad0501",
        "kind": "method",
    },
    {
        "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15.",
        "doi": "10.1186/s13059-017-1382-0",
        "url": "https://doi.org/10.1186/s13059-017-1382-0",
        "kind": "software",
    },
]


def _resource_sha256(s_genes: Sequence[str], g2m_genes: Sequence[str]) -> str:
    payload = "\n".join([*s_genes, *g2m_genes]) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


if (
    len(REGEV_HUMAN_S_GENES) != 43
    or len(REGEV_HUMAN_G2M_GENES) != 54
    or _resource_sha256(REGEV_HUMAN_S_GENES, REGEV_HUMAN_G2M_GENES) != REGEV_HUMAN_RESOURCE_SHA256
):
    raise RuntimeError("Bundled Regev human cell-cycle resource failed its immutable identity check.")


def _strict_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value or value != value.strip():
        raise ValueError(f"{label} must be nonempty and free of surrounding whitespace.")
    return value


def _parse_custom_genes(value: Any, *, label: str) -> tuple[str, ...]:
    if isinstance(value, str):
        parts = value.split(",")
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        parts = list(value)
    else:
        raise TypeError(f"{label} must be a comma-separated string or sequence of strings.")
    genes = []
    for position, item in enumerate(parts):
        if not isinstance(item, str):
            raise TypeError(f"{label} item {position} must be a string.")
        gene = item.strip()
        if not gene:
            raise ValueError(f"{label} contains a blank gene identifier.")
        genes.append(gene)
    if len(genes) != len(set(genes)):
        raise ValueError(f"{label} contains duplicate gene identifiers.")
    return tuple(genes)


def _canonical_axes(adata: AnnData) -> tuple[list[str], list[str]]:
    try:
        if getattr(adata, "isbacked", False):
            raise ValueError("Cell Cycle Score requires an in-memory AnnData; call to_memory() first.")
        obs_names = list(adata.obs_names)
        var_names = list(adata.var_names)
    except AttributeError as error:
        raise TypeError("Cell Cycle Score requires an AnnData input.") from error
    if not obs_names or not var_names:
        raise ValueError("Cell Cycle Score requires nonempty cell and feature axes.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell Cycle Score requires unique cell and feature identifiers.")
    for values, label in ((obs_names, "Cell"), (var_names, "Feature")):
        if any(not isinstance(value, str) or not value or value != value.strip() for value in values):
            raise ValueError(f"{label} identifiers must be canonical nonempty strings.")
    return obs_names, var_names


def _matrix_fingerprint(matrix: Any, *, obs_names: Sequence[str], var_names: Sequence[str], science: Any) -> str:
    digest = hashlib.sha256()
    header = json.dumps(
        {"obs_names": list(obs_names), "var_names": list(var_names)},
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    digest.update(header.encode("utf-8"))
    if science.sparse.issparse(matrix):
        canonical = science.sparse.csr_matrix(matrix, dtype=science.np.float64).copy()
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        digest.update(b"csr-f8")
        digest.update(science.np.asarray(canonical.shape, dtype="<i8").tobytes())
        digest.update(science.np.asarray(canonical.indptr, dtype="<i8").tobytes())
        digest.update(science.np.asarray(canonical.indices, dtype="<i8").tobytes())
        digest.update(science.np.asarray(canonical.data, dtype="<f8").tobytes())
    else:
        canonical = science.np.ascontiguousarray(matrix, dtype="<f8")
        digest.update(b"dense-f8")
        digest.update(science.np.asarray(canonical.shape, dtype="<i8").tobytes())
        digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


def _validate_expression(
    adata: AnnData,
    *,
    source_kind: str,
    layer_name: str | None,
    science: Any,
) -> tuple[Any, str, list[str], str, dict[str, Any], list[str]]:
    if source_kind == "raw":
        if layer_name is not None:
            raise ValueError("Cell Cycle Score Raw source cannot carry a layer name.")
        if adata.raw is None:
            raise ValueError("Cell Cycle Score selected Raw, but this AnnData has no Raw snapshot.")
        if not adata.raw.obs_names.equals(adata.obs_names):
            raise ValueError("Cell Cycle Score Raw observations are not aligned to current cells.")
        matrix = adata.raw.X
        source_var_names = list(adata.raw.var_names)
        source_index = adata.raw.var_names
        source_label = "raw.X"
    elif source_kind == "X":
        if layer_name is not None:
            raise ValueError("Cell Cycle Score X source cannot carry a layer name.")
        matrix = adata.X
        source_var_names = list(adata.var_names)
        source_index = adata.var_names
        source_label = "X"
    elif source_kind == "layer":
        layer_name = _strict_text(layer_name, label="Cell Cycle Score layer name")
        if layer_name not in adata.layers:
            raise ValueError(f"Cell Cycle Score layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        source_var_names = list(adata.var_names)
        source_index = adata.var_names
        source_label = f"layers[{layer_name!r}]"
    else:
        raise ValueError(f"Unsupported Cell Cycle Score expression source: {source_kind!r}.")
    if not source_index.is_unique or any(
        not isinstance(value, str) or not value or value != value.strip() for value in source_var_names
    ):
        raise ValueError("Cell Cycle Score selected source requires unique canonical feature identifiers.")
    expected_shape = (int(adata.n_obs), len(source_var_names))
    if getattr(matrix, "shape", None) != expected_shape:
        raise ValueError("Cell Cycle Score expression source is not aligned to its declared cell/feature axes.")
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    if (
        science.np.issubdtype(values.dtype, science.np.bool_)
        or science.np.issubdtype(values.dtype, science.np.complexfloating)
        or not science.np.issubdtype(values.dtype, science.np.number)
    ):
        raise TypeError("Cell Cycle Score expression must be real numeric values.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError("Cell Cycle Score expression contains NaN or infinity.")
    numeric = values.astype(float, copy=False)
    integer_like = bool(science.np.allclose(numeric, science.np.rint(numeric), rtol=0.0, atol=1e-8))
    contains_negative = bool(numeric.size and (numeric < 0).any())
    total_values = expected_shape[0] * expected_shape[1]
    if science.sparse.issparse(matrix):
        has_implicit_zero = int(matrix.nnz) < total_values
        constant = bool(
            (not numeric.size)
            or (
                not has_implicit_zero
                and float(science.np.ptp(numeric)) == 0.0
            )
            or (has_implicit_zero and bool((numeric == 0).all()))
        )
    else:
        constant = bool((not numeric.size) or float(science.np.ptp(numeric)) == 0.0)
    if constant:
        state = "constant_user_selected"
    elif contains_negative:
        state = "signed_user_selected"
    elif integer_like:
        state = "nonnegative_integer_like_user_selected"
    else:
        state = "nonnegative_noninteger_user_selected"
    evidence = {
        "basis": "explicit_source_and_value_profile",
        "history_used": False,
        "biological_state_verified": False,
        "integer_like": integer_like,
        "contains_negative": contains_negative,
        "constant": constant,
    }
    advisories = [
        "Expression state and full-gene completeness are user-declared; full-gene log-normalized expression is the recommended Cell Cycle Score input."
    ]
    if source_kind == "raw":
        advisories.append("Raw was selected explicitly; the Raw slot itself does not certify counts or transformed state.")
    if integer_like and not contains_negative:
        advisories.append("The selected source is integer-like and may be counts; phase scores can differ from the recommended log-normalized workflow.")
    if contains_negative:
        advisories.append("The selected source contains negative values and may be scaled/residual expression; interpretation is expert-declared.")
    if constant:
        advisories.append("The selected source is constant, so cell-cycle score separation may be degenerate.")
    return matrix, source_label, source_var_names, state, evidence, advisories


def analyze_cell_cycle_score(
    adata: AnnData,
    *,
    source_kind: str,
    layer_name: str | None,
    gene_set_source: str,
    organism: str,
    s_genes: Any,
    g2m_genes: Any,
    output_prefix: str,
    overwrite_existing: bool,
    random_seed: int,
) -> tuple[AnnData, dict[str, Any]]:
    from . import dependencies

    science = dependencies.require_scientific_dependencies()
    obs_names, current_var_names = _canonical_axes(adata)
    matrix, source_label, var_names, source_state, state_evidence, expression_advisories = _validate_expression(
        adata,
        source_kind=source_kind,
        layer_name=layer_name,
        science=science,
    )
    if organism not in {"human", "mouse", "other"}:
        raise ValueError(f"Unsupported Cell Cycle Score organism: {organism!r}.")
    if gene_set_source == "regev_human_97":
        requested_s = REGEV_HUMAN_S_GENES
        requested_g2m = REGEV_HUMAN_G2M_GENES
        resource_sha256 = REGEV_HUMAN_RESOURCE_SHA256
        resource_name = "Regev/Seurat human cell-cycle programs (43 S + 54 G2/M)"
        identifier_namespace = "human_gene_symbol"
        organism_compatibility = "matched_human" if organism == "human" else "nonhuman_expert_override"
    elif gene_set_source == "custom":
        requested_s = _parse_custom_genes(s_genes, label="Custom S-phase genes")
        requested_g2m = _parse_custom_genes(g2m_genes, label="Custom G2/M-phase genes")
        resource_sha256 = _resource_sha256(requested_s, requested_g2m)
        resource_name = "user-supplied exact cell-cycle programs"
        identifier_namespace = "user_declared_exact"
        organism_compatibility = "user_declared_custom"
    else:
        raise ValueError(f"Unsupported Cell Cycle Score gene-set source: {gene_set_source!r}.")
    overlap = sorted(set(requested_s) & set(requested_g2m))

    available = set(var_names)
    found_s = [gene for gene in requested_s if gene in available]
    found_g2m = [gene for gene in requested_g2m if gene in available]
    missing_s = [gene for gene in requested_s if gene not in available]
    missing_g2m = [gene for gene in requested_g2m if gene not in available]
    if not found_s or not found_g2m:
        raise ValueError(
            "Cell Cycle Score requires at least one observed gene in each phase program; "
            f"found S={len(found_s)}, G2/M={len(found_g2m)}."
        )

    output_prefix = _strict_text(output_prefix, label="Cell Cycle Score output prefix")
    if any(character.isspace() for character in output_prefix):
        raise ValueError("Cell Cycle Score output prefix cannot contain whitespace.")
    output_keys = {
        "s_score": f"{output_prefix}_s_score",
        "g2m_score": f"{output_prefix}_g2m_score",
        "phase": f"{output_prefix}_phase",
    }
    if not isinstance(overwrite_existing, bool):
        raise TypeError("overwrite_existing must be a Boolean.")
    if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 0 <= random_seed <= 2**31 - 1:
        raise ValueError("random_seed must be an integer in 0..2^31-1.")
    collisions = [key for key in output_keys.values() if key in adata.obs]
    if collisions and not overwrite_existing:
        raise ValueError(f"Cell Cycle Score output columns already exist: {collisions}.")

    private_obs = science.pd.DataFrame(index=science.pd.Index(obs_names))
    private_var = science.pd.DataFrame(index=science.pd.Index(var_names))
    work = science.ad.AnnData(matrix.copy(), obs=private_obs, var=private_var)
    with _CELL_CYCLE_RNG_LOCK:
        numpy_state = science.np.random.get_state()
        try:
            returned = science.sc.tl.score_genes_cell_cycle(
                work,
                s_genes=found_s,
                g2m_genes=found_g2m,
                gene_pool=var_names,
                n_bins=25,
                ctrl_as_ref=True,
                use_raw=False,
                random_state=random_seed,
                copy=False,
            )
        finally:
            science.np.random.set_state(numpy_state)
    if returned is not None:
        raise RuntimeError("Scanpy score_genes_cell_cycle(copy=False) unexpectedly returned a value.")
    if any(key not in work.obs for key in ("S_score", "G2M_score", "phase")):
        raise RuntimeError("Scanpy did not return the complete cell-cycle score bundle.")
    s_values = science.pd.to_numeric(work.obs["S_score"], errors="raise").to_numpy(dtype=float)
    g2m_values = science.pd.to_numeric(work.obs["G2M_score"], errors="raise").to_numpy(dtype=float)
    phases = work.obs["phase"].astype(str).tolist()
    if (
        s_values.shape != (len(obs_names),)
        or g2m_values.shape != (len(obs_names),)
        or not bool(science.np.isfinite(s_values).all())
        or not bool(science.np.isfinite(g2m_values).all())
        or any(phase not in {"G1", "S", "G2M"} for phase in phases)
    ):
        raise RuntimeError("Scanpy returned malformed cell-cycle scores or phase labels.")

    output = adata.copy()
    output.obs[output_keys["s_score"]] = s_values
    output.obs[output_keys["g2m_score"]] = g2m_values
    output.obs[output_keys["phase"]] = science.pd.Categorical(phases, categories=["G1", "S", "G2M"])
    phase_counts = {phase: int(phases.count(phase)) for phase in ("G1", "S", "G2M")}
    phase_proportions = {phase: count / len(phases) for phase, count in phase_counts.items()}
    control_size = min(len(found_s), len(found_g2m))
    warnings = [
        *expression_advisories,
        "Cell-cycle scores and phase labels are supervised per-cell annotations, not a continuous biological clock.",
        "Cells are not biological replicates; this node performs no Sample-level Condition inference or Technical-batch adjustment.",
        "Regressing cell-cycle scores is a separate modeling decision and was not performed.",
    ]
    if missing_s or missing_g2m:
        warnings.append("Some requested phase genes were absent; scoring used only the explicitly reported observed genes.")
    if gene_set_source == "regev_human_97" and (
        len(found_s) < 20
        or len(found_s) * 2 < len(requested_s)
        or len(found_g2m) < 20
        or len(found_g2m) * 2 < len(requested_g2m)
    ):
        warnings.append(
            "Bundled-program coverage is below the documented practice target of at least 20 and half of each phase; the result is low-coverage exploratory scoring."
        )
    if gene_set_source == "custom" and (len(found_s) < 10 or len(found_g2m) < 10):
        warnings.append(
            "A custom phase program has fewer than 10 observed genes; the result is low-coverage exploratory scoring."
        )
    if overlap:
        warnings.append(
            "The S and G2/M programs overlap; shared genes contribute to both scores and weaken phase-specific interpretation."
        )
    if gene_set_source == "regev_human_97" and organism != "human":
        warnings.append(
            "The literal bundled human-symbol program was applied under a nonhuman organism declaration as an expert override; no orthology mapping or species validity is implied."
        )
    parameters = {
        "source_kind": source_kind,
        "layer_name": layer_name,
        "gene_set_source": gene_set_source,
        "organism": organism,
        "output_prefix": output_prefix,
        "overwrite_existing": overwrite_existing,
        "random_seed": random_seed,
        "n_bins": 25,
        "ctrl_as_ref": True,
        "control_size": control_size,
    }
    results = (
        f"Scored {len(obs_names):,} cells using {len(found_s)}/{len(requested_s)} observed S-phase genes and "
        f"{len(found_g2m)}/{len(requested_g2m)} observed G2/M genes from {resource_name}; phase calls were "
        + ", ".join(f"{phase}={phase_counts[phase]:,}" for phase in ("G1", "S", "G2M"))
        + "."
    )
    summary = {
        "schema_version": CELL_CYCLE_SCHEMA,
        "node_id": CELL_CYCLE_NODE_ID,
        "status": "descriptive_cell_cycle_program_score",
        "methods": (
            "Scanpy score_genes_cell_cycle calculated the mean observed phase-program expression minus "
            "expression-bin-matched control genes (25 bins; Scanpy 1.12 ctrl_as_ref policy), then assigned G1, "
            "S, or G2M from the two scores."
        ),
        "results": results,
        "key_results": {
            "input_cells": len(obs_names),
            "input_features": len(current_var_names),
            "selected_source_features": len(var_names),
            "expression_source": source_label,
            "expression_state": source_state,
            "expression_state_evidence": state_evidence,
            "full_gene_status": "user_declared_not_programmatically_verified",
            "expression_fingerprint_sha256": _matrix_fingerprint(
                matrix, obs_names=obs_names, var_names=var_names, science=science
            ),
            "resource": {
                "name": resource_name,
                "sha256": resource_sha256,
                "organism": organism,
                "organism_compatibility": organism_compatibility,
                "gene_identifier_namespace": identifier_namespace,
                "requested_s_genes": list(requested_s),
                "requested_g2m_genes": list(requested_g2m),
                "found_s_genes": found_s,
                "found_g2m_genes": found_g2m,
                "missing_s_genes": missing_s,
                "missing_g2m_genes": missing_g2m,
                "cross_phase_overlap": overlap,
                "s_coverage": len(found_s) / len(requested_s),
                "g2m_coverage": len(found_g2m) / len(requested_g2m),
            },
            "score_distributions": {
                "s_score": summarize_numeric(s_values),
                "g2m_score": summarize_numeric(g2m_values),
            },
            "phase_counts": phase_counts,
            "phase_proportions": phase_proportions,
            "output_keys": output_keys,
        },
        "parameters": parameters,
        "references": CELL_CYCLE_REFERENCES,
        "software_versions": collect_software_versions(("scanpy", "anndata", "numpy", "pandas", "scipy")),
        "warnings": warnings,
        "limitations": [
            "The gene programs are supervised signatures and may not represent every organism, tissue, disease state, or assay.",
            "Scores depend on the user-selected expression state, observed gene coverage, gene pool, binning, and seed; this node does not certify normalization or full-gene completeness.",
            "Phase is not synchronized time, lineage, proliferation rate, or a hypothesis test.",
        ],
    }
    json.dumps(summary, allow_nan=False, ensure_ascii=False)
    return output, summary


def cell_cycle_code(
    *,
    source_kind: str,
    layer_name: str | None,
    gene_set_source: str,
    organism: str,
    s_genes: Any,
    g2m_genes: Any,
    output_prefix: str,
    overwrite_existing: bool,
    random_seed: int,
) -> str:
    helpers = (
        _resource_sha256,
        _strict_text,
        _parse_custom_genes,
        _canonical_axes,
        _matrix_fingerprint,
        _validate_expression,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    implementation = dedent(inspect.getsource(analyze_cell_cycle_score)).strip()
    dependency_import = "    from . import dependencies\n\n    science = dependencies.require_scientific_dependencies()"
    if dependency_import not in implementation:
        raise RuntimeError("Cell-cycle source extraction could not locate its dependency seam.")
    implementation = implementation.replace(
        dependency_import,
        "    science = _standalone_cell_cycle_science_dependencies()",
        1,
    )
    return f'''from __future__ import annotations

import hashlib
import json
import platform
import threading
from collections.abc import Sequence
from importlib import metadata as importlib_metadata
from types import SimpleNamespace
from typing import Any

CELL_CYCLE_NODE_ID = {CELL_CYCLE_NODE_ID!r}
CELL_CYCLE_SCHEMA = {CELL_CYCLE_SCHEMA!r}
REGEV_HUMAN_S_GENES = {REGEV_HUMAN_S_GENES!r}
REGEV_HUMAN_G2M_GENES = {REGEV_HUMAN_G2M_GENES!r}
REGEV_HUMAN_RESOURCE_SHA256 = {REGEV_HUMAN_RESOURCE_SHA256!r}
CELL_CYCLE_REFERENCES = {CELL_CYCLE_REFERENCES!r}
_CELL_CYCLE_RNG_LOCK = threading.RLock()


def _standalone_cell_cycle_science_dependencies():
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    return SimpleNamespace(ad=ad, np=np, pd=pd, sc=sc, sparse=sparse)


def collect_software_versions(packages):
    versions = {{
        "python": platform.python_version(),
        "openbio-singlecell": {PLUGIN_VERSION!r},
    }}
    for package in dict.fromkeys(packages):
        if package in versions:
            continue
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def summarize_numeric(values):
    import numpy as np

    array = np.asarray(values, dtype=float).ravel()
    finite = array[np.isfinite(array)]
    missing = int(array.size - finite.size)
    if finite.size == 0:
        return {{
            "n": int(array.size),
            "missing": missing,
            "min": None,
            "q1": None,
            "median": None,
            "mean": None,
            "q3": None,
            "max": None,
        }}
    quantiles = np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {{
        "n": int(array.size),
        "missing": missing,
        "min": float(quantiles[0]),
        "q1": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(finite.mean()),
        "q3": float(quantiles[3]),
        "max": float(quantiles[4]),
    }}


{helper_source}


{implementation}


def score_cell_cycle(adata):
    """Return (annotated_adata, strict_summary) without importing OpenBio."""
    return analyze_cell_cycle_score(
        adata,
        source_kind={source_kind!r},
        layer_name={layer_name!r},
        gene_set_source={gene_set_source!r},
        organism={organism!r},
        s_genes={s_genes!r},
        g2m_genes={g2m_genes!r},
        output_prefix={output_prefix!r},
        overwrite_existing={overwrite_existing!r},
        random_seed={random_seed!r},
    )
'''


__all__ = [
    "CELL_CYCLE_NODE_ID",
    "CELL_CYCLE_REFERENCES",
    "CELL_CYCLE_SCHEMA",
    "REGEV_HUMAN_G2M_GENES",
    "REGEV_HUMAN_RESOURCE_SHA256",
    "REGEV_HUMAN_S_GENES",
    "analyze_cell_cycle_score",
    "cell_cycle_code",
]
