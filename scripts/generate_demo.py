from __future__ import annotations

import argparse
import os
from collections.abc import Sequence
from pathlib import Path

CELL_COUNT = 600
GENE_COUNT = 500
RANDOM_SEED = 42
METADATA_SCHEMA_VERSION = 1

MT_GENES = [
    "MT-ND1",
    "MT-ND2",
    "MT-CO1",
    "MT-CO2",
    "MT-ATP6",
    "MT-ATP8",
    "MT-CYB",
    "MT-RNR2",
]

KNOWN_MARKERS = {
    "T cell": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB", "LCK", "MAL", "TCF7", "CCR7"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA", "CD22", "BANK1", "CD19", "CD83"],
    "Monocyte": ["LYZ", "S100A8", "S100A9", "FCN1", "CTSS", "TYMP", "LGALS3", "VCAN", "MNDA", "CTSD"],
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate the local openbio-singlecell demonstration AnnData.")
    parser.add_argument(
        "--comfy-root",
        type=Path,
        help="ComfyUI root. Defaults to OPENBIO_COMFYUI_ROOT, the Manager install root, or a sibling ComfyUI repository.",
    )
    parser.add_argument("--force", action="store_true", help="Replace an existing demo file, including an invalid one.")
    return parser.parse_args(argv)


def resolve_comfy_root(requested: Path | None) -> Path:
    plugin_root = Path(__file__).resolve().parents[1]
    configured = os.environ.get("OPENBIO_COMFYUI_ROOT")
    manager_configured = os.environ.get("COMFYUI_FOLDERS_BASE_PATH")

    if requested is not None:
        candidates = [requested]
        source = "--comfy-root"
    elif configured:
        candidates = [Path(configured)]
        source = "OPENBIO_COMFYUI_ROOT"
    elif manager_configured:
        candidates = [Path(manager_configured)]
        source = "COMFYUI_FOLDERS_BASE_PATH"
    else:
        candidates = []
        if plugin_root.parent.name.lower() == "custom_nodes":
            candidates.append(plugin_root.parent.parent)
        candidates.append(plugin_root.parent / "ComfyUI")
        source = None

    for candidate in candidates:
        candidate = candidate.expanduser()
        if (candidate / "main.py").is_file():
            return candidate.resolve()

    if source is not None:
        raise SystemExit(f"{source} does not point to a ComfyUI root containing main.py: {candidates[0]}")
    raise SystemExit("ComfyUI root was not found. Pass --comfy-root or set OPENBIO_COMFYUI_ROOT.")


def load_dependencies():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy import sparse
    except ImportError as error:
        raise SystemExit("Demo dependencies are missing. Run the openbio-singlecell install script first.") from error
    return ad, np, pd, sparse


def validate_existing(path: Path, ad, sparse) -> tuple[bool, str]:
    try:
        adata = ad.read_h5ad(path)
    except Exception as error:
        return False, str(error)

    required_obs = {"cell_type", "sample", "condition", "batch"}
    expected_markers = {gene for markers in KNOWN_MARKERS.values() for gene in markers}
    expected_group_counts = {group: CELL_COUNT // len(KNOWN_MARKERS) for group in KNOWN_MARKERS}
    metadata = adata.uns.get("openbio_singlecell", {})
    if not isinstance(metadata, dict):
        return False, "OpenBio metadata must be a mapping"
    if metadata.get("schema_version") != METADATA_SCHEMA_VERSION:
        return (
            False,
            "OpenBio metadata has unsupported schema_version "
            f"{metadata.get('schema_version')!r}; expected {METADATA_SCHEMA_VERSION}",
        )
    if not isinstance(metadata.get("analysis_history"), dict):
        return False, "OpenBio metadata analysis_history must be a mapping"
    cell_groups = set(adata.obs["cell_type"].astype(str)) if "cell_type" in adata.obs else set()
    group_counts = (
        adata.obs["cell_type"].astype(str).value_counts().to_dict() if required_obs.issubset(adata.obs.columns) else {}
    )
    sample_names = set(adata.obs["sample"].astype(str)) if "sample" in adata.obs else set()
    condition_names = set(adata.obs["condition"].astype(str)) if "condition" in adata.obs else set()
    batch_names = set(adata.obs["batch"].astype(str)) if "batch" in adata.obs else set()
    integer_counts = sparse.isspmatrix_csr(adata.X) and adata.X.dtype.kind in "iu"
    nonnegative_counts = integer_counts and (adata.X.data.size == 0 or bool((adata.X.data >= 0).all()))
    checks = [
        (adata.shape == (CELL_COUNT, GENE_COUNT), f"unexpected shape {adata.shape}"),
        (sparse.isspmatrix_csr(adata.X), "X is not CSR sparse"),
        (required_obs.issubset(adata.obs.columns), "required obs columns are missing"),
        (cell_groups == set(KNOWN_MARKERS), "cell groups do not match"),
        (group_counts == expected_group_counts, "cell groups do not contain exactly 200 cells each"),
        (sample_names == {"sample_1", "sample_2", "sample_3", "sample_4"}, "samples do not match"),
        (condition_names == {"control", "treated"}, "conditions do not match"),
        (batch_names == {"batch_1", "batch_2"}, "batches do not match"),
        (expected_markers.issubset(set(adata.var_names.astype(str))), "known marker genes are missing"),
        (sum(name.startswith("MT-") for name in adata.var_names) == len(MT_GENES), "MT genes do not match"),
        (integer_counts, "X does not contain integer counts"),
        (nonnegative_counts, "X contains negative counts"),
        (metadata.get("display_name") == "OpenBio single-cell demo", "OpenBio metadata is missing"),
    ]
    for valid, reason in checks:
        if not valid:
            return False, reason
    return True, "valid"


def build_demo(np, pd, sparse, ad):
    rng = np.random.default_rng(RANDOM_SEED)
    marker_genes = [gene for markers in KNOWN_MARKERS.values() for gene in markers]
    other_genes = [f"GENE_{index:04d}" for index in range(1, GENE_COUNT - len(MT_GENES) - len(marker_genes) + 1)]
    gene_names = MT_GENES + marker_genes + other_genes

    groups = np.repeat(np.array(list(KNOWN_MARKERS), dtype=object), CELL_COUNT // len(KNOWN_MARKERS))
    rng.shuffle(groups)
    samples = rng.choice(["sample_1", "sample_2", "sample_3", "sample_4"], size=CELL_COUNT)
    batches = np.where(np.isin(samples, ["sample_1", "sample_2"]), "batch_1", "batch_2")
    conditions = np.where(np.isin(samples, ["sample_1", "sample_3"]), "control", "treated")

    counts = rng.negative_binomial(1, 0.8, size=(CELL_COUNT, GENE_COUNT)).astype(np.int32)
    gene_index = {gene: index for index, gene in enumerate(gene_names)}
    for group, markers in KNOWN_MARKERS.items():
        cell_indices = np.flatnonzero(groups == group)
        marker_indices = [gene_index[gene] for gene in markers]
        counts[np.ix_(cell_indices, marker_indices)] += rng.poisson(
            4.0, size=(len(cell_indices), len(marker_indices))
        ).astype(np.int32)

    high_mt_cells = rng.choice(CELL_COUNT, size=30, replace=False)
    mt_indices = [gene_index[gene] for gene in MT_GENES]
    counts[np.ix_(high_mt_cells, mt_indices)] += rng.poisson(5.0, size=(len(high_mt_cells), len(mt_indices))).astype(
        np.int32
    )

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(groups, categories=list(KNOWN_MARKERS)),
            "sample": pd.Categorical(samples),
            "condition": pd.Categorical(conditions, categories=["control", "treated"]),
            "batch": pd.Categorical(batches),
        },
        index=[f"cell_{index:04d}" for index in range(1, CELL_COUNT + 1)],
    )
    marker_for = {gene: group for group, markers in KNOWN_MARKERS.items() for gene in markers}
    var = pd.DataFrame(
        {
            "feature_types": "Gene Expression",
            "mt": [gene in MT_GENES for gene in gene_names],
            "known_marker_for": [marker_for.get(gene, "") for gene in gene_names],
        },
        index=gene_names,
    )

    adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.uns["openbio_singlecell"] = {
        "schema_version": METADATA_SCHEMA_VERSION,
        "version": "0.2.0",
        "display_name": "OpenBio single-cell demo",
        "source": {"kind": "generated_demo"},
        "random_seed": RANDOM_SEED,
        "warnings": [],
        "analysis_history": {},
        "known_markers": KNOWN_MARKERS,
    }
    return adata


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    comfy_root = resolve_comfy_root(args.comfy_root)
    output = comfy_root / "input" / "openbio-singlecell" / "openbio_singlecell_demo.h5ad"
    ad, np, pd, sparse = load_dependencies()

    if output.exists() and not args.force:
        valid, reason = validate_existing(output, ad, sparse)
        if valid:
            print(f"Demo AnnData already exists and is valid: {output}")
            return 0
        raise SystemExit(f"Existing demo AnnData is invalid ({reason}). Re-run with --force to replace it: {output}")

    output.parent.mkdir(parents=True, exist_ok=True)
    adata = build_demo(np, pd, sparse, ad)
    adata.write_h5ad(output, compression="gzip")
    print(f"Generated {adata.n_obs} cells x {adata.n_vars} genes: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
