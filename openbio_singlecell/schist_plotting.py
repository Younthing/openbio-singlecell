from __future__ import annotations

import inspect
import textwrap

from .schist_analysis import _plain_json, _standalone_hierarchy_evidence_fingerprint


def _standalone_schist_hierarchy_plot(adata, *, key_added="nsbm", _return_details=False):
    import io
    from collections.abc import Mapping

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import LineCollection
    from matplotlib.figure import Figure

    operation = "Schist Hierarchy Plot"
    if not hasattr(adata, "obs") or not hasattr(adata, "obsm") or not hasattr(adata, "uns"):
        raise TypeError(f"{operation} requires an AnnData-like input.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if int(adata.n_obs) < 1 or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires uniquely identified observations.")
    if not isinstance(key_added, str) or not key_added.strip():
        raise ValueError(f"{operation} key_added must be a nonblank string.")
    key_added = key_added.strip()

    openbio = adata.uns.get("openbio_singlecell")
    history = openbio.get("analysis_history") if isinstance(openbio, Mapping) else None
    matching_history = []
    if isinstance(history, Mapping):
        for entry in history.values():
            parameters = entry.get("parameters") if isinstance(entry, Mapping) else None
            if (
                isinstance(parameters, Mapping)
                and entry.get("operation") == "schist_nested_sbm_hierarchy"
                and parameters.get("key_added") == key_added
            ):
                matching_history.append(entry)
    if not matching_history:
        raise ValueError(
            f"{operation} requires matching OpenBio Schist producer history for key_added={key_added!r}."
        )
    producer = matching_history[-1]
    if producer.get("input_cells") != int(adata.n_obs) or producer.get("input_genes") != int(adata.n_vars):
        raise ValueError(f"{operation} current axes do not match the recorded Schist producer axes.")

    schist_root = adata.uns.get("schist")
    bundle = schist_root.get(key_added) if isinstance(schist_root, Mapping) else None
    if not isinstance(bundle, Mapping):
        raise ValueError(f"{operation} requires adata.uns['schist'][{key_added!r}].")
    if not all(isinstance(bundle.get(name), Mapping) for name in ("stats", "blocks", "params")):
        raise ValueError(f"{operation} requires complete stored stats, blocks, and params mappings.")
    stats, blocks, params = bundle["stats"], bundle["blocks"], bundle["params"]
    if bundle.get("openbio_evidence_schema") != "openbio-singlecell/schist-hierarchy-evidence/v1":
        raise ValueError(f"{operation} requires the current stored hierarchy evidence schema.")
    stored_evidence_sha256 = bundle.get("openbio_evidence_sha256")
    if not isinstance(stored_evidence_sha256, str) or len(stored_evidence_sha256) != 64:
        raise ValueError(f"{operation} requires the stored hierarchy evidence fingerprint.")
    if params.get("nested") is not True or params.get("collect_marginals") is not True:
        raise ValueError(f"{operation} requires a nested Schist result with collected marginals.")
    if params.get("key_added") != key_added:
        raise ValueError(f"{operation} stored key_added provenance is inconsistent.")

    try:
        levels = sorted(int(key) for key in blocks)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{operation} stored block keys must be integer strings.") from error
    if any(str(level) not in blocks for level in levels) or levels != list(range(len(levels))) or not levels:
        raise ValueError(f"{operation} requires consecutive string-keyed hierarchy levels starting at zero.")

    memberships = []
    marginal_ranges = []
    cluster_counts = []
    cluster_sizes = []
    parent_maps = []
    for level in levels:
        obs_key = f"{key_added}_level_{level}"
        obsm_key = f"CM_{key_added}_level_{level}"
        if obs_key not in adata.obs or obsm_key not in adata.obsm:
            raise ValueError(f"{operation} requires stored membership {obs_key!r} and marginal {obsm_key!r}.")
        try:
            membership = np.asarray([int(str(value)) for value in adata.obs[obs_key]], dtype=np.int64)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{operation} membership {obs_key!r} must contain integer block labels.") from error
        if membership.shape != (int(adata.n_obs),) or bool((membership < 0).any()):
            raise ValueError(f"{operation} membership {obs_key!r} is not observation-aligned nonnegative labels.")
        unique = np.unique(membership)
        if not np.array_equal(unique, np.arange(unique.size, dtype=np.int64)):
            raise ValueError(f"{operation} membership {obs_key!r} must use consecutive zero-based labels.")
        count = int(unique.size)
        sizes = np.bincount(membership, minlength=count).astype(np.int64)
        marginal = np.asarray(adata.obsm[obsm_key])
        if marginal.shape != (int(adata.n_obs), count):
            raise ValueError(f"{operation} marginal {obsm_key!r} does not match observations and stored blocks.")
        if (
            not np.issubdtype(marginal.dtype, np.number)
            or np.issubdtype(marginal.dtype, np.bool_)
            or np.iscomplexobj(marginal)
        ):
            raise TypeError(f"{operation} marginal {obsm_key!r} must contain real numeric probabilities.")
        marginal = np.asarray(marginal, dtype=float)
        row_sums = marginal.sum(axis=1)
        if (
            not bool(np.isfinite(marginal).all())
            or bool(((marginal < 0) | (marginal > 1)).any())
            or bool((row_sums <= 0).any())
            or bool((row_sums > 1.0 + 1e-8).any())
        ):
            raise ValueError(f"{operation} marginal {obsm_key!r} contains invalid probability mass.")
        memberships.append(membership)
        cluster_counts.append(count)
        cluster_sizes.append(sizes)
        marginal_ranges.append([float(row_sums.min()), float(row_sums.max())])

        block_array = np.asarray(blocks[str(level)])
        expected_length = int(adata.n_obs) if level == 0 else cluster_counts[level - 1]
        if (
            block_array.ndim != 1
            or int(block_array.size) != expected_length
            or not np.issubdtype(block_array.dtype, np.integer)
        ):
            raise ValueError(f"{operation} stored block array for level {level} has invalid shape or dtype.")
        block_array = block_array.astype(np.int64, copy=False)
        if level == 0:
            if not np.array_equal(block_array, membership):
                raise ValueError(f"{operation} finest stored block array does not match cell membership.")
        else:
            if bool(((block_array < 0) | (block_array >= count)).any()):
                raise ValueError(f"{operation} stored parent labels at level {level} are out of range.")
            if not np.array_equal(block_array[memberships[level - 1]], membership):
                raise ValueError(f"{operation} nested parent mapping is inconsistent at level {level}.")
            parent_maps.append(block_array)

    if any(right > left for left, right in zip(cluster_counts, cluster_counts[1:], strict=False)):
        raise ValueError(f"{operation} block counts must be nonincreasing from finest to root.")
    if cluster_counts[-1] != 1:
        raise ValueError(f"{operation} hierarchy must end in one root block.")

    try:
        total_entropy = float(stats.get("entropy"))
        modularity = np.asarray(stats.get("modularity"), dtype=float).ravel()
        level_entropy = np.asarray(stats.get("level_entropy"), dtype=float).ravel()
    except (TypeError, ValueError) as error:
        raise ValueError(f"{operation} stored hierarchy statistics are invalid.") from error
    level_count = len(levels)
    if not np.isfinite(total_entropy) or total_entropy < 0:
        raise ValueError(f"{operation} total entropy must be finite and nonnegative.")
    if modularity.shape != (level_count,) or not bool(np.isfinite(modularity).all()):
        raise ValueError(f"{operation} modularity must be finite and aligned to hierarchy levels.")
    if bool(((modularity < -1) | (modularity > 1)).any()):
        raise ValueError(f"{operation} modularity values must lie within [-1, 1].")
    if level_entropy.shape != (level_count,) or not bool(np.isfinite(level_entropy).all()):
        raise ValueError(f"{operation} level entropy must be finite and aligned to hierarchy levels.")
    observed_evidence_sha256 = _standalone_hierarchy_evidence_fingerprint(
        adata,
        key_added=key_added,
        np=np,
    )
    if observed_evidence_sha256 != stored_evidence_sha256:
        raise ValueError(f"{operation} stored hierarchy evidence fingerprint does not match current content.")

    figure = Figure(figsize=(15.0, 4.5), constrained_layout=True)
    FigureCanvasAgg(figure)
    axes = figure.subplots(1, 4, gridspec_kw={"width_ratios": [1.0, 1.0, 1.0, 2.2]})
    x_values = np.asarray(levels, dtype=int)
    axes[0].plot(x_values, cluster_counts, marker="o", color="#3274a1")
    axes[0].set_title("Blocks by level")
    axes[0].set_ylabel("Blocks")
    axes[1].plot(x_values, modularity, marker="o", color="#55a868")
    axes[1].axhline(0.0, color="#777777", linewidth=0.7)
    axes[1].set_title("Modularity")
    axes[1].set_ylabel("Modularity")
    axes[2].plot(x_values, level_entropy, marker="o", color="#c44e52")
    axes[2].set_title("Level entropy")
    axes[2].set_ylabel("Entropy")
    for axis in axes[:3]:
        axis.set_xlabel("Hierarchy level")
        axis.set_xticks(x_values)
        axis.grid(alpha=0.2)

    total_blocks = int(sum(cluster_counts))
    render_limit = 500
    selected = []
    if total_blocks <= render_limit:
        selected = [np.arange(count, dtype=np.int64) for count in cluster_counts]
    else:
        quota = max(1, render_limit // level_count)
        for sizes in cluster_sizes:
            selected.append(np.argsort(-sizes, kind="stable")[:quota])
    positions = []
    for level, chosen in enumerate(selected):
        full_count = cluster_counts[level]
        y_all = np.linspace(0.0, 1.0, full_count) if full_count > 1 else np.asarray([0.5])
        positions.append({int(block): float(y_all[block]) for block in chosen})
    segments = []
    for level, parent_map in enumerate(parent_maps):
        for child in selected[level]:
            parent = int(parent_map[int(child)])
            if parent in positions[level + 1]:
                segments.append(
                    [(float(level), positions[level][int(child)]), (float(level + 1), positions[level + 1][parent])]
                )
    if segments:
        axes[3].add_collection(LineCollection(segments, colors="#9a9a9a", linewidths=0.8, alpha=0.55))
    palette = ("#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#64b5cd")
    for level, chosen in enumerate(selected):
        y_values = [positions[level][int(block)] for block in chosen]
        point_sizes = [20.0 + 280.0 * cluster_sizes[level][int(block)] / int(adata.n_obs) for block in chosen]
        axes[3].scatter(
            np.full(len(chosen), level),
            y_values,
            s=point_sizes,
            color=palette[level % len(palette)],
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )
    axes[3].set_xlim(-0.2, level_count - 0.8 if level_count > 1 else 0.8)
    axes[3].set_ylim(-0.05, 1.05)
    axes[3].set_xticks(x_values)
    axes[3].set_yticks([])
    axes[3].set_xlabel("Hierarchy level")
    axes[3].set_title("Nested block hierarchy")
    axes[3].grid(axis="x", alpha=0.15)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")

    rendered_counts = [int(len(chosen)) for chosen in selected]
    details = {
        "key_added": key_added,
        "hierarchy_levels": level_count,
        "cluster_counts_finest_to_root": cluster_counts,
        "cluster_sizes_finest_to_root": [sizes[:50].astype(int).tolist() for sizes in cluster_sizes],
        "cluster_sizes_record_limit": 50,
        "cluster_sizes_truncated": [len(sizes) > 50 for sizes in cluster_sizes],
        "parent_edges": int(sum(array.size for array in parent_maps)),
        "rendered_parent_edges": len(segments),
        "rendered_block_counts": rendered_counts,
        "hierarchy_truncated": rendered_counts != cluster_counts,
        "modularity": modularity.tolist(),
        "level_entropy": level_entropy.tolist(),
        "total_entropy": total_entropy,
        "marginal_row_sum_ranges": marginal_ranges,
        "axis_validation": "matched OpenBio Schist producer dimensions",
        "hierarchy_evidence_sha256": observed_evidence_sha256,
    }
    return (png, details) if _return_details else png


def schist_hierarchy_plot_code(*, key_added: str) -> str:
    implementation = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (_plain_json, _standalone_hierarchy_evidence_fingerprint, _standalone_schist_hierarchy_plot)
    )
    return f"""{implementation}


def plot_schist_hierarchy(adata):
    return _standalone_schist_hierarchy_plot(adata, key_added={key_added!r})
"""


__all__ = ["_standalone_schist_hierarchy_plot", "schist_hierarchy_plot_code"]
