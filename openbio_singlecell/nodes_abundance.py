from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_table_result
from .node_types import AnnDataType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData


ABUNDANCE_CATEGORY = "openbio/single-cell/differential-abundance"


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError(
            "Pertpy differential-abundance nodes require pertpy and the relevant optional extras."
        ) from error
    return pertpy


def _require_schist() -> Any:
    try:
        import schist
    except (ImportError, OSError) as error:
        raise RuntimeError("Schist Nested Model requires the schist package.") from error
    return schist


def _comma_separated_columns(value: str) -> list[str]:
    return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))


def _table_with_index(frame: Any, index_name: str, science: dependencies.ScientificDependencies) -> Any:
    table = science.pd.DataFrame(frame).copy()
    if isinstance(table.index, science.pd.MultiIndex):
        table.index = table.index.set_names(
            [name or f"{index_name}_{position + 1}" for position, name in enumerate(table.index.names)]
        )
    elif table.index.name is None:
        table.index.name = index_name
    return table.reset_index()


def _sample_composition_table(
    adata: AnnData,
    sample_key: str,
    group_key: str,
    annotation_key: str,
) -> tuple[Any, int]:
    science = dependencies.require_scientific_dependencies()
    keys = [sample_key, group_key, annotation_key]
    if len(set(keys)) != len(keys):
        raise ValueError("Composition sample, group, and annotation columns must be different.")
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise ValueError(f"Composition observation columns not found: {missing}")

    frame = adata.obs[keys].copy()
    valid = frame.notna().all(axis=1)
    dropped_cells = int((~valid).sum())
    frame = frame.loc[valid]
    if frame.empty:
        raise ValueError("No observations remain after excluding missing composition metadata.")
    for key in keys:
        frame[key] = frame[key].astype(str)

    groups_per_sample = frame.groupby(sample_key, observed=True)[group_key].nunique()
    ambiguous_samples = groups_per_sample[groups_per_sample > 1]
    if not ambiguous_samples.empty:
        raise ValueError(f"Each sample must map to one group; conflicting samples: {ambiguous_samples.index.tolist()}")

    sample_metadata = frame[[sample_key, group_key]].drop_duplicates(sample_key)
    samples = sample_metadata[sample_key].tolist()
    annotations = list(dict.fromkeys(frame[annotation_key].tolist()))
    index = science.pd.MultiIndex.from_product(
        [samples, annotations],
        names=[sample_key, annotation_key],
    )
    counts = frame.groupby([sample_key, annotation_key], observed=True).size().reindex(index, fill_value=0)
    table = counts.rename("count").reset_index()
    table = table.merge(sample_metadata, on=sample_key, how="left", validate="many_to_one")
    table["count"] = table["count"].astype(int)
    totals = table.groupby(sample_key, observed=True)["count"].transform("sum")
    table["proportion"] = table["count"] / totals
    table = table.rename(
        columns={
            sample_key: "sample",
            group_key: "group",
            annotation_key: "annotation",
        }
    )
    return table[["sample", "group", "annotation", "count", "proportion"]], dropped_cells


def _adjust_bh(table: Any, science: dependencies.ScientificDependencies) -> Any:
    from statsmodels.stats.multitest import multipletests

    table = table.copy()
    table["p_adj"] = science.np.nan
    finite = science.np.isfinite(table["p"].to_numpy(dtype=float))
    if finite.any():
        table.loc[finite, "p_adj"] = multipletests(table.loc[finite, "p"], method="fdr_bh")[1]
    return table


class OpenBioSingleCellSampleCompositionSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSampleCompositionSummary",
            display_name="Sample Composition Summary",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("group_key", default="group"),
                io.String.Input("annotation_key", default="cell_type"),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        group_key: str = "group",
        annotation_key: str = "cell_type",
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, dropped_cells = _sample_composition_table(adata, sample_key, group_key, annotation_key)
        parameters = {
            "sample_key": sample_key,
            "group_key": group_key,
            "annotation_key": annotation_key,
        }
        warnings = []
        if dropped_cells:
            warnings.append(f"Excluded {dropped_cells} cells with missing composition metadata.")
        result = make_table_result(
            title=f"Sample composition by {annotation_key}",
            operation="sample_composition_summary",
            parameters=parameters,
            description="Per-sample cell counts and proportions for each annotation category.",
            warnings=warnings,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellDifferentialCompositionTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDifferentialCompositionTest",
            display_name="Differential Composition Test",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("group_key", default="group"),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("control_group", default=""),
                io.String.Input("comparison_groups", default=""),
                io.Float.Input("pseudocount", default=0.001, min=1e-12, step=0.001, advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        group_key: str = "group",
        annotation_key: str = "cell_type",
        control_group: str = "",
        comparison_groups: str = "",
        pseudocount: float = 0.001,
    ) -> io.NodeOutput:
        from scipy import stats

        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        table, dropped_cells = _sample_composition_table(adata, sample_key, group_key, annotation_key)
        control_group = control_group.strip()
        if not control_group:
            raise ValueError("Differential composition control group cannot be empty.")

        present_groups = list(dict.fromkeys(table["group"].astype(str)))
        if control_group not in present_groups:
            raise ValueError(f"Differential composition control group not found: {control_group!r}")
        comparisons = _comma_separated_columns(comparison_groups)
        if not comparisons:
            comparisons = [group for group in present_groups if group != control_group]
        missing_groups = [group for group in comparisons if group not in present_groups]
        if missing_groups:
            raise ValueError(f"Differential composition comparison groups not found: {missing_groups}")
        comparisons = [group for group in comparisons if group != control_group]
        if not comparisons:
            raise ValueError("Differential composition requires at least one non-control comparison group.")

        selected_groups = [control_group, *comparisons]
        analysis_table = table[table["group"].isin(selected_groups)]
        global_rows = []
        pairwise_rows = []
        for annotation, frame in analysis_table.groupby("annotation", sort=False, observed=True):
            values_by_group = {
                group: frame.loc[frame["group"] == group, "proportion"].to_numpy(dtype=float)
                for group in selected_groups
            }
            try:
                statistic, p_value = stats.kruskal(*(values_by_group[group] for group in selected_groups))
            except ValueError:
                statistic, p_value = science.np.nan, science.np.nan
            global_rows.append(
                {
                    "scope": "global",
                    "annotation": str(annotation),
                    "test": "Kruskal-Wallis",
                    "comparison": "all selected groups",
                    "group": "",
                    "control": "",
                    "n_group": int(sum(len(values_by_group[group]) for group in selected_groups)),
                    "n_control": science.np.nan,
                    "mean_group": science.np.nan,
                    "mean_control": science.np.nan,
                    "difference": science.np.nan,
                    "log2_fold_change": science.np.nan,
                    "statistic": statistic,
                    "p": p_value,
                }
            )

            control_values = values_by_group[control_group]
            control_mean = float(science.np.mean(control_values))
            for group in comparisons:
                group_values = values_by_group[group]
                group_mean = float(science.np.mean(group_values))
                statistic, p_value = stats.mannwhitneyu(group_values, control_values, alternative="two-sided")
                pairwise_rows.append(
                    {
                        "scope": "pairwise",
                        "annotation": str(annotation),
                        "test": "Mann-Whitney U",
                        "comparison": f"{group} vs {control_group}",
                        "group": group,
                        "control": control_group,
                        "n_group": int(len(group_values)),
                        "n_control": int(len(control_values)),
                        "mean_group": group_mean,
                        "mean_control": control_mean,
                        "difference": group_mean - control_mean,
                        "log2_fold_change": float(
                            science.np.log2((group_mean + pseudocount) / (control_mean + pseudocount))
                        ),
                        "statistic": statistic,
                        "p": p_value,
                    }
                )

        global_table = _adjust_bh(science.pd.DataFrame.from_records(global_rows), science)
        pairwise_table = _adjust_bh(science.pd.DataFrame.from_records(pairwise_rows), science)
        result_table = science.pd.concat([global_table, pairwise_table], ignore_index=True)
        parameters = {
            "sample_key": sample_key,
            "group_key": group_key,
            "annotation_key": annotation_key,
            "control_group": control_group,
            "comparison_groups": comparisons,
            "pseudocount": pseudocount,
            "p_adjust": "fdr_bh",
        }
        warnings = []
        if dropped_cells:
            warnings.append(f"Excluded {dropped_cells} cells with missing composition metadata.")
        result = make_table_result(
            title=f"Differential composition by {annotation_key}",
            operation="differential_composition_test",
            parameters=parameters,
            description="Sample-level Kruskal-Wallis tests and pairwise Mann-Whitney tests versus the control group.",
            warnings=warnings,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=result_table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellSchistNestedModel(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSchistNestedModel",
            display_name="Schist Nested Model",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, random_seed: int = 123) -> io.NodeOutput:
        schist = _require_schist()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        schist.inference.nested_model(output, random_seed=random_seed)

        parameters = {"random_seed": random_seed}
        finish_adata(
            output,
            "schist_nested_model",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellMiloDifferentialAbundance(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMiloDifferentialAbundance",
            display_name="Milo Differential Abundance",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("design", default="~ group"),
                io.String.Input("contrast", default=""),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("use_rep", default="X_pca"),
                io.Int.Input("n_neighbors", default=30, min=2, max=2**31 - 1),
                io.Float.Input(
                    "neighborhood_proportion",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "mixed_annotation_threshold",
                    default=0.6,
                    min=0.0,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.String.Input("neighbors_key", default="openbio_milo", advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        design: str = "~ group",
        contrast: str = "",
        annotation_key: str = "cell_type",
        use_rep: str = "X_pca",
        n_neighbors: int = 30,
        neighborhood_proportion: float = 0.1,
        mixed_annotation_threshold: float = 0.6,
        neighbors_key: str = "openbio_milo",
        random_seed: int = 123,
    ) -> io.NodeOutput:
        if sample_key not in adata.obs:
            raise ValueError(f"Milo sample column not found in obs: {sample_key!r}")
        if annotation_key not in adata.obs:
            raise ValueError(f"Milo annotation column not found in obs: {annotation_key!r}")
        if use_rep not in adata.obsm:
            raise ValueError(f"Milo representation not found in obsm: {use_rep!r}")
        if not design.strip():
            raise ValueError("Milo design cannot be empty.")
        if not neighbors_key.strip():
            raise ValueError("Milo neighbors_key cannot be empty.")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        milo = pertpy.tl.Milo()
        mdata = milo.load(adata.copy())
        rna = mdata["rna"]
        science.sc.pp.neighbors(
            rna,
            n_neighbors=n_neighbors,
            use_rep=use_rep,
            key_added=neighbors_key,
            random_state=random_seed,
        )
        milo.make_nhoods(mdata, neighbors_key=neighbors_key, prop=neighborhood_proportion)
        counted = milo.count_nhoods(mdata, sample_col=sample_key)
        if counted is not None:
            mdata = counted

        da_options = {"design": design.strip()}
        if contrast.strip():
            da_options["model_contrasts"] = contrast.strip()
        milo.da_nhoods(mdata, **da_options)
        milo.annotate_nhoods(mdata, anno_col=annotation_key)

        abundance = mdata["milo"].var.copy()
        if "nhood_annotation" in abundance and "nhood_annotation_frac" in abundance:
            mixed = abundance["nhood_annotation_frac"] < mixed_annotation_threshold
            abundance["nhood_annotation"] = abundance["nhood_annotation"].astype(object)
            abundance.loc[mixed, "nhood_annotation"] = "Mixed"
        table = _table_with_index(abundance, "neighborhood", science)

        parameters = {
            "sample_key": sample_key,
            "design": design.strip(),
            "contrast": contrast.strip(),
            "annotation_key": annotation_key,
            "use_rep": use_rep,
            "n_neighbors": n_neighbors,
            "neighborhood_proportion": neighborhood_proportion,
            "mixed_annotation_threshold": mixed_annotation_threshold,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
        result = make_table_result(
            title="Milo differential abundance",
            operation="milo_differential_abundance",
            parameters=parameters,
            description="Differential abundance statistics and majority annotation for each cell neighborhood.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellSccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSccodaDifferentialComposition",
            display_name="scCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("covariate_keys", default="group"),
                io.String.Input("formula", default="group"),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("estimated_fdr", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        annotation_key: str = "cell_type",
        covariate_keys: str = "group",
        formula: str = "group",
        reference_cell_type: str = "automatic",
        estimated_fdr: float = 0.05,
        num_samples: int = 10000,
        num_warmup: int = 1000,
        random_seed: int = 123,
    ) -> io.NodeOutput:
        covariates = _comma_separated_columns(covariate_keys)
        if not covariates:
            raise ValueError("scCODA requires at least one covariate column.")
        required_obs = [sample_key, annotation_key, *covariates]
        missing_obs = [key for key in required_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"scCODA observation columns not found: {missing_obs}")
        if not formula.strip():
            raise ValueError("scCODA formula cannot be empty.")
        reference = reference_cell_type.strip() or "automatic"

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        model = pertpy.tl.Sccoda()
        mdata = model.load(
            adata.copy(),
            type="cell_level",
            generate_sample_level=True,
            cell_type_identifier=annotation_key,
            sample_identifier=sample_key,
            covariate_obs=covariates,
        )
        prepared = model.prepare(
            mdata,
            modality_key="coda",
            formula=formula.strip(),
            reference_cell_type=reference,
        )
        if prepared is not None:
            mdata = prepared
        model.run_nuts(
            mdata,
            modality_key="coda",
            rng_key=random_seed,
            num_samples=num_samples,
            num_warmup=num_warmup,
        )
        model.set_fdr(mdata, estimated_fdr, modality_key="coda")

        summaries = model.summary_prepare(mdata["coda"], est_fdr=estimated_fdr)
        sections = ("intercept", "effect", "node_effect")
        tables = []
        for section, frame in zip(sections, summaries, strict=False):
            if frame is None:
                continue
            table = _table_with_index(frame, "parameter", science)
            table.insert(0, "summary_section", section)
            tables.append(table)
        if not tables:
            raise RuntimeError("scCODA did not return summary tables.")
        table = science.pd.concat(tables, ignore_index=True, sort=False)

        parameters = {
            "sample_key": sample_key,
            "annotation_key": annotation_key,
            "covariate_keys": covariates,
            "formula": formula.strip(),
            "reference_cell_type": reference,
            "estimated_fdr": estimated_fdr,
            "num_samples": num_samples,
            "num_warmup": num_warmup,
            "random_seed": random_seed,
        }
        result = make_table_result(
            title="scCODA differential composition",
            operation="sccoda_differential_composition",
            parameters=parameters,
            description="Intercept and cell-type effect summaries from the scCODA compositional model.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellTasccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellTasccodaDifferentialComposition",
            display_name="tascCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("covariate_keys", default="group"),
                io.String.Input("formula", default="group"),
                io.String.Input(
                    "hierarchy_levels",
                    default="",
                ),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("phi", default=0.0, min=0.0, step=0.1, advanced=True),
                io.Float.Input("lambda_1", default=3.5, min=0.0, step=0.1, advanced=True),
                io.Float.Input("estimated_fdr", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=1234, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        annotation_key: str = "cell_type",
        covariate_keys: str = "group",
        formula: str = "group",
        hierarchy_levels: str = "",
        reference_cell_type: str = "automatic",
        phi: float = 0.0,
        lambda_1: float = 3.5,
        estimated_fdr: float = 0.05,
        num_samples: int = 10000,
        num_warmup: int = 1000,
        random_seed: int = 1234,
    ) -> io.NodeOutput:
        covariates = _comma_separated_columns(covariate_keys)
        levels = _comma_separated_columns(hierarchy_levels)
        required_obs = [sample_key, annotation_key, *covariates, *levels]
        missing_obs = [key for key in required_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"tascCODA observation columns not found: {missing_obs}")
        if not covariates:
            raise ValueError("tascCODA requires at least one covariate column.")
        if not levels:
            raise ValueError("tascCODA requires at least one hierarchy level.")
        if not formula.strip():
            raise ValueError("tascCODA formula cannot be empty.")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        model = pertpy.tl.Tasccoda()
        mdata = model.load(
            adata.copy(),
            type="cell_level",
            cell_type_identifier=annotation_key,
            sample_identifier=sample_key,
            covariate_obs=covariates,
            levels_orig=levels,
            add_level_name=True,
        )
        prepared = model.prepare(
            mdata,
            modality_key="coda",
            reference_cell_type=reference_cell_type or "automatic",
            formula=formula.strip(),
            pen_args={"phi": phi, "lambda_1": lambda_1},
            tree_key="tree",
        )
        if prepared is not None:
            mdata = prepared
        model.run_nuts(
            mdata,
            modality_key="coda",
            rng_key=random_seed,
            num_samples=num_samples,
            num_warmup=num_warmup,
        )

        summaries = model.summary_prepare(mdata["coda"], est_fdr=estimated_fdr)
        sections = ("intercept", "effect", "node_effect")
        tables = []
        for section, frame in zip(sections, summaries, strict=False):
            table = _table_with_index(frame, "parameter", science)
            table.insert(0, "summary_section", section)
            tables.append(table)
        table = science.pd.concat(tables, ignore_index=True, sort=False)

        parameters = {
            "sample_key": sample_key,
            "annotation_key": annotation_key,
            "covariate_keys": covariates,
            "formula": formula.strip(),
            "hierarchy_levels": levels,
            "reference_cell_type": reference_cell_type or "automatic",
            "phi": phi,
            "lambda_1": lambda_1,
            "estimated_fdr": estimated_fdr,
            "num_samples": num_samples,
            "num_warmup": num_warmup,
            "random_seed": random_seed,
        }
        result = make_table_result(
            title="tascCODA differential composition",
            operation="tasccoda_differential_composition",
            parameters=parameters,
            description="Intercept, cell-type effect, and hierarchy-node effect summaries from tascCODA.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


ABUNDANCE_NODE_CLASSES = [
    OpenBioSingleCellSampleCompositionSummary,
    OpenBioSingleCellDifferentialCompositionTest,
    OpenBioSingleCellSchistNestedModel,
    OpenBioSingleCellMiloDifferentialAbundance,
    OpenBioSingleCellSccodaDifferentialComposition,
    OpenBioSingleCellTasccodaDifferentialComposition,
]


__all__ = [
    "ABUNDANCE_NODE_CLASSES",
    "OpenBioSingleCellDifferentialCompositionTest",
    "OpenBioSingleCellMiloDifferentialAbundance",
    "OpenBioSingleCellSampleCompositionSummary",
    "OpenBioSingleCellSchistNestedModel",
    "OpenBioSingleCellSccodaDifferentialComposition",
    "OpenBioSingleCellTasccodaDifferentialComposition",
]
