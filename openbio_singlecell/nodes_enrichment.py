from __future__ import annotations

import importlib
import os
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_result
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData

    from .contracts import SingleCellResult


CATEGORY = "openbio/single-cell/enrichment"
GENE_SET_EXTENSIONS = (".csv", ".tsv", ".gmt")


def _require_optional_dependency(name: str) -> Any:
    try:
        return importlib.import_module(name)
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"The {name!r} package is required for this enrichment node ({exc}).") from exc


def _required_name(value: str, description: str) -> str:
    name = value.strip()
    if not name:
        raise ValueError(f"{description} cannot be empty.")
    return name


def _read_gene_sets(
    relative_path: str,
    source_column: str,
    target_column: str,
    science: dependencies.ScientificDependencies,
) -> Any:
    path = resolve_input_path(relative_path, extensions=GENE_SET_EXTENSIONS)
    if os.path.splitext(path)[1].lower() == ".gmt":
        records: list[dict[str, str]] = []
        with open(path, encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\n\r").split("\t")
                if len(fields) < 3:
                    continue
                records.extend({source_column: fields[0], target_column: gene} for gene in fields[2:] if gene)
        network = science.pd.DataFrame.from_records(records, columns=[source_column, target_column])
    else:
        separator = "\t" if path.lower().endswith(".tsv") else ","
        network = science.pd.read_csv(path, sep=separator)

    missing = [column for column in (source_column, target_column) if column not in network]
    if missing:
        raise ValueError(f"Gene-set file is missing columns: {missing}")
    network = network[[source_column, target_column]].dropna().drop_duplicates().copy()
    network[source_column] = network[source_column].astype(str)
    network[target_column] = network[target_column].astype(str)
    if network.empty:
        raise ValueError("Gene-set file contains no source-target pairs.")
    return network


def _expression_adata(adata: AnnData, source: str, layer_name: str, *, dense: bool = False) -> AnnData:
    science = dependencies.require_scientific_dependencies()
    if source == "raw":
        if adata.raw is None:
            raise ValueError("Expression source 'raw' was selected, but adata.raw is unavailable.")
        work = adata.raw.to_adata()
    else:
        work = adata.copy()
        if source == "layer":
            layer_name = _required_name(layer_name, "Expression layer")
            if layer_name not in adata.layers:
                raise ValueError(f"Expression layer not found: {layer_name!r}")
            work.X = adata.layers[layer_name].copy()

    if science.sparse.issparse(work.X):
        work.X = work.X.toarray() if dense else work.X.tocsr()
    return work


def _copy_decoupler_scores(output: AnnData, work: AnnData, generated_key: str, output_key: str) -> None:
    if generated_key not in work.obsm:
        raise RuntimeError(f"decoupler did not create obsm[{generated_key!r}].")
    output.obsm[_required_name(output_key, "Score output key")] = work.obsm[generated_key].copy()


def _rank_for_enrichment(adata: AnnData, groupby: str, source: str, layer_name: str) -> AnnData:
    science = dependencies.require_scientific_dependencies()
    groupby = _required_name(groupby, "Enrichment groupby column")
    if groupby not in adata.obs:
        raise ValueError(f"Enrichment groupby column not found in obs: {groupby!r}")
    if source == "raw" and adata.raw is None:
        raise ValueError("Expression source 'raw' was selected, but adata.raw is unavailable.")
    if source == "layer":
        layer_name = _required_name(layer_name, "Expression layer")
        if layer_name not in adata.layers:
            raise ValueError(f"Expression layer not found: {layer_name!r}")

    work = adata.copy()
    science.sc.tl.rank_genes_groups(
        work,
        groupby=groupby,
        method="wilcoxon",
        use_raw=source == "raw",
        layer=layer_name if source == "layer" else None,
    )
    return work


def _pertpy_table(values: Any, science: dependencies.ScientificDependencies) -> Any:
    if not isinstance(values, dict) or not values:
        raise RuntimeError("Pertpy enrichment did not return grouped result tables.")

    tables = []
    for group, frame in values.items():
        if not isinstance(frame, science.pd.DataFrame):
            raise RuntimeError(f"Pertpy returned an unsupported result for group {group!r}.")
        table = frame.reset_index()
        table = table.rename(columns={table.columns[0]: "term"})
        table.insert(0, "group", str(group))
        tables.append(table)
    return science.pd.concat(tables, ignore_index=True)


class OpenBioSingleCellAUCellScores(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAUCellScores",
            display_name="AUCell Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="counts"),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.Int.Input("min_n", default=5, min=1, max=2**31 - 1, advanced=True),
                io.String.Input("output_key", default="aucell_estimate", advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        source: str = "raw",
        layer_name: str = "counts",
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        min_n: int = 5,
        output_key: str = "aucell_estimate",
        random_seed: int = 123,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        network = _read_gene_sets(gene_sets_file, source_column, target_column, science)
        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = _expression_adata(output, source, layer_name)
        decoupler.run_aucell(
            mat=work,
            net=network,
            source=source_column,
            target=target_column,
            use_raw=False,
            min_n=min_n,
            seed=random_seed,
            verbose=False,
        )
        _copy_decoupler_scores(output, work, "aucell_estimate", output_key)
        parameters = {
            "gene_sets_file": gene_sets_file,
            "source": source,
            "layer_name": layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "min_n": min_n,
            "output_key": output_key,
            "random_seed": random_seed,
        }
        finish_adata(output, "aucell_scores", parameters, cells, genes, started_at, random_seed=random_seed)
        return io.NodeOutput(output)


class OpenBioSingleCellGSVAScores(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGSVAScores",
            display_name="GSVA Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="counts"),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.Int.Input("min_n", default=10, min=1, max=2**31 - 1, advanced=True),
                io.Boolean.Input("mx_diff", default=True, advanced=True),
                io.Boolean.Input("abs_rnk", default=True, advanced=True),
                io.String.Input("output_key", default="gsva_estimate", advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        source: str = "raw",
        layer_name: str = "counts",
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        min_n: int = 10,
        mx_diff: bool = True,
        abs_rnk: bool = True,
        output_key: str = "gsva_estimate",
        random_seed: int = 123,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        network = _read_gene_sets(gene_sets_file, source_column, target_column, science)
        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = _expression_adata(output, source, layer_name, dense=True)
        decoupler.run_gsva(
            mat=work,
            net=network,
            source=source_column,
            target=target_column,
            use_raw=False,
            min_n=min_n,
            mx_diff=mx_diff,
            abs_rnk=abs_rnk,
            seed=random_seed,
            verbose=False,
        )
        _copy_decoupler_scores(output, work, "gsva_estimate", output_key)
        parameters = {
            "gene_sets_file": gene_sets_file,
            "source": source,
            "layer_name": layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "min_n": min_n,
            "mx_diff": mx_diff,
            "abs_rnk": abs_rnk,
            "output_key": output_key,
            "random_seed": random_seed,
        }
        finish_adata(output, "gsva_scores", parameters, cells, genes, started_at, random_seed=random_seed)
        return io.NodeOutput(output)


class OpenBioSingleCellGenePanelScores(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGenePanelScores",
            display_name="Gene Panel Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="counts"),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.String.Input("output_prefix", default="score_", advanced=True),
                io.Int.Input("ctrl_size", default=50, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_bins", default=25, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        source: str = "raw",
        layer_name: str = "counts",
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        output_prefix: str = "score_",
        ctrl_size: int = 50,
        n_bins: int = 25,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        network = _read_gene_sets(gene_sets_file, source_column, target_column, science)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = _expression_adata(output, source, layer_name)
        warnings: list[str] = []
        score_columns: list[str] = []

        for panel, rows in network.groupby(source_column, sort=False):
            panel_genes = [gene for gene in rows[target_column].astype(str) if gene in work.var_names]
            if not panel_genes:
                warnings.append(f"Gene panel {panel!r} has no genes in the selected expression source.")
                continue
            score_name = f"{output_prefix}{panel}"
            science.sc.tl.score_genes(
                work,
                gene_list=panel_genes,
                score_name=score_name,
                ctrl_size=ctrl_size,
                n_bins=n_bins,
                random_state=random_seed,
                use_raw=False,
            )
            output.obs[score_name] = work.obs[score_name].reindex(output.obs_names)
            score_columns.append(score_name)

        if not score_columns:
            raise ValueError("No gene panel could be scored against the selected expression source.")
        parameters = {
            "gene_sets_file": gene_sets_file,
            "source": source,
            "layer_name": layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "output_prefix": output_prefix,
            "ctrl_size": ctrl_size,
            "n_bins": n_bins,
            "score_columns": score_columns,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "gene_panel_scores",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=warnings,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellPathwayScoreTTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPathwayScoreTTest",
            display_name="Pathway Score T-Test",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("cell_type_column", default="cell_type"),
                io.String.Input("group_column", default="group"),
                io.String.Input("group_a", default=""),
                io.String.Input("group_b", default=""),
                io.String.Input("score_key", default="aucell_estimate", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        cell_type_column: str = "cell_type",
        group_column: str = "group",
        group_a: str = "",
        group_b: str = "",
        score_key: str = "aucell_estimate",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        from scipy.stats import ttest_ind
        from statsmodels.stats.multitest import multipletests

        cell_type_column = _required_name(cell_type_column, "Cell-type column")
        group_column = _required_name(group_column, "Group column")
        group_a = _required_name(group_a, "First group")
        group_b = _required_name(group_b, "Second group")
        if group_a == group_b:
            raise ValueError("Pathway score groups must be different.")
        score_key = _required_name(score_key, "Pathway score key")
        missing_obs = [column for column in (cell_type_column, group_column) if column not in adata.obs]
        if missing_obs:
            raise ValueError(f"Pathway score test observation columns not found: {missing_obs}")
        if score_key not in adata.obsm:
            raise ValueError(f"Pathway scores not found in obsm: {score_key!r}")
        scores = adata.obsm[score_key]
        if not isinstance(scores, science.pd.DataFrame):
            raise ValueError(f"obsm[{score_key!r}] must be a named pathway score DataFrame.")

        started_at = time.perf_counter()
        frame = adata.obs[[cell_type_column, group_column]].join(scores)
        rows = []
        for cell_type in frame[cell_type_column].dropna().unique():
            cell_frame = frame[frame[cell_type_column] == cell_type]
            for pathway in scores.columns:
                values_a = science.pd.to_numeric(
                    cell_frame.loc[cell_frame[group_column].astype(str) == group_a, pathway], errors="coerce"
                )
                values_b = science.pd.to_numeric(
                    cell_frame.loc[cell_frame[group_column].astype(str) == group_b, pathway], errors="coerce"
                )
                statistic, p_value = ttest_ind(values_a, values_b, nan_policy="omit")
                rows.append(
                    {
                        "cell_type": str(cell_type),
                        "pathway": str(pathway),
                        "group_a": group_a,
                        "group_b": group_b,
                        "t_stat": statistic,
                        "p": p_value,
                    }
                )
        table = science.pd.DataFrame.from_records(rows)
        if table.empty:
            raise ValueError("Pathway score test produced no cell-type/pathway comparisons.")
        finite = science.np.isfinite(table["p"].to_numpy(dtype=float))
        table["p_adj"] = science.np.nan
        if finite.any():
            table.loc[finite, "p_adj"] = multipletests(table.loc[finite, "p"], method="fdr_bh")[1]

        parameters = {
            "cell_type_column": cell_type_column,
            "group_column": group_column,
            "group_a": group_a,
            "group_b": group_b,
            "score_key": score_key,
        }
        result = make_result(
            kind="table",
            title=f"{group_a} vs {group_b} pathway scores",
            operation="pathway_score_ttest",
            parameters=parameters,
            description="Per-cell-type independent t-tests with Benjamini-Hochberg correction.",
            warnings=[] if finite.all() else ["Some comparisons produced non-finite p-values."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellRankedGSEA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRankedGSEA",
            display_name="Ranked GSEA",
            category=CATEGORY,
            inputs=[
                SingleCellResultType.Input("marker_result"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.String.Input("group", default=""),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.String.Input("score_column", default="score", advanced=True),
                io.Int.Input("min_genes", default=16, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_genes", default=499, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        marker_result: SingleCellResult,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        group: str = "",
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        score_column: str = "score",
        min_genes: int = 16,
        max_genes: int = 499,
        random_seed: int = 123,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        group = _required_name(group, "Marker group")
        if min_genes > max_genes:
            raise ValueError("GSEA min_genes cannot exceed max_genes.")
        table = marker_result.table
        if not isinstance(table, science.pd.DataFrame):
            raise ValueError("Ranked GSEA requires a Marker Genes table result.")
        required_columns = ["group", "gene", score_column]
        missing = [column for column in required_columns if column not in table]
        if missing:
            raise ValueError(f"Marker result is missing columns required by GSEA: {missing}")
        ranked = table.loc[table["group"].astype(str) == group, ["gene", score_column]].dropna().copy()
        if ranked.empty:
            raise ValueError(f"Marker result contains no rows for group {group!r}.")
        ranked["gene"] = ranked["gene"].astype(str)
        ranked = ranked.drop_duplicates("gene").set_index("gene")[[score_column]].T
        ranked.index = [group]

        network = _read_gene_sets(gene_sets_file, source_column, target_column, science)
        sizes = network.groupby(source_column, observed=True).size()
        selected_sets = sizes.index[(sizes >= min_genes) & (sizes <= max_genes)]
        network = network[network[source_column].isin(selected_sets)]
        if network.empty:
            raise ValueError("No gene sets remain after applying the requested size limits.")

        decoupler = _require_optional_dependency("decoupler")
        started_at = time.perf_counter()
        scores, normalized, p_values = decoupler.run_gsea(
            ranked,
            network,
            source=source_column,
            target=target_column,
            seed=random_seed,
        )
        result_table = science.pd.concat(
            [
                scores.iloc[0].rename("score"),
                normalized.iloc[0].rename("normalized_score"),
                p_values.iloc[0].rename("p"),
            ],
            axis=1,
        )
        result_table = result_table.rename_axis("gene_set").reset_index().sort_values("score", ascending=False)
        parameters = {
            "group": group,
            "gene_sets_file": gene_sets_file,
            "source_column": source_column,
            "target_column": target_column,
            "score_column": score_column,
            "min_genes": min_genes,
            "max_genes": max_genes,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title=f"GSEA for {group}",
            operation="ranked_gsea",
            parameters=parameters,
            description="GSEA applied to a selected group from a Marker Genes result.",
            warnings=[],
            input_cells=marker_result.input_cells,
            input_genes=marker_result.input_genes,
            started_at=started_at,
            random_seed=random_seed,
            table=result_table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellGeneSetOverrepresentation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGeneSetOverrepresentation",
            display_name="Gene Set Overrepresentation",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                SingleCellResultType.Input("result"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.String.Input("group", default=""),
                io.String.Input("selection_column", default="p_adj"),
                io.Combo.Input("selection_operator", options=["<=", ">=", "abs>="], default="<="),
                io.Float.Input("threshold", default=0.05, step=0.01),
                io.Combo.Input("universe_source", options=["X", "raw"], default="X"),
                io.String.Input("gene_column", default="gene", advanced=True),
                io.String.Input("group_column", default="group", advanced=True),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        result: SingleCellResult,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        group: str = "",
        selection_column: str = "p_adj",
        selection_operator: str = "<=",
        threshold: float = 0.05,
        universe_source: str = "X",
        gene_column: str = "gene",
        group_column: str = "group",
        source_column: str = "geneset",
        target_column: str = "genesymbol",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if result.kind != "table" or not isinstance(result.table, science.pd.DataFrame):
            raise ValueError("Gene Set Overrepresentation requires a differential table result.")
        gene_column = _required_name(gene_column, "Differential-result gene column")
        if gene_column not in result.table:
            raise ValueError(f"Differential-result gene column not found: {gene_column!r}")
        table = result.table.copy()
        group = group.strip()
        if group:
            group_column = _required_name(group_column, "Differential-result group column")
            if group_column not in table:
                raise ValueError(f"Differential-result group column not found: {group_column!r}")
            table = table[table[group_column].astype(str) == group]

        selection_column = selection_column.strip()
        if selection_column:
            if selection_column not in table:
                raise ValueError(f"Differential-result selection column not found: {selection_column!r}")
            values = science.pd.to_numeric(table[selection_column], errors="coerce")
            if selection_operator == "<=":
                table = table[values <= threshold]
            elif selection_operator == ">=":
                table = table[values >= threshold]
            else:
                table = table[values.abs() >= threshold]
        selected_genes = set(table[gene_column].dropna().astype(str))
        if not selected_genes:
            raise ValueError("No genes passed the requested differential-result selection.")

        if universe_source == "raw":
            if adata.raw is None:
                raise ValueError("Gene-set universe source 'raw' was selected, but adata.raw is unavailable.")
            universe = set(adata.raw.var_names.astype(str))
        else:
            universe = set(adata.var_names.astype(str))
        selected_genes.intersection_update(universe)
        if not selected_genes:
            raise ValueError("Selected differential genes do not overlap the AnnData gene universe.")

        network = _read_gene_sets(gene_sets_file, source_column, target_column, science)
        from scipy.stats import hypergeom
        from statsmodels.stats.multitest import multipletests

        started_at = time.perf_counter()
        rows = []
        for gene_set, frame in network.groupby(source_column, sort=False):
            gene_set_genes = set(frame[target_column].astype(str)).intersection(universe)
            if not gene_set_genes:
                continue
            overlap = sorted(gene_set_genes.intersection(selected_genes))
            p_value = hypergeom.sf(
                len(overlap) - 1,
                len(universe),
                len(gene_set_genes),
                len(selected_genes),
            )
            rows.append(
                {
                    "gene_set": str(gene_set),
                    "intersection_size": len(overlap),
                    "gene_set_size": len(gene_set_genes),
                    "selected_gene_count": len(selected_genes),
                    "universe_size": len(universe),
                    "overlap_genes": ",".join(overlap),
                    "p": p_value,
                }
            )
        enrichment = science.pd.DataFrame.from_records(rows)
        if enrichment.empty:
            raise ValueError("No gene sets overlap the AnnData gene universe.")
        enrichment["p_adj"] = multipletests(enrichment["p"], method="fdr_bh")[1]
        enrichment = enrichment.sort_values(["p_adj", "p", "gene_set"]).reset_index(drop=True)
        parameters = {
            "gene_sets_file": gene_sets_file,
            "group": group,
            "selection_column": selection_column,
            "selection_operator": selection_operator,
            "threshold": threshold,
            "universe_source": universe_source,
            "gene_column": gene_column,
            "group_column": group_column,
            "source_column": source_column,
            "target_column": target_column,
        }
        enriched = make_result(
            kind="table",
            title=f"Gene-set overrepresentation{f' for {group}' if group else ''}",
            operation="gene_set_overrepresentation",
            parameters=parameters,
            description="Hypergeometric enrichment using AnnData genes as the tested universe.",
            warnings=[]
            if bool((enrichment["intersection_size"] > 0).any())
            else ["No gene set contains a selected gene."],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=enrichment,
        )
        return io.NodeOutput(enriched)


class OpenBioSingleCellDGIdbAnnotation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDGIdbAnnotation",
            display_name="DGIdb Drug Annotation",
            category=CATEGORY,
            inputs=[AnnDataType.Input("adata")],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        pertpy = _require_optional_dependency("pertpy")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        pertpy.md.Drug().annotate(output, source="dgidb")
        finish_adata(output, "dgidb_drug_annotation", {"source": "dgidb"}, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellDrugScores(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugScores",
            display_name="DGIdb Drug Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("source", options=["X", "layer"], default="X"),
                io.String.Input("layer_name", default="counts"),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, source: str = "X", layer_name: str = "counts") -> io.NodeOutput:
        pertpy = _require_optional_dependency("pertpy")
        if source == "layer":
            layer_name = _required_name(layer_name, "Drug-score expression layer")
            if layer_name not in adata.layers:
                raise ValueError(f"Drug-score expression layer not found: {layer_name!r}")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        drug = pertpy.md.Drug()
        pertpy.tl.Enrichment().score(
            output,
            targets=drug.dgidb.dictionary,
            layer=layer_name if source == "layer" else None,
        )
        parameters = {"targets": "dgidb", "source": source, "layer_name": layer_name}
        finish_adata(output, "dgidb_drug_scores", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellDrugHypergeometric(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugHypergeometric",
            display_name="Drug Hypergeometric Enrichment",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="counts"),
                io.Float.Input("marker_padj_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        source: str = "raw",
        layer_name: str = "counts",
        marker_padj_threshold: float = 0.05,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        pertpy = _require_optional_dependency("pertpy")
        started_at = time.perf_counter()
        work = _rank_for_enrichment(adata, groupby, source, layer_name)
        drug = pertpy.md.Drug()
        values = pertpy.tl.Enrichment().hypergeometric(
            work,
            targets=drug.dgidb.dictionary,
            padj_threshold=marker_padj_threshold,
        )
        table = _pertpy_table(values, science)
        parameters = {
            "groupby": groupby,
            "source": source,
            "layer_name": layer_name,
            "marker_padj_threshold": marker_padj_threshold,
            "targets": "dgidb",
            "rank_method": "wilcoxon",
        }
        result = make_result(
            kind="table",
            title=f"Drug overrepresentation by {groupby}",
            operation="drug_hypergeometric",
            parameters=parameters,
            description="DGIdb target overrepresentation from Wilcoxon marker genes.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellDrugGSEA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugGSEA",
            display_name="Drug GSEA",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="counts"),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        source: str = "raw",
        layer_name: str = "counts",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        pertpy = _require_optional_dependency("pertpy")
        started_at = time.perf_counter()
        work = _rank_for_enrichment(adata, groupby, source, layer_name)
        drug = pertpy.md.Drug()
        values = pertpy.tl.Enrichment().gsea(work, targets=drug.dgidb.dictionary)
        table = _pertpy_table(values, science)
        parameters = {
            "groupby": groupby,
            "source": source,
            "layer_name": layer_name,
            "targets": "dgidb",
            "rank_method": "wilcoxon",
        }
        result = make_result(
            kind="table",
            title=f"Drug GSEA by {groupby}",
            operation="drug_gsea",
            parameters=parameters,
            description="DGIdb target GSEA from Wilcoxon marker scores.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


ENRICHMENT_NODE_CLASSES = [
    OpenBioSingleCellAUCellScores,
    OpenBioSingleCellGSVAScores,
    OpenBioSingleCellGenePanelScores,
    OpenBioSingleCellPathwayScoreTTest,
    OpenBioSingleCellRankedGSEA,
    OpenBioSingleCellGeneSetOverrepresentation,
    OpenBioSingleCellDGIdbAnnotation,
    OpenBioSingleCellDrugScores,
    OpenBioSingleCellDrugHypergeometric,
    OpenBioSingleCellDrugGSEA,
]


__all__ = [
    "ENRICHMENT_NODE_CLASSES",
    "OpenBioSingleCellAUCellScores",
    "OpenBioSingleCellDGIdbAnnotation",
    "OpenBioSingleCellDrugGSEA",
    "OpenBioSingleCellDrugHypergeometric",
    "OpenBioSingleCellDrugScores",
    "OpenBioSingleCellGenePanelScores",
    "OpenBioSingleCellGeneSetOverrepresentation",
    "OpenBioSingleCellGSVAScores",
    "OpenBioSingleCellPathwayScoreTTest",
    "OpenBioSingleCellRankedGSEA",
]
