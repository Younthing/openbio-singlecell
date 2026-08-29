from __future__ import annotations

import gc
import importlib.util
import json
import os
import pickle
import sys
import tempfile
import uuid
import warnings
import weakref
from dataclasses import replace
from importlib import metadata as distribution_metadata
from pathlib import Path
from types import ModuleType

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

import openbio_singlecell.cnmf_standalone as standalone
import openbio_singlecell.operations_factorization as operations_factorization
from openbio_singlecell.artifact_codecs import (
    ANNDATA_PAYLOAD,
    read_anndata,
    read_table,
    write_anndata,
)
from openbio_singlecell.artifact_envelope import summary_from_metadata, table_from_metadata
from openbio_singlecell.artifact_runtime import NATIVE_AFFINITY_CODECS
from openbio_singlecell.cnmf_native_codec import (
    CNMF_NATIVE_CODEC,
    checkout_cnmf_run,
    write_cnmf_run,
)
from openbio_singlecell.cnmf_run import CNMFRun
from openbio_singlecell.nodes_factorization import (
    OpenBioSingleCellCNMF as CNMFNode,
)
from openbio_singlecell.nodes_factorization import (
    OpenBioSingleCellCNMFRankSurvey as CNMFRankSurveyNode,
)
from openbio_singlecell.operations_factorization import cnmf, cnmf_rank_survey
from openbio_singlecell.operations_input import ANNDATA_CODEC, ANNDATA_KIND
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError, registered_operation_ids


def _save_frame(frame: pd.DataFrame, filename: str | os.PathLike[str]) -> None:
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        data=frame.values,
        index=frame.index.values,
        columns=frame.columns.values,
    )


def _load_frame(filename: str | os.PathLike[str]) -> pd.DataFrame:
    with np.load(filename, allow_pickle=True) as payload:
        return pd.DataFrame(**payload)


class FakeCNMF:
    instances: list[FakeCNMF] = []
    fail_prepare = False
    fail_factorize = False
    omit_restart = False
    tamper_private_counts = False
    bad_result_axis = False
    bad_program_labels = False
    bad_density_cache = False
    omit_result_file = False
    warning_mode: str | None = None

    @classmethod
    def reset(cls) -> None:
        cls.instances = []
        cls.fail_prepare = False
        cls.fail_factorize = False
        cls.omit_restart = False
        cls.tamper_private_counts = False
        cls.bad_result_axis = False
        cls.bad_program_labels = False
        cls.bad_density_cache = False
        cls.omit_result_file = False
        cls.warning_mode = None

    def __init__(self, output_dir=".", name=None):
        self.output_dir = output_dir
        self.name = name
        Path(output_dir, name, "cnmf_tmp").mkdir(parents=True, exist_ok=True)
        self.paths = standalone._expected_paths(Path(output_dir), str(name))
        self.calls: list[tuple[str, dict[str, object]]] = []
        self.results: dict[tuple[int, float], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}
        replicate_path = Path(self.paths["nmf_replicate_parameters"])
        genes_path = Path(self.paths["nmf_genes_list"])
        if replicate_path.is_file() and genes_path.is_file():
            replicate = _load_frame(replicate_path)
            self.components = tuple(int(value) for value in replicate["n_components"].drop_duplicates())
            self.n_iter = int(replicate.groupby("n_components").size().iloc[0])
            self.genes = genes_path.read_text(encoding="utf-8").splitlines()
        type(self).instances.append(self)

    def prepare(
        self,
        counts_fn,
        components,
        n_iter=100,
        densify=False,
        tpm_fn=None,
        seed=None,
        beta_loss="frobenius",
        num_highvar_genes=2000,
        genes_file=None,
        alpha_usage=0.0,
        alpha_spectra=0.0,
        init="random",
        max_NMF_iter=1000,
    ):
        self.calls.append(
            (
                "prepare",
                {
                    "counts_fn": counts_fn,
                    "components": tuple(int(value) for value in components),
                    "n_iter": n_iter,
                    "densify": densify,
                    "tpm_fn": tpm_fn,
                    "seed": seed,
                    "beta_loss": beta_loss,
                    "num_highvar_genes": num_highvar_genes,
                    "genes_file": genes_file,
                    "alpha_usage": alpha_usage,
                    "alpha_spectra": alpha_spectra,
                    "init": init,
                    "max_NMF_iter": max_NMF_iter,
                },
            )
        )
        if self.fail_prepare:
            raise RuntimeError("fake prepare failure")
        counts = ad.read_h5ad(counts_fn)
        if self.tamper_private_counts:
            counts.X = np.asarray(counts.X) + 1
            counts.write_h5ad(counts_fn)
        self.components = tuple(int(value) for value in components)
        self.n_iter = int(n_iter)
        self.seed = int(seed)
        self.genes = [str(value) for value in counts.var_names[: min(int(num_highvar_genes), counts.n_vars)]]
        Path(self.paths["nmf_genes_list"]).write_text("\n".join(self.genes), encoding="utf-8")

        values = counts.X.toarray() if hasattr(counts.X, "toarray") else np.asarray(counts.X)
        highvar = values[:, : len(self.genes)].astype(float)
        std = highvar.std(axis=0)
        std[std == 0] = 1.0
        norm = ad.AnnData(
            X=highvar / std,
            obs=pd.DataFrame(index=counts.obs_names.copy()),
            var=pd.DataFrame(index=pd.Index(self.genes)),
        )
        norm.write_h5ad(self.paths["normalized_counts"])
        totals = values.sum(axis=1, keepdims=True)
        tpm_values = values / totals * 1_000_000
        tpm = ad.AnnData(
            X=tpm_values,
            obs=pd.DataFrame(index=counts.obs_names.copy()),
            var=pd.DataFrame(index=counts.var_names.copy()),
        )
        tpm.write_h5ad(self.paths["tpm"])
        tpm_stats = pd.DataFrame(
            {"__mean": tpm_values.mean(axis=0), "__std": tpm_values.std(axis=0)},
            index=counts.var_names,
        )
        _save_frame(tpm_stats, self.paths["tpm_stats"])
        np.random.seed(seed)
        rows = []
        for k in self.components:
            for iteration in range(self.n_iter):
                rows.append((k, iteration, int(np.random.randint(1, 2**31 - 1)), False))
        replicate = pd.DataFrame(rows, columns=["n_components", "iter", "nmf_seed", "completed"])
        _save_frame(replicate, self.paths["nmf_replicate_parameters"])
        run_parameters = {
            "alpha_H": 0.0,
            "alpha_W": 0.0,
            "beta_loss": "frobenius",
            "init": "random",
            "l1_ratio": 0.0,
            "max_iter": 1000,
            "solver": "cd",
            "tol": 0.0001,
        }
        Path(self.paths["nmf_run_parameters"]).write_text(
            yaml.safe_dump(run_parameters, sort_keys=True), encoding="utf-8"
        )

    def factorize(self, worker_i=0, total_workers=1, skip_completed_runs=False):
        self.calls.append(
            (
                "factorize",
                {
                    "worker_i": worker_i,
                    "total_workers": total_workers,
                    "skip_completed_runs": skip_completed_runs,
                },
            )
        )
        if self.warning_mode == "unknown":
            warnings.warn("unexpected fake backend warning", RuntimeWarning, stacklevel=1)
        if self.warning_mode == "known":
            warnings.warn_explicit(
                "unclosed file <_io.TextIOWrapper name='private.nmf_idvrun_params.yaml' mode='r' encoding='utf-8'>",
                ResourceWarning,
                "C:/fake/cnmf/cnmf.py",
                727,
            )
        for k in self.components:
            for iteration in range(self.n_iter):
                if self.omit_restart and k == self.components[-1] and iteration == self.n_iter - 1:
                    continue
                values = np.full((k, len(self.genes)), 0.05, dtype=float)
                for topic in range(k):
                    values[topic, topic % len(self.genes)] = 2.0 + 0.02 * iteration
                    values[topic, :] += 0.001 * (iteration + 1) * np.arange(1, len(self.genes) + 1)
                frame = pd.DataFrame(values, index=np.arange(1, k + 1), columns=self.genes)
                _save_frame(frame, self.paths["iter_spectra"] % (k, iteration))
                if self.fail_factorize:
                    raise RuntimeError("fake factorize failure")

    def combine(self, components=None, skip_missing_files=False):
        self.calls.append(
            (
                "combine",
                {"components": tuple(int(value) for value in components), "skip_missing_files": skip_missing_files},
            )
        )
        for k in components:
            frames = []
            for iteration in range(self.n_iter):
                frame = _load_frame(self.paths["iter_spectra"] % (int(k), iteration))
                frame.index = [f"iter{iteration}_topic{topic}" for topic in range(1, int(k) + 1)]
                frames.append(frame)
            _save_frame(pd.concat(frames), self.paths["merged_spectra"] % int(k))

    def consensus(
        self,
        k,
        density_threshold=0.5,
        local_neighborhood_size=0.30,
        show_clustering=True,
        build_ref=True,
        skip_density_and_return_after_stats=False,
        close_clustergram_fig=False,
        refit_usage=True,
        normalize_tpm_spectra=False,
        norm_counts=None,
    ):
        self.calls.append(
            (
                "consensus",
                {
                    "k": k,
                    "density_threshold": density_threshold,
                    "local_neighborhood_size": local_neighborhood_size,
                    "show_clustering": show_clustering,
                    "build_ref": build_ref,
                    "skip_density_and_return_after_stats": skip_density_and_return_after_stats,
                    "close_clustergram_fig": close_clustergram_fig,
                    "refit_usage": refit_usage,
                    "normalize_tpm_spectra": normalize_tpm_spectra,
                    "norm_counts": norm_counts,
                },
            )
        )
        if skip_density_and_return_after_stats:
            return pd.DataFrame(
                [float(k), float(density_threshold), 0.90 - 0.01 * int(k), 100.0 / int(k)],
                index=["k", "local_density_threshold", "silhouette", "prediction_error"],
                columns=["stats"],
            )
        merged = _load_frame(self.paths["merged_spectra"] % int(k))
        density, _ = standalone._independent_density(
            merged,
            selected_k=int(k),
            n_iter=self.n_iter,
            density_threshold=float(density_threshold),
            local_neighborhood_size=float(local_neighborhood_size),
        )
        density_frame = density.to_frame()
        if self.bad_density_cache:
            density_frame.iloc[0, 0] += 0.25
        _save_frame(density_frame, self.paths["local_density_cache"] % int(k))
        obs = [str(value) for value in norm_counts.obs_names]
        genes = [f"gene_{index}" for index in range(self._input_gene_count())]
        usage_values = np.fromfunction(lambda row, column: (row + column + 2) % (int(k) + 2) + 1, (len(obs), int(k)))
        usage = pd.DataFrame(usage_values, index=obs, columns=np.arange(1, int(k) + 1))
        score_values = np.fromfunction(
            lambda gene, program: (program + 1) * 10 - gene,
            (len(genes), int(k)),
            dtype=float,
        )
        scores = pd.DataFrame(score_values, index=genes, columns=np.arange(1, int(k) + 1))
        spectra = pd.DataFrame(np.abs(score_values) + 1, index=genes, columns=np.arange(1, int(k) + 1))
        top = pd.DataFrame(
            {program: scores[program].sort_values(ascending=False).index.tolist() for program in scores.columns}
        )
        self.results[(int(k), float(density_threshold))] = (usage, scores, spectra, top)
        token = str(float(density_threshold)).replace(".", "_")
        consensus_spectra = pd.DataFrame(
            np.ones((int(k), len(self.genes))), index=np.arange(1, int(k) + 1), columns=self.genes
        )
        outputs = {
            "consensus_spectra": consensus_spectra,
            "consensus_usages": usage,
            "gene_spectra_score": scores.T,
            "gene_spectra_tpm": spectra.T,
        }
        for key, frame in outputs.items():
            _save_frame(frame, self.paths[key] % (int(k), token))
            text_key = f"{key}__txt"
            text_path = self.paths[text_key] % (int(k), token)
            if self.omit_result_file and key == "gene_spectra_tpm":
                continue
            frame.to_csv(text_path, sep="\t")
        return None

    def _input_gene_count(self) -> int:
        return int(ad.read_h5ad(self.paths["tpm"]).n_vars)

    def load_results(self, K, density_threshold, n_top_genes=100, norm_usage=True):
        self.calls.append(
            (
                "load_results",
                {
                    "K": K,
                    "density_threshold": density_threshold,
                    "n_top_genes": n_top_genes,
                    "norm_usage": norm_usage,
                },
            )
        )
        usage, scores, spectra, top = (frame.copy() for frame in self.results[(int(K), float(density_threshold))])
        if norm_usage:
            usage = usage.div(usage.sum(axis=1), axis=0)
        top = top.iloc[: min(int(n_top_genes), len(top)), :]
        if self.bad_result_axis:
            scores = scores.iloc[:-1, :]
        if self.bad_program_labels:
            usage.columns = ["01", *list(usage.columns[1:])]
        return usage, scores, spectra, top


@pytest.fixture(autouse=True)
def fake_cnmf_module(monkeypatch, tmp_path):
    global _TEST_ROOT
    FakeCNMF.reset()
    _TEST_ROOT = tmp_path
    _ACTIVE_CHECKOUTS.clear()
    _RUN_ARTIFACTS.clear()
    module = ModuleType("cnmf")
    module.__version__ = "1.7.1"
    module.cNMF = FakeCNMF
    module.load_df_from_npz = _load_frame
    monkeypatch.setitem(sys.modules, "cnmf", module)
    yield module
    while _ACTIVE_CHECKOUTS:
        _ACTIVE_CHECKOUTS.pop().__exit__(None, None, None)
    _RUN_ARTIFACTS.clear()
    for instance in FakeCNMF.instances:
        del instance
    gc.collect()


def _adata(*, cells: int = 16, genes: int = 8, noninteger: bool = False) -> ad.AnnData:
    rng = np.random.default_rng(91)
    values = rng.poisson(3.0, size=(cells, genes)).astype(float) + 1.0
    values[:, 0] += np.arange(cells) % 3
    if noninteger:
        values += 0.25
    result = ad.AnnData(X=values.copy())
    result.obs_names = [f"cell_{index}" for index in range(cells)]
    result.var_names = [f"gene_{index}" for index in range(genes)]
    result.layers["counts"] = values.copy()
    return result


_TEST_ROOT: Path
_ACTIVE_CHECKOUTS: list[object] = []
_RUN_ARTIFACTS: dict[int, Path] = {}


def _operation_staging(prefix: str) -> Path:
    root = Path(tempfile.mkdtemp(prefix=prefix, dir=_TEST_ROOT))
    staging = root / "run.partial"
    staging.mkdir()
    return staging


def _run_survey(adata: ad.AnnData | None = None, **kwargs):
    if adata is None:
        adata = _adata()
    arguments = _survey_parameters()
    arguments.update(kwargs)
    staging = _operation_staging("survey-")
    input_root = staging.parent / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    records = cnmf_rank_survey(
        OperationContext(staging, str(uuid.uuid4())),
        {"adata": _artifact_descriptor(input_root, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)},
        arguments,
    )
    by_name = {record["name"]: record for record in records}
    run_root = staging / by_name["run"]["payload"]
    checkout = checkout_cnmf_run(run_root, checkout_parent=staging.parent)
    run = checkout.__enter__()
    _ACTIVE_CHECKOUTS.append(checkout)
    _RUN_ARTIFACTS[id(run)] = run_root
    table, table_metadata = read_table(staging / by_name["k_metrics"]["payload"])
    return (
        run,
        table_from_metadata(table_metadata, table),
        summary_from_metadata(by_name["summary"]["value"]),
        by_name["code"]["value"],
    )


def _run_live_survey(adata: ad.AnnData | None = None, **kwargs):
    if adata is None:
        adata = _adata()
    arguments = {
        "source": "layer:counts",
        "components_min": 2,
        "components_max": 3,
        "n_iter": 4,
        "num_highvar_genes": 6,
        "random_seed": 17,
    }
    arguments.update(kwargs)
    return standalone.cnmf_rank_survey(adata, **arguments)


def _run_consensus(run: CNMFRun, **kwargs):
    parameters = {
        "selected_k": 2,
        "density_threshold": 2.0,
        "local_neighborhood_size": 0.3,
        "n_top_genes": 3,
        "overwrite_existing": False,
    }
    parameters.update(kwargs)
    staging = _operation_staging("consensus-")
    records = cnmf(
        OperationContext(staging, str(uuid.uuid4())),
        {
            "run": _artifact_descriptor(
                _RUN_ARTIFACTS[id(run)],
                kind="OPENBIO_CNMF_RUN",
                codec=CNMF_NATIVE_CODEC,
            )
        },
        parameters,
    )
    by_name = {record["name"]: record for record in records}
    return (
        read_anndata(staging / by_name["adata"]["payload"]),
        summary_from_metadata(by_name["summary"]["value"]),
        by_name["code"]["value"],
    )


def _artifact_descriptor(root: Path, *, kind: str, codec: str) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": kind,
        "codec": codec,
        "path": str(root.resolve()),
    }


def _survey_parameters() -> dict[str, object]:
    return {
        "source": {"source": "layer", "layer_name": "counts"},
        "components_min": 2,
        "components_max": 3,
        "n_iter": 4,
        "num_highvar_genes": 6,
        "random_seed": 17,
    }


def test_rank_survey_takes_ownership_of_worker_private_anndata():
    owned = _adata()

    run, _ = _run_live_survey(owned)

    assert object.__getattribute__(run, "_base_adata") is owned
    run.close()


def test_private_count_h5ad_writer_reuses_worker_owned_matrix(tmp_path, monkeypatch):
    standalone._load_science()
    matrix = np.arange(12, dtype=float).reshape(4, 3)
    obs_names = pd.Index([f"cell_{index}" for index in range(4)])
    var_names = pd.Index([f"gene_{index}" for index in range(3)])
    fingerprint = standalone._matrix_fingerprint(matrix, obs_names, var_names)
    original = standalone.ad.AnnData
    captured: dict[str, object] = {}

    def capture_anndata(*, X, obs, var):
        captured["X"] = X
        return original(X=X, obs=obs, var=var)

    monkeypatch.setattr(standalone.ad, "AnnData", capture_anndata)

    standalone._write_private_counts(tmp_path, matrix, obs_names, var_names, fingerprint)

    assert captured["X"] is matrix


def test_rank_survey_operation_publishes_closed_native_run_and_jsonl_metrics(tmp_path):
    adata = _adata()
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_payload = input_root / ANNDATA_PAYLOAD
    before = input_payload.read_bytes()
    staging = tmp_path / "survey.partial"
    staging.mkdir()
    context = OperationContext(staging, "00000000-0000-4000-8000-000000000010")

    records = cnmf_rank_survey(
        context,
        {"adata": _artifact_descriptor(input_root, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)},
        _survey_parameters(),
    )

    assert [(record["type"], record["name"]) for record in records] == [
        ("artifact", "run"),
        ("artifact", "k_metrics"),
        ("summary", "summary"),
        ("string", "code"),
    ]
    assert records[0] == {
        "type": "artifact",
        "name": "run",
        "kind": "OPENBIO_CNMF_RUN",
        "codec": CNMF_NATIVE_CODEC,
        "payload": "outputs/run",
    }
    assert records[1] == {
        "type": "artifact",
        "name": "k_metrics",
        "kind": "OPENBIO_SINGLE_CELL_TABLE",
        "codec": "table-jsonl-v1",
        "payload": "outputs/k_metrics",
    }
    table, metadata = read_table(staging / records[1]["payload"])
    assert table["k"].tolist() == [2, 3]
    assert metadata["kind"] == "table"
    assert input_payload.read_bytes() == before
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()
    assert CNMF_NATIVE_CODEC in NATIVE_AFFINITY_CODECS
    json.dumps(records, allow_nan=False)

    with pytest.raises(ProtocolError, match="parameters"):
        cnmf_rank_survey(
            context,
            {"adata": _artifact_descriptor(input_root, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)},
            _survey_parameters() | {"unexpected": True},
        )


def test_consensus_operation_uses_private_checkout_and_preserves_native_input(tmp_path):
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _adata())
    survey_staging = tmp_path / "survey.partial"
    survey_staging.mkdir()
    survey_records = cnmf_rank_survey(
        OperationContext(survey_staging, "00000000-0000-4000-8000-000000000011"),
        {"adata": _artifact_descriptor(input_root, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)},
        _survey_parameters(),
    )
    run_root = survey_staging / survey_records[0]["payload"]
    inventory = {
        path.relative_to(run_root).as_posix(): path.read_bytes() for path in run_root.rglob("*") if path.is_file()
    }
    consensus_staging = tmp_path / "consensus.partial"
    consensus_staging.mkdir()

    records = cnmf(
        OperationContext(consensus_staging, "00000000-0000-4000-8000-000000000012"),
        {"run": _artifact_descriptor(run_root, kind="OPENBIO_CNMF_RUN", codec=CNMF_NATIVE_CODEC)},
        {
            "selected_k": 2,
            "density_threshold": 2.0,
            "local_neighborhood_size": 0.5,
            "n_top_genes": 3,
            "overwrite_existing": False,
        },
    )

    assert [(record["type"], record["name"]) for record in records] == [
        ("artifact", "adata"),
        ("summary", "summary"),
        ("string", "code"),
    ]
    output = read_anndata(consensus_staging / records[0]["payload"])
    assert output.obsm["X_cnmf_usage"].shape == (16, 2)
    assert {
        path.relative_to(run_root).as_posix(): path.read_bytes() for path in run_root.rglob("*") if path.is_file()
    } == inventory
    checkout_root = Path(FakeCNMF.instances[-1].output_dir)
    assert checkout_root != run_root
    assert not checkout_root.exists()
    json.dumps(records, allow_nan=False)


def test_native_cnmf_artifact_uses_a_private_writable_checkout(tmp_path):
    run, _ = _run_live_survey()
    artifact = tmp_path / "artifact"
    artifact.mkdir()
    write_cnmf_run(artifact, run)
    source_inventory = {
        path.relative_to(artifact).as_posix(): path.read_bytes() for path in artifact.rglob("*") if path.is_file()
    }

    with checkout_cnmf_run(artifact, checkout_parent=tmp_path) as checkout:
        checkout_root = Path(checkout.private_root)
        assert (
            replace(
                checkout.metadata,
                backend_paths_fingerprint=run.metadata.backend_paths_fingerprint,
            )
            == run.metadata
        )
        assert checkout.metadata.backend_paths_fingerprint != run.metadata.backend_paths_fingerprint
        assert checkout.metrics == run.metrics
        assert checkout_root != artifact
        assert checkout_root.is_dir()
        (checkout_root / "consumer-owned.tmp").write_text("mutable", encoding="utf-8")

    assert not checkout_root.exists()
    assert {
        path.relative_to(artifact).as_posix(): path.read_bytes() for path in artifact.rglob("*") if path.is_file()
    } == source_inventory
    run.close()


def test_native_cnmf_codec_rejects_linked_members(tmp_path):
    run, _ = _run_live_survey()
    outside = tmp_path / "outside.txt"
    outside.write_text("outside", encoding="utf-8")
    linked = Path(run.private_root) / "linked.txt"
    try:
        linked.symlink_to(outside)
    except OSError as error:
        run.close()
        pytest.skip(f"symlink creation unavailable: {error}")
    artifact = tmp_path / "artifact"
    artifact.mkdir()

    with pytest.raises(RuntimeError, match="link or reparse point"):
        write_cnmf_run(artifact, run)

    run.close()


def test_public_schemas_and_typed_staged_outputs_are_preserved():
    survey = CNMFRankSurveyNode.define_schema()
    consensus = CNMFNode.define_schema()
    assert "execute" not in vars(CNMFRankSurveyNode)
    assert "execute" not in vars(CNMFNode)
    assert {"openbio.node.cnmfranksurvey", "openbio.node.cnmf"}.issubset(registered_operation_ids())
    assert [item.id for item in survey.inputs] == [
        "adata",
        "source",
        "components_min",
        "components_max",
        "n_iter",
        "num_highvar_genes",
        "random_seed",
    ]
    assert [item.io_type for item in survey.outputs] == [
        "OPENBIO_CNMF_RUN",
        "OPENBIO_SINGLE_CELL_TABLE",
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]
    assert [item.id for item in consensus.inputs] == [
        "run",
        "selected_k",
        "density_threshold",
        "local_neighborhood_size",
        "n_top_genes",
        "overwrite_existing",
    ]
    assert [item.io_type for item in consensus.outputs] == [
        "OPENBIO_ANNDATA",
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]
    survey_inputs = {item.id: item for item in survey.inputs}
    assert survey_inputs["components_min"].min == 2 and survey_inputs["components_min"].max is None
    assert survey_inputs["components_max"].min == 2 and survey_inputs["components_max"].max is None
    assert survey_inputs["n_iter"].min == 2 and survey_inputs["n_iter"].max is None
    assert survey_inputs["num_highvar_genes"].min == 1 and survey_inputs["num_highvar_genes"].max is None
    consensus_inputs = {item.id: item for item in consensus.inputs}
    assert consensus_inputs["selected_k"].min == 2 and consensus_inputs["selected_k"].max is None


def test_public_raw_source_uses_full_raw_axis_and_generated_code_has_parity():
    full = _adata(genes=10)
    full.raw = full
    adata = full[:, full.var_names[:5]].copy()
    raw_names = list(adata.raw.var_names)
    source_input = next(item for item in CNMFRankSurveyNode.define_schema().inputs if item.id == "source")
    assert [option.key for option in source_input.options] == ["layer", "X", "raw"]

    run, table_result, report, survey_code = _run_survey(
        adata,
        source={"source": "raw"},
        num_highvar_genes=7,
    )
    assert run.metadata.source_kind == "raw"
    assert run.metadata.input_genes == 10
    assert run.metadata.current_features == 5
    assert report.summary["key_results"]["source"] == "raw"
    assert report.summary["key_results"]["source_features"] == 10
    assert report.summary["key_results"]["current_features"] == 5
    assert table_result.table["k"].tolist() == [2, 3]
    producer = FakeCNMF.instances[-2]
    private_counts = ad.read_h5ad(
        _RUN_ARTIFACTS[id(run)] / Path(producer.calls[0][1]["counts_fn"]).relative_to(producer.output_dir)
    )
    assert list(private_counts.var_names) == raw_names

    runtime_output, runtime_report, consensus_code = _run_consensus(
        run,
        selected_k=2,
        n_top_genes=3,
    )
    assert list(runtime_output.var_names) == raw_names
    assert runtime_output.varm["cnmf_gep_scores"].shape == (10, 2)
    assert runtime_output.uns["cnmf"]["survey"]["source_features"] == 10
    assert runtime_output.uns["cnmf"]["survey"]["current_features"] == 5
    assert runtime_report.summary["key_results"]["source"] == "raw"
    assert runtime_report.summary["key_results"]["source_features"] == 10
    assert runtime_report.summary["key_results"]["current_features"] == 5

    survey_namespace: dict[str, object] = {}
    exec(survey_code, survey_namespace)
    generated_run, generated_metrics = survey_namespace["cnmf_rank_survey"](adata)
    assert generated_run.metadata.source_kind == "raw"
    assert generated_run.metadata.input_genes == 10
    assert generated_run.metadata.current_features == 5
    pd.testing.assert_frame_equal(generated_metrics, table_result.table)
    consensus_namespace: dict[str, object] = {}
    exec(consensus_code, consensus_namespace)
    generated_output = consensus_namespace["cnmf_consensus_programs"](generated_run)
    assert list(generated_output.var_names) == raw_names
    assert generated_output.uns["cnmf"]["source"] == runtime_output.uns["cnmf"]["source"]
    assert generated_output.uns["cnmf"]["survey"]["source_features"] == 10
    assert generated_output.uns["cnmf"]["survey"]["current_features"] == 5
    assert {
        key: runtime_report.summary["key_results"][key] for key in ("source", "source_features", "current_features")
    } == {
        "source": generated_output.uns["cnmf"]["source"]["kind"],
        "source_features": generated_output.uns["cnmf"]["survey"]["source_features"],
        "current_features": generated_output.uns["cnmf"]["survey"]["current_features"],
    }
    pd.testing.assert_frame_equal(runtime_output.varm["cnmf_gep_scores"], generated_output.varm["cnmf_gep_scores"])
    pd.testing.assert_frame_equal(runtime_output.varm["cnmf_gep_tpm"], generated_output.varm["cnmf_gep_tpm"])
    run.close()
    generated_run.close()


def test_rank_survey_uses_exact_file_backed_official_chain_and_reports_integrity():
    adata = _adata()
    before = adata.copy()
    legacy_rng_before = np.random.get_state()
    run, table_result, report, code = _run_survey(adata)
    legacy_rng_after = np.random.get_state()
    model = FakeCNMF.instances[-2]
    assert isinstance(run, CNMFRun)
    assert Path(run.private_root).is_dir()
    assert Path(model.output_dir) != Path(run.private_root)
    assert not Path(model.output_dir).exists()
    assert [name for name, _ in model.calls] == [
        "prepare",
        "factorize",
        "combine",
        "consensus",
        "consensus",
    ]
    prepare = model.calls[0][1]
    assert Path(str(prepare["counts_fn"])).parent == Path(model.output_dir)
    assert prepare == {
        "counts_fn": prepare["counts_fn"],
        "components": (2, 3),
        "n_iter": 4,
        "densify": False,
        "tpm_fn": None,
        "seed": 17,
        "beta_loss": "frobenius",
        "num_highvar_genes": 6,
        "genes_file": None,
        "alpha_usage": 0.0,
        "alpha_spectra": 0.0,
        "init": "random",
        "max_NMF_iter": 1000,
    }
    assert model.calls[1][1] == {"worker_i": 0, "total_workers": 1, "skip_completed_runs": False}
    assert model.calls[2][1] == {"components": (2, 3), "skip_missing_files": False}
    assert table_result.table["k"].tolist() == [2, 3]
    assert table_result.table["completed_restarts"].tolist() == [4, 4]
    assert run.metadata.backend_name == "cnmf.cNMF"
    assert run.metadata.backend_version == "1.7.1"
    assert run.metadata.installed_distribution_attested is False
    assert run.metadata.audited_pypi_wheel_sha256 == standalone.CNMF_AUDITED_PYPI_WHEEL_SHA256
    assert len(run.metadata.artifact_hashes) == 7 + 8 + 2
    key_results = report.summary["key_results"]
    assert key_results["artifact_manifest_sha256"].startswith("sha256:")
    assert key_results["execution_mode"] == "CPU, single-worker, owned private file-backed run"
    assert key_results["cleanup_ownership"]["live_run"] == "closed after native artifact encoding"
    assert "Classic cache" in key_results["cleanup_ownership"]["node_graph"]
    assert key_results["cleanup_ownership"]["cached_directory_may_persist"] is True
    assert "declared UMI count source" not in report.summary["methods"]
    assert report.summary["software_versions"]["cnmf"] == "1.7.1"
    assert all("OmicVerse" not in item["citation"] for item in report.summary["references"])
    assert "from openbio_singlecell" not in code
    assert "import openbio_singlecell" not in code
    assert "TemporaryDirectory" in code
    assert 'CNMF_REQUIRED_VERSION = "1.7.1"' in code
    compile(code, "<generated-cnmf-survey>", "exec")
    np.testing.assert_array_equal(adata.X, before.X)
    np.testing.assert_array_equal(adata.layers["counts"], before.layers["counts"])
    assert legacy_rng_before[0] == legacy_rng_after[0]
    np.testing.assert_array_equal(legacy_rng_before[1], legacy_rng_after[1])
    assert legacy_rng_before[2:] == legacy_rng_after[2:]
    run.close()


def test_owned_directory_lives_until_explicit_idempotent_close():
    run, _ = _run_live_survey()
    root = Path(run.private_root)
    assert root.exists()
    run.close()
    assert run.closed
    assert not root.exists()
    run.close()
    with pytest.raises(RuntimeError, match="closed"):
        run.copy_base_adata()
    with pytest.raises(RuntimeError, match="closed"):
        standalone.cnmf_consensus_programs(run, selected_k=2, n_top_genes=3)


def test_owned_directory_finalizer_is_cleanup_fallback():
    run, _ = _run_live_survey()
    root = Path(run.private_root)
    reference = weakref.ref(run)
    del run
    gc.collect()
    assert reference() is None
    assert not root.exists()


def test_cleanup_failure_keeps_retry_and_finalizer_paths_reachable(monkeypatch):
    run, _ = _run_live_survey()
    root = Path(run.private_root)
    original_cleanup = tempfile.TemporaryDirectory.cleanup
    attempts = 0

    def fail_once(directory):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError("simulated cleanup failure")
        return original_cleanup(directory)

    monkeypatch.setattr(tempfile.TemporaryDirectory, "cleanup", fail_once)
    with pytest.raises(OSError, match="simulated cleanup failure"):
        run.close()
    assert run.closed
    assert root.exists()
    assert run._finalizer.alive
    with pytest.raises(RuntimeError, match="closed"):
        run.copy_base_adata()
    run.close()
    assert attempts == 2
    assert not root.exists()
    assert not run._finalizer.alive


def test_construction_cleanup_failure_cannot_mask_primary_analysis_error(monkeypatch):
    FakeCNMF.fail_prepare = True
    original_cleanup = tempfile.TemporaryDirectory.cleanup
    attempts = 0

    def fail_once(directory):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError("construction cleanup failed")
        return original_cleanup(directory)

    monkeypatch.setattr(tempfile.TemporaryDirectory, "cleanup", fail_once)
    with pytest.raises(RuntimeError, match="fake prepare failure") as captured:
        _run_survey()
    notes = tuple(captured.value.__notes__)
    assert any("construction cleanup failed" in note for note in notes)
    del captured
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Implicitly cleaning up <TemporaryDirectory.*>",
            category=ResourceWarning,
        )
        gc.collect()
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()


def test_context_body_error_remains_primary_when_close_cleanup_fails(monkeypatch):
    run, _ = _run_live_survey()
    root = Path(run.private_root)
    original_cleanup = tempfile.TemporaryDirectory.cleanup
    attempts = 0

    def fail_once(directory):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError("context cleanup failed")
        return original_cleanup(directory)

    monkeypatch.setattr(tempfile.TemporaryDirectory, "cleanup", fail_once)
    with pytest.raises(ValueError, match="context body failed") as captured:
        with run:
            raise ValueError("context body failed")
    assert any("context cleanup failed" in note for note in captured.value.__notes__)
    assert run.closed and root.exists()
    run.close()
    assert not root.exists()


def test_node_report_error_remains_primary_when_close_cleanup_fails(monkeypatch):
    original_cleanup = tempfile.TemporaryDirectory.cleanup
    attempts = 0

    def fail_once(directory):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError("report cleanup failed")
        return original_cleanup(directory)

    def fail_report(**kwargs):
        del kwargs
        raise LookupError("summary construction failed")

    monkeypatch.setattr(tempfile.TemporaryDirectory, "cleanup", fail_once)
    monkeypatch.setattr(operations_factorization, "make_analysis_report", fail_report)
    with pytest.raises(LookupError, match="summary construction failed") as captured:
        _run_survey()
    notes = tuple(captured.value.__notes__)
    assert any("report cleanup failed" in note for note in notes)
    del captured
    gc.collect()
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()


@pytest.mark.parametrize("stage", ["prepare", "factorize"])
def test_exception_paths_cleanup_private_directory(stage):
    if stage == "prepare":
        FakeCNMF.fail_prepare = True
    else:
        FakeCNMF.fail_factorize = True
    with pytest.raises(RuntimeError, match=f"fake {stage} failure"):
        _run_survey()
    root = Path(FakeCNMF.instances[-1].output_dir)
    assert not root.exists()


def test_private_count_input_is_reopened_and_backend_mutation_rejected():
    FakeCNMF.tamper_private_counts = True
    with pytest.raises(RuntimeError, match="count input changed during backend preparation"):
        _run_survey()
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()


def test_missing_restart_is_never_combined_or_returned():
    FakeCNMF.omit_restart = True
    with pytest.raises(RuntimeError, match=r"spectra\.k_3\.iter_3.*missing|iter_spectra file is missing"):
        _run_survey()
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()


def test_unknown_warning_is_not_swallowed_and_cleans_run():
    FakeCNMF.warning_mode = "unknown"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(RuntimeWarning, match="unexpected fake backend warning"):
            _run_survey()
    assert not Path(FakeCNMF.instances[-1].output_dir).exists()


def test_exact_known_warning_is_isolated_and_disclosed_under_werror():
    FakeCNMF.warning_mode = "known"
    run, _, report, _ = _run_survey()
    assert run.metadata.known_upstream_warnings == (("cnmf-1.7.1-unclosed-yaml-reader", 1),)
    assert any("unclosed YAML-reader" in warning for warning in report.summary["warnings"])
    run.close()


def test_noninteger_and_logged_expert_source_is_advisory_not_gate():
    adata = _adata(noninteger=True)
    adata.uns["log1p"] = {"base": None}
    run, _, report, _ = _run_survey(adata, source={"source": "X"})
    assert len(run.metadata.input_advisories) == 2
    assert any("not integer-like" in warning for warning in report.summary["warnings"])
    assert any("proven to be logged" in warning for warning in report.summary["warnings"])
    output, final_report, _ = _run_consensus(run, selected_k=2, n_top_genes=3)
    assert output.uns["cnmf"]["survey"]["source_state"] == "logged"
    assert output.uns["cnmf"]["survey"]["source_state_evidence"] == "AnnData uns['log1p'] marker"
    assert list(output.uns["cnmf"]["survey"]["input_advisories"]) == list(run.metadata.input_advisories)
    assert final_report.summary["key_results"]["input_advisories"] == list(run.metadata.input_advisories)
    assert any("not integer-like" in warning for warning in final_report.summary["warnings"])
    assert any("proven to be logged" in warning for warning in final_report.summary["warnings"])
    run.close()


def test_scale_to_layer_history_is_classified_consistently_in_runtime_and_code():
    adata = _adata()
    adata.layers["scaled"] = adata.layers["counts"] / 10.0
    adata.uns["openbio_singlecell"] = {
        "analysis_history": {
            "000000": {
                "operation": "scale_to_layer",
                "parameters": {"output_layer": "scaled"},
            }
        }
    }

    run, _, report, code = _run_survey(
        adata,
        source={"source": "layer", "layer_name": "scaled"},
    )
    assert run.metadata.source_state == "scaled"
    assert run.metadata.source_state_evidence == "OpenBio Scale history"
    assert report.summary["key_results"]["source_state"] == "scaled"
    assert any("proven to be scaled" in warning for warning in report.summary["warnings"])

    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_run, _ = namespace["cnmf_rank_survey"](adata)
    assert generated_run.metadata.source_state == "scaled"
    assert generated_run.metadata.source_state_evidence == "OpenBio Scale history"
    generated_run.close()
    run.close()


def test_unknown_noninteger_source_is_not_described_as_integer_valued_or_umi():
    run, _, report, _ = _run_survey(_adata(noninteger=True))
    provenance_warnings = [
        warning for warning in report.summary["warnings"] if warning.startswith("OpenBio provenance cannot establish")
    ]
    assert len(provenance_warnings) == 1
    assert "integer" not in provenance_warnings[0].lower()
    assert "umi" not in provenance_warnings[0].lower()
    assert "input_advisories" in report.summary["key_results"]
    assert report.summary["key_results"]["input_advisories"] == list(run.metadata.input_advisories)
    run.close()


def test_conservative_resource_envelope_is_advisory_for_experts(monkeypatch):
    monkeypatch.setattr(standalone, "CNMF_RESOURCE_BUDGET_BYTES", 1)
    run, _, report, _ = _run_survey()
    assert run.metadata.requested_resource.estimated_peak_bytes > 1
    assert any("exceeds 2 GiB" in warning for warning in report.summary["warnings"])
    run.close()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda value: value.__setitem__((0, 0), -1), "negative"),
        (lambda value: value.__setitem__((0, 0), np.nan), "finite numeric"),
        (lambda value: value.__setitem__((0, slice(None)), 0), "zero total counts"),
    ],
)
def test_backend_required_count_invariants_remain_hard_gates(mutation, message):
    adata = _adata()
    mutation(adata.layers["counts"])
    with pytest.raises(ValueError, match=message):
        _run_survey(adata)


def test_duplicate_and_unsafe_axes_are_rejected_but_overcomplete_expert_dimensions_are_advisory():
    duplicate = _adata()
    duplicate.obs_names = ["same"] * duplicate.n_obs
    with pytest.raises(ValueError, match="unique obs_names"):
        _run_survey(duplicate)
    run, _, report, _ = _run_survey(
        _adata(cells=4),
        components_min=129,
        components_max=129,
        n_iter=2,
        num_highvar_genes=100_001,
    )
    warnings_list = report.summary["warnings"]
    assert any("exceeds the 4 input cells" in warning for warning in warnings_list)
    assert any("exceeds the 8 realized HVGs" in warning for warning in warnings_list)
    assert any("Requested num_highvar_genes=100001" in warning for warning in warnings_list)
    run.close()
    with pytest.raises(ValueError, match="n_iter must be at least 2"):
        _run_survey(n_iter=1)
    unsafe = _adata()
    unsafe.var_names = ["bad\tgene", *[f"gene_{index}" for index in range(1, unsafe.n_vars)]]
    with pytest.raises(ValueError, match="tab-delimited result axes"):
        _run_survey(unsafe)


def test_zero_total_genes_are_advisory_when_backend_can_exclude_them():
    adata = _adata()
    adata.layers["counts"][:, 0] = 0
    run, _, report, _ = _run_survey(adata)
    assert any("zero-total genes" in warning for warning in report.summary["warnings"])
    run.close()


def test_consensus_loads_official_tuple_and_writes_only_canonical_continuous_state(tmp_path):
    adata = _adata()
    run, _, _, _ = _run_survey(adata)
    output, report, code = _run_consensus(
        run,
        selected_k=2,
        density_threshold=2.0,
        local_neighborhood_size=0.5,
        n_top_genes=3,
    )
    assert run.closed is False
    assert Path(run.private_root).exists()
    assert output is not adata
    assert output.obsm["X_cnmf_usage"].shape == (adata.n_obs, 2)
    assert output.varm["cnmf_gep_scores"].shape == (adata.n_vars, 2)
    assert output.varm["cnmf_gep_tpm"].shape == (adata.n_vars, 2)
    assert list(output.obsm["X_cnmf_usage"].columns) == ["cNMF_1", "cNMF_2"]
    assert "cNMF_cluster" not in output.obs
    assert not any(column.startswith("cNMF_") for column in output.obs.columns)
    state = output.uns["cnmf"]
    assert state["components_before_filtering"] == 8
    assert state["components_after_filtering"] == 8
    assert state["cluster_component_counts"] == {"cNMF_1": 4, "cNMF_2": 4}
    assert (
        state["survey"]["immutable_manifest_before_consensus"] == state["survey"]["immutable_manifest_after_consensus"]
    )
    assert report.summary["key_results"]["selected_k"] == 2
    assert report.summary["key_results"]["installed_distribution_attested"] is False
    assert report.summary["key_results"]["source_state"] == "unknown"
    assert report.summary["key_results"]["cleanup_ownership"]["cached_directory_may_persist"] is True
    assert "tamper-evident" not in report.summary["methods"]
    assert any("private checkout" in warning for warning in report.summary["warnings"])
    assert "load_results" in code and "build_ref=False" in code
    assert "from openbio_singlecell" not in code
    compile(code, "<generated-cnmf-consensus>", "exec")
    output_path = tmp_path / "cnmf-output.h5ad"
    output.write_h5ad(output_path)
    reopened = ad.read_h5ad(output_path)
    assert reopened.obsm["X_cnmf_usage"].shape == (adata.n_obs, 2)
    np.testing.assert_array_equal(adata.X, _adata().X)
    run.close()


def test_repeated_consensus_recomputes_density_cache_for_each_neighborhood():
    run, _, _, _ = _run_survey()
    first, _, _ = _run_consensus(
        run, selected_k=2, density_threshold=2.0, local_neighborhood_size=0.25, n_top_genes=3
    )
    second, _, _ = _run_consensus(
        run, selected_k=2, density_threshold=2.0, local_neighborhood_size=0.5, n_top_genes=3
    )
    assert first.uns["cnmf"]["density_neighbors"] == 1
    assert second.uns["cnmf"]["density_neighbors"] == 2
    assert first.uns["cnmf"]["local_density_summary"] != second.uns["cnmf"]["local_density_summary"]
    final_calls = [
        call
        for model in FakeCNMF.instances[-2:]
        for call in model.calls
        if call[0] == "consensus"
    ]
    assert [call[1]["local_neighborhood_size"] for call in final_calls] == [0.25, 0.5]
    run.close()


def test_result_axis_loss_and_density_mismatch_fail_before_annotation():
    run, _, _, _ = _run_survey()
    FakeCNMF.bad_result_axis = True
    with pytest.raises(RuntimeError, match="axis does not exactly match"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    FakeCNMF.bad_result_axis = False
    FakeCNMF.bad_density_cache = True
    with pytest.raises(RuntimeError, match="disagrees with the independent"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    run.close()


def test_noncanonical_program_labels_and_missing_result_file_fail_closed():
    run, _, _, _ = _run_survey()
    FakeCNMF.bad_program_labels = True
    with pytest.raises(RuntimeError, match="invalid program label"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    FakeCNMF.bad_program_labels = False
    FakeCNMF.omit_result_file = True
    with pytest.raises(RuntimeError, match="gene_spectra_tpm__txt.*missing"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    run.close()


def test_collision_policy_is_transactional_and_overwrite_is_explicit():
    adata = _adata()
    adata.uns["cnmf"] = {"old": True}
    run, _, _, _ = _run_survey(adata)
    with pytest.raises(ValueError, match="result keys already exist"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    output, report, _ = _run_consensus(run, selected_k=2, n_top_genes=3, overwrite_existing=True)
    assert output.uns["cnmf"]["overwrote_existing"] is True
    assert report.summary["key_results"]["overwrote_existing"] is True
    assert adata.uns["cnmf"] == {"old": True}
    run.close()


def test_artifact_tamper_and_backend_path_escape_are_rejected(tmp_path):
    run, _, _, _ = _run_survey()
    model = FakeCNMF.instances[-1]
    merged = _RUN_ARTIFACTS[id(run)] / Path(model.paths["merged_spectra"] % 2).relative_to(run.private_root)
    with merged.open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(RuntimeError, match="immutable SHA-256"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    run.close()

    run, _ = _run_live_survey()
    model = FakeCNMF.instances[-1]
    model.paths["merged_spectra"] = str(tmp_path / "escape-k_%d.npz")
    with pytest.raises(RuntimeError, match="escaped or changed"):
        standalone.cnmf_consensus_programs(run, selected_k=2, n_top_genes=3)
    run.close()

    run, _ = _run_live_survey()
    model = FakeCNMF.instances[-1]
    canonical = Path(model.paths["merged_spectra"])
    model.paths["merged_spectra"] = str(canonical.parent / "nested" / ".." / canonical.name)
    with pytest.raises(RuntimeError, match="escaped or changed|unsafe component"):
        standalone.cnmf_consensus_programs(run, selected_k=2, n_top_genes=3)
    run.close()


def test_symlinked_artifact_is_rejected_when_platform_allows_symlinks(tmp_path):
    run, _, _, _ = _run_survey()
    model = FakeCNMF.instances[-1]
    target = _RUN_ARTIFACTS[id(run)] / Path(model.paths["merged_spectra"] % 2).relative_to(run.private_root)
    outside = tmp_path / "outside.npz"
    outside.write_bytes(target.read_bytes())
    target.unlink()
    try:
        target.symlink_to(outside)
    except OSError as exc:
        run.close()
        pytest.skip(f"symlink creation unavailable: {exc}")
    with pytest.raises(RuntimeError, match="link or reparse point"):
        _run_consensus(run, selected_k=2, n_top_genes=3)
    run.close()
    assert outside.exists()


def test_run_cannot_be_pickled():
    run, _ = _run_live_survey()
    with pytest.raises(TypeError, match="process-local"):
        pickle.dumps(run)
    run.close()


def test_metadata_rebinding_and_low_level_content_change_are_rejected():
    run, _ = _run_live_survey()
    original = run.metadata
    with pytest.raises(AttributeError):
        run.metadata = replace(original, source_state="counts")
    forged = replace(
        original,
        source_state="counts",
        source_state_evidence="forged evidence",
        input_advisories=("forged advisory",),
    )
    object.__setattr__(run, "_metadata", forged)
    with pytest.raises(RuntimeError, match="changed after Survey construction|construction snapshot"):
        standalone.cnmf_consensus_programs(run, selected_k=2, n_top_genes=3)
    object.__setattr__(run, "_metadata", original)
    run.close()


def test_metrics_rebinding_and_forged_stability_are_rejected():
    run, _ = _run_live_survey()
    original = run.metrics
    forged = (replace(original[0], stability=0.123456), *original[1:])
    with pytest.raises(AttributeError):
        run.metrics = forged
    object.__setattr__(run, "_metrics", forged)
    with pytest.raises(RuntimeError, match="changed after Survey construction|construction snapshot"):
        standalone.cnmf_consensus_programs(run, selected_k=2, n_top_genes=3)
    object.__setattr__(run, "_metrics", original)
    run.close()


def test_runtime_and_generated_consensus_reject_duck_and_instance_marker_proxies():
    run, _, _, _ = _run_survey()
    _, _, consensus_code = _run_consensus(run, selected_k=2, n_top_genes=3)

    class DuckProxy:
        def __init__(self, target):
            self._target = target
            self._OPENBIO_PORTABLE_ARTIFACT = standalone._RUN_ARTIFACT_MARKER

        def __getattr__(self, name):
            return getattr(self._target, name)

    proxy = DuckProxy(run)
    with pytest.raises(TypeError, match="concrete or explicitly portable"):
        standalone.cnmf_consensus_programs(proxy, selected_k=2, n_top_genes=3)
    generated_namespace: dict[str, object] = {}
    exec(consensus_code, generated_namespace)
    with pytest.raises(TypeError, match="concrete or explicitly portable"):
        generated_namespace["cnmf_consensus_programs"](proxy)
    run.close()


def test_exact_version_and_explicit_signature_guards(monkeypatch):
    monkeypatch.setattr(standalone.importlib_metadata, "version", lambda name: "1.7.2")
    with pytest.raises(RuntimeError, match="exact cnmf==1.7.1"):
        _run_survey()

    monkeypatch.setattr(standalone.importlib_metadata, "version", lambda name: "1.7.1")

    class IncompatibleCNMF:
        def __init__(self, output_dir, name):
            del output_dir, name

        def prepare(self, counts_fn, components, **kwargs):
            del counts_fn, components, kwargs

        def factorize(self, **kwargs):
            del kwargs

        def combine(self, **kwargs):
            del kwargs

        def consensus(self, **kwargs):
            del kwargs

        def load_results(self, **kwargs):
            del kwargs

    module = ModuleType("cnmf")
    module.__version__ = "1.7.1"
    module.cNMF = IncompatibleCNMF
    module.load_df_from_npz = _load_frame
    monkeypatch.setitem(sys.modules, "cnmf", module)
    with pytest.raises(RuntimeError, match="missing explicit parameters"):
        _run_survey()


def _normalized_cnmf_state(state: dict[str, object]) -> dict[str, object]:
    copied = json.loads(json.dumps(state, allow_nan=False, default=lambda value: value.tolist()))
    survey = copied["survey"]
    for key in (
        "backend_paths_fingerprint",
        "prepare_manifest_sha256",
        "factorization_manifest_sha256",
        "artifact_manifest_sha256",
        "immutable_manifest_before_consensus",
        "immutable_manifest_after_consensus",
    ):
        survey.pop(key, None)
    return copied


def test_generated_staged_code_is_standalone_and_scientifically_equivalent():
    adata = _adata()
    runtime_run, _, _, survey_code = _run_survey(adata)
    runtime_output, _, consensus_code = _run_consensus(
        runtime_run, selected_k=2, density_threshold=2.0, local_neighborhood_size=0.5, n_top_genes=3
    )
    survey_namespace: dict[str, object] = {}
    exec(survey_code, survey_namespace)
    generated_run, generated_metrics = survey_namespace["cnmf_rank_survey"](adata)
    pd.testing.assert_frame_equal(generated_metrics, pd.DataFrame([metric.as_dict() for metric in runtime_run.metrics]))
    consensus_namespace: dict[str, object] = {}
    exec(consensus_code, consensus_namespace)
    generated_output = consensus_namespace["cnmf_consensus_programs"](generated_run)
    pd.testing.assert_frame_equal(runtime_output.obsm["X_cnmf_usage"], generated_output.obsm["X_cnmf_usage"])
    pd.testing.assert_frame_equal(runtime_output.varm["cnmf_gep_scores"], generated_output.varm["cnmf_gep_scores"])
    pd.testing.assert_frame_equal(runtime_output.varm["cnmf_gep_tpm"], generated_output.varm["cnmf_gep_tpm"])
    assert _normalized_cnmf_state(runtime_output.uns["cnmf"]) == _normalized_cnmf_state(generated_output.uns["cnmf"])
    runtime_run.close()
    generated_run.close()


def test_real_cnmf_171_tiny_survey_to_consensus_smoke(monkeypatch):
    monkeypatch.delitem(sys.modules, "cnmf", raising=False)
    if importlib.util.find_spec("cnmf") is None:
        pytest.skip("real cnmf==1.7.1 smoke skipped: optional package is not installed")
    try:
        version = distribution_metadata.version("cnmf")
    except distribution_metadata.PackageNotFoundError:
        pytest.skip("real cnmf==1.7.1 smoke skipped: distribution metadata is unavailable")
    if version != "1.7.1":
        pytest.skip(f"real cnmf==1.7.1 smoke skipped: installed version is {version}")
    counts = np.random.default_rng(7).poisson(3.0, size=(24, 12)).astype(float) + 1
    adata = ad.AnnData(X=counts.copy())
    adata.obs_names = [f"cell_{index}" for index in range(24)]
    adata.var_names = [f"gene_{index}" for index in range(12)]
    adata.layers["counts"] = counts.copy()
    run, metrics = standalone.cnmf_rank_survey(
        adata,
        source="layer:counts",
        components_min=2,
        components_max=2,
        n_iter=4,
        num_highvar_genes=8,
        random_seed=11,
    )
    output = standalone.cnmf_consensus_programs(
        run,
        selected_k=2,
        density_threshold=2.0,
        local_neighborhood_size=0.3,
        n_top_genes=3,
    )
    assert metrics["k"].tolist() == [2]
    assert output.obsm["X_cnmf_usage"].shape == (24, 2)
    assert output.varm["cnmf_gep_scores"].shape == (12, 2)
    root = Path(run.private_root)
    run.close()
    assert not root.exists()
