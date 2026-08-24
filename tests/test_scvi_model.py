from __future__ import annotations

import sys
import types

import pytest

from openbio_singlecell.node_types import AnnDataType, SCVIModelType, TableResultType
from openbio_singlecell.nodes_differential import OpenBioSingleCellSCVIDifferentialExpression
from openbio_singlecell.nodes_integration import OpenBioSingleCellSCVIIntegration
from openbio_singlecell.scvi_model import SCVIModel


def output_values(node_output):
    return node_output.result


@pytest.fixture
def adata(science):
    obs = science.pd.DataFrame(
        {"batch": science.pd.Categorical(["one", "one", "one", "two", "two", "two"])},
        index=[f"cell_{index}" for index in range(6)],
    )
    var = science.pd.DataFrame(index=[f"gene_{index}" for index in range(4)])
    value = science.ad.AnnData(science.np.arange(24, dtype=float).reshape(6, 4), obs=obs, var=var)
    value.layers["counts"] = value.X.copy()
    return value


class FakeTrainedSCVI:
    def __init__(self, registered_adata, science):
        self.adata = registered_adata
        self.is_trained = True
        self.science = science
        self.differential_expression_calls = []
        self.deregister_manager_calls = []
        self.differential_expression_error = None
        self.deregister_manager_error = None
        self.train_calls = 0

    def train(self, **kwargs):
        self.train_calls += 1
        raise AssertionError("Differential expression must not retrain the supplied scVI model.")

    def differential_expression(self, **kwargs):
        analysis_adata = kwargs["adata"]
        self.differential_expression_calls.append(kwargs)
        analysis_adata.obs["_fake_internal_mutation"] = True
        if self.differential_expression_error is not None:
            raise self.differential_expression_error
        return self.science.pd.DataFrame(
            {"lfc_mean": [1.5, -0.5], "proba_de": [0.95, 0.8]},
            index=["gene_0", "gene_1"],
        )

    def deregister_manager(self, analysis_adata=None):
        self.deregister_manager_calls.append(analysis_adata)
        if self.deregister_manager_error is not None:
            raise self.deregister_manager_error


def test_scvi_schemas_use_the_concrete_model_wire():
    integration = OpenBioSingleCellSCVIIntegration.GET_SCHEMA()
    integration_inputs = {input_.id: input_ for input_ in integration.inputs}
    assert [(output.display_name, output.io_type) for output in integration.outputs] == [
        ("adata", AnnDataType.io_type),
        ("model", SCVIModelType.io_type),
    ]
    assert integration_inputs["source"].default == "layer"
    assert integration_inputs["source"].advanced is not True
    assert integration_inputs["counts_layer"].advanced is not True
    assert integration_inputs["batch_key"].default == ""
    assert integration_inputs["size_factor_key"].advanced is True
    assert integration_inputs["compute_mde"].default is False
    assert integration_inputs["compute_mde"].advanced is True
    assert integration_inputs["store_latent_distribution"].default is False
    assert integration_inputs["store_latent_distribution"].advanced is True

    differential = OpenBioSingleCellSCVIDifferentialExpression.GET_SCHEMA()
    differential_inputs = {input_.id: input_ for input_ in differential.inputs}
    assert [input_.id for input_ in differential.inputs] == [
        "adata",
        "model",
        "groupby",
        "group1",
        "group2",
        "subset_column",
        "subset_value",
        "mode",
        "delta",
    ]
    assert differential.inputs[0].io_type == AnnDataType.io_type
    assert differential.inputs[1].io_type == SCVIModelType.io_type
    assert [(output.display_name, output.io_type) for output in differential.outputs] == [
        ("table", TableResultType.io_type)
    ]
    assert differential_inputs["subset_column"].advanced is not True
    assert differential_inputs["subset_value"].advanced is not True
    assert differential_inputs["mode"].default == "vanilla"
    assert differential_inputs["mode"].advanced is not True
    assert differential_inputs["delta"].advanced is True
    assert "change" in differential_inputs["delta"].tooltip


def test_scvi_integration_returns_the_trained_model_wrapper(adata, science, monkeypatch):
    class FakeSCVI:
        setup_calls = []
        instances = []

        @classmethod
        def setup_anndata(cls, registered_adata, **kwargs):
            cls.setup_calls.append((registered_adata, kwargs))

        def __init__(self, registered_adata, **kwargs):
            self.adata = registered_adata
            self.is_trained = False
            self.init_parameters = kwargs
            self.train_parameters = None
            self.n_latent = kwargs["n_latent"]
            self.latent_calls = []
            type(self).instances.append(self)

        def train(self, **kwargs):
            self.train_parameters = kwargs
            self.is_trained = True

        def get_latent_representation(self, *, give_mean=True, return_dist=False):
            self.latent_calls.append((give_mean, return_dist))
            mean = science.np.ones((self.adata.n_obs, self.n_latent))
            if not give_mean and return_dist:
                return mean * 2, mean * 3
            return mean

        def differential_expression(self, **kwargs):
            raise AssertionError("Integration must not run differential expression.")

        def deregister_manager(self, analysis_adata):
            raise AssertionError("Integration must not create temporary AnnData managers.")

    mde_calls = []
    fake_scvi = types.ModuleType("scvi")
    fake_scvi.settings = types.SimpleNamespace(seed=None)
    fake_scvi.model = types.SimpleNamespace(
        SCVI=FakeSCVI,
        utils=types.SimpleNamespace(mde=lambda values: mde_calls.append(values) or values[:, :2]),
    )
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata.obs["total_counts"] = science.np.asarray(adata.X.sum(axis=1)).ravel()
    output, trained_model = output_values(
        OpenBioSingleCellSCVIIntegration.execute(
            adata,
            source="X",
            counts_layer="missing_is_unused_for_x",
            batch_key="",
            size_factor_key="total_counts",
            n_latent=2,
            max_epochs=3,
            random_seed=17,
        )
    )

    raw_model = FakeSCVI.instances[0]
    assert isinstance(trained_model, SCVIModel)
    assert trained_model.registered_adata is output
    assert raw_model.adata is output
    assert raw_model.is_trained is True
    assert raw_model.train_parameters == {"early_stopping": True, "max_epochs": 3}
    setup_adata, setup_kwargs = FakeSCVI.setup_calls[0]
    assert setup_adata is output
    assert setup_kwargs == {
        "layer": None,
        "batch_key": None,
        "size_factor_key": "total_counts",
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
    }
    assert raw_model.latent_calls == [(True, False)]
    assert trained_model.training_parameters["n_latent"] == 2
    assert trained_model.training_parameters["random_seed"] == 17
    assert trained_model.training_parameters["source"] == "X"
    assert trained_model.training_parameters["size_factor_key"] == "total_counts"
    assert list(trained_model.obs_names) == list(output.obs_names)
    assert list(trained_model.var_names) == list(output.var_names)
    assert output.obsm["X_scVI"].shape == (adata.n_obs, 2)
    assert "X_latent_qzm" not in output.obsm
    assert "X_latent_qzv" not in output.obsm
    assert "X_mde" not in output.obsm
    assert "X_scVI" not in adata.obsm

    mde_output, _ = output_values(
        OpenBioSingleCellSCVIIntegration.execute(
            adata,
            source="layer",
            counts_layer="counts",
            compute_mde=True,
            store_latent_distribution=False,
            n_latent=2,
        )
    )
    mde_model = FakeSCVI.instances[1]
    assert mde_model.latent_calls == [(True, False), (False, True)]
    assert "X_mde" in mde_output.obsm
    assert "X_latent_qzm" not in mde_output.obsm
    assert "X_latent_qzv" not in mde_output.obsm
    assert len(mde_calls) == 1

    distribution_output, _ = output_values(
        OpenBioSingleCellSCVIIntegration.execute(
            adata,
            source="layer",
            counts_layer="counts",
            compute_mde=False,
            store_latent_distribution=True,
            n_latent=2,
        )
    )
    distribution_model = FakeSCVI.instances[2]
    assert distribution_model.latent_calls == [(True, False), (False, True)]
    assert "X_mde" not in distribution_output.obsm
    assert distribution_output.obsm["X_latent_qzm"].shape == (adata.n_obs, 2)
    assert distribution_output.obsm["X_latent_qzv"].shape == (adata.n_obs, 2)


def test_scvi_differential_expression_reuses_model_and_aligns_group_subset(adata, science):
    registered = adata.copy()
    raw_model = FakeTrainedSCVI(registered, science)
    model = SCVIModel(raw_model, registered, {"n_latent": 2, "random_seed": 17})
    downstream = registered[["cell_3", "cell_0", "cell_2", "cell_1", "cell_5", "cell_4"], :].copy()
    downstream.obs["group"] = science.pd.Categorical(["B", "A", "B", "A", "B", "A"])
    downstream.obs["cohort"] = science.pd.Categorical(["keep", "keep", "keep", "keep", "drop", "drop"])
    downstream_obs = downstream.obs.copy(deep=True)
    registered_obs = registered.obs.copy(deep=True)

    result = output_values(
        OpenBioSingleCellSCVIDifferentialExpression.execute(
            downstream,
            model,
            groupby="group",
            group1="A",
            group2="B",
            subset_column="cohort",
            subset_value="keep",
            delta=0.5,
        )
    )[0]

    assert raw_model.train_calls == 0
    assert len(raw_model.differential_expression_calls) == 1
    call = raw_model.differential_expression_calls[0]
    assert list(call["adata"].obs_names) == ["cell_3", "cell_0", "cell_2", "cell_1"]
    assert call["adata"].obs["group"].astype(str).tolist() == ["B", "A", "B", "A"]
    assert call["groupby"] == "group"
    assert call["group1"] == "A"
    assert call["group2"] == "B"
    assert call["mode"] == "vanilla"
    assert call["delta"] == 0.5
    assert raw_model.deregister_manager_calls == [None]
    science.pd.testing.assert_frame_equal(downstream.obs, downstream_obs)
    science.pd.testing.assert_frame_equal(registered.obs, registered_obs)
    assert "_fake_internal_mutation" not in downstream.obs
    assert "_fake_internal_mutation" not in registered.obs
    assert list(result.table.columns) == ["gene", "lfc_mean", "proba_de"]
    assert result.input_cells == 4
    assert result.input_genes == 4
    assert result.parameters["mode"] == "vanilla"


def test_scvi_differential_expression_preserves_primary_error_when_cleanup_fails(adata, science):
    registered = adata.copy()
    raw_model = FakeTrainedSCVI(registered, science)
    raw_model.differential_expression_error = ValueError("primary differential-expression failure")
    raw_model.deregister_manager_error = RuntimeError("temporary manager cleanup failure")
    model = SCVIModel(raw_model, registered, {})
    downstream = registered.copy()
    downstream.obs["group"] = science.pd.Categorical(["A", "A", "A", "B", "B", "B"])

    with pytest.raises(ValueError, match="primary differential-expression failure"):
        OpenBioSingleCellSCVIDifferentialExpression.execute(downstream, model, "group", "A", "B")

    assert raw_model.deregister_manager_calls == [None]


def test_scvi_differential_expression_does_not_fail_after_successful_analysis(adata, science):
    class SCVI15TrainedModel(FakeTrainedSCVI):
        def deregister_manager(self, analysis_adata=None):
            self.deregister_manager_calls.append(analysis_adata)
            if analysis_adata is not None:
                raise ValueError("AnnData object was setup with a different model.")

    registered = adata.copy()
    raw_model = SCVI15TrainedModel(registered, science)
    model = SCVIModel(raw_model, registered, {})
    downstream = registered.copy()
    downstream.obs["group"] = science.pd.Categorical(["A", "A", "A", "B", "B", "B"])

    result = model.differential_expression(
        downstream,
        groupby="group",
        group1="A",
        group2="B",
    )

    assert result.index.tolist() == ["gene_0", "gene_1"]
    assert raw_model.deregister_manager_calls == [None]


def test_scvi_model_rejects_incompatible_observation_identity(adata, science):
    registered = adata.copy()
    model = SCVIModel(FakeTrainedSCVI(registered, science), registered, {})
    incompatible = registered.copy()
    incompatible.obs_names = [*incompatible.obs_names[:-1], "unregistered_cell"]
    incompatible.obs["group"] = science.pd.Categorical(["A", "A", "A", "B", "B", "B"])

    with pytest.raises(ValueError, match="obs_names are incompatible"):
        OpenBioSingleCellSCVIDifferentialExpression.execute(incompatible, model, "group", "A", "B")


def test_scvi_model_rejects_incompatible_variable_identity(adata, science):
    registered = adata.copy()
    model = SCVIModel(FakeTrainedSCVI(registered, science), registered, {})
    incompatible = registered[:, list(reversed(registered.var_names))].copy()
    incompatible.obs["group"] = science.pd.Categorical(["A", "A", "A", "B", "B", "B"])

    with pytest.raises(ValueError, match="var_names are incompatible"):
        OpenBioSingleCellSCVIDifferentialExpression.execute(incompatible, model, "group", "A", "B")


def test_scvi_model_rejects_untrained_or_detached_models(adata, science):
    untrained = FakeTrainedSCVI(adata, science)
    untrained.is_trained = False
    with pytest.raises(ValueError, match="requires a trained"):
        SCVIModel(untrained, adata, {})

    detached = FakeTrainedSCVI(adata.copy(), science)
    with pytest.raises(ValueError, match="attached to the registered"):
        SCVIModel(detached, adata, {})


def test_scvi_model_snapshots_training_parameters(adata, science):
    parameters = {"n_latent": 2, "categorical_covariates": ["batch"]}
    model = SCVIModel(FakeTrainedSCVI(adata, science), adata, parameters)
    parameters["categorical_covariates"].append("condition")
    visible_parameters = model.training_parameters
    visible_parameters["categorical_covariates"].append("donor")

    assert model.training_parameters["categorical_covariates"] == ["batch"]
