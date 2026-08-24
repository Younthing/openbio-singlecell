from __future__ import annotations

import copy
import sys
from collections.abc import Mapping, Sequence
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


class SCVIModel:
    """A trained scVI model bound to the AnnData identity used for training."""

    __slots__ = (
        "_model",
        "_model_class",
        "_obs_name_set",
        "_obs_names",
        "_registered_adata",
        "_training_parameters",
        "_var_names",
    )

    def __init__(
        self,
        model: Any,
        registered_adata: AnnData,
        training_parameters: Mapping[str, Any],
    ) -> None:
        if not callable(getattr(model, "differential_expression", None)):
            raise TypeError("SCVIModel requires a scVI model with differential-expression support.")
        if not callable(getattr(model, "deregister_manager", None)):
            raise TypeError("SCVIModel requires a scVI model with AnnData manager lifecycle support.")
        if getattr(model, "adata", None) is not registered_adata:
            raise ValueError("The scVI model must be attached to the registered AnnData object.")
        if not bool(getattr(model, "is_trained", False)):
            raise ValueError("SCVIModel requires a trained scVI model.")
        if not isinstance(training_parameters, Mapping):
            raise TypeError("SCVIModel training_parameters must be a mapping.")

        try:
            obs_names = tuple(registered_adata.obs_names)
            var_names = tuple(registered_adata.var_names)
            obs_names_are_unique = bool(registered_adata.obs_names.is_unique)
            var_names_are_unique = bool(registered_adata.var_names.is_unique)
        except AttributeError as error:
            raise TypeError("SCVIModel requires a registered AnnData object.") from error
        if not obs_names or not var_names:
            raise ValueError("The registered scVI AnnData must contain observations and variables.")
        if not obs_names_are_unique or not var_names_are_unique:
            raise ValueError("The registered scVI AnnData requires unique obs_names and var_names.")

        self._model = model
        self._registered_adata = registered_adata
        self._training_parameters = copy.deepcopy(dict(training_parameters))
        self._obs_names = obs_names
        self._obs_name_set = frozenset(obs_names)
        self._var_names = var_names
        model_type = type(model)
        self._model_class = f"{model_type.__module__}.{model_type.__qualname__}"

    @property
    def registered_adata(self) -> AnnData:
        return self._registered_adata

    @property
    def training_parameters(self) -> Mapping[str, Any]:
        return MappingProxyType(copy.deepcopy(self._training_parameters))

    @property
    def model_class(self) -> str:
        return self._model_class

    @property
    def obs_names(self) -> tuple[str, ...]:
        return self._obs_names

    @property
    def var_names(self) -> tuple[str, ...]:
        return self._var_names

    def differential_expression(
        self,
        adata: AnnData,
        *,
        groupby: str,
        group1: str,
        group2: str,
        subset_column: str = "",
        subset_value: str = "",
        mode: str = "vanilla",
        delta: float = 0.25,
    ) -> DataFrame:
        groupby = groupby.strip()
        group1 = group1.strip()
        group2 = group2.strip()
        subset_column = subset_column.strip()
        subset_value = subset_value.strip()
        mode = mode.strip()
        if not groupby:
            raise ValueError("scVI differential-expression groupby cannot be empty.")
        if not group1 or not group2:
            raise ValueError("scVI differential-expression group1 and group2 values are required.")
        if group1 == group2:
            raise ValueError("scVI differential-expression group1 and group2 must be different.")
        if bool(subset_column) != bool(subset_value):
            raise ValueError("scVI subset_column and subset_value must be provided together.")
        if mode not in {"vanilla", "change"}:
            raise ValueError(f"Unsupported scVI differential-expression mode: {mode!r}")

        self._validate_registered_state()
        required_columns = list(dict.fromkeys([groupby, subset_column] if subset_column else [groupby]))
        analysis_adata = self._copy_compatible_adata(adata, required_columns)
        if subset_column:
            mask = analysis_adata.obs[subset_column].astype(str) == subset_value
            analysis_adata = analysis_adata[mask].copy()
            if analysis_adata.n_obs == 0:
                raise ValueError(f"scVI subset contains no observations: {subset_column}={subset_value!r}")

        observed_groups = set(analysis_adata.obs[groupby].dropna().astype(str))
        missing_groups = [group for group in (group1, group2) if group not in observed_groups]
        if missing_groups:
            raise ValueError(f"scVI differential-expression groups not found after subsetting: {missing_groups}")

        try:
            return self._model.differential_expression(
                adata=analysis_adata,
                groupby=groupby,
                group1=group1,
                group2=group2,
                mode=mode,
                delta=delta,
            )
        finally:
            differential_error = sys.exception()
            try:
                self._model.deregister_manager()
            except Exception:
                if differential_error is None:
                    raise

    def _copy_compatible_adata(self, adata: AnnData, required_obs_columns: Sequence[str]) -> AnnData:
        try:
            downstream_obs_names = tuple(adata.obs_names)
            downstream_var_names = tuple(adata.var_names)
            obs_names_are_unique = bool(adata.obs_names.is_unique)
            var_names_are_unique = bool(adata.var_names.is_unique)
        except AttributeError as error:
            raise TypeError("SCVIModel requires a downstream AnnData object.") from error
        if not downstream_obs_names:
            raise ValueError("Downstream AnnData must contain at least one observation.")
        if not obs_names_are_unique or not var_names_are_unique:
            raise ValueError("Downstream AnnData requires unique obs_names and var_names.")
        if downstream_var_names != self._var_names:
            raise ValueError("Downstream AnnData var_names are incompatible with the trained scVI model.")
        unknown_obs_names = [name for name in downstream_obs_names if name not in self._obs_name_set]
        if unknown_obs_names:
            raise ValueError(
                "Downstream AnnData obs_names are incompatible with the trained scVI model; "
                f"first unknown observation: {unknown_obs_names[0]!r}."
            )
        missing_columns = [column for column in required_obs_columns if column not in adata.obs]
        if missing_columns:
            raise ValueError(f"scVI differential-expression observation columns not found: {missing_columns}")

        analysis_adata = self._registered_adata[list(downstream_obs_names), :].copy()
        for column in required_obs_columns:
            analysis_adata.obs[column] = adata.obs[column].reindex(analysis_adata.obs_names).copy()
        return analysis_adata

    def _validate_registered_state(self) -> None:
        if getattr(self._model, "adata", None) is not self._registered_adata:
            raise RuntimeError("The scVI model is no longer attached to its registered AnnData object.")
        if not bool(getattr(self._model, "is_trained", False)):
            raise RuntimeError("The scVI model is no longer in a trained state.")
        if tuple(self._registered_adata.obs_names) != self._obs_names:
            raise RuntimeError("The registered scVI AnnData obs_names changed after training.")
        if tuple(self._registered_adata.var_names) != self._var_names:
            raise RuntimeError("The registered scVI AnnData var_names changed after training.")


__all__ = ["SCVIModel"]
