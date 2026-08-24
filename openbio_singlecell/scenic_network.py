from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pandas import DataFrame
else:
    DataFrame = Any


REQUIRED_ADJACENCY_COLUMNS = ("TF", "target", "importance")


@dataclass(frozen=True, slots=True)
class ScenicNetwork:
    adjacency: DataFrame
    gene_names: tuple[str, ...]
    grn_method: str
    provenance: dict[str, Any]

    def __post_init__(self) -> None:
        from pandas import DataFrame, to_numeric

        if not isinstance(self.adjacency, DataFrame):
            raise TypeError("ScenicNetwork adjacency must be a pandas DataFrame.")
        if self.adjacency.empty:
            raise ValueError("ScenicNetwork adjacency cannot be empty.")

        missing = [column for column in REQUIRED_ADJACENCY_COLUMNS if column not in self.adjacency.columns]
        if missing:
            raise ValueError(f"ScenicNetwork adjacency is missing required columns: {missing}")

        if isinstance(self.gene_names, str):
            raise TypeError("ScenicNetwork gene names must be an ordered collection of strings.")
        try:
            gene_names = tuple(self.gene_names)
        except TypeError as exc:
            raise TypeError("ScenicNetwork gene names must be an ordered collection of strings.") from exc
        if not gene_names:
            raise ValueError("ScenicNetwork gene names cannot be empty.")
        if any(not isinstance(name, str) or not name.strip() for name in gene_names):
            raise ValueError("ScenicNetwork gene names must be non-empty strings.")
        gene_names = tuple(name.strip() for name in gene_names)
        if len(set(gene_names)) != len(gene_names):
            raise ValueError("ScenicNetwork gene names must be unique.")

        adjacency = self.adjacency.copy(deep=True)
        for column in ("TF", "target"):
            labels = adjacency[column].astype("string")
            if labels.isna().any() or labels.str.strip().eq("").any():
                raise ValueError(f"ScenicNetwork adjacency column {column!r} cannot contain empty values.")
            adjacency[column] = labels.str.strip().astype(str)

        inference_genes = set(gene_names)
        unknown_genes = sorted(set(adjacency["TF"]).union(adjacency["target"]).difference(inference_genes))
        if unknown_genes:
            raise ValueError(
                f"ScenicNetwork adjacency references genes outside the inference gene set: {unknown_genes}"
            )

        importance = to_numeric(adjacency["importance"], errors="coerce")
        if importance.isna().any() or not importance.map(math.isfinite).all():
            raise ValueError("ScenicNetwork adjacency importance values must be finite numbers.")
        adjacency["importance"] = importance

        if not isinstance(self.grn_method, str) or not self.grn_method.strip():
            raise ValueError("ScenicNetwork GRN method cannot be empty.")
        if not isinstance(self.provenance, Mapping) or not self.provenance:
            raise ValueError("ScenicNetwork provenance cannot be empty.")
        if any(not isinstance(key, str) or not key.strip() for key in self.provenance):
            raise ValueError("ScenicNetwork provenance keys must be non-empty strings.")

        object.__setattr__(self, "adjacency", adjacency)
        object.__setattr__(self, "gene_names", gene_names)
        object.__setattr__(self, "grn_method", self.grn_method.strip())
        object.__setattr__(self, "provenance", copy.deepcopy(dict(self.provenance)))


__all__ = ["REQUIRED_ADJACENCY_COLUMNS", "ScenicNetwork"]
