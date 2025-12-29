"""Pseudobulking step.

Groups cells by one or more `.obs` columns and sums counts to create a pooled
"pseudobulk" observation per group. `.obs` metadata columns can be aggregated
per group via a safe, serializable specification.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Annotated, Any, Literal, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class _ObsAgg(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    agg: str
    output_column: str | None = Field(
        default=None,
        description="If provided, writes this aggregation under a different output name.",
    )

    def get_output_column(self, input_column: str) -> str:
        return self.output_column or input_column

    def run(self, series: pd.Series) -> Any:  # pragma: no cover - exercised via subclasses
        raise NotImplementedError


class DropAgg(_ObsAgg):
    agg: Literal["drop"] = "drop"

    def run(self, series: pd.Series) -> Any:
        raise RuntimeError("drop aggregation should not be executed")


class ListAgg(_ObsAgg):
    agg: Literal["list"] = "list"
    unique: bool = Field(
        default=False,
        description="If True, keeps unique values (preserving first-seen order).",
    )

    def run(self, series: pd.Series) -> Any:
        if self.unique:
            seen: set[tuple[str, Any]] = set()
            out: list[Any] = []
            for v in series.tolist():
                try:
                    hash(v)
                    key = ("hashable", v)
                except TypeError:
                    key = ("repr", repr(v))
                if key not in seen:
                    seen.add(key)
                    out.append(v)
            return out
        return series.tolist()


class SumAgg(_ObsAgg):
    agg: Literal["sum"] = "sum"

    def run(self, series: pd.Series) -> Any:
        numeric = pd.to_numeric(series, errors="raise")
        return float(numeric.sum(min_count=1))


class MeanAgg(_ObsAgg):
    agg: Literal["mean"] = "mean"

    def run(self, series: pd.Series) -> Any:
        numeric = pd.to_numeric(series, errors="raise")
        return float(numeric.mean())


class AssertConstantAgg(_ObsAgg):
    agg: Literal["assert_constant"] = "assert_constant"
    allow_null: bool = Field(
        default=True,
        description=(
            "If True, allows all-null groups and returns null. If False, all-null "
            "groups raise."
        ),
    )

    def run(self, series: pd.Series) -> Any:
        # Treat NA consistently
        non_null = series.dropna()
        if len(non_null) == 0:
            if self.allow_null:
                return None
            raise ValueError("assert_constant: group contains only null values")

        first = non_null.iloc[0]
        if not (non_null == first).all():
            examples = non_null.unique().tolist()
            raise ValueError(
                f"assert_constant failed: expected a single value, got {examples}"
            )
        return first


class ProportionsAgg(_ObsAgg):
    agg: Literal["proportions"] = "proportions"
    key_as_str: bool = Field(
        default=True,
        description="If True, converts unique values to string keys in the output dict.",
    )
    include_null: bool = Field(
        default=False,
        description="If True, includes null values as a key in the output dict.",
    )

    def run(self, series: pd.Series) -> Any:
        if self.include_null:
            vc = series.value_counts(dropna=False)
        else:
            vc = series.value_counts(dropna=True)
        total = float(vc.sum())
        out: dict[Any, float] = {}
        for k, v in vc.items():
            key = str(k) if self.key_as_str else k
            out[key] = float(v) / total if total > 0 else 0.0
        return out


ObsAgg = Annotated[
    Union[
        DropAgg,
        ListAgg,
        SumAgg,
        MeanAgg,
        AssertConstantAgg,
        ProportionsAgg,
    ],
    Field(discriminator="agg"),
]


def _sum_by_group(X: Any, group_ids: np.ndarray, n_groups: int) -> Any:
    n_obs = group_ids.shape[0]
    if sp.issparse(X):
        X_csr = X.tocsr()
        rows = group_ids.astype(np.int64, copy=False)
        cols = np.arange(n_obs, dtype=np.int64)
        data = np.ones(n_obs, dtype=np.int8)
        G = sp.csr_matrix((data, (rows, cols)), shape=(n_groups, n_obs))
        return G @ X_csr

    X_np = np.asarray(X)
    if X_np.ndim != 2:
        raise ValueError(f"Expected 2D matrix for X, got shape {X_np.shape}")
    out = np.zeros((n_groups, X_np.shape[1]), dtype=X_np.dtype)
    np.add.at(out, group_ids.astype(np.int64, copy=False), X_np)
    return out


class PseudobulkConfig(FilterStepConfig):
    """Group cells and sum counts to create pseudobulk samples."""

    groupby: list[str] = Field(
        min_length=1,
        description="`.obs` columns whose unique combinations define pooled samples.",
    )
    layer: str | None = Field(
        default=None,
        description=(
            "If provided, sums counts from `.layers[layer]` instead of `.X`."
        ),
    )
    dropna: bool = Field(
        default=False,
        description="If True, drop groups with NA in any groupby key.",
    )
    n_cells_column: str = Field(
        default="n_cells_pooled",
        description="Output `.obs` column recording number of pooled cells per group.",
    )
    obs_aggregations: dict[str, ObsAgg] = Field(
        default_factory=dict,
        description=(
            "Mapping from input `.obs` column name -> aggregation specification. "
            "Unspecified columns are dropped."
        ),
    )
    require_explicit_obs_aggregation: bool = Field(
        default=False,
        description=(
            "If True, requires every non-groupby `.obs` column to appear in "
            "obs_aggregations."
        ),
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        logger.info(f"Applying pseudobulk: groupby={self.groupby}")
        validate_h5ad_columns(adata, self.groupby, location="obs")

        if self.layer is not None and self.layer not in adata.layers:
            raise ValueError(
                f"Layer '{self.layer}' not found in AnnData.layers. "
                f"Available layers: {list(adata.layers.keys())}"
            )

        adata = _ensure_in_memory(adata)

        obs = adata.obs
        assert isinstance(obs, pd.DataFrame)

        missing_specs: list[str] = []
        if self.require_explicit_obs_aggregation:
            for col in obs.columns:
                if col in self.groupby:
                    continue
                if col not in self.obs_aggregations:
                    missing_specs.append(col)
        if missing_specs:
            raise ValueError(
                "Missing obs_aggregations for columns: "
                f"{missing_specs}. Set require_explicit_obs_aggregation=False to drop "
                "unspecified columns."
            )

        dropna = True if self.dropna else False
        grouped = obs.groupby(self.groupby, sort=False, dropna=dropna, observed=True)
        group_ids = grouped.ngroup().to_numpy(dtype=np.int64, copy=False)
        n_groups = int(group_ids.max() + 1) if len(group_ids) > 0 else 0

        if n_groups == 0:
            raise ValueError("Pseudobulk produced 0 groups (check groupby/dropna).")

        X_src = adata.layers[self.layer] if self.layer is not None else adata.X
        X_sum = _sum_by_group(X_src, group_ids, n_groups)

        # Start output obs with group keys (from group index) + pooled cell counts.
        group_sizes = grouped.size()
        out_obs = group_sizes.index.to_frame(index=False)
        out_obs[self.n_cells_column] = group_sizes.to_numpy(dtype=np.int64, copy=False)

        # Aggregate requested metadata columns
        group_names = list(group_sizes.index)

        def per_group_apply(fn: Callable[[pd.Series], Any], col: str) -> list[Any]:
            return [fn(obs.iloc[grouped.indices[name]][col]) for name in group_names]

        for input_col, agg_spec in self.obs_aggregations.items():
            if input_col in self.groupby:
                continue
            if agg_spec.agg == "drop":
                continue
            if input_col not in obs.columns:
                raise ValueError(
                    f"obs_aggregations references missing .obs column '{input_col}'"
                )
            output_col = agg_spec.get_output_column(input_col)
            if output_col in out_obs.columns:
                raise ValueError(
                    f"Duplicate output .obs column '{output_col}' in pseudobulk output."
                )
            out_obs[output_col] = per_group_apply(agg_spec.run, input_col)

        out_obs.index = [f"pb_{i}" for i in range(n_groups)]
        out_var = adata.var.copy()
        out = AnnData(X=X_sum, obs=out_obs, var=out_var, uns=adata.uns.copy())

        stats = calculate_filter_stats(
            adata_before=adata,
            adata_after=out,
            filter_name="Pseudobulk",
        )
        stats["groupby"] = list(self.groupby)
        stats["layer"] = self.layer
        stats["n_groups"] = int(n_groups)
        stats["n_cells_column"] = self.n_cells_column
        stats["dropna"] = self.dropna
        stats["aggregated_obs_columns"] = {
            k: v.get_output_column(k)
            for k, v in self.obs_aggregations.items()
            if v.agg != "drop"
        }

        return out, stats
