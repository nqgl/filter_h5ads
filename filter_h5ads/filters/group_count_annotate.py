"""Group count annotation step.

Adds a `.obs` column that records the number of cells in each group defined by
a set of `.obs` keys (e.g., for each (cell_line, drug, dosage) condition).
"""

from __future__ import annotations

from typing import Any, Literal

import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import ConfigDict, Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class GroupCountAnnotateConfig(FilterStepConfig):
    """Annotate each cell with its group size."""

    groupby: list[str] = Field(
        min_length=1,
        description="`.obs` columns defining each condition/group.",
    )
    output_column: str = Field(
        default="n_cells_in_group",
        description="Output `.obs` column to write group sizes into.",
    )
    overwrite: bool = Field(
        default=False,
        description="If True, allows overwriting an existing output_column.",
    )
    dropna: bool = Field(
        default=False,
        description=(
            "If True, rows with NA in any groupby key are not counted and get NA in "
            "the output column."
        ),
    )
    fillna: int | None = Field(
        default=None,
        description="If set, fills NA group sizes with this value after computation.",
    )
    dtype: Literal["Int64", "int64"] = Field(
        default="Int64",
        description=(
            "Pandas integer dtype for the output. Use 'Int64' to allow NA values."
        ),
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        logger.info(
            f"Applying group count annotation: groupby={self.groupby} "
            f"-> {self.output_column}"
        )
        validate_h5ad_columns(adata, self.groupby, location="obs")
        adata = _ensure_in_memory(adata)
        if getattr(adata, "is_view", False):
            adata = adata.copy()

        if self.output_column in adata.obs.columns and not self.overwrite:
            raise ValueError(
                f"Output column '{self.output_column}' already exists in .obs. "
                "Set overwrite=True to replace it."
            )

        obs = adata.obs
        assert isinstance(obs, pd.DataFrame)

        grouped = obs.groupby(self.groupby, sort=False, dropna=self.dropna, observed=True)
        # Use any series to compute per-row group sizes.
        base_col = self.groupby[0]
        group_sizes = grouped[base_col].transform("size")
        assert isinstance(group_sizes, pd.Series)

        if self.fillna is not None:
            group_sizes = group_sizes.fillna(self.fillna)

        if self.dtype == "Int64":
            group_sizes = group_sizes.astype("Int64")
        else:
            group_sizes = group_sizes.astype("int64")

        adata.obs[self.output_column] = group_sizes

        stats = calculate_filter_stats(
            adata_before=adata,
            adata_after=adata,
            filter_name="Group Count Annotate",
        )
        stats["groupby"] = list(self.groupby)
        stats["output_column"] = self.output_column
        stats["overwrite"] = self.overwrite
        stats["dropna"] = self.dropna
        stats["fillna"] = self.fillna
        stats["dtype"] = self.dtype

        return adata, stats

