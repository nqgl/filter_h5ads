"""Obs value filter implementation.

Filters cells based on whether an .obs column has one of a set of values.
This is useful for including or excluding specific cell types, batches, etc.
"""

from typing import Any

import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import AliasChoices, Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class ObsValueFilterConfig(FilterStepConfig):
    """Configuration for filtering cells by an .obs column's value."""

    key: str = Field(description="Name of the column in .obs to check")
    values: list[str | int | bool] = Field(
        min_length=1,
        description="List of values to include or exclude",
    )
    exclude: bool = Field(
        default=True,
        validation_alias=AliasChoices("exclude", "exclue"),
        description=(
            "If True (default), remove cells where key is in values. "
            "If False, keep only cells where key is in values."
        ),
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply .obs value filter to include or exclude cells."""
        action = "excluding" if self.exclude else "including"
        logger.info(f"Applying obs value filter: {action} {self.key} in {self.values}")

        validate_h5ad_columns(adata, [self.key], location="obs")
        adata = _ensure_in_memory(adata)

        adata_before = adata

        series = adata.obs[self.key]
        assert isinstance(series, pd.Series)
        matches = series.isin(self.values)

        if self.exclude:
            keep_mask = ~matches
        else:
            keep_mask = matches

        adata_filtered = adata[keep_mask].copy()

        stats = calculate_filter_stats(
            adata_before=adata_before,
            adata_after=adata_filtered,
            filter_name="Obs Value Filter",
        )

        del adata_before, adata

        stats["key"] = self.key
        stats["values"] = list(self.values)
        stats["exclude"] = self.exclude
        stats["n_matching_before"] = int(matches.sum())

        logger.info(
            f"Obs value filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells "
            f"retained ({stats['pct_retained']:.1f}%)",
        )

        return adata_filtered, stats

