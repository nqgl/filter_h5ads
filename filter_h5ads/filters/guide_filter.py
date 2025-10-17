"""Guide quality filter implementation.

Filters cells based on whether they passed guide pairing QC in CRISPR screens.
"""

from typing import Any

import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class GuideFilterConfig(FilterStepConfig):
    """Configuration for guide quality filtering.

    Filters cells based on whether they passed guide pairing QC.
    """

    guide_column: str = Field(
        default="pass_guide_filter",
        description="Name of the column in .obs containing guide filter status",
    )
    required_value: int | bool = Field(
        default=1,
        description=(
            "Required value in guide_column for cells to pass (typically 1 or True)"
        ),
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply guide quality filter to keep only correctly perturbed cells.

        Args:
            adata: Input AnnData object

        Returns:
            Tuple of (filtered AnnData, statistics dictionary)

        Raises:
            ValueError: If required column is missing

        """
        logger.info(f"Applying guide filter: {self.guide_column}={self.required_value}")

        # Validate required columns
        validate_h5ad_columns(adata, [self.guide_column], location="obs")

        # Ensure data is in memory for processing
        adata = _ensure_in_memory(adata)

        # Store original for stats
        adata_before = adata

        # Apply filter - handle both boolean and integer columns
        guide_col = adata.obs[self.guide_column]
        assert isinstance(guide_col, pd.Series)

        if self.required_value in (True, 1):
            # Be flexible with boolean/int comparison
            mask = guide_col.astype(bool) | (guide_col == 1)
        elif self.required_value in (False, 0):
            mask = (~guide_col.astype(bool)) | (guide_col == 0)
        else:
            mask = guide_col == self.required_value

        adata_filtered = adata[mask].copy()

        # Calculate statistics
        stats = calculate_filter_stats(
            adata_before=adata_before,
            adata_after=adata_filtered,
            filter_name="Guide Quality Filter",
        )

        # Explicitly delete reference to help garbage collection
        del adata_before, adata

        # Add filter-specific info
        stats["guide_column"] = self.guide_column
        stats["required_value"] = self.required_value

        logger.info(
            f"Guide filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells "
            f"retained ({stats['pct_retained']:.1f}%)",
        )

        return adata_filtered, stats
