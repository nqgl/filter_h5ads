"""Mitochondrial content filter implementation.

Filters cells with high mitochondrial percentage, which typically indicates
dead or dying cells in single-cell experiments.
"""

import logging
from typing import Any

from anndata import AnnData
from pydantic import Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns

logger = logging.getLogger(__name__)


class MitochondrialFilterConfig(FilterStepConfig):
    """Configuration for mitochondrial content filtering.

    Filters cells with high mitochondrial percentage (dead/dying cells).
    """

    max_pct_mt: float = Field(
        default=20.0,
        ge=0.0,
        le=100.0,
        description="Maximum percentage of mitochondrial UMIs allowed",
    )
    mt_column: str = Field(
        default="pct_counts_mt",
        description="Name of the column in .obs containing MT percentage",
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply mitochondrial content filter to remove dead/dying cells.

        Args:
            adata: Input AnnData object

        Returns:
            Tuple of (filtered AnnData, statistics dictionary)

        Raises:
            ValueError: If required column is missing

        """
        logger.info(f"Applying mitochondrial filter: max_pct_mt={self.max_pct_mt}%")

        # Validate required columns
        validate_h5ad_columns(adata, [self.mt_column], location="obs")

        # Ensure data is in memory for processing
        adata = _ensure_in_memory(adata)

        # Store original for stats
        adata_before = adata

        # Apply filter
        mask = adata.obs[self.mt_column] <= self.max_pct_mt
        adata_filtered = adata[mask].copy()

        # Calculate statistics
        stats = calculate_filter_stats(
            adata_before=adata_before,
            adata_after=adata_filtered,
            filter_name="Mitochondrial Content Filter",
        )

        # Explicitly delete reference to help garbage collection
        del adata_before, adata

        # Add filter-specific info
        stats["max_pct_mt_threshold"] = self.max_pct_mt
        stats["mt_column"] = self.mt_column

        logger.info(
            f"Mito filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells "
            f"retained ({stats['pct_retained']:.1f}%)",
        )

        return adata_filtered, stats
