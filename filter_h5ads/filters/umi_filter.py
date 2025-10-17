"""UMI count filter for removing low-quality cells."""

from typing import Any

import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import validate_h5ad_columns


class UMIFilterConfig(FilterStepConfig):
    """Configuration for UMI count filtering.

    Filters cells based on total UMI counts to remove low-quality cells.
    """

    min_counts: int = Field(
        default=15000,
        ge=0,
        description="Minimum total UMI counts required per cell",
    )
    count_column: str = Field(
        default="total_counts",
        description="Name of the column in .obs containing UMI counts",
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply UMI count filter to remove low-quality cells."""
        logger.info(f"Applying UMI filter: min_counts={self.min_counts}")

        # Validate required columns
        validate_h5ad_columns(adata, [self.count_column], location="obs")

        # Ensure data is in memory for processing
        adata = _ensure_in_memory(adata)

        # Capture stats before filtering
        n_cells_before = adata.n_obs
        if "gene_target" in adata.obs.columns:
            gene_target_col = adata.obs["gene_target"]
            assert isinstance(gene_target_col, pd.Series)
            breakdown_before = gene_target_col.value_counts().to_dict()
        else:
            breakdown_before = None

        # Apply filter
        mask = adata.obs[self.count_column] >= self.min_counts
        adata_filtered = adata[mask].copy()

        # Explicitly delete reference to help garbage collection
        del adata

        # Calculate statistics
        n_cells_after = adata_filtered.n_obs
        n_removed = n_cells_before - n_cells_after
        pct_retained = (
            (n_cells_after / n_cells_before * 100) if n_cells_before > 0 else 0.0
        )

        stats = {
            "filter_name": "UMI Count Filter",
            "n_cells_before": n_cells_before,
            "n_cells_after": n_cells_after,
            "n_cells_removed": n_removed,
            "pct_retained": round(pct_retained, 2),
            "min_counts_threshold": self.min_counts,
            "count_column": self.count_column,
        }

        if breakdown_before is not None:
            gene_target_col_after = adata_filtered.obs["gene_target"]
            assert isinstance(gene_target_col_after, pd.Series)
            breakdown_after = gene_target_col_after.value_counts().to_dict()
            breakdown_removed = {
                key: breakdown_before.get(key, 0) - breakdown_after.get(key, 0)
                for key in breakdown_before.keys()
            }
            stats["breakdown_before"] = breakdown_before
            stats["breakdown_after"] = breakdown_after
            stats["breakdown_removed"] = breakdown_removed

        logger.info(
            f"UMI filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells "
            f"retained ({stats['pct_retained']:.1f}%)",
        )

        return adata_filtered, stats
