"""Gene detection filter implementation.

Filters cells with too few detected genes, which typically indicates
low-quality cells in single-cell experiments.
"""

import logging
from typing import Any

from anndata import AnnData
from pydantic import Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns

logger = logging.getLogger(__name__)


class GeneDetectionFilterConfig(FilterStepConfig):
    """Configuration for gene detection filtering.

    Filters cells with too few detected genes (low-quality cells).
    """

    min_genes: int = Field(
        default=200,
        ge=0,
        description="Minimum number of genes with non-zero counts required",
    )
    gene_column: str = Field(
        default="n_genes_by_counts",
        description="Name of the column in .obs containing gene counts",
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply gene detection filter to remove low-quality cells.

        Args:
            adata: Input AnnData object

        Returns:
            Tuple of (filtered AnnData, statistics dictionary)

        Raises:
            ValueError: If required column is missing

        """
        logger.info(f"Applying gene detection filter: min_genes={self.min_genes}")

        # Validate required columns
        validate_h5ad_columns(adata, [self.gene_column], location="obs")

        # Ensure data is in memory for processing
        adata = _ensure_in_memory(adata)

        # Store original for stats
        adata_before = adata

        # Apply filter
        mask = adata.obs[self.gene_column] >= self.min_genes
        adata_filtered = adata[mask].copy()

        # Calculate statistics
        stats = calculate_filter_stats(
            adata_before=adata_before,
            adata_after=adata_filtered,
            filter_name="Gene Detection Filter",
        )

        # Explicitly delete reference to help garbage collection
        del adata_before, adata

        # Add filter-specific info
        stats["min_genes_threshold"] = self.min_genes
        stats["gene_column"] = self.gene_column

        logger.info(
            f"Gene filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells "
            f"retained ({stats['pct_retained']:.1f}%)",
        )

        return adata_filtered, stats
