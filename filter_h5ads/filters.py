"""
Individual filter implementations.

Each function applies a specific filtering step to an AnnData object and returns
the filtered data along with statistics about the filtering operation.
"""

import logging
from typing import TYPE_CHECKING, Any

from anndata import AnnData

from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns

if TYPE_CHECKING:
    from filter_h5ads.config import (
        GeneDetectionFilterConfig,
        GuideFilterConfig,
        MitochondrialFilterConfig,
        UMIFilterConfig,
    )

logger = logging.getLogger(__name__)


def _ensure_in_memory(adata: AnnData) -> AnnData:
    """Ensure AnnData is in memory (not backed).

    Args:
        adata: AnnData object (possibly backed)

    Returns:
        AnnData object in memory
    """
    if adata.isbacked:
        logger.debug("Converting backed AnnData to memory")
        return adata.to_memory()
    return adata


def apply_umi_filter(
    adata: AnnData, config: "UMIFilterConfig"
) -> tuple[AnnData, dict[str, Any]]:
    """Apply UMI count filter to remove low-quality cells.

    Args:
        adata: Input AnnData object
        config: UMI filter configuration

    Returns:
        Tuple of (filtered AnnData, statistics dictionary)

    Raises:
        ValueError: If required column is missing
    """
    logger.info(f"Applying UMI filter: min_counts={config.min_counts}")

    # Validate required columns
    validate_h5ad_columns(adata, [config.count_column], location="obs")

    # Ensure data is in memory for processing
    adata = _ensure_in_memory(adata)

    # Store original for stats
    adata_before = adata

    # Apply filter
    mask = adata.obs[config.count_column] >= config.min_counts
    adata_filtered = adata[mask].copy()

    # Calculate statistics
    stats = calculate_filter_stats(
        adata_before=adata_before, adata_after=adata_filtered, filter_name="UMI Count Filter"
    )

    # Explicitly delete reference to help garbage collection
    del adata_before, adata

    # Add filter-specific info
    stats["min_counts_threshold"] = config.min_counts
    stats["count_column"] = config.count_column

    logger.info(
        f"UMI filter: {stats['n_cells_after']:,} / {stats['n_cells_before']:,} cells "
        f"retained ({stats['pct_retained']:.1f}%)"
    )

    return adata_filtered, stats


def apply_guide_filter(
    adata: AnnData, config: "GuideFilterConfig"
) -> tuple[AnnData, dict[str, Any]]:
    """Apply guide quality filter to keep only correctly perturbed cells.

    Args:
        adata: Input AnnData object
        config: Guide filter configuration

    Returns:
        Tuple of (filtered AnnData, statistics dictionary)

    Raises:
        ValueError: If required column is missing
    """
    logger.info(f"Applying guide filter: {config.guide_column}={config.required_value}")

    # Validate required columns
    validate_h5ad_columns(adata, [config.guide_column], location="obs")

    # Ensure data is in memory for processing
    adata = _ensure_in_memory(adata)

    # Store original for stats
    adata_before = adata

    # Apply filter - handle both boolean and integer columns
    guide_col = adata.obs[config.guide_column]

    if config.required_value in (True, 1):
        # Be flexible with boolean/int comparison
        mask = guide_col.astype(bool) | (guide_col == 1)
    elif config.required_value in (False, 0):
        mask = (~guide_col.astype(bool)) | (guide_col == 0)
    else:
        mask = guide_col == config.required_value

    adata_filtered = adata[mask].copy()

    # Calculate statistics
    stats = calculate_filter_stats(
        adata_before=adata_before, adata_after=adata_filtered, filter_name="Guide Quality Filter"
    )

    # Explicitly delete reference to help garbage collection
    del adata_before, adata

    # Add filter-specific info
    stats["guide_column"] = config.guide_column
    stats["required_value"] = config.required_value

    logger.info(
        f"Guide filter: {stats['n_cells_after']:,} / {stats['n_cells_before']:,} cells "
        f"retained ({stats['pct_retained']:.1f}%)"
    )

    return adata_filtered, stats


def apply_mito_filter(
    adata: AnnData, config: "MitochondrialFilterConfig"
) -> tuple[AnnData, dict[str, Any]]:
    """Apply mitochondrial content filter to remove dead/dying cells.

    Args:
        adata: Input AnnData object
        config: Mitochondrial filter configuration

    Returns:
        Tuple of (filtered AnnData, statistics dictionary)

    Raises:
        ValueError: If required column is missing
    """
    logger.info(f"Applying mitochondrial filter: max_pct_mt={config.max_pct_mt}%")

    # Validate required columns
    validate_h5ad_columns(adata, [config.mt_column], location="obs")

    # Ensure data is in memory for processing
    adata = _ensure_in_memory(adata)

    # Store original for stats
    adata_before = adata

    # Apply filter
    mask = adata.obs[config.mt_column] <= config.max_pct_mt
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
    stats["max_pct_mt_threshold"] = config.max_pct_mt
    stats["mt_column"] = config.mt_column

    logger.info(
        f"Mito filter: {stats['n_cells_after']:,} / {stats['n_cells_before']:,} cells "
        f"retained ({stats['pct_retained']:.1f}%)"
    )

    return adata_filtered, stats


def apply_gene_filter(
    adata: AnnData, config: "GeneDetectionFilterConfig"
) -> tuple[AnnData, dict[str, Any]]:
    """Apply gene detection filter to remove low-quality cells.

    Args:
        adata: Input AnnData object
        config: Gene detection filter configuration

    Returns:
        Tuple of (filtered AnnData, statistics dictionary)

    Raises:
        ValueError: If required column is missing
    """
    logger.info(f"Applying gene detection filter: min_genes={config.min_genes}")

    # Validate required columns
    validate_h5ad_columns(adata, [config.gene_column], location="obs")

    # Ensure data is in memory for processing
    adata = _ensure_in_memory(adata)

    # Store original for stats
    adata_before = adata

    # Apply filter
    mask = adata.obs[config.gene_column] >= config.min_genes
    adata_filtered = adata[mask].copy()

    # Calculate statistics
    stats = calculate_filter_stats(
        adata_before=adata_before, adata_after=adata_filtered, filter_name="Gene Detection Filter"
    )

    # Explicitly delete reference to help garbage collection
    del adata_before, adata

    # Add filter-specific info
    stats["min_genes_threshold"] = config.min_genes
    stats["gene_column"] = config.gene_column

    logger.info(
        f"Gene filter: {stats['n_cells_after']:,} / {stats['n_cells_before']:,} cells "
        f"retained ({stats['pct_retained']:.1f}%)"
    )

    return adata_filtered, stats
