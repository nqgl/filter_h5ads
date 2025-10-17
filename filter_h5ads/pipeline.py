"""Pipeline orchestration for running multiple filters in sequence.

Provides functions to execute a complete filtering pipeline and track
provenance information.
"""

import logging
from datetime import datetime
from typing import Any

from anndata import AnnData

from filter_h5ads.config import FilterPipelineConfig

logger = logging.getLogger(__name__)


def run_pipeline(
    adata: AnnData,
    config: FilterPipelineConfig,
) -> tuple[AnnData, list[dict[str, Any]]]:
    """Execute complete filtering pipeline on an AnnData object.

    The pipeline applies filters in the following order:
    1. UMI count filter
    2. Guide quality filter
    3. Mitochondrial content filter
    4. Gene detection filter

    Only enabled filters are applied.

    Args:
        adata: Input AnnData object
        config: Pipeline configuration

    Returns:
        Tuple of (filtered AnnData, list of statistics from each step)

    """
    logger.info(f"Starting pipeline: {config.pipeline_name}")
    logger.info(f"Initial cell count: {adata.n_obs:,}")

    # Get enabled filters in order
    enabled_filters = config.get_enabled_filters()

    if not enabled_filters:
        logger.warning("No filters enabled in pipeline configuration")
        return adata.copy(), []

    # Track statistics from each step
    all_stats: list[dict[str, Any]] = []

    # Apply filters sequentially
    adata_current = adata
    for filter_name, filter_config in enabled_filters:
        logger.info(f"Applying {filter_name}...")

        # Apply the filter
        adata_current, stats = filter_config.apply(adata_current)

        # Store statistics
        all_stats.append(stats)

        # Log progress
        logger.info(
            f"  → {stats['n_cells_after']:,} cells remaining "
            f"({stats['pct_retained']:.1f}% retained)",
        )

    # Add provenance information
    adata_final = add_provenance(adata_current, config, all_stats)

    logger.info(
        f"Pipeline complete: {adata_final.n_obs:,} / {adata.n_obs:,} cells retained"
    )

    return adata_final, all_stats


def add_provenance(
    adata: AnnData,
    config: FilterPipelineConfig,
    stats: list[dict[str, Any]],
) -> AnnData:
    """Add filtering provenance information to AnnData.uns.

    Stores complete information about the filtering pipeline, including
    configuration, statistics, and timestamp.

    Args:
        adata: AnnData object to annotate
        config: Pipeline configuration that was applied
        stats: List of statistics from each filter step

    Returns:
        AnnData object with provenance added to .uns

    """
    provenance: dict[str, Any] = {
        "pipeline_name": config.pipeline_name,
        "pipeline_hash": config.get_hash(),
        "timestamp": datetime.now().isoformat(),
        "config": config.model_dump(),
        "filter_stats": stats,
        "initial_cells": stats[0]["n_cells_before"] if stats else adata.n_obs,
        "final_cells": adata.n_obs,
        "total_removed": (stats[0]["n_cells_before"] - adata.n_obs) if stats else 0,
    }

    # Store in uns
    if "filtering_provenance" not in adata.uns:
        adata.uns["filtering_provenance"] = []

    # Add this pipeline run to the provenance list
    if isinstance(adata.uns["filtering_provenance"], list):
        adata.uns["filtering_provenance"].append(provenance)
    else:
        # Handle case where uns["filtering_provenance"] exists but isn't a list
        adata.uns["filtering_provenance"] = [provenance]

    logger.info("Added provenance information to adata.uns['filtering_provenance']")

    return adata


def get_pipeline_summary(adata: AnnData) -> dict[str, Any] | None:
    """Extract filtering pipeline provenance from an AnnData object.

    Args:
        adata: AnnData object with filtering provenance

    Returns:
        Most recent provenance dictionary, or None if no provenance exists

    """
    if "filtering_provenance" not in adata.uns:
        return None

    provenance = adata.uns["filtering_provenance"]

    if isinstance(provenance, list) and len(provenance) > 0:
        return provenance[-1]  # Return most recent
    if isinstance(provenance, dict):
        return provenance

    return None


def validate_pipeline_config(config: FilterPipelineConfig, adata: AnnData) -> list[str]:
    """Validate that pipeline configuration is compatible with AnnData object.

    Checks that all required columns referenced in the config exist in the data.

    Args:
        config: Pipeline configuration to validate
        adata: AnnData object to validate against

    Returns:
        List of validation error messages (empty if valid)

    """
    errors: list[str] = []
    obs_columns = set(adata.obs.columns)

    # Check each enabled filter
    enabled_filters = config.get_enabled_filters()

    for filter_name, filter_config in enabled_filters:
        # Get required column based on filter type
        if hasattr(filter_config, "count_column"):
            col = filter_config.count_column
            if col not in obs_columns:
                errors.append(f"{filter_name}: Missing column '{col}' in .obs")

        if hasattr(filter_config, "guide_column"):
            col = filter_config.guide_column
            if col not in obs_columns:
                errors.append(f"{filter_name}: Missing column '{col}' in .obs")

        if hasattr(filter_config, "mt_column"):
            col = filter_config.mt_column
            if col not in obs_columns:
                errors.append(f"{filter_name}: Missing column '{col}' in .obs")

        if hasattr(filter_config, "gene_column"):
            col = filter_config.gene_column
            if col not in obs_columns:
                errors.append(f"{filter_name}: Missing column '{col}' in .obs")

    return errors
