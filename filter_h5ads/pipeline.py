"""Pipeline orchestration for running multiple filters in sequence.

Provides functions to execute a complete filtering pipeline and track
provenance information.
"""

import time
from datetime import datetime
from typing import Any

from anndata import AnnData
from loguru import logger

from filter_h5ads.config import FilterPipelineConfig


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
    pipeline_start = time.time()

    logger.info(f"🚀 Starting pipeline: {config.pipeline_name}")
    logger.info(
        f"📊 Initial dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes"
    )

    # Get enabled filters in order
    enabled_filters = config.get_enabled_filters()

    if not enabled_filters:
        logger.warning("⚠️  No filters enabled in pipeline configuration")
        return adata.copy(), []

    logger.info(f"🔧 {len(enabled_filters)} filters enabled")

    # Track statistics from each step
    all_stats: list[dict[str, Any]] = []

    # Apply filters sequentially
    adata_current = adata
    for i, (filter_name, filter_config) in enumerate(enabled_filters, 1):
        step_start = time.time()
        logger.info(f"\n{'=' * 60}")
        logger.info(f"📍 Step {i}/{len(enabled_filters)}: {filter_name}")
        logger.info(f"{'=' * 60}")

        # Apply the filter
        adata_current, stats = filter_config.apply(adata_current)

        # Store statistics
        all_stats.append(stats)

        step_time = time.time() - step_start

        # Log progress with timing
        logger.success(f"✅ {filter_name} complete in {step_time:.2f}s")
        logger.info(
            f"   Cells: {stats['n_cells_before']:,} → {stats['n_cells_after']:,} "
            f"({stats['pct_retained']:.1f}% retained, "
            f"{stats['n_cells_removed']:,} removed)"
        )

    # Add provenance information
    logger.info("\n📝 Adding provenance information...")
    adata_current.uns["filtering_provenance"] = config.model_dump_json()

    pipeline_time = time.time() - pipeline_start

    logger.success(f"\n{'=' * 60}")
    logger.success(
        f"🎉 Pipeline complete in {pipeline_time:.2f}s ({pipeline_time / 60:.1f} min)"
    )
    logger.success(f"{'=' * 60}")
    logger.info(
        f"📊 Final result: {adata_current.n_obs:,} / {adata.n_obs:,} cells retained "
        f"({adata_current.n_obs / adata.n_obs * 100:.1f}%)"
    )

    return adata_current, all_stats


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
