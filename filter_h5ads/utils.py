"""Utility functions for validation, logging, and helper operations."""

import logging
from typing import Any

import pandas as pd
from anndata import AnnData


def setup_logging(level: str = "INFO") -> None:
    """Configure logging for the package.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    """
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def validate_h5ad_columns(
    adata: AnnData, required_cols: list[str], location: str = "obs"
) -> None:
    """Validate that required columns exist in AnnData object.

    Args:
        adata: AnnData object to validate
        required_cols: List of required column names
        location: Where to look for columns ('obs' or 'var')

    Raises:
        ValueError: If any required columns are missing

    """
    if location == "obs":
        available_cols = set(adata.obs.columns)
        data_type = "cell metadata (.obs)"
    elif location == "var":
        available_cols = set(adata.var.columns)
        data_type = "gene metadata (.var)"
    else:
        raise ValueError(f"Invalid location: {location}. Must be 'obs' or 'var'")

    missing_cols = [col for col in required_cols if col not in available_cols]

    if missing_cols:
        raise ValueError(
            f"Missing required columns in {data_type}: {missing_cols}\n"
            f"Available columns: {sorted(available_cols)}",
        )


def get_column_breakdown(
    adata: AnnData, groupby_col: str = "gene_target"
) -> dict[str, int]:
    """Get cell count breakdown by a grouping column.

    Args:
        adata: AnnData object
        groupby_col: Column name to group by

    Returns:
        Dictionary mapping group values to cell counts

    """
    if groupby_col not in adata.obs.columns:
        return {}

    col_series = adata.obs[groupby_col]
    assert isinstance(col_series, pd.Series)
    value_counts = col_series.value_counts()
    return value_counts.to_dict()


def calculate_filter_stats(
    adata_before: AnnData,
    adata_after: AnnData,
    filter_name: str,
    groupby_col: str = "gene_target",
) -> dict[str, Any]:
    """Calculate comprehensive statistics for a filtering step.

    Args:
        adata_before: AnnData before filtering
        adata_after: AnnData after filtering
        filter_name: Name of the filter
        groupby_col: Column to break down statistics by

    Returns:
        Dictionary with filtering statistics

    """
    n_before = adata_before.n_obs
    n_after = adata_after.n_obs
    n_removed = n_before - n_after
    pct_retained = (n_after / n_before * 100) if n_before > 0 else 0.0

    stats: dict[str, Any] = {
        "filter_name": filter_name,
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "n_cells_removed": n_removed,
        "pct_retained": round(pct_retained, 2),
    }

    # Add breakdown by grouping column if available
    if groupby_col in adata_before.obs.columns:
        breakdown_before = get_column_breakdown(adata_before, groupby_col)
        breakdown_after = get_column_breakdown(adata_after, groupby_col)

        breakdown_removed = {
            key: breakdown_before.get(key, 0) - breakdown_after.get(key, 0)
            for key in breakdown_before.keys()
        }

        stats["breakdown_before"] = breakdown_before
        stats["breakdown_after"] = breakdown_after
        stats["breakdown_removed"] = breakdown_removed

    return stats


def safe_copy(adata: AnnData) -> AnnData:
    """Create a safe copy of AnnData object.

    Args:
        adata: AnnData object to copy

    Returns:
        Copy of the AnnData object

    """
    return adata.copy()


def get_numeric_summary(series: pd.Series) -> dict[str, float]:
    """Get summary statistics for a numeric series.

    Args:
        series: Pandas series with numeric data

    Returns:
        Dictionary with min, max, mean, median, std

    """
    return {
        "min": float(series.min()),
        "max": float(series.max()),
        "mean": float(series.mean()),
        "median": float(series.median()),
        "std": float(series.std()),
    }


def format_number(num: float) -> str:
    """Format large numbers with thousand separators.

    Args:
        num: Number to format

    Returns:
        Formatted string

    """
    if isinstance(num, float):
        return f"{num:,.2f}"
    return f"{num:,}"


def get_unique_examples(series: pd.Series, n: int = 5) -> list[Any]:
    """Get up to n unique examples from a series.

    Args:
        series: Pandas series
        n: Maximum number of examples

    Returns:
        List of unique values

    """
    unique_vals = series.unique()
    if len(unique_vals) > n:
        return list(unique_vals[:n])
    return list(unique_vals)
