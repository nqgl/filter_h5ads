"""
Inspection and overview utilities for AnnData objects.

Provides functions to display comprehensive information about h5ad files
and filtering results.
"""

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from pydantic import BaseModel, Field

from filter_h5ads.utils import format_number, get_numeric_summary, get_unique_examples


class ColumnInfo(BaseModel):
    """Information about a single column in obs or var."""

    dtype: str
    n_unique: int
    examples: list[Any]
    stats: dict[str, float] | None = None


class H5adSummary(BaseModel):
    """Structured summary of an AnnData object."""

    n_obs: int = Field(description="Number of observations (cells)")
    n_vars: int = Field(description="Number of variables (genes)")
    obs_columns: dict[str, ColumnInfo] = Field(description="Cell metadata columns")
    var_columns: dict[str, ColumnInfo] = Field(description="Gene metadata columns")
    layers: list[str] = Field(default_factory=list, description="Layer names")
    obsm_keys: list[str] = Field(default_factory=list, description="Cell embedding keys")
    varm_keys: list[str] = Field(default_factory=list, description="Gene embedding keys")
    obsp_keys: list[str] = Field(default_factory=list, description="Cell pairwise data keys")
    varp_keys: list[str] = Field(default_factory=list, description="Gene pairwise data keys")
    uns_keys: list[str] = Field(default_factory=list, description="Unstructured annotation keys")
    obsm_shapes: dict[str, tuple[int, ...]] = Field(
        default_factory=dict, description="Shapes of obsm arrays"
    )
    varm_shapes: dict[str, tuple[int, ...]] = Field(
        default_factory=dict, description="Shapes of varm arrays"
    )


def get_h5ad_summary(adata: AnnData) -> H5adSummary:
    """Get structured summary information about an AnnData object.

    Args:
        adata: AnnData object to summarize

    Returns:
        H5adSummary containing structured information
    """
    obs_columns = {}
    for col in adata.obs.columns:
        stats = None
        if pd.api.types.is_numeric_dtype(adata.obs[col]):
            stats = get_numeric_summary(adata.obs[col])

        obs_columns[col] = ColumnInfo(
            dtype=str(adata.obs[col].dtype),
            n_unique=adata.obs[col].nunique(),
            examples=get_unique_examples(adata.obs[col], n=5),
            stats=stats,
        )

    var_columns = {}
    for col in adata.var.columns:
        stats = None
        if pd.api.types.is_numeric_dtype(adata.var[col]):
            stats = get_numeric_summary(adata.var[col])

        var_columns[col] = ColumnInfo(
            dtype=str(adata.var[col].dtype),
            n_unique=adata.var[col].nunique(),
            examples=get_unique_examples(adata.var[col], n=5),
            stats=stats,
        )

    return H5adSummary(
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        obs_columns=obs_columns,
        var_columns=var_columns,
        layers=list(adata.layers.keys()) if adata.layers else [],
        obsm_keys=list(adata.obsm.keys()) if adata.obsm else [],
        varm_keys=list(adata.varm.keys()) if adata.varm else [],
        obsp_keys=list(adata.obsp.keys()) if adata.obsp else [],
        varp_keys=list(adata.varp.keys()) if adata.varp else [],
        uns_keys=list(adata.uns.keys()) if adata.uns else [],
        obsm_shapes={key: adata.obsm[key].shape for key in adata.obsm.keys()} if adata.obsm else {},
        varm_shapes={key: adata.varm[key].shape for key in adata.varm.keys()} if adata.varm else {},
    )


def print_h5ad_overview(adata: AnnData, n_examples: int = 5) -> None:
    """Print comprehensive overview of an AnnData object.

    Args:
        adata: AnnData object to inspect
        n_examples: Number of unique examples to show for each column
    """
    print("=" * 80)
    print("H5AD FILE OVERVIEW")
    print("=" * 80)

    # Dimensions
    print(f"\nDimensions: {format_number(adata.n_obs)} cells × {format_number(adata.n_vars)} genes")

    # Observations (cells) metadata
    print(f"\n{'─' * 80}")
    print("CELL METADATA (.obs)")
    print(f"{'─' * 80}")
    print(f"Number of columns: {len(adata.obs.columns)}\n")

    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        n_unique = adata.obs[col].nunique()

        print(f"  • {col}")
        print(f"    Type: {dtype} | Unique values: {format_number(n_unique)}")

        # Show examples
        examples = get_unique_examples(adata.obs[col], n=n_examples)
        if pd.api.types.is_numeric_dtype(adata.obs[col]):
            # For numeric columns, show statistics
            stats = get_numeric_summary(adata.obs[col])
            print(f"    Range: [{stats['min']:.2f}, {stats['max']:.2f}]")
            mean, median, std = stats["mean"], stats["median"], stats["std"]
            print(f"    Mean: {mean:.2f} | Median: {median:.2f} | Std: {std:.2f}")
        else:
            # For categorical, show examples
            examples_str = ", ".join(str(x) for x in examples)
            if n_unique > n_examples:
                examples_str += ", ..."
            print(f"    Examples: {examples_str}")

        print()

    # Variables (genes) metadata
    print(f"{'─' * 80}")
    print("GENE METADATA (.var)")
    print(f"{'─' * 80}")
    print(f"Number of columns: {len(adata.var.columns)}\n")

    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        n_unique = adata.var[col].nunique()

        print(f"  • {col}")
        print(f"    Type: {dtype} | Unique values: {format_number(n_unique)}")

        # Show examples
        examples = get_unique_examples(adata.var[col], n=n_examples)
        if pd.api.types.is_numeric_dtype(adata.var[col]):
            # For numeric columns, show statistics
            stats = get_numeric_summary(adata.var[col])
            print(f"    Range: [{stats['min']:.2f}, {stats['max']:.2f}]")
            mean, median, std = stats["mean"], stats["median"], stats["std"]
            print(f"    Mean: {mean:.2f} | Median: {median:.2f} | Std: {std:.2f}")
        else:
            # For categorical, show examples
            examples_str = ", ".join(str(x) for x in examples)
            if n_unique > n_examples:
                examples_str += ", ..."
            print(f"    Examples: {examples_str}")

        print()

    # Additional data structures
    if adata.layers:
        print(f"{'─' * 80}")
        print("LAYERS")
        print(f"{'─' * 80}")
        for key in adata.layers.keys():
            print(f"  • {key}: {adata.layers[key].shape}")
        print()

    if adata.obsm:
        print(f"{'─' * 80}")
        print("CELL EMBEDDINGS (.obsm)")
        print(f"{'─' * 80}")
        for key in adata.obsm.keys():
            print(f"  • {key}: {adata.obsm[key].shape}")
        print()

    if adata.varm:
        print(f"{'─' * 80}")
        print("GENE EMBEDDINGS (.varm)")
        print(f"{'─' * 80}")
        for key in adata.varm.keys():
            print(f"  • {key}: {adata.varm[key].shape}")
        print()

    if adata.obsp:
        print(f"{'─' * 80}")
        print("CELL PAIRWISE DATA (.obsp)")
        print(f"{'─' * 80}")
        for key in adata.obsp.keys():
            print(f"  • {key}: {adata.obsp[key].shape}")
        print()

    if adata.varp:
        print(f"{'─' * 80}")
        print("GENE PAIRWISE DATA (.varp)")
        print(f"{'─' * 80}")
        for key in adata.varp.keys():
            print(f"  • {key}: {adata.varp[key].shape}")
        print()

    if adata.uns:
        print(f"{'─' * 80}")
        print("UNSTRUCTURED ANNOTATIONS (.uns)")
        print(f"{'─' * 80}")
        for key in adata.uns.keys():
            value = adata.uns[key]
            if isinstance(value, dict):
                print(f"  • {key}: dict with {len(value)} keys")
            elif isinstance(value, (list, tuple)):
                print(f"  • {key}: {type(value).__name__} with {len(value)} items")
            elif isinstance(value, (np.ndarray, pd.DataFrame)):
                print(f"  • {key}: {type(value).__name__} with shape {value.shape}")
            else:
                print(f"  • {key}: {type(value).__name__}")
        print()

    print("=" * 80)


def print_filtering_summary(stats_list: list[dict[str, Any]]) -> None:
    """Print summary of filtering pipeline results.

    Args:
        stats_list: List of statistics dictionaries from each filter step
    """
    print("\n" + "=" * 80)
    print("FILTERING PIPELINE SUMMARY")
    print("=" * 80)

    if not stats_list:
        print("No filters were applied.")
        return

    # Overall summary
    initial_cells = stats_list[0]["n_cells_before"]
    final_cells = stats_list[-1]["n_cells_after"]
    total_removed = initial_cells - final_cells
    total_pct_retained = (final_cells / initial_cells * 100) if initial_cells > 0 else 0

    print(f"\nInitial cells: {format_number(initial_cells)}")
    print(f"Final cells:   {format_number(final_cells)}")
    print(f"Removed:       {format_number(total_removed)} ({100 - total_pct_retained:.1f}%)")
    print(f"Retained:      {total_pct_retained:.1f}%")

    # Individual filter steps
    print(f"\n{'─' * 80}")
    print("FILTER STEPS")
    print(f"{'─' * 80}\n")

    for i, stats in enumerate(stats_list, 1):
        print(f"Step {i}: {stats['filter_name']}")
        print(f"  Before: {format_number(stats['n_cells_before'])} cells")
        print(f"  After:  {format_number(stats['n_cells_after'])} cells")
        removed = format_number(stats["n_cells_removed"])
        pct_removed = 100 - stats["pct_retained"]
        print(f"  Removed: {removed} cells ({pct_removed:.1f}%)")
        print(f"  Retained: {stats['pct_retained']:.1f}%")

        # Show breakdown by gene target if available
        if "breakdown_after" in stats and stats["breakdown_after"]:
            print("\n  Breakdown by gene_target:")
            breakdown = stats["breakdown_after"]

            # Sort by count (descending)
            sorted_breakdown = sorted(breakdown.items(), key=lambda x: x[1], reverse=True)

            # Show top targets
            for target, count in sorted_breakdown[:10]:
                print(f"    • {target}: {format_number(count)} cells")

            if len(sorted_breakdown) > 10:
                remaining = len(sorted_breakdown) - 10
                print(f"    ... and {remaining} more targets")

        print()

    print("=" * 80)
