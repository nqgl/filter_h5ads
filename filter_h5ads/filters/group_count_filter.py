"""Group count filter.

Filters cells based on the number of cells present in each unique combination
of `.obs` metadata keys (e.g., enforce ≥50 cells per (cell_line, drug, dosage)).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import ConfigDict, Field

from filter_h5ads.filters.base import FilterStepConfig, _ensure_in_memory
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class GroupCountFilterConfig(FilterStepConfig):
    """Filter cells by group size (count) in `.obs`."""

    groupby: list[str] = Field(
        min_length=1,
        description="`.obs` columns defining each condition/group.",
    )
    min_cells: int = Field(
        ge=1,
        description="Minimum number of cells required per group to keep that group.",
    )
    dropna: bool = Field(
        default=False,
        description="If True, rows with NA in any groupby key are removed.",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        logger.info(f"Applying group count filter: groupby={self.groupby}")
        validate_h5ad_columns(adata, self.groupby, location="obs")
        adata = _ensure_in_memory(adata)

        adata_before = adata

        obs = adata.obs
        assert isinstance(obs, pd.DataFrame)

        grouped = obs.groupby(self.groupby, sort=False, dropna=self.dropna, observed=True)
        group_ids = grouped.ngroup().to_numpy(dtype=np.int64, copy=False)

        if self.dropna:
            valid = group_ids >= 0
        else:
            valid = np.ones_like(group_ids, dtype=bool)

        if valid.any():
            group_ids_valid = group_ids[valid]
            n_groups = int(group_ids_valid.max() + 1)
            sizes = np.bincount(group_ids_valid, minlength=n_groups)
            per_cell_sizes = np.zeros_like(group_ids, dtype=np.int64)
            per_cell_sizes[valid] = sizes[group_ids_valid]
        else:
            n_groups = 0
            sizes = np.array([], dtype=np.int64)
            per_cell_sizes = np.zeros_like(group_ids, dtype=np.int64)

        keep_mask = valid & (per_cell_sizes >= self.min_cells)
        adata_after = adata_before[keep_mask].copy()

        stats = calculate_filter_stats(
            adata_before=adata_before,
            adata_after=adata_after,
            filter_name="Group Count Filter",
        )
        stats["groupby"] = list(self.groupby)
        stats["min_cells"] = int(self.min_cells)
        stats["dropna"] = self.dropna
        stats["n_groups_observed"] = int(n_groups)
        stats["n_groups_kept"] = int((sizes >= self.min_cells).sum()) if n_groups > 0 else 0
        stats["n_cells_excluded_na"] = int((~valid).sum()) if self.dropna else 0

        logger.info(
            f"Group count filter: {stats['n_cells_after']:,} / "
            f"{stats['n_cells_before']:,} cells retained "
            f"({stats['pct_retained']:.1f}%)"
        )

        return adata_after, stats

