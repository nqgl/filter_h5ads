"""Base classes and utilities for filter modules.

Contains the FilterStepConfig base class that all filter configurations inherit from,
and shared utility functions used across filter implementations.
"""

from typing import Any

from anndata import AnnData
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field


class FilterStepConfig(BaseModel):
    """Base class for all filter step configurations.

    All filter steps should inherit from this class and implement the apply() method.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    enabled: bool = Field(default=True, description="Whether this filter is enabled")

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply filter and return (filtered_adata, stats_dict).

        Args:
            adata: Input AnnData object

        Returns:
            Tuple of (filtered AnnData, statistics dictionary)

        """
        raise NotImplementedError("Subclasses must implement apply()")


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
