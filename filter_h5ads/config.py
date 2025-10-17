"""
Configuration models for the filtering pipeline.

All filter configurations are immutable Pydantic models that define parameters
for each filtering step. Each config has an apply() method that executes the filter.
"""

import hashlib
import json
from pathlib import Path
from typing import Any

from anndata import AnnData
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


class UMIFilterConfig(FilterStepConfig):
    """Configuration for UMI count filtering.

    Filters cells based on total UMI counts to remove low-quality cells.
    """

    min_counts: int = Field(
        default=15000, ge=0, description="Minimum total UMI counts required per cell"
    )
    count_column: str = Field(
        default="total_counts", description="Name of the column in .obs containing UMI counts"
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply UMI count filter."""
        from filter_h5ads.filters import apply_umi_filter

        return apply_umi_filter(adata, self)


class GuideFilterConfig(FilterStepConfig):
    """Configuration for guide quality filtering.

    Filters cells based on whether they passed guide pairing QC.
    """

    guide_column: str = Field(
        default="pass_guide_filter",
        description="Name of the column in .obs containing guide filter status",
    )
    required_value: int | bool = Field(
        default=1,
        description="Required value in guide_column for cells to pass (typically 1 or True)",
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply guide quality filter."""
        from filter_h5ads.filters import apply_guide_filter

        return apply_guide_filter(adata, self)


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
        default="pct_counts_mt", description="Name of the column in .obs containing MT percentage"
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply mitochondrial content filter."""
        from filter_h5ads.filters import apply_mito_filter

        return apply_mito_filter(adata, self)


class GeneDetectionFilterConfig(FilterStepConfig):
    """Configuration for gene detection filtering.

    Filters cells with too few detected genes (low-quality cells).
    """

    min_genes: int = Field(
        default=200, ge=0, description="Minimum number of genes with non-zero counts required"
    )
    gene_column: str = Field(
        default="n_genes_by_counts", description="Name of the column in .obs containing gene counts"
    )

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply gene detection filter."""
        from filter_h5ads.filters import apply_gene_filter

        return apply_gene_filter(adata, self)


class EnsemblConversionConfig(FilterStepConfig):
    """Configuration for Ensembl ID conversion.

    Converts gene names to Ensembl IDs for both target genes and variable names.
    """

    gene_remapping_ref: Any = Field(
        description="SafeGeneSymbolRemappingRef for gene name to Ensembl ID mapping"
    )
    mask_getter: Any = Field(
        description="Callable that returns list/set of allowed Ensembl IDs"
    )
    source_target_column: str = Field(
        default="gene_target",
        description="Column in .obs containing gene targets to convert",
    )
    target_column_name: str = Field(
        default="target_gene",
        description="Name for new column with converted gene targets",
    )

    model_config = ConfigDict(frozen=True, extra="forbid", arbitrary_types_allowed=True)

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply Ensembl ID conversion."""
        from filter_h5ads.conversions import apply_ensembl_conversion

        # Open the gene remapping reference
        gene_remapping = self.gene_remapping_ref.open_target()

        # Get the mask
        mask = self.mask_getter()

        return apply_ensembl_conversion(
            adata=adata,
            gene_remapping=gene_remapping.gene_symbol_to_ensembl_id,
            mask=mask,
            source_target_column=self.source_target_column,
            target_column_name=self.target_column_name,
        )


class FilterPipelineConfig(BaseModel):
    """Complete pipeline configuration.

    Defines the complete filtering pipeline by composing individual filter configs.
    The pipeline is executed in the order: Ensembl -> UMI -> Guide -> Mito -> Gene.
    """

    model_config = ConfigDict(frozen=True, extra="forbid", arbitrary_types_allowed=True)

    pipeline_name: str = Field(description="Name for this pipeline configuration")
    ensembl_conversion: EnsemblConversionConfig | None = Field(
        default=None, description="Ensembl ID conversion configuration"
    )
    umi_filter: UMIFilterConfig | None = Field(
        default=None, description="UMI count filter configuration"
    )
    guide_filter: GuideFilterConfig | None = Field(
        default=None, description="Guide quality filter configuration"
    )
    mito_filter: MitochondrialFilterConfig | None = Field(
        default=None, description="Mitochondrial content filter configuration"
    )
    gene_filter: GeneDetectionFilterConfig | None = Field(
        default=None, description="Gene detection filter configuration"
    )

    def get_hash(self) -> str:
        """Generate deterministic 8-character hash of this configuration.

        Returns:
            8-character hexadecimal hash string
        """
        # Convert to dict, then to sorted JSON for deterministic hashing
        config_dict = self.model_dump()
        config_json = json.dumps(config_dict, sort_keys=True, indent=None)
        hash_obj = hashlib.sha256(config_json.encode("utf-8"))
        return hash_obj.hexdigest()[:8]

    def save_json(self, path: Path) -> None:
        """Save configuration to JSON file.

        Args:
            path: Path where to save the JSON file
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with path.open("w", encoding="utf-8") as f:
            f.write(self.model_dump_json(indent=2))

    @classmethod
    def load_json(cls, path: Path) -> "FilterPipelineConfig":
        """Load configuration from JSON file.

        Args:
            path: Path to the JSON configuration file

        Returns:
            FilterPipelineConfig instance
        """
        path = Path(path)
        with path.open("r", encoding="utf-8") as f:
            config_dict = json.load(f)
        return cls.model_validate(config_dict)

    def get_enabled_filters(self) -> list[tuple[str, FilterStepConfig]]:
        """Get list of enabled filter configurations in execution order.

        Returns:
            List of (filter_name, filter_config) tuples for enabled filters
        """
        filters: list[tuple[str, FilterStepConfig]] = []

        # Ensembl conversion should happen first (preprocessing)
        if self.ensembl_conversion is not None and self.ensembl_conversion.enabled:
            filters.append(("ensembl_conversion", self.ensembl_conversion))

        if self.umi_filter is not None and self.umi_filter.enabled:
            filters.append(("umi_filter", self.umi_filter))
        if self.guide_filter is not None and self.guide_filter.enabled:
            filters.append(("guide_filter", self.guide_filter))
        if self.mito_filter is not None and self.mito_filter.enabled:
            filters.append(("mito_filter", self.mito_filter))
        if self.gene_filter is not None and self.gene_filter.enabled:
            filters.append(("gene_filter", self.gene_filter))

        return filters
