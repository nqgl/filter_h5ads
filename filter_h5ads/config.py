"""Configuration models for the filtering pipeline.

All filter configurations are immutable Pydantic models that define parameters
for each filtering step. Each config has an apply() method that executes the filter.
"""

import hashlib
import json
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from filter_h5ads.filters import (
    EnsemblConversionConfig,
    FilterStepConfig,
    GeneDetectionFilterConfig,
    GuideFilterConfig,
    MitochondrialFilterConfig,
    ObsColumnTransformConfig,
    ObsValueFilterConfig,
    PseudobulkConfig,
    UMIFilterConfig,
)


class FilterPipelineConfig(BaseModel):
    """Complete pipeline configuration.

    Defines the complete filtering pipeline by composing individual filter configs.
    The pipeline is executed in the order:
    Ensembl -> Obs Column Transform -> Obs Value -> UMI -> Guide -> Mito -> Gene -> Pseudobulk.
    """

    model_config = ConfigDict(frozen=True, extra="forbid", arbitrary_types_allowed=True)

    pipeline_name: str = Field(description="Name for this pipeline configuration")
    ensembl_conversion: EnsemblConversionConfig | None = Field(
        default=None,
        description="Ensembl ID conversion configuration",
    )
    obs_column_transform: ObsColumnTransformConfig | None = Field(
        default=None,
        description="Derive a new .obs column by transforming existing .obs values",
    )
    obs_value_filter: ObsValueFilterConfig | None = Field(
        default=None,
        description="Filter cells by an .obs column value (include/exclude values)",
    )
    umi_filter: UMIFilterConfig | None = Field(
        default=None,
        description="UMI count filter configuration",
    )
    guide_filter: GuideFilterConfig | None = Field(
        default=None,
        description="Guide quality filter configuration",
    )
    mito_filter: MitochondrialFilterConfig | None = Field(
        default=None,
        description="Mitochondrial content filter configuration",
    )
    gene_filter: GeneDetectionFilterConfig | None = Field(
        default=None,
        description="Gene detection filter configuration",
    )
    pseudobulk: PseudobulkConfig | None = Field(
        default=None,
        description="Group cells by .obs keys and sum counts into pooled samples",
    )

    def get_hash(self) -> str:
        """Generate deterministic 8-character hash of this configuration.

        Returns:
            8-character hexadecimal hash string

        """
        # Convert to dict, then to sorted JSON for deterministic hashing
        config_json = self.model_dump_json()
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

        if self.obs_column_transform is not None and self.obs_column_transform.enabled:
            filters.append(("obs_column_transform", self.obs_column_transform))

        if self.obs_value_filter is not None and self.obs_value_filter.enabled:
            filters.append(("obs_value_filter", self.obs_value_filter))

        if self.umi_filter is not None and self.umi_filter.enabled:
            filters.append(("umi_filter", self.umi_filter))
        if self.guide_filter is not None and self.guide_filter.enabled:
            filters.append(("guide_filter", self.guide_filter))
        if self.mito_filter is not None and self.mito_filter.enabled:
            filters.append(("mito_filter", self.mito_filter))
        if self.gene_filter is not None and self.gene_filter.enabled:
            filters.append(("gene_filter", self.gene_filter))

        if self.pseudobulk is not None and self.pseudobulk.enabled:
            filters.append(("pseudobulk", self.pseudobulk))

        return filters
