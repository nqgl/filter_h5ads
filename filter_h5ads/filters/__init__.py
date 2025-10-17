"""Filters subpackage for h5ad filtering operations.

Each filter is implemented as a config class with an apply() method.
"""

from filter_h5ads.filters.base import FilterStepConfig as FilterStepConfig
from filter_h5ads.filters.ensembl_conversion import (
    EnsemblConversionConfig as EnsemblConversionConfig,
)
from filter_h5ads.filters.gene_detection import (
    GeneDetectionFilterConfig as GeneDetectionFilterConfig,
)
from filter_h5ads.filters.guide_filter import GuideFilterConfig as GuideFilterConfig
from filter_h5ads.filters.mitochondrial import (
    MitochondrialFilterConfig as MitochondrialFilterConfig,
)
from filter_h5ads.filters.umi_filter import UMIFilterConfig as UMIFilterConfig
