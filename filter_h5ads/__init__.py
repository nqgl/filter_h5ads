"""filter_h5ads: A robust pipeline for filtering h5ad files from CRISPR screens.

This package provides a configurable, type-safe filtering pipeline for AnnData objects,
with particular focus on CRISPR perturbation screen quality control.
"""

from filter_h5ads.config import (
    EnsemblConversionConfig as EnsemblConversionConfig,
)
from filter_h5ads.config import (
    FilterPipelineConfig as FilterPipelineConfig,
)
from filter_h5ads.config import (
    FilterStepConfig as FilterStepConfig,
)
from filter_h5ads.config import (
    GeneDetectionFilterConfig as GeneDetectionFilterConfig,
)
from filter_h5ads.config import (
    GroupCountAnnotateConfig as GroupCountAnnotateConfig,
)
from filter_h5ads.config import (
    GroupCountFilterConfig as GroupCountFilterConfig,
)
from filter_h5ads.config import (
    GuideFilterConfig as GuideFilterConfig,
)
from filter_h5ads.config import (
    MitochondrialFilterConfig as MitochondrialFilterConfig,
)
from filter_h5ads.config import (
    ObsColumnTransformConfig as ObsColumnTransformConfig,
)
from filter_h5ads.config import (
    PseudobulkConfig as PseudobulkConfig,
)
from filter_h5ads.config import (
    ObsValueFilterConfig as ObsValueFilterConfig,
)
from filter_h5ads.config import (
    UMIFilterConfig as UMIFilterConfig,
)
from filter_h5ads.filters.ensembl_conversion import (
    convert_gene_labels as convert_gene_labels,
)
from filter_h5ads.inspection import (
    ColumnInfo as ColumnInfo,
)
from filter_h5ads.inspection import (
    H5adSummary as H5adSummary,
)
from filter_h5ads.inspection import (
    get_h5ad_summary as get_h5ad_summary,
)
from filter_h5ads.inspection import (
    print_filtering_summary as print_filtering_summary,
)
from filter_h5ads.inspection import (
    print_h5ad_overview as print_h5ad_overview,
)
from filter_h5ads.io import (
    batch_filter_h5ads as batch_filter_h5ads,
)
from filter_h5ads.io import (
    load_h5ad as load_h5ad,
)
from filter_h5ads.io import (
    save_filtered_h5ad as save_filtered_h5ad,
)
from filter_h5ads.pipeline import (
    run_pipeline as run_pipeline,
)

__version__ = "0.1.0"
